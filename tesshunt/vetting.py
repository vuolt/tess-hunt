"""Automated vetting of single-event dips from a full-sector search.

Labels each dip with:
  common_mode     many other stars within CM_RADIUS_DEG dip at the same time
                  (local Poisson test), i.e. a detector systematic
  partial         the box reaches within 1 h of a data gap
  category        common_mode | repeating | secondary | single | single_partial
  rp_rjup         companion radius implied by the depth and the TIC stellar radius
  tier            A: planet-sized, not partial, outside sector-wide pile-up
                  windows; B: other candidates

Phase 3 adds per-candidate checks (see the section further down):
  1 shape_test        limb-darkened transit vs box / ramp / step / flare-decay (BIC)
  2 duration_check    T14 long enough for P > DUR_P_MIN (unless grazing)
  3 edge_check        no dips against data gaps unless resolved on both sides
  4 centroid_test, neighbour_test   difference imaging on TESScut pixels
  5 asteroid_check    SkyBoT known solar-system objects at the dip time
  6 catalogue_check   ExoFOP TOIs/CTOIs, TESS EB catalogue
  7 fpp_check         TRICERATOPS false-positive probability
"""

from __future__ import annotations

import os
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.stats import poisson

CM_BIN_D = 0.25        # time bin for the sector-wide histogram (figure only)
CM_RADIUS_DEG = 2.0    # neighbourhood for the local common-mode test
CM_DT = 0.25           # |dt| (d) for two dips to count as coincident
CM_PVALUE = 1e-3       # Poisson tail probability for a dip to count as common-mode
REPEAT_FACTOR = 2.0    # dips within this factor in duration and depth "repeat"
RSUN_RJUP = 9.731      # solar radius in Jupiter radii
RP_MAX_RJUP = 2.0      # largest plausible planet (inflated hot Jupiters reach ~2 RJ)
BTJD = 2457000.0


def common_mode_bins(t0, bin_d=CM_BIN_D):
    """Time bins holding far more dips (from different stars) than typical.

    Returns (edges, counts, flagged_bin_mask, expected_rate)."""
    edges = np.arange(np.floor(t0.min()), np.ceil(t0.max()) + bin_d, bin_d)
    counts, _ = np.histogram(t0, edges)
    lam = float(np.median(counts[counts > 0]))
    flagged = poisson.sf(counts - 1, lam) < 1e-4
    return edges, counts, flagged, lam


def local_common_mode(d, radius_deg=CM_RADIUS_DEG, dt=CM_DT, pvalue=CM_PVALUE):
    """Flag dips that coincide in time with unusually many dips on nearby stars.

    Systematics (scattered light, thermal settling, pointing jitter) hit stars
    on the same part of a detector at the same time; a real transit does not.
    For each dip, count dips on *other* stars within ``radius_deg`` and
    ``|dt|`` and compare with the count expected from that patch's own dip
    rate over the whole sector (Poisson). Returns (flag, n_coincident, expected).
    """
    from scipy.spatial import cKDTree
    ra, dec = np.radians(d.ra.values), np.radians(d.dec.values)
    xyz = np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]
    tree = cKDTree(xyz)
    chord = 2 * np.sin(np.radians(radius_deg) / 2)
    span = d.t0.max() - d.t0.min()
    t0, tic = d.t0.values, d.tic.values
    n_co = np.zeros(len(d), int)
    lam = np.zeros(len(d))
    for i, nb in enumerate(tree.query_ball_point(xyz, chord)):
        nb = np.asarray(nb)
        nb = nb[tic[nb] != tic[i]]
        n_co[i] = int(np.sum(np.abs(t0[nb] - t0[i]) < dt))
        lam[i] = max(len(nb) * 2 * dt / span, 0.5)
    flag = poisson.sf(n_co - 1, lam) < pvalue
    return flag, n_co, lam


def fetch_tois(tics):
    q = ("select tid,toi,tfopwg_disp,pl_orbper,pl_tranmid,pl_trandurh,pl_trandep from toi")
    url = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
           + urllib.parse.quote(q) + "&format=csv")
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            tois = pd.read_csv(r)
    except Exception as e:  # noqa: BLE001
        print(f"TOI table unavailable ({e}); skipping cross-match")
        return pd.DataFrame(columns=["tid", "toi", "tfopwg_disp", "pl_orbper",
                                     "pl_tranmid", "pl_trandurh", "pl_trandep"])
    return tois[tois.tid.isin(set(tics))]


def match_toi(dips, tois):
    """Name of a TOI whose ephemeris predicts a transit at the dip time."""
    out = np.full(len(dips), "", dtype=object)
    host = np.full(len(dips), "", dtype=object)
    by_tic = {k: g for k, g in tois.groupby("tid")}
    for i, (tic, t0) in enumerate(zip(dips.tic.values, dips.t0.values)):
        g = by_tic.get(tic)
        if g is None:
            continue
        host[i] = ",".join(str(x) for x in g.toi)
        for _, r in g.iterrows():
            if not (np.isfinite(r.pl_orbper) and np.isfinite(r.pl_tranmid)) or r.pl_orbper <= 0:
                continue
            tm = r.pl_tranmid - BTJD
            n = np.round((t0 - tm) / r.pl_orbper)
            tol = max((r.pl_trandurh if np.isfinite(r.pl_trandurh) else 3) / 24, 0.1)
            if abs(t0 - (tm + n * r.pl_orbper)) < tol:
                out[i] = str(r.toi)
                break
    return out, host


def classify(d):
    """Per star, look at the dips that are not common-mode.

    A periodic signal (eclipsing binary, short-period planet) repeats with the
    same shape, so the strongest dip is 'repeating' if another dip on the star
    has a duration and depth within a factor REPEAT_FACTOR of it. Otherwise the
    strongest dip is a single-transit candidate ('single', or 'single_partial'
    if it touches a data gap) and any weaker, dissimilar dips on the same star
    are 'secondary'.
    """
    cat = pd.Series("common_mode", index=d.index, dtype=object)
    for _, g in d[~d.common_mode].groupby("tic"):
        g = g.sort_values("snr", ascending=False)
        p = g.iloc[0]
        rest = g.iloc[1:]
        same = ((rest.duration_h / p.duration_h).between(1 / REPEAT_FACTOR, REPEAT_FACTOR)
                & (rest.depth_ppm / p.depth_ppm).between(1 / REPEAT_FACTOR, REPEAT_FACTOR))
        if same.any():
            cat[g.index] = "secondary"
            cat[[g.index[0]] + list(rest.index[same])] = "repeating"
        else:
            cat[g.index] = "secondary"
            cat[g.index[0]] = "single_partial" if p.partial else "single"
    return cat.values


def toi_recall(d, tois, searched_tics, windows):
    """For each TOI on a searched star with a predicted transit inside the
    sector's data windows, was a dip found at that time?"""
    rows = []
    for _, r in tois[tois.tid.isin(searched_tics)].iterrows():
        if not (np.isfinite(r.pl_orbper) and np.isfinite(r.pl_tranmid)) or r.pl_orbper <= 0:
            continue
        tm, per = r.pl_tranmid - BTJD, r.pl_orbper
        lo, hi = windows[0][0], windows[-1][1]
        n = np.arange(np.ceil((lo - tm) / per), np.floor((hi - tm) / per) + 1)
        pred = tm + n * per
        pred = [t for t in pred if any(a + 0.1 <= t <= b - 0.1 for a, b in windows)]
        if not pred:
            continue
        g = d[(d.tic == r.tid) & (d.toi_match == str(r.toi))]
        rows.append(dict(toi=r.toi, tic=r.tid, disp=r.tfopwg_disp, period_d=per,
                         depth_ppm=r.pl_trandep, dur_h=r.pl_trandurh, n_pred=len(pred),
                         recovered=len(g) > 0,
                         category=g.sort_values("snr").category.iloc[-1] if len(g) else "",
                         snr=g.snr.max() if len(g) else np.nan))
    return pd.DataFrame(rows)


def vet(stars, dips, tois, coords):
    d = dips.merge(stars[["tic", "tmag", "teff", "radius", "window_d", "max_dur_h",
                          "long_limited", "hf_variable", "var_period_d", "var_amp_ppm",
                          "sigma_pt_ppm", "crowdsap"]], on="tic", how="left")
    d = d.merge(coords, on="tic", how="left")
    edges, counts, flagged, lam = common_mode_bins(d.t0.values)
    b = np.clip(np.digitize(d.t0.values, edges) - 1, 0, len(counts) - 1)
    d["crowded_time_bin"] = flagged[b]          # sector-wide pile-up (figure/diagnostic)
    d["common_mode"], d["n_coincident"], d["n_coincident_expected"] = local_common_mode(d)
    # Box reaches within 1 h of a data gap: the dip may be a truncated (partial)
    # transit or an edge ramp. Kept, but labelled.
    d["partial"] = d.edge_dist_h < d.duration_h / 2 + 1
    d["n_dips_star"] = d.groupby("tic").t0.transform("size")
    d["toi_match"], d["toi_host"] = match_toi(d, tois)
    d["category"] = classify(d)
    d["candidate"] = np.isin(d.category, ["single", "single_partial"]) & (d.hf_variable == 0)
    # Companion radius implied by the depth and the TIC stellar radius. Anything
    # above RP_MAX_RJUP is a star (eclipsing binary), not a planet. Unknown
    # radius counts as planet-sized so nothing is silently dropped.
    d["rp_rjup"] = np.sqrt(np.clip(d.depth_ppm, 0, None) * 1e-6) * d.radius * RSUN_RJUP
    d["planet_sized"] = d.rp_rjup.isna() | (d.rp_rjup <= RP_MAX_RJUP)
    # Tier A: the cleanest candidates. Tier B: passed vetting but planet-sized
    # fails, the dip touches a gap, or it falls in a sector-wide pile-up window
    # (orbit starts/ends), where most dips are systematics.
    clean = d.planet_sized & ~d.partial & ~d.crowded_time_bin
    d["tier"] = np.where(d.candidate & clean, "A", np.where(d.candidate, "B", ""))
    d = d.sort_values("snr", ascending=False).reset_index(drop=True)
    d.insert(0, "rank", np.arange(1, len(d) + 1))
    cm = dict(edges=edges, counts=counts, flagged=flagged, lam=lam)
    return d, cm


# =====================================================================
# Phase 3: per-candidate vetting checks
# =====================================================================
#
# Each check returns a dict with at least {"pass": bool, ...details}. The
# light-curve checks (shape, duration, edge) are pure functions of arrays so
# they can be unit-tested; the pixel, asteroid, catalogue and FPP checks wrap
# network services and keep their decision logic in separate pure functions.

LD_U = (0.40, 0.25)        # quadratic limb darkening, TESS band, ~solar (fixed)
SHAPE_P_REF = 100.0        # d; transit shape at fixed T14 is insensitive to P
SHAPE_DBIC_ALT = 10.0      # transit must beat ramp/step/flare-decay by this much
SHAPE_DBIC_BOX = -6.0      # ...and a box must not be preferred by more than this
DUR_P_MIN = 20.0           # d; the orbit we want to remain possible
DUR_FRAC = 0.5             # allow T14 down to this fraction of the central duration
GRAZE_B = 0.9              # impact parameter above which the fit counts as grazing
EDGE_GAP_H = 1.0           # a gap longer than this (h) is a data edge
EDGE_H = 3.0               # dips whose ingress/egress is within this of an edge...
EDGE_RESOLVE_FRAC = 0.8    # ...pass only with this baseline coverage on both sides
EDGE_STEP_SIGMA = 3.0      # ...and pre/post baselines agreeing within this (sigma)
EDGE_STEP_DEPTH = 0.25     # ...or within this fraction of the transit depth
EDGE_BASE_H = (1.0, 3.0)   # baseline window each side: T14/2 clipped to this range (h)
G_CGS = 6.674e-8
MSUN_G = 1.989e33
RSUN_CM = 6.957e10


def _robust_sigma_pt(f):
    d = np.diff(f)
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2))


# ------------------------------------------------------------------ models

def transit_model(t, t0, rp, t14, b, u=LD_U, per=SHAPE_P_REF, exp_time=None):
    """Limb-darkened transit (batman), parameterized by T14 and impact parameter.

    For a circular orbit, sin(pi T14 / P) = sqrt((1 + rp)^2 - b^2) / (a sin i),
    with b = a cos i, which gives a/R* and i from (T14, b, rp).
    """
    import batman
    b = min(b, (1 + rp) * 0.999)
    s = np.sin(np.pi * t14 / per)
    a = np.sqrt(b ** 2 + ((1 + rp) ** 2 - b ** 2) / s ** 2)
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a = t0, per, rp, a
    p.inc = float(np.degrees(np.arccos(b / a)))
    p.ecc, p.w, p.u, p.limb_dark = 0.0, 90.0, list(u), "quadratic"
    if exp_time:
        m = batman.TransitModel(p, t, supersample_factor=5, exp_time=exp_time)
    else:
        m = batman.TransitModel(p, t)
    return m.light_curve(p)


def _linfit(cols, f, w):
    """Weighted linear least squares; returns (coef, chi2)."""
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A * w[:, None], f * w, rcond=None)
    r = (f - A @ coef) * w
    return coef, float(r @ r)


@dataclass
class ModelFit:
    name: str
    k: int                      # free parameters
    chi2: float
    params: dict
    model: np.ndarray = field(repr=False, default=None)

    def bic(self, n):
        return self.chi2 + self.k * np.log(n)


def fit_transit(t, f, sigma, t0, dur, exp_time=None, box=None):
    """Least-squares limb-darkened transit fit (t0, rp, T14, b, linear baseline).

    Multi-start over T14 and impact parameter; if ``box`` (a fitted box
    ModelFit) is given, its centre, duration and depth seed extra starts and
    widen the bounds, so a dip that is longer or offset from the detection box
    is not missed by a local minimum."""
    from scipy.optimize import least_squares
    x = t - t0
    depth0 = max(1 - np.median(f[np.abs(x) < dur / 4]) if np.any(np.abs(x) < dur / 4) else 1e-3,
                 1e-4)
    centres, durs = [t0], [dur]
    if box is not None:
        centres.append(box.params["t0"])
        durs.append(box.params["dur_h"] / 24)
        depth0 = max(depth0, box.params["depth_ppm"] * 1e-6)
    span = max(max(durs), 1 / 24)
    lo = [min(centres) - span, 0.003, 0.25 / 24, 0.0, 0.9, -1.0]
    hi = [max(centres) + span, 0.9, max(3 * max(durs), max(durs) + 4 / 24), 0.99, 1.1, 1.0]

    def resid(p):
        m = transit_model(t, p[0], p[1], p[2], p[3] * (1 + p[1]), exp_time=exp_time)
        return (f - m * (p[4] + p[5] * x)) / sigma

    best = None
    for c, d in zip(centres, durs):
        for tf in (0.8, 1.0, 1.3):
            for beta in (0.1, 0.6, 0.9):
                rp0 = np.clip(np.sqrt(depth0) * (1.2 if beta > 0.8 else 1.0), 0.004, 0.85)
                p0 = np.clip([c, rp0, d * tf, beta, 1.0, 0.0], lo, hi)
                try:
                    r = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac", max_nfev=400)
                except Exception:  # noqa: BLE001
                    continue
                if best is None or r.cost < best.cost:
                    best = r
    p = best.x
    b = p[3] * (1 + p[1])
    model = transit_model(t, p[0], p[1], p[2], b, exp_time=exp_time) * (p[4] + p[5] * x)
    return ModelFit("transit", 6, float(2 * best.cost), dict(
        t0=float(p[0]), rp=float(p[1]), t14_h=float(p[2] * 24), b=float(b),
        depth_ppm=float((1 - transit_model(np.array([p[0]]), p[0], p[1], p[2], b)[0]) * 1e6),
        grazing=bool(b > 1 - p[1] or b >= GRAZE_B)), model)


def fit_box(t, f, sigma, t0, dur, step=10 / 1440):
    """Box with linear baseline: coarse grid in (centre, duration), then a fine
    grid (quarter-cadence steps) around the best coarse solution, so the box
    is not handicapped against the continuous transit model."""
    x = t - t0
    w = np.full(len(t), 1 / sigma)
    one = np.ones_like(t)
    cad = float(np.median(np.diff(t)))

    def search(centres, durs, best=None):
        for d in durs:
            for c in centres:
                inb = (np.abs(t - c) < d / 2).astype(float)
                if inb.sum() < 2:
                    continue
                coef, chi2 = _linfit([one, x, -inb], f, w)
                if best is None or chi2 < best[0]:
                    best = (chi2, c, d, coef)
        return best

    coarse = np.array([0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32]) / 24
    durs = [d for d in coarse if 0.3 * dur <= d <= 2.5 * dur] or [dur]
    span = max(dur, 1 / 24)
    best = search(np.arange(t0 - span, t0 + span + 1e-9, step), durs)
    _, c, d, _ = best
    fine = cad / 4
    best = search(np.arange(c - step, c + step + 1e-9, fine),
                  np.arange(max(d * 0.6, 2 * cad), d * 1.4 + 1e-9, fine), best)
    chi2, c, d, coef = best
    inb = (np.abs(t - c) < d / 2).astype(float)
    return ModelFit("box", 5, chi2, dict(t0=float(c), dur_h=float(d * 24),
                                         depth_ppm=float(coef[2] * 1e6)),
                    coef[0] + coef[1] * x - coef[2] * inb)


def fit_ramp(t, f, sigma, t0):
    x = t - t0
    w = np.full(len(t), 1 / sigma)
    coef, chi2 = _linfit([np.ones_like(t), x, x ** 2], f, w)
    return ModelFit("ramp", 3, chi2, dict(c=[float(c) for c in coef]),
                    coef[0] + coef[1] * x + coef[2] * x ** 2)


def _grid_times(t, step=10 / 1440):
    return np.arange(t.min() + step, t.max() - step, step)


def fit_step(t, f, sigma, t0):
    x = t - t0
    w = np.full(len(t), 1 / sigma)
    one = np.ones_like(t)
    best = None
    for ts in _grid_times(t):
        h = (t >= ts).astype(float)
        if h.sum() < 3 or (1 - h).sum() < 3:
            continue
        coef, chi2 = _linfit([one, x, h], f, w)
        if best is None or chi2 < best[0]:
            best = (chi2, ts, coef)
    chi2, ts, coef = best
    return ModelFit("step", 4, chi2, dict(t_step=float(ts), jump_ppm=float(coef[2] * 1e6)),
                    coef[0] + coef[1] * x + coef[2] * (t >= ts))


def fit_flare_decay(t, f, sigma, t0):
    """Sudden jump at t_s followed by exponential relaxation, either sign.

    Negative amplitude: a sudden drop and slow recovery, e.g. a pixel
    sensitivity dropout after a cosmic ray. Positive: a flare."""
    x = t - t0
    w = np.full(len(t), 1 / sigma)
    one = np.ones_like(t)
    best = None
    for tau in np.geomspace(20 / 1440, 1.0, 9):
        for ts in _grid_times(t):
            e = np.where(t >= ts, np.exp(-(t - ts) / tau), 0.0)
            if (e > 0.05).sum() < 2:
                continue
            coef, chi2 = _linfit([one, x, e], f, w)
            if best is None or chi2 < best[0]:
                best = (chi2, ts, tau, coef)
    chi2, ts, tau, coef = best
    e = np.where(t >= ts, np.exp(-(t - ts) / tau), 0.0)
    return ModelFit("flare_decay", 5, chi2, dict(t_start=float(ts), tau_h=float(tau * 24),
                                                 amp_ppm=float(coef[2] * 1e6)),
                    coef[0] + coef[1] * x + coef[2] * e)


def shape_passes(dbic_alt, dbic_box, dbic_box_min=None):
    """Shape decision from the stored BIC differences (see shape_test)."""
    dbic_box_min = SHAPE_DBIC_BOX if dbic_box_min is None else dbic_box_min
    return bool(dbic_alt >= SHAPE_DBIC_ALT and dbic_box >= dbic_box_min)


def shape_window(time, flux, t0, dur):
    """Data used for the shape fits: +/- max(2.5 durations, 0.5 d) around t0."""
    half = max(2.5 * dur, 0.5)
    sel = np.abs(time - t0) < half
    return time[sel], flux[sel]


def shape_test(t, f, t0, dur, exp_time=None):
    """Fit transit, box, ramp, step and flare-decay; compare by BIC.

    Pass if the transit model has the lowest BIC among transit/ramp/step/
    flare-decay by at least SHAPE_DBIC_ALT, and a box is not preferred over
    the transit by more than |SHAPE_DBIC_BOX|. (At low SNR a limb-darkened
    transit and a box are indistinguishable, so demanding a clear win over the
    box would reject most real shallow transits; a strongly box-preferred dip
    has unphysically sharp edges.)
    """
    sigma = _robust_sigma_pt(f)
    n = len(t)
    box = fit_box(t, f, sigma, t0, dur)
    fits = {m.name: m for m in (
        fit_transit(t, f, sigma, t0, dur, exp_time, box=box), box,
        fit_ramp(t, f, sigma, t0), fit_step(t, f, sigma, t0),
        fit_flare_decay(t, f, sigma, t0))}
    bic = {k: m.bic(n) for k, m in fits.items()}
    d_alt = min(bic[k] for k in ("ramp", "step", "flare_decay")) - bic["transit"]
    d_box = bic["box"] - bic["transit"]
    ok = shape_passes(d_alt, d_box)
    return dict(passed=bool(ok), dbic_alt=float(d_alt), dbic_box=float(d_box),
                bic={k: float(v) for k, v in bic.items()}, sigma_ppm=sigma * 1e6, n=n,
                best_alt=min(("ramp", "step", "flare_decay"), key=lambda k: bic[k]),
                transit=fits["transit"].params, fits=fits)


# ------------------------------------------------------------------ duration

def central_duration_h(period_d, rstar, mstar, k=0.0):
    """T14 (h) of a central (b = 0) transit on a circular orbit."""
    a = (G_CGS * mstar * MSUN_G * (period_d * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    ar = a / (rstar * RSUN_CM)
    return float(period_d * 24 / np.pi * np.arcsin(min((1 + k) / ar, 1.0)))


def period_from_duration(t14_h, rstar, mstar, b=0.0, k=0.0):
    """Circular-orbit period (d) giving duration t14 at impact parameter b."""
    from scipy.optimize import brentq
    chord = np.sqrt(max((1 + k) ** 2 - b ** 2, 1e-6))

    def f(p):
        return central_duration_h(p, rstar, mstar, 0.0) * chord - t14_h
    try:
        return float(brentq(f, 0.05, 1e5))
    except ValueError:
        return float("nan")


def duration_check(t14_h, b, rp, rstar, mstar, grazing):
    """Is the fitted duration compatible with an orbit longer than DUR_P_MIN?

    A central transit at P = DUR_P_MIN lasts T_c; non-zero impact parameter or
    eccentricity can shorten it, so durations down to DUR_FRAC * T_c pass.
    Shorter dips pass only if the transit fit itself says the dip is grazing.
    """
    if not (np.isfinite(rstar) and rstar > 0):
        return dict(passed=True, note="no stellar radius", t_min_h=np.nan, p_b0_d=np.nan)
    if not (np.isfinite(mstar) and mstar > 0):
        mstar = rstar            # dwarf approximation, noted
    tc = central_duration_h(DUR_P_MIN, rstar, mstar, rp)
    tmin = DUR_FRAC * tc
    p0 = period_from_duration(t14_h, rstar, mstar, 0.0, rp)
    ok = t14_h >= tmin or grazing
    note = "grazing" if (t14_h < tmin and grazing) else ""
    return dict(passed=bool(ok), t_central_h=tc, t_min_h=tmin, p_b0_d=p0, note=note)


# ------------------------------------------------------------------ edge

def edge_check(time, flux, t0, t14_h, depth=None, gap_h=EDGE_GAP_H, edge_h=EDGE_H):
    """Reject dips close to a data gap unless both sides are resolved.

    Ingress/egress within ``edge_h`` of the edge of the data segment (gaps >
    ``gap_h``) fail unless (a) data cover >= EDGE_RESOLVE_FRAC of a baseline
    window (T14/2 clipped to EDGE_BASE_H) immediately before ingress and after
    egress, and of the transit itself, and (b) the pre- and post-transit
    baselines agree within EDGE_STEP_SIGMA white-noise sigma or
    EDGE_STEP_DEPTH of the depth, whichever is larger. A ramp into a gap
    makes the baselines differ by an amount comparable to the "depth".
    """
    t14 = t14_h / 24
    ti, te = t0 - t14 / 2, t0 + t14 / 2
    gaps = np.where(np.diff(time) > gap_h / 24)[0]
    starts = np.concatenate([[time[0]], time[gaps + 1]])
    ends = np.concatenate([time[gaps], [time[-1]]])
    seg = np.where((starts <= t0) & (ends >= t0))[0]
    if seg.size == 0:
        return dict(passed=False, note="mid-time in a gap", d_before_h=np.nan, d_after_h=np.nan)
    s, e = starts[seg[0]], ends[seg[0]]
    d_before, d_after = (ti - s) * 24, (e - te) * 24
    near = min(d_before, d_after) < edge_h
    cad = np.median(np.diff(time))
    base = np.clip(t14 / 2, EDGE_BASE_H[0] / 24, EDGE_BASE_H[1] / 24)

    def cover(a, b):
        return np.sum((time >= a) & (time < b)) * cad / max(b - a, cad)

    c_pre, c_post, c_in = cover(ti - base, ti), cover(te, te + base), cover(ti, te)
    pre = flux[(time >= ti - base) & (time < ti)]
    post = flux[(time > te) & (time <= te + base)]
    if len(pre) > 2 and len(post) > 2:
        sig = max(_robust_sigma_pt(flux), 1e-9)
        jump = abs(np.median(pre) - np.median(post))
        step = jump / (sig * np.sqrt(1 / len(pre) + 1 / len(post)))
        step_ok = step < EDGE_STEP_SIGMA or (depth is not None and jump < EDGE_STEP_DEPTH * depth)
    else:
        step, step_ok = np.inf, False
    resolved = (min(c_pre, c_post, c_in) >= EDGE_RESOLVE_FRAC) and step_ok
    ok = (not near) or resolved
    return dict(passed=bool(ok), near_edge=bool(near), resolved=bool(resolved),
                d_before_h=float(d_before), d_after_h=float(d_after),
                cover_pre=float(c_pre), cover_post=float(c_post), cover_in=float(c_in),
                baseline_step_sigma=float(step))


# ------------------------------------------------------------------ pixels

PIX_ARCSEC = 21.0          # TESS pixel scale
CENTROID_MAX_PX = 1.0      # difference-image centroid offset allowed regardless of sigma
CENTROID_SIGMA = 3.0       # offsets beyond CENTROID_MAX_PX fail only if this significant
DIFF_MIN_SNR = 3.0         # below this the difference image cannot test anything
NEIGHBOUR_PX = 2.0         # neighbours within this many pixels are tested
NEIGHBOUR_SIGMA = 3.0
NEIGHBOUR_RATIO = 1.5


def tesscut_cutout(ra, dec, sector, size=15):
    """Download a TESScut FFI cutout; returns dict of arrays (nothing kept on disk)."""
    import warnings
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS
    from astroquery.mast import Tesscut
    from .ffi import DEFAULT_BITMASK
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hdul = Tesscut.get_cutouts(coordinates=SkyCoord(ra, dec, unit="deg"), size=size,
                                   sector=sector)[0]
        d = hdul[1].data
        time = np.asarray(d["TIME"], float)
        cube = np.asarray(d["FLUX"], float)
        q = np.asarray(d["QUALITY"], int)
        good = (np.isfinite(time) & ((q & DEFAULT_BITMASK) == 0)
                & np.all(np.isfinite(cube), axis=(1, 2)))
        wcs = WCS(hdul[2].header)
        x, y = wcs.all_world2pix(ra, dec, 0)
        out = dict(time=time[good], cube=cube[good], wcs=wcs, x=float(x), y=float(y),
                   col0=int(hdul[1].header["1CRV4P"]), row0=int(hdul[1].header["2CRV4P"]))
        hdul.close()
    # Local background per frame: median of the faintest 30 % of pixels.
    med = np.median(out["cube"], axis=0)
    faint = med <= np.percentile(med, 30)
    out["cube"] = out["cube"] - np.median(out["cube"][:, faint], axis=1)[:, None, None]
    return out


def difference_image(time, cube, t0, t14):
    """Out-of-transit minus in-transit image, with a per-pixel linear baseline.

    Returns (diff, noise, oot_image, n_in, n_oot). ``diff`` is positive where
    flux dropped during the transit."""
    dt = np.abs(time - t0)
    inn = dt < 0.4 * t14
    if inn.sum() < 2:
        inn = dt < t14 / 2
    gap = 1 / 24
    oot = (dt > t14 / 2 + gap) & (dt < t14 / 2 + gap + max(t14, 3 / 24))
    if inn.sum() < 1 or oot.sum() < 4:
        return None
    x = time - t0
    A = np.column_stack([np.ones(oot.sum()), x[oot]])
    Y = cube[oot].reshape(oot.sum(), -1)
    coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
    resid = Y - A @ coef
    sig = resid.std(axis=0, ddof=2)
    pred_in = (np.column_stack([np.ones(inn.sum()), x[inn]]) @ coef).mean(axis=0)
    obs_in = cube[inn].reshape(inn.sum(), -1).mean(axis=0)
    shape = cube.shape[1:]
    diff = (pred_in - obs_in).reshape(shape)
    noise = (sig * np.sqrt(1 / inn.sum() + 1 / oot.sum())).reshape(shape)
    return diff, noise, cube[oot].mean(axis=0), int(inn.sum()), int(oot.sum())


def centroid_test(diff, noise, x0, y0, radius=2.5):
    """Flux-weighted centroid of the positive difference image near the target.

    Pass if the centroid is within CENTROID_MAX_PX of the target, or the
    offset is not significant (< CENTROID_SIGMA). Also reports where the
    difference image peaks (in SNR) within 5 px."""
    ny, nx = diff.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = np.hypot(xx - x0, yy - y0)
    sel = (r <= radius) & (diff > 0)
    snr_map = diff / np.where(noise > 0, noise, np.inf)
    tot_snr = diff[r <= radius].sum() / np.sqrt((noise[r <= radius] ** 2).sum())
    if sel.sum() == 0 or tot_snr < DIFF_MIN_SNR:
        return dict(passed=True, inconclusive=True, diff_snr=float(tot_snr),
                    offset_px=np.nan, offset_sigma=np.nan, xc=np.nan, yc=np.nan,
                    peak_px=np.nan)
    w = diff[sel]
    xc, yc = (w * xx[sel]).sum() / w.sum(), (w * yy[sel]).sum() / w.sum()
    sx = np.sqrt(((xx[sel] - xc) ** 2 * noise[sel] ** 2).sum()) / w.sum()
    sy = np.sqrt(((yy[sel] - yc) ** 2 * noise[sel] ** 2).sum()) / w.sum()
    off = float(np.hypot(xc - x0, yc - y0))
    sig = float(off / max(np.hypot(sx, sy), 1e-3))
    near = r <= 5
    k = np.argmax(np.where(near, snr_map, -np.inf))
    peak = float(r.flat[k])
    ok = off < CENTROID_MAX_PX or sig < CENTROID_SIGMA
    return dict(passed=bool(ok), inconclusive=False, diff_snr=float(tot_snr), offset_px=off,
                offset_sigma=sig, xc=float(xc), yc=float(yc), peak_px=peak,
                peak_snr=float(snr_map.flat[k]))


def neighbour_test(diff, noise, oot, x0, y0, neighbours, depth):
    """Does a neighbour within NEIGHBOUR_PX show the dip more strongly than the target?

    ``neighbours``: iterable of (tic, x, y, dTmag). Only stars bright enough to
    produce the observed depth (dTmag <= -2.5 log10(depth) + 0.5) are tested.
    Compares the fractional depth in the pixel nearest the neighbour with the
    pixel nearest the target. Neighbours sharing the target's pixel cannot be
    separated here and are reported as unresolved (left to TRICERATOPS)."""
    def frac(ix, iy):
        f = oot[iy, ix]
        if f <= 0:
            return np.nan, np.inf
        return diff[iy, ix] / f, noise[iy, ix] / f

    ny, nx = diff.shape
    tx, ty = int(round(x0)), int(round(y0))
    ft, st = frac(tx, ty)
    dmax = -2.5 * np.log10(max(depth, 1e-6)) + 0.5
    tested, failing, unresolved = [], [], []
    for tic, x, y, dm in neighbours:
        sep = np.hypot(x - x0, y - y0)
        if sep > NEIGHBOUR_PX or dm > dmax:
            continue
        ix, iy = int(round(x)), int(round(y))
        if not (0 <= ix < nx and 0 <= iy < ny):
            continue
        if (ix, iy) == (tx, ty):
            unresolved.append(int(tic))
            continue
        fn, sn = frac(ix, iy)
        z = (fn - ft) / np.hypot(sn, st)
        tested.append(dict(tic=int(tic), sep_px=float(sep), dTmag=float(dm),
                           frac_nb=float(fn), frac_target=float(ft), z=float(z)))
        if z > NEIGHBOUR_SIGMA and fn > NEIGHBOUR_RATIO * ft:
            failing.append(int(tic))
    return dict(passed=not failing, failing=failing, tested=tested, unresolved=unresolved,
                frac_target=float(ft))


def tic_neighbours(ra, dec, wcs, tmag, radius_px=7.5):
    """TIC stars around the target, in cutout pixel coordinates."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astroquery.mast import Catalogs
    t = Catalogs.query_region(SkyCoord(ra, dec, unit="deg"),
                              radius=radius_px * PIX_ARCSEC * u.arcsec, catalog="TIC")
    t = t[np.isfinite(np.asarray(t["Tmag"], float))]
    x, y = wcs.all_world2pix(np.asarray(t["ra"], float), np.asarray(t["dec"], float), 0)
    rows = [(int(i), float(xx), float(yy), float(tm) - tmag, float(d))
            for i, xx, yy, tm, d in zip(t["ID"], x, y, t["Tmag"], t["dstArcSec"])]
    target = min(rows, key=lambda r: r[4])
    return [r[:4] for r in rows if r is not target]


def aperture_mask(oot, x0, y0, radius=2.0, frac=0.2):
    """Pixels within ``radius`` of the target brighter than ``frac`` of the local peak."""
    ny, nx = oot.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = np.hypot(xx - x0, yy - y0)
    near = r <= radius
    peak = oot[near].max()
    ap = near & (oot >= frac * peak)
    ap[int(round(y0)), int(round(x0))] = True
    return ap


# ------------------------------------------------------------------ asteroids

ASTEROID_V_MAX = 19.0      # fainter objects cannot mimic a >~0.1 % dip at T <= 13
ASTEROID_R_ARCSEC = 60.0   # ~3 px: aperture plus the pixels used for its background


def skybot_query(ra, dec, jd, radius_deg=0.2, loc="C57", timeout=60):
    """Known solar-system objects near (ra, dec) at JD (SkyBoT, TESS = C57)."""
    q = {"-ep": f"{jd:.5f}", "-ra": ra, "-dec": dec, "-rd": radius_deg, "-mime": "text",
         "-output": "all", "-loc": loc, "-filter": 0, "-objFilter": "111",
         "-refsys": "EQJ2000"}
    url = ("https://ssp.imcce.fr/webservices/skybot/api/conesearch.php?"
           + urllib.parse.urlencode(q))
    with urllib.request.urlopen(url, timeout=timeout) as r:
        text = r.read().decode()
    return parse_skybot(text)


def parse_skybot(text):
    from astropy.coordinates import Angle
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = [x.strip() for x in line.split("|")]
        if len(f) < 10:
            continue
        rows.append(dict(name=f[1], cls=f[4], ra=Angle(f[2], unit="hourangle").deg,
                         dec=Angle(f[3], unit="deg").deg, v=float(f[5]),
                         sep_arcsec=float(f[7]), dra_arcsec_h=float(f[8]),
                         ddec_arcsec_h=float(f[9])))
    return rows


def asteroid_check(objects, ra, dec, t14_h, v_max=ASTEROID_V_MAX, r_max=ASTEROID_R_ARCSEC):
    """Fail if a known object brighter than v_max passes within r_max of the
    target during the dip (+/- 1 h), using its linear sky motion."""
    hits = []
    for o in objects:
        dx0 = (o["ra"] - ra) * np.cos(np.radians(dec)) * 3600
        dy0 = (o["dec"] - dec) * 3600
        tt = np.linspace(-(t14_h / 2 + 1), t14_h / 2 + 1, 50)
        d = np.hypot(dx0 + o["dra_arcsec_h"] * tt, dy0 + o["ddec_arcsec_h"] * tt).min()
        o = dict(o, min_sep_arcsec=float(d))
        if d < r_max and o["v"] <= v_max:
            hits.append(o)
    return dict(passed=not hits, hits=hits, n_objects=len(objects))


# ------------------------------------------------------------------ catalogues

EXOFOP_TOI = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
EXOFOP_CTOI = "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv"
TESS_EBS = ("https://archive.stsci.edu/hlsps/tess-ebs/"
            "hlsp_tess-ebs_tess_lcf-ffi_s0001-s0026_tess_v1.0_cat.csv")
VILLANOVA_EBS = "https://tessebs.villanova.edu/?order_by=tic&page={page}"


def load_catalogues(cache_dir):
    """ExoFOP TOIs and CTOIs and the TESS EB catalogue (Prsa et al. 2022; the
    MAST HLSP table and the Villanova web catalogue, merged), downloaded once
    into ``cache_dir``."""
    os.makedirs(cache_dir, exist_ok=True)
    out = {}
    for name, url in (("exofop_toi", EXOFOP_TOI), ("exofop_ctoi", EXOFOP_CTOI),
                      ("tess_ebs", TESS_EBS)):
        path = os.path.join(cache_dir, name + ".csv")
        if not os.path.exists(path):
            with urllib.request.urlopen(url, timeout=300) as r, open(path, "wb") as fh:
                fh.write(r.read())
        out[name] = pd.read_csv(path, low_memory=False)
    vpath = os.path.join(cache_dir, "villanova_tics.txt")
    if not os.path.exists(vpath):
        try:
            _fetch_villanova(vpath)
        except Exception:  # noqa: BLE001
            pass
    out["villanova_ebs"] = (set(int(x) for x in open(vpath) if x.strip())
                            if os.path.exists(vpath) else set())
    return out


def _fetch_villanova(path, max_pages=200):
    """TIC IDs listed in the Villanova TESS EB catalogue web pages."""
    import re
    import time
    tics = set()
    for page in range(1, max_pages + 1):
        with urllib.request.urlopen(VILLANOVA_EBS.format(page=page), timeout=60) as r:
            html = r.read().decode()
        found = set(re.findall(r'href="(\d{10})"', html))
        tics |= found
        if not found or f"page={page + 1}" not in html:
            break
        time.sleep(0.5)
    with open(path, "w") as fh:
        fh.write("\n".join(str(int(t)) for t in sorted(tics)))


def catalogue_check(tic, cats):
    """Known TOI/CTOI -> 'known' (passes, but is not new); TESS EB or a TOI
    dispositioned false positive/false alarm -> fail."""
    toi = cats["exofop_toi"]
    ctoi = cats["exofop_ctoi"]
    ebs = cats["tess_ebs"]
    t = toi[toi["TIC ID"] == tic]
    c = ctoi[ctoi["TIC ID"] == tic]
    eb = bool((ebs["tess_id"] == tic).any()) or tic in cats.get("villanova_ebs", set())
    disp = [str(x) for x in t["TFOPWG Disposition"].fillna("")]
    fp = any(d in ("FP", "FA") for d in disp)
    return dict(passed=not (eb or fp), eb=eb, toi_fp=fp,
                tois=[str(x) for x in t["TOI"]], toi_disp=disp,
                ctois=[str(x) for x in c["CTOI"]], known=bool(len(t) or len(c)))


# ------------------------------------------------------------------ TRICERATOPS

FPP_MAX = 0.5              # reject when a false positive is more likely than not
NFPP_MAX = 0.1             # ...or when a nearby-star false positive is plausible
TIC_BG_FIELD = 0.1         # deg^2 of TIC stars used as the background population


def tic_background_file(ra, dec, path, field=TIC_BG_FIELD):
    """Background-star population for TRICERATOPS from the TIC (itself built
    on Gaia DR2), in the CSV format of triceratops.funcs.query_gaia_background.
    Fallback for when the Gaia archive cannot be queried."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astroquery.mast import Catalogs
    radius = np.sqrt(field / np.pi)
    t = Catalogs.query_region(SkyCoord(ra, dec, unit="deg"), radius=radius * u.deg,
                              catalog="TIC").to_pandas()
    t = t[(t.objType == "STAR")]
    df = pd.DataFrame(dict(Mact=t.mass, logg=t.logg, logTe=np.log10(t.Teff), **{"[M/H]": 0.0},
                           TESS=t.Tmag, J=t.Jmag, H=t.Hmag, Ks=t.Kmag)).dropna()
    df.to_csv(path, index=False)
    return path, len(df)


class _CoordCatalogs:
    """Stand-in for astroquery's Catalogs inside triceratops: resolves the
    target by coordinates, since the MAST name resolver is unreachable here."""

    def __init__(self, ra, dec):
        from astroquery.mast import Catalogs
        self._c, self.ra, self.dec = Catalogs, ra, dec

    def query_object(self, name, radius, catalog):
        from astropy.coordinates import SkyCoord
        t = self._c.query_region(SkyCoord(self.ra, self.dec, unit="deg"), radius=radius,
                                 catalog=catalog)
        t.sort("dstArcSec")
        return t

    def query_region(self, *a, **k):
        return self._c.query_region(*a, **k)


GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"


class _GaiaHTTPS:
    """Minimal stand-in for astroquery.gaia.Gaia used by
    triceratops.funcs.query_gaia_background: runs the ADQL on the HTTPS
    synchronous TAP endpoint (astroquery's client also makes plain-HTTP
    requests, which the egress proxy here refuses)."""
    ROW_LIMIT = -1

    class _Job:
        def __init__(self, table):
            self._t = table

        def get_results(self):
            return self._t

    @classmethod
    def launch_job(cls, adql, verbose=False):
        from astropy.io import ascii as asc
        data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
                                       "QUERY": adql}).encode()
        with urllib.request.urlopen(GAIA_TAP_SYNC, data=data, timeout=300) as r:
            text = r.read().decode()
        return cls._Job(asc.read(text, format="csv", fill_values=[("", "nan")]))

    launch_job_async = launch_job


class _TesscutSearch:
    """Stand-in for lightkurve.search_tesscut(...) inside triceratops, which
    would route coordinates through the (unreachable) MAST name resolver."""

    def __init__(self, target, sector):
        self.target, self.sector = target, sector

    def download_all(self, cutout_size):
        from types import SimpleNamespace
        from astroquery.mast import Tesscut
        size = cutout_size[0] if isinstance(cutout_size, (tuple, list)) else cutout_size
        hduls = Tesscut.get_cutouts(coordinates=self.target, size=size, sector=self.sector)
        return [SimpleNamespace(hdu=h) for h in hduls]


def triceratops_fpp(tic, ra, dec, sector, time_from_t0, flux, flux_err, depth, ap_abs,
                    period_range, workdir, n_draws=200_000, seed=0):
    """Run TRICERATOPS; returns dict(fpp, nfpp, scenario probabilities).

    The background-star population comes from Gaia DR3 (TRICERATOPS default);
    if that query fails, a TIC-based population (tic_background_file) is used
    instead, and ``background`` in the result says which."""
    import contextlib
    import io
    import warnings
    import astroquery.gaia
    import triceratops.triceratops as trm
    os.makedirs(workdir, exist_ok=True)
    np.random.seed(seed)
    old = trm.Catalogs, trm.lightkurve.search_tesscut, astroquery.gaia.Gaia
    astroquery.gaia.Gaia = _GaiaHTTPS
    trm.Catalogs = _CoordCatalogs(ra, dec)
    trm.lightkurve.search_tesscut = lambda target, sector=None, **k: _TesscutSearch(target, sector)
    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            tg = trm.target(ID=tic, sectors=np.array([sector]))
            background = "gaia_dr3"
            if tg.trilegal_fname is None:
                path, _ = tic_background_file(ra, dec, f"{tic}_tic_background.csv")
                tg.trilegal_fname = path
                background = "tic"
            n_bg = sum(1 for _ in open(tg.trilegal_fname)) - 1
            if not tg.pix_coords:
                raise RuntimeError("triceratops could not get a TESScut image")
            tg.calc_depths(tdepth=depth, all_ap_pixels=[np.asarray(ap_abs)])
            tg.calc_probs(time=time_from_t0, flux_0=flux, flux_err_0=float(flux_err),
                          P_orb=list(period_range), N=n_draws, parallel=False, verbose=0,
                          exptime=10 / 1440)
        probs = tg.probs[["scenario", "prob"]].groupby("scenario").prob.sum().to_dict()
        return dict(fpp=float(tg.FPP), nfpp=float(tg.NFPP), background=background,
                    n_background=n_bg, n_stars=int(len(tg.stars)),
                    probs={k: float(v) for k, v in probs.items()})
    finally:
        trm.Catalogs, trm.lightkurve.search_tesscut, astroquery.gaia.Gaia = old
        os.chdir(cwd)


def fpp_check(res):
    ok = res["fpp"] < FPP_MAX and res["nfpp"] < NFPP_MAX
    return dict(passed=bool(ok), **res)
