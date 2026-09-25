"""Automated vetting of single-event dips from a full-sector search.

Labels each dip with:
  common_mode     many other stars within CM_RADIUS_DEG dip at the same time
                  (local Poisson test), i.e. a detector systematic
  partial         the box reaches within 1 h of a data gap
  category        common_mode | repeating | secondary | single | single_partial
  rp_rjup         companion radius implied by the depth and the TIC stellar radius
  tier            A: planet-sized, not partial, outside sector-wide pile-up
                  windows; B: other candidates
"""

from __future__ import annotations

import urllib.parse
import urllib.request

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
