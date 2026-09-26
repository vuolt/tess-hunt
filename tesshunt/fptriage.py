"""Phase 6: triage of existing TESS candidates (TOIs and CTOIs) for false positives.

  load_tables        ExoFOP TOI and CTOI tables (bulk download, cached) -> one table
                     with a group per row: unresolved / planet (CP, KP) / fp (FP, FA)
  gaia_ids           Gaia DR3 source ids and TIC masses (TIC bulk query)
  gaia_orbits        Gaia DR3 nss_two_body_orbit + vari_eclipsing_binary, batched
  companion_mass     companion mass from an SB1 / SB2 / astrometric solution
  period_match       transit period vs Gaia period (x1, x2, x3, 1/2, 1/3)
  lc_checks          all-sector light-curve checks for one candidate: odd/even
                     depths, secondary eclipse at any phase, centroid shift in transit
"""

from __future__ import annotations

import io
import os
import tempfile

import numpy as np
import pandas as pd

from . import net
from .detrend import detrend
from .vetting import EXOFOP_CTOI, EXOFOP_TOI, local_dip_snr

BTJD0 = 2457000.0
M_H_BURN = 0.08                 # M_sun: below this a companion is not a star
MATCH_FACTORS = (1.0, 2.0, 3.0, 0.5, 1 / 3)
MATCH_SIGMA = 3.0
MATCH_FLOOR = 0.001             # fractional period tolerance floor (0.1 %)
PX_ARCSEC = 21.0


# ------------------------------------------------------------------ tables

def load_tables() -> pd.DataFrame:
    """TOIs and CTOIs in one table (ExoFOP bulk CSVs, downloaded once, cached)."""
    toi = pd.read_csv(io.BytesIO(net.get(EXOFOP_TOI, net.service_of(EXOFOP_TOI))), low_memory=False)
    ctoi = pd.read_csv(io.BytesIO(net.get(EXOFOP_CTOI, net.service_of(EXOFOP_CTOI))), low_memory=False)
    disp = toi["TFOPWG Disposition"].fillna("")
    group = np.select([disp.isin(["CP", "KP"]), disp.isin(["FP", "FA"])], ["planet", "fp"], "unresolved")
    t = pd.DataFrame(dict(
        name="TOI " + toi["TOI"].astype(str), kind="TOI", tic=toi["TIC ID"].astype(int),
        disposition=disp, group=group,
        period=pd.to_numeric(toi["Period (days)"], errors="coerce"),
        period_err=pd.to_numeric(toi["Period (days) err"], errors="coerce"),
        epoch_btjd=pd.to_numeric(toi["Epoch (BJD)"], errors="coerce") - BTJD0,
        duration_h=pd.to_numeric(toi["Duration (hours)"], errors="coerce"),
        depth_ppm=pd.to_numeric(toi["Depth (ppm)"], errors="coerce"),
        rstar=pd.to_numeric(toi["Stellar Radius (R_Sun)"], errors="coerce"),
        mstar=pd.to_numeric(toi["Stellar Mass (M_Sun)"], errors="coerce"),
        tmag=pd.to_numeric(toi["TESS Mag"], errors="coerce")))
    t["ra"], t["dec"] = _radec(toi["RA"], toi["Dec"])
    promoted = ctoi["Promoted to TOI"].notna()
    cdisp = ctoi["TFOPWG Disposition"].fillna(ctoi["User Disposition"].fillna(""))
    keep = ~promoted
    c = ctoi[keep]
    cd = cdisp[keep]
    cgroup = np.select([cd.isin(["CP", "KP"]), cd.isin(["FP", "FA"])], ["planet", "fp"], "unresolved")
    cc = pd.DataFrame(dict(
        name="CTOI " + c["CTOI"].astype(str), kind="CTOI", tic=c["TIC ID"].astype(int),
        disposition=cd.values, group=cgroup,
        period=pd.to_numeric(c["Period (days)"], errors="coerce"),
        period_err=pd.to_numeric(c["Period (days) Error"], errors="coerce"),
        epoch_btjd=pd.to_numeric(c["Transit Epoch (BJD)"], errors="coerce") - BTJD0,
        duration_h=pd.to_numeric(c["Duration (hrs)"], errors="coerce"),
        depth_ppm=pd.to_numeric(c["Depth ppm"], errors="coerce"),
        rstar=pd.to_numeric(c["Stellar Radius (R_Sun)"], errors="coerce"),
        mstar=np.nan, tmag=pd.to_numeric(c["TESS Mag"], errors="coerce")))
    cc["ra"], cc["dec"] = _radec(c["RA"], c["Dec"])
    out = pd.concat([t, cc], ignore_index=True)
    out["periodic"] = (out.period > 0) & np.isfinite(out.epoch_btjd) & (out.duration_h > 0) & (out.depth_ppm > 0)
    return out


def _radec(ra, dec):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    ra, dec = ra.astype(str), dec.astype(str)
    if ra.str.contains(":").any():
        c = SkyCoord(ra.values, dec.values, unit=(u.hourangle, u.deg))
        return c.ra.deg, c.dec.deg
    return pd.to_numeric(ra, errors="coerce").values, pd.to_numeric(dec, errors="coerce").values


# ------------------------------------------------------------------ Gaia

def gaia_ids(tics) -> pd.DataFrame:
    """TIC -> Gaia DR3 source id, TIC mass and radius (bulk TIC query, cached)."""
    from . import tic as tic_mod
    tt = tic_mod.query_ids(sorted(set(int(x) for x in tics)))
    df = tt.to_pandas() if hasattr(tt, "to_pandas") else pd.DataFrame(tt)
    df = df.rename(columns={"ID": "tic"})
    df["tic"] = df["tic"].astype(str).str.strip().astype("int64")
    # exact integer parsing: via float64 a 19-digit source_id would be rounded
    df["gaia_id"] = pd.array([int(g) if str(g).strip().isdigit() else pd.NA
                              for g in df["GAIA"].astype(str)], dtype="Int64")
    for c in ("mass", "rad", "Teff", "Tmag"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df[["tic", "gaia_id", "mass", "rad", "Teff", "Tmag"]]


NSS_COLS = ("source_id, nss_solution_type, period, period_error, eccentricity, parallax, "
            "a_thiele_innes, b_thiele_innes, f_thiele_innes, g_thiele_innes, "
            "semi_amplitude_primary, semi_amplitude_secondary, mass_ratio, inclination, "
            "goodness_of_fit, significance")


def _gaia_batched(template: str, ids, batch: int) -> pd.DataFrame:
    from .binarity import _gaia
    ids = sorted({int(i) for i in ids if pd.notna(i)})
    out = []
    for k in range(0, len(ids), batch):
        part = ",".join(str(i) for i in ids[k:k + batch])
        out.append(_gaia(template.format(ids=part)))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def gaia_orbits(source_ids, batch: int = 400) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Gaia DR3 two-body orbits and photometric eclipsing-binary candidates for
    the given sources, queried in batches (the Gaia archive is rate-limited by net)."""
    nss = _gaia_batched(f"SELECT {NSS_COLS} FROM gaiadr3.nss_two_body_orbit WHERE source_id IN ({{ids}})",
                        source_ids, batch)
    eb = _gaia_batched("SELECT source_id, frequency, frequency_error, global_ranking, model_type "
                       "FROM gaiadr3.vari_eclipsing_binary WHERE source_id IN ({ids})", source_ids, batch)
    for df in (nss, eb):
        if len(df):
            df["source_id"] = df["source_id"].astype("int64")
    return nss, eb


def _solve_m2(f, m1):
    """Smallest M2 with M2^3 / (M1 + M2)^2 = f (Newton iterations, M_sun)."""
    if not (f > 0 and m1 > 0):
        return np.nan
    m2 = max(f ** (1 / 3) * m1 ** (2 / 3), 1e-4)
    for _ in range(60):
        g = m2 ** 3 / (m1 + m2) ** 2 - f
        dg = (3 * m2 ** 2 * (m1 + m2) - 2 * m2 ** 3) / (m1 + m2) ** 3
        m2 = max(m2 - g / dg, 1e-6)
    return float(m2)


def thiele_innes_a0(A, B, F, G):
    """Photocentre semi-major axis (same unit as A..G) from Thiele-Innes elements."""
    u = (A ** 2 + B ** 2 + F ** 2 + G ** 2) / 2
    v = A * G - B * F
    return float(np.sqrt(u + np.sqrt(max((u + v) * (u - v), 0.0))))


def companion_mass(sol: dict, m1: float) -> tuple[float, str]:
    """Companion mass (M_sun) and how it was obtained.

    SB1: mass function with sin i = 1, appropriate only when the companion
         eclipses (i.e. the transit period matches the orbit).
    SB2: M2 = q M1 from the two semi-amplitudes (mass_ratio).
    astrometric: photocentre orbit, assuming a dark companion (a lower limit
         if the companion contributes light)."""
    P = sol.get("period")
    e = sol.get("eccentricity")
    e = 0.0 if e is None or not np.isfinite(e) else e
    if not (P and np.isfinite(P) and m1 and np.isfinite(m1)):
        return np.nan, ""
    q = sol.get("mass_ratio")
    K1, K2 = sol.get("semi_amplitude_primary"), sol.get("semi_amplitude_secondary")
    if not (q is not None and np.isfinite(q) and q > 0) and all(
            x is not None and np.isfinite(x) and x > 0 for x in (K1, K2)):
        q = K1 / K2                                     # Gaia DR3 SB2 rows leave mass_ratio empty
    if q is not None and np.isfinite(q) and q > 0:
        return float(q * m1), "SB2 mass ratio (K1/K2)"
    K = sol.get("semi_amplitude_primary")
    if K is not None and np.isfinite(K) and K > 0:
        f = 1.0361e-7 * (1 - e ** 2) ** 1.5 * K ** 3 * P
        return _solve_m2(f, m1), "SB1 mass function (sin i = 1)"
    A, B, F, G = (sol.get(k) for k in ("a_thiele_innes", "b_thiele_innes", "f_thiele_innes", "g_thiele_innes"))
    plx = sol.get("parallax")
    if all(x is not None and np.isfinite(x) for x in (A, B, F, G)) and plx and plx > 0:
        a0_au = thiele_innes_a0(A, B, F, G) / plx
        f = a0_au ** 3 / (P / 365.25) ** 2
        return _solve_m2(f, m1), "astrometric orbit (dark companion; lower limit)"
    return np.nan, ""


def period_match(P, sP, Pg, sPg, factors=MATCH_FACTORS):
    """Best factor k with P ~ k * Pg; returns (k, n_sigma, rel_diff) or None."""
    if not (P and Pg and np.isfinite(P) and np.isfinite(Pg)):
        return None
    sP = sP if sP and np.isfinite(sP) else 0.0
    sPg = sPg if sPg and np.isfinite(sPg) else 0.0
    best = None
    for k in factors:
        tol = np.sqrt(sP ** 2 + (k * sPg) ** 2 + (MATCH_FLOOR * P) ** 2)
        d = abs(P - k * Pg)
        ns = d / tol
        if ns <= MATCH_SIGMA and (best is None or ns < best[1]):
            best = (k, float(ns), float(d / P))
    return best


# ------------------------------------------------------------------ light curves

def _tic_path(tic):
    t = f"{int(tic):016d}"
    return f"{t[0:4]}/{t[4:8]}/{t[8:12]}/{t[12:16]}"


def lc_urls(tic, sector):
    """SPOC 2-min (S3 tid) first, then the TESS-SPOC FFI HLSP (S3, then MAST).
    Both carry PDCSAP flux and flux-weighted centroids."""
    urls = []
    keys = net.s3_list(f"tess/public/tid/s{sector:04d}/{_tic_path(tic)}/")
    urls += [f"{net.S3}/{k}" for k in keys if k.endswith("-s_lc.fits")]
    from . import ffi
    urls += ffi.lc_urls(int(tic), int(sector))
    return urls


def read_lc(path):
    from astropy.io import fits
    from .ffi import DEFAULT_BITMASK
    with fits.open(path, memmap=False) as h:
        d = h[1].data
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d["PDCSAP_FLUX"], float)
        cx = np.asarray(d["MOM_CENTR1"], float)
        cy = np.asarray(d["MOM_CENTR2"], float)
        q = np.asarray(d["QUALITY"], int)
    ok = np.isfinite(t) & np.isfinite(f) & (f > 0) & ((q & DEFAULT_BITMASK) == 0) & np.isfinite(cx) & np.isfinite(cy)
    t, f, cx, cy = t[ok], f[ok], cx[ok], cy[ok]
    return t, f / np.median(f), cx, cy


def _transit_mask(t, P, t0, dur, pad=1.0):
    if not (P and P > 0):
        return np.zeros(len(t), bool)
    ph = (t - t0 + 0.5 * P) % P - 0.5 * P
    return np.abs(ph) < pad * dur / 2 + 1 / 48


def _event_centroid(t, c, t0, dur):
    """In-transit minus out-of-transit centroid, with a linear local baseline."""
    dt = t - t0
    inn = np.abs(dt) < 0.4 * dur
    lo, hi = dur / 2 + 1 / 24, dur / 2 + 1 / 24 + max(dur, 3 / 24)
    oot = (np.abs(dt) > lo) & (np.abs(dt) < hi)
    if inn.sum() < 2 or oot.sum() < 6:
        return None
    cf = np.polyfit(dt[oot], c[oot], 1)
    r = c - np.polyval(cf, dt)
    sig = 1.4826 * np.median(np.abs(r[oot] - np.median(r[oot])))
    return float(r[inn].mean()), float(sig * np.sqrt(1 / inn.sum() + 1 / oot.sum()))


def sector_measure(t, f, cx, cy, P, t0, dur, depth, others):
    """Per-transit depths and centroid shifts in one sector, and the detrended
    out-of-transit flux (for the secondary-eclipse scan)."""
    ex = _transit_mask(t, P, t0, dur, pad=1.5)
    for (Po, t0o, do) in others:                        # other planets on the star
        ex |= _transit_mask(t, Po, t0o, do, pad=1.5)
    window = max(3 * dur, 0.5)
    mask, flat, _ = detrend(t, f, window=window, exclude=ex)
    t, flat, cx, cy = t[mask], flat[mask], cx[mask], cy[mask]
    n0, n1 = int(np.ceil((t.min() - t0) / P)), int(np.floor((t.max() - t0) / P))
    events = []
    for n in range(n0, n1 + 1):
        tc = t0 + n * P
        d, snr = local_dip_snr(t, flat, tc, dur)
        if not np.isfinite(d) or not np.isfinite(snr) or snr == 0:
            continue
        ev = dict(n=n, depth=float(d), err=float(abs(d / snr)))
        cxs, cys = _event_centroid(t, cx, tc, dur), _event_centroid(t, cy, tc, dur)
        if cxs and cys:
            ev.update(dcx=cxs[0], ecx=cxs[1], dcy=cys[0], ecy=cys[1])
        events.append(ev)
    masked = _transit_mask(t, P, t0, dur, pad=1.5)
    for (Po, t0o, do) in others:
        masked |= _transit_mask(t, Po, t0o, do, pad=1.5)
    return events, t[~masked], flat[~masked]


def secondary_scan(t, f, P, t0, dur, min_frac=0.4):
    """Box search for a second dip of the same duration at any phase away from the
    transit (eccentric orbits put it anywhere). Returns best (phase, depth, snr)."""
    if len(t) < 50:
        return None
    ph = ((t - t0) / P) % 1.0
    o = np.argsort(ph)
    ph, f = ph[o], f[o]
    w = dur / P
    step = max(w / 4, 1e-4)
    sig = 1.4826 * np.median(np.abs(np.diff(f))) / np.sqrt(2)
    cad = np.median(np.diff(np.sort(t)))
    n_full = max(dur / cad, 1)
    best = None
    cs = np.concatenate([[0.0], np.cumsum(1 - f)])
    for c in np.arange(0.0, 1.0, step):
        if min(c, 1 - c) < 1.5 * w:                     # away from the primary transit
            continue
        lo, hi = np.searchsorted(ph, c - 0.4 * w), np.searchsorted(ph, c + 0.4 * w)
        n = hi - lo
        if n < min_frac * 0.8 * n_full:              # need most of at least one event
            continue
        d = (cs[hi] - cs[lo]) / n
        snr = d / (sig / np.sqrt(n))
        if best is None or snr > best[2]:
            best = (float(c), float(d), float(snr), int(n))
    return best


def combine(events):
    """Odd/even depth difference and a combined centroid-shift significance."""
    ev = [e for e in events if np.isfinite(e["err"]) and e["err"] > 0]
    out = dict(n_events=len(ev))
    if not ev:
        return out
    d = np.array([e["depth"] for e in ev])
    s = np.array([e["err"] for e in ev])
    w = 1 / s ** 2
    out["depth"] = float((w * d).sum() / w.sum())
    out["depth_err"] = float(1 / np.sqrt(w.sum()))
    odd = np.array([e["n"] % 2 == 1 for e in ev])
    if odd.sum() >= 1 and (~odd).sum() >= 1:
        do = (w[odd] * d[odd]).sum() / w[odd].sum()
        de = (w[~odd] * d[~odd]).sum() / w[~odd].sum()
        so, se = 1 / np.sqrt(w[odd].sum()), 1 / np.sqrt(w[~odd].sum())
        out.update(depth_odd=float(do), depth_even=float(de),
                   oddeven_sigma=float(abs(do - de) / np.hypot(so, se)),
                   oddeven_frac=float(abs(do - de) / max(abs(do + de) / 2, 1e-9)))
    cen = [e for e in ev if "dcx" in e and e["ecx"] > 0 and e["ecy"] > 0]
    if cen:
        # centroid shifts are combined per sector-orientation-free chi^2 per event
        chi2 = sum((e["dcx"] / e["ecx"]) ** 2 + (e["dcy"] / e["ecy"]) ** 2 for e in cen)
        dof = 2 * len(cen)
        from scipy.stats import chi2 as chi2d, norm
        p = chi2d.sf(chi2, dof)
        out["centroid_sigma"] = max(float(norm.isf(max(p, 1e-300))), 0.0)
        wx = np.array([1 / e["ecx"] ** 2 for e in cen])
        wy = np.array([1 / e["ecy"] ** 2 for e in cen])
        mx = float((wx * [e["dcx"] for e in cen]).sum() / wx.sum())
        my = float((wy * [e["dcy"] for e in cen]).sum() / wy.sum())
        out["centroid_shift_px"] = float(np.hypot(mx, my))
        dep = max(out["depth"], 1e-6)
        out["source_offset_px"] = float(np.hypot(mx, my) / dep)   # |dc| / depth
    return out


def lc_checks(tic, ra, dec, P, t0, dur_h, depth_ppm, others=(), last_sector=108, tmp_dir=None):
    """All-sector checks for one candidate. Light curves are downloaded (S3 first)
    to a temporary folder and deleted after reading."""
    from tess_stars2px import tess_stars2px_function_entry as t2p
    dur = dur_h / 24.0
    out = t2p(int(tic), float(ra), float(dec))
    sectors = sorted({int(s) for s in out[3] if 0 < s <= last_sector})
    events, st, sf, used = [], [], [], []
    tmp_dir = tmp_dir or tempfile.mkdtemp(prefix="p6_")
    for s in sectors:
        try:
            urls = lc_urls(tic, s)
            if not urls:
                continue
            path = net.download(urls, net.service_of, dest=os.path.join(tmp_dir, f"{tic}_{s}.fits"),
                                cache=False)
        except (FileNotFoundError, net.ServiceError, OSError):
            continue
        try:
            t, f, cx, cy = read_lc(path)
        finally:
            if os.path.exists(path) and path.startswith(tmp_dir):
                os.remove(path)
        if len(t) < 200:
            continue
        ev, tt, ff = sector_measure(t, f, cx, cy, P, t0, dur, depth_ppm * 1e-6, others)
        for e in ev:
            e["sector"] = s
        events += ev
        st.append(tt)
        sf.append(ff)
        used.append(s)
    res = dict(tic=int(tic), sectors=used, **combine(events))
    if st:
        sec = secondary_scan(np.concatenate(st), np.concatenate(sf), P, t0, dur)
        if sec:
            res.update(sec_phase=sec[0], sec_depth=sec[1], sec_snr=sec[2], sec_n=sec[3])
            if res.get("depth"):
                res["sec_ratio"] = sec[1] / max(res["depth"], 1e-9)
    return res
