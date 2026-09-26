"""Everything TESS has on one star: other sectors, second dips, allowed periods.

  observed_sectors   sectors a star fell on silicon (tess-point, offline),
                     limited to sectors whose data are public
  best_lc            the best light curve per sector: SPOC 2-min (S3), then
                     TESS-SPOC FFI, then QLP FFI (S3 mirror first, MAST second),
                     then our own TESScut aperture photometry
  candidate_detrend  detrending that preserves a transit of the candidate's duration
  matched_dips       dips of the candidate's duration in one sector (duotransit search)
  ruled_out_map      where the data exclude a transit of the candidate's depth
  period_scan        which orbital periods survive all of TESS's coverage
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import ffi, net
from .detect import box_search
from .detrend import detrend
from .lightcurves import SectorLC

TID_PREFIX = "tess/public/tid"
QLP_S3 = f"{net.S3}/mast/hlsp/qlp"
QLP_MAST = "https://archive.stsci.edu/hlsps/qlp"

DUO_SNR = 5.0                 # matched-filter SNR for a possible second transit
DUO_DEPTH = (0.5, 2.0)        # depth ratio to the original for "consistent"
EXCLUDE_SIGMA = 3.0           # a transit is ruled out if the box depth there is
EXCLUDE_FRAC = 0.5            # < depth - 3 sigma AND < 0.5 x depth
MIN_COVER = 0.5               # fraction of in-transit cadences needed to rule out


def tic_path(tic: int) -> str:
    t = f"{int(tic):016d}"
    return f"{t[0:4]}/{t[4:8]}/{t[8:12]}/{t[12:16]}"


def latest_sector() -> int:
    """Most recent sector with FFIs on the S3 mirror."""
    prefixes = net.s3_list("tess/public/ffi/", delimiter="/")
    secs = [int(p.rstrip("/").split("/s")[-1]) for p in prefixes if "/s" in p]
    return max(s for s in secs if s < 1000)


def observed_sectors(tic: int, ra: float, dec: float, last: int | None = None) -> list[int]:
    from tess_stars2px import tess_stars2px_function_entry as t2p
    last = last or latest_sector()
    out = t2p(int(tic), float(ra), float(dec))
    return sorted({int(s) for s in out[3] if 0 < s <= last})


# ------------------------------------------------------------------ products

@dataclass
class Product:
    kind: str          # spoc2min | tess-spoc | qlp | tesscut
    sector: int
    lc: SectorLC
    cadence_min: float


def _spoc_2min(tic, sector):
    keys = net.s3_list(f"{TID_PREFIX}/s{sector:04d}/{tic_path(tic)}/")
    lcs = [k for k in keys if k.endswith("-s_lc.fits")]
    if not lcs:
        raise FileNotFoundError("no SPOC 2-min light curve")
    path = net.download([f"{net.S3}/{lcs[0]}"], net.service_of)
    return ffi.read(path)[0]


def _tess_spoc(tic, sector):
    return ffi.read(ffi.download(tic, sector, cache=True))[0]


def _qlp(tic, sector):
    from astropy.io import fits
    t16 = f"{int(tic):016d}"
    name = f"s{sector:04d}/{tic_path(tic)}/hlsp_qlp_tess_ffi_s{sector:04d}-{t16}_tess_v01_llc.fits"
    path = net.download([f"{QLP_S3}/{name}", f"{QLP_MAST}/{name}"], net.service_of)
    with fits.open(path, memmap=False) as h:
        d = h[1].data
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d["SAP_FLUX"], float)
        q = np.asarray(d["QUALITY"], int)
    good = np.isfinite(t) & np.isfinite(f) & (q == 0) & (f > 0)
    t, f = t[good], f[good]
    med = np.median(f)
    return SectorLC(int(tic), sector, t, f / med, np.full(len(t), np.nan))


def _tesscut(tic, sector, ra, dec, tmag):
    from .vetting import aperture_mask, tesscut_cutout
    cut = tesscut_cutout(ra, dec, sector)
    med = np.median(cut["cube"], axis=0)
    ap = aperture_mask(med, cut["x"], cut["y"])
    f = cut["cube"][:, ap].sum(axis=1)
    good = np.isfinite(f) & (f > 0)
    t, f = cut["time"][good], f[good]
    return SectorLC(int(tic), sector, t, f / np.median(f), np.full(len(t), np.nan))


def best_lc(tic: int, ra: float, dec: float, tmag: float, sector: int,
            kinds=("spoc2min", "tess-spoc", "qlp", "tesscut")) -> Product | None:
    """First available light curve for (star, sector) in order of preference."""
    for kind in kinds:
        try:
            if kind == "spoc2min":
                lc = _spoc_2min(tic, sector)
            elif kind == "tess-spoc":
                lc = _tess_spoc(tic, sector)
            elif kind == "qlp":
                lc = _qlp(tic, sector)
            else:
                lc = _tesscut(tic, sector, ra, dec, tmag)
        except (FileNotFoundError, KeyError, OSError):
            continue
        if len(lc.time) < 200:
            continue
        cad = float(np.median(np.diff(lc.time)) * 1440)
        return Product(kind, sector, lc, cad)
    return None


# ------------------------------------------------------------------ search

def candidate_detrend(lc: SectorLC, t14: float):
    """Detrend with a window >= 3x the candidate's duration (and >= 1 d), mask
    the strongest candidate-like dips, detrend again. Returns (time, flat)."""
    w = max(3 * t14, 1.0)
    mask, flat, trend = detrend(lc.time, lc.flux, window=w)
    t, f = lc.time[mask], flat[mask]
    res = box_search(t, f, durations_h=(t14 * 24,), threshold=DUO_SNR, max_events=3,
                     step_min=_step_min(t14, t))
    if res.events:
        ex = np.zeros(len(lc.time), bool)
        for e in res.events:
            ex |= np.abs(lc.time - e.t0) < t14
        mask, flat, trend = detrend(lc.time, lc.flux, window=w, exclude=ex)
        t, f = lc.time[mask], flat[mask]
    return t, f


def _step_min(t14, t):
    cad = np.median(np.diff(t)) * 1440
    return float(np.clip(t14 * 1440 / 4, cad, 10.0))


def matched_dips(t, f, t14, depth, snr_min=DUO_SNR, max_events=5):
    """Box search at the candidate's own duration; returns events with a
    'consistent' flag (depth within DUO_DEPTH of the original)."""
    res = box_search(t, f, durations_h=(t14 * 24,), threshold=snr_min, max_events=max_events,
                     step_min=_step_min(t14, t))
    out = []
    for e in res.events:
        r = e.depth / depth if depth > 0 else np.nan
        out.append(dict(t0=e.t0, depth_ppm=e.depth * 1e6, snr=e.snr, depth_ratio=r,
                        consistent=bool(DUO_DEPTH[0] <= r <= DUO_DEPTH[1])))
    return out, float(res.sigma[0]) if len(res.sigma) else np.nan


# ------------------------------------------------------------------ periods

def ruled_out_map(t, f, t14, depth, sigma_box=None):
    """For box centres on a t14/4 grid across the data: 0 = no usable data,
    1 = data present and a transit of this depth is allowed (or present),
    2 = data exclude a transit of this depth. Returns (grid, status, measured depth)."""
    order = np.argsort(t)
    t, f = t[order], f[order]
    cad = np.median(np.diff(t))
    step = t14 / 4
    grid = np.arange(t[0], t[-1] + step, step)
    lo = np.searchsorted(t, grid - t14 / 2)
    hi = np.searchsorted(t, grid + t14 / 2)
    n = hi - lo
    cs = np.concatenate([[0.0], np.cumsum(1 - f)])
    d = np.where(n > 0, (cs[hi] - cs[lo]) / np.maximum(n, 1), np.nan)
    cover = n * cad / t14
    if sigma_box is None or not np.isfinite(sigma_box):
        sig_pt = 1.4826 * np.median(np.abs(np.diff(f))) / np.sqrt(2)
        sigma_box = sig_pt / np.sqrt(max(t14 / cad, 1))
    n_typ = max(t14 / cad, 1)
    sig = sigma_box * np.sqrt(n_typ / np.maximum(n, 1))
    status = np.zeros(len(grid), np.int8)
    have = cover >= MIN_COVER
    status[have] = 1
    excl = have & (d < depth - EXCLUDE_SIGMA * sig) & (d < EXCLUDE_FRAC * depth)
    status[excl] = 2
    return grid, status, d


def period_scan(t0, t14, maps, pmin=1.0, pmax=None, chunk=2000):
    """Which periods are allowed by all the coverage in ``maps``.

    ``maps`` is a list of (grid, status) from ruled_out_map. A period is
    excluded if any predicted transit t0 + nP (n != 0) lands on a grid point
    whose status is 2. Returns dict(P, allowed, n_observed) where n_observed
    counts predicted transits that fell on usable data (status >= 1)."""
    grid = np.concatenate([m[0] for m in maps])
    status = np.concatenate([m[1] for m in maps])
    o = np.argsort(grid)
    grid, status = grid[o], status[o]
    step = t14 / 4
    span = max(abs(grid[0] - t0), abs(grid[-1] - t0))
    pmax = pmax or span * 1.05
    dlnp = min(step / (2 * span), 1e-3)
    P = np.exp(np.arange(np.log(pmin), np.log(pmax), dlnp))
    allowed = np.ones(len(P), bool)
    n_obs = np.zeros(len(P), int)
    tmin, tmax = grid[0], grid[-1]
    for i in range(0, len(P), chunk):
        p = P[i:i + chunk]
        nlo = np.ceil((tmin - t0) / p).astype(int)
        nhi = np.floor((tmax - t0) / p).astype(int)
        kmax = int((nhi - nlo).max()) + 1
        ks = np.arange(kmax)
        n = nlo[:, None] + ks[None, :]
        valid = (n <= nhi[:, None]) & (n != 0)
        ep = t0 + n * p[:, None]
        j = np.clip(np.searchsorted(grid, ep), 1, len(grid) - 1)
        jj = np.where(np.abs(grid[j] - ep) < np.abs(grid[j - 1] - ep), j, j - 1)
        near = np.abs(grid[jj] - ep) <= step / 2 + 1e-9
        st = np.where(valid & near, status[jj], 0)
        allowed[i:i + chunk] = ~(st == 2).any(axis=1)
        n_obs[i:i + chunk] = (st >= 1).sum(axis=1)
    return dict(P=P, allowed=allowed, n_observed=n_obs, pmax=pmax, span=span)


def allowed_intervals(scan, min_width_frac=0.0):
    """Contiguous allowed period ranges [(lo, hi), ...]."""
    P, a = scan["P"], scan["allowed"]
    out, start = [], None
    for i, ok in enumerate(a):
        if ok and start is None:
            start = i
        if (not ok or i == len(a) - 1) and start is not None:
            end = i if ok else i - 1
            out.append((float(P[start]), float(P[end])))
            start = None
    return [(lo, hi) for lo, hi in out if hi >= lo * (1 + min_width_frac)]


def alias_periods(t0, t1, pmin=1.0):
    """Periods consistent with two transits at t0 and t1: |t1 - t0| / n."""
    dt = abs(t1 - t0)
    n = np.arange(1, int(dt / pmin) + 1)
    return dt / n
