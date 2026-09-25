"""Injection-recovery of synthetic single transits into real light curves."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .lightcurves import SectorLC
from .pipeline import SearchConfig, search_sector


def trapezoid(time: np.ndarray, t0: float, depth: float, duration: float,
              ingress_frac: float = 0.1) -> np.ndarray:
    """Multiplicative trapezoid transit model (1 out of transit).

    duration is first-to-fourth contact; ingress and egress each take
    ``ingress_frac * duration``.
    """
    dt = np.abs(time - t0)
    half = duration / 2
    ing = max(ingress_frac * duration, 1e-6)
    shape = np.clip((half - dt) / ing, 0.0, 1.0)
    return 1.0 - depth * shape


def inject(lc: SectorLC, t0: float, depth: float, duration: float,
           ingress_frac: float = 0.1) -> SectorLC:
    out = lc.copy()
    out.flux = out.flux * trapezoid(out.time, t0, depth, duration, ingress_frac)
    return out


@dataclass
class Trial:
    sector: int
    t0: float
    depth_ppm: float
    duration_h: float
    expected_snr: float     # mean in-box depth / empirical noise for this duration
    measured_snr: float     # best box SNR within tolerance of t0 (NaN if no valid box)
    recovered: bool         # an event above threshold within tolerance of t0
    top_ranked: bool        # ...and it is the strongest event in the sector
    det_t0: float
    det_depth_ppm: float
    det_duration_h: float
    n_events: int


def run_trial(lc: SectorLC, t0: float, depth: float, duration: float,
              cfg: SearchConfig, noise: dict[float, float],
              ingress_frac: float = 0.1) -> Trial:
    """Inject one transit, run the full pipeline, score the outcome.

    ``noise`` maps trial duration (days) -> empirical box noise measured on the
    un-injected sector; used only to compute ``expected_snr``.
    """
    ss = search_sector(inject(lc, t0, depth, duration, ingress_frac), cfg)
    res = ss.search
    tol = max(duration / 2, 0.5 / 24)

    near = np.abs(res.grid - t0) < tol
    msnr = np.nanmax(res.max_snr[near]) if near.any() and np.isfinite(res.max_snr[near]).any() else np.nan

    hits = [e for e in res.events if abs(e.t0 - t0) < tol]
    best = max(hits, key=lambda e: e.snr) if hits else None
    top = bool(best is not None and best is res.events[0])

    durs = np.array(sorted(noise))
    nearest = durs[np.argmin(np.abs(durs - duration))]
    esnr = depth * (1 - ingress_frac) / noise[nearest] if noise[nearest] > 0 else np.nan

    return Trial(
        sector=lc.sector, t0=t0, depth_ppm=round(depth * 1e6, 6),
        duration_h=round(duration * 24, 6),
        expected_snr=float(esnr), measured_snr=float(msnr),
        recovered=best is not None, top_ranked=top,
        det_t0=best.t0 if best else np.nan,
        det_depth_ppm=best.depth * 1e6 if best else np.nan,
        det_duration_h=best.duration_h if best else np.nan,
        n_events=len(res.events),
    )


def trial_row(t: Trial) -> dict:
    return asdict(t)
