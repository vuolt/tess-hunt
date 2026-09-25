"""Single-event box search.

For each trial duration a box is slid across the detrended light curve on a
fixed time grid. The box depth is the mean flux deficit inside it. The noise
for that duration is measured *empirically* from the spread of box depths
across the whole sector (robust MAD), so red noise on the transit timescale
is included automatically rather than assumed white. SNR = depth / noise.

A box only counts if the data cover enough of it, and if there is
out-of-transit baseline on both sides. The two-sided baseline rule throws away
most false dips at the edges of data gaps (scattered-light ramps, thermal
settling after downlinks), at the price of missing transits that are cut in
half by a gap.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

DEFAULT_DURATIONS_H = (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0)


@dataclass
class Event:
    t0: float        # BTJD, box centre
    duration: float  # days
    depth: float     # fractional
    snr: float
    sector: int | None = None

    @property
    def duration_h(self) -> float:
        return self.duration * 24.0


@dataclass
class BoxSearch:
    grid: np.ndarray            # box centres, BTJD
    durations: np.ndarray       # days
    snr: np.ndarray             # (n_dur, n_grid), NaN where box invalid
    depth: np.ndarray           # (n_dur, n_grid)
    sigma: np.ndarray           # (n_dur,) empirical noise per duration
    events: list[Event] = field(default_factory=list)

    @property
    def max_snr(self) -> np.ndarray:
        """Best SNR over all durations at each grid time."""
        with np.errstate(all="ignore"):
            out = np.nanmax(np.where(np.isfinite(self.snr), self.snr, -np.inf), axis=0)
        out[~np.isfinite(out)] = np.nan
        return out


def _count(t: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    return np.searchsorted(t, hi) - np.searchsorted(t, lo)


def box_search(time: np.ndarray, flux: np.ndarray,
               durations_h=DEFAULT_DURATIONS_H, step_min: float = 10.0,
               min_coverage: float = 0.6, baseline_coverage: float = 0.5,
               min_baseline_h: float = 6.0,
               threshold: float = 7.0, max_events: int = 5) -> BoxSearch:
    """Search a detrended light curve (baseline ~1) for single box-shaped dips."""
    order = np.argsort(time)
    t = time[order]
    x = 1.0 - flux[order]          # positive = dimming
    cs = np.concatenate([[0.0], np.cumsum(x)])
    cadence = np.median(np.diff(t))
    grid = np.arange(t[0], t[-1], step_min / 1440.0)
    durations = np.asarray(durations_h, float) / 24.0

    snr = np.full((len(durations), len(grid)), np.nan)
    depth = np.full_like(snr, np.nan)
    sigma = np.full(len(durations), np.nan)

    for i, dur in enumerate(durations):
        lo = np.searchsorted(t, grid - dur / 2)
        hi = np.searchsorted(t, grid + dur / 2)
        n = hi - lo
        cover = n * cadence / dur
        side = max(dur, min_baseline_h / 24.0)  # baseline window on each side
        n_left = _count(t, grid - dur / 2 - side, grid - dur / 2)
        n_right = _count(t, grid + dur / 2, grid + dur / 2 + side)
        valid = ((cover >= min_coverage)
                 & (n_left * cadence / side >= baseline_coverage)
                 & (n_right * cadence / side >= baseline_coverage))
        if valid.sum() < 10:
            continue
        d = np.where(n > 0, (cs[hi] - cs[lo]) / np.maximum(n, 1), np.nan)
        d = d - np.median(d[valid])
        s = 1.4826 * np.median(np.abs(d[valid]))
        # Scale for boxes that are only partly filled with data.
        n_typ = np.median(n[valid])
        s_box = s * np.sqrt(n_typ / np.maximum(n, 1))
        sigma[i] = s
        depth[i, valid] = d[valid]
        snr[i, valid] = d[valid] / s_box[valid]

    result = BoxSearch(grid, durations, snr, depth, sigma)
    result.events = extract_events(result, threshold, max_events)
    return result


def extract_events(res: BoxSearch, threshold: float, max_events: int) -> list[Event]:
    """Greedy peak picking: take the best box, mask around it, repeat."""
    best = res.max_snr.copy()
    with np.errstate(all="ignore"):
        best_idx = np.nanargmax(np.where(np.isfinite(res.snr), res.snr, -np.inf), axis=0)
    events = []
    while len(events) < max_events and np.any(np.isfinite(best)):
        k = int(np.nanargmax(best))
        if best[k] < threshold:
            break
        i = best_idx[k]
        dur = res.durations[i]
        events.append(Event(t0=float(res.grid[k]), duration=float(dur),
                            depth=float(res.depth[i, k]), snr=float(best[k])))
        # Blank every box (of any trial duration) that overlaps this event.
        # Allow for the true transit being up to 2x longer than the best box.
        reach = dur + res.durations.max() / 2
        best[np.abs(res.grid - res.grid[k]) < reach] = np.nan
    return events
