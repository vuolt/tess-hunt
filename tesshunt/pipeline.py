"""Detrend + search, per sector."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .detect import DEFAULT_DURATIONS_H, BoxSearch, box_search
from .detrend import detrend
from .lightcurves import SectorLC


@dataclass
class SearchConfig:
    window: float = 3.0            # biweight window, days (>= 3x longest duration)
    durations_h: tuple = DEFAULT_DURATIONS_H
    step_min: float = 10.0
    threshold: float = 7.0
    max_events: int = 5
    baseline_coverage: float = 0.5   # data coverage required on EACH side of a box
    min_baseline_h: float = 6.0      # minimum width of each baseline window


@dataclass
class SectorSearch:
    lc: SectorLC
    mask: np.ndarray       # points kept after upward-outlier clipping
    trend: np.ndarray
    detrended: np.ndarray
    search: BoxSearch

    @property
    def time(self) -> np.ndarray:
        return self.lc.time[self.mask]

    @property
    def flat(self) -> np.ndarray:
        return self.detrended[self.mask]


def _run(lc: SectorLC, cfg: SearchConfig, exclude=None):
    mask, flat, trend = detrend(lc.time, lc.flux, window=cfg.window, exclude=exclude)
    res = box_search(lc.time[mask], flat[mask], durations_h=cfg.durations_h,
                     step_min=cfg.step_min, threshold=cfg.threshold,
                     max_events=cfg.max_events,
                     baseline_coverage=cfg.baseline_coverage,
                     min_baseline_h=cfg.min_baseline_h)
    return mask, flat, trend, res


def search_sector(lc: SectorLC, cfg: SearchConfig | None = None) -> SectorSearch:
    """Detrend and search one sector, in two passes.

    Even a robust filter partly absorbs a transit that is a sizeable fraction
    of its window, which makes the dip shallower and shorter and leaves bumps
    on either side. So if pass 1 finds anything, the events (+/- one box
    duration, generous because the box underestimates the true duration) are
    excluded from the trend fit, the trend is interpolated across them, and
    the search is repeated.
    """
    cfg = cfg or SearchConfig()
    mask, flat, trend, res = _run(lc, cfg)
    if res.events:
        exclude = np.zeros(len(lc.time), bool)
        for ev in res.events:
            exclude |= np.abs(lc.time - ev.t0) < ev.duration
        mask, flat, trend, res = _run(lc, cfg, exclude)
    for ev in res.events:
        ev.sector = lc.sector
    return SectorSearch(lc, mask, trend, flat, res)
