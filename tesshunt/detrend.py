"""Robust detrending that preserves long, single transits.

A time-windowed Tukey biweight location (as in Hippke et al. 2019, "wotan") is
evaluated on a coarse grid inside each continuous data segment and linearly
interpolated back onto the cadences. The biweight down-weights points far from
the local location, so a transit shorter than roughly a third of the window is
largely ignored by the trend instead of being fitted away.
"""

from __future__ import annotations

import numpy as np


def split_segments(time: np.ndarray, gap: float = 0.5) -> list[slice]:
    """Split a sorted time array wherever consecutive points are > ``gap`` days apart."""
    breaks = np.where(np.diff(time) > gap)[0] + 1
    edges = np.concatenate([[0], breaks, [len(time)]])
    return [slice(a, b) for a, b in zip(edges[:-1], edges[1:]) if b > a]


def biweight_location(x: np.ndarray, c: float = 5.0, n_iter: int = 5) -> float:
    loc = np.median(x)
    for _ in range(n_iter):
        mad = np.median(np.abs(x - loc))
        if mad == 0:
            return loc
        u = (x - loc) / (c * mad)
        w = (1 - u**2) ** 2
        w[np.abs(u) >= 1] = 0
        new = np.sum(w * x) / np.sum(w)
        if abs(new - loc) < 1e-10:
            return new
        loc = new
    return loc


def biweight_trend(time: np.ndarray, flux: np.ndarray, window: float = 3.0,
                   grid_step: float | None = None, gap: float = 0.5,
                   min_points: int = 20, exclude: np.ndarray | None = None) -> np.ndarray:
    """Return a smooth trend for ``flux`` (same length).

    window     full width of the sliding window, days. Use >= 3x the longest
               transit duration you want to keep.
    grid_step  spacing of trend evaluation points (default window / 30).
    exclude    boolean mask of points to ignore when estimating the trend
               (e.g. a detected transit); the trend is interpolated across them.
    """
    step = grid_step or window / 30.0
    trend = np.full_like(flux, np.nan)
    half = window / 2.0
    for seg in split_segments(time, gap):
        t, f = time[seg], flux[seg]
        if exclude is not None:
            use = ~exclude[seg]
            if use.sum() < min_points:
                use[:] = True
            t_fit, f_fit = t[use], f[use]
        else:
            t_fit, f_fit = t, f
        n_grid = max(2, int(np.ceil((t[-1] - t[0]) / step)) + 1)
        grid = np.linspace(t[0], t[-1], n_grid)
        lo = np.searchsorted(t_fit, grid - half)
        hi = np.searchsorted(t_fit, grid + half)
        vals = np.full(n_grid, np.nan)
        for k in range(n_grid):
            if hi[k] - lo[k] >= min_points:
                vals[k] = biweight_location(f_fit[lo[k]:hi[k]])
        ok = np.isfinite(vals)
        if ok.sum() == 0:
            vals[:] = np.median(f_fit)
            ok[:] = True
        trend[seg] = np.interp(t, grid[ok], vals[ok])
    return trend


def detrend(time: np.ndarray, flux: np.ndarray, window: float = 3.0,
            clip_upper: float = 4.0, **kw) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Detrend and remove upward outliers (flares, cosmic rays).

    Returns (mask, detrended_flux, trend) where ``mask`` selects the points
    kept. Downward outliers are never clipped: they might be the transit.
    """
    trend = biweight_trend(time, flux, window=window, **kw)
    resid = flux / trend
    sigma = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    mask = resid < 1 + clip_upper * sigma
    return mask, resid, trend
