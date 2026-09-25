"""Synthetic-data tests (no network): detrending keeps transits, search finds them."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tesshunt.detect import box_search  # noqa: E402
from tesshunt.detrend import biweight_trend, split_segments  # noqa: E402
from tesshunt.injection import trapezoid  # noqa: E402
from tesshunt.lightcurves import SectorLC  # noqa: E402
from tesshunt.pipeline import search_sector  # noqa: E402


def fake_sector(seed=0, sigma=500e-6, gap=(13.0, 14.0)):
    rng = np.random.default_rng(seed)
    t = np.arange(0, 27, 2 / 1440)
    t = t[(t < gap[0]) | (t > gap[1])]
    # Slow (15 d) variability: a 3 d biweight window tracks it well. Much faster
    # variability needs a shorter window (see RESULTS.md, limitations).
    variability = 1 + 1e-3 * np.sin(2 * np.pi * t / 15.0)
    f = variability * (1 + sigma * rng.standard_normal(len(t)))
    return SectorLC(tic=1, sector=1, time=t + 2000, flux=f, flux_err=np.full(len(t), sigma))


def test_split_segments():
    t = np.array([0.0, 0.1, 0.2, 10.0, 10.1])
    assert [(s.start, s.stop) for s in split_segments(t, gap=0.5)] == [(0, 3), (3, 5)]


def test_trapezoid_depth_and_width():
    t = np.linspace(-1, 1, 20001)
    m = trapezoid(t, 0.0, 1e-3, 0.5, ingress_frac=0.1)
    assert np.isclose(m.min(), 1 - 1e-3)
    assert np.all(m[np.abs(t) > 0.25] == 1)


def test_trend_ignores_transit():
    # 12 h transit, 3 d window: most of the depth survives a single pass...
    lc = fake_sector(seed=3, sigma=200e-6)
    lc.flux = lc.flux * trapezoid(lc.time, 2006.0, 5e-3, 0.5)
    in_tr = np.abs(lc.time - 2006.0) < 0.2
    trend = biweight_trend(lc.time, lc.flux, window=3.0)
    assert 1 - (lc.flux / trend)[in_tr].mean() > 0.8 * 5e-3
    # ...and essentially all of it after the pipeline's masked second pass.
    ss = search_sector(lc)
    assert 1 - ss.detrended[in_tr].mean() > 0.95 * 5e-3


def test_detects_injected_single_transit():
    lc = fake_sector(seed=1)
    lc.flux = lc.flux * trapezoid(lc.time, 2020.0, 2e-3, 8 / 24)
    ss = search_sector(lc)
    assert ss.search.events, "no events"
    ev = ss.search.events[0]
    assert abs(ev.t0 - 2020.0) < 4 / 24
    assert 1e-3 < ev.depth < 3e-3


def test_quiet_noise_has_no_strong_events():
    for seed in range(3):
        ss = search_sector(fake_sector(seed=10 + seed))
        assert not ss.search.events


def test_no_events_flush_against_gap():
    # A dip right at the start of the data has no pre-transit baseline.
    lc = fake_sector(seed=2)
    lc.flux = lc.flux * trapezoid(lc.time, 2000.05, 3e-3, 2 / 24)
    res = box_search(lc.time, lc.flux / biweight_trend(lc.time, lc.flux))
    assert all(abs(e.t0 - 2000.05) > 0.2 for e in res.events)
