"""Phase 2 unit tests: variability-adaptive window, FFI URLs, resumable store (no network)."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tesshunt import ffi  # noqa: E402
from tesshunt.injection import trapezoid  # noqa: E402
from tesshunt.lightcurves import SectorLC  # noqa: E402
from tesshunt.survey import Store, analyse, edge_distance  # noqa: E402
from tesshunt.variability import measure, residual_fraction  # noqa: E402


def ffi_like(seed=0, sigma=1000e-6, amp=0.0, period=5.0):
    rng = np.random.default_rng(seed)
    t = np.arange(0, 26, 10 / 1440)
    t = t[(t < 12.5) | (t > 13.5)]
    f = 1 + amp * np.sin(2 * np.pi * t / period) + sigma * rng.standard_normal(len(t))
    return SectorLC(tic=1, sector=48, time=t + 2600, flux=f, flux_err=np.full(len(t), sigma))


def test_lc_url():
    assert ffi.lc_url(298663873, 48) == (
        "https://archive.stsci.edu/hlsps/tess-spoc/s0048/target/0000/0002/9866/3873/"
        "hlsp_tess-spoc_tess_phot_0000000298663873-s0048_tess_v1_lc.fits")


def test_residual_fraction():
    assert residual_fraction(3.0, 1000.0) < 1e-4       # slow trend: fully followed
    assert residual_fraction(3.0, 3.0) == 1.0           # window = period: nothing followed


def test_quiet_star_keeps_full_window():
    v = measure(*_tf(ffi_like()))
    assert v.window == 3.0 and not v.long_limited and not v.hf_variable


def test_spotted_star_gets_short_window_and_flag():
    v = measure(*_tf(ffi_like(amp=0.01, period=3.0)))
    assert np.isclose(v.period, 3.0, rtol=0.05)
    assert v.window < 3.0 and v.long_limited and not v.hf_variable
    assert v.max_duration_h == v.window * 8


def test_fast_pulsator_flagged_hf():
    v = measure(*_tf(ffi_like(amp=0.01, period=0.1)))
    assert v.hf_variable and v.window == 3.0


def test_transit_does_not_masquerade_as_variability():
    lc = ffi_like(seed=3, sigma=300e-6)
    lc.flux = lc.flux * trapezoid(lc.time, 2606.0, 5e-3, 1.0)
    var, cfg, ss = analyse(lc)
    assert var.window == 3.0
    assert ss.search.events and abs(ss.search.events[0].t0 - 2606.0) < 0.2


def test_edge_distance():
    t = np.concatenate([np.arange(0, 5, 0.01), np.arange(6, 10, 0.01)])
    assert np.isclose(edge_distance(t, 1.0), 24.0)
    assert edge_distance(t, 5.5) < 0          # in the gap


def test_store_is_atomic_and_resumable(tmp_path):
    db = str(tmp_path / "x.sqlite")
    s = Store(db)
    s.write({"star": {"tic": 5, "status": "ok", "n_dips": 1},
             "dips": [{"tic": 5, "rank_in_star": 1, "t0": 1.0, "duration_h": 2.0,
                       "depth_ppm": 100.0, "snr": 8.0, "edge_dist_h": 5.0, "n_in": 12}]})
    # rewriting the same star replaces rather than duplicates its rows
    s.write({"star": {"tic": 5, "status": "ok", "n_dips": 1},
             "dips": [{"tic": 5, "rank_in_star": 1, "t0": 1.0, "duration_h": 2.0,
                       "depth_ppm": 100.0, "snr": 8.0, "edge_dist_h": 5.0, "n_in": 12}]})
    assert Store(db).done() == {5}
    assert s.db.execute("SELECT count(*) FROM dips").fetchone()[0] == 1


def _tf(lc):
    return lc.time, lc.flux
