"""Transit prediction helpers."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tesshunt import ephemeris as ep  # noqa: E402


def test_predict_error_propagation():
    pr = ep.predict(0.0, 0.01, 80.0, 0.01, 1, 100.0, 400.0)
    assert np.allclose(pr["tc"], [160.0, 240.0, 320.0, 400.0])
    k = pr["k"]
    assert np.allclose(pr["sigma"], np.sqrt((1 + k) ** 2 + k ** 2) * 0.01)
    half = ep.predict(0.0, 0.01, 80.0, 0.01, 2, 100.0, 400.0)
    assert np.isclose(half["P"], 40.0) and len(half["tc"]) == 8
    # the n = 2 ephemeris contains every n = 1 transit, with the same error
    common = np.isin(half["tc"], pr["tc"])
    assert np.allclose(half["sigma"][common], pr["sigma"])


def test_sector_bounds_are_contiguous():
    b = ep.sector_bounds()
    assert b[48][0] < 2621.1 < b[48][1]          # the S48 transit lies inside S48
    assert abs(b[48][1] - b[49][0]) < 1e-6
    # most sectors are ~27 d; S97-S98 have mid-times 40-56 d apart in tess-point
    assert all(20 < e - s < 60 for s, e in b.values())
