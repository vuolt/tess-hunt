"""Unit tests for the Phase 3 vetting checks (synthetic data, no network)."""

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tesshunt import vetting as v  # noqa: E402

CAD = 10 / 1440


def window(seed=0, sigma=500e-6, half=1.0):
    rng = np.random.default_rng(seed)
    t = np.arange(-half, half, CAD) + 2620.0
    return t, 1 + sigma * rng.standard_normal(len(t)), rng


# ---------------------------------------------------------------- 1 shape

def test_shape_transit_wins_for_a_transit():
    t, f, _ = window()
    f = f * v.transit_model(t, 2620.0, 0.1, 6 / 24, 0.3, exp_time=CAD)
    r = v.shape_test(t, f, 2620.0, 6 / 24, exp_time=CAD)
    assert r["passed"], r
    assert abs(r["transit"]["t14_h"] - 6) < 1
    assert abs(r["transit"]["depth_ppm"] - 1e4) < 2e3


@pytest.mark.parametrize("kind", ["step", "ramp", "flare_decay"])
def test_shape_rejects_systematics(kind):
    t, f, _ = window(seed=1)
    x = t - 2620.0
    if kind == "step":
        f = f - 3e-3 * (x > 0)
    elif kind == "ramp":
        f = f - 4e-3 * (1 - (x / 1.0) ** 2)            # broad curved depression
    else:
        f = f - 5e-3 * np.where(x > 0, np.exp(-x / 0.1), 0)
    r = v.shape_test(t, f, 2620.0, 4 / 24, exp_time=CAD)
    assert not r["passed"], r
    assert r["best_alt"] == kind or r["dbic_alt"] < v.SHAPE_DBIC_ALT


def test_shape_rejects_sharp_box_when_resolved():
    # A deep, instantaneous box (e.g. a data-exclusion artefact) at high SNR
    t, f, _ = window(seed=2, sigma=100e-6)
    f = f - 5e-3 * (np.abs(t - 2620.0013) < 0.1234)
    r = v.shape_test(t, f, 2620.0, 6 / 24, exp_time=CAD)
    assert r["dbic_box"] < v.SHAPE_DBIC_BOX and not r["passed"]


# ---------------------------------------------------------------- 2 duration

def test_central_duration_sun_20d():
    # Sun-like star, P = 20 d: T14 ~ 4.9 h
    assert 4.5 < v.central_duration_h(20, 1.0, 1.0) < 5.3
    assert np.isclose(v.period_from_duration(v.central_duration_h(50, 1, 1), 1, 1), 50, rtol=1e-3)


def test_duration_check():
    assert v.duration_check(5.0, 0.2, 0.05, 1.0, 1.0, False)["passed"]
    assert not v.duration_check(1.0, 0.2, 0.05, 1.0, 1.0, False)["passed"]
    g = v.duration_check(1.0, 0.95, 0.05, 1.0, 1.0, True)
    assert g["passed"] and g["note"] == "grazing"
    assert v.duration_check(1.2, 0.2, 0.05, 0.25, 0.25, False)["passed"]   # M dwarf: T_c ~ 2 h


# ---------------------------------------------------------------- 3 edge

def test_edge_check():
    t = np.concatenate([np.arange(0, 5, CAD), np.arange(6, 10, CAD)])
    f = 1 + 1e-4 * np.random.default_rng(1).standard_normal(len(t))
    assert v.edge_check(t, f, 2.5, 4)["passed"]                          # mid-segment
    r = v.edge_check(t, f, 4.95, 2)                                     # runs into the gap
    assert r["near_edge"] and not r["passed"]
    r = v.edge_check(t, f, 4.85, 2)                                     # 2.6 h from edge, resolved
    assert r["near_edge"] and r["resolved"] and r["passed"]
    ramp = f - 0.01 * np.clip((t - 4.7) / 0.3, 0, None) * (t < 5)       # ramp into gap
    assert not v.edge_check(t, ramp, 4.85, 2, depth=0.005)["passed"]
    assert not v.edge_check(t, f, 5.5, 2)["passed"]                     # inside the gap


# ---------------------------------------------------------------- 4 pixels

def _psf(nx, x, y, flux, sig=0.75):
    yy, xx = np.mgrid[0:nx, 0:nx]
    return flux * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sig ** 2))


def _cube(source_xy, others, depth=0.01, seed=0, n=600):
    rng = np.random.default_rng(seed)
    t = np.arange(n) * CAD
    t0, t14 = t[n // 2], 4 / 24
    frames = []
    for ti in t:
        img = np.zeros((15, 15))
        for (x, y, fl) in [source_xy] + others:
            dip = depth if ((x, y, fl) == source_xy and abs(ti - t0) < t14 / 2) else 0
            img += _psf(15, x, y, fl * (1 - dip))
        frames.append(img + rng.normal(0, 5, img.shape))
    return t, np.array(frames), t0, t14


def test_centroid_on_target():
    tgt = (7.0, 7.0, 5000.0)
    t, cube, t0, t14 = _cube(tgt, [(9.0, 7.0, 3000.0)])
    diff, noise, oot, *_ = v.difference_image(t, cube, t0, t14)
    r = v.centroid_test(diff, noise, 7.0, 7.0)
    assert r["passed"] and r["offset_px"] < 0.5


def test_centroid_and_neighbour_catch_offset_source():
    nb = (9.0, 7.0, 3000.0)
    t, cube, t0, t14 = _cube(nb, [(7.0, 7.0, 5000.0)], depth=0.03)
    diff, noise, oot, *_ = v.difference_image(t, cube, t0, t14)
    c = v.centroid_test(diff, noise, 7.0, 7.0)
    assert not c["passed"] and c["offset_px"] > 1
    n = v.neighbour_test(diff, noise, oot, 7.0, 7.0, [(123, 9.0, 7.0, 0.5)], 0.01)
    assert not n["passed"] and n["failing"] == [123]


def test_neighbour_passes_when_target_is_source():
    tgt = (7.0, 7.0, 5000.0)
    t, cube, t0, t14 = _cube(tgt, [(9.0, 7.0, 3000.0)])
    diff, noise, oot, *_ = v.difference_image(t, cube, t0, t14)
    n = v.neighbour_test(diff, noise, oot, 7.0, 7.0, [(123, 9.0, 7.0, 0.5)], 0.01)
    assert n["passed"]
    n = v.neighbour_test(diff, noise, oot, 7.0, 7.0, [(5, 7.2, 7.1, 1.0)], 0.01)
    assert n["unresolved"] == [5]


# ---------------------------------------------------------------- 5 asteroids

SKYBOT = """# Flag: 1
# Num | Name | RA(h) | DE(deg) | Class | Mv | Err(arcsec) | d(arcsec) | dRA(arcsec/h) | dDEC(arcsec/h) | Dg(ua)
 - | 2013 RL195 | 10 00 00.00 | +15 00 30.0 | MB>Cybele | 17.5 | 0.2 | 30.0 | -20.0 | 0.0 | 2.7
 459640 | 2013 JE32 | 10 00 30.00 | +15 10 00.0 | MB>Inner | 16.0 | 0.2 | 700.0 | 5.0 | 5.0 | 1.4
"""


def test_parse_skybot_and_asteroid_check():
    objs = v.parse_skybot(SKYBOT)
    assert len(objs) == 2 and objs[0]["v"] == 17.5
    r = v.asteroid_check(objs, 150.0, 15.0, 4.0)
    assert not r["passed"] and r["hits"][0]["name"] == "2013 RL195"
    faint = [dict(o, v=21.0) for o in objs]
    assert v.asteroid_check(faint, 150.0, 15.0, 4.0)["passed"]


# ---------------------------------------------------------------- 6 catalogues

def test_catalogue_check():
    cats = dict(
        exofop_toi=pd.DataFrame({"TIC ID": [1, 2], "TOI": [100.01, 200.01],
                                 "TFOPWG Disposition": ["PC", "FP"]}),
        exofop_ctoi=pd.DataFrame({"TIC ID": [3], "CTOI": ["3.01"]}),
        tess_ebs=pd.DataFrame({"tess_id": [4]}), villanova_ebs={6})
    assert v.catalogue_check(1, cats)["passed"] and v.catalogue_check(1, cats)["known"]
    assert not v.catalogue_check(2, cats)["passed"]                    # TOI marked FP
    assert v.catalogue_check(3, cats)["known"]
    assert not v.catalogue_check(4, cats)["passed"]                    # known EB
    assert not v.catalogue_check(6, cats)["passed"]                    # Villanova-only EB
    r = v.catalogue_check(5, cats)
    assert r["passed"] and not r["known"]


# ---------------------------------------------------------------- 7 FPP

def test_fpp_check_thresholds():
    assert v.fpp_check(dict(fpp=0.1, nfpp=0.01))["passed"]
    assert not v.fpp_check(dict(fpp=0.7, nfpp=0.01))["passed"]
    assert not v.fpp_check(dict(fpp=0.1, nfpp=0.3))["passed"]
