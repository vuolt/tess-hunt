"""Physical checks and the planet score (no network)."""

import os
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from tesshunt import planetscore as ps  # noqa: E402


def test_a_over_r_earth_sun():
    assert ps.a_over_r(365.25, 1.0) == pytest.approx(215, rel=0.01)


def test_max_duration_hot_jupiter():
    # WASP-18 b-like: P = 0.94 d, a/R* ~ 3.5 -> a central transit lasts ~2.2 h
    ar = 3.5
    assert ps.max_duration_d(0.9415, 0.0095, ar) * 24 == pytest.approx(2.3, abs=0.3)


def test_density_from_mass_or_logg():
    assert ps.stellar_density(1.0, 1.0, np.nan) == pytest.approx(1.0)
    assert ps.stellar_density(np.nan, 1.0, 4.438) == pytest.approx(1.0)
    assert ps.stellar_density(np.nan, 10.0, 2.438) == pytest.approx(0.001)   # a giant


def test_planet_secondary_limit():
    # WASP-18 b: measured TESS secondary ~340 ppm must be allowed
    assert ps.max_planet_secondary(0.0095, 3.5, 6400) > 340e-6
    # a warm Neptune around a Sun-like star cannot make a 1000 ppm secondary
    assert ps.max_planet_secondary(0.001, 30, 5800) < 1e-5


def test_table_checks_flag_a_long_transit():
    t = pd.DataFrame(dict(period=[10.0, 10.0], depth_ppm=[1000.0, 1000.0], duration_h=[4.0, 20.0],
                          mstar=[1.0, 1.0], rstar=[1.0, 1.0], logg=np.nan, rp_re=[3.5, 3.5]))
    q = ps.table_quantities(t)
    assert q.duration_ratio[0] < 1.5 < q.duration_ratio[1]
    assert q.rp_rj[0] == pytest.approx(3.5 / 11.21)


def test_score_learns_simple_signal():
    rng = np.random.default_rng(0)
    n = 400
    y = (rng.random(n) < 0.5).astype(int)
    t = pd.DataFrame(dict(name=[f"T{i}" for i in range(n)], period=10 ** rng.uniform(0, 2, n),
                          depth_ppm=np.where(y, 2000, 20000) * rng.uniform(0.5, 1.5, n),
                          duration_h=3.0, mstar=1.0, rstar=1.0, logg=4.4, teff=5800.0,
                          rp_re=np.where(y, 5, 25) * rng.uniform(0.7, 1.3, n), n_on_star=1))
    f = ps.features(t, None, set())
    model, p = ps.fit_score(f, y)
    assert ps.auc(y, p) > 0.9
