"""Planet Hunters TESS-style plot (no network)."""

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from tesshunt.phtplot import _binned, pht_plot  # noqa: E402


def test_plot_is_written(tmp_path):
    rng = np.random.default_rng(0)
    t = 2600 + np.arange(0, 27, 10 / 1440)
    f = 1 + 1e-3 * rng.standard_normal(len(t))
    f[np.abs(t - 2610) < 0.1] -= 3e-3
    p = pht_plot(t, f, 2610.0, 0.2, 3e-3, "TIC 1", str(tmp_path / "a" / "tic1_pht.png"), sector=48)
    assert os.path.getsize(p) > 20_000


def test_binned_medians():
    t = np.arange(0, 1, 0.01)
    tb, fb = _binned(t, np.ones_like(t), 0.1)
    assert len(tb) == 10 and np.allclose(fb, 1)


def test_app_finds_the_plot():
    import pht_plots
    from app import data
    p = pht_plots.plot_path(48, 95747180)
    if os.path.exists(p):
        assert data.sheets(48, 95747180)["pht"] == p
    assert pht_plots.plot_path(21, 1).endswith(os.path.join("phase4", "s0021", "pht", "tic1_pht.png"))
