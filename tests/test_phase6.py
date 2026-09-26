"""Phase 6: false-positive triage logic (no network)."""

import os
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from tesshunt import fptriage as fp  # noqa: E402
from tesshunt.injection import trapezoid  # noqa: E402


def test_period_match_factors_and_tolerance():
    assert fp.period_match(10.0, 1e-4, 10.005, 0.01)[0] == 1.0
    assert fp.period_match(5.0, 1e-4, 10.0, 0.01)[0] == 0.5        # EB found at half period
    assert fp.period_match(20.0, 1e-4, 10.0, 0.01)[0] == 2.0
    assert fp.period_match(10.0, 1e-4, 11.0, 0.01) is None
    assert fp.period_match(np.nan, 1e-4, 10.0, 0.01) is None


def test_companion_mass_sb1_wasp18():
    """WASP-18 b: K = 1.8 km/s, P = 0.941 d, M1 = 1.2 M_sun -> ~10 M_Jup (sin i ~ 1)."""
    m2, how = fp.companion_mass(dict(nss_solution_type="SB1", period=0.9415, eccentricity=0.0,
                                     semi_amplitude_primary=1.82), 1.2)
    assert 9 < m2 * 1047.6 < 12 and "SB1" in how


def test_companion_mass_sb2_uses_k_ratio():
    m2, how = fp.companion_mass(dict(nss_solution_type="SB2", period=5.0, semi_amplitude_primary=40.0,
                                     semi_amplitude_secondary=50.0, mass_ratio=np.nan), 1.0)
    assert m2 == pytest.approx(0.8) and "SB2" in how


def test_thiele_innes_face_on_circle():
    # face-on circular orbit of radius a: A = G = a, B = F = 0
    assert fp.thiele_innes_a0(2.0, 0.0, 0.0, 2.0) == pytest.approx(2.0)


def test_solve_m2_roundtrip():
    m1, m2 = 1.0, 0.3
    f = m2 ** 3 / (m1 + m2) ** 2
    assert fp._solve_m2(f, m1) == pytest.approx(m2, rel=1e-6)


def _eb(P=3.0, t0=1.0, dur=0.12, d1=0.02, d2=0.01, odd_even=0.0, days=27, cad=10 / 1440, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(0, days, cad)
    f = 1 + 3e-4 * rng.standard_normal(len(t))
    for n in range(int(days / P) + 1):
        dd = d1 * (1 + (odd_even if n % 2 else -odd_even))
        f = f * trapezoid(t, t0 + n * P, dd, dur)
        if d2:
            f = f * trapezoid(t, t0 + (n + 0.5) * P, d2, dur)
    return t, f


def test_secondary_scan_finds_the_secondary():
    t, f = _eb()
    keep = ~fp._transit_mask(t, 3.0, 1.0, 0.12, pad=1.5)
    ph, depth, snr, n = fp.secondary_scan(t[keep], f[keep], 3.0, 1.0, 0.12)
    assert abs(ph - 0.5) < 0.02 and depth == pytest.approx(0.01, rel=0.3) and snr > 20
    t2, f2 = _eb(d2=0.0)                                   # no secondary
    keep = ~fp._transit_mask(t2, 3.0, 1.0, 0.12, pad=1.5)
    assert fp.secondary_scan(t2[keep], f2[keep], 3.0, 1.0, 0.12)[2] < 5


def test_odd_even_detected_in_events():
    t, f = _eb(d2=0.0, odd_even=0.25)
    ev, _, _ = fp.sector_measure(t, f, np.zeros_like(t), np.zeros_like(t), 3.0, 1.0, 0.12, 0.02, [])
    r = fp.combine(ev)
    assert r["oddeven_sigma"] > 10 and r["oddeven_frac"] > 0.3
    t, f = _eb(d2=0.0, odd_even=0.0, seed=2)
    ev, _, _ = fp.sector_measure(t, f, np.zeros_like(t), np.zeros_like(t), 3.0, 1.0, 0.12, 0.02, [])
    assert fp.combine(ev)["oddeven_sigma"] < 4


def test_other_planets_are_masked_from_secondary_scan():
    """A second planet on the star must not look like a secondary eclipse."""
    t, f = _eb(d2=0.0)
    f = f * np.prod([trapezoid(t, 2.2 + k * 7.0, 0.01, 0.12) for k in range(4)], axis=0)
    others = [(7.0, 2.2, 0.12)]
    _, tt, ff = fp.sector_measure(t, f, np.zeros_like(t), np.zeros_like(t), 3.0, 1.0, 0.12, 0.02, others)
    assert fp.secondary_scan(tt, ff, 3.0, 1.0, 0.12)[2] < 5


def test_gaia_flag_rules():
    from phase6_report import gaia_flags
    m = pd.DataFrame([
        dict(name="A", tic=1, group="unresolved", source="Gaia NSS SB1", gaia_period=5.0, factor=1.0,
             n_sigma=0.1, m2=0.3, m2_method="SB1", eclipsing_solution=False),
        dict(name="B", tic=2, group="planet", source="Gaia NSS SB1", gaia_period=0.94, factor=1.0,
             n_sigma=0.1, m2=0.010, m2_method="SB1", eclipsing_solution=False),       # WASP-18 b-like
        dict(name="C", tic=3, group="unresolved", source="Gaia eclipsing binary (photometric)",
             gaia_period=2.0, factor=0.5, n_sigma=0.3, m2=np.nan, m2_method="", eclipsing_solution=True),
        dict(name="D", tic=4, group="unresolved", source="Gaia NSS SB1", gaia_period=9.0, factor=1.0,
             n_sigma=0.1, m2=0.09, m2_method="SB1", eclipsing_solution=False),
    ])
    g = gaia_flags(m).set_index("name")
    assert g.gaia_flag["A"] and g.gaia_conf["A"] == "high"
    assert not g.gaia_flag["B"] and g.gaia_substellar["B"]
    assert g.gaia_flag["C"] and g.gaia_conf["C"] == "high"
    assert g.gaia_flag["D"] and g.gaia_conf["D"] == "medium"


def test_variant_choice_rejects_checks_that_flag_planets():
    from phase6_report import choose_variants
    import phase6 as p6
    rows = []
    for g, n_bad in (("planet", 20), ("fp", 60), ("unresolved", 10)):
        for i in range(100):
            bad = i < n_bad
            rows.append(dict(sample="random", group=g, n_events=5,
                             oddeven_sigma=20.0 if bad else 1.0, oddeven_frac=0.5 if (bad and g == "fp") else 0.05,
                             sec_snr=0.0, sec_ratio=0.0, centroid_sigma=0.0, source_offset_px=0.0))
    var, chosen = choose_variants(pd.DataFrame(rows))
    # plain sigma cuts flag 20 % of planets; the ">= 20 % different" variant flags none
    assert chosen["odd_even"][1] == p6.VARIANTS["odd_even"][2][0]
    assert chosen["secondary"] is not None and chosen["centroid"] is not None
