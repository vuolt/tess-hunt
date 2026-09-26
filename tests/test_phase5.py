"""Phase 5 unit tests: the decision logic of the expert checks (no network)."""

import os
import sys

import numpy as np
import pytest
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from tesshunt import expert as ex  # noqa: E402


def test_ecc_range_and_e_min():
    lo, hi = ex.ecc_factor_range(0.5)
    assert np.isclose(lo, (0.5 / np.sqrt(0.75)) ** 3) and np.isclose(hi, (1.5 / np.sqrt(0.75)) ** 3)
    from phase5 import e_min_from_ratio
    assert np.isclose(e_min_from_ratio(1.0), 0.0)
    # rho ratio g^3 with g = sqrt((1+e)/(1-e)) for w = 90 deg -> recovers e
    for e in (0.2, 0.5, 0.8):
        g = np.sqrt((1 + e) / (1 - e))
        assert np.isclose(e_min_from_ratio(g ** 3), e)
        assert np.isclose(e_min_from_ratio(g ** -3), e)


def test_period_consistency():
    P = np.full(1000, 100.0)
    pc = ex.period_consistency(P, [(90, 110)], None, None)
    assert pc["frac_allowed"] == 1.0 and pc["reachable_with_ecc"]
    pc = ex.period_consistency(P, [(5000, 6000)], None, None)
    assert pc["frac_allowed"] == 0.0 and not pc["reachable_with_ecc"]
    pc = ex.period_consistency(P, [(300, 400)], None, None)       # needs e ~0.5
    assert pc["frac_allowed"] == 0.0 and pc["reachable_with_ecc"]
    pc = ex.period_consistency(P, [], 80.0, 900)                   # above the floor
    assert pc["frac_allowed"] == 1.0


def test_stellar_check_classes():
    dwarf = ex.stellar_check(dict(evolstage_flame=250, radius_flame=1.0), 1.1, 1.0, 0.01)
    assert dwarf["cls"] == "dwarf" and not dwarf["flag"]
    assert np.isclose(dwarf["rp_rj"], 0.1 * 1.0 * ex.RSUN_RJUP)
    sg = ex.stellar_check(dict(evolstage_flame=450, radius_flame=2.6), 1.2, 1.2, 0.01)
    assert sg["cls"] == "subgiant" and sg["flag"] and sg["serious"]        # 2.5 R_J
    giant = ex.stellar_check(dict(logg_gspphot=2.8, radius_gspphot=8.0), 8.0, 1.5, 1e-4)
    assert giant["cls"] == "giant" and not giant["serious"]                 # small companion
    tic_only = ex.stellar_check(None, 1.0, 1.0, 0.01)             # falls back to the TIC radius
    assert tic_only["cls"] == "dwarf" and tic_only["radius_source"] == "TIC"
    assert "no evolutionary" in tic_only["basis"]
    nothing = ex.stellar_check(None, None, None, 0.01)
    assert nothing["cls"] == "unknown"


def test_binarity_extras_thresholds():
    ok = ex.binarity_extras(dict(ruwe=1.0, ipd_frac_multi_peak=0, ipd_gof_harmonic_amplitude=0.02,
                                 rv_nb_transits=20, rv_chisq_pvalue=0.5, rv_renormalised_gof=1.0,
                                 radial_velocity_error=1.0))
    assert not ok["serious"] and ok["flags"] == []
    rv = ex.binarity_extras(dict(ruwe=1.0, rv_nb_transits=20, rv_chisq_pvalue=1e-5, rv_renormalised_gof=8,
                                 rv_amplitude_robust=30))
    assert rv["serious"] and "rv_variable" in rv["flags"]
    # a resolved/partly resolved companion only dilutes (TOI-3837 b is a confirmed planet
    # with ipd_frac_multi_peak = 42): flagged, but not serious
    mp = ex.binarity_extras(dict(ruwe=1.1, ipd_frac_multi_peak=25))
    assert "ipd_multi_peak" in mp["flags"] and not mp["serious"]
    ruwe_only = ex.binarity_extras(dict(ruwe=2.0, ipd_gof_harmonic_amplitude=0.05))
    assert not ruwe_only["serious"] and ruwe_only["flags"] == ["ruwe"]
    harm = ex.binarity_extras(dict(ruwe=2.0, ipd_gof_harmonic_amplitude=0.2))
    assert "ipd_harmonic" in harm["flags"] and not harm["serious"]


def test_variability_check():
    quiet = ex.variability_check(dict(phot_variable_flag="NOT_AVAILABLE"), pd.DataFrame())
    assert quiet["flags"] == [] and not quiet["serious"]
    rot = ex.variability_check(dict(phot_variable_flag="VARIABLE", best_class_name="SOLAR_LIKE"), None)
    assert "gaia_variable" in rot["flags"] and not rot["serious"]
    eb = ex.variability_check(None, pd.DataFrame(dict(Name=["X Cnc"], Type=["EA"], Period=["2.1"], _r=[3.0])))
    assert eb["serious"]
    eb2 = ex.variability_check(None, pd.DataFrame(dict(Name=["Y"], Type=["EA/RS"], Period=[""], _r=[1.0])))
    assert eb2["serious"]
    # VSX type EP = a known transiting planet (e.g. HD 191939): not variability
    ep = ex.variability_check(None, pd.DataFrame(dict(Name=["HD 191939"], Type=["EP"], Period=[""], _r=[0.5])))
    assert ep["flags"] == [] and not ep["serious"] and "transiting planet" in ep["notes"][0]
    rotv = ex.variability_check(None, pd.DataFrame(dict(Name=["Z"], Type=["ROT"], Period=[""], _r=[0.5])))
    assert rotv["flags"] == ["vsx"] and not rotv["serious"]


def test_psf_image_normalised():
    tot, ims = ex._psf_image((21, 21), [(10.2, 9.7, 100.0)], 0.8, 0.3, 3.0)
    assert abs(ims[0].sum() - 100.0) < 1.0          # flux conserved on a large stamp


def _t(**kw):
    t = dict(role="candidate", category="maybe", reasons="S1 FPP 0.12 >= 0.1", fpp=0.12)
    t.update(kw)
    return t


def test_rerank_rules():
    from phase5 import rerank
    strong = dict(verdict="strong", serious=[], minor=[])
    doubt = dict(verdict="doubtful", serious=["aperture"], minor=[])
    plaus = dict(verdict="plausible", serious=[], minor=["ruwe"])
    assert rerank(_t(), strong)[0] == "submit"                             # FPP-only maybe promoted
    assert rerank(_t(fpp=0.35, reasons="S1 FPP 0.35 >= 0.1"), strong)[0] == "maybe"
    assert rerank(_t(reasons="S2 host binarity: RUWE 2"), strong)[0] == "maybe"
    assert rerank(_t(), doubt)[0] == "drop"
    assert rerank(_t(category="submit", reasons=""), doubt)[0] == "maybe"
    assert rerank(_t(category="submit", reasons=""), plaus)[0] == "submit"
    assert rerank(_t(role="validation"), strong)[0] is None


def test_matching_period_picks_the_right_toi():
    from phase5 import _matching_period
    tois = pd.DataFrame(dict(tid=[1, 1, 1, 2], toi=["1.01", "1.02", "1.03", "2.01"],
                             pl_orbper=[8.88, 28.58, 38.35, float("nan")],
                             pl_tranmid=[2458715.36, 2458726.05, 2458743.55, 2458700.0]))
    t0 = 2458726.05 - 2457000 + 20 * 28.58          # a transit of the 28.58 d planet only
    assert _matching_period(tois, 1, t0) == 28.58
    assert _matching_period(tois, 1, t0 + 5.0) is None
    assert _matching_period(tois, 2, 1700.0) is None  # no period


def test_needed_density_scales_linearly_with_period():
    """Fixed duration and impact parameter: a/R* ~ P, so rho = 3 pi (a/R*)^3 / (G P^2) ~ P
    (this is how phase5 scales one fit to every alias of a duotransit)."""
    t14, b, rp = 4.6 / 24, 0.5, 0.086

    def rho_needed(P):
        s = np.sin(np.pi * t14 / P)
        a = np.sqrt(b ** 2 + ((1 + rp) ** 2 - b ** 2) / s ** 2)
        return 3 * np.pi * a ** 3 / (ex.G_CGS * (P * 86400) ** 2)
    assert rho_needed(80.48) / rho_needed(40.24) == pytest.approx(2.0, rel=1e-3)


def test_other_transits_are_cut_out():
    """Another planet's transit near the candidate must be removed before fitting."""
    t = np.arange(0, 10, 0.01)
    keep = ex._keep_mask(t, [(3.0, 0.1)])
    assert not keep[np.abs(t - 3.0) < 0.1].any()
    assert keep[np.abs(t - 3.0) > 0.2].all()
    assert ex._keep_mask(t, None).all()
