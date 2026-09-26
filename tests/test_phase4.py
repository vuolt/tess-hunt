"""Phase 4 unit tests: polite networking, period constraints, binarity, ranking."""

import http.server
import os
import sys
import threading
import time

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from tesshunt import binarity as bn  # noqa: E402
from tesshunt import multisector as ms  # noqa: E402
from tesshunt import net  # noqa: E402
from tesshunt.injection import trapezoid  # noqa: E402


# ---------------------------------------------------------------- net

class _Handler(http.server.BaseHTTPRequestHandler):
    hits = {}

    def do_GET(self):  # noqa: N802
        _Handler.hits[self.path] = _Handler.hits.get(self.path, 0) + 1
        if self.path.startswith("/missing"):
            self.send_response(404)
            self.end_headers()
            return
        if self.path.startswith("/flaky") and _Handler.hits[self.path] < 3:
            self.send_response(503)
            self.end_headers()
            return
        self.send_response(200)
        self.end_headers()
        self.wfile.write(b"hello " + self.path.encode())

    def log_message(self, *a):
        pass


@pytest.fixture()
def server(tmp_path, monkeypatch):
    monkeypatch.setattr(net, "CACHE", str(tmp_path / "cache"))
    monkeypatch.setitem(net.SERVICES, "test", 0.2)
    real_sleep = time.sleep
    monkeypatch.setattr(net.time, "sleep", lambda s: None if s > 1 else real_sleep(s))
    _Handler.hits = {}
    srv = http.server.HTTPServer(("127.0.0.1", 0), _Handler)
    th = threading.Thread(target=srv.serve_forever, daemon=True)
    th.start()
    yield f"http://127.0.0.1:{srv.server_port}"
    srv.shutdown()


def test_net_caches_responses(server):
    assert net.get(f"{server}/a", "test") == b"hello /a"
    assert net.get(f"{server}/a", "test") == b"hello /a"
    assert _Handler.hits["/a"] == 1                    # second call served from disk


def test_net_caches_404(server):
    for _ in range(2):
        with pytest.raises(FileNotFoundError):
            net.get(f"{server}/missing", "test")
    assert _Handler.hits["/missing"] == 1


def test_net_backs_off_on_503(server):
    assert net.get(f"{server}/flaky", "test") == b"hello /flaky"
    assert _Handler.hits["/flaky"] == 3


def test_service_slot_spacing(tmp_path, monkeypatch):
    monkeypatch.setattr(net, "CACHE", str(tmp_path))
    monkeypatch.setitem(net.SERVICES, "spaced", 0.15)
    ts = []
    for _ in range(4):
        with net.service_slot("spaced"):
            ts.append(time.time())
    assert np.diff(ts).min() >= 0.14


def test_cached_call(tmp_path, monkeypatch):
    monkeypatch.setattr(net, "CACHE", str(tmp_path))
    calls = []
    f = lambda: calls.append(1) or {"x": 1}   # noqa: E731
    assert net.cached_call("test", ("k",), f) == {"x": 1}
    assert net.cached_call("test", ("k",), f) == {"x": 1}
    assert len(calls) == 1


# ---------------------------------------------------------------- periods

def _sector(t_start, days=25, cad=10 / 1440, sigma=3e-4, seed=0, transits=(), depth=5e-3,
            t14=0.25):
    rng = np.random.default_rng(seed)
    t = np.arange(t_start, t_start + days, cad)
    f = 1 + sigma * rng.standard_normal(len(t))
    for tc in transits:
        f = f * trapezoid(t, tc, depth, t14)
    return t, f


def test_ruled_out_map_marks_flat_data():
    t, f = _sector(0)
    grid, status, d = ms.ruled_out_map(t, f, 0.25, 5e-3)
    assert (status == 2).mean() > 0.95                 # a 5 ppt dip would be obvious anywhere


def test_ruled_out_map_keeps_real_transit_allowed():
    t, f = _sector(0, transits=(10.0,))
    grid, status, d = ms.ruled_out_map(t, f, 0.25, 5e-3)
    near = np.abs(grid - 10.0) < 0.05
    assert (status[near] == 1).all()


def test_period_scan_recovers_true_period():
    t0, P = 100.0, 37.3
    maps = []
    for k, start in enumerate([90.0, 200.0, 310.0]):
        tr = [tc for tc in t0 + P * np.arange(-3, 10) if start <= tc <= start + 25]
        t, f = _sector(start, transits=tr, seed=k)
        grid, status, _ = ms.ruled_out_map(t, f, 0.25, 5e-3)
        maps.append((grid, status))
    scan = ms.period_scan(t0, 0.25, maps, pmin=5.0)
    i = np.argmin(np.abs(scan["P"] - P))
    assert scan["allowed"][i - 2:i + 3].any()
    assert scan["allowed"].mean() < 0.5                # most periods are excluded
    assert all(lo <= hi for lo, hi in ms.allowed_intervals(scan))


def test_alias_periods():
    a = ms.alias_periods(0.0, 90.0, pmin=10)
    assert np.allclose(a[:3], [90, 45, 30]) and a.min() >= 10


def test_matched_dips_flags_consistency():
    t, f = _sector(0, transits=(8.0, 17.0), depth=5e-3)
    dips, sig = ms.matched_dips(t, f, 0.25, 5e-3)
    found = [d for d in dips if d["consistent"]]
    assert len(found) >= 2 and np.isfinite(sig)


# ---------------------------------------------------------------- binarity

def _target(**kw):
    base = dict(source_id=1, ruwe=1.0, non_single_star=0, ipd_frac_multi_peak=0,
                rv_chisq_pvalue=np.nan, rv_nb_transits=0, parallax=10.0, parallax_error=0.05,
                pmra=50.0, pmdec=-20.0, phot_g_mean_mag=11.0)
    base.update(kw)
    return base


def test_binarity_clean_star():
    r = bn.assess(_target(), [], pd.DataFrame(), pd.DataFrame(), pd.DataFrame())
    assert not r["host_binary"] and r["flags"] == []


def test_binarity_ruwe_and_close_wds():
    r = bn.assess(_target(ruwe=2.1), [], pd.DataFrame(), pd.DataFrame(), pd.DataFrame())
    assert r["host_binary"] and "ruwe" in r["flags"]
    wds = pd.DataFrame(dict(WDS=["x"], Comp=[""], sep2=["0.4"], mag1=["11.3"], mag2=["11.4"]))
    r = bn.assess(_target(), [], pd.DataFrame(), pd.DataFrame(), wds)
    assert r["host_binary"] and "wds_close" in r["flags"]


def test_binarity_comoving_is_informational():
    nb = pd.DataFrame(dict(source_id=[1, 2], phot_g_mean_mag=[11.0, 14.0], parallax=[10.0, 10.02],
                           parallax_error=[0.05, 0.05], pmra=[50.0, 50.5], pmdec=[-20.0, -20.3],
                           ruwe=[1.0, 1.0], sep_arcsec=[0.0, 30.0]))
    r = bn.assess(_target(), [], nb, pd.DataFrame(), pd.DataFrame())
    assert "comoving_companion" in r["flags"] and not r["host_binary"]


# ---------------------------------------------------------------- ranking

def _g(**kw):
    g = dict(tic=1, role="candidate", review_notes="", fpp=0.01, rp_rj=1.0, rp_rj_diluted=1.0,
             binarity=dict(flags=[], host_binary=False, summary="none", ruwe=1.0),
             confirm=[dict(kind="qlp", available=True, expected_snr=20, depth_ratio=0.9, snr=18)],
             duos=[], joint_periods=None, recurrent_failing=[], other_dips=[])
    g.update(kw)
    return g


def test_rank_rules():
    from phase4_report import classify
    assert classify(_g())[0] == "submit"
    assert classify(_g(confirm=[dict(kind="qlp", available=True, expected_snr=20,
                                     depth_ratio=0.1, snr=1)]))[0] == "drop"          # D1
    assert classify(_g(review_notes="very long dip: x; grazing/V-shaped: y"))[0] == "drop"  # D2
    bin_host = dict(flags=["ruwe"], host_binary=True, summary="RUWE 2", ruwe=2.0)
    assert classify(_g(binarity=bin_host, rp_rj_diluted=2.5))[0] == "drop"          # D3
    assert classify(_g(binarity=bin_host))[0] == "maybe"                             # S2
    assert classify(_g(fpp=0.3))[0] == "maybe"                                       # S1
    assert classify(_g(rp_rj_diluted=1.9))[0] == "maybe"                             # S3
    assert classify(_g(review_notes="difference image too faint for a centroid"))[0] == "submit"


# ---------------------------------------------------------------- step test

def test_one_sided_check_transit_on_slope_is_two_sided():
    from tesshunt import vetting as v
    t, f = _sector(0, days=4, sigma=5e-4, transits=(2.0,), depth=5e-3, t14=0.25)
    f = f * (1 + 4e-3 * (t - 2.0))                    # steep linear trend
    r = v.one_sided_check(t, f, 2.0, 0.25, 5e-3)
    assert r["one_sided"] is False and r["pre"] > 0.75 and r["post"] > 0.75


def test_one_sided_check_flags_step():
    from tesshunt import vetting as v
    t, f = _sector(0, days=4, sigma=5e-4)
    f = np.where(t > 2.0 + 0.125, f * 1.01, f)         # flux jumps up at the "egress"
    r = v.one_sided_check(t, f, 2.0, 0.25, 5e-3)
    assert r["one_sided"] is True and abs(r["pre"]) < 0.25 and r["post"] > 1.5


def test_rank_rule_d6_gaia_orbit():
    from phase4_report import classify
    nss = [dict(nss_solution_type="Orbital", period=202.2, period_error=1.9)]
    bin_host = dict(flags=["gaia_nss"], host_binary=True, summary="NSS", ruwe=1.6, nss=nss)
    d = dict(sector=56, t0=0.0, partial=False, snr=10.0, kind="qlp")
    assert classify(_g(binarity=bin_host, duos=[d], joint_periods=[202.66]))[0] == "drop"
    assert classify(_g(binarity=bin_host, duos=[d], joint_periods=[101.3]))[0] == "maybe"
