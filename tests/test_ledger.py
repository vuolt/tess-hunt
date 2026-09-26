"""Sector ledger and skip-list (no network)."""

import json
import os
import sys
import types

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from tesshunt import ledger  # noqa: E402
from tesshunt.survey import Store  # noqa: E402


@pytest.fixture
def repo(tmp_path, monkeypatch):
    """A fake results tree: S48 finished (top-level folders), S21 up to Phase 3."""
    monkeypatch.setattr(ledger, "ROOT", str(tmp_path))
    monkeypatch.setattr(ledger, "LEDGER", str(tmp_path / "results" / "sectors_searched.csv"))
    p2 = tmp_path / "results" / "phase2"
    p2.mkdir(parents=True)
    for s, n in ((48, 3), (21, 2)):
        (p2 / f"s{s:04d}_summary.json").write_text(json.dumps(dict(
            searched=n, candidates=1, status=[dict(status="ok", n=n), dict(status="error", n=1)])))
        pd.DataFrame(dict(tic=[1, 2, 3, 4][:n + 1], status=["ok"] * (n - 1) + ["skipped", "error"],
                          reason=[""] * (n + 1), tmag=12.0, n_dips=0, injected=0)
                     ).to_csv(p2 / f"s{s:04d}_stars.csv.gz", index=False)
        pd.DataFrame(dict(tic=[1, 1], rank_in_star=[1, 2], t0=[10.0, 11.0], duration_h=3.0,
                          depth_ppm=500.0, snr=8.0, edge_dist_h=5.0, n_in=6, rank=[1, 2], candidate=True)
                     ).to_csv(p2 / f"s{s:04d}_dips.csv", index=False)
        pd.DataFrame(dict(tic=[2], trial=[0], t0=[12.0], depth_ppm=[900.0], recovered=[1])
                     ).to_csv(p2 / f"s{s:04d}_injections.csv", index=False)
    for sub in ("", "s0021"):
        d = tmp_path / "results" / "phase3" / sub
        d.mkdir(parents=True, exist_ok=True)
        (d / "summary.json").write_text(json.dumps(dict(survivors=5, survivors_new=3)))
    (tmp_path / "results" / "phase4").mkdir()
    (tmp_path / "results" / "phase4" / "phase4_summary.json").write_text(
        json.dumps(dict(categories=dict(submit=1, maybe=2, drop=4))))
    return tmp_path


def test_ledger_rows_follow_the_committed_results(repo):
    df = ledger.rebuild().set_index("sector")
    assert list(df.index) == [21, 48]
    assert df.stage[48] == "phase4" and df.phase4_submit[48] == 1 and df.phase4_maybe[48] == 2
    assert df.stage[21] == "phase3" and pd.isna(df.phase4_submit[21]) and df.phase3_new[21] == 3
    assert df.stars_failed[48] == 1
    assert os.path.exists(ledger.LEDGER)
    assert "2 submit" not in ledger.describe(ledger.entry(21))
    assert "1 submit and 2 maybe" in ledger.describe(ledger.entry(48))
    assert ledger.entry(49) is None


def test_searched_stars_leave_out_errors(repo):
    assert ledger.searched_stars(48) == {1, 2, 3}          # star 4 failed: retried
    assert ledger.searched_stars(49) == set()
    assert ledger.searched_elsewhere([1, 2, 9], 48) == {21: 2}


def test_seed_store_copies_searched_stars_once(repo, tmp_path):
    st = Store(str(tmp_path / "work" / "s0048.sqlite"))
    st.write(dict(star=dict(tic=3, status="ok", reason="local", tmag=11.0)))   # already local
    assert ledger.seed_store(st, 48) == 2                   # stars 1 and 2; not 3 (local), not 4 (error)
    assert ledger.seed_store(st, 48) == 0                   # idempotent
    assert st.done() == {1, 2, 3}
    assert st.db.execute("SELECT reason FROM stars WHERE tic=3").fetchone()[0] == "local"
    assert st.db.execute("SELECT COUNT(*) FROM dips WHERE tic=1").fetchone()[0] == 2
    assert st.db.execute("SELECT recovered FROM injections WHERE tic=2").fetchone()[0] == 1


def test_run_sector_stops_for_a_finished_sector(repo, monkeypatch, capsys):
    import run_sector

    def boom(*a, **k):
        raise AssertionError("no step may run")
    monkeypatch.setattr(run_sector, "subprocess", types.SimpleNamespace(run=boom))
    monkeypatch.setattr(sys, "argv", ["run_sector.py", "--sector", "48"])
    run_sector.main()
    out = capsys.readouterr().out
    assert "already searched" in out and "Nothing to do" in out

    # a sector stopped after Phase 3 resumes; --force reruns a finished one
    for argv in (["--sector", "21"], ["--sector", "48", "--force"]):
        monkeypatch.setattr(sys, "argv", ["run_sector.py", *argv, "--dry-run"])
        run_sector.main()
        assert "=== [" in capsys.readouterr().out
