"""Plain-English run summaries and commit messages (no network, no git writes)."""

import json
import os
import sys

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from tesshunt import findings, ledger  # noqa: E402


@pytest.fixture
def repo(tmp_path, monkeypatch):
    for mod in (ledger, findings):
        monkeypatch.setattr(mod, "ROOT", str(tmp_path))
    r = tmp_path / "results"
    (r / "phase2").mkdir(parents=True)
    (r / "phase2" / "s0022_summary.json").write_text(json.dumps(dict(searched=1000, candidates=20, status=[])))
    pd.DataFrame(dict(tic=[1], status=["ok"])).to_csv(r / "phase2" / "s0022_stars.csv.gz", index=False)
    (r / "phase3" / "s0022").mkdir(parents=True)
    (r / "phase3" / "s0022" / "summary.json").write_text(json.dumps(dict(survivors=3, survivors_new=2)))
    (r / "phase4" / "s0022").mkdir(parents=True)
    (r / "phase4" / "s0022" / "phase4_summary.json").write_text(
        json.dumps(dict(categories=dict(submit=1, maybe=1, drop=0))))
    pd.DataFrame(dict(tic=[11, 12], role="candidate", category=["submit", "maybe"],
                      reasons=[None, "S1 review notes: neighbour inside the target pixel"])
                 ).to_csv(r / "phase4" / "s0022" / "followup.csv", index=False)
    (r / "phase5").mkdir()
    (r / "phase6").mkdir()
    return tmp_path


def _p5(repo, rows):
    pd.DataFrame(rows).to_csv(repo / "results" / "phase5" / "phase5_checks.csv", index=False)


def test_sector_with_a_candidate(repo):
    _p5(repo, [dict(sector=22, tic=11, role="candidate", phase5="submit", verdict="strong", minor=None,
                    rp_rj=0.4, star_class="dwarf"),
               dict(sector=22, tic=12, role="candidate", phase5="maybe", verdict="strong", minor=None,
                    rp_rj=0.8, star_class="dwarf")])
    f = findings.sector_findings(22)
    assert f["found"] and f["title"] == "Sector 22: 1 candidate(s) worth submitting"
    assert "TIC 11: strong; ~0.40 Jupiter radii around a dwarf star; every check passed" in f["text"]
    assert "TIC 12" in f["text"] and "neighbour inside the target pixel" in f["text"]
    assert "**" not in f["body"]


def test_sector_with_nothing(repo):
    _p5(repo, [dict(sector=22, tic=12, role="candidate", phase5="maybe", verdict="plausible", minor="gp_weak",
                    rp_rj=0.8, star_class="dwarf")])
    f = findings.sector_findings(22)
    assert not f["found"] and "nothing new worth submitting" in f["title"]
    assert "should still be committed" in f["text"]
    assert findings.sector_findings(23)["done"] is False


def test_fp_findings_compare_with_the_last_commit(repo, monkeypatch):
    out = repo / "results" / "phase6"
    pd.DataFrame(dict(name=["TOI 1.01", "TOI 2.01"], tic=[1, 2], confidence=["high", "low"],
                      evidence=["Gaia orbit | odd/even", "centroid"])).to_csv(out / "likely_false_positives.csv",
                                                                               index=False)
    pd.DataFrame(dict(name=["TOI 1.01", "TOI 2.01", "TOI 3.01"])).to_csv(out / "lc_checks.csv.gz", index=False)
    committed = {"likely_false_positives.csv": pd.DataFrame(dict(name=["TOI 1.01", "TOI 9.01"])),
                 "lc_checks.csv.gz": pd.DataFrame(dict(name=["TOI 1.01"]))}
    monkeypatch.setattr(findings, "committed_csv", lambda p: committed.get(os.path.basename(p)))
    f = findings.fp_findings()
    assert f["found"] and "1 new likely false positive" in f["title"]
    assert "TOI 2.01 (TIC 2, low): centroid" in f["text"] and "TOI 1.01 (" not in f["text"]
    assert "2 more candidates" in f["text"] and "TOI 9.01" in f["text"]      # dropped flag reported


def test_commit_commands():
    c = findings.commit_commands("Title", "Body line")
    assert c.startswith("git add results plots\ngit commit -F - <<'EOF'\nTitle\n\nBody line\nEOF")
    assert c.endswith("git push")
