"""Tests for the Streamlit app's workflow logic, run-log parsing and ExoFOP export."""

import json
import os
import sys
from datetime import date, timedelta

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app import runner, texts  # noqa: E402
from app import workflow as wf  # noqa: E402


def _through_feedback(s, fb="positive"):
    s.mark_post_prepared()
    s.record_feedback("https://www.zooniverse.org/x", fb, "notes")


# ------------------------------------------------------------------ order

def test_steps_must_go_in_order():
    s = wf.CandidateStatus("s0048_1")
    assert s.current_step() == "feedback"
    assert s.unlocked("feedback") and not s.unlocked("submission")
    assert not s.unlocked("recorded") and not s.unlocked("followup")
    with pytest.raises(wf.WorkflowError):
        s.mark_submission_prepared()
    with pytest.raises(wf.WorkflowError):
        s.record_submission(date.today())
    with pytest.raises(wf.WorkflowError):
        s.mark_followup_done()
    _through_feedback(s)
    assert s.current_step() == "submission" and s.unlocked("submission")
    with pytest.raises(wf.WorkflowError):
        s.mark_followup_done()                     # still needs b and c
    s.mark_submission_prepared()
    s.record_submission(date.today(), "")
    s.mark_followup_done()
    assert s.current_step() is None and s.progress() == 1.0


def test_feedback_needs_prepared_post_and_valid_input():
    s = wf.CandidateStatus("k")
    with pytest.raises(wf.WorkflowError):
        s.record_feedback("https://x", "positive")          # post not prepared
    s.mark_post_prepared()
    with pytest.raises(wf.WorkflowError):
        s.record_feedback("not a link", "positive")
    with pytest.raises(wf.WorkflowError):
        s.record_feedback("https://x", "great")
    s.record_feedback(" https://x ", "unclear", " n ")
    assert s.post_url == "https://x" and s.feedback_notes == "n" and s.feedback_done


# ------------------------------------------------------------------ unlock rules

@pytest.mark.parametrize("fb", ["negative", "unclear"])
def test_submission_locked_unless_positive(fb):
    s = wf.CandidateStatus("k")
    _through_feedback(s, fb)
    assert not s.unlocked("submission")
    assert "override" in s.lock_reason("submission")
    with pytest.raises(wf.WorkflowError):
        s.mark_submission_prepared()


def test_override_requires_a_reason():
    s = wf.CandidateStatus("k")
    with pytest.raises(wf.WorkflowError):
        s.override("because I think it's real")            # feedback not recorded yet
    _through_feedback(s, "negative")
    for bad in ("", "   ", "short"):
        with pytest.raises(wf.WorkflowError):
            s.override(bad)
    assert not s.unlocked("submission")
    s.override("second transit seen in new data")
    assert s.unlocked("submission")
    s.mark_submission_prepared()
    assert s.done("submission")


def test_record_submission_date_rules_and_late_ctoi():
    s = wf.CandidateStatus("k")
    _through_feedback(s)
    s.mark_submission_prepared()
    with pytest.raises(wf.WorkflowError):
        s.record_submission("26/09/2026")
    with pytest.raises(wf.WorkflowError):
        s.record_submission(date.today() + timedelta(days=2))
    s.record_submission("2026-09-20")
    assert s.recorded_done and s.ctoi == ""
    s.set_ctoi("TIC 1.01")
    assert s.ctoi == "TIC 1.01"


def test_set_ctoi_before_recording_fails():
    with pytest.raises(wf.WorkflowError):
        wf.CandidateStatus("k").set_ctoi("TIC 1.01")


def test_reset_clears_progress_but_keeps_key():
    s = wf.CandidateStatus("s0021_9")
    _through_feedback(s)
    s.reset()
    assert s.key == "s0021_9" and s.current_step() == "feedback" and not s.post_prepared
    assert s.history[-1]["what"] == "reset"


# ------------------------------------------------------------------ persistence

def test_save_and_load_roundtrip(tmp_path):
    p = str(tmp_path / "status.json")
    assert wf.load(p) == {}
    st = {}
    a = wf.get(st, "s0048_1")
    _through_feedback(a, "negative")
    a.override("seen again by another survey")
    b = wf.get(st, "s0021_2")
    b.mark_post_prepared()
    wf.save(st, p)
    back = wf.load(p)
    assert set(back) == {"s0048_1", "s0021_2"}
    assert back["s0048_1"].feedback == "negative" and back["s0048_1"].unlocked("submission")
    assert back["s0021_2"].post_prepared and not back["s0021_2"].feedback_done
    assert len(back["s0048_1"].history) == 3
    assert not any(f.endswith(".tmp") for f in os.listdir(tmp_path))   # atomic write cleaned up


def test_load_ignores_unknown_fields_and_rejects_corrupt_files(tmp_path):
    p = tmp_path / "status.json"
    p.write_text(json.dumps({"candidates": {"k": {"key": "k", "post_prepared": True, "future_field": 1}}}))
    assert wf.load(str(p))["k"].post_prepared
    p.write_text("{not json")
    with pytest.raises(ValueError):
        wf.load(str(p))
    assert p.read_text() == "{not json"            # never overwritten silently


# ------------------------------------------------------------------ run log

LOG = """
=== [07:42:54] select: phase2.py --sector 21 select
=== select done in 3.7 min

=== [07:46:38] search: phase2.py --sector 21 search
[07:49:41] batch 1: 500/2000 {'ok': 500}  2.7 stars/s
[07:52:52] batch 4: 1500/2000 {'ok': 500}  5.4 stars/s
"""


def test_progress_parsing():
    p = runner.progress(LOG)
    assert p["done"] == ["select"] and p["current"] == "search"
    assert p["fraction"] == pytest.approx(0.75)
    assert p["overall"] == pytest.approx((1 + 0.75) / len(runner.STEP_NAMES))
    assert not p["finished"] and p["failed"] is None


def test_progress_failure_and_resume():
    log = LOG + "step 'search' failed (exit 1); fix and rerun with --from search\n"
    assert runner.progress(log)["failed"] == "search"
    log += "\n=== [09:00:00] started from the app: run_sector.py --sector 21\n"
    assert runner.progress(log)["failed"] is None


def test_progress_finished_with_skipped_step():
    steps = [s for s in runner.STEP_NAMES if s != "injection-vetting"]
    log = "".join(f"=== [00:00:00] {s}: x\n=== {s} done in 1 min\n" for s in steps)
    assert runner.progress(log, steps)["finished"]
    assert not runner.progress(log)["finished"]


# ------------------------------------------------------------------ ExoFOP / text

def test_exofop_row_matches_template():
    f = dict(TIC=95747180, sector=48, epoch_bjd=2459621.1344, epoch_unc=0.0033, depth_ppm=7447,
             depth_unc=570, duration_h=4.55, duration_unc=0.7, radius_re=12.3, impact=0.5,
             period=None, period_note="")
    row = texts.exofop_row(f, "20260926_me_tesshunt", "https://doi.org/x")
    cols = row.split("|")
    assert len(cols) == len(texts.EXOFOP_COLUMNS) == 46
    d = dict(zip(texts.EXOFOP_COLUMNS, cols))
    assert d["target"] == "TIC95747180.01" and d["flag"] == "newctoi" and d["disp"] == "PC"
    assert d["epoch"] == "2459621.1344" and d["depth"] == "7447" and d["period"] == ""
    assert d["prop_period"] == "0" and d["paper"] == "https://doi.org/x"
    assert len(d["notes"]) <= 120
    assert texts.exofop_file([row]).splitlines()[0] == "|".join(texts.EXOFOP_COLUMNS)


def test_glossary_tooltip():
    html = texts.term("FPP")
    assert "<abbr" in html and "false-positive" in html.lower()
    assert texts.term("unknown word") == "unknown word"


def test_exofop_columns_match_official_template():
    """Column order copied from ExoFOP's template (fetched 2026-09-26 from
    texts.EXOFOP_TEMPLATE and stored in tests/fixtures)."""
    p = os.path.join(os.path.dirname(__file__), "fixtures", "exofop_params_planet_template.txt")
    header = next(ln for ln in open(p).read().splitlines() if ln.startswith("target|"))
    assert header.split("|") == list(texts.EXOFOP_COLUMNS)
