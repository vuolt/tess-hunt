"""End-to-end test of the app's workflow page (Streamlit AppTest, no browser).

Runs only if the repository's Phase 4 results exist (they are committed)."""

import os
import sys
from datetime import date

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

st_testing = pytest.importorskip("streamlit.testing.v1")
pytestmark = pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "results", "phase4", "followup.csv")),
                                reason="no Phase 4 results")
APP = os.path.join(ROOT, "app", "main.py")


def _app(page, status_file):
    os.environ["TESSHUNT_WORKFLOW_FILE"] = str(status_file)
    at = st_testing.AppTest.from_file(APP, default_timeout=120)
    at.session_state["page"] = page
    return at


def _by_label(widgets, text):
    return next(w for w in widgets if text in (w.label or ""))


@pytest.mark.parametrize("page", ["Candidates", "Workflow", "Results", "Run", "Glossary"])
def test_pages_render(page, tmp_path):
    at = _app(page, tmp_path / "s.json")
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    assert not (tmp_path / "s.json").exists()            # viewing never writes


def test_workflow_end_to_end(tmp_path):
    from app import workflow as wf
    sf = tmp_path / "status.json"
    at = _app("Workflow", sf)
    at.run()
    key = at.session_state["wf_key"]
    # a. feedback: prepare, then record negative feedback
    at.button(key="prep").click().run()
    assert wf.load(str(sf))[key].post_prepared
    _by_label(at.text_input, "Link to your forum post").input("https://www.zooniverse.org/talk/1")
    _by_label(at.radio, "community say").set_value("negative")
    next(b for b in at.button if b.label == "Save feedback").click().run()
    s = wf.load(str(sf))[key]
    assert s.feedback_done and s.feedback == "negative" and not s.unlocked("submission")
    # b. locked; a too-short override is refused, a real reason unlocks it
    _by_label(at.text_input, "Override anyway").input("short")
    next(b for b in at.button if b.label == "Override").click().run()
    assert any("at least" in e.value for e in at.error)
    assert not wf.load(str(sf))[key].unlocked("submission")
    _by_label(at.text_input, "Override anyway").input("new photometry shows a second dip")
    next(b for b in at.button if b.label == "Override").click().run()
    assert wf.load(str(sf))[key].unlocked("submission")
    assert any("refereed" in w.value for w in at.warning)          # ExoFOP rule shown
    next(b for b in at.button if b.label.startswith("Submission prepared")).click().run()
    assert wf.load(str(sf))[key].submission_done
    # c. record (date defaults to today)
    next(b for b in at.button if b.label == "Save").click().run()
    s = wf.load(str(sf))[key]
    assert s.recorded_done and s.submitted_on == date.today().isoformat()
    # d. follow-up
    next(b for b in at.button if b.label == "Mark follow-up reviewed").click().run()
    s = wf.load(str(sf))[key]
    assert s.followup_done and s.progress() == 1.0
    assert not at.exception


def _snapshot():
    out = {}
    for top in ("results", "plots", "scripts", "tesshunt"):
        for dp, _, fs in os.walk(os.path.join(ROOT, top)):
            for f in fs:
                p = os.path.join(dp, f)
                st = os.stat(p)
                out[p] = (st.st_size, st.st_mtime_ns)
    return out


def test_app_never_touches_results_or_code(tmp_path):
    before = _snapshot()
    for page in ("Candidates", "Workflow", "Results", "Run", "Glossary"):
        _app(page, tmp_path / "s.json").run()
    test_workflow_end_to_end(tmp_path)
    assert _snapshot() == before


def test_runner_start_progress_stop(tmp_path, monkeypatch):
    import time
    from app import runner
    monkeypatch.setattr(runner, "RUN_DIR", str(tmp_path))
    monkeypatch.setattr(runner, "STATE", str(tmp_path / "current.json"))
    child = "import time; print('=== [00:00:00] select: x', flush=True); time.sleep(60)"
    # the child's command line contains 'run_sector.py' so the PID-reuse guard accepts it
    st = runner.start(99, _cmd=[sys.executable, "-c", child, "run_sector.py"])
    time.sleep(1.0)
    cur = runner.current()
    assert cur["running"] and cur["sector"] == 99
    with pytest.raises(RuntimeError):
        runner.start(99, _cmd=[sys.executable, "-c", child, "run_sector.py"])   # one run at a time
    assert runner.progress(runner.log_text(st["log"]), cur["steps"])["current"] == "select"
    assert runner.stop()
    for _ in range(50):
        if not runner.current()["running"]:
            break
        time.sleep(0.1)
    assert not runner.current()["running"]
    assert "stopped from the app" in runner.log_text(st["log"])


def test_run_page_warns_about_a_searched_sector(tmp_path):
    at = _app("Run", tmp_path / "s.json")
    at.run()
    _by_label(at.number_input, "Sector").set_value(48).run()
    assert any("already searched" in w.value for w in at.warning)
    start = _by_label(at.button, "Start")
    assert start.disabled                                   # finished: needs "run again anyway"
    _by_label(at.checkbox, "Run this sector again").check().run()
    assert not _by_label(at.button, "Start").disabled
    assert not (tmp_path / "s.json").exists()


def test_false_positive_run_is_followed_to_the_end(tmp_path, monkeypatch):
    import time
    from app import runner
    monkeypatch.setattr(runner, "RUN_DIR", str(tmp_path))
    monkeypatch.setattr(runner, "STATE", str(tmp_path / "current.json"))
    lines = "".join(f"print('=== [00:00:00] {s}: x'); print('=== {s} done in 0.1 min');"
                    for s in runner.FP_STEP_NAMES)
    st = runner.start_fp(_cmd=[sys.executable, "-c", lines, "run_fp_triage.py"])
    for _ in range(50):
        if not runner.current()["running"]:
            break
        time.sleep(0.1)
    cur = runner.current()
    assert cur["kind"] == "fp" and runner.label(cur) == "False-positive check"
    assert runner.progress(runner.log_text(st["log"]), cur["steps"])["finished"]


def test_run_page_offers_both_kinds_of_run(tmp_path):
    at = _app("Run", tmp_path / "s.json")
    at.run()
    assert not at.exception
    labels = [b.label for b in at.button]
    assert any("Start / resume the search" in x for x in labels)
    assert any("Start the check" in x for x in labels)
    assert not (tmp_path / "s.json").exists()
