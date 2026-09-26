"""Start, watch and stop `scripts/run_sector.py` in the background.

The run's state (PID, sector, log file) lives in work/app_runs/, the pipeline's
git-ignored scratch area, never in results/. A run started from the app keeps
going if the app is closed; reopening the app finds it again.

Resuming: every pipeline step caches its work and skips what is already done,
so starting the same sector again continues an interrupted run.
"""

from __future__ import annotations

import json
import os
import re
import signal
import subprocess
import sys
import time

from . import data

RUN_DIR = data.path("work", "app_runs")
STATE = os.path.join(RUN_DIR, "current.json")
STEP_NAMES = ["select", "search", "phase2-report", "phase3-lc", "phase3-pixels", "phase3-fpp",
              "phase3-report", "phase4", "pht-plots", "expert-checks", "expert-report",
              "injection-vetting"]
FP_STEP_NAMES = ["fp-tables", "fp-lc", "fp-report"]
STEP_WORDS = {
    "select": "choosing the stars to search",
    "search": "searching every light curve for dips",
    "phase2-report": "summarising the search",
    "phase3-lc": "vetting: dip shape, duration and data gaps",
    "phase3-pixels": "vetting: checking the raw pixels and neighbouring stars",
    "phase3-fpp": "vetting: statistical false-positive test",
    "phase3-report": "writing the vetting report",
    "phase4": "follow-up: other sectors, orbits and companions",
    "pht-plots": "drawing forum-style light curves of the best candidates",
    "expert-checks": "expert checks: pixels, noise model, Gaia star and companions, orbit shape",
    "expert-report": "writing the verdicts and plain-English summaries",
    "fp-tables": "loading the TOI / CTOI lists and cross-matching Gaia orbits",
    "fp-lc": "checking light curves for eclipsing-binary signs",
    "fp-report": "comparing with known planets and writing the list of likely false positives",
    "injection-vetting": "measuring how many real planets the checks would keep",
}


def _read_state() -> dict | None:
    try:
        with open(STATE) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except (ProcessLookupError, PermissionError, OSError):
        return False
    cmd = f"/proc/{pid}/cmdline"          # guard against PID reuse where /proc exists
    if os.path.exists(cmd):
        try:
            line = open(cmd, "rb").read()
            return b"run_sector.py" in line or b"run_fp_triage.py" in line
        except OSError:
            return False
    return True


def current() -> dict | None:
    st = _read_state()
    if st:
        st["running"] = _alive(int(st["pid"]))
    return st


def label(st: dict) -> str:
    return "False-positive check" if st.get("kind") == "fp" else f"Sector {st['sector']}"


def _launch(cmd: list, log: str, steps: list, **info) -> dict:
    st = current()
    if st and st["running"]:
        raise RuntimeError(f"{label(st)} is already running; only one run at a time.")
    os.makedirs(RUN_DIR, exist_ok=True)
    with open(log, "a") as fh:
        fh.write(f"\n=== [{time.strftime('%H:%M:%S')}] started from the app: {' '.join(cmd[1:])}\n")
        p = subprocess.Popen(cmd, cwd=data.ROOT, stdout=fh, stderr=subprocess.STDOUT,
                             stdin=subprocess.DEVNULL, start_new_session=True)
    st = dict(pid=p.pid, log=log, started=time.time(), cmd=cmd, steps=steps, **info)
    with open(STATE, "w") as fh:
        json.dump(st, fh)
    return st


def start(sector: int, procs: int = 3, skip_injection: bool = False, force: bool = False,
          _cmd: list | None = None) -> dict:
    """Start run_sector.py in its own process group (``force`` reruns a sector
    that was already searched; ``_cmd`` replaces the command in tests)."""
    cmd = [sys.executable, data.path("scripts", "run_sector.py"), "--sector", str(int(sector)),
           "--procs", str(int(procs))]
    if skip_injection:
        cmd += ["--skip", "injection-vetting"]
    if force:
        cmd += ["--force"]
    steps = [x for x in STEP_NAMES if not (skip_injection and x == "injection-vetting")]
    return _launch(_cmd or cmd, os.path.join(RUN_DIR, f"run_s{int(sector):04d}.log"), steps,
                   kind="sector", sector=int(sector))


def start_fp(procs: int = 2, refresh: bool = False, limit: int = 500, _cmd: list | None = None) -> dict:
    """Start run_fp_triage.py: check other people's candidates for false positives."""
    cmd = [sys.executable, data.path("scripts", "run_fp_triage.py"), "--procs", str(int(procs)),
           "--limit", str(int(limit))]
    if refresh:
        cmd += ["--refresh"]
    return _launch(_cmd or cmd, os.path.join(RUN_DIR, "run_fp_triage.log"), list(FP_STEP_NAMES),
                   kind="fp", sector=None)


def stop() -> bool:
    """Stop the run and every worker it started (they share a process group)."""
    st = current()
    if not st or not st["running"]:
        return False
    pid = int(st["pid"])
    try:
        os.killpg(os.getpgid(pid), signal.SIGTERM)
    except (ProcessLookupError, PermissionError, AttributeError):
        os.kill(pid, signal.SIGTERM)
    with open(st["log"], "a") as fh:
        fh.write(f"\n=== [{time.strftime('%H:%M:%S')}] stopped from the app\n")
    return True


_START = re.compile(r"^=== \[[\d:]+\] ([\w-]+):")
_DONE = re.compile(r"^=== ([\w-]+) done")
_FAIL = re.compile(r"^step '([\w-]+)' failed")
_FRAC = re.compile(r"(\d+)/(\d+)")


def progress(log_text: str, steps: list[str] | None = None) -> dict:
    """Parse a run log: which steps are done, the current step, and how far
    through it the last "n/N" counter says the run is. ``steps`` are the steps
    this run was asked to do (skipped ones print nothing). Pure, for tests.
    Only the part of the log after the last "started from the app" line counts,
    except that steps finished by earlier attempts stay done (they are cached)."""
    steps = steps or STEP_NAMES
    known = set(STEP_NAMES) | set(FP_STEP_NAMES) | set(steps)
    done, current_step, failed, frac = [], None, None, None
    for line in log_text.splitlines():
        m = _START.match(line)
        if m and m.group(1) in known:
            current_step, frac = m.group(1), None
            continue
        m = _DONE.match(line)
        if m and m.group(1) in known:
            if m.group(1) not in done:
                done.append(m.group(1))
            if current_step == m.group(1):
                current_step, frac = None, None
            continue
        m = _FAIL.match(line)
        if m:
            failed, current_step, frac = m.group(1), None, None
            continue
        if "started from the app" in line:
            failed = None
        if current_step:
            ms = _FRAC.findall(line)
            if ms:
                a, b = map(int, ms[-1])
                if 0 < b and a <= b:
                    frac = a / b
    done = [x for x in done if x in steps]
    overall = (len(done) + (frac or 0.0)) / len(steps)
    return dict(done=done, current=current_step, fraction=frac, overall=min(overall, 1.0),
                failed=failed, finished=all(x in done for x in steps))


def tail(path: str, n: int = 15) -> str:
    try:
        with open(path) as fh:
            return "".join(fh.readlines()[-n:])
    except OSError:
        return ""


def log_text(path: str) -> str:
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return ""
