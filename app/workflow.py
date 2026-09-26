"""Per-candidate workflow: the four steps, their order, unlock rules and the
status file. Pure logic (no Streamlit), so it can be unit-tested.

Steps, which can only be completed in this order:
  feedback     prepare a Planet Hunters TESS forum post, then record the post
               link and the community's feedback (positive / negative / unclear)
  submission   prepare the ExoFOP CTOI fields; unlocked only by positive feedback,
               or by an override that states a reason
  recorded     record the ExoFOP submission date (and the CTOI number once assigned)
  followup     look at predicted transits and future TESS coverage

The status of every candidate is kept in one JSON file (default
app/workflow_status.json). It is the only file the app writes.
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import asdict, dataclass, field
from datetime import date, datetime

STEPS = ("feedback", "submission", "recorded", "followup")
STEP_TITLES = {
    "feedback": "a. Get feedback",
    "submission": "b. Prepare submission",
    "recorded": "c. Record submission",
    "followup": "d. Follow-up",
}
FEEDBACK_VALUES = ("positive", "negative", "unclear")
MIN_OVERRIDE_REASON = 10          # characters
DEFAULT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "workflow_status.json")
ENV_PATH = "TESSHUNT_WORKFLOW_FILE"      # overrides DEFAULT_PATH (used by the tests)


def status_path() -> str:
    return os.environ.get(ENV_PATH) or DEFAULT_PATH


class WorkflowError(ValueError):
    """A step was attempted out of order or with invalid input."""


@dataclass
class CandidateStatus:
    key: str                                   # e.g. "s0048_95747180"
    post_prepared: bool = False
    post_url: str = ""
    feedback: str = ""                         # positive | negative | unclear
    feedback_notes: str = ""
    feedback_done: bool = False
    override_reason: str = ""
    submission_prepared: bool = False
    submission_done: bool = False
    submitted_on: str = ""                     # ISO date
    ctoi: str = ""                             # e.g. "TIC 95747180.01" once assigned
    recorded_done: bool = False
    followup_done: bool = False
    history: list = field(default_factory=list)

    # ------------------------------------------------------------ queries
    def done(self, step: str) -> bool:
        return {"feedback": self.feedback_done, "submission": self.submission_done,
                "recorded": self.recorded_done, "followup": self.followup_done}[step]

    def current_step(self) -> str | None:
        for s in STEPS:
            if not self.done(s):
                return s
        return None

    def unlocked(self, step: str) -> bool:
        """A step is open when every earlier step is complete, and (for the
        submission step) the feedback was positive or an override was given."""
        i = STEPS.index(step)
        if not all(self.done(s) for s in STEPS[:i]):
            return False
        if step == "submission":
            return self.feedback == "positive" or bool(self.override_reason)
        return True

    def lock_reason(self, step: str) -> str:
        i = STEPS.index(step)
        missing = [STEP_TITLES[s] for s in STEPS[:i] if not self.done(s)]
        if missing:
            return "Finish " + ", ".join(missing) + " first."
        if step == "submission" and not self.unlocked(step):
            return ("Unlocked only when the feedback is positive. You can override this, "
                    "but you must type a reason.")
        return ""

    def progress(self) -> float:
        return sum(self.done(s) for s in STEPS) / len(STEPS)

    # ------------------------------------------------------------ actions
    def _log(self, what: str):
        self.history.append(dict(at=datetime.now().isoformat(timespec="seconds"), what=what))

    def _require(self, step: str):
        if not self.unlocked(step):
            raise WorkflowError(self.lock_reason(step) or f"{step} is locked")

    def mark_post_prepared(self):
        self._require("feedback")
        self.post_prepared = True
        self._log("forum post prepared")

    def record_feedback(self, post_url: str, feedback: str, notes: str = ""):
        self._require("feedback")
        if not self.post_prepared:
            raise WorkflowError("Prepare the forum post first.")
        post_url = (post_url or "").strip()
        if not post_url.startswith(("http://", "https://")):
            raise WorkflowError("Enter the link to your forum post (starting with https://).")
        if feedback not in FEEDBACK_VALUES:
            raise WorkflowError(f"Feedback must be one of {', '.join(FEEDBACK_VALUES)}.")
        self.post_url, self.feedback, self.feedback_notes = post_url, feedback, notes.strip()
        self.feedback_done = True
        self._log(f"feedback recorded: {feedback}")

    def override(self, reason: str):
        """Unlock the submission step despite non-positive feedback."""
        if not self.feedback_done:
            raise WorkflowError("Record the feedback first.")
        reason = (reason or "").strip()
        if len(reason) < MIN_OVERRIDE_REASON:
            raise WorkflowError(f"Type a reason of at least {MIN_OVERRIDE_REASON} characters.")
        self.override_reason = reason
        self._log(f"override: {reason}")

    def mark_submission_prepared(self):
        self._require("submission")
        self.submission_prepared = True
        self.submission_done = True
        self._log("submission prepared")

    def record_submission(self, submitted_on: str | date, ctoi: str = ""):
        self._require("recorded")
        d = submitted_on.isoformat() if isinstance(submitted_on, date) else str(submitted_on).strip()
        try:
            parsed = date.fromisoformat(d)
        except ValueError as e:
            raise WorkflowError("Enter the submission date as YYYY-MM-DD.") from e
        if parsed > date.today():
            raise WorkflowError("The submission date cannot be in the future.")
        self.submitted_on, self.ctoi = parsed.isoformat(), (ctoi or "").strip()
        self.recorded_done = True
        self._log(f"submission recorded {self.submitted_on} {self.ctoi}".strip())

    def set_ctoi(self, ctoi: str):
        """Add the CTOI number later, once ExoFOP assigns it."""
        if not self.recorded_done:
            raise WorkflowError("Record the submission first.")
        self.ctoi = (ctoi or "").strip()
        self._log(f"CTOI number: {self.ctoi}")

    def mark_followup_done(self):
        self._require("followup")
        self.followup_done = True
        self._log("follow-up reviewed")

    def reset(self):
        key = self.key
        self.__init__(key)
        self._log("reset")


# ------------------------------------------------------------------ storage

def load(path: str | None = None) -> dict[str, CandidateStatus]:
    """All statuses; a missing file means nothing started yet. A corrupt file is
    not overwritten silently: it raises, so the user can fix or move it."""
    path = path or status_path()
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        raw = json.load(fh)
    out = {}
    known = set(CandidateStatus.__dataclass_fields__)
    for key, d in raw.get("candidates", {}).items():
        out[key] = CandidateStatus(**{k: v for k, v in d.items() if k in known})
    return out


def save(statuses: dict[str, CandidateStatus], path: str | None = None):
    """Atomic write (temp file + rename), so a crash never leaves half a file."""
    path = path or status_path()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    data = dict(version=1, saved=datetime.now().isoformat(timespec="seconds"),
                candidates={k: asdict(v) for k, v in sorted(statuses.items())})
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), suffix=".tmp")
    with os.fdopen(fd, "w") as fh:
        json.dump(data, fh, indent=1)
    os.replace(tmp, path)


def get(statuses: dict[str, CandidateStatus], key: str) -> CandidateStatus:
    if key not in statuses:
        statuses[key] = CandidateStatus(key)
    return statuses[key]
