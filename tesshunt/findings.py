"""What a finished run found, in plain English, with a ready commit message.

    python -m tesshunt.findings sector 22     # after run_sector.py --sector 22
    python -m tesshunt.findings fp            # after run_fp_triage.py

Read-only: it looks at the result files and at git (to compare with what
is already committed) and never writes or commits anything. The app shows
the same text when a run finishes.
"""

from __future__ import annotations

import io
import os
import subprocess
import sys
import textwrap

import pandas as pd

from . import ledger

ROOT = ledger.ROOT
COMMIT_PATHS = ["results", "plots"]


# ------------------------------------------------------------------ git (read-only)

def _git(*args) -> str | None:
    try:
        r = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout if r.returncode == 0 else None


def changed_files(paths=COMMIT_PATHS) -> list[str]:
    """Files under ``paths`` that differ from the last commit (new or modified)."""
    out = _git("status", "--porcelain", "--untracked-files=all", "--", *paths)
    if out is None:
        return []
    return [line[3:].strip().strip('"') for line in out.splitlines() if line.strip()]


def committed_csv(path: str) -> pd.DataFrame | None:
    """The last committed version of a CSV (None if it was never committed)."""
    rel = os.path.relpath(path, ROOT).replace(os.sep, "/")
    try:
        r = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=ROOT, capture_output=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    if r.returncode != 0 or not r.stdout:
        return None
    comp = "gzip" if rel.endswith(".gz") else None
    return pd.read_csv(io.BytesIO(r.stdout), compression=comp)


def _wrap(body: str) -> str:
    """Wrap a commit body at 72 characters, keeping list items as items."""
    out = []
    for line in body.split("\n"):
        if line.startswith("- "):
            out.append(textwrap.fill(line, 72, subsequent_indent="  "))
        else:
            out.append(textwrap.fill(line, 72) if line else "")
    return "\n".join(out)


def commit_commands(title: str, body: str, paths=COMMIT_PATHS) -> str:
    """Commands to paste into a terminal in the project folder (macOS / Linux)."""
    msg = title + ("\n\n" + _wrap(body) if body else "")
    return (f"git add {' '.join(paths)}\n"
            f"git commit -F - <<'EOF'\n{msg}\nEOF\n"
            "git push")


# ------------------------------------------------------------------ sector runs

def _read(path):
    return pd.read_csv(path) if os.path.exists(path) else None


def _why(r) -> str:
    """One line on a candidate from its Phase 5 row."""
    size = f"~{r['rp_rj']:.2f} Jupiter radii" if pd.notna(r.get("rp_rj")) else "size unknown"
    star = r.get("star_class") if isinstance(r.get("star_class"), str) else "star type unknown"
    notes = []
    minor = r.get("minor") if isinstance(r.get("minor"), str) else ""
    words = dict(gp_weak="the independent reanalysis finds it only moderately (5-7 sigma)",
                 evolved="the star is slightly evolved", eccentric="needs an eccentric orbit",
                 star_unverified="Gaia cannot confirm the star's size",
                 aperture_marginal="the pixel test is borderline", ruwe="Gaia sees a noisy star",
                 variable="the star is variable", companion="Gaia hints at a companion")
    notes = [words.get(m, m) for m in minor.split(";") if m]
    s = f"{r['verdict']}; {size} around a {star} star"
    return s + ("; " + ", ".join(notes) if notes else "; every check passed")


def sector_findings(sector: int) -> dict:
    """Plain-English summary of a sector run and a commit message."""
    row = ledger.row(sector)
    if not row:
        return dict(done=False, text=f"Sector {sector} has no results yet.", title="", body="", files=[])
    sub = "" if int(sector) == 48 else f"s{int(sector):04d}"
    fu = _read(os.path.join(ROOT, "results", "phase4", sub, "followup.csv"))
    p5 = _read(os.path.join(ROOT, "results", "phase5", "phase5_checks.csv"))
    lines = [ledger.describe(row).replace(" was already searched", " searched")]
    submit, strong_held = [], []
    if p5 is not None:
        c = p5[(p5.sector == int(sector)) & (p5.role == "candidate")]
        for r in c.to_dict("records"):
            if r["phase5"] == "submit":
                submit.append(f"TIC {r['tic']}: {_why(r)}")
            elif r["verdict"] == "strong":
                reason = ""
                if fu is not None:
                    m = fu[fu.tic == r["tic"]]
                    if len(m) and isinstance(m.iloc[0].get("reasons"), str):
                        reason = f" (held back: {m.iloc[0]['reasons']})"
                strong_held.append(f"TIC {r['tic']}: strong but still 'maybe'{reason}")
        n_maybe = int((c.phase5 == "maybe").sum())
    elif fu is not None:
        c = fu[fu.role == "candidate"]
        submit = [f"TIC {r.tic}: Phase 4 submit (expert checks not run yet)"
                  for r in c[c.category == "submit"].itertuples()]
        n_maybe = int((c.category == "maybe").sum())
    else:
        n_maybe = 0
    found = bool(submit or strong_held)
    if submit:
        lines.append(f"**{len(submit)} candidate(s) worth submitting:**\n" + "\n".join(f"- {x}" for x in submit))
    if strong_held:
        lines.append("**Strong, but held back for a manual look:**\n" + "\n".join(f"- {x}" for x in strong_held))
    if not found:
        lines.append("**Nothing new worth submitting in this sector.** The record of searched "
                     "stars and sectors should still be committed, so nobody searches them again.")
    if n_maybe:
        lines.append(f"{n_maybe} more are 'maybe' (see the Candidates page).")
    title = (f"Sector {sector}: {len(submit)} candidate(s) worth submitting" if submit else
             f"Sector {sector}: strong candidates to check by eye" if strong_held else
             f"Sector {sector} searched: nothing new worth submitting")
    body = "\n\n".join(x.replace("**", "") for x in lines)
    return dict(done=True, found=found, text="\n\n".join(lines), title=title, body=body,
                files=changed_files())


# ------------------------------------------------------------------ false-positive runs

def fp_findings() -> dict:
    """New false-positive flags since the last commit, and a commit message."""
    out = os.path.join(ROOT, "results", "phase6")
    now = _read(os.path.join(out, "likely_false_positives.csv"))
    if now is None:
        return dict(done=False, text="No false-positive check has been run yet.", title="", body="",
                    files=[])
    before = committed_csv(os.path.join(out, "likely_false_positives.csv"))
    led_now = pd.read_csv(os.path.join(out, "lc_checks.csv.gz")) \
        if os.path.exists(os.path.join(out, "lc_checks.csv.gz")) else pd.DataFrame(columns=["name"])
    led_before = committed_csv(os.path.join(out, "lc_checks.csv.gz"))
    n_checked = len(set(led_now.name) - set(led_before.name if led_before is not None else []))
    old = set(before.name) if before is not None else set()
    new = now[~now.name.isin(old)]
    gone = sorted(old - set(now.name))
    lines = [f"{len(now)} candidates are flagged as likely false positives in total "
             f"({(now.confidence == 'high').sum()} high confidence). This run checked "
             f"{n_checked} more candidates' light curves."]
    if len(new):
        items = [f"- {r.name} (TIC {r.tic}, {r.confidence}): {str(r.evidence).split(' | ')[0]}"
                 for r in new.head(40).itertuples()]
        if len(new) > 40:
            items.append(f"- ... and {len(new) - 40} more (results/phase6/likely_false_positives.csv)")
        lines.append(f"**{len(new)} new flag(s):**\n" + "\n".join(items))
    else:
        lines.append("**No new false positives found.** The record of what was checked should "
                     "still be committed, so nobody checks the same candidates again.")
    if gone:
        lines.append(f"{len(gone)} earlier flag(s) dropped (the candidate was resolved by TFOPWG "
                     "or removed from the tables): " + ", ".join(gone[:10])
                     + (" ..." if len(gone) > 10 else ""))
    title = (f"False-positive check: {len(new)} new likely false positive(s)" if len(new) else
             f"False-positive check: {n_checked} more candidates checked, nothing new")
    body = "\n\n".join(x.replace("**", "") for x in lines)
    return dict(done=True, found=bool(len(new)), text="\n\n".join(lines), title=title, body=body,
                files=changed_files())


def main(argv):
    if len(argv) >= 2 and argv[0] == "sector":
        f = sector_findings(int(argv[1]))
    elif argv[:1] == ["fp"]:
        f = fp_findings()
    else:
        sys.exit(__doc__)
    print(f["text"])
    if f["done"]:
        print(f"\n{len(f['files'])} new or changed result files. To save them:\n")
        print(commit_commands(f["title"], f["body"]))


if __name__ == "__main__":
    main(sys.argv[1:])
