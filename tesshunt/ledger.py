"""Which sectors and stars have already been searched.

    python -m tesshunt.ledger          # rebuild results/sectors_searched.csv and print it

The ledger is rebuilt from the committed results, so it cannot drift from
them: a sector is listed once results/phase2/sXXXX_summary.json exists, and
its later columns fill in as the Phase 3 and Phase 4 summaries appear.

Stars: results/phase2/sXXXX_stars.csv.gz lists every star searched in a
sector. `seed_store` copies those stars (and their dips and injection
trials) into the local search database, so a fresh clone or a rerun never
downloads or searches them again. Stars that failed with an error are not
copied and are retried. A star searched in one sector is still searched in
another: every sector is new data.
"""

from __future__ import annotations

import datetime
import json
import os
import subprocess

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEDGER = os.path.join(ROOT, "results", "sectors_searched.csv")
DONE_STATUS = ("ok", "skipped")
COLUMNS = ["sector", "stage", "stars_searched", "stars_failed", "phase2_candidates",
           "phase3_survivors", "phase3_new", "phase4_submit", "phase4_maybe", "results_date",
           "stars_file"]


def _tag(sector: int) -> str:
    return f"s{int(sector):04d}"


def _results(phase: int, sector: int, name: str) -> str:
    """Sector 48 predates multi-sector support and uses the top-level folders."""
    sub = "" if int(sector) == 48 else _tag(sector)
    return os.path.join(ROOT, "results", f"phase{phase}", sub, name)


def stars_file(sector: int) -> str:
    return os.path.join(ROOT, "results", "phase2", f"{_tag(sector)}_stars.csv.gz")


def _json(path: str) -> dict | None:
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _commit_date(path: str) -> str:
    try:
        out = subprocess.run(["git", "log", "-1", "--format=%cs", "--", path], cwd=ROOT,
                             capture_output=True, text=True, timeout=20).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        out = ""
    return out or datetime.date.today().isoformat()


def sectors_with_results() -> list[int]:
    d = os.path.join(ROOT, "results", "phase2")
    if not os.path.isdir(d):
        return []
    return sorted(int(f[1:5]) for f in os.listdir(d)
                  if f.startswith("s") and f.endswith("_summary.json") and f[1:5].isdigit())


def row(sector: int) -> dict | None:
    """One ledger row from the committed results of a sector (None if not searched)."""
    p2 = _json(os.path.join(ROOT, "results", "phase2", f"{_tag(sector)}_summary.json"))
    if not p2:
        return None
    failed = sum(s["n"] for s in p2.get("status", []) if s.get("status") not in DONE_STATUS)
    r = dict(sector=int(sector), stage="phase2", stars_searched=p2.get("searched"), stars_failed=failed,
             phase2_candidates=p2.get("candidates"), phase3_survivors=None, phase3_new=None,
             phase4_submit=None, phase4_maybe=None, results_date=_commit_date(stars_file(sector)),
             stars_file=os.path.relpath(stars_file(sector), ROOT))
    p3 = _json(_results(3, sector, "summary.json"))
    if p3:
        r.update(stage="phase3", phase3_survivors=p3.get("survivors"), phase3_new=p3.get("survivors_new"))
    p4 = _json(_results(4, sector, "phase4_summary.json"))
    if p4:
        cat = p4.get("categories", {})
        r.update(stage="phase4", phase4_submit=cat.get("submit", 0), phase4_maybe=cat.get("maybe", 0))
    return r


def rebuild(write: bool = True) -> pd.DataFrame:
    rows = [r for r in (row(s) for s in sectors_with_results()) if r]
    df = pd.DataFrame(rows, columns=COLUMNS)
    for c in COLUMNS[2:9]:
        df[c] = df[c].astype("Int64")
    if write:
        os.makedirs(os.path.dirname(LEDGER), exist_ok=True)
        df.to_csv(LEDGER + ".tmp", index=False)
        os.replace(LEDGER + ".tmp", LEDGER)
    return df


def entry(sector: int) -> dict | None:
    """The ledger row of a sector, from the committed results."""
    return row(sector)


def describe(e: dict) -> str:
    """One plain sentence about an already-searched sector."""
    s = (f"Sector {e['sector']} was already searched (results committed {e['results_date']}): "
         f"{e['stars_searched']:,} stars, {e['phase2_candidates']:,} candidate dips")
    if e.get("phase3_survivors") is not None:
        s += f", {e['phase3_survivors']} passed vetting ({e['phase3_new']} new)"
    if e.get("phase4_submit") is not None:
        s += f", {e['phase4_submit']} submit and {e['phase4_maybe']} maybe after follow-up"
    return s + "."


def _stars(sector: int) -> pd.DataFrame | None:
    p = stars_file(sector)
    if not os.path.exists(p):
        return None
    return pd.read_csv(p)


def searched_stars(sector: int) -> set[int]:
    """TIC IDs already searched in this sector (status ok or skipped)."""
    s = _stars(sector)
    if s is None:
        return set()
    return set(int(t) for t in s.tic[s.status.isin(DONE_STATUS)])


def searched_elsewhere(tics, sector: int) -> dict[int, int]:
    """For each other searched sector, how many of `tics` were searched there."""
    tics = set(int(t) for t in tics)
    return {s: len(tics & searched_stars(s)) for s in sectors_with_results() if s != int(sector)}


def seed_store(store, sector: int) -> int:
    """Copy stars already searched in this sector (from the committed results)
    into the local search database, with their dips and injection trials, so
    they are not downloaded or searched again. Returns the number copied."""
    s = _stars(sector)
    if s is None:
        return 0
    done = store.done()
    new = s[s.status.isin(DONE_STATUS) & ~s.tic.isin(done)]
    if not len(new):
        return 0
    cols = {t: [r[1] for r in store.db.execute(f"PRAGMA table_info({t})")]
            for t in ("stars", "dips", "injections")}
    keep = set(int(t) for t in new.tic)
    tables = {"stars": new}
    for t in ("dips", "injections"):
        p = os.path.join(ROOT, "results", "phase2", f"{_tag(sector)}_{t}.csv")
        df = pd.read_csv(p) if os.path.exists(p) else pd.DataFrame(columns=["tic"])
        tables[t] = df[df.tic.isin(keep)]
    with store.db:
        for t, df in tables.items():
            c = [x for x in cols[t] if x in df.columns]
            if not len(df):
                continue
            vals = df[c].astype(object).where(df[c].notna(), None).values.tolist()
            store.db.executemany(f"INSERT {'OR REPLACE ' if t == 'stars' else ''}INTO {t} "
                                 f"({','.join(c)}) VALUES ({','.join('?' * len(c))})", vals)
    return len(new)


if __name__ == "__main__":
    print(rebuild().to_string(index=False))
