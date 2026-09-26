"""Phase 6: find likely false positives among existing TESS candidates (TOIs, CTOIs).

    python scripts/phase6.py gaia [--refresh]         # tables + Gaia orbit / EB cross-match (batched)
    python scripts/phase6.py lc [--procs 2]           # light-curve checks on the selected samples
    python scripts/phase6.py lc --all --limit 500     # ... and on 500 more unchecked candidates
    python scripts/phase6.py report                   # validation, flags, table, summary
    python scripts/run_fp_triage.py                   # all three, as the app runs them

--refresh downloads today's TOI and CTOI tables (once) so new submissions are
included. Every light-curve check is recorded in results/phase6/lc_checks.csv.gz
(committed), so no candidate is checked twice, on this machine or any other.

Groups: "unresolved" (TFOPWG PC / APC / none; CTOIs not promoted to TOI),
"planet" (CP, KP) and "fp" (FP, FA). Planets and false positives are the
validation sets: every check is run on them too, and a check that flags more
than MAX_PLANET_FLAG_RATE of the planets is not used to flag candidates.

External services follow CLAUDE.md: the TOI/CTOI tables are downloaded once in
bulk; TIC and Gaia queries are batched (thousands of ids per query) through
tesshunt.net (cache, rate limit, backoff); light curves come from the S3
mirror and are deleted after reading.
"""

import argparse
import datetime
import glob
import json
import os
import sys
import tempfile
import time
import traceback
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from tesshunt import fptriage as fp  # noqa: E402

WORK = os.path.join(ROOT, "work", "phase6")
OUT = os.path.join(ROOT, "results", "phase6")
PLOTS = os.path.join(ROOT, "plots", "phase6")
LEDGER = os.path.join(OUT, "lc_checks.csv.gz")      # every light-curve check done so far
SAMPLE = 300
SEED = 6
MAX_PLANET_FLAG_RATE = 0.05          # a check flagging more planets than this is not used
PREFERRED_PLANET_RATE = 0.02         # among a check's variants, prefer the most sensitive below this

# Light-curve check variants, fixed before looking at the validation results.
# report() keeps, per check, the most sensitive variant whose planet flag rate
# is <= PREFERRED_PLANET_RATE (else <= MAX_PLANET_FLAG_RATE, else the check is dropped).
VARIANTS = {
    "odd_even": [
        ("3 sigma", lambda r: r.get("oddeven_sigma", 0) > 3),
        ("5 sigma", lambda r: r.get("oddeven_sigma", 0) > 5),
        ("5 sigma and >= 20 % different", lambda r: r.get("oddeven_sigma", 0) > 5 and r.get("oddeven_frac", 0) >= 0.2),
    ],
    "secondary": [
        ("SNR >= 7", lambda r: r.get("sec_snr", 0) >= 7),
        ("SNR >= 7 and >= 10 % of the transit depth", lambda r: r.get("sec_snr", 0) >= 7 and r.get("sec_ratio", 0) >= 0.1),
        ("SNR >= 10 and >= 10 % of the transit depth", lambda r: r.get("sec_snr", 0) >= 10 and r.get("sec_ratio", 0) >= 0.1),
    ],
    "centroid": [
        ("3 sigma", lambda r: r.get("centroid_sigma", 0) > 3),
        ("5 sigma", lambda r: r.get("centroid_sigma", 0) > 5),
        ("5 sigma and source >= 1 px away", lambda r: r.get("centroid_sigma", 0) > 5 and r.get("source_offset_px", 0) >= 1.0),
    ],
}


def jdump(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    os.replace(path + ".tmp", path)


# ------------------------------------------------------------------ gaia

def stage_gaia(args):
    os.makedirs(WORK, exist_ok=True)
    t = fp.load_tables(refresh=args.refresh)
    print(f"TOI/CTOI tables downloaded {fp.tables_date()}")
    ids = fp.gaia_ids(t.tic.unique())
    t = t.merge(ids, on="tic", how="left")
    t["m1"] = t.mstar.where(t.mstar > 0, t["mass"])
    print(f"{len(t)} candidates, {t.gaia_id.notna().sum()} with a Gaia DR3 id")
    nss, eb = fp.gaia_orbits(t.gaia_id.dropna().astype("int64").unique())
    print(f"Gaia NSS two-body orbits: {len(nss)} solutions on {nss.source_id.nunique() if len(nss) else 0} "
          f"sources; vari_eclipsing_binary: {len(eb)}")
    t.to_csv(os.path.join(WORK, "candidates.csv"), index=False)
    nss.to_csv(os.path.join(WORK, "nss.csv"), index=False)
    eb.to_csv(os.path.join(WORK, "gaia_eb.csv"), index=False)
    gaia_match(t, nss, eb).to_csv(os.path.join(WORK, "gaia_matches.csv"), index=False)


def gaia_match(t, nss, eb):
    """One row per (candidate, Gaia solution) whose period matches."""
    rows = []
    nss_by = {k: g.to_dict("records") for k, g in nss.groupby("source_id")} if len(nss) else {}
    eb_by = {k: g.to_dict("records") for k, g in eb.groupby("source_id")} if len(eb) else {}
    for r in t[t.periodic & t.gaia_id.notna()].itertuples():
        sid = int(r.gaia_id)
        for sol in nss_by.get(sid, []):
            m = fp.period_match(r.period, r.period_err, sol["period"], sol.get("period_error"))
            if not m:
                continue
            m2, how = fp.companion_mass(sol, r.m1)
            rows.append(dict(name=r.name, tic=r.tic, group=r.group, source="Gaia NSS " + str(sol["nss_solution_type"]),
                             gaia_period=sol["period"], gaia_period_err=sol.get("period_error"),
                             factor=m[0], n_sigma=m[1], m1=r.m1, m2=m2, m2_method=how,
                             eclipsing_solution="Eclipsing" in str(sol["nss_solution_type"])))
        for e in eb_by.get(sid, []):
            f = e.get("frequency")
            if not (f and f > 0):
                continue
            Pg = 1.0 / f
            sPg = (e.get("frequency_error") or 0) / f ** 2
            m = fp.period_match(r.period, r.period_err, Pg, sPg)
            if m:
                rows.append(dict(name=r.name, tic=r.tic, group=r.group, source="Gaia eclipsing binary (photometric)",
                                 gaia_period=Pg, gaia_period_err=sPg, factor=m[0], n_sigma=m[1],
                                 m1=r.m1, m2=np.nan, m2_method="", eclipsing_solution=True))
    return pd.DataFrame(rows)


def chance_matches(t, nss, eb, n_perm=200, seed=1):
    """Expected number of matches if Gaia periods were unrelated to the transits:
    Gaia periods are shuffled among the candidates that have a Gaia solution."""
    rng = np.random.default_rng(seed)
    with_nss = t[t.periodic & t.gaia_id.isin(set(nss.source_id) | set(eb.source_id))]
    pers = []
    for sid in with_nss.gaia_id.astype("int64"):
        p = list(nss.period[nss.source_id == sid]) + list(1 / eb.frequency[(eb.source_id == sid) & (eb.frequency > 0)])
        pers.append(p)
    counts = []
    for _ in range(n_perm):
        perm = rng.permutation(len(pers))
        n = 0
        for (i, r), j in zip(enumerate(with_nss.itertuples()), perm):
            if any(fp.period_match(r.period, r.period_err, pg, None) for pg in pers[j]):
                n += 1
        counts.append(n)
    return float(np.mean(counts)), float(np.std(counts)), len(with_nss)


# ------------------------------------------------------------------ light curves

def select_lc(t, matches, prior=None, full=False):
    """Candidates for the light-curve checks, with the sample each belongs to.

    Random samples of SAMPLE per group (the validation sets) are drawn once
    and then kept: ``prior`` (name -> sample, from the checks already done)
    fixes them even after the tables are refreshed. Gaia-matched candidates
    are always included. ``full`` adds every other unresolved candidate."""
    prior = prior or {}
    per = t[t.periodic]
    flagged = set(matches.name) if len(matches) else set()
    if any(v == "random" for v in prior.values()):
        pick = {n for n, v in prior.items() if v == "random"}
    else:
        rng = np.random.default_rng(SEED)
        pick = set()
        for g in ("unresolved", "planet", "fp"):
            pool = per[(per.group == g) & ~per.name.isin(flagged)]
            n = min(SAMPLE, len(pool))
            pick |= set(pool.name.values[rng.choice(len(pool), n, replace=False)])
    keep = per.name.isin(pick | flagged | set(prior))
    if full:
        keep |= per.group == "unresolved"
    sel = per[keep].copy()
    sel["sample"] = [prior.get(n) or ("gaia_flagged" if n in flagged else "random" if n in pick else "full")
                     for n in sel.name]
    return sel


def _json_path(name):
    return os.path.join(WORK, "lc", name.replace(" ", "_") + ".json")


def seed_from_ledger():
    """Write the committed light-curve results back into work/phase6/lc/, so
    candidates checked on any machine are not checked again. Returns count."""
    if not os.path.exists(LEDGER):
        return 0
    led = pd.read_csv(LEDGER)
    n = 0
    for r in led.to_dict("records"):
        p = _json_path(r["name"])
        if os.path.exists(p):
            continue
        r = {k: v for k, v in r.items() if not (isinstance(v, float) and np.isnan(v))}
        r["sectors"] = [int(x) for x in str(r.get("sectors", "")).split(";") if x.strip()]
        jdump(p, r)
        n += 1
    return n


def done_checks():
    """name -> result of every light-curve check on this machine."""
    out = {}
    for p in glob.glob(os.path.join(WORK, "lc", "*.json")):
        with open(p) as fh:
            r = json.load(fh)
        out[r["name"]] = r
    return out


def write_ledger(current_names=None):
    """Save every light-curve check to results/phase6/lc_checks.csv.gz."""
    rows = []
    for r in done_checks().values():
        r = dict(r)
        r["sectors"] = ";".join(map(str, r.get("sectors") or []))
        rows.append(r)
    if not rows:
        return 0
    df = pd.DataFrame(rows)
    if current_names is not None:
        df["in_current_tables"] = df.name.isin(current_names)
    first = ["name", "tic", "group", "sample", "checked", "tables_date", "n_events", "sectors"]
    df = df[[c for c in first if c in df] + [c for c in df if c not in first]].sort_values("name")
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(LEDGER + ".tmp", index=False, float_format="%.6g", compression="gzip")
    os.replace(LEDGER + ".tmp", LEDGER)
    return len(df)


def _others(t, r):
    """Ephemerides of the other candidates on the same star (masked in the scans)."""
    o = t[(t.tic == r["tic"]) & (t.name != r["name"]) & t.periodic]
    return [(float(x.period), float(x.epoch_btjd), float(x.duration_h) / 24) for x in o.itertuples()]


_T = None


def _init():
    global _T
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    _T = pd.read_csv(os.path.join(WORK, "candidates.csv"))


def lc_job(r):
    out = _json_path(r["name"])
    if os.path.exists(out):
        return r["name"], "cached"
    try:
        tmp = tempfile.mkdtemp(prefix="p6_", dir=os.path.join(ROOT, "work", "tmp"))
        res = fp.lc_checks(r["tic"], r["ra"], r["dec"], r["period"], r["epoch_btjd"], r["duration_h"],
                           r["depth_ppm"], others=_others(_T, r), tmp_dir=tmp)
        os.rmdir(tmp) if not os.listdir(tmp) else None
        res.update(name=r["name"], group=r["group"], sample=r["sample"],
                   checked=datetime.date.today().isoformat(), tables_date=r.get("tables_date"))
        jdump(out, res)
        return r["name"], "ok" if res.get("n_events") else "no data"
    except Exception as e:  # noqa: BLE001
        jdump(os.path.join(WORK, "errors", r["name"].replace(" ", "_") + ".json"),
              dict(error=f"{type(e).__name__}: {e}", tb=traceback.format_exc()))
        return r["name"], "error"


def stage_lc(args):
    n = seed_from_ledger()
    if n:
        print(f"{n} light-curve checks restored from {os.path.relpath(LEDGER, ROOT)}")
    t = pd.read_csv(os.path.join(WORK, "candidates.csv"))
    m = pd.read_csv(os.path.join(WORK, "gaia_matches.csv")) if os.path.getsize(
        os.path.join(WORK, "gaia_matches.csv")) > 1 else pd.DataFrame(columns=["name"])
    done = done_checks()
    sel = select_lc(t, m, prior={k: v.get("sample") for k, v in done.items()}, full=args.all)
    sel["tables_date"] = fp.tables_date()
    sel.to_csv(os.path.join(WORK, "lc_selection.csv"), index=False)
    todo = sel[~sel.name.isin(done)]
    # most useful first: Gaia-matched, then the validation samples, then the rest (fixed random order)
    rank = todo["sample"].map({"gaia_flagged": 0, "random": 1, "full": 2}).fillna(3)
    shuffle = np.random.default_rng(SEED).permutation(len(todo))
    todo = todo.assign(_r=rank.values, _s=shuffle).sort_values(["_r", "_s"]).drop(columns=["_r", "_s"])
    unres = sel[sel.group == "unresolved"]
    print(f"light-curve checks: {len(sel) - len(todo)} done, {len(todo)} to do "
          f"({unres.name.isin(done).sum()} of {int((t.periodic & (t.group == 'unresolved')).sum())} "
          "unresolved candidates checked so far)", flush=True)
    os.makedirs(os.path.join(ROOT, "work", "tmp"), exist_ok=True)
    recs = todo.to_dict("records")
    if args.limit:
        recs = recs[:args.limit]
    if not recs:
        print("nothing new to check")
        return
    t0, counts = time.time(), {}
    with Pool(args.procs, initializer=_init) as pool:
        for i, (k, s) in enumerate(pool.imap_unordered(lc_job, recs), 1):
            counts[s] = counts.get(s, 0) + 1
            if i % 25 == 0 or i == len(recs):
                print(f"[{time.strftime('%H:%M:%S')}] lc {i}/{len(recs)} {counts} "
                      f"{i / (time.time() - t0):.2f}/s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gaia")
    g.add_argument("--refresh", action="store_true", help="download today's TOI/CTOI tables")
    s = sub.add_parser("lc")
    s.add_argument("--procs", type=int, default=2)
    s.add_argument("--limit", type=int, default=0, help="check at most this many new candidates")
    s.add_argument("--all", action="store_true",
                   help="also check unresolved candidates outside the samples, not yet checked")
    sub.add_parser("report")
    args = ap.parse_args()
    if args.cmd == "report":
        from phase6_report import report
        report()
    else:
        {"gaia": stage_gaia, "lc": stage_lc}[args.cmd](args)


if __name__ == "__main__":
    main()
