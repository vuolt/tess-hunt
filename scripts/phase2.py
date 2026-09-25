"""Phase 2: full-sector single-transit search on TESS-SPOC FFI light curves.

    python scripts/phase2.py select              # target list + TIC -> dwarf sample
    python scripts/phase2.py search              # resumable; Ctrl-C and rerun any time
    python scripts/phase2.py search --limit 200  # quick test on the first 200 pending stars
    python scripts/phase2_report.py              # CSVs, rankings, plots, summary numbers

State lives in work/ (git-ignored): the TIC table, the sample, and
work/s00XX.sqlite, which records every processed star. Light curves are
downloaded to work/tmp/ and deleted as soon as each star is analysed.
"""

import argparse
import csv
import os
import sys
import time
import urllib.request
from multiprocessing import Pool

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from tesshunt import ffi, tic  # noqa: E402
from tesshunt.survey import Store, process_star  # noqa: E402

WORK = os.path.join(ROOT, "work")


def paths(sector):
    s = f"s{sector:04d}"
    return dict(targets=os.path.join(WORK, f"{s}_targets.csv"),
                tic=os.path.join(WORK, f"{s}_tic.csv"),
                sample=os.path.join(WORK, f"{s}_sample.csv"),
                db=os.path.join(WORK, f"{s}.sqlite"),
                tmp=os.path.join(WORK, "tmp"))


def cmd_select(args):
    p = paths(args.sector)
    os.makedirs(WORK, exist_ok=True)
    if not os.path.exists(p["targets"]):
        urllib.request.urlretrieve(ffi.target_list_url(args.sector), p["targets"])
    ids = np.loadtxt(p["targets"], delimiter=",", skiprows=1, usecols=0, dtype=np.int64)
    print(f"S{args.sector}: {len(ids)} TESS-SPOC targets")

    if not os.path.exists(p["tic"]):
        t0 = time.time()
        table = tic.query_ids(ids, progress=lambda k, n: print(
            f"  TIC {k}/{n} ({time.time() - t0:.0f}s)", end="\r", flush=True))
        print()
        table.write(p["tic"], format="csv", overwrite=True)
    from astropy.table import Table
    table = Table.read(p["tic"], format="csv")
    sel = tic.select_dwarfs(table, args.tmag_max)
    lum = np.asarray(table["lumclass"].filled("") if hasattr(table["lumclass"], "filled")
                     else table["lumclass"])
    tm = np.asarray(table["Tmag"], float)
    counts = dict(targets=len(ids), tic_rows=len(table),
                  dwarf=int(np.sum(lum == "DWARF")),
                  dwarf_bright=int(np.sum((lum == "DWARF") & (tm < args.tmag_max))),
                  selected=int(sel.sum()))
    print(counts)
    table[sel].write(p["sample"], format="csv", overwrite=True)
    with open(os.path.join(WORK, f"s{args.sector:04d}_selection_counts.csv"), "w") as fh:
        w = csv.writer(fh)
        w.writerows(counts.items())


def cmd_search(args):
    p = paths(args.sector)
    os.makedirs(p["tmp"], exist_ok=True)
    sample = np.loadtxt(p["sample"], delimiter=",", skiprows=1, usecols=0, dtype=np.int64)
    sample = sorted(int(x) for x in sample)

    # Injection subset: fixed by seed, independent of processing order.
    rng = np.random.default_rng(args.seed)
    inj = set(int(x) for x in rng.choice(sample, size=min(args.inject_stars, len(sample)),
                                         replace=False))

    store = Store(p["db"])
    if args.retry_errors:
        with store.db:
            n = store.db.execute("DELETE FROM stars WHERE status='error'").rowcount
        print(f"retrying {n} errored stars")
    done = store.done()
    # Process injection stars first so the sensitivity result is available early.
    pending = [t for t in sample if t not in done]
    pending.sort(key=lambda t: (t not in inj, t))
    if args.limit:
        pending = pending[:args.limit]
    print(f"{len(sample)} stars in sample, {len(done)} done, {len(pending)} to process "
          f"({len(inj & set(pending))} injection stars pending)")

    t_start, n_done = time.time(), 0
    jobs = [(t, args.sector, p["tmp"], args.n_inject if t in inj else 0, args.seed + t)
            for t in pending]
    with Pool(args.procs, maxtasksperchild=500) as pool:
        for b in range(0, len(jobs), args.batch):
            batch = jobs[b:b + args.batch]
            status = {}
            for res in pool.imap_unordered(_run, batch, chunksize=4):
                store.write(res)
                s = res["star"]["status"]
                status[s] = status.get(s, 0) + 1
                if "traceback" in res and args.verbose:
                    print(res["traceback"])
            n_done += len(batch)
            rate = n_done / (time.time() - t_start)
            eta = (len(jobs) - n_done) / rate / 3600
            print(f"[{time.strftime('%H:%M:%S')}] batch {b // args.batch + 1}: "
                  f"{n_done}/{len(jobs)} {status}  {rate:.1f} stars/s  ETA {eta:.1f} h",
                  flush=True)


def _run(job):
    return process_star(*job)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=48)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--tmag-max", type=float, default=13.0)
    s = sub.add_parser("search")
    s.add_argument("--procs", type=int, default=12)
    s.add_argument("--batch", type=int, default=500)
    s.add_argument("--limit", type=int, default=0)
    s.add_argument("--inject-stars", type=int, default=200)
    s.add_argument("--n-inject", type=int, default=20, help="trials per injection star")
    s.add_argument("--seed", type=int, default=2024)
    s.add_argument("--retry-errors", action="store_true")
    s.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    {"select": cmd_select, "search": cmd_search}[args.cmd](args)


if __name__ == "__main__":
    main()
