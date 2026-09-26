"""Check other people's TESS candidates (TOIs, CTOIs) for signs of false positives.

    python scripts/run_fp_triage.py                      # 500 more candidates
    python scripts/run_fp_triage.py --refresh --limit 0  # today's tables, every candidate left

Steps (each resumable; the log format matches run_sector.py so the app can
follow it):

  fp-tables    TOI and CTOI tables (downloaded once; --refresh fetches today's)
               and the Gaia DR3 orbit / eclipsing-binary cross-match, in batches
  fp-lc        light-curve checks (odd/even, secondary, centroid) on candidates
               not checked before; results/phase6/lc_checks.csv.gz lists every
               check ever done, so nothing is checked twice
  fp-report    validation on known planets and false positives, flags, summary

Afterwards `python -m tesshunt.findings fp` says what is new and gives a
commit message.
"""

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def steps(args):
    lc = ["phase6.py", "lc", "--all", "--procs", str(args.procs)]
    if args.limit:
        lc += ["--limit", str(args.limit)]
    return [("fp-tables", ["phase6.py", "gaia"] + (["--refresh"] if args.refresh else [])),
            ("fp-lc", lc),
            ("fp-report", ["phase6.py", "report"])]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--refresh", action="store_true", help="download today's TOI and CTOI tables")
    ap.add_argument("--limit", type=int, default=500,
                    help="new candidates to check in this run (0 = all that are left)")
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    env = dict(os.environ, OMP_NUM_THREADS="1", PYTHONWARNINGS="ignore")
    for name, cmd in steps(args):
        argv = [sys.executable, os.path.join(HERE, cmd[0])] + cmd[1:]
        print(f"\n=== [{time.strftime('%H:%M:%S')}] {name}: {' '.join(argv[1:])}", flush=True)
        if args.dry_run:
            continue
        t0 = time.time()
        r = subprocess.run(argv, env=env)
        if r.returncode != 0:
            print(f"step '{name}' failed (exit {r.returncode}); fix and run again (finished work "
                  "is kept)", file=sys.stderr)
            sys.exit(r.returncode)
        print(f"=== {name} done in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
