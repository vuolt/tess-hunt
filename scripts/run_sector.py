"""Run the whole single-transit pipeline for one TESS sector with one command.

    python scripts/run_sector.py --sector 49
    python scripts/run_sector.py --sector 49 --from phase3-lc     # resume from a step
    python scripts/run_sector.py --sector 49 --skip injection-vetting
    python scripts/run_sector.py --sector 48 --force               # run a searched sector again
    python -m tesshunt.findings sector 49                          # what was found + commit message

Steps (each is resumable and caches its work, so rerunning the command after
an interruption continues where it stopped):

  select               TESS-SPOC target list + TIC -> dwarf sample (Tmag < 13)
  search               Phase 2 search of every star (+ injection-recovery on 200)
  phase2-report        vetting labels, ranking, candidate list, figures
  phase3-lc            vetting checks 1-3 (shape, duration, edge)
  phase3-pixels        checks 4-6 (TESScut pixels, SkyBoT, catalogues)
  phase3-fpp           check 7 (TRICERATOPS)
  phase3-report        funnel, shortlist, vetting sheets
  phase4               other sectors, second transits, allowed periods, binarity,
                       submit / maybe / drop, CTOI summaries, follow-up sheets
  pht-plots            light curves of the submit / maybe candidates in the style of
                       Planet Hunters TESS, for the forum
  expert-checks        Phase 5 on every searched sector's submit / maybe candidates
                       (cached, so earlier sectors cost little): aperture test, GP,
                       Gaia star / binarity / variability, density-prior fit
  expert-report        verdicts, updated ranking and CTOI summaries (results/phase5)
  injection-vetting    the Phase 2 injections through all seven checks

Sectors already searched are listed in results/sectors_searched.csv (rebuilt
from the committed results after every run). Starting a finished one (Phase 4
done) again stops with a note unless --force or --from is given; even then, stars already in
results/phase2/sXXXX_stars.csv.gz are not downloaded or searched again.

External services are used politely (tesshunt/net.py): bulk light curves from
the AWS S3 mirror first, every response cached in work/cache, small services
(SkyBoT, ExoFOP, Gaia, VizieR, MAST catalogue and TESScut queries) serialized
at <= ~1.7 requests/s with exponential backoff.
"""

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from tesshunt import ledger  # noqa: E402

STEPS = [
    ("select", ["phase2.py", "--sector", "{s}", "select", "--tmag-max", "{tmag}"]),
    ("search", ["phase2.py", "--sector", "{s}", "search", "--procs", "{procs}"]),
    ("phase2-report", ["phase2_report.py", "--sector", "{s}"]),
    ("phase3-lc", ["phase3_vet.py", "lc", "--sector", "{s}", "--procs", "{procs}"]),
    ("phase3-pixels", ["phase3_vet.py", "pixels", "--sector", "{s}", "--procs", "2"]),
    ("phase3-fpp", ["phase3_vet.py", "fpp", "--sector", "{s}", "--procs", "{fpp_procs}"]),
    ("phase3-report", ["phase3_vet.py", "report", "--sector", "{s}"]),
    ("phase4", ["phase4.py", "all", "--sector", "{s}", "--procs", "{procs}"]),
    ("pht-plots", ["pht_plots.py", "--sector", "{s}"]),
    ("expert-checks", ["phase5.py", "run", "--all-sectors", "--procs", "{fpp_procs}"]),
    ("expert-report", ["phase5.py", "report"]),
    ("injection-vetting", ["phase4_injection_vetting.py", "all", "--sector", "{s}",
                           "--procs", "{procs}"]),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--tmag-max", type=float, default=13.0)
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--fpp-procs", type=int, default=2)
    ap.add_argument("--from", dest="start", choices=[s for s, _ in STEPS])
    ap.add_argument("--skip", nargs="*", default=[], choices=[s for s, _ in STEPS])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true", help="run a sector that was already searched")
    args = ap.parse_args()
    done = ledger.entry(args.sector)
    if done:
        print(ledger.describe(done), flush=True)
        if done["stage"] == "phase4" and not (args.force or args.start):
            print("Nothing to do. Results are in results/phase2-4; use --force to run it again "
                  "(already-searched stars are still not searched again), or --from STEP to redo "
                  "later steps.", flush=True)
            return
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               PYTHONWARNINGS="ignore")
    started = args.start is None
    for name, cmd in STEPS:
        started = started or name == args.start
        if not started or name in args.skip:
            continue
        argv = [sys.executable, os.path.join(HERE, cmd[0])] + [
            c.format(s=args.sector, tmag=args.tmag_max, procs=args.procs,
                     fpp_procs=args.fpp_procs) for c in cmd[1:]]
        print(f"\n=== [{time.strftime('%H:%M:%S')}] {name}: {' '.join(argv[1:])}", flush=True)
        if args.dry_run:
            continue
        t0 = time.time()
        r = subprocess.run(argv, env=env)
        if r.returncode != 0:
            print(f"step '{name}' failed (exit {r.returncode}); fix and rerun with "
                  f"--from {name}", file=sys.stderr)
            sys.exit(r.returncode)
        print(f"=== {name} done in {(time.time() - t0) / 60:.1f} min", flush=True)
    if not args.dry_run:
        ledger.rebuild()
        print(f"ledger updated: {os.path.relpath(ledger.LEDGER, os.path.dirname(HERE))}", flush=True)


if __name__ == "__main__":
    main()
