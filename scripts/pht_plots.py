"""Planet Hunters TESS-style light-curve plots for the submit / maybe candidates.

    python scripts/pht_plots.py --sector 48          # skips plots that exist
    python scripts/pht_plots.py --sector 48 --redo

Reads results/phase4[/sXXXX]/followup.csv and writes
plots/phase4[/sXXXX]/pht/tic<TIC>_pht.png: the whole sector as PHT shows it
(normalised PDCSAP brightness, not detrended) plus a zoom on the dip. Light
curves come from the S3 mirror first (tesshunt/net.py) and stay cached in
work/cache, so reruns download nothing.
"""

import argparse
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from tesshunt import ffi  # noqa: E402
from tesshunt.phtplot import pht_plot  # noqa: E402

CATEGORIES = ("submit", "maybe")


def dirs(sector):
    sub = "" if sector == 48 else f"s{sector:04d}"
    return (os.path.join(ROOT, "results", "phase4", sub, "followup.csv"),
            os.path.join(ROOT, "plots", "phase4", sub, "pht"))


def plot_path(sector, tic):
    return os.path.join(dirs(sector)[1], f"tic{int(tic)}_pht.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--redo", action="store_true")
    args = ap.parse_args()
    src, out = dirs(args.sector)
    if not os.path.exists(src):
        sys.exit(f"no {os.path.relpath(src, ROOT)}: run Phase 4 first")
    f = pd.read_csv(src)
    f = f[f.category.isin(CATEGORIES) & (f.role == "candidate")]
    made = 0
    for r in f.itertuples():
        p = plot_path(args.sector, r.tic)
        if os.path.exists(p) and not args.redo:
            continue
        try:
            lc, _ = ffi.read(ffi.download(int(r.tic), args.sector, cache=True))
        except Exception as e:  # noqa: BLE001
            print(f"TIC {r.tic}: no light curve ({type(e).__name__}: {e})")
            continue
        title = (f"TIC {r.tic} · Sector {args.sector} · TESS mag {r.tmag:.1f} · "
                 f"pipeline verdict: {r.category}")
        pht_plot(lc.time, lc.flux, r.t0_btjd, r.t14_h / 24, r.depth_ppm * 1e-6, title, p,
                 sector=args.sector)
        made += 1
    print(f"{made} PHT-style plots written to {os.path.relpath(out, ROOT)} "
          f"({len(f)} submit/maybe candidates)")


if __name__ == "__main__":
    main()
