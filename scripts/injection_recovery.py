"""Injection-recovery test for the single-transit detector.

Synthetic trapezoid transits are injected into the *raw* (pre-detrending)
PDCSAP light curves of TOI-2180, one per trial, and the full pipeline
(detrend -> search -> re-detrend -> search) is run. Sectors that contain a
known TOI-2180 b transit are left out, so the only real signal present is
whatever noise and systematics the star and spacecraft provide.

    python scripts/injection_recovery.py            # full grid
    python scripts/injection_recovery.py --n 5      # quick smoke run
"""

import argparse
import csv
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tesshunt.injection import inject, run_trial, trial_row  # noqa: E402
from tesshunt.lightcurves import load_sectors  # noqa: E402
from tesshunt.pipeline import SearchConfig, search_sector  # noqa: E402
from tesshunt.plotting import event_figure  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TIC = 298663873
SECTORS_WITH_KNOWN_TRANSIT = {19, 48, 57}   # TOI-2180 b, see run_known_target.py

DEPTHS_PPM = [50, 100, 200, 300, 500, 1000, 2000, 5000]
DURATIONS_H = [1, 2, 4, 8, 16, 24]

CFG = SearchConfig()
OUT = os.path.join(ROOT, "results")
_LCS = {}
_NOISE = {}


def _trial(args):
    sector, t0, depth, dur = args
    return trial_row(run_trial(_LCS[sector], t0, depth, dur, CFG, _NOISE[sector]))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=60, help="trials per (depth, duration) cell")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--procs", type=int, default=os.cpu_count())
    p.add_argument("--baseline-coverage", type=float, default=CFG.baseline_coverage)
    p.add_argument("--out", default=os.path.join(ROOT, "results"),
                   help="directory for the CSV and summary")
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args()
    CFG.baseline_coverage = args.baseline_coverage
    global OUT
    OUT = args.out

    lcs = [lc for lc in load_sectors(TIC) if lc.sector not in SECTORS_WITH_KNOWN_TRANSIT]
    print(f"{len(lcs)} clean sectors")

    # Null run: the un-injected sectors. Gives the false-alarm rate and the
    # empirical per-duration noise used for expected SNR.
    null_events, null_max = [], []
    for lc in lcs:
        ss = search_sector(lc, CFG)
        _LCS[lc.sector] = lc
        _NOISE[lc.sector] = dict(zip(ss.search.durations, ss.search.sigma))
        null_events += ss.search.events
        null_max.append(np.nanmax(ss.search.max_snr))
    print(f"null: {len(null_events)} events above SNR {CFG.threshold} "
          f"in {len(lcs)} sectors; max SNR per sector median "
          f"{np.median(null_max):.1f}, max {np.max(null_max):.1f}")

    rng = np.random.default_rng(args.seed)
    sectors = sorted(_LCS)
    jobs = []
    for depth in DEPTHS_PPM:
        for dur in DURATIONS_H:
            for _ in range(args.n):
                s = sectors[rng.integers(len(sectors))]
                t = _LCS[s].time
                t0 = float(t[rng.integers(len(t))])   # centre on a real cadence
                jobs.append((s, t0, depth * 1e-6, dur / 24))

    t_start = time.time()
    with Pool(args.procs) as pool:
        rows = pool.map(_trial, jobs, chunksize=8)
    print(f"{len(rows)} trials in {time.time() - t_start:.0f}s")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "injection_recovery.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    summarize(rows, null_events, null_max, len(lcs), plots=not args.no_plots)
    if not args.no_plots:
        example_plot(rng)


def grid_fraction(rows, key="recovered"):
    frac = np.zeros((len(DEPTHS_PPM), len(DURATIONS_H)))
    for i, d in enumerate(DEPTHS_PPM):
        for j, h in enumerate(DURATIONS_H):
            cell = [r[key] for r in rows if r["depth_ppm"] == d and r["duration_h"] == h]
            frac[i, j] = np.mean(cell)
    return frac


def summarize(rows, null_events, null_max, n_sectors, plots=True):
    frac = grid_fraction(rows)
    top = grid_fraction(rows, "top_ranked")
    n = len(rows) // (len(DEPTHS_PPM) * len(DURATIONS_H))

    # Table for RESULTS.md
    lines = ["| depth \\ duration | " + " | ".join(f"{h} h" for h in DURATIONS_H) + " |",
             "|---|" + "---|" * len(DURATIONS_H)]
    for i, d in enumerate(DEPTHS_PPM):
        lines.append(f"| {d} ppm | " + " | ".join(f"{frac[i, j]:.0%}" for j in range(len(DURATIONS_H))) + " |")
    table = "\n".join(lines)

    esnr = np.array([r["expected_snr"] for r in rows])
    rec = np.array([r["recovered"] for r in rows])
    bins = np.array([0, 2, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 50, 1e9])
    idx = np.digitize(esnr, bins) - 1
    snr_lines = ["| expected SNR | trials | recovered |", "|---|---|---|"]
    centers, fr = [], []
    for k in range(len(bins) - 1):
        sel = idx == k
        if sel.sum() == 0:
            continue
        hi = "inf" if bins[k + 1] > 1e8 else f"{bins[k+1]:g}"
        snr_lines.append(f"| {bins[k]:g}-{hi} | {sel.sum()} | {rec[sel].mean():.0%} |")
        centers.append(min(np.median(esnr[sel]), 80))
        fr.append(rec[sel].mean())

    with open(os.path.join(OUT, "injection_summary.md"), "w") as fh:
        fh.write(f"Baseline coverage required per side: {CFG.baseline_coverage}; "
                 f"SNR threshold {CFG.threshold}\n")
        fh.write(f"Trials per cell: {n}; total {len(rows)}\n\n")
        fh.write("Recovery fraction (event above threshold within half a duration of the injected T0):\n\n")
        fh.write(table + "\n\nRecovery vs expected SNR:\n\n" + "\n".join(snr_lines) + "\n\n")
        fh.write(f"Null (no injection): {len(null_events)} events above threshold in "
                 f"{n_sectors} sectors; per-sector max SNR median {np.median(null_max):.1f}, "
                 f"max {np.max(null_max):.1f}\n")
        fh.write("Null events: " + ", ".join(f"S{e.sector} t0={e.t0:.2f} SNR={e.snr:.1f} "
                                             f"{e.depth*1e6:.0f}ppm {e.duration_h:.0f}h"
                                             for e in null_events) + "\n")
    print(table)
    print("\n".join(snr_lines))
    if not plots:
        return

    # Figures
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, f, ttl in [(axes[0], frac, "Recovered (above threshold, within T/2)"),
                       (axes[1], top, "Recovered AND strongest event in sector")]:
        im = ax.imshow(f, origin="lower", cmap="viridis", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(DURATIONS_H)), [f"{h}" for h in DURATIONS_H])
        ax.set_yticks(range(len(DEPTHS_PPM)), [f"{d}" for d in DEPTHS_PPM])
        ax.set_xlabel("Injected duration [h]")
        ax.set_ylabel("Injected depth [ppm]")
        ax.set_title(ttl)
        for i in range(len(DEPTHS_PPM)):
            for j in range(len(DURATIONS_H)):
                ax.text(j, i, f"{f[i, j]:.0%}", ha="center", va="center",
                        color="w" if f[i, j] < 0.6 else "k", fontsize=8)
    fig.colorbar(im, ax=axes, shrink=0.8, label="fraction")
    fig.suptitle(f"Single-transit injection-recovery, TOI-2180 (T=8.6), {n_sectors} "
                 f"sectors, {n} trials/cell, SNR threshold {CFG.threshold:g}")
    fig.savefig(os.path.join(ROOT, "plots", "injection_recovery_grid.png"), dpi=130,
                bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    ax = axes[0]
    for j, h in enumerate(DURATIONS_H):
        ax.plot(DEPTHS_PPM, frac[:, j], "o-", label=f"{h} h")
    ax.set_xscale("log")
    ax.set_xlabel("Injected depth [ppm]")
    ax.set_ylabel("Recovery fraction")
    ax.set_ylim(-0.03, 1.03)
    ax.legend(title="duration", fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    ax.plot(centers, fr, "o-", color="k")
    ax.axvline(CFG.threshold, color="C3", ls="--", lw=1, label="detection threshold")
    ax.set_xscale("log")
    ax.set_xlabel("Expected SNR (depth / empirical noise at that duration)")
    ax.set_ylabel("Recovery fraction")
    ax.set_ylim(-0.03, 1.03)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "plots", "injection_recovery_curves.png"), dpi=130)
    plt.close(fig)


def example_plot(rng):
    """Diagnostic figure for one representative injection (500 ppm, 8 h)."""
    lc = _LCS[sorted(_LCS)[5]]
    t0 = float(np.median(lc.time))
    t0 = float(lc.time[np.argmin(np.abs(lc.time - t0))])
    ss = search_sector(inject(lc, t0, 500e-6, 8 / 24), CFG)
    if ss.search.events:
        event_figure(ss, ss.search.events[0],
                     f"Injected 500 ppm / 8 h transit at {t0:.2f} (dotted), sector {lc.sector}",
                     os.path.join(ROOT, "plots", "injection_example.png"),
                     CFG.threshold, [t0], known_label="injected T0")


if __name__ == "__main__":
    main()
