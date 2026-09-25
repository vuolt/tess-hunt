"""Phase 1 demo: blind single-event search on a star with a known long-period planet.

Default target is TOI-2180 b (TIC 298663873; Dalba et al. 2022), a ~260 d
Jupiter that TESS sees as isolated ~24 h transits. The search is blind: the
published ephemeris is only used afterwards to score the detections.

    python scripts/run_known_target.py
"""

import argparse
import csv
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tesshunt.lightcurves import load_sectors  # noqa: E402
from tesshunt.pipeline import SearchConfig, search_sector  # noqa: E402
from tesshunt.plotting import event_figure, overview_figure  # noqa: E402

# TOI-2180 b ephemeris, Dalba et al. 2024 via the NASA Exoplanet Archive.
# BTJD = BJD - 2457000.
TOI2180 = dict(tic=298663873, name="TOI-2180", period=260.1706, t0=1830.7561,
               duration_h=24.1, depth=4750e-6)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def known_transits(t0, period, tmin, tmax):
    n = np.arange(np.ceil((tmin - t0) / period), np.floor((tmax - t0) / period) + 1)
    return list(t0 + n * period)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tic", type=int, default=TOI2180["tic"])
    p.add_argument("--name", default=TOI2180["name"])
    p.add_argument("--sectors", type=int, nargs="*")
    p.add_argument("--threshold", type=float, default=SearchConfig.threshold)
    args = p.parse_args()

    cfg = SearchConfig(threshold=args.threshold)
    lcs = load_sectors(args.tic, sectors=args.sectors)
    print(f"{args.name}: {len(lcs)} sectors")

    is_toi2180 = args.tic == TOI2180["tic"]
    searches, rows = [], []
    t_start = time.time()
    for lc in lcs:
        ss = search_sector(lc, cfg)
        searches.append(ss)
        for ev in ss.search.events:
            rows.append(ev)
        print(f"  S{lc.sector:02d}: {len(lc.time):6d} pts, "
              f"{len(ss.search.events)} events"
              + "".join(f"  [t0={e.t0:.2f} SNR={e.snr:.1f} "
                        f"{e.depth*1e6:.0f}ppm {e.duration_h:.0f}h]"
                        for e in ss.search.events))
    print(f"search time {time.time() - t_start:.1f}s")

    tmin = min(lc.time.min() for lc in lcs)
    tmax = max(lc.time.max() for lc in lcs)
    known = (known_transits(TOI2180["t0"], TOI2180["period"], tmin, tmax)
             if is_toi2180 else [])
    observed = [t for t in known
                if any(lc.time.min() <= t <= lc.time.max()
                       and np.min(np.abs(lc.time - t)) < 0.1 for lc in lcs)]

    slug = args.name.lower().replace(" ", "").replace("-", "")
    os.makedirs(os.path.join(ROOT, "plots"), exist_ok=True)
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)

    csv_path = os.path.join(ROOT, "results", f"{slug}_events.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sector", "t0_btjd", "duration_h", "depth_ppm", "snr",
                    "matches_known_transit", "offset_h"])
        for ev in rows:
            off = min((abs(ev.t0 - t) for t in known), default=np.inf)
            w.writerow([ev.sector, f"{ev.t0:.4f}", f"{ev.duration_h:.1f}",
                        f"{ev.depth*1e6:.0f}", f"{ev.snr:.1f}",
                        off < TOI2180["duration_h"] / 48, f"{off*24:.2f}"
                        if np.isfinite(off) else ""])

    overview_figure(searches, f"{args.name} (TIC {args.tic}): blind single-event "
                    f"search, all sectors (dotted = published transit times)",
                    os.path.join(ROOT, "plots", f"{slug}_overview.png"),
                    cfg.threshold, known)

    # Diagnostic for the strongest event in each sector that has one above threshold.
    for ss in searches:
        if not ss.search.events:
            continue
        ev = ss.search.events[0]
        event_figure(ss, ev, f"{args.name} (TIC {args.tic}) sector {ss.lc.sector}: "
                     f"strongest single event",
                     os.path.join(ROOT, "plots", f"{slug}_s{ss.lc.sector:02d}_event.png"),
                     cfg.threshold, known)

    if known:
        print(f"\nPublished transits inside observed data: "
              f"{[round(t, 2) for t in observed]}")
        for t in observed:
            hits = [e for e in rows if abs(e.t0 - t) < TOI2180["duration_h"] / 48]
            if hits:
                e = max(hits, key=lambda e: e.snr)
                print(f"  {t:.2f}: RECOVERED in S{e.sector} t0={e.t0:.3f} "
                      f"(offset {(e.t0 - t)*24:+.1f} h) depth={e.depth*1e6:.0f} ppm "
                      f"dur={e.duration_h:.0f} h SNR={e.snr:.1f}")
            else:
                print(f"  {t:.2f}: missed")
        others = [e for e in rows
                  if min(abs(e.t0 - t) for t in known) >= TOI2180["duration_h"] / 48]
        print(f"Other events above threshold: {len(others)}")


if __name__ == "__main__":
    main()
