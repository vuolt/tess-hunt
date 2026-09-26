"""Stars with single-transit dips in two searched sectors (possible duotransits).

    python scripts/compare_sectors.py --sectors 21 48

Cross-matches the two sectors by TIC at four levels:

  dips        every dip Phase 2 recorded (results/phase2/sXXXX_dips.csv)
  candidates  Phase 2 candidates (automated vetting)
  vetted      events passing Phase 3 checks 1-3 (shape, duration, edge)
  shortlist   events passing all seven Phase 3 checks (shortlist.csv, known included)

For each star in both it pairs the best (highest-SNR) dip of each sector and
checks whether the two could be the same object: depth ratio within 0.5-2 and
duration ratio within 0.5-2 (Phase 2 durations are on a coarse grid). It lists
the implied periods dt/n down to 10 d, and the chance expectation (stars
dipping in both sectors if dips were independent of the star):
N_a * N_b / N_common, with N_a, N_b the dipping stars among the N_common stars
searched in both sectors.

Writes results/phase4/s<A>/duotransits_with_s<B>.md and .csv.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from phase3_vet import sector_dirs  # noqa: E402

RATIO = (0.5, 2.0)
P_MIN_ALIAS = 10.0


def phase2(sector, what):
    return pd.read_csv(os.path.join(ROOT, "results", "phase2", f"s{sector:04d}_{what}.csv"))


def searched(sector):
    s = pd.read_csv(os.path.join(ROOT, "results", "phase2", f"s{sector:04d}_stars.csv.gz"),
                    usecols=["tic", "status"])
    return set(s.tic[s.status == "ok"])


def phase3(sector):
    _, res, _ = sector_dirs(sector, 3)
    va = pd.read_csv(os.path.join(res, "vetting_all.csv"))
    va = va[va.is_candidate & (va.lc_status == "ok")]
    vetted = va[va["shape"].eq(True) & va["duration"].eq(True) & va["edge"].eq(True)]
    sl = pd.read_csv(os.path.join(res, "shortlist.csv"))
    return vetted, sl


def levels(sector):
    d = phase2(sector, "dips")
    vetted, sl = phase3(sector)
    std = lambda df, t0, dur, dep, snr: pd.DataFrame(dict(  # noqa: E731
        tic=df.tic.astype(int), t0=df[t0], dur_h=df[dur], depth_ppm=df[dep], snr=df[snr]))
    out = {
        "dips": std(d, "t0", "duration_h", "depth_ppm", "snr"),
        "candidates": std(d[d.candidate], "t0", "duration_h", "depth_ppm", "snr"),
        "vetted": std(vetted, "fit_t0", "fit_t14_h", "fit_depth_ppm", "snr"),
        "shortlist": std(sl, "fit_t0", "fit_t14_h", "fit_depth_ppm", "snr"),
    }
    extra = {"category": d.set_index(["tic", "t0"]).category.to_dict() if "category" in d else {},
             "shortlist_new": sl.set_index("tic").new.to_dict(),
             "notes": sl.set_index("tic").review_notes.to_dict()}
    return out, extra


def phase4_status(sector):
    _, res, _ = sector_dirs(sector, 4)
    p = os.path.join(res, "followup.csv")
    if not os.path.exists(p):
        return {}
    f = pd.read_csv(p)
    return f.set_index("tic").category.to_dict()


def pair(a, b):
    """Best dip per star in each sector, matched by TIC."""
    ba = a.sort_values("snr", ascending=False).drop_duplicates("tic")
    bb = b.sort_values("snr", ascending=False).drop_duplicates("tic")
    m = ba.merge(bb, on="tic", suffixes=("_a", "_b"))
    m["depth_ratio"] = m.depth_ppm_b / m.depth_ppm_a
    m["dur_ratio"] = m.dur_h_b / m.dur_h_a
    m["consistent"] = (m.depth_ratio.between(*RATIO) & m.dur_ratio.between(*RATIO))
    m["dt_d"] = (m.t0_b - m.t0_a).abs()
    m["n_alias_ge10d"] = np.floor(m.dt_d / P_MIN_ALIAS).astype(int)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sectors", type=int, nargs=2, required=True)
    args = ap.parse_args()
    sa, sb = args.sectors
    common = searched(sa) & searched(sb)
    la, xa = levels(sa)
    lb, xb = levels(sb)
    st4 = {sa: phase4_status(sa), sb: phase4_status(sb)}
    rows, summary = [], []
    for lev in ("dips", "candidates", "vetted", "shortlist"):
        a = la[lev][la[lev].tic.isin(common)]
        b = lb[lev][lb[lev].tic.isin(common)]
        m = pair(a, b)
        na, nb = a.tic.nunique(), b.tic.nunique()
        exp = na * nb / max(len(common), 1)
        summary.append(dict(level=lev, stars_a=na, stars_b=nb, both=len(m), expected_by_chance=exp,
                            consistent=int(m.consistent.sum())))
        m.insert(0, "level", lev)
        rows.append(m)
    allm = pd.concat(rows, ignore_index=True)
    for s, lab in ((sa, "a"), (sb, "b")):
        allm[f"phase4_{lab}"] = allm.tic.map(st4[s])
    allm[f"new_in_s{sa}"] = allm.tic.map(xa["shortlist_new"])
    allm[f"new_in_s{sb}"] = allm.tic.map(xb["shortlist_new"])
    _, out, _ = sector_dirs(sa, 4)
    os.makedirs(out, exist_ok=True)
    stem = os.path.join(out, f"duotransits_with_s{sb:04d}")
    allm.to_csv(stem + ".csv", index=False, float_format="%.6g")
    sm = pd.DataFrame(summary)

    L = [f"# Stars with dips in both Sector {sa} and Sector {sb}\n",
         f"`python scripts/compare_sectors.py --sectors {sa} {sb}`. "
         f"{len(common):,} stars were searched in both sectors. Full table: `{os.path.relpath(stem, ROOT)}.csv`.\n",
         "A match means the same TIC has a dip at that level in both sectors. \"Consistent\" "
         "means the best dips agree in depth and duration within a factor of 2. \"By chance\" "
         "is the number of shared stars expected if dips fell on stars independently of each "
         "other.\n",
         "| level | stars with dips, S%d | stars with dips, S%d | in both | by chance | consistent |" % (sa, sb),
         "|---|---|---|---|---|---|"]
    for r in summary:
        L.append(f"| {r['level']} | {r['stars_a']:,} | {r['stars_b']:,} | **{r['both']}** | "
                 f"{r['expected_by_chance']:.1f} | {r['consistent']} |")
    L.append("")
    L.append("Two things push the \"dips\" level well above chance: some stars dip in every "
             "sector (variable stars, eclipsing binaries, persistent systematics), and some "
             "pixels are affected in both sectors. The later levels are the meaningful ones.\n")
    for lev in ("candidates", "vetted", "shortlist"):
        m = allm[allm.level == lev].sort_values(["consistent", "snr_a"], ascending=[False, False])
        L.append(f"## {lev.capitalize()} in both sectors ({len(m)})\n")
        if m.empty:
            L.append("None.\n")
            continue
        L.append(f"| TIC | S{sa} t0 | S{sa} depth (ppm) / dur (h) / SNR | S{sb} t0 | "
                 f"S{sb} depth / dur / SNR | depth ratio | dur ratio | consistent | Δt (d) | "
                 f"Phase 4 (S{sa} / S{sb}) |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for _, r in m.head(60).iterrows():
            p4 = f"{r.phase4_a if isinstance(r.phase4_a, str) else '–'} / " \
                 f"{r.phase4_b if isinstance(r.phase4_b, str) else '–'}"
            L.append(f"| {r.tic} | {r.t0_a:.2f} | {r.depth_ppm_a:.0f} / {r.dur_h_a:.1f} / {r.snr_a:.1f} | "
                     f"{r.t0_b:.2f} | {r.depth_ppm_b:.0f} / {r.dur_h_b:.1f} / {r.snr_b:.1f} | "
                     f"{r.depth_ratio:.2f} | {r.dur_ratio:.2f} | {'**yes**' if r.consistent else 'no'} | "
                     f"{r.dt_d:.1f} | {p4} |")
        if len(m) > 60:
            L.append(f"\n({len(m) - 60} more in the CSV.)")
        L.append("")
    with open(stem + ".md", "w") as fh:
        fh.write("\n".join(L))
    print("\n".join(L[:12]))
    print(sm.to_string(index=False))


if __name__ == "__main__":
    main()
