"""Predict future transits of a duotransit candidate for both allowed periods.

    python scripts/predict_transits.py --tic 95747180 --events 45:2540.66 48:2621.1 \
        --n 1 2 --years 2

Fits both TESS transits jointly (shared shape, one mid-time each), propagates
the mid-time errors to every transit in the next ``--years`` years for
P = (t_b - t_a) / n, checks whether TESS observes the star then (tess-point
pointings, Sectors up to 134), and computes night-time visibility from the LCO
1-m network sites. Writes results/phase4/tic<TIC>_predictions.md (+ .csv) and
plots/phase4/tic<TIC>_duo_fit.png.
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from astropy.utils import iers  # noqa: E402

iers.conf.auto_download = False       # offline: IERS-B/erfa defaults are ample for minutes
iers.conf.auto_max_age = None

from astropy.time import Time  # noqa: E402

from tesshunt import ephemeris as ep  # noqa: E402
from tesshunt import multisector as ms  # noqa: E402

OUT = os.path.join(ROOT, "results", "phase4")
PLOTS = os.path.join(ROOT, "plots", "phase4")
VISIBLE_FRAC = 0.5         # a site sees at least this much of the transit + baseline
EDGE_WARN_D = 1.5          # predicted windows this close to an estimated sector edge are flagged


def utc(btjd, fmt="%Y-%m-%d %H:%M"):
    return Time(btjd + 2457000.0, format="jd", scale="tdb").utc.datetime.strftime(fmt)


def load_event(tic, star, sector, t0, t14):
    p = ms.best_lc(tic, star.ra, star.dec, star.Tmag, sector)
    t, f = ms.candidate_detrend(p.lc, t14)
    half = max(2.5 * t14, 0.5)
    sel = np.abs(t - t0) < half
    return dict(t=t[sel], f=f[sel], t0=t0, exp_time=p.cadence_min / 1440, sector=sector,
                kind=p.kind)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tic", type=int, required=True)
    ap.add_argument("--events", nargs=2, required=True, help="sector:t0 (BTJD) for both transits")
    ap.add_argument("--n", type=int, nargs="+", default=[1, 2], help="P = dt / n for each n")
    ap.add_argument("--years", type=float, default=2.0)
    ap.add_argument("--start", default=None, help="UTC start (default: now)")
    ap.add_argument("--sample", default=os.path.join(ROOT, "work", "s0048_sample.csv"))
    args = ap.parse_args()
    warnings.simplefilter("ignore")

    sm = pd.read_csv(args.sample)
    star = sm[sm.ID == args.tic].iloc[0]
    fu = pd.read_csv(os.path.join(OUT, "followup.csv"))
    row = fu[fu.tic == args.tic].iloc[0]
    t14_0, depth0 = row.t14_h / 24, row.depth_ppm * 1e-6
    evs = [load_event(args.tic, star, int(s), float(t), t14_0)
           for s, t in (e.split(":") for e in args.events)]
    fit = ep.joint_fit(evs, t14_0, np.sqrt(depth0))
    (ea, eb) = fit["events"]
    ta, sa, tb, sb = ea["t0"], ea["sigma"], eb["t0"], eb["sigma"]

    t_start = (Time(args.start) if args.start else Time.now()).tdb.jd - 2457000.0
    t_end = t_start + 365.25 * args.years
    bounds = ep.sector_bounds()
    secs = ep.tess_sectors(args.tic, star.ra, star.dec)
    last_known = max(bounds)

    rows = []
    for n in args.n:
        pr = ep.predict(ta, sa, tb, sb, n, t_start, t_end)
        for k, tc, s in zip(pr["k"], pr["tc"], pr["sigma"]):
            w_lo, w_hi = tc - 3 * s - fit["t14"] / 2, tc + 3 * s + fit["t14"] / 2
            tess = [x for x in secs if bounds[x][0] < w_hi and bounds[x][1] > w_lo]
            edge = [x for x in tess if min(abs(w_lo - bounds[x][0]), abs(w_hi - bounds[x][1]),
                                          abs(w_lo - bounds[x][1]), abs(w_hi - bounds[x][0])) < EDGE_WARN_D]
            r = dict(P_d=pr["P"], sigma_P_d=pr["sigma_P"], n=n, k=int(k), tc_btjd=tc,
                     tc_utc=utc(tc), sigma_h=s * 24, win_start_utc=utc(w_lo), win_end_utc=utc(w_hi),
                     win_h=(w_hi - w_lo) * 24,
                     tess_sector=";".join(map(str, tess)) if tess else "",
                     tess_near_edge=bool(edge), beyond_tess_plan=bool(w_lo > bounds[last_known][1]),
                     also_other_period=False)
            # nominal transit plus 1 h baseline, and the full +/-3 sigma window
            n_lo = tc - fit["t14"] / 2 - ep.BASELINE_H / 24
            n_hi = tc + fit["t14"] / 2 + ep.BASELINE_H / 24
            any_ok = None
            for hemi, sites in ep.SITES.items():
                best, union = (0.0, ""), None
                for name, site in sites.items():
                    fn, _, _ = ep.night_coverage(star.ra, star.dec, n_lo, n_hi, site)
                    _, ok, _ = ep.night_coverage(star.ra, star.dec, w_lo, w_hi, site)
                    union = ok if union is None else union | ok
                    if fn > best[0]:
                        best = (fn, name)
                any_ok = union if any_ok is None else any_ok | union
                r[f"{hemi}_site"] = best[1]
                r[f"{hemi}_transit_frac"] = best[0]
                r[f"{hemi}_window_frac"] = float(union.mean())
            r["network_window_frac"] = float(any_ok.mean())
            rows.append(r)
    df = pd.DataFrame(rows)
    # transits shared by both ephemerides (every other 40 d transit)
    if len(args.n) > 1:
        for i, r in df.iterrows():
            others = df[(df.n != r.n)]
            d = np.abs(others.tc_btjd - r.tc_btjd)
            df.loc[i, "also_other_period"] = bool((d < 3 * np.hypot(r.sigma_h, others.sigma_h) / 24).any())
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(os.path.join(OUT, f"tic{args.tic}_predictions.csv"), index=False, float_format="%.5f")
    _plot(fit, evs, args.tic)
    _markdown(args, star, row, fit, evs, df, secs, bounds, t_start, t_end)


def _plot(fit, evs, tic):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(evs), figsize=(6 * len(evs), 3.6), sharey=True)
    for ax, e, m, fe in zip(axes, evs, fit["models"], fit["events"]):
        x = (e["t"] - fe["t0"]) * 24
        ax.plot(x, e["f"], ".", color="0.6", ms=3)
        ax.plot(x, m, "C3", lw=1.5)
        ax.set_title(f"S{e['sector']} ({e['kind']}): t0 = {fe['t0']:.4f} ± {fe['sigma'] * 1440:.1f} min",
                     fontsize=9)
        ax.set_xlabel("hours from fitted mid-time")
    axes[0].set_ylabel("detrended flux")
    fig.suptitle(f"TIC {tic}: joint fit, depth {fit['depth'] * 1e6:.0f} ppm, "
                 f"T14 {fit['t14'] * 24:.2f} h, b {fit['b']:.2f}", fontsize=10)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(os.path.join(PLOTS, f"tic{tic}_duo_fit.png"), dpi=110)
    plt.close(fig)


def _yn(frac):
    return "full" if frac >= 0.999 else (f"{frac:.0%}" if frac > 0 else "–")


def _markdown(args, star, row, fit, evs, df, secs, bounds, t_start, t_end):
    ea, eb = fit["events"]
    L = []
    L.append(f"# TIC {args.tic}: predicted transits, {utc(t_start, '%Y-%m-%d')} to "
             f"{utc(t_end, '%Y-%m-%d')}\n")
    L.append(f"Generated by `python scripts/predict_transits.py --tic {args.tic} --events "
             f"{' '.join(args.events)} --n {' '.join(map(str, args.n))} --years {args.years:g}"
             + (f" --start {args.start}" if args.start else "") + "`.\n")
    L.append(f"**Star.** RA {star.ra:.5f}, Dec {star.dec:+.5f}, Tmag {star.Tmag:.2f}, "
             f"Teff {star.Teff:.0f} K, R★ {star.rad:.2f} R☉ (TIC).\n")
    L.append("## Transits in hand\n")
    L.append("Both events were fitted jointly with a limb-darkened transit model. The fit "
             "shares R_p/R★, T14 and b between the events, and gives each event its own "
             "mid-time and linear baseline. The mid-time errors are the larger of two "
             "estimates: the covariance-matrix error, and the scatter of the mid-times "
             "when the residuals are shifted cyclically and the model refitted (the "
             "\"prayer-bead\" method), which captures correlated noise.\n")
    L.append("| sector | product | mid-time (BTJD) | UTC | σ (formal) | σ (prayer bead) | σ used |")
    L.append("|---|---|---|---|---|---|---|")
    for e, fe in zip(evs, fit["events"]):
        L.append(f"| {e['sector']} | {e['kind']} | {fe['t0']:.5f} | {utc(fe['t0'])} | "
                 f"{fe['sigma_formal'] * 1440:.1f} min | {fe['sigma_prayer'] * 1440:.1f} min | "
                 f"**{fe['sigma'] * 1440:.1f} min** |")
    L.append("")
    L.append(f"Shared shape: depth {fit['depth'] * 1e6:.0f} ppm, R_p/R★ {fit['rp']:.3f} "
             f"(R_p ≈ {fit['rp'] * star.rad * 9.731:.2f} R_J), T14 {fit['t14'] * 24:.2f} h, "
             f"b {fit['b']:.2f}. Reduced χ² is {fit['chi2_red']:.2f}. The plot is "
             f"`plots/phase4/tic{args.tic}_duo_fit.png`.\n")
    L.append("## Ephemerides\n")
    L.append("Phase 4's period scan left only these two periods for the separation "
             "between the transits. All shorter aliases (Δt/n, n ≥ 3) are ruled out by TESS "
             "data.\n")
    L.append("| n | period (d) | σ_P |")
    L.append("|---|---|---|")
    for n in args.n:
        d = df[df.n == n].iloc[0]
        L.append(f"| {n} | {d.P_d:.4f} | {d.sigma_P_d * 1440:.1f} min |")
    L.append("")
    L.append("Every second transit of the 40 d ephemeris coincides with an 80 d transit. "
             "The others (**\"decisive\"** below) happen only if P ≈ 40 d.\n")
    L.append("- Seeing a transit at a decisive time settles P ≈ 40 d.")
    L.append("- Not seeing one only favours P ≈ 80 d. To rule out 40 d, the whole ±3σ "
             "window has to be covered, and even six sites together fall short of that "
             "(see the tables).\n")
    L.append("## TESS\n")
    fut = [s for s in secs if bounds[s][1] > t_start]
    L.append(f"tess-point (pointings to Sector {max(bounds)}, which ends around "
             f"{utc(bounds[max(bounds)][1], '%Y-%m-%d')}) puts the star on silicon in Sectors "
             f"{', '.join(map(str, secs))}. Of those, "
             + (f"Sectors {', '.join(map(str, fut))} are in the future." if fut else
                "none are in the future.") + " ")
    hit = df[df.tess_sector != ""]
    if len(hit):
        L.append("Predicted windows inside a future TESS sector:\n")
        for _, r in hit.iterrows():
            L.append(f"- P = {r.P_d:.2f} d, {r.tc_utc} UTC ± {r.sigma_h:.1f} h: Sector "
                     f"{r.tess_sector}" + (" (near an estimated sector edge)" if r.tess_near_edge else ""))
        L.append("")
    else:
        L.append("**No predicted window falls inside a TESS sector in this period, for either "
                 "period.** The period has to be settled from the ground.\n")
    L.append(f"Sector start and end times are estimated from tess-point's sector mid-times to "
             f"about ±1 d. Mid-sector downlink gaps are not modelled. "
             f"The last {max(0.0, (t_end - bounds[max(bounds)][1])):.0f} d of the interval lie "
             f"beyond the published pointing plan.\n")
    L.append("## Ground-based visibility\n")
    L.append("The criteria are the Sun below −12° and the star above 30° altitude "
             "(airmass < 2). The sites are the LCO 1-m network:\n")
    L.append("- **North:** Teide, McDonald, Haleakala.")
    L.append("- **South:** CTIO, SAAO, Siding Spring.\n")
    L.append("Two numbers are given per hemisphere. **transit + 1 h** is the fraction of "
             "the nominal transit, plus 1 h of baseline on each side, that the best single "
             "site can observe. **3σ window** is the fraction of the whole uncertainty "
             "window (t_c ± 3σ ± T14/2) that the hemisphere's three sites cover together. "
             "\"full\" means 100 %. With σ ≈ 3 h, no single night covers a 3σ window, so a "
             "campaign would combine sites (the last column combines all six).\n")
    for n in args.n:
        d = df[df.n == n]
        P = d.P_d.iloc[0]
        L.append(f"### P = {P:.3f} d ({len(d)} transits)\n")
        L.append("| # | mid-time UTC | ±1σ | 3σ window (UTC) | TESS | north: transit + 1 h | north: 3σ window | "
                 "south: transit + 1 h | south: 3σ window | all 6 sites: 3σ window | note |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for _, r in d.iterrows():
            note = []
            if n > min(args.n) and not r.also_other_period:
                note.append("**decisive**")
            if r.beyond_tess_plan:
                note.append("after S134")
            ns = r.north_site.split(" (")[0] if r.north_transit_frac > 0 else ""
            ss = r.south_site.split(" (")[0] if r.south_transit_frac > 0 else ""
            L.append(f"| {r.k} | {r.tc_utc} | {r.sigma_h:.1f} h | {r.win_start_utc[5:]} → "
                     f"{r.win_end_utc[5:]} | {r.tess_sector or '–'} | "
                     f"{_yn(r.north_transit_frac)}{' ' + ns if ns else ''} | {_yn(r.north_window_frac)} | "
                     f"{_yn(r.south_transit_frac)}{' ' + ss if ss else ''} | {_yn(r.south_window_frac)} | "
                     f"{_yn(r.network_window_frac)} | {', '.join(note)} |")
        L.append("")
    df["decisive"] = (df.n > min(args.n)) & ~df.also_other_period
    uniq = df.sort_values("n").drop_duplicates("tc_utc").sort_values("tc_btjd")
    good = uniq[(uniq.north_transit_frac >= VISIBLE_FRAC) | (uniq.south_transit_frac >= VISIBLE_FRAC)]
    L.append("## Windows visible at night\n")
    L.append(f"Each item below is a predicted transit where at least one site can observe "
             f"≥ {VISIBLE_FRAC:.0%} of the nominal transit plus baseline. Transits shared by "
             f"both periods are listed once, as \"both\".\n")
    for _, r in good.iterrows():
        per = f"P = {r.P_d:.2f} d only — **decisive**" if r.decisive else (
            "both periods" if r.also_other_period else f"P = {r.P_d:.2f} d")
        where = [f"{h}: {r[f'{h}_site']} {_yn(r[f'{h}_transit_frac'])} of transit + 1 h"
                 for h in ("north", "south") if r[f"{h}_transit_frac"] >= VISIBLE_FRAC]
        L.append(f"- **{r.tc_utc} UTC** ± {r.sigma_h:.1f} h ({per}). "
                 + "; ".join(where)
                 + f". The six sites together cover {_yn(r.network_window_frac)} of the 3σ window.")
    if not len(good):
        L.append("- none")
    L.append("")
    top = good[good.decisive].sort_values("network_window_frac", ascending=False).head(2)
    if len(top):
        L.append("**Reading this.** The best opportunities are the 40 d-only (\"decisive\") "
                 "transits of " + " and ".join(
                     f"**{r.tc_utc[:10]}** ({r.north_site if r.north_transit_frac >= r.south_transit_frac else r.south_site})"
                     for _, r in top.iterrows())
                 + ". A transit detected there confirms P ≈ 40 d. Covering the rest of the "
                 "window from other sites strengthens a non-detection.")
    L.append("")
    L.append("The mid-time uncertainty (about 3 h) comes from having only one cycle "
             "between the TESS transits. A single ground-based detection would pin the "
             "ephemeris to minutes. The transit is deep enough for a 1 m telescope: "
             f"{fit['depth'] * 1e6:.0f} ppm on a T = {star.Tmag:.1f} star.\n")
    L.append(f"The full table, including partial windows, is in "
             f"`results/phase4/tic{args.tic}_predictions.csv`.\n")
    path = os.path.join(OUT, f"tic{args.tic}_predictions.md")
    keep = ""
    if os.path.exists(path):                  # keep the archival section (archival_search.py)
        old = open(path).read()
        a, z = "<!-- archival:start -->", "<!-- archival:end -->"
        if a in old and z in old:
            keep = "\n" + old[old.index(a):old.index(z) + len(z)] + "\n"
    with open(path, "w") as fh:
        fh.write("\n".join(L) + keep)
    print("\n".join(L))


if __name__ == "__main__":
    main()
