"""Phase 2 report: export, vet, rank and plot the results of scripts/phase2.py.

    python scripts/phase2_report.py              # everything
    python scripts/phase2_report.py --no-top     # skip re-downloading for the top-N plots

Outputs
  results/phase2/s00XX_stars.csv.gz       every star: status, variability, window, flags
  results/phase2/s00XX_dips.csv           every dip, ranked by SNR, with vetting flags
  results/phase2/s00XX_candidates.csv     dips that survive the automated vetting
  results/phase2/s00XX_injections.csv     injection-recovery trials
  results/phase2/s00XX_summary.json       numbers quoted in RESULTS.md
  plots/phase2/*.png                      summary figures
  plots/phase2/top100/*.png               diagnostics for the 100 highest-SNR dips
"""

import argparse
import json
import os
import sqlite3
import sys
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import poisson

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tesshunt import ffi  # noqa: E402
from tesshunt.detect import Event  # noqa: E402
from tesshunt.plotting import binned, event_figure  # noqa: E402
from tesshunt.survey import analyse  # noqa: E402

CM_BIN_D = 0.25        # common-mode time bin
CM_PVALUE = 1e-4       # Poisson tail probability for a bin to count as common-mode
BTJD = 2457000.0


def load(db_path):
    db = sqlite3.connect(db_path)
    return (pd.read_sql("SELECT * FROM stars", db), pd.read_sql("SELECT * FROM dips", db),
            pd.read_sql("SELECT * FROM injections", db))


def common_mode_bins(t0, bin_d=CM_BIN_D, pvalue=CM_PVALUE):
    """Time bins holding far more dips (from different stars) than typical.

    Returns (edges, counts, flagged_bin_mask, expected_rate)."""
    edges = np.arange(np.floor(t0.min()), np.ceil(t0.max()) + bin_d, bin_d)
    counts, _ = np.histogram(t0, edges)
    lam = float(np.median(counts[counts > 0]))
    flagged = poisson.sf(counts - 1, lam) < pvalue
    return edges, counts, flagged, lam


def fetch_tois(tics):
    q = ("select tid,toi,tfopwg_disp,pl_orbper,pl_tranmid,pl_trandurh,pl_trandep from toi")
    url = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
           + urllib.parse.quote(q) + "&format=csv")
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            tois = pd.read_csv(r)
    except Exception as e:  # noqa: BLE001
        print(f"TOI table unavailable ({e}); skipping cross-match")
        return pd.DataFrame(columns=["tid", "toi", "tfopwg_disp", "pl_orbper",
                                     "pl_tranmid", "pl_trandurh", "pl_trandep"])
    return tois[tois.tid.isin(set(tics))]


def match_toi(dips, tois):
    """Name of a TOI whose ephemeris predicts a transit at the dip time."""
    out = np.full(len(dips), "", dtype=object)
    host = np.full(len(dips), "", dtype=object)
    by_tic = {k: g for k, g in tois.groupby("tid")}
    for i, (tic, t0) in enumerate(zip(dips.tic.values, dips.t0.values)):
        g = by_tic.get(tic)
        if g is None:
            continue
        host[i] = ",".join(str(x) for x in g.toi)
        for _, r in g.iterrows():
            if not (np.isfinite(r.pl_orbper) and np.isfinite(r.pl_tranmid)) or r.pl_orbper <= 0:
                continue
            tm = r.pl_tranmid - BTJD
            n = np.round((t0 - tm) / r.pl_orbper)
            tol = max((r.pl_trandurh if np.isfinite(r.pl_trandurh) else 3) / 24, 0.1)
            if abs(t0 - (tm + n * r.pl_orbper)) < tol:
                out[i] = str(r.toi)
                break
    return out, host


def vet(stars, dips, tois):
    d = dips.merge(stars[["tic", "tmag", "teff", "radius", "window_d", "max_dur_h",
                          "long_limited", "hf_variable", "var_period_d", "var_amp_ppm",
                          "sigma_pt_ppm", "crowdsap"]], on="tic", how="left")
    edges, counts, flagged, lam = common_mode_bins(d.t0.values)
    b = np.clip(np.digitize(d.t0.values, edges) - 1, 0, len(counts) - 1)
    d["common_mode"] = flagged[b]
    d["near_edge"] = (d.edge_dist_h < d.duration_h)
    clean = ~d.common_mode & ~d.near_edge
    d["n_dips_star"] = d.groupby("tic").t0.transform("size")
    d["n_clean_star"] = clean.groupby(d.tic).transform("sum")
    d["toi_match"], d["toi_host"] = match_toi(d, tois)

    cat = np.where(d.common_mode, "common_mode",
          np.where(d.near_edge, "edge",
          np.where(d.n_clean_star >= 2, "multi", "single")))
    d["category"] = cat
    d["candidate"] = (d.category == "single") & (d.hf_variable == 0)
    d = d.sort_values("snr", ascending=False).reset_index(drop=True)
    d.insert(0, "rank", np.arange(1, len(d) + 1))
    cm = dict(edges=edges, counts=counts, flagged=flagged, lam=lam)
    return d, cm


# ---------------------------------------------------------------- figures

def fig_times(d, cm, path):
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    e, c, f = cm["edges"], cm["counts"], cm["flagged"]
    ax = axes[0]
    ax.bar(e[:-1], c, width=np.diff(e), align="edge", color=np.where(f, "C3", "0.5"))
    ax.axhline(cm["lam"], color="k", ls=":", lw=1, label=f"median {cm['lam']:.0f} per bin")
    ax.set_yscale("log")
    ax.set_ylabel(f"dips per {CM_BIN_D:g} d")
    ax.set_title("All dips: bins in red are common-mode (Poisson p < 1e-4), i.e. systematics")
    ax.legend(fontsize=8)
    ax = axes[1]
    cand = d[d.candidate]
    ax.hist(cand.t0, bins=e, color="C0")
    ax.set_ylabel("candidates")
    ax.set_xlabel("Time [BTJD]")
    ax.set_title("Single-dip candidates after vetting")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def fig_dips(d, path):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    ax = axes[0]
    bins = np.logspace(np.log10(7), np.log10(max(d.snr.max(), 8)), 40)
    for c in ["common_mode", "edge", "multi", "single"]:
        ax.hist(d.snr[d.category == c], bins=bins, histtype="step", lw=1.5,
                label=f"{c} ({(d.category == c).sum()})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("SNR")
    ax.set_ylabel("dips")
    ax.legend(fontsize=8)
    ax = axes[1]
    cand = d[d.candidate]
    jit = np.exp(np.random.default_rng(0).normal(0, 0.06, len(cand)))
    sc = ax.scatter(cand.duration_h * jit, cand.depth_ppm, c=np.log10(cand.snr), s=6,
                    cmap="viridis")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("box duration [h] (jittered)")
    ax.set_ylabel("depth [ppm]")
    ax.set_title("Candidates")
    fig.colorbar(sc, ax=ax, label="log10 SNR")
    ax = axes[2]
    ax.scatter(cand.tmag, cand.depth_ppm, s=4, c="C0", alpha=0.5, label="candidates")
    ax.set_yscale("log")
    ax.set_xlabel("Tmag")
    ax.set_ylabel("depth [ppm]")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def fig_variability(s, path):
    ok = s[s.status == "ok"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    ax = axes[0]
    per = ok[np.isfinite(ok.var_period_d)]
    col = np.where(per.hf_variable == 1, "C3", np.where(per.long_limited == 1, "C1", "0.6"))
    ax.scatter(per.var_period_d, per.var_amp_ppm, s=1, c=col, alpha=0.5, rasterized=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("variability period [d]")
    ax.set_ylabel("semi-amplitude [ppm]")
    ax.set_title("grey: full window, orange: long-limited, red: fast variable", fontsize=9)
    ax = axes[1]
    w, n = np.unique(ok.window_d, return_counts=True)
    ax.bar([f"{x:g}" for x in w], n, color="C0")
    ax.set_yscale("log")
    ax.set_xlabel("detrending window [d]")
    ax.set_ylabel("stars")
    for i, v in enumerate(n):
        ax.text(i, v, f"{v}", ha="center", va="bottom", fontsize=8)
    ax = axes[2]
    tb = np.arange(np.floor(ok.tmag.min()), 13.01, 0.5)
    idx = np.digitize(ok.tmag, tb)
    x, fl, fh = [], [], []
    for k in range(1, len(tb)):
        sel = idx == k
        if sel.sum() > 20:
            x.append((tb[k - 1] + tb[k]) / 2)
            fl.append(ok.long_limited[sel].mean())
            fh.append(ok.hf_variable[sel].mean())
    ax.plot(x, fl, "o-", label="long-limited (window < 3 d)")
    ax.plot(x, fh, "s-", label="fast variable (unfilterable)")
    ax.set_xlabel("Tmag")
    ax.set_ylabel("fraction of stars")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


SNR_BINS = np.array([0, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 45, 70, 1e9])


def recovery_curve(inj, key="expected_snr", bins=SNR_BINS):
    idx = np.digitize(inj[key], bins) - 1
    rows = []
    for k in range(len(bins) - 1):
        sel = idx == k
        if sel.sum():
            rows.append((bins[k], bins[k + 1], int(sel.sum()), float(inj.recovered[sel].mean()),
                         float(np.median(inj[key][sel]))))
    return pd.DataFrame(rows, columns=["lo", "hi", "n", "frac", "median"])


def _snr_axis(ax):
    from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator([2, 3, 5, 7, 10, 15, 20, 30, 50, 100]))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_locator(NullLocator())


def fig_sensitivity(inj, path):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    groups = [("all", inj),
              ("full window", inj[(inj.long_limited == 0) & (inj.hf_variable == 0)]),
              ("long-limited", inj[inj.long_limited == 1]),
              ("fast variable", inj[inj.hf_variable == 1])]
    ax = axes[0]
    for name, g in groups:
        if len(g) < 20:
            continue
        rc = recovery_curve(g)
        ax.plot(np.minimum(rc["median"], 100), rc.frac, "o-", label=f"{name} ({len(g)})")
    ax.axvline(7, color="C3", ls="--", lw=1)
    _snr_axis(ax)
    ax.set_xlabel("expected SNR (empirical noise)")
    ax.set_ylabel("recovery fraction")
    ax.set_ylim(-0.03, 1.03)
    ax.legend(fontsize=8, title="star type")
    ax.grid(alpha=0.3)
    ax = axes[1]
    base = groups[1][1]
    for h in sorted(base.duration_h.unique()):
        g = base[base.duration_h == h]
        rc = recovery_curve(g)
        ax.plot(np.minimum(rc["median"], 100), rc.frac, "o-", label=f"{h:g} h")
    ax.axvline(7, color="C3", ls="--", lw=1)
    _snr_axis(ax)
    ax.set_xlabel("expected SNR (empirical noise)")
    ax.set_title("full-window stars, by duration", fontsize=9)
    ax.set_ylim(-0.03, 1.03)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[2]
    rc = recovery_curve(inj, "expected_snr_white")
    ax.plot(np.minimum(rc["median"], 100), rc.frac, "o-", color="k", label="white-noise SNR")
    rc = recovery_curve(inj)
    ax.plot(np.minimum(rc["median"], 100), rc.frac, "o-", color="C0", label="empirical SNR")
    ax.axvline(7, color="C3", ls="--", lw=1)
    _snr_axis(ax)
    ax.set_xlabel("expected SNR")
    ax.set_title("all injections: red noise shifts the curve", fontsize=9)
    ax.set_ylim(-0.03, 1.03)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------- top-N plots

def star_search(tic, sector, tmp):
    path = os.path.join(tmp, f"plot_{tic}.fits")
    try:
        ffi.download(tic, sector, path)
        lc, _ = ffi.read(path)
    finally:
        if os.path.exists(path):
            os.remove(path)
    var, cfg, ss = analyse(lc)
    return var, cfg, ss


def plot_top(d, sector, outdir, tmp, n=100):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(tmp, exist_ok=True)
    top = d.head(n)
    for tic, g in top.groupby("tic", sort=False):
        try:
            var, cfg, ss = star_search(int(tic), sector, tmp)
        except Exception as e:  # noqa: BLE001
            print(f"  TIC {tic}: {e}")
            continue
        for _, r in g.iterrows():
            ev = Event(r.t0, r.duration_h / 24, r.depth_ppm * 1e-6, r.snr, sector)
            flags = [r.category] + (["TOI " + r.toi_match] if r.toi_match else []) \
                + (["long-limited"] if r.long_limited else []) \
                + (["fast-var"] if r.hf_variable else [])
            title = (f"#{r['rank']}  TIC {tic}  T={r.tmag:.1f}  S{sector}  "
                     f"window {var.window:g} d  [{', '.join(flags)}]")
            event_figure(ss, ev, title, os.path.join(outdir, f"r{r['rank']:03d}_tic{tic}.png"),
                         cfg.threshold, dpi=65)
    print(f"  top-{n} plots -> {outdir}")


def contact_sheet(d, sector, path, tmp, n=50):
    cand = d[d.candidate].head(n)
    cols = 5
    rows = int(np.ceil(len(cand) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(3.2 * cols, 2.3 * rows))
    axes = np.atleast_2d(axes)
    for ax, (_, r) in zip(axes.flat, cand.iterrows()):
        try:
            _, _, ss = star_search(int(r.tic), sector, tmp)
        except Exception:  # noqa: BLE001
            ax.set_axis_off()
            continue
        dur = r.duration_h / 24
        w = max(3 * dur, 0.5)
        sel = np.abs(ss.time - r.t0) < w
        ax.plot((ss.time[sel] - r.t0) * 24, ss.flat[sel], ".", ms=1.5, color="0.6")
        if sel.sum() > 10:
            tb, fb = binned(ss.time[sel], ss.flat[sel], max(dur / 6, 0.5 / 24))
            ax.plot((tb - r.t0) * 24, fb, "o", ms=2.5, color="C0")
        tt = np.linspace(-w, w, 500)
        ax.plot(tt * 24, np.where(np.abs(tt) < dur / 2, 1 - r.depth_ppm * 1e-6, 1), "C3", lw=1)
        ax.set_title(f"#{r['rank']} TIC {r.tic}\nSNR {r.snr:.0f}, {r.depth_ppm:.0f} ppm, "
                     f"{r.duration_h:g} h{' TOI ' + r.toi_match if r.toi_match else ''}",
                     fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in list(axes.flat)[len(cand):]:
        ax.set_axis_off()
    fig.suptitle(f"Top {len(cand)} vetted single-dip candidates, S{sector} (hours from dip centre)")
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=48)
    ap.add_argument("--no-top", action="store_true")
    ap.add_argument("--top", type=int, default=100)
    args = ap.parse_args()
    tag = f"s{args.sector:04d}"
    work = os.path.join(ROOT, "work")
    out = os.path.join(ROOT, "results", "phase2")
    plots = os.path.join(ROOT, "plots", "phase2")
    os.makedirs(out, exist_ok=True)
    os.makedirs(plots, exist_ok=True)

    stars, dips, inj = load(os.path.join(work, f"{tag}.sqlite"))
    ok = stars[stars.status == "ok"]
    print(f"{len(stars)} stars recorded, {len(ok)} searched, {len(dips)} dips, "
          f"{len(inj)} injections")

    tois = fetch_tois(stars.tic.values)
    d, cm = vet(stars, dips, tois)

    stars.sort_values("tic").to_csv(os.path.join(out, f"{tag}_stars.csv.gz"), index=False,
                                    float_format="%.6g")
    d.to_csv(os.path.join(out, f"{tag}_dips.csv"), index=False, float_format="%.6g")
    d[d.candidate].to_csv(os.path.join(out, f"{tag}_candidates.csv"), index=False,
                          float_format="%.6g")
    inj.to_csv(os.path.join(out, f"{tag}_injections.csv"), index=False, float_format="%.6g")

    fig_times(d, cm, os.path.join(plots, "detection_times.png"))
    fig_dips(d, os.path.join(plots, "dip_distributions.png"))
    fig_variability(stars, os.path.join(plots, "variability_windows.png"))
    if len(inj):
        fig_sensitivity(inj, os.path.join(plots, "sensitivity.png"))

    # Summary numbers
    sel_counts = {}
    sc_path = os.path.join(work, f"{tag}_selection_counts.csv")
    if os.path.exists(sc_path):
        sel_counts = dict(pd.read_csv(sc_path, header=None).values.tolist())
    tois_in_sample = tois[tois.tid.isin(ok.tic)]
    summ = dict(
        selection=sel_counts,
        recorded=len(stars),
        status=stars.groupby(["status", "reason"]).size().rename("n").reset_index()
        .to_dict("records"),
        searched=len(ok),
        long_limited=int(ok.long_limited.sum()),
        hf_variable=int(ok.hf_variable.sum()),
        both_flags=int(((ok.long_limited == 1) & (ok.hf_variable == 1)).sum()),
        periodic=int(np.isfinite(ok.var_period_d).sum()),
        window_counts={f"{k:g}": int(v) for k, v in ok.window_d.value_counts().sort_index().items()},
        median_sigma_pt_ppm=float(ok.sigma_pt_ppm.median()),
        median_red_noise_4h=float(ok.red_noise_4h.median()),
        stars_with_dips=int((ok.n_dips > 0).sum()),
        dips=len(d),
        category_counts=d.category.value_counts().to_dict(),
        candidates=int(d.candidate.sum()),
        candidate_stars=int(d[d.candidate].tic.nunique()),
        common_mode_bins=int(cm["flagged"].sum()),
        common_mode_rate_per_bin=cm["lam"],
        tois_in_sample=int(tois_in_sample.tid.nunique()),
        toi_matched_dips=int((d.toi_match != "").sum()),
        toi_matched_candidates=int(((d.toi_match != "") & d.candidate).sum()),
        toi2180=d[d.tic == 298663873][["rank", "t0", "depth_ppm", "duration_h", "snr",
                                        "category"]].to_dict("records"),
    )
    if len(inj):
        base = inj[(inj.long_limited == 0) & (inj.hf_variable == 0)]
        summ["injection"] = dict(
            stars=int(inj.tic.nunique()), trials=len(inj),
            curve_all=recovery_curve(inj).to_dict("records"),
            curve_full_window=recovery_curve(base).to_dict("records"),
            curve_white=recovery_curve(inj, "expected_snr_white").to_dict("records"),
            by_group={name: dict(n=len(g), frac_snr_gt_15=float(g[g.expected_snr > 15].recovered.mean())
                                 if (g.expected_snr > 15).any() else None)
                      for name, g in [("full_window", base),
                                      ("long_limited", inj[inj.long_limited == 1]),
                                      ("hf_variable", inj[inj.hf_variable == 1])]},
            long_limited_by_duration=inj[inj.long_limited == 1].groupby("duration_h")
            .recovered.agg(["size", "mean"]).reset_index().to_dict("records"),
            full_window_by_duration_snr_gt_15=base[base.expected_snr > 15].groupby("duration_h")
            .recovered.agg(["size", "mean"]).reset_index().to_dict("records"),
            false_detections_elsewhere=None,
        )
    with open(os.path.join(out, f"{tag}_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1, default=float)
    print(json.dumps({k: v for k, v in summ.items() if k != "injection"}, indent=1,
                     default=float))

    if not args.no_top:
        tmp = os.path.join(work, "tmp")
        plot_top(d, args.sector, os.path.join(plots, "top100"), tmp, n=args.top)
        contact_sheet(d, args.sector, os.path.join(plots, "top_candidates_sheet.png"), tmp)


if __name__ == "__main__":
    main()
