"""Phase 3 report: funnel, validation table, shortlist and vetting sheets.

Called by ``python scripts/phase3_vet.py report``.
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from phase3_vet import CHECKS, OUT, PLOTS, jpath, load_json  # noqa: E402
from tesshunt import vetting as v  # noqa: E402
from tesshunt.plotting import binned  # noqa: E402

LABEL = {"shape": "1 Shape", "duration": "2 Duration", "edge": "3 Edge",
         "pixels": "4 Pixels", "asteroid": "5 Asteroids", "catalogue": "6 Catalogues",
         "fpp": "7 FPP"}


def _get(path):
    return load_json(path) if os.path.exists(path) else None


def collect(tg):
    rows = []
    for r in tg.to_dict("records"):
        k = r["key"]
        lc, px, fp = _get(jpath("lc", k)), _get(jpath("pix", k)), _get(jpath("fpp", k))
        o = dict(key=k, tic=r["tic"], t0=r["t0"], tier=r.get("tier") or "", rank=r["rank"],
                 snr=r["snr"], depth_ppm=r["depth_ppm"], duration_h=r["duration_h"],
                 tmag=r["Tmag"], rstar=r["rad"], mstar=r["mass"], ra=r["ra"], dec=r["dec"],
                 is_candidate=r["is_candidate"], is_validation=r["is_validation"],
                 toi_phase2=r.get("toi_match") if isinstance(r.get("toi_match"), str) else
                 (f"{r['toi_match']:.2f}" if pd.notna(r.get("toi_match")) else ""))
        o["lc_status"] = lc["status"] if lc else "not run"
        if lc and lc["status"] == "ok":
            sh, du, ed, tr = lc["shape"], lc["duration"], lc["edge"], lc["shape"]["transit"]
            o.update(shape=sh["passed"], dbic_alt=sh["dbic_alt"], dbic_box=sh["dbic_box"],
                     best_alt=sh["best_alt"], fit_t0=tr["t0"], fit_rp=tr["rp"],
                     fit_t14_h=tr["t14_h"], fit_b=tr["b"], fit_depth_ppm=tr["depth_ppm"],
                     grazing=tr["grazing"], duration=du["passed"], t_min_h=du.get("t_min_h"),
                     p_b0_d=du.get("p_b0_d"), duration_note=du.get("note", ""),
                     edge=ed["passed"], near_edge=ed.get("near_edge"),
                     edge_resolved=ed.get("resolved"))
        o["pix_status"] = px["status"] if px else "not run"
        if px and px["status"] == "ok":
            p, c, n = px["pixels"], px["pixels"]["centroid"], px["pixels"]["neighbours"]
            o.update(pixels=p["passed"], centroid=c["passed"], centroid_offset_px=c["offset_px"],
                     centroid_sigma=c["offset_sigma"], diff_snr=c["diff_snr"],
                     centroid_inconclusive=c["inconclusive"], neighbours=n["passed"],
                     neighbours_failing=";".join(map(str, n["failing"])),
                     n_unresolved=len(n["unresolved"]), ap_snr=p["ap_snr"],
                     asteroid=px["asteroid"]["passed"],
                     asteroid_hits=";".join(h["name"] for h in px["asteroid"]["hits"]),
                     catalogue=px["catalogue"]["passed"], known=px["catalogue"]["known"],
                     eb=px["catalogue"]["eb"], toi_fp=px["catalogue"]["toi_fp"],
                     tois=";".join(px["catalogue"]["tois"]),
                     ctois=";".join(px["catalogue"]["ctois"]))
        o["fpp_status"] = fp["status"] if fp else "not run"
        if fp and fp["status"] == "ok":
            o.update(fpp=fp["fpp"]["passed"], fpp_value=fp["fpp"]["fpp"],
                     nfpp_value=fp["fpp"]["nfpp"], p_lo=fp["fpp"]["p_range"][0],
                     p_hi=fp["fpp"]["p_range"][1])
        rows.append(o)
    return pd.DataFrame(rows)


def funnel(df):
    """Sequential funnel over CHECKS; errors are counted separately."""
    alive = df.copy()
    out = [dict(step="Phase 2 candidates", remaining=len(alive), removed=0, errors=0)]
    for c in CHECKS:
        if c not in alive:
            break
        err = alive[c].isna()
        fail = alive[c].eq(False)
        out.append(dict(step=LABEL[c], removed=int(fail.sum()), errors=int(err.sum()),
                        remaining=int((alive[c] == True).sum())))  # noqa: E712
        alive = alive[alive[c] == True]  # noqa: E712
    return pd.DataFrame(out), alive


def independent(df):
    """How many fail each check when it is applied on its own (to the set the
    check was run on: 1-3 on every candidate, 4-6 on survivors of 1-3, 7 on
    survivors of 1-6)."""
    res = {}
    for c in CHECKS:
        if c in df:
            ran = df[c].notna()
            res[c] = dict(ran=int(ran.sum()), failed=int((df[c] == False).sum()))  # noqa: E712
    return res


# ---------------------------------------------------------------- sheets

def sheet(row, path):
    k = row["key"]
    lcz = np.load(jpath("lc", k, "npz"))
    pix = np.load(jpath("pix", k, "npz")) if os.path.exists(jpath("pix", k, "npz")) else None
    t0, t14 = row["fit_t0"], row["fit_t14_h"] / 24
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.8, 1.2, 1.3], hspace=0.35, wspace=0.25)

    ax = fig.add_subplot(gs[0, :])
    t, f = lcz["t"], lcz["f"]
    ax.plot(t, f, ",", color="0.6", rasterized=True)
    tb, fb = binned(t, f, 1 / 24)
    ax.plot(tb, fb, ".", ms=2, color="C0")
    ax.axvspan(t0 - t14 / 2, t0 + t14 / 2, color="C3", alpha=0.25)
    lo, hi = np.nanpercentile(f, [0.3, 99.7])
    ax.set_ylim(min(lo, 1 - 1.3 * row["fit_depth_ppm"] * 1e-6), hi)
    ax.set_xlabel("BTJD")
    ax.set_ylabel("detrended flux")
    known = f" | known: TOI {row['tois']}" if isinstance(row.get("tois"), str) and row["tois"] else ""
    ax.set_title(f"TIC {row['tic']}  T={row['tmag']:.2f}  R*={row['rstar']:.2f}  "
                 f"M*={row['mstar']:.2f}  | S48 dip T0={t0:.3f}  depth={row['fit_depth_ppm']:.0f} ppm  "
                 f"T14={row['fit_t14_h']:.1f} h  b={row['fit_b']:.2f}  SNR={row['snr']:.1f}  "
                 f"tier {row['tier'] or '-'}{known}", fontsize=10)

    ax = fig.add_subplot(gs[1, :2])
    tw, fw = lcz["tw"], lcz["fw"]
    x = (tw - t0) * 24
    ax.plot(x, fw, ".", ms=2, color="0.65")
    tb, fb = binned(tw, fw.astype(float), max(t14 / 8, 0.5 / 24))
    ax.plot((tb - t0) * 24, fb, "o", ms=3, color="C0", label="binned")
    styles = dict(transit=("C3", "-", 2.2), box=("C1", "--", 1.2), ramp=("C2", ":", 1.5),
                  step=("C4", "-.", 1.2), flare_decay=("C5", "--", 1.2))
    bic = {m: row.get(f"bic_{m}") for m in styles}
    for m, (c, ls, lw) in styles.items():
        if f"m_{m}" in lcz:
            lab = m if bic[m] is None else f"{m} (BIC {bic[m]:.0f})"
            ax.plot(x, lcz[f"m_{m}"], ls, color=c, lw=lw, label=lab)
    ax.set_xlabel(f"hours from {t0:.3f}")
    ax.set_ylabel("detrended flux")
    ax.legend(fontsize=7, ncol=3, loc="lower left")
    ax.set_title(f"Shape: dBIC(best alt: {row['best_alt']}) = {row['dbic_alt']:.1f}, "
                 f"dBIC(box) = {row['dbic_box']:.1f}", fontsize=9)

    ax = fig.add_subplot(gs[1, 2])
    if pix is not None:
        ax.plot((pix["tl"] - t0) * 24, pix["fl"], ".", ms=2, color="0.6")
        tb, fb = binned(pix["tl"], pix["fl"], max(t14 / 8, 0.5 / 24))
        ax.plot((tb - t0) * 24, fb, "o", ms=3, color="C0")
        ax.axvspan(-t14 * 12, t14 * 12, color="C3", alpha=0.15)
        ax.set_title(f"TESScut aperture ({int(pix['ap'].sum())} px): dip SNR "
                     f"{row.get('ap_snr', np.nan):.1f}", fontsize=9)
    ax.set_xlabel("hours from T0")

    if pix is not None:
        nbs = pix["nbs"]
        for j, (img, ttl) in enumerate([(pix["oot"], "out-of-transit image"),
                                        (pix["diff"] / pix["noise"], "difference image (SNR)")]):
            ax = fig.add_subplot(gs[2, j])
            if j == 0:
                im = ax.imshow(np.log10(np.clip(img, np.nanpercentile(img, 5), None)),
                               origin="lower", cmap="gray")
                yy, xx = np.nonzero(pix["ap"])
                for a, b in zip(xx, yy):
                    ax.add_patch(plt.Rectangle((a - .5, b - .5), 1, 1, fill=False, ec="C1", lw=1))
            else:
                vmax = np.nanmax(np.abs(img))
                im = ax.imshow(img, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
                if np.isfinite(row.get("centroid_offset_px", np.nan)):
                    c = load_json(jpath("pix", row["key"]))["pixels"]["centroid"]
                    ax.plot(c["xc"], c["yc"], "o", mfc="none", mec="k", ms=10, mew=1.5,
                            label="diff centroid")
            fig.colorbar(im, ax=ax, shrink=0.8)
            ax.plot(pix["x"], pix["y"], "+", color="C3", ms=14, mew=2, label="target")
            ax.add_patch(plt.Circle((pix["x"], pix["y"]), v.NEIGHBOUR_PX, fill=False,
                                    ec="C3", ls=":"))
            for _, nx, ny, dm in nbs:
                if dm < 6:
                    ax.plot(nx, ny, "x", color="yellow" if j == 0 else "k", ms=6)
                    ax.text(nx + 0.3, ny + 0.3, f"{dm:+.1f}", fontsize=6,
                            color="yellow" if j == 0 else "k")
            ax.set_xlim(-0.5, img.shape[1] - 0.5)
            ax.set_ylim(-0.5, img.shape[0] - 0.5)
            ax.set_title(ttl, fontsize=9)
            if j == 1:
                ax.legend(fontsize=7, loc="upper right")

    ax = fig.add_subplot(gs[2, 2])
    ax.set_axis_off()
    lines = check_lines(row)
    for i, (name, ok, detail) in enumerate(lines):
        y = 0.95 - i * 0.105
        col = {True: "C2", False: "C3", None: "0.5"}[ok]
        mark = {True: "PASS", False: "FAIL", None: "n/a"}[ok]
        ax.text(0.0, y, name, fontsize=9, transform=ax.transAxes, weight="bold")
        ax.text(0.42, y, mark, fontsize=9, transform=ax.transAxes, color=col, weight="bold")
        ax.text(0.0, y - 0.045, detail, fontsize=7, transform=ax.transAxes, color="0.3")
    fig.savefig(path, dpi=80, bbox_inches="tight")
    plt.close(fig)


def _b(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else bool(x)


def check_lines(r):
    f = lambda x, fmt: ("-" if x is None or (isinstance(x, float) and np.isnan(x))  # noqa: E731
                        else format(x, fmt))
    return [
        ("1 Shape", _b(r.get("shape")),
         f"dBIC vs {r.get('best_alt')} {f(r.get('dbic_alt'), '.1f')}, vs box {f(r.get('dbic_box'), '.1f')}"),
        ("2 Duration", _b(r.get("duration")),
         f"T14 {f(r.get('fit_t14_h'), '.1f')} h vs min {f(r.get('t_min_h'), '.1f')} h (P>20 d) "
         f"{r.get('duration_note') or ''}"),
        ("3 Edge", _b(r.get("edge")),
         f"near edge: {r.get('near_edge')}, resolved: {r.get('edge_resolved')}"),
        ("4a Centroid", _b(r.get("centroid")),
         f"offset {f(r.get('centroid_offset_px'), '.2f')} px ({f(r.get('centroid_sigma'), '.1f')} sig), "
         f"diff SNR {f(r.get('diff_snr'), '.0f')}"
         + (" [inconclusive]" if r.get("centroid_inconclusive") else "")),
        ("4b Neighbours", _b(r.get("neighbours")),
         f"stronger on: {r.get('neighbours_failing') or 'none'}; unresolved: {r.get('n_unresolved', '-')}"),
        ("5 Asteroids", _b(r.get("asteroid")), f"hits: {r.get('asteroid_hits') or 'none'}"),
        ("6 Catalogues", _b(r.get("catalogue")),
         f"TOI {r.get('tois') or '-'} CTOI {r.get('ctois') or '-'} EB {r.get('eb')}"),
        ("7 TRICERATOPS", _b(r.get("fpp")),
         f"FPP {f(r.get('fpp_value'), '.3f')}  NFPP {f(r.get('nfpp_value'), '.3f')}  "
         f"P {f(r.get('p_lo'), '.0f')}-{f(r.get('p_hi'), '.0f')} d"),
    ]


def funnel_figure(fn, path):
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.bar(fn.step, fn.remaining, color="C0")
    for i, (r, d) in enumerate(zip(fn.remaining, fn.removed)):
        ax.text(i, r, f"{r}" + (f"\n(-{d})" if d else ""), ha="center", va="bottom", fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel("candidates remaining")
    ax.set_title("Phase 3 vetting funnel (Sector 48 single-dip candidates)")
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def report(tg):
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(PLOTS, "sheets"), exist_ok=True)
    os.makedirs(os.path.join(PLOTS, "validation"), exist_ok=True)
    df = collect(tg)
    # BICs for the sheets
    for i, r in df.iterrows():
        lc = _get(jpath("lc", r["key"]))
        if lc and lc["status"] == "ok":
            for m, b in lc["shape"]["bic"].items():
                df.loc[i, f"bic_{m}"] = b - lc["shape"]["bic"]["transit"]
    df.to_csv(os.path.join(OUT, "vetting_all.csv"), index=False, float_format="%.5g")

    cand = df[df.is_candidate]
    fn, surv = funnel(cand)
    fn_tier = {t: funnel(cand[cand.tier == t])[0] for t in ("A", "B")}
    indep = independent(cand)
    funnel_figure(fn, os.path.join(PLOTS, "funnel.png"))

    surv = surv.copy()
    surv["new"] = ~surv.known.astype(bool)
    surv = surv.sort_values(["new", "fpp_value", "snr"], ascending=[False, True, False])
    surv.insert(0, "shortlist_rank", np.arange(1, len(surv) + 1))
    cols = ["shortlist_rank", "tic", "tier", "tmag", "rstar", "fit_t0", "fit_depth_ppm",
            "fit_t14_h", "fit_b", "snr", "dbic_alt", "dbic_box", "p_b0_d", "fpp_value",
            "nfpp_value", "p_lo", "p_hi", "centroid_offset_px", "n_unresolved", "tois", "ctois",
            "new"]
    surv[cols].to_csv(os.path.join(OUT, "shortlist.csv"), index=False, float_format="%.5g")
    for _, r in surv.iterrows():
        sheet(r, os.path.join(PLOTS, "sheets", f"{int(r.shortlist_rank):02d}_tic{r.tic}.png"))

    val = df[df.is_validation].copy()
    for _, r in val.iterrows():
        if r.lc_status == "ok" and r.pix_status == "ok":
            toi = r.toi_phase2 or "x"
            sheet(r, os.path.join(PLOTS, "validation", f"toi{toi}_tic{r.tic}.png"))
    vt = val[["tic", "toi_phase2", "tmag", "snr", "fit_depth_ppm", "fit_t14_h"] +
             [c for c in CHECKS if c in val] + ["fpp_value", "nfpp_value"]]
    vt.to_csv(os.path.join(OUT, "validation.csv"), index=False, float_format="%.4g")

    summ = dict(
        candidates=int(len(cand)),
        funnel=fn.to_dict("records"),
        funnel_by_tier={t: f.to_dict("records") for t, f in fn_tier.items()},
        independent=indep,
        survivors=int(len(surv)), survivors_new=int(surv.new.sum()),
        survivors_known=int((~surv.new).sum()),
        validation=dict(
            n=int(len(val)),
            per_check={c: dict(ran=int(val[c].notna().sum()),
                               passed=int((val[c] == True).sum()))  # noqa: E712
                       for c in CHECKS if c in val},
            all_passed=int(val[[c for c in CHECKS if c in val]].eq(True).all(axis=1).sum()),
            failures={r.toi_phase2: [c for c in CHECKS if c in val and r[c] == False]  # noqa: E712
                      for _, r in val.iterrows()}),
        errors={s: df[f"{s}_status"].value_counts().to_dict() for s in ("lc", "pix", "fpp")},
    )
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(fn.to_string(index=False))
    print(json.dumps({k: summ[k] for k in ("independent", "survivors", "survivors_new",
                                           "survivors_known", "validation", "errors")},
                     indent=1, default=str))
