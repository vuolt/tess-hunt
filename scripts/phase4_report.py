"""Phase 4 ranking, CTOI summaries, follow-up sheets and summary tables.

Called by ``python scripts/phase4.py rank|report``.

Ranking rules (fixed before the Phase 4 results were looked at):

DROP if any of
  D1 independent photometry of the same sector (SPOC 2-min or QLP) should
     have seen the dip at SNR >= 5 but measures < 30 % of its depth;
  D2 Phase 3 flagged the dip as very long (T14 > 30 h) AND V-shaped or with
     a significant centroid offset: the signature of variability/systematics;
  D3 the host shows binarity (RUWE > 1.4, Gaia NSS, image doubling, RV
     variability or a close similar-brightness WDS pair) AND the companion
     radius, corrected for that dilution, exceeds 2 R_J;
  D4 vetted second dips exist but no period is consistent with all of them;
  D5 the "dip" is a flux step: in the undetrended TESScut pixels AND in every
     other light curve covering it, the in-transit level continues the
     baseline on one side (vetting.one_sided_check). Added after rules D1-D4
     and S1-S5 had been applied, when TIC 27068699 turned out to be a
     single-pixel jump; checked against the validation planets before use
     (results/phase4/confirm_validation.csv).

SUBMIT (as CTOI) if not dropped and all of
  S1 Phase 3 FPP < 0.1 and no Phase 3 review notes other than
     "difference image too faint for a centroid";
  S2 no binarity flag on the host (resolved wide companions are fine);
  S3 implied radius (dilution-corrected) <= 1.8 R_J;
  S4 the original dip is confirmed in independent photometry where that
     photometry is sensitive enough (expected SNR >= 5), or cannot be tested;
  S5 no other TESS sector shows a consistent-depth dip that fails vetting
     (a recurrent systematic on this star).

MAYBE otherwise.
"""

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from phase4 import jdump, jload, paths  # noqa: E402
from tesshunt import multisector as ms  # noqa: E402
from tesshunt.plotting import binned  # noqa: E402

RSUN_RJUP = 9.731
RJUP_REARTH = 11.209
CONFIRM_SNR = 5.0
CONFIRM_MIN_RATIO = 0.3
R_DROP = 2.0
R_SUBMIT = 1.8
FPP_SUBMIT = 0.1
BENIGN_NOTES = ("difference image too faint for a centroid",)


def _load(P, stage, tic, ext="json"):
    p = os.path.join(P["work"], stage, f"{tic}.{ext}")
    if not os.path.exists(p):
        return None
    return jload(p) if ext == "json" else np.load(p)


def joint_periods(periods, t0, t14):
    """Periods allowed by the exclusion map AND consistent with every vetted
    second dip (each at an integer number of cycles from t0)."""
    duos = periods.get("duos") or []
    if not duos:
        return None
    cands = [a["period"] for a in periods.get("allowed_aliases", []) if a["sector"] == duos[0]["sector"]]
    out = []
    for p in cands:
        ok = True
        for d in duos[1:]:
            n = (d["t0"] - t0) / p
            if abs(n - round(n)) * p > max(t14, 0.2):
                ok = False
        if ok:
            out.append(p)
    return out


def gather(c, P, sector):
    tic = c["tic"]
    cov = _load(P, "coverage", tic)
    srch = _load(P, "search", tic)
    per = _load(P, "periods", tic)
    conf = _load(P, "confirm", tic)
    bn = _load(P, "binarity", tic)
    g = dict(c)
    g["sectors"] = cov["sectors"] if cov else []
    g["n_sectors"] = len(g["sectors"])
    kinds = [p["kind"] for p in (cov or {}).get("products", [])]
    g["products"] = {k: kinds.count(k) for k in ("spoc2min", "tess-spoc", "qlp", "tesscut")}
    g["n_no_data"] = kinds.count(None)
    # other-sector dips
    other = []
    if srch:
        for r in srch["sectors"]:
            for d in r.get("matched_dips", []):
                if d.get("is_original") or not d.get("consistent"):
                    continue
                vt = d.get("vetting") or {}
                other.append(dict(sector=r["sector"], kind=r["kind"], t0=d["t0"], snr=d["snr"],
                                  depth_ratio=d["depth_ratio"], passed_1_4=vt.get("passed_1_4"),
                                  partial=vt.get("partial", False),
                                  consistent_shape=vt.get("consistent_shape"),
                                  shape=(vt.get("shape") or {}).get("passed"),
                                  edge=vt.get("edge_ok", (vt.get("edge") or {}).get("passed")),
                                  pixels=(vt.get("pixels") or {}).get("passed")))
    g["other_dips"] = other
    g["duos"] = [d for d in other if (d["passed_1_4"] or d["partial"]) and d["consistent_shape"]]
    g["recurrent_failing"] = [d for d in other if not d["passed_1_4"] and not d["partial"]]
    g["periods"] = per
    g["joint_periods"] = joint_periods(per, c["t0"], c["t14"]) if per else None
    g["confirm"] = conf["measurements"] if conf else []
    g["binarity"] = bn
    # radius, with dilution by a close similar-brightness companion if flagged
    rp = np.sqrt(max(c["depth"], 0)) * c["rad"] * RSUN_RJUP if c["rad"] == c["rad"] else np.nan
    dil = 1.0
    if bn and "wds_close" in bn.get("flags", []):
        dil = 2.0          # equal-brightness pair: the transited star gives ~half the light
    g["rp_rj"], g["rp_rj_diluted"] = rp, rp * np.sqrt(dil)
    return g


def classify(g):
    reasons_drop, reasons_maybe, good = [], [], []
    notes = g.get("review_notes") or ""
    bn = g.get("binarity") or {}
    flags = bn.get("flags", [])
    host_bin = bn.get("host_binary", False)
    # D1 / S4: independent confirmation
    tested = [m for m in g["confirm"] if m.get("available") and m["kind"] in ("spoc2min", "qlp")
              and (m.get("expected_snr") or 0) >= CONFIRM_SNR]
    if tested:
        best = max(tested, key=lambda m: m.get("expected_snr") or 0)
        r = best.get("depth_ratio")
        if r is not None and r < CONFIRM_MIN_RATIO:
            reasons_drop.append(f"D1 not seen in {best['kind']} (depth ratio {r:.2f}, "
                                f"expected SNR {best['expected_snr']:.0f})")
        else:
            good.append(f"confirmed in {best['kind']} (depth ratio {r:.2f}, SNR {best['snr']:.1f})")
    # D5: a step, not a dip
    stepped = [m for m in g["confirm"] if m.get("available") and (m.get("step") or {}).get("one_sided") is not None]
    if (stepped and any(m["kind"] == "tesscut" for m in stepped) and len(stepped) >= 2
            and all(m["step"]["one_sided"] for m in stepped)):
        reasons_drop.append("D5 flux step, not a dip: in-transit level continues one side's baseline in "
                            + ", ".join(f"{m['kind']} ({m['step']['pre']:+.2f}/{m['step']['post']:+.2f} of "
                                        "the depth before/after)" for m in stepped))
    # D2
    if "very long dip" in notes and ("grazing" in notes or "centroid offset" in notes):
        reasons_drop.append("D2 very long and V-shaped/off-centre: variability or systematics signature")
    # D3
    if host_bin and g["rp_rj_diluted"] > R_DROP:
        reasons_drop.append(f"D3 binary host ({', '.join(flags)}) and Rp {g['rp_rj_diluted']:.1f} R_J "
                            f"> {R_DROP:g} after dilution")
    # D4
    if g["duos"] and g["joint_periods"] is not None and len(g["joint_periods"]) == 0:
        reasons_drop.append("D4 second dip(s) inconsistent with every allowed period")
    # submit criteria
    fpp = g.get("fpp")
    other_notes = [n.strip() for n in notes.split(";") if n.strip() and n.strip() not in BENIGN_NOTES]
    if fpp is None or fpp >= FPP_SUBMIT:
        reasons_maybe.append(f"S1 FPP {fpp:.2f} >= {FPP_SUBMIT:g}" if fpp is not None else "S1 no FPP")
    if other_notes:
        reasons_maybe.append("S1 review notes: " + "; ".join(other_notes))
    if host_bin:
        reasons_maybe.append("S2 host binarity: " + bn.get("summary", ""))
    if not (g["rp_rj_diluted"] <= R_SUBMIT):
        reasons_maybe.append(f"S3 Rp {g['rp_rj_diluted']:.2f} R_J > {R_SUBMIT:g}")
    if g["recurrent_failing"]:
        reasons_maybe.append("S5 consistent-depth dips failing vetting in sector(s) "
                             + ", ".join(str(d["sector"]) for d in g["recurrent_failing"]))
    if g["duos"]:
        jp = g["joint_periods"] or []
        good.append(f"second transit candidate(s) in sector(s) "
                    f"{', '.join(str(d['sector']) + (' (partial)' if d['partial'] else '') for d in g['duos'])}"
                    + (f"; periods consistent with all: {', '.join(f'{p:.2f}' for p in jp[:6])} d"
                       if jp else ""))
    if reasons_drop:
        cat = "drop"
    elif not reasons_maybe:
        cat = "submit"
    else:
        cat = "maybe"
    return cat, reasons_drop, [r for r in reasons_maybe if r], good


def epoch_uncertainty(t14_d, rp, snr):
    """Carter et al. (2008): sigma_tc ~ T/Q * sqrt(tau / 2T), tau = ingress time."""
    tau = t14_d * rp / (1 + rp)
    return float(t14_d / max(snr, 1) * np.sqrt(max(tau, 1e-4) / (2 * t14_d)))


def period_summary(g):
    per = g["periods"]
    if not per:
        return "no period constraints (search not run)"
    floor = per.get("longest_excluded_below")
    iv = per.get("intervals", [])
    short = [(lo, hi) for lo, hi in iv if hi < (floor or np.inf)]
    parts = [f"all P > {floor:.0f} d allowed (longer than TESS's span around the event)"
             if floor else "no fully allowed long-period range"]
    if short:
        parts.append(f"{len(short)} allowed windows below that, e.g. " +
                     ", ".join(f"{lo:.1f}-{hi:.1f}" for lo, hi in short[:5]) + " d")
    if g.get("p_min_circ") and g["p_min_circ"] == g["p_min_circ"]:
        parts.append(f"duration implies P >= {g['p_min_circ']:.0f} d for a central circular orbit")
    if g.get("joint_periods"):
        parts.append("periods consistent with the second dip(s): "
                     + ", ".join(f"{p:.2f}" for p in g["joint_periods"][:8]) + " d")
    return "; ".join(parts)


def rank(cands, sector):
    P = paths(sector)
    os.makedirs(P["out"], exist_ok=True)
    rows = []
    for c in cands.to_dict("records"):
        g = gather(c, P, sector)
        cat, drop, maybe, good = classify(g)
        g.update(category=cat, reasons_drop=drop, reasons_maybe=maybe, positives=good,
                 period_summary=period_summary(g))
        rows.append(g)
    order = {"submit": 0, "maybe": 1, "drop": 2}
    rows.sort(key=lambda g: (g["role"] != "candidate", order[g["category"]],
                             g.get("fpp") if g.get("fpp") is not None else 1, -g["snr"]))
    jdump(os.path.join(P["work"], "ranked.json"), rows)
    for g in rows:
        print(f"{g['role']:9s} TIC {g['tic']:>10} {g['category']:6s} "
              f"{'; '.join(g['reasons_drop'] or g['reasons_maybe'] or g['positives'])[:150]}")
    return rows


# ---------------------------------------------------------------- report

def ctoi_rows(rows, sector):
    out, md = [], []
    for g in rows:
        if g["role"] != "candidate" or g["category"] != "submit":
            continue
        et = epoch_uncertainty(g["t14"], g["rp"], g["snr"])
        dep_ppm = g["depth"] * 1e6
        rp_re = g["rp_rj_diluted"] * RJUP_REARTH
        bn = g["binarity"] or {}
        steps = [m for m in g["confirm"] if m.get("available") and (m.get("step") or {}).get("one_sided") is False]
        vet = (f"Passes 7-step vetting (tess-hunt Phase 3): transit preferred over the best "
               f"non-transit model ({g.get('best_alt')}) by dBIC {g.get('dbic_alt', float('nan')):.0f}; "
               f"centroid on target; dip present in raw pixels; no known asteroid; not a known "
               f"EB/TOI; TRICERATOPS FPP {g['fpp']:.3f}. Not a flux step: the in-transit level "
               f"is below both local baselines in {', '.join(m['kind'] for m in steps)} (Phase 4).")
        if g["snr"] < 10 or g["fpp"] > 0.05:
            vet += (f" Caveat: marginal (SNR {g['snr']:.1f}, FPP {g['fpp']:.3f}); "
                    "follow-up photometry should precede heavy investment.")
        conf = "; ".join(g["positives"])
        per = g["period_summary"]
        out.append({
            "TIC ID": g["tic"], "Sector": sector,
            "Transit Epoch (BJD)": round(g["t0"] + 2457000.0, 5),
            "Transit Epoch error (d)": round(et, 5),
            "Period (days)": 0, "Period note": per,
            "Depth (ppm)": round(dep_ppm), "Depth error (ppm)": round(dep_ppm / g["snr"]),
            "Duration (hrs)": round(g["t14"] * 24, 2),
            "Duration error (hrs)": round(g["t14"] * 24 / g["snr"] * 2, 2),
            "Planet Radius (R_Earth)": round(rp_re, 1),
            "Impact parameter": round(g["b"], 2),
            "Stellar radius (R_Sun, TIC)": round(g["rad"], 3), "Tmag": round(g["Tmag"], 2),
            "RUWE": bn.get("ruwe"), "Notes": f"Single transit. {vet} {conf}".strip(),
        })
        md.append(
            f"### TIC {g['tic']}  (Tmag {g['Tmag']:.2f}, R* {g['rad']:.2f} Rsun, Teff "
            f"{g['Teff']:.0f} K)\n\n"
            f"- **Epoch:** BJD_TDB {g['t0'] + 2457000:.4f} ± {et:.4f} (TESS S{sector})\n"
            f"- **Depth:** {dep_ppm:.0f} ± {dep_ppm / g['snr']:.0f} ppm; **duration:** "
            f"{g['t14'] * 24:.2f} h; b ≈ {g['b']:.2f}; **Rp ≈ {rp_re:.1f} R⊕ "
            f"({g['rp_rj_diluted']:.2f} R_J)**\n"
            f"- **Period constraints:** {per}\n"
            f"- **Vetting:** {vet}\n"
            f"- **Binarity:** RUWE {bn.get('ruwe', float('nan')):.2f}; "
            f"{bn.get('summary', 'n/a')}\n"
            f"- **Other TESS data:** {g['n_sectors']} sectors observed "
            f"({', '.join(f'{k} {n}' for k, n in g['products'].items() if n)}); "
            f"{conf or 'no second transit found'}\n")
    return pd.DataFrame(out), "\n".join(md)


def period_plot(ax, g, P):
    z = _load(P, "periods", g["tic"], "npz")
    if z is None:
        ax.text(0.5, 0.5, "no period scan", transform=ax.transAxes, ha="center")
        return
    Pp = z["P"]
    per = g["periods"] or {}
    # Drawn from the full-resolution intervals, widened to stay visible on a
    # log axis: many allowed windows are narrower than a pixel.
    ax.axvspan(Pp[0], Pp[-1] * 1.3, color="0.8", lw=0, label="excluded by TESS data")
    for k, (lo, hi) in enumerate(per.get("intervals") or []):
        mid = np.sqrt(lo * hi)
        w = max(hi / lo, 1.004) ** 0.5
        ax.axvspan(mid / w, mid * w, color="C2", alpha=0.8, lw=0, label="allowed" if k == 0 else None)
    floor = per.get("longest_excluded_below")
    if floor:
        ax.axvline(floor, color="k", ls=":", lw=1)
        ax.text(floor, 1.02, f" all P > {floor:.0f} d allowed", fontsize=7, va="bottom")
    if g.get("p_min_circ") and g["p_min_circ"] == g["p_min_circ"] and g["p_min_circ"] > Pp[0]:
        ax.axvline(g["p_min_circ"], color="C1", ls="--", lw=1.2,
                   label=f"min P if circular (duration, b=0): {g['p_min_circ']:.0f} d")
    for p in (g.get("joint_periods") or [])[:10]:
        ax.axvline(p, color="C3", lw=1.5, label="period from second transit(s)" if p == g["joint_periods"][0] else None)
    ax.set_xscale("log")
    ax.set_xlim(Pp[0], Pp[-1] * 1.3)
    ax.set_ylim(0, 1.15)
    ax.set_yticks([])
    ax.set_xlabel("orbital period [d]")
    ax.legend(fontsize=7, loc="upper left")


def followup_sheet(g, P, sector, path):
    fig = plt.figure(figsize=(15, 10.5))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.55, 1.1, 1.0], hspace=0.45, wspace=0.25)
    # coverage timeline
    ax = fig.add_subplot(gs[0, :])
    srch = _load(P, "search", g["tic"]) or {"sectors": []}
    colors = {"spoc2min": "C0", "tess-spoc": "C2", "qlp": "C4", "tesscut": "C1", None: "0.85"}
    for r in srch["sectors"]:
        if r.get("t_start") is None:
            continue
        ax.axvspan(r["t_start"], r["t_end"], ymin=0.2, ymax=0.8, color=colors[r["kind"]], alpha=0.6)
        ax.text((r["t_start"] + r["t_end"]) / 2, 0.86, str(r["sector"]), ha="center", fontsize=6,
                transform=ax.get_xaxis_transform())
    ax.axvline(g["t0"], color="C3", lw=2)
    for d in g["other_dips"]:
        good = (d["passed_1_4"] or d["partial"]) and d["consistent_shape"]
        ax.axvline(d["t0"], color="k" if good else "0.6", lw=1.5, ls="-" if good else ":")
    for k, col in colors.items():
        if k:
            ax.plot([], [], "s", color=col, label=k)
    ax.plot([], [], color="C3", lw=2, label="original dip")
    ax.plot([], [], color="k", lw=1.5, label="vetted 2nd dip")
    ax.plot([], [], color="0.6", ls=":", label="consistent-depth dip failing vetting")
    ax.set_yticks([])
    ax.set_ylim(0, 1)
    ax.set_xlabel("BTJD")
    ax.legend(fontsize=7, ncol=7, loc="lower center", bbox_to_anchor=(0.5, 1.05))
    ax.set_title(f"TIC {g['tic']}: {g['category'].upper()}  |  {g['n_sectors']} sectors observed",
                 fontsize=11, pad=22)
    # original dip in independent photometry
    ax = fig.add_subplot(gs[1, 0])
    shown = 0
    for kind, col in (("tess-spoc", "C2"), ("qlp", "C4"), ("spoc2min", "C0")):
        p = ms.best_lc(g["tic"], g["ra"], g["dec"], g["Tmag"], sector, kinds=(kind,))
        if p is None:
            continue
        t, f = ms.candidate_detrend(p.lc, g["t14"])
        w = max(3 * g["t14"], 0.5)
        sel = np.abs(t - g["t0"]) < w
        if sel.sum() > 5:
            tb, fb = binned(t[sel], f[sel], max(g["t14"] / 8, 0.5 / 24))
            ax.plot((tb - g["t0"]) * 24, fb - 0.6 * shown * g["depth"] * 3, "o-", ms=3,
                    color=col, label=f"{kind}" + (" (offset)" if shown else ""))
            shown += 1
    ax.axvspan(-g["t14"] * 12, g["t14"] * 12, color="C3", alpha=0.1)
    ax.set_xlabel("hours from T0")
    ax.set_title(f"S{sector} dip in each product", fontsize=9)
    ax.legend(fontsize=7)
    # second dips
    duo_axes = [fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])]
    others = sorted(g["other_dips"], key=lambda d: (not (d["passed_1_4"] or d["partial"]),
                                                    -d["snr"]))[:2]
    for ax, d in zip(duo_axes, others):
        p = ms.best_lc(g["tic"], g["ra"], g["dec"], g["Tmag"], d["sector"], kinds=(d["kind"],))
        t, f = ms.candidate_detrend(p.lc, g["t14"])
        w = max(3 * g["t14"], 0.5)
        sel = np.abs(t - d["t0"]) < w
        ax.plot((t[sel] - d["t0"]) * 24, f[sel], ".", ms=2, color="0.6")
        tb, fb = binned(t[sel], f[sel], max(g["t14"] / 8, 0.5 / 24))
        ax.plot((tb - d["t0"]) * 24, fb, "o", ms=3, color="C0")
        ax.axhline(1 - g["depth"], color="C3", ls="--", lw=1, label="original depth")
        ax.set_title(f"S{d['sector']} ({d['kind']}) dip at {d['t0']:.2f}: SNR {d['snr']:.1f}, "
                     f"depth x{d['depth_ratio']:.2f}\nchecks 1-4 "
                     f"{'PASS' if d['passed_1_4'] else ('PARTIAL' if d['partial'] else 'FAIL')}"
                     f" (shape {d['shape']}, edge "
                     f"{d['edge']}, pixels {d['pixels']})", fontsize=8)
        ax.set_xlabel("hours from dip")
        ax.legend(fontsize=7)
    for ax in duo_axes[len(others):]:
        ax.set_axis_off()
        ax.text(0.5, 0.5, "no other consistent-depth dip\nin any TESS sector", ha="center",
                va="center", transform=ax.transAxes, fontsize=10, color="0.4")
    # periods
    ax = fig.add_subplot(gs[2, :2])
    period_plot(ax, g, P)
    ax.set_title("Orbital periods still allowed by all TESS coverage", fontsize=9)
    # text
    ax = fig.add_subplot(gs[2, 2])
    ax.set_axis_off()
    bn = g["binarity"] or {}
    lines = [f"Decision: {g['category'].upper()}", ""]
    lines += [f"- {r}" for r in g["reasons_drop"] + g["reasons_maybe"]]
    lines += [f"+ {r}" for r in g["positives"]]
    lines += ["", f"Binarity: RUWE {bn.get('ruwe') or float('nan'):.2f}",
              bn.get("summary", "n/a"), "",
              f"Rp ~ {g['rp_rj_diluted']:.2f} R_J; Phase 3 FPP "
              f"{g['fpp'] if g['fpp'] is not None else float('nan'):.3f}"]
    import textwrap
    txt = "\n".join(textwrap.fill(ln, 60) if ln else "" for ln in lines)
    ax.text(0, 1, txt, va="top", fontsize=7.5, family="monospace", transform=ax.transAxes)
    fig.savefig(path, dpi=80, bbox_inches="tight")
    plt.close(fig)


def report(cands, sector):
    P = paths(sector)
    rows = rank(cands, sector)
    os.makedirs(os.path.join(P["plots"], "sheets"), exist_ok=True)
    os.makedirs(os.path.join(P["plots"], "periods"), exist_ok=True)
    table = []
    for i, g in enumerate(rows, 1):
        tag = f"{i:02d}_tic{g['tic']}" if g["role"] == "candidate" else f"control_tic{g['tic']}"
        followup_sheet(g, P, sector, os.path.join(P["plots"], "sheets", f"{tag}_followup.png"))
        vetting_sheet(g, P, sector, os.path.join(P["plots"], "sheets", f"{tag}_vetting.png"))
        fig, ax = plt.subplots(figsize=(9, 2.8))
        period_plot(ax, g, P)
        ax.set_title(f"TIC {g['tic']}: allowed periods ({g['category']})", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(P["plots"], "periods", f"{tag}.png"), dpi=100)
        plt.close(fig)
        bn = g["binarity"] or {}
        per = g["periods"] or {}
        table.append(dict(
            order=i, tic=g["tic"], role=g["role"], category=g["category"], rank_phase3=g.get("rank3"),
            tier=g["tier"], tmag=g["Tmag"], rstar=g["rad"], t0_btjd=g["t0"],
            depth_ppm=g["depth"] * 1e6, t14_h=g["t14"] * 24, snr=g["snr"], fpp=g["fpp"],
            rp_rj=g["rp_rj"], rp_rj_diluted=g["rp_rj_diluted"], ruwe=bn.get("ruwe"),
            binarity_flags=";".join(bn.get("flags", [])), n_sectors=g["n_sectors"],
            **{f"n_{k}": n for k, n in g["products"].items()},
            n_second_dips=len(g["duos"]), n_failing_dips=len(g["recurrent_failing"]),
            p_floor_d=per.get("longest_excluded_below"),
            allowed_fraction_logP=per.get("allowed_fraction_log"),
            joint_periods=";".join(f"{p:.3f}" for p in (g["joint_periods"] or [])),
            confirm=";".join(f"{m['kind']}:{m.get('depth_ratio') or float('nan'):.2f}/"
                             f"{m.get('snr') or float('nan'):.1f}"
                             for m in g["confirm"] if m.get("available") and m["kind"] != "tess-spoc"),
            step=";".join(f"{m['kind']}:{m['step']['pre']:+.2f}/{m['step']['post']:+.2f}"
                          f"{'*' if m['step']['one_sided'] else ''}"
                          for m in g["confirm"] if m.get("available") and m.get("step")
                          and m["step"]["one_sided"] is not None),
            reasons="; ".join(g["reasons_drop"] + g["reasons_maybe"]),
            positives="; ".join(g["positives"]), period_summary=g["period_summary"]))
    df = pd.DataFrame(table)
    df.to_csv(os.path.join(P["out"], "followup.csv"), index=False, float_format="%.8g")
    ctoi, md = ctoi_rows(rows, sector)
    ctoi.to_csv(os.path.join(P["out"], "ctoi_candidates.csv"), index=False)
    with open(os.path.join(P["out"], "ctoi_summaries.md"), "w") as fh:
        fh.write(f"# CTOI-ready summaries (TESS Sector {sector} single transits)\n\n"
                 "Prepared for submission to ExoFOP as community TOIs; not submitted.\n\n" + md)
    cand = df[df.role == "candidate"]
    summ = dict(candidates=len(cand), categories=cand.category.value_counts().to_dict(),
                with_second_dip=int((cand.n_second_dips > 0).sum()),
                with_failing_dips=int((cand.n_failing_dips > 0).sum()),
                binarity_flagged=int(cand.binarity_flags.str.len().gt(0).sum()),
                controls=df[df.role == "control"].to_dict("records"))
    with open(os.path.join(P["out"], "phase4_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in summ.items() if k != "controls"}, indent=1))


def vetting_sheet(g, P, sector, path):
    """The Phase 3 vetting sheet, with Gaia DR3 binarity added to the check table."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import phase3_report as p3
    va = pd.read_csv(os.path.join(P["p3"], "vetting_all.csv"))
    row = va[va.key == g["key"]]
    if row.empty:
        return
    r = row.iloc[0].to_dict()
    for m in ("transit", "box", "ramp", "step", "flare_decay"):
        r.setdefault(f"bic_{m}", None)
    bn = g["binarity"] or {}
    r["binarity_line"] = (f"RUWE {bn.get('ruwe') or float('nan'):.2f}; "
                          + bn.get("summary", "n/a"))
    r["binarity_ok"] = not bn.get("host_binary", False)
    p3.sheet(r, path)
