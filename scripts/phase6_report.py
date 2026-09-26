"""Phase 6 report: validation of each check, flags, table and plain-English summary.

Called by `python scripts/phase6.py report`.

Gaia orbit check. A candidate is flagged when its transit period matches a
Gaia DR3 period (x1, x2, x3, 1/2, 1/3 within 3 sigma) AND either Gaia classifies
the system as an eclipsing binary, or the companion Gaia measures is a star
(>= 0.08 M_sun). A match with a substellar companion (e.g. WASP-18 b, TOI-503 b)
confirms the transiting object's orbit instead and is listed separately.

Light-curve checks (odd/even depths, secondary eclipse, centroid shift). Their
thresholds were fixed before the validation (phase6.VARIANTS). Rates are measured
on the random samples only (the Gaia-flagged candidates were selected by another
check). Per check, the most sensitive variant flagging <= 2 % of the confirmed/known
planets is used (else <= 5 %); a check above 5 % in every variant is not used.

Confidence:
  high    Gaia eclipsing-binary solution with a matching period; or a Gaia stellar
          companion >= 0.10 M_sun with a matching period; or >= 2 independent checks
  medium  Gaia companion 0.08-0.10 M_sun with a matching period; or one light-curve
          check far beyond its threshold (odd/even or centroid > 10 sigma, secondary
          SNR > 15 and >= 20 % of the transit depth)
  low     one light-curve check just beyond its threshold
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import phase6 as p6  # noqa: E402
from tesshunt import planetscore as ps  # noqa: E402

M_STAR = 0.08
M_HIGH = 0.10


def _load():
    W = p6.WORK
    t = pd.read_csv(os.path.join(W, "candidates.csv"), dtype={"gaia_id": "Int64"})
    m = pd.read_csv(os.path.join(W, "gaia_matches.csv"))
    nss = pd.read_csv(os.path.join(W, "nss.csv"))
    eb = pd.read_csv(os.path.join(W, "gaia_eb.csv"))
    p6.seed_from_ledger()
    lc = pd.DataFrame(list(p6.done_checks().values()))
    if len(lc):
        # use today's disposition (a candidate may have been confirmed or refuted since its check)
        grp = dict(zip(t.name, t.group))
        lc["group_at_check"] = lc.group
        lc = lc[lc.name.isin(grp)].copy()
        lc["group"] = lc.name.map(grp)
    t = pd.concat([t, ps.table_quantities(t)], axis=1)
    if len(lc):
        rows = t.set_index("name").reindex(lc.name).reset_index()
        lc = lc.reset_index(drop=True)
        lc["sec_excess"] = ps.secondary_excess(lc, rows).values
    return t, m, nss, eb, lc


def gaia_flags(m):
    """Per candidate: flagged / substellar / evidence text."""
    rows = []
    for name, g in m.groupby("name"):
        ev, flag, conf, sub = [], False, None, False
        for r in g.itertuples():
            k = {1.0: "", 2.0: " (transit period = 2x Gaia's)", 0.5: " (transit period = half Gaia's)",
                 3.0: " (3x)"}.get(round(r.factor, 3), f" (x{r.factor:.2f})")
            if r.eclipsing_solution:
                flag, conf = True, "high"
                ev.append(f"{r.source}: P = {r.gaia_period:.4f} d{k}")
            elif np.isfinite(r.m2) and r.m2 >= M_STAR:
                flag = True
                c = "high" if r.m2 >= M_HIGH else "medium"
                conf = "high" if conf == "high" else c
                ev.append(f"{r.source} orbit P = {r.gaia_period:.3f} d{k}; companion {r.m2:.2f} M_sun "
                          f"({r.m2_method})")
            elif np.isfinite(r.m2):
                sub = True
                ev.append(f"{r.source} orbit P = {r.gaia_period:.3f} d{k}; companion {r.m2 * 1047.6:.0f} M_Jup "
                          "(below 0.08 M_sun: Gaia may see the transiting object's own orbit)")
            else:
                ev.append(f"{r.source} orbit P = {r.gaia_period:.3f} d{k}; no mass estimate")
                if not flag:
                    flag, conf = True, "medium"
        rows.append(dict(name=name, tic=g.tic.iloc[0], group=g.group.iloc[0], gaia_flag=flag,
                         gaia_conf=conf if flag else None, gaia_substellar=sub and not flag,
                         gaia_evidence="; ".join(ev)))
    return pd.DataFrame(rows)


def _pick(rows, check):
    best = None
    for r in rows:
        if r["check"] != check:
            continue
        pr, fr = r["planet_rate"], r["fp_rate"]
        for cap, tier in ((p6.PREFERRED_PLANET_RATE, 0), (p6.MAX_PLANET_FLAG_RATE, 1)):
            if pr <= cap:
                key = (tier, -fr)
                if best is None or key < best[0]:
                    best = (key, r["variant"], r["_fn"], pr, fr)
                break
    for r in rows:
        if r["check"] == check:
            r["used"] = bool(best and r["variant"] == best[1])
    return best


def choose_table_variants(per):
    """The checks on catalogue numbers, validated on every known planet and
    false positive with a period."""
    rows, chosen = [], {}
    recs = {g: per[per.group == g].to_dict("records") for g in ("planet", "fp", "unresolved")}
    for check, variants in p6.TABLE_VARIANTS.items():
        for label, fn in variants:
            k = {g: int(sum(bool(fn(r)) for r in recs[g])) for g in recs}
            n = {g: len(recs[g]) for g in recs}
            rows.append(dict(check=check, variant=label, planets_flagged=k["planet"], planets=n["planet"],
                             planet_rate=k["planet"] / max(n["planet"], 1), fps_flagged=k["fp"], fps=n["fp"],
                             fp_rate=k["fp"] / max(n["fp"], 1), unresolved_flagged=k["unresolved"],
                             unresolved=n["unresolved"], _fn=fn))
        chosen[check] = _pick(rows, check)
    out = pd.DataFrame(rows).drop(columns="_fn")
    return out, chosen


def choose_variants(lc):
    """Flag rates of every variant on the random samples, and the variant used.
    The held-out samples are only reported, never used to choose."""
    rnd = lc[(lc["sample"] == "random") & (lc.n_events.fillna(0) > 0)]
    hold = lc[(lc["sample"] == "holdout") & (lc.n_events.fillna(0) > 0)] if len(lc) else lc
    rows, chosen = [], {}
    for check, variants in p6.VARIANTS.items():
        best = None
        for label, fn in variants:
            rate = {}
            for g in ("planet", "fp", "unresolved"):
                sub = rnd[rnd.group == g]
                n = int(sum(bool(fn(r)) for r in sub.to_dict("records")))
                rate[g] = (n, len(sub))
            pr = rate["planet"][0] / max(rate["planet"][1], 1)
            fr = rate["fp"][0] / max(rate["fp"][1], 1)
            hr = {}
            for g in ("planet", "fp"):
                sub = hold[hold.group == g] if len(hold) else hold
                hr[g] = (int(sum(bool(fn(r)) for r in sub.to_dict("records"))), len(sub))
            rows.append(dict(check=check, variant=label,
                             planets_flagged=rate["planet"][0], planets=rate["planet"][1], planet_rate=pr,
                             fps_flagged=rate["fp"][0], fps=rate["fp"][1], fp_rate=fr,
                             unresolved_flagged=rate["unresolved"][0], unresolved=rate["unresolved"][1],
                             holdout_planets_flagged=hr["planet"][0], holdout_planets=hr["planet"][1],
                             holdout_fps_flagged=hr["fp"][0], holdout_fps=hr["fp"][1]))
            for cap, tier in ((p6.PREFERRED_PLANET_RATE, 0), (p6.MAX_PLANET_FLAG_RATE, 1)):
                if pr <= cap:
                    key = (tier, -fr)
                    if best is None or key < best[0]:
                        best = (key, label, fn, pr, fr)
                    break
        chosen[check] = best
        for r in rows:
            if r["check"] == check:
                r["used"] = bool(best and r["variant"] == best[1])
    return pd.DataFrame(rows), chosen


def table_evidence(r, tchosen):
    flags, ev = [], []
    if tchosen.get("size") and tchosen["size"][2](r):
        flags.append("size")
        ev.append(f"implied size {r['rp_rj']:.1f} Jupiter radii: too big for a planet"
                  + (f" on a {r['period']:.0f}-day orbit" if r["period"] > 10 else ""))
    if tchosen.get("density") and tchosen["density"][2](r):
        flags.append("density")
        ev.append(f"the transit lasts {r['duration_ratio']:.1f}x longer than any orbit around a star of "
                  "the catalogued density allows (a blend, or a misclassified star)")
    return flags, ev


def planet_scores(t, lc, gf, flags):
    """Combined planet probability for every candidate with light-curve checks.

    Trained on the known planets and false positives of the random and
    held-out samples (never on the Gaia-selected sample, whose selection
    depends on one of the features); the quoted performance is from 5-fold
    cross-validation."""
    have = lc[lc.n_events.fillna(0) > 0] if len(lc) else lc
    if not len(have):
        return pd.DataFrame(), {}
    gflag = set(gf[gf.gaia_flag].name) if len(gf) else set()
    tt = t.set_index("name")
    train_names = have[have["sample"].isin(["random", "holdout"]) & have.group.isin(["planet", "fp"])].name
    tr = tt.loc[train_names].reset_index()
    y = (tr.group == "planet").astype(int).values
    if len(set(y)) < 2 or min(np.bincount(y)) < 20:
        return pd.DataFrame(), {}
    ftr = ps.features(tr, have, gflag)
    model, p_cv = ps.fit_score(ftr, y)
    order = np.argsort(-p_cv)
    top = order[: max(1, len(order) // 10)]
    info = dict(n_train=len(y), n_planets=int(y.sum()), auc=ps.auc(y, p_cv),
                top_decile_planet_frac=float(y[top].mean()),
                p90_planet_frac=float(y[p_cv >= 0.9].mean()) if (p_cv >= 0.9).any() else float("nan"),
                n_p90=int((p_cv >= 0.9).sum()),
                coef=dict(zip(ps.FEATURES, model[-1].coef_[0].round(2))))
    un = tt.loc[have[have.group == "unresolved"].name].reset_index()
    fu = ps.features(un, have, gflag)
    un["p_planet"] = model.predict_proba(fu[ps.FEATURES].values)[:, 1]
    flagged = set(flags.name) if len(flags) else set()
    un["flagged"] = un.name.isin(flagged)
    cols = ["name", "tic", "disposition", "period", "depth_ppm", "rp_rj", "duration_ratio", "n_on_star",
            "tmag", "p_planet", "flagged"]
    return un[cols].sort_values("p_planet", ascending=False), info


def lc_evidence(r, chosen):
    flags, ev, strong = [], [], []
    if chosen.get("odd_even") and chosen["odd_even"][2](r):
        flags.append("odd_even")
        ev.append(f"odd and even transits differ by {r['oddeven_sigma']:.1f} sigma "
                  f"({r['depth_odd'] * 1e6:.0f} vs {r['depth_even'] * 1e6:.0f} ppm)")
        if r["oddeven_sigma"] > 10:
            strong.append("odd_even")
    if chosen.get("secondary") and chosen["secondary"][2](r):
        flags.append("secondary")
        ev.append(f"secondary eclipse at phase {r['sec_phase']:.3f}: {r['sec_depth'] * 1e6:.0f} ppm "
                  f"({r['sec_ratio']:.0%} of the transit depth), SNR {r['sec_snr']:.1f}, "
                  f"{r.get('sec_excess', float('nan')):.0f} sigma deeper than a planet could make it")
        if r["sec_snr"] > 15 and r.get("sec_ratio", 0) >= 0.2:
            strong.append("secondary")
    if chosen.get("centroid") and chosen["centroid"][2](r):
        flags.append("centroid")
        ev.append(f"the star's image shifts during transit ({r['centroid_sigma']:.1f} sigma; the dimming "
                  f"source would be ~{r['source_offset_px'] * p6.fp.PX_ARCSEC:.0f} arcsec away)")
        if r["centroid_sigma"] > 10:
            strong.append("centroid")
    return flags, ev, strong


def report():
    t, m, nss, eb, lc = _load()
    os.makedirs(p6.OUT, exist_ok=True)
    per = t[t.periodic]
    gf = gaia_flags(m) if len(m) else pd.DataFrame(columns=["name", "gaia_flag"])
    # Gaia check validation (all periodic candidates of each group)
    gv = []
    for g in ("planet", "fp", "unresolved"):
        n = int((per.group == g).sum())
        k = int(gf[(gf.group == g) & gf.gaia_flag].shape[0]) if len(gf) else 0
        s = int(gf[(gf.group == g) & gf.gaia_substellar].shape[0]) if len(gf) else 0
        gv.append(dict(check="gaia_orbit", variant="period match + stellar companion or EB",
                       group=g, flagged=k, total=n, rate=k / max(n, 1), substellar_matches=s))
    gv = pd.DataFrame(gv)
    chance = p6.chance_matches(t, nss, eb, n_perm=200)
    var, chosen = choose_variants(lc)
    var.to_csv(os.path.join(p6.OUT, "validation_lightcurve_checks.csv"), index=False, float_format="%.4f")
    tvar, tchosen = choose_table_variants(per)
    tvar.to_csv(os.path.join(p6.OUT, "validation_catalogue_checks.csv"), index=False, float_format="%.4f")
    gv.to_csv(os.path.join(p6.OUT, "validation_gaia.csv"), index=False, float_format="%.4f")

    # flags for candidates that are still unresolved
    lcd = {r["name"]: r for r in lc.to_dict("records")} if len(lc) else {}
    gfd = {r["name"]: r for r in gf.to_dict("records")} if len(gf) else {}
    rows = []
    for r in per[per.group == "unresolved"].itertuples():
        checks, ev, conf = [], [], None
        g = gfd.get(r.name)
        if g and g["gaia_flag"]:
            checks.append("gaia_orbit")
            ev.append(g["gaia_evidence"])
            conf = g["gaia_conf"]
        L = lcd.get(r.name)
        if L is not None and L.get("n_events"):
            f, e, strong = lc_evidence(L, chosen)
            checks += f
            ev += e
            if len(checks) >= 2:
                conf = "high"
            elif f and conf is None:
                conf = "medium" if strong else "low"
        f, e = table_evidence(r._asdict(), tchosen)
        if f:
            checks += f
            ev += e
            conf = "high" if len(checks) >= 2 else (conf or "medium")
        if checks:
            rows.append(dict(name=r.name, tic=r.tic, disposition=r.disposition or "(none)",
                             period_d=r.period, depth_ppm=r.depth_ppm, checks=";".join(checks),
                             confidence=conf, evidence=" | ".join(ev),
                             in_lc_sample=r.name in lcd,
                             sectors=",".join(map(str, lcd[r.name]["sectors"])) if r.name in lcd else ""))
    flags = pd.DataFrame(rows)
    order = {"high": 0, "medium": 1, "low": 2}
    if len(flags):
        flags = flags.sort_values(["confidence", "name"], key=lambda c: c.map(order) if c.name == "confidence" else c)
    flags.to_csv(os.path.join(p6.OUT, "likely_false_positives.csv"), index=False, float_format="%.6g")
    sub = gf[gf.gaia_substellar & (gf.group == "unresolved")] if len(gf) else gf
    sub.to_csv(os.path.join(p6.OUT, "gaia_substellar_companions.csv"), index=False)
    scores, score_info = planet_scores(t, lc, gf, flags)
    scores.to_csv(os.path.join(p6.OUT, "planet_scores.csv"), index=False, float_format="%.4g")
    _plot(var, gv)
    _summary(t, per, gv, var, chosen, flags, sub, lc, chance, tvar, scores, score_info)
    n = p6.write_ledger(current_names=set(t.name))
    print(f"{n} light-curve checks recorded in {os.path.relpath(p6.LEDGER, p6.ROOT)}")
    print(open(os.path.join(p6.OUT, "summary.md")).read()[:6000])


def _plot(var, gv):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(p6.PLOTS, exist_ok=True)
    rows = [("Gaia orbit", gv.set_index("group").rate.get("planet", 0), gv.set_index("group").rate.get("fp", 0), True)]
    for r in var.itertuples():
        rows.append((f"{r.check}: {r.variant}", r.planet_rate, r.fp_rate, r.used))
    fig, ax = plt.subplots(figsize=(10, 0.45 * len(rows) + 1.2))
    y = np.arange(len(rows))[::-1]
    ax.barh(y + 0.18, [x[1] * 100 for x in rows], 0.35, color="C2", label="confirmed / known planets")
    ax.barh(y - 0.18, [x[2] * 100 for x in rows], 0.35, color="C3", label="known false positives")
    ax.axvline(p6.MAX_PLANET_FLAG_RATE * 100, color="k", ls="--", lw=1, label="5 % limit for planets")
    ax.set_yticks(y, [("✓ " if x[3] else "   ") + x[0] for x in rows], fontsize=8)
    ax.set_xlabel("% of the group flagged")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("Phase 6: how often each check flags known planets vs known false positives (✓ = used)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(p6.PLOTS, "validation_rates.png"), dpi=110)
    plt.close(fig)


def _summary(t, per, gv, var, chosen, flags, sub, lc, chance, tvar=None, scores=None, score_info=None):
    L = ["# Phase 6: likely false positives among existing TESS candidates\n",
         "`python scripts/run_fp_triage.py` (or the app's Run page). Tables: `likely_false_positives.csv`, "
         "`validation_gaia.csv`, `validation_lightcurve_checks.csv`, `gaia_substellar_companions.csv`; "
         "plot: `plots/phase6/validation_rates.png`.\n"]
    nu = int((t.group == "unresolved").sum())
    nup = int((per.group == "unresolved").sum())
    L.append("## What was checked\n")
    L.append(f"- **Candidates:** {nu:,} TOIs and CTOIs without a final disposition. {nup:,} of them "
             "have a period, which the checks need. The validation sets come from the same "
             f"tables: {int((per.group == 'planet').sum()):,} confirmed or known planets (CP, KP) and "
             f"{int((per.group == 'fp').sum()):,} known false positives (FP, FA).")
    L.append("- **Gaia orbit check:** applied to every candidate with a period. It used "
             f"{t.gaia_id.notna().sum():,} Gaia DR3 source IDs, queried in batches.")
    n_lc = lc[lc.n_events.fillna(0) > 0]
    lc_unres = lc[lc.group == "unresolved"] if len(lc) else lc
    lc_unres_data = n_lc[n_lc.group == "unresolved"]
    tdate = p6.fp.tables_date()
    L.append(f"- **Light-curve checks:** {len(n_lc):,} candidates had usable TESS light curves in at "
             "least one sector: all the Gaia-matched candidates, random samples of about "
             f"{p6.SAMPLE} from each group, and any others checked in later runs. So far "
             f"{len(lc_unres):,} of the {nup:,} unresolved candidates with a period have been checked "
             f"({len(lc_unres_data):,} had usable data); each is recorded in `lc_checks.csv.gz` and "
             "never checked twice.")
    if tdate:
        L.append(f"- **Tables:** ExoFOP TOI and CTOI lists as downloaded on {tdate}.")
    L.append("")
    L.append("## Validation: how often each check flags known planets\n")
    L.append("Light-curve checks are chosen on the first random samples. The held-out samples, "
             "checked separately and never used to choose, show whether the choice holds up.\n")
    L.append("| check | planets flagged | known false positives flagged | held-out planets | "
             "held-out false positives | used? |")
    L.append("|---|---|---|---|---|---|")
    g = gv.set_index("group")
    L.append(f"| Gaia orbit (period match + stellar companion or eclipsing binary) | "
             f"{g.flagged['planet']}/{g.total['planet']} ({g.rate['planet']:.1%}) | "
             f"{g.flagged['fp']}/{g.total['fp']} ({g.rate['fp']:.1%}) | – | – | yes |")
    if tvar is not None:
        for r in tvar.itertuples():
            L.append(f"| {r.check} (catalogue): {r.variant} | {r.planets_flagged}/{r.planets} "
                     f"({r.planet_rate:.1%}) | {r.fps_flagged}/{r.fps} ({r.fp_rate:.1%}) | – | – | "
                     f"{'**yes**' if r.used else 'no'} |")

    def _h(k, n):
        return f"{k}/{n} ({k / n:.1%})" if n else "–"
    for r in var.itertuples():
        hp = _h(getattr(r, "holdout_planets_flagged", 0), getattr(r, "holdout_planets", 0))
        hf = _h(getattr(r, "holdout_fps_flagged", 0), getattr(r, "holdout_fps", 0))
        L.append(f"| {r.check}: {r.variant} | {r.planets_flagged}/{r.planets} ({r.planet_rate:.1%}) | "
                 f"{r.fps_flagged}/{r.fps} ({r.fp_rate:.1%}) | {hp} | {hf} | {'**yes**' if r.used else 'no'} |")
    L.append("")
    used = var[var.used & (var.get("holdout_planets", 0) > 0)] if "holdout_planets" in var else var.iloc[0:0]
    if len(used):
        hr = {r.check: r.holdout_planets_flagged / r.holdout_planets for r in used.itertuples()}
        worse = [k for k, v in hr.items() if v > p6.PREFERRED_PLANET_RATE]
        L.append("**On the held-out known planets** the light-curve checks in use flag "
                 + ", ".join(f"{v:.1%} ({k.replace('_', '/')})" for k, v in hr.items())
                 + (f". All stay within the {p6.MAX_PLANET_FLAG_RATE:.0%} limit"
                    if all(v <= p6.MAX_PLANET_FLAG_RATE for v in hr.values())
                    else f". Some exceed the {p6.MAX_PLANET_FLAG_RATE:.0%} limit")
                 + (f", but {', '.join(w.replace('_', '/') for w in worse)} flag more planets than on the "
                    "sample they were chosen on, so a candidate flagged by one of those alone could still "
                    "be a planet." if worse else ".") + "\n")
    L.append(f"The Gaia check found the Gaia orbit of the transiting object itself, with a substellar "
             f"companion mass, on {g.substellar_matches['planet']} known planets (WASP-18 b, the brown dwarf "
             f"TOI-503 b and others). Such matches are *not* counted as false positives. If Gaia periods "
             f"were unrelated to the transits, shuffling them among the {chance[2]} candidates that have a "
             f"Gaia solution gives {chance[0]:.1f} ± {chance[1]:.1f} matches by chance.\n")
    dropped = [k for k, v in chosen.items() if v is None]
    if dropped:
        L.append("**Not used, because they flag too many known planets:** " + ", ".join(dropped) + ".\n")
    L.append("## Result\n")
    if len(flags):
        c = flags.confidence.value_counts()
        by = flags.checks.str.split(";").explode().value_counts()
        L.append(f"**{len(flags)} unresolved candidates show evidence of being false positives:** "
                 f"{c.get('high', 0)} high, {c.get('medium', 0)} medium and {c.get('low', 0)} low confidence.\n")
        L.append("How they were flagged (a candidate can be flagged by more than one check):\n")
        words = dict(gaia_orbit="Gaia sees a companion star (or an eclipsing binary) on the transit's period",
                     odd_even="alternate transits have different depths: two stars eclipsing at twice the period",
                     secondary="a second eclipse half an orbit later, deeper than any planet's glow: "
                               "the companion gives off light, so it is a star",
                     centroid="the light dims off-centre: the eclipse is on a neighbouring star",
                     size="the object would be too big to be a planet",
                     density="the transit lasts longer than any orbit around this star allows")
        for k, n in by.items():
            label = dict(gaia_orbit="the Gaia orbit check", odd_even="the odd/even check").get(k, f"the {k} check")
            L.append(f"- {n} by {label}: {words.get(k, k)}.")
        L.append("")
        rnd = n_lc[(n_lc["sample"] == "random") & (n_lc.group == "unresolved")]
        rf = flags[flags.name.isin(rnd.name)]
        if len(rnd):
            share = len(rf) / len(rnd)
            L.append(f"In the random sample of {len(rnd)} unresolved candidates, {len(rf)} ({share:.0%}) were "
                     "flagged by a light-curve check. If the sample is representative, that is about "
                     f"{share * nup:.0f} of the {nup:,} unresolved candidates with periods. "
                     f"{max(nup - len(lc_unres), 0):,} have not had their light curves checked yet; "
                     "each later run checks more of them.\n")
        exp = [f"{v[3]:.1%} by the {k.replace('_', '/')} check" for k, v in chosen.items() if v]
        L.append("**How reliable the flags are.** Known planets were flagged " + ", ".join(exp)
                 + f" and {g.rate['planet']:.2%} by the Gaia check. A candidate flagged only by a light-curve "
                 "check could still be a planet, so a low-confidence flag means *look again*, not "
                 "*false positive*.\n")
        pl = n_lc[n_lc.group == "planet"]
        for k, v in chosen.items():
            if v and k == "odd_even":
                hit = sorted(r["name"] for r in pl.to_dict("records") if v[2](r))
                if hit:
                    L.append(f"Known planets flagged by the odd/even check: {', '.join(hit)}. Some are known "
                             "to have transit-timing variations or young, spotted host stars (e.g. TOI-1136, "
                             "TOI-2076, TOI-451), where a fixed ephemeris catches some transits only partly. "
                             "A candidate flagged only by this check should be checked for timing variations "
                             "first.\n")
        if score_info:
            si = score_info
            L.append("## The most likely planets among the undecided candidates\n")
            L.append(f"All checks are combined into one planet probability by a logistic regression trained "
                     f"on {si['n_train']} known cases ({si['n_planets']} confirmed planets, the rest known false "
                     "positives) using physical features only: size, orbit period, depth, transit duration "
                     "against the star's density, several candidates on one star, and the Gaia, odd/even, "
                     "centroid and secondary-eclipse evidence. Brightness and distance are left out on purpose, "
                     "because they mostly reflect which stars were followed up.\n")
            L.append(f"- **How well it separates known cases** (5-fold cross-validation): AUC {si['auc']:.2f} "
                     "(1 = perfect, 0.5 = guessing). Among the 10 % of known cases it ranked most planet-like, "
                     f"{si['top_decile_planet_frac']:.0%} were real planets"
                     + (f"; of the {si['n_p90']} it gave 90 % or more, {si['p90_planet_frac']:.0%} were."
                        if si["n_p90"] else "."))
            L.append("- **What the probability means:** the chance of being a planet for a candidate drawn "
                     "from a mix like the training set (about half planets). It ranks candidates well, but "
                     "it is not a validation: that needs follow-up observations or a full statistical "
                     "validation (e.g. TRICERATOPS with high-resolution imaging).\n")
            if scores is not None and len(scores):
                best = scores[~scores.flagged].head(15)
                L.append(f"Top {len(best)} of the {len(scores):,} undecided candidates checked so far, with no "
                         "red flag from any check (`planet_scores.csv` has all):\n")
                L.append("| candidate | TIC | disposition | period (d) | size (R_J) | candidates on the star | "
                         "planet probability |")
                L.append("|---|---|---|---|---|---|---|")
                for r in best.itertuples():
                    L.append(f"| {r.name} | {r.tic} | {r.disposition or '(none)'} | {r.period:.3f} | "
                             f"{r.rp_rj:.2f} | {r.n_on_star} | {r.p_planet:.0%} |")
                L.append("")
        L.append("## The high-confidence flags\n")
        L.append("| candidate | TIC | disposition | period (d) | checks | evidence |")
        L.append("|---|---|---|---|---|---|")
        for r in flags[flags.confidence == "high"].itertuples():
            L.append(f"| {r.name} | {r.tic} | {r.disposition} | {r.period_d:.4f} | {r.checks.replace(';', ', ')} | {r.evidence.replace('|', ';')} |")
        L.append("")
    if len(sub):
        L.append(f"**Also useful:** {len(sub)} unresolved candidates have a Gaia orbit on the transit period "
                 "with a companion below 0.08 M_sun. Gaia may be seeing the transiting object itself: a massive "
                 "planet, a brown dwarf or, near 75-80 M_Jup, one of the lowest-mass stars "
                 "(`gaia_substellar_companions.csv`).\n")
    L.append("## Caveats\n")
    L.append("- A flag is evidence, not a verdict. TFOPWG makes dispositions from all the data, including "
             "follow-up observations that these checks don't use.")
    L.append("- **Masses.** The companion mass assumes the TIC mass for the primary star. The SB1 masses "
             "assume sin i = 1, which is appropriate when the companion eclipses. Masses from astrometry "
             "assume a dark companion, so they are lower limits.")
    L.append("- **Centroids** here are the light curves' flux-weighted centroids, not difference images. "
             "Crowded fields and saturated stars can shift them.")
    L.append("- **Secondary eclipses:** hot Jupiters show real, shallow ones (WASP-18 b: 3.6 % of the "
             "transit depth). That's why the check requires a secondary that is deep relative to the "
             "transit.")
    L.append("- The light-curve checks cover only a sample of the unresolved candidates, not all of them.")
    with open(os.path.join(p6.OUT, "summary.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
