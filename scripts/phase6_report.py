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

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import phase6 as p6  # noqa: E402

M_STAR = 0.08
M_HIGH = 0.10


def _load():
    W = p6.WORK
    t = pd.read_csv(os.path.join(W, "candidates.csv"), dtype={"gaia_id": "Int64"})
    m = pd.read_csv(os.path.join(W, "gaia_matches.csv"))
    nss = pd.read_csv(os.path.join(W, "nss.csv"))
    eb = pd.read_csv(os.path.join(W, "gaia_eb.csv"))
    lc = []
    for p in glob.glob(os.path.join(W, "lc", "*.json")):
        with open(p) as fh:
            lc.append(json.load(fh))
    lc = pd.DataFrame(lc)
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
                          "(substellar: Gaia sees the transiting object's own orbit)")
            else:
                ev.append(f"{r.source} orbit P = {r.gaia_period:.3f} d{k}; no mass estimate")
                if not flag:
                    flag, conf = True, "medium"
        rows.append(dict(name=name, tic=g.tic.iloc[0], group=g.group.iloc[0], gaia_flag=flag,
                         gaia_conf=conf if flag else None, gaia_substellar=sub and not flag,
                         gaia_evidence="; ".join(ev)))
    return pd.DataFrame(rows)


def choose_variants(lc):
    """Flag rates of every variant on the random samples, and the variant used."""
    rnd = lc[(lc["sample"] == "random") & (lc.n_events.fillna(0) > 0)]
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
            rows.append(dict(check=check, variant=label,
                             planets_flagged=rate["planet"][0], planets=rate["planet"][1], planet_rate=pr,
                             fps_flagged=rate["fp"][0], fps=rate["fp"][1], fp_rate=fr,
                             unresolved_flagged=rate["unresolved"][0], unresolved=rate["unresolved"][1]))
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
                  f"({r['sec_ratio']:.0%} of the transit depth), SNR {r['sec_snr']:.1f}")
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
    _plot(var, gv)
    _summary(t, per, gv, var, chosen, flags, sub, lc, chance)
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


def _summary(t, per, gv, var, chosen, flags, sub, lc, chance):
    L = ["# Phase 6: likely false positives among existing TESS candidates\n",
         "`python scripts/phase6.py gaia | lc | report`. Tables: `likely_false_positives.csv`, "
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
    L.append(f"- **Light-curve checks:** {len(n_lc):,} candidates had usable TESS light curves in at "
             "least one sector. These were all the Gaia-matched candidates plus random samples of "
             f"about {p6.SAMPLE} from each group.")
    L.append("")
    L.append("## Validation: how often each check flags known planets\n")
    L.append("| check | planets flagged | known false positives flagged | used? |")
    L.append("|---|---|---|---|")
    g = gv.set_index("group")
    L.append(f"| Gaia orbit (period match + stellar companion or eclipsing binary) | "
             f"{g.flagged['planet']}/{g.total['planet']} ({g.rate['planet']:.1%}) | "
             f"{g.flagged['fp']}/{g.total['fp']} ({g.rate['fp']:.1%}) | yes |")
    for r in var.itertuples():
        L.append(f"| {r.check}: {r.variant} | {r.planets_flagged}/{r.planets} ({r.planet_rate:.1%}) | "
                 f"{r.fps_flagged}/{r.fps} ({r.fp_rate:.1%}) | {'**yes**' if r.used else 'no'} |")
    L.append("")
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
                     secondary="a second, shallower eclipse: the 'planet' also gives off light, so it is a star",
                     centroid="the light dims off-centre: the eclipse is on a neighbouring star")
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
                     f"{share * nup:.0f} of the {nup:,} unresolved candidates with periods; only "
                     f"{len(n_lc[(n_lc['sample'] == 'random') & (n_lc.group == 'unresolved')]):,} were "
                     "examined here.\n")
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
        L.append("## The high-confidence flags\n")
        L.append("| candidate | TIC | disposition | period (d) | checks | evidence |")
        L.append("|---|---|---|---|---|---|")
        for r in flags[flags.confidence == "high"].itertuples():
            L.append(f"| {r.name} | {r.tic} | {r.disposition} | {r.period_d:.4f} | {r.checks.replace(';', ', ')} | {r.evidence.replace('|', ';')} |")
        L.append("")
    if len(sub):
        L.append(f"**Also useful:** {len(sub)} unresolved candidates have a Gaia orbit on the transit period "
                 "with a *substellar* companion. Gaia may be seeing the transiting object itself (a brown "
                 "dwarf or massive planet): `gaia_substellar_companions.csv`.\n")
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
