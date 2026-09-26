"""Phase 5: expert-level checks on the "submit" and "maybe" candidates.

    python scripts/phase5.py run      [--sectors 48 21] [--procs 2]   # all checks, cached
    python scripts/phase5.py report                                    # verdicts, ranking, CTOIs

Targets: every Phase 4 "submit" or "maybe" candidate of the listed sectors,
the Phase 4 control stars, and the Phase 3 validation set of each sector (known
TOIs, including the known false positives) as a calibration of each check's
false-flag rate. Checks (tesshunt/expert.py):

  1 aperture   depth in 1 / 1.5 / 2 / 3 px TESScut apertures vs a PSF model of the
               target and its TIC neighbours
  2 gp         joint celerite2 GP + transit fit of the undetrended light curve
  3 star       dwarf status and radius from Gaia DR3 (FLAME / GSP-Phot)
  4 binarity   ipd_frac_multi_peak, ipd_gof_harmonic_amplitude with RUWE, RV scatter
  5 variable   Gaia DR3 variability + AAVSO VSX
  6 physical   emcee transit fit with a stellar-density prior: implied circular period
               vs the periods TESS allows; for a known period, the eccentricity needed

Verdict per candidate:
  doubtful   any serious flag: aperture test fails; GP SNR < 7 or GP depth < 0.5x
             Phase 3; evolved host making R_p > 2 R_J; image doubling, a partly
             resolved companion or RV variability; a catalogued eclipsing binary;
             shape needing e > 0.5 (or R_p/R* > 0.3)
  plausible  no serious flag but at least one minor one: evolved host, RUWE > 1.4,
             other catalogued variability, eccentricity 0.3-0.5 needed, GP SNR 7-9
  strong     no flag at all

Ranking update (from Phase 4):
  submit -> submit unless doubtful (then maybe)
  maybe  -> drop if doubtful; submit if strong and every Phase 4 reason was only an
            FPP in 0.1-0.3; otherwise maybe
"""

import argparse
import json
import os
import sys
import time
import traceback
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tesshunt import expert as ex  # noqa: E402

WORK = os.path.join(ROOT, "work", "phase5")
OUT = os.path.join(ROOT, "results", "phase5")
PLOTS = os.path.join(ROOT, "plots", "phase5")
RHO_ERR_FRAC = 0.25            # density uncertainty when Gaia FLAME has no mass/radius range
GP_SNR_GOOD = 9.0
E_MINOR, E_SERIOUS = 0.3, 0.5
RP_RATIO_MAX = 0.3


def jdump(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh, default=_js, indent=1)
    os.replace(path + ".tmp", path)


def _js(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def jload(path):
    with open(path) as fh:
        return json.load(fh)


# ------------------------------------------------------------------ targets

def _toi_periods():
    """The full TOI table, downloaded once in bulk (cached)."""
    import io
    import urllib.parse
    from tesshunt import net
    q = "select tid,toi,tfopwg_disp,pl_orbper from toi"
    url = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=" + urllib.parse.quote(q)
           + "&format=csv")
    return pd.read_csv(io.StringIO(net.get_text(url, "exoplanet_archive")))


def targets(sectors):
    import phase4
    from phase3_vet import sector_dirs
    tois = _toi_periods()
    rows = []
    for s in sectors:
        _, res4, _ = sector_dirs(s, 4)
        fpath = os.path.join(res4, "followup.csv")
        if not os.path.exists(fpath):
            print(f"S{s}: no Phase 4 results yet, skipped")
            continue
        fu = pd.read_csv(fpath)
        cands = phase4.candidates(s)
        cmap = cands.set_index("tic").to_dict("index")
        for _, r in fu.iterrows():
            if r.role == "candidate" and r.category not in ("submit", "maybe"):
                continue
            c = cmap.get(r.tic)
            if c is None:
                continue
            per = None
            if isinstance(r.joint_periods, str) and r.joint_periods:
                per = [float(x) for x in r.joint_periods.split(";")][:2]
            if r.role == "control":
                tp = tois[tois.tid == r.tic].pl_orbper.dropna()
                per = [float(tp.iloc[0])] if len(tp) else per
            rows.append(dict(key=f"s{s:04d}_{r.tic}", tic=int(r.tic), sector=s, role=r.role,
                             category=r.category, reasons=r.reasons if isinstance(r.reasons, str) else "",
                             t0=c["t0"], t14=c["t14"], depth=c["depth"], snr=c["snr"], fpp=c["fpp"],
                             ra=c["ra"], dec=c["dec"], Tmag=c["Tmag"], rad=c["rad"], mass=c["mass"],
                             periods=per, binarity_flags=r.binarity_flags if isinstance(r.binarity_flags, str) else ""))
        val = phase4.validation_planets(s)
        for c in val.to_dict("records"):
            disp = str(c.get("tfop_disp") or "")
            tp = tois[tois.tid == c["tic"]].pl_orbper.dropna()
            rows.append(dict(key=f"s{s:04d}_val_{c['key']}", tic=int(c["tic"]), sector=s,
                             role="false_positive" if ("FP" in disp or "FA" in disp) else "validation",
                             category="", reasons="", toi=c.get("toi"), tfop_disp=disp,
                             t0=c["t0"], t14=c["t14"], depth=c["depth"], snr=c["snr"], fpp=c["fpp"],
                             ra=c["ra"], dec=c["dec"], Tmag=c["Tmag"], rad=c["rad"], mass=c["mass"],
                             periods=[float(tp.iloc[0])] if len(tp) else None, binarity_flags=""))
    df = pd.DataFrame(rows)
    # Gaia DR3 ids from the TIC (cached, one bulk query)
    from tesshunt import tic as tic_mod
    tt = tic_mod.query_ids(sorted(set(df.tic)))
    gmap = {int(i): (int(g) if str(g).isdigit() else None)
            for i, g in zip(tt["ID"], [str(x) for x in tt["GAIA"]])}
    df["gaia_id"] = df.tic.map(gmap)
    return df


# ------------------------------------------------------------------ checks

def _lc(tic, sector):
    from tesshunt import ffi
    lc, _ = ffi.read(ffi.download(int(tic), int(sector), cache=True))
    return lc


def _cached(key, name, fn):
    p = os.path.join(WORK, key, f"{name}.json")
    if os.path.exists(p):
        return jload(p)
    try:
        out = fn()
    except Exception as e:  # noqa: BLE001
        out = dict(error=f"{type(e).__name__}: {e}", tb=traceback.format_exc())
    jdump(p, out)
    return out


def rho_prior(star, gs):
    """Stellar density and uncertainty (cgs) from Gaia FLAME mass/radius, else TIC."""
    m, r = star.get("mass"), star.get("radius")
    if gs and ex._num(gs.get("mass_flame")) and ex._num(gs.get("radius_flame")):
        m, r = ex._num(gs["mass_flame"]), ex._num(gs["radius_flame"])
        dm = (ex._num(gs.get("mass_flame_upper")) or m * 1.1) - (ex._num(gs.get("mass_flame_lower")) or m * 0.9)
        dr = (ex._num(gs.get("radius_flame_upper")) or r * 1.05) - (ex._num(gs.get("radius_flame_lower")) or r * 0.95)
        frac = float(np.hypot(dm / 2 / m, 3 * dr / 2 / r))
        src = "Gaia FLAME"
    else:
        frac, src = RHO_ERR_FRAC, "TIC (±25 %)"
    if not (m and r):
        return None, None, "none"
    rho = m / r ** 3 * ex.RHO_SUN
    return float(rho), float(max(frac, 0.1) * rho), src


def e_min_from_ratio(ratio):
    """Smallest eccentricity that makes a transit shorter (ratio > 1) or longer
    (ratio < 1) than circular by the density ratio rho_needed / rho_star."""
    g2 = np.asarray(ratio, float) ** (2 / 3)
    return np.abs(g2 - 1) / (g2 + 1)


def run_target(t):
    key = t["key"]
    from phase3_vet import sector_dirs
    s = int(t["sector"])
    res = {}
    res["aperture"] = _cached(key, "aperture", lambda: ex.aperture_test(
        t["ra"], t["dec"], t["Tmag"], s, t["t0"], t["t14"]))
    lc = _lc(t["tic"], s)
    cad = float(np.median(np.diff(lc.time)))

    def gp():
        g = ex.gp_reprocess(lc.time, lc.flux, t["t0"], t["t14"], t["depth"], cad)
        g.pop("arrays", None)
        return g
    res["gp"] = _cached(key, "gp", gp)
    res["gaia"] = _cached(key, "gaia", lambda: ex.gaia_star(t["gaia_id"]) if t.get("gaia_id") else None)
    gs = res["gaia"] if isinstance(res["gaia"], dict) and "error" not in res["gaia"] else None
    dil = 2.0 if "wds_close" in (t.get("binarity_flags") or "") else 1.0
    res["star"] = ex.stellar_check(gs, t["rad"], t["mass"], t["depth"], dil)
    res["binarity"] = ex.binarity_extras(gs)
    vs = _cached(key, "vsx", lambda: ex.vsx(t["ra"], t["dec"]).to_dict("records"))
    res["variability"] = ex.variability_check(gs, pd.DataFrame(vs) if isinstance(vs, list) else None)

    def physical():
        from tesshunt import multisector as ms
        tt, ff = ms.candidate_detrend(lc, t["t14"])
        w = np.abs(tt - t["t0"]) < max(2.5 * t["t14"], 0.5)
        rho, drho, src = rho_prior(dict(mass=t["mass"], radius=res["star"]["radius"] or t["rad"]), gs)
        if rho is None:
            return dict(error="no stellar density")
        out = dict(rho_star=rho, rho_err=drho, rho_source=src)
        rng = np.random.default_rng(3)
        if t.get("periods"):
            out["known"] = []
            for P in t["periods"]:
                o = ex.density_fit(tt[w], ff[w], t["t0"], t["t14"], t["depth"], cad, rho, drho,
                                   known_period=P, nsteps=2000, burn=800)
                rs = np.array(o.pop("rho_samples"))
                ratio = rs / np.clip(rng.normal(rho, drho, len(rs)), 0.05 * rho, None)
                em = e_min_from_ratio(ratio)
                o.pop("model", None)
                out["known"].append(dict(P=P, rho_needed=o["rho"], ratio=[float(x) for x in np.percentile(ratio, [16, 50, 84])],
                                         e_min=[float(x) for x in np.percentile(em, [16, 50, 84])],
                                         rp=o["rp"], b=o["b"], grazing_prob=o["grazing_prob"]))
        o = ex.density_fit(tt[w], ff[w], t["t0"], t["t14"], t["depth"], cad, rho, drho,
                           nsteps=2500, burn=1000)
        out["model"] = o.pop("model")
        out["t_fit"] = tt[w].tolist()
        out["f_fit"] = ff[w].tolist()
        per_path = os.path.join(sector_dirs(s, 4)[0], "periods", f"{t['tic']}.json")
        if os.path.exists(per_path) and t["role"] in ("candidate", "control"):
            per = jload(per_path)
            out["periods"] = ex.period_consistency(o["P_samples"], per.get("intervals"),
                                                   per.get("longest_excluded_below"), per.get("pmax"))
            out["periods"]["intervals"] = per.get("intervals")
            out["periods"]["floor"] = per.get("longest_excluded_below")
        out["P_hist"] = np.histogram(np.log10(o.pop("P_samples")), bins=60, range=(0, 4))[0].tolist()
        out["single"] = o
        return out
    res["physical"] = _cached(key, "physical", physical)
    return res


def _run(t):
    t0 = time.time()
    try:
        r = run_target(t)
        errs = [k for k, v in r.items() if isinstance(v, dict) and "error" in v]
        return t["key"], ("ok" if not errs else f"errors in {errs}"), time.time() - t0
    except Exception as e:  # noqa: BLE001
        return t["key"], f"{type(e).__name__}: {e}", time.time() - t0


def cmd_run(args):
    tg = targets(args.sectors)
    os.makedirs(WORK, exist_ok=True)
    tg.to_json(os.path.join(WORK, "targets.json"), orient="records", indent=1)
    print(tg.groupby(["sector", "role"]).size().to_string())
    recs = tg.to_dict("records")
    if args.only:
        recs = [x for x in recs if x["key"] in args.only]
    with Pool(args.procs) as pool:
        for i, (k, st, dt) in enumerate(pool.imap_unordered(_run, recs), 1):
            print(f"[{time.strftime('%H:%M:%S')}] {i}/{len(recs)} {k}: {st} ({dt:.0f}s)", flush=True)


# ------------------------------------------------------------------ assessment

def assess(t, r):
    """Plain-English sentence and flags per check; overall verdict."""
    lines, serious, minor = {}, [], []
    ap = r["aperture"]
    if "error" in ap:
        lines["aperture"] = f"Aperture test could not run ({ap['error']})."
    else:
        dep = ", ".join(f"{x['depth'] * 1e6:.0f}" for x in ap["rows"])
        if ap["flag"]:
            serious.append("aperture")
            lines["aperture"] = (f"**The depth changes with aperture size** ({dep} ppm in 1, 1.5, 2 "
                                 f"and 3 px apertures; worst {ap['worst_resid_sigma']:.1f}σ from what a "
                                 "signal on this star would give)"
                                 + (f", and fits neighbour TIC {ap['best_other']['source']} "
                                    f"{ap['best_other']['sep_px']:.1f} px away better"
                                    if ap.get("prefers_other") else "")
                                 + ". The dip may not come from this star, or may be instrumental.")
        else:
            un = (f" (TIC {', '.join(ap['unresolved'])} lies within 1 px and cannot be separated "
                  "this way)") if ap.get("unresolved") else ""
            lines["aperture"] = (f"The depth is the same, within the noise, in small and large "
                                 f"apertures ({dep} ppm), as expected if the dip comes from this star{un}.")
            if ap.get("worst_resid_sigma") and ap["worst_resid_sigma"] > 3:
                minor.append("aperture_marginal")
    gp = r["gp"]
    if "error" in gp:
        lines["gp"] = f"GP re-analysis failed ({gp['error']})."
    else:
        ratio = gp["depth"] / t["depth"] if t["depth"] else np.nan
        txt = (f"Re-analysed with a different noise model (a Gaussian process fitted together "
               f"with the transit): depth {gp['depth'] * 1e6:.0f} ppm ({ratio:.2f}× Phase 3), "
               f"significance {gp['snr']:.1f}σ (Phase 3: {t['snr']:.1f}).")
        if gp["snr"] < ex.GP_SNR_MIN or ratio < ex.GP_DEPTH_MIN:
            serious.append("gp")
            txt += " **The dip is not robust to the choice of noise model.**"
        elif gp["snr"] < GP_SNR_GOOD:
            minor.append("gp_weak")
            txt += " It survives, but only moderately."
        lines["gp"] = txt
    st = r["star"]
    if st["cls"] == "dwarf" and st["radius_source"] == "TIC":
        lines["star"] = (f"Gaia DR3 has no evolutionary parameters for this star; its TIC radius "
                         f"({st['radius']:.2f} R☉) is consistent with a dwarf, giving a companion of "
                         f"about {st['rp_rj']:.2f} R_J ({st['rp_rj'] * 11.2:.0f} R⊕). Not independently confirmed.")
    elif st["cls"] == "dwarf":
        lines["star"] = (f"Gaia confirms a main-sequence (dwarf) star ({st['basis']}); radius "
                         f"{st['radius']:.2f} R☉ ({st['radius_source']}), so the companion is about "
                         f"{st['rp_rj']:.2f} R_J ({st['rp_rj'] * 11.2:.0f} R⊕).")
    elif st["cls"] == "unknown":
        lines["star"] = (f"Gaia DR3 has no evolutionary parameters for this star; with the TIC "
                         f"radius ({t['rad']:.2f} R☉) the companion is about {st['rp_rj_tic']:.2f} R_J.")
    else:
        lines["star"] = (f"The star is **{st['cls']}** ({st['basis']}; radius {st['radius']:.2f} R☉ "
                         f"from {st['radius_source']}). The companion is then about {st['rp_rj']:.2f} R_J"
                         f" (TIC radius gave {st['rp_rj_tic']:.2f}).")
        (serious if st["serious"] else minor).append("evolved")
    bx = r["binarity"]
    if bx["serious"]:
        serious.append("binarity")
        lines["binarity"] = "**Gaia sees signs of a close companion:** " + "; ".join(bx["notes"]) + "."
    elif "ruwe" in bx["flags"]:
        minor.append("ruwe")
        lines["binarity"] = "Gaia's astrometry is noisier than for a single star (" + "; ".join(bx["notes"]) + \
            "), but the image shape and velocities show no companion."
    else:
        lines["binarity"] = "No sign of a close companion in Gaia's image shape or velocities (" + \
            "; ".join(bx["notes"]) + ")."
    vv = r["variability"]
    if vv["serious"]:
        serious.append("eclipsing")
        lines["variable"] = "**Catalogued as an eclipsing binary:** " + "; ".join(vv["notes"]) + "."
    elif vv["flags"]:
        minor.append("variable")
        lines["variable"] = "Catalogued as variable: " + "; ".join(vv["notes"]) + "."
    else:
        lines["variable"] = "Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue."
    ph = r["physical"]
    if "error" in ph:
        lines["physical"] = f"Physical-consistency fit failed ({ph['error']})."
    else:
        sg = ph["single"]
        rp = sg["rp"][1]
        txt = ""
        if ph.get("known"):
            parts = []
            worst_e = []
            for k in ph["known"]:
                e16, e50, _ = k["e_min"]
                parts.append(f"for P = {k['P']:.2f} d the shape needs {k['ratio'][1]:.1f}× the star's "
                             f"density on a circular orbit, i.e. an eccentricity of at least "
                             f"~{e50:.2f} (≥ {e16:.2f} at 1σ)")
                worst_e.append(e16)
            txt = "Fitted with this star's density (" + ph["rho_source"] + "): " + "; ".join(parts) + "."
            best = min(worst_e)
            if best > E_SERIOUS:
                serious.append("physical")
                txt += " **No listed period fits a moderately eccentric orbit.**"
            elif best > E_MINOR:
                minor.append("eccentric")
                txt += " This is possible but needs a fairly eccentric orbit."
            elif len(ph["known"]) > 1:
                es = [k["e_min"][1] for k in ph["known"]]
                pref = ph["known"][int(np.argmin(es))]["P"]
                txt += f" The shape mildly favours P = {pref:.2f} d."
        if "periods" in ph and not ph.get("known"):
            pc = ph["periods"]
            txt = (f"Fitted with this star's density ({ph['rho_source']}), the dip's shape implies a "
                   f"circular-orbit period of ~{pc['P_circ'][1]:.0f} d ({pc['P_circ'][0]:.0f}–"
                   f"{pc['P_circ'][2]:.0f}); {pc['frac_allowed']:.0%} of that range is still "
                   "allowed by TESS's other observations.")
            if pc["frac_allowed"] < 0.05:
                if pc["reachable_with_ecc"]:
                    minor.append("eccentric")
                    txt += " An allowed period is reachable with a moderately eccentric orbit (e ≤ 0.5)."
                else:
                    serious.append("physical")
                    txt += " **No allowed period is reachable even with e ≤ 0.5.**"
        elif not ph.get("known"):
            txt = (f"Fitted with this star's density ({ph['rho_source']}), the shape implies a "
                   f"circular-orbit period of ~{sg['P'][1]:.0f} d ({sg['P'][0]:.0f}–{sg['P'][2]:.0f}).")
        if rp > RP_RATIO_MAX:
            serious.append("too_large")
            txt += f" **R_p/R★ = {rp:.2f} is too large for a planet.**"
        if sg["grazing_prob"] > 0.5:
            txt += f" The fit is probably grazing (p = {sg['grazing_prob']:.2f}), so the size is uncertain."
        lines["physical"] = txt
    verdict = "doubtful" if serious else ("plausible" if minor else "strong")
    return dict(lines=lines, serious=serious, minor=minor, verdict=verdict)


def rerank(t, a):
    cat = t["category"]
    if t["role"] != "candidate":
        return None, ""
    v = a["verdict"]
    reasons = [x.strip() for x in (t.get("reasons") or "").split(";") if x.strip()]
    only_fpp = bool(reasons) and all(x.startswith("S1 FPP") for x in reasons) and (t["fpp"] or 1) < 0.3
    if cat == "submit":
        return ("maybe", "Phase 5 doubtful: " + ", ".join(a["serious"])) if v == "doubtful" else ("submit", "")
    if v == "doubtful":
        return "drop", "Phase 5 doubtful: " + ", ".join(a["serious"])
    if v == "strong" and only_fpp:
        return "submit", f"Phase 5 strong; Phase 4 held it back only for FPP {t['fpp']:.2f} < 0.3"
    return "maybe", ""


def cmd_report(args):
    tg = pd.read_json(os.path.join(WORK, "targets.json"), orient="records")
    rows, md, cal = [], [], []
    for t in tg.to_dict("records"):
        try:
            r = run_target(t)            # all cached after `run`
        except Exception as e:  # noqa: BLE001
            print("skip", t["key"], e)
            continue
        a = assess(t, r)
        new_cat, why = rerank(t, a)
        ph = r["physical"] if "error" not in r["physical"] else {}
        rows.append(dict(key=t["key"], tic=t["tic"], sector=t["sector"], role=t["role"],
                         phase4=t["category"], phase5=new_cat, rank_change=why, verdict=a["verdict"],
                         serious=";".join(a["serious"]), minor=";".join(a["minor"]),
                         ap_worst_sigma=r["aperture"].get("worst_resid_sigma"),
                         gp_depth_ppm=r["gp"].get("depth", np.nan) * 1e6 if "depth" in r["gp"] else np.nan,
                         gp_snr=r["gp"].get("snr"), phase3_depth_ppm=t["depth"] * 1e6, phase3_snr=t["snr"],
                         star_class=r["star"]["cls"], radius_rsun=r["star"]["radius"],
                         radius_source=r["star"]["radius_source"], rp_rj=r["star"]["rp_rj"] or r["star"]["rp_rj_tic"],
                         ruwe=r["binarity"].get("ruwe"), binarity_flags=";".join(r["binarity"]["flags"]),
                         variability_flags=";".join(r["variability"]["flags"]),
                         P_circ_med=(ph.get("single") or {}).get("P", [None, None])[1],
                         frac_P_allowed=(ph.get("periods") or {}).get("frac_allowed"),
                         e_min_known=";".join(f"{k['P']:.2f}:{k['e_min'][1]:.2f}" for k in ph.get("known", []))))
        if t["role"] in ("candidate", "control"):
            md.append(_card(t, a, new_cat, why))
        _plot(t, r, a)
    df = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(os.path.join(OUT, "phase5_checks.csv"), index=False, float_format="%.5g")
    # calibration on the validation planets
    val = df[df.role == "validation"]
    fps = df[df.role == "false_positive"]
    checks = ["aperture", "gp", "evolved", "binarity", "eclipsing", "physical", "too_large"]
    for c in checks:
        cal.append(dict(check=c, planets=len(val), flagged_serious=int(val.serious.str.contains(c).sum()),
                        false_positives=len(fps), fp_flagged=int(fps.serious.str.contains(c).sum())))
    cal = pd.DataFrame(cal)
    cal.to_csv(os.path.join(OUT, "calibration.csv"), index=False)
    ver = val.verdict.value_counts().to_dict()
    with open(os.path.join(OUT, "candidates.md"), "w") as fh:
        fh.write("# Phase 5: expert checks, candidate by candidate\n\n"
                 "Produced by `python scripts/phase5.py report`. Numbers: `phase5_checks.csv`; "
                 "calibration on known planets: `calibration.csv`; plots: `plots/phase5/`.\n\n"
                 f"**Calibration.** The same checks on {len(val)} known planets from the "
                 f"validation sets: {ver.get('strong', 0)} strong, {ver.get('plausible', 0)} plausible, "
                 f"{ver.get('doubtful', 0)} doubtful. Serious flags per check: "
                 + ", ".join(f"{r.check} {r.flagged_serious}" for r in cal.itertuples())
                 + f". On the {len(fps)} known false positives: "
                 + ", ".join(f"{v} {n}" for v, n in fps.verdict.value_counts().items()) + ".\n\n"
                 + "\n".join(md))
    ctoi(df, tg)
    print(df[df.role.isin(["candidate", "control"])][["key", "phase4", "phase5", "verdict", "serious", "minor"]].to_string(index=False))
    print(cal.to_string(index=False))
    print("validation verdicts:", ver, " false positives:", fps.verdict.value_counts().to_dict())


def _card(t, a, new_cat, why):
    names = dict(aperture="Aperture test", gp="Independent reprocessing", star="Stellar check",
                 binarity="Gaia binarity", variable="Variability", physical="Physical consistency")
    head = (f"## TIC {t['tic']} (Sector {t['sector']}"
            + (", control" if t["role"] == "control" else "") + f") — **{a['verdict']}**")
    L = [head, ""]
    if t["role"] == "candidate":
        L.append(f"Phase 4: **{t['category']}** → Phase 5: **{new_cat}**" + (f" ({why})" if why else "") + "\n")
    for k in ("aperture", "gp", "star", "binarity", "variable", "physical"):
        L.append(f"- **{names[k]}.** {a['lines'].get(k, '–')}")
    L.append("")
    why_v = {"strong": "every check passed.",
             "plausible": "no check failed, but some raise minor concerns: " + ", ".join(a["minor"]) + ".",
             "doubtful": "failed: " + ", ".join(a["serious"]) + "."}[a["verdict"]]
    L.append(f"**Verdict: {a['verdict']}** — {why_v}")
    L.append(f"Plot: `plots/phase5/{t['key']}.png`\n")
    return "\n".join(L)


def _plot(t, r, a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(PLOTS, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.8))
    ax = axes[0]
    ap = r["aperture"]
    if "rows" in ap:
        x = [q["radius"] for q in ap["rows"]]
        ax.errorbar(x, [q["depth"] * 1e6 for q in ap["rows"]], [q["err"] * 1e6 for q in ap["rows"]],
                    fmt="o", color="k", label="measured")
        ax.plot(x, np.array(ap["hypotheses"][0]["expected"]) * 1e6, "C0-", label="expected: on target")
        if ap.get("best_other"):
            ax.plot(x, np.array(ap["best_other"]["expected"]) * 1e6, "C3--",
                    label=f"on TIC {ap['best_other']['source']}")
        ax.set_xlabel("aperture radius [px]", fontsize=8)
        ax.set_ylabel("depth [ppm]", fontsize=8)
        ax.legend(fontsize=7)
    ax.set_title("1 aperture test", fontsize=9)
    ax = axes[1]
    ph = r["physical"]
    if "t_fit" in ph:
        tt = np.array(ph["t_fit"])
        ax.plot((tt - t["t0"]) * 24, ph["f_fit"], ".", color="0.6", ms=2)
        ax.plot((tt - t["t0"]) * 24, ph["model"], "C3", lw=1.5, label="density-prior fit")
        ax.set_xlabel("hours from T0", fontsize=8)
        gp = r["gp"]
        if "depth" in gp:
            ax.axhline(1 - gp["depth"], color="C0", ls=":", label=f"GP depth ({gp['snr']:.1f}σ)")
        ax.legend(fontsize=7)
    ax.set_title("2/6 transit fits", fontsize=9)
    ax = axes[2]
    if "P_hist" in ph:
        edges = np.logspace(0, 4, 61)
        ax.stairs(ph["P_hist"], edges, fill=True, alpha=0.5, label="circular-orbit P from shape")
        if ph.get("periods"):
            for lo, hi in ph["periods"].get("intervals") or []:
                ax.axvspan(lo, hi, color="C2", alpha=0.25, lw=0)
            if ph["periods"].get("floor"):
                ax.axvspan(ph["periods"]["floor"], 1e4, color="C2", alpha=0.25, lw=0,
                           label="allowed by TESS")
        for k in ph.get("known", []):
            ax.axvline(k["P"], color="C3", label=f"P = {k['P']:.1f} d (e ≥ {k['e_min'][1]:.2f})")
        ax.set_xscale("log")
        ax.set_xlabel("period [d]", fontsize=8)
        ax.legend(fontsize=7)
    ax.set_title("6 period vs shape", fontsize=9)
    fig.suptitle(f"TIC {t['tic']} (S{t['sector']}, {t['role']}): {a['verdict']}"
                 + (f" — {', '.join(a['serious'])}" if a["serious"] else ""), fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, f"{t['key']}.png"), dpi=90)
    plt.close(fig)


def ctoi(df, tg):
    """CTOI summaries for the Phase 5 "submit" list, built with the Phase 4 code
    and extended with the Phase 5 verdict."""
    import phase4
    from phase4_report import ctoi_rows, gather
    sub = df[(df.role == "candidate") & (df.phase5 == "submit")]
    rows_all, md_all = [], []
    for s in sorted(sub.sector.unique()):
        cands = phase4.candidates(s)
        P = phase4.paths(s)
        gs = []
        for tic in sub[sub.sector == s].tic:
            c = cands[cands.tic == tic].iloc[0].to_dict()
            g = gather(c, P, s)
            g["category"], g["role"] = "submit", "candidate"
            g["positives"] = g.get("positives") or []
            gs.append(g)
        for g in gs:                        # one at a time so each gets its own Phase 5 line
            rows, md = ctoi_rows([g], s)
            v = df[(df.tic == g["tic"]) & (df.sector == s)].iloc[0]
            summary = v.verdict + (f"; minor concerns: {v.minor}" if v.minor else "; all six checks passed")
            rows["Notes"] = rows["Notes"] + f" Phase 5 expert checks: {summary}."
            rows.insert(1, "Sector searched", s)
            md = md.replace("- **Other TESS data:**", f"- **Phase 5 expert checks:** {summary} "
                            f"(details: `results/phase5/candidates.md`)\n- **Other TESS data:**", 1)
            rows_all.append(rows)
            md_all.append(md)
    if rows_all:
        pd.concat(rows_all).to_csv(os.path.join(OUT, "ctoi_candidates.csv"), index=False)
    with open(os.path.join(OUT, "ctoi_summaries.md"), "w") as fh:
        fh.write("# CTOI-ready summaries after Phase 5\n\nPrepared for submission to ExoFOP as "
                 "community TOIs; **not submitted**.\n\n" + "\n".join(md_all))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--sectors", type=int, nargs="+", default=[48, 21])
    r.add_argument("--procs", type=int, default=2)
    r.add_argument("--only", nargs="*", help="restrict to these target keys (tests)")
    sub.add_parser("report")
    args = ap.parse_args()
    {"run": cmd_run, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    main()
