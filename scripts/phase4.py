"""Phase 4: strengthen the shortlist with all of TESS's data on each star.

    python scripts/phase4.py coverage   # sectors + best light curve per sector (S3 first)
    python scripts/phase4.py search     # Phase 2 search + duotransit search + checks 1-4
    python scripts/phase4.py periods    # allowed-period maps
    python scripts/phase4.py confirm    # the S48 dip in QLP / SPOC 2-min / TESScut + step test
    python scripts/phase4.py binarity   # Gaia DR3 RUWE, NSS, companions, WDS, El-Badry
    python scripts/phase4.py validate   # rules D1 and D5 on the Phase 3 validation planets
    python scripts/phase4.py rank       # submit / maybe / drop + CTOI summaries
    python scripts/phase4.py report     # sheets, CSVs, summary
    python scripts/phase4.py all

Candidates are the new ones in results/phase{3}/shortlist.csv, plus control
stars (default: TOI-2180 b, whose other transits and 260 d period are known)
that go through exactly the same steps. Per-candidate results are cached in
work/phase4/ and each stage skips what is already done.
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

from tesshunt import binarity as bn  # noqa: E402
from tesshunt import multisector as ms  # noqa: E402
from tesshunt import vetting as v  # noqa: E402
from tesshunt.survey import analyse  # noqa: E402

CONTROLS = {298663873: "TOI-2180 b (P = 260.17 d, transits in S19, S48, S57)"}


def paths(sector):
    from phase3_vet import sector_dirs
    tag = f"s{sector:04d}"
    p3work, p3, _ = sector_dirs(sector, 3)
    _, out, plots = sector_dirs(sector, 4)
    return dict(work=os.path.join(ROOT, "work", "phase4", tag), p3work=p3work, p3=p3,
                out=out, plots=plots, sample=os.path.join(ROOT, "work", f"{tag}_sample.csv"))


def jdump(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=_js, indent=1)
    os.replace(tmp, path)


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


# ---------------------------------------------------------------- candidates

def candidates(sector, controls=True):
    """New shortlist candidates (+ controls) with precise Phase 3 fit parameters."""
    P = paths(sector)
    sl = pd.read_csv(os.path.join(P["p3"], "shortlist.csv"))
    va = pd.read_csv(os.path.join(P["p3"], "vetting_all.csv"))
    rows = []
    for _, r in sl[sl.new].iterrows():
        rows.append(dict(tic=int(r.tic), role="candidate", rank3=int(r.shortlist_rank)))
    if controls:
        for tic, label in CONTROLS.items():
            rows.append(dict(tic=tic, role="control", rank3=None, label=label))
    return _with_fits(P, sl, va, rows)


def validation_planets(sector):
    """The Phase 3 validation set (known TOIs with one predicted transit in the
    sector, plus TOI-2180 b), with the same fit parameters as candidates."""
    P = paths(sector)
    sl = pd.read_csv(os.path.join(P["p3"], "shortlist.csv"))
    va = pd.read_csv(os.path.join(P["p3"], "vetting_all.csv"))
    val = va[va.is_validation & (va.lc_status == "ok")]
    rows = [dict(tic=int(t), role="validation", rank3=None) for t in sorted(set(val.tic))]
    df = _with_fits(P, sl, va, rows, validation=True)
    info = val.sort_values("snr", ascending=False).drop_duplicates("tic").set_index("tic")
    df["toi"] = [info.tois.get(t) for t in df.tic]
    df["tfop_disp"] = [info.toi_disp.get(t) for t in df.tic]
    return df


def _with_fits(P, sl, va, rows, validation=False):
    out = []
    for r in rows:
        g = va[(va.tic == r["tic"]) & (va.lc_status == "ok")]
        if validation:
            g = g[g.is_validation]
        elif r["role"] == "candidate":
            g = g[g.is_candidate]
        g = g.sort_values("snr", ascending=False)
        if g.empty:
            continue
        k = g.iloc[0]
        lc = jload(os.path.join(P["p3work"], "lc", f"{k.key}.json"))
        tr = lc["shape"]["transit"]
        r.update(key=k.key, t0=tr["t0"], t14=tr["t14_h"] / 24, depth=tr["depth_ppm"] * 1e-6,
                 b=tr["b"], rp=tr["rp"], snr=float(k.snr), tier=k.tier if isinstance(k.tier, str) else "",
                 fpp=float(k.fpp_value) if k.fpp_value == k.fpp_value else None,
                 p_min_circ=lc["duration"].get("p_b0_d"), review_notes="",
                 dbic_alt=lc["shape"].get("dbic_alt"), best_alt=lc["shape"].get("best_alt"))
        out.append(r)
    df = pd.DataFrame(out)
    sample = pd.read_csv(P["sample"], usecols=["ID", "ra", "dec", "Tmag", "rad", "mass", "Teff"])
    df = df.merge(sample.rename(columns={"ID": "tic"}), on="tic", how="left")
    notes = sl.set_index("tic").review_notes.to_dict()
    df["review_notes"] = [notes.get(t) if isinstance(notes.get(t), str) else "" for t in df.tic]
    return df


# ---------------------------------------------------------------- stages

def stage_coverage(c, P, sector):
    out = os.path.join(P["work"], "coverage", f"{c['tic']}.json")
    if os.path.exists(out):
        return "cached"
    secs = ms.observed_sectors(c["tic"], c["ra"], c["dec"])
    prods = []
    for s in secs:
        p = ms.best_lc(c["tic"], c["ra"], c["dec"], c["Tmag"], s)
        prods.append(dict(sector=s, kind=p.kind if p else None,
                          n=len(p.lc.time) if p else 0,
                          cadence_min=p.cadence_min if p else None))
    jdump(out, dict(tic=c["tic"], sectors=secs, products=prods))
    return "ok"


def _product(c, s, kind):
    return ms.best_lc(c["tic"], c["ra"], c["dec"], c["Tmag"], s, kinds=(kind,))


def stage_search(c, P, sector):
    out = os.path.join(P["work"], "search", f"{c['tic']}.json")
    if os.path.exists(out):
        return "cached"
    cov = jload(os.path.join(P["work"], "coverage", f"{c['tic']}.json"))
    t14, depth, t0 = c["t14"], c["depth"], c["t0"]
    per_sector, maps = [], {}
    for pr in cov["products"]:
        s, kind = pr["sector"], pr["kind"]
        if kind is None:
            per_sector.append(dict(sector=s, kind=None))
            continue
        p = _product(c, s, kind)
        rec = dict(sector=s, kind=kind, cadence_min=p.cadence_min, n=len(p.lc.time),
                   t_start=float(p.lc.time.min()), t_end=float(p.lc.time.max()))
        # (a) the unmodified Phase 2 search
        try:
            var, cfg, ss = analyse(p.lc)
            rec["phase2_events"] = [dict(t0=e.t0, duration_h=e.duration_h,
                                         depth_ppm=e.depth * 1e6, snr=e.snr)
                                    for e in ss.search.events]
            rec["window_d"] = var.window
        except Exception as e:  # noqa: BLE001
            rec["phase2_error"] = f"{type(e).__name__}: {e}"
        # (b) matched search at the candidate's duration + exclusion map
        t, f = ms.candidate_detrend(p.lc, t14)
        dips, sig = ms.matched_dips(t, f, t14, depth)
        grid, status, _ = ms.ruled_out_map(t, f, t14, depth, sig)
        maps[s] = (grid, status)
        rec["sigma_box_ppm"] = sig * 1e6
        rec["n_excluded_grid"] = int((status == 2).sum())
        rec["n_usable_grid"] = int((status >= 1).sum())
        for d in dips:
            d["is_original"] = bool(abs(d["t0"] - t0) < max(t14, 0.25))
            if d["consistent"] and not d["is_original"]:
                d["vetting"] = vet_dip(c, s, t, f, d, p.cadence_min / 1440)
        rec["matched_dips"] = dips
        if s == sector:
            rec["original_here"] = original_confirmation(c, s, t, f, kind)
        per_sector.append(rec)
    np.savez_compressed(os.path.join(P["work"], "search", f"{c['tic']}_maps.npz"),
                        **{f"g{s}": m[0] for s, m in maps.items()},
                        **{f"s{s}": m[1] for s, m in maps.items()})
    jdump(out, dict(tic=c["tic"], sectors=per_sector))
    return "ok"


def vet_dip(c, sector, t, f, d, cad):
    """Checks 1-4 on a possible second transit."""
    res = {}
    tw, fw = v.shape_window(t, f, d["t0"], c["t14"])
    sh = v.shape_test(tw, fw, d["t0"], c["t14"], exp_time=cad)
    tr = sh["transit"]
    res["shape"] = dict(passed=sh["passed"], dbic_alt=sh["dbic_alt"], dbic_box=sh["dbic_box"],
                        best_alt=sh["best_alt"], transit=tr)
    du = v.duration_check(tr["t14_h"], tr["b"], tr["rp"], c["rad"], c["mass"], tr["grazing"])
    res["duration"] = dict(du, passed=v.duration_passes(du))
    res["edge"] = v.edge_check(t, f, tr["t0"], tr["t14_h"], tr["depth_ppm"] * 1e-6)
    ratio_t14 = tr["t14_h"] / (c["t14"] * 24)
    res["t14_ratio"] = ratio_t14
    res["depth_ratio_fit"] = tr["depth_ppm"] * 1e-6 / c["depth"]
    res["edge_pixels"] = None
    try:
        px = v.pixel_checks(c["ra"], c["dec"], c["Tmag"], sector, tr["t0"], tr["t14_h"] / 24,
                            max(tr["depth_ppm"], 1) * 1e-6, d["snr"])
        arr = px.pop("arrays", None)
        res["pixels"] = px
        if arr is not None and px["confirm"]["passed"]:
            # same TESScut route for the edge check as in Phase 3
            res["edge_pixels"] = bool(v.edge_check(arr["tl"], arr["fl"], tr["t0"], tr["t14_h"],
                                                   tr["depth_ppm"] * 1e-6)["passed"])
    except Exception as e:  # noqa: BLE001
        res["pixels"] = dict(passed=None, error=f"{type(e).__name__}: {e}")
    edge_ok = bool(res["edge"]["passed"] or res["edge_pixels"])
    res["edge_ok"] = edge_ok
    res["passed_1_4"] = bool(res["shape"]["passed"] and res["duration"]["passed"]
                             and edge_ok and res["pixels"].get("passed"))
    # Everything but the edge passes: a real transit cut short by a data gap.
    res["partial"] = bool(not edge_ok and res["shape"]["passed"] and res["duration"]["passed"]
                          and res["pixels"].get("passed"))
    res["consistent_shape"] = bool(0.5 <= ratio_t14 <= 2.0 and
                                   0.5 <= res["depth_ratio_fit"] <= 2.0)
    return res


def original_confirmation(c, sector, t, f, kind):
    """Is the original dip also in this product (2-min / QLP / TESScut)?"""
    dep, snr = v.local_dip_snr(t, f, c["t0"], c["t14"])
    return dict(kind=kind, depth_ppm=dep * 1e6, snr=snr,
                depth_ratio=dep / c["depth"] if c["depth"] else None)


def confirm_measurements(c, sector):
    """Measure the original dip in independently processed photometry of the
    same sector: SPOC 2-min (if any) and QLP (not the TESS-SPOC light curve
    the dip was found in), plus our TESScut aperture. Every product, the
    discovery light curve included, also gets the step test on its
    undetrended flux (vetting.one_sided_check)."""
    res = []
    for kind in ("spoc2min", "qlp", "tesscut", "tess-spoc"):
        p = _product(c, sector, kind)
        if p is None:
            res.append(dict(kind=kind, available=False))
            continue
        step = v.one_sided_check(p.lc.time, p.lc.flux, c["t0"], c["t14"], c["depth"])
        t, f = ms.candidate_detrend(p.lc, c["t14"])
        dep, snr = v.local_dip_snr(t, f, c["t0"], c["t14"])
        sig_pt = 1.4826 * np.median(np.abs(np.diff(f))) / np.sqrt(2)
        n_in = max(0.8 * c["t14"] / (p.cadence_min / 1440), 1)
        exp_snr = c["depth"] / (sig_pt / np.sqrt(n_in)) if sig_pt > 0 else np.nan
        res.append(dict(kind=kind, available=True, cadence_min=p.cadence_min,
                        depth_ppm=dep * 1e6, snr=snr, expected_snr=exp_snr,
                        depth_ratio=dep / c["depth"] if c["depth"] else None,
                        step=step))
    return res


def stage_confirm(c, P, sector):
    out = os.path.join(P["work"], "confirm", f"{c['tic']}.json")
    if os.path.exists(out):
        return "cached"
    jdump(out, dict(tic=c["tic"], sector=sector, measurements=confirm_measurements(c, sector)))
    return "ok"


def stage_periods(c, P, sector):
    out = os.path.join(P["work"], "periods", f"{c['tic']}.json")
    if os.path.exists(out):
        return "cached"
    z = np.load(os.path.join(P["work"], "search", f"{c['tic']}_maps.npz"))
    secs = sorted(int(k[1:]) for k in z.files if k.startswith("g"))
    maps = [(z[f"g{s}"], z[f"s{s}"]) for s in secs]
    scan = ms.period_scan(c["t0"], c["t14"], maps)
    srch = jload(os.path.join(P["work"], "search", f"{c['tic']}.json"))
    duos = [dict(sector=r["sector"], partial=d["vetting"].get("partial", False), **d)
            for r in srch["sectors"] for d in r.get("matched_dips", [])
            if d.get("vetting") and d["vetting"]["consistent_shape"]
            and (d["vetting"]["passed_1_4"] or d["vetting"].get("partial"))]
    alias = []
    for d in duos:
        for p in ms.alias_periods(c["t0"], d["t0"]):
            i = np.searchsorted(scan["P"], p)
            if 0 < i < len(scan["P"]):
                alias.append(dict(sector=d["sector"], t1=d["t0"], period=float(p),
                                  n=int(round(abs(d["t0"] - c["t0"]) / p)),
                                  allowed=bool(scan["allowed"][i - 1] and scan["allowed"][i])))
    iv = ms.allowed_intervals(scan)
    frac = float(np.mean(scan["allowed"]))
    # downsample for plotting / storage
    k = max(1, len(scan["P"]) // 20000)
    np.savez_compressed(os.path.join(P["work"], "periods", f"{c['tic']}.npz"),
                        P=scan["P"][::k], allowed=scan["allowed"].reshape(-1)[::k])
    jdump(out, dict(tic=c["tic"], allowed_fraction_log=frac, span_d=scan["span"],
                    pmax=scan["pmax"], intervals=iv, n_intervals=len(iv),
                    longest_excluded_below=_first_allowed_run(scan), duos=duos, aliases=alias,
                    allowed_aliases=[a for a in alias if a["allowed"]]))
    return "ok"


def _first_allowed_run(scan):
    """Smallest P above which every period is allowed (the 'safe' long-period floor)."""
    a = scan["allowed"]
    bad = np.where(~a)[0]
    return float(scan["P"][bad[-1] + 1]) if len(bad) and bad[-1] + 1 < len(a) else (
        float(scan["P"][0]) if not len(bad) else None)


def stage_binarity(c, P, sector):
    out = os.path.join(P["work"], "binarity", f"{c['tic']}.json")
    if os.path.exists(out):
        return "cached"
    from tesshunt import tic as tic_mod
    row = tic_mod.query_ids([c["tic"]]).to_pandas().iloc[0].to_dict()
    res = bn.binarity(row)
    jdump(out, res)
    return "ok"


STAGES = {"coverage": stage_coverage, "search": stage_search, "confirm": stage_confirm,
          "periods": stage_periods, "binarity": stage_binarity}


def stage_validate(sector):
    """Rules D1 and D5 applied to known planets: neither may drop them."""
    P = paths(sector)
    rows = []
    for c in validation_planets(sector).to_dict("records"):
        cache = os.path.join(P["work"], "validate", f"{c['tic']}.json")
        if not os.path.exists(cache):
            jdump(cache, dict(tic=c["tic"], measurements=confirm_measurements(c, sector)))
        for m in jload(cache)["measurements"]:
            st = m.get("step") or {}
            rows.append(dict(tic=c["tic"], toi=c["toi"], tfop_disp=c["tfop_disp"],
                             depth_ppm=round(c["depth"] * 1e6), t14_h=round(c["t14"] * 24, 1),
                             product=m["kind"], available=m.get("available"),
                             expected_snr=m.get("expected_snr"), depth_ratio=m.get("depth_ratio"),
                             snr=m.get("snr"), step_pre=st.get("pre"), step_post=st.get("post"),
                             one_sided=st.get("one_sided")))
    df = pd.DataFrame(rows)
    os.makedirs(P["out"], exist_ok=True)
    df.to_csv(os.path.join(P["out"], "confirm_validation.csv"), index=False, float_format="%.4g")
    from phase4_report import CONFIRM_MIN_RATIO, CONFIRM_SNR
    a = df[df.available.fillna(False).astype(bool)]
    d1 = a[a["product"].isin(["spoc2min", "qlp"]) & (a.expected_snr >= CONFIRM_SNR)]
    d1_best = d1.sort_values("expected_snr").groupby("tic").tail(1)
    st = a[a.one_sided.notna()]
    d5 = st.groupby("tic").apply(lambda g: bool((g["product"] == "tesscut").any() and len(g) >= 2
                                                and g.one_sided.astype(bool).all()))
    summ = dict(planets=int(df.tic.nunique()),
                d1_testable=int(d1_best.tic.nunique()),
                d1_false_drops=int((d1_best.depth_ratio < CONFIRM_MIN_RATIO).sum()),
                d1_depth_ratio_range=[float(d1_best.depth_ratio.min()), float(d1_best.depth_ratio.max())],
                d5_testable=int(len(d5)), d5_false_drops=int(d5.sum()),
                one_sided_products=int(st.one_sided.astype(bool).sum()), products_tested=int(len(st)))
    jdump(os.path.join(P["out"], "confirm_validation_summary.json"), summ)
    print(json.dumps(summ, indent=1))


def _run(args):
    name, c, sector = args
    P = paths(sector)
    os.makedirs(os.path.join(P["work"], name), exist_ok=True)
    try:
        return c["tic"], STAGES[name](c, P, sector)
    except Exception as e:  # noqa: BLE001
        jdump(os.path.join(P["work"], "errors", f"{name}_{c['tic']}.json"),
              dict(error=f"{type(e).__name__}: {e}", tb=traceback.format_exc()))
        return c["tic"], "error"


def run_stage(name, cands, sector, procs):
    t_start = time.time()
    jobs = [(name, c, sector) for c in cands.to_dict("records")]
    counts = {}
    # Network-bound stages are serialized by the rate limiter anyway; CPU-bound
    # 'search' benefits from a few processes.
    with Pool(procs if name == "search" else 1) as pool:
        for tic, status in pool.imap_unordered(_run, jobs):
            counts[status] = counts.get(status, 0) + 1
            print(f"[{time.strftime('%H:%M:%S')}] {name} TIC {tic}: {status}", flush=True)
    print(f"{name}: {counts} in {time.time() - t_start:.0f}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=list(STAGES) + ["validate", "rank", "report", "all"])
    ap.add_argument("--sector", type=int, default=48)
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--no-controls", action="store_true")
    ap.add_argument("--only", type=int, nargs="*", help="restrict to these TIC IDs")
    args = ap.parse_args()
    cands = candidates(args.sector, controls=not args.no_controls)
    if args.only:
        cands = cands[cands.tic.isin(args.only)]
    os.makedirs(paths(args.sector)["work"], exist_ok=True)
    stages = list(STAGES) + ["validate", "rank", "report"] if args.stage == "all" else [args.stage]
    for st in stages:
        if st in STAGES:
            run_stage(st, cands, args.sector, args.procs)
        elif st == "validate":
            stage_validate(args.sector)
        else:
            from phase4_report import rank, report
            (rank if st == "rank" else report)(cands, args.sector)


if __name__ == "__main__":
    main()
