"""Pipeline check: run the pipeline on chosen known planets and false positives.

Give it TOIs or TICs whose nature scientists have already settled (confirmed or
known planets, and false positives), or ask for a random sample of each. For
every entry it picks a TESS sector where a transit was observed, then runs the
real pipeline code on that one star, in its own workspace (work/benchmark/NAME/,
never touching sector results):

  search      is the dip found, and is it a single-dip candidate?
  checks 1-7  shape, duration, data gaps, pixels, asteroid, catalogue, FPP
  expert      Phase 5 expert checks and verdict (strong / plausible / doubtful)

The pipeline is not told the answer: check 6 runs "blind" (the star's own
entries are removed from the TOI and CTOI lists, so a known false positive is
not rejected just because ExoFOP already says so), and the known period is not
passed on. Phase 4 (other sectors) and the sector-wide pile-up test need a
whole sector and are not run.

    python scripts/benchmark.py TOI-2099.01 "TIC 377367194" 1135.01
    python scripts/benchmark.py --planets 10 --false-positives 10
    python scripts/benchmark.py --planets 5 --false-positives 5 --sector 75

Results: results/benchmark/NAME/results.csv and summary.md.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import traceback

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path[:0] = [ROOT, HERE]

from tesshunt import ffi, vetting as v  # noqa: E402
from tesshunt.survey import MIN_POINTS, analyse, edge_distance  # noqa: E402

BTJD = 2457000.0
PLANET = ("CP", "KP")
FALSE_POS = ("FP",)
MAX_SECTOR_TRIES = 6
CHECKS = [("shape", "1. dip shape"), ("duration", "2. duration fits the star"),
          ("edge", "3. not a data-gap artefact"), ("pixels", "4. on the target star (pixels)"),
          ("asteroid", "5. not an asteroid"), ("catalogue", "6. not a known eclipsing binary"),
          ("fpp", "7. false-positive probability")]


# ------------------------------------------------------------------ entries

def toi_table():
    cats = v.load_catalogues(os.path.join(ROOT, "work", "phase3"))
    return cats, cats["exofop_toi"]


def parse_entries(texts, toi):
    """'TOI-2099.01', '2099.01', 'TOI 2099' (all its planets), 'TIC 230386397' or a
    bare TIC number (every TOI on that star). Returns TOI-table rows."""
    rows, unknown = [], []
    for raw in texts:
        s = raw.strip()
        if not s:
            continue
        m_tic = re.fullmatch(r"(?i)tic[\s-]*(\d+)", s)
        num = re.sub(r"(?i)^(toi|ctoi)[\s-]*", "", s)
        if m_tic or (num.isdigit() and int(num) > 100000):
            t = toi[toi["TIC ID"] == int(m_tic.group(1) if m_tic else num)]
        elif re.fullmatch(r"\d+\.\d+", num):
            t = toi[np.isclose(toi["TOI"], float(num))]
        elif num.isdigit():
            t = toi[np.floor(toi["TOI"]) == int(num)]
        else:
            t = toi.iloc[0:0]
        if len(t):
            rows += [r for _, r in t.iterrows()]
        else:
            unknown.append(s)
    return rows, unknown


def sample(toi, n_planets, n_fp, min_period, seed):
    """Random known planets (TFOPWG CP/KP) and false positives (FP) with a usable
    ephemeris, long enough periods (this pipeline hunts single transits) and
    stars bright enough for TESS-SPOC light curves."""
    ok = (toi["Period (days)"].fillna(0) >= min_period) & toi["Epoch (BJD)"].notna() \
        & (toi["TESS Mag"].fillna(99) <= 13.0) & toi["Sectors"].notna()
    rng = np.random.default_rng(seed)
    out = []
    for disp, n in ((PLANET, n_planets), (FALSE_POS, n_fp)):
        g = toi[ok & toi["TFOPWG Disposition"].isin(disp)]
        idx = rng.choice(len(g), size=min(n, len(g)), replace=False) if n else []
        out += [g.iloc[i] for i in idx]
    return out


def truth(disp):
    return "planet" if disp in PLANET else ("false positive" if disp in FALSE_POS or disp == "FA"
                                            else "undecided")


# ------------------------------------------------------------------ sector choice

def predicted(t0_btjd, period, lo, hi):
    if not (np.isfinite(period) and period > 0):
        return [t0_btjd] if lo <= t0_btjd <= hi else []
    n = np.arange(np.ceil((lo - t0_btjd) / period), np.floor((hi - t0_btjd) / period) + 1)
    return list(t0_btjd + n * period)


def choose_sector(r, sector=None):
    """A sector with a TESS-SPOC light curve covering a predicted transit; prefer
    one with a single transit (what this pipeline is built for), newest first."""
    tic = int(r["TIC ID"])
    t0, per = float(r["Epoch (BJD)"]) - BTJD, float(r["Period (days)"] or np.nan)
    dur = float(r["Duration (hours)"]) / 24 if np.isfinite(r["Duration (hours)"]) else 0.2
    secs = [sector] if sector else sorted({int(x) for x in str(r["Sectors"]).split(",")
                                           if x.strip().isdigit()}, reverse=True)
    best, tried = None, 0
    for s in secs:
        if tried >= MAX_SECTOR_TRIES:
            break
        tried += 1
        try:
            path = ffi.download(tic, s, cache=True)
        except FileNotFoundError:
            continue
        lc, info = ffi.read(path)
        if len(lc.time) < MIN_POINTS:
            continue
        tt = lc.time
        cov = [p for p in predicted(t0, per, tt.min(), tt.max())
               if np.sum(np.abs(tt - p) < dur / 2) >= 3]
        if not cov:
            continue
        cand = (len(cov) != 1, s, lc, info, cov)
        if best is None or cand[0] < best[0]:
            best = cand
        if len(cov) == 1:
            break
    return (None, None, None, []) if best is None else best[1:]


# ------------------------------------------------------------------ one entry

def blind(cats, tic):
    """The catalogues without this star's own TOI/CTOI entries."""
    out = dict(cats)
    out["exofop_toi"] = cats["exofop_toi"][cats["exofop_toi"]["TIC ID"] != tic]
    out["exofop_ctoi"] = cats["exofop_ctoi"][cats["exofop_ctoi"]["TIC ID"] != tic]
    return out


def star_info(tic):
    from tesshunt import tic as tic_mod
    t = tic_mod.query_ids([tic])
    if not len(t):
        return {}
    r = t[0]
    g = str(r["GAIA"])

    def f(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan
    return dict(ra=f(r["ra"]), dec=f(r["dec"]), Tmag=f(r["Tmag"]), rad=f(r["rad"]),
                mass=f(r["mass"]), gaia_id=int(g) if g.isdigit() else None)


def run_entry(r, cats, work, sector=None):
    import phase3_vet as p3
    import phase5 as p5
    tic, disp = int(r["TIC ID"]), str(r["TFOPWG Disposition"] if pd.notna(r["TFOPWG Disposition"]) else "")
    out = dict(toi=f"{r['TOI']:.2f}", tic=tic, disp=disp, truth=truth(disp),
               period_d=r["Period (days)"], depth_ppm=r["Depth (ppm)"], duration_h=r["Duration (hours)"],
               tmag=r["TESS Mag"])
    s, lc, info, cov = choose_sector(r, sector)
    if s is None:
        return dict(out, stage="no data", outcome="no TESS-SPOC light curve with a transit in it "
                    f"(tried up to {MAX_SECTOR_TRIES} sectors)")
    out.update(sector=s, n_transits_in_sector=len(cov), t_expected=round(cov[0], 4))
    # --- search (as in survey.process_star)
    var, cfg, ss = analyse(lc)
    ev = ss.search.events
    dur = float(r["Duration (hours)"]) if np.isfinite(r["Duration (hours)"]) else 3.0
    tol = max(dur / 24, 0.1)
    d = pd.DataFrame([dict(tic=tic, t0=e.t0, duration_h=e.duration_h, depth_ppm=e.depth * 1e6,
                           snr=e.snr, edge_dist_h=edge_distance(ss.time, e.t0)) for e in ev])
    if d.empty:
        return dict(out, stage="search", outcome="missed: no dip found in this sector")
    d["common_mode"] = False
    d["partial"] = d.edge_dist_h < d.duration_h / 2 + 1
    d["category"] = v.classify(d)
    hit = d[np.min(np.abs(d.t0.values[:, None] - np.array(cov)[None, :]), axis=1) < tol]
    if hit.empty:
        return dict(out, stage="search", outcome="missed: no dip at the expected time",
                    n_dips=len(d))
    h = hit.sort_values("snr").iloc[-1]
    candidate = h.category in ("single", "single_partial") and not var.hf_variable
    out.update(found=True, snr=round(float(h.snr), 1), found_depth_ppm=round(float(h.depth_ppm)),
               category=h.category, candidate=bool(candidate))
    # --- vetting checks 1-7 (Phase 3 code, own workspace, blind catalogue)
    star = star_info(tic)
    p3.configure(s)
    p3.WORK = os.path.join(work, "phase3", f"s{s:04d}")
    p3._CATS = blind(cats, tic)
    row = dict(key=p3.key(tic, float(h.t0)), tic=tic, t0=float(h.t0), duration_h=float(h.duration_h),
               snr=float(h.snr), is_validation=False, **{k: star.get(k) for k in
                                                         ("ra", "dec", "Tmag", "rad", "mass")})
    k = row["key"]
    p3.lc_job(row)
    lcr = p3.load_json(p3.jpath("lc", k))
    if lcr["status"] != "ok":
        return dict(out, stage="checks", outcome=f"check 1 could not run ({lcr.get('error')})")
    res = {"shape": bool(v.shape_passes(lcr["shape"]["dbic_alt"], lcr["shape"]["dbic_box"])),
           "duration": bool(v.duration_passes(lcr["duration"]))}
    p3.pixel_job(row)
    pix = p3.load_json(p3.jpath("pix", k))
    lc_edge, pix_edge = p3.edge_results(k, row["snr"])
    res["edge"] = bool(lc_edge or pix_edge)
    if pix["status"] == "ok":
        ap_depth, ap_snr, lc_depth = p3.pixel_dip(k)
        pixd = dict(pix["pixels"], neighbours=p3.neighbours_now(k, pix["pixels"]))
        res["pixels"] = bool(v.pixel_passes(pixd, row["snr"], ap_snr, ap_depth, lc_depth))
        res["asteroid"] = pix["asteroid"].get("passed")
        res["catalogue"] = bool(v.catalogue_check_passes(pix["catalogue"]))
    p3.fpp_job(row)
    fr = p3.load_json(p3.jpath("fpp", k))
    if fr["status"] == "ok":
        res["fpp"] = bool(fr["fpp"]["passed"])
        out["fpp"] = round(float(fr["fpp"]["fpp"]), 3)
    for c, _ in CHECKS:
        out[f"check_{c}"] = res.get(c)
    # --- expert checks (Phase 5 code, own workspace, no known period)
    tr = lcr["shape"]["transit"]
    p5.WORK = os.path.join(work, "phase5")
    t = dict(key=f"bench_{k}", tic=tic, sector=s, role="benchmark", category="", reasons="",
             t0=tr["t0"], t14=tr["t14_h"] / 24, depth=max(tr["depth_ppm"], 1.0) * 1e-6, snr=row["snr"],
             fpp=out.get("fpp"), periods=None, binarity_flags="", **star)
    try:
        a = p5.assess(t, p5.run_target(t))
        out.update(verdict=a["verdict"], serious=";".join(a["serious"]), minor=";".join(a["minor"]))
    except Exception as e:  # noqa: BLE001
        out.update(verdict=None, expert_error=f"{type(e).__name__}: {e}")
    # --- the pipeline's decision, in its own order
    if not candidate:
        why = "the star is too variable" if var.hf_variable else f"it looks like a {h.category} signal"
        out.update(stage="search", outcome=f"found (SNR {h.snr:.0f}), but not kept as a candidate: {why}")
        return out
    failed = [name for c, name in CHECKS if res.get(c) is False]
    if failed:
        out.update(stage="checks", outcome=f"rejected at check {failed[0]}")
    elif any(res.get(c) is None for c, _ in CHECKS):
        missing = [name for c, name in CHECKS if res.get(c) is None]
        out.update(stage="checks", outcome=f"passed the checks that ran; could not run {', '.join(missing)}")
    else:
        out.update(stage="expert", outcome=f"passed checks 1-7; expert verdict {out.get('verdict') or '?'}")
    return out


def kept(o):
    """Would the pipeline put this forward as a planet candidate?"""
    return (o.get("candidate") is True and o.get("stage") == "expert"
            and o.get("verdict") in ("strong", "plausible"))


def scorecard(df):
    """Plain-English score lines: for planets, how many were kept; for false
    positives, how many the checks rejected (one the search never found is
    counted separately: nothing was tested)."""
    lines = []
    for label in ("planet", "false positive"):
        g = df[df.truth == label]
        if g.empty:
            continue
        nodata = int((g.stage == "no data").sum())
        tested = g[g.stage != "no data"]
        found = tested[tested.get("found", pd.Series(False, index=tested.index)).fillna(False).astype(bool)]
        k = int(tested.apply(kept, axis=1).astype(bool).sum()) if len(tested) else 0
        missed = len(tested) - len(found)
        extra = [f"{missed} not found by the search"] if missed else []
        if nodata:
            extra.append(f"{nodata} with no usable data")
        if label == "planet":
            head = f"{k} of {len(tested)} kept as candidates"
        else:
            head = f"{len(found) - k} of {len(found)} found were rejected or rated doubtful"
        lines.append(f"- **Known {label}s:** {head}" + (f" ({'; '.join(extra)})" if extra else ""))
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("entries", nargs="*", help="TOIs (e.g. TOI-2099.01) or TICs (e.g. 'TIC 230386397')")
    ap.add_argument("--planets", type=int, default=0, help="random known planets (CP/KP) to add")
    ap.add_argument("--false-positives", type=int, default=0, help="random known false positives to add")
    ap.add_argument("--min-period", type=float, default=15.0,
                    help="days; for the random sample (this pipeline hunts single transits)")
    ap.add_argument("--sector", type=int, default=None, help="use this sector instead of choosing one")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--name", default=time.strftime("%Y%m%d-%H%M%S"))
    args = ap.parse_args()
    cats, toi = toi_table()
    rows, unknown = parse_entries(args.entries, toi)
    rows += sample(toi, args.planets, args.false_positives, args.min_period, args.seed)
    for u in unknown:
        print(f"not in the TOI list, skipped: {u}", flush=True)
    if not rows:
        sys.exit("nothing to check: give TOIs/TICs, or --planets / --false-positives")
    work = os.path.join(ROOT, "work", "benchmark", args.name)
    outdir = os.path.join(ROOT, "results", "benchmark", args.name)
    os.makedirs(outdir, exist_ok=True)
    print(f"\n=== [{time.strftime('%H:%M:%S')}] benchmark: {len(rows)} entries -> "
          f"{os.path.relpath(outdir, ROOT)}", flush=True)
    results = []
    for i, r in enumerate(rows, 1):
        t_start = time.time()
        try:
            o = run_entry(r, cats, work, args.sector)
        except Exception as e:  # noqa: BLE001
            from tesshunt import net
            if isinstance(e, net.ServiceHalt):
                raise
            o = dict(toi=f"{r['TOI']:.2f}", tic=int(r["TIC ID"]), disp=str(r["TFOPWG Disposition"]),
                     truth=truth(str(r["TFOPWG Disposition"])), stage="error",
                     outcome=f"{type(e).__name__}: {e}", tb=traceback.format_exc())
        results.append(o)
        pd.DataFrame(results).to_csv(os.path.join(outdir, "results.csv"), index=False)
        print(f"[{time.strftime('%H:%M:%S')}] {i}/{len(rows)} TOI-{o['toi']} ({o['truth']}): "
              f"{o['outcome']} ({time.time() - t_start:.0f}s)", flush=True)
    df = pd.DataFrame(results)
    md = [f"# Pipeline check {args.name}", "", *scorecard(df), "",
          "| TOI | TIC | truth | sector | outcome |", "|---|---|---|---|---|"]
    md += [f"| {o['toi']} | {o['tic']} | {o['truth']} | {o.get('sector', '')} | {o['outcome']} |"
           for o in results]
    with open(os.path.join(outdir, "summary.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    with open(os.path.join(outdir, "settings.json"), "w") as fh:
        json.dump(vars(args), fh, indent=1)
    print("\n".join(scorecard(df)), flush=True)
    print("=== benchmark done", flush=True)


if __name__ == "__main__":
    main()
