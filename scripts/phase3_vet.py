"""Phase 3: vet the Phase 2 single-dip candidates.

    python scripts/phase3_vet.py lc        # checks 1-3 on every candidate (+ validation TOIs)
    python scripts/phase3_vet.py pixels    # checks 4-6 on survivors of 1-3 (+ validation)
    python scripts/phase3_vet.py fpp       # check 7 (TRICERATOPS) on survivors of 1-6 (+ validation)
    python scripts/phase3_vet.py report    # funnel, shortlist, vetting sheets, summary

Each stage is resumable: per-candidate results are cached as JSON/NPZ in
work/phase3/ (git-ignored) and skipped on rerun. Light curves and pixel
cutouts are fetched, reduced, and discarded; only derived arrays are cached.

Checks (tesshunt/vetting.py): 1 shape, 2 duration, 3 edge, 4 pixels
(centroid + neighbours), 5 asteroids, 6 catalogues, 7 FPP.
The validation set is TOI-2180 b plus every TOI with a single predicted
transit in the sector that Phase 2 recovered; it goes through every check
regardless of earlier results, to measure how often each check would reject
a real planet.
"""

import argparse
import json
import os
import sys
import time
import traceback
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

# Run as a script, this module is "__main__"; phase3_report does `import phase3_vet`,
# which would load a second, unconfigured copy (still pointing at Sector 48's
# folders). Register this module under its own name so both share one state.
if __name__ == "__main__":
    sys.modules.setdefault("phase3_vet", sys.modules[__name__])

from tesshunt import ffi, vetting as v  # noqa: E402
from tesshunt.survey import analyse  # noqa: E402

P2 = os.path.join(ROOT, "results", "phase2")


def sector_dirs(sector, phase):
    """(work, results, plots) for a phase and sector. Sector 48 predates
    multi-sector support and keeps the top-level folders; other sectors get
    an sXXXX subfolder."""
    sub = "" if sector == 48 else f"s{sector:04d}"
    return (os.path.join(ROOT, "work", f"phase{phase}", sub),
            os.path.join(ROOT, "results", f"phase{phase}", sub),
            os.path.join(ROOT, "plots", f"phase{phase}", sub))


def configure(sector):
    global SECTOR, TAG, WORK, OUT, PLOTS
    SECTOR, TAG = sector, f"s{sector:04d}"
    WORK, OUT, PLOTS = sector_dirs(sector, 3)


configure(48)
CHECKS = ["shape", "duration", "edge", "pixels", "asteroid", "catalogue", "fpp"]


def key(tic, t0):
    return f"{int(tic)}_{t0:.4f}"


def jpath(stage, k, ext="json"):
    d = os.path.join(WORK, stage)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{k}.{ext}")


def save_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=_jsonable)
    os.replace(tmp, path)


def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


# ---------------------------------------------------------------- targets

def targets():
    """Candidates (tier A first, then B) plus the validation TOIs."""
    cand = pd.read_csv(os.path.join(P2, f"{TAG}_candidates.csv"))
    dips = pd.read_csv(os.path.join(P2, f"{TAG}_dips.csv"))
    rec = pd.read_csv(os.path.join(P2, f"{TAG}_toi_recall.csv"))
    rec = rec[(rec.n_pred == 1) & rec.recovered]
    val = []
    for _, r in rec.iterrows():
        g = dips[(dips.tic == r.tic) & (dips.toi_match.astype(str) == f"{r.toi:.2f}")]
        if len(g):
            val.append(g.sort_values("snr").iloc[-1])
    val = pd.DataFrame(val)
    cand["is_candidate"] = True
    val["is_validation"] = True
    both = pd.concat([cand, val]).drop_duplicates(["tic", "t0"], keep="first")
    both = both.merge(val[["tic", "t0", "is_validation"]], on=["tic", "t0"], how="left",
                      suffixes=("_x", ""))
    both = both.drop(columns=[c for c in both.columns if c.endswith("_x")])
    both["is_candidate"] = both.is_candidate.fillna(False).astype(bool)
    both["is_validation"] = both.is_validation.fillna(False).astype(bool)
    sample = pd.read_csv(os.path.join(ROOT, "work", f"{TAG}_sample.csv"),
                         usecols=["ID", "ra", "dec", "mass", "rad", "Tmag"])
    both = both.drop(columns=[c for c in ("ra", "dec") if c in both.columns])
    both = both.merge(sample.rename(columns={"ID": "tic"}), on="tic", how="left")
    both["key"] = [key(a, b) for a, b in zip(both.tic, both.t0)]
    order = {"A": 0, "B": 1}
    both["order"] = [order.get(t, 2) for t in both.tier.fillna("")]
    return both.sort_values(["order", "rank"]).reset_index(drop=True)


# ---------------------------------------------------------------- stage: lc

def lc_job(row):
    k = row["key"]
    out = jpath("lc", k)
    if os.path.exists(out):
        return k, "cached"
    try:
        lc, _ = ffi.read(ffi.download(int(row["tic"]), SECTOR, cache=True))
        var, cfg, ss = analyse(lc)
        t, f = ss.time, ss.flat
        cad = float(np.median(np.diff(t)))
        dur = row["duration_h"] / 24
        tw, fw = v.shape_window(t, f, row["t0"], dur)
        sh = v.shape_test(tw, fw, row["t0"], dur, exp_time=cad)
        tr = sh["transit"]
        du = v.duration_check(tr["t14_h"], tr["b"], tr["rp"], row["rad"], row["mass"],
                              tr["grazing"])
        ed = v.edge_check(t, f, tr["t0"], tr["t14_h"], tr["depth_ppm"] * 1e-6)
        np.savez_compressed(jpath("lc", k, "npz"), t=t.astype("f8"), f=f.astype("f4"),
                            tw=tw, fw=fw,
                            **{f"m_{n}": m.model.astype("f4") for n, m in sh["fits"].items()})
        sh = {kk: vv for kk, vv in sh.items() if kk != "fits"}
        save_json(out, dict(key=k, tic=int(row["tic"]), t0_det=row["t0"], shape=sh,
                            duration=du, edge=ed, cadence=cad, t_first=float(t[0]),
                            t_last=float(t[-1]), window_d=var.window, status="ok"))
        return k, "ok"
    except Exception as e:  # noqa: BLE001
        save_json(out, dict(key=k, tic=int(row["tic"]), status="error",
                            error=f"{type(e).__name__}: {e}", tb=traceback.format_exc()))
        return k, "error"


def lc_ok(k):
    """Checks 1-2 (shape, duration) from the cached light-curve result."""
    p = jpath("lc", k)
    if not os.path.exists(p):
        return None
    r = load_json(p)
    if r["status"] != "ok":
        return None
    return v.shape_passes(r["shape"]["dbic_alt"], r["shape"]["dbic_box"]) and \
        v.duration_passes(r["duration"])


def pixel_dip(k):
    """(depth, SNR, light-curve depth) of the dip in the cached TESScut
    aperture light curve, measured against local baselines
    (vetting.local_dip_snr); the last value is the transit-fit depth."""
    tr = load_json(jpath("lc", k))["shape"]["transit"]
    lc_depth = tr["depth_ppm"] * 1e-6
    pz = jpath("pix", k, "npz")
    if not os.path.exists(pz):
        return np.nan, np.nan, lc_depth
    d = np.load(pz)
    dep, snr = v.local_dip_snr(d["tl"], d["fl"], tr["t0"], tr["t14_h"] / 24)
    return dep, snr, lc_depth


def neighbours_now(k, pix):
    """Re-run the neighbour test from the cached difference image, so rule
    changes (e.g. NEIGHBOUR_MIN_PX) apply without re-downloading."""
    pz = jpath("pix", k, "npz")
    if not os.path.exists(pz):
        return pix["neighbours"]
    d = np.load(pz)
    dep = max(load_json(jpath("lc", k))["shape"]["transit"]["depth_ppm"], 1) * 1e-6
    return v.neighbour_test(d["diff"], d["noise"], d["oot"], float(d["x"]), float(d["y"]),
                            [tuple(r) for r in d["nbs"]], dep)


def edge_results(k, lc_snr=None):
    """Check 3 on the PDCSAP light curve and, if a pixel cutout exists, on the
    TESScut aperture light curve (which keeps cadences PDCSAP blanks near
    orbit boundaries). The pixel route only counts if the dip is actually
    confirmed in those pixels (check 4c). Returns (lc_edge, pixel_edge or None)."""
    r = load_json(jpath("lc", k))
    lc_edge = r["edge"]["passed"]
    pz = jpath("pix", k, "npz")
    if not os.path.exists(pz):
        return lc_edge, None
    d = np.load(pz)
    tr = r["shape"]["transit"]
    pe = v.edge_check(d["tl"], d["fl"], tr["t0"], tr["t14_h"], tr["depth_ppm"] * 1e-6)
    ap_depth, ap_snr, lc_depth = pixel_dip(k)
    snr = r.get("lc_snr") if lc_snr is None else lc_snr
    confirmed = snr is not None and v.pixel_confirm_check(ap_snr, snr, ap_depth,
                                                          lc_depth)["passed"]
    return lc_edge, bool(pe["passed"] and confirmed)


def passed_lc(k, lc_snr=None):
    ok = lc_ok(k)
    if not ok:
        return ok
    lc_edge, pix_edge = edge_results(k, lc_snr)
    return bool(lc_edge or pix_edge)


# ---------------------------------------------------------------- stage: pixels

_CATS = None


def pixel_job(row):
    global _CATS
    k = row["key"]
    out = jpath("pix", k)
    if os.path.exists(out):
        r = load_json(out)
        if r["status"] == "ok" and (r["asteroid"].get("passed") is None
                                    or r["asteroid"].get("error")):
            return k, retry_asteroid(row, r, out)
        return k, "cached"
    try:
        lcr = load_json(jpath("lc", k))
        tr = lcr["shape"]["transit"]
        t0, t14 = tr["t0"], tr["t14_h"] / 24
        depth = max(tr["depth_ppm"], 1.0) * 1e-6
        cut = v.tesscut_cutout(row["ra"], row["dec"], SECTOR)
        di = v.difference_image(cut["time"], cut["cube"], t0, t14)
        if di is None:
            raise RuntimeError("not enough cadences for a difference image")
        diff, noise, oot, n_in, n_oot = di
        nbs = v.tic_neighbours(row["ra"], row["dec"], cut["wcs"], row["Tmag"])
        ct = v.centroid_test(diff, noise, cut["x"], cut["y"])
        nt = v.neighbour_test(diff, noise, oot, cut["x"], cut["y"], nbs, depth)
        ap = v.aperture_mask(oot, cut["x"], cut["y"])
        # aperture light curve around the dip, normalized by a linear OOT baseline
        sel = np.abs(cut["time"] - t0) < max(2.5 * t14, 0.5)
        tl = cut["time"][sel]
        fl = cut["cube"][sel][:, ap].sum(axis=1)
        oo = np.abs(tl - t0) > t14 / 2 + 1 / 24
        c = np.polyfit(tl[oo] - t0, fl[oo], 1) if oo.sum() > 3 else [0, np.median(fl)]
        fl = fl / np.polyval(c, tl - t0)
        ap_depth, ap_snr = v.local_dip_snr(tl, fl, t0, t14)
        pix_depth_args = (ap_depth, tr["depth_ppm"] * 1e-6)
        pix = dict(centroid=ct, neighbours=nt, n_in=n_in, n_oot=n_oot,
                   ap_depth_ppm=ap_depth * 1e6, ap_snr=ap_snr, n_ap=int(ap.sum()))
        pix["confirm"] = v.pixel_confirm_check(ap_snr, row["snr"], *pix_depth_args)
        pix["passed"] = v.pixel_passes(pix, row["snr"], ap_snr, *pix_depth_args)
        ap_abs = [[cut["col0"] + int(ix), cut["row0"] + int(iy)]
                  for iy, ix in zip(*np.nonzero(ap))]
        np.savez_compressed(jpath("pix", k, "npz"), diff=diff, noise=noise, oot=oot, ap=ap,
                            x=cut["x"], y=cut["y"], tl=tl, fl=fl,
                            nbs=np.array([n for n in nbs], float).reshape(-1, 4))
        # 5: asteroids
        try:
            objs = v.skybot_query(row["ra"], row["dec"], t0 + 2457000.0)
            ast = v.asteroid_check(objs, row["ra"], row["dec"], t14 * 24)
        except Exception as e:  # noqa: BLE001
            ast = dict(passed=None, error=f"{type(e).__name__}: {e}", hits=[], n_objects=None)
        # 6: catalogues
        if _CATS is None:
            _CATS = v.load_catalogues(WORK)
        cat = v.catalogue_check(int(row["tic"]), _CATS)
        save_json(out, dict(key=k, status="ok", pixels=pix, asteroid=ast, catalogue=cat,
                            ap_abs=ap_abs))
        return k, "ok"
    except Exception as e:  # noqa: BLE001
        save_json(out, dict(key=k, status="error", error=f"{type(e).__name__}: {e}",
                            tb=traceback.format_exc()))
        return k, "error"


def retry_asteroid(row, r, out):
    lcr = load_json(jpath("lc", row["key"]))
    tr = lcr["shape"]["transit"]
    try:
        objs = v.skybot_query(row["ra"], row["dec"], tr["t0"] + 2457000.0)
        r["asteroid"] = v.asteroid_check(objs, row["ra"], row["dec"], tr["t14_h"])
    except Exception as e:  # noqa: BLE001
        r["asteroid"].update(passed=None, error=f"{type(e).__name__}: {e}")
        save_json(out, r)
        return "error"
    save_json(out, r)
    return "asteroid retried"


def passed_pix(k, lc_snr=None):
    p = jpath("pix", k)
    if not os.path.exists(p):
        return None
    r = load_json(p)
    if r["status"] != "ok" or r["asteroid"].get("passed") is None or r["asteroid"].get("error"):
        return None
    ap_depth, ap_snr, lc_depth = pixel_dip(k)
    pixd = dict(r["pixels"], neighbours=neighbours_now(k, r["pixels"]))
    pix = (v.pixel_passes(pixd, lc_snr, ap_snr, ap_depth, lc_depth)
           if lc_snr is not None else r["pixels"]["passed"])
    return pix and r["asteroid"]["passed"] and v.catalogue_check_passes(r["catalogue"])


# ---------------------------------------------------------------- stage: fpp

def fpp_job(row):
    k = row["key"]
    out = jpath("fpp", k)
    if os.path.exists(out):
        return k, "cached"
    try:
        lcr = load_json(jpath("lc", k))
        pix = load_json(jpath("pix", k))
        tr = lcr["shape"]["transit"]
        t0 = tr["t0"]
        d = np.load(jpath("lc", k, "npz"))
        tw, fw = d["tw"], d["fw"].astype(float)
        err = v._robust_sigma_pt(fw)
        # Single transit: no second transit inside the data bounds the period
        # from below; the duration (b = 0, circular) bounds it from above.
        p_lo = max(t0 - lcr["t_first"], lcr["t_last"] - t0)
        p_b0 = lcr["duration"].get("p_b0_d", np.nan)
        p_hi = max(p_b0 if np.isfinite(p_b0) else 3 * p_lo, 1.5 * p_lo)
        res = v.triceratops_fpp(int(row["tic"]), row["ra"], row["dec"], SECTOR, tw - t0, fw,
                                err, max(tr["depth_ppm"], 1.0) * 1e-6, pix["ap_abs"],
                                (p_lo, p_hi), os.path.join(WORK, "tri", k))
        res["p_range"] = [p_lo, p_hi]
        save_json(out, dict(key=k, status="ok", fpp=v.fpp_check(res)))
        return k, "ok"
    except Exception as e:  # noqa: BLE001
        save_json(out, dict(key=k, status="error", error=f"{type(e).__name__}: {e}",
                            tb=traceback.format_exc()))
        return k, "error"


def _init_worker():
    os.environ.setdefault("OMP_NUM_THREADS", "1")


def run_stage(fn, rows, procs, threads=False, label=""):
    jobs = [r for r in rows]
    if not jobs:
        print(f"{label}: nothing to do")
        return
    t_start, n, counts = time.time(), 0, {}
    pool = ThreadPool(procs) if threads else Pool(procs, initializer=_init_worker)
    with pool:
        for k, status in pool.imap_unordered(fn, jobs):
            n += 1
            counts[status] = counts.get(status, 0) + 1
            if n % 50 == 0 or n == len(jobs):
                rate = n / (time.time() - t_start)
                print(f"[{time.strftime('%H:%M:%S')}] {label} {n}/{len(jobs)} {counts} "
                      f"{rate:.2f}/s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["lc", "pixels", "fpp", "report"])
    ap.add_argument("--sector", type=int, default=48)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--retry-errors", action="store_true")
    args = ap.parse_args()
    configure(args.sector)
    tg = targets()
    rows = tg.to_dict("records")
    if args.stage == "report":
        from phase3_report import report
        report(tg)
        return
    stage = {"lc": "lc", "pixels": "pix", "fpp": "fpp"}[args.stage]
    if args.retry_errors:
        for r in rows:
            p = jpath(stage, r["key"])
            if os.path.exists(p) and load_json(p)["status"] == "error":
                os.remove(p)
    if args.stage == "lc":
        todo = rows
        fn, threads = lc_job, False
    elif args.stage == "pixels":
        # Pixels are needed for everything passing checks 1-2: the edge check
        # can also be satisfied by the TESScut aperture light curve.
        todo = [r for r in rows if lc_ok(r["key"]) or
                (r["is_validation"] and lc_ok(r["key"]) is not None)]
        fn, threads = pixel_job, True
    else:
        todo = [r for r in rows if (passed_lc(r["key"], r["snr"]) and passed_pix(r["key"], r["snr"])) or
                (r["is_validation"] and passed_pix(r["key"], r["snr"]) is not None)]
        fn, threads = fpp_job, False
    if args.stage != "pixels":
        todo = [r for r in todo if not os.path.exists(jpath(stage, r["key"]))]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{args.stage}: {len(tg)} targets ({int(tg.is_candidate.sum())} candidates, "
          f"{int(tg.is_validation.sum())} validation), {len(todo)} to process")
    run_stage(fn, todo, args.procs, threads, args.stage)


if __name__ == "__main__":
    main()
