"""Phase 4 item 4: run the Phase 2 injected transits through all seven vetting checks.

    python scripts/phase4_injection_vetting.py prep     # per-star downloads (cached, polite)
    python scripts/phase4_injection_vetting.py checks   # checks 1-4 (+6) on every recovered injection
    python scripts/phase4_injection_vetting.py sky      # check 5 (SkyBoT) for survivors of 1-4
    python scripts/phase4_injection_vetting.py fpp      # check 7 on an SNR-stratified subset
    python scripts/phase4_injection_vetting.py report

This is an independent test of rules that were tuned on the validation TOIs:
the injections were made in Phase 2, before any vetting rule existed.

Injections are trapezoids (10 % ingress) put into the raw PDCSAP flux exactly
as in Phase 2. For the pixel checks the same transit is injected into the
TESScut cutout as a PSF-shaped signal at the target position,
cube - (1 - m(t)) * F * G(x, y), with G a 0.75 px Gaussian and F the target
flux fitted to the out-of-transit image, so a real on-target transit is what
the pixel checks see.
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

from tesshunt import ffi, vetting as v  # noqa: E402
from tesshunt.injection import inject, trapezoid  # noqa: E402
from tesshunt.survey import INGRESS_FRAC, analyse  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from phase3_vet import sector_dirs  # noqa: E402


def configure(sector):
    global SECTOR, WORK, OUT, PLOTS, P3WORK
    SECTOR = sector
    WORK = os.path.join(ROOT, "work", "phase4", f"s{sector:04d}", "injvet")
    _, OUT, PLOTS = sector_dirs(sector, 4)
    P3WORK = sector_dirs(sector, 3)[0]


configure(48)
FPP_SUBSET = 150
SNR_BINS = [0, 7, 8, 9, 10, 12, 15, 20, 30, 50, 1e9]


def jdump(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    os.replace(tmp, path)


def jload(path):
    with open(path) as fh:
        return json.load(fh)


def trials():
    inj = pd.read_csv(os.path.join(ROOT, "results", "phase2", f"s{SECTOR:04d}_injections.csv"))
    sample = pd.read_csv(os.path.join(ROOT, "work", f"s{SECTOR:04d}_sample.csv"),
                         usecols=["ID", "ra", "dec", "Tmag", "rad", "mass"])
    inj = inj.merge(sample.rename(columns={"ID": "tic"}), on="tic", how="left")
    inj["key"] = [f"{t}_{k}" for t, k in zip(inj.tic, inj.trial)]
    # physically plausible single transits: long enough for P >= 20 d, short
    # enough for P <= 100 yr (the duration check's own limits)
    plaus = []
    for r in inj.itertuples():
        du = v.duration_check(r.duration_h, 0.0, np.sqrt(r.depth_ppm * 1e-6), r.rad, r.mass, False)
        plaus.append(bool(v.duration_passes(du)))
    inj["plausible"] = plaus
    return inj


# ---------------------------------------------------------------- prep

def prep(inj):
    """Sequential, cached downloads per star: light curve (S3 first), TESScut
    cutout, TIC neighbours. Network access is rate-limited in tesshunt.net."""
    stars = inj.drop_duplicates("tic")
    for i, r in enumerate(stars.itertuples(), 1):
        try:
            ffi.download(int(r.tic), SECTOR, cache=True)
            cut = v.tesscut_cutout(r.ra, r.dec, SECTOR)
            v.tic_neighbours(r.ra, r.dec, cut["wcs"], r.Tmag)
            status = "ok"
        except Exception as e:  # noqa: BLE001
            status = f"{type(e).__name__}: {e}"
        print(f"[{time.strftime('%H:%M:%S')}] prep {i}/{len(stars)} TIC {r.tic}: {status}",
              flush=True)
    v.load_catalogues(P3WORK)


# ---------------------------------------------------------------- checks 1-4 (+6)

def psf_inject(cut, t0, depth, dur):
    """Inject the transit into the cutout as a PSF at the target position."""
    ny, nx = cut["cube"].shape[1:]
    yy, xx = np.mgrid[0:ny, 0:nx]
    g = np.exp(-((xx - cut["x"]) ** 2 + (yy - cut["y"]) ** 2) / (2 * 0.75 ** 2))
    g /= g.sum()
    img = np.median(cut["cube"], axis=0)
    F = float((img * g).sum() / (g * g).sum())
    m = trapezoid(cut["time"], t0, depth, dur, INGRESS_FRAC)
    return cut["cube"] - (1 - m)[:, None, None] * F * g[None, :, :]


_CATS = None


def check_trial(r):
    global _CATS
    out = os.path.join(WORK, "checks", f"{r['key']}.json")
    if os.path.exists(out):
        return r["key"], "cached"
    try:
        tic, t0 = int(r["tic"]), r["t0"]
        depth, dur = r["depth_ppm"] * 1e-6, r["duration_h"] / 24
        lc, _ = ffi.read(ffi.download(tic, SECTOR, cache=True))
        var, cfg, ss = analyse(inject(lc, t0, depth, dur, INGRESS_FRAC))
        tol = max(dur / 2, 0.5 / 24)
        hits = [e for e in ss.search.events if abs(e.t0 - t0) < tol]
        if not hits:
            jdump(out, dict(key=r["key"], detected=False))
            return r["key"], "not detected"
        e = max(hits, key=lambda e: e.snr)
        t, f = ss.time, ss.flat
        cad = float(np.median(np.diff(t)))
        tw, fw = v.shape_window(t, f, e.t0, e.duration)
        sh = v.shape_test(tw, fw, e.t0, e.duration, exp_time=cad)
        tr = sh["transit"]
        du = v.duration_check(tr["t14_h"], tr["b"], tr["rp"], r["rad"], r["mass"], tr["grazing"])
        ed = v.edge_check(t, f, tr["t0"], tr["t14_h"], tr["depth_ppm"] * 1e-6)
        res = dict(key=r["key"], detected=True, det_snr=e.snr, det_t0=e.t0,
                   shape=bool(v.shape_passes(sh["dbic_alt"], sh["dbic_box"])),
                   shape_strict_box=bool(v.shape_passes(sh["dbic_alt"], sh["dbic_box"], -6.0)),
                   dbic_alt=sh["dbic_alt"], dbic_box=sh["dbic_box"], best_alt=sh["best_alt"],
                   fit=tr, duration=bool(v.duration_passes(du)), edge_lc=bool(ed["passed"]))
        # 4: pixels, with the transit injected into the cutout
        cut = v.tesscut_cutout(r["ra"], r["dec"], SECTOR)
        cube = psf_inject(cut, t0, depth, dur)
        t14 = tr["t14_h"] / 24
        fdep = max(tr["depth_ppm"], 1) * 1e-6
        di = v.difference_image(cut["time"], cube, tr["t0"], t14)
        if di is None:
            res.update(pixels=None, edge=res["edge_lc"])
        else:
            diff, noise, oot, *_ = di
            nbs = v.tic_neighbours(r["ra"], r["dec"], cut["wcs"], r["Tmag"])
            ct = v.centroid_test(diff, noise, cut["x"], cut["y"])
            nt = v.neighbour_test(diff, noise, oot, cut["x"], cut["y"], nbs, fdep)
            ap = v.aperture_mask(oot, cut["x"], cut["y"])
            sel = np.abs(cut["time"] - tr["t0"]) < max(2.5 * t14, 0.5)
            tl = cut["time"][sel]
            fl = cube[sel][:, ap].sum(axis=1)
            fl = fl / np.median(fl)
            ap_depth, ap_snr = v.local_dip_snr(tl, fl, tr["t0"], t14)
            conf = v.pixel_confirm_check(ap_snr, e.snr, ap_depth, fdep)
            pe = v.edge_check(tl, fl, tr["t0"], tr["t14_h"], fdep)["passed"] and conf["passed"]
            res.update(centroid=bool(ct["passed"]), neighbours=bool(nt["passed"]),
                       pixel_confirm=bool(conf["passed"]), pixel_snr_ratio=conf["ratio"],
                       pixel_depth_ratio=conf["depth_ratio"],
                       pixels=bool(ct["passed"] and nt["passed"] and conf["passed"]),
                       edge_pixels=bool(pe), edge=bool(res["edge_lc"] or pe),
                       ap_abs=[[cut["col0"] + int(ix), cut["row0"] + int(iy)]
                               for iy, ix in zip(*np.nonzero(ap))])
        # 6: catalogues (star level)
        if _CATS is None:
            _CATS = v.load_catalogues(P3WORK)
        cat = v.catalogue_check(tic, _CATS)
        res["catalogue"] = bool(v.catalogue_check_passes(cat))
        jdump(out, res)
        return r["key"], "ok"
    except Exception as e:  # noqa: BLE001
        jdump(os.path.join(WORK, "errors", f"{r['key']}.json"),
              dict(error=f"{type(e).__name__}: {e}", tb=traceback.format_exc()))
        return r["key"], "error"


def passes_1_4(res):
    return bool(res.get("detected") and res.get("shape") and res.get("duration")
                and res.get("edge") and res.get("pixels"))


# ---------------------------------------------------------------- check 5

def sky(inj):
    n = 0
    for r in inj.to_dict("records"):
        p = os.path.join(WORK, "checks", f"{r['key']}.json")
        if not os.path.exists(p):
            continue
        res = jload(p)
        if not passes_1_4(res) or "asteroid" in res:
            continue
        objs = v.skybot_query(r["ra"], r["dec"], res["fit"]["t0"] + 2457000.0)
        res["asteroid"] = bool(v.asteroid_check(objs, r["ra"], r["dec"],
                                                res["fit"]["t14_h"])["passed"])
        jdump(p, res)
        n += 1
        if n % 50 == 0:
            print(f"[{time.strftime('%H:%M:%S')}] sky {n}", flush=True)
    print(f"sky: {n} queried")


# ---------------------------------------------------------------- check 7

def fpp_subset(inj, n=FPP_SUBSET, seed=7):
    rows = []
    for r in inj.to_dict("records"):
        p = os.path.join(WORK, "checks", f"{r['key']}.json")
        if os.path.exists(p):
            res = jload(p)
            if passes_1_4(res) and res.get("asteroid") and res.get("catalogue"):
                rows.append(dict(r, det_snr=res["det_snr"]))
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["bin"] = pd.cut(df.expected_snr, SNR_BINS)
    per = max(1, n // df["bin"].nunique())
    return (df.groupby("bin", observed=True, group_keys=False)
              .apply(lambda g: g.sample(min(len(g), per), random_state=seed)))


def fpp_one(r):
    p = os.path.join(WORK, "checks", f"{r['key']}.json")
    res = jload(p)
    if "fpp" in res:
        return r["key"], "cached"
    try:
        tic, t0 = int(r["tic"]), r["t0"]
        lc, _ = ffi.read(ffi.download(tic, SECTOR, cache=True))
        var, cfg, ss = analyse(inject(lc, t0, r["depth_ppm"] * 1e-6, r["duration_h"] / 24,
                                      INGRESS_FRAC))
        tr = res["fit"]
        tw, fw = v.shape_window(ss.time, ss.flat, tr["t0"], tr["t14_h"] / 24)
        err = v._robust_sigma_pt(fw)
        p_lo = max(tr["t0"] - ss.time[0], ss.time[-1] - tr["t0"])
        du = v.duration_check(tr["t14_h"], tr["b"], tr["rp"], r["rad"], r["mass"], tr["grazing"])
        p_hi = max(du["p_b0_d"] if np.isfinite(du["p_b0_d"]) else 3 * p_lo, 1.5 * p_lo)
        fr = v.triceratops_fpp(tic, r["ra"], r["dec"], SECTOR, tw - tr["t0"], fw, err,
                               max(tr["depth_ppm"], 1) * 1e-6, res["ap_abs"], (p_lo, p_hi),
                               os.path.join(WORK, "tri", r["key"]))
        res["fpp"] = bool(v.fpp_check(fr)["passed"])
        res["fpp_value"], res["nfpp_value"] = fr["fpp"], fr["nfpp"]
        jdump(p, res)
        return r["key"], "ok"
    except Exception as e:  # noqa: BLE001
        jdump(os.path.join(WORK, "errors", f"fpp_{r['key']}.json"),
              dict(error=f"{type(e).__name__}: {e}", tb=traceback.format_exc()))
        return r["key"], "error"


# ---------------------------------------------------------------- report

CHECK_ORDER = ["shape", "duration", "edge", "pixels", "asteroid", "catalogue", "fpp"]


def collect(inj):
    rows = []
    for r in inj.to_dict("records"):
        p = os.path.join(WORK, "checks", f"{r['key']}.json")
        res = jload(p) if os.path.exists(p) else {}
        o = {k: r[k] for k in ("key", "tic", "trial", "t0", "depth_ppm", "duration_h",
                               "expected_snr", "recovered", "plausible")}
        for k in ("detected", "det_snr", "dbic_alt", "dbic_box", "best_alt", "shape_strict_box",
                  "edge_lc", "edge_pixels", "centroid", "neighbours", "pixel_confirm",
                  "pixel_snr_ratio", "pixel_depth_ratio", "fpp_value", "nfpp_value") + tuple(CHECK_ORDER):
            o[k] = res.get(k)
        rows.append(o)
    return pd.DataFrame(rows)


def recall_table(df, snr_col="expected_snr"):
    """Per SNR bin: n injected, recovered, passing 1-6, and the 1-7 estimate
    (1-6 pass rate x FPP pass rate measured on the subset in that bin)."""
    out = []
    df = df.copy()
    df["bin"] = pd.cut(df[snr_col], SNR_BINS)
    for b, g in df.groupby("bin", observed=True):
        rec = g[g.recovered == 1]
        ok16 = rec[[c for c in CHECK_ORDER[:-1]]].eq(True).all(axis=1)
        sub = rec[ok16 & rec.fpp.notna()]
        fpp_rate = float(sub.fpp.eq(True).mean()) if len(sub) else np.nan
        out.append(dict(snr_bin=str(b), n=len(g), recovered=float((g.recovered == 1).mean()),
                        n_recovered=len(rec),
                        vet_1_6=float(ok16.mean()) if len(rec) else np.nan,
                        n_fpp=len(sub), fpp_pass=fpp_rate,
                        vet_1_7=float(ok16.mean() * fpp_rate) if len(rec) and len(sub) else np.nan,
                        **{f"pass_{c}": float(rec[c].eq(True).mean()) if len(rec) else np.nan
                           for c in CHECK_ORDER[:-1]}))
    return pd.DataFrame(out)


def report(inj):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(PLOTS, exist_ok=True)
    df = collect(inj)
    df.to_csv(os.path.join(OUT, "injection_vetting.csv"), index=False, float_format="%.6g")
    tabs = {"all": recall_table(df), "plausible": recall_table(df[df.plausible])}
    for k, t in tabs.items():
        t.to_csv(os.path.join(OUT, f"injection_vetting_recall_{k}.csv"), index=False,
                 float_format="%.3f")
    rec = df[df.recovered == 1]
    summ = dict(
        trials=len(df), recovered=int(len(rec)), plausible=int(df.plausible.sum()),
        recovered_plausible=int((rec.plausible).sum()),
        detected_again=int(rec.detected.eq(True).sum()),
        pass_by_check={c: dict(ran=int(rec[c].notna().sum()), passed=int(rec[c].eq(True).sum()))
                       for c in CHECK_ORDER},
        pass_by_check_plausible={c: dict(ran=int(rec[rec.plausible][c].notna().sum()),
                                         passed=int(rec[rec.plausible][c].eq(True).sum()))
                                 for c in CHECK_ORDER},
        shape_box_rule_extra_fail=int((rec["shape"].eq(True) & rec.shape_strict_box.eq(False)).sum()),
        fpp_subset=int(rec.fpp.notna().sum()),
        errors=len(os.listdir(os.path.join(WORK, "errors")))
        if os.path.exists(os.path.join(WORK, "errors")) else 0,
        recall_all=tabs["all"].to_dict("records"),
        recall_plausible=tabs["plausible"].to_dict("records"))
    with open(os.path.join(OUT, "injection_vetting_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1, default=float)

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
    for ax, (name, t) in zip(axes[:2], tabs.items()):
        x = np.arange(len(t))
        ax.plot(x, t.recovered, "o-", color="0.5", label="detected (Phase 2)")
        ax.plot(x, t.recovered * t.vet_1_6, "s-", color="C0", label="detected + checks 1-6")
        ax.plot(x, t.recovered * t.vet_1_7, "D-", color="C3",
                label="detected + all 7 (FPP from subset)")
        ax.set_xticks(x, [s.replace("(", "").replace("]", "").replace(", ", "-")
                          .replace("1000000000.0", "inf") for s in t.snr_bin], rotation=30)
        ax.set_xlabel("expected SNR")
        ax.set_ylabel("fraction of injected transits")
        ax.set_ylim(-0.03, 1.03)
        ax.grid(alpha=0.3)
        ax.set_title(f"{'all injections' if name == 'all' else 'physically plausible durations'}"
                     f" (n={int(t.n.sum())})", fontsize=10)
        ax.legend(fontsize=8)
    ax = axes[2]
    t = tabs["plausible"]
    x = np.arange(len(t))
    for c in CHECK_ORDER[:-1]:
        ax.plot(x, t[f"pass_{c}"], "o-", label=c)
    ax.plot(x, t.fpp_pass, "k--", label="fpp (subset)")
    ax.set_xticks(x, [s.replace("(", "").replace("]", "").replace(", ", "-")
                      .replace("1000000000.0", "inf") for s in t.snr_bin], rotation=30)
    ax.set_xlabel("expected SNR")
    ax.set_ylabel("pass rate among detected")
    ax.set_ylim(-0.03, 1.03)
    ax.grid(alpha=0.3)
    ax.set_title("per-check pass rate, plausible durations", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "injection_vetting_recall.png"), dpi=110)
    plt.close(fig)
    print(json.dumps({k: v_ for k, v_ in summ.items() if not k.startswith("recall")}, indent=1,
                     default=float))
    print(tabs["plausible"].round(2).to_string(index=False))


def _init():
    os.environ.setdefault("OMP_NUM_THREADS", "1")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["prep", "checks", "sky", "fpp", "report", "all"])
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--sector", type=int, default=48)
    args = ap.parse_args()
    configure(args.sector)
    inj = trials()
    todo = inj[inj.recovered == 1]
    if args.limit:
        todo = todo.head(args.limit)
    stages = ["prep", "checks", "sky", "fpp", "report"] if args.stage == "all" else [args.stage]
    for st in stages:
        if st == "prep":
            prep(todo)
        elif st == "checks":
            t_start, counts = time.time(), {}
            with Pool(args.procs, initializer=_init) as pool:
                for i, (k, s) in enumerate(pool.imap_unordered(check_trial, todo.to_dict("records")), 1):
                    counts[s] = counts.get(s, 0) + 1
                    if i % 100 == 0 or i == len(todo):
                        print(f"[{time.strftime('%H:%M:%S')}] checks {i}/{len(todo)} {counts} "
                              f"{i / (time.time() - t_start):.1f}/s", flush=True)
        elif st == "sky":
            sky(todo)
        elif st == "fpp":
            sub = fpp_subset(todo)
            print(f"fpp subset: {len(sub)}")
            with Pool(min(args.procs, 3), initializer=_init) as pool:
                for i, (k, s) in enumerate(pool.imap_unordered(fpp_one, sub.to_dict("records")), 1):
                    print(f"[{time.strftime('%H:%M:%S')}] fpp {i}/{len(sub)} {k}: {s}", flush=True)
        else:
            report(inj)


if __name__ == "__main__":
    main()
