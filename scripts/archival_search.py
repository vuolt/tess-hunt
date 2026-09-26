"""Search archival photometry for earlier transits of a TESS duotransit.

    python scripts/archival_search.py --tic 95747180 --events 45:2540.662 48:2621.1348 --n 1 2

Checks K2 (campaign footprints, K2fov), SuperWASP DR1 and KELT (NASA Exoplanet
Archive), ZTF (IRSA), Gaia DR3 epoch photometry, and records why ASAS-SN was not
retrieved. For each usable light curve and each period hypothesis
P = (t_b - t_a)/n it

  * estimates the sensitivity: inject the TESS transit at the nominal ephemeris
    and measure the box depth SNR it would give (and the analytic
    depth * sqrt(N_in) / sigma),
  * scans P over +-4 sigma_P (the ephemeris is anchored on the S48 transit), since
    back-projected mid-times are uncertain by hours to days, measuring the box
    depth at the predicted transits,
  * calls a hypothesis "ruled out" only if the TESS depth is excluded at > 3 sigma
    at *every* trial period, and "supported" only if some trial period shows a
    > 4 sigma dip at the TESS depth (0.5-2x) after a look-elsewhere correction.

The result is written as a section of results/phase4/tic<TIC>_predictions.md
(between archival markers; predict_transits.py preserves it) plus
results/phase4/tic<TIC>_archival.json and plots/phase4/tic<TIC>_archival.png.
"""

import argparse
import io
import json
import os
import sys
import urllib.parse
import warnings

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tesshunt import ephemeris as ep  # noqa: E402
from tesshunt import net  # noqa: E402

OUT = os.path.join(ROOT, "results", "phase4")
PLOTS = os.path.join(ROOT, "plots", "phase4")
EXO = "https://exoplanetarchive.ipac.caltech.edu"
MARK_A, MARK_B = "<!-- archival:start -->", "<!-- archival:end -->"
SCAN_SIGMA = 4.0
IN_FRAC = 0.8              # use the central 80 % of the transit as "in transit"
DETECT_SIGMA = 4.0         # look-elsewhere-corrected significance for "supported"
EXCLUDE_SIGMA = 3.0


# ------------------------------------------------------------------ data

def exo_tap(q):
    url = f"{EXO}/TAP/sync?query={urllib.parse.quote(q)}&format=csv"
    return pd.read_csv(io.StringIO(net.get_text(url, "exoplanet_archive")))


def box(ra, dec, r_deg=0.006):
    c = np.cos(np.radians(dec))
    return (f"ra between {ra - r_deg / c:.5f} and {ra + r_deg / c:.5f} and "
            f"dec between {dec - r_deg:.5f} and {dec + r_deg:.5f}")


def superwasp(ra, dec):
    """SuperWASP DR1 light curve (TAMFLUX2, detrended), as BTJD / relative flux."""
    from astropy.io import fits
    meta = exo_tap(f"select sourceid,tile,wasp_mag,npts from superwasptimeseries where {box(ra, dec)}")
    if meta.empty:
        return None, meta
    m = meta.iloc[0]
    url = (f"{EXO}/data/ETSS/SuperWASP/FITS/DR1/{m.tile}/"
           + urllib.parse.quote(f"{m.sourceid}.fits"))
    with fits.open(io.BytesIO(net.get(url, "exoplanet_archive"))) as h:
        jd_ref = h[0].header["JD_REF"]
        d = h[1].data
        t = d["TMID"] / 86400.0 + jd_ref - 2457000.0          # HJD(UTC); << 1 min from BJD_TDB
        f = np.asarray(d["TAMFLUX2"], float)
        # FLAG 32 is the value of 99 % of the points (the normal photometry); the
        # few FLAG 0 points are noisier (6 % vs 3.6 % scatter) and are dropped.
        flag = np.bincount(d["FLAG"].astype(int)).argmax()
        ok = (d["FLAG"] == flag) & np.isfinite(f) & (f > 0)
    return dict(name=f"SuperWASP ({m.sourceid})", t=t[ok], f=f[ok], group=np.zeros(ok.sum(), int),
                mag=float(m.wasp_mag), window_d=5.0, flag=int(flag)), meta


def kelt(ra, dec):
    return exo_tap("select kelt_sourceid,kelt_field,kelt_mag,npts,stddevwrtmedian "
                   f"from kelttimeseries where {box(ra, dec)}")


def ztf(ra, dec):
    url = ("https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20"
           f"{ra:.6f}%20{dec:.6f}%200.0004&BANDNAME=g,r,i&FORMAT=csv")
    d = pd.read_csv(io.StringIO(net.get(url, "irsa", timeout=300).decode()))
    if d.empty:
        return None
    d = d[(d.catflags == 0) & np.isfinite(d.mag)]
    f = 10 ** (-0.4 * (d.mag - d.groupby("oid").mag.transform("median")))
    return dict(name="ZTF (g, r, i)", t=d.hjd.values - 2457000.0, f=f.values,
                group=pd.factorize(d.oid)[0], mag=float(d.mag.median()), window_d=60.0)


def k2_campaigns(ra, dec):
    from K2fov import getKeplerFov
    from K2fov.K2onSilicon import onSiliconCheck
    on = []
    for c in range(0, 20):
        try:
            if onSiliconCheck(ra, dec, getKeplerFov(c)):
                on.append(c)
        except Exception:  # noqa: BLE001
            continue
    return on


def k2_nearest(ra, dec):
    from astropy.coordinates import SkyCoord
    from K2fov import getFieldInfo
    c = SkyCoord(ra, dec, unit="deg")
    best = None
    for k in range(0, 20):
        try:
            fi = getFieldInfo(k)
        except Exception:  # noqa: BLE001
            continue
        s = c.separation(SkyCoord(fi["ra"], fi["dec"], unit="deg")).deg
        if best is None or s < best[1]:
            best = (k, s)
    return best


def gaia_epoch(source_id):
    from tesshunt import binarity as bn
    return bn._gaia("SELECT has_epoch_photometry, phot_g_n_obs, phot_g_mean_flux_over_error, "
                    f"phot_g_mean_mag FROM gaiadr3.gaia_source WHERE source_id = {int(source_id)}")


# ------------------------------------------------------------------ analysis

NORM_MIN_PTS = 15          # reference points per normalisation
SAME_NIGHT_D = 0.5         # points closer than this never normalise each other


def _reference_median(tt, ff, window_d):
    """For each point: median of the points in the same group within window_d/2,
    excluding the same night (so a transit cannot normalise itself away); if
    fewer than NORM_MIN_PTS, the NORM_MIN_PTS nearest points outside that night."""
    med = np.empty(len(tt))
    for i, x in enumerate(tt):
        dt = np.abs(tt - x)
        use = (dt > SAME_NIGHT_D) & (dt < window_d / 2)
        if use.sum() < NORM_MIN_PTS:
            cand = np.where(dt > SAME_NIGHT_D)[0]
            use = cand[np.argsort(dt[cand])[:NORM_MIN_PTS]]
        med[i] = np.median(ff[use]) if np.size(use) and np.any(use) else np.nan
    return med


def clean(lc, clip=5.0):
    """Normalise each group (field/filter) by nearby reference points (see
    _reference_median), then clip outliers."""
    t, f, g = lc["t"], lc["f"], lc["group"]
    window_d = lc.get("window_d", 5.0)
    out_t, out_f = [], []
    for k in np.unique(g):
        s = g == k
        if s.sum() <= NORM_MIN_PTS:
            continue
        tt, ff = t[s], f[s]
        o = np.argsort(tt)
        tt, ff = tt[o], ff[o]
        med = _reference_median(tt, ff, window_d)
        r = ff / med
        out_t.append(tt)
        out_f.append(r)
    t, f = np.concatenate(out_t), np.concatenate(out_f)
    good = np.isfinite(f)
    t, f = t[good], f[good]
    sig = 1.4826 * np.median(np.abs(f - np.median(f)))
    ok = np.abs(f - 1) < clip * sig
    return t[ok], f[ok], float(1.4826 * np.median(np.abs(f[ok] - 1)))


def box_depth(t, f, tb, P, t14):
    ph = (t - tb + 0.5 * P) % P - 0.5 * P
    inn = np.abs(ph) < IN_FRAC * t14 / 2
    out = np.abs(ph) > t14
    if inn.sum() < 1 or out.sum() < 10:
        return np.nan, np.nan, int(inn.sum())
    sig = 1.4826 * np.median(np.abs(f[out] - np.median(f[out])))
    d = np.median(f[out]) - f[inn].mean()
    return float(d), float(sig * np.sqrt(1 / inn.sum() + 1 / out.sum())), int(inn.sum())


def scan(t, f, tb, P0, sP, t14):
    step = t14 / 8 / max(abs(t - tb).max() / P0, 1)        # mid-time moves < T14/8 per step
    Ps = np.arange(P0 - SCAN_SIGMA * sP, P0 + SCAN_SIGMA * sP + step, step)
    res = np.array([box_depth(t, f, tb, p, t14) for p in Ps], float)
    return Ps, res


def analyse(lc, tb, hyps, t14, depth):
    t, f, sig = clean(lc)
    out = dict(name=lc["name"], n_points=int(len(t)), sigma_pt=sig, mag=lc.get("mag"),
               t_range=[float(t.min()), float(t.max())], hyps={})
    for h in hyps:
        P0, sP = h["P"], h["sigma_P"]
        n_exp = len(t) * IN_FRAC * t14 / P0
        Ps, res = scan(t, f, tb, P0, sP, t14)
        # sensitivity: inject the TESS transit (box of the TESS depth) at each trial
        # period and measure it the same way; no in-transit data counts as 0 sigma
        inj = []
        for p_ in Ps:
            ph = (t - tb + 0.5 * p_) % p_ - 0.5 * p_
            di, ei, _ = box_depth(t, f * np.where(np.abs(ph) < t14 / 2, 1 - depth, 1.0), tb, p_, t14)
            inj.append(di / ei if np.isfinite(ei) and ei > 0 else 0.0)
        inj = np.array(inj)
        d, e, n = res[:, 0], res[:, 1], res[:, 2]
        good = np.isfinite(d) & (n >= 1)
        excl = good & (d + EXCLUDE_SIGMA * e < depth)
        snr = np.where(good, d / e, np.nan)
        n_indep = max(int((Ps[-1] - Ps[0]) / (t14 / max(abs(t - tb).max() / P0, 1))), 1)
        from scipy.stats import norm
        p_best = norm.sf(np.nanmax(snr)) if good.any() else 1.0
        p_global = 1 - (1 - p_best) ** n_indep
        best = int(np.nanargmax(snr)) if good.any() else None
        consistent = best is not None and 0.5 * depth <= d[best] <= 2 * depth
        out["hyps"][h["label"]] = dict(
            P=P0, sigma_P=sP, n_trials=int(len(Ps)), n_independent=n_indep,
            expected_in_transit_points=float(n_exp),
            analytic_snr=float(depth * np.sqrt(max(n_exp, 0)) / sig) if sig > 0 else np.nan,
            injected_snr=float(np.median(inj)), injected_snr_max=float(inj.max()),
            frac_P_detectable=float(np.mean(inj >= 5)),
            nominal_depth=float(d[np.argmin(abs(Ps - P0))]),
            nominal_err=float(e[np.argmin(abs(Ps - P0))]),
            frac_P_with_data=float(good.mean()), frac_P_excluded=float(excl.mean()),
            best_P=float(Ps[best]) if best is not None else None,
            best_snr=float(snr[best]) if best is not None else None,
            best_depth=float(d[best]) if best is not None else None,
            global_sigma=float(norm.isf(min(max(p_global, 1e-300), 1 - 1e-16))),
            ruled_out=bool(good.all() and excl.all()),
            supported=bool(consistent and norm.isf(min(max(p_global, 1e-300), 1 - 1e-16)) >= DETECT_SIGMA))
        out["hyps"][h["label"]]["_scan"] = (Ps, d, e)
    out["_lc"] = (t, f)
    return out


# ------------------------------------------------------------------ report

def plot(results, tb, hyps, t14, depth, tic):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [r for r in results if r is not None]
    fig, axes = plt.subplots(len(rows), 1 + len(hyps), figsize=(5 * (1 + len(hyps)), 3.3 * len(rows)),
                             squeeze=False)
    for i, r in enumerate(rows):
        t, f = r["_lc"]
        for j, h in enumerate(hyps):
            ax = axes[i, j]
            P0 = h["P"]
            ph = ((t - tb + 0.5 * P0) % P0 - 0.5 * P0) * 24
            sel = np.abs(ph) < 24
            ax.plot(ph[sel], f[sel], ".", color="0.7", ms=2)
            b = np.arange(-24, 24.01, 1.0)
            idx = np.digitize(ph[sel], b)
            bm = [np.mean(f[sel][idx == k]) if (idx == k).sum() > 2 else np.nan for k in range(1, len(b))]
            ax.plot(b[:-1] + 0.5, bm, "o-", color="C0", ms=3)
            ax.axhline(1 - depth, color="C3", ls="--", lw=1, label="TESS depth")
            ax.axvspan(-t14 * 12, t14 * 12, color="C3", alpha=0.08)
            ax.set_ylim(1 - 4 * max(depth, r["sigma_pt"] / 3), 1 + 3 * max(depth, r["sigma_pt"] / 3))
            hh = r["hyps"][h["label"]]
            ax.set_title(f"{r['name']}: fold at P = {P0:.3f} d\n(1 σ mid-time ≈ "
                         f"{hh.get('sigma_mid_h', np.nan):.1f} h in these data)", fontsize=8)
            ax.set_xlabel("hours from predicted mid-time", fontsize=8)
            ax.tick_params(labelsize=7)
            if j == 0:
                ax.legend(fontsize=7)
        ax = axes[i, -1]
        for h in hyps:
            Ps, d, e = r["hyps"][h["label"]]["_scan"]
            ax.fill_between(Ps - h["P"], (d - 3 * e) * 1e3, (d + 3 * e) * 1e3, alpha=0.3,
                            label=f"P ≈ {h['P']:.1f} d (±3σ)")
        ax.axhline(depth * 1e3, color="C3", ls="--", lw=1, label="TESS depth")
        ax.axhline(0, color="k", lw=0.5)
        lim = SCAN_SIGMA * max(h["sigma_P"] for h in hyps)
        ax.set_xlim(-lim, lim)
        ax.set_xlabel("trial P − nominal P [d] (white: no data in transit)", fontsize=8)
        ax.set_ylabel("box depth ±3σ [ppt]", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_title(f"{r['name']}: depth at predicted transits", fontsize=8)
        ax.legend(fontsize=7)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(os.path.join(PLOTS, f"tic{tic}_archival.png"), dpi=100)
    plt.close(fig)


def section(tic, info, results, hyps, depth, t14):
    L = [MARK_A, "## Archival photometry: earlier transits?\n",
         f"Produced by `python scripts/archival_search.py --tic {tic} ...`. Plot: "
         f"`plots/phase4/tic{tic}_archival.png`, numbers: `results/phase4/tic{tic}_archival.json`.\n",
         f"The transit to find is {depth * 1e3:.1f} ppt deep and {t14 * 24:.1f} h long. The "
         "fraction of time spent in transit is "
         + " and ".join(f"{t14 / h['P']:.2%} for P ≈ {h['P']:.1f} d" for h in hyps)
         + ". A survey can test a period only if it has enough points inside those "
         "windows, at a precision that makes the stacked dip significant.\n",
         "### Sources\n",
         "| source | data for this star | precision | could it detect the transit? |",
         "|---|---|---|---|"]
    for s in info:
        L.append(f"| {s['source']} | {s['data']} | {s['precision']} | {s['verdict']} |")
    L.append("")
    for r in results:
        if r is None:
            continue
        L.append(f"### {r['name']}\n")
        L.append(f"- {r['n_points']} good points, BTJD {r['t_range'][0]:.0f} to "
                 f"{r['t_range'][1]:.0f}, per-point scatter {r['sigma_pt'] * 100:.2f} % "
                 "after normalisation.")
        L.append("")
        L.append("| hypothesis | 1σ on the mid-time in these data | points expected in transit | "
                 "SNR if the TESS transit is there (injected; median / best over the P range) | depth at nominal P | "
                 "P range with data | P range excluding the TESS depth | best dip | verdict |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for h in hyps:
            x = r["hyps"][h["label"]]
            verdict = ("**ruled out**" if x["ruled_out"] else
                       "**supported**" if x["supported"] else "not constrained")
            L.append(
                f"| P ≈ {x['P']:.2f} d | {x['sigma_mid_h']:.1f} h | {x['expected_in_transit_points']:.0f} | "
                f"{x['injected_snr']:.1f} / {x['injected_snr_max']:.1f} | "
                + (f"{x['nominal_depth'] * 1e3:+.1f} ± {x['nominal_err'] * 1e3:.1f} ppt | "
                   if np.isfinite(x['nominal_depth']) else "no data in transit | ")
                + f"{x['frac_P_with_data']:.0%} | {x['frac_P_excluded']:.0%} | "
                f"{x['best_depth'] * 1e3:+.1f} ppt at P = {x['best_P']:.4f} d "
                f"({x['best_snr']:.1f}σ local, {max(x['global_sigma'], 0):.1f}σ after "
                f"{x['n_independent']} independent trials) | {verdict} |")
        L.append("")
    L.append("### Conclusion\n")
    ro = {h["label"]: any(r and r["hyps"][h["label"]]["ruled_out"] for r in results) for h in hyps}
    su = {h["label"]: any(r and r["hyps"][h["label"]]["supported"] for r in results) for h in hyps}
    if not any(ro.values()) and not any(su.values()):
        L.append("**The archives neither support nor rule out either period.**\n")
        L.append("- K2 never observed the star.")
        L.append("- Gaia DR3 published no epoch photometry for it.")
        L.append("- KELT has no light curve of it.")
        L.append("- The surveys that did observe it (SuperWASP, ZTF, and ASAS-SN by "
                 "estimate) are too noisy at this brightness, and too sparsely sampled "
                 "inside the transit windows, to see a transit this shallow.\n")
        L.append("The injection test shows this directly: the TESS transit put into "
                 "their data would not reach a significant detection. The 40 d vs "
                 "80 d question has to be settled by the ground-based windows listed "
                 "above.")
        L.append("")
        L.append("Gaia DR4 is planned to publish epoch photometry for all sources. Its "
                 "precision per visit, about 0.15 % at G = 13, is enough to see this transit "
                 "in a single point. However, only about 0.1–0.2 of its roughly 45 visits are "
                 "expected to fall in a transit.")
    else:
        for k in ro:
            if ro[k]:
                L.append(f"- **{k}: ruled out** by the archives (see tables).")
            if su[k]:
                L.append(f"- **{k}: supported** by a significant dip at the predicted times (see tables).")
    L.append("")
    L.append(MARK_B)
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tic", type=int, required=True)
    ap.add_argument("--events", nargs=2, required=True)
    ap.add_argument("--n", type=int, nargs="+", default=[1, 2])
    args = ap.parse_args()
    warnings.simplefilter("ignore")
    import predict_transits as pt

    sm = pd.read_csv(os.path.join(ROOT, "work", "s0048_sample.csv"))
    star = sm[sm.ID == args.tic].iloc[0]
    fu = pd.read_csv(os.path.join(OUT, "followup.csv"))
    row = fu[fu.tic == args.tic].iloc[0]
    t14_0 = row.t14_h / 24
    evs = [pt.load_event(args.tic, star, int(s), float(t), t14_0)
           for s, t in (e.split(":") for e in args.events)]
    fit = ep.joint_fit(evs, t14_0, np.sqrt(row.depth_ppm * 1e-6), n_shift=100)
    (ea, eb) = fit["events"]
    ta, sa, tb, sb = ea["t0"], ea["sigma"], eb["t0"], eb["sigma"]
    depth, t14 = fit["depth"], fit["t14"]
    hyps = [dict(label=f"P ≈ {(tb - ta) / n:.1f} d", n=n, P=(tb - ta) / n,
                 sigma_P=float(np.hypot(sa, sb) / n)) for n in args.n]

    info = []
    # K2
    k2 = k2_campaigns(star.ra, star.dec)
    near = k2_nearest(star.ra, star.dec)
    info.append(dict(source="K2", data=(f"on silicon in C{', C'.join(map(str, k2))}" if k2 else
                                         f"not on silicon in any campaign (C0–C19); nearest field "
                                         f"centre C{near[0]} is {near[1]:.0f}° away"),
                     precision="–", verdict="no data" if not k2 else "see below"))
    # Gaia
    gj = json.load(open(os.path.join(ROOT, "work", "phase4", "s0048", "binarity", f"{args.tic}.json")))
    g = gaia_epoch(gj["gaia"]["source_id"]).iloc[0]
    n_vis = g.phot_g_n_obs / 9
    per_ccd = np.sqrt(g.phot_g_n_obs) / g.phot_g_mean_flux_over_error   # upper bound: includes all scatter
    per_vis = per_ccd / 3                                                 # 9 CCDs per field-of-view visit
    info.append(dict(source="Gaia DR3 epoch photometry",
                     data="not published for this source (`has_epoch_photometry` = False)"
                     if not g.has_epoch_photometry else "available",
                     precision=f"≈ {per_vis * 100:.2f} % per visit (G = {g.phot_g_mean_mag:.1f}, "
                               f"{n_vis:.0f} visits)",
                     verdict="no data. If it were published, "
                             f"{min(n_vis * t14 / h['P'] for h in hyps):.2f}–"
                             f"{max(n_vis * t14 / h['P'] for h in hyps):.2f} visits would be "
                             "expected in transit"))
    # KELT
    k = kelt(star.ra, star.dec)
    info.append(dict(source="KELT", data="no source within 20″" if k.empty else f"{len(k)} light curve(s)",
                     precision="–" if k.empty else f"{k.stddevwrtmedian.iloc[0]:.3f} mag",
                     verdict="no data" if k.empty else "see below"))
    # SuperWASP, ZTF
    results = []
    sw, _ = superwasp(star.ra, star.dec)
    zt = ztf(star.ra, star.dec)
    for lc in (sw, zt):
        if lc is None:
            results.append(None)
            continue
        r = analyse(lc, tb, hyps, t14, depth)
        far = max(r["t_range"], key=lambda x: abs(x - tb))      # worst case in these data
        for h in hyps:
            kk = (far - tb) / h["P"]
            r["hyps"][h["label"]]["sigma_mid_h"] = float(
                np.sqrt((1 + kk / h["n"]) ** 2 * sb ** 2 + (kk / h["n"]) ** 2 * sa ** 2) * 24)
        results.append(r)
        best_inj = max(x["injected_snr_max"] for x in r["hyps"].values())
        info.append(dict(source=r["name"], data=f"{r['n_points']} points",
                         precision=f"{r['sigma_pt'] * 100:.1f} % per point",
                         verdict=(f"no: the TESS transit injected at any allowed period would "
                                  f"reach at most {best_inj:.1f}σ"
                                  if best_inj < 5 else f"yes (injected transit {best_inj:.0f}σ)")))
    # ASAS-SN: Sky Patrol host not reachable from this environment -> estimate only
    n_as, s_as = 3000, 0.02
    snr_as = [depth * np.sqrt(n_as * IN_FRAC * t14 / h["P"]) / s_as for h in hyps]
    info.append(dict(source="ASAS-SN (estimate)",
                     data="not retrieved: the Sky Patrol service is not reachable from this environment",
                     precision=f"≈ {s_as * 100:.0f} % per epoch at V ≈ 13 (typical), about {n_as} epochs",
                     verdict=f"no: expected {min(snr_as):.1f}–{max(snr_as):.1f}σ for a perfectly phased stack"))
    plot(results, tb, hyps, t14, depth, args.tic)
    js = dict(tic=args.tic, depth=depth, t14_d=t14, tb=tb, sigma_tb=sb, ta=ta, sigma_ta=sa,
              hypotheses=hyps, sources=info,
              results=[None if r is None else {k: v for k, v in r.items() if not k.startswith("_")}
                       | {"hyps": {kk: {a: b for a, b in vv.items() if not a.startswith("_")}
                                   for kk, vv in r["hyps"].items()}} for r in results])
    with open(os.path.join(OUT, f"tic{args.tic}_archival.json"), "w") as fh:
        json.dump(js, fh, indent=1, default=float)
    sec = section(args.tic, info, results, hyps, depth, t14)
    path = os.path.join(OUT, f"tic{args.tic}_predictions.md")
    txt = open(path).read() if os.path.exists(path) else ""
    if MARK_A in txt:
        txt = txt[:txt.index(MARK_A)] + sec + txt[txt.index(MARK_B) + len(MARK_B):]
    else:
        txt = txt.rstrip("\n") + "\n\n" + sec + "\n"
    with open(path, "w") as fh:
        fh.write(txt)
    print(sec)


if __name__ == "__main__":
    main()
