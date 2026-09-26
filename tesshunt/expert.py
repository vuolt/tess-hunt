"""Phase 5: expert-level checks on a single-transit candidate.

  aperture_test      depth in 1 / 1.5 / 2 / 3 px TESScut apertures vs the depth a
                     PSF model predicts for a signal on the target (and on each
                     TIC neighbour); flags aperture-dependent depths
  gp_reprocess       joint celerite2 GP + transit fit on the undetrended light
                     curve; depth and SNR (from the likelihood gain) vs Phase 3
  gaia_star          one Gaia DR3 query: astrometry, photometry, GSP-Phot/FLAME
                     astrophysics, IPD/RV binarity indicators, variability
  stellar_check      dwarf / subgiant / giant from FLAME stage, radius, log g
  binarity_extras    ipd_frac_multi_peak, ipd_gof_harmonic_amplitude + RUWE,
                     RV chi^2 / renormalised GOF (Katz et al. 2023)
  variability_check  Gaia DR3 variability flag/classes + AAVSO VSX (VizieR)
  density_fit        emcee transit fit with a stellar-density prior: the
                     circular-orbit period the shape implies, compared with the
                     periods TESS still allows (or, for a known period, the
                     density the shape needs vs the star's)
"""

from __future__ import annotations

import numpy as np

from . import net
from .vetting import (_robust_sigma_pt, local_dip_snr, tesscut_cutout,
                      tic_neighbours)

G_CGS = 6.674e-8
RSUN_CM, MSUN_G = 6.957e10, 1.989e33
RHO_SUN = MSUN_G / (4 / 3 * np.pi * RSUN_CM ** 3)       # 1.41 g/cm^3
RSUN_RJUP = 9.731
LD_U = (0.40, 0.25)

APERTURE_RADII = (1.0, 1.5, 2.0, 3.0)
APERTURE_SIGMA = 3.0          # an aperture "deviates" beyond this ...
APERTURE_SINGLE_SIGMA = 5.0   # ... flag if >= 2 deviate, or one by this much
NULL_EPOCHS = 40
GP_SNR_MIN = 7.0              # flag if the GP likelihood gain gives less than this ...
GP_DEPTH_MIN = 0.5            # ... or the GP depth is below this fraction of Phase 3's
E_MAX = 0.5                     # eccentricity allowed when judging physical consistency


# ------------------------------------------------------------------ 1. apertures

def _psf_image(shape, stars, sigma, halo_frac=0.0, halo_scale=3.0):
    """Core + halo double-Gaussian PSF (TESS PSFs have broad wings) summed over
    stars = [(x, y, flux)]; returns (total, per-star list). Sampled at pixel
    centres, 3x3 sub-sampled so a sub-pixel core is integrated properly."""
    ny, nx = shape
    sub = (np.arange(3) - 1) / 3.0
    ims = []
    for x, y, f in stars:
        im = np.zeros(shape)
        for dy in sub:
            for dx in sub:
                yy, xx = np.mgrid[0:ny, 0:nx]
                r2 = (xx + dx - x) ** 2 + (yy + dy - y) ** 2
                core = np.exp(-r2 / (2 * sigma ** 2)) / (2 * np.pi * sigma ** 2)
                hs = sigma * halo_scale
                halo = np.exp(-r2 / (2 * hs ** 2)) / (2 * np.pi * hs ** 2)
                im += (1 - halo_frac) * core + halo_frac * halo
        ims.append(f * im / 9.0)
    return np.sum(ims, axis=0), ims


def _fit_psf(oot, x0, y0, nbs):
    """Fit target flux, PSF width and a flat background to the out-of-transit
    image (neighbour fluxes fixed relative to the target by their TIC dTmag)."""
    from scipy.optimize import least_squares
    ny, nx = oot.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    use = np.hypot(xx - x0, yy - y0) <= 4.5
    rel = [(x, y, 10 ** (-0.4 * dm)) for _, x, y, dm in nbs if dm < 6]
    dist = [float(np.hypot(x - x0, y - y0)) for _, x, y, dm in nbs if dm < 6]

    def model(p):
        a, s, hf, hsc, bg = p
        tot, _ = _psf_image(oot.shape, [(x0, y0, a)] + [(x, y, a * r) for x, y, r in rel], s, hf, hsc)
        return tot + bg

    # weight by Poisson-like noise so the bright core does not dominate the wings
    w = 1 / np.sqrt(np.abs(oot) + np.median(np.abs(oot)))
    a0 = oot[int(round(y0)), int(round(x0))] * 2 * np.pi * 0.8 ** 2
    r = least_squares(lambda p: ((model(p) - oot) * w)[use], [max(a0, 1.0), 0.7, 0.2, 3.0, 0.0],
                      bounds=([0, 0.3, 0.0, 2.0, -np.inf], [np.inf, 2.0, 0.5, 4.0, np.inf]))
    a, s, hf, hsc, bg = r.x
    return dict(flux=float(a), sigma=float(s), halo_frac=float(hf), halo_scale=float(hsc), bg=float(bg),
                stars=[(x0, y0, a)] + [(x, y, a * rr) for x, y, rr in rel],
                ids=["target"] + [str(i) for i, x, y, dm in nbs if dm < 6],
                dist=[0.0] + dist)


def _keep_mask(t, other_transits):
    """False within other known transits on the star: (centre, duration) pairs,
    padded by half a duration plus 30 min on each side."""
    keep = np.ones(len(t), bool)
    for tc, dur in other_transits or []:
        keep &= np.abs(t - tc) > dur + 1 / 48
    return keep


def aperture_test(ra, dec, tmag, sector, t0, t14, radii=APERTURE_RADII, other_transits=()):
    cut = tesscut_cutout(ra, dec, sector)
    k = _keep_mask(cut["time"], other_transits)
    cut = dict(cut, time=cut["time"][k], cube=cut["cube"][k])
    x0, y0 = cut["x"], cut["y"]
    t = cut["time"]
    oot_sel = np.abs(t - t0) > max(t14, 0.25)
    oot = np.median(cut["cube"][oot_sel], axis=0)
    nbs = tic_neighbours(ra, dec, cut["wcs"], tmag)
    psf = _fit_psf(oot, x0, y0, nbs)
    _, ims = _psf_image(oot.shape, psf["stars"], psf["sigma"], psf["halo_frac"], psf["halo_scale"])
    ny, nx = oot.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = np.hypot(xx - x0, yy - y0)
    half = max(2.5 * t14, 0.5)
    sel = np.abs(t - t0) < half
    # Empirical depth errors: the same estimator at random epochs away from the
    # dip captures red noise and pointing jitter, which matter most in small apertures.
    rng = np.random.default_rng(0)
    cand = t[(np.abs(t - t0) > 2 * t14 + half) & (t > t.min() + half) & (t < t.max() - half)]
    null_t = rng.choice(cand, size=min(NULL_EPOCHS, len(cand)), replace=False) if len(cand) else []
    rows = []
    for rad in radii:
        ap = r <= rad
        ap[int(round(y0)), int(round(x0))] = True
        full = cut["cube"][:, ap].sum(axis=1)
        full = full / np.median(full)
        d, snr = local_dip_snr(t[sel], full[sel], t0, t14)
        err_formal = abs(d / snr) if snr and np.isfinite(snr) and snr != 0 else np.nan
        nulls = []
        for tn in null_t:
            w = np.abs(t - tn) < half
            dn, _ = local_dip_snr(t[w], full[w], tn, t14)
            if np.isfinite(dn):
                nulls.append(dn)
        err_emp = float(1.4826 * np.median(np.abs(np.array(nulls) - np.median(nulls)))) \
            if len(nulls) >= 10 else np.nan
        err = float(np.nanmax([err_formal, err_emp]))
        tot = sum(im[ap].sum() for im in ims) + psf["bg"] * ap.sum()
        frac = [im[ap].sum() / tot for im in ims]          # share of aperture flux per star
        rows.append(dict(radius=rad, npix=int(ap.sum()), depth=float(d), err=float(err),
                         err_formal=float(err_formal), err_empirical=err_emp, n_null=len(nulls),
                         snr=float(snr), frac=frac))
    depth = np.array([x["depth"] for x in rows])
    err = np.array([x["err"] for x in rows])
    ok = np.isfinite(depth) & np.isfinite(err) & (err > 0)
    hyp = []
    for j, name in enumerate(psf["ids"]):
        c = np.array([x["frac"][j] for x in rows])
        w = 1 / err[ok] ** 2
        amp = float((w * c[ok] * depth[ok]).sum() / (w * c[ok] ** 2).sum()) if ok.any() else np.nan
        chi2 = float((((depth[ok] - amp * c[ok]) / err[ok]) ** 2).sum()) if ok.any() else np.nan
        hyp.append(dict(source=name, sep_px=psf["dist"][j], intrinsic_depth=amp, chi2=chi2,
                        physical=bool(amp <= 1.0), expected=(amp * c).tolist()))
    tgt = hyp[0]
    resid = (depth - np.array(tgt["expected"])) / err
    worst = float(np.nanmax(np.abs(resid[ok]))) if ok.any() else np.nan
    n_dev = int(np.sum(np.abs(resid[ok]) > APERTURE_SIGMA))
    # neighbours closer than 1 px cannot be told apart from the target by aperture size
    alt = [h for h in hyp[1:] if h["physical"] and np.isfinite(h["chi2"]) and h["sep_px"] >= 1.0]
    unresolved = [h["source"] for h in hyp[1:] if h["sep_px"] < 1.0]
    best_alt = min(alt, key=lambda h: h["chi2"]) if alt else None
    flag = bool(ok.sum() >= 3 and (n_dev >= 2 or worst >= APERTURE_SINGLE_SIGMA))
    prefers_other = bool(best_alt is not None and best_alt["chi2"] + 9 < tgt["chi2"])
    return dict(rows=rows, psf_sigma_px=psf["sigma"], psf_halo_frac=psf["halo_frac"],
                psf_halo_scale=psf["halo_scale"], hypotheses=hyp,
                worst_resid_sigma=worst, n_deviating=n_dev, unresolved=unresolved,
                chi2_target=tgt["chi2"], n_ap=int(ok.sum()),
                best_other=best_alt, flag=flag or prefers_other,
                prefers_other=prefers_other)


# ------------------------------------------------------------------ 2. GP

def transit_lc(t, t0, per, rp, a, inc, exp_time=None, u=LD_U):
    import batman
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a, p.inc = t0, per, rp, a, inc
    p.ecc, p.w, p.u, p.limb_dark = 0.0, 90.0, list(u), "quadratic"
    kw = dict(supersample_factor=7, exp_time=exp_time) if exp_time else {}
    return batman.TransitModel(p, t, **kw).light_curve(p)


def _shape_model(t, t0, rp, t14, b, exp_time, per=100.0):
    b = min(b, (1 + rp) * 0.999)
    s = np.sin(np.pi * t14 / per)
    a = np.sqrt(b ** 2 + ((1 + rp) ** 2 - b ** 2) / s ** 2)
    return transit_lc(t, t0, per, rp, a, np.degrees(np.arccos(b / a)), exp_time)


def gp_reprocess(t, f, t0, t14, depth0, exp_time, window_d=4.0, other_transits=()):
    """Joint GP (celerite2 SHO, Q = 1/sqrt 2) + transit fit, and the same GP
    without a transit. The GP timescale is kept >= 3 T14 so it cannot absorb
    the transit. Returns depth, SNR = sqrt(2 dlnL), and the GP-detrended flux."""
    import celerite2
    from celerite2 import terms
    from scipy.optimize import minimize
    sel = (np.abs(t - t0) < window_d) & _keep_mask(t, other_transits)
    t, f = t[sel], f[sel] / np.median(f[sel])
    sig = _robust_sigma_pt(f)
    rho_min = 3 * t14

    def gp_for(p):
        ls, lr, lj = p
        k = terms.SHOTerm(sigma=np.exp(ls), rho=np.exp(lr), Q=1 / np.sqrt(2))
        gp = celerite2.GaussianProcess(k, mean=0.0)
        gp.compute(t, yerr=np.sqrt(sig ** 2 + np.exp(2 * lj)))
        return gp

    def nll_null(p):
        m = p[3]
        try:
            return -gp_for(p[:3]).log_likelihood(f - m)
        except Exception:  # noqa: BLE001
            return 1e25

    def nll_tr(p):
        m, lt0, rp, t14_, b = p[3:]
        try:
            mod = _shape_model(t, t0 + lt0, rp, t14_, b, exp_time) * m
            return -gp_for(p[:3]).log_likelihood(f - mod)
        except Exception:  # noqa: BLE001
            return 1e25

    b_gp = [(np.log(1e-5), np.log(0.1)), (np.log(rho_min), np.log(30.0)), (np.log(1e-6), np.log(0.05))]
    p0 = [np.log(max(np.std(f), 1e-4)), np.log(max(rho_min, 1.0)), np.log(sig / 3)]
    r0 = minimize(nll_null, p0 + [1.0], method="L-BFGS-B", bounds=b_gp + [(0.9, 1.1)])
    best = None
    for bb in (0.2, 0.6, 0.85):
        x0 = list(r0.x) + [0.0, np.sqrt(max(depth0, 1e-5)), t14, bb]
        bnds = b_gp + [(0.9, 1.1), (-t14 / 2, t14 / 2), (0.003, 0.5), (0.3 * t14, 3 * t14), (0.0, 1.3)]
        r = minimize(nll_tr, x0, method="L-BFGS-B", bounds=bnds)
        if best is None or r.fun < best.fun:
            best = r
    dlnl = float(r0.fun - best.fun)
    ls, lr, lj, m, lt0, rp, t14f, b = best.x
    tc = t0 + lt0
    dep = float(1 - _shape_model(np.array([tc]), tc, rp, t14f, b, None)[0])
    snr = float(np.sqrt(max(2 * dlnl, 0)))
    gp = gp_for(best.x[:3])
    mod = _shape_model(t, tc, rp, t14f, b, exp_time) * m
    gp.compute(t, yerr=np.sqrt(sig ** 2 + np.exp(2 * lj)))
    trend = gp.predict(f - mod, t=t, return_cov=False)
    return dict(depth=dep, snr=snr, t0=float(tc), t14=float(t14f), b=float(b), rp=float(rp),
                gp_rho_d=float(np.exp(lr)), gp_sigma=float(np.exp(ls)), dlnl=dlnl,
                arrays=dict(t=t, f=f, trend=trend + m, model=mod))


# ------------------------------------------------------------------ 3-5. Gaia, VSX

GAIA_COLS = (
    "g.source_id, g.parallax, g.parallax_error, g.phot_g_mean_mag, g.bp_rp, g.ruwe, "
    "g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.ipd_frac_odd_win, "
    "g.radial_velocity, g.radial_velocity_error, g.rv_nb_transits, g.rv_chisq_pvalue, "
    "g.rv_renormalised_gof, g.rv_amplitude_robust, g.grvs_mag, g.non_single_star, "
    "g.phot_variable_flag, g.has_epoch_photometry, "
    "ap.teff_gspphot, ap.logg_gspphot, ap.radius_gspphot, ap.ag_gspphot, "
    "ap.radius_flame, ap.radius_flame_lower, ap.radius_flame_upper, ap.lum_flame, "
    "ap.mass_flame, ap.mass_flame_lower, ap.mass_flame_upper, ap.evolstage_flame, ap.age_flame, "
    "vc.best_class_name, vc.best_class_score, eb.global_ranking AS eb_global_ranking, "
    "eb.frequency AS eb_frequency")


def gaia_star(source_id):
    from .binarity import _gaia
    q = (f"SELECT {GAIA_COLS} FROM gaiadr3.gaia_source AS g "
         "LEFT OUTER JOIN gaiadr3.astrophysical_parameters AS ap ON ap.source_id = g.source_id "
         "LEFT OUTER JOIN gaiadr3.vari_classifier_result AS vc ON vc.source_id = g.source_id "
         "LEFT OUTER JOIN gaiadr3.vari_eclipsing_binary AS eb ON eb.source_id = g.source_id "
         f"WHERE g.source_id = {int(source_id)}")
    df = _gaia(q)
    if df.empty:
        return None
    rec = df.iloc[0].to_dict()
    rec["source_id"] = int(df.source_id.iloc[0])
    return {k: (None if isinstance(v, float) and not np.isfinite(v) else v) for k, v in rec.items()}


def _num(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def stellar_check(gs, tic_rad, tic_mass, depth, dilution=1.0):
    """Dwarf status and the companion radius with Gaia's stellar radius."""
    stage = _num(gs.get("evolstage_flame")) if gs else None
    logg = _num(gs.get("logg_gspphot")) if gs else None
    r_fl = _num(gs.get("radius_flame")) if gs else None
    r_gs = _num(gs.get("radius_gspphot")) if gs else None
    plx = _num(gs.get("parallax")) if gs else None
    gmag = _num(gs.get("phot_g_mean_mag")) if gs else None
    ag = _num(gs.get("ag_gspphot")) if gs else None
    radius, source = (r_fl, "Gaia FLAME") if r_fl else ((r_gs, "Gaia GSP-Phot") if r_gs else (tic_rad, "TIC"))
    mg = (gmag + 5 * np.log10(plx / 100.0) - (ag or 0.0)) if (gmag and plx and plx > 0) else None
    if stage is not None:
        cls = "dwarf" if stage <= 360 else ("subgiant" if stage <= 490 else "giant")
        basis = f"FLAME evolutionary stage {stage:.0f}"
    elif logg is not None:
        cls = "dwarf" if logg >= 3.9 else ("subgiant" if logg >= 3.5 else "giant")
        basis = f"GSP-Phot log g {logg:.2f}"
    elif radius:
        cls = "dwarf" if radius <= 1.6 else ("subgiant" if radius <= 3 else "giant")
        basis = (f"{source} radius {radius:.2f} R☉ only; Gaia DR3 has no evolutionary "
                 "parameters for this star")
    else:
        cls, basis = "unknown", "no Gaia astrophysical parameters"
    if cls == "dwarf" and radius and radius > 2.0:
        cls, basis = "subgiant", basis + f", but radius {radius:.2f} R☉"
    rp = float(np.sqrt(max(depth, 0) * dilution) * radius * RSUN_RJUP) if radius else None
    rp_tic = float(np.sqrt(max(depth, 0) * dilution) * tic_rad * RSUN_RJUP) if tic_rad else None
    mass = _num(gs.get("mass_flame")) if gs else None
    mass = mass or tic_mass
    return dict(cls=cls, basis=basis, radius=radius, radius_source=source, mass=mass,
                abs_g=mg, logg=logg, teff=_num(gs.get("teff_gspphot")) if gs else None,
                rp_rj=rp, rp_rj_tic=rp_tic,
                flag=bool(cls in ("subgiant", "giant")),
                serious=bool(cls in ("subgiant", "giant") and rp is not None and rp > 2.0))


def binarity_extras(gs):
    if not gs:
        return dict(flags=[], notes=["no Gaia DR3 source"], serious=False)
    flags, notes = [], []
    ruwe = _num(gs.get("ruwe"))
    mp = _num(gs.get("ipd_frac_multi_peak"))
    gh = _num(gs.get("ipd_gof_harmonic_amplitude"))
    nrv = _num(gs.get("rv_nb_transits"))
    pv = _num(gs.get("rv_chisq_pvalue"))
    gof = _num(gs.get("rv_renormalised_gof"))
    if mp is not None and mp > 10:
        flags.append("ipd_multi_peak")
        notes.append(f"{mp:.0f} % of Gaia scans see a double image (ipd_frac_multi_peak)")
    if gh is not None and gh > 0.1 and ruwe is not None and ruwe > 1.4:
        flags.append("ipd_harmonic")
        notes.append(f"elongated image (ipd_gof_harmonic_amplitude {gh:.2f}) with RUWE {ruwe:.2f}: "
                     "a partly resolved companion")
    if nrv is not None and nrv >= 10 and pv is not None and gof is not None and pv < 0.01 and gof > 4:
        flags.append("rv_variable")
        notes.append(f"radial velocity varies (p = {pv:.1e}, renormalised GOF {gof:.1f}, "
                     f"amplitude {_num(gs.get('rv_amplitude_robust')) or float('nan'):.1f} km/s)")
    elif gs.get("radial_velocity_error") is not None and nrv:
        notes.append(f"RV error {_num(gs.get('radial_velocity_error')):.1f} km/s over {nrv:.0f} "
                     "transits: no significant scatter")
    else:
        notes.append("no Gaia RV")
    if ruwe is not None and ruwe > 1.4:
        flags.append("ruwe")
        notes.append(f"RUWE {ruwe:.2f}")
    # Image doubling or elongation means a resolved or partly resolved companion:
    # planets do orbit such stars (TOI-3837 b), so it only dilutes the depth (minor).
    # A large RV swing means a stellar-mass companion on a short orbit (serious).
    return dict(flags=flags, notes=notes, ruwe=ruwe, multi_peak=mp, harmonic=gh,
                rv_pvalue=pv, rv_gof=gof, rv_n=nrv,
                serious=bool("rv_variable" in flags))


VSX_EB = ("EA", "EB", "EW", "ELL", "E")      # VSX eclipsing types; "EP" is a planetary transit


def _vsx_eclipsing(typ: str) -> bool:
    import re as _re
    tokens = [x for x in _re.split(r"[+/|:() ]+", typ.upper()) if x]
    return any(t in VSX_EB or t.startswith(("EA", "EB", "EW")) for t in tokens)


def vsx(ra, dec, radius_arcsec=30):
    import io
    import urllib.parse
    import pandas as pd
    url = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urllib.parse.urlencode(
        {"-source": "B/vsx/vsx", "-c": f"{ra} {dec}", "-c.rs": radius_arcsec,
         "-out": "Name,Type,Period,max,min,_r", "-out.max": 20}))
    txt = net.get_text(url, "vizier")
    lines = [ln for ln in txt.splitlines() if ln and not ln.startswith("#")]
    if len(lines) < 3:
        return pd.DataFrame()
    return pd.read_csv(io.StringIO("\n".join([lines[0]] + lines[3:])), sep="\t")


def variability_check(gs, vsx_df):
    flags, notes = [], []
    if gs and gs.get("phot_variable_flag") == "VARIABLE":
        cls = gs.get("best_class_name")
        notes.append(f"Gaia DR3 lists it as variable ({cls or 'unclassified'}"
                     + (f", score {_num(gs.get('best_class_score')):.2f}" if gs.get("best_class_score") else "")
                     + ")")
        flags.append("gaia_variable")
        if cls and str(cls).startswith("ECL"):
            flags.append("gaia_eclipsing")
    if gs and _num(gs.get("eb_global_ranking")) is not None:
        flags.append("gaia_eclipsing")
        f = _num(gs.get("eb_frequency"))
        notes.append("in the Gaia DR3 eclipsing-binary table"
                     + (f" (P = {1 / f:.3f} d)" if f else ""))
    if vsx_df is not None and not vsx_df.empty:
        for _, r in vsx_df.iterrows():
            typ = str(r.get("Type", "")).strip()
            per = str(r.get("Period", "")).strip()
            dist = str(r.get("_r", "")).strip()
            note = f"VSX {str(r.get('Name', '')).strip()} ({typ}"
            if per and per.lower() != "nan":
                note += f", P = {per} d"
            note += ")"
            if dist and dist.lower() != "nan":
                note += f" at {float(dist):.0f}″"
            if typ.upper().startswith("EP"):
                notes.append(note.replace("(EP", "(EP: a known transiting planet"))
                continue
            notes.append(note)
            flags.append("vsx")
            if _vsx_eclipsing(typ):
                flags.append("vsx_eclipsing")
    if not notes:
        notes.append("not variable in Gaia DR3; no VSX entry within 30″")
    return dict(flags=sorted(set(flags)), notes=notes,
                serious=bool({"gaia_eclipsing", "vsx_eclipsing"} & set(flags)))


# ------------------------------------------------------------------ 6. density

def a_over_r(rho_cgs, per_d):
    return (G_CGS * rho_cgs * (per_d * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)


def density_fit(t, f, t0, t14, depth, exp_time, rho_star, rho_err, known_period=None,
                nwalkers=32, nsteps=2500, burn=1000, seed=1):
    """emcee fit. Single transit: parameters t0, rp, b, log10 P, rho (Gaussian prior
    rho_star +- rho_err), baseline; circular orbit. With ``known_period`` P is fixed
    and rho gets a broad prior, so the posterior rho is what the shape needs."""
    import emcee
    rng = np.random.default_rng(seed)
    sig = _robust_sigma_pt(f)
    fixedP = known_period is not None
    lo_rho = max(rho_star - 5 * rho_err, 1e-3)

    def unpack(p):
        if fixedP:
            tc, rp, b, lrho, c = p
            per = known_period
        else:
            tc, rp, b, lp, lrho, c = p
            per = 10 ** lp
        return tc, rp, b, per, 10 ** lrho, c

    def lnprob(p):
        tc, rp, b, per, rho, c = unpack(p)
        if not (0.003 < rp < 0.6 and 0 <= b < 1 + rp and abs(tc - t0) < t14 and 0.9 < c < 1.1):
            return -np.inf
        if not fixedP and not (0.5 < per < 20000):
            return -np.inf
        if fixedP:
            if not (1e-3 < rho < 100):
                return -np.inf
            lp = 0.0
        else:
            lp = -0.5 * ((rho - rho_star) / rho_err) ** 2 if rho > lo_rho else -np.inf
        a = a_over_r(rho, per)
        if a <= 1 + rp or b >= a:
            return -np.inf
        inc = np.degrees(np.arccos(b / a))
        m = transit_lc(t, tc, per, rp, a, inc, exp_time) * c
        return lp - 0.5 * np.sum(((f - m) / sig) ** 2)

    # start: period from the circular duration relation at the star's density
    rp0 = np.sqrt(max(depth, 1e-5))
    if not fixedP:
        # a/R* ~ P / (pi T14) (b = 0)  and  a/R* = (G rho P^2 / 3 pi)^(1/3)  ->  P
        P0 = float(np.clip((np.pi * t14) ** 3 * G_CGS * rho_star * 86400 ** 2 / (3 * np.pi), 1, 10000))
        p0 = np.array([t0, rp0, 0.3, np.log10(P0), np.log10(rho_star), 1.0])
        scale = np.array([t14 / 20, rp0 / 10, 0.2, 0.2, 0.05, 1e-4])
    else:
        rho0 = 3 * np.pi * (known_period / (np.pi * t14)) ** 3 / (G_CGS * (known_period * 86400) ** 2)
        p0 = np.array([t0, rp0, 0.3, np.log10(np.clip(rho0, 0.01, 50)), 1.0])
        scale = np.array([t14 / 20, rp0 / 10, 0.2, 0.2, 1e-4])
    ndim = len(p0)
    start = []
    while len(start) < nwalkers:
        q = p0 + scale * rng.standard_normal(ndim)
        q[2] = abs(q[2])
        if np.isfinite(lnprob(q)):
            start.append(q)
        elif len(start) == 0 and rng.random() < 0.01:
            p0[2] = 0.1
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
    sampler.run_mcmc(np.array(start), nsteps, progress=False)
    ch = sampler.get_chain(discard=burn, flat=True, thin=5)
    lnp = sampler.get_log_prob(discard=burn, flat=True, thin=5)
    out = dict(n=len(ch), acceptance=float(np.mean(sampler.acceptance_fraction)))
    q = lambda x: [float(v) for v in np.percentile(x, [16, 50, 84])]  # noqa: E731
    out["t0"] = q(ch[:, 0])
    out["rp"] = q(ch[:, 1])
    out["b"] = q(ch[:, 2])
    if fixedP:
        out["rho"] = q(10 ** ch[:, 3])
        out["rho_samples"] = (10 ** ch[:, 3]).tolist()[::10]
    else:
        out["P"] = q(10 ** ch[:, 3])
        out["rho"] = q(10 ** ch[:, 4])
        out["P_samples"] = (10 ** ch[:, 3]).tolist()
    best = ch[np.argmax(lnp)]
    tc, rp, b, per, rho, c = unpack(best)
    a = a_over_r(rho, per)
    out["best"] = dict(t0=float(tc), rp=float(rp), b=float(b), P=float(per), rho=float(rho),
                       depth=float(1 - transit_lc(np.array([tc]), tc, per, rp, a,
                                                  np.degrees(np.arccos(b / a)))[0]))
    out["model"] = (transit_lc(t, tc, per, rp, a, np.degrees(np.arccos(b / a)), exp_time) * c).tolist()
    out["grazing_prob"] = float(np.mean(ch[:, 2] > 1 - ch[:, 1]))
    return out


def ecc_factor_range(e_max=E_MAX):
    """Range of P_true / P_circ for the same transit shape and stellar density:
    the duration scales with g = (1 + e sin w) / sqrt(1 - e^2), P with g^3."""
    g_lo = (1 - e_max) / np.sqrt(1 - e_max ** 2)
    g_hi = (1 + e_max) / np.sqrt(1 - e_max ** 2)
    return float(g_lo ** 3), float(g_hi ** 3)


def period_consistency(P_samples, intervals, floor, pmax, e_max=E_MAX):
    """Posterior mass of the circular-orbit period inside the periods TESS allows,
    and whether any allowed period is within the eccentricity range of it."""
    P = np.asarray(P_samples, float)
    allowed = np.zeros(len(P), bool)
    for lo, hi in intervals or []:
        allowed |= (P >= lo) & (P <= hi)
    if floor:
        allowed |= P >= floor
    f_lo, f_hi = ecc_factor_range(e_max)
    p16, p50, p84 = np.percentile(P, [16, 50, 84])
    lo_e, hi_e = p16 * f_lo, p84 * f_hi
    reach = any(hi >= lo_e and lo <= hi_e for lo, hi in (intervals or [])) or \
        (floor is not None and hi_e >= floor)
    return dict(frac_allowed=float(allowed.mean()), P_circ=[float(p16), float(p50), float(p84)],
                ecc_range=[float(lo_e), float(hi_e)], reachable_with_ecc=bool(reach))
