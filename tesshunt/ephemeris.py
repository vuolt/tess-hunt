"""Transit predictions for a duotransit: joint fit, future windows, TESS and
ground-based visibility.

  joint_fit         shared-shape transit fit to several events, one t0 each,
                    with prayer-bead (cyclic residual shift) t0 uncertainties
  predict           future mid-times for P = (t_b - t_a) / n, with errors
  tess_windows      which predicted windows fall inside future TESS sectors
  night_coverage    fraction of a window observable at night from a site
"""

from __future__ import annotations

import numpy as np

from .vetting import _robust_sigma_pt, transit_model

# LCO 1-m network sites (a standard TESS follow-up network), lat/lon/height.
SITES = {
    "north": {"Teide (Tenerife)": (28.300, -16.511, 2390),
              "McDonald (Texas)": (30.680, -104.015, 2070),
              "Haleakala (Hawaii)": (20.707, -156.257, 3055)},
    "south": {"CTIO (Chile)": (-30.167, -70.805, 2207),
              "SAAO (South Africa)": (-32.380, 20.810, 1760),
              "Siding Spring (Australia)": (-31.273, 149.071, 1165)},
}
SUN_MAX_ALT = -12.0      # deg: nautical twilight or darker
TARGET_MIN_ALT = 30.0    # deg: airmass <= 2
BASELINE_H = 1.0         # out-of-transit baseline wanted on each side


def joint_fit(events, t14_guess, rp_guess, b_guess=0.5, n_shift=200):
    """events: list of dict(t, f, t0, exp_time). Returns dict with shared
    rp, t14 (d), b and, per event, t0, sigma_formal, sigma_prayer, sigma."""
    from scipy.optimize import least_squares
    sig = [_robust_sigma_pt(e["f"]) for e in events]
    ne = len(events)

    def unpack(p):
        rp, t14, beta = p[:3]
        return rp, t14, beta * (1 + rp), p[3:3 + ne], p[3 + ne:].reshape(ne, 2)

    def models(p, evs):
        rp, t14, b, t0s, base = unpack(p)
        return [transit_model(e["t"], t0s[i], rp, t14, b, exp_time=e["exp_time"])
                * (base[i, 0] + base[i, 1] * (e["t"] - e["t0"])) for i, e in enumerate(evs)]

    def resid(p, evs):
        return np.concatenate([(e["f"] - m) / sig[i] for i, (e, m) in enumerate(zip(evs, models(p, evs)))])

    p0 = np.r_[rp_guess, t14_guess, b_guess / (1 + rp_guess), [e["t0"] for e in events],
               np.tile([1.0, 0.0], ne)]
    lo = np.r_[0.005, 0.3 * t14_guess, 0.0, [e["t0"] - t14_guess / 2 for e in events],
               np.tile([0.9, -1.0], ne)]
    hi = np.r_[0.5, 3 * t14_guess, 0.99, [e["t0"] + t14_guess / 2 for e in events],
               np.tile([1.1, 1.0], ne)]
    best = None
    for beta in (0.1, 0.5, 0.85):
        p0[2] = beta
        r = least_squares(resid, np.clip(p0, lo, hi), bounds=(lo, hi), args=(events,), x_scale="jac")
        if best is None or r.cost < best.cost:
            best = r
    p = best.x
    J = best.jac
    dof = max(len(best.fun) - len(p), 1)
    cov = np.linalg.pinv(J.T @ J) * max(2 * best.cost / dof, 1.0)
    perr = np.sqrt(np.clip(np.diag(cov), 0, None))

    # Prayer bead: shift each event's residuals cyclically, refit t0s.
    mods = models(p, events)
    shifts = []
    rng = np.linspace(0, 1, n_shift, endpoint=False)
    for frac in rng:
        evs = []
        for e, m in zip(events, mods):
            r_ = e["f"] - m
            k = int(frac * len(r_))
            evs.append(dict(e, f=m + np.roll(r_, k)))
        try:
            rr = least_squares(resid, p, bounds=(lo, hi), args=(evs,), x_scale="jac", max_nfev=200)
            shifts.append(rr.x[3:3 + ne])
        except Exception:  # noqa: BLE001
            continue
    shifts = np.array(shifts)
    rp, t14, b, t0s, base = unpack(p)
    out = dict(rp=float(rp), t14=float(t14), b=float(b), depth=float(rp ** 2),
               sigma_pt=sig, events=[], chi2_red=float(2 * best.cost / dof),
               models=mods)
    for i, e in enumerate(events):
        sf = float(perr[3 + i])
        sp = float(np.std(shifts[:, i])) if len(shifts) > 10 else np.nan
        out["events"].append(dict(t0=float(t0s[i]), sigma_formal=sf, sigma_prayer=sp,
                                  sigma=float(np.nanmax([sf, sp]))))
    return out


def predict(ta, sa, tb, sb, n, t_start, t_end):
    """Mid-times t_b + k P for P = (t_b - t_a)/n within [t_start, t_end], with
    sigma_k^2 = (1 + k/n)^2 sb^2 + (k/n)^2 sa^2 (independent t_a, t_b)."""
    P = (tb - ta) / n
    sP = np.hypot(sa, sb) / n
    k = np.arange(int(np.ceil((t_start - tb) / P)), int(np.floor((t_end - tb) / P)) + 1)
    tc = tb + k * P
    s = np.sqrt((1 + k / n) ** 2 * sb ** 2 + (k / n) ** 2 * sa ** 2)
    return dict(P=P, sigma_P=sP, k=k, tc=tc, sigma=s)


def sector_bounds(max_sector=None):
    """Approximate sector start/end (BTJD) from tess-point's sector mid-times:
    halfway to the neighbouring sectors' mid-times. Edges good to ~1 d."""
    import tess_stars2px as t2
    sc = t2.TESS_Spacecraft_Pointing_Data()
    s, m = np.array(sc.sectors), np.array(sc.midtimes) - 2457000.0
    keep = s < 1000
    s, m = s[keep], m[keep]
    o = np.argsort(s)
    s, m = s[o], m[o]
    half = np.diff(m) / 2
    start = np.r_[m[0] - half[0], m[1:] - half]
    end = np.r_[m[:-1] + half, m[-1] + half[-1]]
    return {int(a): (float(b), float(c)) for a, b, c in zip(s, start, end)}


def tess_sectors(tic, ra, dec):
    from tess_stars2px import tess_stars2px_function_entry as t2p
    out = t2p(int(tic), float(ra), float(dec))
    return sorted({int(x) for x in out[3] if 0 < x < 1000})


def night_coverage(ra, dec, t_lo, t_hi, site, step_min=5.0):
    """Fraction of [t_lo, t_hi] (BTJD) with the Sun below SUN_MAX_ALT and the
    target above TARGET_MIN_ALT at ``site`` = (lat, lon, height_m)."""
    import astropy.units as u
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
    from astropy.time import Time
    lat, lon, h = site
    loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=h * u.m)
    n = max(int((t_hi - t_lo) * 1440 / step_min), 2)
    tt = Time(np.linspace(t_lo, t_hi, n) + 2457000.0, format="jd", scale="tdb")
    frame = AltAz(obstime=tt, location=loc)
    alt_t = SkyCoord(ra * u.deg, dec * u.deg).transform_to(frame).alt.deg
    alt_s = get_sun(tt).transform_to(frame).alt.deg
    ok = (alt_s < SUN_MAX_ALT) & (alt_t > TARGET_MIN_ALT)
    return float(ok.mean()), ok, tt
