"""Binarity indicators for a candidate host star (Gaia DR3 + catalogues).

All queries go through tesshunt.net (cached, rate-limited, one at a time).

  gaia_target        Gaia DR3 record of the host (by source_id from the TIC):
                     RUWE, non_single_star, image-doubling indicators, RV scatter
  gaia_nss           Gaia DR3 non-single-star orbit, if flagged
  gaia_neighbours    Gaia DR3 sources within 60": co-moving companions and
                     blends inside one TESS pixel
  elbadry_pairs      El-Badry, Rix & Heintz (2021) wide binaries (VizieR)
  wds_pairs          Washington Double Star catalogue entries (VizieR)
  assess             combine into flags and a one-line summary
"""

from __future__ import annotations

import io

import numpy as np
import pandas as pd

from . import net

GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
VIZIER_TSV = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
RUWE_MAX = 1.4                # above: astrometric excess noise, likely unresolved binary
IPD_MULTI_PEAK_MAX = 10       # %: fraction of scans with a double-peaked image
BLEND_ARCSEC = 21.0           # one TESS pixel
BLEND_DMAG = 6.0              # fainter neighbours cannot mimic a >~0.4 % dip
COMOVE_RADIUS = 60.0          # arcsec searched for co-moving companions
WDS_CLOSE_ARCSEC = 2.0        # WDS pairs this close are unresolved by Gaia/TESS...
WDS_CLOSE_DMAG = 3.0          # ...and this similar in brightness dilute/confuse the host


def _gaia(adql: str) -> pd.DataFrame:
    txt = net.get_text(GAIA_TAP, "gaia", data={"REQUEST": "doQuery", "LANG": "ADQL",
                                              "FORMAT": "csv", "QUERY": adql})
    return pd.read_csv(io.StringIO(txt))


def _vizier(params: dict) -> pd.DataFrame:
    import urllib.parse
    txt = net.get_text(f"{VIZIER_TSV}?{urllib.parse.urlencode(params)}", "vizier")
    lines = [ln for ln in txt.splitlines() if ln and not ln.startswith("#")]
    if len(lines) < 3:
        return pd.DataFrame()
    header, rows = lines[0].split("\t"), [ln.split("\t") for ln in lines[3:]]
    return pd.DataFrame(rows, columns=header).apply(lambda c: c.str.strip())


TARGET_COLS = ("source_id, ra, dec, phot_g_mean_mag, bp_rp, parallax, parallax_error, pmra, "
               "pmdec, ruwe, non_single_star, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude, "
               "radial_velocity, radial_velocity_error, rv_amplitude_robust, rv_chisq_pvalue, "
               "rv_nb_transits")


def gaia_target(source_id: int) -> dict | None:
    df = _gaia(f"SELECT {TARGET_COLS} FROM gaiadr3.gaia_source WHERE source_id = {int(source_id)}")
    if df.empty:   # DR2 id that changed in DR3: follow the DR2->DR3 crossmatch
        m = _gaia("SELECT dr3_source_id FROM gaiadr3.dr2_neighbourhood WHERE dr2_source_id = "
                  f"{int(source_id)} ORDER BY angular_distance")
        if m.empty:
            return None
        df = _gaia(f"SELECT {TARGET_COLS} FROM gaiadr3.gaia_source WHERE source_id = "
                   f"{int(m.dr3_source_id.iloc[0])}")
    return None if df.empty else df.iloc[0].to_dict()


def gaia_nss(source_id: int) -> list[dict]:
    df = _gaia("SELECT nss_solution_type, period, period_error, eccentricity FROM "
               f"gaiadr3.nss_two_body_orbit WHERE source_id = {int(source_id)}")
    return df.to_dict("records")


def gaia_neighbours(ra: float, dec: float, radius_arcsec: float = COMOVE_RADIUS) -> pd.DataFrame:
    r = radius_arcsec / 3600
    return _gaia(
        "SELECT source_id, phot_g_mean_mag, parallax, parallax_error, pmra, pmdec, ruwe, "
        f"DISTANCE(POINT(ra, dec), POINT({ra}, {dec})) * 3600 AS sep_arcsec "
        f"FROM gaiadr3.gaia_source WHERE 1 = CONTAINS(POINT(ra, dec), CIRCLE({ra}, {dec}, {r}))")


def elbadry_pairs(source_id: int) -> pd.DataFrame:
    out = []
    for col in ("Source1", "Source2"):
        df = _vizier({"-source": "J/MNRAS/506/2269/catalog", col: int(source_id),
                      "-out": "Source1,Source2,theta,sepAU,BinType,R"})
        out.append(df)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def wds_pairs(ra: float, dec: float, radius_arcsec: float = 60) -> pd.DataFrame:
    return _vizier({"-source": "B/wds/wds", "-c": f"{ra} {dec}", "-c.rs": radius_arcsec,
                    "-out": "WDS,Comp,sep2,mag1,mag2,Obs2"})


def comoving(target: dict, nb: pd.DataFrame) -> pd.DataFrame:
    """Neighbours with parallax and proper motion consistent with the target."""
    if nb.empty or not np.isfinite(target.get("parallax", np.nan)):
        return nb.iloc[0:0]
    nb = nb[nb.source_id != target["source_id"]]
    good = (nb.parallax / nb.parallax_error > 5) & (target["parallax"] / target["parallax_error"] > 5)
    dplx = np.abs(nb.parallax - target["parallax"]) / np.hypot(nb.parallax_error,
                                                               target["parallax_error"])
    mu_t = np.hypot(target["pmra"], target["pmdec"])
    dmu = np.hypot(nb.pmra - target["pmra"], nb.pmdec - target["pmdec"])
    return nb[good & (dplx < 3) & (dmu < max(0.1 * mu_t, 2.0))]


def assess(target: dict | None, nss: list, nb: pd.DataFrame, eb: pd.DataFrame,
           wds: pd.DataFrame) -> dict:
    """Binarity flags from the individual indicators (pure function)."""
    flags, notes = [], []
    ruwe = np.nan if target is None else target.get("ruwe", np.nan)
    if target is None:
        notes.append("no Gaia DR3 match")
    else:
        if np.isfinite(ruwe) and ruwe > RUWE_MAX:
            flags.append("ruwe")
            notes.append(f"RUWE {ruwe:.2f} > {RUWE_MAX}")
        if (target.get("non_single_star") or 0) > 0:
            flags.append("gaia_nss")
            per = ", ".join(f"{o['nss_solution_type']} P={o['period']:.1f} d"
                            for o in nss if o.get("period") == o.get("period"))
            notes.append("Gaia non-single-star" + (f" ({per})" if per else ""))
        if (target.get("ipd_frac_multi_peak") or 0) > IPD_MULTI_PEAK_MAX:
            flags.append("ipd_multi_peak")
            notes.append(f"double-peaked images in {target['ipd_frac_multi_peak']:.0f} % of scans")
        rvp = target.get("rv_chisq_pvalue")
        if rvp is not None and np.isfinite(rvp) and rvp < 0.01 and (target.get("rv_nb_transits") or 0) >= 10:
            flags.append("rv_variable")
            notes.append(f"Gaia RV variable (p={rvp:.1e}, amp {target.get('rv_amplitude_robust', np.nan):.1f} km/s)")
    cm = comoving(target, nb) if target is not None else nb.iloc[0:0]
    if len(cm):
        flags.append("comoving_companion")
        notes.append("co-moving Gaia companion(s) at " + ", ".join(
            f"{s:.1f}\" (G={g:.1f})" for s, g in zip(cm.sep_arcsec, cm.phot_g_mean_mag)))
    if target is not None and not nb.empty:
        g0 = target["phot_g_mean_mag"]
        close = nb[(nb.source_id != target["source_id"]) & (nb.sep_arcsec < BLEND_ARCSEC)
                   & (nb.phot_g_mean_mag < g0 + BLEND_DMAG)]
        if len(close):
            flags.append("blend_1px")
            notes.append("Gaia neighbour(s) within 1 px: " + ", ".join(
                f"{s:.1f}\" dG={g - g0:+.1f}" for s, g in zip(close.sep_arcsec, close.phot_g_mean_mag)))
    if not eb.empty:
        flags.append("wide_binary_catalogue")
        notes.append("El-Badry+21 wide binary: " + ", ".join(
            f"{float(t) * 3600:.0f}\" ({float(a):.0f} AU)"   # theta is in degrees
            for t, a in zip(eb.theta, eb.sepAU)))
    if not wds.empty:
        flags.append("wds")
        sep = pd.to_numeric(wds.sep2, errors="coerce")
        dm = (pd.to_numeric(wds.mag2, errors="coerce") - pd.to_numeric(wds.mag1, errors="coerce")).abs()
        if ((sep < WDS_CLOSE_ARCSEC) & (dm < WDS_CLOSE_DMAG)).any():
            flags.append("wds_close")
            notes.append(f"close WDS pair (< {WDS_CLOSE_ARCSEC:g}\", dmag < {WDS_CLOSE_DMAG:g}): "
                         "unresolved companion dilutes the transit")
        notes.append("WDS: " + ", ".join(f"{w}{c} sep {s}\" mags {m1}/{m2}"
                                         for w, c, s, m1, m2 in zip(wds.WDS, wds.Comp, wds.sep2,
                                                                    wds.mag1, wds.mag2)))
    # Only indicators of an unresolved companion to the host itself count
    # against the candidate; resolved wide companions are informational.
    host_binary = any(f in flags for f in ("ruwe", "gaia_nss", "ipd_multi_peak", "rv_variable",
                                           "wds_close"))
    return dict(ruwe=float(ruwe) if ruwe == ruwe else None, flags=flags, notes=notes,
                host_binary=host_binary, n_comoving=int(len(cm)),
                summary="; ".join(notes) if notes else "no binarity indicators")


def binarity(tic_row: dict) -> dict:
    """Run every query for one star; ``tic_row`` needs GAIA (source id), ra, dec."""
    sid = tic_row.get("GAIA")
    target = gaia_target(int(sid)) if sid == sid and sid is not None else None
    ra, dec = ((target["ra"], target["dec"]) if target else (tic_row["ra"], tic_row["dec"]))
    nss = gaia_nss(target["source_id"]) if target and (target.get("non_single_star") or 0) > 0 else []
    nb = gaia_neighbours(ra, dec)
    eb = elbadry_pairs(target["source_id"]) if target else pd.DataFrame()
    wds = wds_pairs(tic_row["ra"], tic_row["dec"])
    res = assess(target, nss, nb, eb, wds)
    res["gaia"] = target
    return res
