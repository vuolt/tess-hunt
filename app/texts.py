"""Plain-English text: glossary, forum post, ExoFOP CTOI fields, TESS coverage."""

from __future__ import annotations

import html
import io
import os
import zipfile

import numpy as np
import pandas as pd

from . import data

PHT_FORUM = "https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/talk"
EXOFOP = "https://exofop.ipac.caltech.edu/tess/"
EXOFOP_CTOI_LIST = "https://exofop.ipac.caltech.edu/tess/view_ctoi.php"
EXOFOP_BULK_PARAMS = "https://exofop.ipac.caltech.edu/tess/add_new.php?param=bulkparams"
EXOFOP_GUIDELINES = "https://exofop.ipac.caltech.edu/tess/candidate_help.php"
EXOFOP_REQUEST = "https://exofop.ipac.caltech.edu/tess/pub_candidate_upload_request.php"
EXOFOP_TEMPLATE = "https://exofop.ipac.caltech.edu/tess/templates/params_planet_YYYYMMDD_001.txt"

GLOSSARY = {
    "TESS": "NASA's Transiting Exoplanet Survey Satellite. It watches large patches of sky for about "
            "27 days at a time, measuring the brightness of hundreds of thousands of stars.",
    "sector": "One ~27-day stretch of TESS observations of one patch of sky. Sectors are numbered in order.",
    "TIC": "TESS Input Catalog: the star catalogue TESS uses. Every star has a TIC number, e.g. TIC 95747180.",
    "light curve": "A star's brightness measured over time. A planet passing in front of the star shows "
                   "up as a small, short dip.",
    "transit": "A planet passing in front of its star, blocking a little of its light.",
    "dip": "A short drop in a star's brightness. It may be a transit, or noise, starspots, another star, "
           "or an instrument effect.",
    "single transit": "Only one dip seen. The planet's orbit is then too long for TESS to have caught a "
                      "second one, which makes these long-period planets rare and valuable.",
    "duotransit": "Two dips of the same star, seen in different sectors. The time between them gives a "
                  "short list of possible orbital periods.",
    "SNR": "Signal-to-noise ratio: how much the dip stands out from the star's normal flicker. Above ~7 "
           "is a clear detection; below that it may be noise.",
    "ppm": "Parts per million. A 1000 ppm dip means the star got 0.1 % fainter.",
    "depth": "How much fainter the star gets during the dip. A Jupiter-sized planet in front of a "
             "Sun-like star gives about 1 % (10,000 ppm).",
    "duration": "How long the dip lasts. Planets on longer orbits move more slowly and give longer dips.",
    "period": "The time a planet takes to orbit its star once.",
    "alias": "One of several periods that fit two dips equally well (the gap divided by 1, 2, 3, …).",
    "epoch": "The time of the middle of the dip, given as a Julian Date so observers can predict the next one.",
    "BTJD": "Barycentric TESS Julian Date: a date in days used by TESS (BJD − 2,457,000).",
    "FPP": "False-positive probability: the chance, from the TRICERATOPS statistical test, that the dip "
           "is caused by something other than a planet around this star (e.g. an eclipsing binary).",
    "false positive": "A dip that looks like a planet but is not: often a pair of stars eclipsing each "
                      "other, a nearby star's variability, or an instrument effect.",
    "eclipsing binary": "Two stars orbiting each other so that one passes in front of the other. The most "
                        "common impostor of a planet.",
    "vetting": "The checks that try to rule out every non-planet explanation for a dip.",
    "TOI": "TESS Object of Interest: a planet candidate found by the TESS team's own pipelines.",
    "CTOI": "Community TOI: a planet candidate found by someone outside the TESS team and submitted to ExoFOP.",
    "ExoFOP": "The Exoplanet Follow-up Observing Program website, where TESS candidates are listed so "
              "astronomers can observe them.",
    "TFOP": "TESS Follow-up Observing Program: the network of astronomers who observe TESS candidates from "
            "the ground to confirm or rule them out.",
    "Planet Hunters TESS": "A Zooniverse citizen-science project whose forum discusses possible planets "
                           "found in TESS data.",
    "Gaia": "ESA's space telescope that measures the positions, distances and motions of over a billion stars.",
    "RUWE": "A Gaia number that is close to 1 for single stars; above ~1.4 the star's motion looks wobbly, "
            "often because of an unseen companion star.",
    "dwarf": "A normal, main-sequence star like the Sun. 'Subgiant' and 'giant' stars are older and bigger, "
             "which makes any object in front of them bigger too.",
    "aperture": "The group of camera pixels added up to measure a star's brightness.",
    "centroid": "The light-weighted centre of a star's image. If it moves during a dip, the dip may come "
                "from a neighbouring star.",
    "GP": "Gaussian process: a flexible statistical model of a star's natural flicker, used to separate it "
          "from a transit.",
    "R⊕": "Earth radius (6,371 km). Neptune is 3.9 R⊕, Jupiter 11.2 R⊕.",
    "R_J": "Jupiter radius (71,492 km), about 11.2 Earth radii.",
}


def term(word: str, label: str | None = None) -> str:
    """Inline tooltip (HTML <abbr>) for a glossary term; plain text if unknown."""
    key = next((k for k in GLOSSARY if k.lower() == word.lower()), None)
    shown = html.escape(label or word)
    if key is None:
        return shown
    return (f'<abbr title="{html.escape(GLOSSARY[key])}" '
            f'style="text-decoration: underline dotted; cursor: help;">{shown}</abbr>')


# ------------------------------------------------------------------ forum post

def forum_post(c: dict, followup_row) -> str:
    """Ready-to-paste Planet Hunters TESS Talk post (markdown)."""
    btjd = c.get("epoch_btjd")
    bjd = f"{btjd + 2457000:.4f}" if btjd is not None and np.isfinite(btjd) else "?"
    size = c["size"].split(":")[0]
    lines = [
        f"**Possible long-period planet: single transit on TIC {c['tic']} (TESS Sector {c['sector']})**",
        "",
        f"Hi all, I found one transit-like dip on TIC {c['tic']} and would value your opinion.",
        "",
        f"- Star: {c['star']}, T = {c['tmag']:.1f}" if c.get("tmag") else f"- Star: {c['star']}",
        f"- Mid-transit: BJD {bjd} (Sector {c['sector']})",
        f"- Depth: {c['depth_ppm']:.0f} ppm; duration: {c['duration_h']:.1f} h",
        f"- Implied size: {size}",
        f"- Orbit: {c['orbit']}",
    ]
    if c.get("fpp") is not None and np.isfinite(c["fpp"]):
        lines.append(f"- TRICERATOPS false-positive probability: {c['fpp']:.3f}")
    if c.get("verdict"):
        lines.append(f"- Expert checks (aperture, reprocessing, Gaia, variability, stellar density): {c['verdict']}")
    lines += [
        "",
        "It passed our automated vetting (shape, duration, data gaps, raw pixels, centroid, "
        "neighbouring stars, asteroids, known binaries, false-positive probability) and a search "
        "of every other TESS sector. Plots attached: the light curve, the vetting sheet and the "
        "follow-up sheet.",
        "",
        "Does anyone see a reason this is not a planet (e.g. a nearby eclipsing binary or a "
        "systematic in this sector), or know of other data on this star?",
        "",
        "(Found with the open tess-hunt single-transit pipeline.)",
    ]
    return "\n".join(lines)


def plot_bundle(plots: dict) -> bytes:
    """The candidate's key plots zipped in memory (nothing written to disk)."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, p in plots.items():
            if os.path.exists(p):
                z.write(p, arcname=f"{name}_{os.path.basename(p)}")
    return buf.getvalue()


# ------------------------------------------------------------------ ExoFOP

EXOFOP_COLUMNS = (
    "target|flag|disp|discovery|candname|period|period_unc|epoch|epoch_unc|depth|depth_unc|"
    "duration|duration_unc|inc|inc_unc|imp|imp_unc|r_planet|r_planet_err|ar_star|ar_star_err|"
    "radius|radius_unc|mass|mass_unc|temp|temp_unc|insol|insol_unc|dens|dens_unc|sma|sma_unc|"
    "ecc|ecc_unc|arg_peri|arg_peri_err|time_peri|time_peri_unc|vsa|vsa_unc|tag|group|"
    "prop_period|paper|notes").split("|")


def _epoch_unc(t14_d, depth, snr):
    """Carter et al. (2008): sigma_tc ~ (T / Q) sqrt(tau / 2T), tau = ingress time."""
    rp = np.sqrt(max(depth, 1e-8))
    tau = t14_d * rp / (1 + rp)
    return float(t14_d / max(snr, 1) * np.sqrt(tau / (2 * t14_d)))


def ctoi_fields(c: dict, followup_row) -> dict:
    """CTOI parameters from the results (the ExoFOP-ready CSV when the pipeline
    made one for this star, otherwise computed from the follow-up table)."""
    ready = None
    for p in (data.path("results", "phase5", "ctoi_candidates.csv"),
              data.path("results", "phase4", data._sub(c["sector"]), "ctoi_candidates.csv")):
        if os.path.exists(p):
            df = pd.read_csv(p)
            m = df[df["TIC ID"] == c["tic"]]
            if len(m):
                ready = m.iloc[0]
                break
    depth = float(c["depth_ppm"])
    t14_h = float(c["duration_h"])
    snr = float(c["snr"]) if c.get("snr") is not None and np.isfinite(c["snr"]) else 10.0
    epoch = float(c["epoch_btjd"]) + 2457000.0
    f = dict(
        TIC=c["tic"], sector=c["sector"],
        epoch_bjd=round(float(ready["Transit Epoch (BJD)"]), 5) if ready is not None else round(epoch, 5),
        epoch_unc=round(float(ready["Transit Epoch error (d)"]), 5) if ready is not None
        else round(_epoch_unc(t14_h / 24, depth * 1e-6, snr), 5),
        depth_ppm=round(depth), depth_unc=round(depth / snr),
        duration_h=round(t14_h, 2), duration_unc=round(2 * t14_h / snr, 2),
        radius_re=round(float(c["rp_rj"]) * data.RJ_RE, 1) if c.get("rp_rj") is not None else None,
        impact=round(float(ready["Impact parameter"]), 2) if ready is not None else None,
        period=None, period_note=c["orbit"],
    )
    jp = followup_row.get("joint_periods")
    if isinstance(jp, str) and jp.strip() and len(jp.split(";")) == 1:
        f["period"] = round(float(jp), 4)
    f["vetting"] = ("Single transit; passes 7-step vetting and multi-sector search (tess-hunt). "
                    + (f"FPP {c['fpp']:.3f}. " if c.get("fpp") is not None and np.isfinite(c["fpp"]) else "")
                    + (f"Expert checks: {c['verdict']}." if c.get("verdict") else ""))
    return f


def exofop_row(f: dict, tag: str, paper: str, target_suffix: str = "01") -> str:
    """One line of ExoFOP's planet-parameter bulk-upload file (pipe-delimited,
    column order from ExoFOP's template). New CTOIs need a paper URL and a tag."""
    notes = (f"TESS S{f['sector']} single transit; " + (f"P={f['period']} d" if f["period"]
             else "period unconstrained") + "; tess-hunt")[:120]
    vals = dict.fromkeys(EXOFOP_COLUMNS, "")
    vals.update(target=f"TIC{f['TIC']}.{target_suffix}", flag="newctoi", disp="PC", discovery="TESS",
                period=f["period"] or "", epoch=f["epoch_bjd"], epoch_unc=f["epoch_unc"],
                depth=f["depth_ppm"], depth_unc=f["depth_unc"], duration=f["duration_h"],
                duration_unc=f["duration_unc"], imp=f["impact"] if f["impact"] is not None else "",
                radius=f["radius_re"] if f["radius_re"] is not None else "",
                tag=tag, prop_period=0, paper=paper, notes=notes)
    return "|".join(str(vals[c]) for c in EXOFOP_COLUMNS)


def exofop_file(rows: list[str]) -> str:
    return "|".join(EXOFOP_COLUMNS) + "\n" + "\n".join(rows) + "\n"


# ------------------------------------------------------------------ TESS coverage

def future_tess(tic: int, ra: float, dec: float, after_btjd: float) -> list[dict]:
    """Future sectors that put the star on a TESS camera (tess-point's pointing
    plan, offline), with approximate dates."""
    try:
        import tess_stars2px as t2
        from tess_stars2px import tess_stars2px_function_entry as t2p
    except ImportError:
        return []
    sc = t2.TESS_Spacecraft_Pointing_Data()
    mids = {int(s): float(m) - 2457000.0 for s, m in zip(sc.sectors, sc.midtimes) if 0 < s < 1000}
    out = t2p(int(tic), float(ra), float(dec))
    secs = sorted({int(s) for s in out[3] if 0 < s < 1000})
    res = []
    for s in secs:
        mid = mids.get(s)
        if mid is not None and mid > after_btjd:
            from astropy.time import Time
            res.append(dict(sector=s, approx_mid=Time(mid + 2457000.0, format="jd").iso[:10]))
    return res
