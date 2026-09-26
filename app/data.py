"""Read-only access to the pipeline's results for the app. Nothing here writes.

Layout (Sector 48 predates multi-sector support and uses the top-level folders):
  results/phase2/sXXXX_summary.json                    search numbers
  results/phase3[/sXXXX]/summary.json, shortlist.csv   vetting funnel, shortlist
  results/phase4[/sXXXX]/followup.csv                  submit / maybe / drop
  plots/phase4[/sXXXX]/sheets/NN_ticT_{followup,vetting}.png
  plots/phase4[/sXXXX]/pht/ticT_pht.png                light curve in the Planet Hunters TESS style
  results/phase5/phase5_checks.csv, candidates.md      expert checks (optional)
  results/phase4/ticT_predictions.csv/.md              predicted transits (optional)
  results/sectors_searched.csv                         sectors already searched (built
                                                       from the files above; not written here)
"""

from __future__ import annotations

import glob
import json
import os
import re

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RJ_RE = 11.21            # Jupiter radius in Earth radii
RN_RE = 3.88             # Neptune radius in Earth radii


def _sub(sector: int) -> str:
    return "" if int(sector) == 48 else f"s{int(sector):04d}"


def path(*parts) -> str:
    return os.path.join(ROOT, *parts)


def sectors() -> list[int]:
    """Sectors with at least Phase 2 results."""
    out = []
    for p in glob.glob(path("results", "phase2", "s*_summary.json")):
        m = re.search(r"s(\d{4})_summary", p)
        if m:
            out.append(int(m.group(1)))
    return sorted(out)


def searched() -> pd.DataFrame:
    """Sectors already searched, read from the committed results (never written here)."""
    from tesshunt import ledger
    return ledger.rebuild(write=False)


def searched_sector(sector: int) -> dict | None:
    from tesshunt import ledger
    return ledger.entry(sector)


def fp_status() -> dict:
    """How far the false-positive check of other people's candidates has got."""
    out = dict(checked=0, flagged=0, high=0, tables_date=None)
    led = path("results", "phase6", "lc_checks.csv.gz")
    if os.path.exists(led):
        d = pd.read_csv(led, usecols=lambda c: c in ("group", "in_current_tables", "tables_date"))
        cur = d[d.get("in_current_tables", True) == True]  # noqa: E712
        out["checked"] = int((cur.group == "unresolved").sum())
        if "tables_date" in d and d.tables_date.notna().any():
            out["tables_date"] = str(d.tables_date.dropna().max())
    fl = path("results", "phase6", "likely_false_positives.csv")
    if os.path.exists(fl):
        f = pd.read_csv(fl, usecols=["confidence"])
        out["flagged"], out["high"] = len(f), int((f.confidence == "high").sum())
    return out


def _json(p):
    try:
        with open(p) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def funnel(sector: int) -> list[dict]:
    """stars -> stars with dips -> candidates -> vetted -> shortlist -> Phase 4 decision,
    each with a plain-English sentence about what that step removed."""
    s2 = _json(path("results", "phase2", f"s{sector:04d}_summary.json")) or {}
    s3 = _json(path("results", "phase3", _sub(sector), "summary.json")) or {}
    steps = []
    if s2:
        steps.append(dict(step="Stars searched", n=s2.get("searched"),
                          text="Every bright, Sun-like or smaller star with a TESS light curve in "
                               "this sector was searched for a single dip."))
        steps.append(dict(step="Stars with dips", n=s2.get("stars_with_dips"),
                          text="Stars where the search found at least one dip that stands out from "
                               "the noise. Most of these are noise, starspots or instrument effects."))
        steps.append(dict(step="Candidates", n=s2.get("candidates"),
                          text="Automatic filters removed dips that repeat regularly (known variable "
                               "or binary stars) and dips that happen on many stars at once "
                               "(instrument problems)."))
    if s3:
        f = s3.get("funnel") or []
        vetted = next((x["remaining"] for x in f if x["step"].startswith("7")), None)
        steps.append(dict(step="Vetted", n=vetted,
                          text="Seven checks removed dips with the wrong shape, impossible durations, "
                               "dips at data gaps, dips not in the raw pixels or coming from a "
                               "neighbouring star, passing asteroids, known binaries, and dips a "
                               "statistical test says are probably not planets."))
        steps.append(dict(step="Shortlist (new)", n=s3.get("survivors_new"),
                          text="Of those, the ones not already known as planet candidates (TOIs or "
                               "CTOIs): the new finds."))
    fu = followup(sector)
    if fu is not None:
        c = fu[fu.role == "candidate"].category.value_counts()
        steps.append(dict(step="Worth submitting", n=int(c.get("submit", 0)),
                          text=f"After looking at every other TESS observation of these stars, Gaia "
                               f"data on companions and independent light curves: "
                               f"{int(c.get('submit', 0))} to submit, {int(c.get('maybe', 0))} maybe, "
                               f"{int(c.get('drop', 0))} dropped."))
    return steps


def followup(sector: int) -> pd.DataFrame | None:
    p = path("results", "phase4", _sub(sector), "followup.csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def phase5() -> pd.DataFrame | None:
    p = path("results", "phase5", "phase5_checks.csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def phase2_candidates(sector: int) -> pd.DataFrame | None:
    p = path("results", "phase2", f"s{sector:04d}_candidates.csv")
    return pd.read_csv(p, usecols=lambda c: c in ("tic", "teff", "radius", "ra", "dec", "tmag", "snr")) \
        if os.path.exists(p) else None


def sheets(sector: int, tic: int) -> dict:
    d = path("plots", "phase4", _sub(sector), "sheets")
    out = {}
    for kind in ("followup", "vetting"):
        m = sorted(glob.glob(os.path.join(d, f"*_tic{tic}_{kind}.png")))
        if m:
            out[kind] = m[0]
    pht = path("plots", "phase4", _sub(sector), "pht", f"tic{tic}_pht.png")
    if os.path.exists(pht):
        out["pht"] = pht
    p5 = path("plots", "phase5", f"s{sector:04d}_{tic}.png")
    if os.path.exists(p5):
        out["phase5"] = p5
    for p in glob.glob(path("plots", "phase4", f"tic{tic}_*.png")):
        out[os.path.basename(p).split("_", 1)[1].replace(".png", "")] = p
    return out


def predictions(tic: int) -> tuple[pd.DataFrame | None, str | None]:
    csv = path("results", "phase4", f"tic{tic}_predictions.csv")
    md = path("results", "phase4", f"tic{tic}_predictions.md")
    return (pd.read_csv(csv) if os.path.exists(csv) else None,
            open(md).read() if os.path.exists(md) else None)


# ------------------------------------------------------------------ plain English

def star_type(teff: float | None) -> str:
    if teff is None or not np.isfinite(teff):
        return "a star of unknown type"
    if teff >= 7500:
        return "a hot, white star (A-type)"
    if teff >= 6000:
        return "a star a bit hotter than the Sun (F-type)"
    if teff >= 5200:
        return "a Sun-like star (G-type)"
    if teff >= 3900:
        return "an orange star, cooler than the Sun (K-type)"
    return "a small red star (M-type)"


def size_words(rp_rj: float | None) -> str:
    if rp_rj is None or not np.isfinite(rp_rj):
        return "unknown size"
    re_ = rp_rj * RJ_RE
    if re_ < 1.5:
        cmp = "about Earth's size"
    elif re_ < 2.5:
        cmp = "a 'super-Earth', up to 2.5× Earth"
    elif re_ < 6:
        cmp = f"Neptune-like ({re_ / RN_RE:.1f}× Neptune)"
    elif re_ < 10:
        cmp = f"between Neptune and Jupiter, roughly Saturn-sized ({rp_rj:.1f}× Jupiter)"
    elif re_ < 16:
        cmp = f"Jupiter-like ({rp_rj:.1f}× Jupiter)"
    else:
        cmp = f"larger than Jupiter ({rp_rj:.1f}×) — more likely a small star than a planet"
    return f"{re_:.1f}× Earth's radius: {cmp}"


REASON_WORDS = {
    "D1": "an independent analysis of the same TESS data does not see the dip",
    "D2": "the dip is very long and oddly shaped, typical of stellar variability or instrument effects",
    "D3": "the star has a close companion star, and correcting for it makes the object too large for a planet",
    "D4": "a second dip in another sector does not fit any possible orbit",
    "D5": "the 'dip' is really a sudden jump in brightness (an instrument artefact)",
    "D6": "the dips repeat on the orbit of a known companion star: it is an eclipsing star, not a planet",
    "S1": "the statistical false-positive test or a reviewer note leaves some doubt",
    "S2": "Gaia suggests the star may have a close companion star",
    "S3": "the object would be larger than most planets",
    "S4": "an independent check could not confirm the dip",
    "S5": "a similar dip in another sector failed our checks, hinting at a recurring instrument effect",
}


def reason_words(reasons: str | float) -> list[str]:
    if not isinstance(reasons, str) or not reasons.strip():
        return []
    out = []
    for part in reasons.split(";"):
        m = re.match(r"\s*([DS]\d)\b", part)
        if m and REASON_WORDS.get(m.group(1)) and REASON_WORDS[m.group(1)] not in out:
            out.append(REASON_WORDS[m.group(1)])
    return out


def orbit_words(row) -> str:
    jp = row.get("joint_periods")
    if isinstance(jp, str) and jp.strip():
        ps = [float(x) for x in jp.split(";")]
        if len(ps) <= 3:
            return ("A second dip in another sector means the orbit takes "
                    + " or ".join(f"{p:.1f}" for p in ps) + " days.")
        return f"A second dip in another sector allows {len(ps)} possible orbits ({ps[-1]:.0f}–{ps[0]:.0f} days)."
    floor = row.get("p_floor_d")
    if floor and np.isfinite(floor):
        return (f"Only one dip seen, so the orbit is unknown. TESS's other observations rule out "
                f"most orbits shorter than ~{floor:.0f} days (a few narrow windows remain).")
    return "Only one dip seen, so the orbit is unknown."


def card(sector: int, row, p5: pd.DataFrame | None, p2: pd.DataFrame | None) -> dict:
    """Everything a plain-English candidate card needs."""
    tic = int(row["tic"])
    teff = distance = None
    if p2 is not None:
        m = p2[p2.tic == tic]
        if len(m):
            teff = float(m.teff.iloc[0]) if np.isfinite(m.teff.iloc[0]) else None
    rp = row.get("rp_rj_diluted", row.get("rp_rj"))
    status, why = row["category"], reason_words(row.get("reasons"))
    verdict = None
    if p5 is not None:
        m = p5[(p5.tic == tic) & (p5.sector == sector)]
        if len(m):
            r = m.iloc[0]
            verdict = r.get("verdict")
            if isinstance(r.get("phase5"), str):
                status = r["phase5"]
            if isinstance(r.get("rank_change"), str) and r["rank_change"]:
                why = [f"expert checks: {r['rank_change']}"] + why
            if "distance_pc" in r and np.isfinite(r["distance_pc"]):
                distance = float(r["distance_pc"])
            if np.isfinite(r.get("rp_rj", np.nan)):
                rp = r["rp_rj"]
    return dict(tic=tic, sector=sector, status=status, reasons=why, verdict=verdict,
                star=star_type(teff), teff=teff, distance_pc=distance,
                radius_rsun=row.get("rstar"), tmag=row.get("tmag"),
                depth_ppm=row.get("depth_ppm"), duration_h=row.get("t14_h"),
                size=size_words(rp), rp_rj=rp, orbit=orbit_words(row),
                epoch_btjd=row.get("t0_btjd"), fpp=row.get("fpp"), snr=row.get("snr"))
