"""Physical plausibility checks and a combined planet score for TOIs / CTOIs.

Checks (each validated on known planets and known false positives before
use; see scripts/phase6_report.py):

- size: the implied companion is too big for a planet. Real inflated hot
  Jupiters reach ~2 R_J, but only on short orbits.
- density: the transit lasts longer than any orbit around a star of the
  catalogued density allows (duration / maximum central duration). Usually
  a blend or a misclassified star.
- secondary (physical): a dip near phase 0.5 deeper than the hottest
  possible planet could produce by thermal emission plus reflected light.

The planet score is a logistic regression on physically motivated features
only (no brightness or distance, which mostly reflect which targets got
followed up), trained on TFOPWG-confirmed planets vs. known false positives
and checked by cross-validation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

RJ_RE = 11.21
A_R_SUN = 4.206            # a/R_sun for P = 1 d around a star of solar density
LAMBDA_TESS = 0.8e-6       # m, middle of the TESS band


# ------------------------------------------------------------------ physics

def stellar_density(mstar, rstar, logg) -> np.ndarray:
    """Density in solar units: M/R^3 when the mass is known, else from log g."""
    m, r, g = (np.asarray(x, float) for x in (mstar, rstar, logg))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(m > 0, m / r ** 3, 10 ** (g - 4.438) / r)


def a_over_r(period_d, rho_solar) -> np.ndarray:
    with np.errstate(invalid="ignore"):
        return A_R_SUN * np.asarray(rho_solar, float) ** (1 / 3) * np.asarray(period_d, float) ** (2 / 3)


def max_duration_d(period_d, depth, a_r) -> np.ndarray:
    """Longest possible transit (central, circular orbit), in days."""
    k = np.sqrt(np.clip(np.asarray(depth, float), 0, None))
    with np.errstate(invalid="ignore"):
        return np.asarray(period_d, float) / np.pi * np.arcsin(np.clip((1 + k) / a_r, 0, 1))


def _planck(T):
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        return 1.0 / np.expm1(6.626e-34 * 2.998e8 / (LAMBDA_TESS * 1.381e-23 * np.asarray(T, float)))


def max_planet_secondary(depth, a_r, teff) -> np.ndarray:
    """Deepest secondary eclipse a planet could give: day side at the hottest
    possible temperature (no heat redistribution, zero albedo) plus full
    reflection (geometric albedo 1)."""
    k2 = np.clip(np.asarray(depth, float), 0, None)
    a_r = np.asarray(a_r, float)
    tday = np.asarray(teff, float) * np.sqrt(1 / a_r) * (2 / 3) ** 0.25
    with np.errstate(invalid="ignore"):
        return k2 * _planck(tday) / _planck(teff) + k2 / a_r ** 2


# ------------------------------------------------------------------ checks

def table_quantities(t: pd.DataFrame) -> pd.DataFrame:
    """Per candidate: radius in R_J, a/R*, duration ratio (from the tables)."""
    rho = stellar_density(t.mstar, t.rstar, t.logg)
    ar = a_over_r(t.period, rho)
    tmax = max_duration_d(t.period, t.depth_ppm * 1e-6, ar)
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = (t.duration_h / 24) / tmax
    return pd.DataFrame(dict(rp_rj=t.rp_re / RJ_RE, a_r=ar, duration_ratio=ratio), index=t.index)


def secondary_excess(lc: pd.DataFrame, t: pd.DataFrame) -> pd.Series:
    """How many sigma the measured secondary exceeds the deepest planet
    secondary (negative: a planet could make it). ``lc`` and ``t`` aligned."""
    q = table_quantities(t)
    mx = pd.Series(max_planet_secondary(t.depth_ppm * 1e-6, q.a_r, t.teff), index=t.index).fillna(0.0)
    err = lc.sec_depth / lc.sec_snr.where(lc.sec_snr > 0)
    return (lc.sec_depth - mx) / err


# ------------------------------------------------------------------ score

FEATURES = ["log_rp", "big_long", "log_duration_ratio", "log_depth", "log_period", "multi",
            "gaia_flag", "oddeven", "centroid", "secondary"]


def features(t: pd.DataFrame, lc: pd.DataFrame | None, gaia_flagged: set) -> pd.DataFrame:
    """Physically motivated features, with missing values set to neutral."""
    q = table_quantities(t)
    f = pd.DataFrame(index=t.index)
    f["log_rp"] = np.log10(q.rp_rj.clip(0.02, 10)).fillna(np.log10(0.3))
    f["big_long"] = ((q.rp_rj > 1.5) & (t.period > 10)).astype(float)
    f["log_duration_ratio"] = np.log10(q.duration_ratio.clip(0.05, 20)).fillna(0.0)
    f["log_depth"] = np.log10(t.depth_ppm.clip(10, 5e5)).fillna(3.0)
    f["log_period"] = np.log10(t.period.clip(0.2, 1000)).fillna(1.0)
    f["multi"] = (t.n_on_star > 1).astype(float)
    f["gaia_flag"] = t.name.isin(gaia_flagged).astype(float)
    f[["oddeven", "centroid", "secondary"]] = 0.0
    if lc is not None and len(lc):
        L = lc.set_index("name")
        L = L[L.n_events.fillna(0) > 0]
        m = t.name.isin(L.index)
        r = L.reindex(t.name[m])
        f.loc[m, "oddeven"] = (np.clip(r.oddeven_sigma.fillna(0).values, 0, 30)
                               * (r.oddeven_frac.fillna(0).values >= 0.2))
        f.loc[m, "centroid"] = (np.clip(r.centroid_sigma.fillna(0).values, 0, 30)
                                * (r.source_offset_px.fillna(0).values >= 1.0))
        exc = secondary_excess(r.reset_index(drop=True), t[m].reset_index(drop=True)).values
        near = (np.abs(r.sec_phase.fillna(0).values - 0.5) < 0.1)
        f.loc[m, "secondary"] = np.clip(np.nan_to_num(exc), 0, 30) * near
    return f


def fit_score(f: pd.DataFrame, y: np.ndarray, seed: int = 1):
    """Logistic regression (planet = 1) with 5-fold cross-validated
    probabilities for the training set. Returns (model, cv_probabilities)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    model = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000))
    cv = StratifiedKFold(5, shuffle=True, random_state=seed)
    p_cv = cross_val_predict(model, f[FEATURES].values, y, cv=cv, method="predict_proba")[:, 1]
    model.fit(f[FEATURES].values, y)
    return model, p_cv


def auc(y, p) -> float:
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, p))
