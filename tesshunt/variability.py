"""Stellar variability timescale and the detrending window it allows.

The dominant variability period P and semi-amplitude A come from a
Lomb-Scargle periodogram. A sliding location filter of full width W leaves a
residual of roughly

    A * |1 - sinc(W / P)|          (np.sinc(x) = sin(pi x) / (pi x))

on a sinusoid: long windows cannot follow fast variability. The window is set
to the longest value in WINDOWS whose predicted residual stays below

    tolerance = max(3 * sigma_24h, 500 ppm)

Shortening the window does not trade residual against long-transit
sensitivity; it removes long transits entirely (a transit survives only if it
is shorter than about W/3, so the longest duration searched is W/3). So the
window is only shortened when the residual would swamp a giant-planet-like
day-long transit (0.5-1 %: SNR < ~7-14 against a 500 ppm residual), or when it
exceeds a few times the white noise for faint stars. Smaller residuals are left
in and simply raise the empirical noise the SNR is measured against.

Flags
  long_limited  W < 3 d: day-long single transits are not recoverable (they
                would be absorbed by the trend), and the search is capped at W/3.
  hf_variable   no allowed window brings the residual below tolerance:
                variability faster than ~0.5 d (pulsators, contact binaries)
                that no transit-preserving filter can remove.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

WINDOWS = (3.0, 2.5, 2.0, 1.5, 1.25, 1.0, 0.75, 0.5)   # days, longest first
FULL_WINDOW = 3.0          # keeps 24 h transits (3x rule)
TOL_FACTOR = 3.0           # tolerance = max(TOL_FACTOR * sigma_24h, TOL_FLOOR)
TOL_FLOOR = 500e-6         # residual that would swamp a ~0.5-1 % day-long transit
FAP_MAX = 1e-3             # periodogram peak must be at least this significant


@dataclass
class Variability:
    period: float          # d, NaN if no significant periodicity
    amplitude: float       # fractional semi-amplitude of the best-fit sinusoid
    power: float
    fap: float
    sigma_pt: float        # point-to-point white noise per cadence
    sigma_24h: float       # white noise of a 24 h box mean
    tolerance: float
    window: float          # chosen biweight window, d
    residual: float        # predicted residual at the chosen window
    max_duration_h: float  # longest transit duration searched
    long_limited: bool
    hf_variable: bool

    def as_dict(self) -> dict:
        return asdict(self)


def point_to_point(flux: np.ndarray) -> float:
    d = np.diff(flux)
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2))


def residual_fraction(window: float, period: float) -> float:
    return float(abs(1.0 - np.sinc(window / period)))


def measure(time: np.ndarray, flux: np.ndarray, pmin: float = 0.04,
            pmax: float | None = None) -> Variability:
    from astropy.timeseries import LombScargle

    cadence = float(np.median(np.diff(time)))
    sigma_pt = point_to_point(flux)
    sigma_24h = sigma_pt / np.sqrt(1.0 / cadence)
    tol = max(TOL_FACTOR * sigma_24h, TOL_FLOOR)

    # Clip gross outliers (both directions) so they do not dominate the periodogram.
    med = np.median(flux)
    mad = 1.4826 * np.median(np.abs(flux - med))
    keep = np.abs(flux - med) < 8 * max(mad, sigma_pt)
    t, f = time[keep], flux[keep]
    pmax = pmax or (t[-1] - t[0]) / 2

    ls = LombScargle(t, f)
    freq, power = ls.autopower(minimum_frequency=1 / pmax, maximum_frequency=1 / pmin,
                               samples_per_peak=5)
    k = int(np.argmax(power))
    fbest, pbest = float(freq[k]), float(power[k])
    fap = float(ls.false_alarm_probability(pbest, method="baluev",
                                           minimum_frequency=1 / pmax,
                                           maximum_frequency=1 / pmin))
    model = ls.model(np.linspace(0, 1 / fbest, 200), fbest)
    amp = float((model.max() - model.min()) / 2)
    period = 1 / fbest

    periodic = fap < FAP_MAX
    if not periodic:
        period = np.nan
    if not (periodic and amp > tol):
        # No variability that a 3 d window would fail to follow.
        resid = amp * residual_fraction(FULL_WINDOW, period) if periodic else 0.0
        return Variability(period, amp, pbest, fap, sigma_pt, sigma_24h, tol,
                           FULL_WINDOW, resid, FULL_WINDOW * 24 / 3, False, False)

    resid = np.array([amp * residual_fraction(w, period) for w in WINDOWS])
    ok = np.where(resid <= tol)[0]
    if ok.size:
        i = int(ok[0])                       # longest window within tolerance
    else:
        # Nothing works: take the longest window whose residual is within 20 %
        # of the best achievable (shrinking the window does not help fast variables).
        i = int(np.where(resid <= 1.2 * resid.min())[0][0])
    w = WINDOWS[i]
    return Variability(period, amp, pbest, fap, sigma_pt, sigma_24h, tol, w,
                       float(resid[i]), w * 24 / 3, w < FULL_WINDOW,
                       bool(resid[i] > tol))
