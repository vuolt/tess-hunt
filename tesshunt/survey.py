"""Per-star processing for a full-sector search, with a resumable SQLite store.

Each star is: download FFI light curve -> measure variability -> choose the
detrending window -> two-pass detrend + box search -> (optionally) injection-
recovery trials -> delete the file. All rows for a star are written in one
transaction, so a star is either fully recorded or not at all, and a restart
skips every star already in the ``stars`` table.
"""

from __future__ import annotations

import os
import sqlite3
import time
import traceback

import numpy as np

from . import ffi
from .detect import DEFAULT_DURATIONS_H, box_search
from .detrend import split_segments
from .injection import inject
from .pipeline import SearchConfig, search_sector
from .variability import measure

MIN_POINTS = 1000          # ~7 d of 10-min data
INJ_DURATIONS_H = (1, 2, 4, 8, 16, 24)
INJ_SNR_RANGE = (3.0, 60.0)
INGRESS_FRAC = 0.1
GAP_FOR_EDGE = 0.1         # d; gaps longer than this count as data edges for vetting

SCHEMA = """
CREATE TABLE IF NOT EXISTS stars (
  tic INTEGER PRIMARY KEY, status TEXT, reason TEXT, tmag REAL, teff REAL,
  radius REAL, crowdsap REAL, n_points INTEGER, span_d REAL, sigma_pt_ppm REAL,
  sigma_24h_ppm REAL, var_period_d REAL, var_amp_ppm REAL, var_fap REAL,
  window_d REAL, resid_ppm REAL, tol_ppm REAL, max_dur_h REAL,
  long_limited INTEGER, hf_variable INTEGER, red_noise_4h REAL,
  n_dips INTEGER, injected INTEGER, runtime_s REAL
);
CREATE TABLE IF NOT EXISTS dips (
  tic INTEGER, rank_in_star INTEGER, t0 REAL, duration_h REAL, depth_ppm REAL,
  snr REAL, edge_dist_h REAL, n_in INTEGER
);
CREATE TABLE IF NOT EXISTS injections (
  tic INTEGER, trial INTEGER, t0 REAL, depth_ppm REAL, duration_h REAL,
  expected_snr REAL, expected_snr_white REAL, window_d REAL, max_dur_h REAL,
  long_limited INTEGER, hf_variable INTEGER, recovered INTEGER, top_ranked INTEGER,
  det_t0 REAL, det_depth_ppm REAL, det_duration_h REAL, measured_snr REAL
);
CREATE INDEX IF NOT EXISTS dips_tic ON dips(tic);
"""


class Store:
    def __init__(self, path: str):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        self.db = sqlite3.connect(path)
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.executescript(SCHEMA)

    def done(self) -> set[int]:
        return {r[0] for r in self.db.execute("SELECT tic FROM stars")}

    def write(self, res: dict) -> None:
        with self.db:   # one transaction per star
            s = res["star"]
            self.db.execute(f"INSERT OR REPLACE INTO stars ({','.join(s)}) "
                            f"VALUES ({','.join('?' * len(s))})", list(s.values()))
            self.db.execute("DELETE FROM dips WHERE tic=?", (s["tic"],))
            self.db.execute("DELETE FROM injections WHERE tic=?", (s["tic"],))
            for table in ("dips", "injections"):
                rows = res.get(table) or []
                if rows:
                    cols = list(rows[0])
                    self.db.executemany(
                        f"INSERT INTO {table} ({','.join(cols)}) VALUES ({','.join('?' * len(cols))})",
                        [list(r.values()) for r in rows])


N_MASK_FOR_VAR = 3


def analyse(lc):
    """Variability -> window -> search. Returns (Variability, SearchConfig, SectorSearch).

    A deep day-long transit has plenty of low-frequency power, and on its own
    it would make the periodogram report "variability" and shrink the window
    until the transit is filtered away. So a preliminary 3 d-window search runs
    first, and its strongest events are masked before variability is measured.
    """
    full = SearchConfig(max_events=10)
    pre = search_sector(lc, full)
    keep = np.ones(len(lc.time), bool)
    for ev in pre.search.events[:N_MASK_FOR_VAR]:
        keep &= np.abs(lc.time - ev.t0) >= ev.duration
    var = measure(lc.time[keep], lc.flux[keep])
    if var.window == full.window:
        return var, full, pre        # same config: the preliminary search is the answer
    durs = tuple(d for d in DEFAULT_DURATIONS_H if d <= var.max_duration_h + 1e-9)
    cfg = SearchConfig(window=var.window, durations_h=durs, max_events=10)
    return var, cfg, search_sector(lc, cfg)


def edge_distance(time: np.ndarray, t0: float) -> float:
    """Hours from t0 to the nearest edge of the data segment containing it."""
    best = np.inf
    for seg in split_segments(time, GAP_FOR_EDGE):
        a, b = time[seg][0], time[seg][-1]
        if a <= t0 <= b:
            return float(min(t0 - a, b - t0) * 24)
        best = min(best, abs(t0 - a), abs(t0 - b))
    return float(-best * 24)   # negative: t0 falls in a gap


def box_noise(time, flat, cadence):
    """Empirical box noise at every standard duration, plus the white-noise value."""
    res = box_search(time, flat, durations_h=DEFAULT_DURATIONS_H, threshold=np.inf)
    sig_pt = 1.4826 * np.median(np.abs(np.diff(flat))) / np.sqrt(2)
    white = {d: sig_pt / np.sqrt(d / 24 / cadence) for d in DEFAULT_DURATIONS_H}
    emp = dict(zip(DEFAULT_DURATIONS_H, res.sigma))
    return emp, white


def _star_row(tic, status, reason="", info=None, lc=None, var=None, extra=None):
    row = dict(tic=int(tic), status=status, reason=reason)
    info = info or {}
    row.update(tmag=info.get("tessmag"), teff=info.get("teff"), radius=info.get("radius"),
               crowdsap=info.get("crowdsap"))
    if lc is not None and len(lc.time):
        row.update(n_points=len(lc.time), span_d=float(lc.time[-1] - lc.time[0]))
    if var is not None:
        row.update(sigma_pt_ppm=var.sigma_pt * 1e6, sigma_24h_ppm=var.sigma_24h * 1e6,
                   var_period_d=var.period, var_amp_ppm=var.amplitude * 1e6,
                   var_fap=var.fap, window_d=var.window, resid_ppm=var.residual * 1e6,
                   tol_ppm=var.tolerance * 1e6, max_dur_h=var.max_duration_h,
                   long_limited=int(var.long_limited), hf_variable=int(var.hf_variable))
    row.update(extra or {})
    return row


def process_star(tic: int, sector: int, tmpdir: str, n_inject: int = 0,
                 seed: int = 0) -> dict:
    t_start = time.time()
    path = os.path.join(tmpdir, f"{tic}.fits")
    info = {}
    try:
        try:
            ffi.download(tic, sector, path)
            lc, info = ffi.read(path)
        finally:
            if os.path.exists(path):
                os.remove(path)
        if len(lc.time) < MIN_POINTS:
            return {"star": _star_row(tic, "skipped", "too_few_points", info, lc,
                                      extra=dict(runtime_s=time.time() - t_start))}

        var, cfg, ss = analyse(lc)
        cadence = float(np.median(np.diff(lc.time)))
        emp, white = box_noise(ss.time, ss.flat, cadence)
        dips = [dict(tic=int(tic), rank_in_star=i + 1, t0=e.t0, duration_h=e.duration_h,
                     depth_ppm=e.depth * 1e6, snr=e.snr,
                     edge_dist_h=edge_distance(ss.time, e.t0),
                     n_in=int(np.sum(np.abs(ss.time - e.t0) < e.duration / 2)))
                for i, e in enumerate(ss.search.events)]

        injections = []
        if n_inject:
            rng = np.random.default_rng(seed)
            for k in range(n_inject):
                dur_h = float(rng.choice(INJ_DURATIONS_H))
                snr_target = float(np.exp(rng.uniform(*np.log(INJ_SNR_RANGE))))
                sig = emp[dur_h] if np.isfinite(emp[dur_h]) else white[dur_h]
                depth = min(snr_target * sig / (1 - INGRESS_FRAC), 0.05)
                t0 = float(lc.time[rng.integers(len(lc.time))])
                v2, c2, s2 = analyse(inject(lc, t0, depth, dur_h / 24, INGRESS_FRAC))
                tol = max(dur_h / 48, 0.5 / 24)
                hits = [e for e in s2.search.events if abs(e.t0 - t0) < tol]
                best = max(hits, key=lambda e: e.snr) if hits else None
                near = np.abs(s2.search.grid - t0) < tol
                msnr = s2.search.max_snr[near] if near.any() else np.array([np.nan])
                injections.append(dict(
                    tic=int(tic), trial=k, t0=t0, depth_ppm=depth * 1e6, duration_h=dur_h,
                    expected_snr=depth * (1 - INGRESS_FRAC) / sig,
                    expected_snr_white=depth * (1 - INGRESS_FRAC) / white[dur_h],
                    window_d=v2.window, max_dur_h=v2.max_duration_h,
                    long_limited=int(v2.long_limited), hf_variable=int(v2.hf_variable),
                    recovered=int(best is not None),
                    top_ranked=int(best is not None and best is s2.search.events[0]),
                    det_t0=best.t0 if best else None,
                    det_depth_ppm=best.depth * 1e6 if best else None,
                    det_duration_h=best.duration_h if best else None,
                    measured_snr=float(np.nanmax(msnr)) if np.isfinite(msnr).any() else None))

        red4 = emp[4.0] / white[4.0] if np.isfinite(emp[4.0]) else None
        star = _star_row(tic, "ok", "", info, lc, var, dict(
            red_noise_4h=red4, n_dips=len(dips), injected=len(injections),
            runtime_s=time.time() - t_start))
        return {"star": star, "dips": dips, "injections": injections}
    except FileNotFoundError:
        return {"star": _star_row(tic, "skipped", "no_lc_file", info,
                                  extra=dict(runtime_s=time.time() - t_start))}
    except Exception as e:  # noqa: BLE001
        msg = f"{type(e).__name__}: {e}"[:300]
        return {"star": _star_row(tic, "error", msg, info,
                                  extra=dict(runtime_s=time.time() - t_start)),
                "traceback": traceback.format_exc()}
