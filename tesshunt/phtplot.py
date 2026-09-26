"""Light-curve plot in the style of Planet Hunters TESS (Zooniverse).

PHT shows the whole sector as it comes from the SPOC pipeline: normalised
PDCSAP brightness (not detrended) against days since the start of the
sector, one dot per cadence. The same view makes a candidate easy to discuss
on the PHT forum. A second panel zooms in on the dip, with 30-min bins.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

DOT = "#2b2d42"
BIN = "#ef476f"
MARK = "#ffd166"


def _binned(t, f, width):
    if len(t) == 0:
        return t, f
    edges = np.arange(t.min(), t.max() + width, width)
    idx = np.digitize(t, edges)
    tb, fb = [], []
    for i in np.unique(idx):
        m = idx == i
        if m.sum() >= 2:
            tb.append(t[m].mean())
            fb.append(np.median(f[m]))
    return np.array(tb), np.array(fb)


def pht_plot(time, flux, t0: float, t14_d: float, depth: float, title: str, path: str,
             sector: int | None = None, zoom_d: float = 1.0) -> str:
    """Write the plot to `path` and return it. `time` in BTJD, `flux` normalised."""
    time, flux = np.asarray(time, float), np.asarray(flux, float)
    ok = np.isfinite(time) & np.isfinite(flux)
    time, flux = time[ok], flux[ok]
    start = time.min()
    x, x0 = time - start, t0 - start
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 6.8), gridspec_kw=dict(height_ratios=[1.25, 1]))
    fig.patch.set_facecolor("white")

    a1.plot(x, flux, ".", ms=1.6, color=DOT, rasterized=True)
    a1.axvspan(x0 - max(t14_d, 0.1), x0 + max(t14_d, 0.1), color=MARK, alpha=0.45, lw=0, zorder=0)
    lo, hi = np.nanpercentile(flux, [0.2, 99.8])
    pad = 0.15 * (hi - lo)
    a1.set_ylim(min(lo, 1 - 1.4 * depth) - pad, hi + 3 * pad)
    a1.set_xlim(x.min() - 0.3, x.max() + 0.3)
    a1.set_xlabel("Days" + (f" (from the start of Sector {sector})" if sector else ""))
    a1.set_ylabel("Brightness (normalised)")
    a1.set_title(title, fontsize=11, loc="left")
    a1.annotate("dip", (x0, hi + 0.8 * pad), xytext=(x0, hi + 2.3 * pad), ha="center", va="center",
                fontsize=9, arrowprops=dict(arrowstyle="->", color="0.3"), color="0.3")

    w = max(zoom_d, 3 * t14_d)
    m = np.abs(time - t0) < w
    a2.plot((time[m] - t0) * 24, flux[m], ".", ms=2.5, color=DOT, alpha=0.6, label="each measurement")
    tb, fb = _binned(time[m], flux[m], 0.5 / 24)
    a2.plot((tb - t0) * 24, fb, "o-", ms=4, lw=1, color=BIN, label="30-min average")
    a2.axvspan(-t14_d * 12, t14_d * 12, color=MARK, alpha=0.45, lw=0, zorder=0)
    a2.axhline(1, color="0.6", lw=0.8, ls=":")
    if m.any():
        l2, h2 = np.nanpercentile(flux[m], [0.5, 99.5])
        p2 = 0.2 * (h2 - l2)
        a2.set_ylim(min(l2, 1 - 1.3 * depth) - p2, h2 + p2)
    a2.set_xlabel(f"Hours from the middle of the dip (day {x0:.2f}; BTJD {t0:.3f})")
    a2.set_ylabel("Brightness (normalised)")
    a2.legend(fontsize=8, loc="lower left", frameon=False)
    a2.text(0.99, 0.04, f"depth ≈ {depth * 1e6:.0f} ppm ({depth * 100:.2f} %), "
                        f"duration ≈ {t14_d * 24:.1f} h", transform=a2.transAxes, ha="right", fontsize=9)
    for a in (a1, a2):
        a.grid(alpha=0.2)
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path
