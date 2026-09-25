"""Diagnostic figures."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .detect import Event  # noqa: E402
from .pipeline import SectorSearch  # noqa: E402


def binned(t: np.ndarray, f: np.ndarray, width: float):
    edges = np.arange(t.min(), t.max() + width, width)
    idx = np.digitize(t, edges)
    tb, fb = [], []
    for k in np.unique(idx):
        sel = idx == k
        if sel.sum() >= 3:
            tb.append(t[sel].mean())
            fb.append(np.median(f[sel]))
    return np.array(tb), np.array(fb)


def event_figure(ss: SectorSearch, ev: Event, title: str, path: str,
                 threshold: float, known_t0: list[float] | None = None,
                 known_label: str = "published T0") -> None:
    """Four-panel diagnostic for one detected single-transit event."""
    lc, res = ss.lc, ss.search
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 0.8, 1.3])

    ax = fig.add_subplot(gs[0, :])
    ax.plot(lc.time, lc.flux, ",", color="0.5", rasterized=True)
    ax.plot(lc.time, ss.trend, "-", color="C1", lw=1.2, label="biweight trend")
    ax.set_ylabel("PDCSAP flux")
    ax.set_title(title)
    ax.legend(loc="lower left", fontsize=8)

    ax = fig.add_subplot(gs[1, :])
    ax.plot(ss.time, ss.flat, ",", color="0.4", rasterized=True)
    tb, fb = binned(ss.time, ss.flat, 1 / 24)
    ax.plot(tb, fb, ".", color="C0", ms=3, label="1 h bins")
    for e in res.events:
        ax.axvspan(e.t0 - e.duration / 2, e.t0 + e.duration / 2,
                   color="C3" if e is ev else "C2", alpha=0.25)
    shown = [t for t in known_t0 or [] if lc.time.min() <= t <= lc.time.max()]
    for n, t0 in enumerate(shown):
        ax.axvline(t0, color="k", ls=":", lw=1, label=known_label if n == 0 else None)
    ax.set_ylabel("Detrended flux")
    ax.legend(loc="lower left", fontsize=8)
    lo, hi = np.percentile(ss.flat, [0.2, 99.8])
    ax.set_ylim(min(lo, 1 - 1.5 * ev.depth), hi)

    ax = fig.add_subplot(gs[2, :], sharex=fig.axes[0])
    ax.plot(res.grid, res.max_snr, "-", color="k", lw=0.8)
    ax.axhline(threshold, color="C3", ls="--", lw=1, label=f"threshold {threshold:g}")
    ax.set_ylabel("Box SNR\n(best duration)")
    ax.set_xlabel("Time [BTJD]")
    ax.legend(loc="upper left", fontsize=8)

    # Zoom on the event with the best-fitting box.
    ax = fig.add_subplot(gs[3, 0])
    w = max(3 * ev.duration, 0.5)
    sel = np.abs(ss.time - ev.t0) < w
    ax.plot((ss.time[sel] - ev.t0) * 24, ss.flat[sel], ".", color="0.6", ms=2)
    tb, fb = binned(ss.time[sel], ss.flat[sel], 0.5 / 24)
    ax.plot((tb - ev.t0) * 24, fb, "o", color="C0", ms=3, label="30 min bins")
    tt = np.linspace(-w, w, 2000)
    model = np.where(np.abs(tt) < ev.duration / 2, 1 - ev.depth, 1.0)
    ax.plot(tt * 24, model, "-", color="C3", lw=1.5, label="best box")
    ax.set_xlabel(f"Hours from {ev.t0:.3f} BTJD")
    ax.set_ylabel("Detrended flux")
    ax.legend(loc="lower left", fontsize=8)

    # SNR vs trial duration at the event time.
    ax = fig.add_subplot(gs[3, 1])
    k = int(np.argmin(np.abs(res.grid - ev.t0)))
    ax.plot(res.durations * 24, res.snr[:, k], "o-", color="k")
    ax.set_xscale("log")
    ax.set_xlabel("Trial duration [h]")
    ax.set_ylabel("SNR at event")
    ax.text(0.03, 0.95,
            f"Sector {lc.sector}\nT0 = {ev.t0:.3f} BTJD\n"
            f"depth = {ev.depth * 1e6:.0f} ppm\nduration = {ev.duration_h:.1f} h\n"
            f"SNR = {ev.snr:.1f}",
            transform=ax.transAxes, va="top", fontsize=9, family="monospace")

    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def overview_figure(searches: list[SectorSearch], title: str, path: str,
                    threshold: float, known_t0: list[float] | None = None) -> None:
    """Best box SNR across every sector, one strip per observing season."""
    fig, ax = plt.subplots(figsize=(12, 3.5))
    for ss in searches:
        ax.plot(ss.search.grid, ss.search.max_snr, "-", color="k", lw=0.5)
        for e in ss.search.events:
            ax.plot(e.t0, e.snr, "v", color="C3", ms=7)
    for t0 in known_t0 or []:
        ax.axvline(t0, color="C0", ls=":", lw=1)
    ax.axhline(threshold, color="C3", ls="--", lw=1)
    ax.set_ylim(-8, None)
    ax.set_xlabel("Time [BTJD]")
    ax.set_ylabel("Box SNR")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
