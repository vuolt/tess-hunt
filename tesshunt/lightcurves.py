"""Download and load TESS SPOC light curves.

Targets are addressed by TIC ID rather than by name: the MAST name resolver is
not always reachable, but TIC-based queries go straight to the archive.
"""

from __future__ import annotations

import glob
import os
import warnings
from dataclasses import dataclass

import numpy as np

DEFAULT_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")


@dataclass
class SectorLC:
    """One sector of normalized PDCSAP photometry."""

    tic: int
    sector: int
    time: np.ndarray  # BTJD = BJD - 2457000
    flux: np.ndarray  # normalized to median 1
    flux_err: np.ndarray
    tessmag: float | None = None

    def copy(self) -> "SectorLC":
        return SectorLC(self.tic, self.sector, self.time.copy(), self.flux.copy(),
                        self.flux_err.copy(), self.tessmag)


def _lk():
    # lightkurve emits noisy import-time warnings about optional extras.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import lightkurve
    return lightkurve


def download(tic: int, sectors=None, exptime: int = 120, author: str = "SPOC",
             data_dir: str = DEFAULT_DATA_DIR) -> list[str]:
    """Download SPOC light-curve files for a TIC ID into ``data_dir``.

    Files already present are reused. Returns the list of local FITS paths.
    """
    lk = _lk()
    result = lk.search_lightcurve(f"TIC {tic}", author=author, exptime=exptime)
    if sectors is not None:
        wanted = {int(s) for s in sectors}
        keep = [i for i, m in enumerate(result.mission)
                if int(m.split()[-1]) in wanted]
        result = result[keep]
    if len(result) == 0:
        raise ValueError(f"No {author} {exptime}s light curves for TIC {tic}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result.download_all(download_dir=data_dir, quality_bitmask="default")
    return local_files(tic, data_dir)


def local_files(tic: int, data_dir: str = DEFAULT_DATA_DIR) -> list[str]:
    pattern = os.path.join(data_dir, "mastDownload", "TESS", f"*-{tic:016d}-*",
                           "*_lc.fits")
    return sorted(glob.glob(pattern))


def read_sector(path: str) -> SectorLC:
    """Read one SPOC LC file: PDCSAP flux, default quality mask, NaNs removed."""
    lk = _lk()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lc = lk.read(path, flux_column="pdcsap_flux", quality_bitmask="default")
    time = np.asarray(lc.time.value, float)
    flux = np.asarray(lc.flux.value, float)
    err = np.asarray(lc.flux_err.value, float)
    good = np.isfinite(time) & np.isfinite(flux) & np.isfinite(err)
    time, flux, err = time[good], flux[good], err[good]
    med = np.median(flux)
    return SectorLC(
        tic=int(lc.meta["TICID"]),
        sector=int(lc.meta["SECTOR"]),
        time=time,
        flux=flux / med,
        flux_err=err / med,
        tessmag=lc.meta.get("TESSMAG"),
    )


def load_sectors(tic: int, sectors=None, data_dir: str = DEFAULT_DATA_DIR,
                 fetch: bool = True) -> list[SectorLC]:
    """Load all (or selected) sectors for a TIC ID, downloading if needed."""
    files = local_files(tic, data_dir)
    if fetch and not files:
        files = download(tic, sectors=sectors, data_dir=data_dir)
    lcs = [read_sector(f) for f in files]
    if sectors is not None:
        wanted = {int(s) for s in sectors}
        missing = wanted - {lc.sector for lc in lcs}
        if missing and fetch:
            download(tic, sectors=sorted(missing), data_dir=data_dir)
            lcs = [read_sector(f) for f in local_files(tic, data_dir)]
        lcs = [lc for lc in lcs if lc.sector in wanted]
    return sorted(lcs, key=lambda lc: lc.sector)
