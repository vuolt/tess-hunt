"""TESS-SPOC full-frame-image light curves (MAST HLSP), fetched by direct URL.

Direct URLs avoid a MAST search per star, which matters when processing a
whole sector. Files are downloaded to a scratch path and deleted by the caller.
"""

from __future__ import annotations

import numpy as np
from astropy.io import fits

from . import net
from .lightcurves import SectorLC

BASE = "https://archive.stsci.edu/hlsps/tess-spoc"

# Same as lightkurve's TessQualityFlags.DEFAULT_BITMASK (17087): attitude
# tweak, safe mode, coarse/Earth point, argabrightening, desaturation, manual
# exclude, impulsive outlier, bad calibration. Scattered-light cadences are
# already NaN in TESS-SPOC PDCSAP flux.
DEFAULT_BITMASK = 17087


def target_list_url(sector: int) -> str:
    return f"{BASE}/target_lists/s{sector:04d}.csv"


def lc_path(tic: int, sector: int) -> str:
    t = f"{tic:016d}"
    return (f"s{sector:04d}/target/{t[0:4]}/{t[4:8]}/{t[8:12]}/{t[12:16]}/"
            f"hlsp_tess-spoc_tess_phot_{t}-s{sector:04d}_tess_v1_lc.fits")


def lc_url(tic: int, sector: int) -> str:
    """MAST archive URL of a TESS-SPOC FFI light curve."""
    return f"{BASE}/{lc_path(tic, sector)}"


def lc_urls(tic: int, sector: int) -> list[str]:
    """Where to look, in order: the AWS S3 mirror, then MAST."""
    return [f"{net.S3}/mast/hlsp/tess-spoc/{lc_path(tic, sector)}", lc_url(tic, sector)]


def download(tic: int, sector: int, dest: str | None = None, cache: bool = False) -> str:
    """Fetch one TESS-SPOC light curve (S3 mirror first, MAST second).

    cache=False writes to ``dest`` for the caller to delete (bulk sector
    processing, where keeping ~50 GB per sector is not an option); cache=True
    keeps the file in work/cache and returns that path (do not delete it)."""
    return net.download(lc_urls(tic, sector), net.service_of, dest=dest, cache=cache)


def read(path: str, bitmask: int = DEFAULT_BITMASK) -> tuple[SectorLC, dict]:
    """Read PDCSAP flux; returns (normalized light curve, header info)."""
    with fits.open(path, memmap=False) as h:
        hdr = h[0].header
        d = h[1].data
        time_ = np.asarray(d["TIME"], float)
        flux = np.asarray(d["PDCSAP_FLUX"], float)
        err = np.asarray(d["PDCSAP_FLUX_ERR"], float)
        q = np.asarray(d["QUALITY"], int)
        info = dict(n_cadences=len(time_), tessmag=hdr.get("TESSMAG"),
                    teff=hdr.get("TEFF"), radius=hdr.get("RADIUS"),
                    sector=hdr.get("SECTOR"), crowdsap=h[1].header.get("CROWDSAP"))
    good = (np.isfinite(time_) & np.isfinite(flux) & np.isfinite(err)
            & ((q & bitmask) == 0) & (flux > 0))
    time_, flux, err = time_[good], flux[good], err[good]
    med = np.median(flux) if good.any() else 1.0
    lc = SectorLC(tic=int(hdr["TICID"]), sector=int(hdr["SECTOR"]), time=time_,
                  flux=flux / med, flux_err=err / med, tessmag=info["tessmag"])
    return lc, info
