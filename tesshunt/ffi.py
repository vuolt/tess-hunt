"""TESS-SPOC full-frame-image light curves (MAST HLSP), fetched by direct URL.

Direct URLs avoid a MAST search per star, which matters when processing a
whole sector. Files are downloaded to a scratch path and deleted by the caller.
"""

from __future__ import annotations

import os
import time
import urllib.request

import numpy as np
from astropy.io import fits

from .lightcurves import SectorLC

BASE = "https://archive.stsci.edu/hlsps/tess-spoc"

# Same as lightkurve's TessQualityFlags.DEFAULT_BITMASK (17087): attitude
# tweak, safe mode, coarse/Earth point, argabrightening, desaturation, manual
# exclude, impulsive outlier, bad calibration. Scattered-light cadences are
# already NaN in TESS-SPOC PDCSAP flux.
DEFAULT_BITMASK = 17087


def target_list_url(sector: int) -> str:
    return f"{BASE}/target_lists/s{sector:04d}.csv"


def lc_url(tic: int, sector: int) -> str:
    t = f"{tic:016d}"
    return (f"{BASE}/s{sector:04d}/target/{t[0:4]}/{t[4:8]}/{t[8:12]}/{t[12:16]}/"
            f"hlsp_tess-spoc_tess_phot_{t}-s{sector:04d}_tess_v1_lc.fits")


def download(tic: int, sector: int, dest: str, retries: int = 4) -> str:
    """Download one light curve to ``dest``; retries with backoff on network errors."""
    url = lc_url(tic, sector)
    for attempt in range(retries):
        try:
            tmp = dest + ".part"
            with urllib.request.urlopen(url, timeout=60) as r, open(tmp, "wb") as fh:
                fh.write(r.read())
            os.replace(tmp, dest)
            return dest
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise FileNotFoundError(url) from e
            err = e
        except Exception as e:  # noqa: BLE001  (timeouts, resets)
            err = e
        time.sleep(2 ** (attempt + 1))
    raise RuntimeError(f"download failed after {retries} tries: {url}: {err}")


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
