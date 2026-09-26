"""TESS Input Catalog lookups for target selection."""

from __future__ import annotations

import warnings

import numpy as np

COLUMNS = ["ID", "Tmag", "Teff", "logg", "rad", "mass", "lumclass", "objType",
           "disposition", "contratio", "ra", "dec", "GAIA", "GAIAmag", "pmRA", "pmDEC"]


def query_ids(ids, chunk: int = 2000, retries: int = 4, progress=None):
    """Return an astropy Table of TIC rows for ``ids`` (queried in chunks)."""
    from astropy.table import vstack
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from astroquery.mast import Catalogs
    from . import net
    ids = [int(i) for i in ids]
    parts = []
    for k in range(0, len(ids), chunk):
        part = ids[k:k + chunk]
        t = net.cached_call("mast_catalog", ("tic_ids", tuple(part)),
                            lambda part=part: Catalogs.query_criteria(catalog="Tic", ID=part))
        parts.append(t[[c for c in COLUMNS if c in t.colnames]])
        if progress:
            progress(min(k + chunk, len(ids)), len(ids))
    return vstack(parts)


def select_dwarfs(table, tmag_max: float = 13.0):
    """Boolean mask: TIC luminosity class DWARF, Tmag < tmag_max, real stars.

    TIC v8 ``lumclass`` is DWARF/GIANT/SUBGIANT from Teff and radius (Stassun
    et al. 2019). Rows flagged as artifacts, duplicates or split sources are
    dropped.
    """
    def col(name, fill):
        c = table[name]
        return np.asarray(c.filled(fill) if hasattr(c, "filled") else c)

    lum = col("lumclass", "")
    obj = col("objType", "")
    disp = col("disposition", "")
    tmag = col("Tmag", np.nan).astype(float)
    return ((lum == "DWARF") & (obj == "STAR") & (tmag < tmag_max)
            & ~np.isin(disp, ["ARTIFACT", "DUPLICATE", "SPLIT"]))
