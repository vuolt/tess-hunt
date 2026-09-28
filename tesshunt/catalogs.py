"""Local copies of star catalogues, downloaded once in bulk (see CLAUDE.md).

Instead of one VizieR query per star, the whole catalogue is streamed once from
the CDS archive (cdsarc.cds.unistra.fr, free and public; one request, through
the rate-limited, logged net layer), reduced to the columns the checks use and
kept in work/catalogs/ as Parquet. Lookups are then local cone searches.

  VSX  (B/vsx, AAVSO International Variable Star Index, ~10 M stars)  2.2 GB, ~35 min
  WDS  (B/wds, Washington Double Star catalogue, 158 k pairs)       6.5 MB -> 4.4 MB

    python -m tesshunt.catalogs ensure     # download what is missing or stale
    python -m tesshunt.catalogs status
"""

from __future__ import annotations

import gzip
import os
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

from . import net

CDSARC = "https://cdsarc.cds.unistra.fr/ftp"
DIR = os.path.join(os.path.dirname(net.CACHE), "catalogs")
MAX_AGE_DAYS = 180          # refresh a local copy after this long
CHUNK = 1 << 20


def _path(name: str) -> str:
    return os.path.join(DIR, f"{name}.parquet")


def _stream_lines(url: str, gz: bool):
    """Yield the lines of a (possibly gzipped) file, streamed in one request."""
    def open_stream():
        req = urllib.request.Request(url, headers={"User-Agent": "tess-hunt (research; bulk catalogue)"})
        return urllib.request.urlopen(req, timeout=600)
    # the slot spaces this request like any other small-service call; retries,
    # backoff and the halt rule apply to opening the stream
    r = net.cached_call("cdsarc", ("stream", url), open_stream, cache=False)
    t0, n, mark = time.time(), 0, 200 * CHUNK
    try:
        src = gzip.GzipFile(fileobj=r) if gz else r
        for raw in src:
            n += len(raw)
            if n >= mark:
                mark += 200 * CHUNK
                print(f"  {url.rsplit('/', 1)[1]}: {n / 1e6:.0f} MB read ({time.time() - t0:.0f} s)",
                      flush=True)
            yield raw.decode("latin-1").rstrip("\n")
    finally:
        r.close()


def _num(s: str):
    s = s.strip()
    try:
        return float(s) if s else np.nan
    except ValueError:
        return np.nan


def _build_vsx() -> pd.DataFrame:
    rows = []
    for ln in _stream_lines(f"{CDSARC}/B/vsx/vsx.dat", gz=False):
        ra, dec = _num(ln[42:51]), _num(ln[52:61])
        if np.isnan(ra) or np.isnan(dec):
            continue
        rows.append((ln[9:39].strip(), ra, dec, ln[62:92].strip(), _num(ln[95:102]),
                     _num(ln[120:127]), _num(ln[158:177])))
    return pd.DataFrame(rows, columns=["Name", "ra", "dec", "Type", "max", "min", "Period"])


def _build_wds() -> pd.DataFrame:
    rows = []
    for ln in _stream_lines(f"{CDSARC}/B/wds/wds.dat.gz", gz=True):
        rah, ram, ras = _num(ln[112:114]), _num(ln[114:116]), _num(ln[116:121])
        ded, dem, des = _num(ln[122:124]), _num(ln[124:126]), _num(ln[126:130])
        if np.isnan(rah) or np.isnan(ded):
            continue
        ra = 15 * (rah + np.nan_to_num(ram) / 60 + np.nan_to_num(ras) / 3600)
        dec = ded + np.nan_to_num(dem) / 60 + np.nan_to_num(des) / 3600
        if ln[121:122] == "-":
            dec = -dec
        rows.append((ln[0:10].strip(), ln[17:22].strip(), ra, dec, _num(ln[52:58]),
                     _num(ln[58:64]), _num(ln[64:69]), _num(ln[28:32])))
    return pd.DataFrame(rows, columns=["WDS", "Comp", "ra", "dec", "sep2", "mag1", "mag2", "Obs2"])


BUILDERS = {"vsx": _build_vsx, "wds": _build_wds}


def age_days(name: str) -> float | None:
    p = _path(name)
    return (time.time() - os.path.getmtime(p)) / 86400 if os.path.exists(p) else None


def ensure(names=tuple(BUILDERS), force: bool = False):
    """Download (once) the catalogues that are missing or older than MAX_AGE_DAYS."""
    os.makedirs(DIR, exist_ok=True)
    for name in names:
        age = age_days(name)
        if not force and age is not None and age < MAX_AGE_DAYS:
            continue
        print(f"downloading {name} from the CDS archive (one request)...", flush=True)
        df = BUILDERS[name]()
        df = df.sort_values("dec", ignore_index=True)
        tmp = _path(name) + ".part"
        df.to_parquet(tmp, index=False)
        os.replace(tmp, _path(name))
        _TABLES.pop(name, None)
        print(f"  {name}: {len(df):,} rows, {os.path.getsize(_path(name)) / 1e6:.1f} MB", flush=True)


_TABLES: dict[str, pd.DataFrame] = {}


def available(name: str) -> bool:
    return os.path.exists(_path(name))


def cone(name: str, ra: float, dec: float, radius_arcsec: float) -> pd.DataFrame:
    """Rows of a local catalogue within radius_arcsec of (ra, dec), nearest first,
    with the separation in arcseconds as '_r'."""
    if name not in _TABLES:
        _TABLES[name] = pd.read_parquet(_path(name))
    t = _TABLES[name]
    r = radius_arcsec / 3600
    d = t["dec"].to_numpy()
    lo, hi = np.searchsorted(d, dec - r), np.searchsorted(d, dec + r, side="right")
    sub = t.iloc[lo:hi]
    if sub.empty:
        return sub.assign(_r=[])
    ra1, de1 = np.radians(sub["ra"].to_numpy()), np.radians(sub["dec"].to_numpy())
    ra0, de0 = np.radians(ra), np.radians(dec)
    sep = 2 * np.degrees(np.arcsin(np.sqrt(np.sin((de1 - de0) / 2) ** 2 + np.cos(de0) * np.cos(de1)
                                           * np.sin((ra1 - ra0) / 2) ** 2))) * 3600
    out = sub.assign(_r=np.round(sep, 2))
    return out[out["_r"] <= radius_arcsec].sort_values("_r", ignore_index=True)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    if cmd == "ensure":
        ensure(force="--force" in sys.argv)
    for n in BUILDERS:
        a = age_days(n)
        print(f"{n}: " + ("not downloaded" if a is None else
                          f"{os.path.getsize(_path(n)) / 1e6:.1f} MB, {a:.0f} days old"))
