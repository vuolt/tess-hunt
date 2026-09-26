"""Polite, cached network access for every external service (see CLAUDE.md).

* Every response is cached on disk under work/cache/ and never re-fetched.
* "Small" services (SkyBoT, ExoFOP, Gaia archive, VizieR, MAST catalogue and
  TESScut queries) are serialized *across processes* with a lock file and
  spaced by at least SMALL_INTERVAL seconds, i.e. at most ~1.7 requests/s in
  total with no parallel requests.
* HTTP 429/503 (and dropped connections) are retried with exponential
  backoff; after MAX_RETRIES the call raises ServiceError so the caller can
  stop and report instead of hammering the service.
* Bulk files come from the AWS S3 mirror (stpubdata) first, MAST second.
"""

from __future__ import annotations

import contextlib
import fcntl
import hashlib
import os
import pickle
import random
import time
import urllib.error
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.environ.get("TESSHUNT_CACHE", os.path.join(ROOT, "work", "cache"))

SMALL_INTERVAL = 0.6      # s between requests to one small service (~1.7/s)
MAX_RETRIES = 6
TIMEOUT = 120

# service -> minimum spacing (s); None = bulk service (S3), not serialized
SERVICES = {
    "skybot": SMALL_INTERVAL, "exofop": SMALL_INTERVAL, "gaia": SMALL_INTERVAL,
    "vizier": SMALL_INTERVAL, "mast_catalog": SMALL_INTERVAL, "tesscut": SMALL_INTERVAL,
    "mast_api": SMALL_INTERVAL, "mast_files": 0.25, "exoplanet_archive": SMALL_INTERVAL,
    "villanova": SMALL_INTERVAL, "irsa": SMALL_INTERVAL,
    "s3": None,
}

S3 = "https://stpubdata.s3.amazonaws.com"
MAST_ARCHIVE = "https://archive.stsci.edu"


class ServiceError(RuntimeError):
    """A service kept failing after MAX_RETRIES backoff attempts."""


def _key(*parts) -> str:
    h = hashlib.sha1()
    for p in parts:
        h.update(repr(p).encode())
    return h.hexdigest()


def _cache_path(service: str, key: str, ext: str) -> str:
    d = os.path.join(CACHE, service, key[:2])
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, key + ext)


@contextlib.contextmanager
def service_slot(service: str):
    """Exclusive, rate-limited access to a small service across processes."""
    interval = SERVICES.get(service, SMALL_INTERVAL)
    if interval is None:
        yield
        return
    os.makedirs(CACHE, exist_ok=True)
    lock_path = os.path.join(CACHE, f".{service}.lock")
    with open(lock_path, "a+") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            fh.seek(0)
            txt = fh.read().strip()
            last = float(txt) if txt else 0.0
            wait = last + interval - time.time()
            if wait > 0:
                time.sleep(wait)
            yield
        finally:
            fh.seek(0)
            fh.truncate()
            fh.write(f"{time.time():.3f}")
            fh.flush()
            fcntl.flock(fh, fcntl.LOCK_UN)


def _retrying(service, fn):
    """Run fn() inside the service slot, backing off on 429/503/connection errors."""
    for attempt in range(MAX_RETRIES):
        try:
            with service_slot(service):
                return fn()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise FileNotFoundError(getattr(e, "url", "")) from e
            if e.code not in (429, 500, 502, 503, 504):
                raise
            err = e
        except FileNotFoundError:
            raise
        except (urllib.error.URLError, ConnectionError, TimeoutError, OSError) as e:
            err = e
        time.sleep(min(2 ** (attempt + 1), 60) * (1 + 0.25 * random.random()))
    raise ServiceError(f"{service}: giving up after {MAX_RETRIES} attempts: {err}")


def get(url: str, service: str, data: dict | None = None, cache: bool = True,
        timeout: int = TIMEOUT) -> bytes:
    """GET (or POST form ``data``) with disk cache, rate limit and backoff.

    404s are cached too (as FileNotFoundError), so a missing product is only
    asked for once."""
    body = urllib.parse.urlencode(data).encode() if data else None
    k = _key(url, body)
    path, miss = _cache_path(service, k, ".bin"), _cache_path(service, k, ".404")
    if cache and os.path.exists(path):
        with open(path, "rb") as fh:
            return fh.read()
    if cache and os.path.exists(miss):
        raise FileNotFoundError(url)

    def fetch():
        req = urllib.request.Request(url, data=body,
                                     headers={"User-Agent": "tess-hunt (research; cached)"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    try:
        out = _retrying(service, fetch)
    except FileNotFoundError:
        if cache:
            open(miss, "w").close()
        raise
    if cache:
        tmp = path + ".part"
        with open(tmp, "wb") as fh:
            fh.write(out)
        os.replace(tmp, path)
    return out


def get_text(url, service, data=None, cache=True) -> str:
    return get(url, service, data, cache).decode()


def download(urls, service_for_url, dest: str | None = None, cache: bool = True) -> str:
    """Fetch the first available of ``urls`` to a local file and return its path.

    ``service_for_url(url)`` names the service for each URL. With cache=True
    the file lives in the cache (kept); with cache=False it is written to
    ``dest`` for the caller to delete (bulk sector processing)."""
    last = None
    for url in urls:
        service = service_for_url(url)
        k = _key(url)
        cpath = _cache_path("files", k, os.path.splitext(url)[1] or ".dat")
        miss = cpath + ".404"
        if cache and os.path.exists(cpath):
            return cpath
        if os.path.exists(miss):
            last = FileNotFoundError(url)
            continue
        try:
            out = get(url, service, cache=False)
        except FileNotFoundError as e:
            open(miss, "w").close()
            last = e
            continue
        path = cpath if cache else dest
        tmp = path + ".part"
        with open(tmp, "wb") as fh:
            fh.write(out)
        os.replace(tmp, path)
        return path
    raise last if last else FileNotFoundError(urls)


def service_of(url: str) -> str:
    host = urllib.parse.urlparse(url).netloc
    if "amazonaws.com" in host:
        return "s3"
    if "tesscut" in url:
        return "tesscut"
    if host == "archive.stsci.edu":
        return "mast_files"
    if "mast.stsci.edu" in host:
        return "mast_api"
    if "exofop" in host:
        return "exofop"
    if "imcce" in host:
        return "skybot"
    if "esac.esa.int" in host:
        return "gaia"
    if "irsa.ipac.caltech.edu" in host:
        return "irsa"
    if "cds" in host:
        return "vizier"
    if "exoplanetarchive" in host:
        return "exoplanet_archive"
    return "mast_files"


def cached_call(service: str, key_parts, fn, cache: bool = True):
    """Cache the (picklable) result of a library call such as an astroquery
    query, and run it inside the service slot with backoff on failure."""
    k = _key(*key_parts)
    path = _cache_path(service, k, ".pkl")
    if cache and os.path.exists(path):
        with open(path, "rb") as fh:
            return pickle.load(fh)
    for attempt in range(MAX_RETRIES):
        try:
            with service_slot(service):
                out = fn()
            break
        except Exception as e:  # noqa: BLE001  (astroquery raises many types)
            msg = str(e)
            transient = any(s in msg for s in ("429", "503", "502", "504", "timed out",
                                               "Connection", "reset", "Temporary"))
            if not transient or attempt == MAX_RETRIES - 1:
                raise ServiceError(f"{service}: {type(e).__name__}: {e}") from e
            time.sleep(min(2 ** (attempt + 1), 60) * (1 + 0.25 * random.random()))
    if cache:
        tmp = path + ".part"
        with open(tmp, "wb") as fh:
            pickle.dump(out, fh)
        os.replace(tmp, path)
    return out


def s3_list(prefix: str, delimiter: str | None = None, cache: bool = True) -> list[str]:
    """Keys (or common prefixes, with a delimiter) under an S3 prefix of stpubdata."""
    import re
    q = {"list-type": 2, "prefix": prefix, "max-keys": 1000}
    if delimiter:
        q["delimiter"] = delimiter
    out, token = [], None
    while True:
        if token:
            q["continuation-token"] = token
        xml = get_text(f"{S3}/?{urllib.parse.urlencode(q)}", "s3", cache=cache)
        tag = "Prefix" if delimiter else "Key"
        found = re.findall(rf"<{tag}>([^<]+)</{tag}>", xml)
        out += [f for f in found if f != prefix]
        m = re.search(r"<NextContinuationToken>([^<]+)</", xml)
        if not m:
            return out
        token = m.group(1)
