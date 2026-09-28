# Be a good citizen to external data services:

For bulk TESS downloads, prefer the AWS S3 mirror (stpubdata) over MAST's web endpoints.
Cache every response in work/ and never re-download what's already cached.
For small services (SkyBoT, ExoFOP, Gaia archive, VizieR, TRILEGAL, MAST catalogue queries), make at most 1–2 requests per second in total, with no parallel requests.
On HTTP 429 or 503, back off exponentially and retry; after repeated failures, stop and report instead of retrying endlessly.
Download catalogue tables once in bulk (e.g. the full TOI/CTOI lists) instead of querying per star.
If a service answers 429 (too many requests) or sends Retry-After, treat it as a signal to slow down for the rest of the run, not just to retry that one call.

# Use only free, public services:

Use only services that are free and meant for public use: MAST, the AWS S3 mirror stpubdata (free through the AWS Open Data programme; access it anonymously), the Gaia archive, CDS (VizieR, SkyBoT), NASA Exoplanet Archive, ExoFOP and similar.
Never use paid services, paid APIs, paid cloud compute, requester-pays buckets or anyone's credentials or API keys.
If a paid or private service appears reachable by accident (a misconfigured bucket, an open endpoint, a leaked key), do not use it; stop and tell the user instead.
