# Be a good citizen to external data services:

For bulk TESS downloads, prefer the AWS S3 mirror (stpubdata) over MAST's web endpoints.
Cache every response in work/ and never re-download what's already cached.
For small services (SkyBoT, ExoFOP, Gaia archive, VizieR, TRILEGAL, MAST catalogue queries), make at most 1–2 requests per second in total, with no parallel requests.
On HTTP 429 or 503, back off exponentially and retry; after repeated failures, stop and report instead of retrying endlessly.
Download catalogue tables once in bulk (e.g. the full TOI/CTOI lists) instead of querying per star.
