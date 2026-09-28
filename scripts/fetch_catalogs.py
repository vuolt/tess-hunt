"""Download the bulk star catalogues (VSX, WDS) once, if missing or stale.

Run as the first step of run_sector.py; see tesshunt/catalogs.py."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tesshunt import catalogs  # noqa: E402

if __name__ == "__main__":
    catalogs.ensure(force="--force" in sys.argv)
