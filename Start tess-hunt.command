#!/bin/bash
# Double-click this file (macOS Finder) to start the tess-hunt app in your browser.
# First start: creates a Python environment in .venv and installs the requirements
# (several minutes). Later starts reuse it and open the app in a few seconds.
# Close the Terminal window (or press Ctrl+C in it) to stop the app.

cd "$(dirname "$0")" || exit 1

fail() {
    echo
    echo "ERROR: $1"
    echo
    read -r -p "Press Return to close this window. " _
    exit 1
}

# 1. Python 3.10 or newer
PY=""
for c in python3.13 python3.12 python3.11 python3.10 python3; do
    if command -v "$c" >/dev/null 2>&1 && \
       "$c" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' 2>/dev/null; then
        PY="$c"
        break
    fi
done
[ -n "$PY" ] || fail "Python 3.10 or newer is needed. Install it from https://www.python.org/downloads/ and double-click this file again."

# 2. The environment (created once; reinstalled only when requirements.txt changes)
if [ ! -x .venv/bin/python ]; then
    echo "First start: creating the Python environment in .venv ..."
    "$PY" -m venv .venv || fail "Could not create the Python environment (.venv)."
fi
REQ_HASH=$(shasum requirements.txt 2>/dev/null || sha1sum requirements.txt)
REQ_HASH=${REQ_HASH%% *}
if [ "$(cat .venv/.tesshunt-requirements 2>/dev/null)" != "$REQ_HASH" ]; then
    echo "Installing the requirements (this takes a few minutes the first time) ..."
    .venv/bin/python -m pip install --upgrade pip >/dev/null
    .venv/bin/python -m pip install -r requirements.txt || fail "Installing the requirements failed (see the messages above). Check your internet connection and try again."
    echo "$REQ_HASH" > .venv/.tesshunt-requirements
fi

# 3. The app (Streamlit opens it in your default browser)
echo
echo "Starting tess-hunt. It opens in your browser; if not, go to http://localhost:8501"
echo "To stop it, close this window or press Ctrl+C."
echo
exec .venv/bin/python -m streamlit run app/main.py --browser.gatherUsageStats false
