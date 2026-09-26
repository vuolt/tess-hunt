@echo off
rem Double-click this file (Windows Explorer) to start the tess-hunt app in your browser.
rem First start: creates a Python environment in .venv and installs the requirements
rem (several minutes). Later starts reuse it. Close this window to stop the app.

cd /d "%~dp0"

where py >nul 2>nul && (set "PY=py -3") || (set "PY=python")
%PY% -c "import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)" >nul 2>nul
if errorlevel 1 (
    echo Python 3.10 or newer is needed. Install it from https://www.python.org/downloads/
    echo ^(tick "Add python.exe to PATH"^) and double-click this file again.
    pause
    exit /b 1
)

if not exist ".venv\Scripts\python.exe" (
    echo First start: creating the Python environment in .venv ...
    %PY% -m venv .venv || (echo Could not create .venv & pause & exit /b 1)
)

for /f "skip=1 tokens=* delims=" %%h in ('certutil -hashfile requirements.txt SHA1') do if not defined REQHASH set "REQHASH=%%h"
set /p OLDHASH=<".venv\.tesshunt-requirements" 2>nul
if not "%REQHASH%"=="%OLDHASH%" (
    echo Installing the requirements ^(this takes a few minutes the first time^) ...
    ".venv\Scripts\python.exe" -m pip install --upgrade pip >nul
    ".venv\Scripts\python.exe" -m pip install -r requirements.txt || (echo Installing the requirements failed. & pause & exit /b 1)
    > ".venv\.tesshunt-requirements" echo %REQHASH%
)

echo.
echo Starting tess-hunt. It opens in your browser; if not, go to http://localhost:8501
echo To stop it, close this window.
echo.
".venv\Scripts\python.exe" -m streamlit run app\main.py --browser.gatherUsageStats false
pause
