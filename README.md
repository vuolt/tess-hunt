# tess-hunt
Searching public TESS data for long-period exoplanets that transit only once, the ones automated pipelines tend to miss.

## Phase 1: single-event detection

```
tesshunt/        pipeline package
  lightcurves.py   download / load SPOC 2-min light curves by TIC ID (lightkurve)
  detrend.py       segment-wise time-windowed biweight trend
  detect.py        single-event box search with empirical (red-noise) SNR
  pipeline.py      two-pass detrend + search per sector
  injection.py     trapezoid injection + recovery scoring
  plotting.py      diagnostic figures
scripts/
  run_known_target.py    blind search on TOI-2180 b, a known ~260 d planet
  injection_recovery.py  injection-recovery grid on real TOI-2180 light curves
tests/             synthetic-data unit tests (no network)
plots/, results/   committed outputs
```

Setup and run:

```
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
python scripts/run_known_target.py        # downloads ~75 MB into data/ (git-ignored)
python scripts/injection_recovery.py      # ~1.5 min on 4 cores
python -m pytest tests
```

See [RESULTS.md](RESULTS.md) for what it found.
