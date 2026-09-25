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

## Phase 2: full-sector search on FFI light curves

```
tesshunt/
  ffi.py           TESS-SPOC FFI light curves by direct MAST URL
  tic.py           TIC lookup; dwarf, Tmag < 13 selection
  variability.py   variability timescale -> detrending window, long-limited / fast-variable flags
  survey.py        per-star download -> analyse -> delete; injections; resumable SQLite store
  vetting.py       common-mode, repeating/single, partial, implied radius, tiers, TOI match
scripts/
  phase2.py          select | search   (resumable: rerun to continue; --retry-errors)
  phase2_report.py   CSVs, rankings, summary JSON, figures, top-100 diagnostics
results/phase2/    stars (gz), dips ranked by SNR, candidates, injections, TOI recall, summary
plots/phase2/      summary figures, top100/, top_candidates_sheet.png, validation plot
```

```
export OMP_NUM_THREADS=1                  # one thread per worker process
python scripts/phase2.py select           # ~4 min: target list + TIC -> work/s0048_sample.csv
python scripts/phase2.py search           # ~2 h on 4 cores for 128k stars; safe to interrupt
python scripts/phase2_report.py           # ~5 min, re-downloads ~150 stars for plots
```

Working state (TIC table, sample, `work/s0048.sqlite`, temporary FITS files)
lives in `work/`, which is git-ignored. Each light curve is deleted as soon as
its star is analysed.
