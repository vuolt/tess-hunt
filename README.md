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

## Phase 3: vetting

```
tesshunt/vetting.py      checks 1-7: shape (batman transit vs box/ramp/step/flare-decay, BIC),
                         duration vs P > 20 d, edge, TESScut difference imaging + centroid +
                         neighbours + in-pixel confirmation, SkyBoT asteroids, ExoFOP TOI/CTOI +
                         TESS EB catalogue, TRICERATOPS FPP
scripts/phase3_vet.py    resumable stages: lc | pixels | fpp | report
scripts/phase3_report.py funnel, validation table, shortlist, vetting sheets
results/phase3/          vetting_all.csv, shortlist.csv, validation.csv, summary.json
plots/phase3/            funnel.png, sheets/ (shortlist), validation/ (known TOIs)
```

```
python scripts/phase3_vet.py lc        # ~15 min: re-downloads each candidate's light curve
python scripts/phase3_vet.py pixels    # ~20 min: TESScut, SkyBoT, catalogues
python scripts/phase3_vet.py fpp       # ~1.5 h: TRICERATOPS (Gaia DR3 over HTTPS)
python scripts/phase3_vet.py report
```

Cutouts and light curves are reduced and discarded; cached per-candidate
results live in `work/phase3/` (git-ignored). TRICERATOPS needs
`setuptools<81` (its `pytransit` dependency imports `pkg_resources`).

## Phase 4: follow-up of the shortlist, and the whole pipeline in one command

```
tesshunt/net.py          polite networking: disk cache for every response, per-service rate
                         limits shared across processes, backoff on 429/5xx (see CLAUDE.md)
tesshunt/multisector.py  observed sectors (tess-point), best light curve per sector
                         (SPOC 2-min -> TESS-SPOC -> QLP -> own TESScut photometry; S3 first),
                         duotransit search, exclusion maps, allowed-period scans
tesshunt/binarity.py     Gaia DR3 RUWE / NSS / image doubling / RV scatter, co-moving
                         companions, 1-px blends, El-Badry wide binaries, WDS
scripts/phase4.py        coverage | search | confirm | periods | binarity | rank | report | all
scripts/phase4_report.py submit / maybe / drop rules, CTOI summaries, follow-up sheets
scripts/phase4_injection_vetting.py   Phase 2 injections through all seven checks
scripts/run_sector.py    everything, for any sector, in one command
results/phase4/          followup.csv, ctoi_candidates.csv, ctoi_summaries.md,
                         injection_vetting*.csv/json, phase4_summary.json
plots/phase4/            sheets/ (vetting + follow-up page per candidate), periods/,
                         injection_vetting_recall.png
```

### Run any sector

```
python -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
python scripts/run_sector.py --sector 49              # select -> search -> vet -> follow up
python scripts/run_sector.py --sector 49 --from phase3-lc
python scripts/run_sector.py --sector 49 --skip injection-vetting
python scripts/run_sector.py --sector 49 --dry-run    # list the steps
```

Every step caches its work and skips what is done, so rerunning the command after
an interruption continues where it stopped. Outputs for sectors other than 48 go
to `sXXXX/` subfolders of `results/phaseN`, `plots/phaseN` and `work/phaseN`
(Sector 48 predates this and uses the top-level folders).

Rough cost for a sector of ~130k stars on 4 cores: selection a few minutes, the
search ~2 h when TESS-SPOC light curves are on the S3 mirror (S1-S81 at the time
of writing), vetting ~2-3 h (dominated by TRICERATOPS), follow-up ~30 min.
Sectors newer than the mirror fall back to MAST at 4 files/s, which is slow by
design.

### External services

All network access goes through `tesshunt/net.py`:

- bulk light curves come from the AWS S3 mirror (`stpubdata`: SPOC 2-min under
  `tess/public/tid`, TESS-SPOC and QLP under `mast/hlsp`), with MAST as fallback;
- every response is cached under `work/cache/` and never fetched twice (the
  bulk Phase 2 search is the exception: its ~50 GB of light curves per sector
  are streamed and deleted, and its resumable database stops re-downloads);
- small services (SkyBoT, ExoFOP, Gaia archive, VizieR, MAST catalogue and
  TESScut queries) are serialized across processes at <= ~1.7 requests/s;
- HTTP 429/5xx and dropped connections back off exponentially, and after 6
  failures the call raises `ServiceError` instead of retrying forever;
- catalogue tables (ExoFOP TOIs/CTOIs, TESS EB catalogue) are downloaded once
  in bulk.
