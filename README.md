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
scripts/phase4.py        coverage | search | confirm | periods | binarity | validate | rank |
                         report | all  (validate: rules D1/D5 on known planets)
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

## The app: explore and follow up candidates without the command line

A local web app for someone with no astronomy background. It reads the
pipeline's results and walks you through getting feedback on a candidate and
preparing an ExoFOP submission.

### Start it with a double-click

- **macOS:** double-click **`Start tess-hunt.command`** in Finder.
- **Windows:** double-click **`Start tess-hunt.bat`**.

A terminal window opens, and then the app opens in your browser at
http://localhost:8501.

- **First start:** it creates a Python environment in `.venv` and installs the
  requirements, which takes several minutes.
- **Later starts:** the app opens in a few seconds. The requirements are
  reinstalled only when `requirements.txt` changes.
- **To stop the app:** close the terminal window.
- **Python:** 3.10 or newer is needed. If it's missing, the window says so and
  links to python.org.
- **macOS, if the repository was downloaded as a ZIP:** macOS may refuse to run
  a downloaded script the first time. Right-click the file, choose **Open**, then
  confirm. After that, a double-click works. A `git clone` doesn't have this
  problem.
- **macOS, "permission denied":** run `chmod +x "Start tess-hunt.command"` once.

### Or from a terminal

```
python -m venv .venv && . .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app/main.py
```

### Pages

- **Candidates.** One card per shortlisted candidate, in plain words:
  - what kind of star it is, and how far away (if the Phase 5 checks have run);
  - how deep and how long the dip is, and how big the object would be compared
    with Earth, Neptune and Jupiter;
  - what orbits are still possible;
  - its status (submit / maybe / drop) and the reason.

  The follow-up sheets and vetting plots are next to it. Words with a dotted
  underline show an explanation when you hover over them.
- **Workflow.** Four steps per candidate, which can only be done in order:
  - **a. Get feedback.** A ready-to-paste Planet Hunters TESS forum post, the key
    plots as a zip to attach, and a link to the forum. After posting, record
    the post link and whether the feedback was positive, negative or unclear,
    with notes.
  - **b. Prepare submission.** Unlocked only by positive feedback, or by an
    "override anyway" with a typed reason. It gives:
    - the ExoFOP CTOI fields to copy (TIC, epoch, period constraints, depth,
      duration, size, vetting summary);
    - a file in ExoFOP's planet-parameter bulk-upload format;
    - links to ExoFOP's guidelines, its upload-access request page and its
      bulk-upload page.

    **ExoFOP accepts community candidates only once they are published in a
    refereed journal**, and uploading needs approved access. The app says so,
    and asks for the paper URL and your ExoFOP data tag.
  - **c. Record submission.** The date you submitted, and the CTOI number once
    ExoFOP assigns one (you can add it later).
  - **d. Follow-up.** Predicted future transits, where a period is known
    (for example TIC 95747180), and whether TESS will observe the star again,
    from TESS's published pointing plan.
- **Results.** For each searched sector, the funnel from stars searched to
  dips, candidates, vetted, the new shortlist, and "worth submitting". One
  sentence per step says what it removed.
- **Run.** Choose a sector and start the whole pipeline
  (`scripts/run_sector.py`) in the background. There is a live progress bar and
  a stop button. A sector takes several hours; you can close the app and come
  back. Pressing *Start / resume* on an interrupted sector continues where it
  stopped, because every step keeps its finished work.
- **Glossary.** Plain-English explanations of TIC, SNR, FPP, duotransit, CTOI,
  ExoFOP, TFOP and more.

### What the app does and does not do

- It **never logs in to or submits anything to an external website.** You post
  on the forum and upload to ExoFOP yourself, using the text and files it
  prepares.
- It **only reads** `results/` and `plots/`, and never changes results or
  pipeline code. A test checks this (`tests/test_app_ui.py`).
- It writes only:
  - `app/workflow_status.json`, your workflow progress, which you can commit if
    you want to keep it in the repository;
  - `work/app_runs/`, the log of a search started from the Run page, in the
    git-ignored scratch area;
  - the pipeline's own outputs, but only while a search you started is running.
- Tests: `pytest tests/test_app_workflow.py tests/test_app_ui.py`. They cover
  the step order, the unlock and override rules, saving and loading the status
  file, run-log parsing, start and stop, the ExoFOP file format, and a
  click-through of the whole workflow.
