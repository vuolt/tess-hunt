# tess-hunt
Searching public TESS data for long-period exoplanets that transit only once, the ones automated pipelines tend to miss.

## Quick start (no programming needed)

1. **Get the project:** `git clone https://github.com/vuolt/tess-hunt.git`, or
   download it as a ZIP from GitHub and unzip it.
2. **Start the app:** double-click **`Start tess-hunt.command`** (macOS) or
   **`Start tess-hunt.bat`** (Windows). The first start installs everything
   (several minutes); after that it opens in your browser in seconds. You
   need Python 3.10 or newer.
3. **Run page:** pick one of two jobs.
   - **Search a new sector for planets.** Pick a TESS sector nobody has
     searched yet; the page lists the searched ones. It runs for 5–8 hours
     in the background.
   - **Check other people's candidates.** This looks for false positives
     among the TOIs and CTOIs that nobody has confirmed or ruled out yet.
     Each run checks the next batch (500 by default, about half an hour).
     Tick "Download today's lists" now and then to include new
     submissions.
4. **When a run finishes,** the page says in plain words what it found, or
   that nothing was found. It gives you a ready commit message and the
   commands to paste into a terminal (or what to type into GitHub Desktop).
   Committing matters even when nothing was found. The record of searched
   sectors and stars (`results/sectors_searched.csv`) and of checked
   candidates (`results/phase6/lc_checks.csv.gz`) is how the project, and
   anyone who uses it after you, avoids doing the same work twice.
5. **Candidates and Workflow pages:** everything about a candidate in plain
   English, and the steps to get it looked at by the Planet Hunters TESS
   community and submitted to ExoFOP.

If you don't have permission to push to this repository, push to your own
fork and open a pull request.

The same jobs from a terminal:

```
python scripts/run_sector.py --sector 22          # search a sector
python -m tesshunt.findings sector 22             # what it found + commit message
python scripts/run_fp_triage.py [--refresh]       # check 500 more TOIs/CTOIs
python -m tesshunt.findings fp                    # what it found + commit message
```

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
scripts/pht_plots.py     light curve of each submit / maybe candidate in the style of
                         Planet Hunters TESS (plots/phase4[/sXXXX]/pht/), for the forum
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
python scripts/run_sector.py --sector 48 --force      # a sector already searched, again
```

Every step caches its work and skips what is done, so rerunning the command after
an interruption continues where it stopped. Outputs for sectors other than 48 go
to `sXXXX/` subfolders of `results/phaseN`, `plots/phaseN` and `work/phaseN`
(Sector 48 predates this and uses the top-level folders).

### Sectors and stars already searched

`results/sectors_searched.csv` lists every sector searched so far, how far it
got (phase2 / phase3 / phase4) and its headline numbers. It is rebuilt from the
committed results at the end of every run (or with `python -m tesshunt.ledger`),
so it always matches them.

- **Finished sector:** `run_sector.py` prints what was found and stops.
  `--force` runs it again; `--from STEP` redoes later steps.
- **Stars:** `results/phase2/sXXXX_stars.csv.gz` lists every star searched
  in that sector. On a fresh clone, the search copies those stars and their
  dips into the local database instead of downloading their light curves
  again. Stars that failed with an error are retried. `--search-again` on
  `phase2.py search` turns this off.
- **Other sectors:** a star searched in one sector is still searched in
  another, because each sector is new data. The search prints how many
  such stars there are.
- **The app's Run page** shows the same list and warns before rerunning a
  finished sector.

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

## Phase 5: expert checks on the submit and maybe candidates

```
tesshunt/expert.py   aperture test (TESScut, PSF model of Gaia neighbours), GP + transit
                     refit (celerite2), Gaia dwarf/subgiant check, Gaia binarity extras,
                     Gaia variability + AAVSO VSX, density-prior MCMC fit
scripts/phase5.py    run [--sectors 48 21] | report
results/phase5/      candidates.md (plain-English write-up per candidate),
                     phase5_checks.csv, calibration.csv, ctoi_summaries.md
```

Each candidate gets a verdict: strong, plausible or doubtful. A doubtful
"submit" becomes "maybe" and a doubtful "maybe" is dropped. The checks are
calibrated on known planets in the same sectors (`calibration.csv`).

## Phase 6: likely false positives among existing TOIs and CTOIs

```
python scripts/run_fp_triage.py                  # all three steps, 500 new candidates
python scripts/phase6.py gaia [--refresh]        # TOI/CTOI tables (--refresh: today's), Gaia orbits
python scripts/phase6.py lc --all --limit 500    # odd/even, secondary, centroid on S3 light curves
python scripts/phase6.py report                  # validation, flags, results/phase6/summary.md
```

Every light-curve check is recorded in `results/phase6/lc_checks.csv.gz`
(committed). A new clone restores the checks from that record instead of
downloading the light curves again, and each run checks only candidates
that haven't been checked yet. The validation samples (about 300 each of
unresolved candidates, planets and known false positives) are drawn once
and kept, even after the tables are refreshed. The TOI/CTOI tables are
cached until `--refresh` downloads today's copies, once per run.

Every check is first run on confirmed planets and known false positives from
the same tables. A check variant that flags more than 5 % of known planets is
not used (2 % preferred). Results: `results/phase6/`.

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

  Next to it is the light curve drawn the way Planet Hunters TESS shows it
  (whole sector, then a zoom on the dip), then the follow-up sheets and vetting plots. Words with a dotted
  underline show an explanation when you hover over them.
- **Workflow.** Four steps per candidate, which can only be done in order:
  - **a. Get feedback.** A ready-to-paste Planet Hunters TESS forum post, the
    forum-style light-curve plot (PNG) and the key plots as a zip to attach, and a link to the forum. After posting, record
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
- **Run.** Two background jobs, with a live progress bar and a stop button.
  You can close the app and come back.
  - **Search a sector** (`scripts/run_sector.py`): the whole pipeline,
    including the Phase 5 expert checks. It takes several hours. Pressing
    *Start / resume* on an interrupted sector continues where it stopped,
    and a finished sector asks before running again.
  - **Check other people's candidates** (`scripts/run_fp_triage.py`): the
    Phase 6 false-positive checks on the next batch of TOIs/CTOIs not yet
    checked, optionally after downloading today's lists.

  When a run finishes, the page lists what was found (candidates worth
  submitting, strong ones held back, or new likely false positives), or
  says nothing was found. It then gives a ready commit message with the
  commands to save the results. The app never commits anything itself.
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
