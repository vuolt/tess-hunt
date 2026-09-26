# Results

- [Phase 1](#phase-1-single-transit-detection-on-one-star-2-min-data): the pipeline
  on one star with a known long-period planet, plus injection-recovery.
- [Phase 2](#phase-2-full-sector-search-tess-spoc-ffi-sector-48): a full-sector
  search of 128,455 dwarfs in TESS-SPOC FFI light curves.
- [Phase 3](#phase-3-vetting-the-sector-48-candidates): seven vetting checks
  cut the 2,627 single-dip candidates to a hand-reviewable shortlist of
  **28 new candidates** (plus 26 already-known TOIs/CTOIs).

## Phase 1: single-transit detection on one star (2-min data)

Everything here can be reproduced with `scripts/run_known_target.py` and
`scripts/injection_recovery.py` (see README). Raw FITS files go to `data/`,
which git ignores.

### Pipeline

1. **Data.** SPOC 2-min PDCSAP light curves from MAST, queried by TIC ID,
   using the default quality mask. Each sector is normalized to its median.
2. **Detrending.** A time-windowed Tukey biweight (3 d window, c = 5),
   computed separately for each continuous data segment. Pass 1 searches;
   any event found (±1 box duration) is then masked out of the trend fit and
   the sector is detrended and searched again. Upward outliers above 4σ are
   clipped. Downward points are never clipped.
3. **Detection.** A single-event box search on a 10-min grid with 10 trial
   durations from 1 to 24 h. The noise at each duration is measured
   empirically as the robust spread of box depths across the sector, so red
   noise is included. A box counts only if it is ≥ 60 % covered by data and
   has ≥ 50 % coverage in a baseline window on each side (width
   max(duration, 6 h)). Threshold: **SNR ≥ 7**.

### 1. Known planet: TOI-2180 b

TOI-2180 b (TIC 298663873, T = 8.6) is a ~260 d warm Jupiter. Published
ephemeris from Dalba et al. 2024, via the NASA Exoplanet Archive:
P = 260.17 d, T0 = BTJD 1830.756, T14 ≈ 24.1 h, depth ≈ 0.47 %. The search
covered all 40 SPOC 2-min sectors (S14–S86) blind. The ephemeris was used
only afterwards, to score the results.

| Published transit (BTJD) | Sector | Detected T0 | Offset | Box depth | Box dur. | SNR |
|---|---|---|---|---|---|---|
| 1830.76 | 19 | 1830.785 | +0.7 h | 4311 ppm | 16 h | 40.6 |
| 2611.27 | 48 | 2611.259 | −0.2 h | 3740 ppm | 16 h | 55.7 |
| 2871.44 | 57 | 2871.487 | +1.2 h | 4676 ppm | 16 h | 83.9 |

**All 3 transits that fall in the data were recovered, each as the
strongest event in its sector.** The best box (16 h) is shorter than T14
because a limb-darkened, U-shaped transit concentrates its SNR in the core.
The S48 depth is low because the transit sits at the start of a data
segment, where the trend has least support. Plots:
`plots/toi2180_overview.png` and `plots/toi2180_s{19,48,57}_event.png`.

Two other events cross the threshold in the remaining 37 sectors. Both
look like systematics or stellar noise, not planets:

- **S75, BTJD 3341.14** (SNR 11.1, 1.8 ppt, 4 h). A drop into the end of a
  data segment just before a gap: a classic downlink or scattered-light
  artifact (`plots/toi2180_s75_event.png`).
- **S81, BTJD 3514.12** (SNR 10.2, 190 ppm, 24 h). A shallow, day-long
  depression, probably low-frequency stellar or instrumental noise
  (`plots/toi2180_s81_event.png`).

The empirical false-alarm rate is therefore **~2 per 37 sectors (≈ 0.05
per star-sector)** at SNR 7, on a bright, quiet star. That is too many for a
blind search over thousands of stars without vetting (Phase 2).

### 2. Injection-recovery

**Setup.** Trapezoid transits (10 % ingress and egress) were injected
multiplicatively into the raw PDCSAP flux of the 37 TOI-2180 sectors with
no known transit, and the full two-pass pipeline was run on each. Each
cell of 8 depths × 6 durations got 60 trials (2,880 total, seed 42). Each
trial used a random sector and a mid-time placed on a real cadence, so some
injections land next to gaps, just as real transits would. An injection
counts as **recovered** if an event with SNR ≥ 7 is found within half the
injected duration of the injected T0.

Recovery fraction (`plots/injection_recovery_grid.png`):

| depth \ duration | 1 h | 2 h | 4 h | 8 h | 16 h | 24 h |
|---|---|---|---|---|---|---|
| 50 ppm | 0% | 0% | 0% | 0% | 2% | 0% |
| 100 ppm | 0% | 0% | 0% | 0% | 0% | 0% |
| 200 ppm | 0% | 0% | 0% | 2% | 2% | 0% |
| 300 ppm | 0% | 0% | 0% | 7% | 18% | 15% |
| 500 ppm | 2% | 2% | 47% | 65% | 70% | 58% |
| 1000 ppm | 73% | 95% | 87% | 97% | 90% | 92% |
| 2000 ppm | 92% | 97% | 100% | 100% | 100% | 100% |
| 5000 ppm | 93% | 98% | 100% | 98% | 100% | 100% |

Recovery as a function of *expected* SNR (injected mean box depth divided by
the empirical noise for that sector and duration) carries over to other
stars better than depth does (`plots/injection_recovery_curves.png`):

| expected SNR | 5–6 | 6–7 | 7–8 | 8–9 | 9–10 | 10–12 | 12–15 | 15–20 | ≥ 20 |
|---|---|---|---|---|---|---|---|---|---|
| recovered | 3% | 19% | 33% | 60% | 59% | 80% | 84% | 97% | 98–99% |

**Takeaways**

- For this T = 8.6 star, recovery is ≥ 90 % for single transits of
  **≥ 1 ppt lasting ≥ 2 h** and ≈ 50–70 % for **500 ppm lasting 4–24 h**.
  In practice that means Jupiters and Saturns around Sun-like stars, and
  Neptunes only when they are long-duration.
- 50 % recovery sits at expected SNR ≈ 8, and ≥ 95 % needs SNR ≳ 15. The
  transition is broader than for pure white noise because both the injected
  signal and the noise estimate vary from sector to sector.
- The ~2–7 % of high-SNR injections that are still missed land within a few
  hours of a data gap. Gaps include the 5 h downlink gaps in recent sectors
  and quality-flag dropouts, which are shorter than the 0.5 d
  segment-split threshold. The two-sided baseline rule rejects those boxes
  by design, since it is what suppresses edge systematics. 1 h transits
  suffer most because their 6 h baseline windows are long compared with
  the transit.
- Recovered events are accurate: the median |ΔT0| is 0.13 h (3 % of the
  duration), the median depth ratio is 0.96 (16–84 %: 0.88–1.04), and the
  box duration comes out slightly short (e.g. 6 h for an 8 h trapezoid).
- Ablation: requiring 75 % baseline coverage instead of 50 % gave the same
  false-alarm count (2) but lower recovery near gaps (e.g. 2000 ppm × 4 h:
  95 % vs 100 %), so 50 % is the default.

### Limitations and next steps

- **One bright, quiet star.** The absolute depth limits apply to T ≈ 8.6.
  For fainter stars use the expected-SNR curve. Injecting into a sample of
  stars across magnitude and variability is the obvious next step.
- **Fast variability.** A location-only biweight with a 3 d window cannot
  follow variability faster than ~10 d. In a synthetic test, a 5.3 d, 0.2 %
  sinusoid left ~900 ppm residuals. Rapid rotators need a shorter window,
  which costs sensitivity to long transits, or a GP- or spline-based
  detrender.
- **Gaps.** Transits that are partly in a gap are mostly lost. A dedicated
  partial-transit search would fix that.
- **Vetting.** There is no centroid, odd/even, background, or
  momentum-dump check, and no shape test (box vs. transit vs. ramp). The
  two null events above would need to be rejected by exactly these.
- Only 2-min SPOC data. FFI light curves (TESS-SPOC / QLP) would add far
  more stars.

## Phase 2: full-sector search (TESS-SPOC FFI, Sector 48)

Reproduce with `scripts/phase2.py select`, `scripts/phase2.py search` and
`scripts/phase2_report.py` (see README). Every number below comes from
`results/phase2/s0048_summary.json`.

### Sample and processing

| step | stars |
|---|---|
| TESS-SPOC FFI targets in S48 | 159,585 |
| TIC luminosity class DWARF | 156,647 |
| … and Tmag < 13 (objType STAR, not artifact/duplicate) | **128,486** |
| searched | **128,455** |
| skipped: < 1,000 good cadences | 30 |
| skipped: no light-curve file on MAST | 1 |
| errors | 0 |

- **Data.** S48 (BTJD 2610–2636, Jan–Feb 2022) was chosen because TOI-2180 b
  transits in it, which gives an end-to-end check against a known
  long-period planet. The light curves are TESS-SPOC 10-min PDCSAP, fetched
  by direct URL and filtered with lightkurve's default quality mask.
- **Batching and resuming.** Stars are processed in batches of 1,000 with
  10 worker processes. Each star is downloaded, analysed, and its file
  deleted, and its results are written to SQLite in a single transaction. A
  restart skips every star already recorded. This was exercised once
  mid-run: after 88 stars the search was restarted and resumed with none
  lost. The run took about 2 h of wall time on 4 cores. Peak disk use was
  the ~56 MB results database; no light curve is kept.
- **Noise.** The median 10-min point-to-point scatter is 1,418 ppm:
  241 ppm for T < 9, 389 for T 9–10, 633 for T 10–11, 1,091 for T 11–12,
  and 1,917 for T 12–13. The median red-noise factor at 4 h (empirical box
  noise / white) is 1.09.

### Variability and adaptive detrending

For each star, a preliminary 3 d-window search masks its three strongest dips
first, so that a deep transit is not mistaken for variability. A
Lomb-Scargle periodogram then gives the dominant period P and
semi-amplitude A. The detrending window W is the longest of
{3, 2.5, 2, 1.5, 1.25, 1, 0.75, 0.5} d whose predicted residual,
A·|1 − sinc(W/P)|, stays below max(3 σ₂₄ₕ, 500 ppm). The longest duration
searched is W/3.

| | stars | fraction |
|---|---|---|
| significant periodicity (FAP < 10⁻³) | 55,274 | 43.0 % |
| full 3 d window | 120,430 | 93.8 % |
| **long-limited**: W < 3 d, day-long transits not searched | **8,025** | 6.2 % |
| **fast variable**: no allowed window removes the variability | **3,965** | 3.1 % |
| either flag | 10,560 | 8.2 % |

The shortened windows break down as 2.5 d: 1,537; 2 d: 1,492; 1.5 d: 1,469;
1.25 d: 646; 1 d: 536; 0.75 d: 486; 0.5 d: 1,859. Flagged stars are
spotted rotators with P ≈ 1–13 d and eclipsing binaries or pulsators
(`plots/phase2/variability_windows.png`). The flagged fraction is higher
for bright stars, because the 500 ppm tolerance is fixed while their noise
is lower, and because hot pulsators are brighter.

Injection-recovery (below) shows what the long-limited flag means in
practice. Only 9 flagged stars were in the injection set, so these numbers
are rough. At expected SNR > 15:

- transits within the shortened search range: 94 % recovered (48/51);
- transits longer than W/3: 45 % recovered (10/22).

So "day-long transits unrecoverable" is closer to "sensitivity to day-long
transits roughly halved". Longer transits are sometimes still caught by a
shorter box on their core.

### Detections

**10,912 dips above SNR 7 on 5,078 stars (4.0 % of those searched).** Each
dip is recorded in `results/phase2/s0048_dips.csv` with its TIC, T0,
duration, depth, SNR and rank, and vetted automatically
(`tesshunt/vetting.py`):

| category | dips | meaning |
|---|---|---|
| repeating | 5,583 | the star's strongest dip has a partner within ×2 in duration and depth: eclipsing binary, short-period planet, or periodic variability |
| common_mode | 1,738 | stars within 2° dip within ±0.25 d far more often than that patch's own rate predicts (Poisson p < 10⁻³): detector systematics |
| secondary | 906 | weaker, dissimilar dips on a star whose strongest dip is single |
| **single** | **2,560** | the strongest non-systematic dip on a star with no repeat |
| **single_partial** | **125** | as `single`, but the box reaches within 1 h of a data gap |

- **Candidates** are `single` or `single_partial` on stars that are not
  fast-variable: **2,627 candidates on 2,627 stars (2.0 %)**, listed in
  `s0048_candidates.csv`.
  - 2,459 are planet-sized: √depth · R★(TIC) ≤ 2 R_J.
  - 1,464 of those have SNR > 10, and 581 have SNR > 15.
  - Their box durations are spread fairly evenly from 1 to 24 h.
- **Tier A** (planet-sized, not partial, outside the sector-wide pile-up
  windows): **991 candidates**, of which 548 have SNR > 10 and 256 have
  SNR > 15. 24 of them match known TOIs. The top 50 are shown in
  `plots/phase2/top_candidates_sheet.png`.
- **The candidate list is still dominated by systematics.** A real
  transit population would be flat in time. Instead, **55 % of candidates
  (1,443) fall in the 4.25 d of sector-wide pile-up windows**
  (`plots/phase2/detection_times.png`), which are only ~19 % of the 22 d of
  data: the orbit starts, the orbit ends, and one mid-orbit-1 episode. The
  local common-mode test removes 1,681 dips in those windows, but not all
  of them. The contact sheet also still contains ramps and slopes from
  imperfect detrending near segment edges. Tier A excludes the pile-up
  windows, but it needs shape and centroid vetting before anyone looks at
  it as a planet list.
- **Top 100 by SNR** (`plots/phase2/top100/`, one diagnostic per dip): all
  100 are eclipsing-binary depths (≥ 12.9 %, implied companion > 2 R_J) on
  40 stars. By category: 82 repeating, 15 single (lone eclipses of
  longer-period binaries), 2 common-mode, 1 secondary.

### Known planets as an external check

**TOI-2180 b** (P = 260 d) was recovered blind at T0 = 2611.433 BTJD
(TOI-table T0 2611.446), with SNR 27, a 24 h box and depth 2.0 ppt. The
depth is underestimated (true 4.7 ppt) because the transit starts on the
orbit-1 ramp. It is labelled `single_partial`, tier B, since it is next to
the gap and inside an orbit-start pile-up window
(`plots/phase2/validation_tic298663873.png`). Its window stayed at 3 d
because the transit was masked before variability was measured. Without
that step the transit's own low-frequency power shrank the window to
1.25 d and the transit was filtered out. That bug was caught on this star
before the full run.

**TOI recall.** The sample contains 381 TOI hosts, and 363 TOIs have a
predicted transit inside the data windows (NASA Exoplanet Archive
ephemerides; data windows approximated from where dips occur). **129
(36 %) were recovered.** By TOI depth:

| TOI depth | TOIs | recovered |
|---|---|---|
| < 1 ppt | 150 | 7 % |
| 1–3 ppt | 95 | 16 % |
| 3–10 ppt | 77 | 82 % |
| > 10 ppt | 41 | 98 % |

This matches the injection sensitivity below. Of the 62 TOIs with exactly
one predicted transit in the sector, 28 were recovered, including 88 %
(15/17) of those deeper than 3 ppt. 24 of the 28 were labelled `single` or
`single_partial`: for example TOI-5718.01 (P = 717 d), TOI-1859.01 (63 d),
TOI-5722.01 (61 d) and TOI-2180.01 (260 d). The remaining 4 were labelled
repeating (2), secondary (1) and common-mode (1). Most multi-transit TOIs
came out as `repeating` (97 of the 129 recoveries).

### Sensitivity across real stars (injection-recovery)

**Setup.** 200 stars were drawn at random (seed 2024) from the sample, with
20 injected single transits each (4,000 trials). Each trial used:

- a trapezoid shape with 10 % ingress;
- a duration drawn from {1, 2, 4, 8, 16, 24} h;
- a depth set so the expected SNR (mean box depth over that star's
  empirical box noise at that duration) is log-uniform from 3 to 60;
- a mid-time on a random real cadence.

The full pipeline was re-run on every injected light curve, including the
variability measurement and window choice. An injection counts as recovered
if an SNR ≥ 7 dip is found within half a duration of the injected T0. The
200 stars span T = 7.4–13.0 and include 9 long-limited and 5 fast-variable
stars.

**Recovery vs expected SNR**, full-window stars (3,764 trials;
`plots/phase2/sensitivity.png`):

| expected SNR | 5–6 | 6–7 | 7–8 | 8–9 | 9–10 | 10–12 | 12–15 | 15–20 | 20–30 | ≥ 30 |
|---|---|---|---|---|---|---|---|---|---|---|
| recovered | 5 % | 19 % | 44 % | 62 % | 77 % | 80 % | 92 % | 97 % | 99 % | 99 % |

- 50 % recovery is at expected SNR ≈ 8.
- Above SNR 15, recovery is 97–100 % at every duration from 1 to 24 h.
- Long transits (16–24 h) need about 1.5× more SNR than short ones to reach
  the same recovery. The likely reasons are that a 3 d window partly
  absorbs a day-long dip in the first pass, and that day-long box noise is
  estimated from few independent boxes.
- Against *white-noise* SNR the curve is shifted right, mostly for long
  transits. At white-noise SNR 15–30, 2–8 h transits are recovered 92–99 %
  of the time, but 24 h transits only 80–87 %. That is red noise on day
  timescales, which the empirical SNR already accounts for.
- By star type at SNR > 15: full window 98.5 %, long-limited 79 %,
  fast variable 82 %.

**Depth needed for ≥ 90 % recovery** (full-window stars; bins are
lower–upper depth; `s0048_sensitivity_by_tmag.csv`):

| Tmag | 2 h transit | 8 h transit | 24 h transit |
|---|---|---|---|
| 10–11 | 1.8–3.2 ppt (95 %) | 1.0–1.8 ppt (89 %) | 1.8–3.2 ppt (100 %) |
| 11–12 | 3.2–5.6 ppt (95 %) | 1.8–3.2 ppt (90 %) | 3.2–5.6 ppt (93 %) |
| 12–13 | 5.6–10 ppt (94 %) | 3.2–5.6 ppt (94 %) | 3.2–5.6 ppt (98 %) |

The T < 10 row is too sparse (11 stars) for a clean threshold. In planet
terms:

- Jupiter- and Saturn-size single transits (≈ 0.6–1 % around a Sun-like
  dwarf; the 5.6–10 ppt bin) are recovered ≥ 94 % of the time down to
  T = 13 at every duration from 2 to 24 h.
- 3–6 ppt transits are recovered ≥ 93 % down to T ≈ 12, and at T 12–13
  only if they last ≥ 8 h.
- Neptune-size transits (≈ 0.1 %) are only partly reachable, and only for
  T ≲ 11 (1–1.8 ppt: 71–89 % at 8–24 h).

### Limitations and next steps (Phase 3)

- **Vetting is the bottleneck.** Needed steps:
  - a shape test that separates a transit from a ramp or a linear trend
    (a box vs. linear-trend Δχ²);
  - centroid and background checks;
  - common-mode tests by camera/CCD rather than sky distance;
  - checks against other sectors, since a single transit should not
    repeat.
- **Eclipsing binaries.** The ≤ 2 R_J cut relies on TIC radii. Grazing or
  diluted binaries still pass.
- **The variability rule is analytic** (a sinusoid through a running
  location filter) and its tolerance is a judgment call. The flagged-star
  injection sample is small (9 + 5 stars).
- **One sector at 10-min cadence.** 1 h boxes contain only 6 cadences.
  The TOI-recall data windows are approximate.

## Phase 3: vetting the Sector 48 candidates

Reproduce with `scripts/phase3_vet.py lc | pixels | fpp | report`; the checks
live in `tesshunt/vetting.py` and have unit tests in
`tests/test_vetting_checks.py`. Every number below comes from
`results/phase3/summary.json`.

### The checks

| # | check | a candidate fails if |
|---|---|---|
| 1 | **Shape** | a limb-darkened `batman` transit model does not beat every non-transit alternative (quadratic ramp, step, and a sudden drop with exponential recovery, i.e. "flare-decay") by ΔBIC ≥ 10. A box is also fitted and its ΔBIC reported. |
| 2 | **Duration** | the fitted T14 is shorter than half the duration of a central transit at P = 20 d (from TIC R★ and M★), unless the fit is grazing (b ≥ 0.9 or b > 1 − Rp/R★). It also fails if even a central transit would need P > 100 yr. |
| 3 | **Edge** | ingress or egress lies within 3 h of a data gap (> 1 h), unless the dip is resolved on both sides. "Resolved" means ≥ 80 % coverage in a baseline window before and after, and in-transit, with pre/post baselines agreeing within 3σ or 25 % of the depth. This can be satisfied in either the PDCSAP light curve or the TESScut pixel light curve. |
| 4a | **Centroid** | the difference-image centroid (TESScut 15×15, out-of-transit minus in-transit) is > 1 px from the target *and* the offset is ≥ 3σ |
| 4b | **Neighbours** | a TIC star within 2 px that is bright enough to cause the dip shows a fractional dip in its pixel that is > 3σ and > 1.5× stronger than the target's |
| 4c | **In pixels** | in a simple TESScut aperture light curve (local, sigma-clipped baselines, no PDC), the dip's SNR is < 0.3× its light-curve SNR, or its depth is outside 0.4–2.5× the light-curve depth |
| 5 | **Asteroids** | SkyBoT (observer code C57 = TESS) finds a known object with V ≤ 19 passing within 60″ during the dip ± 1 h |
| 6 | **Catalogues** | the star is in the TESS EB catalogue (MAST HLSP ∪ Villanova, identical 4,578-star S1–26 lists), unless it hosts a CP/KP planet; or a TOI on the star is TFOPWG FP/FA. A TOI or CTOI match marks the candidate as *known*, not failed. |
| 7 | **FPP** | TRICERATOPS gives FPP ≥ 0.5 or NFPP ≥ 0.1. Inputs: Gaia DR3 background population; period range from "no second transit in the sector" up to the minimum circular period from the duration. |

### What each check removed

Candidates (2,627; tier A first as requested) through the checks in order:

| step | removed | remaining | tier A left | tier B left |
|---|---|---|---|---|
| Phase 2 candidates | – | 2,627 | 991 | 1,636 |
| 1 Shape | 1,704 | 923 | 451 | 472 |
| 2 Duration | 170 | 753 | 369 | 384 |
| 3 Edge | 250 | 503 | 307 | 196 |
| 4 Pixels | 301 | 202 | 83 | 119 |
| 5 Asteroids | 1 | 201 | 82 | 119 |
| 6 Catalogues | 13 (12 EBs, 1 TOI FP) | 188 | 79 | 109 |
| 7 FPP | 134 | **54** | 43 | 11 |

**54 survivors: 28 new and 26 already known as TOIs or CTOIs.** There were
no errors in the final run; transient MAST/SkyBoT failures were retried.
`plots/phase3/funnel.png` shows the same funnel.

Each check applied on its own, to the set it ran on:

| check | ran on | failed |
|---|---|---|
| shape | 2,627 | 1,704 |
| duration | 2,627 | 488 |
| edge | 2,627 | 1,241 |
| centroid (4a) | 767 passing 1–2 | 60 |
| neighbours (4b) | 767 | 76 |
| in pixels (4c) | 767 | 523 |
| asteroids | 767 | 47 (e.g. (146) Lucina, V = 12.4, 9″ from TIC 139415069) |
| catalogues | 767 | 17 |
| FPP | 285 passing 1–6 | 200 |

The shape test and the in-pixels check do most of the work. Most shape
failures are best fit by the sudden-drop/recovery or step models. Of the
dips that reach the pixel stage, 68 % are *not* in the raw pixels at the
depth PDCSAP says. Those are systematics or light-curve processing, not
astrophysics. Asteroids removed only one candidate in sequence because most of the 47
asteroid-affected dips had already failed the edge or pixel checks.

### Validation: does the vetting reject real planets?

Every check was run on TOI-2180 b plus the 27 other TOIs with exactly one
predicted transit in S48 that Phase 2 recovered, regardless of earlier
results. Two of these 28 are TFOPWG false positives (TOI-5621.01 FP,
TOI-1609.01 FA), so the planet set is **26 TOIs (8 CP, 16 PC, 2 APC)**:

| check | planets passing |
|---|---|
| 1 Shape | 21 / 26 |
| 2 Duration | 24 / 26 |
| 3 Edge | 26 / 26 |
| 4 Pixels (4a / 4b / 4c) | 25 / 26 (26 / 26 / 25) |
| 5 Asteroids | 26 / 26 |
| 6 Catalogues | 24 / 26 |
| 7 FPP | 23 / 26 |
| **all seven** | **18 / 26** |

**All 8 confirmed planets (CP) pass every check**: TOI-2180 b,
TOI-1339.02 (HD 191939 c), TOI-1136.01, TOI-1751.01, TOI-1742.01,
TOI-1824.01, TOI-1859.01 and TOI-3837.01. The eight losses, all PC/APC:

- **Shape, SNR ≤ 9.7 (5):** TOI-5635.01, 5710.01, 5718.01, 4109.01 and
  7357.01. At these SNRs the flare-decay model fits about as well as a
  transit. This is the main cost of the vetting: **below SNR ~10, about
  half of real single transits fail the shape test.**
- **Duration (2):** TOI-2271.01 (APC) and TOI-5200.01 have fitted T14 of
  2.4 h and 2.6 h, shorter than half the P = 20 d central duration. TOI-2271.01
  is also in the EB catalogue and has FPP 0.99. TOI-1642.01 (APC) is likewise
  EB-listed with FPP 0.996. The two APCs may well be binaries.
- **Pixels (1):** TOI-5710.01, SNR 7.8, has a pixel SNR of only 1.95 (it
  also fails shape).
- The two known false positives are both rejected: TOI-5621.01 by the
  catalogue check, and TOI-1609.01 by five checks.

### Rule changes made because of the validation set

The spec was followed except where the validation planets showed a check
was rejecting real planets without discriminating. Each change is
recorded in code comments and tests:

- **Box vs. transit is reported but does not reject.** Requiring "a box
  is not strongly preferred" (ΔBIC > −6) rejected 2 of 26 validation
  planets (TOI-1339.02 = HD 191939 c, and TOI-3832.01). It removed only
  35 candidates, a similar rate. At 10-min cadence a giant planet's ingress
  spans 1–2 cadences, so real transits can look box-like. The stricter
  rule is available as `shape_passes(..., dbic_box_min=-6)`.
- **Transit fit seeded from the box fit.** Without it, the fit fell into
  local minima for dips longer than, or offset from, the detection box.
  TOI-1339.02 went from ΔBIC(box) −86 to −51.
- **Edge check accepts TESScut data.** The raw FFI pixels keep cadences that
  PDCSAP blanks at orbit starts. TOI-2180 b's transit is resolved on both
  sides in TESScut but not in PDCSAP. The TESScut route counts only if the
  dip is also confirmed in those pixels (4c). Without that condition,
  orbit-start ramps were rescued.
- **Pixel dip measured against local, sigma-clipped baselines, with a
  depth-consistency condition.** A single global baseline let a
  scattered-light spike fake a pixel "dip" (TIC 17238617). Validation
  planets have pixel/light-curve depth ratios of 0.65–1.74; the
  orbit-start artefacts had 0.13–0.26 or 3.3–6.
- **EB catalogue matches do not reject CP/KP planet hosts.** The EB
  catalogue lists HD 191939.
- **Long-duration limit added** (P > 100 yr for a central transit). One
  shortlisted dip needed a 158-yr orbit.

### Shortlist

The 28 new candidates, ranked by FPP then SNR. Each has a vetting sheet
with a light-curve zoom, every model fit, the TESScut aperture light
curve, the out-of-transit and difference images, and a pass/fail table.
R_p = √depth · R★(TIC). P_min is the minimum circular period consistent
with the duration.

| # | TIC | tier | T | R★ | T0 (BTJD) | depth (ppm) | T14 (h) | b | SNR | R_p (R_J) | P_min (d) | FPP | review notes | sheet |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 239198203 | A | 12.4 | 0.74 | 2628.90 | 45868 | 4.1 | 0.66 | 75.8 | 1.54 | 12 | 0.000 | deep (>2 %): check for stellar companion | [sheet](plots/phase3/sheets/01_tic239198203.png) |
| 2 | 376848623 | A | 9.6 | 3.57 | 2633.20 | 1002 | 42.8 | 1.24 | 14.7 | 1.10 | 195 | 0.000 | very long dip: compare with stellar variability; grazing/V-shaped: EB-like geometry | [sheet](plots/phase3/sheets/02_tic376848623.png) |
| 3 | 147971475 | A | 12.4 | 2.44 | 2634.00 | 2602 | 56.0 | 0.47 | 18.2 | 1.21 | 1700 | 0.000 | very long dip: compare with stellar variability; centroid offset 0.76 px at 22 sigma | [sheet](plots/phase3/sheets/03_tic147971475.png) |
| 4 | 252795236 | A | 12.1 | 0.81 | 2614.60 | 7938 | 13.9 | 0.01 | 15.0 | 0.70 | 675 | 0.000 | difference image too faint for a centroid | [sheet](plots/phase3/sheets/04_tic252795236.png) |
| 5 | 142905733 | B | 12.4 | 0.61 | 2627.00 | 2416 | 34.4 | 0.00 | 18.0 | 0.29 | 15852 | 0.000 | very long dip: compare with stellar variability | [sheet](plots/phase3/sheets/05_tic142905733.png) |
| 6 | 459793183 | A | 10.1 | 1.62 | 2619.70 | 3932 | 8.6 | 0.82 | 26.9 | 0.99 | 23 | 0.000 |  | [sheet](plots/phase3/sheets/06_tic459793183.png) |
| 7 | 298092113 | A | 11.9 | 0.88 | 2614.50 | 4432 | 19.7 | 0.00 | 13.3 | 0.57 | 1461 | 0.000 |  | [sheet](plots/phase3/sheets/07_tic298092113.png) |
| 8 | 284993804 | B | 12.6 | 1.16 | 2614.40 | 6897 | 25.7 | 0.76 | 11.4 | 0.94 | 1354 | 0.000 | difference image too faint for a centroid | [sheet](plots/phase3/sheets/08_tic284993804.png) |
| 9 | 16222047 | A | 11.3 | 1.78 | 2613.10 | 1816 | 11.2 | 0.00 | 12.4 | 0.74 | 39 | 0.000 |  | [sheet](plots/phase3/sheets/09_tic16222047.png) |
| 10 | 335661035 | B | 12.9 | 0.69 | 2627.10 | 9408 | 60.4 | 1.26 | 34.2 | 0.65 | 28204 | 0.000 | very long dip: compare with stellar variability; grazing/V-shaped: EB-like geometry; centroid offset 0.68 px at 21 sigma | [sheet](plots/phase3/sheets/10_tic335661035.png) |
| 11 | 88704042 | A | 7.1 | 1.05 | 2628.80 | 212 | 59.7 | 0.99 | 13.5 | 0.15 | 28869 | 0.000 | very long dip: compare with stellar variability; grazing/V-shaped: EB-like geometry; centroid offset 0.78 px at 51 sigma | [sheet](plots/phase3/sheets/11_tic88704042.png) |
| 12 | 155873261 | B | 10.2 | 2.40 | 2625.60 | 5229 | 9.9 | 0.88 | 11.1 | 1.69 | 9 | 0.004 |  | [sheet](plots/phase3/sheets/12_tic155873261.png) |
| 13 | 165685135 | A | 11.7 | 1.20 | 2618.70 | 2737 | 11.8 | 0.80 | 16.5 | 0.61 | 147 | 0.012 |  | [sheet](plots/phase3/sheets/13_tic165685135.png) |
| 14 | 188506574 | B | 11.5 | 1.87 | 2614.40 | 1860 | 23.8 | 0.01 | 12.4 | 0.78 | 367 | 0.013 |  | [sheet](plots/phase3/sheets/14_tic188506574.png) |
| 15 | 27068699 | A | 12.9 | 0.67 | 2618.70 | 19892 | 5.6 | 0.60 | 20.7 | 0.92 | 45 | 0.014 |  | [sheet](plots/phase3/sheets/15_tic27068699.png) |
| 16 | 156372726 | A | 12.4 | 0.84 | 2620.30 | 2691 | 46.2 | 0.90 | 13.7 | 0.43 | 18468 | 0.018 | very long dip: compare with stellar variability; difference image too faint for a centroid | [sheet](plots/phase3/sheets/16_tic156372726.png) |
| 17 | 154565237 | A | 10.5 | 0.95 | 2613.70 | 1214 | 4.4 | 0.00 | 8.4 | 0.32 | 17 | 0.026 | neighbour inside the target pixel | [sheet](plots/phase3/sheets/17_tic154565237.png) |
| 18 | 237109179 | A | 10.5 | 1.20 | 2611.80 | 1374 | 5.6 | 0.33 | 8.7 | 0.43 | 15 | 0.094 |  | [sheet](plots/phase3/sheets/18_tic237109179.png) |
| 19 | 165410329 | A | 12.8 | 0.64 | 2619.00 | 5708 | 4.9 | 1.27 | 8.1 | 0.47 | 19 | 0.109 | grazing/V-shaped: EB-like geometry | [sheet](plots/phase3/sheets/19_tic165410329.png) |
| 20 | 157264264 | A | 10.2 | 1.99 | 2620.30 | 3547 | 5.3 | 0.00 | 28.8 | 1.16 | 3 | 0.116 |  | [sheet](plots/phase3/sheets/20_tic157264264.png) |
| 21 | 159540437 | B | 11.6 | 1.85 | 2611.50 | 4804 | 11.6 | 0.92 | 8.4 | 1.25 | 45 | 0.118 |  | [sheet](plots/phase3/sheets/21_tic159540437.png) |
| 22 | 85911158 | A | 12.8 | 1.11 | 2628.20 | 4916 | 31.0 | 1.26 | 9.7 | 0.76 | 1409 | 0.137 | very long dip: compare with stellar variability; grazing/V-shaped: EB-like geometry; centroid offset 0.70 px at 16 sigma | [sheet](plots/phase3/sheets/22_tic85911158.png) |
| 23 | 159159589 | A | 11.5 | 0.95 | 2619.70 | 14185 | 4.2 | 0.85 | 41.0 | 1.10 | 9 | 0.140 |  | [sheet](plots/phase3/sheets/23_tic159159589.png) |
| 24 | 233575173 | A | 12.1 | 2.12 | 2631.80 | 2972 | 29.0 | 0.93 | 14.8 | 1.12 | 457 | 0.143 |  | [sheet](plots/phase3/sheets/24_tic233575173.png) |
| 25 | 55779787 | B | 12.5 | 0.96 | 2611.00 | 4455 | 10.3 | 1.31 | 7.9 | 0.63 | 80 | 0.194 | grazing/V-shaped: EB-like geometry | [sheet](plots/phase3/sheets/25_tic55779787.png) |
| 26 | 29235065 | B | 12.2 | 1.70 | 2610.40 | 5021 | 18.7 | 1.27 | 8.1 | 1.17 | 114 | 0.227 | grazing/V-shaped: EB-like geometry | [sheet](plots/phase3/sheets/26_tic29235065.png) |
| 27 | 95747180 | A | 12.6 | 1.22 | 2621.10 | 8537 | 4.6 | 0.75 | 13.1 | 1.10 | 8 | 0.278 |  | [sheet](plots/phase3/sheets/27_tic95747180.png) |
| 28 | 198206622 | A | 12.1 | 1.74 | 2626.10 | 8836 | 7.5 | 0.72 | 19.1 | 1.59 | 10 | 0.362 |  | [sheet](plots/phase3/sheets/28_tic198206622.png) |

**How to read the list.**

- **Look first at the 14 candidates with no review notes:** #6, 7, 9,
  12–15, 18, 20, 21, 23, 24, 27 and 28. After them come #1 (deep), #4 and #8
  (difference image too faint for a centroid). #6 (TIC 459793183: 3.9 ppt,
  8.6 h, ~1 R_J on a T = 10.1 star) and #1 (TIC 239198203: 4.6 %, 4.1 h,
  ~1.5 R_J) are the cleanest examples.
- **Long, grazing dips with sub-pixel centroid offsets** are more likely
  stellar variability, systematics or blends, despite FPP ≈ 0. These carry
  the notes "very long dip" (T14 > 30 h) and "centroid offset … at N sigma".
- **Depths of 2–5 %** give R_p ≈ 1.1–1.7 R_J and could equally be low-mass
  stars. Radial velocities would decide.
- **FPP values are indicative only.** TRICERATOPS is designed for
  phase-folded periodic transits. For a single transit the period is only
  bounded, and FPP ≈ 0 mostly reflects the lack of a better-fitting
  eclipsing-binary scenario over that range. No candidate here is
  "validated".

**Known candidates that survived all checks** (a positive control for
the checks):

| # | TIC | known as | depth (ppm) | T14 (h) | SNR | FPP |
|---|---|---|---|---|---|---|
| 29 | 298663873 | TOI 2180.01;CTOI 298663873.01 | 2559 | 23.3 | 27.0 | 0.000 |
| 30 | 289718679 | CTOI 289718679.01 | 718 | 8.5 | 14.9 | 0.000 |
| 31 | 258799245 | CTOI 258799245.01 | 2826 | 9.1 | 24.5 | 0.001 |
| 32 | 287080092 | TOI 1751.01 | 585 | 8.7 | 10.5 | 0.001 |
| 33 | 229742722 | TOI 1859.01;CTOI 229742722.01 | 5312 | 4.8 | 37.5 | 0.003 |
| 34 | 142387023 | TOI 1824.01;CTOI 142387023.01 | 1135 | 4.2 | 7.7 | 0.006 |
| 35 | 219857012 | TOI 1742.01 | 486 | 7.5 | 8.4 | 0.012 |
| 36 | 288636342 | TOI 1692.01;1692.02;CTOI 288636342.02 | 2874 | 7.7 | 16.8 | 0.017 |
| 37 | 420112217 | TOI 4163.01 | 4958 | 6.3 | 9.5 | 0.019 |
| 38 | 459969957 | TOI 1274.01 | 12857 | 4.3 | 43.1 | 0.024 |
| 39 | 53481079 | TOI 5717.01 | 1432 | 4.1 | 8.2 | 0.024 |
| 40 | 85293053 | TOI 1772.01;1772.02;CTOI 85293053.02 | 3904 | 5.5 | 20.0 | 0.027 |
| 41 | 224299081 | TOI 5619.01;CTOI 224299081.01;224299081.02 | 2833 | 6.0 | 14.3 | 0.031 |
| 42 | 99919551 | CTOI 99919551.01 | 7432 | 4.1 | 11.9 | 0.031 |
| 43 | 302728777 | CTOI 302728777.01 | 3677 | 5.9 | 11.3 | 0.033 |
| 44 | 115861501 | CTOI 115861501.01 | 2870 | 6.3 | 10.8 | 0.048 |
| 45 | 353807936 | TOI 5650.01 | 4258 | 7.6 | 11.6 | 0.052 |
| 46 | 356736743 | TOI 3832.01 | 7847 | 4.3 | 11.5 | 0.054 |
| 47 | 198189972 | TOI 4117.01 | 6585 | 5.8 | 14.3 | 0.075 |
| 48 | 458481611 | TOI 5722.01;CTOI 458481611.01 | 7669 | 6.0 | 14.8 | 0.087 |
| 49 | 288132261 | TOI 1258.01 | 844 | 4.6 | 8.1 | 0.100 |
| 50 | 147898114 | CTOI 147898114.01 | 9260 | 4.0 | 15.5 | 0.253 |
| 51 | 229539379 | CTOI 229539379.01 | 6586 | 8.5 | 10.9 | 0.332 |
| 52 | 233496435 | CTOI 233496435.01 | 7131 | 4.1 | 10.9 | 0.376 |
| 53 | 156512187 | TOI 7018.01 | 36202 | 5.4 | 48.3 | 0.418 |
| 54 | 308073360 | TOI 3837.01;CTOI 308073360.01 | 7280 | 5.4 | 20.0 | 0.466 |

### Limitations

- **One sector, 10-min FFIs.** Every candidate needs a second transit or
  RV follow-up. Checking other TESS sectors for these stars (no repeat
  allowed for P below their span) is the natural next vetting step.
- **Shape tests lose real planets below SNR ~10** (5 of 11 validation
  planets with SNR < 10 failed). A candidate list that needs completeness
  at low SNR should rank by ΔBIC rather than cut on it.
- **The pixel analysis is simple:** a 15×15 TESScut cutout, a
  flux-threshold aperture, a per-frame median background, and a
  flux-weighted difference-image centroid, not PRF fitting. Neighbours
  within the same pixel are not tested (they are listed and left to
  TRICERATOPS).
- **The EB catalogue covers only S1–26 stars.** Stars first observed later
  get no EB check.
- **TRICERATOPS uses the TESScut aperture defined here, not the SPOC
  aperture.** The MAST name resolver is unreachable from this environment,
  so TRICERATOPS target lookup and TESScut searches were patched to use
  coordinates. Its Gaia DR3 query runs over the archive's HTTPS TAP
  endpoint.
- The validation set is small (26 planets). The per-check pass rates carry
  ±10–20 % binomial uncertainty.
