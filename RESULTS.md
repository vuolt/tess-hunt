# Results

- [Phase 1](#phase-1-single-transit-detection-on-one-star-2-min-data): the pipeline
  on one star with a known long-period planet, plus injection-recovery.
- [Phase 2](#phase-2-full-sector-search-tess-spoc-ffi-sector-48): a full-sector
  search of 128,455 dwarfs in TESS-SPOC FFI light curves.

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
