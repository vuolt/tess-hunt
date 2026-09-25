# Phase 1 results: single-transit detection in TESS

Everything here can be reproduced with `scripts/run_known_target.py` and
`scripts/injection_recovery.py` (see README). Raw FITS files go to `data/`,
which git ignores.

## Pipeline

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

## 1. Known planet: TOI-2180 b

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

## 2. Injection-recovery

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

## Limitations and next steps

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
