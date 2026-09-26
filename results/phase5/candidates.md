# Phase 5: expert checks, candidate by candidate

Produced by `python scripts/phase5.py report`. Numbers: `phase5_checks.csv`; calibration on known planets: `calibration.csv`; plots: `plots/phase5/`.

**Calibration.** The same checks on 26 known planets from the validation sets: 7 strong, 15 plausible, 4 doubtful. Serious flags per check: aperture 0, gp 1, evolved 0, binarity 3, eclipsing 0, physical 0, too_large 0. On the 2 known false positives: plausible 1, doubtful 1.

## TIC 165685135 (Sector 48) — **plausible**

Phase 4: **submit** → Phase 5: **submit**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (2162, 2356, 2479, 2543 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 2719 ppm (0.99× Phase 3), significance 8.1σ (Phase 3: 16.5).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 364; radius 1.20 R☉ from Gaia FLAME). The companion is then about 0.61 R_J (TIC radius gave 0.61).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 2.0 km/s over 10 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the shape implies a circular-orbit period of ~217 d (105–776).

**Verdict: plausible** — no check failed, but some raise minor concerns: evolved.
Plot: `plots/phase5/s0048_165685135.png`

## TIC 237109179 (Sector 48) — **plausible**

Phase 4: **submit** → Phase 5: **submit**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (1165, 1454, 2763, 5170 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1355 ppm (0.99× Phase 3), significance 6.2σ (Phase 3: 8.7). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 399; radius 1.23 R☉ from Gaia FLAME). The companion is then about 0.44 R_J (TIC radius gave 0.43).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.3 km/s over 16 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the shape implies a circular-orbit period of ~22 d (14–84).

**Verdict: plausible** — no check failed, but some raise minor concerns: aperture_marginal, gp_weak, evolved.
Plot: `plots/phase5/s0048_237109179.png`

## TIC 239198203 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (41560, 42192, 41469, 40942 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 45828 ppm (1.00× Phase 3), significance 22.9σ (Phase 3: 75.9).
- **Stellar check.** Gaia DR3 has no evolutionary parameters for this star; its TIC radius (0.74 R☉) is consistent with a dwarf, giving a companion of about 1.54 R_J (17 R⊕). Not independently confirmed.
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.9 km/s over 29 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)): of the 40 periods allowed by the second dip, 27 fit a near-circular orbit; the best is P = 22.7 d (eccentricity ≥ ~0.05).

**Verdict: plausible** — no check failed, but some raise minor concerns: star_unverified.
Plot: `plots/phase5/s0048_239198203.png`

## TIC 142905733 (Sector 48) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: gp)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (987, 1115, 1766, 2011 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1791 ppm (0.74× Phase 3), significance 3.6σ (Phase 3: 18.0). **The dip is not robust to the choice of noise model.**
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 334); radius 0.68 R☉ (Gaia FLAME), so the companion is about 0.32 R_J (4 R⊕).
- **Gaia binarity.** Gaia's astrometry is noisier than for a single star (RV error 1.2 km/s over 28 transits: no significant scatter; RUWE 1.75), but the image shape and velocities show no companion.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the shape implies a circular-orbit period of ~13540 d (10702–17038).

**Verdict: doubtful** — failed: gp.
Plot: `plots/phase5/s0048_142905733.png`

## TIC 459793183 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (3414, 3548, 3669, 3818 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 3900 ppm (0.99× Phase 3), significance 13.5σ (Phase 3: 26.9).
- **Stellar check.** The star is **subgiant** (TIC radius 1.62 R☉ only; Gaia DR3 has no evolutionary parameters for this star; radius 1.62 R☉ from TIC). The companion is then about 1.40 R_J (TIC radius gave 1.40).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (no Gaia DR3 source).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the shape implies a circular-orbit period of ~65 d (35–99).

**Verdict: plausible** — no check failed, but some raise minor concerns: evolved.
Plot: `plots/phase5/s0048_459793183.png`

## TIC 16222047 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (1260, 1260, 1656, 1605 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1892 ppm (1.04× Phase 3), significance 6.1σ (Phase 3: 12.4). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 448; radius 1.75 R☉ from Gaia FLAME). The companion is then about 0.73 R_J (TIC radius gave 0.74).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.5 km/s over 34 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the shape implies a circular-orbit period of ~59 d (43–134).

**Verdict: plausible** — no check failed, but some raise minor concerns: gp_weak, evolved.
Plot: `plots/phase5/s0048_16222047.png`

## TIC 155873261 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (4448, 4071, 4005, 3992 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 7563 ppm (1.45× Phase 3), significance 10.2σ (Phase 3: 11.1).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 458; radius 2.65 R☉ from Gaia FLAME). The companion is then about 1.87 R_J (TIC radius gave 1.69).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (no Gaia RV).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the shape implies a circular-orbit period of ~1 d (1–3).

**Verdict: plausible** — no check failed, but some raise minor concerns: evolved.
Plot: `plots/phase5/s0048_155873261.png`

## TIC 154565237 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (869, 1063, 1092, 1180 ppm), as expected if the dip comes from this star (TIC 950815449 lies within 1 px and cannot be separated this way).
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1231 ppm (1.01× Phase 3), significance 8.3σ (Phase 3: 8.4).
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 113); radius 0.94 R☉ (Gaia FLAME), so the companion is about 0.32 R_J (4 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.3 km/s over 14 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME): of the 6 periods allowed by the second dip, 0 fit a near-circular orbit; the best is P = 64.3 d (eccentricity ≥ ~0.34). Every allowed period needs a fairly eccentric orbit.

**Verdict: plausible** — no check failed, but some raise minor concerns: eccentric.
Plot: `plots/phase5/s0048_154565237.png`

## TIC 157264264 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (3655, 3459, 3252, 3385 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 3566 ppm (1.01× Phase 3), significance 13.7σ (Phase 3: 28.8).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 457; radius 1.88 R☉ from Gaia FLAME). The companion is then about 1.09 R_J (TIC radius gave 1.16).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.2 km/s over 18 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the shape implies a circular-orbit period of ~4 d (3–6).

**Verdict: plausible** — no check failed, but some raise minor concerns: evolved.
Plot: `plots/phase5/s0048_157264264.png`

## TIC 159540437 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (2237, 2331, 2952, 3776 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 4293 ppm (0.89× Phase 3), significance 5.0σ (Phase 3: 8.4). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (TIC radius 1.85 R☉ only; Gaia DR3 has no evolutionary parameters for this star; radius 1.85 R☉ from TIC). The companion is then about 1.25 R_J (TIC radius gave 1.25).
- **Gaia binarity.** Gaia's astrometry is noisier than for a single star (RV error 1.5 km/s over 16 transits: no significant scatter; RUWE 2.84), but the image shape and velocities show no companion.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the shape implies a circular-orbit period of ~235 d (92–393). The fit is probably grazing (p = 0.58), so the size is uncertain.

**Verdict: plausible** — no check failed, but some raise minor concerns: gp_weak, evolved, ruwe.
Plot: `plots/phase5/s0048_159540437.png`

## TIC 159159589 (Sector 48) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: binarity)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (12362, 11898, 11897, 11737 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 14184 ppm (1.00× Phase 3), significance 16.8σ (Phase 3: 41.0).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 363; radius 1.11 R☉ from Gaia FLAME). The companion is then about 1.29 R_J (TIC radius gave 1.10).
- **Gaia binarity.** **Gaia's radial velocities swing by tens of km/s: a companion star on a short orbit, which could itself cause the dip.** radial velocity varies (p = 9.0e-06, renormalised GOF 4.0, amplitude 16.7 km/s); RUWE 3.88.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the shape implies a circular-orbit period of ~18 d (12–25).

**Verdict: doubtful** — failed: binarity.
Plot: `plots/phase5/s0048_159159589.png`

## TIC 29235065 (Sector 48) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: aperture)

- **Aperture test.** **The depth changes with aperture size** (2418, 5044, 8896, 17847 ppm in 1, 1.5, 2 and 3 px apertures; worst 9.6σ from what a signal on this star would give), and fits neighbour TIC 29235068 3.3 px away better. The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 5801 ppm (1.16× Phase 3), significance 5.9σ (Phase 3: 8.1). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (TIC radius 1.70 R☉ only; Gaia DR3 has no evolutionary parameters for this star; radius 1.70 R☉ from TIC). The companion is then about 1.17 R_J (TIC radius gave 1.17).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (no Gaia RV).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the shape implies a circular-orbit period of ~30 d (16–191).

**Verdict: doubtful** — failed: aperture.
Plot: `plots/phase5/s0048_29235065.png`

## TIC 95747180 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (9265, 8276, 7772, 8387 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 8508 ppm (1.00× Phase 3), significance 9.1σ (Phase 3: 13.1).
- **Stellar check.** Gaia DR3 has no evolutionary parameters for this star; its TIC radius (1.22 R☉) is consistent with a dwarf, giving a companion of about 1.10 R_J (12 R⊕). Not independently confirmed.
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.6 km/s over 28 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)): for P = 80.47 d the shape needs 9.1× the star's density on a circular orbit, i.e. an eccentricity of at least ~0.63 (≥ 0.40 at 1σ); for P = 40.24 d the shape needs 4.5× the star's density on a circular orbit, i.e. an eccentricity of at least ~0.47 (≥ 0.19 at 1σ). This is possible but most likely needs a fairly eccentric orbit.

**Verdict: plausible** — no check failed, but some raise minor concerns: star_unverified, eccentric.
Plot: `plots/phase5/s0048_95747180.png`

## TIC 298663873 (Sector 48, control) — **plausible**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (3977, 4638, 4704, 4788 ppm), as expected if the dip comes from this star (TIC 298663875 lies within 1 px and cannot be separated this way).
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 5136 ppm (2.01× Phase 3), significance 15.9σ (Phase 3: 27.0).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 447; radius 1.63 R☉ from Gaia FLAME). The companion is then about 0.80 R_J (TIC radius gave 0.80).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.1 km/s over 17 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME): for P = 260.17 d the shape needs 0.6× the star's density on a circular orbit, i.e. an eccentricity of at least ~0.18 (≥ 0.15 at 1σ).

**Verdict: plausible** — no check failed, but some raise minor concerns: aperture_marginal, evolved.
Plot: `plots/phase5/s0048_298663873.png`
