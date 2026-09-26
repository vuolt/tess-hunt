# Phase 5: expert checks, candidate by candidate

Produced by `python scripts/phase5.py report`. Numbers: `phase5_checks.csv`; calibration on known planets: `calibration.csv`; plots: `plots/phase5/`.

**Calibration.** The same checks on 52 known planets from the validation sets: 17 strong, 27 plausible, 8 doubtful. Serious flags per check: aperture 0, gp 0, evolved 0, binarity 8, eclipsing 0, physical 0, too_large 0, deep_eclipse 0. On the 6 known false positives: plausible 3, doubtful 3.

## TIC 165685135 (Sector 48) — **plausible**

Phase 4: **submit** → Phase 5: **submit**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (2162, 2356, 2479, 2543 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 2719 ppm (0.99× Phase 3), significance 8.1σ (Phase 3: 16.5).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 364; radius 1.20 R☉ from Gaia FLAME). The companion is then about 0.61 R_J (TIC radius gave 0.61).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 2.0 km/s over 10 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~170 d (104–484); 87% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: evolved.
Plot: `plots/phase5/s0048_165685135.png`

## TIC 237109179 (Sector 48) — **plausible**

Phase 4: **submit** → Phase 5: **submit**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (1165, 1454, 2763, 5170 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1355 ppm (0.99× Phase 3), significance 6.2σ (Phase 3: 8.7). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 399; radius 1.23 R☉ from Gaia FLAME). The companion is then about 0.44 R_J (TIC radius gave 0.43).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.3 km/s over 16 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~20 d (13–56); 1% of that range is still allowed by TESS's other observations. An allowed period is reachable with a moderately eccentric orbit (e ≤ 0.5).
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: aperture_marginal, gp_weak, evolved, eccentric.
Plot: `plots/phase5/s0048_237109179.png`

## TIC 239198203 (Sector 48) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: deep_eclipse)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (41560, 42192, 41469, 40942 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 45828 ppm (1.00× Phase 3), significance 22.9σ (Phase 3: 75.9).
- **Stellar check.** Gaia DR3 has no evolutionary parameters for this star; its TIC radius (0.74 R☉) is consistent with a dwarf, giving a companion of about 1.54 R_J (17 R⊕). Not independently confirmed.
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.9 km/s over 29 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)): of the 40 periods allowed by the second dip, 27 fit a near-circular orbit; the best is P = 22.7 d (eccentricity ≥ ~0.05).
- **Other TESS sectors.** **Another TESS sector shows a much deeper eclipse on this star**: 21.1 % deep in Sector 21 (BTJD 1891.00, SNR 595), against 4.59 % for this dip. The star is an eclipsing binary; this dip is probably its shallower (secondary) eclipse.

**Verdict: doubtful** — failed: deep_eclipse.
Plot: `plots/phase5/s0048_239198203.png`

## TIC 142905733 (Sector 48) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: gp)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (987, 1115, 1766, 2011 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1791 ppm (0.74× Phase 3), significance 3.6σ (Phase 3: 18.0). **The dip is not robust to the choice of noise model.**
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 334); radius 0.68 R☉ (Gaia FLAME), so the companion is about 0.32 R_J (4 R⊕).
- **Gaia binarity.** Gaia's astrometry is noisier than for a single star (RV error 1.2 km/s over 28 transits: no significant scatter; RUWE 1.75), but the image shape and velocities show no companion.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~13821 d (10929–17433); 100% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: gp.
Plot: `plots/phase5/s0048_142905733.png`

## TIC 459793183 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (3414, 3548, 3669, 3818 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 3900 ppm (0.99× Phase 3), significance 13.5σ (Phase 3: 26.9).
- **Stellar check.** The star is **subgiant** (TIC radius 1.62 R☉ only; Gaia DR3 has no evolutionary parameters for this star; radius 1.62 R☉ from TIC). The companion is then about 1.40 R_J (TIC radius gave 1.40).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (no Gaia DR3 source).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the dip's shape implies a circular-orbit period of ~65 d (37–100); 28% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: evolved.
Plot: `plots/phase5/s0048_459793183.png`

## TIC 16222047 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (1260, 1260, 1656, 1605 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1892 ppm (1.04× Phase 3), significance 6.1σ (Phase 3: 12.4). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 448; radius 1.75 R☉ from Gaia FLAME). The companion is then about 0.73 R_J (TIC radius gave 0.74).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.5 km/s over 34 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~58 d (43–140); 63% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: gp_weak, evolved.
Plot: `plots/phase5/s0048_16222047.png`

## TIC 155873261 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (4448, 4071, 4005, 3992 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 7563 ppm (1.45× Phase 3), significance 10.2σ (Phase 3: 11.1).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 458; radius 2.65 R☉ from Gaia FLAME). The companion is then about 1.87 R_J (TIC radius gave 1.69).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (no Gaia RV).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~1 d (1–856); 22% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: evolved.
Plot: `plots/phase5/s0048_155873261.png`

## TIC 154565237 (Sector 48) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: deep_eclipse)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (869, 1063, 1092, 1180 ppm), as expected if the dip comes from this star (TIC 950815449 lies within 1 px and cannot be separated this way).
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1231 ppm (1.01× Phase 3), significance 8.3σ (Phase 3: 8.4).
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 113); radius 0.94 R☉ (Gaia FLAME), so the companion is about 0.32 R_J (4 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.3 km/s over 14 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME): of the 6 periods allowed by the second dip, 0 fit a near-circular orbit; the best is P = 64.3 d (eccentricity ≥ ~0.34). Every allowed period needs a fairly eccentric orbit.
- **Other TESS sectors.** **Another TESS sector shows a much deeper eclipse on this star**: 5.6 % deep in Sector 47 (BTJD 2594.32, SNR 196), against 0.12 % for this dip. The star is an eclipsing binary; this dip is probably its shallower (secondary) eclipse.

**Verdict: doubtful** — failed: deep_eclipse.
Plot: `plots/phase5/s0048_154565237.png`

## TIC 157264264 (Sector 48) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: physical)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (3655, 3459, 3252, 3385 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 3566 ppm (1.01× Phase 3), significance 13.7σ (Phase 3: 28.8).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 457; radius 1.88 R☉ from Gaia FLAME). The companion is then about 1.09 R_J (TIC radius gave 1.16).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.2 km/s over 18 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~4 d (3–6); 0% of that range is still allowed by TESS's other observations. **No allowed period is reachable even with e ≤ 0.5.**
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: physical.
Plot: `plots/phase5/s0048_157264264.png`

## TIC 159540437 (Sector 48) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: deep_eclipse)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (2237, 2331, 2952, 3776 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 4293 ppm (0.89× Phase 3), significance 5.0σ (Phase 3: 8.4). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (TIC radius 1.85 R☉ only; Gaia DR3 has no evolutionary parameters for this star; radius 1.85 R☉ from TIC). The companion is then about 1.25 R_J (TIC radius gave 1.25).
- **Gaia binarity.** Gaia's astrometry is noisier than for a single star (RV error 1.5 km/s over 16 transits: no significant scatter; RUWE 2.84), but the image shape and velocities show no companion.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the dip's shape implies a circular-orbit period of ~268 d (130–435); 19% of that range is still allowed by TESS's other observations. The fit is probably grazing (p = 0.67), so the size is uncertain.
- **Other TESS sectors.** **Another TESS sector shows a much deeper eclipse on this star**: 4.1 % deep in Sector 60 (BTJD 2939.60, SNR 127), against 0.48 % for this dip. The star is an eclipsing binary; this dip is probably its shallower (secondary) eclipse.

**Verdict: doubtful** — failed: deep_eclipse.
Plot: `plots/phase5/s0048_159540437.png`

## TIC 159159589 (Sector 48) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: binarity)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (12362, 11898, 11897, 11737 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 14184 ppm (1.00× Phase 3), significance 16.8σ (Phase 3: 41.0).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 363; radius 1.11 R☉ from Gaia FLAME). The companion is then about 1.29 R_J (TIC radius gave 1.10).
- **Gaia binarity.** **Gaia's radial velocities swing by tens of km/s: a companion star on a short orbit, which could itself cause the dip.** radial velocity varies (p = 9.0e-06, renormalised GOF 4.0, amplitude 16.7 km/s); RUWE 3.88.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the dip's shape implies a circular-orbit period of ~17 d (11–25); 1% of that range is still allowed by TESS's other observations. An allowed period is reachable with a moderately eccentric orbit (e ≤ 0.5).
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: binarity.
Plot: `plots/phase5/s0048_159159589.png`

## TIC 29235065 (Sector 48) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: aperture)

- **Aperture test.** **The depth changes with aperture size** (2418, 5044, 8896, 17847 ppm in 1, 1.5, 2 and 3 px apertures; worst 9.6σ from what a signal on this star would give), and fits neighbour TIC 29235068 3.3 px away better. The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 5801 ppm (1.16× Phase 3), significance 5.9σ (Phase 3: 8.1). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (TIC radius 1.70 R☉ only; Gaia DR3 has no evolutionary parameters for this star; radius 1.70 R☉ from TIC). The companion is then about 1.17 R_J (TIC radius gave 1.17).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (no Gaia RV).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the dip's shape implies a circular-orbit period of ~29 d (12–107); 30% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture.
Plot: `plots/phase5/s0048_29235065.png`

## TIC 95747180 (Sector 48) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (9265, 8276, 7772, 8387 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 8508 ppm (1.00× Phase 3), significance 9.1σ (Phase 3: 13.1).
- **Stellar check.** Gaia DR3 has no evolutionary parameters for this star; its TIC radius (1.22 R☉) is consistent with a dwarf, giving a companion of about 1.10 R_J (12 R⊕). Not independently confirmed.
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.6 km/s over 28 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)): for P = 80.47 d the shape needs 8.4× the star's density on a circular orbit, i.e. an eccentricity of at least ~0.61 (≥ 0.40 at 1σ); for P = 40.24 d the shape needs 4.2× the star's density on a circular orbit, i.e. an eccentricity of at least ~0.44 (≥ 0.20 at 1σ). This is possible but most likely needs a fairly eccentric orbit.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: star_unverified, eccentric.
Plot: `plots/phase5/s0048_95747180.png`

## TIC 298663873 (Sector 48, control) — **plausible**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (3977, 4638, 4704, 4788 ppm), as expected if the dip comes from this star (TIC 298663875 lies within 1 px and cannot be separated this way).
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 5136 ppm (2.01× Phase 3), significance 15.9σ (Phase 3: 27.0).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 447; radius 1.63 R☉ from Gaia FLAME). The companion is then about 0.80 R_J (TIC radius gave 0.80).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.1 km/s over 17 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME): for P = 260.17 d the shape needs 0.6× the star's density on a circular orbit, i.e. an eccentricity of at least ~0.19 (≥ 0.15 at 1σ).
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: aperture_marginal, evolved.
Plot: `plots/phase5/s0048_298663873.png`

## TIC 361977189 (Sector 21) — **doubtful**

Phase 4: **submit** → Phase 5: **maybe** (Phase 5 doubtful: aperture)

- **Aperture test.** **The depth changes with aperture size** (8775, 2605, 2348, 1641 ppm in 1, 1.5, 2 and 3 px apertures; worst 5.8σ from what a signal on this star would give). The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1432 ppm (0.84× Phase 3), significance 5.3σ (Phase 3: 11.2). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (TIC radius 2.93 R☉ only; Gaia DR3 has no evolutionary parameters for this star; radius 2.93 R☉ from TIC). The companion is then about 1.18 R_J (TIC radius gave 1.18).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.2 km/s over 23 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the dip's shape implies a circular-orbit period of ~94 d (60–200); 42% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture.
Plot: `plots/phase5/s0021_361977189.png`

## TIC 285082902 (Sector 21) — **plausible**

Phase 4: **submit** → Phase 5: **submit**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (3555, 3217, 2722, 4183 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1842 ppm (0.90× Phase 3), significance 5.5σ (Phase 3: 11.0). It survives, but only moderately.
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 318); radius 1.07 R☉ (Gaia FLAME), so the companion is about 0.47 R_J (5 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.5 km/s over 18 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~716 d (462–2316); 97% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: gp_weak.
Plot: `plots/phase5/s0021_285082902.png`

## TIC 219116234 (Sector 21) — **doubtful**

Phase 4: **submit** → Phase 5: **maybe** (Phase 5 doubtful: deep_eclipse)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (7511, -1145, -507, 1055 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 2160 ppm (0.93× Phase 3), significance 5.8σ (Phase 3: 14.0). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 405; radius 0.84 R☉ from Gaia FLAME). The companion is then about 0.40 R_J (TIC radius gave 0.40).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.6 km/s over 164 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)): for P = 1695.09 d the shape needs 0.4× the star's density on a circular orbit, i.e. an eccentricity of at least ~0.27 (≥ 0.13 at 1σ); for P = 339.02 d the shape needs 0.1× the star's density on a circular orbit, i.e. an eccentricity of at least ~0.67 (≥ 0.58 at 1σ). The shape mildly favours P = 1695.09 d.
- **Other TESS sectors.** **Another TESS sector shows a much deeper eclipse on this star**: 1.3 % deep in Sector 78 (BTJD 3436.34, SNR 23), against 0.23 % for this dip. The star is an eclipsing binary; this dip is probably its shallower (secondary) eclipse.

**Verdict: doubtful** — failed: deep_eclipse.
Plot: `plots/phase5/s0021_219116234.png`

## TIC 284625888 (Sector 21) — **doubtful**

Phase 4: **submit** → Phase 5: **maybe** (Phase 5 doubtful: gp)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (3716, 4196, 4047, 5934 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 2018 ppm (0.77× Phase 3), significance 4.7σ (Phase 3: 10.4). **The dip is not robust to the choice of noise model.**
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 344); radius 0.92 R☉ (Gaia FLAME), so the companion is about 0.46 R_J (5 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.9 km/s over 29 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~1820 d (1006–5726); 98% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: gp.
Plot: `plots/phase5/s0021_284625888.png`

## TIC 1044288 (Sector 21) — **strong**

Phase 4: **submit** → Phase 5: **submit**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (1220, 1994, 2190, 2377 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 4624 ppm (0.95× Phase 3), significance 7.2σ (Phase 3: 7.6).
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 104); radius 0.63 R☉ (Gaia FLAME), so the companion is about 0.43 R_J (5 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.6 km/s over 18 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~34 d (16–161); 41% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: strong** — every check passed.
Plot: `plots/phase5/s0021_1044288.png`

## TIC 308015761 (Sector 21) — **doubtful**

Phase 4: **submit** → Phase 5: **maybe** (Phase 5 doubtful: aperture, gp)

- **Aperture test.** **The depth changes with aperture size** (3309, 3514, 5023, 9851 ppm in 1, 1.5, 2 and 3 px apertures; worst 5.2σ from what a signal on this star would give), and fits neighbour TIC 308015767 1.4 px away better. The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 9116 ppm (1.92× Phase 3), significance 2.1σ (Phase 3: 10.3). **The dip is not robust to the choice of noise model.**
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 123); radius 0.97 R☉ (Gaia FLAME), so the companion is about 0.65 R_J (7 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 2.2 km/s over 16 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~21 d (12–56); 49% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture, gp.
Plot: `plots/phase5/s0021_308015761.png`

## TIC 21442437 (Sector 21) — **plausible**

Phase 4: **submit** → Phase 5: **submit**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (2886, 2005, 808, 595 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1963 ppm (1.02× Phase 3), significance 6.8σ (Phase 3: 7.6). It survives, but only moderately.
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 195); radius 1.07 R☉ (Gaia FLAME), so the companion is about 0.46 R_J (5 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.5 km/s over 20 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~11 d (6–36); 15% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: gp_weak.
Plot: `plots/phase5/s0021_21442437.png`

## TIC 288676782 (Sector 21) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (1636, 1092, 1462, -393 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1554 ppm (1.09× Phase 3), significance 6.6σ (Phase 3: 9.7). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 462; radius 1.96 R☉ from Gaia FLAME). The companion is then about 0.72 R_J (TIC radius gave 0.73).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.5 km/s over 16 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~861 d (564–1858); 54% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: gp_weak, evolved.
Plot: `plots/phase5/s0021_288676782.png`

## TIC 82203007 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: aperture)

- **Aperture test.** **The depth changes with aperture size** (666, 1053, 1683, 2678 ppm in 1, 1.5, 2 and 3 px apertures; worst 11.0σ from what a signal on this star would give). The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1550 ppm (1.03× Phase 3), significance 8.6σ (Phase 3: 12.0).
- **Stellar check.** Gaia DR3 has no evolutionary parameters for this star; its TIC radius (1.60 R☉) is consistent with a dwarf, giving a companion of about 0.60 R_J (7 R⊕). Not independently confirmed.
- **Gaia binarity.** Gaia's astrometry is noisier than for a single star (RV error 0.6 km/s over 28 transits: no significant scatter; RUWE 2.46), but the image shape and velocities show no companion.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the dip's shape implies a circular-orbit period of ~427 d (285–724); 77% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture.
Plot: `plots/phase5/s0021_82203007.png`

## TIC 24004731 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: aperture, binarity)

- **Aperture test.** **The depth changes with aperture size** (5816, 1475, 415, -1565 ppm in 1, 1.5, 2 and 3 px apertures; worst 10.4σ from what a signal on this star would give). The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 2343 ppm (1.00× Phase 3), significance 7.9σ (Phase 3: 10.0).
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 359); radius 1.11 R☉ (Gaia FLAME), so the companion is about 0.52 R_J (6 R⊕).
- **Gaia binarity.** **Gaia's radial velocities swing by tens of km/s: a companion star on a short orbit, which could itself cause the dip.** radial velocity varies (p = 0.0e+00, renormalised GOF 8.3, amplitude 11.9 km/s); RUWE 5.13.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~713 d (571–1322); 95% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture, binarity.
Plot: `plots/phase5/s0021_24004731.png`

## TIC 91844695 (Sector 21) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (2989, 394, -25, -404 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1225 ppm (0.95× Phase 3), significance 7.8σ (Phase 3: 10.1).
- **Stellar check.** Gaia DR3 has no evolutionary parameters for this star; its TIC radius (0.87 R☉) is consistent with a dwarf, giving a companion of about 0.30 R_J (3 R⊕). Not independently confirmed.
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (no Gaia DR3 source).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the dip's shape implies a circular-orbit period of ~263 d (170–524); 60% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: aperture_marginal, star_unverified.
Plot: `plots/phase5/s0021_91844695.png`

## TIC 458423861 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: aperture)

- **Aperture test.** **The depth changes with aperture size** (1615, 2656, 4902, 9037 ppm in 1, 1.5, 2 and 3 px apertures; worst 3.4σ from what a signal on this star would give), and fits neighbour TIC 458423860 2.0 px away better. The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 3401 ppm (0.97× Phase 3), significance 6.1σ (Phase 3: 11.7). It survives, but only moderately.
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 292); radius 1.46 R☉ (Gaia FLAME), so the companion is about 0.84 R_J (9 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 2.2 km/s over 36 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~129 d (88–242); 76% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture.
Plot: `plots/phase5/s0021_458423861.png`

## TIC 230062978 (Sector 21) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (-1613, 1472, 2947, 2948 ppm), as expected if the dip comes from this star (TIC 1102535496 lies within 1 px and cannot be separated this way).
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 2767 ppm (0.98× Phase 3), significance 5.7σ (Phase 3: 9.2). It survives, but only moderately.
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 276); radius 0.61 R☉ (Gaia FLAME), so the companion is about 0.32 R_J (4 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.8 km/s over 19 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the dip's shape implies a circular-orbit period of ~2270 d (1055–8046); 92% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: aperture_marginal, gp_weak.
Plot: `plots/phase5/s0021_230062978.png`

## TIC 233163483 (Sector 21) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (2165, 1552, 1551, 1257 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1424 ppm (0.92× Phase 3), significance 7.0σ (Phase 3: 9.0). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 452; radius 2.01 R☉ from Gaia FLAME). The companion is then about 0.77 R_J (TIC radius gave 0.75).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.7 km/s over 26 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME): of the 10 periods allowed by the second dip, 1 fit a near-circular orbit; the best is P = 17.3 d (eccentricity ≥ ~0.26).
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: gp_weak, evolved.
Plot: `plots/phase5/s0021_233163483.png`

## TIC 353932738 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: aperture)

- **Aperture test.** **The depth changes with aperture size** (6375, 4397, 538, 1226 ppm in 1, 1.5, 2 and 3 px apertures; worst 3.7σ from what a signal on this star would give). The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 2178 ppm (1.01× Phase 3), significance 6.1σ (Phase 3: 11.0). It survives, but only moderately.
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 247); radius 1.62 R☉ (Gaia FLAME), so the companion is about 0.73 R_J (8 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.7 km/s over 26 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~856 d (387–8329); 100% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture.
Plot: `plots/phase5/s0021_353932738.png`

## TIC 285007702 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: aperture)

- **Aperture test.** **The depth changes with aperture size** (3045, 2001, 1696, 5179 ppm in 1, 1.5, 2 and 3 px apertures; worst 8.4σ from what a signal on this star would give), and fits neighbour TIC 285007703 3.1 px away better. The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1745 ppm (0.98× Phase 3), significance 8.8σ (Phase 3: 17.7).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 432; radius 1.48 R☉ from Gaia FLAME). The companion is then about 0.61 R_J (TIC radius gave 0.65).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.4 km/s over 20 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~117 d (92–212); 73% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture.
Plot: `plots/phase5/s0021_285007702.png`

## TIC 441735282 (Sector 21) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (5390, 845, 1102, 733 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 470 ppm (1.05× Phase 3), significance 6.4σ (Phase 3: 11.6). It survives, but only moderately.
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 195); radius 0.90 R☉ (Gaia FLAME), so the companion is about 0.18 R_J (2 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.2 km/s over 12 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~3689 d (2610–8823); 99% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: gp_weak.
Plot: `plots/phase5/s0021_441735282.png`

## TIC 86177767 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: aperture)

- **Aperture test.** **The depth changes with aperture size** (2849, 2922, 971, 836 ppm in 1, 1.5, 2 and 3 px apertures; worst 4.9σ from what a signal on this star would give). The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 710 ppm (1.00× Phase 3), significance 5.7σ (Phase 3: 11.5). It survives, but only moderately.
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 316); radius 1.97 R☉ (Gaia FLAME), so the companion is about 0.51 R_J (6 R⊕).
- **Gaia binarity.** Gaia's astrometry is noisier than for a single star (RV error 0.3 km/s over 22 transits: no significant scatter; RUWE 7.35), but the image shape and velocities show no companion.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~2433 d (1109–10030); 99% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture.
Plot: `plots/phase5/s0021_86177767.png`

## TIC 385453267 (Sector 21) — **strong**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (11130, 11462, 11468, 11214 ppm), as expected if the dip comes from this star (TIC 900390889 lies within 1 px and cannot be separated this way).
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 13786 ppm (1.00× Phase 3), significance 16.2σ (Phase 3: 29.7).
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 189); radius 0.67 R☉ (Gaia FLAME), so the companion is about 0.76 R_J (9 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.1 km/s over 22 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~65 d (43–119); 56% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: strong** — every check passed.
Plot: `plots/phase5/s0021_385453267.png`

## TIC 219860619 (Sector 21) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (5933, 871, 944, 753 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 669 ppm (0.93× Phase 3), significance 5.6σ (Phase 3: 9.6). It survives, but only moderately.
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 265); radius 1.71 R☉ (Gaia FLAME), so the companion is about 0.44 R_J (5 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.4 km/s over 11 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~947 d (579–3551); 56% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: aperture_marginal, gp_weak.
Plot: `plots/phase5/s0021_219860619.png`

## TIC 125552729 (Sector 21) — **strong**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (5607, 5434, 5082, 4066 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 5810 ppm (0.96× Phase 3), significance 9.3σ (Phase 3: 13.6).
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 319); radius 1.05 R☉ (Gaia FLAME), so the companion is about 0.80 R_J (9 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.8 km/s over 27 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~65 d (52–103); 69% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: strong** — every check passed.
Plot: `plots/phase5/s0021_125552729.png`

## TIC 160333759 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: deep_eclipse)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (3964, 2665, 2482, 2306 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 3322 ppm (0.97× Phase 3), significance 7.5σ (Phase 3: 10.7).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 367; radius 1.20 R☉ from Gaia FLAME). The companion is then about 0.68 R_J (TIC radius gave 0.74).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 2.5 km/s over 19 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~195 d (126–581); 89% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** **Another TESS sector shows a much deeper eclipse on this star**: 1.4 % deep in Sector 48 (BTJD 2619.96, SNR 39), against 0.34 % for this dip. The star is an eclipsing binary; this dip is probably its shallower (secondary) eclipse.

**Verdict: doubtful** — failed: deep_eclipse.
Plot: `plots/phase5/s0021_160333759.png`

## TIC 286750182 (Sector 21) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (3310, 2926, 2253, 1740 ppm), as expected if the dip comes from this star (TIC 286750183 lies within 1 px and cannot be separated this way).
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 4847 ppm (0.96× Phase 3), significance 8.0σ (Phase 3: 14.5).
- **Stellar check.** The star is **subgiant** (FLAME evolutionary stage 325, but radius 2.08 R☉; radius 2.08 R☉ from Gaia FLAME). The companion is then about 1.44 R_J (TIC radius gave 1.19).
- **Gaia binarity.** Gaia's astrometry is noisier than for a single star (RV error 0.8 km/s over 32 transits: no significant scatter; RUWE 3.23), but the image shape and velocities show no companion.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~12 d (7–39); 2% of that range is still allowed by TESS's other observations. An allowed period is reachable with a moderately eccentric orbit (e ≤ 0.5).
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: evolved, ruwe, eccentric.
Plot: `plots/phase5/s0021_286750182.png`

## TIC 86379954 (Sector 21) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (5385, 6360, 9220, 11220 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 4548 ppm (1.10× Phase 3), significance 5.1σ (Phase 3: 12.0). It survives, but only moderately.
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 303); radius 1.76 R☉ (Gaia FLAME), so the companion is about 1.10 R_J (12 R⊕).
- **Gaia binarity.** Gaia's astrometry is noisier than for a single star (RV error 1.9 km/s over 29 transits: no significant scatter; RUWE 2.57), but the image shape and velocities show no companion.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~733 d (280–2552); 91% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: aperture_marginal, gp_weak, ruwe.
Plot: `plots/phase5/s0021_86379954.png`

## TIC 117409179 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: binarity, deep_eclipse)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (5735, 5712, 5251, 5144 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 5867 ppm (0.96× Phase 3), significance 11.3σ (Phase 3: 17.8).
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 295); radius 1.74 R☉ (Gaia FLAME), so the companion is about 1.32 R_J (15 R⊕).
- **Gaia binarity.** **Gaia's radial velocities swing by tens of km/s: a companion star on a short orbit, which could itself cause the dip.** radial velocity varies (p = 3.7e-14, renormalised GOF 6.4, amplitude 51.1 km/s); RUWE 2.04.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME): of the 16 periods allowed by the second dip, 1 fit a near-circular orbit; the best is P = 39.7 d (eccentricity ≥ ~0.26).
- **Other TESS sectors.** **Another TESS sector shows a much deeper eclipse on this star**: 5.4 % deep in Sector 72 (BTJD 3275.79, SNR 125), against 0.61 % for this dip. The star is an eclipsing binary; this dip is probably its shallower (secondary) eclipse.

**Verdict: doubtful** — failed: binarity, deep_eclipse.
Plot: `plots/phase5/s0021_117409179.png`

## TIC 417913356 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: gp)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (3595, 3046, 3008, 3069 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 3460 ppm (0.77× Phase 3), significance 3.4σ (Phase 3: 11.9). **The dip is not robust to the choice of noise model.**
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 339); radius 1.03 R☉ (Gaia FLAME), so the companion is about 0.68 R_J (8 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.5 km/s over 30 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~1963 d (735–6390); 95% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: gp.
Plot: `plots/phase5/s0021_417913356.png`

## TIC 88115338 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: aperture)

- **Aperture test.** **The depth changes with aperture size** (4902, 2632, 472, -359 ppm in 1, 1.5, 2 and 3 px apertures; worst 4.3σ from what a signal on this star would give). The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 2314 ppm (0.98× Phase 3), significance 7.2σ (Phase 3: 8.7).
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 325); radius 1.22 R☉ (Gaia FLAME), so the companion is about 0.58 R_J (6 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.9 km/s over 16 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~46 d (25–202); 36% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture.
Plot: `plots/phase5/s0021_88115338.png`

## TIC 302330626 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: aperture)

- **Aperture test.** **The depth changes with aperture size** (4013, 8155, 10117, 16938 ppm in 1, 1.5, 2 and 3 px apertures; worst 6.2σ from what a signal on this star would give). The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 3277 ppm (0.94× Phase 3), significance 5.4σ (Phase 3: 9.9). It survives, but only moderately.
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 129); radius 0.96 R☉ (Gaia FLAME), so the companion is about 0.55 R_J (6 R⊕).
- **Gaia binarity.** Gaia sees a close, partly resolved companion star (13 % of Gaia scans see a double image (ipd_frac_multi_peak); RV error 1.2 km/s over 38 transits: no significant scatter). Planets do orbit such stars, but its light dilutes the dip, so the object may be larger than estimated.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~211 d (84–1817); 63% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture.
Plot: `plots/phase5/s0021_302330626.png`

## TIC 147560583 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: deep_eclipse)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (12263, 11580, 11124, 10406 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 13043 ppm (1.00× Phase 3), significance 13.1σ (Phase 3: 27.4).
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 268); radius 0.99 R☉ (Gaia FLAME), so the companion is about 1.10 R_J (12 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 2.8 km/s over 28 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME): of the 19 periods allowed by the second dip, 12 fit a near-circular orbit; the best is P = 74.3 d (eccentricity ≥ ~0.08).
- **Other TESS sectors.** **Another TESS sector shows a much deeper eclipse on this star**: 8.9 % deep in Sector 14 (BTJD 1702.45, SNR 176), against 1.31 % for this dip. The star is an eclipsing binary; this dip is probably its shallower (secondary) eclipse.

**Verdict: doubtful** — failed: deep_eclipse.
Plot: `plots/phase5/s0021_147560583.png`

## TIC 236874984 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: deep_eclipse)

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (953, 1298, 1340, 1142 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 1853 ppm (0.95× Phase 3), significance 7.1σ (Phase 3: 9.7).
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (FLAME evolutionary stage 276); radius 1.06 R☉ (Gaia FLAME), so the companion is about 0.46 R_J (5 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 0.5 km/s over 29 transits: no significant scatter).
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (Gaia FLAME), the dip's shape implies a circular-orbit period of ~436 d (136–4967); 44% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** **Another TESS sector shows a much deeper eclipse on this star**: 5.0 % deep in Sector 78 (BTJD 3434.89, SNR 254), against 0.20 % for this dip. The star is an eclipsing binary; this dip is probably its shallower (secondary) eclipse.

**Verdict: doubtful** — failed: deep_eclipse.
Plot: `plots/phase5/s0021_236874984.png`

## TIC 165986822 (Sector 21) — **doubtful**

Phase 4: **maybe** → Phase 5: **drop** (Phase 5 doubtful: aperture)

- **Aperture test.** **The depth changes with aperture size** (5709, 11531, 14155, 12475 ppm in 1, 1.5, 2 and 3 px apertures; worst 4.6σ from what a signal on this star would give), and fits neighbour TIC 165986821 1.3 px away better. The dip may not come from this star, or may be instrumental.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 13466 ppm (1.97× Phase 3), significance 6.7σ (Phase 3: 9.7). It survives, but only moderately.
- **Stellar check.** The star is **subgiant** (TIC radius 2.26 R☉ only; Gaia DR3 has no evolutionary parameters for this star; radius 2.26 R☉ from TIC). The companion is then about 1.82 R_J (TIC radius gave 1.82).
- **Gaia binarity.** Gaia's astrometry is noisier than for a single star (RV error 3.4 km/s over 21 transits: no significant scatter; RUWE 1.43), but the image shape and velocities show no companion.
- **Variability.** Not listed as variable by Gaia DR3 or the AAVSO VSX catalogue.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the dip's shape implies a circular-orbit period of ~22 d (11–65); 50% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: doubtful** — failed: aperture.
Plot: `plots/phase5/s0021_165986822.png`

## TIC 202441153 (Sector 21) — **plausible**

Phase 4: **maybe** → Phase 5: **maybe**

- **Aperture test.** The depth is the same, within the noise, in small and large apertures (4116, 4605, 4621, 4696 ppm), as expected if the dip comes from this star.
- **Independent reprocessing.** Re-analysed with a different noise model (a Gaussian process fitted together with the transit): depth 5846 ppm (1.03× Phase 3), significance 8.7σ (Phase 3: 15.7).
- **Stellar check.** Gaia confirms a main-sequence (dwarf) star (GSP-Phot log g 4.64); radius 0.69 R☉ (Gaia FLAME), so the companion is about 0.51 R_J (6 R⊕).
- **Gaia binarity.** No sign of a close companion in Gaia's image shape or velocities (RV error 1.5 km/s over 25 transits: no significant scatter).
- **Variability.** Catalogued as variable: Gaia DR3 lists it as variable (SOLAR_LIKE, score 0.46); VSX Gaia DR3 1644940730266100992 (ROT) at 0″.
- **Physical consistency.** Fitted with this star's density (TIC (±25 %)), the dip's shape implies a circular-orbit period of ~839 d (331–3368); 78% of that range is still allowed by TESS's other observations.
- **Other TESS sectors.** No much deeper eclipse in any other TESS sector of this star.

**Verdict: plausible** — no check failed, but some raise minor concerns: variable.
Plot: `plots/phase5/s0021_202441153.png`
