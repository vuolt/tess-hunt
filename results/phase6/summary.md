# Phase 6: likely false positives among existing TESS candidates

`python scripts/phase6.py gaia | lc | report`. Tables: `likely_false_positives.csv`, `validation_gaia.csv`, `validation_lightcurve_checks.csv`, `gaia_substellar_companions.csv`; plot: `plots/phase6/validation_rates.png`.

## What was checked

- **Candidates:** 9,615 TOIs and CTOIs without a final disposition. 8,333 of them have a period, which the checks need. The validation sets come from the same tables: 1,403 confirmed or known planets (CP, KP) and 1,410 known false positives (FP, FA).
- **Gaia orbit check:** applied to every candidate with a period. It used 12,432 Gaia DR3 source IDs, queried in batches.
- **Light-curve checks:** 1,090 candidates had usable TESS light curves in at least one sector. These were all the Gaia-matched candidates plus random samples of about 300 from each group.

## Validation: how often each check flags known planets

| check | planets flagged | known false positives flagged | used? |
|---|---|---|---|
| Gaia orbit (period match + stellar companion or eclipsing binary) | 1/1403 (0.1%) | 68/1410 (4.8%) | yes |
| odd_even: 3 sigma | 56/298 (18.8%) | 75/273 (27.5%) | no |
| odd_even: 5 sigma | 18/298 (6.0%) | 43/273 (15.8%) | no |
| odd_even: 5 sigma and >= 20 % different | 10/298 (3.4%) | 33/273 (12.1%) | **yes** |
| secondary: SNR >= 7 | 33/298 (11.1%) | 90/273 (33.0%) | no |
| secondary: SNR >= 7 and >= 10 % of the transit depth | 25/298 (8.4%) | 81/273 (29.7%) | no |
| secondary: SNR >= 10 and >= 10 % of the transit depth | 15/298 (5.0%) | 64/273 (23.4%) | no |
| centroid: 3 sigma | 116/298 (38.9%) | 166/273 (60.8%) | no |
| centroid: 5 sigma | 87/298 (29.2%) | 150/273 (54.9%) | no |
| centroid: 5 sigma and source >= 1 px away | 0/298 (0.0%) | 21/273 (7.7%) | **yes** |

The Gaia check found the Gaia orbit of the transiting object itself, with a substellar companion mass, on 3 known planets (WASP-18 b, the brown dwarf TOI-503 b and others). Such matches are *not* counted as false positives. If Gaia periods were unrelated to the transits, shuffling them among the 472 candidates that have a Gaia solution gives 4.7 ± 2.1 matches by chance.

**Not used, because they flag too many known planets:** secondary.

## Result

**198 unresolved candidates show evidence of being false positives:** 175 high, 14 medium and 9 low confidence.

How they were flagged (a candidate can be flagged by more than one check):

- 181 by the Gaia orbit check: Gaia sees a companion star (or an eclipsing binary) on the transit's period.
- 43 by the odd/even check: alternate transits have different depths: two stars eclipsing at twice the period.
- 10 by the centroid check: the light dims off-centre: the eclipse is on a neighbouring star.

In the random sample of 269 unresolved candidates, 17 (6%) were flagged by a light-curve check. If the sample is representative, that is about 527 of the 8,333 unresolved candidates with periods; only 269 were examined here.

**How reliable the flags are.** Known planets were flagged 3.4% by the odd/even check, 0.0% by the centroid check and 0.07% by the Gaia check. A candidate flagged only by a light-curve check could still be a planet, so a low-confidence flag means *look again*, not *false positive*.

Known planets flagged by the odd/even check: TOI 1136.03, TOI 2076.01, TOI 2207.01, TOI 2449.01, TOI 3353.01, TOI 396.01, TOI 4127.01, TOI 451.02, TOI 6551.01, TOI 6647.01. Some are known to have transit-timing variations or young, spotted host stars (e.g. TOI-1136, TOI-2076, TOI-451), where a fixed ephemeris catches some transits only partly. A candidate flagged only by this check should be checked for timing variations first.

## The high-confidence flags

| candidate | TIC | disposition | period (d) | checks | evidence |
|---|---|---|---|---|---|
| CTOI 102625324.01 | 102625324 | PC | 9.1764 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 9.176 d; companion 0.32 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 57.0 sigma (1704 vs 13526 ppm) |
| CTOI 1030389830.01 | 1030389830 | PC | 0.6877 | gaia_orbit | Gaia eclipsing binary (photometric): P = 0.6877 d |
| CTOI 115115136.01 | 115115136 | PC | 3.9708 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 3.961 d; companion 0.24 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 32.3 sigma (7202 vs 11185 ppm) |
| CTOI 11799874.01 | 11799874 | PC | 4.2581 | gaia_orbit | Gaia NSS SB1 orbit P = 4.258 d; companion 0.14 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 118389230.01 | 118389230 | PC | 3.7575 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 3.758 d; companion 0.16 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 13.0 sigma (9455 vs 5871 ppm) |
| CTOI 121420805.01 | 121420805 | PC | 1.7201 | gaia_orbit | Gaia eclipsing binary (photometric): P = 1.7169 d |
| CTOI 121420805.02 | 121420805 | PC | 1.7181 | gaia_orbit | Gaia eclipsing binary (photometric): P = 1.7169 d |
| CTOI 126945917.02 | 126945917 | PC | 1.3424 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 1.341 d; companion 0.16 M_sun (SB1 mass function (sin i = 1)); Gaia eclipsing binary (photometric): P = 1.3411 d ; odd and even transits differ by 367.1 sigma (81552 vs 131348 ppm) |
| CTOI 135230723.01 | 135230723 | PC | 6.7760 | gaia_orbit | Gaia NSS SB1 orbit P = 6.779 d; companion 0.16 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 136886352.01 | 136886352 | PC | 2.7223 | gaia_orbit | Gaia NSS SB1 orbit P = 2.722 d; companion 0.56 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 138854818.01 | 138854818 | PC | 5.7890 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 5.788 d; companion 0.15 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 10.1 sigma (2638 vs 4704 ppm) |
| CTOI 138941853.01 | 138941853 | PC | 11.5902 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 11.602 d; companion 0.20 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 7.1 sigma (8161 vs -1208 ppm) |
| CTOI 140736015.01 | 140736015 | PC | 7.3135 | gaia_orbit, odd_even | Gaia NSS SB2 orbit P = 7.313 d; companion 1.45 M_sun (SB2 mass ratio (K1/K2)) ; odd and even transits differ by 71.0 sigma (11712 vs 6918 ppm) |
| CTOI 141268467.01 | 141268467 | PC | 4.1806 | gaia_orbit | Gaia NSS SB2 orbit P = 4.181 d; companion 1.17 M_sun (SB2 mass ratio (K1/K2)) |
| CTOI 141268467.02 | 141268467 | PC | 4.1750 | gaia_orbit | Gaia NSS SB2 orbit P = 4.181 d; companion 1.17 M_sun (SB2 mass ratio (K1/K2)) |
| CTOI 144462697.01 | 144462697 | PC | 1.2333 | gaia_orbit | Gaia eclipsing binary (photometric): P = 1.2340 d |
| CTOI 144462697.02 | 144462697 | PC | 1.2354 | gaia_orbit | Gaia eclipsing binary (photometric): P = 1.2340 d |
| CTOI 152223725.01 | 152223725 | PC | 4.4319 | gaia_orbit | Gaia eclipsing binary (photometric): P = 4.4365 d |
| CTOI 152223725.02 | 152223725 | PC | 4.4403 | gaia_orbit | Gaia eclipsing binary (photometric): P = 4.4365 d |
| CTOI 153734545.01 | 153734545 | PC | 3.2196 | gaia_orbit | Gaia NSS SB2 orbit P = 6.438 d (transit period = half Gaia's); companion 1.69 M_sun (SB2 mass ratio (K1/K2)) |
| CTOI 153735144.01 | 153735144 | PC | 3.7785 | gaia_orbit | Gaia NSS SB1 orbit P = 3.779 d; companion 0.29 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 155916287.01 | 155916287 | PC | 6.0531 | gaia_orbit | Gaia NSS SB1 orbit P = 6.053 d; companion 0.71 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 167005326.01 | 167005326 | PC | 0.7204 | gaia_orbit, centroid | Gaia eclipsing binary (photometric): P = 0.7204 d ; the star's image shifts during transit (15.6 sigma; the dimming source would be ~4131 arcsec away) |
| CTOI 176871438.01 | 176871438 | PC | 15.1383 | gaia_orbit, centroid | Gaia NSS SB1 orbit P = 15.135 d; companion 0.22 M_sun (SB1 mass function (sin i = 1)) ; the star's image shifts during transit (15.0 sigma; the dimming source would be ~2194 arcsec away) |
| CTOI 190736332.01 | 190736332 | PC | 2.4628 | gaia_orbit | Gaia NSS SB1 orbit P = 2.463 d; companion 0.19 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 192831602.01 | 192831602 | PC | 19.5799 | gaia_orbit | Gaia NSS SB1 orbit P = 39.144 d (transit period = half Gaia's); companion 0.39 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 198358825.01 | 198358825 | PC | 23.4995 | gaia_orbit | Gaia NSS SB1 orbit P = 23.519 d; companion 0.31 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 204497617.01 | 204497617 | PC | 1.5903 | gaia_orbit | Gaia eclipsing binary (photometric): P = 0.7947 d (transit period = 2x Gaia's) |
| CTOI 207080350.01 | 207080350 | PC | 9.1056 | gaia_orbit, odd_even | Gaia NSS SB2 orbit P = 9.108 d; companion 1.46 M_sun (SB2 mass ratio (K1/K2)); Gaia eclipsing binary (photometric): P = 9.1070 d ; odd and even transits differ by 1964.1 sigma (131555 vs 6630 ppm) |
| CTOI 214299966.01 | 214299966 | PC | 3.0681 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 3.067 d; companion 0.27 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 22.5 sigma (3821 vs 5626 ppm) |
| CTOI 21818238.01 | 21818238 | PC | 10.4429 | gaia_orbit | Gaia NSS SB1 orbit P = 10.448 d; companion 0.26 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 219201178.01 | 219201178 | PC | 2.4168 | gaia_orbit, centroid | Gaia NSS SB1 orbit P = 2.417 d; companion 0.57 M_sun (SB1 mass function (sin i = 1)) ; the star's image shifts during transit (20.6 sigma; the dimming source would be ~63 arcsec away) |
| CTOI 219308281.01 | 219308281 | PC | 2.8604 | gaia_orbit | Gaia eclipsing binary (photometric): P = 2.8574 d |
| CTOI 219308281.02 | 219308281 | PC | 2.8514 | gaia_orbit | Gaia eclipsing binary (photometric): P = 2.8574 d |
| CTOI 219322317.01 | 219322317 | PC | 8.5460 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 8.544 d; companion 0.33 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 15.7 sigma (10643 vs 7422 ppm) |
| CTOI 220523550.01 | 220523550 | PC | 13.5236 | gaia_orbit | Gaia NSS SB1 orbit P = 13.491 d; companion 0.28 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 220523550.02 | 220523550 | PC | 13.4951 | gaia_orbit | Gaia NSS SB1 orbit P = 13.491 d; companion 0.28 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 220567625.01 | 220567625 | PC | 19.4444 | gaia_orbit | Gaia NSS SB1 orbit P = 19.427 d; companion 0.50 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 225336886.01 | 225336886 | PC | 2.1708 | gaia_orbit | Gaia eclipsing binary (photometric): P = 2.1712 d |
| CTOI 231714759.02 | 231714759 | PC | 5.4597 | gaia_orbit, odd_even | Gaia eclipsing binary (photometric): P = 5.4612 d ; odd and even transits differ by 490.9 sigma (246391 vs 148198 ppm) |
| CTOI 233126177.01 | 233126177 | PC | 65.9450 | gaia_orbit | Gaia NSS SB1 orbit P = 66.441 d; companion 0.17 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 235009317.01 | 235009317 | PC | 7.4583 | gaia_orbit | Gaia NSS SB1 orbit P = 7.457 d; companion 0.59 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 235009317.02 | 235009317 | PC | 7.4528 | gaia_orbit | Gaia NSS SB1 orbit P = 7.457 d; companion 0.59 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 237332597.01 | 237332597 | PC | 2.9299 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 2.927 d; companion 0.50 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 410.9 sigma (22866 vs 12983 ppm) |
| CTOI 237332597.02 | 237332597 | PC | 2.9306 | gaia_orbit | Gaia NSS SB1 orbit P = 2.927 d; companion 0.50 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 237342298.01 | 237342298 | PC | 2.6921 | gaia_orbit | Gaia NSS SB1 orbit P = 2.692 d; companion 0.23 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 25226885.01 | 25226885 | PC | 1.6313 | gaia_orbit, odd_even | Gaia eclipsing binary (photometric): P = 1.6321 d ; odd and even transits differ by 106.7 sigma (71876 vs 51885 ppm) |
| CTOI 25226885.02 | 25226885 | PC | 1.6306 | gaia_orbit, odd_even | Gaia eclipsing binary (photometric): P = 1.6321 d ; odd and even transits differ by 23.4 sigma (12435 vs 7641 ppm) |
| CTOI 257483992.01 | 257483992 | PC | 8.7867 | gaia_orbit | Gaia NSS SB1 orbit P = 8.788 d; companion 0.21 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 264485594.01 | 264485594 | PC | 7.1700 | gaia_orbit | Gaia NSS SB1 orbit P = 7.168 d; companion 0.65 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 265489698.01 | 265489698 | PC | 1.1562 | gaia_orbit, centroid | Gaia eclipsing binary (photometric): P = 1.1555 d ; the star's image shifts during transit (5.2 sigma; the dimming source would be ~10117 arcsec away) |
| CTOI 265489698.02 | 265489698 | PC | 1.1535 | gaia_orbit | Gaia eclipsing binary (photometric): P = 1.1555 d |
| CTOI 266996791.01 | 266996791 | PC | 2.1264 | gaia_orbit, odd_even | Gaia eclipsing binary (photometric): P = 2.1237 d ; odd and even transits differ by 1129.2 sigma (224239 vs 486544 ppm) |
| CTOI 266996791.02 | 266996791 | PC | 2.1229 | gaia_orbit | Gaia eclipsing binary (photometric): P = 2.1237 d |
| CTOI 268290940.01 | 268290940 | PC | 9.3067 | gaia_orbit, odd_even | Gaia NSS SB2 orbit P = 9.305 d; companion 1.24 M_sun (SB2 mass ratio (K1/K2)) ; odd and even transits differ by 275.2 sigma (3855 vs 16163 ppm) |
| CTOI 272357134.01 | 272357134 | PC | 4.1875 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 4.197 d; companion 0.24 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 5.8 sigma (177 vs 131 ppm) |
| CTOI 272357134.02 | 272357134 | PC | 4.1958 | gaia_orbit | Gaia NSS SB1 orbit P = 4.197 d; companion 0.24 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 273792220.01 | 273792220 | PC | 1.9868 | gaia_orbit | Gaia NSS SB2C orbit P = 1.988 d; companion 1.26 M_sun (SB2 mass ratio (K1/K2)) |
| CTOI 273792220.02 | 273792220 | PC | 1.9868 | gaia_orbit | Gaia NSS SB2C orbit P = 1.988 d; companion 1.26 M_sun (SB2 mass ratio (K1/K2)) |
| CTOI 277649549.01 | 277649549 | PC | 14.4882 | gaia_orbit | Gaia NSS SB2 orbit P = 14.491 d; companion 1.38 M_sun (SB2 mass ratio (K1/K2)) |
| CTOI 277649549.02 | 277649549 | PC | 14.4903 | gaia_orbit | Gaia NSS SB2 orbit P = 14.491 d; companion 1.38 M_sun (SB2 mass ratio (K1/K2)) |
| CTOI 278683641.01 | 278683641 | PC | 4.2521 | gaia_orbit | Gaia NSS SB1 orbit P = 4.242 d; companion 0.55 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 279322914.01 | 279322914 | PC | 18.8689 | gaia_orbit | Gaia NSS SB1 orbit P = 18.874 d; companion 0.24 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 282502866.01 | 282502866 | PC | 7.4007 | gaia_orbit | Gaia NSS SB2 orbit P = 7.397 d; companion 1.35 M_sun (SB2 mass ratio (K1/K2)) |
| CTOI 284990199.01 | 284990199 | PC | 4.0931 | odd_even, centroid | odd and even transits differ by 11.1 sigma (-14202 vs -1771 ppm) ; the star's image shifts during transit (37.0 sigma; the dimming source would be ~68469 arcsec away) |
| CTOI 287351976.01 | 287351976 | PC | 1.9243 | gaia_orbit | Gaia NSS EclipsingBinary: P = 1.9219 d; Gaia eclipsing binary (photometric): P = 1.9219 d |
| CTOI 287351976.02 | 287351976 | PC | 1.9236 | gaia_orbit | Gaia NSS EclipsingBinary: P = 1.9219 d; Gaia eclipsing binary (photometric): P = 1.9219 d |
| CTOI 295176393.01 | 295176393 | PC | 27.2519 | gaia_orbit | Gaia NSS SB1 orbit P = 27.249 d; companion 0.53 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 300903537.01 | 300903537 | PC | 5.8236 | gaia_orbit | Gaia NSS SB1 orbit P = 5.822 d; companion 0.20 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 301980639.01 | 301980639 | PC | 4.7826 | gaia_orbit | Gaia NSS SB1 orbit P = 4.783 d; companion 0.30 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 306635088.01 | 306635088 | PC | 17.0004 | gaia_orbit | Gaia NSS SB1 orbit P = 17.006 d; companion 0.23 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 318357083.01 | 318357083 | PC | 3.2400 | gaia_orbit | Gaia NSS SB1 orbit P = 3.244 d; companion 0.15 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 318707367.01 | 318707367 | PC | 13.9757 | gaia_orbit | Gaia NSS SB1 orbit P = 6.987 d (transit period = 2x Gaia's); companion 0.43 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 32606889.01 | 32606889 | PC | 4.6869 | gaia_orbit | Gaia NSS SB1 orbit P = 4.687 d; companion 0.21 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 339502646.01 | 339502646 | PC | 15.6250 | odd_even, centroid | odd and even transits differ by 6.8 sigma (-616 vs 2694 ppm) ; the star's image shifts during transit (7.5 sigma; the dimming source would be ~49 arcsec away) |
| CTOI 344087362.01 | 344087362 | PC | 14.4959 | gaia_orbit | Gaia NSS SB1 orbit P = 14.496 d; companion 0.23 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 350479101.02 | 350479101 | PC | 11.9229 | gaia_orbit, odd_even | Gaia NSS SB2C orbit P = 11.918 d; companion 1.40 M_sun (SB2 mass ratio (K1/K2)) ; odd and even transits differ by 42.9 sigma (3990 vs 3114 ppm) |
| CTOI 350480660.01 | 350480660 | PC | 4.4630 | gaia_orbit | Gaia NSS SB1 orbit P = 4.464 d; companion 0.16 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 356435413.01 | 356435413 | PC | 27.4796 | gaia_orbit | Gaia NSS SB1 orbit P = 27.479 d; companion 0.68 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 368803791.01 | 368803791 | PC | 15.5260 | gaia_orbit | Gaia NSS SB1 orbit P = 7.761 d (transit period = 2x Gaia's); companion 0.19 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 372909935.01 | 372909935 | PC | 3.6125 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 3.611 d; companion 0.31 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 9.5 sigma (1338 vs 960 ppm) |
| CTOI 375477302.01 | 375477302 | PC | 3.0786 | gaia_orbit | Gaia NSS SB1 orbit P = 3.072 d; companion 0.16 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 381543300.01 | 381543300 | PC | 3.8557 | gaia_orbit | Gaia NSS SB1 orbit P = 3.855 d; companion 0.26 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 38850860.01 | 38850860 | PC | 8.4171 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 8.413 d; companion 0.14 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 56.4 sigma (15612 vs 3481 ppm) |
| CTOI 39018208.01 | 39018208 | PC | 6.6458 | gaia_orbit | Gaia NSS SB1 orbit P = 6.650 d; companion 0.43 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 39018208.02 | 39018208 | PC | 6.6604 | gaia_orbit | Gaia NSS SB1 orbit P = 6.650 d; companion 0.43 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 406526291.01 | 406526291 | PC | 3.5100 | gaia_orbit | Gaia NSS SB1 orbit P = 3.510 d; companion 0.14 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 416455230.01 | 416455230 | PC | 7.7925 | gaia_orbit | Gaia NSS SB1 orbit P = 7.799 d; companion 0.13 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 41896812.01 | 41896812 | PC | 1.2590 | gaia_orbit | Gaia eclipsing binary (photometric): P = 2.5216 d (transit period = half Gaia's) |
| CTOI 41896812.02 | 41896812 | PC | 1.2632 | gaia_orbit | Gaia eclipsing binary (photometric): P = 2.5216 d (transit period = half Gaia's) |
| CTOI 441422220.01 | 441422220 | PC | 3.7000 | gaia_orbit | Gaia NSS SB1 orbit P = 3.700 d; companion 0.67 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 441422220.02 | 441422220 | PC | 3.6993 | gaia_orbit | Gaia NSS SB1 orbit P = 3.700 d; companion 0.67 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 467113954.01 | 467113954 | PC | 7.5708 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 7.567 d; companion 0.31 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 67.7 sigma (30610 vs 18184 ppm) |
| CTOI 51238317.01 | 51238317 | PC | 14.2610 | gaia_orbit | Gaia NSS SB1 orbit P = 14.249 d; companion 0.20 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 52079677.01 | 52079677 | PC | 9.0347 | gaia_orbit | Gaia NSS SB1 orbit P = 9.028 d; companion 0.69 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 52079677.02 | 52079677 | PC | 9.0417 | gaia_orbit | Gaia NSS SB1 orbit P = 9.028 d; companion 0.69 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 55497281.01 | 55497281 | PC | 14.9028 | gaia_orbit | Gaia NSS SB2 orbit P = 14.901 d; companion 1.16 M_sun (SB2 mass ratio (K1/K2)) |
| CTOI 62573638.01 | 62573638 | PC | 4.7225 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 4.716 d; companion 0.34 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 58.7 sigma (1406 vs 2988 ppm) |
| CTOI 65557265.01 | 65557265 | PC | 1.6828 | gaia_orbit | Gaia eclipsing binary (photometric): P = 3.3657 d (transit period = half Gaia's) |
| CTOI 70663347.01 | 70663347 | PC | 8.7336 | gaia_orbit | Gaia NSS SB1 orbit P = 8.731 d; companion 0.19 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 90083037.01 | 90083037 | PC | 9.9415 | gaia_orbit | Gaia NSS SB1 orbit P = 9.940 d; companion 0.17 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 91369561.02 | 91369561 | PC | 3.9931 | gaia_orbit, odd_even | Gaia eclipsing binary (photometric): P = 3.9830 d ; odd and even transits differ by 14.8 sigma (39037 vs 48927 ppm) |
| CTOI 91369561.03 | 91369561 | PC | 3.9924 | gaia_orbit | Gaia eclipsing binary (photometric): P = 3.9830 d |
| CTOI 92349924.01 | 92349924 | PC | 8.8007 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 8.799 d; companion 0.63 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 11.4 sigma (4171 vs 5833 ppm) |
| CTOI 92594505.01 | 92594505 | PC | 22.8687 | gaia_orbit, odd_even | Gaia NSS SB2 orbit P = 11.437 d (transit period = 2x Gaia's); companion 1.31 M_sun (SB2 mass ratio (K1/K2)) ; odd and even transits differ by 689.0 sigma (58415 vs 28264 ppm) |
| CTOI 92866460.01 | 92866460 | PC | 19.1605 | gaia_orbit | Gaia NSS SB1 orbit P = 19.195 d; companion 0.54 M_sun (SB1 mass function (sin i = 1)) |
| CTOI 98672702.01 | 98672702 | PC | 15.1288 | gaia_orbit | Gaia NSS SB1 orbit P = 15.167 d; companion 0.35 M_sun (SB1 mass function (sin i = 1)) |
| TOI 1059.01 | 380783252 | PC | 9.4497 | gaia_orbit | Gaia NSS SB1 orbit P = 9.450 d; companion 0.23 M_sun (SB1 mass function (sin i = 1)) |
| TOI 1115.01 | 379286801 | APC | 4.4520 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 4.453 d; companion 0.60 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 14.8 sigma (1928 vs 716 ppm) |
| TOI 1119.01 | 261369656 | APC | 10.9566 | gaia_orbit | Gaia NSS SB1 orbit P = 10.959 d; companion 0.60 M_sun (SB1 mass function (sin i = 1)) |
| TOI 1153.01 | 154840461 | APC | 6.0361 | gaia_orbit | Gaia NSS SB1 orbit P = 6.036 d; companion 0.66 M_sun (SB1 mass function (sin i = 1)) |
| TOI 121.01 | 207081058 | APC | 14.7741 | gaia_orbit | Gaia NSS SB1 orbit P = 14.781 d; companion 0.11 M_sun (SB1 mass function (sin i = 1)) |
| TOI 1371.01 | 162130054 | APC | 5.0095 | gaia_orbit | Gaia NSS SB1 orbit P = 5.010 d; companion 0.18 M_sun (SB1 mass function (sin i = 1)) |
| TOI 1455.01 | 387259626 | APC | 3.6231 | gaia_orbit | Gaia NSS SB1 orbit P = 3.623 d; companion 0.15 M_sun (SB1 mass function (sin i = 1)) |
| TOI 1461.01 | 44631965 | APC | 3.5686 | gaia_orbit | Gaia NSS SB1 orbit P = 3.568 d; companion 0.17 M_sun (SB1 mass function (sin i = 1)) |
| TOI 1553.01 | 327341626 | APC | 5.9382 | gaia_orbit | Gaia NSS SB1 orbit P = 5.938 d; companion 0.38 M_sun (SB1 mass function (sin i = 1)) |
| TOI 162.01 | 99493790 | APC | 7.7647 | gaia_orbit | Gaia NSS SB1 orbit P = 7.765 d; companion 0.13 M_sun (SB1 mass function (sin i = 1)) |
| TOI 1951.01 | 356590734 | PC | 4.8294 | gaia_orbit | Gaia NSS SB1 orbit P = 4.830 d; companion 0.13 M_sun (SB1 mass function (sin i = 1)) |
| TOI 1974.01 | 295562522 | APC | 5.4382 | gaia_orbit | Gaia NSS SB1 orbit P = 5.438 d; companion 0.44 M_sun (SB1 mass function (sin i = 1)) |
| TOI 2038.01 | 235937532 | APC | 3.7321 | gaia_orbit | Gaia NSS SB1 orbit P = 3.732 d; companion 0.12 M_sun (SB1 mass function (sin i = 1)) |
| TOI 2159.01 | 270515566 | APC | 10.0510 | gaia_orbit | Gaia NSS SB1 orbit P = 10.051 d; companion 0.12 M_sun (SB1 mass function (sin i = 1)) |
| TOI 2349.01 | 405452527 | APC | 11.5738 | gaia_orbit | Gaia NSS SB1 orbit P = 11.562 d; companion 0.16 M_sun (SB1 mass function (sin i = 1)) |
| TOI 2366.01 | 39218269 | APC | 17.1870 | gaia_orbit | Gaia NSS SB1 orbit P = 17.175 d; companion 0.34 M_sun (SB1 mass function (sin i = 1)) |
| TOI 2606.01 | 355800238 | PC | 1.4509 | gaia_orbit | Gaia NSS SB1 orbit P = 1.451 d; companion 0.12 M_sun (SB1 mass function (sin i = 1)) |
| TOI 2819.01 | 387275908 | APC | 4.3906 | gaia_orbit | Gaia NSS SB1 orbit P = 4.391 d; companion 0.56 M_sun (SB1 mass function (sin i = 1)) |
| TOI 3006.01 | 386161311 | PC | 8.6135 | gaia_orbit | Gaia NSS SB1 orbit P = 8.611 d; companion 0.83 M_sun (SB1 mass function (sin i = 1)) |
| TOI 3010.01 | 357202877 | PC | 9.5538 | gaia_orbit | Gaia NSS SB1 orbit P = 9.543 d; companion 0.20 M_sun (SB1 mass function (sin i = 1)) |
| TOI 3142.01 | 387791274 | PC | 1.1537 | gaia_orbit | Gaia eclipsing binary (photometric): P = 2.3073 d (transit period = half Gaia's) |
| TOI 3148.01 | 323980895 | PC | 4.2814 | gaia_orbit | Gaia NSS SB1 orbit P = 4.281 d; companion 0.19 M_sun (SB1 mass function (sin i = 1)) |
| TOI 3191.01 | 342375120 | PC | 4.4368 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 8.879 d (transit period = half Gaia's); companion 1.11 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 5.1 sigma (6116 vs 8744 ppm) |
| TOI 3241.01 | 241326434 | PC | 6.8381 | gaia_orbit | Gaia NSS SB1 orbit P = 6.836 d; companion 0.16 M_sun (SB1 mass function (sin i = 1)) |
| TOI 3260.01 | 370138999 | PC | 11.4520 | gaia_orbit | Gaia NSS SB1 orbit P = 11.451 d; companion 0.32 M_sun (SB1 mass function (sin i = 1)) |
| TOI 3298.01 | 407772585 | PC | 3.0538 | gaia_orbit | Gaia NSS SB1 orbit P = 3.053 d; companion 0.37 M_sun (SB1 mass function (sin i = 1)) |
| TOI 3501.01 | 98957720 | APC | 15.3503 | gaia_orbit | Gaia NSS SB1 orbit P = 15.346 d; companion 0.19 M_sun (SB1 mass function (sin i = 1)) |
| TOI 3502.01 | 452281203 | PC | 3.4468 | gaia_orbit | Gaia NSS SB1 orbit P = 3.447 d; companion 0.27 M_sun (SB1 mass function (sin i = 1)) |
| TOI 3756.01 | 67300566 | APC | 4.4237 | gaia_orbit | Gaia NSS SB1 orbit P = 4.424 d; companion 0.62 M_sun (SB1 mass function (sin i = 1)) |
| TOI 3904.01 | 233462817 | APC | 10.1615 | gaia_orbit | Gaia NSS SB1 orbit P = 10.160 d; companion 0.54 M_sun (SB1 mass function (sin i = 1)) |
| TOI 4338.01 | 253434221 | APC | 7.4997 | gaia_orbit | Gaia NSS SB1 orbit P = 7.499 d; companion 0.29 M_sun (SB1 mass function (sin i = 1)) |
| TOI 4400.01 | 231630147 | PC | 7.0642 | gaia_orbit | Gaia NSS SB1 orbit P = 7.061 d; companion 0.14 M_sun (SB1 mass function (sin i = 1)) |
| TOI 446.01 | 1449640 | APC | 3.5018 | gaia_orbit | Gaia NSS SB1 orbit P = 3.502 d; companion 0.16 M_sun (SB1 mass function (sin i = 1)) |
| TOI 4657.01 | 434110695 | PC | 8.1157 | gaia_orbit | Gaia NSS SB1 orbit P = 8.111 d; companion 0.60 M_sun (SB1 mass function (sin i = 1)) |
| TOI 4792.01 | 265350422 | APC | 8.4482 | gaia_orbit | Gaia NSS SB1 orbit P = 8.448 d; companion 0.40 M_sun (SB1 mass function (sin i = 1)) |
| TOI 4807.01 | 174546820 | PC | 17.6869 | gaia_orbit | Gaia NSS SB1 orbit P = 17.689 d; companion 0.27 M_sun (SB1 mass function (sin i = 1)) |
| TOI 4973.01 | 397522142 | PC | 7.9486 | gaia_orbit | Gaia NSS SB1 orbit P = 7.948 d; companion 0.26 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5103.01 | 366443576 | APC | 4.9899 | gaia_orbit | Gaia NSS SB1 orbit P = 4.990 d; companion 0.74 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5149.01 | 286094277 | PC | 27.3715 | gaia_orbit | Gaia NSS SB1 orbit P = 27.388 d; companion 0.13 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5156.01 | 468889418 | PC | 22.8524 | gaia_orbit | Gaia NSS SB1 orbit P = 22.853 d; companion 0.20 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5197.01 | 219431828 | APC | 20.5773 | gaia_orbit | Gaia NSS SB1 orbit P = 20.573 d; companion 0.39 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5200.01 | 21101556 | PC | 17.4373 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 17.436 d; companion 0.52 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 5.1 sigma (1802 vs 1443 ppm) |
| TOI 5307.01 | 369477998 | APC | 5.2976 | gaia_orbit | Gaia NSS SB1 orbit P = 5.298 d; companion 0.62 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5372.01 | 371470493 | PC | 11.6271 | gaia_orbit | Gaia NSS SB1 orbit P = 11.630 d; companion 0.79 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5379.01 | 328325110 | APC | 12.7427 | gaia_orbit | Gaia NSS SB1 orbit P = 12.744 d; companion 0.50 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5394.01 | 61109252 | APC | 15.1935 | gaia_orbit | Gaia NSS SB1 orbit P = 15.195 d; companion 0.40 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5475.01 | 78601408 | APC | 10.0532 | gaia_orbit | Gaia NSS SB1 orbit P = 10.052 d; companion 0.84 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5677.01 | 287959034 | APC | 3.8481 | gaia_orbit | Gaia NSS SB1 orbit P = 3.848 d; companion 0.30 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5761.01 | 428703516 | APC | 19.5422 | gaia_orbit | Gaia NSS SB1 orbit P = 19.546 d; companion 0.58 M_sun (SB1 mass function (sin i = 1)) |
| TOI 5860.01 | 462372373 | APC | 5.2499 | gaia_orbit | Gaia NSS SB1 orbit P = 5.248 d; companion 0.27 M_sun (SB1 mass function (sin i = 1)) |
| TOI 6111.01 | 189322727 | APC | 5.8340 | gaia_orbit | Gaia NSS SB1 orbit P = 5.833 d; companion 0.70 M_sun (SB1 mass function (sin i = 1)) |
| TOI 6568.01 | 253126207 | PC | 17.4368 | gaia_orbit | Gaia NSS SB1 orbit P = 17.432 d; companion 0.13 M_sun (SB1 mass function (sin i = 1)) |
| TOI 6611.01 | 403368444 | PC | 6.5569 | gaia_orbit | Gaia NSS SB1 orbit P = 6.569 d; companion 0.12 M_sun (SB1 mass function (sin i = 1)) |
| TOI 668.01 | 102195674 | APC | 4.3787 | gaia_orbit | Gaia NSS SB1 orbit P = 4.378 d; companion 0.17 M_sun (SB1 mass function (sin i = 1)) |
| TOI 6806.01 | 144327080 | APC | 23.2491 | gaia_orbit | Gaia NSS SB1 orbit P = 23.229 d; companion 0.16 M_sun (SB1 mass function (sin i = 1)) |
| TOI 6936.01 | 63082902 | PC | 10.2986 | gaia_orbit | Gaia NSS SB1 orbit P = 10.298 d; companion 1.06 M_sun (SB1 mass function (sin i = 1)) |
| TOI 7093.01 | 62867537 | APC | 2.8372 | gaia_orbit | Gaia NSS SB1 orbit P = 2.837 d; companion 0.23 M_sun (SB1 mass function (sin i = 1)) |
| TOI 759.01 | 152147232 | APC | 4.2136 | gaia_orbit | Gaia NSS SB1 orbit P = 4.213 d; companion 0.11 M_sun (SB1 mass function (sin i = 1)) |
| TOI 7590.01 | 148322608 | PC | 7.3364 | gaia_orbit | Gaia NSS SB1 orbit P = 7.336 d; companion 0.80 M_sun (SB1 mass function (sin i = 1)) |
| TOI 7599.01 | 142757891 | APC | 2.3429 | gaia_orbit | Gaia NSS SB1 orbit P = 4.687 d (transit period = half Gaia's); companion 1.22 M_sun (SB1 mass function (sin i = 1)) |
| TOI 764.01 | 181159386 | APC | 5.6317 | gaia_orbit | Gaia NSS SB1 orbit P = 5.631 d; companion 0.26 M_sun (SB1 mass function (sin i = 1)) |
| TOI 7749.01 | 457542392 | PC | 20.6448 | gaia_orbit | Gaia NSS SB1 orbit P = 20.667 d; companion 0.27 M_sun (SB1 mass function (sin i = 1)) |
| TOI 7761.01 | 45636452 | PC | 26.6556 | gaia_orbit | Gaia NSS SB1 orbit P = 26.668 d; companion 0.60 M_sun (SB1 mass function (sin i = 1)) |
| TOI 7863.01 | 464286189 | PC | 18.2993 | gaia_orbit, odd_even | Gaia NSS SB1 orbit P = 18.294 d; companion 0.38 M_sun (SB1 mass function (sin i = 1)) ; odd and even transits differ by 5.2 sigma (2722 vs 1301 ppm) |
| TOI 865.01 | 44797824 | APC | 0.7456 | gaia_orbit, odd_even | Gaia NSS SB2 orbit P = 1.491 d (transit period = half Gaia's); companion 0.92 M_sun (SB2 mass ratio (K1/K2)) ; odd and even transits differ by 7.8 sigma (711 vs 1245 ppm) |
| TOI 901.01 | 214361331 | PC | 3.2311 | gaia_orbit | Gaia NSS SB1 orbit P = 3.231 d; companion 0.50 M_sun (SB1 mass function (sin i = 1)) |
| TOI 924.01 | 382068562 | PC | 12.1273 | gaia_orbit | Gaia NSS SB1 orbit P = 12.122 d; companion 0.29 M_sun (SB1 mass function (sin i = 1)) |
| TOI 948.01 | 146438872 | APC | 12.7009 | gaia_orbit | Gaia NSS SB1 orbit P = 12.711 d; companion 0.45 M_sun (SB1 mass function (sin i = 1)) |

**Also useful:** 5 unresolved candidates have a Gaia orbit on the transit period with a companion below 0.08 M_sun. Gaia may be seeing the transiting object itself: a massive planet, a brown dwarf or, near 75-80 M_Jup, one of the lowest-mass stars (`gaia_substellar_companions.csv`).

## Caveats

- A flag is evidence, not a verdict. TFOPWG makes dispositions from all the data, including follow-up observations that these checks don't use.
- **Masses.** The companion mass assumes the TIC mass for the primary star. The SB1 masses assume sin i = 1, which is appropriate when the companion eclipses. Masses from astrometry assume a dark companion, so they are lower limits.
- **Centroids** here are the light curves' flux-weighted centroids, not difference images. Crowded fields and saturated stars can shift them.
- **Secondary eclipses:** hot Jupiters show real, shallow ones (WASP-18 b: 3.6 % of the transit depth). That's why the check requires a secondary that is deep relative to the transit.
- The light-curve checks cover only a sample of the unresolved candidates, not all of them.
