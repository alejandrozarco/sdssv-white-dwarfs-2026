# Short-period white dwarfs with irradiated companions

`tables/irradiated_companions.csv` lists four white dwarfs with photometric periods of 81-130 min. In each:
- the amplitude is larger in red than in blue light;
- the WISE photometry (and VISTA photometry, where available) exceeds the white-dwarf prediction.

`tables/irradiated_companions_periods.csv` gives the highest periodogram peak and the semi-amplitude at the adopted period for each data set. `tables/desi_halpha_6914922055508553984.csv` gives the H-alpha emission fit of WDJ205249.27−032419.53. The method is in [METHODS.md](../METHODS.md#methods).

| star | Gaia DR3 | G | distance (pc) | GF21 Teff, mass (H) | P (min) | semi-amplitude | W1, W2 / white-dwarf model | companion M_W1 |
|---|---|---|---|---|---|---|---|---|
| WDJ205249.27−032419.53 | 6914922055508553984 | 17.47 | 334 | 16.2 kK, 0.30 Msun | 97.703 | ZTF g 2.9%, r 7.0%; TESS 10.1-10.7% | 3.7, 4.4 | 9.25 |
| WDJ212738.67+593755.72 | 2191618770599895296 | 16.87 | 241 | 13.8 kK, 0.26 Msun | 130.182 | ZTF g 1.5%, r 4.2%; TESS 6.1-8.7% | 3.6, 3.5 | 9.28 |
| WDJ070106.16−534811.37 | 5503429908930455808 | 18.43 | 602 | 16.2 kK, 0.26 Msun | 81.493 | Gaia G 8.8%, BP 6.6 ± 2.5%, RP 14 ± 5% | 2.2, 2.1 | 9.82 |
| WDJ040444.35−395043.1 | 4844023064578952320 | 17.66 | 632 | 20.6 kK, 0.24 Msun | 117.319 | Gaia G 9.8%, BP 3.9 ± 1.9%, RP 22 ± 3% | 2.0, 1.3 | 9.33 |

- **Temperatures and masses:** Gentile Fusillo et al. (2021) H-atmosphere fits to Gaia photometry.
- **Distance:** 1/parallax.
- **Infrared ratios:** observed flux divided by the pure-H model flux at that Teff and log g, scaled to Gaia G.
- **Companion M_W1:** the W1 excess converted to an absolute magnitude.
- **Limits of the photometry:**
  - W1 of WDJ070106.16−534811.37 is 18.07; its W2 and the W2 of WDJ040444.35−395043.1 have errors of 0.28-0.31 mag, close to the CatWISE limit.
  - TESS full-frame-image amplitudes are fractions of the total aperture flux, not corrected for other stars; only their periods are used.

## WDJ205249.27−032419.53 (Gaia DR3 6914922055508553984)
- **Period:** P = 97.703 min in ZTF, Gaia DR3 and TESS (sectors 55 and 81, 120 s). The nearest alias has Δχ² = 4918.
- **Existing classifications:**
  - VSX and Chen et al. (2020) list type DSCT with P = 0.0678498 d.
  - DESI DR1: DAe (Amorim et al. 2026) and WD+MS (Swan et al. 2026).
  - Kilic et al. (2026): DA, 18.4 kK, log g 7.26.
- **DESI DR1 spectrum** (one 494-s exposure, tile 20836, petal 8, MJD 59358.46):
  - Balmer absorption and a narrow H-alpha emission line at +189 ± 8 km/s (σ = 1.9 Å), at photometric phase 0.82 (phase 0 = maximum light).
  - The median residual of 28 galaxy and QSO spectra on the same tile and petal peaks at 6561.6 Å (+18%). This is the flux-calibration feature described by Swan et al. (2026).
  - A fit with that residual as a free multiplicative term needs the emission line: Δχ² = 211 for three parameters.

<img src="../figures/irradiated_companions/6914922055508553984.png" width="800">

<img src="../figures/irradiated_companions/6914922055508553984_desi_halpha.png" width="800">

## WDJ212738.67+593755.72 (Gaia DR3 2191618770599895296)
- **Period:** P = 130.182 min in ZTF, Gaia DR3 and TESS (sectors 76, 77, 83 and 84, 120 s). The nearest alias has Δχ² = 4206.
- **Existing classifications:**
  - VSX and Chen et al. (2020) list type DSCT with P = 0.0904050 d.
  - Jestin et al. (2026) list it as periodic.
- **Spectra:** none found in SDSS, DESI DR1, SDSS-V DR20 or LAMOST DR11.

<img src="../figures/irradiated_companions/2191618770599895296.png" width="800">

## WDJ070106.16−534811.37 (Gaia DR3 5503429908930455808)
- **Period:** P = 81.493 min in Gaia DR3 and TESS full-frame images (sectors 88, 89, 93, 96 and 98). The nearest alias has Δχ² = 84.
- **Existing classifications:**
  - VSX lists Gaia DR3 type VAR with P = 0.0565929 d.
  - The period is also in Ranaivomanana et al. (2025).
  - Gaia XP class DB (Vincent et al. 2024).
- **Spectra:** none found.

<img src="../figures/irradiated_companions/5503429908930455808.png" width="800">

## WDJ040444.35−395043.1 (Gaia DR3 4844023064578952320)
- **Period:** P = 117.319 min in Gaia DR3 and TESS full-frame images (sectors 106 and 107). The nearest alias has Δχ² = 1366.
- **Existing classifications:** VSX lists type WD with the Gaia period (0.0814713 d).
- **Spectrum:** SDSS-V DR20 has one visit (S/N 14), classified DA by SnowWhite (sdss_id 93071386).

<img src="../figures/irradiated_companions/4844023064578952320.png" width="800">
