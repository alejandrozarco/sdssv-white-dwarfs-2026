# SDSS-V DR20 white dwarfs: measurements

Measurements of white dwarfs from public SDSS-V DR20 spectra (Astra 0.8.1), TESS, HST/COS, GALEX, ATLAS, ZTF and Gaia DR3 epoch photometry. Each table gives the measured quantities with existing SIMBAD, MWDD (snapshot 2026-08-05) and SDSS-V SnowWhite classifications. The scripts in `scripts/` download the public data and recompute every table and figure. Data were retrieved 2026-09-23 to 2026-09-25. Methods are in [METHODS.md](METHODS.md).

| topic | table(s) | objects |
|---|---|---|
| [Zeeman splitting](#zeeman-splitting) | `magnetic_zeeman.csv` | 30 |
| [Carbon lines](#carbon-lines) | `carbon_white_dwarfs.csv`, `carbon_*_optical_CII_features.csv`, `carbon_5208047381438507520_cos_features.csv` | 4 |
| [TESS amplitude spectra](#tess-amplitude-spectra) | `zz_ceti_objects.csv`, `zz_ceti_tess_sectors.csv` | 5 |
| [Eclipse](#eclipse) | `eclipsing_4731701084150029824.csv` | 1 |
| [Balmer emission](#balmer-emission) | `balmer_emission.csv` | 2 |
| [Photometric periods](#photometric-periods) | `periodic_6021870154194477312.csv`, `periodic_white_dwarfs.csv` | 8 |

## Zeeman splitting
H-alpha and H-beta are fitted with three Gaussian components. B_split is the linear-Zeeman field implied by the separation of the two outer components. Red lines in the figure mark the fitted centres.

<img src="figures/magnetic_zeeman_halpha_hbeta.png" width="600">

## Carbon lines
Line cross-correlation (C II, C I, H I, He I) for four SDSS-V white dwarfs, Gaussian fits to optical C II lines, HST/COS G130M equivalent widths, and GALEX FUV−NUV colours.

The first figure shows the SDSS-V spectra of Gaia DR3 5208047381438507520 and Gaia DR3 6886051830805052288 (middle two). Above them is an SDSS-V DA of similar colour, and below them is the hot DQ SDSS J234843.30−094245.3 (Dufour et al. 2008). Orange lines mark C II positions; blue dotted lines mark Balmer positions. A three-spectrum version with only 5208047381438507520 is `figures/hot_dq_comparison_5208047381438507520.png`.

<img src="figures/hot_dq_comparison_sdssv.png" width="700">

<img src="figures/carbon_optical_spectra.png" width="700">
<img src="figures/carbon_5208047381438507520_cos.png" width="700">
<img src="figures/galex_fuv_nuv.png" width="450">

## TESS amplitude spectra
The highest TESS peak per light curve for five SDSS-V white dwarfs, with pixel-level fits to locate the signal.

<img src="figures/zz_ceti_tess_amplitude_spectra.png" width="700">

## Eclipse
Eclipse ephemeris of Gaia DR3 4731701084150029824 from ATLAS forced photometry.

<img src="figures/eclipsing_4731701084150029824_atlas_phase.png" width="600">

## Balmer emission
H-alpha and H-beta emission equivalent widths and double-peak separations.

<img src="figures/balmer_emission.png" width="600">

## Photometric periods
Gaia DR3 6021870154194477312 is shown in the first figure:
- top: its SDSS-V spectrum;
- middle: the ATLAS periodogram;
- bottom: ATLAS, Gaia and TESS light curves folded on one frequency.

`periodic_6021870154194477312.csv` also gives fits for the three Gaia sources within 13″.

<img src="figures/periodic_6021870154194477312.png" width="650">

`periodic_white_dwarfs.csv` covers seven white dwarfs selected by their Gaia DR3 GLS frequency. Each frequency was recovered in ATLAS or ZTF, with the Gaia and TESS amplitudes and times of maximum at that frequency.

<img src="figures/periodic_white_dwarfs.png" width="650">

## Reproduction
```
pip install -r requirements.txt
cd scripts
python zeeman_split.py ../data/zeeman_input.csv
python carbon_lines.py 95077848 5208047381438507520
python carbon_lines.py 110600288 6466745168812781568
python carbon_lines.py 102600838 5836110898905253760
python carbon_lines.py 114554634 6886051830805052288
python cos_lines.py
python galex_colours.py 5208047381438507520
python galex_colours.py 6886051830805052288
python tess_periodogram.py 2055170284 102 120
python tess_pixel_test.py 6492083311194727168 2055170284 102 82.34 120
python j0353_eclipse.py
python cv_balmer.py 65701864
python periodic_6021870154194477312.py
python tess_periodogram.py 1251484163 65 120 0.5 50
python periodic_white_dwarfs.py
python figures.py
```
Arguments for the other TESS light curves and pixel tests are in the table columns (TIC, sector, cadence, frequency).

## Data sources
SDSS-V DR20 and SDSS DR17; Gaia DR3 (ESA/Gaia/DPAC), including epoch photometry (VizieR I/355/epphot) and the Gaia Synthetic Photometry Catalogue (VizieR J/A+A/674/A33); TESS SPOC and HST/COS program 17420 (MAST); GALEX GUVcat AIS (Bianchi et al. 2017); ATLAS forced photometry (Tonry et al. 2018; Shingles et al. 2021); ZTF public data releases (IRSA); NIST Atomic Spectra Database; Montreal White Dwarf Database; SIMBAD.
