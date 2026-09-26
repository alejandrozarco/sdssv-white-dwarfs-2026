"""Aperture photometry of TESS full-frame images (TESScut, 11 x 11 pixels) at a position, and the signal at a given frequency.
Usage: python tess_ffi_photometry.py <gaia_dr3> <ra> <dec> <frequency_cd> <Tmag> <sector> [<sector> ...]
Per sector: QUALITY == 0; aperture = 3 x 3 pixels centred on the target pixel; background = per-cadence median of the pixels outside a
5 x 5 box, times 9; 1-day running median subtracted; 5-sigma clip. Lomb-Scargle over 5-30 c/d (highest peak, Baluev false-alarm
probability) and the power rank of the given frequency (+-0.01 c/d). Amplitude of a sinusoid at the given frequency as a fraction of
the aperture flux, and scaled to the target's expected flux (15,000 e/s x 10^(-0.4 (Tmag - 10))), which assumes all other flux in the
aperture is constant. Writes ../tables/tess_ffi_<gaia_dr3>.csv."""
import sys, os, numpy as np, pandas as pd
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
from astropy.timeseries import LombScargle
from astropy.wcs import WCS
from scipy.ndimage import median_filter

gid, ra, dec, f0, tmag = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]); sectors = [int(s) for s in sys.argv[6:]]
c = SkyCoord(ra, dec, unit="deg"); expected = 15000 * 10 ** (-0.4 * (tmag - 10)); rows = []
for sec in sectors:
    h = Tesscut.get_cutouts(coordinates=c, size=11, sector=sec)[0]
    d = h[1].data; x, y = WCS(h[2].header).world_to_pixel(c); xi, yi = int(round(float(x))), int(round(float(y)))
    q = d["QUALITY"] == 0; fl = d["FLUX"][q]; t = d["TIME"][q] + 2457000.0
    yy, xx = np.mgrid[0:fl.shape[1], 0:fl.shape[2]]; bg = np.nanmedian(fl[:, (np.abs(xx - xi) > 2) | (np.abs(yy - yi) > 2)], axis=1)
    ap = (np.abs(xx - xi) <= 1) & (np.abs(yy - yi) <= 1); lc = np.nansum(fl[:, ap], axis=1) - bg * ap.sum(); ok = np.isfinite(lc) & np.isfinite(t)
    t, lc = t[ok], lc[ok]; o = np.argsort(t); t, lc = t[o], lc[o]; med = np.median(lc); r = lc / med - 1
    r = r - median_filter(r, size=int(1 / np.median(np.diff(t))) | 1, mode="nearest"); s = 1.4826 * np.median(np.abs(r)); m = np.abs(r) < 5 * s; t, r = t[m], r[m]
    fr = np.arange(5, 30, 0.001); ls = LombScargle(t, r); P = ls.power(fr); k = int(np.argmax(P)); near = np.abs(fr - f0) < 0.01
    X = np.vstack([np.ones_like(t), np.cos(2 * np.pi * f0 * t), np.sin(2 * np.pi * f0 * t)]).T; p = np.linalg.lstsq(X, r, rcond=None)[0]
    amp = float(np.hypot(p[1], p[2])); err = float(np.std(r - X @ p) * np.sqrt(2 / len(t)))
    rows.append(dict(gaia_dr3=gid, sector=sec, n=len(t), cadence_min=round(float(np.median(np.diff(t))) * 1440, 1), aperture_median_e_s=round(float(med), 1),
                     peak_cd=round(float(fr[k]), 4), peak_fap=float(f"{ls.false_alarm_probability(P[k], minimum_frequency=5, maximum_frequency=30):.2g}"),
                     power_rank_at_frequency=int(np.sum(P > P[near].max())), frequency_cd=f0, amplitude_aperture_frac=round(amp, 4), e_amplitude_aperture_frac=round(err, 4),
                     expected_target_e_s=round(expected, 1), amplitude_target_frac=round(amp * med / expected, 3)))
out = pd.DataFrame(rows); out.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tables", f"tess_ffi_{gid}.csv"), index=False); print(out.to_string())
