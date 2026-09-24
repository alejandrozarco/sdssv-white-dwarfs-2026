"""Lomb-Scargle amplitude spectrum of a TESS SPOC light curve. Usage: python tess_periodogram.py <TIC> <sector> [120|20] [fmin] [fmax]
PDCSAP flux, QUALITY == 0, normalised, 5-sigma (MAD) clip; frequency step 0.001 c/d; S/N = amplitude / mean amplitude;
Baluev false-alarm probability of the highest peak."""
import sys, os, subprocess, numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle
from astroquery.mast import Observations
tic, sec = sys.argv[1], int(sys.argv[2]); cad = sys.argv[3] if len(sys.argv) > 3 else "120"
fmin = float(sys.argv[4]) if len(sys.argv) > 4 else 20.0; fmax = float(sys.argv[5]) if len(sys.argv) > 5 else (359.0 if cad == "120" else 2159.0)
obs = Observations.query_criteria(obs_collection="TESS", target_name=tic, sequence_number=sec, dataproduct_type="timeseries", provenance_name="SPOC")
names = [n for n in Observations.get_product_list(obs)["productFilename"] if n.endswith("fast-lc.fits" if cad == "20" else "s_lc.fits")]
fn = names[0]; path = os.path.join("..", "data", "cache", fn); os.makedirs(os.path.dirname(path), exist_ok=True)
if not os.path.exists(path):
    subprocess.run(["curl", "-sL", "-m", "600", "-o", path, f"https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS/product/{fn}"], check=True)
h = fits.open(path); d = h[1].data
m = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]) & np.isfinite(d["TIME"]); t = d["TIME"][m]; y = d["PDCSAP_FLUX"][m] / np.median(d["PDCSAP_FLUX"][m]) - 1
mad = 1.4826 * np.median(np.abs(y - np.median(y))); ok = np.abs(y - np.median(y)) < 5 * mad; t, y = t[ok], y[ok]
fr = np.arange(fmin, fmax, 0.001); amp = np.sqrt(4 * LombScargle(t, y, normalization="psd").power(fr, method="fast") / len(t)); mean = amp.mean()
peaks = []
for k in np.argsort(amp)[::-1]:
    if all(abs(fr[k] - q[0]) > 0.05 for q in peaks):
        peaks.append((fr[k], amp[k]))
    if len(peaks) == 5:
        break
ls = LombScargle(t, y); fap = ls.false_alarm_probability(ls.power(np.array([peaks[0][0]]))[0], minimum_frequency=fmin, maximum_frequency=fmax, method="baluev")
print(f"TIC {tic} S{sec} {cad}-s: n {len(t)}, CROWDSAP {h[1].header.get('CROWDSAP')}, mean amplitude {mean*1e3:.2f} ppt, FAP(top) {fap:.2g}")
for f, a in peaks:
    print(f"  {f:10.4f} c/d  P = {86400/f:8.1f} s  amplitude {a*1e3:6.2f} ppt  S/N {a/mean:5.1f}")
