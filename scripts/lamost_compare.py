"""LAMOST LRS spectrum of Gaia DR3 883885440381808000 (obsid 8002245, 2011-11-24) against its SDSS-V spectrum (sdss_id 57623143).
Downloads the LAMOST DR10 v2.0 FITS file, applies the carbon_screen.py matched filter and plots both spectra around C I and C II
lines. Usage: python lamost_compare.py"""
import os, requests, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import percentile_filter, gaussian_filter1d
from sdssv import visits, C, CACHE
import carbon_screen as cs
path = os.path.join(CACHE, "lamost_dr10_8002245.fits.gz")
if not os.path.exists(path):
    r = requests.get("https://www.lamost.org/dr10/v2.0/spectrum/fits/8002245", timeout=300); r.raise_for_status(); open(path, "wb").write(r.content)
d = fits.open(path)[1].data[0]
w = np.array(d["WAVELENGTH"], float); f = np.array(d["FLUX"], float); iv = np.array(d["IVAR"], float); ok = (iv > 0) & np.isfinite(f) & (d["ANDMASK"] == 0)
fl = np.interp(cs.GRID, w[ok], f[ok], left=np.nan, right=np.nan); ivl = np.interp(cs.GRID, w[ok], iv[ok], left=0, right=0)
p = cs.prep(fl, ivl)
for sp in ("C", "C I", "C II", "He I"):
    c = cs.contrast(*p, sp); k = int(np.argmax(np.where(cs.INW, c, -99))); kk = int(np.argmax(c))
    print(f"LAMOST {sp}: contrast {c[k]:.2f} at {cs.VELS[k]:+.0f} km/s (-200..+400 window); maximum over -1500..+1500: {c[kk]:.2f} at {cs.VELS[kk]:+.0f} km/s")
v = visits("57623143")[0]; okv = (v["ivar"] > 0) & np.isfinite(v["flux"]); fs = np.interp(cs.GRID, v["wave"][okv], v["flux"][okv], left=np.nan, right=np.nan)


def norm(f):
    m = np.isfinite(f); ff = np.interp(np.arange(len(cs.GRID)), np.where(m)[0], f[m])
    n = int(round(np.log10(1 + 150 / 5000) / 6e-5)); return gaussian_filter1d(ff / gaussian_filter1d(percentile_filter(ff, 85, size=n), n / 3), 1.2)


ns, nl = norm(fs), norm(fl); vref = 110
fig, axs = plt.subplots(2, 4, figsize=(15, 7))
for ax, (a, b) in zip(axs.flat, [(3880, 4130), (4200, 4420), (4700, 5100), (5100, 5420), (5850, 6100), (6480, 6680), (7050, 7300), (8250, 8400)]):
    m = (cs.GRID > a) & (cs.GRID < b)
    ax.plot(cs.GRID[m], ns[m], "k-", lw=0.8, label="SDSS-V MJD 60669"); ax.plot(cs.GRID[m], nl[m] + 0.3, color="tab:purple", lw=0.8, label="LAMOST 2011-11-24 (+0.3)")
    for sp, col in (("C I", "tab:red"), ("C II", "tab:blue")):
        for lam in cs.LINES[sp]:
            if a < lam * (1 + vref / C) < b:
                ax.axvline(lam * (1 + vref / C), color=col, lw=0.7)
    for lam in cs.BALMER:
        if a < lam < b:
            ax.axvline(lam, color="tab:green", ls=":", lw=0.8)
    ax.set_xlim(a, b); ax.set_xlabel("vacuum wavelength (A)", fontsize=7)
axs.flat[0].legend(fontsize=7)
fig.suptitle(f"Gaia DR3 883885440381808000: C I (red) and C II (blue) at {vref:+d} km/s; Balmer rest wavelengths (green dotted)", fontsize=9)
fig.tight_layout(); fig.savefig("../figures/carbon_883885440381808000_lamost.png", dpi=85)
