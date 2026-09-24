"""Pixel-level localisation of a TESS periodic signal. Usage: python tess_pixel_test.py <gaia_dr3> <TIC> <sector> <freq_c_per_d> [120|20]
Target pixel file (SPOC); WCS offset and Gaussian PSF width fitted to the median image with Gaia DR3 sources (G < 19.5, proper motion
applied); sinusoid at the frequency fitted to every pixel; signed in-phase amplitude map (phase of the aperture sum) fitted with a
single PSF at each Gaia source within the stamp; chi2 per source (same number of parameters)."""
import sys, os, json, subprocess, numpy as np, requests
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from astroquery.mast import Observations
gaia, tic, sec, f0 = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]); cad = sys.argv[5] if len(sys.argv) > 5 else "120"
obs = Observations.query_criteria(obs_collection="TESS", target_name=tic, sequence_number=sec, dataproduct_type="timeseries", provenance_name="SPOC")
fn = [n for n in Observations.get_product_list(obs)["productFilename"] if n.endswith("fast-tp.fits" if cad == "20" else "s_tp.fits")][0]
path = os.path.join("..", "data", "cache", fn); os.makedirs(os.path.dirname(path), exist_ok=True)
if not os.path.exists(path):
    subprocess.run(["curl", "-sL", "-m", "1200", "-o", path, f"https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS/product/{fn}"], check=True)
th = fits.open(path); T = th[1].data; wcs = WCS(th[2].header); ap = th[2].data
q = (T["QUALITY"] == 0) & np.isfinite(T["TIME"]); t = T["TIME"][q]; FL = T["FLUX"][q]; good = np.all(np.isfinite(FL), axis=(1, 2)); t, FL = t[good], FL[good]
ny, nx = FL.shape[1:]; med = np.median(FL, axis=0); yy, xx = np.mgrid[0:ny, 0:nx]
ra0, de0 = wcs.all_pix2world([[nx / 2 - 0.5, ny / 2 - 0.5]], 0)[0]
TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
def tap(qs): return requests.get(TAP, params=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="json", QUERY=qs), timeout=120).json()["data"]
G = tap(f"SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE({ra0},{de0},{0.006*max(nx,ny)})) AND phot_g_mean_mag < 19.5")
G += [r for r in tap(f"SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE source_id={gaia}") if str(r[0]) not in {str(x[0]) for x in G}]
yr = 2016.0 + (2457000 + np.median(t) - 2457389.0) / 365.25; src = []
for sid, ra, de, pmra, pmde, gm in G:
    px, py = wcs.all_world2pix([[ra + (pmra or 0) * (yr - 2016) / 3.6e6 / np.cos(np.radians(de)), de + (pmde or 0) * (yr - 2016) / 3.6e6]], 0)[0]
    if -2 < px < nx + 1 and -2 < py < ny + 1:
        src.append((str(sid), gm, px, py))
psf = lambda x0, y0, s: np.exp(-0.5 * ((xx - x0) ** 2 + (yy - y0) ** 2) / s ** 2)
fl0 = np.array([10 ** (-0.4 * (g - 16.0)) for _, g, _, _ in src])
model = lambda p: p[3] * sum(f * psf(x + p[0], y + p[1], p[2]) for f, (_, _, x, y) in zip(fl0, src)) + p[4]
dx, dy, s = least_squares(lambda p: (model(p) - med).ravel(), [0, 0, 1.0, med.max(), 0], bounds=([-3, -3, 0.4, 0, -np.inf], [3, 3, 3, np.inf, np.inf])).x[:3]
X = np.vstack([np.ones_like(t), np.sin(2 * np.pi * f0 * t), np.cos(2 * np.pi * f0 * t)]).T; XtXi = np.linalg.inv(X.T @ X)
B = np.zeros((ny, nx, 3)); Sg = np.zeros((ny, nx))
for iy in range(ny):
    for ix in range(nx):
        v = FL[:, iy, ix] - np.median(FL[:, iy, ix]); b = XtXi @ (X.T @ v); B[iy, ix] = b; Sg[iy, ix] = np.sqrt(np.sum((v - X @ b) ** 2) / (len(v) - 3))
apm = ap & 2 > 0; phi = np.arctan2(B[..., 2][apm].sum(), B[..., 1][apm].sum())
a = B[..., 1] * np.cos(phi) + B[..., 2] * np.sin(phi); e = Sg * np.sqrt(XtXi[1, 1])
tg = [z for z in src if z[0] == gaia][0]; m = np.hypot(xx - (tg[2] + dx), yy - (tg[3] + dy)) < 4; out = []
for sid, gm, x, y in src:
    P = psf(x + dx, y + dy, s)[m]
    if P.sum() < 0.3:
        continue
    A = np.sum(a[m] * P / e[m] ** 2) / np.sum(P ** 2 / e[m] ** 2); eA = 1 / np.sqrt(np.sum(P ** 2 / e[m] ** 2))
    out.append(dict(chi2=round(float(np.sum(((a[m] - A * P) / e[m]) ** 2)), 1), gaia=sid, G=round(float(gm), 2), amp=round(float(A), 4), e_amp=round(float(eA), 4),
                    sep_arcsec=round(float(np.hypot(x - tg[2], y - tg[3]) * 21), 1)))
out.sort(key=lambda r: r["chi2"])
res = dict(gaia=gaia, tic=tic, sector=sec, freq=f0, cadence=cad, no_signal_chi2=round(float(np.sum((a[m] / e[m]) ** 2)), 1), n_pix=int(m.sum()), sources=out)
json.dump(res, open(f"pixel_{gaia}_S{sec}.json", "w"), indent=1)
for r in out[:6]:
    print(("TARGET " if r["gaia"] == gaia else "       ") + json.dumps(r))
