"""HST/COS G130M spectrum lfac0z010 (program 17420) of Gaia DR3 5208047381438507520: download from MAST, equivalent widths
of C/Si/N features against local median continua (several sideband choices), and core positions of C II 1334.53/1335.71."""
import numpy as np, pandas as pd
from astropy.io import fits
from astroquery.mast import Observations
t = Observations.query_criteria(obs_id="lfac0z010")
p = Observations.get_product_list(t); p = p[[f == "lfac0z010_x1dsum.fits" for f in p["productFilename"]]]
path = Observations.download_products(p, download_dir="../data/cache")["Local Path"][0]
d = fits.open(path)[1].data
W = np.concatenate([r["WAVELENGTH"] for r in d]); F = np.concatenate([r["FLUX"] for r in d]); E = np.concatenate([r["ERROR"] for r in d])
Q = np.concatenate([r["DQ_WGT"] for r in d]); o = np.argsort(W); W, F, E, Q = W[o], F[o], E[o], Q[o]; ok = (Q > 0) & (E > 0); W, F, E = W[ok], F[ok], E[ok]
dw = np.median(np.diff(W)); C = 299792.458


def ew(center, half, side):
    c = np.median(np.concatenate([F[(W > center + a) & (W < center + b)] for a, b in side])); s = np.abs(W - center) < half
    return np.sum(1 - F[s] / c) * dw, np.sqrt(np.sum((E[s] / c) ** 2)) * dw


rows = []
for name, cen, half, sides in [("C II 1335", 1335.1, 2.0, [((-6, -3), (3, 6)), ((-9, -4), (4, 9)), ((-12, -6), (6, 12)), ((-15, -8), (8, 15)), ((-20, -12), (12, 20))]),
                               ("C III 1175", 1175.6, 2.0, [((-8, -4), (4, 8))]), ("C III 1247", 1247.38, 0.8, [((-4, -2), (2, 4))]),
                               ("C II 1324", 1323.95, 0.8, [((-4, -2), (2, 4))]), ("Si II 1260", 1260.42, 0.8, [((-4, -2), (2, 4))]),
                               ("Si IV 1394", 1393.76, 0.8, [((-4, -2), (2, 4))]), ("N V 1239", 1238.82, 0.8, [((-4, -2), (2, 4)), ((-10, -6), (6, 10))])]:
    for sd in sides:
        v, e = ew(cen, half, sd); rows.append(dict(feature=name, center_A=cen, half_width_A=half, sidebands_A=str(sd), EW_A=round(v, 3), e_EW_A=round(e, 3)))
        print(f"{name:11s} sidebands {sd}: EW {v:+.3f} +- {e:.3f} A")
for lab in (1334.532, 1335.708):
    s = np.abs(W - lab) < 0.6; k = np.argmin(np.convolve(F[s], np.ones(3) / 3, mode="same")[1:-1]) + 1
    print(f"core of C II {lab}: {W[s][k]:.3f} A, {(W[s][k] / lab - 1) * C:+.0f} km/s")
    rows.append(dict(feature=f"C II {lab} core", center_A=round(W[s][k], 3), velocity_kms=round((W[s][k] / lab - 1) * C)))
pd.DataFrame(rows).to_csv("cos_features_5208047381438507520.csv", index=False)
