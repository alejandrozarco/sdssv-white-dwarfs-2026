"""Zeeman splitting of H-alpha and H-beta in SDSS-V spectra.
Model per line: linear continuum x (1 - broad concentric Gaussian - three Gaussian absorption components with free centres,
widths and depths). B_split = half-separation of the outer components / (4.67e-13 lambda0^2 [A/G]), i.e. 20.13 A/MG (H-alpha)
and 11.04 A/MG (H-beta).
Usage: python zeeman_split.py tables/magnetic_zeeman.csv  (columns gaia_dr3, sdss_id, spectrum = coadd or visit MJD)"""
import sys, numpy as np, pandas as pd
from scipy.optimize import least_squares
from sdssv import star_spectrum, visits
L0 = {"Ha": 6564.61, "Hb": 4862.68}; K = {k: 4.67e-13 * v ** 2 * 1e6 for k, v in L0.items()}
WIN = {"Ha": (6300, 6830), "Hb": (4660, 5060)}


def model(p, x, l0):
    c0, c1, dw, sw, d0, s0, dm, sm, dp, sp, lm, lp, sh = p
    g = lambda mu, s: np.exp(-0.5 * ((x - mu) / s) ** 2)
    return (c0 + c1 * (x - l0) / 100) * (1 - dw * g(l0 + sh, sw) - d0 * g(l0 + sh, s0) - dm * g(l0 + sh - lm, sm) - dp * g(l0 + sh + lp, sp))


def fit(x, y, iv, l0, dl0):
    m = (iv > 0) & np.isfinite(y); x, y, w = x[m], y[m], np.sqrt(iv[m])
    med = np.median(y[np.abs(x - l0) > 0.8 * (x.max() - l0)]); best = None
    for dl in dl0 * np.array([0.6, 0.8, 1.0, 1.25, 1.6]):
        for sw in (40, 80):
            p0 = [med, 0, 0.15, sw, 0.2, 6, 0.15, 8, 0.15, 8, dl, dl, 0]
            lo = [0, -np.inf, 0, 15, 0, 1.5, 0, 1.5, 0, 1.5, 5, 5, -15]; hi = [np.inf, np.inf, 1, 300, 1, 40, 1, 60, 1, 60, 400, 400, 15]
            r = least_squares(lambda p: (model(p, x, l0) - y) * w, p0, bounds=(lo, hi), max_nfev=4000)
            if best is None or r.cost < best.cost:
                best = r
    chi2r = 2 * best.cost / (len(x) - len(best.x)); C = np.linalg.inv(best.jac.T @ best.jac) * max(chi2r, 1)
    return best.x, C


def measure(sdss_id, spectrum="coadd"):
    if spectrum == "coadd":
        lam, f, iv, snr = star_spectrum(sdss_id)
    else:
        v = [x for x in visits(sdss_id) if str(x["mjd"]) == str(spectrum)][0]; lam, f, iv, snr = v["wave"], v["flux"], v["ivar"], v["snr"]
    out = dict(snr=round(snr, 1))
    for k in ("Ha", "Hb"):
        m = (lam > WIN[k][0]) & (lam < WIN[k][1]); p, C = fit(lam[m], f[m], iv[m], L0[k], 6.0 * K[k])
        out[f"B_split_{k}_MG"] = round((p[10] + p[11]) / 2 / K[k], 2)
        out[f"e_B_split_{k}_MG"] = round(0.5 * np.sqrt(C[10, 10] + C[11, 11] + 2 * C[10, 11]) / K[k], 2)
    return out


if __name__ == "__main__":
    t = pd.read_csv(sys.argv[1], dtype={"gaia_dr3": str, "sdss_id": str})
    rows = []
    for _, r in t.iterrows():
        res = measure(r["sdss_id"], str(r["spectrum"])); rows.append(dict(gaia_dr3=r["gaia_dr3"], **res))
        print(r["gaia_dr3"], res, flush=True)
    pd.DataFrame(rows).to_csv("zeeman_split_measured.csv", index=False)
