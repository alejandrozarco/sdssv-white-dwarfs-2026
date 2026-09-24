"""H-alpha and H-beta emission of SDSS-V spectra. Usage: python cv_balmer.py <sdss_id>
Coadd of all visits (sdssv.coadd). Equivalent width (positive = emission) within +-35 A (H-alpha) / +-30 A (H-beta) against a
linear continuum from the sidebands 6480-6520 and 6610-6650 A (H-alpha), 4780-4815 and 4950-4990 A (H-beta). Peak separation from
a fit of two Gaussians of equal width on a linear continuum."""
import sys, numpy as np
from scipy.optimize import curve_fit
from sdssv import visits, coadd, C
w, f, iv = coadd(visits(sys.argv[1]))
for name, l0, half, sb in (("H-alpha", 6564.61, 35, [(6480, 6520), (6610, 6650)]), ("H-beta", 4862.68, 30, [(4780, 4815), (4950, 4990)])):
    s = np.zeros_like(w, bool)
    for a, b in sb:
        s |= (w > a) & (w < b)
    s &= np.isfinite(f) & (iv > 0); cp = np.polyfit(w[s], f[s], 1)
    win = (np.abs(w - l0) < half) & np.isfinite(f) & (iv > 0); cont = np.polyval(cp, w[win])
    dw = np.gradient(w[win]); ew = np.sum((f[win] / cont - 1) * dw); eew = np.sqrt(np.sum((dw / np.sqrt(iv[win]) / cont) ** 2))
    x = w[win] - l0; y = f[win] / cont - 1; e = 1 / np.sqrt(iv[win]) / cont
    g2 = lambda x, a1, a2, m, sep, sg: a1 * np.exp(-0.5 * ((x - m + sep / 2) / sg) ** 2) + a2 * np.exp(-0.5 * ((x - m - sep / 2) / sg) ** 2)
    p, cv = curve_fit(g2, x, y, p0=[y.max(), y.max(), 0, 20, 6], sigma=e, absolute_sigma=True, bounds=([0, 0, -20, 2, 1], [np.inf, np.inf, 20, 60, 30]), maxfev=20000)
    print(f"{name}: EW {ew:.1f} +- {eew:.1f} A; peak separation {p[3] / l0 * C:.0f} +- {np.sqrt(cv[3, 3]) / l0 * C:.0f} km/s; peak width (sigma) {p[4] / l0 * C:.0f} km/s")
