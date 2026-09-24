"""Eclipse ephemeris of Gaia DR3 4731701084150029824 from ATLAS forced photometry (data/atlas_forced_photometry_4731701084150029824.txt;
difference fluxes uJy, bands o and c). Cuts: duJy > 0, err == 0, chi/N < 10, duJy < 3 x median; per-season median subtracted;
5-sigma clip; times converted to BJD_TDB. Model per band: offset + trapezoid (common centre, total width, flat fraction) with free
depth, offset and depth solved linearly on a grid; period from the chi2 minimum, error from delta chi2 = 1 after scaling chi2_r to 1."""
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
geo = EarthLocation.from_geocentric(0, 0, 0, unit="m"); c0 = SkyCoord(58.30184 * u.deg, -55.04399 * u.deg)
L = [l for l in open("../data/atlas_forced_photometry_4731701084150029824.txt").read().splitlines() if l.strip()]
hdr = L[0].lstrip("#").split(); R = [dict(zip(hdr, l.split())) for l in L[1:]]
ok = [x for x in R if float(x["duJy"]) > 0 and float(x["err"]) == 0 and float(x["chi/N"]) < 10]
D = {}
for b in ("o", "c"):
    s = [x for x in ok if x["F"] == b]; med = np.median([float(x["duJy"]) for x in s]); s = [x for x in s if float(x["duJy"]) < 3 * med]
    mjd = np.array([float(x["MJD"]) for x in s]); f = np.array([float(x["uJy"]) for x in s]); e = np.array([float(x["duJy"]) for x in s])
    season = np.floor((mjd - 57000) / 365.25 + 0.3).astype(int)
    for sv in np.unique(season):
        f[season == sv] -= np.median(f[season == sv])
    clip = np.abs(f) < 5 * 1.4826 * np.median(np.abs(f)) + 3 * np.median(e)
    t = Time(mjd[clip], format="mjd", scale="utc", location=geo); D[b] = dict(t=(t.tdb + t.light_travel_time(c0)).jd, f=f[clip], e=e[clip])
Tref = 2460670.38343


def trap(x, T, fr):
    a = np.abs(x); tf = fr * T; m = np.zeros_like(a); m[a <= tf / 2] = 1
    r = (a > tf / 2) & (a < T / 2); m[r] = (T / 2 - a[r]) / max(T / 2 - tf / 2, 1e-12); return m


CEN = np.arange(-0.02, 0.02001, 0.0005); TW = np.arange(0.016, 0.0801, 0.002); FR = (0.0, 0.25, 0.5, 0.75, 0.95)


def scan(P):
    best = (np.inf, None); ph = {b: ((d["t"] - Tref) / P + 0.5) % 1 - 0.5 for b, d in D.items()}
    for T in TW:
        for fr in FR:
            for c in CEN:
                chi = 0; par = []
                for b, d in D.items():
                    m = trap(ph[b] - c, T, fr); w = 1 / d["e"] ** 2
                    S, Sm, Smm, Sy, Sym = np.sum(w), np.sum(w * m), np.sum(w * m * m), np.sum(w * d["f"]), np.sum(w * d["f"] * m)
                    det = S * Smm - Sm ** 2; off = (Smm * Sy - Sm * Sym) / det; dep = -(S * Sym - Sm * Sy) / det
                    chi += np.sum(w * (d["f"] - off + dep * m) ** 2); par += [off, dep]
                if chi < best[0]:
                    best = (chi, (c, T, fr, *par))
    return best


Ps = 0.1478696 + np.linspace(-1.5e-6, 1.5e-6, 31); chis = np.array([scan(P)[0] for P in Ps]); j = int(np.argmin(chis))
s2 = chis[j] / (sum(len(d["t"]) for d in D.values()) - 7); sl = slice(max(j - 6, 0), j + 7)
a, b, _ = np.polyfit(Ps[sl] - Ps[j], chis[sl], 2); P = Ps[j] - b / (2 * a); sigP = np.sqrt(s2 / a)
chi, (c, T, fr, oo, do, oc, dc) = scan(P)
print(f"P = {P:.8f} +- {sigP:.8f} d; chi2_r {s2:.2f}; T0 = BJD_TDB {Tref + c * P:.5f}; total duration {T * P * 1440:.1f} min ({T * 100:.1f}% of P), "
      f"flat fraction {fr}; depth o {do:.1f} uJy, c {dc:.1f} uJy; n(o) {len(D['o']['t'])}, n(c) {len(D['c']['t'])}")
