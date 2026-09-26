"""Gaia DR3 3107374277060584064 (WDJ064438.09-004550.51): periodic modulation and moving emission lines.

Data:
- CoRoT faint-star light curves of CoRoT 102743730 (runs IRa01, LRa01, LRa06; CDS B/corot, downloaded to ../data/cache/). CoRoT 102743730
  is Gaia DR3 3107374272762041856 (G = 16.17), 4.1 arcsec away; the CoRoT photometric mask also contains this star. BAR extension,
  STATUS == 0, DATEBARTT + 2400000.0 = BJD(TT). Per run: flux divided by its median, 3-day running median subtracted, 5-sigma clip.
- ZTF DR light curves of this star and of the neighbour (../data/ztf_<gaia_dr3>.csv; IRSA light-curve service, 2 arcsec):
  catflags == 0, magerr < 0.25; fractional flux about the median of each filter. HJD is used as BJD (difference < 10 s).
- Gaia DR3 epoch photometry (../data/gaia_dr3_epoch_photometry_3107374277060584064.csv; VizieR I/355/epphot): G transits with
  GrVFlag == 0; TimeG + 2455197.5 = BJD.
- SDSS-V DR20 visit spectra (sdss_id 74709777; sdssv.visits: XCSAO shift removed for in_stack visits).
Per data set: generalised Lomb-Scargle over 0.05-20 c/d (highest peak, Baluev false-alarm probability) and a sinusoid plus first
harmonic at the adopted frequency (amplitude and time of maximum after BJD 2459300.0).
Adopted frequency: joint fit of the three CoRoT runs (30-min bins, errors from the scatter within each bin) and ZTF g and r, with a
common frequency and phase and one amplitude and offset per data set, on a grid of 2e-7 c/d over 1.6862-1.6872 c/d. The delta chi2
of the nearest cycle-count aliases is printed; chi2 is dominated by the CoRoT bins.
Emission lines: per visit, Gaussian fits with a quadratic baseline to H-alpha (+-1800 km/s) and to the Ca II triplet (8500.35,
8544.44, 8664.52 A; +-1200 km/s; common velocity and width, one amplitude per line); visit time = mean of the TAI start and end
(no barycentric correction, < 0.01 in phase); phase 0 = maximum of the ZTF r fit.
Usage: python reflection_3107374277060584064.py (writes ../tables/reflection_3107374277060584064.csv and _visits.csv)."""
import os, subprocess, numpy as np, pandas as pd
from astropy.io import fits
from astropy.timeseries import LombScargle
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit
from sdssv import CACHE, fetch, visits

GID, NEIGHBOUR, SDSS_ID = "3107374277060584064", "3107374272762041856", "74709777"
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"); TAB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tables")
COROT = {"IRa01": "N2-4.4/2007/02/03/EN2_STAR_MON_0102743730_20070203T130553_20070402T070158.fits",
         "LRa01": "N2-4.4/2007/10/23/EN2_STAR_MON_0102743730_20071023T223035_20080303T093534.fits",
         "LRa06": "N2-4.4/2012/01/12/EN2_STAR_MON_0102743730_20120112T183055_20120329T092714.fits"}
T0 = 2459300.0; C = 299792.458; CAT = [8500.35, 8544.44, 8664.52]


def corot(run):
    path = os.path.join(CACHE, os.path.basename(COROT[run])); os.makedirs(CACHE, exist_ok=True)
    if not (os.path.exists(path) and os.path.getsize(path) > 10000):
        subprocess.run(["curl", "-sL", "-m", "600", "-o", path, f"https://cdsarc.cds.unistra.fr/ftp/B/corot/files/{COROT[run]}"], check=True)
    d = fits.open(path)["BAR"].data; ok = (d["STATUS"] == 0) & (d["WHITEFLUX"] > 0)
    t = d["DATEBARTT"][ok].astype(float) + 2400000.0; y = d["WHITEFLUX"][ok].astype(float); o = np.argsort(t); t, y = t[o], y[o] / np.median(y[o]) - 1
    r = y - median_filter(y, size=int(3 / np.median(np.diff(t))) | 1, mode="nearest"); s = 1.4826 * np.median(np.abs(r)); k = np.abs(r) < 5 * s
    return t[k], r[k]


def binned(t, y, width=1 / 48):
    b = np.floor((t - t.min()) / width); g = pd.DataFrame(dict(t=t, y=y, b=b)).groupby("b").agg(t=("t", "mean"), y=("y", "mean"), n=("y", "size"), s=("y", "std"))
    g = g[g.n >= 3]; return g.t.values, g.y.values, (g.s / np.sqrt(g.n)).values


def ztf(gid):
    d = pd.read_csv(os.path.join(D, f"ztf_{gid}.csv")); d = d[(d.catflags == 0) & (d.magerr < 0.25)]; out = {}
    for band in ("zg", "zr"):
        s = d[d.filtercode == band]
        if len(s) >= 20:
            f = 10 ** (-0.4 * (s.mag.values - np.median(s.mag))) - 1; out[f"ZTF {band}"] = (s.hjd.values, f, 0.921 * s.magerr.values * (f + 1))
    return out


def gaia():
    d = pd.read_csv(os.path.join(D, f"gaia_dr3_epoch_photometry_{GID}.csv")); d = d[(d.GrVFlag == 0) & np.isfinite(d.FG) & np.isfinite(d.TimeG)]
    return d.TimeG.values + 2455197.5, d.FG.values / np.median(d.FG) - 1, d.e_FG.values / np.median(d.FG)


def fit(t, y, e, f):
    x = 2 * np.pi * f * (t - T0); X = np.vstack([np.ones_like(t), np.cos(x), np.sin(x), np.cos(2 * x), np.sin(2 * x)]).T
    w = 1 / e ** 2; A = X.T @ (X * w[:, None]); p = np.linalg.solve(A, X.T @ (w * y)); chi2 = np.sum(w * (y - X @ p) ** 2)
    Cv = np.linalg.inv(A) * max(chi2 / (len(t) - 5), 1); a1 = np.hypot(p[1], p[2])
    return dict(amplitude_frac=a1, e_amplitude_frac=np.sqrt((Cv[1, 1] + Cv[2, 2]) / 2), harmonic2_frac=np.hypot(p[3], p[4]),
                t_max_bjd=T0 + (np.arctan2(p[2], p[1]) % (2 * np.pi)) / (2 * np.pi * f))


def datasets():
    sets = {}
    for run in COROT:
        t, y = corot(run); sets[f"CoRoT {run} (CoRoT 102743730, blend)"] = (t, y, np.full(len(t), 1.4826 * np.median(np.abs(y))))
    sets.update(ztf(GID)); sets["Gaia DR3 epoch photometry G"] = gaia()
    sets.update({f"{k} (neighbour Gaia DR3 {NEIGHBOUR})": v for k, v in ztf(NEIGHBOUR).items()})
    return sets


def joint(sets):
    J = [binned(*sets[k][:2]) for k in sets if k.startswith("CoRoT")] + [sets[k] for k in ("ZTF zg", "ZTF zr")]
    def chi(f):
        tot, ph = 0.0, []
        for t, y, e in J:
            X = np.vstack([np.ones_like(t), np.cos(2 * np.pi * f * (t - T0)), np.sin(2 * np.pi * f * (t - T0))]).T
            p = np.linalg.lstsq(X / e[:, None], y / e, rcond=None)[0]; ph.append((np.hypot(p[1], p[2]), np.arctan2(p[2], p[1]), np.sum(1 / e ** 2)))
        w = np.array([a ** 2 * s for a, _, s in ph]); phi = np.angle(np.sum(w * np.exp(1j * np.array([q for _, q, _ in ph]))))
        for t, y, e in J:
            X = np.vstack([np.ones_like(t), np.cos(2 * np.pi * f * (t - T0) - phi)]).T; p = np.linalg.lstsq(X / e[:, None], y / e, rcond=None)[0]
            tot += np.sum(((y - X @ p) / e) ** 2)
        return tot
    fr = np.arange(1.6862, 1.6872, 2e-7); cs = np.array([chi(f) for f in fr]); k = int(np.argmin(cs))
    mins = [i for i in range(1, len(cs) - 1) if cs[i] < cs[i - 1] and cs[i] < cs[i + 1] and abs(fr[i] - fr[k]) > 1e-4]
    alias = min(mins, key=lambda i: cs[i]) if mins else None
    return fr[k], (fr[alias], cs[alias] - cs[k]) if alias is not None else (np.nan, np.nan), cs[k] / sum(len(j[0]) for j in J)


def model(lines):
    def m(x, *p):
        c = x.mean(); y = p[0] + p[1] * (x - c) + p[2] * (x - c) ** 2
        for k, l in enumerate(lines):
            y = y + p[3 + k] * np.exp(-0.5 * ((x - l * (1 + p[-2] / C)) / (l * p[-1] / C)) ** 2)
        return y
    return m


def emission(f, t_max):
    vs = visits(SDSS_ID); rows = []
    with fits.open(fetch(SDSS_ID, "visit")) as h:
        tm = {int(r["mjd"]): (r["tai_beg"] + r["tai_end"]) / 2 / 86400 + 2400000.5 for i in (1, 2) if h[i].data is not None and len(h[i].data) for r in h[i].data}
    for v in vs:
        row = dict(mjd=v["mjd"], jd_mid=round(tm[v["mjd"]], 4), phase=round(((tm[v["mjd"]] - t_max) * f) % 1, 3), snr=round(v["snr"], 1), in_stack=v["in_stack"])
        for name, lines, win in (("halpha", [6564.61], 1800), ("caii", CAT, 1200)):
            m = np.zeros(len(v["wave"]), bool)
            for l in lines:
                m |= np.abs(v["wave"] / l - 1) * C < win
            x, y, iv = v["wave"][m], v["flux"][m], v["ivar"][m]; ok = (iv > 0) & np.isfinite(y); x, y, iv = x[ok], y[ok], iv[ok]; nl = len(lines)
            p0 = [np.median(y), 0, 0] + [0.3 * np.median(y)] * nl + [0, 150]
            try:
                p, cv = curve_fit(model(lines), x, y, p0=p0, sigma=1 / np.sqrt(iv), absolute_sigma=True, maxfev=20000,
                                  bounds=([-np.inf] * (3 + nl) + [-600, 30], [np.inf] * (3 + nl) + [600, 800]))
                pe = np.sqrt(np.diag(cv)); a, ea = p[3:3 + nl], pe[3:3 + nl]
                row.update({f"{name}_v_kms": round(p[-2]), f"{name}_e_v_kms": round(pe[-2]), f"{name}_sigma_kms": round(p[-1]),
                            f"{name}_amp_snr": round(float(np.sum(a) / np.sqrt(np.sum(ea ** 2))), 1)})
            except Exception:
                row[f"{name}_v_kms"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    sets = datasets(); f, (fa, dchi), chi2r = joint(sets); rows = []
    print(f"adopted frequency {f:.7f} c/d (P = {1 / f:.8f} d = {24 / f:.5f} h); nearest cycle-count alias {fa:.7f} c/d, delta chi2 {dchi:.0f} (chi2 per point {chi2r:.1f})")
    for k, (t, y, e) in sets.items():
        fr = np.arange(0.05, 20, 0.1 / (t.max() - t.min())); ls = LombScargle(t, y, e if not k.startswith("CoRoT") else None); P = ls.power(fr); j = int(np.argmax(P))
        r = fit(t, y, e, f)
        rows.append(dict(dataset=k, n=len(t), bjd_first=round(t.min(), 3), bjd_last=round(t.max(), 3), peak_cd=round(fr[j], 5),
                         peak_fap=float(f"{ls.false_alarm_probability(P[j], minimum_frequency=0.05, maximum_frequency=20, method='baluev'):.2g}"),
                         power_at_adopted=round(float(ls.power(f)), 3), frequency_cd=round(f, 7), period_d=round(1 / f, 8),
                         amplitude_frac=round(r["amplitude_frac"], 4), e_amplitude_frac=round(r["e_amplitude_frac"], 4), harmonic2_frac=round(r["harmonic2_frac"], 4),
                         t_max_bjd=round(r["t_max_bjd"], 4)))
    out = pd.DataFrame(rows); out.to_csv(os.path.join(TAB, f"reflection_{GID}.csv"), index=False)
    t_max = float(out[out.dataset == "ZTF zr"].t_max_bjd.iloc[0])
    vis = emission(f, t_max); vis.to_csv(os.path.join(TAB, f"reflection_{GID}_visits.csv"), index=False)
    pd.set_option("display.width", 250); print(out.to_string()); print(vis.to_string())


if __name__ == "__main__":
    main()
