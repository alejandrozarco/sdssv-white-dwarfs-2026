"""Eclipse ephemeris of Gaia DR3 4851800979770492544 (WDJ030317.61-420658.71) from TESS 2-min SPOC light curves (TIC 651059545,
sectors 97, 105, 106; SAP flux, QUALITY == 0, divided by its median and scaled to the star's share of the aperture with the SPOC
CROWDSAP value) and ATLAS forced photometry (data/atlas_forced_photometry_4851800979770492544.txt; cuts duJy > 0, err == 0, chi/N < 10,
duJy < 3 x median; per-season median subtracted; fractional flux relative to the Gaia synthetic SDSS magnitudes c = (g + r)/2,
o = (r + i)/2, J/A+A/674/A33; times to BJD_TDB).
Period search: a box of 2.5-min half-width is slid over a grid of frequency and phase; the frequency and reference time that maximise the
summed eclipse depth are refined by timing every cycle with at least four points within 6 min (box centre on a 0.05-min grid, only
cycles whose in-box deficit exceeds 0.8 in units of the star's flux); a linear ephemeris is fitted to the mid-times, with the
uncertainties from the scatter of the residuals. Profile: SAP flux folded on the ephemeris in 1-min bins. ATLAS: mean fractional flux
within 2 min of the predicted mid-eclipse against a scan of ephemeris shifts. Usage: python eclipse_4851800979770492544.py (writes
../tables/eclipsing_4851800979770492544.csv and ../tables/eclipsing_4851800979770492544_profile.csv)."""
import os, subprocess, numpy as np, pandas as pd
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import LombScargle
from astropy.coordinates import SkyCoord, EarthLocation
from astroquery.mast import Observations
import astropy.units as u
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"); T = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tables")
GID, TIC, SECTORS = "4851800979770492544", "651059545", (97, 105, 106); RA, DEC = 45.82364258397, -42.11636797807
G, GS, RS, IS = 17.946447, 17.815, 18.1213, 18.3554
# --- TESS
os.makedirs(os.path.join(D, "cache"), exist_ok=True); tt, yy, ss = [], [], []
obs = Observations.query_criteria(obs_collection="TESS", target_name=TIC, dataproduct_type="timeseries", provenance_name="SPOC")
prods = Observations.get_product_list(obs); prods = prods[[str(x).endswith("s_lc.fits") for x in prods["productFilename"]]]
for fn in prods["productFilename"]:
    sec = int(str(fn).split("-s")[1][:4])
    if sec not in SECTORS: continue
    path = os.path.join(D, "cache", str(fn))
    if not os.path.exists(path): subprocess.run(["curl", "-sL", "-m", "600", "-o", path, f"https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS/product/{fn}"], check=True)
    with fits.open(path) as h:
        d = h[1].data; q = (d["QUALITY"] == 0) & np.isfinite(d["SAP_FLUX"]); cs = float(h[1].header["CROWDSAP"])
        t = d["TIME"][q] + 2457000.0; f = d["SAP_FLUX"][q] / np.nanmedian(d["SAP_FLUX"][q]); tt.append(t); yy.append((f - 1) / cs + 1); ss.append(np.full(len(t), sec))
t = np.concatenate(tt); y = np.concatenate(yy); sec = np.concatenate(ss)
def score(f, t0, w=2.5 / 1440):
    ph = ((t - t0) * f + 0.5) % 1 - 0.5; m = np.abs(ph / f) < w
    return (1 - y[m]).sum() / np.sqrt(m.sum()) if m.sum() > 5 else 0.0
fr = np.linspace(2, 60, 580001); ls = LombScargle(t, y); pw = ls.power(fr); f_ls = fr[np.argmax(pw)]
ph = ((t - t.min()) * f_ls) % 1; nb = 100; idx = np.digitize(ph, np.linspace(0, 1, nb + 1)) - 1
prof = np.array([np.nanmean(y[idx == k]) if (idx == k).sum() else np.nan for k in range(nb)]); t0 = t.min() + (np.nanargmin(prof) + 0.5) / nb / f_ls
best = (0.0, f_ls, t0)
for f in np.arange(f_ls - 0.0004, f_ls + 0.0004, 0.000002):
    for dt in np.arange(-3, 3.01, 0.25) / 1440:
        v = score(f, t0 + dt)
        if v > best[0]: best = (v, f, t0 + dt)
_, fb, t0b = best; Pb = 1 / fb; n = np.round((t - t0b) / Pb); times = []
for c in np.unique(n):
    m = (n == c) & (np.abs(t - (t0b + c * Pb)) < 6 / 1440)
    if m.sum() < 4: continue
    grid = np.arange(-3, 3, 0.05) / 1440; sc = []
    for dt in grid:
        ine = np.abs(t[m] - (t0b + c * Pb + dt)) < 2.0 / 1440; sc.append(-(1 - y[m][ine]).sum() if ine.sum() >= 1 else 0.0)
    if min(sc) < -0.8: times.append((c, t0b + c * Pb + grid[int(np.argmin(sc))], int(sec[m][0])))
E = np.array([x[0] for x in times]); TM = np.array([x[1] for x in times])
A = np.vstack([np.ones_like(E), E]).T; coef = np.linalg.lstsq(A, TM, rcond=None)[0]; resid = TM - A @ coef; sig = resid.std(); cov = np.linalg.inv(A.T @ A) * sig ** 2
T0, P = coef; eT0, eP = np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1])
ph = ((t - T0) / P + 0.5) % 1 - 0.5; dtm = ph * P * 1440; bins = np.arange(-10, 10.1, 1.0); idx = np.digitize(dtm, bins) - 1
profile = [dict(minutes_from_mid=round(0.5 * (bins[k] + bins[k + 1]), 1), flux_fraction=round(float(np.nanmedian(y[(idx == k) & (np.abs(dtm) < 10)])), 3), n=int(((idx == k) & (np.abs(dtm) < 10)).sum())) for k in range(len(bins) - 1)]
pf = pd.DataFrame(profile); depth = 1 - pf.flux_fraction.min(); below_half = int((pf.flux_fraction < 0.5).sum()); below_015 = int((pf.flux_fraction < 0.15).sum())
# --- ATLAS
L = [l for l in open(os.path.join(D, f"atlas_forced_photometry_{GID}.txt")).read().splitlines() if l.strip()]; hdr = L[0].lstrip("#").split(); R = [dict(zip(hdr, l.split())) for l in L[1:]]
ok = [x for x in R if float(x["duJy"]) > 0 and float(x["err"]) == 0 and float(x["chi/N"]) < 10]; c0 = SkyCoord(RA * u.deg, DEC * u.deg); geo = EarthLocation.from_geocentric(0, 0, 0, unit="m")
ref = {"c": 3631e6 * 10 ** (-0.4 * (GS + RS) / 2), "o": 3631e6 * 10 ** (-0.4 * (RS + IS) / 2)}; atl = {}
for b in ("c", "o"):
    x = [r for r in ok if r["F"] == b]; med = np.median([float(r["duJy"]) for r in x]); x = [r for r in x if float(r["duJy"]) < 3 * med]
    mjd = np.array([float(r["MJD"]) for r in x]); fl = np.array([float(r["uJy"]) for r in x]); e = np.array([float(r["duJy"]) for r in x])
    season = np.floor((mjd - 57000) / 365.25 + 0.3).astype(int)
    for sv in np.unique(season): fl[season == sv] -= np.median(fl[season == sv])
    tb = Time(mjd, format="mjd", scale="utc", location=geo); bjd = (tb.tdb + tb.light_travel_time(c0)).jd; fr_ = fl / ref[b]
    pha = ((bjd - T0) / P + 0.5) % 1 - 0.5; dm = pha * P * 1440; ine = np.abs(dm) < 2.0; out = np.abs(dm) > 8
    shifts = np.arange(-10, 10.1, 1.0); scan = [float(np.mean(fr_[np.abs(dm - sh) < 2.0])) if (np.abs(dm - sh) < 2.0).sum() > 3 else np.nan for sh in shifts]
    atl[b] = dict(n=len(x), mjd_first=round(mjd.min(), 1), mjd_last=round(mjd.max(), 1), in_eclipse_n=int(ine.sum()), in_eclipse_mean=round(float(np.mean(fr_[ine])), 3), in_eclipse_err=round(float(np.std(fr_[ine]) / np.sqrt(max(ine.sum(), 1))), 3),
                  out_of_eclipse_median=round(float(np.median(fr_[out])), 3), scan_min_shift_min=float(shifts[int(np.nanargmin(scan))]), scan_min=round(float(np.nanmin(scan)), 3))
row = dict(gaia_dr3=GID, name="WDJ030317.61-420658.71", ra_deg=RA, dec_deg=DEC, G=G, tic=TIC, sectors=" ".join(map(str, SECTORS)), n_points=len(t), ls_peak_cd=round(float(f_ls), 5),
           n_eclipses_timed=len(E), cycle_first=int(E.min()), cycle_last=int(E.max()), T0_bjd_tdb=round(T0, 5), e_T0_d=float(f"{eT0:.1g}"), period_d=round(P, 8), e_period_d=float(f"{eP:.1g}"), period_min=round(P * 1440, 5),
           o_minus_c_rms_min=round(sig * 1440, 2), depth_fraction=round(float(depth), 3), minutes_below_half=below_half, minutes_below_0p15=below_015,
           **{f"atlas_{b}_{k}": v for b in atl for k, v in atl[b].items()})
pd.DataFrame([row]).to_csv(os.path.join(T, f"eclipsing_{GID}.csv"), index=False); pf.to_csv(os.path.join(T, f"eclipsing_{GID}_profile.csv"), index=False)
print(pd.Series(row).to_string()); print(pf.to_string(index=False))
