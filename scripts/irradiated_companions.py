"""Short-period white dwarfs with a red-rising photometric modulation and an infrared excess (data/irradiated_companions_sources.csv).

Light curves:
- ZTF DR (data/ztf_<gaia_dr3>.csv; IRSA light-curve service, 1.5 arcsec): catflags == 0, magerr < 0.25; each ZTF object/filter light
  curve converted to fractional flux about its median; light curves with fewer than 15 points dropped.
- Gaia DR3 epoch photometry (data/gaia_dr3_epoch_photometry_<gaia_dr3>.csv; VizieR I/355/epphot): G, BP and RP transits without the
  variability-rejection flag; fractional flux about the median; TimeG/BP/RP + 2455197.5 used as BJD.
- ATLAS forced photometry (data/atlas_forced_photometry_<gaia_dr3>.txt; c and o bands) where it exists: duJy > 0, err == 0, chi/N < 10,
  duJy below three times the band median; per-season (365.25-d) median flux subtracted; fractional flux relative to the G-band flux
  3631e6 x 10^(-0.4 G) uJy; MJD (UTC) converted to BJD_TDB at the star's position.
- TESS: SPOC 120-s PDCSAP light curves (QUALITY == 0) where they exist; otherwise TESScut full-frame-image cutouts (7 x 7 pixels):
  3 x 3-pixel aperture on the target pixel, per-cadence background = median of the outer ring of pixels, times 9. A 1-day running median
  is subtracted and points beyond 5 sigma are clipped. PDCSAP fractions include the SPOC crowding correction; FFI fractions are relative
  to the total aperture flux and are not corrected for other stars in the aperture.
Per data set: the highest Lomb-Scargle peak (0.5-50 c/d; 2-60 c/d for TESS) with its Baluev false-alarm probability, and the semi-amplitude
of the fundamental of a sinusoid-plus-first-harmonic fit at the adopted frequency.
Adopted frequency: common-phase sinusoid over ZTF, TESS and Gaia G; each data set is divided by its own semi-amplitude (independent fit at
the highest combined Lomb-Scargle peak, 2-50 c/d) and has its own offset. Candidate frequencies are the six highest Lomb-Scargle maxima of
the scaled data within +-0.02 c/d; each is refined on a 1e-7 c/d grid and the lowest chi2 is adopted. The uncertainty is the half-range with chi2 <= chi2_min + chi2_r; next_alias_delta_chi2 is
the chi2 difference to the best other candidate. t_max: first maximum of the joint
fundamental after BJD 2459000.0.
Infrared: CatWISE2020 W1/W2 (VizieR II/365, 2 arcsec) and VHS DR5 J/Ks (VizieR II/367, 1.5 arcsec), Vega magnitudes. The white-dwarf
prediction is the Montreal pure-H synthetic photometry (Holberg & Bergeron 2006; Bedard et al. 2020; Table_DA) interpolated at the Gentile
Fusillo et al. (2021) H-atmosphere Teff and log g (log g clipped to 7.0-9.0), scaled to the Gaia G magnitude: m_pred = G - (G3 - m)_model.
Zero points (Jy): W1 309.54, W2 171.787, J 1594, Ks 666.7 (2MASS values used for VISTA). The excess flux is converted to an absolute
magnitude of the companion with d = 1/parallax.
Usage: python irradiated_companions.py (writes ../tables/irradiated_companions.csv and ../tables/irradiated_companions_periods.csv)."""
import os, io, numpy as np, pandas as pd, requests, warnings
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.timeseries import LombScargle
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
from astroquery.mast import Observations, Tesscut
from scipy.ndimage import median_filter

warnings.filterwarnings("ignore")
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"); CACHE = os.path.join(D, "cache"); os.makedirs(CACHE, exist_ok=True)
SRC = pd.read_csv(os.path.join(D, "irradiated_companions_sources.csv"), dtype={"gaia_dr3": str, "tic": str}).set_index("gaia_dr3")
TABLE_DA = "https://www.astro.umontreal.ca/~bergeron/CoolingModels/Tables/Table_DA"
ZP = {"W1": 309.54, "W2": 171.787, "J": 1594.0, "Ks": 666.7}
V = Vizier(columns=["**"], row_limit=-1); V.TIMEOUT = 300


def gaia_info(gid):
    g = V.query_constraints(catalog="I/355/gaiadr3", Source=gid)[0][0]
    w = V.query_constraints(catalog="J/MNRAS/508/3877/maincat", GaiaEDR3=gid)[0][0]
    return dict(ra=float(g["RA_ICRS"]), dec=float(g["DE_ICRS"]), G=float(g["Gmag"]), bp_rp=float(g["BP-RP"]), plx=float(g["Plx"]), e_plx=float(g["e_Plx"]),
                wdj=str(w["WDJname"]), teff=float(w["TeffH"]), logg=float(w["loggH"]), mass=float(w["MassH"]))


def gaia_epochs(gid):
    p = os.path.join(D, f"gaia_dr3_epoch_photometry_{gid}.csv")
    if not os.path.exists(p):
        q = V.query_constraints(catalog="I/355/epphot", Source=gid)
        if not q:
            return {}
        q[0].to_pandas().to_csv(p, index=False)
    e = pd.read_csv(p); out = {}
    for b, tc, fc, ec, fl in (("G", "TimeG", "FG", "e_FG", "GrVFlag"), ("BP", "TimeBP", "FBP", "e_FBP", "BPrVFlag"), ("RP", "TimeRP", "FRP", "e_FRP", "RPrVFlag")):
        m = np.isfinite(e[tc]) & np.isfinite(e[fc]) & (e[fc] > 0) & (e[fl] == 0)
        f = e[fc][m].values; med = np.median(f); out[f"Gaia {b}"] = (e[tc][m].values + 2455197.5, f / med - 1, e[ec][m].values / med)
    return out


def ztf(gid):
    p = os.path.join(D, f"ztf_{gid}.csv")
    if not os.path.exists(p):
        return {}
    d = pd.read_csv(p); d = d[(d.catflags == 0) & (d.magerr < 0.25)]; out = {}
    for band in ("zg", "zr"):
        t, y, e = [], [], []
        for _, s in d[d.filtercode == band].groupby("oid"):
            if len(s) < 15:
                continue
            fl = 10 ** (-0.4 * (s.mag.values - np.median(s.mag.values))) - 1
            t += list(s.hjd.values); y += list(fl - np.mean(fl)); e += list(0.921 * s.magerr.values * (fl + 1))
        if t:
            out[f"ZTF {band}"] = tuple(map(np.array, (t, y, e)))
    return out


def atlas(gid, ra, dec, G):
    p = os.path.join(D, f"atlas_forced_photometry_{gid}.txt")
    if not os.path.exists(p):
        return {}
    L = [l for l in open(p).read().splitlines() if l.strip()]; hdr = L[0].lstrip("#").split(); R = [dict(zip(hdr, l.split())) for l in L[1:]]
    ok = [x for x in R if float(x["duJy"]) > 0 and float(x["err"]) == 0 and float(x["chi/N"]) < 10]; ref = 3631e6 * 10 ** (-0.4 * G); out = {}
    c0 = SkyCoord(ra * u.deg, dec * u.deg); geo = EarthLocation.from_geodetic(lon=-156.2569 * u.deg, lat=20.7075 * u.deg, height=3055 * u.m)
    for b in ("c", "o"):
        x = [r for r in ok if r["F"] == b]
        if len(x) < 15:
            continue
        med = np.median([float(r["duJy"]) for r in x]); x = [r for r in x if float(r["duJy"]) < 3 * med]
        mjd = np.array([float(r["MJD"]) for r in x]); f = np.array([float(r["uJy"]) for r in x]); e = np.array([float(r["duJy"]) for r in x])
        season = np.floor((mjd - 57000) / 365.25 + 0.3).astype(int)
        for sv in np.unique(season):
            f[season == sv] -= np.median(f[season == sv])
        clip = np.abs(f) < 5 * 1.4826 * np.median(np.abs(f)) + 3 * np.median(e)
        t = Time(mjd[clip], format="mjd", scale="utc", location=geo); bjd = (t.tdb + t.light_travel_time(c0)).jd
        out[f"ATLAS {b}"] = (np.array(bjd), f[clip] / ref, e[clip] / ref)
    return out


def detrend(t, y):
    o = np.argsort(t); t, y = t[o], y[o]; w = max(int(1.0 / np.median(np.diff(t))) | 1, 3); r = y - median_filter(y, w, mode="nearest")
    s = 1.4826 * np.median(np.abs(r)); m = np.abs(r) < 5 * s
    return t[m], r[m], np.full(m.sum(), s)


def tess_spoc(tic, sectors):
    out = {}; obs = Observations.query_criteria(target_name=tic, obs_collection="TESS", dataproduct_type="timeseries")
    for sec in sectors:
        o = obs[(obs["sequence_number"] == sec) & (obs["t_exptime"] == 120)]
        if not len(o):
            continue
        pl = Observations.get_product_list(o); pl = pl[[x.endswith("_lc.fits") for x in pl["productFilename"]]]
        path = Observations.download_products(pl, download_dir=CACHE)["Local Path"][0]
        d = fits.open(path)[1].data; q = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"])
        f = d["PDCSAP_FLUX"][q]; out[f"TESS S{sec}"] = detrend(d["TIME"][q] + 2457000.0, f / np.median(f) - 1)
    return out


def tess_ffi(ra, dec, sectors):
    out = {}; c = SkyCoord(ra, dec, unit="deg")
    for sec in sectors:
        p = os.path.join(CACHE, f"tesscut_{ra:.5f}_{dec:.5f}_s{sec}_7x7.fits")
        if not os.path.exists(p):
            Tesscut.get_cutouts(coordinates=c, size=7, sector=sec)[0].writeto(p, overwrite=True)
        h = fits.open(p); d = h[1].data
        x, y = WCS(h[2].header).world_to_pixel(c); xi, yi = int(round(float(x))), int(round(float(y)))
        q = (d["QUALITY"] == 0) & np.isfinite(d["TIME"]); fl = d["FLUX"][q]; t = d["TIME"][q] + 2457000.0
        yy, xx = np.mgrid[0:fl.shape[1], 0:fl.shape[2]]; ring = (np.abs(xx - xi) > 1) | (np.abs(yy - yi) > 1)
        ap = (np.abs(xx - xi) <= 1) & (np.abs(yy - yi) <= 1); lc = np.nansum(fl[:, ap], axis=1) - np.nanmedian(fl[:, ring], axis=1) * ap.sum()
        ok = np.isfinite(lc); out[f"TESS S{sec} FFI"] = detrend(t[ok], lc[ok] / np.median(lc[ok]) - 1)
    return out


def harmonic_fit(t, y, e, f, t0=2459000.0):
    X = np.vstack([np.ones_like(t)] + [fn(2 * np.pi * k * f * (t - t0)) for k in (1, 2) for fn in (np.cos, np.sin)]).T
    W = X / e[:, None]; c = np.linalg.lstsq(W, y / e, rcond=None)[0]; chi = np.sum(((y - X @ c) / e) ** 2) / max(len(t) - 5, 1)
    cov = np.linalg.inv(W.T @ W) * max(chi, 1); a = np.hypot(c[1], c[2])
    ea = np.sqrt((c[1] ** 2 * cov[1, 1] + c[2] ** 2 * cov[2, 2]) / a ** 2) if a > 0 else np.nan
    return a, ea


def adopted_frequency(sets, f_guess):
    """Common-phase fit: each data set is scaled by its own semi-amplitude (independent fit at f_guess) and offset; the scaled data share one
    sinusoid. Candidates are the six highest Lomb-Scargle maxima of the scaled data within +-0.02 c/d; each is refined on a 1e-7 c/d grid."""
    names = [k for k in sets if k not in ("Gaia BP", "Gaia RP")]; t0 = 2459000.0; T, Y, E, I = [], [], [], []
    for i, k in enumerate(names):
        t, y, e = sets[k]; a, _ = harmonic_fit(t, y, e, f_guess)
        if not np.isfinite(a) or a <= 0:
            continue
        T.append(t); Y.append((y - np.median(y)) / a); E.append(e / a); I.append(np.full(len(t), i))
    t, y, e, gi = map(np.concatenate, (T, Y, E, I)); ids = np.unique(gi); G = (gi[:, None] == ids[None, :]).astype(float)

    def chi2(f):
        X = np.hstack([G, np.cos(2 * np.pi * f * (t - t0))[:, None], np.sin(2 * np.pi * f * (t - t0))[:, None]])
        c = np.linalg.lstsq(X / e[:, None], y / e, rcond=None)[0]; return np.sum(((y - X @ c) / e) ** 2), c
    fr = np.arange(f_guess - 0.02, f_guess + 0.02, 1e-6); pw = LombScargle(t, y, e).power(fr)
    loc = sorted([i for i in range(1, len(fr) - 1) if pw[i] >= pw[i - 1] and pw[i] >= pw[i + 1]], key=lambda i: -pw[i]); cand = []
    for i in loc:
        if all(abs(fr[i] - c) > 2e-5 for c in cand):
            cand.append(fr[i])
        if len(cand) == 6:
            break
    res = []
    for fc in cand:
        grid = np.arange(fc - 1e-5, fc + 1e-5, 1e-7); ch = np.array([chi2(f)[0] for f in grid]); k = int(np.argmin(ch)); res.append((ch[k], grid[k], grid, ch))
    res.sort(key=lambda x: x[0]); c2, f0, grid, ch = res[0]; chir = c2 / (len(t) - len(ids) - 2); ok = grid[ch <= c2 + chir]; err = max((ok.max() - ok.min()) / 2, 1e-7)
    alias_f, alias_dchi2 = (res[1][1], res[1][0] - c2) if len(res) > 1 else (np.nan, np.nan)
    _, c = chi2(f0); ph = np.arctan2(c[-1], c[-2]); tmax = t0 + (ph / (2 * np.pi)) / f0; tmax += np.ceil((t0 - tmax) * f0) / f0
    return f0, err, tmax, alias_f, alias_dchi2, chir


def wd_colours():
    p = os.path.join(CACHE, "Table_DA")
    if not os.path.exists(p):
        open(p, "w").write(requests.get(TABLE_DA, timeout=120).text)
    L = open(p).read().splitlines(); hdr = L[1].replace("log g", "logg").split()
    M = np.array([[float(x) for x in l.split()] for l in L[2:] if len(l.split()) == len(hdr)]); ix = {k: hdr.index(k) for k in ("G3", "W1", "W2", "Ks")}; ix["J"] = hdr.index("J")

    def colour(teff, logg, band):
        lg = min(max(logg, 7.0), 9.0); lo = min(np.floor(lg * 2) / 2, 8.5); hi = lo + 0.5; v = []
        for l in (lo, hi):
            q = M[np.isclose(M[:, 1], l)]; q = q[np.argsort(q[:, 0])]; v.append(np.interp(teff, q[:, 0], q[:, ix["G3"]] - q[:, ix[band]]))
        return float(np.interp(lg, [lo, hi], v))
    return colour


def infrared(info, colour):
    c = SkyCoord(info["ra"], info["dec"], unit="deg"); obs = {}
    w = V.query_region(c, radius=2 * u.arcsec, catalog="II/365/catwise")
    if len(w):
        r = w[0][0]; obs.update(W1=(float(r["W1mproPM"]), float(r["e_W1mproPM"])), W2=(float(r["W2mproPM"]), float(r["e_W2mproPM"])))
    v = V.query_region(c, radius=1.5 * u.arcsec, catalog="II/367/vhs_dr5")
    if len(v):
        r = v[0][0]
        for b, col in (("J", "Jpmag"), ("Ks", "Kspmag")):
            if not np.ma.is_masked(r[col]) and np.isfinite(float(r[col])):
                obs[b] = (float(r[col]), float(r["e_" + col]))
    dm = 5 * np.log10(1000 / info["plx"] / 10); out = {}
    for b, (m, e) in obs.items():
        pred = info["G"] - colour(info["teff"], info["logg"], b); fo = ZP[b] * 10 ** (-0.4 * m) * 1e6; fp = ZP[b] * 10 ** (-0.4 * pred) * 1e6; ex = fo - fp
        out.update({f"{b}_obs": round(m, 3), f"e_{b}_obs": round(e, 3), f"{b}_wd_model": round(pred, 3), f"{b}_ratio": round(fo / fp, 2),
                    f"{b}_excess_uJy": round(ex, 1), f"e_{b}_excess_uJy": round(fo * e / 1.0857, 1),
                    f"M_{b}_companion": round(-2.5 * np.log10(ex * 1e-6 / ZP[b]) - dm, 2) if ex > 0 else np.nan})
    return out


def main():
    colour = wd_colours(); rows, prow = [], []
    for gid, s in SRC.iterrows():
        info = gaia_info(gid); sets = {}
        sets.update(ztf(gid)); sets.update(atlas(gid, info["ra"], info["dec"], info["G"]))
        secs = [int(x) for x in str(s.tess_sectors).split()]
        sets.update(tess_spoc(s.tic, secs) if s.tess_mode == "spoc" else tess_ffi(info["ra"], info["dec"], secs))
        gaia = gaia_epochs(gid)
        tt = np.concatenate([sets[k][0] for k in sets]); yy = np.concatenate([sets[k][1] for k in sets]); ee = np.concatenate([sets[k][2] for k in sets])
        fr = np.linspace(2, 50, 480001); f_guess = fr[np.argmax(LombScargle(tt, yy, ee).power(fr))]
        f0, ef, tmax, af, adchi, chir = adopted_frequency({**sets, **({"Gaia G": gaia["Gaia G"]} if "Gaia G" in gaia else {})}, f_guess); allsets = {**sets, **gaia}
        for k, (t, y, e) in allsets.items():
            lo = 2 if k.startswith("TESS") else 0.5; fq = np.linspace(lo, 50, 400001); ls = LombScargle(t, y, e); pw = ls.power(fq); j = np.argmax(pw)
            a, ea = harmonic_fit(t, y, e, f0)
            prow.append(dict(gaia_dr3=gid, dataset=k, n=len(t), peak_cd=round(float(fq[j]), 5), peak_fap=float(f"{ls.false_alarm_probability(pw[j], minimum_frequency=lo, maximum_frequency=50):.2g}"),
                             semi_amplitude_pct=round(100 * a, 2), e_semi_amplitude_pct=round(100 * ea, 2)))
        amp = {r["dataset"]: r["semi_amplitude_pct"] for r in prow if r["gaia_dr3"] == gid}
        red = amp.get("ZTF zr", np.nan) / amp.get("ZTF zg", np.nan) if "ZTF zr" in amp else (amp.get("Gaia RP", np.nan) / amp.get("Gaia BP", np.nan) if "Gaia RP" in amp else amp.get("ATLAS o", np.nan) / amp.get("ATLAS c", np.nan))
        rows.append(dict(gaia_dr3=gid, name=s["name"], ra_deg=round(info["ra"], 6), dec_deg=round(info["dec"], 6), G=round(info["G"], 3), bp_rp=round(info["bp_rp"], 3),
                         parallax_mas=round(info["plx"], 3), e_parallax_mas=round(info["e_plx"], 3), distance_pc=round(1000 / info["plx"]), M_G=round(info["G"] + 5 * np.log10(info["plx"] / 100), 2),
                         gf21_teff_H=round(info["teff"]), gf21_logg_H=round(info["logg"], 2), gf21_mass_H=round(info["mass"], 3),
                         frequency_cd=round(f0, 7), e_frequency_cd=float(f"{ef:.2g}"), period_min=round(1440 / f0, 3), t_max_bjd=round(tmax, 5),
                         next_alias_cd=round(af, 6), next_alias_delta_chi2=round(adchi, 1), chi2_r=round(chir, 2),
                         red_to_blue_amplitude=round(red, 2), **infrared(info, colour)))
        print(gid, "done", flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(D, "..", "tables", "irradiated_companions.csv"), index=False)
    pd.DataFrame(prow).to_csv(os.path.join(D, "..", "tables", "irradiated_companions_periods.csv"), index=False)


if __name__ == "__main__":
    main()
