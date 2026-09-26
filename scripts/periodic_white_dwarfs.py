"""Photometric periods of white dwarfs with a Gaia DR3 GLS frequency (data/periodic_white_dwarfs_sources.csv).

Ground-based light curves:
- ATLAS forced photometry (data/atlas_forced_photometry_<gaia_dr3>.txt; positions propagated to 2020.5): cuts duJy > 0, err == 0,
  chi/N < 10, duJy < 3 x median; per-season median subtracted; 5-sigma clip; fractional flux relative to the Gaia synthetic SDSS
  magnitudes (VizieR J/A+A/674/A33, white-dwarf table), c = (g + r)/2 and o = (r + i)/2 in flux; where the star is not in that table,
  the Gaia G flux is used for both bands.
- ZTF DR light curves (data/ztf_<gaia_dr3>.csv; IRSA light-curve service, 2 arcsec): catflags == 0; each ZTF object/filter light curve
  converted to fractional flux about its median; light curves with fewer than 20 points dropped.
Frequency: generalised Lomb-Scargle over 0.05-50 c/d (Baluev false-alarm probability of the highest peak), then a least-squares
sinusoid with one offset per light curve on a fine grid; the uncertainty is the half-range where chi2 <= chi2_min + chi2_r.
Amplitudes: sinusoid plus first harmonic at the adopted frequency; for ZTF also per filter (rows "ZTF zg", "ZTF zr").
Gaia DR3 epoch photometry (VizieR I/355/epphot; data/gaia_dr3_epoch_photometry_<gaia_dr3>.csv): G transits without a rejection flag;
TimeG + 2455197.5 used as BJD. TESS (where SPOC light curves exist): PDCSAP, QUALITY == 0, 5-sigma clip, highest peak over 0.2-50 c/d.
t_max: first maximum of the fitted fundamental after T0 = BJD 2458000.0 at the adopted frequency.
Usage: python periodic_white_dwarfs.py (writes ../tables/periodic_white_dwarfs.csv)."""
import os, subprocess, numpy as np, pandas as pd
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import LombScargle
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
T0 = 2458000.0; GEO = EarthLocation.from_geocentric(0, 0, 0, unit="m")
SRC = pd.read_csv(os.path.join(D, "periodic_white_dwarfs_sources.csv"), dtype={"gaia_dr3": str}).set_index("gaia_dr3")
TESS = {"2883364038621038208": ("705345754", (87, 98)), "2888030331609338240": ("705508671", (98,)), "6639666736903611136": ("201655627", (27, 67, 94, 103, 104)),
        "974895286283420160": ("407569944", (20, 47, 60)), "2795150147707769728": ("611449439", (57,))}
fl = lambda m: 3631e6 * 10 ** (-0.4 * m)


def load_atlas(gid):
    s = SRC.loc[gid]; c0 = SkyCoord(s.ra_deg * u.deg, s.dec_deg * u.deg)
    ref = {"c": (fl(s.g_sdss_syn) + fl(s.r_sdss_syn)) / 2, "o": (fl(s.r_sdss_syn) + fl(s.i_sdss_syn)) / 2} if np.isfinite(s.g_sdss_syn) else {"c": fl(s.G), "o": fl(s.G)}
    L = [l for l in open(os.path.join(D, f"atlas_forced_photometry_{gid}.txt")).read().splitlines() if l.strip()]
    hdr = L[0].lstrip("#").split(); R = [dict(zip(hdr, l.split())) for l in L[1:]]
    ok = [x for x in R if float(x["duJy"]) > 0 and float(x["err"]) == 0 and float(x["chi/N"]) < 10]; T, Y, E, G = [], [], [], []
    for b in ("c", "o"):
        x = [r for r in ok if r["F"] == b]; med = np.median([float(r["duJy"]) for r in x]); x = [r for r in x if float(r["duJy"]) < 3 * med]
        mjd = np.array([float(r["MJD"]) for r in x]); f = np.array([float(r["uJy"]) for r in x]); e = np.array([float(r["duJy"]) for r in x])
        season = np.floor((mjd - 57000) / 365.25 + 0.3).astype(int)
        for sv in np.unique(season):
            f[season == sv] -= np.median(f[season == sv])
        clip = np.abs(f) < 5 * 1.4826 * np.median(np.abs(f)) + 3 * np.median(e)
        t = Time(mjd[clip], format="mjd", scale="utc", location=GEO)
        T += list((t.tdb + t.light_travel_time(c0)).jd); Y += list(f[clip] / ref[b]); E += list(e[clip] / ref[b]); G += [f"ATLAS {b}"] * int(clip.sum())
    return np.array(T), np.array(Y), np.array(E), np.array(G)


def load_ztf(gid):
    d = pd.read_csv(os.path.join(D, f"ztf_{gid}.csv")); d = d[d.catflags == 0]; T, Y, E, G = [], [], [], []
    s = SRC.loc[gid]; c0 = SkyCoord(s.ra_deg * u.deg, s.dec_deg * u.deg)
    for (oid, band), x in d.groupby(["oid", "filtercode"]):
        if len(x) < 20:
            continue
        t = Time(x.mjd.values, format="mjd", scale="utc", location=GEO)
        T += list((t.tdb + t.light_travel_time(c0)).jd); Y += list(10 ** (-0.4 * (x.mag.values - np.median(x.mag))) - 1)
        E += list(0.921 * x.magerr.values); G += [f"ZTF {band} {oid}"] * len(x)
    return np.array(T), np.array(Y), np.array(E), np.array(G)


def load_gaia(gid):
    d = pd.read_csv(os.path.join(D, f"gaia_dr3_epoch_photometry_{gid}.csv")); d = d[(d.GrVFlag == 0) & np.isfinite(d.FG) & np.isfinite(d.TimeG)]
    return d.TimeG.values + 2455197.5, d.FG.values / np.median(d.FG) - 1, d.e_FG.values / np.median(d.FG)


def sinefit(t, y, e, f, groups=None, harm=2):
    x = 2 * np.pi * f * (t - T0); cols = [np.ones_like(t)]
    if groups is not None:
        cols += [(groups == g).astype(float) for g in np.unique(groups)[1:]]
    k0 = len(cols)
    for h in range(1, harm + 1):
        cols += [np.sin(h * x), np.cos(h * x)]
    X = np.vstack(cols).T; w = 1 / e ** 2; A = X.T @ (X * w[:, None]); p = np.linalg.solve(A, X.T @ (w * y))
    chi2 = float(np.sum(w * (y - X @ p) ** 2)); C = np.linalg.inv(A) * max(chi2 / (len(t) - len(p)), 1.0); res = dict(chi2=chi2)
    for h in range(1, harm + 1):
        i = k0 + 2 * (h - 1); a, b = p[i], p[i + 1]
        res[f"amp{h}"] = float(np.hypot(a, b)); res[f"e_amp{h}"] = float(np.sqrt((C[i, i] + C[i + 1, i + 1]) / 2))
        if h == 1:
            res["t_max"] = T0 + (np.arctan2(a, b) % (2 * np.pi)) / (2 * np.pi * f); res["e_t_max"] = res["e_amp1"] / res["amp1"] / (2 * np.pi * f)
    return res


def adopted_frequency(t, y, e, g):
    fr = np.arange(0.05, 50, 0.2 / (t.max() - t.min())); ls = LombScargle(t, y, e); P = ls.power(fr); k = int(np.argmax(P))
    fap = float(ls.false_alarm_probability(P[k], minimum_frequency=0.05, maximum_frequency=50, method="baluev"))
    step = 1 / (t.max() - t.min()); fg = np.arange(fr[k] - 2 * step, fr[k] + 2 * step, step / 200)
    chi = np.array([sinefit(t, y, e, f, g, harm=1)["chi2"] for f in fg]); j = int(np.argmin(chi))
    s2 = chi[j] / (len(t) - len(np.unique(g)) - 2); inside = fg[chi <= chi[j] + s2]
    return dict(f=float(fg[j]), e_f=float((inside.max() - inside.min()) / 2), gls_f=float(fr[k]), gls_fap=fap), (fr, P)


def tess(gid, f):
    if gid not in TESS:
        return []
    from astroquery.mast import Observations
    tic, sectors = TESS[gid]; rows = []
    obs = Observations.query_criteria(obs_collection="TESS", target_name=tic, dataproduct_type="timeseries", provenance_name="SPOC")
    for fn in sorted({n for n in Observations.get_product_list(obs)["productFilename"] if n.endswith("s_lc.fits")}):
        path = os.path.join(D, "cache", fn); os.makedirs(os.path.dirname(path), exist_ok=True)
        if not os.path.exists(path):
            subprocess.run(["curl", "-sL", "-m", "600", "-o", path, f"https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS/product/{fn}"], check=True)
        h = fits.open(path); d = h[1].data; sec = h[0].header["SECTOR"]
        if sec not in sectors:
            continue
        m = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]); t = d["TIME"][m] + 2457000.0; y = d["PDCSAP_FLUX"][m] / np.median(d["PDCSAP_FLUX"][m]) - 1
        e = d["PDCSAP_FLUX_ERR"][m] / np.median(d["PDCSAP_FLUX"][m]); mad = 1.4826 * np.median(np.abs(y - np.median(y))); ok = np.abs(y - np.median(y)) < 5 * mad
        t, y, e = t[ok], y[ok], e[ok]; fr = np.arange(0.2, 50, 0.001); ls = LombScargle(t, y); P = ls.power(fr); k = int(np.argmax(P)); r = sinefit(t, y, e, f, harm=1)
        rows.append(dict(dataset=f"TESS S{sec} PDCSAP (TIC {tic}, CROWDSAP {h[1].header['CROWDSAP']:.3f})", n=len(t), bjd_first=round(t.min(), 3), bjd_last=round(t.max(), 3),
                         peak_cd=round(float(fr[k]), 4), peak_fap=float(f"{ls.false_alarm_probability(P[k], minimum_frequency=0.2, maximum_frequency=50, method='baluev'):.2g}"),
                         amplitude_frac=round(r["amp1"], 4), e_amplitude_frac=round(r["e_amp1"], 4), t_max_bjd=round(r["t_max"], 5), e_t_max_min=round(r["e_t_max"] * 1440, 1)))
    return rows


def main():
    rows = []
    for gid, s in SRC.iterrows():
        t, y, e, g = load_atlas(gid) if s.ground == "ATLAS" else load_ztf(gid)
        F, _ = adopted_frequency(t, y, e, g); f = F["f"]; r = sinefit(t, y, e, f, g)
        base = dict(gaia_dr3=gid, name=s["name"], frequency_cd=round(f, 7), e_frequency_cd=round(F["e_f"], 7), period_h=round(24 / f, 5), e_period_h=round(24 * F["e_f"] / f ** 2, 5))
        rows.append(dict(base, dataset=s.ground, n=len(t), bjd_first=round(t.min(), 3), bjd_last=round(t.max(), 3), peak_cd=round(F["gls_f"], 5), peak_fap=float(f"{F['gls_fap']:.2g}"),
                         amplitude_frac=round(r["amp1"], 4), e_amplitude_frac=round(r["e_amp1"], 4), harmonic2_frac=round(r["amp2"], 4), e_harmonic2_frac=round(r["e_amp2"], 4),
                         t_max_bjd=round(r["t_max"], 5), e_t_max_min=round(r["e_t_max"] * 1440, 1)))
        if s.ground == "ZTF":
            for band in ("zg", "zr"):
                m = np.array([x.startswith(f"ZTF {band}") for x in g])
                if m.sum() < 20:
                    continue
                rb = sinefit(t[m], y[m], e[m], f, g[m])
                rows.append(dict(base, dataset=f"ZTF {band}", n=int(m.sum()), bjd_first=round(t[m].min(), 3), bjd_last=round(t[m].max(), 3), amplitude_frac=round(rb["amp1"], 4),
                                 e_amplitude_frac=round(rb["e_amp1"], 4), harmonic2_frac=round(rb["amp2"], 4), e_harmonic2_frac=round(rb["e_amp2"], 4),
                                 t_max_bjd=round(rb["t_max"], 5), e_t_max_min=round(rb["e_t_max"] * 1440, 1)))
        tg, yg, eg = load_gaia(gid); rg = sinefit(tg, yg, eg, f, harm=1)
        rows.append(dict(base, dataset="Gaia DR3 epoch photometry G", n=len(tg), bjd_first=round(tg.min(), 3), bjd_last=round(tg.max(), 3), peak_cd=round(s.gaia_gls_freq_cd, 5),
                         peak_fap=float(f"{s.gaia_gls_fap:.2g}"), amplitude_frac=round(rg["amp1"], 4), e_amplitude_frac=round(rg["e_amp1"], 4), t_max_bjd=round(rg["t_max"], 5),
                         e_t_max_min=round(rg["e_t_max"] * 1440, 1)))
        for tr in tess(gid, f):
            rows.append(dict(base, **tr))
    out = pd.DataFrame(rows); out.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tables", "periodic_white_dwarfs.csv"), index=False)
    pd.set_option("display.width", 250); print(out.to_string())


if __name__ == "__main__":
    main()
