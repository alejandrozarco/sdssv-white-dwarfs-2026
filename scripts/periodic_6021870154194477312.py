"""Periodic signal of Gaia DR3 6021870154194477312 (GALEX J161854.1-355427) in ATLAS, Gaia DR3 and TESS photometry.

ATLAS forced photometry (data/atlas_forced_photometry_<gaia_dr3>.txt) at the Gaia DR3 positions of the star and of the three Gaia DR3
sources within 13 arcsec (data/periodic_6021870154194477312_sources.csv). Cuts: duJy > 0, err == 0, chi/N < 10, duJy < 3 x median;
per-season median subtracted; 5-sigma clip; times converted to BJD_TDB.
Frequency: generalised Lomb-Scargle of the combined c and o fluxes of the star (each band divided by its reference flux) over
0.5-50 c/d, then a least-squares sinusoid with separate band offsets on a fine grid; the 1-sigma range is where chi2 <= chi2_min + chi2_r.
Amplitudes: sinusoid plus first harmonic at the adopted frequency, per band. Fractional amplitudes use reference fluxes from the
Gaia-synthesised SDSS magnitudes of the Gaia DR3 white dwarf catalogue (J/A+A/674/A33: g 17.7409, r 18.1185, i 18.5434):
c = (g + r)/2 and o = (r + i)/2 in flux.
Gaia DR3 epoch photometry (VizieR I/355/epphot; data/gaia_dr3_epoch_photometry_6021870154194477312.csv): transits with a rejection
flag removed; TimeG + 2455197.5 is used as BJD (the TCB-TDB offset of ~19 s is ignored).
TESS: SPOC PDCSAP light curve of TIC 1251484163, sector 65 (120 s), QUALITY == 0, 5-sigma clip. PDCSAP includes the SPOC crowding
correction (CROWDSAP), so its amplitude depends on that correction.
Phase: t_max is the first maximum of the fitted fundamental after T0 = BJD 2458000.0, using the adopted frequency for all data
sets. Usage: python periodic_6021870154194477312.py (writes ../tables/periodic_6021870154194477312.csv)."""
import os, subprocess, numpy as np, pandas as pd
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import LombScargle
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

GID = "6021870154194477312"; TIC, SECTOR = "1251484163", 65; T0 = 2458000.0
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
GEO = EarthLocation.from_geocentric(0, 0, 0, unit="m")
fl = lambda m: 3631e6 * 10 ** (-0.4 * m)
REF = {"c": (fl(17.7409) + fl(18.1185)) / 2, "o": (fl(18.1185) + fl(18.5434)) / 2}
SRC = pd.read_csv(os.path.join(D, f"periodic_{GID}_sources.csv"), dtype={"gaia_dr3": str}).set_index("gaia_dr3")


def load_atlas(gid):
    ra, dec = SRC.loc[gid, "ra_deg"], SRC.loc[gid, "dec_deg"]; c0 = SkyCoord(ra * u.deg, dec * u.deg)
    L = [l for l in open(os.path.join(D, f"atlas_forced_photometry_{gid}.txt")).read().splitlines() if l.strip()]
    hdr = L[0].lstrip("#").split(); R = [dict(zip(hdr, l.split())) for l in L[1:]]
    ok = [x for x in R if float(x["duJy"]) > 0 and float(x["err"]) == 0 and float(x["chi/N"]) < 10]; out = {}
    for b in ("c", "o"):
        s = [x for x in ok if x["F"] == b]; med = np.median([float(x["duJy"]) for x in s]); s = [x for x in s if float(x["duJy"]) < 3 * med]
        mjd = np.array([float(x["MJD"]) for x in s]); f = np.array([float(x["uJy"]) for x in s]); e = np.array([float(x["duJy"]) for x in s])
        season = np.floor((mjd - 57000) / 365.25 + 0.3).astype(int)
        for sv in np.unique(season):
            f[season == sv] -= np.median(f[season == sv])
        clip = np.abs(f) < 5 * 1.4826 * np.median(np.abs(f)) + 3 * np.median(e)
        t = Time(mjd[clip], format="mjd", scale="utc", location=GEO); out[b] = dict(t=(t.tdb + t.light_travel_time(c0)).jd, f=f[clip], e=e[clip])
    return out


def load_gaia():
    d = pd.read_csv(os.path.join(D, f"gaia_dr3_epoch_photometry_{GID}.csv")); out = {}
    for b, tc, fc, ec, flag in (("G", "TimeG", "FG", "e_FG", "GrVFlag"), ("BP", "TimeBP", "FBP", "e_FBP", "BPrVFlag"), ("RP", "TimeRP", "FRP", "e_FRP", "RPrVFlag")):
        s = d[(d[flag] == 0) & np.isfinite(d[fc]) & np.isfinite(d[tc])]; med = np.median(s[fc])
        out[b] = dict(t=s[tc].values + 2455197.5, f=s[fc].values / med - 1, e=s[ec].values / med)
    return out


def load_tess():
    from astroquery.mast import Observations
    obs = Observations.query_criteria(obs_collection="TESS", target_name=TIC, sequence_number=SECTOR, dataproduct_type="timeseries", provenance_name="SPOC")
    fn = [n for n in Observations.get_product_list(obs)["productFilename"] if n.endswith("s_lc.fits")][0]
    path = os.path.join(D, "cache", fn); os.makedirs(os.path.dirname(path), exist_ok=True)
    if not os.path.exists(path):
        subprocess.run(["curl", "-sL", "-m", "600", "-o", path, f"https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS/product/{fn}"], check=True)
    h = fits.open(path); d = h[1].data
    m = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]) & np.isfinite(d["TIME"]); t = d["TIME"][m]; y = d["PDCSAP_FLUX"][m] / np.median(d["PDCSAP_FLUX"][m]) - 1
    e = d["PDCSAP_FLUX_ERR"][m] / np.median(d["PDCSAP_FLUX"][m]); mad = 1.4826 * np.median(np.abs(y - np.median(y))); ok = np.abs(y - np.median(y)) < 5 * mad
    return dict(t=t[ok] + 2457000.0, f=y[ok], e=e[ok]), float(h[1].header.get("CROWDSAP"))


def sinefit(t, y, e, f, harm=2, groups=None):
    x = 2 * np.pi * f * (t - T0); cols = [np.ones_like(t)]
    if groups is not None:
        cols += [(groups == g).astype(float) for g in np.unique(groups)[1:]]
    k0 = len(cols)
    for h in range(1, harm + 1):
        cols += [np.sin(h * x), np.cos(h * x)]
    X = np.vstack(cols).T; w = 1 / e ** 2; A = X.T @ (X * w[:, None]); p = np.linalg.solve(A, X.T @ (w * y)); C = np.linalg.inv(A)
    chi2 = float(np.sum(w * (y - X @ p) ** 2)); s2 = max(chi2 / (len(t) - len(p)), 1.0); C = C * s2
    res = dict(chi2=chi2, dof=len(t) - len(p))
    for h in range(1, harm + 1):
        a, b = p[k0 + 2 * (h - 1)], p[k0 + 2 * (h - 1) + 1]; ca, cb = C[k0 + 2 * (h - 1), k0 + 2 * (h - 1)], C[k0 + 2 * (h - 1) + 1, k0 + 2 * (h - 1) + 1]
        res[f"amp{h}"] = float(np.hypot(a, b)); res[f"e_amp{h}"] = float(np.sqrt((ca + cb) / 2))
        if h == 1:
            phi = np.arctan2(a, b) % (2 * np.pi); res["t_max"] = T0 + phi / (2 * np.pi * f)
            res["e_t_max"] = float(np.sqrt((ca + cb) / 2) / np.hypot(a, b) / (2 * np.pi * f))
    return res


def adopted_frequency(at):
    t = np.concatenate([at[b]["t"] for b in ("c", "o")]); y = np.concatenate([at[b]["f"] / REF[b] for b in ("c", "o")])
    e = np.concatenate([at[b]["e"] / REF[b] for b in ("c", "o")]); g = np.concatenate([[b] * len(at[b]["t"]) for b in ("c", "o")])
    fr = np.arange(0.5, 50, 0.2 / (t.max() - t.min())); ls = LombScargle(t, y, e); P = ls.power(fr); k = int(np.argmax(P))
    fap = float(ls.false_alarm_probability(P[k], minimum_frequency=0.5, maximum_frequency=50, method="baluev"))
    fg = np.arange(fr[k] - 2e-4, fr[k] + 2e-4, 2e-7); chi = np.array([sinefit(t, y, e, f, harm=1, groups=g)["chi2"] for f in fg])
    j = int(np.argmin(chi)); s2 = chi[j] / (len(t) - 4); inside = fg[chi <= chi[j] + s2]
    return dict(f=float(fg[j]), e_f=float((inside.max() - inside.min()) / 2), gls_f=float(fr[k]), gls_power=float(P[k]), gls_fap=fap), (fr, P)


def main():
    at = load_atlas(GID); F, _ = adopted_frequency(at); f = F["f"]; rows = []
    print(f"adopted frequency {f:.7f} +- {F['e_f']:.7f} c/d, P = {1440 / f:.4f} min (GLS peak {F['gls_f']:.5f} c/d, power {F['gls_power']:.3f}, Baluev FAP {F['gls_fap']:.1e})")
    for gid in SRC.index:
        d = at if gid == GID else load_atlas(gid)
        sep = SkyCoord(SRC.loc[gid, "ra_deg"] * u.deg, SRC.loc[gid, "dec_deg"] * u.deg).separation(SkyCoord(SRC.loc[GID, "ra_deg"] * u.deg, SRC.loc[GID, "dec_deg"] * u.deg)).arcsec
        for b in ("c", "o"):
            r = sinefit(d[b]["t"], d[b]["f"], d[b]["e"], f)
            row = dict(gaia_dr3=gid, separation_arcsec=round(sep, 1), G=SRC.loc[gid, "G"], dataset="ATLAS", band=b, n=len(d[b]["t"]),
                       bjd_first=round(d[b]["t"].min(), 3), bjd_last=round(d[b]["t"].max(), 3), frequency_cd=round(f, 7), e_frequency_cd=round(F["e_f"], 7),
                       amplitude_uJy=round(r["amp1"], 2), e_amplitude_uJy=round(r["e_amp1"], 2), harmonic2_uJy=round(r["amp2"], 2), e_harmonic2_uJy=round(r["e_amp2"], 2),
                       t_max_bjd=round(r["t_max"], 5), e_t_max_min=round(r["e_t_max"] * 1440, 1))
            if gid == GID:
                row.update(amplitude_frac=round(r["amp1"] / REF[b], 4), e_amplitude_frac=round(r["e_amp1"] / REF[b], 4),
                           harmonic2_frac=round(r["amp2"] / REF[b], 4), e_harmonic2_frac=round(r["e_amp2"] / REF[b], 4))
            rows.append(row)
    ga = load_gaia()
    for b, d in ga.items():
        r = sinefit(d["t"], d["f"], d["e"], f, harm=1); p1 = LombScargle(d["t"], d["f"], d["e"]).power(np.array([f]))[0]
        rows.append(dict(gaia_dr3=GID, separation_arcsec=0.0, G=SRC.loc[GID, "G"], dataset="Gaia DR3 epoch photometry", band=b, n=len(d["t"]),
                         bjd_first=round(d["t"].min(), 3), bjd_last=round(d["t"].max(), 3), frequency_cd=round(f, 7), amplitude_frac=round(r["amp1"], 4),
                         e_amplitude_frac=round(r["e_amp1"], 4), t_max_bjd=round(r["t_max"], 5), e_t_max_min=round(r["e_t_max"] * 1440, 1), gls_power_at_f=round(float(p1), 3)))
    tt, crowd = load_tess(); r = sinefit(tt["t"], tt["f"], tt["e"], f, harm=1)
    fr = np.arange(0.5, 50, 0.001); ls = LombScargle(tt["t"], tt["f"]); P = ls.power(fr); k = int(np.argmax(P))
    rows.append(dict(gaia_dr3=GID, separation_arcsec=0.0, G=SRC.loc[GID, "G"], dataset=f"TESS S{SECTOR} PDCSAP (TIC {TIC}, CROWDSAP {crowd:.3f})", band="TESS",
                     n=len(tt["t"]), bjd_first=round(tt["t"].min(), 3), bjd_last=round(tt["t"].max(), 3), frequency_cd=round(f, 7), amplitude_frac=round(r["amp1"], 4),
                     e_amplitude_frac=round(r["e_amp1"], 4), t_max_bjd=round(r["t_max"], 5), e_t_max_min=round(r["e_t_max"] * 1440, 1),
                     tess_gls_peak_cd=round(float(fr[k]), 4),
                     tess_gls_peak_fap=float(f"{ls.false_alarm_probability(P[k], minimum_frequency=0.5, maximum_frequency=50, method='baluev'):.2g}")))
    out = pd.DataFrame(rows); out.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tables", f"periodic_{GID}.csv"), index=False)
    pd.set_option("display.width", 250); print(out.to_string())


if __name__ == "__main__":
    main()
