"""Carbon screen of SDSS-V visit spectra with a weighted matched filter.
Usage: python carbon_screen.py <sdss_id> [<sdss_id> ...]      (prints one row per spectrum)
       python carbon_screen.py --sample <out.csv>             (selects the sample below from the SnowWhite table and screens it)
       python carbon_screen.py --table                        (screens the objects in ../data/carbon_screen_objects.csv -> ../tables/carbon_screen.csv)
Sample: SnowWhite classification containing DA and not MS; parallax/error > 5; S/N > 5; and (log g >= 8.5 with 10-40 kK, or
M_G > 11.05 + 3.3 (BP-RP + 0.3) + 0.2 with -0.5 < BP-RP < 0.35, or log g >= 9.4).
Per spectrum: inverse-variance coadd of all visits (sdssv.visits, XCSAO shift removed for in-stack visits) on a log grid
3850-9250 A (step 6e-5 dex); depth = 1 - f/cont with cont = running 80th percentile over 25 A, Gaussian-smoothed; Balmer cores
(+-35 A), 5574-5582, 6297-6304, 6860-6960, 7590-7700 and 8940-8990 A masked; weights = ivar * cont^2 capped at their 90th
percentile. Matched filter z(v) = sum(w d T_v) / sqrt(sum(w T_v^2)) for Gaussian templates (FWHM 5 A) of the lines below,
v = -1500..+1500 km/s in 10 km/s steps. Contrast = (z - median) / (1.4826 MAD) with median and MAD over v < -600 or v > +800 km/s;
the reported peak is the maximum inside -200..+400 km/s. 'C' is the combined C I + C II template."""
import sys, numpy as np, pandas as pd
from scipy.ndimage import percentile_filter, gaussian_filter1d
from sdssv import visits, cas, C
GRID = 10 ** np.arange(np.log10(3850), np.log10(9250), 6e-5)
LINES = {"C I": [4773.1, 4933.4, 5053.6, 5381.8, 6014.8, 7115.1, 7117.1, 7118.9, 8337.4, 9063.9, 9091.0, 9097.3, 9114.3],
         "C II": [3921.8, 4077.0, 4268.4, 4375.5, 4620.5, 5146.6, 6579.9, 6584.7, 7233.3, 7238.4],
         "He I": [4027.3, 4389.2, 4472.7, 4714.5, 4923.3, 5017.1, 5877.3, 6679.99, 7067.1]}
LINES["C"] = LINES["C I"] + LINES["C II"]
BALMER = [3890.2, 3971.2, 4102.9, 4341.7, 4862.7, 6564.6]
MASK = [(b - 35, b + 35) for b in BALMER] + [(5574, 5582), (6297, 6304), (6860, 6960), (7590, 7700), (8940, 8990)]
VELS = np.arange(-1500, 1501, 10.0)
INW = (VELS >= -200) & (VELS <= 400); FAR = (VELS < -600) | (VELS > 800)
TT = {}
for sp, L in LINES.items():
    T = np.zeros((len(VELS), len(GRID)))
    for lam in L:
        T += np.exp(-0.5 * ((GRID[None, :] - lam * (1 + VELS / C)[:, None]) / (5 / 2.3548)) ** 2)
    TT[sp] = T


def on_grid(v):
    ok = (v["ivar"] > 0) & np.isfinite(v["flux"])
    return np.interp(GRID, v["wave"][ok], v["flux"][ok], left=np.nan, right=np.nan), np.interp(GRID, v["wave"][ok], v["ivar"][ok], left=0, right=0)


def prep(f, iv):
    m = np.isfinite(f) & (iv > 0)
    if m.sum() < 2000:
        return None
    ff = np.interp(np.arange(len(GRID)), np.where(m)[0], f[m])
    n = int(round(np.log10(1 + 25 / 5000) / 6e-5))
    cont = gaussian_filter1d(percentile_filter(ff, 80, size=n), n / 3)
    d = 1 - ff / cont; w = iv * cont ** 2; good = m.copy()
    for a, b in MASK:
        good &= ~((GRID > a) & (GRID < b))
    w = np.where(good, w, 0); w = np.minimum(w, np.percentile(w[good], 90)); d[~good] = 0
    return d, w


def contrast(d, w, sp):
    z = (TT[sp] @ (w * d)) / np.sqrt((TT[sp] ** 2) @ w)
    med = np.median(z[FAR]); mad = 1.4826 * np.median(np.abs(z[FAR] - med))
    return (z - med) / mad


def screen(sdss_id):
    vs = visits(sdss_id); g = [on_grid(v) for v in vs]
    num = np.nansum([f * iv for f, iv in g], axis=0); den = np.sum([iv for f, iv in g], axis=0)
    p = prep(np.where(den > 0, num / np.where(den > 0, den, 1), np.nan), den)
    if p is None:
        return dict(sdss_id=sdss_id, status="too few pixels")
    row = dict(sdss_id=sdss_id, n_visits=len(vs), snr_max=round(max(v["snr"] for v in vs), 1))
    for sp in ("C", "C I", "C II", "He I"):
        c = contrast(*p, sp); k = int(np.argmax(np.where(INW, c, -99)))
        row[f"{sp.replace(' ', '')}_contrast"] = round(float(c[k]), 2); row[f"{sp.replace(' ', '')}_v_kms"] = int(VELS[k])
    k0 = int(np.argmin(np.abs(VELS - row["C_v_kms"]))); per = []
    for f, iv in g:
        pv = prep(f, iv)
        if pv is not None:
            per.append(f"{contrast(*pv, 'C')[k0]:.1f}")
    row["C_contrast_per_visit"] = "/".join(per); row["status"] = "ok"
    return row


if __name__ == "__main__":
    if sys.argv[1] == "--sample":
        sw = cas("SELECT sdss_id, gaia_dr3_source_id, classification, teff, logg, snr, plx, e_plx, g_mag, bp_mag, rp_mag FROM snow_white_boss_star "
                 "WHERE classification LIKE '%DA%' AND classification NOT LIKE '%MS%' AND plx > 5 * e_plx AND snr > 5")
        bprp = sw.bp_mag - sw.rp_mag; mg = sw.g_mag + 5 * np.log10(sw.plx / 100)
        sel = ((sw.logg >= 8.5) & (sw.teff >= 10000) & (sw.teff <= 40000)) | ((mg > 11.05 + 3.3 * (bprp + 0.3) + 0.2) & (bprp > -0.5) & (bprp < 0.35)) | (sw.logg >= 9.4)
        rows = []
        for sid in sw[sel].sdss_id.astype("int64").astype(str):
            try:
                rows.append(screen(sid))
            except Exception as e:
                rows.append(dict(sdss_id=sid, status=f"hole: {type(e).__name__}"))
        pd.DataFrame(rows).to_csv(sys.argv[2], index=False); print(len(rows), "spectra screened")
    elif sys.argv[1] == "--table":
        obj = pd.read_csv("../data/carbon_screen_objects.csv", dtype={"sdss_id": str, "gaia_dr3": str}).fillna("")
        sw = cas("SELECT sdss_id, classification, g_mag, bp_mag, rp_mag, plx FROM snow_white_boss_star WHERE sdss_id IN (" + ",".join(obj.sdss_id) + ")")
        sw["sdss_id"] = sw.sdss_id.astype("int64").astype(str)
        res = pd.DataFrame([screen(s) for s in obj.sdss_id])
        t = obj.merge(sw, on="sdss_id", how="left").merge(res, on="sdss_id")
        t["bp_rp"] = (t.bp_mag - t.rp_mag).round(3); t["M_G"] = (t.g_mag + 5 * np.log10(t.plx / 100)).round(2); t["G"] = t.g_mag.round(3)
        t = t.rename(columns={"classification": "snowwhite_class"})
        cols = ["gaia_dr3", "sdss_id", "G", "bp_rp", "M_G", "snowwhite_class", "n_visits", "snr_max", "C_contrast", "C_v_kms", "CI_contrast", "CI_v_kms",
                "CII_contrast", "CII_v_kms", "HeI_contrast", "C_contrast_per_visit", "existing_classification", "reference"]
        t.sort_values("C_contrast", ascending=False)[cols].to_csv("../tables/carbon_screen.csv", index=False); print(len(t), "rows")
    else:
        print(pd.DataFrame([screen(s) for s in sys.argv[1:]]).to_string(index=False))
