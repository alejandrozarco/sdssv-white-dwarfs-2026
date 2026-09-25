"""Ca II triplet emission of one white dwarf in every available spectrum: SDSS-V visits and coadd, SDSS/BOSS/DESI spectra
from SPARCL (NOIRLab Astro Data Lab) and ESO X-shooter VIS spectra (phase 3 products).
Usage: python gas_disc_epochs.py <gaia_dr3> <sdss_id> <ra> <dec>   -> ../tables/gas_disc_epochs_<gaia_dr3>.csv
Each spectrum is normalised by a quadratic continuum over 8330-8830 A fitted outside +-1100 km/s of the Ca II lines
(vacuum 8500.35, 8544.44, 8664.52 A) and O I 8448.7 A (3 iterations of 3-sigma clipping) and resampled onto the SDSS-V
log grid (ivar scaled by the pixel-size ratio).
ew_A: summed equivalent width of the emission (positive) within +-900 km/s of the three lines, with its formal error.
amp_rel_sdssv: least-squares scale of the SDSS-V coadd profile (normalised flux - 1 within +-1100 km/s of the lines, zero
elsewhere, smoothed by 1 pixel) fitted over 8350-8830 A; 1 = the SDSS-V coadd strength.
v_blue_kms, v_red_kms: velocities of the highest point of the profile (1-pixel smoothing) in -700..-100 and +100..+700 km/s,
averaged over the three lines; given only where ew_A / ew_err_A > 10 and the Gaussian FWHM exceeds 600 km/s (double-peaked).
v_gauss_kms, fwhm_gauss_kms: common centroid and FWHM of a fit of one Gaussian per line (shared velocity and width, free
amplitudes, linear continuum) over 8350-8830 A; given where ew_A / ew_err_A > 5. Not corrected for instrumental broadening.
X-shooter wavelengths are converted from air to vacuum (Morton 2000) and are topocentric; SDSS, BOSS, DESI and SDSS-V are
heliocentric vacuum wavelengths."""
import sys, os, json, subprocess, numpy as np, pandas as pd, requests, io
from astropy.io import fits
from astropy.time import Time
from scipy.ndimage import gaussian_filter1d
from sdssv import visits, CACHE, C
from gas_disc_screen import GRID, CAT, on_grid
M = (GRID > 8250) & (GRID < 8950); W = GRID[M]
WIN = np.zeros(len(W), bool)
for _l in CAT:
    WIN |= np.abs(W / _l - 1) * C < 1100


def air_to_vac(w):
    s2 = (1e4 / w) ** 2
    return w * (1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))


def norm(w, f, iv):
    w = np.asarray(w, float); f = np.asarray(f, float); iv = np.asarray(iv, float)
    sel = (w > 8330) & (w < 8830) & (iv > 0) & np.isfinite(f); use = sel.copy()
    for l in list(CAT) + [8448.7]:
        use &= np.abs(w / l - 1) * C > 1100
    for _ in range(3):
        p = np.polyfit(w[use], f[use], 2, w=np.sqrt(iv[use])); P = np.polyval(p, w); use &= np.abs((f - P) * np.sqrt(iv)) < 3
    n = np.interp(W, w[sel], (f / P)[sel], left=np.nan, right=np.nan); v = np.interp(W, w[sel], (iv * P ** 2)[sel], left=0, right=0)
    v *= np.median(np.diff(W)) / np.median(np.diff(w[sel]))
    return n, np.where(np.isfinite(n), v, 0)


def measure(n, v, T):
    dl = np.gradient(W); win = np.zeros(len(W), bool)
    for l in CAT:
        win |= np.abs(W / l - 1) * C < 900
    ok = win & (v > 0) & np.isfinite(n); ew = np.sum(((n - 1) * dl)[ok]); e = np.sqrt(np.sum((dl ** 2 / v)[ok]))
    ok = (v > 0) & np.isfinite(n) & (np.abs(W - 8590) < 240)
    A = np.sum(v[ok] * T[ok] * (n[ok] - 1)) / np.sum(v[ok] * T[ok] ** 2); eA = 1 / np.sqrt(np.sum(v[ok] * T[ok] ** 2))
    return round(float(ew), 2), round(float(e), 2), round(float(A), 3), round(float(eA), 3)


def peaks(n):
    vb, vr = [], []
    for l in CAT:
        v = (W / l - 1) * C; y = gaussian_filter1d(np.nan_to_num(n, nan=1.0), 1)
        b = (v > -700) & (v < -100); r = (v > 100) & (v < 700); vb.append(v[b][np.argmax(y[b])]); vr.append(v[r][np.argmax(y[r])])
    return int(round(np.mean(vb), -1)), int(round(np.mean(vr), -1))


def gauss_fit(n, v):
    from scipy.optimize import curve_fit
    ok = (v > 0) & np.isfinite(n) & (W > 8350) & (W < 8830)
    def g3(w, a1, a2, a3, vel, sg, c0, c1):
        y = c0 + c1 * (w - 8600) / 300
        for a, l in zip((a1, a2, a3), CAT):
            y = y + a * np.exp(-0.5 * ((w - l * (1 + vel / C)) / (l * sg / C)) ** 2)
        return y
    p, cv = curve_fit(g3, W[ok], n[ok], p0=[0.3, 0.4, 0.3, 0, 150, 1, 0], sigma=1 / np.sqrt(v[ok]), maxfev=20000)
    return int(round(p[3])), int(round(2.3548 * abs(p[4])))


def sparcl_spectra(ra, dec, r=3.0):
    from sparcl.client import SparclClient
    d = r / 3600
    for k in range(5):
        try:
            cl = SparclClient(); break
        except Exception:
            if k == 4:
                raise
    f = cl.find(outfields=["sparcl_id"], constraints={"ra": [ra - d / np.cos(np.radians(dec)), ra + d / np.cos(np.radians(dec))], "dec": [dec - d, dec + d]}, limit=50)
    ids = [x["sparcl_id"] for x in f.records]
    if not ids:
        return []
    recs = cl.retrieve(uuid_list=ids, include=["flux", "ivar", "wavelength", "data_release", "specid", "dateobs"]).records
    out = []
    for x in recs:
        dt = x["dateobs"]
        if isinstance(dt, str) and dt.startswith("["):
            dt = json.loads(dt)
        dt = str(dt[0] if isinstance(dt, (list, tuple)) else dt)
        out.append(dict(dataset=x["data_release"], identifier=str(x["specid"]), date=str(dt)[:10], mjd=round(Time(str(dt)[:19].replace(" ", "T")).mjd, 3),
                        w=np.array(x["wavelength"]), f=np.array(x["flux"]), iv=np.array(x["ivar"])))
    return out


def eso_xshooter(ra, dec, arm="VIS"):
    q = ("SELECT dp_id, t_min, em_min FROM ivoa.ObsCore WHERE instrument_name='XSHOOTER' AND "
         f"CONTAINS(POINT('ICRS',s_ra,s_dec),CIRCLE('ICRS',{ra},{dec},0.003))=1")
    t = pd.read_csv(io.StringIO(requests.post("https://archive.eso.org/tap_obs/sync", data=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q), timeout=300).text))
    lo = {"UVB": (2.9e-7, 3.1e-7), "VIS": (5.2e-7, 5.5e-7)}[arm]; t = t[(t.em_min > lo[0]) & (t.em_min < lo[1])]
    out = []
    for r in t.itertuples():
        path = os.path.join(CACHE, f"{r.dp_id}.fits")
        if not os.path.exists(path):
            subprocess.run(["curl", "-sL", "-m", "600", "-o", path, f"https://dataportal.eso.org/dataportal_new/file/{r.dp_id}"], check=True)
        d = fits.open(path)[1].data; e = np.asarray(d["ERR"][0], float)
        out.append(dict(dataset=f"ESO XSHOOTER {arm}", identifier=r.dp_id, date=Time(r.t_min, format="mjd").iso[:10], mjd=round(float(r.t_min), 3),
                        w=air_to_vac(np.asarray(d["WAVE"][0], float) * 10), f=np.asarray(d["FLUX"][0], float), iv=np.where(e > 0, 1 / np.where(e > 0, e, 1) ** 2, 0)))
    return sorted(out, key=lambda x: x["mjd"])


def sdssv_spectra(sdss_id):
    vs = [v for v in visits(sdss_id) if ((v["ivar"] > 0) & np.isfinite(v["flux"])).sum() >= 1000]
    g = [on_grid(v) for v in vs]; num = np.nansum([f * iv for f, iv in g], axis=0); den = np.sum([iv for _, iv in g], axis=0)
    co = dict(dataset="SDSS-V coadd", identifier=f"sdss_id {sdss_id}, {len(vs)} visits", date="", mjd=np.nan, w=GRID, f=np.where(den > 0, num / np.where(den > 0, den, 1), np.nan), iv=den)
    per = [dict(dataset="SDSS-V visit", identifier=f"sdss_id {sdss_id}", date=Time(v["mjd"], format="mjd").iso[:10], mjd=v["mjd"], w=GRID, f=f, iv=iv) for v, (f, iv) in zip(vs, g)]
    return co, per


def all_spectra(gaia, sdss_id, ra, dec):
    co, per = sdssv_spectra(sdss_id)
    return co, per, sparcl_spectra(ra, dec), eso_xshooter(ra, dec)


if __name__ == "__main__":
    gaia, sdss_id, ra, dec = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
    co, per, sp, xs = all_spectra(gaia, sdss_id, ra, dec)
    nc, vc = norm(co["w"], co["f"], co["iv"]); T = np.where(WIN, gaussian_filter1d(np.nan_to_num(nc - 1), 1), 0.0)
    rows = []
    for s in [co] + per + sp + xs:
        n, v = norm(s["w"], s["f"], s["iv"]); ew, e, A, eA = measure(n, v, T)
        vg, fw = gauss_fit(n, v) if ew > 5 * e else ("", "")
        vb, vr = peaks(n) if (ew > 10 * e and fw != "" and fw > 600) else ("", "")
        rows.append(dict(gaia_dr3=gaia, dataset=s["dataset"], identifier=s["identifier"], date_utc=s["date"], mjd=s["mjd"], ew_A=ew, ew_err_A=e,
                         amp_rel_sdssv=A, amp_err=eA, v_blue_kms=vb, v_red_kms=vr, v_gauss_kms=vg, fwhm_gauss_kms=fw))
    t = pd.DataFrame(rows); t.to_csv(f"../tables/gas_disc_epochs_{gaia}.csv", index=False); print(t.to_string())
