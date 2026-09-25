"""SDSS-V DR20 helpers: Gaia DR3 -> sdss_id lookup (SnowWhite table), Astra 0.8.1 mwmStar/mwmVisit download and reading.
Visit spectra: the XCSAO velocity shift is removed for visits with in_stack = True (wavelength = grid * (1 + v/c));
visits with in_stack = False are used as delivered."""
import io, os, time, requests, numpy as np, pandas as pd
from astropy.io import fits
C = 299792.458
SAS = "https://data.sdss.org/sas/dr20/spectro/astra/0.8.1/spectra"
CAS = "https://skyserver.sdss.org/dr20/SkyServerWS/SearchTools/SqlSearch"
CACHE = os.path.join(os.path.dirname(__file__), "..", "data", "cache")


def cas(sql, dtype=None):
    r = requests.get(CAS, params=dict(cmd=sql, format="csv"), timeout=300)
    r.raise_for_status()
    t = r.text.split("\n", 1)[1] if r.text.startswith("#Table") else r.text
    return pd.read_csv(io.StringIO(t), dtype=dtype)


def snowwhite(gaia_ids):
    ids = ",".join(str(int(g)) for g in gaia_ids)
    return cas("SELECT sdss_id, gaia_dr3_source_id, ra, dec, g_mag, plx, classification, p_da, p_dah, teff, logg, snr, telescope "
               f"FROM snow_white_boss_star WHERE gaia_dr3_source_id IN ({ids})")


def fetch(sdss_id, kind="star"):
    sid = str(sdss_id); name = f"mwm{'Star' if kind == 'star' else 'Visit'}-0.8.1-{sid}.fits"
    path = os.path.join(CACHE, name); os.makedirs(CACHE, exist_ok=True)
    if not (os.path.exists(path) and os.path.getsize(path) > 1000):
        url = f"{SAS}/{kind}/{sid[-4:-2]}/{sid[-2:]}/{name}"
        for k in range(3):
            r = requests.get(url, timeout=180)
            if r.status_code == 200 and len(r.content) > 1000:
                open(path, "wb").write(r.content); break
            time.sleep(3)
        else:
            raise IOError(f"download failed: {url}")
    return path


def star_spectrum(sdss_id):
    """Highest-S/N coadd (APO or LCO HDU): wavelength, flux, ivar, snr."""
    best = None
    with fits.open(fetch(sdss_id, "star")) as h:
        for i in (1, 2):
            if h[i].data is None or len(h[i].data) == 0:
                continue
            d = h[i].data[0]
            if best is None or d["snr"] > best[3]:
                best = (np.array(d["wavelength"], float), np.array(d["flux"], float), np.array(d["ivar"], float), float(d["snr"]))
    return best


def visits(sdss_id):
    out = []
    with fits.open(fetch(sdss_id, "visit")) as h:
        for i in (1, 2):
            if h[i].data is None or len(h[i].data) == 0:
                continue
            hd = h[i].header; wg = 10 ** (hd["CRVAL"] + hd["CDELT"] * np.arange(hd["NPIXELS"]))
            for r in h[i].data:
                v = float(r["xcsao_v_rad"]); ins = bool(r["in_stack"])
                w = wg * (1 + v / C) if (ins and np.isfinite(v)) else wg
                out.append(dict(mjd=int(r["mjd"]), xcsao_v=v, in_stack=ins, snr=float(r["snr"]), wave=w,
                                flux=np.array(r["flux"], float), ivar=np.array(r["ivar"], float)))
    return out


def coadd(vs, lo=3700, hi=9300, dlog=6e-5):
    grid = 10 ** np.arange(np.log10(lo), np.log10(hi), dlog); num = np.zeros_like(grid); den = np.zeros_like(grid)
    for v in vs:
        ok = (v["ivar"] > 0) & np.isfinite(v["flux"])
        f = np.interp(grid, v["wave"][ok], v["flux"][ok], left=np.nan, right=np.nan)
        iv = np.interp(grid, v["wave"][ok], v["ivar"][ok], left=0, right=0)
        m = np.isfinite(f); num[m] += f[m] * iv[m]; den[m] += iv[m]
    return grid, np.where(den > 0, num / np.where(den > 0, den, 1), np.nan), den
