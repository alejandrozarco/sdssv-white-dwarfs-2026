"""H-alpha emission in the DESI DR1 spectrum of WDJ205249.27-032419.53 (Gaia DR3 6914922055508553984, TARGETID 39627705435554416) and a test
against the flux-calibration residual at H-alpha noted for DESI spectra by Swan et al. (2026).

- Target spectrum: SPARCL (DESI-DR1). The exposure is a single 494-s bright-time exposure on tile 20836, petal 8 (night 20210523).
- Residual template: all spectra on the same tile and petal (Astro Data Lab desi_dr1.fiberassign and desi_dr1.ztile; SPARCL), keeping
  GALAXY and QSO spectra with no strong line redshifted into 6530-6610 A and a continuum S/N >= 8 in 6500-6550 and 6585-6640 A; each is
  divided by a linear continuum fitted to those sidebands, and the template is the pixel median of the fractional residuals.
- Fit (6520-6610 A, inverse-variance weights): [a + b (lambda - 6565) - Gaussian absorption] x (1 + s x template), with and without an added
  Gaussian emission line; delta chi2 between the two fits.
- Photometric phase of the exposure: ephemeris from ../tables/irradiated_companions.csv (phase 0 = maximum light); the DESI MEAN_MJD is
  converted to BJD_TDB at Kitt Peak.
Writes ../tables/desi_halpha_6914922055508553984.csv and data/cache/desi_halpha_6914922055508553984.npz (spectrum and template, for the figure)."""
import os, io, numpy as np, pandas as pd, requests
from scipy.optimize import least_squares
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

H = os.path.dirname(os.path.abspath(__file__)); CACHE = os.path.join(H, "..", "data", "cache"); os.makedirs(CACHE, exist_ok=True)
TID, GID, TILE, PETAL, C = 39627705435554416, "6914922055508553984", 20836, 8, 299792.458
LINES = [3727.4, 4102.9, 4341.7, 4862.7, 4960.3, 5008.2, 6549.9, 6564.6, 6585.3, 6718.3, 6732.7]


def tap(q):
    r = requests.get("https://datalab.noirlab.edu/tap/sync", params=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q), timeout=600)
    return pd.read_csv(io.StringIO(r.text))


def main():
    from sparcl.client import SparclClient
    cl = SparclClient(connect_timeout=60, read_timeout=600)
    z = tap(f"SELECT targetid, mean_mjd, coadd_exptime, coadd_numexp FROM desi_dr1.zpix WHERE targetid={TID}").iloc[0]
    mates = tap(f"SELECT f.targetid, z.spectype, z.z FROM desi_dr1.fiberassign AS f JOIN desi_dr1.ztile AS z ON z.targetid=f.targetid AND z.tileid=f.tileid "
                f"WHERE f.tileid={TILE} AND f.petal_loc={PETAL}")
    ids = [int(x) for x in mates.targetid]; spec = {}
    for k in range(0, len(ids), 100):
        r = cl.retrieve_by_specid(specid_list=ids[k:k + 100], include=["specid", "wavelength", "flux", "ivar", "redshift", "spectype"], dataset_list=["DESI-DR1"])
        for x in r.records:
            spec[int(x.specid)] = x
    tg = spec[TID]; w0, f0, i0 = np.array(tg.wavelength), np.array(tg.flux), np.array(tg.ivar); m = (w0 > 6520) & (w0 < 6610) & (i0 > 0); w, f, iv = w0[m], f0[m], i0[m]
    R = []
    for sid, x in spec.items():
        if sid == TID or x.spectype not in ("GALAXY", "QSO") or any(6530 < l * (1 + x.redshift) < 6610 for l in LINES):
            continue
        ww, ff, ii = np.array(x.wavelength), np.array(x.flux), np.array(x.ivar); ok = ii > 0; side = ok & (((ww > 6500) & (ww < 6550)) | ((ww > 6585) & (ww < 6640)))
        if side.sum() < 40:
            continue
        p = np.polyfit(ww[side], ff[side], 1); c = np.polyval(p, ww); snr = np.median(c[side]) / np.std(ff[side] - c[side])
        if snr >= 8:
            R.append(np.interp(w, ww, ff / c - 1))
    T = np.median(np.array(R), 0)

    def model(p, em):
        out = (p[0] + p[1] * (w - 6565) - p[2] * np.exp(-0.5 * ((w - p[3]) / p[4]) ** 2)) * (1 + p[5] * T)
        return out + p[6] * np.exp(-0.5 * ((w - p[7]) / p[8]) ** 2) if em else out
    res = lambda p, em: (f - model(p, em)) * np.sqrt(iv)
    p0 = [21, 0, 4, 6563, 12, 1.0]; a = least_squares(res, p0, args=(False,)); b = least_squares(res, p0 + [8, 6568.7, 1.7], args=(True,))
    J = b.jac; cov = np.linalg.inv(J.T @ J) * np.sum(b.fun ** 2) / (len(w) - len(b.x)); e = np.sqrt(np.diag(cov))
    eph = pd.read_csv(os.path.join(H, "..", "tables", "irradiated_companions.csv"), dtype={"gaia_dr3": str}).set_index("gaia_dr3").loc[GID]
    kp = EarthLocation.of_site("kitt peak"); tm = Time(float(z.mean_mjd), format="mjd", scale="utc", location=kp); sc = SkyCoord(eph.ra_deg, eph.dec_deg, unit="deg")
    bjd = (tm.tdb + tm.light_travel_time(sc)).jd; ph = ((bjd - eph.t_max_bjd) * eph.frequency_cd) % 1
    out = dict(gaia_dr3=GID, targetid=TID, tile=TILE, petal=PETAL, mean_mjd=float(z.mean_mjd), exptime_s=float(z.coadd_exptime), n_exp=int(z.coadd_numexp), bjd_tdb=round(bjd, 5),
               phase_from_max=round(ph, 3), n_template_spectra=len(R), template_peak_A=round(float(w[np.argmax(T)]), 2), template_peak_frac=round(float(T.max()), 3),
               chi2_template_only=round(float(np.sum(a.fun ** 2)), 1), chi2_template_plus_emission=round(float(np.sum(b.fun ** 2)), 1), n_pix=len(w),
               template_scale=round(float(b.x[5]), 2), emission_centre_A_vac=round(float(b.x[7]), 2), e_emission_centre_A=round(float(e[7]), 2),
               emission_velocity_kms=round(float((b.x[7] - 6564.61) / 6564.61 * C), 0), e_emission_velocity_kms=round(float(e[7] / 6564.61 * C), 0),
               emission_sigma_A=round(float(abs(b.x[8])), 2), emission_peak_over_continuum=round(float(b.x[6] / b.x[0]), 2))
    pd.DataFrame([out]).to_csv(os.path.join(H, "..", "tables", f"desi_halpha_{GID}.csv"), index=False)
    np.savez(os.path.join(CACHE, f"desi_halpha_{GID}.npz"), w=w, f=f, iv=iv, T=T, m_template=model(a.x, False), m_full=model(b.x, True), w0=w0, f0=f0, i0=i0)
    print(out)


if __name__ == "__main__":
    main()
