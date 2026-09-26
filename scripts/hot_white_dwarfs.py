"""Effective temperatures of hot SDSS-V white dwarfs from fits of TMAP H+He NLTE model spectra (TheoSSA) to their SDSS-V DR20 coadds.
Models: the grid in ../data/theossa_hhe_grid.csv (Teff 60-200 kK, log g 6.5-8.0, He mass fraction 0-1; 280 spectra, downloaded to
../data/cache/theossa/). The model parameters are read from each file (text header, or VOTable parameters for files served as VOTable).
Normalised flux, air wavelengths converted to vacuum (Morton 2000), convolved with a Gaussian of FWHM = lambda/1800.
Observed: coadd of the in-stack visits (sdssv.py). Line windows (vacuum centre, half-width in A):
  "H+He":   He I 4027, H-delta, H-gamma, He I 4472, He II 4542, He II 4686, H-beta, He I 4922, He II 5412, He I 5876, H-alpha;
  "He only": the He I and He II windows.
In each window the spectrum is fitted as (a + b x) * model (weighted linear least squares); chi2 is summed over the windows.
Velocity: the 20 models with the lowest chi2 at v = 0 are fitted on a grid of -150..+150 km/s (10 km/s); the best velocity is then used
for the whole grid. Range: models with chi2 <= chi2_min + 6 * max(chi2_r, 1).
Stars and literature values: ../data/hot_white_dwarfs_sources.csv. Writes ../tables/hot_white_dwarfs.csv (one row per star and line set).
python hot_white_dwarfs.py"""
import os, sys, pickle, warnings, numpy as np, pandas as pd, requests
from scipy.ndimage import gaussian_filter1d
from sdssv import visits, coadd, CACHE, C
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "..", "data"); TAB = os.path.join(HERE, "..", "tables")
MDIR = os.path.join(CACHE, "theossa"); VS = np.arange(-150, 151, 10.0); R = 1800
SETS = {"H+He": [(4027.3, 15), (4102.9, 60), (4341.7, 60), (4472.7, 15), (4542.8, 20), (4687.0, 35), (4862.7, 70), (4923.3, 12), (5413.0, 30), (5877.3, 20), (6564.6, 60)],
        "He only": [(4027.3, 15), (4472.7, 15), (4542.8, 20), (4687.0, 35), (4923.3, 12), (5413.0, 30), (5877.3, 20)]}


def read_model(path):
    with open(path) as fh:
        first = fh.readline()
    if first.startswith("<?xml"):
        from astropy.io.votable import parse
        warnings.filterwarnings("ignore"); v = parse(path)
        p = {x.name: x.value for x in v.iter_fields_and_params() if getattr(x, "value", None) is not None}
        t = v.get_first_table().to_table()
        return float(p["t_eff"]), float(p["log_g"]), 1 - float(p["w_H"]), np.asarray(t["spectral"], float), np.asarray(t["flux_norm"], float)
    hd = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("*"):
                break
            if "=" in line:
                k, v = line[1:].strip().split("=", 1); hd[k.strip()] = v.strip()
    d = np.loadtxt(path, comments="*", usecols=(0, 2))
    return float(hd["Teff"]), float(hd["logg"]), float(hd["Mass.Fraction.HE"]), d[:, 0], d[:, 1]


def models():
    pk = os.path.join(CACHE, "theossa_hhe_vacuum_R1800.pkl")
    if os.path.exists(pk):
        return pickle.load(open(pk, "rb"))
    os.makedirs(MDIR, exist_ok=True); M = []
    for _, g in pd.read_csv(os.path.join(DATA, "theossa_hhe_grid.csv")).iterrows():
        path = os.path.join(MDIR, g.file)
        if not os.path.exists(path):
            r = requests.get(g.url, timeout=600); r.raise_for_status(); open(path, "wb").write(r.content)
        teff, logg, xhe, w, n = read_model(path)
        keep = (w > 3900) & (w < 6800); w, n = w[keep], n[keep]
        s2 = (1e4 / w) ** 2; w = w * (1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))
        lg = np.arange(np.log(w[0]), np.log(w[-1]), 2e-5); nn = np.interp(lg, np.log(w), n)
        M.append(dict(teff=teff, logg=logg, xhe=xhe, lg=lg, n=gaussian_filter1d(nn, (1 / R / 2.355) / 2e-5)))
    pickle.dump(M, open(pk, "wb"))
    return M


def spectrum(sdss_id):
    vs = visits(str(sdss_id)); vs = [v for v in vs if v["in_stack"]] or vs
    return coadd(vs)[:3]


def chi2(w, f, iv, m, v, lines):
    tot = 0.0; lw = np.log(w / (1 + v / C))
    for l, hw in lines:
        s = (np.abs(w - l) < hw) & (iv > 0) & np.isfinite(f)
        if s.sum() < 10:
            continue
        mod = np.interp(lw[s], m["lg"], m["n"]); x = (w[s] - l) / hw; sw = np.sqrt(iv[s])
        A = np.vstack([mod, mod * x]).T * sw[:, None]; c, *_ = np.linalg.lstsq(A, f[s] * sw, rcond=None); r = f[s] * sw - A @ c; tot += float(r @ r)
    return tot


def npix(w, f, iv, lines):
    return int(sum(((np.abs(w - l) < hw) & (iv > 0) & np.isfinite(f)).sum() for l, hw in lines if ((np.abs(w - l) < hw) & (iv > 0) & np.isfinite(f)).sum() >= 10))


def fit(w, f, iv, M, lines):
    pre = sorted(M, key=lambda m: chi2(w, f, iv, m, 0.0, lines))[:20]
    v0 = min(((chi2(w, f, iv, m, v, lines), v) for m in pre for v in VS))[1]
    res = sorted(((chi2(w, f, iv, m, v0, lines), i) for i, m in enumerate(M)))
    c0, i0 = res[0]; cr = c0 / npix(w, f, iv, lines); ok = [M[i] for c, i in res if c <= c0 + 6 * max(cr, 1)]
    return M[i0], v0, cr, ok


if __name__ == "__main__":
    M = models(); S = pd.read_csv(os.path.join(DATA, "hot_white_dwarfs_sources.csv"), dtype={"sdss_id": str, "gaia_dr3": str}); rows = []
    for _, s in S.iterrows():
        w, f, iv = spectrum(s.sdss_id)
        for name, lines in SETS.items():
            b, v0, cr, ok = fit(w, f, iv, M, lines); T = [m["teff"] for m in ok]; G = [m["logg"] for m in ok]; X = [m["xhe"] for m in ok]
            rows.append(dict(gaia_dr3=s.gaia_dr3, sdss_id=s.sdss_id, name=s["name"], role=s.role, line_set=name, teff_kK=b["teff"] / 1e3, teff_min_kK=min(T) / 1e3,
                             teff_max_kK=max(T) / 1e3, logg=b["logg"], logg_min=min(G), logg_max=max(G), he_mass_fraction=round(b["xhe"], 2),
                             he_min=round(min(X), 2), he_max=round(max(X), 2), v_kms=v0, chi2_r=round(cr, 2), lit_teff_kK=s.lit_teff_kK,
                             fit_over_lit=round(b["teff"] / 1e3 / s.lit_teff_kK, 2) if np.isfinite(s.lit_teff_kK) else np.nan))
            print(rows[-1], flush=True)
    out = pd.DataFrame(rows); out.to_csv(os.path.join(TAB, "hot_white_dwarfs.csv"), index=False)
    c = out[(out.role == "control")]
    for name in SETS:
        q = c[c.line_set == name].merge(S[["sdss_id", "snr"]], on="sdss_id"); q = q[q.snr >= 19]
        print(f"{name}: controls with S/N >= 19: n = {len(q)}, median fit/lit = {q.fit_over_lit.median():.2f}, 16-84% = {q.fit_over_lit.quantile(.16):.2f}-{q.fit_over_lit.quantile(.84):.2f}")
