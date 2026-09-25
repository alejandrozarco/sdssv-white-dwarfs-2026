"""Ca II triplet emission (gaseous debris disc) screen of SDSS-V DR20 white dwarf spectra.
Usage: python gas_disc_screen.py --sample              pass 1: every object with a SnowWhite classification (about 51,000;
                                                       each mwmVisit file is downloaded, reduced to the regions below and deleted)
       python gas_disc_screen.py --pass2 <out.csv>     pass 2: template subtraction and matched filter for every stored object
       python gas_disc_screen.py --table <pass2.csv>   writes ../tables/gas_disc_screen.csv for ../data/gas_disc_objects.csv
Pass 1: visits with the XCSAO shift removed for in-stack visits (sdssv.visits), inverse-variance coadd on a log grid
(3850-9250 A, 6e-5 dex). The coadd over REG and each visit over VREG are stored in ../data/cache/gas_disc_store/.
Pass 2, Ca II triplet region 8250-8950 A (vacuum lines 8500.35, 8544.44, 8664.52 A):
  n = f / P2, a quadratic continuum fitted with ivar weights outside +-1100 km/s of the Ca II lines and O I 8448.7 A,
  3 iterations of 3-sigma clipping.
  Sample: parallax/error > 3 and M_G > 8.5, plus all objects with parallax/error <= 3 (M_G = 11.2 + 3.3 (BP-RP) for the
  neighbour search). Groups: H (SnowWhite class starting DA, without MS), nonH (other white-dwarf classes), MS, CV.
  Template = pixel median of n over the 100 nearest neighbours in (BP-RP / 0.04, M_G / 0.2) of the same group with
  coadd S/N > 10 in the region. Model n = b + a (template - 1) fitted outside the masked windows; r = n - model;
  weights w = ivar P2^2 / s^2 with s = max(1, 1.4826 MAD of r sqrt(ivar P2^2) outside the windows).
  Matched filter z = sum(w r T) / sqrt(sum(w T^2)) for emission templates summed over the three lines: single Gaussians
  (sigma 60, 150 km/s) and double-peaked profiles (half-separation vp = 150, 250, 350, 500 km/s, peak sigma max(70, 0.45 vp)),
  centre velocity -400..+400 km/s in 20 km/s steps. z_cat is the maximum; also per-line z, the summed equivalent width over
  +-700 km/s windows, z_fake (same templates with the line set shifted by -2900, +2900, +7500 km/s) and z_abs (maximum of -z).
  z per visit at the coadd best template and velocity."""
import sys, os, time, glob, requests, numpy as np, pandas as pd, warnings
from concurrent.futures import ThreadPoolExecutor
from scipy.spatial import cKDTree
from sdssv import visits, cas, CACHE, C
warnings.filterwarnings("ignore")
GRID = 10 ** np.arange(np.log10(3850), np.log10(9250), 6e-5)
CAT = np.array([8500.35, 8544.44, 8664.52]); OI = 8448.7
REG = [(3880, 4000), (4780, 4940), (4990, 5200), (5820, 5940), (6400, 6750), (7740, 7810), (8250, 8950)]
VREG = [(6400, 6750), (8250, 8950)]
STORE = os.path.join(CACHE, "gas_disc_store")
VEL = np.arange(-400, 401, 20.0)
PROF = [("n60", 0, 60), ("n150", 0, 150)] + [(f"dp{vp}", vp, max(70, 0.45 * vp)) for vp in (150, 250, 350, 500)]
MC = (GRID > 8250) & (GRID < 8950); WC = GRID[MC]; XC = (WC - 8600) / 350
WIN = np.zeros(len(WC), bool)
for _lam in list(CAT) + [OI]:
    WIN |= np.abs(WC / _lam - 1) * C < 1100


def tmpl(lines, vp, sig, wave):
    T = np.zeros((len(VEL), len(wave)))
    for lam in lines:
        for s in ((0,) if vp == 0 else (-vp, vp)):
            mu = lam * (1 + (VEL + s) / C)[:, None]; sg = lam * sig / C
            T += np.exp(-0.5 * ((wave[None, :] - mu) / sg) ** 2)
    return T


TT = {n: tmpl(CAT, vp, sg, WC) for n, vp, sg in PROF}
TL = {n: [tmpl([l], vp, sg, WC) for l in CAT] for n, vp, sg in PROF}
FAKE = {d: {n: tmpl(CAT * (1 + d / C), vp, sg, WC) for n, vp, sg in PROF} for d in (-2900, 2900, 7500)}


def on_grid(v):
    ok = (v["ivar"] > 0) & np.isfinite(v["flux"])
    return np.interp(GRID, v["wave"][ok], v["flux"][ok], left=np.nan, right=np.nan), np.interp(GRID, v["wave"][ok], v["ivar"][ok], left=0, right=0)


def store(sid):
    vs = [on_grid(v) + (v["mjd"], v["snr"]) for v in visits(sid) if ((v["ivar"] > 0) & np.isfinite(v["flux"])).sum() >= 1000]
    if not vs:
        return "no usable visit"
    num = np.nansum([f * iv for f, iv, _, _ in vs], axis=0); den = np.sum([iv for _, iv, _, _ in vs], axis=0)
    f = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan); d = {}
    for a, b in REG:
        m = (GRID > a) & (GRID < b); d[f"f_{a}"] = f[m].astype(np.float32); d[f"iv_{a}"] = den[m].astype(np.float32)
    for a, b in VREG:
        m = (GRID > a) & (GRID < b)
        d[f"vf_{a}"] = np.array([x[0][m] for x in vs], np.float32); d[f"viv_{a}"] = np.array([x[1][m] for x in vs], np.float32)
    d["mjd"] = np.array([x[2] for x in vs]); d["vsnr"] = np.array([x[3] for x in vs], np.float32)
    os.makedirs(os.path.join(STORE, sid[-2:]), exist_ok=True); np.savez_compressed(os.path.join(STORE, sid[-2:], f"{sid}.npz"), **d)
    return "ok"


def norm(f, iv, x, mask):
    ok = np.isfinite(f) & (iv > 0); use = ok & ~mask
    for _ in range(3):
        if use.sum() < 60:
            return None
        p = np.polyfit(x[use], f[use], 2, w=np.sqrt(iv[use])); P = np.polyval(p, x)
        use = use & (np.abs((f - P) * np.sqrt(iv)) < 3)
    if np.median(P) <= 0:
        return None
    return np.where(ok, f / P, np.nan), np.where(ok, iv * P ** 2, 0.0)


def z(T, r, w):
    den = np.sqrt((T ** 2) @ w); den[den == 0] = np.inf
    return (T @ (w * r)) / den


def resid(n, w, t, mask):
    ok = (w > 0) & np.isfinite(n) & np.isfinite(t); use = ok & ~mask
    if use.sum() < 60:
        return None
    A = np.vstack([np.ones(use.sum()), t[use] - 1]).T * np.sqrt(w[use])[:, None]
    b, a = np.linalg.lstsq(A, n[use] * np.sqrt(w[use]), rcond=None)[0]
    r = np.where(ok, n - (b + a * (t - 1)), 0.0); ww = np.where(ok, w, 0.0)
    q = r[use] * np.sqrt(w[use]); s = max(1.0, 1.4826 * np.median(np.abs(q - np.median(q))))
    return r, ww / s ** 2, s, a


def snowwhite_all():
    return cas("SELECT sdss_id, gaia_dr3_source_id, ra, dec, classification, teff, logg, snr, n_boss_visits, g_mag, bp_mag, rp_mag, plx, e_plx "
               "FROM snow_white_boss_star WHERE classification IS NOT NULL AND classification <> ''", dtype=str)


def pass2(outf):
    sw = snowwhite_all(); have = {os.path.basename(p)[:-4] for p in glob.glob(os.path.join(STORE, "*", "*.npz"))}
    t = sw[sw.sdss_id.isin(have)].drop_duplicates("sdss_id").reset_index(drop=True)
    t["bprp"] = t.bp_mag.astype(float) - t.rp_mag.astype(float); plx = t.plx.astype(float)
    t["MG"] = np.where(plx > 0, t.g_mag.astype(float) + 5 * np.log10(np.where(plx > 0, plx, 1) / 100), np.nan)
    good = plx / t.e_plx.astype(float) > 3
    t = t[(good & (t.MG > 8.5)) | ~good].reset_index(drop=True); good = t.plx.astype(float) / t.e_plx.astype(float) > 3
    t["MGf"] = np.where(good, t.MG, 11.2 + 3.3 * t.bprp); t["plx_ok"] = good
    c = t.classification.fillna("")
    t["grp"] = np.where(c.str.contains("CV"), "CV", np.where(c.str.contains("MS"), "MS", np.where(c.str.startswith("DA"), "H", "nonH")))
    N = len(t); NC = np.full((N, len(WC)), np.nan, np.float32); WCw = np.zeros((N, len(WC)), np.float32); vis = [None] * N
    for i, sid in enumerate(t.sdss_id):
        d = np.load(os.path.join(STORE, sid[-2:], f"{sid}.npz"))
        a = norm(d["f_8250"].astype(float), d["iv_8250"].astype(float), XC, WIN)
        if a is not None:
            NC[i], WCw[i] = a
        vis[i] = (d["vf_8250"].astype(float), d["viv_8250"].astype(float), d["mjd"])
    snc = np.sqrt(np.nanmedian(np.where(WCw > 0, WCw, np.nan), axis=1)); t["snr_cat"] = np.round(snc, 1)
    feat = np.vstack([t.bprp / 0.04, t.MGf / 0.2]).T; rows = []
    for g in ("H", "nonH", "MS", "CV"):
        gi = np.where(t.grp == g)[0]; ref = gi[(snc[gi] > 10) & np.isfinite(feat[gi]).all(1)]
        if len(ref) < 20:
            ref = gi
        tree = cKDTree(feat[ref]); k = min(101, len(ref))
        for i in gi:
            _, nb = tree.query(feat[i], k=k); nb = ref[np.atleast_1d(nb)]; nb = nb[nb != i][:100]
            row = dict(sdss_id=t.sdss_id.iat[i], gaia_dr3=t.gaia_dr3_source_id.iat[i], snowwhite_class=t.classification.iat[i], group=g,
                       plx_ok=bool(t.plx_ok.iat[i]), G=round(float(t.g_mag.iat[i]), 2), bp_rp=round(float(t.bprp.iat[i]), 3),
                       M_G=round(float(t.MG.iat[i]), 2), snr_cat=t.snr_cat.iat[i])
            tc = np.nanmedian(np.where(WCw[nb] > 0, NC[nb], np.nan), axis=0)
            rr = resid(NC[i].astype(float), WCw[i].astype(float), tc, WIN) if np.isfinite(NC[i]).sum() > 100 else None
            if rr is not None:
                r, w, s, a = rr; best = (-1e9, None, 0)
                for n in TT:
                    zz = z(TT[n], r, w); kk = int(np.argmax(zz))
                    if zz[kk] > best[0]:
                        best = (float(zz[kk]), n, kk)
                zb, nbest, kb = best
                win = np.zeros(len(WC), bool)
                for lam in CAT:
                    win |= np.abs(WC / (lam * (1 + VEL[kb] / C)) - 1) * C < 700
                row.update(z_cat=round(zb, 2), template=nbest, v_kms=int(VEL[kb]), noise_scale=round(s, 2),
                           z_lines="/".join(f"{z(T[kb:kb + 1], r, w)[0]:.1f}" for T in TL[nbest]),
                           z_abs=round(float(max((-z(TT[n], r, w)).max() for n in TT)), 2),
                           z_fake=round(float(max(z(FAKE[d][n], r, w).max() for d in FAKE for n in TT)), 2),
                           ew_A=round(float(np.sum((r * np.gradient(WC))[win & (w > 0)])), 2))
                pv = []
                if vis[i][0].shape[0] > 1:
                    for vf, viv in zip(vis[i][0], vis[i][1]):
                        q = norm(vf, viv, XC, WIN); qq = None if q is None else resid(q[0], q[1], tc, WIN)
                        pv.append("nan" if qq is None else f"{z(TT[nbest][kb:kb + 1], qq[0], qq[1])[0]:.1f}")
                row["z_visits"] = "/".join(pv); row["visit_mjds"] = "/".join(str(m) for m in vis[i][2])
            rows.append(row)
    pd.DataFrame(rows).to_csv(outf, index=False); print(len(rows), "objects")


if __name__ == "__main__":
    if sys.argv[1] == "--sample":
        sw = snowwhite_all(); have = {os.path.basename(p)[:-4] for p in glob.glob(os.path.join(STORE, "*", "*.npz"))}
        todo = [s for s in sw.sdss_id.drop_duplicates() if s not in have]; print(len(todo), "objects to store")

        def one(sid):
            for k in range(3):
                try:
                    st = store(sid); break
                except Exception as e:
                    st = f"hole: {type(e).__name__}"; time.sleep(3)
            p = os.path.join(CACHE, f"mwmVisit-0.8.1-{sid}.fits")
            if os.path.exists(p):
                os.remove(p)
            return sid, st
        with ThreadPoolExecutor(4) as ex:
            res = list(ex.map(one, todo))
        bad = [r for r in res if r[1] != "ok"]; print(len(res) - len(bad), "stored;", len(bad), "not stored:", bad[:20])
    elif sys.argv[1] == "--pass2":
        pass2(sys.argv[2])
    elif sys.argv[1] == "--table":
        p = pd.read_csv(sys.argv[2], dtype={"sdss_id": str, "gaia_dr3": str})
        obj = pd.read_csv("../data/gas_disc_objects.csv", dtype={"sdss_id": str, "gaia_dr3": str})
        out = obj.merge(p.drop(columns=["gaia_dr3"]), on="sdss_id", how="left")
        out["v_kms"] = out.v_kms.astype("Int64"); out.to_csv("../tables/gas_disc_screen.csv", index=False); print(out.to_string())
