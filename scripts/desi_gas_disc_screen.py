"""Ca II triplet emission screen of DESI DR1 white-dwarf spectra, with the method of gas_disc_screen.py (pass 2).
Usage: python desi_gas_disc_screen.py --sample              store the spectra of all TARGETIDs in the Amorim et al. (2026) DESI DR1
                                                            class file (SPARCL retrieve_by_specid, DESI-DR1, batches of 200)
       python desi_gas_disc_screen.py --pass2 <out.csv>     template subtraction and matched filter for every stored spectrum
       python desi_gas_disc_screen.py --table <pass2.csv>   writes ../tables/desi_gas_disc_screen.csv for ../data/desi_gas_disc_objects.csv
Storage: flux and ivar (mask != 0 -> ivar 0) interpolated onto the SDSS-V screen grid (3850-9250 A, 6e-5 dex; ivar scaled by the
pixel-size ratio) over 3880-4000, 6400-6750, 7740-7810 and 8250-8950 A, in ../data/cache/desi_gas_store/.
Pass 2: as gas_disc_screen.py, with neighbours in (BP-RP/0.04, M_G/0.2) from the class file; groups from its CLASS: H (starting DA,
no '+'), nonH (other classes without '+', not CV/AMCVn), MS (contains '+'), CV. Objects without M_G use 11.2 + 3.3 (BP-RP).
The table adds the DESI DR1 spectral class of Swan et al. (2026), matched by WDJ name."""
import sys, os, glob, time, subprocess, numpy as np, pandas as pd, warnings
from scipy.spatial import cKDTree
from sdssv import CACHE
from gas_disc_screen import norm, z, resid, TT, TL, FAKE, WC, XC, WIN, CAT, VEL, C, GRID
warnings.filterwarnings("ignore")
STORE = os.path.join(CACHE, "desi_gas_store")
REG = [(3880, 4000), (6400, 6750), (7740, 7810), (8250, 8950)]
AMORIM = "https://raw.githubusercontent.com/AstroWeljr/DESI-DR1-Class/main/DESI_CLASS_FINAL.txt"
SWAN = "https://cygnus.astro.warwick.ac.uk/phsdaj/WhiteDwarfs/Tables/DESI_DR1_catalogue_1.0.csv.gz"


def cached(url):
    path = os.path.join(CACHE, os.path.basename(url)); os.makedirs(CACHE, exist_ok=True)
    if not os.path.exists(path):
        subprocess.run(["curl", "-sL", "-m", "900", "-o", path, url], check=True)
    return path


def classes():
    return pd.read_csv(cached(AMORIM), sep=r"\s+", dtype={"DESIID": str, "edr3id": str})


def save(rec):
    t = str(rec["specid"]); w = np.array(rec["wavelength"], float); f = np.array(rec["flux"], float)
    iv = np.array(rec["ivar"], float) * (np.array(rec["mask"]) == 0); ok = (iv > 0) & np.isfinite(f)
    if ok.sum() < 1000:
        return "bad spectrum"
    d = {}
    for lo, hi in REG:
        m = (GRID > lo) & (GRID < hi)
        d[f"f_{lo}"] = np.interp(GRID[m], w[ok], f[ok], left=np.nan, right=np.nan).astype(np.float32)
        d[f"iv_{lo}"] = np.interp(GRID[m], w[ok], iv[ok], left=0, right=0).astype(np.float32) * (np.median(np.diff(GRID[m])) / np.median(np.diff(w)))
    os.makedirs(os.path.join(STORE, t[-2:]), exist_ok=True); np.savez_compressed(os.path.join(STORE, t[-2:], f"{t}.npz"), **d)
    return "ok"


def sample():
    from sparcl.client import SparclClient
    c = SparclClient(read_timeout=900); a = classes()
    have = {os.path.basename(p)[:-4] for p in glob.glob(os.path.join(STORE, "*", "*.npz"))}
    todo = [i for i in a.DESIID if i not in have]; stat = {}; holes = []
    for i in range(0, len(todo), 200):
        chunk = todo[i:i + 200]
        for k in range(3):
            try:
                r = c.retrieve_by_specid([int(x) for x in chunk], include=["specid", "flux", "ivar", "wavelength", "mask"], dataset_list=["DESI-DR1"], limit=2000); break
            except Exception:
                r = None; time.sleep(30 * (k + 1))
        if r is None:
            holes.append(i); continue
        got = {str(x["specid"]): x for x in r.records}
        for t in chunk:
            s = save(got[t]) if t in got else "missing"; stat[s] = stat.get(s, 0) + 1
    print(stat, "failed batches (start index):", holes)


def pass2(outf):
    a = classes(); have = {os.path.basename(p)[:-4] for p in glob.glob(os.path.join(STORE, "*", "*.npz"))}
    t = a[a.DESIID.isin(have)].drop_duplicates("DESIID").reset_index(drop=True)
    t["bprp"] = pd.to_numeric(t.bp_rp, errors="coerce"); t["MGv"] = pd.to_numeric(t.MG, errors="coerce")
    t["MGf"] = np.where(np.isfinite(t.MGv), t.MGv, 11.2 + 3.3 * t.bprp); cl = t.CLASS.astype(str)
    t["grp"] = np.where(cl.str.contains("CV|AMCVn"), "CV", np.where(cl.str.contains(r"\+"), "MS", np.where(cl.str.startswith("DA"), "H", "nonH")))
    N = len(t); NC = np.full((N, len(WC)), np.nan, np.float32); WCw = np.zeros((N, len(WC)), np.float32)
    for i, tid in enumerate(t.DESIID):
        d = np.load(os.path.join(STORE, tid[-2:], f"{tid}.npz")); q = norm(d["f_8250"].astype(float), d["iv_8250"].astype(float), XC, WIN)
        if q is not None:
            NC[i], WCw[i] = q
    snc = np.sqrt(np.nanmedian(np.where(WCw > 0, WCw, np.nan), axis=1)); t["snr_cat"] = np.round(snc, 1)
    feat = np.vstack([t.bprp.fillna(0) / 0.04, t.MGf.fillna(12) / 0.2]).T; rows = []
    for g in ("H", "nonH", "MS", "CV"):
        gi = np.where(t.grp == g)[0]; ref = gi[(snc[gi] > 10) & np.isfinite(feat[gi]).all(1)]
        if len(ref) < 20:
            ref = gi
        tree = cKDTree(feat[ref]); k = min(101, len(ref))
        for i in gi:
            _, nb = tree.query(feat[i], k=k); nb = ref[np.atleast_1d(nb)]; nb = nb[nb != i][:100]
            row = dict(targetid=t.DESIID.iat[i], wdj_name=t["#Name"].iat[i], ra=t["RA(deg)"].iat[i], dec=t["DEC(deg)"].iat[i], amorim_class=t.CLASS.iat[i],
                       group=g, G=t.G.iat[i], bp_rp=t.bp_rp.iat[i], M_G=t.MG.iat[i], snr_cat=t.snr_cat.iat[i])
            tc = np.nanmedian(np.where(WCw[nb] > 0, NC[nb], np.nan), axis=0)
            rr = resid(NC[i].astype(float), WCw[i].astype(float), tc, WIN) if np.isfinite(NC[i]).sum() > 100 else None
            if rr is not None:
                r, w, s, aa = rr; best = (-1e9, None, 0)
                for n in TT:
                    zz = z(TT[n], r, w); kk = int(np.argmax(zz))
                    if zz[kk] > best[0]:
                        best = (float(zz[kk]), n, kk)
                zb, nbest, kb = best; win = np.zeros(len(WC), bool)
                for lam in CAT:
                    win |= np.abs(WC / (lam * (1 + VEL[kb] / C)) - 1) * C < 700
                row.update(z_cat=round(zb, 2), template=nbest, v_kms=int(VEL[kb]), noise_scale=round(s, 2),
                           z_lines="/".join(f"{z(T[kb:kb + 1], r, w)[0]:.1f}" for T in TL[nbest]),
                           z_abs=round(float(max((-z(TT[n], r, w)).max() for n in TT)), 2),
                           z_fake=round(float(max(z(FAKE[d][n], r, w).max() for d in FAKE for n in TT)), 2),
                           ew_A=round(float(np.sum((r * np.gradient(WC))[win & (w > 0)])), 2))
            rows.append(row)
    pd.DataFrame(rows).to_csv(outf, index=False); print(len(rows), "spectra")


if __name__ == "__main__":
    if sys.argv[1] == "--sample":
        sample()
    elif sys.argv[1] == "--pass2":
        pass2(sys.argv[2])
    elif sys.argv[1] == "--table":
        p = pd.read_csv(sys.argv[2], dtype={"targetid": str})
        obj = pd.read_csv("../data/desi_gas_disc_objects.csv", dtype={"targetid": str, "gaia_dr3": str})
        sw = pd.read_csv(cached(SWAN), usecols=["WDJname", "specType"]).rename(columns={"WDJname": "wdj_name", "specType": "swan_class"})
        out = obj.merge(p, on="targetid", how="left").merge(sw.drop_duplicates("wdj_name"), on="wdj_name", how="left")
        out["v_kms"] = out.v_kms.astype("Int64"); out.to_csv("../tables/desi_gas_disc_screen.csv", index=False); print(out.to_string())
