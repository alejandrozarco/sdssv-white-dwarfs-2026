"""Figures in ../figures from the tables in ../tables and public data (run after the measurement scripts)."""
import os, json, subprocess, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d, percentile_filter
from astropy.io import fits
from astropy.timeseries import LombScargle
from astroquery.mast import Observations
from sdssv import visits, coadd, C
T = "../tables"; F = "../figures"; os.makedirs(F, exist_ok=True)


def norm(w, f, iv, width=120, pct=85):
    m = np.isfinite(f) & (iv > 0); w, f = w[m], f[m]
    n = int(round(np.log10(1 + width / 5000) / np.median(np.diff(np.log10(w)))))
    return w, f / gaussian_filter1d(percentile_filter(f, pct, size=n), n / 4)


def spectrum_for(row):
    s = str(row["spectrum"])
    if s == "coadd":
        return coadd([x for x in visits(row["sdss_id"]) if x["in_stack"]])
    v = [x for x in visits(row["sdss_id"]) if str(x["mjd"]) == s.split()[-1]][0]
    return v["wave"], v["flux"], v["ivar"]


def zeeman():
    t = pd.read_csv(f"{T}/magnetic_zeeman.csv", dtype={"gaia_dr3": str, "sdss_id": str})
    fig, ax = plt.subplots(len(t), 2, figsize=(11, 1.25 * len(t)), sharex="col")
    for i, (_, r) in enumerate(t.iterrows()):
        w, f, iv = spectrum_for(r); m = (iv > 0) & np.isfinite(f)
        for j, (ln, l0, lo, hi) in enumerate((("Ha", 6564.61, 6300, 6830), ("Hb", 4862.68, 4660, 5060))):
            s = m & (w > lo) & (w < hi); y = gaussian_filter1d(f[s], 1.5); y = y / np.median(y)
            ax[i, j].plot(w[s], y, "k", lw=0.6); ax[i, j].axvline(l0, color="0.6", lw=0.6, ls=":")
            for x in (r[f"{ln}_sigma_minus_A"], r[f"{ln}_pi_A"], r[f"{ln}_sigma_plus_A"]):
                ax[i, j].axvline(x, color="r", lw=0.8)
            ax[i, j].set_yticks([])
            if j == 0:
                ax[i, j].text(0.01, 0.08, f"{r.gaia_dr3}  S/N {r.snr}", transform=ax[i, j].transAxes, fontsize=6)
    ax[-1, 0].set_xlabel("wavelength (A)"); ax[-1, 1].set_xlabel("wavelength (A)")
    ax[0, 0].set_title("H-alpha; red: fitted component centres; grey dotted: 6564.61 A", fontsize=8); ax[0, 1].set_title("H-beta; red: fitted component centres; grey dotted: 4862.68 A", fontsize=8)
    plt.tight_layout(h_pad=0.1); plt.savefig(f"{F}/magnetic_zeeman_halpha_hbeta.png", dpi=90); plt.close()


def carbon_optical():
    L = json.load(open("../data/nist_vacuum_lines.json"))
    rows = [("5208047381438507520", "95077848"), ("6466745168812781568", "110600288"), ("5836110898905253760", "102600838")]
    fig, ax = plt.subplots(len(rows), 1, figsize=(15, 3.2 * len(rows)))
    for a, (g, sid) in zip(ax, rows):
        d = json.load(open(f"carbon_lines_{g}.json")); v = d["reference_velocity_for_line_depths"]
        sp, col = ("C I", "orange") if d["ccf_coadd"]["C I"]["contrast"] > d["ccf_coadd"]["C II"]["contrast"] else ("C II", "m")
        w, f = norm(*coadd(visits(sid))); s = (w > 3800) & (w < 9200)
        a.plot(w[s], gaussian_filter1d(f[s], 1.0), "k", lw=0.6)
        for lam, I in L[sp]:
            if 3800 < lam < 9200 and I >= (100 if sp == "C II" else 1e5):
                a.axvline(lam * (1 + v / C), color=col, lw=0.5, alpha=0.5)
        for lam in (6564.632, 4862.691):
            a.axvline(lam * (1 + v / C), color="b", ls="--", lw=0.7)
        a.set_xlim(3800, 9200); a.set_ylim(0.45, 1.2)
        a.set_title(f"Gaia DR3 {g} (sdss_id {sid}); {sp} NIST lines ({col}) and H-alpha/H-beta (blue) at {v:+d} km/s", fontsize=8)
    ax[-1].set_xlabel("vacuum wavelength (A)")
    plt.tight_layout(); plt.savefig(f"{F}/carbon_optical_spectra.png", dpi=90); plt.close()


def carbon_cos():
    p = Observations.get_product_list(Observations.query_criteria(obs_id="lfac0z010")); p = p[[x == "lfac0z010_x1dsum.fits" for x in p["productFilename"]]]
    d = fits.open(Observations.download_products(p, download_dir="../data/cache")["Local Path"][0])[1].data
    W = np.concatenate([r["WAVELENGTH"] for r in d]); Fl = np.concatenate([r["FLUX"] for r in d]); Q = np.concatenate([r["DQ_WGT"] for r in d]); o = np.argsort(W)
    W, Fl, Q = W[o], Fl[o], Q[o]; W, Fl = W[Q > 0], Fl[Q > 0]
    lines = {"C III": [1174.93, 1175.26, 1175.59, 1175.71, 1175.99, 1176.37, 1247.38], "C II": [1323.95, 1334.53, 1335.71], "Si II": [1260.42, 1264.74],
             "Si IV": [1393.76, 1402.77], "Ly-alpha": [1215.67], "O I (airglow)": [1302.17, 1304.86, 1306.03]}
    cols = {"C III": "m", "C II": "r", "Si II": "g", "Si IV": "g", "Ly-alpha": "b", "O I (airglow)": "0.5"}
    fig, ax = plt.subplots(2, 1, figsize=(15, 7))
    for a, (lo, hi) in zip(ax, [(1130, 1285), (1285, 1432)]):
        s = (W > lo) & (W < hi); a.plot(W[s], np.convolve(Fl[s], np.ones(7) / 7, mode="same"), "k", lw=0.6)
        for nm, ls in lines.items():
            for l in ls:
                if lo < l < hi:
                    a.axvline(l, color=cols[nm], lw=0.7, alpha=0.7)
        a.set_xlim(lo, hi); a.set_ylim(0, np.percentile(Fl[s], 99.5) * 1.1); a.set_ylabel("flux (erg/s/cm2/A)")
    ax[0].set_title("HST/COS G130M lfac0z010, Gaia DR3 5208047381438507520; magenta C III, red C II, green Si II/IV, blue Ly-alpha, grey O I airglow", fontsize=8)
    ax[1].set_xlabel("wavelength (A)"); plt.tight_layout(); plt.savefig(f"{F}/carbon_5208047381438507520_cos.png", dpi=90); plt.close()


def galex():
    d = pd.read_csv("galex_colours_5208047381438507520.csv", dtype={"source_id": str}); fig, a = plt.subplots(figsize=(7, 5))
    for gname, col, mk in (("DA", "0.6", "."), ("DB", "tab:green", "s"), ("DQ", "tab:purple", "D")):
        s = d[(d.group == gname) & d.fuv_nuv.notna() & (d.e_FUVmag < 0.1)]; a.scatter(s.bp_rp, s.fuv_nuv, s=10 if gname == "DA" else 22, c=col, marker=mk, label=f"{gname} ({len(s)})")
    t = d[d.group == "target"].iloc[0]; a.scatter([t.bp_rp], [t.fuv_nuv], s=150, c="r", marker="*", label="Gaia DR3 5208047381438507520")
    a.set_xlabel("Gaia BP-RP"); a.set_ylabel("GALEX FUV-NUV"); a.legend(fontsize=7)
    a.set_title("GUVcat AIS; SDSS-V SnowWhite DA/DB spectra and MWDD DQ types (Teff >= 15 kK)", fontsize=8)
    plt.tight_layout(); plt.savefig(f"{F}/galex_fuv_nuv.png", dpi=90); plt.close()


def eclipse():
    ns = {}
    exec(open("j0353_eclipse.py").read().split("Tref = 2460670.38343")[0], ns)
    D = ns["D"]; P = 0.14786971; T0 = 2460670.38336
    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for a, b in zip(ax, ("o", "c")):
        ph = ((D[b]["t"] - T0) / P + 0.5) % 1 - 0.5; a.errorbar(ph, D[b]["f"], D[b]["e"], fmt=".", ms=2, color="0.6", elinewidth=0.3, alpha=0.5)
        edges = np.linspace(-0.5, 0.5, 101); c = 0.5 * (edges[1:] + edges[:-1]); idx = np.digitize(ph, edges) - 1
        mb = [np.average(D[b]["f"][idx == i], weights=1 / D[b]["e"][idx == i] ** 2) if np.any(idx == i) else np.nan for i in range(100)]
        a.plot(c, mb, "k-", lw=1); a.set_ylabel(f"ATLAS {b} flux (uJy, difference)")
    ax[1].set_xlabel(f"phase (P = {P} d, T0 = BJD_TDB {T0})"); ax[0].set_title("Gaia DR3 4731701084150029824", fontsize=9)
    plt.tight_layout(); plt.savefig(f"{F}/eclipsing_4731701084150029824_atlas_phase.png", dpi=90); plt.close()


def balmer():
    rows = [("2002597083798483200", "65701864"), ("1977447164064222976", "65312747")]
    fig, ax = plt.subplots(len(rows), 2, figsize=(11, 3.2 * len(rows)))
    for i, (g, sid) in enumerate(rows):
        w, f, iv = coadd(visits(sid))
        for j, (l0, lo, hi) in enumerate(((6564.61, 6450, 6680), (4862.68, 4760, 4970))):
            s = (w > lo) & (w < hi) & np.isfinite(f); ax[i, j].plot((w[s] / l0 - 1) * C, gaussian_filter1d(f[s], 1.0), "k", lw=0.7)
            ax[i, j].set_title(f"Gaia DR3 {g} {'H-alpha' if j == 0 else 'H-beta'}", fontsize=8)
        ax[i, 0].set_ylabel("flux")
    ax[-1, 0].set_xlabel("velocity (km/s)"); ax[-1, 1].set_xlabel("velocity (km/s)")
    plt.tight_layout(); plt.savefig(f"{F}/balmer_emission.png", dpi=90); plt.close()


def tess_amplitude(tic, sec, cad, lo, hi, ax, label):
    obs = Observations.query_criteria(obs_collection="TESS", target_name=str(tic), sequence_number=int(sec), dataproduct_type="timeseries", provenance_name="SPOC")
    fn = [n for n in Observations.get_product_list(obs)["productFilename"] if n.endswith("fast-lc.fits" if int(cad) == 20 else "s_lc.fits")][0]
    path = f"../data/cache/{fn}"
    if not os.path.exists(path):
        subprocess.run(["curl", "-sL", "-m", "600", "-o", path, f"https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS/product/{fn}"], check=True)
    d = fits.open(path)[1].data; m = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]); t = d["TIME"][m]; y = d["PDCSAP_FLUX"][m] / np.median(d["PDCSAP_FLUX"][m]) - 1
    mad = 1.4826 * np.median(np.abs(y - np.median(y))); ok = np.abs(y - np.median(y)) < 5 * mad; t, y = t[ok], y[ok]
    fr = np.arange(lo, hi, 0.001); amp = np.sqrt(4 * LombScargle(t, y, normalization="psd").power(fr, method="fast") / len(t)) * 1e3
    ax.plot(fr, amp, "k", lw=0.5); ax.set_title(label, fontsize=8); ax.set_ylabel("amplitude (ppt)")


def zz():
    z = pd.read_csv(f"{T}/zz_ceti_tess_sectors.csv", dtype={"gaia_dr3": str}); z = z[z.pixel_chi2_target.notna()]
    fig, ax = plt.subplots(len(z), 1, figsize=(10, 2.2 * len(z)))
    for a, (_, r) in zip(ax, z.iterrows()):
        tess_amplitude(r.tic, r.sector, r.cadence_s, 20, 359 if r.cadence_s == 120 else 400, a, f"Gaia DR3 {r.gaia_dr3}  TIC {r.tic}  S{r.sector} {r.cadence_s}-s")
        a.axvline(r.pixel_freq_cd, color="r", lw=0.6, alpha=0.5)
    ax[-1].set_xlabel("frequency (c/d)"); plt.tight_layout(); plt.savefig(f"{F}/zz_ceti_tess_amplitude_spectra.png", dpi=80); plt.close()


if __name__ == "__main__":
    for fn in (zeeman, carbon_optical, carbon_cos, galex, eclipse, balmer, zz):
        fn(); print(fn.__name__, "done", flush=True)
