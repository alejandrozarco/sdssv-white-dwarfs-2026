"""Figures in ../figures from the tables in ../tables and public data (run after the measurement scripts)."""
import os, json, subprocess, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d, percentile_filter
from astropy.io import fits
from astropy.timeseries import LombScargle
from astroquery.mast import Observations
from sdssv import visits, coadd, C
T = "../tables"; F = "../figures"


def out(topic, name):
    os.makedirs(f"{F}/{topic}", exist_ok=True)
    return f"{F}/{topic}/{name}"


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


def zeeman_panel(ax, r, w, f, iv, ln, l0, lo, hi, label=True):
    m = (iv > 0) & np.isfinite(f) & (w > lo) & (w < hi); y = gaussian_filter1d(f[m], 1.5); y = y / np.median(y)
    ax.plot(w[m], y, "k", lw=0.7); ax.axvline(l0, color="0.6", lw=0.6, ls=":")
    for x in (r[f"{ln}_sigma_minus_A"], r[f"{ln}_pi_A"], r[f"{ln}_sigma_plus_A"]):
        ax.axvline(x, color="r", lw=0.8)
    ax.set_yticks([])


def zeeman():
    """One figure per star (H-alpha and H-beta with the fitted component centres) and an overview of the six highest-S/N stars."""
    t = pd.read_csv(f"{T}/magnetic_zeeman.csv", dtype={"gaia_dr3": str, "sdss_id": str}); spec = {}
    for _, r in t.iterrows():
        w, f, iv = spec.setdefault(r.gaia_dr3, spectrum_for(r))
        fig, ax = plt.subplots(1, 2, figsize=(10, 2.4))
        for a, (ln, l0, lo, hi, nm) in zip(ax, (("Ha", 6564.61, 6300, 6830, "H-alpha"), ("Hb", 4862.68, 4660, 5060, "H-beta"))):
            zeeman_panel(a, r, w, f, iv, ln, l0, lo, hi); a.set_xlabel("wavelength (A)")
            a.set_title(f"{nm}: B_split {r[f'B_split_{ln}_MG']} MG", fontsize=8)
        fig.suptitle(f"Gaia DR3 {r.gaia_dr3} ({r['name']}), SDSS-V {r.spectrum}, S/N {r.snr}; red: fitted component centres", fontsize=8)
        plt.tight_layout(); plt.savefig(out("zeeman", f"{r.gaia_dr3}.png"), dpi=90); plt.close()
    top = t.astype({"snr": float}).sort_values("snr", ascending=False).head(6)
    fig, ax = plt.subplots(3, 2, figsize=(10, 6.5))
    for a, (_, r) in zip(ax.ravel(), top.iterrows()):
        w, f, iv = spec[r.gaia_dr3]; zeeman_panel(a, r, w, f, iv, "Ha", 6564.61, 6300, 6830)
        a.set_title(f"Gaia DR3 {r.gaia_dr3}: B_split {r.B_split_Ha_MG} MG", fontsize=8)
    for a in ax[-1]:
        a.set_xlabel("wavelength (A)")
    fig.suptitle("H-alpha of the six highest-S/N stars in magnetic_zeeman.csv; red: fitted component centres; grey dotted: 6564.61 A", fontsize=8)
    plt.tight_layout(); plt.savefig(out("zeeman", "overview.png"), dpi=90); plt.close()


def carbon_optical():
    L = json.load(open("../data/nist_vacuum_lines.json"))
    rows = [("5208047381438507520", "95077848"), ("6466745168812781568", "110600288"), ("5836110898905253760", "102600838")]
    for g, sid in rows:
        fig, a = plt.subplots(figsize=(12, 3.2))
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
        a.set_xlabel("vacuum wavelength (A)"); plt.tight_layout(); plt.savefig(out("carbon", f"{g}_lines.png"), dpi=90); plt.close()


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
    ax[1].set_xlabel("wavelength (A)"); plt.tight_layout(); plt.savefig(out("carbon", "5208047381438507520_cos.png"), dpi=90); plt.close()



def carbon_screen_spectra():
    import carbon_screen as cs
    t = pd.read_csv(f"{T}/carbon_screen.csv", dtype={"gaia_dr3": str, "sdss_id": str}).set_index("gaia_dr3")
    rows = ["883885440381808000", "4847399905305694080", "2076678981825545088", "6465542891501713408", "343958710690034944"]
    for g in rows:
        fig, a = plt.subplots(figsize=(12, 3.2))
        r = t.loc[g]; v = int(r.C_v_kms)
        w, f = norm(*coadd(visits(r.sdss_id))); s = (w > 3800) & (w < 9200)
        a.plot(w[s], gaussian_filter1d(f[s], 1.0), "k", lw=0.6)
        for sp, col in (("C I", "orange"), ("C II", "m")):
            for lam in cs.LINES[sp]:
                a.axvline(lam * (1 + v / C), color=col, lw=0.6, alpha=0.7)
        for lam in (6564.632, 4862.691, 4341.69, 4102.89):
            a.axvline(lam * (1 + v / C), color="b", ls="--", lw=0.7)
        a.set_xlim(3800, 9200); a.set_ylim(0.4, 1.25)
        a.set_title(f"Gaia DR3 {g} (sdss_id {r.sdss_id}, S/N {r.snr_max}); C I (orange), C II (magenta) and Balmer (blue) at {v:+d} km/s; "
                    f"screen contrast C {r.C_contrast}, C I {r.CI_contrast}, C II {r.CII_contrast}", fontsize=8)
        a.set_xlabel("vacuum wavelength (A)"); plt.tight_layout(); plt.savefig(out("carbon", f"{g}_sdssv.png"), dpi=90); plt.close()


def galex():
    d = pd.read_csv("galex_colours_5208047381438507520.csv", dtype={"source_id": str}); fig, a = plt.subplots(figsize=(7, 5))
    for gname, col, mk in (("DA", "0.6", "."), ("DB", "tab:green", "s"), ("DQ", "tab:purple", "D")):
        s = d[(d.group == gname) & d.fuv_nuv.notna() & (d.e_FUVmag < 0.1)]; a.scatter(s.bp_rp, s.fuv_nuv, s=10 if gname == "DA" else 22, c=col, marker=mk, label=f"{gname} ({len(s)})")
    t = d[d.group == "target"].iloc[0]; a.scatter([t.bp_rp], [t.fuv_nuv], s=150, c="r", marker="*", label="Gaia DR3 5208047381438507520")
    a.set_xlabel("Gaia BP-RP"); a.set_ylabel("GALEX FUV-NUV"); a.legend(fontsize=7)
    a.set_title("GUVcat AIS; SDSS-V SnowWhite DA/DB spectra and MWDD DQ types (Teff >= 15 kK)", fontsize=8)
    plt.tight_layout(); plt.savefig(out("carbon", "galex_fuv_nuv.png"), dpi=90); plt.close()


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
    plt.tight_layout(); plt.savefig(out("eclipse", "4731701084150029824_atlas_phase.png"), dpi=90); plt.close()


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
    plt.tight_layout(); plt.savefig(out("balmer_emission", "balmer_emission.png"), dpi=90); plt.close()


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
    for _, r in z.iterrows():
        fig, a = plt.subplots(figsize=(9, 2.6))
        tess_amplitude(r.tic, r.sector, r.cadence_s, 20, 359 if r.cadence_s == 120 else 400, a, f"Gaia DR3 {r.gaia_dr3}  TIC {r.tic}  S{r.sector} {r.cadence_s}-s; red: pixel-fit frequency")
        a.axvline(r.pixel_freq_cd, color="r", lw=0.6, alpha=0.5); a.set_xlabel("frequency (c/d)")
        plt.tight_layout(); plt.savefig(out("zz_ceti", f"{r.gaia_dr3}_s{r.sector}.png"), dpi=90); plt.close()


def periodic():
    import periodic_6021870154194477312 as P6
    g = P6.GID; t = pd.read_csv(f"{T}/periodic_{g}.csv", dtype={"gaia_dr3": str}); f = float(t.frequency_cd.iloc[0])
    t0 = float(t[(t.gaia_dr3 == g) & (t.dataset == "ATLAS") & (t.band == "o")].t_max_bjd.iloc[0])
    at = P6.load_atlas(g); _, (fr, pw) = P6.adopted_frequency(at); ga = P6.load_gaia(); tt, crowd = P6.load_tess()
    fig, a = plt.subplots(figsize=(11, 3.2)); w, fl, iv = coadd(visits("106653583")); s = np.isfinite(fl) & (iv > 0) & (w > 3750) & (w < 9000)
    a.plot(w[s], gaussian_filter1d(fl[s], 2), "k", lw=0.6)
    for l in (6564.61, 4862.68, 4341.69, 4102.89):
        a.axvline(l, color="0.6", lw=0.6, ls=":")
    a.set_xlabel("wavelength (A)"); a.set_ylabel("flux"); a.set_title(f"Gaia DR3 {g}: SDSS-V spectrum (sdss_id 106653583, all visits); grey dotted: H-alpha to H-delta", fontsize=8)
    plt.tight_layout(); plt.savefig(out("periodic", f"{g}_spectrum.png"), dpi=90); plt.close()
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.6)); a = ax[0]; a.plot(fr, pw, "k", lw=0.4); a.set_xlabel("frequency (c/d)"); a.set_ylabel("GLS power")
    a.set_title(f"ATLAS c + o, fractional flux; adopted f = {f:.7f} c/d (P = {1440 / f:.3f} min)", fontsize=8)
    ph = lambda x: ((x - t0) * f) % 1
    a = ax[1]
    for b, col in (("c", "c"), ("o", "orange")):
        p = ph(at[b]["t"]); y = at[b]["f"] / P6.REF[b]; e = at[b]["e"] / P6.REF[b]; edges = np.linspace(0, 1, 21); c = (edges[1:] + edges[:-1]) / 2
        m = [np.sum(y[(p >= lo) & (p < hi)] / e[(p >= lo) & (p < hi)] ** 2) / np.sum(1 / e[(p >= lo) & (p < hi)] ** 2) for lo, hi in zip(edges[:-1], edges[1:])]
        se = [1 / np.sqrt(np.sum(1 / e[(p >= lo) & (p < hi)] ** 2)) for lo, hi in zip(edges[:-1], edges[1:])]
        for k in (0, 1):
            a.errorbar(c + k, m, se, fmt="o", ms=3, color=col, label=f"ATLAS {b}" if k == 0 else None)
    a.set_xlabel(f"phase (t_max = BJD_TDB {t0:.5f})"); a.set_ylabel("fractional flux (20 bins)"); a.legend(fontsize=7)
    a = ax[2]; d = ga["G"]; p = ph(d["t"])
    for k in (0, 1):
        a.errorbar(p + k, d["f"], d["e"], fmt="o", ms=3, color="k", label="Gaia DR3 G" if k == 0 else None)
    p = ph(tt["t"]); edges = np.linspace(0, 1, 21); c = (edges[1:] + edges[:-1]) / 2
    m = [np.mean(tt["f"][(p >= lo) & (p < hi)]) for lo, hi in zip(edges[:-1], edges[1:])]
    for k in (0, 1):
        a.plot(c + k, m, "s", ms=3, color="r", label=f"TESS S{P6.SECTOR} PDCSAP, CROWDSAP {crowd:.3f} (20 bins)" if k == 0 else None)
    a.set_xlabel("phase"); a.set_ylabel("fractional flux"); a.legend(fontsize=7)
    plt.tight_layout(); plt.savefig(out("periodic", f"{g}_lightcurve.png"), dpi=90); plt.close()


def periodic_white_dwarfs():
    import periodic_white_dwarfs as PW
    t = pd.read_csv(f"{T}/periodic_white_dwarfs.csv", dtype={"gaia_dr3": str}); ids = list(PW.SRC.index)
    for gid in ids:
        fig, ax = plt.subplots(1, 2, figsize=(11, 3.1))
        s = PW.SRC.loc[gid]; r = t[(t.gaia_dr3 == gid) & (t.dataset == s.ground)].iloc[0]; f = float(r.frequency_cd); t0 = float(r.t_max_bjd)
        x, y, e, g = PW.load_atlas(gid) if s.ground == "ATLAS" else PW.load_ztf(gid); _, (fr, pw) = PW.adopted_frequency(x, y, e, g)
        ax[0].plot(fr, pw, "k", lw=0.4); ax[0].set_xscale("log"); ax[0].set_ylabel("GLS power"); ax[0].set_xlabel("frequency (c/d)")
        ax[0].set_title(f"Gaia DR3 {gid} ({s['name']}): {s.ground}, f = {f:.6f} c/d, P = {24 / f:.4f} h", fontsize=8)
        p = ((x - t0) * f) % 1; edges = np.linspace(0, 1, 21); c = (edges[1:] + edges[:-1]) / 2
        m = [np.sum(y[(p >= lo) & (p < hi)] / e[(p >= lo) & (p < hi)] ** 2) / np.sum(1 / e[(p >= lo) & (p < hi)] ** 2) for lo, hi in zip(edges[:-1], edges[1:])]
        se = [1 / np.sqrt(np.sum(1 / e[(p >= lo) & (p < hi)] ** 2)) for lo, hi in zip(edges[:-1], edges[1:])]
        tg, yg, eg = PW.load_gaia(gid); pg = ((tg - t0) * f) % 1
        for k in (0, 1):
            ax[1].errorbar(c + k, m, se, fmt="o", ms=3, color="C0", label=f"{s.ground} (20 bins)" if k == 0 else None)
            ax[1].errorbar(pg + k, yg, eg, fmt=".", ms=3, color="0.4", alpha=0.7, label="Gaia DR3 G" if k == 0 else None)
        ax[1].set_ylabel("fractional flux"); ax[1].legend(fontsize=7); ax[1].set_xlabel("phase (0 = t_max of the ground-based fit)")
        plt.tight_layout(); plt.savefig(out("periodic", f"{gid}.png"), dpi=90); plt.close()


def hot_dq_comparison(extra=False):
    """SDSS-V spectrum of Gaia DR3 5208047381438507520 between an SDSS-V DA of similar colour and the SDSS DR17 spectrum of the hot DQ
    SDSS J234843.30-094245.3 (Dufour et al. 2008). Each spectrum is smoothed (Gaussian, 1.5 pixels) and scaled to its median flux at 4500-4600 A."""
    path = os.path.join("..", "data", "cache", "spec-7166-56602-0536.fits")
    if not os.path.exists(path):
        subprocess.run(["curl", "-sL", "-m", "300", "-o", path, "https://data.sdss.org/sas/dr17/eboss/spectro/redux/v5_13_2/spectra/lite/7166/spec-7166-56602-0536.fits"], check=True)
    d = fits.open(path)[1].data; w_dq, f_dq = 10 ** d["loglam"], d["flux"]
    specs = [("DA white dwarf Gaia DR3 2293913930823813888 (SDSS-V; BP-RP -0.42)", *coadd([v for v in visits("69198817") if v["in_stack"]])[:2]),
             ("Gaia DR3 5208047381438507520 (SDSS-V; BP-RP -0.41)", *coadd([v for v in visits("95077848") if v["in_stack"]])[:2]),
             ("hot DQ SDSS J234843.30-094245.3 (SDSS DR17; Dufour et al. 2008)", w_dq, f_dq)]
    if extra:
        specs.insert(2, ("Gaia DR3 6886051830805052288 (SDSS-V; BP-RP -0.40)", *coadd(visits("114554634"))[:2]))
    top = 1.3 * (len(specs) - 1) + 2.3
    fig, ax = plt.subplots(figsize=(11, 6.5))
    for k, (lab, w, f) in enumerate(specs):
        s = np.isfinite(f) & (w > 3820) & (w < 5000); y = gaussian_filter1d(f[s], 1.5); y = y / np.median(y[(w[s] > 4500) & (w[s] < 4600)])
        n = len(specs) - 1; col = ("0.35", "k", "k", "C3")[k] if extra else ("0.35", "k", "C3")[k]
        off = 1.3 * (n - k); ax.plot(w[s], y + off, color=col, lw=0.9)
        ax.text(4995, off + 1.32, lab, fontsize=10, ha="right", va="bottom", color=col, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    for l in (3920.7, 4075.9, 4267.3, 4372.5, 4619.2):
        ax.axvline(l, color="C1", lw=0.8, alpha=0.6); ax.text(l + 2, top - 0.45, "C II", fontsize=8, color="C1", rotation=90, va="bottom")
    for l, n in ((4862.68, "Hb"), (4341.69, "Hg"), (4102.89, "Hd"), (3971.2, "He")):
        ax.axvline(l, color="C0", lw=0.8, ls=":"); ax.text(l + 2, top - 0.45, n, fontsize=8, color="C0", rotation=90, va="bottom")
    ax.set_xlim(3820, 5000); ax.set_ylim(0.2, top); ax.set_yticks([]); ax.set_xlabel("wavelength (A, vacuum)"); ax.set_ylabel("scaled flux + offset")
    ax.set_title("Orange: C II line positions; blue dotted: hydrogen Balmer lines", fontsize=9)
    plt.tight_layout(); plt.savefig(out("carbon", "hot_dq_comparison_sdssv.png" if extra else "hot_dq_comparison_5208047381438507520.png"), dpi=100); plt.close()


def gas_discs():
    """Ca II triplet emission: WD 0856+048 (Gaia DR3 578709631539357440) by epoch and in X-shooter, and WD J1959+2208 (Gaia DR3
    1827014701883095680) in its two SDSS-V visits. Spectra normalised as in gas_disc_epochs.py."""
    import gas_disc_epochs as GE
    CAT = list(GE.CAT); W = GE.W
    def binned(n, v, lo=8420, hi=8720, step=3.0):
        e = np.arange(lo, hi + step, step); c = 0.5 * (e[1:] + e[:-1]); k = np.digitize(W, e) - 1; ok = (k >= 0) & (k < len(c)) & (v > 0) & np.isfinite(n)
        num = np.bincount(k[ok], (n * v)[ok], len(c)); den = np.bincount(k[ok], v[ok], len(c))
        return c, np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
    # WD 0856+048: spectra by epoch and equivalent width against date
    g, sid, ra, dec = "578709631539357440", "55774610", 134.841202, 4.636784
    co, per, sp, xs = GE.all_spectra(g, sid, ra, dec); ep = pd.read_csv(f"{T}/gas_disc_epochs_{g}.csv")
    boss = [x for x in sp if x["dataset"] == "BOSS-DR17"][0]; desi = [x for x in sp if x["dataset"] == "DESI-DR1"][0]
    rows = [(f"BOSS {boss['date']}", boss, "0.45"), ("SDSS-V coadd 2021-2023 (7 visits)", co, "C0"), (f"DESI {desi['date']}", desi, "C2"), (f"X-shooter {xs[-1]['date']}", xs[-1], "C3")]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw=dict(width_ratios=[1.5, 1]))
    for i, (lab, s_, col) in enumerate(rows):
        c, y = binned(*GE.norm(s_["w"], s_["f"], s_["iv"])); ax[0].plot(c, y + 0.55 * i, color=col, drawstyle="steps-mid", lw=1)
        ax[0].text(8424, 1.3 + 0.55 * i, lab, color=col, fontsize=8)
    for l in CAT:
        ax[0].axvline(l, color="0.5", lw=0.6, ls=":")
    ax[0].set_xlim(8420, 8720); ax[0].set_xlabel("vacuum wavelength (A)"); ax[0].set_ylabel("normalised flux + offset"); ax[0].set_yticks([])
    ax[0].set_title("Ca II triplet (dotted: 8500.35, 8544.44, 8664.52 A), 3 A bins", fontsize=8)
    e = ep[ep.dataset != "SDSS-V coadd"].copy(); e["year"] = 2000 + (e.mjd - 51544.5) / 365.25
    for ds, col in (("SDSS-DR17", "0.45"), ("BOSS-DR17", "0.45"), ("SDSS-V visit", "C0"), ("DESI-DR1", "C2"), ("ESO XSHOOTER VIS", "C3")):
        q = e[e.dataset == ds]; ax[1].errorbar(q.year, q.ew_A, q.ew_err_A, fmt="o", ms=4, capsize=2, color=col, label=ds)
    ax[1].axhline(0, color="0.6", lw=0.5); ax[1].set_xlabel("year"); ax[1].set_ylabel("emission equivalent width, 3 lines (A)"); ax[1].legend(fontsize=7)
    ax[1].set_title("gas_disc_epochs_578709631539357440.csv", fontsize=8)
    fig.suptitle("WD 0856+048 = Gaia DR3 578709631539357440", fontsize=9); plt.tight_layout(); plt.savefig(out("gas_discs", f"{g}_epochs.png"), dpi=90); plt.close()
    # X-shooter 2025: line profiles and photospheric metal lines
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.4), gridspec_kw=dict(width_ratios=[1.6, 1, 1]))
    n, v = GE.norm(xs[-1]["w"], xs[-1]["f"], xs[-1]["iv"])
    for l, col in zip(CAT, ("C0", "C1", "C3")):
        m = np.abs(W / l - 1) * C < 1500; ax[0].plot((W[m] / l - 1) * C, gaussian_filter1d(np.nan_to_num(n[m], nan=1), 1), color=col, lw=0.9, label=f"Ca II {l:.2f}")
    ax[0].axvline(0, color="0.5", lw=0.5); ax[0].legend(fontsize=7); ax[0].set_xlabel("velocity (km/s, topocentric)"); ax[0].set_ylabel("normalised flux")
    ax[0].set_title(f"X-shooter VIS {xs[-1]['date']} ({xs[-1]['identifier']})", fontsize=8)
    uvb = GE.eso_xshooter(ra, dec, arm="UVB")[-1]
    for a_, (l, nm) in zip(ax[1:], ((4482.58, "Mg II 4481"), (3934.78, "Ca II K"))):
        m = np.abs(uvb["w"] / l - 1) * C < 2500; mc = m & (np.abs(uvb["w"] / l - 1) * C > 600); p = np.polyfit(uvb["w"][mc], uvb["f"][mc], 1)
        a_.plot((uvb["w"][m] / l - 1) * C, gaussian_filter1d(uvb["f"][m] / np.polyval(p, uvb["w"][m]), 2), "k", lw=0.8); a_.axvline(0, color="0.5", lw=0.5)
        a_.set_title(f"{nm}, X-shooter UVB {uvb['date']}", fontsize=8); a_.set_xlabel("velocity (km/s, topocentric)")
    plt.tight_layout(); plt.savefig(out("gas_discs", f"{g}_xshooter.png"), dpi=90); plt.close()
    # WD J1959+2208: SDSS-V visits
    g, sid = "1827014701883095680", "63867520"; co, per = GE.sdssv_spectra(sid)
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.4), gridspec_kw=dict(width_ratios=[1.4, 1]))
    c, y = binned(*GE.norm(co["w"], co["f"], co["iv"]), step=2.0); ax[0].plot(c, y, "k", drawstyle="steps-mid", lw=0.9)
    for l in CAT:
        ax[0].axvline(l, color="0.5", lw=0.6, ls=":")
    ax[0].set_xlim(8420, 8720); ax[0].set_xlabel("vacuum wavelength (A)"); ax[0].set_ylabel("normalised flux"); ax[0].set_title(f"SDSS-V coadd (sdss_id {sid}, 2 visits), 2 A bins", fontsize=8)
    for k, s_ in enumerate(per):
        n, v = GE.norm(s_["w"], s_["f"], s_["iv"])
        for l, ls in zip(CAT, ("-", "--", ":")):
            m = np.abs(W / l - 1) * C < 1500; ax[1].plot((W[m] / l - 1) * C, gaussian_filter1d(np.nan_to_num(n[m], nan=1), 1) + 0.6 * k, color=f"C{k}", ls=ls, lw=0.9,
                                                       label=f"{s_['date']} Ca II {l:.0f}")
    ax[1].axvline(0, color="0.5", lw=0.5); ax[1].legend(fontsize=6, ncol=2); ax[1].set_xlabel("velocity (km/s)"); ax[1].set_title("each visit, three lines overlaid (offset by visit)", fontsize=8)
    fig.suptitle("WD J1959+2208 = Gaia DR3 1827014701883095680", fontsize=9); plt.tight_layout(); plt.savefig(out("gas_discs", f"{g}_sdssv.png"), dpi=90); plt.close()


if __name__ == "__main__":
    import sys
    ALL = dict(zeeman=zeeman, carbon_optical=carbon_optical, carbon_screen_spectra=carbon_screen_spectra, carbon_cos=carbon_cos, galex=galex,
               eclipse=eclipse, balmer=balmer, zz=zz, periodic=periodic, periodic_white_dwarfs=periodic_white_dwarfs,
               hot_dq_comparison=hot_dq_comparison, hot_dq_comparison_sdssv=lambda: hot_dq_comparison(extra=True), gas_discs=gas_discs)
    for name in (sys.argv[1:] or ALL):
        ALL[name](); print(name, "done", flush=True)
