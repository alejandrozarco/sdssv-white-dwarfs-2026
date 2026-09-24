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


def periodic():
    import periodic_6021870154194477312 as P6
    g = P6.GID; t = pd.read_csv(f"{T}/periodic_{g}.csv", dtype={"gaia_dr3": str}); f = float(t.frequency_cd.iloc[0])
    t0 = float(t[(t.gaia_dr3 == g) & (t.dataset == "ATLAS") & (t.band == "o")].t_max_bjd.iloc[0])
    at = P6.load_atlas(g); _, (fr, pw) = P6.adopted_frequency(at); ga = P6.load_gaia(); tt, crowd = P6.load_tess()
    fig = plt.figure(figsize=(11, 10)); gs = fig.add_gridspec(3, 2)
    a = fig.add_subplot(gs[0, :]); w, fl, iv = coadd(visits("106653583")); s = np.isfinite(fl) & (iv > 0) & (w > 3750) & (w < 9000)
    a.plot(w[s], gaussian_filter1d(fl[s], 2), "k", lw=0.6)
    for l in (6564.61, 4862.68, 4341.69, 4102.89):
        a.axvline(l, color="0.6", lw=0.6, ls=":")
    a.set_xlabel("wavelength (A)"); a.set_ylabel("flux"); a.set_title(f"Gaia DR3 {g}: SDSS-V spectrum (sdss_id 106653583, all visits); grey dotted: H-alpha to H-delta", fontsize=8)
    a = fig.add_subplot(gs[1, :]); a.plot(fr, pw, "k", lw=0.4); a.set_xlabel("frequency (c/d)"); a.set_ylabel("GLS power")
    a.set_title(f"ATLAS c + o, fractional flux; adopted f = {f:.7f} c/d (P = {1440 / f:.3f} min)", fontsize=8)
    ph = lambda x: ((x - t0) * f) % 1
    a = fig.add_subplot(gs[2, 0])
    for b, col in (("c", "c"), ("o", "orange")):
        p = ph(at[b]["t"]); y = at[b]["f"] / P6.REF[b]; e = at[b]["e"] / P6.REF[b]; edges = np.linspace(0, 1, 21); c = (edges[1:] + edges[:-1]) / 2
        m = [np.sum(y[(p >= lo) & (p < hi)] / e[(p >= lo) & (p < hi)] ** 2) / np.sum(1 / e[(p >= lo) & (p < hi)] ** 2) for lo, hi in zip(edges[:-1], edges[1:])]
        se = [1 / np.sqrt(np.sum(1 / e[(p >= lo) & (p < hi)] ** 2)) for lo, hi in zip(edges[:-1], edges[1:])]
        for k in (0, 1):
            a.errorbar(c + k, m, se, fmt="o", ms=3, color=col, label=f"ATLAS {b}" if k == 0 else None)
    a.set_xlabel(f"phase (t_max = BJD_TDB {t0:.5f})"); a.set_ylabel("fractional flux (20 bins)"); a.legend(fontsize=7)
    a = fig.add_subplot(gs[2, 1]); d = ga["G"]; p = ph(d["t"])
    for k in (0, 1):
        a.errorbar(p + k, d["f"], d["e"], fmt="o", ms=3, color="k", label="Gaia DR3 G" if k == 0 else None)
    p = ph(tt["t"]); edges = np.linspace(0, 1, 21); c = (edges[1:] + edges[:-1]) / 2
    m = [np.mean(tt["f"][(p >= lo) & (p < hi)]) for lo, hi in zip(edges[:-1], edges[1:])]
    for k in (0, 1):
        a.plot(c + k, m, "s", ms=3, color="r", label=f"TESS S{P6.SECTOR} PDCSAP, CROWDSAP {crowd:.3f} (20 bins)" if k == 0 else None)
    a.set_xlabel("phase"); a.set_ylabel("fractional flux"); a.legend(fontsize=7)
    plt.tight_layout(); plt.savefig(f"{F}/periodic_{g}.png", dpi=90); plt.close()


def periodic_white_dwarfs():
    import periodic_white_dwarfs as PW
    t = pd.read_csv(f"{T}/periodic_white_dwarfs.csv", dtype={"gaia_dr3": str}); ids = list(PW.SRC.index)
    fig, ax = plt.subplots(len(ids), 2, figsize=(11, 3.0 * len(ids)))
    for i, gid in enumerate(ids):
        s = PW.SRC.loc[gid]; r = t[(t.gaia_dr3 == gid) & (t.dataset == s.ground)].iloc[0]; f = float(r.frequency_cd); t0 = float(r.t_max_bjd)
        x, y, e, g = PW.load_atlas(gid) if s.ground == "ATLAS" else PW.load_ztf(gid); _, (fr, pw) = PW.adopted_frequency(x, y, e, g)
        ax[i, 0].plot(fr, pw, "k", lw=0.4); ax[i, 0].set_xscale("log"); ax[i, 0].set_ylabel("GLS power")
        ax[i, 0].set_title(f"Gaia DR3 {gid} ({s['name']}): {s.ground}, f = {f:.6f} c/d, P = {24 / f:.4f} h", fontsize=8)
        p = ((x - t0) * f) % 1; edges = np.linspace(0, 1, 21); c = (edges[1:] + edges[:-1]) / 2
        m = [np.sum(y[(p >= lo) & (p < hi)] / e[(p >= lo) & (p < hi)] ** 2) / np.sum(1 / e[(p >= lo) & (p < hi)] ** 2) for lo, hi in zip(edges[:-1], edges[1:])]
        se = [1 / np.sqrt(np.sum(1 / e[(p >= lo) & (p < hi)] ** 2)) for lo, hi in zip(edges[:-1], edges[1:])]
        tg, yg, eg = PW.load_gaia(gid); pg = ((tg - t0) * f) % 1
        for k in (0, 1):
            ax[i, 1].errorbar(c + k, m, se, fmt="o", ms=3, color="C0", label=f"{s.ground} (20 bins)" if k == 0 else None)
            ax[i, 1].errorbar(pg + k, yg, eg, fmt=".", ms=3, color="0.4", alpha=0.7, label="Gaia DR3 G" if k == 0 else None)
        ax[i, 1].set_ylabel("fractional flux"); ax[i, 1].legend(fontsize=7)
    ax[-1, 0].set_xlabel("frequency (c/d)"); ax[-1, 1].set_xlabel("phase (0 = t_max of the ground-based fit)")
    plt.tight_layout(); plt.savefig(f"{F}/periodic_white_dwarfs.png", dpi=90); plt.close()


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
    plt.tight_layout(); plt.savefig(f"{F}/hot_dq_comparison_sdssv.png" if extra else f"{F}/hot_dq_comparison_5208047381438507520.png", dpi=110); plt.close()


if __name__ == "__main__":
    for fn in (zeeman, carbon_optical, carbon_cos, galex, eclipse, balmer, zz, periodic, periodic_white_dwarfs, hot_dq_comparison, lambda: hot_dq_comparison(extra=True)):
        fn(); print(fn.__name__, "done", flush=True)
