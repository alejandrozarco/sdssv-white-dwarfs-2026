"""Line test on SDSS-V visit spectra. Usage: python carbon_lines.py <sdss_id> <label>
Coadd of all visits (sdssv.coadd); pseudo-continuum = running 85th percentile over 120 A, smoothed; depth = 1 - flux/continuum.
Cross-correlation of the depth spectrum with templates of the 15 strongest NIST vacuum lines per species (Gaussian FWHM 4 A),
-3000..+3000 km/s; contrast = (peak - median) / (1.4826 MAD) of the correlation more than 1000 km/s from the peak.
Gaussian fits to nine C II features; velocities relative to intensity-weighted NIST vacuum wavelengths."""
import sys, json, numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import percentile_filter, gaussian_filter1d
from sdssv import visits, coadd, C
sid, label = sys.argv[1], sys.argv[2]
L = json.load(open("../data/nist_vacuum_lines.json")); L["He II"] = [(4687.02, 100.0), (5413.03, 50.0), (6562.0, 40.0), (4542.9, 20.0)]
SPEC = ["C II", "He I", "H I", "C I", "O I", "O II", "C III", "He II", "Mg II", "Si II"]
FEATS = [(3918, 3924), (4074, 4079), (4266, 4271), (4371, 4378), (4618, 4623), (5889, 5895), (6577, 6587), (6780, 6787), (7230, 7241)]


def prep(w, f, iv):
    m = np.isfinite(f) & (iv > 0); w, f, iv = w[m], f[m], iv[m]
    npix = int(round(np.log10(1 + 120 / 5000) / np.median(np.diff(np.log10(w)))))
    cont = gaussian_filter1d(percentile_filter(f, 85, size=npix), npix / 4)
    return w, f, cont, 1 - f / cont, 1 / np.sqrt(iv) / cont


def strongest(sp, n=15, lo=3850, hi=9200, merge=1.5):
    out = []
    for lam, I in sorted([x for x in L[sp] if lo < x[0] < hi], key=lambda x: -x[1]):
        if all(abs(lam - q) > merge for q in out):
            out.append(lam)
        if len(out) == n:
            break
    return out


def ccf(w, depth, species):
    vels = np.arange(-3000, 3001, 10.0); res = {}
    for sp in species:
        lines = strongest(sp)
        if not lines:
            continue
        cc = np.array([np.corrcoef(depth, sum(np.exp(-0.5 * ((w - l * (1 + v / C)) / 1.7) ** 2) for l in lines))[0, 1] for v in vels])
        k = np.argmax(cc); far = np.abs(vels - vels[k]) > 1000; mad = 1.4826 * np.median(np.abs(cc[far] - np.median(cc[far])))
        res[sp] = dict(v_peak=float(vels[k]), r_peak=round(float(cc[k]), 3), contrast=round(float((cc[k] - np.median(cc[far])) / mad), 1))
    return res


vs = visits(sid)
w, f, cont, depth, edep = prep(*coadd(vs))
out = dict(sdss_id=sid, label=label, visits=[dict(mjd=v["mjd"], snr=round(v["snr"], 1), xcsao_v=round(v["xcsao_v"], 1), in_stack=v["in_stack"]) for v in vs],
           ccf_coadd=ccf(w, depth, SPEC), ccf_visits={})
for v in vs:
    wv, fv, cv, dv, ev = prep(v["wave"], v["flux"], v["ivar"]); out["ccf_visits"][v["mjd"]] = ccf(wv, dv, ["C II", "C I", "He I", "H I"])
g = lambda x, a, mu, s, c0, c1: c0 + c1 * (x - mu) - a * np.exp(-0.5 * ((x - mu) / s) ** 2)
vg = out["ccf_coadd"]["C II"]["v_peak"]; out["c_ii_features"] = []
for lo, hi in FEATS:
    comps = [x for x in L["C II"] if lo <= x[0] <= hi]; lab = sum(a * b for a, b in comps) / sum(b for a, b in comps)
    mu0 = lab * (1 + vg / C); s = np.abs(w - mu0) < 25
    try:
        p, cv = curve_fit(g, w[s], f[s] / cont[s], p0=[0.1, mu0, 3.0, 1, 0], sigma=edep[s], absolute_sigma=True,
                          bounds=([0, mu0 - 8, 0.8, 0.5, -1], [1, mu0 + 8, 15, 1.5, 1]), maxfev=20000)
        e = np.sqrt(np.diag(cv))
        out["c_ii_features"].append(dict(lab_vac=round(lab, 2), center=round(float(p[1]), 2), v=round(float((p[1] / lab - 1) * C)), e_v=round(float(e[1] / lab * C)),
                                         depth=round(float(p[0]), 3), depth_snr=round(float(p[0] / e[0]), 1), sigma_A=round(float(p[2]), 2)))
    except RuntimeError:
        out["c_ii_features"].append(dict(lab_vac=round(lab, 2), fit="failed"))
good = [x for x in out["c_ii_features"] if x.get("depth_snr", 0) > 4 and x.get("e_v", 1e9) < 200]
vm = float(np.average([x["v"] for x in good], weights=[1 / x["e_v"] ** 2 for x in good])) if good else vg
out["c_ii_weighted_mean_v"] = round(vm)
if out["ccf_coadd"]["C I"]["contrast"] > out["ccf_coadd"]["C II"]["contrast"]:
    vm = out["ccf_coadd"]["C I"]["v_peak"]
out["reference_velocity_for_line_depths"] = round(vm)
cii = np.array([x[0] for x in L["C II"] if x[1] >= 30]); rng = np.random.default_rng(1); base = []
for mu in [x for x in rng.uniform(3900, 7400, 4000) if np.min(np.abs(cii * (1 + vm / C) - x)) > 8][:400]:
    s = np.abs(w - mu) < 3; base.append(np.sum(depth[s] / edep[s] ** 2) / np.sum(1 / edep[s] ** 2))
out["depth_baseline_random_windows"] = dict(median=round(float(np.median(base)), 3), p99=round(float(np.percentile(base, 99)), 3))
out["depth_at_lines"] = {}
for nm, lam in [("H-alpha", 6564.632), ("H-beta", 4862.691), ("He I 4472", 4472.735), ("He I 5877", 5877.25), ("He I 6680", 6679.99), ("He II 4687", 4687.02)]:
    s = np.abs(w - lam * (1 + vm / C)) < 3; out["depth_at_lines"][nm] = round(float(np.sum(depth[s] / edep[s] ** 2) / np.sum(1 / edep[s] ** 2)), 3)
json.dump(out, open(f"carbon_lines_{label}.json", "w"), indent=1)
print(json.dumps(out["ccf_coadd"], indent=0)); print(json.dumps(out["c_ii_features"]))
