"""Object pages and the README object index.

Builds docs/objects/<gaia_dr3>.md (one page per object: description, measurement, supporting data) from the measurement
tables, and rewrites the README section between "<!-- object-index:start -->" and "<!-- object-index:end -->".
All numbers are taken verbatim from tables/*.csv; descriptions repeat the topic pages.
Usage: python object_pages.py (run from scripts/)."""
import os, glob, collections
import pandas as pd
import numpy as np

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(R, "docs", "objects")
os.makedirs(OUT, exist_ok=True)
T = lambda f: pd.read_csv(os.path.join(R, "tables", f), dtype=str).fillna("")

def pct(x, nd=1):
    try: return f"{float(x) * 100:.{nd}f}%"
    except (TypeError, ValueError): return ""

def cat(r, *cols):
    parts = []
    for c, lab in cols:
        v = str(r.get(c, "")).strip()
        if v and v.lower() not in ("nan", "none", "-", ""): parts.append(f"{lab} {v}")
    return "; ".join(parts) if parts else "none found"

# Descriptions as on the topic pages (docs/periodic.md, docs/hot_wd_periods.md, docs/dae_wd_periods.md).
DESC = {
    "6021870154194477312": "White dwarf with a 103.4-min period; three Gaia sources within 13 arcsec are fitted separately.",
    "3107374277060584064": "Hot white dwarf with He II 4686 absorption (SIMBAD WD* DO:, MWDD DO:, SDSS-V SnowWhite DA:; VSX type WD without a period); H-alpha, H-beta and Ca II emission whose velocity follows the photometric phase.",
    "2883364038621038208": "White dwarf selected by its Gaia DR3 GLS frequency (GALEX J060343.7-380911).",
    "6722639595190126208": "Hot massive DA white dwarf (GF21 H-atmosphere 36.0 kK, 1.20 Msun; Gaia XP fit in the MWDD 63.7 kK, 1.30 Msun), G = 16.79, 120 pc; the period was found in TESS 2-min light curves (sectors 93 and 104) and recovered in ATLAS.",
    "3161618477052648192": "DC white dwarf (GF21 H-atmosphere 7.6 kK, 0.77 Msun), G = 17.40, 59 pc; the period was found in ZTF and recovered in three TESS sectors and in ATLAS; the highest TESS peak (7.52 c/d) belongs to another star in the aperture. Steen et al. (2024) list it as a likely spotted variable at P = 0.4914 h, the one-cycle-per-day alias of this period.",
    "4851800979770492544": "DA white dwarf (GF21 H-atmosphere 13.9 kK, 0.35 Msun; Gaia XP fit in the MWDD 14.5 kK, 0.31 Msun), G = 17.95, 322 pc; no other Gaia source within 28 arcsec.",
    "178685757799822080": "White dwarf selected by its Gaia DR3 GLS frequency (GALEX J043613.3+383720).",
    "6456720612064924928": "White dwarf selected by its Gaia DR3 GLS frequency (GALEX J211204.8-571801).",
    "2888030331609338240": "White dwarf selected by its Gaia DR3 GLS frequency (GALEX J054140.8-362248).",
    "3496637913394359680": "White dwarf selected by its Gaia DR3 GLS frequency (GALEX J124819.8-261413).",
    "437628614520520320": "White dwarf selected by its Gaia DR3 GLS frequency (WDJ025503.24+475833.96).",
    "6639666736903611136": "White dwarf selected by its Gaia DR3 GLS frequency (GALEX J191430.4-572023).",
    "3890059941364406144": "DA with narrow Balmer lines; Gentile Fusillo et al. (2021) H-atmosphere fit 22,100 K, 0.32 Msun.",
    "974895286283420160": "DA at 51 pc; the TESS period is listed by Oliveira da Rosa et al. (2024).",
    "2795150147707769728": "DA; Gentile Fusillo et al. (2021) H-atmosphere fit 27,000 K, 1.08 Msun.",
    "6170660401283991680": "Gaia XP class DO (Vincent et al. 2024).",
    "6136817910121524096": "Hot-subdwarf candidate in Geier et al. (2019); GALEX J132200.5-422412.",
    "3123625093275668736": "White dwarf and M dwarf in Rebassa-Mansergas et al. (2025); SDSS-V SnowWhite DA_MS.",
    "3354819845628139904": "Hot-subdwarf candidate in Geier et al. (2019).",
    "3230486971974872192": "DO white dwarf, Teff 105.6 kK, log g 8.0 (Kilic et al. 2026a, via the MWDD).",
    "1094376947131876352": "DESI DR1 class DOA, Teff 78.2 kK (Swan et al. 2026).",
    "1038176780370360576": "SBSS 0910+584; spectral type DO (MWDD). G = 17.73, 810 pc.",
    "1157401396015448960": "GALEX J151215.7+065156; UHE white dwarf, DOZ (Reindl et al. 2021, where the period is published). G = 17.22, 990 pc.",
    "1879989790567353344": "GALEX J221519.8+253059; spectral type DOZ (MWDD). G = 17.05, 1190 pc.",
    "920621124593362816": "KUV 07523+4017; DOZ / PG 1159 (Reindl et al. 2021, where the period is published). G = 17.80, 1052 pc.",
    "953685015492787456": "Gaia XP class DO (Vincent et al. 2024). G = 17.53, 483 pc.",
    "5157333438398813824": "Spectral type DA with a photometric temperature above 100 kK (MWDD). G = 17.37, 1072 pc.",
    "1415911839725510528": "Spectral type DA, 68.0 kK (Kilic et al. 2026, via the MWDD). G = 17.59, 833 pc.",
    "2311285729210966144": "GALEX J234931.8-353916; SDSS-V DR20 SnowWhite class DA. G = 17.88, 882 pc.",
    "3365371721281530880": "No spectrum found. G = 18.10, 1100 pc.",
    "1420761029600606592": "DESI DR1 class DAe; GF21 15.0 kK, 0.11 Msun; W1, W2 excess 1.8, 1.9 mag. G = 16.26, 403 pc. SDSS J1724+5620, post-common-envelope binary with the orbital period published by Rebassa-Mansergas et al. (2008).",
    "1337970174853051392": "DESI DR1 class DAE; GF21 12.1 kK, 0.17 Msun; W1, W2 excess 2.0, 1.7 mag. G = 18.94, 741 pc.",
    "4409006786607484672": "DESI DR1 class DAE; GF21 22.4 kK, 0.24 Msun; W1, W2 excess 2.1, 2.6 mag. G = 18.95, 1233 pc.",
    "3717349170269867520": "DESI DR1 class DAe; GF21 22.8 kK, 0.34 Msun; W1, W2 excess 1.9, 2.0 mag. G = 18.99, 955 pc.",
    "1769157090045264128": "DESI DR1 class DAE; GF21 16.0 kK, 0.24 Msun; W1, W2 excess 1.9, 1.6 mag. G = 19.45, 1010 pc.",
    "709815329316284928": "DESI DR1 class DAE; GF21 13.2 kK, 0.25 Msun; W1, W2 excess 2.1, 2.4 mag. G = 19.26, 742 pc.",
    "926161868627454976": "DESI DR1 class DAe; GF21 18.9 kK, 0.21 Msun; W1, W2 excess 2.2, 2.9 mag. G = 19.40, 1449 pc.",
}

Entry = collections.namedtuple("Entry", "topic title page desc obs short tables datafiles scripts extra")
objects = collections.defaultdict(lambda: {"name": "", "entries": []})

def add(gaia, name, e):
    o = objects[str(gaia)]
    if name and (not o["name"] or o["name"].startswith("Gaia DR3")): o["name"] = name
    o["entries"].append(e)

def period_topic(csv, topic, title, page, script, descfun):
    d = T(csv)
    for gaia, g in d.groupby("gaia_dr3", sort=False):
        r = g.iloc[0]
        obs = (f"P = {r.period_h} ± {r.e_period_h} h ({r.dataset}, n = {r.n}, BJD {r.bjd_first}-{r.bjd_last}); "
               f"semi-amplitude {pct(r.amplitude_frac)} ± {pct(r.e_amplitude_frac)}")
        rows = ["| data set | n | semi-amplitude | first harmonic |", "|---|---|---|---|"] + [
            f"| {x.dataset} | {x.n} | {pct(x.amplitude_frac, 2)} ± {pct(x.e_amplitude_frac, 2)} | {pct(x.harmonic2_frac, 2) or '-'} |"
            for _, x in g.iterrows()]
        add(gaia, r["name"], Entry(topic, title, page, descfun(gaia, r), obs,
                                   f"P = {r.period_h} h, {pct(r.amplitude_frac)} ({r.dataset})",
                                   [csv], None, [script], "\n".join(rows)))

period_topic("hot_wd_periods.csv", "hot_wd_periods", "Day-scale periods of hot white dwarfs", "hot_wd_periods.md",
             "hot_dae_wd_periods.py", lambda g, r: DESC.get(g, "Hot white dwarf."))
period_topic("dae_wd_periods.csv", "dae_wd_periods", "Periods of DA white dwarfs with emission lines", "dae_wd_periods.md",
             "hot_dae_wd_periods.py", lambda g, r: DESC.get(g, "DA white dwarf with emission lines (DESI DR1)."))
period_topic("periodic_white_dwarfs.csv", "periodic", "Photometric periods", "periodic.md",
             "periodic_white_dwarfs.py", lambda g, r: DESC.get(g, "White dwarf selected by its Gaia DR3 GLS frequency."))

# Photometric periods: the two objects with their own tables.
d = T("periodic_6021870154194477312.csv")
add("6021870154194477312", "Gaia DR3 6021870154194477312", Entry(
    "periodic", "Photometric periods", "periodic.md", DESC["6021870154194477312"],
    "103.4-min period in ATLAS, Gaia DR3 epoch photometry and TESS (frequency, amplitudes and times of maximum per data set in the table).",
    "103.4-min period (ATLAS, Gaia, TESS)", ["periodic_6021870154194477312.csv"], None,
    ["periodic_6021870154194477312.py", "tess_periodogram.py 1251484163 65 120 0.5 50"], None))
d = T("reflection_3107374277060584064.csv")
add("3107374277060584064", "WDJ064438.09-004550.51", Entry(
    "periodic", "Photometric periods", "periodic.md", DESC["3107374277060584064"],
    "P = 0.59288582 d (14.2293 h) in CoRoT (2007-2012), Gaia DR3 epoch photometry and ZTF; semi-amplitude 11.6% (ZTF r), "
    "3.9% (ZTF g), 10.2% (Gaia G); emission velocities +242/-166/+152 km/s (H-alpha) at phases 0.67/0.23/0.90 from maximum light.",
    "P = 14.229 h; emission follows the phase",
    ["reflection_3107374277060584064.csv", "reflection_3107374277060584064_visits.csv"], None,
    ["reflection_3107374277060584064.py"], None))

# Irradiated companions.
d = T("irradiated_companions.csv"); dp = T("irradiated_companions_periods.csv")
NEW_PERIOD = {"4996506979251027584", "303768056000635776", "5467851842959399808", "6177529630243170432", "2482810406432480512"}
def ratio_err(gaia):
    q = dp[dp.gaia_dr3 == gaia].set_index("dataset")
    if "Gaia RP" not in q.index or "Gaia BP" not in q.index: return ""
    rp, bp = q.loc["Gaia RP"], q.loc["Gaia BP"]; ar, er, ab, eb = (float(x) for x in (rp.semi_amplitude_pct, rp.e_semi_amplitude_pct, bp.semi_amplitude_pct, bp.e_semi_amplitude_pct))
    v = ar / ab; return f"{v:.1f} ± {v * ((er / ar) ** 2 + (eb / ab) ** 2) ** 0.5:.1f}"
for _, r in d.iterrows():
    desc = (f"Low-mass white dwarf: GF21 H-atmosphere Teff {r.gf21_teff_H} K, {r.gf21_mass_H} Msun; "
            f"G = {r.G}, {r.distance_pc} pc.")
    w1 = f"; W1 {r.W1_ratio}x the white-dwarf model (companion M_W1 = {r.M_W1_companion})" if r.W1_ratio else ""
    rerr = ratio_err(r.gaia_dr3) or str(r.red_to_blue_amplitude)
    obs = f"P = {r.period_min} min; red-to-blue (Gaia RP/BP) semi-amplitude ratio {rerr}{w1}"
    add(r.gaia_dr3, r["name"], Entry("irradiated_companions", "Short-period white dwarfs with irradiated companions",
        "irradiated_companions.md", desc, obs, (f"**P = {r.period_min} min**; " if r.gaia_dr3 in NEW_PERIOD else f"P = {r.period_min} min; ") + (f"**W1 {r.W1_ratio}x model**" if r.W1_ratio else f"red/blue amplitude {rerr}"),
        ["irradiated_companions.csv", "irradiated_companions_periods.csv"], None, ["irradiated_companions.py"], None))

# Gaseous discs.
d = T("gas_disc_white_dwarfs.csv")
for _, r in d.iterrows():
    desc = f"White dwarf, G = {r.G}; catalogued: {cat(r, ('snowwhite_class', 'SnowWhite'), ('simbad_type', 'SIMBAD'), ('mwdd_spectype', 'MWDD'), ('desi_dr1_class', 'DESI DR1'))}."
    peaks = f", peak velocities {r.v_blue_kms}/{r.v_red_kms} km/s" if r.v_blue_kms else ""
    ew = f"; EW {r.ew_sdssv_coadd_A} ± {r.e_ew_sdssv_coadd_A} A in the SDSS-V coadd" if r.ew_sdssv_coadd_A else ""
    obs = f"Ca II triplet emission ({r.profile}{peaks}){ew}; other spectra: {r.other_spectra or 'none'}"
    add(r.gaia_dr3, r["name"], Entry("gas_discs", "Ca II triplet emission (gaseous discs)", "gas_discs.md", desc, obs,
        f"Ca II emission, EW {r.ew_sdssv_coadd_A} A" if r.ew_sdssv_coadd_A else f"Ca II emission ({r.profile})", ["gas_disc_white_dwarfs.csv", f"gas_disc_epochs_{r.gaia_dr3}.csv"],
        None, ["gas_disc_epochs.py"], None))

# Carbon lines.
d = T("carbon_white_dwarfs.csv")
for _, r in d.iterrows():
    desc = f"White dwarf, G = {r.G}; catalogued: {cat(r, ('simbad_type', 'SIMBAD'), ('mwdd_spectype', 'MWDD'))}."
    cos = "; C II/C III also in HST/COS" if str(r.HST_COS).strip() not in ("", "nan", "no") else ""
    obs = f"Carbon lines in the SDSS-V spectra: CCF contrast C II {r.ccf_CII_contrast}, C I {r.ccf_CI_contrast} (velocities {r.ccf_CII_v_kms}/{r.ccf_CI_v_kms} km/s){cos}"
    add(r.gaia_dr3, r["name"], Entry("carbon", "Carbon lines", "carbon.md", desc, obs,
        f"C II/C I lines (CCF {r.ccf_CII_contrast}/{r.ccf_CI_contrast})",
        ["carbon_white_dwarfs.csv", f"carbon_{r.gaia_dr3}_optical_CII_features.csv", f"carbon_{r.gaia_dr3}_cos_features.csv"],
        None, ["carbon_lines.py"], None))

# Zeeman splitting.
d = T("magnetic_zeeman.csv")
for _, r in d.iterrows():
    desc = f"DA white dwarf, G = {r.G}; catalogued: {cat(r, ('snowwhite_class', 'SnowWhite'), ('simbad_type', 'SIMBAD'), ('mwdd_spectype', 'MWDD'))}." + (f" {r.note[0].upper() + r.note[1:]}{'' if r.note.endswith('.') else '.'}" if r.note else "")
    obs = (f"Zeeman-split Balmer lines: B = {r.B_split_Ha_MG} ± {r.e_B_split_Ha_MG} MG (H-alpha), "
           f"{r.B_split_Hb_MG} ± {r.e_B_split_Hb_MG} MG (H-beta); spectrum S/N {r.snr}")
    inconsistent = abs(float(r.B_split_Ha_MG) - float(r.B_split_Hb_MG)) > 1.0
    if inconsistent: obs += "; the two lines disagree by more than 1 MG, so the field is not established"
    add(r.gaia_dr3, r["name"], Entry("zeeman", "Zeeman splitting", "zeeman.md", desc, obs,
        f"B = {r.B_split_Ha_MG} MG (H-alpha) vs {r.B_split_Hb_MG} MG (H-beta): inconsistent" if inconsistent else f"B = {r.B_split_Ha_MG} MG",
        ["magnetic_zeeman.csv"], None, ["zeeman_split.py"], None))

# TESS amplitude spectra (ZZ Ceti candidates).
zo, zs = T("zz_ceti_objects.csv"), T("zz_ceti_tess_sectors.csv")
for _, r in zo.iterrows():
    desc = (f"White dwarf, G = {r.G}; SnowWhite {r.snowwhite_class} (Teff {r.snowwhite_teff} K, log g {r.snowwhite_logg}); "
            f"catalogued: {cat(r, ('simbad_type', 'SIMBAD'), ('mwdd_spectype', 'MWDD'))}.")
    g = zs[zs.gaia_dr3 == r.gaia_dr3]
    per = "; ".join(f"sector {x.sector}: {x.period_s} s at {x.amp_ppt} ppt (FAP {x.fap_baluev})" for _, x in g.iterrows())
    px = g[g.pixel_chi2_target.astype(str).str.strip() != ""]
    pix = "; ".join(f"pixel-level test (sector {x.sector}): chi2 {x.pixel_chi2_target} at the target vs {x.pixel_chi2_best_other} "
                    f"at the best other star (G = {x.pixel_best_other_G}, {x.pixel_best_other_sep_arcsec} arcsec)" for _, x in px.iterrows())
    obs = f"TESS pulsation signal(s): {per}" + (f"; {pix}" if pix else "")
    best = (px if len(px) else g[g.cadence_s == 120]).sort_values("fap_baluev").iloc[0]
    add(r.gaia_dr3, r["name"], Entry("zz_ceti", "TESS amplitude spectra", "zz_ceti.md", desc, obs,
        f"pulsations, {best.period_s} s (sector {best.sector}, FAP {best.fap_baluev})", ["zz_ceti_objects.csv", "zz_ceti_tess_sectors.csv"],
        None, ["tess_periodogram.py", "tess_pixel_test.py"], None))

# Eclipses and Balmer emission.
r = T("eclipsing_4731701084150029824.csv").iloc[0]
add(r.gaia_dr3, r["name"], Entry("eclipse_and_emission", "Eclipse and Balmer emission", "eclipse_and_emission.md",
    f"Star with M-dwarf colours (BP-RP {r.bp_rp}, M_G {float(r.G) - 5 * np.log10(1000 / float(r.parallax_mas)) + 5:.2f}), G = {r.G}; the eclipsed object is not identified; catalogued: {cat(r, ('simbad_type', 'SIMBAD'), ('mwdd_spectype', 'MWDD'))}.",
    f"Eclipses in ATLAS: P = {r.period_d} ± {r.e_period_d} d, total duration {r.total_duration_min} min, depths {r.depth_o_uJy}/{r.depth_c_uJy} uJy (o/c)",
    f"eclipses, P = {r.period_d} d", ["eclipsing_4731701084150029824.csv"], None, ["j0353_eclipse.py"], None))
r = T("eclipsing_4851800979770492544.csv").iloc[0]
add(r.gaia_dr3, r["name"], Entry("eclipse_and_emission", "Eclipse and Balmer emission", "eclipse_and_emission.md",
    DESC.get(r.gaia_dr3, "White dwarf."),
    f"Eclipses in TESS (sectors {r.sectors}; {r.n_eclipses_timed} eclipses timed) and ATLAS: P = {r.period_d} ± {r.e_period_d} d, T0 = BJD_TDB {r.T0_bjd_tdb}, "
    f"depth {r.depth_fraction} of the star's flux in the 1-min TESS profile ({r.minutes_below_half} min below half flux); ATLAS in-eclipse flux {r.atlas_o_in_eclipse_mean} ± {r.atlas_o_in_eclipse_err} (o), {r.atlas_c_in_eclipse_mean} ± {r.atlas_c_in_eclipse_err} (c) of the star",
    f"eclipses, P = {r.period_min} min", ["eclipsing_4851800979770492544.csv", "eclipsing_4851800979770492544_profile.csv"], None, ["eclipse_4851800979770492544.py"], None))
for _, r in T("balmer_emission.csv").iterrows():
    obs = (f"H-alpha emission EW {r.Halpha_EW_A} ± {r.e_Halpha_EW_A} A (peak separation {r.Halpha_peak_sep_kms} km/s); "
           f"H-beta EW {r.Hbeta_EW_A} ± {r.e_Hbeta_EW_A} A")
    add(r.gaia_dr3, r["name"], Entry("eclipse_and_emission", "Eclipse and Balmer emission", "eclipse_and_emission.md",
        f"G = {r.G}; catalogued: {cat(r, ('simbad_type', 'SIMBAD'), ('mwdd_spectype', 'MWDD'))}.", obs,
        f"Balmer emission, EW {r.Halpha_EW_A} A", ["balmer_emission.csv"], None, ["cv_balmer.py"], None))

# Hot white dwarfs with He II lines (targets; the table also carries the control stars).
d = T("hot_white_dwarfs.csv")
for c in ("teff_kK", "teff_min_kK", "teff_max_kK"): d[c] = d[c].astype(float)
TEFF_NOT_CONSTRAINED = {"6365804611201098368", "5671975077144346112"}  # as judged on docs/hot_white_dwarfs.md (grid ceiling; S/N)
for gaia, g in d[d.role == "target"].groupby("gaia_dr3", sort=False):
    he = g[g.line_set == "He only"].iloc[0]
    fits = "; ".join(f"{x.line_set}: Teff {x.teff_kK} kK ({x.teff_min_kK}-{x.teff_max_kK}), log g {x.logg}" for _, x in g.iterrows())
    lit = g.iloc[0].lit_teff_kK
    obs = f"He II 4686 and Balmer absorption (DAO); TMAP model fits - {fits}" + (f"; literature Teff {lit} kK" if str(lit).strip() else "")
    add(gaia, g.iloc[0]["name"], Entry("hot_white_dwarfs", "Hot white dwarfs with He II lines", "hot_white_dwarfs.md",
        "Hot white dwarf; no earlier spectrum found.", obs,
        "DAO (He II 4686 + Balmer); " + ("Teff not constrained" if gaia in TEFF_NOT_CONSTRAINED else f"Teff {he.teff_kK:g} kK (He-only TMAP fit" + (f", {he.teff_min_kK:g}-{he.teff_max_kK:g}" if he.teff_min_kK != he.teff_max_kK else "") + ")"),
        ["hot_white_dwarfs.csv"], None, ["hot_white_dwarfs.py"], None))

# Write the object pages.
PAGE_LEAD = {
    "2249098833310553728": "A catalogued but unclassified 63.7-minute variable: a coherent modulation on an over-luminous "
    "low-mass white dwarf, 1.6 ± 0.3 times larger in Gaia RP than in BP. With the GF21 parameters (16.8 kK, 0.14 Msun) the "
    "15.6% G-band semi-amplitude is more than a passive companion inside its Roche lobe can reflect, so either the companion "
    "is at or beyond its Roche lobe or the primary is hotter than the photometric fit; no spectrum and no infrared "
    "measurement (a red source 3.7 arcsec away dominates WISE) exist, so the companion's nature is open.",
}
ORDER = ["gas_discs", "carbon", "zeeman", "periodic", "zz_ceti", "eclipse_and_emission", "hot_white_dwarfs",
         "irradiated_companions", "hot_wd_periods", "dae_wd_periods"]
for gaia, o in objects.items():
    name = o["name"] or f"Gaia DR3 {gaia}"
    head = name if gaia in name else f"{name} (Gaia DR3 {gaia})"
    L = [f"# {head}", ""]
    if gaia in PAGE_LEAD:
        L += [PAGE_LEAD[gaia], ""]
    figs = sorted(glob.glob(os.path.join(R, "figures", "*", f"{gaia}*.png")))
    datafiles = sorted(glob.glob(os.path.join(R, "data", f"*{gaia}*")))
    for e in sorted(o["entries"], key=lambda e: ORDER.index(e.topic)):
        L += [f"## {e.title}", "", e.desc, "", f"**Measurement.** {e.obs}.", ""]
        if e.extra: L += [e.extra, ""]
        L += [f"Method and context: [{e.page.replace('.md', '')}](../{e.page})", ""]
        tabs = ", ".join(f"[{t}](../../tables/{t})" for t in e.tables if os.path.exists(os.path.join(R, "tables", t)))
        L += [f"Tables: {tabs}. Scripts: {', '.join('`scripts/' + s + '`' for s in e.scripts)} "
              f"(commands in the [reproduction list](../../README.md#reproduction)).", ""]
    if datafiles:
        L += ["## Input data", ""] + [f"- [data/{os.path.basename(f)}](../../data/{os.path.basename(f)})" for f in datafiles] + [""]
    if figs:
        L += ["## Figures", ""] + [f'<img src="../../figures/{os.path.relpath(f, os.path.join(R, "figures"))}" width="700">' for f in figs] + [""]
    open(os.path.join(OUT, f"{gaia}.md"), "w").write("\n".join(L))

# README index between the markers.
TITLES = {"gas_discs": "Ca II triplet emission (gaseous discs)", "carbon": "Carbon lines", "zeeman": "Zeeman splitting",
          "periodic": "Photometric periods", "zz_ceti": "TESS amplitude spectra (pulsation candidates)",
          "eclipse_and_emission": "Eclipse and Balmer emission", "hot_white_dwarfs": "Hot white dwarfs with He II lines",
          "irradiated_companions": "Short-period white dwarfs with irradiated companions",
          "hot_wd_periods": "Day-scale periods of hot white dwarfs", "dae_wd_periods": "Periods of DA white dwarfs with emission lines"}
MARK = {
    "4851800979770492544": 4, "6722639595190126208": 4, "3161618477052648192": 3, "4996506979251027584": 3, "303768056000635776": 3, "5467851842959399808": 3, "6177529630243170432": 3, "2482810406432480512": 3,
    "2249098833310553728": 5,
    "578709631539357440": 4, "6914922055508553984": 4, "6021870154194477312": 4, "1980205739970324224": 4,
    "3107374277060584064": 4, "2191618770599895296": 4, "5503429908930455808": 4, "4844023064578952320": 4,
    "1827014701883095680": 3, "2527617665632689024": 3, "1764314497240770176": 3, "1283510882895711872": 3,
    "5208047381438507520": 3, "6886051830805052288": 3, "883885440381808000": 3, "4847399905305694080": 3,
    "2076678981825545088": 3, "5836110898905253760": 3,
    "3890059941364406144": 3, "2795150147707769728": 3, "6170660401283991680": 3, "3230486971974872192": 3,
    "1094376947131876352": 3, "6492083311194727168": 3, "6558472750993181568": 3, "4731701084150029824": 3,
    "2995107164834343680": 3, "5570041179495992704": 3, "3090786872841030016": 3, "4036084504408126976": 3,
    "6365804611201098368": 3, "5671975077144346112": 3, "4705562733524591232": 3,
    "1038176780370360576": 3, "1879989790567353344": 3, "953685015492787456": 3, "5157333438398813824": 3,
    "1420761029600606592": 3, "1337970174853051392": 3, "4409006786607484672": 3, "3717349170269867520": 3,
    "1769157090045264128": 3, "709815329316284928": 3, "926161868627454976": 3,
    "1379988076130545536": 2, "6465542891501713408": 2, "343958710690034944": 2, "6466745168812781568": 2,
    "2883364038621038208": 2, "178685757799822080": 2, "6456720612064924928": 2, "2888030331609338240": 2,
    "3496637913394359680": 2, "437628614520520320": 2, "6639666736903611136": 2, "6136817910121524096": 2,
    "3123625093275668736": 2, "3354819845628139904": 2,
    "2908195134345338496": 2, "4729763229265811328": 2, "6472153670805656832": 2,
    "2002597083798483200": 2, "1977447164064222976": 2,
    "1415911839725510528": 2, "2311285729210966144": 2, "3365371721281530880": 2,
    "5660016586818424832": 2, "6499095244738784128": 2, "4307667617377160704": 2, "5680077137810906624": 2,
    "4198738558061020928": 2, "5807585134758743040": 2, "2883030508640778752": 2, "6412133010376770560": 2,
    "6430762242043644032": 2, "5848754492268362624": 2, "428300220431345536": 2, "4241409569220727424": 2,
    "50526755482496000": 2, "1199447137276626048": 2, "4236646794083461120": 2, "2833867392391927936": 2,
    "3115382600062991872": 2, "2069622487994113408": 2, "657989024107287168": 2, "4680104332056706304": 2,
    "4800536829945050752": 2, "4187365308538865792": 2, "375892788968158848": 2, "4833309388918586624": 2,
    "6465559379883395328": 2, "5665371272869153152": 2, "3327361677328480256": 2, "4679463733391272448": 2,
    "3290180587821828480": 2,
    "974895286283420160": 1, "1157401396015448960": 1, "920621124593362816": 1,
}

NOBOLD = {"974895286283420160", "1157401396015448960", "920621124593362816", "1420761029600606592"}

idx = []
for topic in ORDER:
    rows = [(o["name"] or f"Gaia DR3 {g}", g, e) for g, o in objects.items() for e in o["entries"] if e.topic == topic]
    if not rows: continue
    idx += [f"### {TITLES[topic]} ({len(rows)})", "", "| object | description | measurement |", "|---|---|---|"]
    for name, g, e in rows:
        m = MARK.get(g, 0)
        marks = " " + "\\*" * m if m else ""
        s = e.short if g in NOBOLD or "**" in e.short or "inconsistent" in e.short or "red/blue amplitude" in e.short else f"**{e.short}**"
        idx.append(f"| [{name}](docs/objects/{g}.md){marks} | {e.desc} | {s} |")
    idx.append("")
p = os.path.join(R, "README.md"); t = open(p).read()
a, b = "<!-- object-index:start -->", "<!-- object-index:end -->"
t = t[:t.index(a) + len(a)] + "\n" + "\n".join(idx) + t[t.index(b):]
open(p, "w").write(t)
print(len(objects), "object pages;", sum(len(o['entries']) for o in objects.values()), "entries")
