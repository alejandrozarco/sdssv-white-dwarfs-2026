"""GALEX FUV-NUV versus Gaia BP-RP for a target and comparison samples.
Comparison samples: SDSS-V SnowWhite DA and DB/DBA spectra (plx > 4 mas, -0.55 < BP-RP < -0.25, S/N > 10); MWDD DQ-type white
dwarfs with MWDD Teff >= 15 kK (data/mwdd_dq_teff_ge15kK.csv, MWDD table of 2026-08-05). Gaia DR3 positions propagated to 2007.0;
GUVcat AIS (VizieR II/335/galex_ais) via CDS XMatch within 4 arcsec. Usage: python galex_colours.py <gaia_dr3_id>
Options: --window LO HI (BP-RP range of the comparison samples, default -0.55 -0.25);
--target-galex FUV eFUV NUV eNUV (target photometry from another GALEX catalogue, e.g. GR6/7 GII, replacing the AIS match)."""
import sys, io, numpy as np, pandas as pd, requests
import astropy.units as u
from astropy.table import Table
from astroquery.xmatch import XMatch
from sdssv import cas
target = sys.argv[1]; args = sys.argv[2:]
lo, hi = (float(args[args.index("--window") + 1]), float(args[args.index("--window") + 2])) if "--window" in args else (-0.55, -0.25)
tg = [float(v) for v in args[args.index("--target-galex") + 1:args.index("--target-galex") + 5]] if "--target-galex" in args else None
sw = cas(f"SELECT gaia_dr3_source_id, classification FROM snow_white_boss_star WHERE plx > 4 AND (bp_mag - rp_mag) BETWEEN {lo} AND {hi} "
         "AND snr > 10 AND classification IN ('DA','DB','DBA')")
grp = {str(int(g)): ("DA" if c == "DA" else "DB") for g, c in zip(sw.gaia_dr3_source_id, sw.classification) if g > 0}
for g in pd.read_csv("../data/mwdd_dq_teff_ge15kK.csv", dtype={"gaia_dr3": str}).gaia_dr3:
    grp[g] = "DQ"
grp[target] = "target"; ids = list(grp); rows = []
for i in range(0, len(ids), 300):
    q = ("SELECT source_id, ra, dec, pmra, pmdec, parallax, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag FROM gaiadr3.gaia_source "
         "WHERE source_id IN (" + ",".join(ids[i:i + 300]) + ")")
    r = requests.post("https://gea.esac.esa.int/tap-server/tap/sync", data=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q), timeout=300)
    rows.append(pd.read_csv(io.StringIO(r.text), dtype={"source_id": str}))
g = pd.concat(rows); g["group"] = g.source_id.map(grp)
g["ra07"] = g.ra + g.pmra.fillna(0) * -9 / 3.6e6 / np.cos(np.radians(g.dec)); g["de07"] = g.dec + g.pmdec.fillna(0) * -9 / 3.6e6
x = XMatch.query(cat1=Table.from_pandas(g[["source_id", "ra07", "de07"]]), cat2="vizier:II/335/galex_ais", max_distance=4 * u.arcsec, colRA1="ra07", colDec1="de07").to_pandas()
x["source_id"] = x.source_id.astype(str); x = x.sort_values("angDist").drop_duplicates("source_id")
d = g.merge(x[["source_id", "FUVmag", "e_FUVmag", "NUVmag", "e_NUVmag"]], on="source_id", how="left")
if tg:
    d.loc[d.group == "target", ["FUVmag", "e_FUVmag", "NUVmag", "e_NUVmag"]] = tg
d["bp_rp"] = d.phot_bp_mean_mag - d.phot_rp_mean_mag; d["fuv_nuv"] = d.FUVmag - d.NUVmag; d["M_G"] = d.phot_g_mean_mag + 5 * np.log10(d.parallax / 100)
d.to_csv(f"galex_colours_{target}.csv", index=False)
T = d[d.group == "target"].iloc[0]
print(f"target {target}: BP-RP {T.bp_rp:.3f}, M_G {T.M_G:.2f}, FUV {T.FUVmag:.3f} +- {T.e_FUVmag:.3f}, NUV {T.NUVmag:.3f} +- {T.e_NUVmag:.3f}, FUV-NUV {T.fuv_nuv:+.3f}")
for gname in ("DA", "DB", "DQ"):
    s = d[(d.group == gname) & d.fuv_nuv.notna() & (d.e_FUVmag < 0.1)]
    for dc, mg in ((0.06, None), (0.10, T.M_G - 0.35)):
        n = s[(np.abs(s.bp_rp - T.bp_rp) < dc) & ((s.M_G > mg) if mg is not None else True)]
        if len(n):
            print(f"{gname}: |dBP-RP| < {dc}{'' if mg is None else f', M_G > {mg:.2f}'}: n = {len(n)}, FUV-NUV median {n.fuv_nuv.median():+.2f}, "
                  f"min {n.fuv_nuv.min():+.2f}, max {n.fuv_nuv.max():+.2f}, n >= target {int((n.fuv_nuv >= T.fuv_nuv).sum())}")
