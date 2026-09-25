"""WISE W1/W2 of a white dwarf from the Legacy Surveys DR10 Tractor catalogue (forced, deblended WISE photometry), compared with
a Rayleigh-Jeans extrapolation of its z-band flux: f(W) = f(z) (lambda_z / lambda_W)^2, lambda_z = 0.92, lambda_W1 = 3.37,
lambda_W2 = 4.62 micron. The extrapolation ignores the Balmer/Paschen opacity and is an approximate photospheric level.
Usage: python wise_excess.py <gaia_dr3> <ra> <dec>   (Astro Data Lab TAP; fluxes in nanomaggies, AB)
Prints the matched source (within 1.5 arcsec) and all Tractor sources within 8 arcsec with W1 flux."""
import sys, io, requests, numpy as np, pandas as pd
g, ra, dec = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]); d = 8 / 3600
q = ("SELECT ra, dec, type, flux_g, flux_r, flux_z, flux_ivar_z, flux_w1, flux_ivar_w1, flux_w2, flux_ivar_w2, fracflux_w1, fracflux_w2 "
     f"FROM ls_dr10.tractor WHERE ra BETWEEN {ra - d / np.cos(np.radians(dec))} AND {ra + d / np.cos(np.radians(dec))} AND dec BETWEEN {dec - d} AND {dec + d}")
r = requests.get("https://datalab.noirlab.edu/tap/sync", params=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q), timeout=300)
t = pd.read_csv(io.StringIO(r.text))
if not len(t) or "ra" not in t.columns:
    sys.exit("no Legacy Surveys DR10 sources (outside the footprint or query failed)")
t["sep"] = np.hypot((t.ra - ra) * np.cos(np.radians(dec)), t.dec - dec) * 3600; t = t.sort_values("sep")
print(t[["sep", "type", "flux_g", "flux_r", "flux_z", "flux_w1", "flux_w2", "fracflux_w1"]].round(3).to_string(index=False))
m = t.iloc[0]
if m.sep > 1.5:
    sys.exit("no source within 1.5 arcsec")
e1, e2 = 1 / np.sqrt(m.flux_ivar_w1), 1 / np.sqrt(m.flux_ivar_w2)
rj1, rj2 = m.flux_z * (0.92 / 3.37) ** 2, m.flux_z * (0.92 / 4.62) ** 2
print(f"Gaia DR3 {g}: W1 {m.flux_w1:.2f} +- {e1:.2f} nMgy, W2 {m.flux_w2:.2f} +- {e2:.2f} nMgy; Rayleigh-Jeans from z ({m.flux_z:.2f} nMgy): "
      f"W1 {rj1:.2f}, W2 {rj2:.2f} nMgy; ratio W1 {m.flux_w1 / rj1:.1f}, W2 {m.flux_w2 / rj2:.1f}; excess significance "
      f"W1 {(m.flux_w1 - rj1) / e1:.1f} sigma, W2 {(m.flux_w2 - rj2) / e2:.1f} sigma; fracflux_w1 {m.fracflux_w1:.2f}")
