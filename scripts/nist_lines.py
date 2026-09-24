"""NIST ASD vacuum wavelengths and relative intensities (3700-9300 A) for the species used in carbon_lines.py."""
import re, json
import astropy.units as u
from astroquery.nist import Nist
out = {}
for sp in ["C II", "C I", "He I", "H I", "O I", "O II", "C III", "Mg II", "Si II"]:
    t = Nist.query(3700 * u.AA, 9300 * u.AA, linename=sp, wavelength_type="vacuum"); rows = []
    for r in t:
        try:
            w = float(str(r["Ritz"]).strip() or str(r["Observed"]).strip())
        except ValueError:
            continue
        m = re.match(r"\s*(\d+)", str(r["Rel."]))
        if m:
            rows.append((w, float(m.group(1))))
    out[sp] = rows
json.dump(out, open("../data/nist_vacuum_lines.json", "w"))
