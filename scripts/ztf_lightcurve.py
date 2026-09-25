"""ZTF light-curve summary at a position: python ztf_lightcurve.py <gaia_dr3> <ra> <dec>
IRSA ZTF light-curve service, 1.8 arcsec radius, catflags = 0, magerr < 0.25; data sets (ZTF object and filter) with at
least 20 points. Per data set: median magnitude, rms and median error of the fractional flux, seasonal medians (season =
year starting in October) with errors 1.25 x median error / sqrt(n), and the number of points below -3 sigma.
Nights with more than 50 points: rms of 10-minute bins and the lowest bin. Generalised Lomb-Scargle of all data sets
together (0.02-300 c/d): the five highest distinct peaks, and for peaks below 50 c/d their rank in each data set separately. Writes ../data/ztf_<gaia_dr3>.csv."""
import sys, io, requests, numpy as np, pandas as pd
from astropy.timeseries import LombScargle
g, ra, dec = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
for k in range(3):
    r = requests.get("https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves", params=dict(POS=f"CIRCLE {ra} {dec} 0.0005", BANDNAME="g,r,i", FORMAT="csv"), timeout=300)
    if r.status_code == 200 and r.text.startswith("oid"):
        break
else:
    sys.exit("IRSA query failed")
df = pd.read_csv(io.StringIO(r.text)); df.to_csv(f"../data/ztf_{g}.csv", index=False)
df = df[(df.catflags == 0) & (df.magerr < 0.25)]; parts = []
for (o, f), s in df.groupby(["oid", "filtercode"]):
    if len(s) < 20:
        continue
    med = np.median(s.mag); fl = 10 ** (-0.4 * (s.mag - med)); e = 0.921 * s.magerr * fl; zz = (fl - 1) / e
    seas = np.floor(2000 + (s.mjd - 51544.5) / 365.25 + 0.25)
    sm = pd.DataFrame(dict(fl=fl, e=e, seas=seas)).groupby("seas").agg(n=("fl", "size"), med=("fl", "median"), err=("e", lambda x: 1.2533 * np.median(x) / np.sqrt(len(x))))
    print(f"{f} {o}: n {len(s)}, median mag {med:.3f}, rms {100 * np.std(fl):.2f}%, median error {100 * np.median(e):.2f}%, points < -3 sigma {int((zz < -3).sum())}")
    print("   seasonal median (%, error, n):", {int(k): (round(100 * (v.med - 1), 2), round(100 * v.err, 2), int(v.n)) for k, v in sm.iterrows()})
    for night, x in s.assign(fl=fl).groupby(np.floor(s.mjd)):
        if len(x) > 50:
            b = x.groupby(np.floor((x.mjd - x.mjd.min()) * 144)).fl.mean() - 1
            print(f"   night MJD {int(night)}: {len(x)} points, 10-min bin rms {100 * b.std():.2f}%, lowest bin {100 * b.min():.2f}%")
    parts.append((f, o, pd.DataFrame(dict(t=s.hjd.values, f=fl.values - 1, e=e.values))))
d = pd.concat([p[2] for p in parts]); fr = np.arange(0.02, 300, 0.1 / (d.t.max() - d.t.min())); ls = LombScargle(d.t, d.f, d.e); p = ls.power(fr, method="fast")
peaks = []
for k in np.argsort(p)[::-1]:
    if all(abs(fr[k] - q) > 0.01 for q in peaks):
        peaks.append(fr[k])
    if len(peaks) == 5:
        break
for f0 in peaks:
    pw = ls.power(f0)
    print(f"all data: peak {f0:.5f} c/d (P {24 / f0:.4f} h), amplitude ~{100 * np.sqrt(4 * pw * np.var(d.f)):.2f}%, "
          f"Baluev FAP {ls.false_alarm_probability(pw, minimum_frequency=0.02, maximum_frequency=300, method='baluev'):.1e}")
    if f0 > 50:
        continue
    for f, o, x in parts:
        if len(x) < 80:
            continue
        l1 = LombScargle(x.t, x.f, x.e); fr1 = np.arange(0.05, 50, 0.1 / (x.t.max() - x.t.min())); p1 = l1.power(fr1, method="fast")
        near = np.abs(fr1 - f0) < 0.002
        print(f"   {f} {o}: rank of the highest power within 0.002 c/d: {int((p1 > p1[near].max()).sum()) + 1} of {len(fr1)}")
