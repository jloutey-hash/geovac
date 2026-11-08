"""GY: Cauchy-radius check: Taylor coefficients from two circles (rho = 1 and rho = 2) must agree; the kite's
integrand is analytic for |z| < 4 (the density pair's M vanishes first at z = -4; the Gaussian pair has no finite point).
usage: python gy_radius_check.py <json1> <json2> [nmax] [key]"""
import json, sys
from pathlib import Path
import mpmath as mp
a = json.loads(Path(sys.argv[1]).read_text()); b = json.loads(Path(sys.argv[2]).read_text())
nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 100
key = sys.argv[4] if len(sys.argv) > 4 else "phi_n"
mp.mp.dps = 290
worst = mp.mpf(0)
for n in range(min(nmax, len(a[key]), len(b[key]))):
    x, y = mp.mpf(a[key][n]), mp.mpf(b[key][n])
    if abs(y) < mp.mpf(10) ** -200:      # an exact zero (phi_b,0 = 0) shows only its noise floor
        continue
    worst = max(worst, abs(x - y) / abs(y))
msg = f"{key}: rho={a['rho']} vs rho={b['rho']}, {len(a[key])}/{len(b[key])}-coefficient runs, worst relative difference over n<{nmax}: {mp.nstr(worst, 3)}"
print(msg, flush=True)
with open(Path(__file__).parent / "data" / "gy_radius_check.log", "a") as f:
    f.write(msg + "\n")
