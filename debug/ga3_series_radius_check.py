"""GA3: Cauchy-radius check of the lin kite series at gamma = 2: the Taylor coefficients from circles of radius
rho = 1 and rho = 2 (short runs, ga3_4_series.py 2 256 120 <nodes> <rho> lin) must agree; at gamma = 2 the kite
integrand is analytic for |z| < 4 (M = sqrt(1 + u(1-u) z) vanishes first at z = -4; N = sqrt(4 + v(1-v) z) only
at |z| >= 16), so both circles lie inside the disc of convergence.  The pseudo-threshold x = (1 - gamma)^2 = 1
(t = 1) is not on the principal sheet.
usage: python ga3_series_radius_check.py <rho1 json> <rho2 json> [nmax]"""
import json, sys
from pathlib import Path
import mpmath as mp

a = json.loads(Path(sys.argv[1]).read_text())
b = json.loads(Path(sys.argv[2]).read_text())
nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 100
mp.mp.dps = 130
worst = mp.mpf(0)
for n in range(min(nmax, len(a["phi_n"]), len(b["phi_n"]))):
    x, y = mp.mpf(a["phi_n"][n]), mp.mpf(b["phi_n"][n])
    worst = max(worst, abs(x - y) / abs(y))
msg = f"rho={a['rho']} vs rho={b['rho']}, {len(a['phi_n'])}-coefficient runs, worst relative difference over n<{nmax}: {mp.nstr(worst, 3)}"
print(msg, flush=True)
(Path(__file__).parent / "data" / "ga3_series_radius_check.log").write_text(msg + "\n")
