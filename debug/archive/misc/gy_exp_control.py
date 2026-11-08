"""GY positive control of the Cauchy-series driver: gy_4_series.py in mode 'exp' (the exponential pair of record QG,
ring_ext_QG_common.phi_gem_z, gamma = 3/2) must reproduce the stored gamma = 3/2 series of record QG coefficient by coefficient.
usage: python gy_exp_control.py <new json> <stored json> [nmax]"""
import json, sys
from pathlib import Path
import mpmath as mp
a = json.loads(Path(sys.argv[1]).read_text()); b = json.loads(Path(sys.argv[2]).read_text())
nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 100
mp.mp.dps = 120
worst = mp.mpf(0)
for n in range(min(nmax, len(a["phi_n"]), len(b["phi_n"]))):
    x, y = mp.mpf(a["phi_n"][n]), mp.mpf(b["phi_n"][n])
    worst = max(worst, abs(x - y) / abs(y))
msg = f"exp control: new driver (M={a['Mp']}, dps={a['dps']}, rho={a['rho']}) vs stored QG gamma=3/2 series (M={b['Mp']}, dps={b['dps']}): worst relative difference over n<{nmax}: {mp.nstr(worst, 3)}"
print(msg, flush=True)
(Path(__file__).parent / "data" / "gy_exp_control.log").write_text(msg + "\n")
