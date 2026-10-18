"""QR step 3a: theta_t^j f(t), j <= J, at log-spaced t, float64, for the ring of four S (t on d),
the open chain of four (t on the end d, edge da removed), and the controls: Coulomb triangle G
(t on one density) and the chain of three C (t on the middle density).  Partial-wave evaluator
debug/ring_ext_QR_pw.py.  Output debug/data/ring_ext_QR_data.npz"""
import sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QR_pw as pw

J = 10
ts = np.exp(np.linspace(np.log(0.25), np.log(4.0), 64))
g = pw.Grid(xmin=-30, xmax=np.log(240), H=0.03)
out = {"ts": ts}
t0 = time.time()
out["C3"], _ = pw.batch(g, 3, 1, ts, J, 0, closed=False); print("C3", time.time() - t0, flush=True)
out["C4"], _ = pw.batch(g, 4, 3, ts, J, 0, closed=False); print("C4", time.time() - t0, flush=True)
out["G"], _ = pw.batch(g, 3, 2, ts, J, 110); print("G", time.time() - t0, flush=True)
out["S"], _ = pw.batch(g, 4, 3, ts, J, 70); print("S", time.time() - t0, flush=True)
np.savez(Path(__file__).parent / "data" / "ring_ext_QR_data.npz", **out)
