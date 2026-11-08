"""Independent cross-check of a linear-method single-geminal optimum by DIRECT VMC:
sample |g Psi_CI|^2 with g = 1 + c (F - Fbar) and average the standard local energy
(lih_vmc.vmc_linear).  --nsweep TOTAL incl. --nburn.  Writes data/mg_direct_<tag>.npz."""
import os
import sys
import types
import argparse

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
if "geovac" not in sys.modules:
    _p = types.ModuleType("geovac"); _p.__path__ = [os.path.join(os.path.dirname(HERE), "geovac")]
    sys.modules["geovac"] = _p
import lih_vmc
from lih_pi_m16_vmc import _build_wf_from_no

ap = argparse.ArgumentParser()
ap.add_argument('--geminal', default='exp'); ap.add_argument('--gamma', type=float, default=2.0)
ap.add_argument('--c', type=float, required=True); ap.add_argument('--Fbar', type=float, default=0.0)
ap.add_argument('--nwalk', type=int, default=1000); ap.add_argument('--nsweep', type=int, default=8000)
ap.add_argument('--nburn', type=int, default=4000); ap.add_argument('--seed', type=int, default=0)
a = ap.parse_args()
z = np.load(os.path.join(HERE, "data", "lih_pi_m16_energy_trunc1e-2.npz"), allow_pickle=True)
wf = _build_wf_from_no(z['ct'], z['T_no'], int(z['Mk']), float(z['E_trunc']))
lj = lih_vmc.LinearGeminalJastrow(a.geminal, a.gamma, a.c, a.Fbar, wf.na, wf.nb)
E, err, acc = lih_vmc.vmc_linear(wf, lj, nwalk=a.nwalk, nsweep=a.nsweep, nburn=a.nburn, seed=a.seed)
np.savez(os.path.join(HERE, "data", f"mg_direct_{a.geminal}{a.gamma}_s{a.seed}.npz"), E=E, err=err)
print(f"[direct] {a.geminal}{a.gamma} c={a.c} seed={a.seed} E={E:.6f} +/- {err:.6f}", flush=True)
