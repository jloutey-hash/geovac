r"""PHASE C -- GATE G2-pi: independent VMC cross-check of the M=16 pi marriage energy.

Standalone driver (launch AFTER debug/lih_pi_m16_energy.py has written its .npz):
    python debug/lih_pi_m16_vmc.py [--geminal linexp|exp] [--nwalk 2000 --nsweep 4000]

Reconstructs the SAME truncated M=16 pi natural-orbital reference the marriage used (so gate-6,
c=0, reproduces E_trunc), then evaluates E_VMC[(1 + c(F - Fbar)) Psi_CI] with the STANDARD local
energy (NOT the gradient form; lih_vmc.local_energy_linear) reading c, Fbar from the energy .npz.
G2-pi PASS: |E_VMC - E_R12| < 1 mHa AND gate-6 reproduces E_trunc within stats.
Writes debug/data/lih_pi_m16_vmc.log.
"""
from __future__ import annotations

import os
import sys
import argparse
from itertools import combinations

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, "data")

import lih_pi_m16 as M
import lih_marriage_phase0 as P0
import lih_vmc
ENERGY_NPZ = os.path.join(DATA, "lih_pi_m16_energy.npz")


def _build_wf_from_no(ct: np.ndarray, T_no: np.ndarray, Mk: int, E_trunc: float):
    """LiHWavefunction on the TRUNCATED M=16 pi natural-orbital reference (matches the marriage)."""
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    prims = ref['prims']
    na = nb = M.NELEC // 2
    apairs = list(combinations(range(Mk), na)); bpairs = list(combinations(range(Mk), nb))
    CS, _, _ = P0.fold_cs(ct, Mk)                       # grouped-order signed coefficient matrix
    Vnn = ref['Vnn']
    return lih_vmc.LiHWavefunction(prims, T_no, ct, apairs, bpairs, CS, ref['R'], 3.0, 1.0, Vnn,
                                   E_trunc - Vnn, Mk, na, nb)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--geminal', default='linexp', choices=['linexp', 'exp'])
    ap.add_argument('--nwalk', type=int, default=2000)
    ap.add_argument('--nsweep', type=int, default=4000)
    ap.add_argument('--nburn', type=int, default=800)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()

    z = np.load(ENERGY_NPZ, allow_pickle=True)
    g = args.geminal
    if f"{g}_E_R12" not in z.files:
        raise SystemExit(f"{g} not in {ENERGY_NPZ}; run lih_pi_m16_energy.py --geminal {g}")
    E_R12 = float(z[f"{g}_E_R12"]); c = float(z[f"{g}_c"]); Fbar = float(z[f"{g}_Fbar"])
    E_trunc = float(z['E_trunc']); Mk = int(z['Mk']); gamma = float(z['GAM'])
    T_no = z['T_no']; ct = z['ct']
    print("=" * 90)
    print(f"GATE G2-pi -- VMC cross-check of the M=16 pi marriage ({g})")
    print(f"  E_R12 = {E_R12:.6f}  c = {c:+.5f}  Fbar = {Fbar:.5f}  gamma = {gamma}  E_trunc = {E_trunc:.6f}")
    print("=" * 90, flush=True)

    wf = _build_wf_from_no(ct, T_no, Mk, E_trunc)

    # gate-6: J=1 VMC must reproduce E_trunc
    E0, e0, a0 = lih_vmc.vmc(wf, nwalk=args.nwalk, nsweep=args.nsweep, nburn=args.nburn,
                             jas=None, seed=args.seed)
    dev0 = (E0 - E_trunc) / max(e0, 1e-9)
    print(f"  gate-6 (c=0): E_VMC = {E0:.5f} +/- {e0:.5f}  vs E_trunc {E_trunc:.6f}  "
          f"dev = {dev0:+.1f} sigma  {'PASS' if abs(dev0) < 4 else 'FAIL'}", flush=True)

    # correlated: (1 + c(F-Fbar)) Psi_CI
    linjas = lih_vmc.LinearGeminalJastrow(g, gamma, c, Fbar, wf.na, wf.nb)
    E, err, acc = lih_vmc.vmc_linear(wf, linjas, nwalk=args.nwalk, nsweep=args.nsweep,
                                     nburn=args.nburn, seed=args.seed)
    dER = (E - E_R12) * 1e3
    ok = abs(dER) < 1.0 and abs(dev0) < 4
    print(f"  E_VMC(linear) = {E:.6f} +/- {err:.6f}  vs E_R12 {E_R12:.6f}  d = {dER:+.3f} mHa  "
          f"[{'PASS' if abs(dER) < 1.0 else 'FAIL'}]")
    print(f"\n  GATE G2-pi: {'PASS' if ok else 'FAIL'}", flush=True)


if __name__ == "__main__":
    main()
