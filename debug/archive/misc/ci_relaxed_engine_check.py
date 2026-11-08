r"""Grid-vs-engine check of the CI block used by geovac/lih_r12ci/ci_relaxed.py.

The relaxed eigenproblem uses H_CI on the marriage's grid quadrature (so its frozen restriction is the
banked 2x2 exactly).  This compares it with the exact-integral engine H_CI (h1_no / eri_no via
Slater-Condon) in the same chi_I = dA dB / 2 basis, and re-solves ansatz (A) with the engine CI block
swapped in (a quadrature-sensitivity diagnostic of the relaxed energy).
Usage: python debug/ci_relaxed_engine_check.py 1e-2
"""
import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def main() -> int:
    tr = sys.argv[1] if len(sys.argv) > 1 else "1e-2"
    import geovac.lih_r12ci.ref_tensor as P0
    from geovac.prolate_allelectron.fci_grid import _dets, _matel
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.pi_m16 as M
    import scipy.linalg as sl
    KG.USE_EXACT_NEUMANN = True
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    _, info = M.build_ref_from_prims(ref, "chk", cfg_tag='pi', trunc=float(tr), verbose=False)
    Mk = int(info['Mk'])
    dets = _dets(Mk, 2, 2)
    _, _, sigma = P0.fold_cs(info['ct'], Mk)
    sig = sigma.reshape(-1)
    kept = np.asarray(info['kept'])
    CS, _, _ = P0.fold_cs(info['ct'], Mk)
    ia, ib = np.nonzero(CS)
    n_ap = CS.shape[1]
    order = ia * n_ap + ib                                   # chi_I order of ci_relaxed
    assert set(order.tolist()) == set(kept.tolist())
    N = len(order)
    He = np.empty((N, N))
    for I in range(N):
        for J in range(N):
            He[I, J] = sig[order[I]] * sig[order[J]] * _matel(info['h1_no'], info['eri_no'],
                                                              dets[order[I]], dets[order[J]])
    z = np.load(os.path.join(ROOT, "debug", "data", f"ci_relaxed_{tr}.npz"), allow_pickle=True)
    Hg = z['H_CI']; d = z['d_ref']
    print(f"N={N}  max|H_grid - H_engine| = {np.max(np.abs(Hg - He)):.2e}")
    print(f"d.Hg.d + Vnn = {d @ Hg @ d + M.V_NN:.9f}   d.He.d + Vnn = {d @ He @ d + M.V_NN:.9f} "
          f"(E_trunc {info['E_trunc']:.9f})")
    for I in range(min(N, 12)):
        print(f"  I={I:3d}  (H_grid d)_I={Hg[I] @ d:+.6f}  (H_engine d)_I={He[I] @ d:+.6f}  "
              f"diff {(Hg[I] - He[I]) @ d:+.2e}")
    np.savez(os.path.join(ROOT, "debug", "data", f"ci_relaxed_engine_{tr}.npz"), H_engine=He, d=d)
    H = z['H'].copy(); S = z['S']
    H[:N, :N] = He
    w = sl.eigh(H, S, eigvals_only=True)
    print(f"E(A) grid CI block = {float(z['E_A']):.9f};  E(A) engine CI block = {w[0] + M.V_NN:.9f}")
    print(f"CI-only lowest: grid {np.linalg.eigvalsh(Hg)[0] + M.V_NN:.9f}  engine "
          f"{np.linalg.eigvalsh(He)[0] + M.V_NN:.9f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
