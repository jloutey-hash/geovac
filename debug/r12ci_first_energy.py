"""First H2 energy from the exact mpf r12 engine (debug/prolate_r12_mpf.py).

Builds a mixed p={0,1} prolate-spheroidal James-Coolidge basis, assembles the
exact (mpf) S and H = T + V_ne + V_ee, and solves the generalized eigenproblem.
Small basis: direct scipy eigh.  Large basis: canonical orthogonalization
(re-basing) against S, mirroring Paper 12's prolate_recondition route.

Compares to E_EXACT = -1.174475 Ha (Kolos-Wolniewicz), D_e = 0.174475 Ha.
Run from the project root:  python debug/r12ci_first_energy.py [j_max l_max]
"""
import os
import sys
import numpy as np
import scipy.linalg as sla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geovac import prolate_recondition as pr
import prolate_r12_mpf as m

R = 1.4011
ALPHA = 1.0
E_EXACT = -1.174475
DE_EXACT = 0.174475
NUC = 1.0 / R


def build_basis(j_max, l_max, p_set=(0, 1)):
    """Gerade sigma product basis (mu=0), (l+m) even, x p in p_set."""
    rlist = [(j, k) for j in range(j_max + 1) for k in range(j_max + 1)]
    alist = [(l, mm) for l in range(l_max + 1) for mm in range(l_max + 1)
             if (l + mm) % 2 == 0]
    basis = []
    for p in p_set:
        for (j, k) in rlist:
            for (l, mm) in alist:
                basis.append((pr.ProductFn(j, l, k, mm, 0, ALPHA), p))
    return basis


def solve_direct(S, H):
    """Direct generalized eigensolve (small, well-conditioned basis)."""
    w = sla.eigvalsh(H, S)
    return float(w[0])


def solve_canonical(S, H, thresholds=(1e-14, 1e-13, 1e-12, 1e-11, 1e-10)):
    """Canonical orthogonalization sweep: drop eigvecs of S below tol, solve in
    the retained subspace, return the LOWEST variational energy across the sweep."""
    sval, svec = np.linalg.eigh(S)
    best = None
    for tol in thresholds:
        keep = sval > tol * sval[-1]
        if keep.sum() == 0:
            continue
        X = svec[:, keep] / np.sqrt(sval[keep])
        Ho = X.T @ H @ X
        w = np.linalg.eigvalsh(Ho)
        e = float(w[0])
        if best is None or e < best[0]:
            best = (e, tol, int(keep.sum()))
    return best


def _report(tag, e, n, extra=""):
    et = e + NUC
    de = (-1.0 - et) / DE_EXACT * 100
    print(f"{tag}: E_elec={e:.7f}  E_tot={et:.7f}  "
          f"D_e%={de:.3f}  err={(E_EXACT-et)*1000:+.4f} mHa {extra}")


def run(j_max, l_max, use_mpf=False):
    basis = build_basis(j_max, l_max)
    n = len(basis)
    print(f"basis: j_max={j_max} l_max={l_max}  n={n}  (p=0 and p=1)")
    if use_mpf:
        Smp, Hmp = m.assemble_mixed(basis, R, ALPHA, mpf_out=True)
        Sf = m._tofloat(Smp)
        print(f"cond(S) = {np.linalg.cond(Sf):.3e}")
        bf = solve_canonical(Sf, m._tofloat(Hmp))
        if bf:
            _report("float64  ", bf[0], n, f"(tol={bf[1]:g}, kept {bf[2]}/{n})")
        bm = m.solve_canonical_mpf(Smp, Hmp)
        if bm:
            _report("mpf      ", bm[0], n, f"(tol={bm[1]:g}, kept {bm[2]}/{n})")
        return
    S, H = m.assemble_mixed(basis, R, ALPHA)
    condS = np.linalg.cond(S)
    print(f"cond(S) = {condS:.3e}")
    try:
        e_dir = solve_direct(S, H)
        _report("direct   ", e_dir, n)
    except Exception as ex:
        print("direct   : FAILED", ex)
    best = solve_canonical(S, H)
    if best:
        _report("canonical", best[0], n, f"(tol={best[1]:g}, kept {best[2]}/{n})")


def run_p0_only(j_max, l_max):
    """Control: SAME basis with p=0 only (no r12) -> shows the r12 gain."""
    basis = build_basis(j_max, l_max, p_set=(0,))
    n = len(basis)
    S, H = m.assemble_mixed(basis, R, ALPHA)
    best = solve_canonical(S, H)
    e, tol, k = best
    et = e + NUC
    de = (-1.0 - et) / DE_EXACT * 100
    print(f"[p0-only control] n={n}  E_tot={et:.6f}  D_e%={de:.3f}  "
          f"err={(E_EXACT-et)*1000:+.3f} mHa")


if __name__ == "__main__":
    jm = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    lm = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    use_mpf = "mpf" in sys.argv[3:]
    run_p0_only(jm, lm)
    run(jm, lm, use_mpf=use_mpf)
