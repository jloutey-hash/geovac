"""aha Track 1 / STEP 1 -- the two-functionals-of-one-spectrum THEOREM check.

Claim under test
----------------
For a two-center basis that is orthonormal WITHIN each center, the joint metric is
    S = [[I, C], [C^T, I]],           C = cross-center overlap block,
so with the SVD C = U diag(sigma) V^T,
    (i)   spec(S) = {1 + sigma_k} U {1 - sigma_k}
    (ii)  cond(S) = (1 + sigma_max) / (1 - sigma_max)
    (iii) if C is symmetric, the gerade/ungerade blocks are I +/- C, so
          cond(gerade) = (1+lam_max)/(1+lam_min),  cond(ungerade) = (1-lam_min)/(1-lam_max)
          with lam = EIGENvalues of C  (signed);  |lam|_max = sigma_max.
    (iv)  ||[P_A,P_B]|| = max_k sigma_k sqrt(1-sigma_k^2)   (v4.73.0 convention:
          sigma_k = cos(principal angle theta_k), so this is max_k |sin(2 theta_k)|/2).

So the Paper-60 metric-conditioning wall and the v4.73.0 composition wall are two
EXACT functionals of one object.  Both are checked numerically here on H2+ SW data.

Run:  python debug/aha_t1_theorem.py
"""
from __future__ import annotations

import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from debug.aha_t1_core import (sw_block_trapz, sw_block_gl, assemble_two_center,
                               sigma_spectrum, cond_from_sigma, commutator_from_sigma,
                               commutator_direct)

np.set_printoptions(precision=6, suppress=True)


def main() -> None:
    print("=" * 86)
    print("STEP 1 -- theorem check on H2+ Shibuya-Wulfman data (k=1, R = 2 bohr => s = kR = 2)")
    print("=" * 86)

    # ---------------------------------------------------------------- (0) intra = I?
    print("\n[0] Are the SW intra-center blocks EXACTLY the identity?")
    print("    SW intra = (2/pi) int_0^pi sin(a chi) sin(b chi) dchi = delta_ab  [ANALYTIC]")
    print("    -> identity holds in exact arithmetic; the deviation below is pure quadrature.")
    for nmax in (4, 8, 10):
        Iq = sw_block_trapz(0.0, nmax, "S")
        Ig = sw_block_gl(0.0, nmax, "S")
        print(f"    nmax={nmax:>2}:  max|intra_trapz - I| = {np.abs(Iq-np.eye(nmax)).max():.3e}"
              f"   max|intra_GL - I| = {np.abs(Ig-np.eye(nmax)).max():.3e}")
    print("    => intra blocks are the identity to quadrature accuracy; NO orthonormalization")
    print("       step is needed.  All S below are built with the EXACT identity intra-block.")

    # ------------------------------------------- (0b) independent quadrature cross-check
    print("\n[0b] Cross-center block C: two INDEPENDENT quadratures agree")
    print("     (uniform trapezoid in chi  vs  panelled Gauss-Legendre in v = cot(chi/2))")
    for s in (1.4, 2.0, 4.0, 10.0):
        Ct = sw_block_trapz(s, 10, "S")
        Cg = sw_block_gl(s, 10, "S")
        print(f"     s=kR={s:>5.1f}:  max|C_trapz - C_GL| = {np.abs(Ct-Cg).max():.3e}"
              f"   (max|C| = {np.abs(Ct).max():.4f})")

    # --------------------------------------------------------- (i)-(ii) spectrum + cond
    print("\n[i,ii] spec(S) = 1 +/- sigma_k   and   cond(S) = (1+sigma_max)/(1-sigma_max)")
    print(f"  {'nmax':>4} {'N':>3} | {'cond(S) direct':>15} {'(1+sm)/(1-sm)':>15}"
          f" {'rel.resid':>11} | {'max|spec-(1+-sig)|':>19}")
    print("  " + "-" * 80)
    worst_cond, worst_spec = 0.0, 0.0
    for nmax in range(2, 11):
        C = sw_block_trapz(2.0, nmax, "S")
        S = assemble_two_center(C)
        sig = sigma_spectrum(C)
        c_direct = float(np.linalg.cond(S))
        c_theory = cond_from_sigma(sig)
        rel = abs(c_direct - c_theory) / c_theory
        ev = np.sort(np.linalg.eigvalsh(S))
        pred = np.sort(np.concatenate([1 + sig, 1 - sig]))
        dspec = float(np.abs(ev - pred).max())
        worst_cond = max(worst_cond, rel)
        worst_spec = max(worst_spec, dspec)
        print(f"  {nmax:>4} {2*nmax:>3} | {c_direct:>15.9f} {c_theory:>15.9f}"
              f" {rel:>11.2e} | {dspec:>19.2e}")
    print(f"  WORST relative cond residual = {worst_cond:.2e}   "
          f"WORST spectrum residual = {worst_spec:.2e}")

    # ------------------------------------------------------ (iii) gerade / ungerade
    print("\n[iii] gerade/ungerade blocks are I +/- C  =>  the g/u lever is a SIGN statement")
    print("      about the EIGENvalues of C:  the near-linear-dependent direction has")
    print("      lam ~= +1, so it lands in (I - C) = UNGERADE.  Gerade never sees it.")
    print(f"  {'nmax':>4} | {'lam_min(C)':>11} {'lam_max(C)':>11} {'sigma_max':>10}"
          f" | {'cond(g) dir':>11} {'(1+lmax)/(1+lmin)':>18} | {'cond(u) dir':>11}"
          f" {'(1-lmin)/(1-lmax)':>18}")
    print("  " + "-" * 104)
    for nmax in (2, 4, 6, 8, 10):
        C = sw_block_trapz(2.0, nmax, "S")
        lam = np.linalg.eigvalsh(C)
        g, u = np.eye(nmax) + C, np.eye(nmax) - C
        cg_d, cu_d = float(np.linalg.cond(g)), float(np.linalg.cond(u))
        cg_t = (1 + lam.max()) / (1 + lam.min())
        cu_t = (1 - lam.min()) / (1 - lam.max())
        print(f"  {nmax:>4} | {lam.min():>11.6f} {lam.max():>11.6f} {np.abs(lam).max():>10.6f}"
              f" | {cg_d:>11.5f} {cg_t:>18.5f} | {cu_d:>11.3f} {cu_t:>18.3f}")

    # --------------------------------------------------------------- (iv) commutator
    print("\n[iv] ||[P_A,P_B]|| = max_k sigma_k sqrt(1-sigma_k^2)  (explicit projectors vs formula)")
    print(f"  {'nmax':>4} {'s=kR':>5} | {'direct ||[PA,PB]||':>19} {'formula':>12} {'resid':>10}"
          f" | {'argmax sigma_k':>14}")
    print("  " + "-" * 78)
    worst_comm = 0.0
    for s in (1.4, 2.0, 4.0):
        for nmax in (2, 4, 6, 8, 10):
            C = sw_block_trapz(s, nmax, "S")
            sig = sigma_spectrum(C)
            d = commutator_direct(C)
            f = commutator_from_sigma(sig)
            worst_comm = max(worst_comm, abs(d - f))
            karg = int(np.argmax(sig * np.sqrt(np.maximum(0, 1 - sig ** 2))))
            print(f"  {nmax:>4} {s:>5.1f} | {d:>19.15f} {f:>12.9f} {abs(d-f):>10.2e}"
                  f" | {sig[karg]:>14.6f}")
    print(f"  WORST commutator residual = {worst_comm:.2e}")

    # ---------------------------------- (v) does the ACTUAL driver S (quadrature intra) obey it?
    print("\n[v] Same test on the assembled matrix the Paper-60 driver actually builds")
    print("    (intra = quadrature block, NOT exact I): residual = the intra quadrature error.")
    for nmax in (4, 8, 10):
        intra = sw_block_trapz(0.0, nmax, "S")
        C = sw_block_trapz(2.0, nmax, "S")
        S_drv = np.block([[intra, C], [C.T, intra]])
        c_direct = float(np.linalg.cond(S_drv))
        c_theory = cond_from_sigma(sigma_spectrum(C))
        print(f"    nmax={nmax:>2}: cond(driver S) = {c_direct:.9f}   theory = {c_theory:.9f}"
              f"   rel = {abs(c_direct-c_theory)/c_theory:.2e}"
              f"   (max|intra-I| = {np.abs(intra-np.eye(nmax)).max():.2e})")

    print("\nDONE.")


if __name__ == "__main__":
    main()
