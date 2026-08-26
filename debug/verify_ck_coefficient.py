"""
Verification script: is geovac/sturmian_solver.py::_ck_coefficient a genuine bug?
====================================================================================

Flag under test: `SturmianCI._ck_coefficient` (geovac/sturmian_solver.py:214) uses
    q = mc - ma
in the second Wigner 3j argument, instead of the standard Condon-Shortley
    q = ma - mc
(the sign the function's OWN docstring formula implies, and the sign used
elsewhere in the codebase's "Rule B / global-M_L / accuracy" convention,
e.g. geovac/casimir_ci.py::_gaunt_ck).

This script computes c^k(l=1,m=+1; l=1,m=0) [needed for the claimed example
<2p+1 2p-1|1/r12|2p0 2p0>] three independent ways:

  (a1) 3j-formula, CORRECT sign (q = ma - mc), reusing sturmian_solver's own
       _wigner3j (isolates: "is the bug just the sign flip?").
  (a2) geovac/casimir_ci.py::_gaunt_ck -- an independently written 3j-based
       implementation elsewhere in the codebase (different file, different
       author session, its own _wigner3j).
  (b)  Direct 2D numerical quadrature over the sphere of the DEFINING integral
       c^k(l,m,l',m') = sqrt(4pi/(2k+1)) * Integral[ Y*_lm Y_{k,m-m'} Y_l'm' dOmega ]
       -- no 3j symbol anywhere, using scipy.special.sph_harm on a fine grid.

...and compares all three to the framework's actual (a3) sturmian_solver._ck_coefficient.

Then assembles the full ERI <2p+1 2p-1|1/r12|2p0 2p0> using the corrected c^k's
times the module's own (trusted, angular-independent) radial Slater integral
_slater_rk, and compares to what SturmianCI._build_eri actually produces.

Finally, scopes impact:
  - scans the whole s+p (max_n=2) angular table for how many (a,c) pairs with
    ma != mc are silently zeroed;
  - empirically measures the energy shift on He (n_e=2, Z=2) FCI at max_n=2/3
    with the bug vs. the corrected coefficient (monkeypatch within this
    process only -- no production files touched);
  - does the same for a 4-electron toy system at max_n=2 (forces partial p
    occupation), to bound the effect size for a genuinely p-occupied case.

Read-only w.r.t. geovac/ and papers/: this script imports and monkeypatches
in its own process; it edits nothing on disk under geovac/ or papers/.
"""
import sys
import os
import copy
import numpy as np
from scipy.special import sph_harm
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import geovac.sturmian_solver as ss
import geovac.casimir_ci as cc


# ---------------------------------------------------------------------------
# (a1) Corrected 3j-formula c^k, reusing sturmian_solver's OWN _wigner3j
# ---------------------------------------------------------------------------
def ck_correct_3j(la: int, ma: int, lc: int, mc: int, k: int) -> float:
    """Standard Condon-Shortley c^k, q = ma - mc (the sign the module's own
    docstring formula implies: c^k(l m, l' m') ~ Integral[Y*_lm Y_{k,m-m'} Y_l'm']).
    Uses ss._wigner3j (the SAME 3j engine the framework already has), so this
    isolates whether the bug is purely the q-sign flip."""
    q = ma - mc
    pre = ((-1) ** ma * np.sqrt((2 * la + 1) * (2 * lc + 1)))
    w1 = ss._wigner3j(la, k, lc, 0, 0, 0)
    if abs(w1) < 1e-15:
        return 0.0
    w2 = ss._wigner3j(la, k, lc, -ma, q, mc)
    return pre * w1 * w2


# ---------------------------------------------------------------------------
# (b) Direct numerical quadrature of the DEFINING integral -- no 3j at all
# ---------------------------------------------------------------------------
def ck_numeric_sphere(la: int, ma: int, lc: int, mc: int, k: int,
                       n_theta: int = 96, n_phi: int = 192) -> complex:
    """c^k(l,m,l',m') = sqrt(4pi/(2k+1)) * Integral Y*_lm(Omega) Y_{k,q}(Omega) Y_l'm'(Omega) dOmega
    with q = m - m', evaluated by brute-force 2D grid quadrature on the sphere.
    Gauss-Legendre nodes in cos(theta); uniform (rectangle-rule) nodes in phi
    -- exact for the band-limited trig content here (|m| <= 3)."""
    q = ma - mc
    if abs(q) > k:
        return 0.0
    x, w = np.polynomial.legendre.leggauss(n_theta)  # nodes/weights on [-1,1] for cos(theta)
    theta = np.arccos(x)
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    dphi = 2 * np.pi / n_phi

    TH, PH = np.meshgrid(theta, phi, indexing='ij')
    Wtheta = np.broadcast_to(w[:, None], TH.shape)

    Y_la_ma = sph_harm(ma, la, PH, TH)
    Y_k_q = sph_harm(q, k, PH, TH)
    Y_lc_mc = sph_harm(mc, lc, PH, TH)

    integrand = np.conj(Y_la_ma) * Y_k_q * Y_lc_mc
    integral = np.sum(integrand * Wtheta) * dphi
    return np.sqrt(4 * np.pi / (2 * k + 1)) * integral


def main():
    print("=" * 78)
    print("PART 1 -- c^k(l=1,m=+1; l=1,m=0) and c^k(l=1,m=-1; l=1,m=0), k=0,1,2")
    print("=" * 78)
    cases = [(1, 1, 1, 0), (1, -1, 1, 0), (1, 0, 1, 0), (1, 1, 1, -1)]
    for (la, ma, lc, mc) in cases:
        print(f"\n-- c^k(l={la},m={ma}; l={lc},m={mc}) --  [ma {'==' if ma==mc else '!='} mc]")
        for k in range(0, la + lc + 1):
            fw = ss._ck_coefficient(la, ma, lc, mc, k)
            a1 = ck_correct_3j(la, ma, lc, mc, k)
            a2 = cc._gaunt_ck(la, ma, lc, mc, k)
            b = ck_numeric_sphere(la, ma, lc, mc, k)
            b_re, b_im = b.real, b.imag
            print(f"  k={k}: framework(sturmian)={fw: .8f}  "
                  f"correct-3j={a1: .8f}  casimir._gaunt_ck={a2: .8f}  "
                  f"numeric-sphere={b_re: .8f}{'+%.1eJ'%b_im if abs(b_im)>1e-9 else ''}")

    print()
    print("=" * 78)
    print("PART 2 -- the claimed ERI  <2p+1 2p-1 | 1/r12 | 2p0 2p0>")
    print("=" * 78)
    # a = 2p+1, c = 2p0  (electron-1 pair);  b = 2p-1, d = 2p0  (electron-2 pair)
    la, ma = 1, 1
    lc, mc = 1, 0
    lb, mb = 1, -1
    ld, md = 1, 0
    k_scale = 1.0     # Sturmian scale parameter (Z_eff = n * k_scale)
    n_val = 2
    Z_eff = n_val * k_scale

    ks_allowed = [k for k in range(0, la + lc + 1) if (la + lc + k) % 2 == 0]
    print(f"Allowed k (parity + triangle): {ks_allowed}   (Z_eff={Z_eff}, k_scale={k_scale})")

    rk_vals = {}
    for k in ks_allowed:
        rk_vals[k] = ss._slater_rk(n_val, la, Z_eff, n_val, lb, Z_eff,
                                    n_val, lc, Z_eff, n_val, ld, Z_eff, k)
    print("Radial R^k (module's own trusted _slater_rk, angular-independent):")
    for k, v in rk_vals.items():
        print(f"  R^{k} = {v:.8f}")

    # Framework (buggy) assembly
    eri_framework = 0.0
    for k in ks_allowed:
        c_ac = ss._ck_coefficient(la, ma, lc, mc, k)
        c_bd = ss._ck_coefficient(lb, mb, ld, md, k)
        eri_framework += c_ac * c_bd * rk_vals[k]

    # Method (a1): corrected 3j (same _wigner3j engine, sign fixed)
    eri_a1 = 0.0
    for k in ks_allowed:
        c_ac = ck_correct_3j(la, ma, lc, mc, k)
        c_bd = ck_correct_3j(lb, mb, ld, md, k)
        eri_a1 += c_ac * c_bd * rk_vals[k]

    # Method (a2): casimir_ci's independent _gaunt_ck implementation
    eri_a2 = 0.0
    for k in ks_allowed:
        c_ac = cc._gaunt_ck(la, ma, lc, mc, k)
        c_bd = cc._gaunt_ck(lb, mb, ld, md, k)
        eri_a2 += c_ac * c_bd * rk_vals[k]

    # Method (b): pure numeric-sphere Gaunt (no 3j at all)
    eri_b = 0.0
    for k in ks_allowed:
        c_ac = ck_numeric_sphere(la, ma, lc, mc, k).real
        c_bd = ck_numeric_sphere(lb, mb, ld, md, k).real
        eri_b += c_ac * c_bd * rk_vals[k]

    print(f"\n  Framework (SturmianCI._ck_coefficient, as-shipped): {eri_framework:.8f}")
    print(f"  Method a1 (corrected-sign 3j, same wigner3j engine): {eri_a1:.8f}")
    print(f"  Method a2 (casimir_ci._gaunt_ck, independent file):  {eri_a2:.8f}")
    print(f"  Method b  (direct numeric sphere quadrature, no 3j): {eri_b:.8f}")
    print(f"  Task-flagged reference value: approx -0.0342 (unspecified k_scale/Z_eff in the flag)")

    print()
    print("=" * 78)
    print("PART 3 -- does SturmianCI._build_eri actually drop this entry end-to-end?")
    print("=" * 78)
    solver = ss.SturmianCI(Z=2, n_electrons=2, max_n=2)
    states = solver.states
    print("States (max_n=2):", states)
    idx = {s: i for i, s in enumerate(states)}
    a_i, c_i = idx[(2, 1, 1)], idx[(2, 1, 0)]
    b_i, d_i = idx[(2, 1, -1)], idx[(2, 1, 0)]
    eri_table = solver._build_eri(k=k_scale)
    key = (a_i, b_i, c_i, d_i)
    present = key in eri_table
    print(f"  orbital index map: 2p+1={a_i}, 2p-1={b_i}, 2p0={c_i}(=d_i={d_i})")
    print(f"  (a,b,c,d)=(2p+1,2p-1,2p0,2p0) key {key} present in eri dict? {present}"
          f"{'  value=' + str(eri_table[key]) if present else '  (=> exactly 0.0, dropped)'}")

    print()
    print("=" * 78)
    print("PART 4 -- systemic scope: how many m-changing (a,c) pairs are silently zeroed?")
    print("=" * 78)
    n_sp = solver.n_spatial
    total_pairs = 0
    off_diag_pairs = 0
    fw_zero_when_should_be_nonzero = 0
    for a in range(n_sp):
        la_, ma_ = states[a][1], states[a][2]
        for c in range(n_sp):
            lc_, mc_ = states[c][1], states[c][2]
            total_pairs += 1
            if ma_ == mc_:
                continue
            off_diag_pairs += 1
            # does the CORRECT coefficient say this pair is nonzero for some k?
            any_nonzero_correct = False
            any_nonzero_fw = False
            for k in range(0, la_ + lc_ + 1):
                if abs(ck_correct_3j(la_, ma_, lc_, mc_, k)) > 1e-12:
                    any_nonzero_correct = True
                if abs(ss._ck_coefficient(la_, ma_, lc_, mc_, k)) > 1e-12:
                    any_nonzero_fw = True
            if any_nonzero_correct and not any_nonzero_fw:
                fw_zero_when_should_be_nonzero += 1
    print(f"  Total (a,c) orbital pairs in max_n=2 s+p basis: {total_pairs}")
    print(f"  Pairs with m_a != m_c (off-diagonal-in-m):     {off_diag_pairs}")
    print(f"  Off-diagonal pairs where CORRECT c^k is genuinely")
    print(f"  nonzero for some k, but framework returns 0 for ALL k: "
          f"{fw_zero_when_should_be_nonzero} / {off_diag_pairs}")

    print()
    print("=" * 78)
    print("PART 5 -- empirical energy impact: He (Z=2, n_e=2) FCI, max_n=2 and max_n=3")
    print("=" * 78)
    # Monkeypatch the module-level name used inside _build_eri; restore after.
    orig_ck = ss._ck_coefficient
    try:
        for max_n in (2, 3):
            solver_bug = ss.SturmianCI(Z=2, n_electrons=2, max_n=max_n)
            r_bug = solver_bug.solve(k=1.5)
            E_bug = r_bug['energy']

            ss._ck_coefficient = ck_correct_3j
            solver_fix = ss.SturmianCI(Z=2, n_electrons=2, max_n=max_n)
            r_fix = solver_fix.solve(k=1.5)
            E_fix = r_fix['energy']
            ss._ck_coefficient = orig_ck

            dE_mHa = (E_fix - E_bug) * 1000.0
            print(f"  max_n={max_n}: E(buggy c^k)   = {E_bug:.6f} Ha")
            print(f"           E(corrected c^k) = {E_fix:.6f} Ha")
            print(f"           dE = {dE_mHa:+.4f} mHa   "
                  f"(eri_count buggy={r_bug['eri_count']}, fixed={r_fix['eri_count']})")
    finally:
        ss._ck_coefficient = orig_ck

    print()
    print("=" * 78)
    print("PART 6 -- 4-electron toy system at max_n=2 (forces partial p occupation)")
    print("=" * 78)
    orig_ck = ss._ck_coefficient
    try:
        for Z_test in (4, 6):
            solver_bug = ss.SturmianCI(Z=Z_test, n_electrons=4, max_n=2)
            r_bug = solver_bug.solve(k=1.5)
            E_bug = r_bug['energy']

            ss._ck_coefficient = ck_correct_3j
            solver_fix = ss.SturmianCI(Z=Z_test, n_electrons=4, max_n=2)
            r_fix = solver_fix.solve(k=1.5)
            E_fix = r_fix['energy']
            ss._ck_coefficient = orig_ck

            dE_mHa = (E_fix - E_bug) * 1000.0
            print(f"  Z={Z_test}, n_e=4, max_n=2: E(buggy)={E_bug:.6f} Ha  "
                  f"E(fixed)={E_fix:.6f} Ha  dE={dE_mHa:+.4f} mHa")
    finally:
        ss._ck_coefficient = orig_ck

    print()
    print("=" * 78)
    print("PART 7 -- is this the SAME deliberate 'Rule A' convention as composed_qubit.py?")
    print("=" * 78)
    import geovac.composed_qubit as cq
    same_sign_convention = True
    for (la, ma, lc, mc, k) in [(1, 1, 1, 0, 2), (1, -1, 1, 0, 2), (2, 2, 2, 0, 2)]:
        v_sturm = ss._ck_coefficient(la, ma, lc, mc, k)
        v_cq = cq._ck_coefficient(la, ma, lc, mc, k)
        match = abs(v_sturm - v_cq) < 1e-12
        same_sign_convention = same_sign_convention and match
        print(f"  ({la},{ma},{lc},{mc},k={k}): sturmian_solver={v_sturm:.8f}  "
              f"composed_qubit={v_cq:.8f}  bit-match={match}")
    print(f"  sturmian_solver._ck_coefficient uses the IDENTICAL q=mc-ma convention "
          f"as composed_qubit._ck_coefficient: {same_sign_convention}")
    print("  NOTE: composed_qubit._ck_coefficient / lattice_index._ck_coefficient are the")
    print("  criteria.md-registered, disclosed 'Rule A' (pair-diagonal, sparsity-for-the-QC-")
    print("  product) homes. sturmian_solver.SturmianCI/StandardFCI/GeneralizedSturmianCI are")
    print("  NOT in that registered list -- their stated purpose (compare_he/compare_he_")
    print("  generalized: benchmarking energy ACCURACY against NIST) is exactly the 'B is the")
    print("  point' case per criteria.md's own Dual-rule ERI framing.")


if __name__ == "__main__":
    main()
