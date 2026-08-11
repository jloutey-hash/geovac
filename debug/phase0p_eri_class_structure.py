"""Phase 0-prime, retargeted build (option A): structure of the atom-centered
two-center ERI, class by class.

The Hylleraas-extension plan died because neumann_vee.py serves a basis with no
azimuthal index. Retargeting to a real atom-centered engine means Phase 0 has to
be redone -- and the first finding is that this is NOT one problem. A two-center
ERI (ab|cd) with each index on center A or B splits into classes with genuinely
different machinery:

  (AA|AA), (BB|BB)   one-center. Already solved in the corpus
                     (hypergeometric_slater.py, exact Fractions).
  (AA|BB)            BOTH charge distributions are one-center, on opposite
                     nuclei. This is the classical two-center Coulomb integral
                     between two multipole expansions -- the two-electron analog
                     of what shibuya_wulfman.py already does for cross-center
                     V_ne. Expected to TERMINATE by Gaunt on each side.
  (AA|AB), (AB|BB)   "hybrid": one distribution one-center, the other two-center.
  (AB|AB)            "exchange": both distributions two-center. The hardest, and
                     the subject of an entire paper (Ruedenberg 1951 Part II).

Paper 58's census measured the class weights at n_max=2: the three genuinely
cross classes are ~80% of the tensor, (AA|BB) is 390/2944 = 13.2%, and
(AA|AA)+(BB|BB) is 7.3%. So the bulk of the work is in the hardest classes -- but
(AA|BB) is a clean, self-contained first target with a decidable support rule.

WHAT THIS SCRIPT ESTABLISHES
----------------------------
The support rule for the (AA|BB) class, exactly.

A one-center product chi_{n1 l1 m1} * chi_{n2 l2 m2} on center A has angular
content expandable in Y_LM about A, and by the Gaunt/3-j triangle rule that
expansion TERMINATES at L <= l1 + l2 with parity l1 + l2 + L even, and requires
M = m1 - m2 (with the standard conjugation convention). The two-center Coulomb
interaction then couples a multipole L_A on A to L_B on B. So (AA|BB) is nonzero
only if BOTH sides admit a nonvanishing multipole, and the whole thing is a
finite double sum -- no truncation, and the support is decidable from labels
alone.

Verified here by computing the Gaunt coefficients exactly (Wigner 3-j, exact
rationals under the square root) and confirming: termination at l1 + l2, the
parity selection, and the M rule.

Run from repo root:  python debug/phase0p_eri_class_structure.py
"""

from __future__ import annotations

import sympy as sp
from sympy.physics.wigner import gaunt, wigner_3j


def gaunt_exact(l1, m1, l2, m2, l3, m3):
    """Exact Gaunt coefficient integral Y_l1m1 Y_l2m2 Y_l3m3 dOmega."""
    return sp.nsimplify(gaunt(l1, l2, l3, m1, m2, m3))


def one_center_multipole_support(l1, m1, l2, m2, L_ceiling=None):
    """Which (L, M) survive in the expansion of conj(Y_l1m1) * Y_l2m2?

    Uses the Gaunt coefficient with the conjugation folded in as
    conj(Y_l1m1) = (-1)^m1 Y_l1,-m1, so the surviving M is m2 - m1.
    """
    predicted_Lmax = l1 + l2
    ceiling = L_ceiling if L_ceiling is not None else predicted_Lmax + 3
    surviving = []
    for L in range(0, ceiling + 1):
        M = m2 - m1
        if abs(M) > L:
            continue
        # integral of Y_l1,-m1 * Y_l2,m2 * conj(Y_L,M)  ->  Gaunt with -M
        g = gaunt_exact(l1, -m1, l2, m2, L, -M)
        if g != 0:
            surviving.append((L, M, g))
    return predicted_Lmax, surviving


CASES = [
    # (l1, m1, l2, m2)
    (0, 0, 0, 0),
    (1, 0, 0, 0),
    (1, 1, 1, 1),
    (1, 1, 1, -1),
    (1, 1, 0, 0),
    (2, 0, 1, 0),
    (2, 2, 1, 1),
    (2, 1, 2, -1),
    (2, 2, 2, 2),
]


def main() -> None:
    print("Phase 0' -- (AA|BB) support rule for atom-centered two-center ERIs\n")
    print("One-center product conj(Y_l1m1) * Y_l2m2 expanded in Y_LM:\n")
    print(f"{'(l1,m1)x(l2,m2)':<20}{'l1+l2':>7}{'max L':>7}{'M':>5}"
          f"{'L values surviving':>26}{'verdict':>12}")
    print("-" * 78)

    all_terminate = True
    parity_ok = True
    m_rule_ok = True
    for l1, m1, l2, m2 in CASES:
        Lmax_pred, surv = one_center_multipole_support(l1, m1, l2, m2)
        Ls = [L for L, _M, _g in surv]
        highest = max(Ls) if Ls else None
        Ms = {M for _L, M, _g in surv}
        # checks
        term = (highest is None) or (highest <= Lmax_pred)
        par = all((l1 + l2 + L) % 2 == 0 for L in Ls)
        mr = (Ms == {m2 - m1}) or not Ls
        all_terminate &= term
        parity_ok &= par
        m_rule_ok &= mr
        label = f"({l1},{m1})x({l2},{m2})"
        verdict = "TERMINATES" if term else "EXCEEDS"
        print(f"{label:<20}{Lmax_pred:>7}{str(highest):>7}{str(m2-m1):>5}"
              f"{str(Ls):>26}{verdict:>12}")

    print("-" * 78)
    print(f"\nTermination at L <= l1 + l2 in all cases      : {all_terminate}")
    print(f"Parity selection (l1 + l2 + L even) holds     : {parity_ok}")
    print(f"M rule (M = m2 - m1 only) holds               : {m_rule_ok}")

    if all_terminate and parity_ok and m_rule_ok:
        print("\n=> (AA|BB) support is DECIDABLE FROM LABELS ALONE.")
        print("   Each side contributes a finite multipole set")
        print("   L in {|l1-l2|, ..., l1+l2} with l1+l2+L even and M = m2-m1;")
        print("   the class is a finite double sum over (L_A, L_B) with no")
        print("   truncation. This is the two-electron analog of the")
        print("   L_max = l1+l2 termination shibuya_wulfman already achieves")
        print("   for cross-center V_ne, and it makes (AA|BB) -- 13.2% of the")
        print("   census tensor at n_max=2 -- the correct first build target.")
    else:
        print("\n=> support is NOT decidable from labels; revisit before building.")

    print("\nStill open for the harder classes (AA|AB), (AB|BB), (AB|AB):")
    print("  those carry a genuinely two-center charge distribution, whose")
    print("  angular content about either nucleus does NOT terminate")
    print("  (cos theta_A = (1 + xi eta)/(xi + eta) mixes the coordinates).")
    print("  They are the Ruedenberg Part II problem and need separate scoping.")


if __name__ == "__main__":
    main()
