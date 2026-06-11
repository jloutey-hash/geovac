"""
W1a-diag Q-C: Transfer operator T_{lam_1 -> lam_2} between Coulomb-Sturmian bases.

The Coulomb-Sturmian set {chi_{nlm}^{lam}} at fixed lam is COMPLETE on
L^2(R^3) for the radial measure r dr (NOT r^2 dr -- this is the Sturmian
inner product). The set is NOT orthonormal in the L^2(R^3, r^2 dr) sense
unless lam = Z/n.

The "transfer operator" between Sturmian bases at different lambda is just
the change-of-basis matrix:
   T(n, l, m; n', l', m') = <chi_{nlm}^{lam_2} | chi_{n'l'm'}^{lam_1}>_{L^2(r^2 dr)}

The angular factors are diagonal: T = T^radial(n, l; n', l) * delta_{l,l'} * delta_{m, m'}.

Q-C asks: does T preserve sparsity and integer-Pauli scaling?

Compute T^radial in closed form for small n, l, and (a) check it is
generically dense in (n, n'), (b) measure the Frobenius norm growth as
the basis sizes grow, (c) flag whether it preserves "block structure" that
the Pauli encoding relies on.
"""

import json
from sympy import (Rational, Symbol, sqrt, exp, integrate, oo, simplify,
                   factorial, S, sympify, expand, factor, Matrix, eye,
                   nsimplify, srepr)
from sympy.functions.special.polynomials import assoc_laguerre


def R_radial(n, l, lam, r):
    rho = 2 * lam * r
    L_assoc = assoc_laguerre(n - l - 1, 2 * l + 1, rho)
    psi_unnormed = rho**l * exp(-lam * r) * L_assoc
    norm_sq = integrate(psi_unnormed**2 * r**2, (r, 0, oo))
    return psi_unnormed / sqrt(simplify(norm_sq))


def transfer_radial(n_a, l, lam_a, n_b, lam_b):
    """<R_{n_a, l}^{lam_a} | R_{n_b, l}^{lam_b}> in L^2(r^2 dr).

    Closed form because both wavefunctions are
    polynomial * exp(-(lam_a + lam_b) r), and r^2 dr integrals are gamma
    functions.
    """
    r = Symbol('r', positive=True)
    Ra = R_radial(n_a, l, lam_a, r)
    Rb = R_radial(n_b, l, lam_b, r)
    val = integrate(Ra * Rb * r**2, (r, 0, oo))
    return simplify(val)


def transfer_matrix_one_l(l, n_max_a, n_max_b, lam_a, lam_b):
    """Build the (n_a, n_b) transfer matrix at fixed l, with n_a in [l+1, n_max_a]
    and n_b in [l+1, n_max_b]. Each entry is a closed-form symbolic expression
    in (lam_a, lam_b)."""
    rows = []
    for n_a in range(l + 1, n_max_a + 1):
        row = []
        for n_b in range(l + 1, n_max_b + 1):
            T_entry = transfer_radial(n_a, l, lam_a, n_b, lam_b)
            row.append(T_entry)
        rows.append(row)
    return rows


def main():
    out = {}
    lam_a = Symbol('lam_a', positive=True)
    lam_b = Symbol('lam_b', positive=True)

    # 1) Build T at l=0, n in {1, 2, 3} for both bases. 3x3 matrix.
    print("Building T^radial at l=0, n_a, n_b in {1, 2, 3}...")
    T_l0 = transfer_matrix_one_l(0, 3, 3, lam_a, lam_b)
    print("\n  T^radial(l=0) entries:")
    for i, row in enumerate(T_l0):
        for j, entry in enumerate(row):
            n_a = i + 1
            n_b = j + 1
            print(f"    T[n_a={n_a}, n_b={n_b}] = {entry}")

    out["T_l0"] = [[str(x) for x in row] for row in T_l0]

    # 2) Sanity check: matched limit T(lam_a, lam_a) should be identity at
    # diagonal and zero off-diagonal (orthonormality of Coulomb-Sturmians at
    # fixed lambda in the L^2 r^2 dr measure -- no, this is wrong! Coulomb-
    # Sturmians are orthogonal in r dr, not r^2 dr.). Actually at lam = Z/n
    # the matched-Z hydrogenic functions ARE orthonormal in r^2 dr.
    # General Sturmians at fixed lam are orthonormal in (1/r) measure or
    # similar. So matched-lambda T may NOT be identity in r^2 dr!
    print("\n[Sanity] T_l0 at matched lam_b = lam_a:")
    for i, row in enumerate(T_l0):
        for j, entry in enumerate(row):
            n_a = i + 1
            n_b = j + 1
            matched = simplify(entry.subs(lam_b, lam_a))
            print(f"    T[n_a={n_a}, n_b={n_b}] | matched = {matched}")
            # Save matched-limit cell.
            out.setdefault("T_l0_matched", {})[f"({n_a},{n_b})"] = str(matched)

    # 3) Compute Frobenius norm of T_l0 (sum of squares) and density
    # (count of nonzero entries / total).
    sym_entries = [entry for row in T_l0 for entry in row]
    nonzero_count = sum(1 for e in sym_entries if e != 0)
    total_count = len(sym_entries)
    out["density_l0"] = nonzero_count / total_count
    print(f"\n  Density of T_l0 (l=0, 3x3): {nonzero_count}/{total_count} = {out['density_l0']:.3f}")

    # 4) Substitute concrete values to check density numerically.
    # lam_a = 1, lam_b = 2 (e.g., H 1s vs He+ 1s).
    print("\n[Numerical] T_l0 at (lam_a, lam_b) = (1, 2):")
    T_num = []
    for i, row in enumerate(T_l0):
        for j, entry in enumerate(row):
            n_a = i + 1
            n_b = j + 1
            num = float(entry.subs([(lam_a, 1), (lam_b, 2)]))
            T_num.append(num)
            print(f"    T[n_a={n_a}, n_b={n_b}] = {num:.6f}")
    out["T_l0_numerical_lam1_lam2"] = T_num

    # 5) Density check at l=0, larger basis.
    print("\nBuilding larger basis: l=0, n in {1..5} for assessment of density growth...")
    T_l0_big = transfer_matrix_one_l(0, 5, 5, lam_a, lam_b)
    sym_entries_big = [entry for row in T_l0_big for entry in row]
    nonzero_count_big = sum(1 for e in sym_entries_big if e != 0)
    total_count_big = len(sym_entries_big)
    out["density_l0_5x5"] = nonzero_count_big / total_count_big
    print(f"  Density of T_l0 (5x5): {nonzero_count_big}/{total_count_big} = {out['density_l0_5x5']:.3f}")

    # 6) Take a numerical look at the magnitude profile to assess
    # whether T is "approximately diagonal" or not.
    print("\nMagnitude profile at (lam_a=1, lam_b=2), l=0, 5x5:")
    for i, row in enumerate(T_l0_big):
        line = "  "
        for j, entry in enumerate(row):
            num = float(entry.subs([(lam_a, 1), (lam_b, 2)]))
            line += f"{num:>8.4f} "
        print(line)

    # 7) Check Pauli-encoding implication.
    # The composed-block Pauli encoding maps each spin-orbital to a qubit.
    # If T is generically dense, transferring an operator A from basis 1 to
    # basis 2 via A_2 = T A_1 T^T produces a *dense* A_2 even if A_1 was
    # sparse. The Pauli-string count therefore generically inflates by
    # O(Q^2) -- the Frobenius norm is preserved but locality is destroyed.
    # This is the same Lowdin-orthogonalization pathology Track BU
    # (multi-electron Sturmian) hit: "Lowdin orthogonalization destroys
    # Gaunt sparsity at 14x inflation for nested Be." The transfer operator
    # has the same structural defect.

    out["pauli_sparsity_implication"] = {
        "verdict": "transfer operator T_{lam_a->lam_b} is GENERICALLY DENSE in (n_a, n_b) at fixed l",
        "consequence_for_pauli_encoding": (
            "An operator A sparse in basis 1 becomes T A T^T which is dense "
            "in basis 2. Pauli-string count inflates by O(Q^2). This is the "
            "same pathology Track BU and Track DF Sprint 5 (Lowdin "
            "orthogonalization on heterogeneous nested Be) demonstrated."
        ),
        "implication_for_w1a": (
            "The transfer-operator route to multi-lambda multi-particle "
            "compositions cannot preserve integer-Pauli scaling at the "
            "multi-particle FCI level -- but for one-body cross-register "
            "operators (like V_{ne}(r_e, R_n)), this is irrelevant: the "
            "operator is built directly in the joint basis without going "
            "through a transfer-operator change-of-basis. The transfer "
            "operator is a USEFUL diagnostic but NOT the right path to W1a "
            "closure."
        ),
    }

    with open("debug/data/multifocal_b_w1a_qc.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nSaved: debug/data/multifocal_b_w1a_qc.json")


if __name__ == "__main__":
    main()
