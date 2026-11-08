"""Master Mellin engine on finite spectral triples — A_F probe.

Tests whether 𝓜[Tr(D_F^k · e^{-tD_F²})] for k ∈ {0,1,2} produces
M1/M2/M3-class transcendentals on:
  (a) SM electroweak A_F = ℂ ⊕ ℍ (Sprint H1's slice, lepton sector)
  (b) Full SM A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ) (one generation, lepton + quark)
  (c) Simpler comparator triples (ℂ ⊕ ℂ, M_2 ⊕ M_3)

The structural claim under test: D_F in the Krajewski lineage is off-diagonal
between bimodule sectors, so σ(D_F) is symmetric about 0, so η-class
invariants (k=1) vanish identically.

Outputs go to debug/data/finite_spectral_triple_engine.json.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import sympy as sp


def cc_lepton_dirac_spectrum():
    """One-generation lepton sector of SM A_F.

    Connes-Chamseddine convention: H_F^matter = (nu_L, e_L, nu_R, e_R) ⊗ 1
    plus antimatter doubling. D_F has 4×4 matter block

        M = [[0, Y], [Y†, 0]] with Y = diag(y_nu, y_e)

    Eigenvalues of M: ±y_nu, ±y_e. After matter+antimatter doubling on H_F,
    each eigenvalue has multiplicity 2.
    """
    y_nu, y_e = sp.symbols("y_nu y_e", positive=True, real=True)
    spectrum = [
        (+y_nu, 2),
        (-y_nu, 2),
        (+y_e, 2),
        (-y_e, 2),
    ]
    return spectrum, [y_nu, y_e]


def cc_full_sm_dirac_spectrum():
    """Full SM A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ), one generation.

    Spectrum from Yukawas y_nu, y_e (lepton) + y_u, y_d (quark, color×3).
    Quark sector: Y_quark = diag(y_u, y_d), color triplets the multiplicities.
    """
    y_nu, y_e, y_u, y_d = sp.symbols("y_nu y_e y_u y_d", positive=True, real=True)
    spectrum = [
        # Lepton sector (no color)
        (+y_nu, 2), (-y_nu, 2),
        (+y_e, 2),  (-y_e, 2),
        # Quark sector (×3 color)
        (+y_u, 6), (-y_u, 6),
        (+y_d, 6), (-y_d, 6),
    ]
    return spectrum, [y_nu, y_e, y_u, y_d]


def comparator_c_plus_c_dirac_spectrum():
    """A_F = ℂ ⊕ ℂ with off-diagonal D_F = [[0, m], [m, 0]] on H_F = ℂ²."""
    m = sp.symbols("m", positive=True, real=True)
    spectrum = [(+m, 1), (-m, 1)]
    return spectrum, [m]


def comparator_m2_plus_m3_dirac_spectrum():
    """A_F = M_2(ℂ) ⊕ M_3(ℂ) with off-diagonal D_F:
    H_F = ℂ² ⊕ ℂ³ has dim 5; D_F is 5×5 with off-diagonal block of size 2×3.
    Off-diagonal block carries singular values σ_1, σ_2.
    Spectrum of D_F: ±σ_1, ±σ_2, 0 (one zero eigenvalue from rank deficit).
    """
    s1, s2 = sp.symbols("sigma_1 sigma_2", positive=True, real=True)
    spectrum = [(+s1, 1), (-s1, 1), (+s2, 1), (-s2, 1), (sp.Integer(0), 1)]
    return spectrum, [s1, s2]


def trace_Dk_heat(spectrum, k, t):
    """Tr(D^k · exp(-t D^2)) = Σ multiplicity_i · λ_i^k · exp(-t λ_i^2).

    For finite Hermitian D with eigenvalues {(λ_i, m_i)}.
    """
    return sum(m_i * (lam**k) * sp.exp(-t * lam**2) for (lam, m_i) in spectrum)


def mellin_transform_heat(trace_expr, t, s, eigenvalues):
    """𝓜[Σ c_i e^{-tλ_i²}](s) = Γ(s) · Σ c_i / λ_i^{2s} (Mellin of single
    Gaussian e^{-tλ²} is Γ(s)·λ^{-2s} for Re(s)>0, λ>0).

    For our purposes we extract the symbolic Dirichlet structure:
    write trace_expr = Σ c_i · λ_i^{k} · e^{-tλ_i²}, get Mellin =
    Γ(s) · Σ c_i · λ_i^{k-2s}.
    """
    # Symbolic by-hand: don't try sp.mellin_transform on a multi-exponential
    # sum, it tends to time out. Just structurally describe the output.
    return None  # caller inspects trace_expr directly


def engine_table_for_triple(name, spectrum, params):
    """Compute Tr(D^k · e^{-tD²}) for k ∈ {0,1,2}, simplify, classify."""
    t, s = sp.symbols("t s", positive=True)

    out = {"name": name, "params": [str(p) for p in params]}
    out["dim_H_F"] = sum(m for _, m in spectrum)
    out["spectrum"] = [(str(lam), int(m)) for lam, m in spectrum]

    for k in [0, 1, 2]:
        T = trace_Dk_heat(spectrum, k, t)
        T_simplified = sp.simplify(T)
        # Mellin output structure: Γ(s) · Σ c_i · λ_i^{k-2s}
        mellin_dirichlet = sum(
            m_i * (lam**k) * (lam ** (-2 * s)) for (lam, m_i) in spectrum if lam != 0
        )
        mellin_dirichlet = sp.simplify(mellin_dirichlet)
        out[f"k={k}"] = {
            "trace_expr": str(T_simplified),
            "is_zero": bool(T_simplified == 0),
            "mellin_dirichlet_structure": str(mellin_dirichlet),
        }
    return out


def check_eta_vanishing_general(spectrum):
    """For finite Hermitian off-diagonal D_F, σ(D_F) symmetric about 0
    ⟹ Tr(D_F · e^{-tD_F²}) ≡ 0 ⟹ η(D_F) = 0 ⟹ M3-class invariants vanish.

    Test: is the spectrum symmetric about 0?
    """
    eigvals = []
    for lam, m in spectrum:
        eigvals.extend([lam] * m)

    # Pair λ with -λ
    paired = sorted(eigvals, key=lambda x: sp.simplify(x))
    pos = [x for x in eigvals if sp.simplify(x) != 0 and not (x.is_real and x < 0) if x.is_real and x > 0]
    neg = [-x for x in eigvals if sp.simplify(x) != 0 and x.is_real and x < 0]
    # Multiset equality between pos and neg ⟹ symmetry
    return sorted(pos, key=str) == sorted(neg, key=str)


def factorized_combined_engine_check():
    """Verify that for D = D_GV ⊗ 1_F + γ_GV ⊗ D_F (CC convention),
    D² = D_GV² ⊗ 1_F + 1_GV ⊗ D_F² (cross term vanishes by γ_GV anticommuting
    with D_GV).

    For 2×2 toy example, verify symbolically.
    """
    a, b, c, d = sp.symbols("a b c d", real=True)
    # Toy D_GV = [[0, a], [a, 0]] (off-diagonal Hermitian, anticommutes with γ_GV = σ_3)
    D_GV = sp.Matrix([[0, a], [a, 0]])
    gamma_GV = sp.Matrix([[1, 0], [0, -1]])
    # Toy D_F = [[0, c], [c, 0]]
    D_F = sp.Matrix([[0, c], [c, 0]])
    I_GV = sp.eye(2)
    I_F = sp.eye(2)

    # Tensor products
    def tp(A, B):
        return sp.Matrix(sp.BlockMatrix([
            [A[0, 0] * B, A[0, 1] * B],
            [A[1, 0] * B, A[1, 1] * B],
        ]))

    D_total = tp(D_GV, I_F) + tp(gamma_GV, D_F)
    D_total_sq = D_total * D_total
    D_total_sq = sp.Matrix([[sp.simplify(D_total_sq[i, j])
                             for j in range(4)] for i in range(4)])

    expected = tp(D_GV * D_GV, I_F) + tp(I_GV, D_F * D_F)
    expected = sp.Matrix([[sp.simplify(expected[i, j])
                           for j in range(4)] for i in range(4)])

    diff = D_total_sq - expected
    diff = sp.Matrix([[sp.simplify(diff[i, j])
                       for j in range(4)] for i in range(4)])

    return {
        "D_total_sq": str(D_total_sq.tolist()),
        "expected_factorized": str(expected.tolist()),
        "difference": str(diff.tolist()),
        "factorizes": bool(diff == sp.zeros(4, 4)),
    }


def main():
    results = {}

    # Leg 1: Engine on canonical finite triples
    cases = [
        ("SM_lepton_sector", *cc_lepton_dirac_spectrum()),
        ("SM_full_one_generation", *cc_full_sm_dirac_spectrum()),
        ("comparator_C_plus_C", *comparator_c_plus_c_dirac_spectrum()),
        ("comparator_M2_plus_M3", *comparator_m2_plus_m3_dirac_spectrum()),
    ]

    triple_results = []
    for name, spectrum, params in cases:
        triple_results.append(engine_table_for_triple(name, spectrum, params))
    results["leg_1_engine_on_triples"] = triple_results

    # Leg 2: η-vanishing check on each
    eta_check = []
    for name, spectrum, params in cases:
        spectrum_symbolic = spectrum
        # Sum λ^1 · m over spectrum and check if = 0
        S1 = sum(m * lam for lam, m in spectrum)
        eta_check.append({
            "name": name,
            "sum_lambda_m": str(sp.simplify(S1)),
            "eta_zero": bool(sp.simplify(S1) == 0),
        })
    results["leg_2_eta_vanishing"] = eta_check

    # Leg 3: Combined T_GV ⊗ T_F factorization
    results["leg_3_combined_factorization"] = factorized_combined_engine_check()

    # Output
    out_path = Path("debug/data/finite_spectral_triple_engine.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Wrote {out_path}")

    # Print summary
    print("\n=== LEG 1: Engine on finite triples ===")
    for tr in triple_results:
        print(f"\n{tr['name']} (dim_H_F = {tr['dim_H_F']})")
        for k in [0, 1, 2]:
            entry = tr[f"k={k}"]
            print(f"  k={k}: zero={entry['is_zero']}")
            print(f"    Tr = {entry['trace_expr']}")
            print(f"    Mellin Dirichlet: {entry['mellin_dirichlet_structure']}")

    print("\n=== LEG 2: η-vanishing ===")
    for ec in eta_check:
        print(f"  {ec['name']}: Σλ·m = {ec['sum_lambda_m']}, η-zero = {ec['eta_zero']}")

    print("\n=== LEG 3: D² factorization ===")
    print(f"  D² = D_GV² ⊗ 1 + 1 ⊗ D_F² ?  {results['leg_3_combined_factorization']['factorizes']}")


if __name__ == "__main__":
    main()
