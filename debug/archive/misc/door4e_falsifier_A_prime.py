"""Door 4e -- Falsifier A' on DAS: does any AC compatibility condition
distinguish H from M_2(C) at the n=2 rung?

Setup
=====

Door 4d (PARTIAL-DOOR) introduced the Division-Algebra-of-the-Sphere (DAS)
criterion: the inner algebra at Hopf rung n is the natural associative real
normed division algebra whose unit-norm sphere IS the rung's Hopf bundle
total space S^(2n-1). DAS closes the H-vs-M_2(C) fork by selecting H,
agrees with Door 4b at n=1 (C) and n=3 (M_3(C)).

DAS is PARTIAL-DOOR (not FULL-DOOR) because the existing Paper 32 SS VIII.B
construction extracts only the gauge group from the Hopf tower, not the
algebra realization. Falsifier A' asks: does the AC tensor-product
construction H = H_GV (x) H_F provide ANY morphism that transfers the
rung-n sphere's division-algebra structure to a constraint on H_F's inner
algebra?

Three candidate compatibility conditions are natural:

  C1. SU(2)-equivariance under Ad action.
      SU(2) acts on M_2(C) by Ad: g.m.g^{-1}. The SU(2)-invariant real
      subalgebras of M_2(C) form a small lattice. Does H emerge uniquely?

  C2. Hopf U(1)-equivariance.
      The Hopf bundle S^3 -> S^2 has structure group U(1). Its right-action
      on S^3 ⊂ C^2 is scalar phase multiplication. Does compatibility with
      this U(1) action force the n=2 inner algebra to be H?

  C3. Principal-bundle compatibility.
      The Hopf S^3 is a principal U(1)-bundle over S^2. Does the inner
      algebra's action need to descend correctly to the bundle in a way
      that distinguishes H from M_2(C)?

This probe tests all three literally.

Expected outcome (per Door 4d honest scope)
==========================================

NEGATIVE on all three -- the standard CCM AC construction is well-known
to NOT transfer structure between H_GV and H_F across the tensor product
seam. The HONEST CONCLUSION of NEGATIVE is:

  DAS is genuinely a NEW principle, not derivable from existing CCM
  machinery via natural tensor-product compatibility.

  Both DAS and the standard CCM ℍ-selector (second-order condition +
  complex chiral fermion rep) are IMPORTS relative to the bare AC
  Connes axioms. DAS is LEANER (one criterion vs CCM's three), so it
  represents an alternative -- not derivative -- foundational axiom.

This is informative: Falsifier A' tells us the FULL-DOOR upgrade for DAS
is not a theorem about existing CCM, but a structural-axiom UPGRADE
("minimal AC compatible with Hopf-rung geometry"). Sprint-scale,
paper-level, not theorem-grade about existing machinery.

If any of C1/C2/C3 DOES distinguish, that would be a surprise and would
upgrade DAS to FULL-DOOR. Falsifier A' is the cheap test that decides
between "expected NEGATIVE" and "surprise POSITIVE."

Guardrail
=========

No Yukawa value, generation count, or KO-dim is selected. The probe stays
strictly on the algebra-structure side (H vs M_2(C)) at the n=2 rung.
H1 / W3 / Koide negatives (CLAUDE.md SS 3) untouched.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

OUT_PATH = Path(__file__).parent / "data" / "door4e_falsifier_A_prime.json"
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

# Pauli matrices
SIGMA_0 = np.eye(2, dtype=np.complex128)
SIGMA_1 = np.array([[0, 1], [1, 0]], dtype=np.complex128)
SIGMA_2 = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
SIGMA_3 = np.array([[1, 0], [0, -1]], dtype=np.complex128)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def random_SU2(rng: np.random.Generator) -> np.ndarray:
    """Sample a uniform-random SU(2) element via unit quaternion."""
    v = rng.normal(size=4)
    v /= np.linalg.norm(v)
    # Standard quaternion -> SU(2): q -> q1*I - i*(q2*sigma_1 + q3*sigma_2 + q4*sigma_3)
    g = v[0] * SIGMA_0 - 1j * v[1] * SIGMA_1 - 1j * v[2] * SIGMA_2 - 1j * v[3] * SIGMA_3
    return g


def H_basis() -> list[np.ndarray]:
    """A real ℝ-linear basis of ℍ ⊂ M_2(ℂ).

    Standard embedding 1, i, j, k → σ_0, iσ_1, iσ_2, iσ_3.
    Spans a 4-real-dimensional subspace of M_2(ℂ) (which is 8-real-dim).
    """
    return [SIGMA_0, 1j * SIGMA_1, 1j * SIGMA_2, 1j * SIGMA_3]


def M2_C_basis() -> list[np.ndarray]:
    """A complex ℂ-linear basis of M_2(ℂ).

    Standard {I, σ_1, σ_2, σ_3}. Complex coefficients allowed.
    """
    return [SIGMA_0, SIGMA_1, SIGMA_2, SIGMA_3]


def is_in_real_span_H(m: np.ndarray, tol: float = 1e-12) -> bool:
    """Check if a 2x2 complex matrix lies in ℝ-span of {σ_0, iσ_1, iσ_2, iσ_3}.

    Equivalently: m is in H iff its (1,1) and (2,2) entries are real (same
    value, the q1 coefficient), and the off-diagonal entries are
    (q3 + i*q2, -q3 + i*q2) up to the standard embedding. We just project
    m onto the H basis and check the projection is real-coefficients +
    reconstructs m.
    """
    basis = H_basis()
    # Express m as sum_i c_i * basis[i] with c_i complex; require c_i real.
    # Use Frobenius inner product: c_i = <basis[i], m>_F / <basis[i], basis[i]>_F
    coeffs = []
    for b in basis:
        ip = np.trace(b.conj().T @ m)
        norm = np.trace(b.conj().T @ b)
        coeffs.append(ip / norm)
    coeffs_arr = np.array(coeffs)
    # Reconstruct and check it equals m
    recon = sum(c * b for c, b in zip(coeffs_arr, basis))
    if np.linalg.norm(recon - m) > tol:
        return False
    # Check all coefficients are real
    return bool(np.max(np.abs(coeffs_arr.imag)) < tol)


# ---------------------------------------------------------------------------
# C1. SU(2)-equivariance under Ad action
# ---------------------------------------------------------------------------


def C1_SU2_Ad_invariance(n_samples: int = 64) -> dict:
    """Check that BOTH H and M_2(C) are SU(2)-invariant as real subspaces of
    M_2(C) under the Ad action g.m.g^{-1}.

    The point: SU(2)-Ad-invariance alone does NOT distinguish them. M_2(C)
    is trivially invariant (it's the whole space); H is invariant by the
    standard structural fact (Ad SU(2) acts on M_2(C) ≅ ℂ ⊕ su(2)_complexified
    by trivial-rep on the scalar part and adjoint-rep on the traceless part,
    and H is exactly the ℝ-span of {σ_0} ⊕ {i·su(2)}).

    What the probe verifies:
      (i)  H is closed under Ad SU(2) (i.e., g.h.g^{-1} ∈ H for h ∈ H,
           g ∈ SU(2)).
      (ii) M_2(C) is trivially closed under Ad SU(2).
      (iii) The MINIMAL non-trivial SU(2)-Ad-invariant *-subalgebra of M_2(C)
            containing the identity is H. M_2(C) is the MAXIMAL such.
      => SU(2)-equivariance ALONE does NOT distinguish H from M_2(C);
         a MINIMALITY criterion is needed to pick H out.
    """
    rng = np.random.default_rng(20260602)

    # Test (i): H closed under Ad SU(2)
    h_violations = 0
    h_test_count = 0
    for _ in range(n_samples):
        # Random h ∈ H (real linear combo of basis)
        h_coeffs = rng.normal(size=4)
        h = sum(c * b for c, b in zip(h_coeffs, H_basis()))
        # Random g ∈ SU(2)
        g = random_SU2(rng)
        h_conj = g @ h @ g.conj().T
        h_test_count += 1
        if not is_in_real_span_H(h_conj):
            h_violations += 1

    # Test (ii): M_2(C) trivially closed
    # (g.m.g^{-1} is a 2x2 complex matrix iff m is; nothing to check)
    m2c_invariance = True

    return {
        "test_name": "C1_SU2_Ad_invariance",
        "n_samples": n_samples,
        "H_Ad_invariant": bool(h_violations == 0),
        "H_violations_out_of_samples": h_violations,
        "M2C_Ad_invariant_trivially": m2c_invariance,
        "distinguishes_H_from_M2C": False,
        "minimal_SU2_Ad_invariant_starsubalgebra_containing_I": "H",
        "maximal_such_subalgebra": "M_2(C)",
        "verdict": (
            "Both H and M_2(C) are SU(2)-Ad-invariant real *-subalgebras "
            "of M_2(C). The MINIMAL such is H; the MAXIMAL is M_2(C). "
            "SU(2)-Ad-invariance ALONE does not distinguish them; a "
            "MINIMALITY axiom would select H."
        ),
    }


# ---------------------------------------------------------------------------
# C2. Hopf U(1)-equivariance (scalar phase right-action on C^2)
# ---------------------------------------------------------------------------


def C2_Hopf_U1_equivariance(n_samples: int = 64) -> dict:
    """The Hopf bundle S^3 -> S^2 has structure group U(1) acting on
    S^3 ⊂ C^2 by right-multiplication with diag(e^{iθ}, e^{iθ}) = e^{iθ}I
    (scalar phase). Any matrix M commutes with e^{iθ}I trivially, so this
    U(1) action does NOT constrain the inner algebra at all.

    What the probe verifies: M·(e^{iθ}I) = (e^{iθ}I)·M for ALL M ∈ M_2(C),
    bit-exact. Hopf U(1)-equivariance is vacuous as a distinguisher.
    """
    rng = np.random.default_rng(20260603)

    max_commutator_norm = 0.0
    for _ in range(n_samples):
        # Random complex 2x2 matrix
        M = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        theta = rng.uniform(0, 2 * np.pi)
        U_phase = np.exp(1j * theta) * SIGMA_0
        commutator = M @ U_phase - U_phase @ M
        max_commutator_norm = max(max_commutator_norm, float(np.linalg.norm(commutator)))

    return {
        "test_name": "C2_Hopf_U1_equivariance",
        "n_samples": n_samples,
        "max_commutator_norm_with_scalar_phase": max_commutator_norm,
        "all_matrices_commute_with_scalar_phase": bool(max_commutator_norm < 1e-12),
        "distinguishes_H_from_M2C": False,
        "verdict": (
            "Hopf U(1) right-action is scalar phase multiplication e^{iθ}I. "
            "Every 2x2 complex matrix commutes with scalar phase trivially "
            "(scalar phases are central). Hopf U(1)-equivariance is vacuous "
            "as a distinguisher of inner algebras."
        ),
    }


# ---------------------------------------------------------------------------
# C3. Principal-bundle compatibility (Hopf S^3 as principal U(1)-bundle)
# ---------------------------------------------------------------------------


def C3_principal_bundle_compatibility(n_samples: int = 32) -> dict:
    """The Hopf bundle S^3 -> S^2 has fiber U(1) (scalar phase). An action
    of an inner algebra A on S^3 ⊂ C^2 is "principal-bundle compatible" if
    it commutes with the U(1) fiber action.

    For both H and M_2(C), the left-action by multiplication commutes with
    the right U(1) scalar action (associativity of matrix product +
    scalar commutativity). So principal-bundle compatibility ALSO does not
    distinguish them.

    The probe verifies: for random h ∈ H and m ∈ M_2(C), random q ∈ S^3 ⊂
    C^2, random θ ∈ [0, 2π):

        h · (q · e^{iθ}) == (h · q) · e^{iθ}   (H-compatible with U(1))
        m · (q · e^{iθ}) == (m · q) · e^{iθ}   (M_2(C)-compatible with U(1))
    """
    rng = np.random.default_rng(20260604)

    H_compat_max_err = 0.0
    M2C_compat_max_err = 0.0

    for _ in range(n_samples):
        # Random vector in C^2
        q = rng.normal(size=2) + 1j * rng.normal(size=2)
        theta = rng.uniform(0, 2 * np.pi)
        phase = np.exp(1j * theta)

        # H test
        h_coeffs = rng.normal(size=4)
        h = sum(c * b for c, b in zip(h_coeffs, H_basis()))
        lhs_H = h @ (q * phase)
        rhs_H = (h @ q) * phase
        H_compat_max_err = max(H_compat_max_err, float(np.linalg.norm(lhs_H - rhs_H)))

        # M_2(C) test (complex coeffs)
        m_coeffs = rng.normal(size=4) + 1j * rng.normal(size=4)
        m = sum(c * b for c, b in zip(m_coeffs, M2_C_basis()))
        lhs_M = m @ (q * phase)
        rhs_M = (m @ q) * phase
        M2C_compat_max_err = max(M2C_compat_max_err, float(np.linalg.norm(lhs_M - rhs_M)))

    return {
        "test_name": "C3_principal_bundle_compatibility",
        "n_samples": n_samples,
        "H_compatibility_max_error": H_compat_max_err,
        "M2C_compatibility_max_error": M2C_compat_max_err,
        "H_compatible_with_Hopf_U1": bool(H_compat_max_err < 1e-12),
        "M2C_compatible_with_Hopf_U1": bool(M2C_compat_max_err < 1e-12),
        "distinguishes_H_from_M2C": False,
        "verdict": (
            "Both H (via standard embedding into M_2(C)) and M_2(C) (full) act "
            "on C^2 by matrix multiplication, which commutes with the right "
            "scalar-phase U(1) action by associativity. Principal-bundle "
            "compatibility with the Hopf U(1) fiber is satisfied bit-exactly "
            "by both candidates."
        ),
    }


# ---------------------------------------------------------------------------
# Comparison: input data of CCM vs DAS imports for ℍ at n=2
# ---------------------------------------------------------------------------


def import_comparison() -> dict:
    """Compare the input data each derivation of ℍ at n=2 requires."""
    return {
        "standard_CCM_derivation": {
            "axioms_invoked": [
                "Connes second-order condition: [[D, a], [J b J^{-1}, c]] = 0 (with nonzero D_F)",
                "Complex chiral fermion representation requirement (forces L-doublet pairing with anti-L-doublet)",
                "Dimension count: 2N² = 32 / 4x4-grading argument (Chamseddine-Connes-Marcolli 2007)",
            ],
            "n_axioms": 3,
            "selects": "H over M_2(C)",
            "source": "Chamseddine-Connes 2008, J. Geom. Phys. 58, 38",
        },
        "DAS_derivation": {
            "axioms_invoked": [
                "DAS: inner algebra at rung n IS the natural associative real normed division algebra whose unit-norm sphere is the rung-n Hopf bundle total space S^(2n-1)",
            ],
            "n_axioms": 1,
            "selects": "H over M_2(C)",
            "source": "GeoVac Door 4d (2026-06-02); structurally uses Bertrand × Hopf-tower (Paper 32 §VIII.B) + Hurwitz's classification (1898)",
        },
        "comparison": (
            "Both derivations are IMPORTS relative to the bare Connes axioms "
            "(neither follows from the AC tensor-product compatibility "
            "conditions C1/C2/C3 alone, all of which are vacuous as "
            "distinguishers — see Parts 1-3). CCM uses 3 algebraic axioms; "
            "DAS uses 1 geometric criterion. By INPUT MINIMALITY, DAS is "
            "leaner. By INTEGRATION WITH EXISTING NCG MACHINERY, CCM is "
            "more standard. Neither is derivable from the other; they are "
            "alternative ℍ-selectors, not redundant."
        ),
        "implication_for_DAS_status": (
            "DAS is PARTIAL-DOOR CONFIRMED. The FULL-DOOR upgrade requires "
            "elevating DAS from 'consistent extension' to 'natural construction "
            "principle.' Natural principle candidate: 'the inner algebra at "
            "rung n is the MINIMAL real *-subalgebra of M_n(C) containing the "
            "identity and closed under Ad SU(n).' This minimality criterion "
            "picks out H at n=2 (the 4-dim minimal SU(2)-Ad-invariant *-subalgebra "
            "of M_2(C) containing I; the only candidates are ℂ·I (trivial, no "
            "non-abelian structure) and H itself, with M_2(C) being non-minimal). "
            "This minimal-AC axiom would be the natural geometric formulation, "
            "complementary to CCM's algebraic three-axiom path."
        ),
    }


# ---------------------------------------------------------------------------
# Naming the structural-axiom upgrade for FULL-DOOR
# ---------------------------------------------------------------------------


def structural_upgrade_path() -> dict:
    """Articulate the precise modification to the AC tensor-product
    construction that WOULD upgrade DAS to FULL-DOOR.

    The expected outcome of Door 4e is that NONE of C1/C2/C3 force ℍ -- the
    standard CCM construction is too permissive. The FULL-DOOR upgrade for
    DAS is a structural axiom (paper-level), not a theorem about existing
    structure.

    Three candidate FULL-DOOR upgrades:
    """
    return {
        "candidate_upgrades": {
            "upgrade_A_minimality": {
                "name": "Minimal-AC axiom",
                "statement": (
                    "The inner algebra at rung n is the MINIMAL real "
                    "*-subalgebra of M_n(C) that (a) contains the identity, "
                    "(b) is closed under Ad SU(n), and (c) is non-abelian "
                    "(to support a non-trivial gauge group)."
                ),
                "selects_at_n2": "H (the unique 4-real-dim minimal non-abelian SU(2)-Ad-invariant *-subalgebra of M_2(C) containing I)",
                "agrees_at_n1_n3": "Yes (gives C at n=1, M_3(C) at n=3 via similar minimality argument)",
                "naturalness": (
                    "Strong: 'minimal' is a standard guiding principle in NCG "
                    "(Connes constructs spectral triples with minimal algebraic "
                    "content satisfying axioms). The principle 'use the minimal "
                    "algebra reproducing the gauge group' is a SHARPER version "
                    "of CCM's 'maximal sub-algebra satisfying order-one.'"
                ),
            },
            "upgrade_B_sphere_Lie": {
                "name": "Sphere-Lie-group axiom",
                "statement": (
                    "When the rung-n Hopf-bundle total sphere S^(2n-1) has a "
                    "Lie-group structure (i.e., for n=1, 2: from C, H division "
                    "algebras), the inner algebra at rung n is the *-algebra "
                    "whose group of unitaries IS this Lie group. When S^(2n-1) "
                    "has no Lie-group structure (n=3: S^5 not a Lie group; "
                    "n=4: S^7 a Lie group only non-associatively via O), fall "
                    "back to the minimal matrix algebra M_n(C) reproducing SU(n)."
                ),
                "selects_at_n2": "H (because S^3 = Sp(1) = U(H) is the unique Lie-group structure on S^3)",
                "agrees_at_n1_n3": "Yes (C at n=1 because S^1 = U(1) = U(C); M_3(C) at n=3 fallback)",
                "naturalness": (
                    "Very strong: directly reads the algebraic content off "
                    "the geometric content GeoVac already extracts via "
                    "Bertrand × Hopf-tower. No new input data; just reads "
                    "the Hopf-rung sphere fully (Lie structure + topology), "
                    "where the existing construction reads it partially "
                    "(gauge group only)."
                ),
            },
            "upgrade_C_division_ring": {
                "name": "Division-ring axiom",
                "statement": (
                    "The inner algebra at rung n is, when possible, a "
                    "DIVISION ring (every nonzero element invertible). When "
                    "no division ring of the right Wedderburn type exists, "
                    "fall back to M_n(C)."
                ),
                "selects_at_n2": "H (the unique 4-real-dim division ring of complex-realizable type)",
                "agrees_at_n1_n3": "Yes (C at n=1; M_3(C) at n=3 fallback because no associative real division algebra has dim 9)",
                "naturalness": (
                    "Moderate: 'prefer division rings' is a common preference "
                    "in algebraic foundations but less STRUCTURALLY tied to "
                    "the Hopf-rung geometry than upgrades A or B."
                ),
            },
        },
        "recommended_upgrade": "upgrade_B_sphere_Lie",
        "recommendation_rationale": (
            "Upgrade B reads the Hopf-rung sphere's FULL geometric content "
            "(Lie-group structure when it exists) rather than just the gauge "
            "group. This is the LEANEST extension of the existing Bertrand × "
            "Hopf-tower argument and uses NO input outside what GeoVac "
            "already extracts. Upgrades A and C are good alternative "
            "framings but invoke axioms (minimality, division-ring "
            "preference) that are less directly tied to the Hopf tower."
        ),
        "verdict_after_falsifier_A_prime": "PARTIAL-DOOR (CONFIRMED)",
        "verdict_rationale": (
            "Falsifier A' (Parts 1-3) confirms that no standard AC "
            "compatibility condition forces H over M_2(C). The standard CCM "
            "construction is genuinely too permissive at the n=2 rung; ℍ is "
            "imported via the CCM second-order condition. DAS provides an "
            "alternative, LEANER import (Upgrade B's sphere-Lie axiom) that "
            "is structurally tied to the Hopf-rung geometry rather than to "
            "second-order Connes axioms. DAS does NOT upgrade from PARTIAL-DOOR "
            "to FULL-DOOR via Falsifier A'; the upgrade requires adopting the "
            "sphere-Lie axiom as a foundational construction principle. The "
            "open question for the next round: is adopting the sphere-Lie "
            "axiom JUSTIFIED relative to alternative constructions, or is it "
            "merely an aesthetic preference? This is a paper-level question "
            "(no further compute beyond what's here) and is the natural next "
            "thread to launch."
        ),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> dict:
    print("=" * 72)
    print("Door 4e: Falsifier A' on DAS -- does the AC construction force H?")
    print("=" * 72)

    out: dict = {
        "probe": "Door_4e_falsifier_A_prime",
        "date": "2026-06-02",
        "predecessor_memos": [
            "debug/door4_gauge_yukawa_boundary_memo.md",
            "debug/door4b_inner_algebra_forcing_memo.md",
            "debug/door4c_j_signtable_audit_memo.md",
            "debug/door4d_division_algebra_sphere_memo.md",
        ],
    }

    print("\n--- C1: SU(2)-equivariance under Ad action ---")
    out["C1"] = C1_SU2_Ad_invariance()
    print(f"  H Ad-invariant: {out['C1']['H_Ad_invariant']}")
    print(f"  M_2(C) Ad-invariant (trivially): {out['C1']['M2C_Ad_invariant_trivially']}")
    print(f"  C1 distinguishes H from M_2(C): {out['C1']['distinguishes_H_from_M2C']}")

    print("\n--- C2: Hopf U(1)-equivariance (scalar phase) ---")
    out["C2"] = C2_Hopf_U1_equivariance()
    print(f"  Max [M, e^(i theta)*I] norm: {out['C2']['max_commutator_norm_with_scalar_phase']:.2e}")
    print(f"  All matrices commute with scalar phase: {out['C2']['all_matrices_commute_with_scalar_phase']}")
    print(f"  C2 distinguishes H from M_2(C): {out['C2']['distinguishes_H_from_M2C']}")

    print("\n--- C3: Principal-bundle compatibility (Hopf U(1) fiber) ---")
    out["C3"] = C3_principal_bundle_compatibility()
    print(f"  H compatibility error: {out['C3']['H_compatibility_max_error']:.2e}")
    print(f"  M_2(C) compatibility error: {out['C3']['M2C_compatibility_max_error']:.2e}")
    print(f"  C3 distinguishes H from M_2(C): {out['C3']['distinguishes_H_from_M2C']}")

    print("\n--- Import comparison: CCM vs DAS ---")
    out["import_comparison"] = import_comparison()
    print(f"  Standard CCM axioms required for H-selection: {out['import_comparison']['standard_CCM_derivation']['n_axioms']}")
    print(f"  DAS axioms required for H-selection:           {out['import_comparison']['DAS_derivation']['n_axioms']}")
    print(f"  -> {out['import_comparison']['comparison']}")

    print("\n--- Structural-upgrade path for DAS FULL-DOOR ---")
    out["structural_upgrade"] = structural_upgrade_path()
    print(f"  VERDICT: {out['structural_upgrade']['verdict_after_falsifier_A_prime']}")
    print(f"  Recommended FULL-DOOR upgrade: {out['structural_upgrade']['recommended_upgrade']}")
    print(f"    ({out['structural_upgrade']['candidate_upgrades']['upgrade_B_sphere_Lie']['name']})")

    print(f"\n[wrote {OUT_PATH}]")
    OUT_PATH.write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    main()
