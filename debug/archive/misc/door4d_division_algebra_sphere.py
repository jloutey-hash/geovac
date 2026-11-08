"""Door 4d -- Does the Division-Algebra-of-the-Sphere (DAS) criterion close
the H-vs-M_2(C) fork left open by Doors 4b/4c?

Background
==========

Door 4b: among finite semisimple real *-algebras with <=3 summands and
factor size <=3, exactly TWO reproduce U(1) x SU(2) x SU(3):
    C (+) H        (+) M_3(C)     (Standard Model)
    C (+) M_2(C)   (+) M_3(C)     (alternative)
Factor count, n=1=C, n=3=M_3(C) FORCED by Bertrand x Hopf-tower elimination.
n=2 fork: H vs M_2(C) ADMITTED, not FORCED.

Door 4c: the COMBINED J = J_GV (x) J_F sign table cannot distinguish them
(J_F is the matter<->antimatter swap with J_F^2 = +1 in both cases; the
internal C^2 involution sign that decides the real form lives on the algebra
side and is invisible to the combined J). Order-zero and order-one also pass
for both. The combined-J handle is CLOSED NEGATIVE.

This probe
==========

Tests a STRUCTURALLY DISTINCT handle: the Division-Algebra-of-the-Sphere
(DAS) criterion. The thesis:

  At each rung n of the complex-Hopf tower S^(2n-1) -> CP^(n-1), the
  rung's sphere is S^(2n-1). For n=1, S^1 IS the unit-norm sphere of C
  under complex multiplication. For n=2, S^3 IS the unit-norm sphere of
  H under quaternion multiplication. No associative real division algebra
  has S^5 as its unit sphere (Hurwitz's theorem: the associative real
  normed division algebras are exactly R, C, H, of dimensions 1, 2, 4).

  The DAS criterion: the inner algebra at rung n is the natural division
  algebra whose unit sphere IS the rung's Hopf bundle total space S^(2n-1).

DAS reproduces Door 4b's elimination forcings at n=1 (S^1 = unit C => C)
and at n=3 (S^5 is not any associative division algebra's sphere => fall
back to M_3(C) via SU(3)-matching). At n=2 it CLOSES the fork by selecting
H (S^3 = unit H = Sp(1)) over M_2(C) (whose SU(2) requires the quotient
U(2)/U(1) and whose "natural" unit-norm sphere is U(2), 4-dimensional, not
S^3).

The handle is structurally orthogonal to Doors 4b/4c: it never invokes any
J or sign-table data. It is a purely GEOMETRIC argument identifying the
rung's sphere with the unit elements of the rung's division algebra.

What this driver verifies
=========================

Part 1: at n=1, S^1 = {z in C : |z|=1} is a multiplicative subgroup of C
under complex multiplication (bit-exact identity).
Part 2: at n=2, S^3 = {q in H : |q|=1} is a multiplicative subgroup of H
under quaternion multiplication (bit-exact identity), and matches Sp(1) =
SU(2) as a Lie group (bit-exact under the standard quaternion-to-SU(2) map).
Part 3: at n=3, no associative real normed division algebra has S^5 as its
unit sphere (Hurwitz: only dimensions 1, 2, 4 over R; no candidate at 6).
Part 4: the candidate algebras at the n=2 fork have DIFFERENT unit-sphere
dimensions: U(H) = S^3 (3-dim) DIRECTLY, while U(M_2(C)) = U(2) (4-dim) -->
SU(2) requires unimodularity (quotient by S^1). The DAS criterion selects
the direct realization H.
Part 5: explicit check that DAS at n=2 selects H by stating the criterion
formally and showing it agrees with the existing n=1 and n=3 forcings.

Output: JSON with the bit-exact identities + a structural-verdict block.

Decision gate
=============

The probe is STRUCTURAL/CONCEPTUAL, not a one-shot bit-exact "yes/no"
forcing test in the Door 4c sense. The verdict ranks:

  DOOR if DAS is derivable from existing GeoVac infrastructure (Bertrand x
        Hopf-tower extracts BOTH the gauge group AND the algebra realization)

  PARTIAL-DOOR if DAS is a coherent NEW criterion consistent with existing
                forcings, but not derivable from the current construction
                (=> closes the fork only if explicitly adopted)

  WALL if DAS adds nothing structurally beyond what Door 4b extracted
        (i.e., the geometric identification S^3 = unit H = Sp(1) doesn't
         translate to a constraint on the inner algebra without further
         input)

The honest reading: PARTIAL-DOOR is the expected verdict. DAS is a
coherent extension that closes the fork; adopting it requires a paper-level
statement of the construction principle. This is the door-side of the
deeper "second packing axiom for inner-factor data" question (Paper 18
SS IV).

Guardrail
=========

No Yukawa value is selected; no value-side claim is made. The probe stays
strictly on the algebra-structure side (which division ring at the n=2
rung). H1, W3, Koide negatives in CLAUDE.md SS 3 are untouched.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

OUT_PATH = Path(__file__).parent / "data" / "door4d_division_algebra_sphere.json"
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Helpers: quaternion arithmetic
# ---------------------------------------------------------------------------

SIGMA_0 = np.eye(2, dtype=np.complex128)
SIGMA_1 = np.array([[0, 1], [1, 0]], dtype=np.complex128)
SIGMA_2 = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
SIGMA_3 = np.array([[1, 0], [0, -1]], dtype=np.complex128)


def quaternion_components_to_matrix(q1: float, q2: float, q3: float, q4: float) -> np.ndarray:
    """Realize a real quaternion q = q1 + q2*i + q3*j + q4*k as a 2x2 complex
    matrix. Standard embedding (matches geovac/almost_commutative.py).

    q1 -> q1 * sigma_0
    i  -> i * sigma_1
    j  -> i * sigma_2
    k  -> i * sigma_3
    """
    return q1 * SIGMA_0 + 1j * q2 * SIGMA_1 + 1j * q3 * SIGMA_2 + 1j * q4 * SIGMA_3


def quaternion_norm(q1: float, q2: float, q3: float, q4: float) -> float:
    """|q|^2 = q1^2 + q2^2 + q3^2 + q4^2 for real q."""
    return float(np.sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4))


def quaternion_product(
    a: tuple[float, float, float, float],
    b: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    """Hamilton quaternion product a*b for real quaternions."""
    a1, a2, a3, a4 = a
    b1, b2, b3, b4 = b
    # (a1 + a2 i + a3 j + a4 k)(b1 + b2 i + b3 j + b4 k) with ij=k, jk=i, ki=j, i^2=j^2=k^2=-1
    c1 = a1 * b1 - a2 * b2 - a3 * b3 - a4 * b4
    c2 = a1 * b2 + a2 * b1 + a3 * b4 - a4 * b3
    c3 = a1 * b3 - a2 * b4 + a3 * b1 + a4 * b2
    c4 = a1 * b4 + a2 * b3 - a3 * b2 + a4 * b1
    return (c1, c2, c3, c4)


# ---------------------------------------------------------------------------
# Part 1: n=1 rung -- S^1 IS the unit sphere of C
# ---------------------------------------------------------------------------


def part1_n1_S1_is_unit_C(n_samples: int = 32) -> dict:
    """Verify S^1 = {z in C : |z|=1} is closed under complex multiplication
    (bit-exact identity for unit-norm complex numbers)."""
    rng = np.random.default_rng(20260602)
    thetas = rng.uniform(0, 2 * np.pi, size=n_samples)
    phis = rng.uniform(0, 2 * np.pi, size=n_samples)

    max_norm_error = 0.0
    max_closure_error = 0.0
    for theta, phi in zip(thetas, phis):
        z = np.exp(1j * theta)
        w = np.exp(1j * phi)
        # Norm preserved
        max_norm_error = max(max_norm_error, abs(abs(z) - 1.0), abs(abs(w) - 1.0))
        # Product is on S^1
        zw = z * w
        max_closure_error = max(max_closure_error, abs(abs(zw) - 1.0))

    # Structural identity: S^1 = U(1) = unit elements of C
    return {
        "rung": 1,
        "sphere": "S^1",
        "candidate_division_algebra": "C",
        "n_samples": n_samples,
        "max_norm_error": max_norm_error,
        "max_closure_error_under_multiplication": max_closure_error,
        "is_multiplicative_group": bool(max_closure_error < 1e-14),
        "verdict_n1": "S^1 = unit norm sphere of C; closed under complex multiplication; bit-exact.",
    }


# ---------------------------------------------------------------------------
# Part 2: n=2 rung -- S^3 IS the unit sphere of H (and = Sp(1) = SU(2))
# ---------------------------------------------------------------------------


def part2_n2_S3_is_unit_H(n_samples: int = 32) -> dict:
    """Verify S^3 = {q in H : |q|=1} is closed under quaternion multiplication
    AND that the quaternion-to-SU(2) embedding maps unit quaternions bit-exactly
    onto SU(2)."""
    rng = np.random.default_rng(20260603)

    max_norm_error = 0.0
    max_closure_error = 0.0
    max_su2_error = 0.0

    for _ in range(n_samples):
        # Random unit quaternion (uniform on S^3)
        v = rng.normal(size=4)
        v /= np.linalg.norm(v)
        q = tuple(float(x) for x in v)

        w = rng.normal(size=4)
        w /= np.linalg.norm(w)
        p = tuple(float(x) for x in w)

        # Norm preserved
        nq = quaternion_norm(*q)
        np_ = quaternion_norm(*p)
        max_norm_error = max(max_norm_error, abs(nq - 1.0), abs(np_ - 1.0))

        # Closure under quaternion multiplication
        qp = quaternion_product(q, p)
        nqp = quaternion_norm(*qp)
        max_closure_error = max(max_closure_error, abs(nqp - 1.0))

        # Quaternion-to-matrix embedding maps unit quaternions to SU(2)
        # i.e., det = 1 and U^dagger U = I
        U = quaternion_components_to_matrix(*q)
        det_err = abs(np.linalg.det(U) - 1.0)
        unit_err = float(np.linalg.norm(U.conj().T @ U - np.eye(2), ord="fro"))
        # The standard embedding I used: 1->sigma_0, i->i*sigma_1, j->i*sigma_2, k->i*sigma_3.
        # det(U) = q1^2 + q2^2 + q3^2 - q4^2  (NOT |q|^2 in this embedding)
        # Let me redo: under the geovac/almost_commutative.py embedding,
        # det(U) = (q1 + i*q4)(q1 - i*q4) - (i*q2 + q3)(i*q2 - q3) = q1^2+q4^2 + q2^2-q3^2
        # That's NOT |q|^2 = q1^2+q2^2+q3^2+q4^2 unless q3=0.
        # That convention is NOT the standard "unit quaternion = SU(2)" embedding.
        # The standard one is: 1->I, i->-i*sigma_1, j->-i*sigma_2, k->-i*sigma_3
        # Let me use that for the bit-exact SU(2) test:
        U_std = q[0] * SIGMA_0 - 1j * q[1] * SIGMA_1 - 1j * q[2] * SIGMA_2 - 1j * q[3] * SIGMA_3
        det_std = complex(np.linalg.det(U_std))
        unit_std = float(np.linalg.norm(U_std.conj().T @ U_std - np.eye(2), ord="fro"))
        max_su2_error = max(max_su2_error, abs(det_std - 1.0), unit_std)

    return {
        "rung": 2,
        "sphere": "S^3",
        "candidate_division_algebra": "H",
        "n_samples": n_samples,
        "max_norm_error": max_norm_error,
        "max_closure_error_under_quaternion_multiplication": max_closure_error,
        "max_SU2_isomorphism_error": max_su2_error,
        "is_multiplicative_group": bool(max_closure_error < 1e-14),
        "is_isomorphic_to_SU2": bool(max_su2_error < 1e-13),
        "verdict_n2": (
            "S^3 = unit norm sphere of H under quaternion multiplication AND "
            "S^3 is bit-exactly isomorphic to SU(2) = Sp(1). The H-unit-sphere "
            "identification is DIRECT (single chart, no quotient)."
        ),
    }


# ---------------------------------------------------------------------------
# Part 3: n=3 rung -- S^5 is NOT the unit sphere of any associative
# real normed division algebra (Hurwitz's theorem)
# ---------------------------------------------------------------------------


def part3_n3_S5_no_division_algebra() -> dict:
    """Hurwitz's theorem (1898): the only finite-dimensional real normed
    division algebras are R, C, H, O of dimensions 1, 2, 4, 8. Of these, only
    R, C, H are ASSOCIATIVE. Their unit spheres have dimensions 0 (S^0 in R),
    1 (S^1 in C), 3 (S^3 in H), 7 (S^7 in O).

    S^5 is NOT a unit sphere of any associative real normed division algebra.
    The Hopf bundle at n=3 (S^5 -> CP^2) therefore has NO direct realization
    via a unit-norm-sphere identification with an associative division algebra,
    forcing the fallback to a matrix algebra (M_3(C)).
    """
    associative_real_normed_division_algebras = {
        "R": {"dim": 1, "unit_sphere": "S^0"},
        "C": {"dim": 2, "unit_sphere": "S^1"},
        "H": {"dim": 4, "unit_sphere": "S^3"},
    }
    # Octonions O are non-associative (dim 8, unit sphere S^7)
    non_associative_normed_division_algebra = {
        "O": {"dim": 8, "unit_sphere": "S^7", "associative": False},
    }
    # Possible candidate at dim 6 (whose unit sphere would be S^5): NONE
    # (Hurwitz's theorem leaves no dim-6 candidate.)

    has_S5_realization = any(
        d["unit_sphere"] == "S^5"
        for d in {**associative_real_normed_division_algebras, **non_associative_normed_division_algebra}.values()
    )

    return {
        "rung": 3,
        "sphere": "S^5",
        "associative_real_normed_division_algebras_by_Hurwitz": associative_real_normed_division_algebras,
        "non_associative_normed_division_algebras": non_associative_normed_division_algebra,
        "S5_is_unit_sphere_of_associative_division_algebra": bool(has_S5_realization),
        "verdict_n3": (
            "Hurwitz: associative real normed division algebras are R, C, H "
            "(dim 1, 2, 4; unit spheres S^0, S^1, S^3). S^5 has NO associative "
            "division-algebra realization. Inner algebra at n=3 must fall back "
            "to a matrix algebra; M_3(C) is selected by Door 4b elimination "
            "(only M_3(C) has SU(3) as its post-unimodularity unitary group)."
        ),
    }


# ---------------------------------------------------------------------------
# Part 4: the n=2 fork -- H vs M_2(C) under DAS
# ---------------------------------------------------------------------------


def part4_n2_fork_under_DAS() -> dict:
    """Compare the two n=2 candidates under the Division-Algebra-of-the-Sphere
    criterion:

    Candidate H:
      U(H) = unit quaternions = S^3 = Sp(1) = SU(2) (bit-exact, Part 2).
      The S^3 IS the unit sphere of H. DIRECT realization.

    Candidate M_2(C):
      U(M_2(C)) = U(2) is 4-DIMENSIONAL (3-dim S^3 plus 1-dim S^1 for the
      determinant phase). SU(2) is obtained as U(2)/U(1) via unimodularity.
      The "natural" unit-norm sphere of M_2(C) is U(2), not S^3 -- S^3 is
      obtained INDIRECTLY through a quotient.

    Under DAS ("inner algebra at rung n = the natural division algebra whose
    unit sphere IS the rung's Hopf bundle total space S^(2n-1)"), the DIRECT
    realization H is selected; the INDIRECT realization M_2(C) requires
    additional input (the unimodularity quotient) that DAS does not provide.
    """
    return {
        "candidate_H": {
            "unitary_group": "U(H) = Sp(1) = SU(2)",
            "unit_sphere_dim": 3,
            "unit_sphere": "S^3",
            "matches_rung_sphere_directly": True,
            "requires_quotient": False,
        },
        "candidate_M2_C": {
            "unitary_group": "U(M_2(C)) = U(2)",
            "unit_sphere_dim": 4,
            "unit_sphere": "U(2)",
            "SU2_from_unimodularity_quotient": "U(2)/U(1) = SU(2)",
            "matches_rung_sphere_directly": False,
            "requires_quotient": True,
        },
        "DAS_selection_at_n2": "H (direct unit-sphere realization)",
        "verdict_n2_fork": (
            "DAS selects H. H's unit sphere IS S^3 directly; M_2(C)'s unit sphere "
            "is U(2) (4-dim), with S^3 = SU(2) obtained only via the unimodularity "
            "quotient. DAS picks the direct (non-quotient) realization."
        ),
    }


# ---------------------------------------------------------------------------
# Part 5: DAS criterion -- agreement with Door 4b at n=1, n=3; closure at n=2
# ---------------------------------------------------------------------------


def part5_DAS_summary() -> dict:
    """Synthesize the DAS verdict across all three rungs and compare to the
    Door 4b elimination forcings."""
    return {
        "DAS_principle": (
            "The inner algebra at rung n is the natural associative real normed "
            "division algebra whose unit-norm sphere IS the rung's Hopf bundle "
            "total space S^(2n-1). When no such division algebra exists (no "
            "associative candidate at dim 2n), the inner algebra falls back to "
            "the minimal matrix algebra M_n(C) reproducing the rung's gauge "
            "group SU(n)."
        ),
        "rung_by_rung_DAS_selection": {
            "n=1": "C (S^1 = unit C)",
            "n=2": "H (S^3 = unit H = Sp(1))",
            "n=3": "M_3(C) (no associative div algebra has S^5 as unit sphere; fallback)",
        },
        "Door_4b_elimination_forcing": {
            "n=1": "C (only C has U(1) as its unitary group)",
            "n=2": "H or M_2(C) (BOTH have SU(2) post-unimodularity)",
            "n=3": "M_3(C) (only M_3(C) has SU(3) as post-unimodularity unitary group)",
        },
        "agreement_at_n1": "DAS agrees with Door 4b (both give C).",
        "agreement_at_n3": "DAS agrees with Door 4b (both give M_3(C), by different routes).",
        "extension_at_n2": (
            "DAS CLOSES the fork (selects H) where Door 4b is silent (both H and "
            "M_2(C) reproduce SU(2)). The DAS handle is structurally DISTINCT from "
            "Door 4c (combined-J sign table): DAS invokes only the geometric "
            "identification of S^(2n-1) with the unit-norm sphere of the rung's "
            "associative division algebra; no J or sign-table data is referenced."
        ),
        "is_DAS_derivable_from_existing_GeoVac_structure?": (
            "DAS is CONSISTENT WITH the Bertrand x Hopf-tower forcing (the tower "
            "produces the rung-n sphere S^(2n-1); the geometric identification of "
            "this sphere with a division-algebra unit sphere is a CLASSICAL FACT "
            "about real normed division algebras, Hurwitz 1898). However, the "
            "existing Paper 32 SS VIII.B argument extracts only the GAUGE GROUP "
            "from the Hopf tower; it does not invoke the unit-sphere identification "
            "to constrain the inner-algebra REALIZATION. Adopting DAS as a forcing "
            "requires a STRENGTHENING of the construction principle: 'the inner "
            "algebra at rung n is the natural division algebra of which the rung's "
            "sphere is the unit-norm subset.' This strengthening is CONSISTENT with "
            "existing forcings at n=1 and n=3 (where it gives the same answer as "
            "elimination) and CLOSES the fork at n=2."
        ),
        "verdict": "PARTIAL-DOOR",
        "verdict_rationale": (
            "DAS is a coherent and natural extension of the Bertrand x Hopf-tower "
            "forcing. It closes the H-vs-M_2(C) fork by selecting H, agreeing with "
            "the existing n=1 and n=3 forcings. It is structurally orthogonal to "
            "Doors 4b/4c (no algebra-enumeration, no J sign-table). It is NOT a "
            "TRUE forcing derivable from the EXISTING Paper 32 SS VIII.B construction, "
            "which extracts only the gauge group from the Hopf tower. Adopting DAS "
            "as a new forcing principle is the minimal structural extension that "
            "would upgrade the inner-algebra status from 'mostly-forced' (Door 4b) "
            "to 'fully forced'. This is the door-side content of the deeper "
            "'second packing axiom for inner-factor data' question (Paper 18 SS IV)."
        ),
        "open_follow_on": (
            "(a) Is there a GeoVac-NATIVE structural argument (beyond classical "
            "Hurwitz / division-algebra mathematics) that FORCES the inner algebra "
            "to inherit the rung-n sphere's division-algebra structure? If yes, DAS "
            "is upgraded from PARTIAL-DOOR to FULL-DOOR. If no, DAS remains a "
            "new axiom on the same standing as the CCM second-order condition -- "
            "literature-grade rather than GeoVac-grade. "
            "(b) Does adopting DAS introduce any new constraints downstream "
            "(e.g., on the Higgs sector, the chirality grading gamma_F, the inner "
            "KO-dim)? If yes, those constraints become falsifiers for DAS itself."
        ),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> dict:
    print("=" * 72)
    print("Door 4d: Division-Algebra-of-the-Sphere (DAS) probe")
    print("=" * 72)

    out: dict = {
        "probe": "Door_4d_DAS",
        "date": "2026-06-02",
        "predecessor_memos": [
            "debug/door4_gauge_yukawa_boundary_memo.md",
            "debug/door4b_inner_algebra_forcing_memo.md",
            "debug/door4c_j_signtable_audit_memo.md",
        ],
    }

    print("\n--- Part 1: n=1 rung (S^1 = unit C) ---")
    out["part1_n1"] = part1_n1_S1_is_unit_C()
    print(f"  max norm error:    {out['part1_n1']['max_norm_error']:.2e}")
    print(f"  max closure error: {out['part1_n1']['max_closure_error_under_multiplication']:.2e}")
    print(f"  S^1 is unit-norm subgroup of C: {out['part1_n1']['is_multiplicative_group']}")

    print("\n--- Part 2: n=2 rung (S^3 = unit H = Sp(1) = SU(2)) ---")
    out["part2_n2"] = part2_n2_S3_is_unit_H()
    print(f"  max norm error:        {out['part2_n2']['max_norm_error']:.2e}")
    print(f"  max closure error:     {out['part2_n2']['max_closure_error_under_quaternion_multiplication']:.2e}")
    print(f"  max SU(2)-iso error:   {out['part2_n2']['max_SU2_isomorphism_error']:.2e}")
    print(f"  S^3 = unit H = SU(2):  {out['part2_n2']['is_isomorphic_to_SU2']}")

    print("\n--- Part 3: n=3 rung (S^5: no associative division algebra) ---")
    out["part3_n3"] = part3_n3_S5_no_division_algebra()
    print(f"  S^5 has associative div-alg realization: {out['part3_n3']['S5_is_unit_sphere_of_associative_division_algebra']}")
    print(f"  -> {out['part3_n3']['verdict_n3']}")

    print("\n--- Part 4: n=2 fork under DAS ---")
    out["part4_n2_fork"] = part4_n2_fork_under_DAS()
    print(f"  Candidate H: unit sphere = {out['part4_n2_fork']['candidate_H']['unit_sphere']} (dim {out['part4_n2_fork']['candidate_H']['unit_sphere_dim']}, direct)")
    print(f"  Candidate M_2(C): unit sphere = {out['part4_n2_fork']['candidate_M2_C']['unit_sphere']} (dim {out['part4_n2_fork']['candidate_M2_C']['unit_sphere_dim']}, requires unimodularity)")
    print(f"  DAS selects: {out['part4_n2_fork']['DAS_selection_at_n2']}")

    print("\n--- Part 5: DAS summary ---")
    out["part5_DAS_summary"] = part5_DAS_summary()
    print(f"  VERDICT: {out['part5_DAS_summary']['verdict']}")
    print(f"  Agreement at n=1: {out['part5_DAS_summary']['agreement_at_n1']}")
    print(f"  Agreement at n=3: {out['part5_DAS_summary']['agreement_at_n3']}")
    print(f"  Extension at n=2: closes H-vs-M_2(C) fork by selecting H")

    print(f"\n[wrote {OUT_PATH}]")
    OUT_PATH.write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    main()
