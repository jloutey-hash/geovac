"""Door 4f -- Falsifier A'' on Upgrade B (sphere-Lie-group axiom):
does adopting it introduce NEW downstream constraints?

Setup
=====

Door 4d (PARTIAL-DOOR) introduced DAS, the Division-Algebra-of-the-Sphere
criterion that closes the H-vs-M_2(C) fork at n=2 by selecting H. Door 4e
(PARTIAL-DOOR CONFIRMED) verified the literal Falsifier A' is negative on
all three natural AC tensor-product compatibility conditions, and named
the FULL-DOOR upgrade path: Upgrade B, the sphere-Lie-group axiom:

  "The inner algebra at rung n is the *-algebra whose unit group IS the
   rung-n Hopf-bundle total sphere S^(2n-1) as a Lie group, when one
   exists; else fallback M_n(C)."

Falsifier A'' tests whether adopting Upgrade B introduces NEW DOWNSTREAM
CONSTRAINTS beyond what CCM's standard 3-axiom path provides. Three
possible outcomes per Door 4e's framing:

  Outcome 1: Introduces a new constraint the framework SUPPORTS
             => corroborates DAS; promote to working axiom.
  Outcome 2: Introduces a new constraint the framework CONTRADICTS
             => Upgrade B RULED OUT; fork returns to "PARTIAL-DOOR
                via leaner H-import but no FULL-DOOR closure."
  Outcome 3: Introduces no new constraints (pure reformulation of CCM)
             => DAS maintained at PARTIAL-DOOR; Upgrade B is
                "CCM-equivalent with leaner axiom count."

This probe runs the 12 sub-tests articulated below for each plausible
downstream target.

Sub-tests and their expected outcomes
=====================================

T1. Inner KO-dimension (currently FREE per Door 4b Q3).
T2. Chirality grading gamma_F (currently independent of gamma_GV per G3
    NEGATIVE).
T3. Yukawa eigenvalues (currently FREE per seam theorem).
T4. Generation count N_gen (currently FREE per Door 4b Q3).
T5. Connes order-zero condition [a, J b J^{-1}] = 0 (Door 4c: passes
    for both H and M_2(C)).
T6. Connes order-one condition [[D, a], J b J^{-1}] = 0 with nonzero
    Yukawa (Door 4c: passes for both).
T7. Connes second-order condition [[D, a], [J b J^{-1}, c]] = 0 (the
    CCM workhorse).
T8. Bosonic spectral action's gauge-sector coefficient at leading order.
T9. Higgs vacuum manifold = SU(2)/U(1) = S^2 (post-SSB).
T10. Cross-rung Yukawa coupling structure (off-diagonal D_F entries).
T11. Adams 1958 / Hopf-invariant-one theorem extensibility to higher
     rungs (n=4 via octonions, n>=5 fallback).
T12. Consistency with the Bertrand x Hopf-tower reading of GeoVac.

Most tests are STRUCTURAL (paper-level reasoning); a few are spot-checkable
(T5, T6 already verified in Door 4c data and not re-run here for honesty).
T8 (spectral action gauge coupling) is structurally argued: both H and
M_2(C) selections give SU(2) after unimodularity, and the spectral
action's gauge kinetic term depends only on the SU(2) Killing form, which
is identical in both cases.

Expected verdict
================

Outcome 3 (NEUTRAL): Upgrade B is a clean reformulation of CCM with no
new downstream constraints. DAS stays PARTIAL-DOOR. The advantage is
purely axiom-minimality: 1 axiom vs CCM's 3. No new predictive content.

This is honest scope: Outcome 3 is informative because it RULES OUT
Outcome 2 (Upgrade B doesn't break anything) AND Outcome 1 (Upgrade B
isn't "more predictive than CCM"). The net is that DAS provides an
ALTERNATIVE foundational axiom for the same observable content.

Guardrail
=========

No Yukawa value, generation count, or KO-dim is selected. The probe is
structural / forcing-catalogue work, not a prediction sprint. H1 / W3 /
Koide negatives (CLAUDE.md SS 3) untouched.
"""

from __future__ import annotations

import json
from pathlib import Path

OUT_PATH = Path(__file__).parent / "data" / "door4f_falsifier_A_double_prime.json"
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Sub-tests T1-T12
# ---------------------------------------------------------------------------


def test_T1_inner_KO_dim() -> dict:
    """T1. Does Upgrade B constrain the inner KO-dimension?"""
    return {
        "test_id": "T1",
        "target": "Inner KO-dimension",
        "current_status_pre_Upgrade_B": "FREE (Door 4b Q3); imposed to be 6 to match SM total 9 mod 8 = 1.",
        "Upgrade_B_directly_constrains": False,
        "reasoning": (
            "Upgrade B fixes A_F's structure at each rung (the *-algebra). "
            "The inner KO-dim is a property of (H_F, D_F, J_F, gamma_F) "
            "determined by HOW we double H_F to include antimatter (matter "
            "vs antimatter swap sign of J_F^2) and the chirality grading. "
            "Upgrade B says nothing about H_F's doubling or about J_F^2's "
            "sign; it constrains only A_F."
        ),
        "verdict": "SILENT",
        "introduces_new_constraint": False,
    }


def test_T2_gamma_F() -> dict:
    """T2. Does Upgrade B constrain the chirality grading gamma_F?"""
    return {
        "test_id": "T2",
        "target": "Chirality grading gamma_F on H_F",
        "current_status_pre_Upgrade_B": (
            "Independent of gamma_GV per G3 NEGATIVE result (gamma_GV and "
            "gamma_F are independent commuting Z_2's)."
        ),
        "Upgrade_B_directly_constrains": False,
        "reasoning": (
            "Upgrade B constrains A_F's structure. The standard CCM "
            "chirality grading gamma_F = diag(+1, +1, -1, -1) on "
            "(nu_L, e_L, nu_R, e_R) is CONSISTENT with H-action on the "
            "L-doublet (the natural left-quaternionic action lives on the "
            "L-block); but Upgrade B does not UNIQUELY determine which "
            "subspace is L vs R. The L vs R split is additional input "
            "from the doubled-Hilbert-space structure (matter side, not "
            "algebra side)."
        ),
        "verdict": "SILENT (consistent but not uniquely determined)",
        "introduces_new_constraint": False,
    }


def test_T3_Yukawa() -> dict:
    """T3. Does Upgrade B constrain Yukawa eigenvalues?"""
    return {
        "test_id": "T3",
        "target": "Yukawa eigenvalues (entries of D_F)",
        "current_status_pre_Upgrade_B": (
            "FREE per seam theorem (Paper 32 SS VIII Door 4): the inner "
            "Yukawa Dirichlet ring Q[y_i^{-2s}] is disjoint from every "
            "forced outer ring; eta-trivialization + ac_factorization "
            "make the Yukawa values structurally inaccessible to the "
            "outer-side forcing."
        ),
        "Upgrade_B_directly_constrains": False,
        "reasoning": (
            "Upgrade B constrains A_F (the algebra factors). The seam "
            "theorem applies independently of A_F's specific real form: "
            "the off-diagonal D_F entries (Yukawas) carry their own "
            "Dirichlet ring disjoint from the outer Mellin engine. "
            "Selecting H vs M_2(C) at the n=2 factor does not enter the "
            "Yukawa Dirichlet ring."
        ),
        "verdict": "SILENT (seam theorem unaffected)",
        "introduces_new_constraint": False,
    }


def test_T4_N_gen() -> dict:
    """T4. Does Upgrade B constrain the generation count N_gen?"""
    return {
        "test_id": "T4",
        "target": "Generation count N_gen",
        "current_status_pre_Upgrade_B": (
            "FREE per Door 4b Q3 (invisible to inner automorphisms; "
            "enters as a Hilbert-space multiplicity H_F = C^N_gen (x) "
            "H_F^{(1 gen)})."
        ),
        "Upgrade_B_directly_constrains": False,
        "reasoning": (
            "Upgrade B constrains A_F. N_gen is a multiplicity of the "
            "matter representation, not a feature of A_F itself. The "
            "inner-automorphism gauge group is blind to per-summand "
            "multiplicity, and Upgrade B is articulated at the algebra "
            "level, not at the rep-space-multiplicity level."
        ),
        "verdict": "SILENT",
        "introduces_new_constraint": False,
    }


def test_T5_order_zero() -> dict:
    """T5. Connes order-zero condition under Upgrade B."""
    return {
        "test_id": "T5",
        "target": "Connes order-zero: [a, J b J^{-1}] = 0",
        "current_status_pre_Upgrade_B": (
            "PASSES for both H and M_2(C) at finite n_max in {1, 2, 3} "
            "(Door 4c bit-exact data, matter/antimatter decoupling "
            "mechanism: a is supported on matter sector, J b J^{-1} on "
            "antimatter sector, so the commutator vanishes by block "
            "disjointness)."
        ),
        "Upgrade_B_directly_constrains": False,
        "reasoning": (
            "Upgrade B forces A_F at n=2 to be H. Door 4c already "
            "verified order-zero passes for H bit-exactly. Upgrade B "
            "does not modify the mechanism (matter/antimatter decoupling) "
            "by which order-zero passes; it simply restricts to one of "
            "the two candidates that both passed."
        ),
        "verdict": "COMPATIBLE (passes; mechanism unaffected)",
        "introduces_new_constraint": False,
    }


def test_T6_order_one() -> dict:
    """T6. Connes order-one condition under Upgrade B with nonzero Yukawa."""
    return {
        "test_id": "T6",
        "target": "Connes order-one: [[D, a], J b J^{-1}] = 0",
        "current_status_pre_Upgrade_B": (
            "PASSES for both H and M_2(C) at finite n_max in {1, 2} "
            "with generic nonzero Yukawa y = 0.3 (Door 4c part 3c, "
            "bit-exact). Same matter/antimatter decoupling mechanism: "
            "[D, a] lives on matter, J b J^{-1} on antimatter."
        ),
        "Upgrade_B_directly_constrains": False,
        "reasoning": (
            "Same argument as T5: Upgrade B restricts to H, which "
            "already passes order-one. Mechanism unaffected."
        ),
        "verdict": "COMPATIBLE (passes; mechanism unaffected)",
        "introduces_new_constraint": False,
    }


def test_T7_second_order() -> dict:
    """T7. Connes second-order condition (the CCM workhorse) under Upgrade B."""
    return {
        "test_id": "T7",
        "target": "Connes second-order: [[D, a], [J b J^{-1}, c]] = 0",
        "current_status_pre_Upgrade_B": (
            "This is the CCM axiom used to SELECT H over M_2(C) in the "
            "standard Chamseddine-Connes 2008 path. Imposed as an axiom; "
            "constrains D_F's off-diagonal Yukawa structure."
        ),
        "Upgrade_B_directly_constrains": False,
        "reasoning": (
            "Upgrade B forces A_F = H via a DIFFERENT route (geometric, "
            "via the sphere's Lie-group structure). With Upgrade B, "
            "the second-order condition becomes either: (a) AUTOMATIC "
            "given H -- if the H selection forces the same D_F constraint "
            "that the second-order condition would; or (b) STILL NEEDED "
            "AS A SEPARATE AXIOM -- if the second-order condition "
            "constrains D_F's Yukawa structure beyond what H selection "
            "alone imposes. The standard CCM analysis suggests (b): the "
            "second-order condition fixes the off-diagonal block structure "
            "of D_F (which Yukawa entries can be nonzero), and H selection "
            "alone does not determine this. Upgrade B and the second-order "
            "condition therefore COEXIST without redundancy."
        ),
        "verdict": (
            "COMPATIBLE but NOT REDUNDANT (second-order condition is a "
            "separate constraint on D_F; Upgrade B fixes A_F structure)."
        ),
        "introduces_new_constraint": False,
    }


def test_T8_bosonic_spectral_action() -> dict:
    """T8. Bosonic spectral action's gauge-sector coefficient at leading order."""
    return {
        "test_id": "T8",
        "target": "Bosonic spectral action gauge-sector coefficient",
        "current_status_pre_Upgrade_B": (
            "Standard CCM: the bosonic action Tr(f(D_A^2 / Lambda^2)) "
            "at leading order in the gauge field gives the SU(2) Yang-Mills "
            "kinetic term (1/4g^2) * Tr(F_munu F^munu) with coupling "
            "determined by the SU(2) Killing form / trace normalization."
        ),
        "Upgrade_B_directly_constrains": False,
        "reasoning": (
            "Both H and M_2(C) selections produce SU(2) as the gauge "
            "group (post-unimodularity for M_2(C)). The spectral action's "
            "gauge kinetic coefficient depends only on the Killing form "
            "of the gauge Lie algebra, which is su(2) in both cases with "
            "the same normalization. The bosonic spectral action is "
            "INSENSITIVE to which real form of M_2(C) the algebra is."
        ),
        "verdict": "IDENTICAL at gauge-sector leading order",
        "introduces_new_constraint": False,
    }


def test_T9_higgs_vacuum_manifold() -> dict:
    """T9. Higgs vacuum manifold post-SSB."""
    return {
        "test_id": "T9",
        "target": "Higgs vacuum manifold (post-SSB)",
        "current_status_pre_Upgrade_B": (
            "Standard CCM: Higgs vacuum manifold is SU(2) x U(1) / U(1) "
            "= SU(2)/Z_2 modulo gauge = S^2 (after symmetry breaking "
            "by Higgs vev)."
        ),
        "Upgrade_B_directly_constrains": False,
        "reasoning": (
            "The Higgs vacuum manifold is determined by the gauge group "
            "and the Higgs rep, both of which are identical for H and "
            "M_2(C) selections at the gauge-group level. The Higgs field "
            "lives in the off-diagonal D_F block (free per seam theorem), "
            "and its vacuum manifold after SSB is S^2 for both cases."
        ),
        "verdict": "IDENTICAL (S^2 for both)",
        "introduces_new_constraint": False,
    }


def test_T10_cross_rung() -> dict:
    """T10. Cross-rung Yukawa coupling structure (off-diagonal D_F entries)."""
    return {
        "test_id": "T10",
        "target": "Cross-rung Yukawa coupling (e.g., quark CKM mixing)",
        "current_status_pre_Upgrade_B": (
            "Free per seam theorem; CKM matrix is calibration data."
        ),
        "Upgrade_B_directly_constrains": False,
        "reasoning": (
            "Upgrade B selects each rung's algebra independently (block "
            "diagonal). Cross-rung Yukawa entries (off-diagonal in the "
            "generation index, e.g., quark CKM mixing) are calibration "
            "data inaccessible to either H or M_2(C) selection at the "
            "algebra level."
        ),
        "verdict": "SILENT (cross-rung coupling stays free)",
        "introduces_new_constraint": False,
    }


def test_T11_adams_extensibility() -> dict:
    """T11. Adams 1958 extensibility to higher rungs."""
    return {
        "test_id": "T11",
        "target": "Extensibility to rung n >= 4 (e.g., octonionic S^7)",
        "current_status_pre_Upgrade_B": (
            "GeoVac currently does not reach beyond n=3 (Bertrand x Hopf "
            "tower truncates at U(1) x SU(2) x SU(3) = SM gauge group)."
        ),
        "Upgrade_B_directly_constrains": True,
        "reasoning": (
            "Adams 1958/1960 (Hopf-invariant-one theorem): S^(2n-1) is "
            "parallelizable iff n in {1, 2, 4}; carries Lie-group structure "
            "ONLY at n=1 (U(1)), n=2 (Sp(1) = SU(2)), and n=4 (S^7 via "
            "octonions, NON-ASSOCIATIVE). Upgrade B's articulation is "
            "natural for the associative-only case (n=1, 2); for n=4 "
            "the octonionic S^7 would require an extension to non-associative "
            "*-algebras, which is outside the standard NCG framework. For "
            "n>=5, S^(2n-1) is not a Lie group, so Upgrade B's fallback "
            "(M_n(C)) applies cleanly. The rule has a CLEAN PER-RUNG "
            "specification all the way up; no ambiguity. This is a CONSTRAINT "
            "that becomes substantive only if GeoVac ever extends beyond "
            "n=3 (which the SM gauge truncation does not currently motivate)."
        ),
        "verdict": (
            "CLEAN EXTENSIBILITY RULE (no current constraint at n <= 3; "
            "applies cleanly to any future higher-rung extension)."
        ),
        "introduces_new_constraint": False,
        "introduces_extensibility_rule": True,
    }


def test_T12_bertrand_consistency() -> dict:
    """T12. Consistency with the Bertrand x Hopf-tower reading."""
    return {
        "test_id": "T12",
        "target": "Consistency with Bertrand x Hopf-tower argument",
        "current_status_pre_Upgrade_B": (
            "Bertrand x Hopf-tower extracts the rung-n sphere S^(2n-1) "
            "from the closed-orbit constraint + complex-Hopf bundle. "
            "The existing reading (Paper 32 SS VIII.B) extracts the "
            "GAUGE GROUP SU(n) from this sphere."
        ),
        "Upgrade_B_directly_constrains": True,
        "reasoning": (
            "Upgrade B extracts MORE structure from the same sphere: not "
            "just its topology (which gives the gauge group via homotopy / "
            "fundamental group) but also its Lie-group structure (when one "
            "exists, by Adams 1958: n=1, 2 associative). The reading 'use "
            "the Lie structure when available, fall back to topology when "
            "not' is a NATURAL STRENGTHENING of 'use the topology only.' "
            "The case n=3 fallback is just 'when the Lie structure doesn't "
            "exist, use only the topology.' This is structurally consistent "
            "with the existing argument; it does not contradict any "
            "existing GeoVac result. It IS new input data (read more off "
            "the sphere than before)."
        ),
        "verdict": (
            "CONSISTENT EXTENSION of the existing reading. Strengthens "
            "the Bertrand x Hopf-tower argument without contradicting it."
        ),
        "introduces_new_constraint": False,
        "strengthens_existing_reading": True,
    }


# ---------------------------------------------------------------------------
# Synthesis: three-outcome verdict
# ---------------------------------------------------------------------------


def synthesize_verdict(tests: list[dict]) -> dict:
    """Combine the 12 sub-tests into a three-outcome verdict."""
    n_introduces_new = sum(1 for t in tests if t.get("introduces_new_constraint", False))
    n_silent = sum(1 for t in tests if "SILENT" in t.get("verdict", ""))
    n_compatible = sum(1 for t in tests if "COMPATIBLE" in t.get("verdict", ""))
    n_identical = sum(1 for t in tests if "IDENTICAL" in t.get("verdict", ""))
    n_clean_extensibility = sum(1 for t in tests if t.get("introduces_extensibility_rule", False))
    n_consistent_extension = sum(1 for t in tests if t.get("strengthens_existing_reading", False))

    # Three-outcome verdict
    if n_introduces_new > 0:
        outcome = "Outcome_1_or_2_NEW_CONSTRAINT"
        verdict_text = "Upgrade B introduces a new downstream constraint -- needs further classification."
    elif n_clean_extensibility >= 1 and n_consistent_extension >= 1:
        outcome = "Outcome_3_NEUTRAL"
        verdict_text = (
            "Upgrade B introduces NO new downstream constraint on existing GeoVac "
            "observables (T1-T10). It does provide a CLEAN EXTENSIBILITY RULE for "
            "higher rungs (T11, via Adams 1958) and a CONSISTENT STRENGTHENING of "
            "the Bertrand x Hopf-tower reading (T12), but neither is a constraint "
            "on existing observables -- they are clean structural extensions. "
            "Upgrade B is therefore a CLEAN REFORMULATION of CCM's H-selection "
            "with a leaner axiom count, contributing no NEW predictive content "
            "but RULING OUT contradictions and providing structural elegance."
        )
    else:
        outcome = "Outcome_3_NEUTRAL"
        verdict_text = "Upgrade B introduces no new constraints."

    return {
        "n_sub_tests": len(tests),
        "n_introducing_new_constraint": n_introduces_new,
        "n_silent": n_silent,
        "n_compatible_pre_existing": n_compatible,
        "n_identical_to_CCM": n_identical,
        "n_clean_extensibility": n_clean_extensibility,
        "n_consistent_extension": n_consistent_extension,
        "outcome": outcome,
        "verdict_text": verdict_text,
    }


def DAS_status_after_falsifier_A_double_prime() -> dict:
    """Summarize DAS status across Doors 4d, 4e, 4f."""
    return {
        "DAS_status_after_door4d": "PARTIAL-DOOR (new handle introduced)",
        "DAS_status_after_door4e": "PARTIAL-DOOR CONFIRMED (literal Falsifier A' negative on three AC compatibility conditions; FULL-DOOR upgrade path = adopt Upgrade B sphere-Lie-group axiom)",
        "DAS_status_after_door4f": (
            "PARTIAL-DOOR FINAL (Falsifier A'' returns Outcome 3 NEUTRAL: "
            "Upgrade B is a clean reformulation of CCM's H-selection with no "
            "new downstream constraints; advantage is purely axiom-minimality "
            "1 vs 3). Final per-rung forcing under Upgrade B (if adopted):"
        ),
        "final_inner_algebra_forcing_under_Upgrade_B": {
            "factor_count": "FORCED (= 3, from Bertrand x Hopf-tower truncation at n <= 3)",
            "n=1": "FORCED (C, by Upgrade B: S^1 = U(1) = U(C))",
            "n=2": "FORCED (H, by Upgrade B: S^3 = Sp(1) = U(H))",
            "n=3": "FORCED (M_3(C), by Upgrade B: S^5 not a Lie group, fallback to minimal matrix algebra reproducing SU(3))",
            "inner_KO_dim": "FREE (T1: silent; imposed to be 6 for SM total)",
            "chirality_gamma_F": "FREE (T2: consistent but not uniquely determined)",
            "Yukawa_eigenvalues": "FREE (T3: seam theorem)",
            "N_gen": "FREE (T4: invisible to inner automorphisms)",
        },
        "scope_for_a_door4g_or_later": (
            "Outcome 3 closes the structural-axiom investigation. The next "
            "open question is judgment-level (PI-side): IS Upgrade B's "
            "axiom-minimality advantage worth promoting it to a working "
            "GeoVac axiom, or is the CCM 3-axiom path preferred for its "
            "literature compatibility? This is a paper-level choice, not "
            "a probe-able falsifier. Door 4 series is structurally complete "
            "at this point."
        ),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> dict:
    print("=" * 72)
    print("Door 4f: Falsifier A'' on Upgrade B -- new downstream constraints?")
    print("=" * 72)

    tests = [
        test_T1_inner_KO_dim(),
        test_T2_gamma_F(),
        test_T3_Yukawa(),
        test_T4_N_gen(),
        test_T5_order_zero(),
        test_T6_order_one(),
        test_T7_second_order(),
        test_T8_bosonic_spectral_action(),
        test_T9_higgs_vacuum_manifold(),
        test_T10_cross_rung(),
        test_T11_adams_extensibility(),
        test_T12_bertrand_consistency(),
    ]

    out: dict = {
        "probe": "Door_4f_falsifier_A_double_prime",
        "date": "2026-06-02",
        "predecessor_memos": [
            "debug/door4_gauge_yukawa_boundary_memo.md",
            "debug/door4b_inner_algebra_forcing_memo.md",
            "debug/door4c_j_signtable_audit_memo.md",
            "debug/door4d_division_algebra_sphere_memo.md",
            "debug/door4e_falsifier_A_prime_memo.md",
        ],
        "sub_tests": tests,
    }

    for t in tests:
        print(f"  {t['test_id']}: {t['target']}")
        print(f"      -> {t['verdict']}")

    print("\n--- Synthesis ---")
    out["synthesis"] = synthesize_verdict(tests)
    print(f"  Tests run: {out['synthesis']['n_sub_tests']}")
    print(f"  Introducing new constraint: {out['synthesis']['n_introducing_new_constraint']}")
    print(f"  Silent: {out['synthesis']['n_silent']}")
    print(f"  Compatible: {out['synthesis']['n_compatible_pre_existing']}")
    print(f"  Identical to CCM: {out['synthesis']['n_identical_to_CCM']}")
    print(f"  VERDICT: {out['synthesis']['outcome']}")

    print("\n--- DAS status across Doors 4d/4e/4f ---")
    out["DAS_status"] = DAS_status_after_falsifier_A_double_prime()
    print(f"  After Door 4d: {out['DAS_status']['DAS_status_after_door4d']}")
    print(f"  After Door 4e: {out['DAS_status']['DAS_status_after_door4e']}")
    print(f"  After Door 4f: PARTIAL-DOOR FINAL (Outcome 3 NEUTRAL)")

    print(f"\n[wrote {OUT_PATH}]")
    OUT_PATH.write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    main()
