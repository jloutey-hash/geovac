"""The GeoVac numeric registry — the corpus's declared quantity graph.

WHY THIS EXISTS
---------------
C16 and C17 are *blocklists*: they name values that are known-wrong and
refuse to let them re-surface.  They work, and they have caught real defects
three separate times.  But they can only guard a class after that class has
already produced a defect, and they say nothing about the relationships
*between* numbers.

Deltas 5–7 showed that the surviving defect classes are relational, not
value-level.  A census of the group4 + group6 numeric surface found:

    2,815 numeral occurrences   825 distinct
      352 appear at more than one LOCUS      (twins)
      201 appear in more than one DOCUMENT   (cross-paper coupling)

Each delta pass found 15–20 genuine defects — about 5% of that coupled
surface.  Manual review *samples* the graph; it does not traverse it.  The
three edge types it kept missing:

  * **twins** — the same quantity tabulated in two papers.  A reviewer
    scoped to one cannot see the other.  (P20's relativistic table sat
    entirely pre-correction while P14's twin was re-measured.)
  * **derivations** — λ/Q, N², a ratio against a fixed baseline, an exponent
    fitted from a printed column.  When a base value moves its derivatives
    silently do not.  Catching these needs recomputation, not reading, which
    is why they survived seven passes.
  * **conventions** — identity-included vs non-identity counts; a fit's
    point range.  Both numbers are correct; only the pairing is wrong.
    Value-checking is blind to this by construction.

So this file declares the graph.  It is an *allowlist* with structure: each
quantity carries its canonical value, the convention it is stated in, where
that value came from, and — if it is derived — the expression that produces
it.  ``check_numeric_consistency.py`` (C21) verifies the corpus against it.

MAINTENANCE RULE (the point of the artifact)
--------------------------------------------
This is a living index, maintained like ``docs/claim_test_matrix.md``.

  1. **When a measured value moves, edit it HERE FIRST**, then run C21.  The
     gate reports every locus that still carries the old value, including
     ones in other papers.  Do not hand-sweep for loci; that is the process
     this artifact replaces.
  2. **When a new load-bearing number enters the corpus, register it.**  C20
     reports high-salience unregistered numerals (multi-locus or
     multi-document) precisely so the registry cannot quietly fall behind
     the corpus.  An unregistered multi-document numeral is a maintenance
     debt, not a defect — but it is *visible*.
  3. **Never register a value you have not either measured or cited.**  The
     ``provenance`` field is not decoration: a registry of guesses would
     launder them into apparent authority.  ``MEASURED`` entries name the
     route; ``CITED`` entries name the source; ``DERIVED`` entries name the
     inputs and are recomputed on every run rather than stored.
  4. **Conventions are declared, not inferred.**  Where a quantity is quoted
     in two conventions (identity in/out), both appear, tied by the
     ``convention`` field, so a checker can tell "different quantity" from
     "wrong value".

Companion gate: ``debug/qa/check_numeric_consistency.py``.
Self-test:      ``tests/test_numeric_registry.py``.
"""
from __future__ import annotations

# ---------------------------------------------------------------------------
# Measured quantities.  value / convention / provenance / (optional) aliases.
#
# `aliases` holds OTHER legitimate renderings of the same physical quantity
# in a different convention -- so the gate can distinguish "stated in the
# other convention" from "wrong".
# ---------------------------------------------------------------------------

MEASURED = {
    # ---- two-center decompactification front (Papers 58 / 60, 2026-09-05/06) ----
    # The front R*(n) is the distance at which <ns_A|ns_B> = 1/sqrt2 (projector
    # principal angle 45 deg, ||[P_A,P_B]|| maximal).  It tracks the decay
    # length n/Z, not the mean radius n^2/Z.  Independent quadrature route:
    # tests/test_paper58_decompactification_front.py.
    "decomp_front_ratio": dict(
        value=2.14, convention="front ratio R*(n)/(n/Z): median over n = 2..8, Z = 1, "
                               "hydrogenic ns-ns overlap = 1/sqrt2 (H2+)",
        q=None,
        provenance="MEASURED 2026-09-05, debug/decompactification_R_sweep.py; "
                   "range 2.12-2.19 (n=2..8), n=1 gives 1.565",
        aliases={2.15: "mean over n = 2..8"}),
    "decomp_front_exponent": dict(
        value=0.98, convention="exponent: log-log slope of R*(n) vs n, n = 2..8",
        q=None,
        provenance="MEASURED 2026-09-05, same driver; the rejected sqrt(ZR) "
                   "window would give 2.0 here",
        aliases={1.02: "inverse fit n*(R) vs R", 1.0: "rounded"}),
    "tail_geometric_const": dict(
        value=1.58, convention="R_rel / sqrt(l_A l_B): 1s-1s tail-reach front "
                               "(|S| -> S(0)/sqrt2), charge pairs (1,1),(2,1),(3,1),(2,2)",
        q=None,
        provenance="MEASURED 2026-09-06, debug/decompactification_correlation_ladder.py "
                   "rung0a; max dev 2.0 % (four pairs), 9.1 % (dense ratio scan t in [1,8])",
        aliases={0.79: "as c in R_rel = c * 2 sqrt(l_A l_B)"}),
    "tail_t_c": dict(
        value=2.664, convention="1s exponent ratio above which no absolute |S| = 1/sqrt2 "
                                "front exists; root of (2 sqrt t/(1+t))^3 = 1/sqrt2 "
                                "(S(R) is monotone decreasing for every t, so the "
                                "united-atom value is the maximum)",
        q=None,
        provenance="closed form, 2026-09-06; pinned in "
                   "tests/test_paper58_decompactification_front.py. The ladder "
                   "driver's 2.7456 was the first point of its discrete t-scan past "
                   "the root, not the root -- caught by the independent test route.",
        aliases={2.66: "rounded"}),
    # Rungs 1-2 of the correlation ladder: percent residual of the occupation-
    # weighted principal-angle front vs the 1s-1s law at the EMPIRICAL zeta_eff.
    "front_resid_h2_hf": dict(
        value=-0.7, convention="percent residual of the front vs the 1s-1s law at "
                               "empirical zeta_eff, H2 Hartree-Fock, base 5-zeta basis", q=None,
        provenance="MEASURED 2026-09-06, decompactification_correlation_ladder.py "
                   "(R* 1.283 vs pred 1.292)", aliases={}),
    "front_resid_h2_fci": dict(
        value=-1.2, convention="percent residual of the front, H2 full CI, base basis", q=None,
        provenance="MEASURED 2026-09-06, same driver (R* 1.270 vs pred 1.285)",
        aliases={}),
    "front_resid_heh_hf": dict(
        value=0.0, convention="percent residual of the front, HeH+ Hartree-Fock, base basis", q=None,
        provenance="MEASURED 2026-09-06, same driver (R* 0.757 vs pred 0.757)",
        aliases={}),
    "front_resid_heh_fci": dict(
        value=-0.2, convention="percent residual of the front, HeH+ full CI, base basis", q=None,
        provenance="MEASURED 2026-09-06, same driver (R* 0.757 vs pred 0.758)",
        aliases={}),
    # Signed cross-center coherence (M2) front: inward shift at FCI vs HF.
    "coherence_shift_h2": dict(
        value=23.1, convention="percent inward shift of the M2 front, H2 FCI vs HF "
                               "(0.987 vs 1.285 bohr)", q=None,
        provenance="MEASURED 2026-09-06, same driver", aliases={23: "rounded"}),
    "coherence_shift_heh": dict(
        value=9.0, convention="percent inward shift of the M2 front, HeH+ FCI vs HF "
                              "(0.690 vs 0.758 bohr)", q=None,
        provenance="MEASURED 2026-09-06, same driver", aliases={}),

    # ---- trunk multi-document constants (registered FULL #7, 2026-09-06) ----
    # C21 had ZERO \gvq surface on all six trunk docs (the completeness-critic's
    # highest-value gap; same class as FULL #4).  saturation_c is the exemplar:
    # the lambda_max saturation-rate constant appears as the identical literal
    # "42.7397" in Papers 0, 7 and the group3 synthesis, so its cross-document
    # consistency now has a gate.  Closed form (Paper 0 SecVI); independently
    # re-derived by the FULL #7 code panel via two disjoint routes.
    "saturation_c": dict(
        value=42.7397,
        convention="lambda_max saturation-rate constant C: deficit = (C+o(1))/n_max^2; "
                   "closed form (pi^2/4)(2+2^{1/3})^2(1+2^{-2/3})",
        q=None,
        provenance="DERIVED closed form (Paper 0 SecVI), numerics-pinned over 591 "
                   "cutoffs; loci P0 SecVI / P7 New-Contributions+caveats / group3 "
                   "synthesis. FULL #7 code panel re-derived it two ways.",
        aliases={42.74: "2 dp", 42.739654: "6 dp", 42.6: "n_max=320 finite sample (understates)"}),
    "l2_rate_4_over_pi": dict(
        value=1.2732395447351628,  # 4/pi at full precision (tight symbolic match)
        convention="constant: asymptotic GH-convergence rate 4/pi (n_max*gamma_n/log n_max) "
                   "in the dual-Coxeter (rotation-angle) metric; 2/pi on the unit S^3; "
                   "2*sqrt2/pi under Kac's basic form -- CONVENTION-DEPENDENT, not canonical",
        q=None,
        provenance="DERIVED (Paper 38/40), numerics-pinned (doubling estimator from above); "
                   "loci P38, P40, group3 synthesis. Annotated FULL #8 follow-up 2026-09-06.",
        aliases={1.2732: "4 dp", 1.27324: "5 dp", 0.63662: "2/pi unit-S^3 half"}),
    # Slater F^0(1s,1s) coefficient and the three alpha-decomposition
    # ingredients: multi-document symbolic-fraction trunk constants that C21
    # could not verify until the symbolic-literal parser (2026-09-06).
    "slater_f0_1s": dict(
        value=0.625,
        convention="constant: Slater F^0(1s,1s) coefficient = 5/8 (F^0 = 5Z/8 on S^3); "
                   "loci Paper 7 SecV, group3 synthesis",
        q=None,
        provenance="SYMBOLIC (Paper 7 eq:f0_s3), test_paper7_vee_s3.py; = 5/8 = 0.625",
        aliases={}),
    "delta_dirac": dict(
        value=0.025,
        convention="constant: Dirac boundary-degeneracy Delta = 1/40 = (g_3^Dirac)^-1 "
                   "(alpha-decomposition ingredient, Paper 2); loci Paper 7, Paper 32",
        q=None,
        provenance="OBSERVATION-tier ingredient (Paper 2/32); = 1/40 = 0.025. "
                   "The COMBINATION K=pi(B+F-Delta) stays an Observation (hard rule).",
        aliases={}),
    "b_casimir": dict(
        value=42.0,
        convention="constant: Casimir count B = 42 (alpha-decomposition ingredient, Paper 2); "
                   "loci Paper 2, Paper 32",
        q=None,
        provenance="derived ingredient (Paper 2); integer 42.",
        aliases={}),
    "f_fock_dirichlet": dict(
        value=1.6449340668482264,  # pi^2/6 at full precision
        convention="constant: Fock Dirichlet F = pi^2/6 = zeta(2) at the packing exponent "
                   "(alpha-decomposition ingredient, Paper 2); loci Paper 2, Paper 32",
        q=None,
        provenance="derived ingredient (Paper 2); = pi^2/6 = 1.6449340668.",
        aliases={1.6449: "4 dp"}),

    # ---- TC composed (electronic-only, PK classically partitioned) ------
    # The TC/standard Pauli ratio is 1.616 to three decimals for all three
    # molecules -- the cross-molecule uniformity survived the exact-rule
    # correction; only the constant moved (retired 1.68).
    "lih_tc_pauli": dict(
        value=1353, convention="non-identity Pauli terms", q=30,
        provenance="MEASURED 2026-08-30, build_tc_composed_hamiltonian("
                   "lih_spec(), pk_in_hamiltonian=False)",
        aliases={1354: "with PK in the Hamiltonian"}),
    "beh2_tc_pauli": dict(
        value=2255, convention="non-identity Pauli terms", q=50,
        provenance="MEASURED 2026-08-30, same route", aliases={}),
    "h2o_tc_pauli": dict(
        value=3157, convention="non-identity Pauli terms", q=70,
        provenance="MEASURED 2026-08-30, same route", aliases={}),
    "tc_pauli_ratio": dict(
        value=1.616, convention="TC/standard Pauli ratio (dimensionless)",
        q=None,
        provenance="MEASURED 2026-08-30; identical to 3 dp across "
                   "LiH/BeH2/H2O at n_max=2. Fixed-basis result: it does NOT "
                   "establish basis-axis behaviour.",
        aliases={}),

    # ---- second/third-row balanced library ------------------------------
    # Pauli support IS isostructurally invariant (NaH = KH, HCl = HBr, ...);
    # QWC group counts are NOT, because the greedy grouper sorts by
    # descending |coefficient| before insertion, so identical support with
    # different coefficients gives different group counts.  Do not "correct"
    # one partner's QWC to match the other's.
    "nah_balanced_lambda": dict(
        value=19.582, convention="non-identity 1-norm (Ha), at NaH's own R = 3.566 bohr",
        q=20,
        provenance="MEASURED 2026-08-31, build_balanced_hamiltonian(nah_spec()) "
                   "after the R-default fix. The prior 20.607 was this molecule "
                   "evaluated at R = 3.015 (LiH's bond length), which the builder "
                   "substituted silently; see debug/qa/balanced_lambda_geometry_finding.md.",
        aliases={}),
    "hcl_balanced_lambda": dict(
        value=866.133, convention="non-identity 1-norm (Ha), at HCl's own R = 2.409 bohr",
        q=50,
        provenance="MEASURED 2026-08-31, same route. The prior 869.72 was HCl "
                   "evaluated at R = 3.015 (LiH's bond length).",
        aliases={}),

    # Balanced lambda_ni, the remaining ten rows of the Papers 14/20
    # column.  Registered 2026-08-31: they had no keys, which is why a
    # systematic wrong-geometry error in twelve published cells was
    # invisible to this gate.  Each at its OWN equilibrium bond length.
    "mgh2_balanced_lambda": dict(
        value=111.800, convention="non-identity 1-norm (Ha), at R = 3.261 bohr",
        q=40,
        provenance="MEASURED 2026-08-31, debug/qa/remeasure_balanced_lambda.py",
        aliases={}),
    "h2s_balanced_lambda": dict(
        value=873.632, convention="non-identity 1-norm (Ha), at R = 2.534 bohr",
        q=60,
        provenance="MEASURED 2026-08-31, debug/qa/remeasure_balanced_lambda.py",
        aliases={}),
    "ph3_balanced_lambda": dict(
        value=889.349, convention="non-identity 1-norm (Ha), at R = 2.680 bohr",
        q=70,
        provenance="MEASURED 2026-08-31, debug/qa/remeasure_balanced_lambda.py",
        aliases={}),
    "sih4_balanced_lambda": dict(
        value=909.154, convention="non-identity 1-norm (Ha), at R = 2.800 bohr",
        q=80,
        provenance="MEASURED 2026-08-31, debug/qa/remeasure_balanced_lambda.py",
        aliases={}),
    "kh_balanced_lambda": dict(
        value=28.124, convention="non-identity 1-norm (Ha), at R = 4.243 bohr",
        q=20,
        provenance="MEASURED 2026-08-31, debug/qa/remeasure_balanced_lambda.py",
        aliases={}),
    "cah2_balanced_lambda": dict(
        value=123.773, convention="non-identity 1-norm (Ha), at R = 3.807 bohr",
        q=40,
        provenance="MEASURED 2026-08-31, debug/qa/remeasure_balanced_lambda.py",
        aliases={}),
    "hbr_balanced_lambda": dict(
        value=875.313, convention="non-identity 1-norm (Ha), at R = 2.670 bohr",
        q=50,
        provenance="MEASURED 2026-08-31, debug/qa/remeasure_balanced_lambda.py",
        aliases={}),
    "h2se_balanced_lambda": dict(
        value=883.443, convention="non-identity 1-norm (Ha), at R = 2.760 bohr",
        q=60,
        provenance="MEASURED 2026-08-31, debug/qa/remeasure_balanced_lambda.py",
        aliases={}),
    "ash3_balanced_lambda": dict(
        value=885.248, convention="non-identity 1-norm (Ha), at R = 2.820 bohr",
        q=70,
        provenance="MEASURED 2026-08-31, debug/qa/remeasure_balanced_lambda.py",
        aliases={}),
    "geh4_balanced_lambda": dict(
        value=874.833, convention="non-identity 1-norm (Ha), at R = 2.870 bohr",
        q=80,
        provenance="MEASURED 2026-08-31, debug/qa/remeasure_balanced_lambda.py",
        aliases={}),
    "nah_balanced_qwc": dict(
        value=69, convention="QWC groups, balanced", q=20,
        provenance="MEASURED 2026-08-30; KH ties at 69 (the only "
                   "isostructural pair that does)",
        aliases={}),
    "mgh2_balanced_qwc": dict(
        value=903, convention="QWC groups, balanced", q=40,
        provenance="MEASURED 2026-08-30; CaH2 partner is 898, NOT equal",
        aliases={}),
    "sih4_balanced_qwc": dict(
        value=1551, convention="QWC groups, balanced", q=80,
        provenance="MEASURED 2026-08-30; GeH4 partner is 1562, NOT equal",
        aliases={}),

    # ---- composed / balanced 1-norms, exact global-M_L rule -------------
    # Every entry here has an identity-included twin; the alias records it so
    # a cell in the other convention is reported as a CONVENTION mismatch,
    # not as a wrong value.  Papers 14 and 20 were interleaving the two
    # inside single columns until 2026-08-30.
    "lih_composed_lambda": dict(
        value=27.715, convention="non-identity 1-norm (Ha)", q=30,
        provenance="MEASURED 2026-08-30, ecosystem hamiltonian('LiH')",
        aliases={34.036: "including identity", 34.0: "including identity"}),
    "beh2_composed_lambda": dict(
        value=54.286, convention="non-identity 1-norm (Ha)", q=50,
        provenance="MEASURED 2026-08-30, same route",
        aliases={68.812: "including identity", 68.8: "including identity"}),
    "h2o_composed_lambda": dict(
        value=186.761, convention="non-identity 1-norm (Ha)", q=70,
        provenance="MEASURED 2026-08-30, same route",
        aliases={372.464: "including identity", 372: "including identity"}),

    "lih_balanced_lambda": dict(
        value=75.231, convention="non-identity 1-norm (Ha)", q=30,
        provenance="MEASURED 2026-08-30, build_balanced_hamiltonian(lih_spec())",
        aliases={78.156: "including identity", 78.2: "including identity"}),
    "beh2_balanced_lambda": dict(
        value=286.893,
        convention="non-identity 1-norm (Ha), at BeH2's own R = 2.502 bohr",
        q=50,
        provenance="MEASURED 2026-09-01 after the R-default fix. The prior "
                   "289.581 / 328.11 pair was BeH2 evaluated at R = 3.015 "
                   "(LiH's bond length); see debug/qa/balanced_lambda_geometry_finding.md.",
        aliases={323.427: 'including identity', 323.4: 'including identity'}),
    "h2o_balanced_lambda": dict(
        value=1439.751, convention="non-identity 1-norm (Ha)", q=70,
        provenance="MEASURED 2026-08-30, same route",
        aliases={1612.010: "including identity", 1612: "including identity"}),

    "h2o_balanced_pauli": dict(
        value=19741, convention="non-identity Pauli terms", q=70,
        provenance="MEASURED 2026-08-30, same route",
        aliases={19742: "including identity"}),
    "beh2_balanced_pauli": dict(
        value=8867, convention="non-identity Pauli terms", q=50,
        provenance="MEASURED 2026-08-30, same route",
        aliases={8868: "including identity"}),
    "h2o_composed_pk_lambda": dict(
        value=28061.879,
        convention="1-norm including identity, PK in-Hamiltonian (Ha)",
        q=70,
        provenance="MEASURED 2026-08-30, build_balanced_hamiltonian "
                   "one_norm_composed (no non-identity split exposed)",
        aliases={}),
    "beh2_composed_pk_lambda": dict(
        value=374.866,
        convention="1-norm including identity, PK in-Hamiltonian (Ha)",
        q=50,
        provenance="MEASURED 2026-08-30, same route",
        aliases={}),

    # ---- angular ERI density (Paper 22), global-M_L vs pair-diagonal ----
    # Not a Pauli/1-norm quantity, but the same failure shape: two correct
    # densities under two rules, and the corpus quoted the stricter one as
    # "the" density in P26.  See C16 pairdiag-density-as-the-angular-density.
    "eri_density_lmax2": dict(
        value=8.52, convention="global-M_L Coulomb angular ERI density (%)",
        q=None,
        provenance="MEASURED 2026-08-30, potential_sparsity.angular_zero_count"
                   "(1, 2); exact-sympy pinned in tests/test_paper22_density.py",
        aliases={2.76: "pair-diagonal D_pd (axially-symmetric / m-decoupled "
                       "bases only)"}),

    # ---- He atomic series, exact global-M_L ERI rule --------------------
    "he_n2_pauli": dict(
        value=287, convention="non-identity Pauli terms", q=10,
        provenance="MEASURED 2026-08-30, JordanWignerEncoder on "
                   "LatticeIndex(He, max_n=2, slater_full)",
        aliases={288: "including identity"}),
    "he_n3_pauli": dict(
        value=14078, convention="non-identity Pauli terms", q=28,
        provenance="MEASURED 2026-08-30, same route",
        aliases={14079: "including identity"}),
    "he_n4_pauli": dict(
        value=250402, convention="non-identity Pauli terms", q=60,
        provenance="MEASURED 2026-08-30, same route",
        aliases={250403: "including identity"}),
    "he_n5_pauli": dict(
        value=2434441, convention="non-identity Pauli terms", q=110,
        provenance="MEASURED 2026-08-30, ecosystem hamiltonian('He', max_n=5); "
                   "28 min build",
        aliases={2434442: "including identity"}),

    "he_n2_lambda": dict(
        value=11.175, convention="non-identity 1-norm (Ha)", q=10,
        provenance="MEASURED 2026-08-30"),
    "he_n3_lambda": dict(
        value=74.207, convention="non-identity 1-norm (Ha)", q=28,
        provenance="MEASURED 2026-08-30"),
    "he_n4_lambda": dict(
        value=275.718, convention="non-identity 1-norm (Ha)", q=60,
        provenance="MEASURED 2026-08-30"),
    "he_n5_lambda": dict(
        value=790.007, convention="non-identity 1-norm (Ha)", q=110,
        provenance="MEASURED 2026-08-30"),

    "he_n2_qwc": dict(value=67, convention="greedy QWC groups", q=10,
                      provenance="MEASURED 2026-08-30"),
    "he_n3_qwc": dict(value=5569, convention="greedy QWC groups", q=28,
                      provenance="MEASURED 2026-08-30"),
    "he_n4_qwc": dict(value=86224, convention="greedy QWC groups", q=60,
                      provenance="MEASURED 2026-08-30, ~62 min"),

    # ---- fitted exponents.  The point range is part of the identity. ----
    "exp_pauli_4pt": dict(
        value=3.773, convention="log-log exponent, 4 points Q=10..110",
        provenance="MEASURED 2026-08-30; R^2 1.0000, max|log resid| 0.0064"),
    "exp_lambda_4pt": dict(
        value=1.774, convention="log-log exponent, 4 points Q=10..110",
        provenance="MEASURED 2026-08-30; R^2 0.9997, max|log resid| 0.0421"),
    "exp_qwc_3pt": dict(
        value=4.013, convention="log-log exponent, 3 points Q=10..60",
        provenance="MEASURED 2026-08-30; R^2 0.9976.  Q=110 infeasible: "
                   "2.4M terms against O(N^2) greedy grouping"),
    "exp_pauli_3pt": dict(
        value=3.779, convention="log-log exponent, 3 points Q=10..60",
        provenance="MEASURED 2026-08-30; the subset the slow tests compute"),
    "exp_lambda_3pt": dict(
        value=1.792, convention="log-log exponent, 3 points Q=10..60",
        provenance="MEASURED 2026-08-30"),
    "eri_decay_4pt": dict(
        value=-0.462, convention="ERI-density decay, 4 points M=5..55",
        provenance="MEASURED 2026-08-30 from the tab:eri_density column"),
    "eri_decay_2pt": dict(
        value=-0.489, convention="ERI-density decay, 2-point M=5->30 endpoint",
        provenance="MEASURED 2026-08-30.  This is the frequently-quoted "
                   "-0.49; it is NOT the four-point value"),

    # ---- isolated orbital blocks ----------------------------------------
    "block_sp_n2_pauli": dict(
        value=279, convention="non-identity Pauli, s+p block n_max=2 (M=5)",
        provenance="MEASURED 2026-08-29"),
    "block_sp_n2_eri": dict(
        value=107, convention="nonzero ERI quartets, s+p n_max=2 (M=5)",
        provenance="MEASURED 2026-08-29, two independent routes "
                   "(lattice_index and casimir_ci)"),
    "block_d_eri": dict(
        value=85, convention="nonzero ERI quartets, pure d shell (M=5)",
        provenance="MEASURED 2026-08-30; = sum of squared m-multiplicities"),
    "block_d_pauli": dict(
        value=343, convention="non-identity Pauli, pure d shell (M=5)",
        provenance="MEASURED 2026-08-30"),
    "block_sonly_pauli": dict(
        value=117, convention="non-identity Pauli, s-only M=3",
        provenance="MEASURED 2026-08-30",
        aliases={118: "including identity"}),
    "block_sp_m9_pauli": dict(
        value=2967, convention="non-identity Pauli, s+p M=9",
        provenance="MEASURED 2026-08-30"),

    # ---- composed architecture ------------------------------------------
    "composed_coeff": dict(
        value=27.90, convention="non-identity Pauli per qubit, main-group",
        provenance="MEASURED 2026-08-29; exact (= 279/10), not a fit"),
    "composed_coeff_d": dict(
        value=30.03, convention="non-identity Pauli per qubit, d-block",
        provenance="MEASURED 2026-08-29; DENSER than main-group"),
    "lih_composed_pauli": dict(
        value=837, convention="non-identity Pauli, LiH composed Q=30",
        provenance="MEASURED 2026-08-30",
        aliases={838: "including identity"}),
    "h2o_composed_pauli": dict(
        value=1953, convention="non-identity Pauli, H2O composed Q=70",
        provenance="MEASURED 2026-08-30",
        aliases={1954: "including identity"}),
    "exp_composed_3pt": dict(
        value=3.1685, convention="within-molecule exponent, 3 pts n_max=1..3",
        provenance="MEASURED 2026-08-30; the paper's small-basis 3.17"),
    "exp_composed_2pt": dict(
        value=2.8163, convention="within-molecule exponent, 2 pts n_max=1,2",
        provenance="FORCED: 1 + log_5(27.90/1.5); identical across all six "
                   "molecules, spread 0.000"),

    # ---- relativistic, post factor-order correction ----------------------
    "lih_rel_n2_pauli": dict(
        value=1501, convention="non-identity Pauli, LiH_rel Q=30",
        provenance="MEASURED 2026-08-30 under the corrected Condon-Shortley "
                   "order; retired order gave 1413"),
    "cah_rel_n2_pauli": dict(
        value=998, convention="non-identity Pauli, CaH/SrH/BaH_rel Q=20",
        provenance="MEASURED 2026-08-30; isostructural across all three"),
    "lih_rel_n2_lambda": dict(
        value=39.53, convention="non-identity 1-norm (Ha), LiH_rel Q=30",
        provenance="MEASURED 2026-08-30"),
    "lih_rel_n3_pauli": dict(
        value=90114, convention="non-identity Pauli, LiH_rel Q=84",
        provenance="MEASURED 2026-08-30"),

    # ---- balanced --------------------------------------------------------
    "lih_balanced_pauli": dict(
        value=2726, convention="Pauli terms, LiH balanced Q=30",
        provenance="MEASURED 2026-08-29; retired value 878"),
}

# ---------------------------------------------------------------------------
# Cited external baselines.  NOT ours; never silently re-measured.
# ---------------------------------------------------------------------------

# `identity=None` on a cited entry means UNKNOWN, not "not applicable":
# these are other groups' numbers and none of the sources state whether the
# identity term is counted.  Check E skips unknowns rather than assuming.
CITED = {
    "he_ccpvdz_pauli": dict(value=156, convention="Pauli, He cc-pVDZ Q=10",
                            source="recomputed Gaussian baseline, Paper 14", identity=None),
    "he_ccpvdz_lambda": dict(value=42.95, convention="1-norm (Ha)",
                             source="same", identity=None),
    "he_ccpvtz_pauli": dict(value=21607, convention="Pauli, He cc-pVTZ Q=28",
                            source="same", identity=None),
    "he_ccpvtz_lambda": dict(value=530.47, convention="1-norm (Ha)",
                             source="same", identity=None),
    "lih_sto3g_raw": dict(value=907, convention="raw JW Pauli, LiH STO-3G",
                          source="recomputed from published integrals", identity=None),
    "lih_sto3g_lambda": dict(value=34.3, convention="1-norm (Ha)",
                             source="same", identity=None),
    "h2o_sto3g_pauli": dict(value=551, convention="Pauli, H2O STO-3G Q=12",
                            source="same", identity=None),
    "lih_ccpvdz_pauli": dict(value=63519, convention="Pauli, LiH cc-pVDZ",
                             source="Trenev et al. 2025, Table 5", identity=None),
    "h2o_ccpvdz_pauli": dict(value=107382, convention="Pauli, H2O cc-pVDZ Q=46",
                             source="Trenev et al. 2025, Table 5", identity=None),
    "chawla_rah": dict(value=12556, convention="relativistic Pauli, RaH Q=18",
                       source="Chawla et al. 2024", identity=None),
}

# ---------------------------------------------------------------------------
# Derived quantities.  Stored as EXPRESSIONS, never as values -- so they
# cannot go stale independently of their inputs.  This is the class that
# survived seven review passes, because catching it needs recomputation
# rather than reading.
#
#   (name, expression over registry keys, tolerance, what it appears as)
# ---------------------------------------------------------------------------

DERIVED = {
    "he_n2_lambda_per_q":   ("he_n2_lambda / 10", 0.005, "lambda/Q at Q=10"),
    "he_n3_lambda_per_q":   ("he_n3_lambda / 28", 0.005, "lambda/Q at Q=28"),
    "he_n4_lambda_per_q":   ("he_n4_lambda / 60", 0.005, "lambda/Q at Q=60"),
    "he_n5_lambda_per_q":   ("he_n5_lambda / 110", 0.005, "lambda/Q at Q=110"),

    "ratio_lambda_q10":     ("he_ccpvdz_lambda / he_n2_lambda", 0.05,
                             "equal-qubit 1-norm advantage at Q=10"),
    "ratio_lambda_q28":     ("he_ccpvtz_lambda / he_n3_lambda", 0.05,
                             "equal-qubit 1-norm advantage at Q=28"),
    "ratio_pauli_q10":      ("he_n2_pauli / he_ccpvdz_pauli", 0.05,
                             "equal-qubit Pauli ratio at Q=10 (>1 = GeoVac denser)"),
    "ratio_pauli_q28":      ("he_ccpvtz_pauli / he_n3_pauli", 0.05,
                             "equal-qubit Pauli advantage at Q=28"),

    "ratio_ccpvdz_lih":     ("lih_ccpvdz_pauli / lih_composed_pauli", 1.0,
                             "LiH cc-pVDZ advantage"),
    "ratio_ccpvdz_h2o":     ("h2o_ccpvdz_pauli / h2o_composed_pauli", 1.0,
                             "H2O cc-pVDZ advantage"),

    "chawla_ratio_lih_n2":  ("lih_rel_n2_pauli / chawla_rah", 0.002,
                             "Chawla ratio, LiH_rel n_max=2"),
    "chawla_ratio_cah_n2":  ("cah_rel_n2_pauli / chawla_rah", 0.002,
                             "Chawla ratio, CaH_rel n_max=2"),

    "groups_per_term_n2":   ("he_n2_qwc / he_n2_pauli", 0.005,
                             "QWC groups per term at Q=10"),
    "groups_per_term_n3":   ("he_n3_qwc / he_n3_pauli", 0.005,
                             "QWC groups per term at Q=28"),
    "groups_per_term_n4":   ("he_n4_qwc / he_n4_pauli", 0.005,
                             "QWC groups per term at Q=60"),

    "composed_coeff_check": ("block_sp_n2_pauli / 10", 0.005,
                             "the composed coefficient IS the block count/10"),
    "lih_composed_check":   ("composed_coeff * 30", 1.0,
                             "LiH composed = 27.90 x 30"),
    "h2o_composed_check":   ("composed_coeff * 70", 1.0,
                             "H2O composed = 27.90 x 70"),

    "exponent_gap_lo":      ("3.9 - exp_pauli_4pt", 0.02,
                             "scaling gap vs the Gaussian band, low end"),
    "exponent_gap_hi":      ("4.3 - exp_pauli_4pt", 0.02,
                             "scaling gap vs the Gaussian band, high end"),
}

# ---------------------------------------------------------------------------
# Known-wrong values, kept so the gate can name WHAT a stale locus is.
# (C16/C17 block them; this maps them to their replacement for diagnostics.)
# ---------------------------------------------------------------------------

RETIRED = {
    # value: (replacement key, REQUIRED context, FORBIDDEN context or None)
    #
    # Context is mandatory.  Several of these numerals are common in other
    # roles (Z ranges, qubit counts, page and reference numbers), and a
    # bare-value match reports them all.  The forbidden column exists for
    # the genuinely ambiguous ones.
    120:     ("he_n2_pauli",        r"Pauli|terms",            r"page|Rev\.|\\textbf|N *= *120"),
    2659:    ("he_n3_pauli",        r"Pauli|terms",            None),
    31039:   ("he_n4_pauli",        r"Pauli|terms",            None),
    227338:  ("he_n5_pauli",        r"Pauli|terms",            None),

    # 1-norm / Pauli values the papers printed before the 2026-08-30 pass.
    5798:    ("h2o_balanced_pauli",   r"Pauli|terms",            None),
    1511:    ("h2o_balanced_lambda",  r"1-norm|lambda|\\lambda", None),
    306.4:   ("beh2_balanced_lambda", r"1-norm|lambda|\\lambda", None),
    373.4:   ("beh2_composed_pk_lambda", r"1-norm|lambda|\\lambda", None),
    354.9:   ("beh2_composed_pk_lambda", r"1-norm|lambda|\\lambda", None),
    66.0:    ("beh2_composed_lambda", r"electronic-only|1-norm",  None),
    28055:   ("h2o_composed_pk_lambda", r"1-norm|composed",       None),
    # Balanced lambda at the WRONG geometry: these are the values the
    # 2026-08-31 correction retired, produced by evaluating every
    # molecule at R = 3.015 (LiH's bond length).  The pre-2026-08-30
    # values 19.6 / 866.1 were RIGHT and are canonical again -- they used
    # to sit in this map, which is how a correct value came to be
    # flagged as retired.
    20.6:    ("nah_balanced_lambda",  r"lambda|1-norm|NaH",  r"page|Rev\.|retired|LiH's"),
    328.11:  ("beh2_balanced_lambda", r"lambda|1-norm|BeH",  r"page|Rev\.|retired"),
    328.1:   ("beh2_balanced_lambda", r"lambda|1-norm|BeH",  r"page|Rev\.|retired"),
    289.581: ("beh2_balanced_lambda", r"lambda|1-norm|BeH",  r"page|Rev\.|retired"),
    289.6:   ("beh2_balanced_lambda", r"lambda|1-norm|BeH",  r"page|Rev\.|retired"),
    869.7:   ("hcl_balanced_lambda",  r"lambda|1-norm|HCl",  r"page|Rev\.|retired|LiH's"),
    110.5:   ("mgh2_balanced_lambda", r"lambda|1-norm|MgH",  r"page|Rev\.|retired"),
    879.6:   ("h2s_balanced_lambda",  r"lambda|1-norm|H\$_2\$S",  r"page|Rev\.|retired"),
    895.3:   ("ph3_balanced_lambda",  r"lambda|1-norm|PH",   r"page|Rev\.|retired"),
    914.1:   ("sih4_balanced_lambda", r"lambda|1-norm|SiH",  r"page|Rev\.|retired"),
    # tab:tc_composed, pre-exact-rule (exempted for passes by a marker 5
    # lines below the table -- see stage4_ledger_notes.md).
    334:     ("lih_composed_pauli",  r"Pauli|TC|standard",      r"page|Rev\."),
    562:     ("lih_tc_pauli",        r"Pauli|TC",               r"page|Rev\."),
    934:     ("beh2_tc_pauli",       r"Pauli|TC",               None),
    1307:    ("h2o_tc_pauli",        r"Pauli|TC",               None),
    1.68:    ("tc_pauli_ratio",      r"TC|Pauli ratio",         None),

    11.29:   ("he_n2_lambda",       r"lambda|1-norm",          None),
    78.36:   ("he_n3_lambda",       r"lambda|1-norm",          None),
    261.57:  ("he_n4_lambda",       r"lambda|1-norm",          None),
    657.07:  ("he_n5_lambda",       r"lambda|1-norm",          None),

    25:      ("he_n2_qwc",          r"QWC|measurement group",  None),
    791:     ("he_n3_qwc",          r"QWC|measurement group",  None),
    10199:   ("he_n4_qwc",          r"QWC|measurement group",  None),

    # --- the retired pair-diagonal composed-scaling family (trunk FULL
    # 2026-09-01).  These three survived a certification because none was
    # registered; C21 samples what it knows.  2.5 is a ubiquitous numeral,
    # so its required context is the Pauli-scaling reading specifically;
    # the replacement is not another exponent but the exact linearity
    # N_Pauli = 27.90 x Q (composed_coeff).
    2.5:     ("composed_coeff",     r"Q\^\{2\.5\}|O\(Q|Pauli[- ]?(?:term )?scal", 
                                    r"retired|pair-diagonal (?:rule |count )?gave|superseded"),
    51:      ("composed_coeff",     r"1712|1\{,\}712",
                                    r"retired|pair-diagonal|gave"),
    1712:    ("composed_coeff",     r"advantage|fewer|Gaussian|Pauli",
                                    r"retired|pair-diagonal|gave|arXiv"),
    3.15:    ("exp_pauli_4pt",      r"Q\^|exponent|scaling|alpha", r"Table~3\.15|Table 3\.15"),
    1.69:    ("exp_lambda_4pt",     r"Q\^|exponent|scaling|alpha|lambda",
                                    r"meV|polarizability"),
    3.36:    ("exp_qwc_3pt",        r"Q\^|exponent|scaling|QWC",   None),
    3.147:   ("exp_pauli_4pt",      r"exponent|alpha|scaling",  None),
    1.690:   ("exp_lambda_4pt",     r"Q\^|exponent|scaling",
                                    r"meV|polarizability"),
    1.694:   ("exp_lambda_4pt",     r"exponent|alpha|scaling",  None),
    3.355:   ("exp_qwc_3pt",        r"exponent|alpha|scaling",  None),

    111:     ("block_sp_n2_pauli",  r"Pauli per|per \$s\$\+\$p\$|block.*Pauli", r"Rev\.|\\textbf|page"),
    # \b required: the context search is case-insensitive, so a bare
    # "ERI" matched inside charactERIzes / inhERIt and admitted any
    # 65 in the corpus (4 false positives, 2026-08-31).
    65:      ("block_sp_n2_eri",    r"\bERI\b|nonzero quartet", None),
    56:      ("block_d_pauli",      r"d-only|pure \$d\$.*Pauli", r"Z\s*=|Q\s*=|--56|through Ba"),
    976:     ("block_sp_m9_pauli",  r"Pauli",                   None),

    # The retired LiH-vs-cc-pVDZ Pauli ratio.  Canonical is 76x
    # (lih_ccpvdz_pauli / lih_composed_pauli = 63,519 / 837).  Added
    # 2026-08-31 after the Phase-0 sweep found it live in the field
    # guide (x2) and Paper 58 with NO registry entry to catch it --
    # it surfaced only via the S15 re-read-the-prose rule.  Context is
    # narrow on purpose: 190 is otherwise an unremarkable integer.
    190:     ("ratio_ccpvdz_lih",   r"cc-pVDZ|cc-pVTZ|Gaussian",
                                    r"page|Rev\.|Q = 190|qubits"),
    11.10:   ("composed_coeff",     r"coefficient|per qubit|\\times Q|Pauli/", None),
    9.23:    ("composed_coeff_d",   r"coefficient|per qubit|Pauli/", None),
    333:     ("lih_composed_pauli", r"Pauli",                   None),
    # (334 is defined once, above, with its forbidden context: a second
    #  334 entry here used to overwrite it and drop that guard.)
    778:     ("h2o_composed_pauli", r"Pauli",                   None),

    1413:    ("lih_rel_n2_pauli",   r"Pauli|N_\{\\mathrm\{Pauli|rel",  None),
    942:     ("cah_rel_n2_pauli",   r"Pauli|N_\{\\mathrm\{Pauli|rel",  None),
    40.59:   ("lih_rel_n2_lambda",  r"lambda|1-norm",           None),
    89226:   ("lih_rel_n3_pauli",   r"Pauli|rel",               None),
    878:     ("lih_balanced_pauli", r"Pauli|balanced",          None),
}



# ---------------------------------------------------------------------------
# Convention parsing.
#
# The `convention` strings are written for a human reader, but the identity
# convention has to be machine-checkable: the one defect class that pure
# value-checking cannot see is a table mixing identity-included and
# non-identity counts in the same column, where every individual value is
# correct and only the pairing is wrong.  `tests/test_numeric_registry.py`
# asserts that every registered convention parses, so a new entry phrased
# some other way fails loudly instead of quietly leaving the check.
# ---------------------------------------------------------------------------

_KINDS = (
    # Two-center decompactification FRONT (Papers 58/60, 2026-09-06): the
    # distance at which the two center projectors reach principal angle 45
    # deg, and quantities stated relative to it (ratios to a decay length,
    # percent residuals, coherence-front shifts).  Listed first so a front
    # convention that also mentions a decay length is not read as an exponent.
    ("front", ("front",)),
    # Dimensionless constants (rate constants, saturation constants): no
    # identity-in/out convention applies, like `density`.  Registered FULL #7
    # for the trunk multi-document constants (saturation_c, l2_rate_4_over_pi).
    ("constant", ("constant",)),
    ("pauli", ("pauli",)),
    ("lambda", ("1-norm",)),
    ("qwc", ("qwc",)),
    ("eri", ("eri quartet", "nonzero eri")),
    ("exponent", ("exponent", "decay")),
    # Angular ERI density: a percentage carrying a SELECTION-RULE convention
    # (global-M_L Coulomb vs the stricter pair-diagonal) rather than an
    # identity convention.  Registered 2026-08-30 after P26 quoted the
    # pair-diagonal 2.76% as "the" density; family() leaves identity None
    # here, so check E compares rules instead of identity terms.
    ("density", ("density",)),
)


def family(key: str):
    """(kind, identity) for a registered quantity.

    kind      -- 'pauli' | 'lambda' | 'qwc' | 'eri' | 'exponent' | None
    identity  -- 'out' (non-identity), 'in' (includes identity), or None
                 where the distinction does not apply.
    """
    d = MEASURED.get(key) or CITED.get(key)
    if not d:
        raise KeyError(key)
    conv = (d.get("convention") or "").lower()
    kind = None
    for name, needles in _KINDS:
        if any(n in conv for n in needles):
            kind = name
            break
    # An explicit field wins over the parse.  External baselines usually
    # do not state their convention, and inferring one from silence would
    # assert something the source does not support (registry rule 3).
    if "identity" in d:
        return kind, d["identity"]
    identity = None
    if kind in ("pauli", "lambda"):
        identity = "out" if "non-identity" in conv else "in"
    return kind, identity


def parses(key: str) -> bool:
    """True if the convention string yields a recognised kind."""
    return family(key)[0] is not None


def resolve(key: str, _seen=None) -> float:
    """Value of a measured, cited, or DERIVED quantity.

    DERIVED was added to this lookup 2026-08-31.  Without it a RETIRED
    entry could not name a derived replacement, so retired RATIOS,
    exponents and per-qubit figures had nowhere to point -- and a ratio
    is inherently derived.  The 190x LiH-vs-cc-pVDZ figure sat live in
    three loci with no registry entry able to catch it for exactly this
    reason.
    """
    if key in MEASURED:
        return float(MEASURED[key]["value"])
    if key in CITED:
        return float(CITED[key]["value"])
    if key in DERIVED:
        # Guard against a cyclic expression rather than blowing the stack.
        _seen = set(_seen or ())
        if key in _seen:
            raise KeyError(f"cyclic derivation through {key!r}")
        _seen.add(key)
        return evaluate(DERIVED[key][0], _seen)
    raise KeyError(f"unregistered quantity: {key}")


def evaluate(expr: str, _seen=None) -> float:
    """Evaluate a DERIVED expression over registry keys."""
    import re as _re
    names = sorted(set(_re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expr)),
                   key=len, reverse=True)
    out = expr
    for nm in names:
        out = out.replace(nm, repr(resolve(nm, _seen)))
    return float(eval(out, {"__builtins__": {}}, {}))


def all_canonical_values():
    """Every value the corpus is allowed to state, with its label."""
    out = {}
    for k, d in MEASURED.items():
        out.setdefault(round(float(d["value"]), 6), []).append(k)
        for alias in (d.get("aliases") or {}):
            out.setdefault(round(float(alias), 6), []).append(f"{k} (alias)")
    for k, d in CITED.items():
        out.setdefault(round(float(d["value"]), 6), []).append(k)
    return out
