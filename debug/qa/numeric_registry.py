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
        value=20.607, convention="non-identity 1-norm (Ha)", q=20,
        provenance="MEASURED 2026-08-30, build_balanced_hamiltonian(nah_spec()); "
                   "reproducible across two builds. The 2026-08-29 re-sync "
                   "left 19.6 live in BOTH P14 tab:second_row and P20 "
                   "tab:molecules.",
        aliases={193.440: "including identity"}),
    "hcl_balanced_lambda": dict(
        value=869.720, convention="non-identity 1-norm (Ha)", q=50,
        provenance="MEASURED 2026-08-30, same route; retired 866.1 was live "
                   "in the same two tables",
        aliases={1239.441: "including identity"}),
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
        value=289.581, convention="non-identity 1-norm (Ha)", q=50,
        provenance="MEASURED 2026-08-30, same route",
        aliases={328.110: "including identity", 328.1: "including identity"}),
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
    120:     ("he_n2_pauli",        r"Pauli|terms",            r"page|Rev\.|\\textbf"),
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
    19.6:    ("nah_balanced_lambda",  r"lambda|1-norm|NaH",      r"page|Rev\."),
    866.1:   ("hcl_balanced_lambda",  r"lambda|1-norm|HCl",      None),
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

    3.15:    ("exp_pauli_4pt",      r"Q\^|exponent|scaling|alpha", r"Table~3\.15|Table 3\.15"),
    1.69:    ("exp_lambda_4pt",     r"Q\^|exponent|scaling|alpha|lambda", None),
    3.36:    ("exp_qwc_3pt",        r"Q\^|exponent|scaling|QWC",   None),
    3.147:   ("exp_pauli_4pt",      r"exponent|alpha|scaling",  None),
    1.694:   ("exp_lambda_4pt",     r"exponent|alpha|scaling",  None),
    3.355:   ("exp_qwc_3pt",        r"exponent|alpha|scaling",  None),

    111:     ("block_sp_n2_pauli",  r"Pauli per|per \$s\$\+\$p\$|block.*Pauli", r"Rev\.|\\textbf|page"),
    65:      ("block_sp_n2_eri",    r"ERI|nonzero quartet",     None),
    56:      ("block_d_pauli",      r"d-only|pure \$d\$.*Pauli", r"Z\s*=|Q\s*=|--56|through Ba"),
    976:     ("block_sp_m9_pauli",  r"Pauli",                   None),

    11.10:   ("composed_coeff",     r"coefficient|per qubit|\\times Q|Pauli/", None),
    9.23:    ("composed_coeff_d",   r"coefficient|per qubit|Pauli/", None),
    333:     ("lih_composed_pauli", r"Pauli",                   None),
    334:     ("lih_composed_pauli", r"Pauli",                   None),
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


def resolve(key: str) -> float:
    """Value of a measured or cited quantity."""
    if key in MEASURED:
        return float(MEASURED[key]["value"])
    if key in CITED:
        return float(CITED[key]["value"])
    raise KeyError(f"unregistered quantity: {key}")


def evaluate(expr: str) -> float:
    """Evaluate a DERIVED expression over registry keys."""
    import re as _re
    names = sorted(set(_re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expr)),
                   key=len, reverse=True)
    out = expr
    for nm in names:
        out = out.replace(nm, repr(resolve(nm)))
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
