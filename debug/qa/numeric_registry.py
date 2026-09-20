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
    # ------------------------------------------------------------- Papers 59/61
    "t2_collinear_period": dict(
        value=0.395355765901713964,
        convention="constant: the collinear three-center Bessel-moment period "
                   "T2, certified to 66 digits (stored here as a double, so "
                   "only the leading ~17 are representable -- the full digit "
                   "string is in provenance)",
        q=None,
        provenance="MEASURED/certified (Paper 61 Sec. 3): six parameter-disjoint "
                   "runs, two independent parallel configurations agreeing to "
                   "1.7e-67, identity cross-validated against an independent 2D "
                   "evaluation to 69-96 digits. Full value "
                   "0.395355765901713964325229296804847564260563977867082108"
                   "935234265469... Pinned by tests/test_paper59_t2_value.py "
                   "(identity, K^-7 tail law, and the digit-19 correction of "
                   "the earlier 18-digit anchor).",
        # No alias: the paper's 18-digit display and its 66-digit certified
        # display are the SAME double, so C21's display-rounding acceptance
        # matches both against this one value.  (The 19th digit of the older
        # 18-digit anchor was superseded by the certification -- a fact about
        # the digits beyond double precision, not about this stored value.)
        aliases={}),
    "three_center_eri_truth": dict(
        value=0.204941722,
        convention="constant: ground-truth value of the genuine three-center "
                   "ERI (XY|XZ) over 1s Slater orbitals, Paper 59 validation "
                   "point",
        q=None,
        provenance="MEASURED (Paper 59 Sec. 4), the reference the closed-form "
                   "route is checked against",
        aliases={}),
    # ---------------------------------------------------------------- Paper 60
    # Isoenergetic secular matrix M = diag(Z R_nu) + T'.  The 1-norm exponent is
    # a WINDOW fit and the window is load-bearing: /qa 2026-09-07 graded "grows
    # sublinearly" a LARGE overclaim because the local slope rises steadily past
    # the fitted range and the sublinear part is the nuclear diagonal, not the
    # pure-number block the paper credited.  Both facts are pinned here so a
    # future edit that moves the number is forced back through the prose.
    "p60_l2_inflation_hydrogenic": dict(
        value=1.19, convention="exponent: log-log fit of the naive-Loewdin "
                               "block-encoding 1-norm lambda vs Q=2N for the "
                               "HYDROGENIC (a=Z/n) atomic basis, N=1..5 (Q=2..10), "
                               "N=1 excluded (eq:blowup)",
        q=None,
        provenance="MEASURED 2026-09-13 via "
                   "geovac.sturmian_l2_encoding.fit_lambda_exponent(hydrogenic): "
                   "1.1909. Test-banded (1.0,1.4) in "
                   "tests/test_sturmian_l2_encoding.py::"
                   "test_paper60_l2_encoding_blowup_exponents."),
    "p60_l2_inflation_sturmian": dict(
        value=3.33, convention="exponent: same fit as p60_l2_inflation_hydrogenic "
                               "but the shared-scale STURMIAN basis; eq:blowup",
        q=None,
        provenance="MEASURED 2026-09-13 via fit_lambda_exponent(sturmian): 3.3328. "
                   "Test-banded (3.0,3.6), same test. The shared-scale set inflates "
                   "faster than hydrogenic by >1 in exponent -- the eq:blowup "
                   "mechanism (density of S^-1/2)."),
    "p60_sw_cond_exponent": dict(
        value=1.85, convention="exponent: log-log fit of cond of the momentum-space "
                               "Shibuya-Wulfman two-center metric vs N=2n, WINDOW "
                               "N=12..16 at R=2 bohr. A PRE-ASYMPTOTIC window reading, "
                               "NOT canonical: the full N=4..16 slope is 1.80 and the "
                               "DERIVED asymptote is N^2 (eq:sigma_law)",
        q=None,
        provenance="MEASURED 2026-09-13 via tests/test_paper60_sturmian.py::"
                   "_sw_metric_mom: 1.8518 on N=12..16 (full-range 1.8003). The SW "
                   "exponent is test-banded (1.5,2.2) as ~N^1.8 in "
                   "test_paper60_sw_better_conditioned_than_l2."),
    "p60_l2_overlap_exponent": dict(
        value=1.71, convention="exponent: log-log fit of cond of the shared-scale "
                               "Coulomb-Sturmian L2 overlap vs N=2,3,5,8 (the "
                               "3.0/5.83/13.93/32.16 sequence). Comparator to "
                               "p60_sw_cond_exponent",
        q=None,
        provenance="MEASURED 2026-09-13: 1.7114 on N=2,3,5,8. REGISTERED BECAUSE THE "
                   "PAPER PRINTED 1.70 (an approximate comparator); corrected to the "
                   "measured 1.71 at its one locus. Sequence pinned in "
                   "tests/test_paper60_sturmian.py::"
                   "test_paper60_l2_overlap_condition_number_grows."),
    "p60_water_a1_exponent": dict(
        value=1.96, convention="exponent: log-log fit of cond of the RAW water A_1 "
                               "block vs N=2n over N=12..192, via _water_A1. This is "
                               "the 'N^1.96 here' value of the sec:resource "
                               "preconditioner table; DISTINCT from the sec:molecular "
                               "three-center-SW probe (N^1.97, 19.9->698 over N=6..36, "
                               "a different construction)",
        q=None,
        provenance="MEASURED 2026-09-13 via tests/test_paper60_preconditioner.py::"
                   "_water_A1: 1.9580 over N=12..192 (cond 183/2696/41700). The raw "
                   "exponent is banded (1.85,2.05) in "
                   "test_uniform_banding_alone_is_not_the_lever."),
    "p60_onenorm_exponent": dict(
        value=0.82, convention="exponent: log-log slope of the entrywise norm of M "
                               "vs K, full s+p+d+f, helium (Z=2), over the window "
                               "K = 74..164, on a CONVERGED radial box; NOT "
                               "asymptotic and NOT a stable exponent",
        q=None,
        provenance="MEASURED 2026-09-07 by TWO independent routes agreeing to 4 dp: "
                   "exact grid-free Slater algebra (dps 60-70) and converged "
                   "quadrature under R_MAX >= 3*n_max^2. Both give 0.8193 on this "
                   "window. The local slope FALLS monotonically -- 0.827, 0.818, "
                   "0.810, 0.802, 0.794, 0.787, 0.781, 0.776, 0.771, 0.766 across "
                   "K = 100..514 -- so no window fit is stable. RETIRED: 0.84 and "
                   "0.854 and the 'rises to 0.906' sequence, all measured on a "
                   "60-bohr box whose relative error GREW x1.011 -> x1.205 across "
                   "the fit range; that growth, not the matrix, produced the rise.",
        aliases={0.8058: "global fit extended to K = 340",
                 0.7976: "global fit extended to K = 514",
                 0.766: "local slope at K = 514"}),
    "p60_onenorm_exponent_sonly": dict(
        value=0.73, convention="exponent: same slope, s-sector configurations only, "
                               "helium, converged box, K = 74..164-equivalent "
                               "s-window", q=None,
        provenance="MEASURED 2026-09-07, converged. RETIRED: 0.78, which was "
                   "measured on a DIFFERENT window (K = 3..21, "
                   "debug/sturmian_he_secular.py) and registered here as 'same "
                   "window' -- that provenance was false. On the paper's own "
                   "n_max = 7..10 window the s-only/spdf pair is 0.726 -> 0.819, "
                   "not the printed 0.78 -> 0.84. The DIRECTION of the "
                   "angular-dilution statement survives; both endpoints move.",
        aliases={0.7154: "s-only, matched n_max 7..14 window"}),
    "p60_T0_asymptotic_exponent": dict(
        value=0.5, convention="exponent: ASYMPTOTIC log-log slope of the nuclear "
                              "diagonal ||T0||_1 = Z*sum(R_nu). Not a fitted value "
                              "-- ||T0||_1 is NOT a power law",
        q=None,
        provenance="SYMBOLIC 2026-09-07. ||T0||_1 = Z*sqrt(2K)*ln(K/2) for lmax=3, "
                   "with K = 2N^2-4N+4 exactly. Needs no ERI and no grid -- pure "
                   "configuration combinatorics, so it is box-independent and "
                   "identical on every route. PM verified the slope to K = 498,004: "
                   "0.709 (K=20) -> 0.705 (164) -> 0.699 (340) -> 0.662 (12484) -> "
                   "0.603 (498004), still falling, consistent with 1/2 + O(1/log K). "
                   "RETIRED as a claim: the fitted 'K^0.70', which is simply what a "
                   "log-log fit returns inside K = 74..164.",
        aliases={0.7050: "fitted slope inside the paper's K = 74..164 window",
                 0.7026: "fitted slope over K = 74..340"}),
    "p60_Tprime_exponent": dict(
        value=0.94, convention="exponent: log-log slope of the OFF-diagonal "
                               "entrywise norm of M (= offdiag of T') vs K, helium, "
                               "converged box, K = 74..340",
        q=None,
        provenance="MEASURED 2026-09-07, two routes agreeing to 4 dp (0.9368). "
                   "SUBLINEAR. RETIRED: 1.05, which carried +0.07 of pure box bias. "
                   "NOTE the paper's own split is M = diag(Z R_nu) + T', so the "
                   "block its prose names is the FULL T' -- see "
                   "p60_Tprime_full_exponent -- not this off-diagonal part.",
        aliases={0.9774: "same leg on the K = 74..164 window",
                 0.9156: "extended to K = 514"}),
    "p60_Tprime_full_exponent": dict(
        value=0.88, convention="exponent: log-log slope of the FULL ||T'||_1 (the "
                               "block the paper's prose names), helium, converged "
                               "box, K = 74..340",
        q=None,
        provenance="MEASURED 2026-09-07, two routes agreeing (0.8750). SUBLINEAR "
                   "over the entire range; its local slope crosses 1 at K ~ 37 and "
                   "falls to 0.805 by K = 514. This is the object 'the pure-number "
                   "block T' is superlinear' referred to, and that claim is FALSE. "
                   "The exact additive identity is "
                   "||M||_1 = ||T0||_1 - ||diag T'||_1 + ||T'||^off, and the paper's "
                   "eq:sublinear_split printed the first and third terms while "
                   "dropping the middle one.",
        aliases={0.9033: "K = 74..164 window", 0.8597: "extended to K = 514"}),
    "p60_energy_floor": dict(
        value=6.44, convention="constant: mHa by which the fixed-l_max s+p+d+f "
                               "Goscinskian family saturates ABOVE the exact "
                               "helium ground state; a floor, not a rate",
        q=None,
        provenance="MECHANISM (2026-09-08): the floor is the SCALE LOCK, not angular "
                   "truncation -- freeing lambda over the identical span "
                   "reaches 1.28 mHa at K=130. And it is GROUND-STATE "
                   "specific: 2^1S sits at 1.786 mHa at K=202. "
                   "MEASURED 2026-09-07 on the converged ladder K = 74..514, fit "
                   "gap(K) = floor + b*K^-q (floor 6.4404, b 49.0, q 0.796; max "
                   "residual 0.0019 mHa). Extrapolates to 6.4455 at K = 1e5 and "
                   "6.4412 at K = 1e6 -- genuinely saturating, not slow "
                   "convergence. 7x the basis (K 74 -> 514) closed 16% of the "
                   "deficit. This is 4.0x chemical accuracy (1.594 mHa), so the "
                   "family cannot reach chemical accuracy at ANY K.",
        aliases={6.78: "gap at K = 514, the largest computed"}),
    "p60_span_deficit_spdf": dict(
        value=1.28, convention="constant: mHa above the exact He ground state "
                               "reached by a VARIATIONAL CI over the identical "
                               "Goscinskian span (l_max=3, K=130) with the global "
                               "scale lambda optimized -- the span's own deficit, "
                               "with the isoenergetic scale-lock removed",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_scale_scan.py. Companion "
                   "locked-scale value at the same K is 7.46 mHa. s-sector "
                   "counterpart (against the known exact s-limit -2.879029 Ha) is "
                   "p60_span_deficit_sonly_free vs p60_span_deficit_sonly_locked, "
                   "both at K=136 (this prose said '4.43 locked' until 2026-09-11, "
                   "which is the K=105 row -- a mismatched pair). Pipeline unit-tested at K=1, "
                   "where both postings coincide and return -2.8476562 = "
                   "-(2-5/16)^2 exactly.",
        aliases={1.640: "l_max=3, K=100", 2.232: "l_max=3, K=74"}),
    "p60_span_deficit_sonly_locked": dict(
        value=4.4038, convention="constant: mHa above the exact He s-limit "
                                 "(-2.879029 Ha) reached by the LOCKED metric-free "
                                 "isoenergetic posing, s-only, K=136 -- the partner "
                                 "of p60_span_deficit_sonly_free at the SAME K",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_freescale_resource.py; read back "
                   "from debug/data/p60_freescale_l0.json 2026-09-11. Registered "
                   "because the paper carried 4.43 here -- the K=105 row -- for "
                   "three days as an unregistered literal C21 could not see. The "
                   "two halves of this comparison are a MATCHED PAIR and must move "
                   "together: quoting one at K=136 and the other at K=105 is the "
                   "defect this key exists to prevent.",
        aliases={4.4341: "K=105", 4.4805: "K=78", 4.5565: "K=55"}),
    "p60_span_deficit_sonly_free": dict(
        value=0.1466, convention="constant: mHa above the exact He s-limit "
                                 "(-2.879029 Ha) reached by a VARIATIONAL CI over "
                                 "the identical s-only Goscinskian span with the "
                                 "global scale lambda optimized, K=136",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_freescale_resource.py; read back "
                   "from debug/data/p60_freescale_l0.json 2026-09-11. Partner of "
                   "p60_span_deficit_sonly_locked at the same K.",
        aliases={0.2219: "K=105", 0.3555: "K=78", 0.6135: "K=55"}),
    "p60_posing_cost_ground": dict(
        value=4.21, convention="constant: mHa, E_iso - min_lambda E_var over the "
                               "SAME span, He ground state, s-only, K=105. "
                               "Reference-free: needs no known limit",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_posing_cost_by_state.py, "
                   "lambda by GRID scan (Brent found a local minimum at nmax=4 and "
                   "reported a NEGATIVE cost, which the variational bound forbids). "
                   "Grows with K: 3.52 (K=36), 4.13 (K=78), 4.21 (K=105).",
        aliases={3.524: "K=36", 4.125: "K=78"}),
    "p60_posing_cost_exc": dict(
        value=0.98, convention="constant: mHa, same quantity as "
                               "p60_posing_cost_ground but for He 2^1S (the second "
                               "root of M), s-only, K=105",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_posing_cost_by_state.py. "
                   "4.3x SMALLER than the ground state. The ratio NARROWS with K "
                   "(5.17 at K=36, 5.11 at 55, 4.66 at 78, 4.29 at 105); "
                   "what widens is the absolute separation, 2.84 -> 3.23 mHa, "
                   "and that flattens by K=105. RETIRED: 'the ratio widens "
                   "with K', written 2026-09-08 and refuted by the backing "
                   "test the same day. "
                   "(0.681/3.524 at K=36 -> 0.983/4.212 at K=105). This is the "
                   "measured form of Avery's split-shell mechanism: a Goscinskian "
                   "1s^2 configuration pins both electrons to one exponent, an "
                   "excited configuration gets two free from n_a != n_b.",
        aliases={0.681: "K=36", 0.885: "K=78"}),
    "p60_exc_gap_k452": dict(
        value=1.7163, convention="constant: mHa above the exact He 2^1S energy "
                                 "(-2.145974046 Ha) reached by the METRIC-FREE "
                                 "isoenergetic posing, full s+p+d+f, at the largest basis "
                                 "where BOTH roots were computed, K=452 (n_max=16; the ground-state-only ladder reaches K=514). A MEASURED "
                                 "ladder endpoint, deliberately NOT an extrapolated "
                                 "floor -- see the caveat below",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_excited_ladder.py, extended to "
                   "n_max=16. Companion ground-state value at the SAME K and the "
                   "same ||M||_1 is p60_gnd_gap_k452 = 6.8196 mHa. Ladder: 1.972 "
                   "(K=74), 1.813 (164), 1.786 (202), 1.7485 (290), 1.7163 (452). "
                   "CAVEAT ON THE FLOOR: a free-floor fit is WINDOW-stable (0.4% "
                   "drift over 21 windows, drifting UP, i.e. approaching from "
                   "below) and Shanks brackets it from above at [1.647, 1.676], "
                   "but it is MODEL-family sensitive -- a two-parameter c + b/lnK "
                   "form, rejected at 340x worse RMS, puts the floor at 0.82x "
                   "chemical accuracy, i.e. BELOW it. The paper therefore cites "
                   "this measured endpoint and states that 'above chemical "
                   "accuracy' is a model-selection conclusion.",
        aliases={1.786: "K=202", 1.647: "Shanks lower bracket on the floor"}),
    "chem_accuracy_mha": dict(
        value=1.5936014616, convention="constant: chemical accuracy = 1 kcal/mol "
                                       "expressed in mHa (4.184 kJ/mol divided by "
                                       "2625.4996 kJ/mol per Hartree)",
        q=None,
        provenance="DEFINITION, CODATA-consistent unit conversion. Registered "
                   "2026-09-08 because Paper 60 states three accuracies as "
                   "MULTIPLES of it, and those multiples must be derived from the "
                   "energies rather than typed independently.",
        aliases={1.594: "3 dp, as printed in Sec.4"}),
    "p60_gnd_gap_k452": dict(
        value=6.8196, convention="constant: mHa above the exact He ground state "
                                 "(-2.903724377 Ha) reached by the METRIC-FREE "
                                 "isoenergetic posing, full s+p+d+f, at the largest basis "
                                 "where BOTH roots were computed, K=452 (n_max=16; the ground-state-only ladder reaches K=514)",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_excited_ladder.py. Paired with "
                   "p60_exc_gap_k452 at the SAME K and the same ||M||_1 -- the pair "
                   "is the paper's state-dependence headline, so both are "
                   "registered and both ratios are DERIVED from them rather than "
                   "typed. Shanks brackets this floor from above at [6.47, 6.62].",
        aliases={8.036: "K=74", 7.158: "K=202", 6.9143: "K=340"}),
    "p60_he_chain_s_k55": dict(
        value=-2.8745, convention="constant: Ha, the He ground-state energy from "
                                  "the LOCKED metric-free isoenergetic posing on "
                                  "the s-only (l_max=0) Goscinskian family at "
                                  "n_max=10, K=55 -- the second rung of the "
                                  "convergence chain. Printed rounded to four "
                                  "decimals",
        q=None,
        provenance="MEASURED 2026-09-12 on a converged grid (box 500, 40000 pts): "
                   "-2.874468, i.e. 29.256 mHa above the exact -2.903724377. "
                   "REGISTERED BECAUSE THE PAPER CARRIED -2.873, which is the "
                   "n_max=4 (K=10) value (-2.873219) -- so the chain silently "
                   "MIXED basis sizes, against this paper's own requirement that "
                   "every quantity name its basis-growth family. Found by the "
                   "owed-items recheck after /qa paper_60 FULL 2026-09-12, which "
                   "had confirmed only the K=164 endpoint. Grid note: an apparent "
                   "box drift (-2.8744 -> -2.8739 over boxes 300..1200) is a GRID "
                   "artifact; at 40000 points the value is stable to 2e-5 across "
                   "the same boxes.",
        aliases={-2.873219: "n_max=4, K=10 -- the value the paper had printed",
                 -2.894672: "the +p rung, n_max=10, K=100",
                 -2.847651: "the 1s^2 single-configuration rung"}),
    "p60_he_chain_spdf_k164": dict(
        value=-2.8964, convention="constant: Ha, the He ground-state energy "
                                  "reached by the LOCKED metric-free isoenergetic "
                                  "posing on the full s+p+d+f Goscinskian family "
                                  "at K=164 -- the endpoint of the convergence "
                                  "chain quoted in the abstract and in sec:atomic. "
                                  "Printed ROUNDED to four decimals as -2.8964 "
                                  "(corrected 2026-09-13: this line said THREE "
                                  "decimals, and a 2026-09-12 edit asserted the "
                                  "chain was truncated -- it is rounded, as the "
                                  "exact 1s^2 value -(27/16)^2 = -2.84765625 "
                                  "shows, truncating to -2.8476 where the paper "
                                  "prints -2.8477)",
        q=None,
        provenance="DERIVED 2026-09-12 from the registered gap at the same K: "
                   "the extended-ladder alias gives 7.289 mHa above the exact "
                   "-2.903724377 Ha, so E = -2.896435. REGISTERED BECAUSE THE "
                   "PAPER CARRIED -2.897 -- a retired 60-bohr-domain value -- at "
                   "TWO loci (abstract and sec:atomic) as an unregistered literal "
                   "C21 could not see. It was also self-refuting: K=164 is a "
                   "nested sub-family of the K=244 pool whose own error is 7.06 "
                   "mHa, so eq:no_selection's interlacing forces E(164) >= "
                   "-2.896667, which -2.897 violates. Found by /qa paper_60 FULL "
                   "2026-09-12. The direct grid measurement gives -2.896432 "
                   "against the -2.896435 derived from the gap alias, a 3e-6 "
                   "spread well inside the printed precision. OWED DISCHARGED "
                   "2026-09-12/13: both earlier rungs were remeasured -- the s "
                   "rung WAS stale (-2.873 was the n_max=4 value; it is now "
                   "registered as p60_he_chain_s_k55 = -2.8745) and the +p rung "
                   "was confirmed at -2.894672.",
        aliases={-2.895688: "K=74, from the 8.036 mHa alias"}),
    "p60_window_richardson_pi2": dict(
        value=9.86949, convention="constant: first-order Richardson limit in 1/n "
                                  "of n^2 * min<theta^2>, the minimal mean-square "
                                  "spread in the Fock polar angle theta = pi - chi "
                                  "over span{sin(a chi)}_{a<=n}; target pi^2 = "
                                  "9.8696044",
        q=None,
        provenance="MEASURED 2026-09-12, debug/p60_window_constant_probe.py test F2, "
                   "from n = 20..640. This is the Kac-Murdock-Szego constant c_1 "
                   "recomputed in a SECOND representation -- a band-limited "
                   "concentration problem with NO Bessel function anywhere in it -- "
                   "which is the whole evidential point: it shows the pi^2 of "
                   "eq:sigma_law is a truncation constant, not Bessel content. "
                   "First route was the exact tridiagonal spectrum (v5.10.18, "
                   "tests/test_paper60_kms_attribution.py). Two routes, one "
                   "constant, per the independent-route rule.",
        aliases={9.84287: "raw n=640, unextrapolated",
                 9.81624: "raw n=320", 9.76331: "raw n=160"}),
    "p60_window_rms_richardson_pi": dict(
        value=3.14158, convention="constant: first-order Richardson limit in 1/n of "
                                  "n * rms(theta) for the near-null (top singular) "
                                  "direction of the SW cross block at kR = 2; "
                                  "target pi = 3.1415927",
        q=None,
        provenance="MEASURED 2026-09-12, debug/p60_window_constant_probe.py test F1, "
                   "from n = 20..640. Says the degeneracy direction ACHIEVES the "
                   "band-limited minimum of p60_window_richardson_pi2 (same "
                   "constant, squared), i.e. the near-null direction is the optimal "
                   "concentrator at the p = 0 pole. Paired with F3, which checks "
                   "1 - sigma_max = (kR)^2 <theta^2>/24 on that direction to 0.3% "
                   "at n = 320 across kR = 0.5..4.",
        aliases={3.13733: "raw n=640, unextrapolated", 3.13309: "raw n=320"}),
    "p60_weighted_collapse_control": dict(
        value=0.828, convention="constant: the collapse (1-sigma_max)(n/kR)^2 at "
                                "kR=2, n=160 under the CONTROL weight W = 1+cos(chi), "
                                "which VANISHES at the degeneracy chi = pi -- the "
                                "value that must differ from pi^2/24 = 0.4112 for "
                                "the weight-independence measurement to mean "
                                "anything",
        q=None,
        provenance="MEASURED 2026-09-12, debug/p60_contraction_seam_probe.py test E. "
                   "MATCHED SET -- this control and the four smooth-weight values in "
                   "aliases must move together; quoting the agreement without the "
                   "control would report an insensitivity as a measurement. Fitted "
                   "exponent stays -1.97 here, so the control moves the CONSTANT "
                   "only, not the exponent.",
        aliases={0.40718: "smooth W = 1 (SW reference), n=160",
                 0.41229: "smooth W = 1 + 0.8 cos(chi), n=160",
                 0.41065: "smooth W = 2 + sin(chi), n=160",
                 0.41635: "smooth W = exp(-chi), n=160"}),
    "p60_stateprep_overlap_exc": dict(
        value=0.798, convention="constant: L2-metric overlap between the normalized "
                                "dominant single configuration and the true 2^1S "
                                "root, K=164 -- the state-preparation cost driver "
                                "for an interior root",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_stateprep_overlap.py. The four "
                   "lowest roots give 0.992, 0.798, 0.864, 0.889 -- 2^1S is the "
                   "HARDEST of the four, not the deepest, so the driver is mixing "
                   "at the bottom of the Rydberg series and not spectral depth. "
                   "Costs 1.25x the ground state in rotations, 1.54x in "
                   "repetitions; two configurations reach 0.99.",
        aliases={0.992: "ground state", 0.6371: "overlap squared, 2^1S"}),
    "p60_best102_locked": dict(
        value=7.25, convention="constant: mHa above the exact He ground state "
                               "reached by the BEST 102 configurations (ranked by "
                               "ground-state weight) drawn from the K=244 pool, "
                               "locked-scale posing -- the sharpest selection test "
                               "of Avery's '102 optimized configurations'",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_avery_102_probe.py. ABOVE the "
                   "pool's own 7.057 mHa, as Cauchy interlacing requires: a "
                   "principal submatrix cannot have a larger top eigenvalue. 200 "
                   "random 102-subsets reach 13.9 mHa at best. The cited Avery figure "
                   "is 1.224 mHa, unreachable in this posing at any K.",
        aliases={7.057: "the full K=244 pool",
                 13.862: "BEST of 200 random 102-subsets",
                 903.4: "WORST of 200 random 102-subsets -- the driver's max() over\n                         energies selected the least-bound subset;  mislabelled as\n                         'best' until 2026-09-08"}),
    "p60_freescale_set_sonly": dict(
        value=0.72, convention="exponent: log-log slope of the LOCKED ||M||_1 vs K "
                               "on the S-ONLY free-scale comparison ladder, "
                               "K = 21..136 (n_max 6..16). NOT the headline "
                               "eq:sublinear exponent -- that is 0.82, full "
                               "s+p+d+f, window K=74..164. Different sector, "
                               "different window, different ladder",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_freescale_resource.py. This is "
                   "the FIRST member of a matched SET measured on one ladder, and "
                   "the set must move together: locked ||M||_1 K^0.716, "
                   "||H(lambda*)||_1 K^1.952, whitened ||S^-1/2 H S^-1/2||_1 "
                   "K^2.749, cond(S) K^0.936. Quoting any one against a value from "
                   "another ladder is the defect C17 blocked on 2026-09-08, when "
                   "0.72 was written in the headline ||M||_1 ~ K^p form.",
        aliases={1.95: "||H(lambda*)||_1 exponent, same ladder",
                 2.75: "whitened exponent, same ladder",
                 0.94: "cond(S) exponent, same ladder"}),
    "p60_cond_S_converged": dict(
        value=16.0, convention="constant: condition number of the L2 overlap S at "
                               "N=8 (K=100), helium, CONVERGED radial box",
        q=None,
        provenance="MEASURED 2026-09-07, PM-verified independently at R = 120/240/"
                   "480 (16.01 at all three). RETIRED: 3673, which is a pure "
                   "radial-box artifact -- it switches on exactly where n_max^2 "
                   "first exceeds R_MAX = 60, and every box agrees to 4 digits "
                   "below that point. Converged cond(S) grows mildly, ~0.12*K "
                   "(4.07 at K=10 to 23.49 at K=164 to 56.3 at K=452): an ordinary "
                   "Gram matrix, not an ill-conditioned one.",
        aliases={23.49: "N=10, K=164, converged", 4.07: "N=3, K=10"}),

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
    # ---- Paper 12 azimuthal channels (2026-09-14) ------------------------
    # The sigma-only restriction, not the cusp, is Paper 12's 7.6% residual.
    # Memo: debug/sprint_tmr_method_memo.md.  Backing:
    # tests/test_paper12_azimuthal_channels.py.
    "p12_sigma_only_de_pct": dict(
        value=92.42, convention="% of D_e at (j,l)=(3,3), H2 R=1.4011, sigma only, "
                                "canonical orthogonalisation",
        provenance="MEASURED 2026-09-14; Paper 12's own published 92.45 "
                   "(E=-1.161304) is 58 uHa lower and reflects cond(S)=2.6e14 "
                   "linear dependence rather than physics",
        aliases={92.45: "Paper 12's unconditioned N=72 value",
                 92.25: "(j,l)=(2,2), N=27"}),
    "p12_azimuthal_de_pct": dict(
        value=99.09, convention="% of D_e at (j,l)=(3,3), H2 R=1.4011, |m|<=1, "
                                "alpha=1.25, canonical-orthogonalisation "
                                "threshold 1e-11",
        provenance="MEASURED 2026-09-14 in Paper 12's own basis; reproduced "
                   "to 99.10 by an independent Cartesian-Gaussian full CI. "
                   "NOT a stable fourth digit: at this basis cond(S)=2e16, the "
                   "solver returns NON-VARIATIONAL values at a MAJORITY of "
                   "alpha -- six of nine grid points over [0.90, 1.30] at the "
                   "declared threshold, including alpha=1.00, the natural "
                   "default (0.95 -> -277 Ha, 1.15 -> -4.1, 1.20 -> -8.1) -- "
                   "and the value moves 99.14/99.09/98.99/98.41 across "
                   "thresholds 1e-12/1e-11/1e-10/1e-8. Over the full "
                   "(alpha, threshold) grid the weakest VARIATIONAL value is "
                   "76.11% at (0.90, 1e-8), where the sigma-only solve gives "
                   "91.06% -- i.e. the channels lose ground there; 95.5% is "
                   "the floor only for alpha >= 1.10 and must not be quoted "
                   "as a grid floor (C16: p12-grid-floor-955). Quote 99.1% at "
                   "summary surfaces; the precise value only where alpha and "
                   "threshold are stated. Envelope: 99.0-99.1%.",
        aliases={98.96: "(j,l)=(2,2)", 99.00: "(j,l)=(3,2)",
                 99.10: "independent Gaussian-basis route",
                 99.14: "same point at threshold 1e-12 (measured 99.1446)",
                 99.1: "the 3-s.f. value quoted at summary surfaces"}),
    "p12_azimuthal_gain_mha": dict(
        value=11.64, convention="mHa gained by opening |m|<=1 at (j,l)=(3,3)",
        provenance="MEASURED 2026-09-14; the sigma-axis control over the same "
                   "enlargement is 0.34 mHa"),
    "p12_sigma_growth_mha": dict(
        value=0.34, convention="mHa gained by N=27 -> 72 along the sigma axis",
        provenance="MEASURED; read from Paper 12's own Table tab:convergence "
                   "(92.2 -> 92.4 -> 92.4 % of D_e)"),
    "p12_cond_s_33": dict(
        value=2.6e14, convention="cond(S), (j,l)=(3,3) sigma only, N=72, "
                                 "alpha=1.05 -- the value at which this "
                                 "literal, the 2.0e16 for |m|<=1 and the "
                                 "-79 Ha direct-eigensolve figure all land "
                                 "together; at alpha=1.0 cond is 3.04e14. "
                                 "It varies ~3x over alpha in [0.9,1.4], "
                                 "~15x for |m|<=1",
        provenance="MEASURED 2026-09-14; rises to 2.0e16 at |m|<=1, N=144"),
    # ---- Paper 12 re-conditioning (2026-09-16) ---------------------------
    # The 99.1% monomial cap is a CONDITIONING artifact of the xi^j radial set
    # (a Hankel moment problem), not a structural ceiling.  Re-basing the SAME
    # span to an orthogonal-polynomial family is an exact change of basis and
    # the energy climbs monotonically/variationally to chemical accuracy.
    # Extended precision is required where the MONOMIAL matrices are re-based;
    # since v5.13.8 S and H1 are built directly in the orthogonal basis instead
    # (no congruence applied to them), leaving only V_ee on the mpf path.  The
    # values below are unchanged by that -- both routes agree to the float64
    # downcast floor.
    #
    # TWO DISTINCT QUANTITIES, deliberately kept as separate keys (2026-09-18,
    # v5.13.9, PI-approved).  The `p12_rebased_*` keys are the endpoint of the
    # monotone alpha = 1.0 ladder printed in tab:recondition; the
    # `p12_rebased_*_aopt` keys are the same basis at its VARIATIONAL OPTIMUM
    # alpha = 1.40.  They are not competing values and neither supersedes the
    # other: the ladder's own claim ("every point containing its predecessor and
    # lying below it") holds only at a consistent alpha, so moving one row of it
    # to the optimum would break the monotonicity it asserts.  The abstract and
    # conclusion quote the optimum; the table remains the alpha = 1.0 ladder.
    # "alpha=1.0; best variational point" in the conventions below is therefore
    # true of the DISCARD-THRESHOLD sweep, not of alpha.  See
    # debug/sprint_explicit_correlation_scoping_memo.md Sec. 4b.
    #
    # Memo: debug/sprint_h2_recondition_memo.md.  Backing:
    # tests/test_paper12_recondition.py; module geovac/prolate_recondition.py.
    "p12_rebased_de_pct": dict(
        value=99.77, convention="% of D_e, H2 R=1.4011, re-based (5,5)+delta "
                                "(|m|<=2), alpha=1.0; best variational point",
        provenance="MEASURED 2026-09-16 via geovac/prolate_recondition.py "
                   "(E=-1.1740693, 0.41 mHa, inside chemical accuracy). "
                   "Laguerre x Legendre and mu-adapted Gegenbauer give the "
                   "IDENTICAL energy (same span); Gegenbauer conditions the "
                   "downcast solve ~1000x tighter (9.1e4 vs 9.4e10).",
        aliases={99.7674: "full precision",
                 99.71: "(4,4)+delta, CI-tested anchor (E=-1.1739704, 0.50 mHa)",
                 99.216: "(3,3) pi re-based, monomial-capped at 99.09"}),
    "p12_rebased_err_mha": dict(
        value=0.41, convention="mHa above exact D_e, re-based (5,5)+delta, "
                               "alpha=1.0 (the ladder endpoint; the optimum is "
                               "p12_rebased_err_mha_aopt)",
        provenance="MEASURED 2026-09-16; monomial ceiling was 1.58 mHa (99.09)"),
    # ---- Paper 11 H2+ spectral solver.  ADDED 2026-09-19 (/qa group2 CODE run).
    # C21 had ZERO Paper-11 keys, so the numeric gate had nothing to check on
    # this paper -- which is why a headline wrong by ~7.6 orders survived.
    "p11_h2plus_err_ha": dict(
        value=3.6e-14, convention="|E_total(spectral, n_basis=20, R=2.0) - E_ref| "
                                  "in Ha, against E_ref = -0.6026342144949 Ha "
                                  "(1s-sigma-g, Bates/Ledsham/Stewart lineage). "
                                  "REFERENCE-LIMITED: the residual is at or below "
                                  "the precision at which E_ref is conventionally "
                                  "quoted, so the mantissa is NOT a stable "
                                  "quantity -- quote 'machine precision', not a "
                                  "percentage. Retired value: 0.0002% (= 1.21e-6 "
                                  "Ha), which no route reproduces; see "
                                  "p11_h2plus_de_pct_RETIRED",
        provenance="MEASURED 2026-09-19 via ProlateSpheroidalLattice("
                   "R=2.0, radial_method='spectral', n_basis=20).total_energy(): "
                   "E=-0.602634214494936. Ladder: n_basis=5 -> 2.13e-8 Ha "
                   "(3.5e-6 %), 10 -> 2.53e-12 (4.2e-10 %), 20 -> 3.64e-14 "
                   "(6.0e-12 %). Even n_basis=5 is 57x better than the retired "
                   "0.0002% headline. Cross-checked by an independent "
                   "Slater-basis/closed-form-moment/det-root solver at dps=45 "
                   "(-0.602634214494946), agreeing to 15 digits.",
        aliases={0.0: "reported qualitatively as 'machine precision'"}),
    "p11_h2plus_req_bohr": dict(
        value=1.9973, convention="fitted R_eq (bohr) from a FINE PES grid "
                                 "(0.005 spacing) + fit_spectroscopic_constants, "
                                 "vs R_ref = 1.997. The retired 2.005 / 0.38% is "
                                 "a COARSE-grid (0.05) FIT artifact, not a solver "
                                 "property: the H2+ well is flat to 1e-6 Ha over "
                                 "+-0.005 bohr, so a quadratic apex on a 0.05 grid "
                                 "lands milli-bohr off",
        provenance="MEASURED 2026-09-19: fine grid arange(1.95,2.05,0.005) -> "
                   "R_eq=1.99726092 (+0.0131% vs 1.997), E_min=-0.60263464, "
                   "D_e=0.10263464, k=0.10355766. The paper prints 2.005 / "
                   "0.38%, ~29x worse, and thereby prints its BETTER solver as "
                   "worse than its own FD-8000 control (2.001 / 0.21%).",
        aliases={1.997: "reference value", 0.0131: "percent error vs 1.997"}),
    "p12_rebased_de_pct_aopt": dict(
        value=99.81, convention="% of D_e, H2 R=1.4011, re-based (5,5)+delta "
                                "(|m|<=2), laguerre_legendre, at the VARIATIONAL "
                                "OPTIMUM of the shared radial exponent "
                                "alpha=1.40. Distinct from p12_rebased_de_pct, "
                                "which is the same basis at alpha=1.0",
        provenance="MEASURED 2026-09-18 (v5.13.9) via prolate_recondition: "
                   "E=-1.1741513, 0.324 mHa, variational, all 1944 functions "
                   "kept, cond 4.90e10 (BETTER conditioned than alpha=1.0's "
                   "9.37e10). Bracketed scan alpha=1.00/1.20/1.40/1.50 -> "
                   "0.406/0.346/0.324/0.334 mHa",
        aliases={99.8144: "full precision", 99.814: "3 s.f."}),
    "p12_rebased_err_mha_aopt": dict(
        value=0.32, convention="mHa above exact D_e, re-based (5,5)+delta at the "
                               "variational optimum alpha=1.40",
        provenance="MEASURED 2026-09-18; 0.3237 at full precision. The alpha=1.0 "
                   "value is p12_rebased_err_mha (0.41); the gain is 0.082 mHa, "
                   "20% of that residual",
        aliases={0.324: "3 s.f."}),
    "p12_rebased_alpha_opt": dict(
        value=1.40, convention="variational optimum of the single shared radial "
                               "exponent alpha at (5,5)+delta, bracketed (1.50 "
                               "is worse). Drifts UP with basis size (~1.20-1.25 "
                               "at (4,4) mu<=1 vs 1.40 at (5,5)) -- the signature "
                               "of single-exponent strain",
        provenance="MEASURED 2026-09-18 (v5.13.9)"),
    # ---- Explicit-r12 (James-Coolidge-type) extension of the algebraic engine.
    # ADDED 2026-09-20.  DISTINCT method from the re-based CI above: the basis is
    # multiplied by r12^p (p in {0,1}), and every r12 matrix element is evaluated
    # EXACTLY by the same algebraic Neumann machinery as V_ee (odd r12 powers ->
    # A K0 - B K1 via the general-m X-tables; the p0xp1 kinetic collapses by IBP
    # to -1/2 <r12 g_v Lap g_u>, no grad_r12).  No quadrature.  Memo:
    # debug/sprint_neumann_r12_build_memo.md.  Backing: tests/test_paper12_r12.py;
    # engine debug/prolate_r12_mpf.py (migrates to geovac/ when the p>=2 arc ends).
    "p12_r12_err_mha": dict(
        value=0.053, convention="mHa above exact D_e, H2 R=1.4011, explicit-r12 "
                                "prolate CI, p in {0,1}, sigma (mu=0), (j_max=3, "
                                "l_max=4) n=416, alpha=1.0, mpf-orthogonalized "
                                "solve. The re-based-CI residual for comparison is "
                                "p12_rebased_err_mha_aopt (0.32 mHa) -> 6x smaller",
        provenance="MEASURED 2026-09-20 via debug/prolate_r12_mpf.assemble_mixed "
                   "+ solve_canonical_mpf: E_tot=-1.1744215 (n=416, cond 5.0e18, "
                   "all 416 vectors kept; mpf==float64 to 0.65 uHa). n=400 (4,3) "
                   "gives 0.0528 (cond 2.7e19, mpf keeps 399/400, 16 uHa below "
                   "float64). Matches the Tao-McCurdy-Rescigno grid-prolate 0.05 "
                   "mHa in the same coordinates, with EXACT algebraic integrals. "
                   "Verified: monotone from ABOVE over (2,0)->(2,2)->(3,2), "
                   "-20.66 -> -0.167 -> -0.079 mHa; block engine validated vs "
                   "quadrature to 1.4e-4. PLATEAUS here: p<=1 single-alpha space "
                   "near-spanned (cond 1e19); alpha not the lever (scan optimum "
                   "~1.0-1.2, a=1.40 HURTS); microhartree needs r12 powers p>=2",
        aliases={0.0528: "best point, (4,3) n=400 mpf",
                 0.0535: "(3,4) n=416, all vectors kept"}),
    "p12_r12_heh_err_mha": dict(
        value=0.88, convention="mHa above reference for HeH+ (2e heteronuclear, "
                               "R=1.4632), explicit-r12 engine, (j=2,l=3) n=288, "
                               "a=1.6, p in {0,1}. Probe result (B arc): the "
                               "exact-algebraic r12 GENERALIZES to a 2nd center",
        provenance="MEASURED 2026-09-20 via debug/heh_converge.py "
                   "(E_tot=-2.977808 vs E_ref=-2.97869). Converges -133 (2,1) -> "
                   "-11.7 (2,2) -> -0.88 (2,3); lever is ANGULAR basis (l_max+1 = "
                   "-10.8 mHa) not exponent (two-block +0.3) or radial (+0.25). "
                   "New heteronuclear V_ne validated vs quad 6.1e-5. Reference "
                   "-2.97869 is LOAD-BEARING (HeH+ X1Sigma+ BO near R_e; verify)",
        aliases={0.882: "full precision", 0.9: "1 s.f. (paper display)"}),
    "p12_r12_de_pct": dict(
        value=99.97, convention="% of D_e, H2 R=1.4011, explicit-r12 prolate CI "
                                "p in {0,1}, (3,4) n=416; = 1 - p12_r12_err_mha/"
                                "(1000*0.174475). Beats the re-based CI's "
                                "p12_rebased_de_pct_aopt (99.81%)",
        provenance="MEASURED 2026-09-20; 99.9693 at full precision (n=416), "
                   "99.9697 at n=400. p=0 control at the SAME (2,2) basis is "
                   "92.25% -> r12 cuts the error ~80x at matched radial/angular",
        aliases={99.9693: "full precision n=416", 99.9697: "n=400"}),
    "p12_recond_cond_gain": dict(
        value=326, convention="normalized cond(S) ratio Laguerre/Gegenbauer at "
                              "(3,3), |m|<=1 (3.49e5 / 1.07e3 = 326, correct); "
                              "grows to ~1.03e6x at (5,5)+delta (9.37e10 / "
                              "9.14e4), i.e. SIX orders of magnitude",
        provenance="MEASURED 2026-09-16; same span, same energy to <1 uHa. "
                   "CORRECTED 2026-09-18 (v5.13.9): this convention read "
                   "'~1025x at (5,5)+delta' -- that is 1.025e6 with the e3 "
                   "dropped, a 1000x UNDERSTATEMENT of the corpus's own result, "
                   "and Paper 12 inherited it as 'two to three orders of "
                   "magnitude' (true of the 326x (3,3) case, wrong for the "
                   "(5,5)+delta figures it was paired with). Recomputed: "
                   "9.37e10/9.14e4 = 1.025e6. At the new alpha=1.40 optimum the "
                   "ratio is 4.90e10/3.72e4 = 1.32e6 (6.1 orders)"),

    "lih_composed_pauli": dict(
        value=837, convention="non-identity Pauli, LiH composed Q=30",
        provenance="MEASURED 2026-08-30",
        aliases={838: "including identity"}),
    "h2o_composed_pauli": dict(
        value=1953, convention="non-identity Pauli, H2O composed Q=70",
        provenance="MEASURED 2026-08-30",
        aliases={1954: "including identity"}),
    # Added 2026-09-19 (v5.14.6).  The LiH and H2O composed counts were
    # registered; BeH2 -- the middle row of the same block-topology table --
    # was not, so the retired 556 had no canonical value to be checked
    # against and SCOPE_BOUNDARY.md carried it unguarded.
    "beh2_composed_pauli": dict(
        value=1395, convention="non-identity Pauli, BeH2 composed Q=50",
        provenance="CITED, test-backed: tests/test_general_builder.py:112 "
                   "asserts N_pauli in (1395, 1396) and "
                   "tests/test_spin_ful_composed.py:35 pins 1396 for "
                   "beh2_spec; consistent with 27.90 x 50 = 1395 and with the "
                   "exact-rule regression row in docs/validation_benchmarks.md "
                   "(838 / 1,396 / 1,954). Retired pair-diagonal value: 556.",
        aliases={1396: "including identity"}),
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
    "p60_exc_ratio_k452":   ("p60_exc_gap_k452 / chem_accuracy_mha", 0.005,
                             "He 2^1S error at K=452, in units of chemical accuracy"),
    "p60_gnd_ratio_k452":   ("p60_gnd_gap_k452 / chem_accuracy_mha", 0.005,
                             "He ground-state error at K=452, same units, same K, same ||M||_1"),

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
    # MERGED 2026-09-06.  This key was written TWICE -- as `1.69` and as
    # `1.690` -- which are the same float, so the dict literal silently
    # collapsed them and the second overwrote the first, discarding its
    # `alpha|lambda` detection anchors.  A guard disarmed by a dict literal
    # is exactly the "guard that cannot fail" class, so `test_numeric_registry`
    # now scans the source for duplicate literal keys.
    #
    # `forbid` gains spheroidal|DLMF: Paper 58's DLMF section reports a
    # log-gap local slope of -1.69, and the anchor `exponent` matches inside
    # the word "exponential" in that prose.  A 1-norm scaling exponent and a
    # spheroidal eigenvalue splitting are unrelated quantities.
    1.69:    ("exp_lambda_4pt",     r"Q\^|exponent|scaling|alpha|lambda",
                                    r"meV|polarizability|spheroidal|DLMF"),
    3.36:    ("exp_qwc_3pt",        r"Q\^|exponent|scaling|QWC",   None),
    3.147:   ("exp_pauli_4pt",      r"exponent|alpha|scaling",  None),
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
    # 11.11 added 2026-09-19 (v5.14.6): docs/validation_benchmarks.md carried
    # the ratio as "11.11 +- 0.1", not 11.10, so an exact-value check on the
    # registered retired literal could not have caught it even in scope.
    11.11:   ("composed_coeff",     r"coefficient|per qubit|\\times Q|Pauli/", None),
    9.23:    ("composed_coeff_d",   r"coefficient|per qubit|Pauli/", None),
    # 556 added 2026-09-19 (v5.14.6): the BeH2 block-topology count, retired
    # with 334/778 but never registered alongside them.
    556:     ("beh2_composed_pauli", r"Pauli",                  None),
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
    # Accuracy-side quantities (Paper 12 azimuthal channels, 2026-09-14):
    # a percentage of D_e, an energy gain in mHa, a condition number.  No
    # identity-in/out convention applies -- these count no Pauli terms -- so
    # family() leaves identity None, as for `density` and `constant`.  What
    # the convention string must still carry is the BASIS and truncation the
    # number was measured at, because 92.25 and 92.42 are the same quantity
    # at different (j_max, l_max) and pairing them would be the twin defect
    # this registry exists to catch.
    ("accuracy", ("% of d_e", "mha gained", "cond(s)")),
    # Absolute energy RESIDUALS (Paper 11 H2+, Paper 12 re-based H2;
    # registered 2026-09-19).  |E - E_ref| in Ha, or mHa above exact.  No
    # identity-in/out convention applies.  What the convention string must
    # carry instead is (a) the REFERENCE the residual is measured against and
    # (b) the basis/truncation -- because a residual smaller than the
    # reference's own quoted precision is not a measurement of the method.
    # That is precisely the Paper-11 defect this kind was added for: a
    # 3.6e-14 Ha residual reported as "0.0002%" against a reference given to
    # 13 digits.  Added rather than reworded: the three keys it covers
    # (p11_h2plus_err_ha, p12_rebased_err_mha, p12_rebased_err_mha_aopt) are
    # genuinely residuals, and relabelling them "% of d_e" to satisfy the
    # parser would register a false convention to quiet a gate.
    ("residual", ("in ha, against", "mha above exact", "residual")),
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
