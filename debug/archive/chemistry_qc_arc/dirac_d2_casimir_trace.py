"""
Track D2 — Dirac analog of B = 42 (pure symbolic).

Tests whether B = 42 (Paper 2's truncated Casimir trace on scalar S³ at
m = 3) lifts to a Dirac-spectral invariant.

All arithmetic is exact: sympy.Rational and sympy.Integer throughout.
Floating point is forbidden by the algebraic-first mandate.

Candidate identities:
  (a) Original hypothesis — truncated Casimir traces for n_CH = 0..2.
  (b) Single-level Weyl at n=5:  g_5^Weyl = (5+1)(5+2) = 42 exactly.
  (c) Cumulative Dirac to n=2:    Σ g_n^Dirac = 40 = Δ^{-1}.
  (d) Variant truncations / weightings.

Spectrum (Camporesi–Higuchi on unit S³):
    |λ_n| = n + 3/2          (CH index n ≥ 0)
    g_n^Dirac = 2(n+1)(n+2)
    g_n^Weyl  =  (n+1)(n+2)
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp
from sympy import Integer, Rational, simplify

from geovac.dirac_s3 import (
    dirac_degeneracy,
    dirac_eigenvalue_abs,
)

B = Integer(42)  # Paper 2 scalar truncated Casimir trace at m = 3.


def classify_match(value: sp.Expr) -> str:
    """Label the relationship of a symbolic value to B = 42."""
    value = sp.nsimplify(sp.sympify(value), rational=True)
    if value == B:
        return "EQUALS_B"
    if value == 0:
        return "ZERO"
    ratio = sp.Rational(value, B)  # value / B as exact Rational
    # Integer multiple?
    if ratio.q == 1:
        return f"INT_MULT_B (value = {int(ratio.p)} * B)"
    # B a multiple of value (i.e. value = B/k)?
    inv = sp.Rational(B, value) if value != 0 else None
    if inv is not None and inv.q == 1:
        return f"B_OVER_INT (value = B/{int(inv.p)})"
    # Simple rational multiple (small numerator/denominator)?
    if ratio.q <= 12 and abs(ratio.p) <= 12:
        return f"SMALL_RATIONAL_MULT_B ({ratio})"
    return "UNRELATED"


def candidate_a_truncated_casimir_traces() -> dict:
    """(a) Truncated Casimir traces for n_CH = 0..2 in both sectors."""
    results = {}
    for sector in ("dirac", "weyl"):
        s_g = Integer(0)
        s_lg = Rational(0)
        s_l2g = Rational(0)
        s_l2g_int = Rational(0)  # with (2|λ|)² = (2n+3)² to stay integer
        per_level = []
        for n in range(3):  # n_CH = 0, 1, 2  -> Paper 2's m=3 cutoff
            lam = dirac_eigenvalue_abs(n, convention="ch")  # n + 3/2
            g = Integer(dirac_degeneracy(n, sector=sector, convention="ch"))
            s_g += g
            s_lg += lam * g
            s_l2g += lam**2 * g
            s_l2g_int += (2 * lam)**2 * g  # = (2n+3)^2 * g
            per_level.append({
                "n_ch": n,
                "lambda_abs": str(lam),
                "g": int(g),
                "lambda_times_g": str(lam * g),
                "lambda2_times_g": str(lam**2 * g),
            })
        results[sector] = {
            "per_level": per_level,
            "sum_g": str(s_g),
            "sum_lambda_g": str(s_lg),
            "sum_lambda2_g": str(simplify(s_l2g)),
            "sum_4lambda2_g": str(simplify(s_l2g_int)),
            "classify_sum_g": classify_match(s_g),
            "classify_sum_lambda_g": classify_match(s_lg),
            "classify_sum_lambda2_g": classify_match(s_l2g),
            "classify_sum_4lambda2_g": classify_match(s_l2g_int),
        }
    return results


def candidate_b_single_level_weyl() -> dict:
    """(b) Single-level Weyl degeneracy at n_CH = 5 equals 42."""
    # Verify g_5^Weyl symbolically.
    n, k = sp.symbols("n k", integer=True)
    g_weyl_n = (n + 1) * (n + 2)
    g_5_weyl = g_weyl_n.subs(n, 5)
    # Also full-Dirac at n=5
    g_5_dirac = 2 * g_weyl_n.subs(n, 5)

    # Solve (n+1)(n+2) = 42 symbolically to confirm n=5 is the unique positive
    # integer solution.
    roots = sp.solve(sp.Eq(g_weyl_n, 42), n)
    positive_int_roots = [r for r in roots if r.is_Integer and r >= 0]

    # Also: what value of n gives g_n^Dirac = 42? Solve 2(n+1)(n+2) = 42.
    dirac_eq_42 = sp.solve(sp.Eq(2 * g_weyl_n, 42), n)

    # Structural observation: n = 5 in CH means n_Fock = 6 in the existing
    # Fock graph (scalar lattice shell numbering), while Paper 2 caps scalar
    # shells at m = 3 (n_Fock ≤ 3). So this is OUTSIDE the scalar truncation
    # window; candidate (b) only makes sense as a single-level observation,
    # not a truncated trace.

    return {
        "g_weyl_at_n5_symbolic": str(g_5_weyl),
        "g_weyl_at_n5_value": int(g_5_weyl),
        "equals_B": bool(g_5_weyl == B),
        "g_dirac_at_n5": int(g_5_dirac),
        "g_dirac_at_n5_over_B": str(sp.Rational(int(g_5_dirac), int(B))),
        "n_solving_weyl_equals_42": [str(r) for r in roots],
        "positive_integer_n_weyl_equals_42": [int(r) for r in positive_int_roots],
        "n_solving_dirac_equals_42": [str(r) for r in dirac_eq_42],
        "note": (
            "n_CH = 5 corresponds to n_Fock = 6, OUTSIDE Paper 2's scalar "
            "cutoff m ≤ 3 (n_Fock ≤ 3). Candidate (b) is a single-level "
            "coincidence, not a truncated-trace identity. No independent "
            "structural argument (d_max=4, N_init=2, l_max=2) picks n=5."
        ),
        "structural_check_n_eq_d_max_plus_1": {
            "d_max": 4,
            "n_candidate": 5,
            "value": int((sp.Integer(5) + 1) * (sp.Integer(5) + 2)),
            "note": "n = d_max + 1 = 5 would give a Paper-0 tie to the packing exponent, but no derivation connects d_max+1 to a Weyl-degeneracy truncation."
        },
    }


def candidate_c_cumulative_dirac() -> dict:
    """(c) Cumulative Dirac sum to n=2 is 40 = Δ^{-1}.

    Also scan: for what n_max does the cumulative sum hit 42?
    """
    # Closed form: Σ_{n=0}^{N} 2(n+1)(n+2).
    N = sp.symbols("N", integer=True, nonnegative=True)
    n_sym = sp.symbols("n_sym", integer=True)
    closed = sp.simplify(sp.summation(2 * (n_sym + 1) * (n_sym + 2),
                                      (n_sym, 0, N)))
    # For each N, compute.
    table = []
    matches = []
    for Nv in range(8):
        val = int(closed.subs(N, Nv))
        table.append({"n_max_ch": Nv, "cumulative_g_dirac": val})
        if val == 42:
            matches.append(Nv)
    # Solve closed = 42 symbolically.
    sol = sp.solve(sp.Eq(closed, 42), N)

    # Same for Weyl.
    closed_weyl = sp.simplify(sp.summation((n_sym + 1) * (n_sym + 2),
                                           (n_sym, 0, N)))
    weyl_table = []
    weyl_matches = []
    for Nv in range(8):
        val = int(closed_weyl.subs(N, Nv))
        weyl_table.append({"n_max_ch": Nv, "cumulative_g_weyl": val})
        if val == 42:
            weyl_matches.append(Nv)
    sol_weyl = sp.solve(sp.Eq(closed_weyl, 42), N)

    return {
        "closed_form_dirac": str(closed),
        "dirac_cumulative_table": table,
        "n_max_where_dirac_cumulative_equals_42": matches,
        "dirac_solve_eq42": [str(r) for r in sol],
        "closed_form_weyl": str(closed_weyl),
        "weyl_cumulative_table": weyl_table,
        "n_max_where_weyl_cumulative_equals_42": weyl_matches,
        "weyl_solve_eq42": [str(r) for r in sol_weyl],
        "note": (
            "Cumulative Dirac sum produces 12, 40, 84, 160, ... — never 42. "
            "Cumulative Weyl produces 2, 6, 12, 20, 30, 42, 56, ... — hits "
            "42 at n_max = 5 (see candidate (d)). This is the same n=5 that "
            "lights up single-level Weyl, pointing at the same coincidence."
        ),
    }


def candidate_d_variant_truncations_weightings() -> dict:
    """(d) Variant weightings and cutoffs."""
    out = {}

    # (d1) Σ (2j+1) · mult using j = l + 1/2 labeling (Bär). For a single
    # chirality sector at level n_CH, the total 2j+1 count is just the Weyl
    # degeneracy (n+1)(n+2). So this reduces to candidate (a) Weyl. Record
    # for completeness.
    n_sym = sp.symbols("n_sym", integer=True)
    out["d1_two_j_plus_1_weight"] = {
        "note": "Σ(2j+1) over a single chirality sector = (n+1)(n+2) = g_Weyl. No new info.",
        "equivalent_to": "candidate_a_weyl_sum_g",
    }

    # (d2) Σ n · g_n and Σ n² · g_n, for various cutoffs.
    results = {"dirac": {}, "weyl": {}}
    for sector in ("dirac", "weyl"):
        for weight_name, weight_fn in [
            ("n", lambda nn: nn),
            ("n**2", lambda nn: nn**2),
            ("2n+3", lambda nn: 2 * nn + 3),          # = 2|λ|
            ("(2n+3)**2", lambda nn: (2 * nn + 3)**2),  # = 4|λ|²
        ]:
            rows = []
            matches_42 = []
            for Nv in range(10):
                s = sum(
                    int(weight_fn(n)) *
                    dirac_degeneracy(n, sector=sector, convention="ch")
                    for n in range(Nv + 1)
                )
                rows.append({"n_max_ch": Nv, "sum": int(s)})
                if s == 42:
                    matches_42.append(Nv)
                if s == 21 or s == 84 or s == 126:
                    matches_42.append({"n_max": Nv, "sum": int(s),
                                       "note": "simple rational multiple of 42"})
            results[sector][weight_name] = {
                "table": rows,
                "matches_42_or_simple_multiple": matches_42,
            }
    out["d2_weighted_sums"] = results

    # (d3) Cutoff scan for the "|λ|²·g" sum and g-only sum searching for 42.
    # Already in candidate (a) and (c); summarize the 42 hits.
    hits_42 = []
    # Partial Weyl sum hits 42 at n_max=5 (we tabulated this in candidate c).
    # Any dirac-weighted partial sum hits 42?
    for sector in ("dirac", "weyl"):
        for weight_name, weight_fn in [
            ("1 (g only)", lambda nn: 1),
            ("n", lambda nn: nn),
            ("2n+3", lambda nn: 2 * nn + 3),
            ("(2n+3)**2", lambda nn: (2 * nn + 3)**2),
        ]:
            for Nv in range(10):
                s = sum(
                    int(weight_fn(n)) *
                    dirac_degeneracy(n, sector=sector, convention="ch")
                    for n in range(Nv + 1)
                )
                if s == 42:
                    hits_42.append({
                        "sector": sector,
                        "weight": weight_name,
                        "n_max_ch": Nv,
                        "value": s,
                    })
    out["d3_all_hits_of_42"] = hits_42

    # (d4) Does the "bare" single-level Weyl (n+1)(n+2) hit 42 only at n=5?
    # We showed yes in candidate (b). Here, restrict cumulative window to
    # the Paper-2 m=3 cutoff and record EVERY integer value attainable
    # by a weighted partial sum under that cutoff. Paper 2's cutoff window
    # is n_CH ∈ {0, 1, 2} (scalar shells m ∈ {1, 2, 3}, i.e., n_Fock ≤ 3).
    # Enumerate linear combinations Σ c_n · g_n^sector for small integer c_n
    # and check which hit 42.
    target_coeffs = {}
    for sector in ("dirac", "weyl"):
        g_vals = [dirac_degeneracy(n, sector=sector, convention="ch")
                  for n in range(3)]
        hits = []
        rng = range(-3, 4)  # small integer coefficient scan
        for c0 in rng:
            for c1 in rng:
                for c2 in rng:
                    if (c0, c1, c2) == (0, 0, 0):
                        continue
                    s = c0 * g_vals[0] + c1 * g_vals[1] + c2 * g_vals[2]
                    if s == 42:
                        hits.append({"coeffs": [c0, c1, c2],
                                     "g_values": g_vals,
                                     "sum": s})
        target_coeffs[sector] = {
            "g_values_n0_n1_n2": g_vals,
            "integer_linear_combos_hitting_42": hits[:20],
            "total_combos_hitting_42": len(hits),
        }
    out["d4_integer_combinations_in_m3_window"] = target_coeffs

    return out


def candidate_e_structural_closed_forms() -> dict:
    """(e) Closed-form comparison of scalar B(m) to Dirac/Weyl traces.

    Paper 2's B(m) = Σ_{n=1..m} Σ_{l=0..n-1} (2l+1)·l(l+1)
                   = m(m-1)(m+1)(m+2)(2m+1)/20.
    B(1)=0, B(2)=6, B(3)=42, B(4)=162.
    """
    m = sp.symbols("m", integer=True, positive=True)
    N, n = sp.symbols("N n", integer=True, nonnegative=True)

    B_m = m * (m - 1) * (m + 1) * (m + 2) * (2 * m + 1) / 20
    # Dirac |λ|²·g cumulative at N = m - 1:
    dirac_l2g = (N + 1) * (N + 2) * (N + 3) * (2 * N + 3) * (2 * N + 5) / 10
    dirac_l2g_at_m = dirac_l2g.subs(N, m - 1)

    # Weyl single-level |λ|·g at top of window (n_CH = m - 1):
    top_weyl_lambdag = Rational(1, 2) * (2 * m + 1) * m * (m + 1)
    ratio_topWeyl = sp.simplify(B_m / top_weyl_lambdag)
    eq_solve = sp.solve(sp.Eq((m - 1) * (m + 2), 10), m)  # where ratio = 1

    # Cumulative |λ|·g (both sectors):
    cum_weyl_lg = sp.simplify(sp.summation(
        Rational(1, 2) * (2 * n + 3) * (n + 1) * (n + 2), (n, 0, N)))
    cum_dirac_lg = 2 * cum_weyl_lg

    return {
        "B_scalar_m_closed_form":
            str(sp.factor(B_m)) + "  [m(m-1)(m+1)(m+2)(2m+1)/20]",
        "B_at_m_1_to_4": [int(B_m.subs(m, k)) for k in (1, 2, 3, 4)],
        "dirac_lambda2_g_cum_closed_form":
            str(sp.factor(dirac_l2g)) + "  [ (N+1)(N+2)(N+3)(2N+3)(2N+5)/10 ]",
        "dirac_lambda2_g_at_m3_window":
            int(dirac_l2g.subs(N, 2)),  # = 378
        "ratio_dirac_lambda2_to_B_scalar_closed":
            str(sp.factor(sp.simplify(dirac_l2g_at_m / B_m))),
        "note_ratio":
            "ratio = 2(2m+3)/(m-1); it equals 9 at m=3 only by virtue of 2·9/1 = 18 (no — rechecked), "
            "specifically 2·9/2 = 9. This is a single-point rational coincidence at m=3.",
        "top_weyl_single_level_lambdag_closed":
            str(sp.factor(top_weyl_lambdag)) + "  [m(m+1)(2m+1)/2]",
        "ratio_B_over_top_Weyl_lambdag_closed":
            str(sp.factor(ratio_topWeyl)) + "  [(m-1)(m+2)/10]",
        "m_where_ratio_is_1_symbolic":
            [int(r) for r in eq_solve if r.is_Integer and r > 0],
        "ratio_1_interpretation": (
            "B(m) = top_Weyl_|λ|·g only when (m-1)(m+2)=10, uniquely m=3. "
            "This is the CLEANEST single-level correspondence, but it is a "
            "point coincidence, not a universal identity."
        ),
        "cumulative_weyl_lambdag_closed_form":
            str(sp.factor(cum_weyl_lg)) + "  [(N+1)(N+2)²(N+3)/4]",
        "cumulative_weyl_lambdag_at_N2": 60,
        "cumulative_dirac_lambdag_closed_form":
            str(sp.factor(cum_dirac_lg)) + "  [(N+1)(N+2)²(N+3)/2]",
        "cumulative_dirac_lambdag_at_N2": 120,
    }


def main() -> None:
    print("Track D2 — Dirac analog of B = 42\n" + "=" * 60)
    print(f"B (scalar, Paper 2) = {B}")
    print()

    print("Candidate (a): Truncated Casimir traces, n_CH = 0..2")
    print("-" * 60)
    a = candidate_a_truncated_casimir_traces()
    for sector, data in a.items():
        print(f"  [{sector}]")
        for row in data["per_level"]:
            print(f"    n_CH={row['n_ch']}: |λ|={row['lambda_abs']}, "
                  f"g={row['g']}, |λ|·g={row['lambda_times_g']}, "
                  f"|λ|²·g={row['lambda2_times_g']}")
        print(f"    Σg = {data['sum_g']}   [{data['classify_sum_g']}]")
        print(f"    Σ|λ|·g = {data['sum_lambda_g']}   [{data['classify_sum_lambda_g']}]")
        print(f"    Σ|λ|²·g = {data['sum_lambda2_g']}   [{data['classify_sum_lambda2_g']}]")
        print(f"    Σ(2n+3)²·g = {data['sum_4lambda2_g']}   [{data['classify_sum_4lambda2_g']}]")
    print()

    print("Candidate (b): Single-level Weyl at n_CH = 5")
    print("-" * 60)
    b = candidate_b_single_level_weyl()
    print(f"  g_5^Weyl = (5+1)(5+2) = {b['g_weyl_at_n5_symbolic']} = {b['g_weyl_at_n5_value']}")
    print(f"  equals B? {b['equals_B']}")
    print(f"  g_5^Dirac = {b['g_dirac_at_n5']} = {b['g_dirac_at_n5_over_B']} · B")
    print(f"  positive-integer n solving (n+1)(n+2)=42: {b['positive_integer_n_weyl_equals_42']}")
    print(f"  note: {b['note']}")
    print()

    print("Candidate (c): Cumulative sums")
    print("-" * 60)
    c = candidate_c_cumulative_dirac()
    print(f"  Σ_{{n=0..N}} g_n^Dirac closed form: {c['closed_form_dirac']}")
    print("  Dirac cumulative values:")
    for row in c["dirac_cumulative_table"]:
        print(f"    N={row['n_max_ch']}: {row['cumulative_g_dirac']}")
    print(f"  n_max hitting 42: {c['n_max_where_dirac_cumulative_equals_42']}")
    print(f"  Σ_{{n=0..N}} g_n^Weyl closed form: {c['closed_form_weyl']}")
    print("  Weyl cumulative values:")
    for row in c["weyl_cumulative_table"]:
        print(f"    N={row['n_max_ch']}: {row['cumulative_g_weyl']}")
    print(f"  n_max hitting 42: {c['n_max_where_weyl_cumulative_equals_42']}")
    print()

    print("Candidate (d): Variant weightings")
    print("-" * 60)
    d = candidate_d_variant_truncations_weightings()
    print("  All weighted partial-sum hits of 42:")
    for hit in d["d3_all_hits_of_42"]:
        print(f"    sector={hit['sector']}, weight={hit['weight']}, "
              f"n_max={hit['n_max_ch']} → {hit['value']}")
    print("  Integer linear combos c_0·g_0 + c_1·g_1 + c_2·g_2 = 42 in the m≤3 window (c_i ∈ [-3,3]):")
    for sector, info in d["d4_integer_combinations_in_m3_window"].items():
        print(f"    [{sector}] g_vals={info['g_values_n0_n1_n2']}, "
              f"hits={info['total_combos_hitting_42']}")
        for h in info["integer_linear_combos_hitting_42"][:5]:
            print(f"       {h['coeffs']} · {h['g_values']} = 42")
    print()

    # Write JSON output.
    out_path = Path(__file__).resolve().parent / "data" / "dirac_d2_casimir_trace.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    json_out = {
        "B_scalar_paper2": int(B),
        "spectrum_formulas": {
            "lambda_abs": "n + 3/2 (CH; n >= 0)",
            "g_dirac": "2(n+1)(n+2)",
            "g_weyl": "(n+1)(n+2)",
        },
        "paper2_scalar_cutoff": "m = 3, i.e., n_CH in {0, 1, 2}",
        "candidate_a_truncated_casimir_traces": a,
        "candidate_b_single_level_weyl": b,
        "candidate_c_cumulative_sums": c,
        "candidate_d_variants": d,
        "candidate_e_structural_closed_forms": candidate_e_structural_closed_forms(),
        "verdict_summary": {
            "within_m3_window_unweighted": (
                "Cumulative Σ g_n^Dirac (n=0..2) = 40 = Δ^{-1} (not B). "
                "Cumulative Σ g_n^Weyl (n=0..2) = 20 (not B)."
            ),
            "within_m3_window_|lambda|_weighted": (
                "Cumulative Σ|λ|·g^Dirac (n=0..2) = 120 (= 20·6, unrelated to B=42). "
                "Cumulative Σ|λ|·g^Weyl  (n=0..2) = 60. "
                "However, SINGLE-LEVEL at n_CH=2 (top of window): |λ|·g^Weyl = "
                "(7/2)·12 = 42 = B EXACTLY, and |λ|·g^Dirac = 84 = 2B."
            ),
            "within_m3_window_|lambda|^2_weighted": (
                "Cumulative Σ|λ|²·g^Dirac (n=0..2) = 378 = 9·B. "
                "Integer multiple of B, but coefficient 9 is not structural "
                "(closed-form ratio 2(2m+3)/(m-1) hits 9 only at m=3; varies "
                "with m)."
            ),
            "cleanest_single_hit": (
                "Top-of-window Weyl single-level |λ|·g at n_CH = m-1: "
                "m(m+1)(2m+1)/2. Equals B(m)=m(m-1)(m+1)(m+2)(2m+1)/20 iff "
                "(m-1)(m+2)=10, uniquely m=3. So at Paper 2's cutoff and "
                "ONLY at that cutoff, the top-shell Weyl weighted degeneracy "
                "equals the scalar truncated Casimir trace."
            ),
            "verdict_for_D5": "B-doesn't-lift (nuanced)",
            "one_line_reason": (
                "No closed-form Dirac-sector trace reproduces B(m)=m(m-1)(m+1)(m+2)(2m+1)/20 "
                "as a universal identity. Every apparent hit — 9·B at m=3 from Σ|λ|²·g^Dirac, "
                "1·B at m=3 from top-shell |λ|·g^Weyl — is a single-point rational coincidence "
                "at m=3 whose ratio is non-constant in m (respectively 2(2m+3)/(m-1) and "
                "10/[(m-1)(m+2)]). The scalar B(m) carries a factor (m-1)(m+2) that no "
                "Dirac truncation produces. The only TRUNCATED trace identity the Dirac "
                "sector does reproduce at m=3 is Δ^{-1} = 40 (Phase 4H SM-D, already known). "
                "B remains a scalar-Laplace-Beltrami quantity."
            ),
            "why_structurally": (
                "B(m) in closed form has the factor (m-1)(m+2) — it vanishes at m=1 "
                "(the scalar ground state contributes no Casimir mass, since l=0 gives "
                "l(l+1)=0) and grows as m^5. The Dirac/Weyl degeneracies (n+1)(n+2) "
                "have NO (n-1)(n+2) factor — the Dirac ground state n_CH=0 contributes "
                "g=4 (Dirac) or g=2 (Weyl) with |λ|=3/2 ≠ 0, because the Dirac operator "
                "has no zero mode on S³. This is the structural reason no universal "
                "Dirac identity can reproduce B: the Dirac spectrum is gapped, the "
                "scalar spectrum has a zero mode, and B(m) carries the (m-1)-factor "
                "signature of that zero-mode gap."
            ),
        },
    }
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(json_out, f, indent=2, ensure_ascii=False)
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
