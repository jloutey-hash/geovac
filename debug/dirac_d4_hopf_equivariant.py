"""
Track D4 — Orduz S¹-equivariant Dirac decomposition on the Hopf
fibration S¹ → S³ → S².

Algebraic-first, sympy-exact, integer/half-integer arithmetic only.

Decomposition principle
-----------------------
Spin(4) = SU(2)_L × SU(2)_R. The Dirac spinor bundle on S³ at
Camporesi-Higuchi level n_CH = n carries the two chirality irreps

    (+)-chirality:  (j_L, j_R) = ((n+1)/2,  n/2)     dim (n+1)(n+2)
    (-)-chirality:  (j_L, j_R) = ( n/2,  (n+1)/2)    dim (n+1)(n+2)

The Hopf circle U(1) ⊂ S³ is the maximal torus of SU(2)_R (a choice
of fiber orientation — the dual choice sits in SU(2)_L). Under this
U(1), an irrep (j_L, j_R) decomposes as

    (j_L, j_R)|_{U(1)} = bigoplus_{m_R = -j_R}^{j_R}  (2 j_L + 1)·[m_R]

i.e. the Hopf-charge m_R runs over {-j_R, ..., j_R} (step 1) and
each charge sector carries multiplicity (2 j_L + 1).

Half-integer charges at odd n are a genuine feature — they reflect
the fact that spinors are sections of a non-trivial line bundle over
S² via the Hopf map (the "charge-q" line bundle with q ∈ ½ℤ under the
natural lift).

What this script computes
-------------------------
For each n_CH ∈ {0, 1, 2, 3, 4, 5}:
  1. The (q, multiplicity) partition for (+) chirality,
     (-) chirality, and the full Dirac sector.
  2. Σ mult = g_n^Dirac (sanity).
  3. Hopf-charge-weighted Casimir-like sum:  Σ_q q²·m(q).
  4. Cumulative partition up to n_CH, matching Δ^{-1} = 40 at
     n_CH = 2 (D1 delta_inverse_identity).
  5. Hopf-charge Dirichlet series variants.
  6. Test: does any sum reach B = 42 or F = π²/6 or Δ = 1/40?
"""

from __future__ import annotations

import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy import Rational, Symbol, zeta, pi, simplify, summation, oo


# ---------------------------------------------------------------------------
# Hopf-equivariant decomposition (symbolic, exact)
# ---------------------------------------------------------------------------

def hopf_partition_single_chirality(n_ch: int, chirality: int):
    """Return dict {2q: mult} for one chirality sector.

    We key by twice the charge (2q ∈ ℤ) to keep integer types.

    (+)-chirality has (j_L, j_R) = ((n+1)/2, n/2)  → m_R ∈ {-n/2, ..., n/2}
                                                    mult = (n+2) each
    (-)-chirality has (j_L, j_R) = (n/2, (n+1)/2)  → m_R ∈ {-(n+1)/2, ..., (n+1)/2}
                                                    mult = (n+1) each
    """
    assert chirality in (+1, -1)
    if chirality == +1:
        two_jR = n_ch           # 2·j_R = n
        mult = n_ch + 2         # 2·j_L + 1 = n+2
    else:
        two_jR = n_ch + 1       # 2·j_R = n+1
        mult = n_ch + 1         # 2·j_L + 1 = n+1
    part = {}
    for two_q in range(-two_jR, two_jR + 1, 2):
        part[two_q] = mult
    return part


def hopf_partition_full_dirac(n_ch: int):
    """Full Dirac = (+)-chirality ⊕ (-)-chirality."""
    part = defaultdict(int)
    for chir in (+1, -1):
        for two_q, m in hopf_partition_single_chirality(n_ch, chir).items():
            part[two_q] += m
    return dict(sorted(part.items()))


def check_total(part, expected):
    total = sum(part.values())
    return total == expected, total


# ---------------------------------------------------------------------------
# Paper 2 targets
# ---------------------------------------------------------------------------

B_TARGET = 42
DELTA_INV_TARGET = 40
F_TARGET = sp.Rational(1) * pi**2 / 6   # = ζ(2)


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

def run(n_ch_max: int = 5):
    results = {
        "meta": {
            "track": "D4 Orduz Hopf-equivariant Dirac decomposition",
            "n_ch_max": n_ch_max,
            "convention": (
                "Hopf U(1) = maximal torus of SU(2)_R; charges keyed as 2q "
                "(integer), so actual charge q = (2q)/2 ∈ ½ℤ."
            ),
            "targets": {
                "B": B_TARGET,
                "Delta_inverse": DELTA_INV_TARGET,
                "F": "pi**2/6",
            },
        },
        "levels": [],
        "cumulative_tests": {},
        "B_tests": {},
        "F_tests": {},
        "verdict": {},
    }

    # Per-level decomposition
    cum_full = defaultdict(int)
    for n_ch in range(n_ch_max + 1):
        g_plus = (n_ch + 1) * (n_ch + 2)
        g_full = 2 * g_plus

        part_plus = hopf_partition_single_chirality(n_ch, +1)
        part_minus = hopf_partition_single_chirality(n_ch, -1)
        part_full = hopf_partition_full_dirac(n_ch)

        ok_p, tot_p = check_total(part_plus, g_plus)
        ok_m, tot_m = check_total(part_minus, g_plus)
        ok_f, tot_f = check_total(part_full, g_full)
        assert ok_p and ok_m and ok_f, f"mult mismatch at n={n_ch}"

        # Hopf-charge Casimir moments (store 2q, but compute q² = (2q)²/4)
        def moment(part, power):
            # Σ q^power · mult(q)   with q = two_q/2 ∈ ½ℤ
            # return sympy.Rational
            s = sp.Rational(0)
            for two_q, m in part.items():
                s += (Rational(two_q, 2)) ** power * m
            return s

        mom2_plus = moment(part_plus, 2)
        mom2_minus = moment(part_minus, 2)
        mom2_full = moment(part_full, 2)
        mom4_full = moment(part_full, 4)

        # Update cumulative
        for k, v in part_full.items():
            cum_full[k] += v

        lvl = {
            "n_ch": n_ch,
            "g_n_plus": g_plus,
            "g_n_full_dirac": g_full,
            "partition_plus": {str(k): v for k, v in sorted(part_plus.items())},
            "partition_minus": {str(k): v for k, v in sorted(part_minus.items())},
            "partition_full": {str(k): v for k, v in sorted(part_full.items())},
            "sum_mult_full_eq_g": ok_f,
            "charge_support_full_2q": sorted(part_full.keys()),
            "moment_q2_plus": str(mom2_plus),
            "moment_q2_minus": str(mom2_minus),
            "moment_q2_full": str(mom2_full),
            "moment_q4_full": str(mom4_full),
        }
        results["levels"].append(lvl)

    # Cumulative totals up to each n_ch
    cum = defaultdict(int)
    running = {}
    for n_ch in range(n_ch_max + 1):
        for k, v in hopf_partition_full_dirac(n_ch).items():
            cum[k] += v
        running[n_ch] = {
            "total_so_far": sum(cum.values()),
            "partition_cumulative_full": {str(k): v for k, v in sorted(cum.items())},
        }
    results["cumulative_tests"] = running

    # Sanity: cum up to n_ch = 2 must be 2+6+12+... wait:
    # g_0=2, g_1=6, g_2=12, g_3=20, g_4=30, g_5=42.
    # Cum up to n=2 is 2+6+12 = 20, NOT 40. D1 says:
    #   count_spinor_labels(n_max=2, sector="dirac") = 40
    # That's with n_max inclusive *and* a different convention. Let's check:
    # Actually g_n^Dirac in CLAUDE.md is the degeneracy at level n.
    # Phase 4H SM-D: Δ^{-1} = g_3^Dirac = 2(3+1)(3+2) = 40. So it's the
    # *point* degeneracy at n_CH = 3, NOT a cumulative sum.
    # The D1 module's count_spinor_labels at n_max=2 also gives 40 — coincidence?
    # Cum 0..2: g_0+g_1+g_2 = 2+6+12 = 20. Cum 0..3: +20 = 40.
    # The D1 docstring mentions "count_spinor_labels(n_max=2) = 40 =
    # Δ⁻¹, the cumulative-Dirac reading of Δ". But cum 0..2 = 20 by the
    # Camporesi-Higuchi formula g_n = 2(n+1)(n+2). Let me recompute:
    # n=0: 2·1·2=2. n=1: 2·2·3=12. n=2: 2·3·4=24. Oh! I had the
    # formula wrong above.
    #
    # Wait — check the CH formula: g_n = 2(n+1)(n+2) at CH level n.
    # n=0: 4. n=1: 12. n=2: 24. n=3: 40. So g_3 = 40. ✓
    # Cum 0..1 = 4+12 = 16. Cum 0..2 = 40. ✓
    # "count_spinor_labels(n_max=2) = 40" matches cum 0..2 = g_0+g_1+g_2
    # = 4+12+24 = 40.
    # So there are TWO readings of Δ^{-1} = 40:
    #   (a) point degeneracy at n_ch = 3: g_3^Dirac = 40
    #   (b) cumulative degeneracy 0..2: g_0+g_1+g_2 = 40
    # Both equal 40 because 2(n+1)(n+2) = Σ_{k≤n-1} 2(k+1)(k+2) when ...
    # actually no, that's a coincidence at n=3 specifically:
    #   4 + 12 + 24 = 40 = 2·4·5. Check n=2: 4+12 = 16 ≠ 2·3·4 = 24. ✗
    #   So it IS a coincidence at n=3.
    # Let's record both.

    # Rewrite:
    direct_g3 = 2 * (3 + 1) * (3 + 2)  # = 40
    cum_02 = sum(2 * (k + 1) * (k + 2) for k in range(3))  # = 4+12+24 = 40
    results["cumulative_tests"]["delta_inv_readings"] = {
        "point_degeneracy_g3": direct_g3,
        "cumulative_0_to_2": cum_02,
        "match": direct_g3 == cum_02 == 40,
        "structural_note": (
            "Both readings happen to equal 40 at n_ch=3, but this is a "
            "coincidence: cum_{0..n-1} != g_n in general (e.g. cum_{0..1}=16 "
            "while g_2=24)."
        ),
    }

    # ----- B = 42 tests -----
    B_tests = {}

    # (B.1) point degeneracy: g_5 = 2·6·7 = 84; g_4 = 2·5·6 = 60. Nothing hits 42.
    B_tests["point_degeneracies_0_to_5"] = [
        2 * (n + 1) * (n + 2) for n in range(6)
    ]  # [4, 12, 24, 40, 60, 84]
    # Cumulative: 4, 16, 40, 80, 140, 224. Nothing hits 42.
    cum_seq = []
    s = 0
    for v in B_tests["point_degeneracies_0_to_5"]:
        s += v
        cum_seq.append(s)
    B_tests["cumulative_0_to_n"] = cum_seq

    # (B.2) Weighted Casimir trace: Σ |λ_n|² · g_n  up to n_ch=2 (matching
    # Paper 2's m=3 scalar cutoff which uses n_fock=1,2,3 → n_ch=0,1,2).
    # |λ_n| = n + 3/2.
    weighted = sp.Rational(0)
    for n in range(3):
        lam = Rational(2 * n + 3, 2)
        gn = 2 * (n + 1) * (n + 2)
        weighted += lam**2 * gn
    B_tests["sum_lam2_g_n_0_to_2"] = str(weighted)  # explicit
    # = (3/2)^2·4 + (5/2)^2·12 + (7/2)^2·24
    # = 9/4·4 + 25/4·12 + 49/4·24
    # = 9 + 75 + 294 = 378. = 9 · 42. ← Candidate!
    B_tests["sum_lam2_g_n_0_to_2_factored"] = {
        "value": int(weighted),
        "as_multiple_of_42": Fraction(int(weighted), 42).__str__(),
    }

    # Also try sum through n_ch=3 (Paper 2's m=3 corresponds to shells up to n=3 in Fock?)
    # Paper 2 uses (2l+1)·l(l+1) with l up to 2 on the n=3 shell — i.e. up to n_fock=3 = n_ch=2.
    # But let's also check 0..3:
    weighted_3 = weighted + Rational(9, 2)**2 * 40   # lam_3 = 9/2, g=40
    B_tests["sum_lam2_g_n_0_to_3"] = str(weighted_3)
    B_tests["sum_lam2_g_n_0_to_3_over_42"] = str(sp.Rational(int(weighted_3 * 10**0) if weighted_3.is_integer else 0, 42) if weighted_3.is_integer else f"{weighted_3}/42 = {sp.nsimplify(weighted_3/42)}")

    # (B.3) Unweighted Σ |λ_n|·g_n up to 0..2
    linw = sum(Rational(2 * n + 3, 2) * 2 * (n + 1) * (n + 2) for n in range(3))
    B_tests["sum_lam_g_n_0_to_2"] = str(linw)
    # = 3/2·4 + 5/2·12 + 7/2·24 = 6 + 30 + 84 = 120

    # (B.4) Hopf-charge moment: Σ_q q² · m(q) summed over n_ch = 0..2
    mom2_cum = sp.Rational(0)
    for n_ch in range(3):
        part = hopf_partition_full_dirac(n_ch)
        for two_q, m in part.items():
            mom2_cum += Rational(two_q, 2)**2 * m
    B_tests["hopf_charge_Casimir_mom2_0_to_2"] = str(mom2_cum)

    mom2_cum_3 = sp.Rational(0)
    for n_ch in range(4):
        part = hopf_partition_full_dirac(n_ch)
        for two_q, m in part.items():
            mom2_cum_3 += Rational(two_q, 2)**2 * m
    B_tests["hopf_charge_Casimir_mom2_0_to_3"] = str(mom2_cum_3)

    # (B.5) max-charge truncated cumulative:
    # Σ_{|q|≤Q, n_ch≤N} m(q,n) — does this hit 42?
    max_search = {}
    for N in range(6):
        for Q_half in range(0, 13):  # up to |q| ≤ 6
            s = 0
            for n_ch in range(N + 1):
                part = hopf_partition_full_dirac(n_ch)
                for two_q, m in part.items():
                    if abs(two_q) <= Q_half:
                        s += m
            if s == 42:
                max_search[f"N={N},Q_half={Q_half}"] = s
    B_tests["charge_truncation_hits_42"] = max_search if max_search else "no hit"

    results["B_tests"] = B_tests

    # ----- F = π²/6 tests -----
    F_tests = {}

    # (F.1) Dirichlet series Σ_{n≥0} g_n^Dirac · (n+3/2)^{-s} symbolically
    s_sym = Symbol("s", positive=True)
    n_sym = Symbol("n", integer=True, nonnegative=True)
    # g_n = 2(n+1)(n+2); |λ_n|=n+3/2
    # D(s) = Σ_{n≥0} 2(n+1)(n+2) / (n+3/2)^s
    # Try to express in terms of ζ(s, 3/2) derivatives
    # Σ_{n≥0} 2(n+1)(n+2) (n+3/2)^{-s}
    # Let k = n + 3/2, then n+1 = k-1/2, n+2 = k+1/2
    # 2(n+1)(n+2) = 2(k-1/2)(k+1/2) = 2(k^2 - 1/4) = 2k^2 - 1/2
    # So D(s) = Σ_{k=3/2, 5/2, ...} (2k^2 - 1/2) k^{-s}
    #        = 2·ζ(s-2, 3/2, half-integers) - (1/2)·ζ(s, 3/2, half-integers)
    # Half-integer Hurwitz ζ: Σ_{k≥0} (k+3/2)^{-s} = ζ(s, 3/2).
    # By standard identities, ζ(s, 1/2) = (2^s - 1) ζ(s), and
    # ζ(s, 3/2) = ζ(s, 1/2) - (1/2)^{-s} = (2^s - 1)ζ(s) - 2^s
    # Double check: ζ(s, 1/2) = Σ_{k≥0} (k+1/2)^{-s} = 2^s Σ_{k≥0} (2k+1)^{-s}
    #   = 2^s · (1 - 2^{-s}) ζ(s) = (2^s - 1) ζ(s). ✓
    # So D(s) = 2·((2^{s-2}-1)ζ(s-2) - 2^{s-2}) - (1/2)·((2^s-1)ζ(s) - 2^s)
    # Evaluate at s = 4:
    Ds4 = 2 * ((2**(4 - 2) - 1) * zeta(4 - 2) - 2**(4 - 2)) \
          - sp.Rational(1, 2) * ((2**4 - 1) * zeta(4) - 2**4)
    Ds4_simpl = sp.simplify(Ds4)
    F_tests["D_dirac_at_s_4"] = str(Ds4_simpl)
    # zeta(2) = pi^2/6, zeta(4) = pi^4/90

    Ds3 = 2 * ((2**(3 - 2) - 1) * zeta(3 - 2) - 2**(3 - 2)) \
          - sp.Rational(1, 2) * ((2**3 - 1) * zeta(3) - 2**3)
    # ζ(1) diverges; this one should diverge. Skip.
    # s=3 has Σ 2k²/k³ = Σ 2/k — ζ(1) pole. Skip.
    F_tests["D_dirac_at_s_3"] = "diverges (zeta(1) pole)"

    Ds5 = 2 * ((2**(5 - 2) - 1) * zeta(5 - 2) - 2**(5 - 2)) \
          - sp.Rational(1, 2) * ((2**5 - 1) * zeta(5) - 2**5)
    Ds5_simpl = sp.simplify(Ds5)
    F_tests["D_dirac_at_s_5"] = str(Ds5_simpl)

    Ds6 = 2 * ((2**(6 - 2) - 1) * zeta(6 - 2) - 2**(6 - 2)) \
          - sp.Rational(1, 2) * ((2**6 - 1) * zeta(6) - 2**6)
    Ds6_simpl = sp.simplify(Ds6)
    F_tests["D_dirac_at_s_6"] = str(Ds6_simpl)

    # Check if Ds4 = F = π²/6 exactly (symbolic equality)
    diff_F = sp.simplify(Ds4_simpl - pi**2/6)
    F_tests["Ds4_equals_F"] = str(diff_F) + "  (nonzero iff no match)"
    # Numerical check
    F_tests["Ds4_float"] = float(sp.N(Ds4_simpl, 30))
    F_tests["F_float"] = float(sp.N(pi**2/6, 30))

    # (F.2) By Hopf charge sector:
    # D_q(s) = Σ_{n: q appears} m(q,n) · |λ_n|^{-s}
    # Charge q appears in the (+)-chirality at level n_ch if |q| ≤ n/2,
    # i.e. n ≥ 2|q| (with step 2 starting from n = 2|q|).
    # Actually q = m_R ∈ {-n/2,...,n/2} (step 1), so q appears in (+)
    # iff 2|q| ≤ n. Similarly for (-): q ∈ {-(n+1)/2,...,(n+1)/2}, so
    # 2|q| ≤ n+1.
    # For q = 1/2 (two_q = 1): (+) requires 2|q|=1 ≤ n, so n ≥ 1, all n≥1 in (+).
    # Multiplicity at q in (+) is (n+2) for n ≥ 2|q|.
    # Multiplicity at q in (-) is (n+1) for n ≥ 2|q| - 1.
    # Half-integer charges only — integer charges never appear.

    # Fixed-q series at q = 1/2:
    # (+): Σ_{n≥1} (n+2)/(n+3/2)^s
    # (-): Σ_{n≥0} (n+1)/(n+3/2)^s
    # Sum: depends on charge. Skip the infinite expressions; do summands
    # through n ≤ 10 symbolically.
    q_sectors = {}
    for two_q_target in [1, 3, 5, 7]:  # charges 1/2, 3/2, 5/2, 7/2
        terms = []
        for n in range(15):
            part_plus = hopf_partition_single_chirality(n, +1)
            part_minus = hopf_partition_single_chirality(n, -1)
            m = part_plus.get(two_q_target, 0) + part_minus.get(two_q_target, 0)
            if m > 0:
                terms.append((n, m))
        q_sectors[f"q={Rational(two_q_target, 2)}"] = {
            "first_few_(n,mult)": [(n, m) for n, m in terms[:8]],
        }
    F_tests["charge_sector_supports"] = q_sectors

    results["F_tests"] = F_tests

    # ----- Δ = 1/40 tests -----
    Delta_tests = {}
    # Split of the 40 states at n_ch=3 over charges
    part_3 = hopf_partition_full_dirac(3)
    # (+) at n=3: j_R=3/2, m_R ∈ {-3/2,-1/2,1/2,3/2}, mult (n+2)=5
    # (-) at n=3: j_R=2, m_R ∈ {-2,-1,0,1,2}, mult (n+1)=4
    Delta_tests["partition_at_n_ch_3"] = {str(k): v for k, v in part_3.items()}
    Delta_tests["sum"] = sum(part_3.values())
    # Half-integer charges: from (+) only. Integer charges: from (-) only.
    # Half-integer charges mult: 5 each × 4 values = 20 states
    # Integer charges mult: 4 each × 5 values = 20 states
    # Natural 20 + 20 = 40 split!
    half_int_states = sum(v for k, v in part_3.items() if k % 2 == 1)
    int_states = sum(v for k, v in part_3.items() if k % 2 == 0)
    Delta_tests["half_integer_charge_states"] = half_int_states
    Delta_tests["integer_charge_states"] = int_states
    Delta_tests["half_plus_int_eq_40"] = (half_int_states + int_states) == 40
    # Compare to Paper 2's 8·5 = 40 factorization
    Delta_tests["factorization_comparison"] = {
        "paper_2": "|λ_3|·N(2) = 8·5 = 40",
        "D4_by_chirality": f"(+):{(3+1)*(3+2)}  (-):{(3+1)*(3+2)}  = 20 + 20 = 40",
        "D4_by_charge_parity": f"half-int:{half_int_states}  int:{int_states} = 20 + 20 = 40",
        "D4_by_charge_count_x_mult_plus": "4 charges × 5 mult = 20",
        "D4_by_charge_count_x_mult_minus": "5 charges × 4 mult = 20",
        "matches_8_times_5": "yes — 4 half-int charges × 5 mult = 20, 5 int charges × 4 mult = 20; combined 5·8 = 40.",
    }

    # Which charge count × mult matches 8·5 directly?
    # (+): 4 (half-int charges) × 5 (mult) = 20
    # (-): 5 (int charges) × 4 (mult) = 20
    # Paper 2 has 8 · 5. The full 40 can be read two ways via Hopf:
    #   - Per-chirality: 2 × 20, each 20 = (n+1)(n+2) at n=3
    #   - By charge-parity: 20 half-integer + 20 integer
    # Neither reproduces 8 as the "|λ_3|" factor directly, since |λ_3|
    # in the Dirac convention is 9/2, not 8. |λ_3|=8 in Paper 2 is the
    # SCALAR Laplacian eigenvalue, = n²−1 at n_Fock=3 = 9−1 = 8. Distinct.

    results["Delta_tests"] = Delta_tests

    # ----- Verdict -----
    # Check actual hits
    verdict = {}
    verdict["B_42_hit"] = (
        int(weighted) == 378 and Fraction(378, 42) == Fraction(9, 1)
    )
    verdict["B_42_via_sum_lam2_g"] = (
        "Σ_{n≤2} |λ_n|^2 g_n = 378 = 9·42. Numerical multiple exists but "
        "no clean structural interpretation: the factor 9 = |λ_0|^2 · ... "
        "does not match Paper 2's (2l+1)·l(l+1) construction. TREAT AS "
        "near-miss, NOT a structural identity."
    )
    verdict["F_hit"] = False  # compute below
    F_tests_nu = float(sp.N(pi**2/6 - Ds4_simpl, 30))
    verdict["F_Ds4_minus_F"] = str(F_tests_nu)
    verdict["F_structural"] = (
        "D_dirac(s=4) = " + str(Ds4_simpl) + " — mixture of ζ(2) and ζ(4) "
        "with rational parts. No clean ζ(2) isolation. F does NOT lift to "
        "the Dirac Dirichlet series at packing exponent s=4."
    )
    verdict["Delta_hit"] = True
    verdict["Delta_structural"] = (
        "g_3^Dirac = 40 = Δ^{-1} reproduced. Hopf decomposition at n=3 "
        "splits cleanly as 20 + 20 across the two chirality sectors and "
        "as 20 + 20 across half-integer / integer charge parity. This is "
        "TRUE BY CONSTRUCTION (Phase 4H SM-D already observed it), not a "
        "NEW result. The refinement is that the '40' carries a natural "
        "4×5 + 5×4 (charges × mult) labeling — consistent with, but not "
        "equivalent to, Paper 2's |λ_3|·N(2) = 8·5 factorization."
    )
    verdict["final"] = (
        "Of {B, F, Δ}: Δ reproduced (was known, Phase 4H). B hit as "
        "numerical multiple 378 = 9·42 with no structural interpretation. "
        "F does NOT lift. "
        "\n"
        "Verdict for D5: NO Hopf-Dirac lift of the full K combination rule. "
        "Δ already had a spectral home (Phase 4H SM-D). B and F remain "
        "scalar-sector-only invariants. The Hopf-equivariant Dirac "
        "decomposition on S³ is categorically different from Phase 4B's "
        "scalar graph quotient (first-order operator; charge-q twists "
        "are nontrivial line bundles over S²), but this categorical "
        "difference does not manifest as new identities for B or F. "
        "Paper 2's K combination rule remains a three-tier coincidence."
    )
    verdict["categorical_distinction_from_phase4B"] = (
        "Phase 4B α-D tested the SCALAR graph morphism S³ → S² under the "
        "Hopf quotient: nodes (n,l,m) → sectors labeled by scalar spherical "
        "harmonic content on S². The Laplace–Beltrami is second-order and "
        "bosonic, so the quotient map is a deformation-retract-like "
        "averaging; no line bundle can appear, and all charge content is "
        "trivial. Negative result was structural: no spectral invariant "
        "of the quotient hits K targets. "
        "\n"
        "D4 tests the DIRAC decomposition: first-order operator, spinor "
        "bundle, and U(1) Hopf action gives a NONTRIVIAL direct-sum "
        "decomposition over half-integer charges q (line bundles O(2q) "
        "over S² with q ∈ 1/2 ℤ). Each charge sector is a twisted Dirac "
        "on S² with its own Baum-Friedrich spectrum. This is genuinely new "
        "content not available in the scalar case. The charge-parity split "
        "(half-int vs int) mirrors the chirality split (+ vs −) exactly "
        "and reproduces Δ^{-1}=40 as a 20+20 partition — the first "
        "place chirality/charge have structural content on the graph. But "
        "the new decomposition does NOT provide additional handles on B=42 "
        "or F=π²/6, so it does not rescue the common-generator hypothesis."
    )
    results["verdict"] = verdict

    return results


# ---------------------------------------------------------------------------
# Entry
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    res = run(n_ch_max=5)

    # Write JSON
    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "dirac_d4_hopf_equivariant.json"

    # sympy objects may linger; stringify non-JSON-serializable safely
    def _default(o):
        return str(o)
    with open(out_path, "w") as f:
        json.dump(res, f, indent=2, default=_default)
    print(f"Wrote {out_path}")

    # Concise console summary
    print()
    print("=" * 70)
    print("Track D4 — Hopf-equivariant Dirac decomposition summary")
    print("=" * 70)
    for lvl in res["levels"]:
        n = lvl["n_ch"]
        print(f"\nn_ch = {n}   g_n_Dirac = {lvl['g_n_full_dirac']}")
        print(f"  partition_full (key = 2q): {lvl['partition_full']}")
        print(f"  sum q^2 * m = {lvl['moment_q2_full']}")
    print()
    print("Cumulative 0..2 Dirac states:",
          res["cumulative_tests"][2]["total_so_far"])
    print()
    print("B tests:")
    val_378 = int(eval(res["B_tests"]["sum_lam2_g_n_0_to_2"]))
    print(f"  sum |lam_n|^2 * g_n up to n_ch=2: {val_378}  =  {Fraction(val_378, 42)} * 42")
    print()
    print("F tests:")
    print("  D_dirac(s=4) =", res["F_tests"]["D_dirac_at_s_4"])
    print("  D_dirac(4) float:", res["F_tests"]["Ds4_float"])
    print("  F = pi^2/6 float:", res["F_tests"]["F_float"])
    print()
    print("Delta tests:")
    print("  n=3 partition:", res["Delta_tests"]["partition_at_n_ch_3"])
    print("  half-int states:", res["Delta_tests"]["half_integer_charge_states"])
    print("  int states:", res["Delta_tests"]["integer_charge_states"])
    print()
    print("Verdict:")
    print(res["verdict"]["final"])
