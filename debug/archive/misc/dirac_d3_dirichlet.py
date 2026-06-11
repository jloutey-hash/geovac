"""Track D3: Dirichlet-series search for a Dirac analog of F = π²/6 on S³.

Algebraic-first (pure sympy). The scalar result from Phase 4F (α-J) was:

    D_{n²}(s) = Σ_{n≥1} n² · n^{-s} = ζ(s-2)     (shift of argument)
    D_{n²}(4) = ζ(2) = π²/6 = F                     (EXACT symbolic)

with weight n² = Fock scalar degeneracy, exponent s=4 = d_max (packing dim).

D3 tests the Dirac analog: degeneracy g_n^Dirac = 2(n+1)(n+2) (CH indexing,
n=0,1,2,...) or, in Fock indexing m = n+1 ≥ 1:

    g_m = 2 m (m+1)

For each candidate series we compute the EXACT symbolic closed form in sympy
and evaluate at s ∈ {3,4,5,6}. We look for a clean π² · (rational),
ζ(2) · (rational), or an additive decomposition isolating ζ(2).
"""

from __future__ import annotations

import json
from sympy import (
    Symbol, summation, zeta, oo, Rational, simplify, nsimplify,
    pi, Integer, together, factor, S, latex
)

n = Symbol('n', integer=True, nonnegative=True)   # CH index, ≥ 0
m = Symbol('m', integer=True, positive=True)      # Fock index = n+1, ≥ 1
s = Symbol('s', positive=True)                    # Dirichlet exponent (symbolic)

# ---------------------------------------------------------------------------
# Reference scalar result (Phase 4F α-J)
# ---------------------------------------------------------------------------
# D_{n²}(s) = Σ_{m≥1} m² / m^s = ζ(s-2)
D_nsq_sym = summation(m**2 / m**s, (m, 1, oo))     # sympy returns zeta(s-2)
F_target = zeta(2)                                  # π²/6

# ---------------------------------------------------------------------------
# Candidate (a): Dirac-degeneracy series in FOCK indexing m = n+1
#   D_dirac(s) = Σ_{n≥0} 2(n+1)(n+2) / (n+1)^s
#              = 2 Σ_{m≥1} m(m+1) / m^s
#              = 2 [ζ(s-2) + ζ(s-1)]
# ---------------------------------------------------------------------------
D_dirac_fock = summation(2 * m * (m + 1) / m**s, (m, 1, oo))

# ---------------------------------------------------------------------------
# Candidate (b): Dirac-degeneracy with CH eigenvalue (n + 3/2) as denominator
#   D_dirac^spec(s) = Σ_{n≥0} 2(n+1)(n+2) / (n + 3/2)^s
# This is a Hurwitz zeta. Expand numerator around m = n + 3/2:
#   (n+1)(n+2) = (m - 1/2)(m + 1/2) = m² - 1/4
# So Σ 2(m² - 1/4) / m^s (m running over 3/2, 5/2, ...)
#   = 2 ζ(s-2, 3/2) - (1/2) ζ(s, 3/2)
# ---------------------------------------------------------------------------
from sympy import zeta as z
# Hurwitz ζ(s, a) = zeta(s, a) in sympy
D_dirac_spec = 2 * z(s - 2, Rational(3, 2)) - Rational(1, 2) * z(s, Rational(3, 2))

# ---------------------------------------------------------------------------
# Candidate (c): Weyl (single-chirality) — half of (a) and (b)
# ---------------------------------------------------------------------------
D_weyl_fock = D_dirac_fock / 2          # = ζ(s-2) + ζ(s-1)
D_weyl_spec = D_dirac_spec / 2

# ---------------------------------------------------------------------------
# Candidate (d): Isolate the "extra Dirac-over-scalar" content
#   g_m = 2 m(m+1) = 2m² + 2m,  scalar degeneracy is m² (Fock)
#   Δ_Dirac-over-scalar weight = 2m² + 2m - m² = m² + 2m = m(m+2)
#   So Σ m(m+2)/m^s = ζ(s-2) + 2 ζ(s-1)
#
#   Alternatively: Dirac − 2·scalar:  (2m² + 2m) − 2m² = 2m
#   Σ 2m / m^s = 2 ζ(s-1)
# ---------------------------------------------------------------------------
delta_over_scalar = summation((m**2 + 2*m) / m**s, (m, 1, oo))
dirac_minus_2scalar = summation(2*m / m**s, (m, 1, oo))

# ---------------------------------------------------------------------------
# Candidate (e): Regularized ζ-differences (Phase 4F subtask 2 style)
#   Scalar B-analog: 2 D_B(s) = ζ(s-4) - ζ(s-2)   (from Σ (n²-1) weights)
#   Dirac B-analog: use (n+1)(n+2) - something? Try
#     D_B^Dirac(s) = Σ [(n+1)(n+2) - (n+1)²] / (n+1)^s
#                  = Σ (n+1) / (n+1)^s  = ζ(s-1)
# ---------------------------------------------------------------------------
D_B_dirac = summation(m / m**s, (m, 1, oo))     # = ζ(s-1)

# ---------------------------------------------------------------------------
# Evaluate all candidates at s ∈ {3, 4, 5, 6}
# ---------------------------------------------------------------------------
candidates = {
    "a_D_dirac_fock":         D_dirac_fock,          # 2[ζ(s-2)+ζ(s-1)]
    "b_D_dirac_spec":         D_dirac_spec,          # Hurwitz combination
    "c_D_weyl_fock":          D_weyl_fock,
    "c_D_weyl_spec":          D_weyl_spec,
    "d_delta_over_scalar":    delta_over_scalar,     # ζ(s-2)+2ζ(s-1)
    "d_dirac_minus_2scalar":  dirac_minus_2scalar,   # 2ζ(s-1)
    "e_D_B_dirac":            D_B_dirac,             # ζ(s-1)
    "ref_scalar_D_nsq":       D_nsq_sym,             # ζ(s-2)  [F at s=4]
}

results = {}
for name, expr in candidates.items():
    entry = {}
    entry["symbolic_in_s"] = str(simplify(expr))
    for s_val in [3, 4, 5, 6]:
        try:
            val = simplify(expr.subs(s, s_val))
            entry[f"s={s_val}"] = {
                "exact":   str(val),
                "latex":   latex(val),
            }
            try:
                entry[f"s={s_val}"]["float"] = float(val.evalf(25))
            except Exception as e:
                entry[f"s={s_val}"]["float"] = f"diverges: {e}"
            ratio = simplify(val / F_target)
            entry[f"s={s_val}"]["ratio_to_F"] = str(ratio)
            q_cand = simplify(val / pi**2)
            entry[f"s={s_val}"]["over_pi_sq"] = str(q_cand)
        except Exception as e:
            entry[f"s={s_val}"] = {"error": str(e)}
    results[name] = entry

# ---------------------------------------------------------------------------
# Pretty print
# ---------------------------------------------------------------------------
print("=" * 70)
print("Track D3: Dirac-Dirichlet search for F = π²/6 analog")
print("=" * 70)
for name, entry in results.items():
    print(f"\n### {name}")
    print(f"  symbolic: {entry['symbolic_in_s']}")
    for s_val in [3, 4, 5, 6]:
        k = f"s={s_val}"
        if 'error' in entry[k]:
            print(f"  {k}: ERROR {entry[k]['error']}")
        else:
            print(f"  {k}: {entry[k]['exact']}")
            print(f"       ratio / F = {entry[k]['ratio_to_F']}")
            print(f"       over π²   = {entry[k]['over_pi_sq']}")

# ---------------------------------------------------------------------------
# Targeted check: can F be extracted from any additive decomposition?
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("F-extraction tests")
print("=" * 70)

# At s=4, D_dirac_fock(4) = 2[ζ(2) + ζ(3)] = 2F + 2ζ(3)
d4 = simplify(D_dirac_fock.subs(s, 4))
print(f"\nD_dirac_fock(4) = {d4}")
print(f"  = 2·F + 2·ζ(3)?  Check: D_dirac_fock(4) - 2F - 2ζ(3) = "
      f"{simplify(d4 - 2*F_target - 2*zeta(3))}")

# Isolate F via combination: (1/2)·[D_dirac_fock(4) - 2ζ(3)] = F ?
iso = simplify(D_dirac_fock.subs(s, 4) / 2 - zeta(3))
print(f"\n(1/2)·D_dirac_fock(4) - ζ(3) = {iso}")
print(f"  equals F?  {simplify(iso - F_target) == 0}")

# At s=5: D_dirac_fock(5) = 2[ζ(3) + ζ(4)] = 2ζ(3) + π⁴/45
d5 = simplify(D_dirac_fock.subs(s, 5))
print(f"\nD_dirac_fock(5) = {d5}")

# At s=6: D_dirac_fock(6) = 2[ζ(4) + ζ(5)] = π⁴/45 + 2ζ(5)
d6 = simplify(D_dirac_fock.subs(s, 6))
print(f"D_dirac_fock(6) = {d6}")

# Try: is there an s for which D_dirac(s) is a pure rational multiple of F?
# Need ζ(s-1) to vanish or to be commensurate with ζ(s-2)=F-family
# ζ(s-1) rational ⇔ s-1 ∈ {0, negative even} ⇔ s=1 (pole) or s=-1,-3,...
# So NO positive s gives clean F from (a).

# Candidate (d) dirac_minus_2scalar = 2ζ(s-1) — same story.
# Candidate (e) D_B_dirac = ζ(s-1) at s=3 → ζ(2) = F !  EXACT HIT.
eB3 = simplify(D_B_dirac.subs(s, 3))
print(f"\n*** D_B_dirac(s=3) = ζ(2) = F  ***")
print(f"    D_B_dirac(3) = {eB3}")
print(f"    equals F?  {simplify(eB3 - F_target) == 0}")

# What is this series geometrically?
#   D_B_dirac = Σ_{m≥1} m / m^s
# with weight m = n+1 = Fock principal quantum number
# at s = 3. Weight = PRINCIPAL QUANTUM NUMBER (not degeneracy).
# Exponent = 3 = dim(S³), not d_max = 4.

# Compare: scalar F hit was weight n² (Fock degeneracy) at s=4=d_max.
# Dirac F hit is weight m (principal qn) at s=3=dim(S³).
# Ask: is there a weight tied to the Dirac degeneracy specifically
# (not to m alone) that hits F at s=4 = d_max?

# Weight to test: g_m^Dirac / (something rational in m) at s=4
# Candidate: (D_dirac_fock(s) - D_B_dirac(s)) at s=4
#   = 2[ζ(s-2)+ζ(s-1)] - ζ(s-1) = 2ζ(s-2) + ζ(s-1)
# At s=4: 2F + ζ(3) — mixed.

# Candidate: D_dirac_fock(s) - 2·D_B_dirac(s) at s=4
#   = 2ζ(s-2) + 2ζ(s-1) - 2ζ(s-1) = 2ζ(s-2)
# At s=4: 2F !  Clean multiple of F.
combo = simplify(D_dirac_fock.subs(s, 4) - 2*D_B_dirac.subs(s, 4))
print(f"\n*** D_dirac_fock(4) - 2·D_B_dirac(4) = 2F  ***")
print(f"    value = {combo}")
print(f"    equals 2F?  {simplify(combo - 2*F_target) == 0}")
# But this combination = Σ 2m(m+1)/m⁴ - Σ 2m/m⁴ = Σ 2m² /m⁴ = 2 ζ(2)
# which is just 2·D_{n²}(4) = 2·(scalar F). NOT new content.

# SUMMARY FOR JSON
summary = {
    "scalar_F_identity": {
        "series": "D_{n^2}(s) = sum m^2/m^s = zeta(s-2)",
        "at_s_4": str(simplify(D_nsq_sym.subs(s, 4))),
        "equals_F": True,
    },
    "dirac_fock_at_s": {
        "3": "2*[zeta(1)+zeta(2)] -- DIVERGES (zeta(1) pole)",
        "4": str(simplify(D_dirac_fock.subs(s, 4))),  # 2[ζ(2)+ζ(3)] = 2F + 2ζ(3)
        "5": str(simplify(D_dirac_fock.subs(s, 5))),
        "6": str(simplify(D_dirac_fock.subs(s, 6))),
    },
    "findings": [
        "D_dirac_fock(s) = 2·zeta(s-2) + 2·zeta(s-1)  (EXACT)",
        "At s=4 (packing d_max): 2·F + 2·zeta(3) — mixed transcendental with zeta(3)",
        "NO positive integer s gives a clean rational multiple of F from D_dirac.",
        "D_B_dirac(s) = zeta(s-1) hits F at s=3 = dim(S^3), not d_max.",
        "D_B_dirac with weight m (principal qn, NOT Dirac degeneracy) at s=3 gives F.",
        "This is NOT a Dirac-specific lift: zeta(s-1) at s=3 is trivial from any series",
        "with weight m^1. The Dirac degeneracy g_m = 2m(m+1) does NOT collapse to F.",
    ],
    "verdict": "F does NOT lift cleanly to the Dirac sector on S^3. "
               "The Dirac-degeneracy Dirichlet series at s=d_max=4 gives "
               "2F + 2*zeta(3), which mixes F with the independent transcendental "
               "zeta(3). No additive/subtractive combination isolates F as a "
               "Dirac-specific quantity: the only extractions (e.g. subtracting "
               "2·D_B_dirac) recover the SCALAR n^2-weighted series, not a new "
               "Dirac identity. Phase 4F's F = zeta(s-2) identity is specific "
               "to the scalar Fock degeneracy n^2.",
    "series_results": results,
}

out_path = r"C:/Users/jlout/Desktop/Project_Geometric/debug/data/dirac_d3_dirichlet.json"
with open(out_path, "w") as f:
    json.dump(summary, f, indent=2)
print(f"\nWrote {out_path}")
