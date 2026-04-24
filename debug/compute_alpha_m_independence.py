"""
Sprint A of the alpha-program reframe (CLAUDE.md §1.7 WH5).

Compute K(m) = pi(B(m) + F(m) - Delta(m)) at m = 2, 3, 4, 5, 6, 10
under the canonical structural identifications of Phases 4B/4F/4G/4H.

Outputs:
  debug/data/alpha_m_independence.json
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp
from mpmath import mp, mpf, pi as mp_pi

# 50 dps for numerical floats
mp.dps = 60

# ---------------------------------------------------------------------
# Symbolic forms (from Paper 2, Phases 4B-4H)
# ---------------------------------------------------------------------

m = sp.Symbol("m", positive=True, integer=True)

# B(m): Casimir trace on (n,l) lattice, Paper 2 Eq. (eq:Bclosed).
#   B(m) = sum_{n=1}^{m} n^2 (n^2 - 1) / 2
#        = m (m-1) (m+1) (m+2) (2m+1) / 20
B_sym = m * (m - 1) * (m + 1) * (m + 2) * (2 * m + 1) / 20

# F: infinite Dirichlet series of the Fock degeneracy at the packing
# exponent s = d_max = 4 (Phase 4F).
#   F = D_{n^2}(s = d_max) = sum_n n^2 / n^4 = zeta(2) = pi^2/6
# d_max = 4 is the Paper 0 packing-axiom maximum valence; it is NOT m.
# We treat F as m-INDEPENDENT (the canonical reading).
F_sym = sp.pi**2 / 6

# Delta(m): two canonical forms.
#   Form A (Paper 2, eq:boundary):    1 / [|lambda_m| * N(m-1)]
#                                  =  1 / [(m^2 - 1) * m(m-1)(2m-1)/6]
#   Form B (Phase 4H SM-D):           1 / g_m^Dirac
#                                  =  1 / [2 (m+1) (m+2)]
delta_A_sym = 1 / ((m**2 - 1) * m * (m - 1) * (2 * m - 1) / 6)
delta_B_sym = 1 / (2 * (m + 1) * (m + 2))

# Cubic alpha-from-K closure: alpha^3 - K*alpha + 1 = 0; smallest positive root.
def alpha_inv_from_K(K_value: mpf) -> mpf:
    """
    Solve alpha^3 - K*alpha + 1 = 0 for the smallest positive real root,
    then return 1/alpha.  Cardano's formula or sympy nsolve.
    """
    a = sp.Symbol("a", real=True)
    roots = sp.nroots(a**3 - sp.Float(str(K_value), 60) * a + 1, n=60)
    pos = [r for r in roots if sp.im(r) == 0 and r > 0]
    pos.sort()
    if not pos:
        return None
    smallest = pos[0]
    return mpf(str(1 / smallest))


# ---------------------------------------------------------------------
# CODATA 2022 reference
# ---------------------------------------------------------------------
ALPHA_INV_CODATA_2022 = mpf("137.035999084")  # CODATA 2022 recommended


# ---------------------------------------------------------------------
# Compute K(m) for several m, both Delta forms
# ---------------------------------------------------------------------

def K_of_m(m_val: int, delta_form: str) -> dict:
    """Compute B(m), F, Delta(m), K(m) symbolically and numerically."""
    B_m = sp.simplify(B_sym.subs(m, m_val))
    F_m = F_sym  # m-independent
    if delta_form == "A":
        delta_m = sp.simplify(delta_A_sym.subs(m, m_val))
    elif delta_form == "B":
        delta_m = sp.simplify(delta_B_sym.subs(m, m_val))
    else:
        raise ValueError(delta_form)

    K_sym = sp.pi * (B_m + F_m - delta_m)
    K_sym_simplified = sp.simplify(K_sym)
    K_num = mpf(str(sp.N(K_sym_simplified, 60)))

    delta_K_codata = K_num - ALPHA_INV_CODATA_2022

    # Also compute alpha^-1 from K via the cubic closure
    try:
        alpha_inv_from_cubic = alpha_inv_from_K(K_num)
        delta_alpha_codata = (
            alpha_inv_from_cubic - ALPHA_INV_CODATA_2022
            if alpha_inv_from_cubic is not None
            else None
        )
    except Exception as e:
        alpha_inv_from_cubic = None
        delta_alpha_codata = None

    return {
        "m": m_val,
        "delta_form": delta_form,
        "B_sym": str(B_m),
        "B_num": float(sp.N(B_m, 30)),
        "F_sym": str(F_m),
        "F_num": float(sp.N(F_m, 30)),
        "delta_sym": str(delta_m),
        "delta_num": float(sp.N(delta_m, 30)),
        "K_sym": str(K_sym_simplified),
        "K_num_50dps": mp.nstr(K_num, 50),
        "K_minus_alphainv_codata": mp.nstr(delta_K_codata, 25),
        "alpha_inv_from_cubic_50dps": (
            mp.nstr(alpha_inv_from_cubic, 50)
            if alpha_inv_from_cubic is not None
            else None
        ),
        "alpha_inv_minus_codata": (
            mp.nstr(delta_alpha_codata, 25)
            if delta_alpha_codata is not None
            else None
        ),
    }


# ---------------------------------------------------------------------
# Subtask 1: clarify variable identity
# ---------------------------------------------------------------------
# In Paper 2 §III, B(m) is the cumulative Casimir sum to shell m.
# In Paper 2 §III, Delta uses |lambda_m| (the gap eigenvalue at shell m)
#   and N(m-1) (the cumulative state count up to shell m-1).
# So both B and Delta are indexed by the SAME truncation cutoff m.
# F is fixed by D_{n^2}(s = d_max = 4); d_max is a Paper 0 packing
# constant, NOT m.  So F has no m-dependence in the canonical reading.
#
# The m=3 selection rule fixes m via B(m)/N(m) = dim(S^3) = 3,
# uniquely solved at m=3.  This is the ONLY structural reason Paper 2
# evaluates K at m=3.
#
# Whether Delta uses Form A (canonical, |lambda_m|*N(m-1)) or
# Form B (Phase 4H, g_m^Dirac) is a question we test below: they
# AGREE at m=3 (both give 40^{-1}=1/40) but DISAGREE at other m.

# ---------------------------------------------------------------------
# Subtask 2-3-5: K(m) tables
# ---------------------------------------------------------------------

results: dict = {
    "metadata": {
        "sprint": "Sprint A — alpha-program reframe (WH5)",
        "date": "2026-04-18",
        "purpose": "Compute K(m) = pi(B(m) + F - Delta(m)) at multiple m"
                   " to characterize the m-dependence and test whether"
                   " the m=3 hit is a single-point coincidence,"
                   " convergent series, or growing object.",
        "alpha_inv_codata_2022": str(ALPHA_INV_CODATA_2022),
        "F_form": "F = pi^2/6  (m-independent, fixed by d_max=4)",
        "Delta_forms": {
            "A_canonical_paper2": "Delta(m) = 1 / [|lambda_m| * N(m-1)]"
                                  " = 6 / [(m^2-1) m(m-1)(2m-1)]",
            "B_phase4H_dirac":    "Delta(m) = 1 / g_m^Dirac"
                                  " = 1 / [2(m+1)(m+2)]",
        },
        "agreement_at_m3": "both Forms A and B give 1/40 at m=3",
        "B_form": "B(m) = m(m-1)(m+1)(m+2)(2m+1)/20",
    },
    "K_table_form_A": [],
    "K_table_form_B": [],
    "delta_form_AB_comparison": [],
    "sign_pattern_test_m3": {},
    "sign_pattern_test_m4": {},
}

# Sub-task 2 + 3 + 5: K(m) at m = 2, 3, 4, 5, 6, 10 for both Delta forms
m_values = [2, 3, 4, 5, 6, 10]
for m_val in m_values:
    results["K_table_form_A"].append(K_of_m(m_val, "A"))
    results["K_table_form_B"].append(K_of_m(m_val, "B"))

# Compare Delta(m) Form A vs Form B
for m_val in m_values:
    d_A = sp.simplify(delta_A_sym.subs(m, m_val))
    d_B = sp.simplify(delta_B_sym.subs(m, m_val))
    results["delta_form_AB_comparison"].append({
        "m": m_val,
        "Delta_A": str(d_A),
        "Delta_A_num": float(sp.N(d_A, 30)),
        "Delta_B": str(d_B),
        "Delta_B_num": float(sp.N(d_B, 30)),
        "agree": sp.simplify(d_A - d_B) == 0,
        "ratio_A_over_B_num": float(sp.N(d_A / d_B, 30)),
    })

# Sub-task 5: sign pattern test at m=3 and m=4
def sign_test(m_val: int) -> dict:
    """Try all 8 sign patterns of (B, F, Delta) for K(m) and report which closest hits 1/alpha."""
    B_m = float(sp.N(B_sym.subs(m, m_val), 30))
    F_v = float(sp.N(F_sym, 30))
    d_A = float(sp.N(delta_A_sym.subs(m, m_val), 30))
    target = float(ALPHA_INV_CODATA_2022)
    pi_v = float(sp.N(sp.pi, 30))
    out = []
    for sB in (1, -1):
        for sF in (1, -1):
            for sD in (1, -1):
                K_val = pi_v * (sB * B_m + sF * F_v + sD * d_A)
                rel = abs(K_val - target) / target
                pattern = (
                    ("+" if sB > 0 else "-")
                    + ("+" if sF > 0 else "-")
                    + ("+" if sD > 0 else "-")
                )
                out.append({
                    "pattern_BFΔ": pattern,
                    "K": K_val,
                    "K_minus_alphainv_codata": K_val - target,
                    "rel_err_vs_codata": rel,
                })
    out.sort(key=lambda r: r["rel_err_vs_codata"])
    return {
        "m": m_val,
        "all_patterns_sorted_by_closeness": out,
        "best_pattern": out[0]["pattern_BFΔ"],
        "best_rel_err": out[0]["rel_err_vs_codata"],
    }


results["sign_pattern_test_m3"] = sign_test(3)
results["sign_pattern_test_m4"] = sign_test(4)


# ---------------------------------------------------------------------
# Subtask 4 helpers: monotonicity, convergence character
# ---------------------------------------------------------------------
def characterize(table: list[dict], label: str) -> dict:
    Kvals = [mpf(row["K_num_50dps"]) for row in table]
    diffs_codata = [mpf(row["K_minus_alphainv_codata"]) for row in table]
    growth_ratios = []
    for i in range(1, len(Kvals)):
        if Kvals[i - 1] != 0:
            growth_ratios.append(float(Kvals[i] / Kvals[i - 1]))
    return {
        "label": label,
        "m_values": m_values,
        "K_values_50dps": [mp.nstr(k, 50) for k in Kvals],
        "K_minus_codata_25dps": [mp.nstr(d, 25) for d in diffs_codata],
        "growth_ratios_K(m+1)/K(m)": growth_ratios,
        "monotonicity": (
            "monotonic increasing" if all(Kvals[i] > Kvals[i - 1]
                                           for i in range(1, len(Kvals)))
            else (
                "monotonic decreasing" if all(Kvals[i] < Kvals[i - 1]
                                               for i in range(1, len(Kvals)))
                else "non-monotonic"
            )
        ),
        "min_abs_diff_idx": int(min(range(len(diffs_codata)),
                                     key=lambda i: abs(diffs_codata[i]))),
        "min_abs_diff_m": int(m_values[min(range(len(diffs_codata)),
                                            key=lambda i: abs(diffs_codata[i]))]),
        "min_abs_diff_value": mp.nstr(
            diffs_codata[
                min(range(len(diffs_codata)),
                    key=lambda i: abs(diffs_codata[i]))
            ], 25),
    }


results["characterization_form_A"] = characterize(
    results["K_table_form_A"], "K(m) using Delta = Paper 2 canonical")
results["characterization_form_B"] = characterize(
    results["K_table_form_B"], "K(m) using Delta = 1/g_m^Dirac (Phase 4H)")


# ---------------------------------------------------------------------
# Save and report
# ---------------------------------------------------------------------
out_path = Path("debug/data/alpha_m_independence.json")
out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(json.dumps(results, indent=2, default=str))

# Print summary table
print(f"\nCODATA 2022:  alpha^-1 = {ALPHA_INV_CODATA_2022}")
print()
print("=" * 88)
print("K(m) using Form A  (Delta = 1/[|lambda_m| * N(m-1)]; Paper 2 canonical)")
print("=" * 88)
print(f"{'m':>3}  {'B(m)':>12}  {'F':>14}  {'Delta(m)':>14}  {'K(m)':>22}  {'K - 1/alpha (CODATA)':>28}")
for row in results["K_table_form_A"]:
    print(f"{row['m']:>3}  {row['B_num']:>12g}  {row['F_num']:>14.10g}  "
          f"{row['delta_num']:>14.6g}  {row['K_num_50dps'][:22]:>22}  "
          f"{row['K_minus_alphainv_codata'][:28]:>28}")

print()
print("=" * 88)
print("K(m) using Form B  (Delta = 1/g_m^Dirac; Phase 4H SM-D)")
print("=" * 88)
print(f"{'m':>3}  {'B(m)':>12}  {'F':>14}  {'Delta(m)':>14}  {'K(m)':>22}  {'K - 1/alpha (CODATA)':>28}")
for row in results["K_table_form_B"]:
    print(f"{row['m']:>3}  {row['B_num']:>12g}  {row['F_num']:>14.10g}  "
          f"{row['delta_num']:>14.6g}  {row['K_num_50dps'][:22]:>22}  "
          f"{row['K_minus_alphainv_codata'][:28]:>28}")

print()
print("Form A vs Form B comparison of Delta(m):")
print(f"{'m':>3}  {'Delta_A':>16}  {'Delta_B':>16}  {'ratio A/B':>12}  {'agree?':>8}")
for row in results["delta_form_AB_comparison"]:
    print(f"{row['m']:>3}  {row['Delta_A_num']:>16.10g}  {row['Delta_B_num']:>16.10g}"
          f"  {row['ratio_A_over_B_num']:>12.6g}  {'YES' if row['agree'] else 'NO':>8}")

print()
print("Sign pattern test at m=3 (top 4):")
for r in results["sign_pattern_test_m3"]["all_patterns_sorted_by_closeness"][:4]:
    print(f"  {r['pattern_BFΔ']}: K = {r['K']:>22.15g}  rel_err = {r['rel_err_vs_codata']:.3e}")

print()
print("Sign pattern test at m=4 (top 4):")
for r in results["sign_pattern_test_m4"]["all_patterns_sorted_by_closeness"][:4]:
    print(f"  {r['pattern_BFΔ']}: K = {r['K']:>22.15g}  rel_err = {r['rel_err_vs_codata']:.3e}")

print()
print("Characterization (Form A):")
for k, v in results["characterization_form_A"].items():
    if k != "K_values_50dps" and k != "K_minus_codata_25dps":
        print(f"  {k}: {v}")
print()
print("Characterization (Form B):")
for k, v in results["characterization_form_B"].items():
    if k != "K_values_50dps" and k != "K_minus_codata_25dps":
        print(f"  {k}: {v}")

print()
print(f"Saved JSON to {out_path}")
