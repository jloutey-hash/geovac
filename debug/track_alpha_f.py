"""
Track alpha-F: Delta = 1/(B - N_init) Structural Test
Phase 4C alpha sprint.

Tests whether Delta = 1/(|lambda_{n_max}| * N(n_max - 1)) coincides with
1/(B(n_max) - N_init) where N_init = 2 (Paper 0 Axiom 1), and if so whether
the identity holds polynomially in n_max.

Closed forms (verified against Paper 2):
  |lambda_n| = n^2 - 1                        (S^3 Laplace-Beltrami)
  N(m)       = m(m+1)(2m+1)/6                 (total state count through shell m)
  B(m)       = sum_{n=1}^{m} sum_{l=0}^{n-1} (2l+1) * l*(l+1)
"""
import json
from pathlib import Path
import sympy as sp

OUT_DATA = Path("debug/data/track_alpha_phase4c/track_f_delta_identity.json")
OUT_ANALYSIS = Path("debug/data/track_alpha_phase4c/track_f_analysis.md")

N_INIT = 2  # Paper 0 Axiom 1

# ---------------------------------------------------------------------------
# Closed forms
# ---------------------------------------------------------------------------
m = sp.symbols('m', positive=True, integer=True)
n = sp.symbols('n', positive=True, integer=True)
l = sp.symbols('l', nonnegative=True, integer=True)


def lam_abs(nn):
    """|lambda_n| = n^2 - 1."""
    return nn**2 - 1


def N_count(mm):
    """N(m) = m(m+1)(2m+1)/6, the sum of 1^2 + 2^2 + ... + m^2."""
    return mm * (mm + 1) * (2 * mm + 1) / sp.Integer(6)


def B_casimir(mm):
    """B(m) = sum_{n=1}^{m} sum_{l=0}^{n-1} (2l+1)*l*(l+1).

    Inner sum over l from 0 to n-1 of (2l+1)*l*(l+1).
    """
    inner = sp.summation((2*l + 1) * l * (l + 1), (l, 0, n - 1))
    total = sp.summation(inner, (n, 1, mm))
    return sp.simplify(total)


# Sanity: B(3) should equal 42
B3_check = B_casimir(3)
assert B3_check == 42, f"B(3) expected 42, got {B3_check}"

# Closed-form polynomial for B(m)
B_poly_m = sp.simplify(B_casimir(m))
print("B(m) closed form:", sp.expand(B_poly_m))

# Delta form from Paper 2: Delta = 1 / (|lambda_m| * N(m-1))
def delta_paper2(mm):
    return sp.Rational(1) / (lam_abs(mm) * N_count(mm - 1))


# Candidate: Delta = 1 / (B(m) - N_init)
def delta_candidate(mm):
    return sp.Rational(1) / (B_casimir(mm) - N_INIT)


# ---------------------------------------------------------------------------
# SUBTASK 1: Algebraic equivalence test
# ---------------------------------------------------------------------------
# LHS(m) = |lambda_m| * N(m-1)   (denominator of Paper 2 Delta)
# RHS(m) = B(m) - N_init         (denominator of candidate Delta)
# Identity holds iff LHS(m) == RHS(m).

rows = []
for mm in range(1, 9):
    lhs = lam_abs(mm) * N_count(mm - 1)
    rhs = B_casimir(mm) - N_INIT
    diff = sp.simplify(lhs - rhs)
    delta_p2 = delta_paper2(mm) if lhs != 0 else None
    delta_cand = sp.Rational(1, int(rhs)) if rhs != 0 else None
    rows.append({
        "m": mm,
        "lambda_abs": int(lam_abs(mm)),
        "N_m_minus_1": int(N_count(mm - 1)),
        "LHS": int(lhs),
        "B_m": int(B_casimir(mm)),
        "RHS": int(rhs),
        "LHS_minus_RHS": int(diff),
        "Delta_paper2": str(delta_p2),
        "Delta_candidate": str(delta_cand),
        "agree": bool(diff == 0),
    })

# Polynomial identity in symbolic m
LHS_sym = sp.expand(lam_abs(m) * N_count(m - 1))
RHS_sym = sp.expand(B_poly_m - N_INIT)
diff_sym = sp.simplify(LHS_sym - RHS_sym)

print("\nLHS(m) =", LHS_sym)
print("RHS(m) =", RHS_sym)
print("LHS - RHS =", diff_sym)

identity_holds_poly = (diff_sym == 0)
print("Polynomial identity holds:", identity_holds_poly)

# If not an identity, test at m=3 only
lhs3 = lam_abs(3) * N_count(2)
rhs3 = B_casimir(3) - N_INIT
coincidence_at_3 = (lhs3 == rhs3)

# ---------------------------------------------------------------------------
# SUBTASK 2: Alternative decompositions of Delta = 1/40
# ---------------------------------------------------------------------------
B = 42
F_val = sp.pi**2 / 6  # for reference, not used in rational tests
kappa = sp.Rational(-1, 16)
d_max = 4
lam3 = lam_abs(3)  # = 8
N2 = N_count(2)    # = 5
Delta_target = sp.Rational(1, 40)

checks = []

def check(label, expr):
    val = sp.simplify(expr)
    agree = (val == Delta_target)
    checks.append({"formula": label, "value": str(val), "equals_1/40": bool(agree)})

check("1/(B - N_init) = 1/40", sp.Rational(1) / (B - N_INIT))
check("kappa^2 * B", kappa**2 * B)
check("|kappa|/B", sp.Abs(kappa) / B)
check("B * |kappa|^2", B * kappa**2)
check("1/(|lambda_3| * N(2)) [Paper 2]", sp.Rational(1) / (lam3 * N2))
check("N_init / (B * N(2))", sp.Rational(N_INIT) / (B * N2))
check("1/(d_max^2 * N(2))", sp.Rational(1) / (d_max**2 * 5))
check("1/(8*5)", sp.Rational(1, 8*5))
check("1/(2*20)", sp.Rational(1, 2*20))
check("1/(4*10)", sp.Rational(1, 4*10))
check("1/(d_max * 2 * N(2))", sp.Rational(1) / (d_max * 2 * 5))
check("1/(|lambda_3| * (B - N_init)/|lambda_3|)", sp.Rational(1) / (lam3 * ((B - N_INIT) // 8)))
check("|kappa| * (|lambda_3|/N(2))", sp.Abs(kappa) * sp.Rational(lam3, N2))
check("2*|kappa| * (|lambda_3|/N(2))", 2 * sp.Abs(kappa) * sp.Rational(lam3, N2))

# Systematic search: 1/40 = a/(b * X) where X in {B, N(2), lambda_3, B-N_init, d_max, N_init}
quantities = {
    "B": B,
    "N(2)": int(N2),
    "|lambda_3|": int(lam3),
    "B-N_init": B - N_INIT,
    "d_max": d_max,
    "N_init": N_INIT,
    "B+N_init": B + N_INIT,
}

systematic = []
for qname, qval in quantities.items():
    for a in range(1, 11):
        for b in range(1, 21):
            if a * qval == 40 * b:  # a/(b*q) = 1/40
                systematic.append({
                    "form": f"{a}/({b}*{qname})",
                    "q_value": qval,
                    "equals_1/40": True,
                })

# Two-quantity products
two_q = []
qlist = list(quantities.items())
for i, (n1, v1) in enumerate(qlist):
    for j, (n2, v2) in enumerate(qlist):
        if v1 * v2 == 40:
            two_q.append(f"1/({n1}*{n2}) = 1/{v1*v2}")
        if v1 * v2 == 20 and 2 == 40 // (v1*v2):
            two_q.append(f"2/({n1}*{n2}) = 1/{v1*v2//2}")

# ---------------------------------------------------------------------------
# Save data
# ---------------------------------------------------------------------------
data = {
    "track": "alpha-F",
    "phase": "4C",
    "goal": "Test whether Delta = 1/(B - N_init) is a structural identity",
    "definitions": {
        "lambda_abs_n": "n^2 - 1",
        "N_m": "m(m+1)(2m+1)/6 (sum of squares 1..m)",
        "B_m": "sum_{n=1}^{m} sum_{l=0}^{n-1} (2l+1)*l*(l+1)",
        "B_m_closed_form": str(sp.expand(B_poly_m)),
        "N_init": N_INIT,
        "Delta_paper2": "1/(|lambda_{n_max}| * N(n_max - 1))",
        "Delta_candidate": "1/(B(n_max) - N_init)",
    },
    "subtask_1": {
        "table": rows,
        "LHS_symbolic": str(LHS_sym),
        "RHS_symbolic": str(RHS_sym),
        "LHS_minus_RHS_symbolic": str(diff_sym),
        "polynomial_identity_holds": bool(identity_holds_poly),
        "coincidence_at_m_3": bool(coincidence_at_3),
    },
    "subtask_2": {
        "target": "1/40",
        "manual_checks": checks,
        "systematic_single_quantity": systematic,
        "two_quantity_products": two_q,
    },
}

OUT_DATA.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_DATA, "w") as f:
    json.dump(data, f, indent=2)

print(f"\nData saved to {OUT_DATA}")

# ---------------------------------------------------------------------------
# Analysis report
# ---------------------------------------------------------------------------
def verdict_text():
    if identity_holds_poly:
        return "IDENTITY (polynomial in m)"
    elif coincidence_at_3:
        failing = [r["m"] for r in rows if not r["agree"]]
        return f"COINCIDENCE AT m=3 ONLY (fails at m = {failing})"
    else:
        return "DISAGREEMENT AT m=3 (contradicts Phase 4B observation)"

md = f"""# Track alpha-F Analysis: Delta = 1/(B - N_init)?

## Setup
Tests whether Paper 2's Delta = 1/(|lambda_{{n_max}}| * N(n_max - 1))
coincides algebraically with 1/(B(n_max) - N_init), where N_init = 2
(Paper 0 Axiom 1), and whether any coincidence is polynomial in n_max.

Closed forms:
- |lambda_n| = n^2 - 1
- N(m) = m(m+1)(2m+1)/6
- B(m) = sum_{{n=1}}^{{m}} sum_{{l=0}}^{{n-1}} (2l+1) l(l+1)
- B(m) closed form: {sp.expand(B_poly_m)}
- B(3) = {int(B_casimir(3))} (verified against Paper 2 Eq. box line 226)

## Subtask 1: Algebraic equivalence

Compare LHS(m) = |lambda_m| * N(m-1) against RHS(m) = B(m) - N_init.

| m | \\|lambda_m\\| | N(m-1) | LHS | B(m) | RHS = B(m) - 2 | LHS - RHS | agree |
|---|-----|--------|-----|------|------|------------|-------|
"""
for r in rows:
    md += f"| {r['m']} | {r['lambda_abs']} | {r['N_m_minus_1']} | {r['LHS']} | {r['B_m']} | {r['RHS']} | {r['LHS_minus_RHS']} | {'YES' if r['agree'] else 'no'} |\n"

md += f"""
### Symbolic polynomial check

- LHS(m) = {LHS_sym}
- RHS(m) = {RHS_sym}
- LHS(m) - RHS(m) = {diff_sym}

**Polynomial identity in m: {"HOLDS" if identity_holds_poly else "DOES NOT HOLD"}**

**Verdict for Subtask 1: {verdict_text()}**

## Subtask 2: Alternative decompositions of Delta = 1/40

### Manual candidates
| Formula | Value | Equals 1/40 |
|---------|-------|-------------|
"""
for c in checks:
    md += f"| {c['formula']} | {c['value']} | {'YES' if c['equals_1/40'] else 'no'} |\n"

md += f"""
### Systematic single-quantity search (a/(b*Q) = 1/40)
Found {len(systematic)} hits (including trivial reorderings):

"""
for s in systematic[:40]:
    md += f"- {s['form']} (Q = {s['q_value']})\n"

md += f"""
### Two-quantity products with value 40
"""
for t in set(two_q):
    md += f"- {t}\n"

# Summary
n_clean_alts = sum(1 for c in checks if c["equals_1/40"])
md += f"""

## Verdict

**Subtask 1:** {'Identity holds as a polynomial in m' if identity_holds_poly else 'NOT a polynomial identity; see m-table for failing rows'}
**Subtask 2:** {n_clean_alts} manual decompositions equal 1/40.
The forms 1/(|lambda_3| * N(2)) and 1/(B - N_init) both equal 1/40, but for distinct reasons unless the m-polynomial identity holds.

**Classification of Delta:**
"""
if identity_holds_poly:
    md += """
Delta is STRUCTURALLY DETERMINED by B and N_init (and by the geometric quantities |lambda| and N via an equivalent polynomial identity). The two Paper 2 inputs
B and Delta are NOT independent -- Delta is a derived quantity.

**Recommendation: FLAG Paper 2 for update.** The three-ingredient form
K = pi(B + F - Delta) reduces to a two-ingredient form in the Casimir/N_init
combination. Paper 2 Section III should be revised to either (a) derive Delta
from B and N_init explicitly, or (b) present Delta in the LHS form as a
reminder that it tracks the highest-shell Laplace-Beltrami weight times the
cumulative state count.
"""
else:
    md += f"""
At m = 3 both quantities equal 40, but {"they disagree elsewhere" if not identity_holds_poly else ""}.
This means the Phase 4B observation Delta = 1/(B - N_init) is a single-point
coincidence, not a structural identity. Delta remains an independent ingredient
in Paper 2's K formula.

Failing rows (m != 3): {[r['m'] for r in rows if not r['agree']]}

**Recommendation: CLOSE this line of investigation.** Delta is not reducible to
B and N_init except at the specific cutoff m = 3. Phase 4B's observation was
real at the data point but does not generalize.
"""

with open(OUT_ANALYSIS, "w") as f:
    f.write(md)

print(f"Analysis saved to {OUT_ANALYSIS}")
print("\n=== VERDICT ===")
print(verdict_text())
