"""S^(3) symbolic closed-form assembly (modulo W10).

Assembles S^(3) = S3_linear + product_terms entirely in sympy Rational
arithmetic, substitutes every identified t-value closed form, reduces
Hurwitz zeta via zeta(s,5/2) = (2^s-1)zeta(s) - 2^s - (2/3)^s,
and verifies the result numerically at 200 dps.

Strategy for unidentified t3 pairs: the linear table contains pairs whose
SUM is identified via stuffle relations (R2, R4, R5, R6) without needing
individual values:
  coeff(2,1,3)*t3(2,1,3) + coeff(2,3,1)*t3(2,3,1)
  = -2048*(t3(2,1,3)+t3(2,3,1))  [R2: sum = lam(3)*t2(2,1) - t2(5,1) - t2(2,4) - t3(3,2,1)]
  coeff(4,1,3)*t3(4,1,3) + coeff(4,3,1)*t3(4,3,1)
  = 2048*(t3(4,1,3)+t3(4,3,1))   [R6: sum = lam(3)*t2(4,1) - t2(7,1) - t2(4,4) - t3(3,4,1)]
  coeff(2,2,2)*t3(2,2,2) = -1024*t3(2,2,2)  [R1: = (lam(2)*t2(2,2)-t2(4,2)-t2(2,4))/3]
  coeff(2,3,3)*t3(2,3,3) = 2048*t3(2,3,3)  [R5]
  coeff(2,4,2)*t3(2,4,2) + coeff(4,2,2)*t3(4,2,2) [R4 + stageB]
  => t3(2,4,2) = lam(2)*t2(4,2) - 2*t3(4,2,2) - t2(6,2) - t2(4,4); t3(4,2,2) identified

Unidentified (W10) set: t3(4,3,3) and t3(4,4,2) — weight-10 depth-3
level-2 objects, not yet identified. Their table coefficients are:
    t3(4,3,3): -2048
    t3(4,4,2): -1024
W10 is evaluated numerically from the 220-dps cache for the gate check.

Decision gate: |symbolic_part + W10_numerical - canonical| <= 1e-180.
Individual substitution verification: each t-value residual <= 1e-50 at 60 dps.

Output: debug/data/s3_symbolic_assembly.json
"""
from __future__ import annotations

import ast
import json
from fractions import Fraction
from pathlib import Path

import mpmath
import sympy as sp

HERE = Path(__file__).parent
OUT: dict = {}

# ---------------------------------------------------------------------------
# 0. Precision setup
# ---------------------------------------------------------------------------
DPS_GATE = 200
DPS_INDIV = 60
mpmath.mp.dps = DPS_GATE

CANONICAL_STR = (
    "31.57256120751202275476192954506898337195758446957719412870053697663883"
    "457068544809094473302375788237499843246296775122195468382644692707154045"
    "096529427150951831230805045804466710135470254516592227073"
)
CANONICAL = mpmath.mpf(CANONICAL_STR)

CACHE = json.loads((HERE / "data" / "s3_pslq_cache.json").read_text())


def cval_mp(key: str) -> mpmath.mpf:
    """Read a cache value (prefer 220-dps over 340-dps)."""
    for tag in [f"{key}@220", f"{key}@340"]:
        if tag in CACHE:
            return mpmath.mpf(CACHE[tag])
    raise KeyError(f"cache miss: {key}")


def lam_mp(s: int) -> mpmath.mpf:
    return (1 - mpmath.power(2, -s)) * mpmath.zeta(s)


# ---------------------------------------------------------------------------
# 1. Sympy symbols (exact rational arithmetic everywhere)
# ---------------------------------------------------------------------------
pi   = sp.Symbol("pi")
ln2  = sp.Symbol("ln2")
z3   = sp.Symbol("z3")
z5   = sp.Symbol("z5")
z7   = sp.Symbol("z7")
z9   = sp.Symbol("z9")
z11  = sp.Symbol("z11")
Li4  = sp.Symbol("Li4")   # Li_4(1/2)
t51  = sp.Symbol("t51")   # t2(5,1) — depth-2 generator (irreducible vs products)
t53  = sp.Symbol("t53")   # t2(5,3) — depth-2 generator
t71  = sp.Symbol("t71")   # t2(7,1) — depth-2 generator

pi2 = pi**2; pi4 = pi**4; pi6 = pi**6; pi8 = pi**8

# ---------------------------------------------------------------------------
# 2. λ(s) = (1-2^{-s})ζ(s) expressed symbolically
# ---------------------------------------------------------------------------
ZETA_EVEN = {
    2:  pi2 / 6,
    4:  pi4 / 90,
    6:  pi6 / 945,
    8:  pi8 / 9450,
    10: sp.Rational(1, 93555) * pi**10,
    12: sp.Rational(691, 638512875) * pi**12,
}
ZETA_ODD = {3: z3, 5: z5, 7: z7, 9: z9, 11: z11}


def lam_sym(s: int) -> sp.Expr:
    z = ZETA_EVEN[s] if s % 2 == 0 else ZETA_ODD[s]
    return (1 - sp.Rational(1, 2**s)) * z


# ---------------------------------------------------------------------------
# 3. Hurwitz reduction: ζ(s,5/2) = (2^s-1)ζ(s) - 2^s - (2/3)^s
# ---------------------------------------------------------------------------
def hz_sym(s: int) -> sp.Expr:
    z = ZETA_EVEN[s] if s % 2 == 0 else ZETA_ODD[s]
    return (2**s - 1) * z - 2**s - sp.Rational(2, 3)**s


# P, Q, R3 from the full relation (after Hurwitz reduction)
P_sym = pi2 - pi4 / 12 - sp.Rational(64, 81)
Q_sym = 4*hz_sym(4) - 2*hz_sym(6) + sp.Rational(1, 4)*hz_sym(8)
R3_sym = 8*hz_sym(6) - 6*hz_sym(8) + sp.Rational(3, 2)*hz_sym(10) - sp.Rational(1, 8)*hz_sym(12)

# ---------------------------------------------------------------------------
# 4. S_min symbolic closed form (Paper 28, verified 1.7e-129)
# ---------------------------------------------------------------------------
S_min_sym = (8*pi2*ln2 - sp.Rational(2,3)*pi4*ln2
             - 3*pi2*z3 + sp.Rational(1,2)*pi4*z3
             - sp.Rational(5,2)*pi2*z5
             + pi6/4 - sp.Rational(3,2)*pi4 - pi8/96)

# ---------------------------------------------------------------------------
# 5. t2 substitution table (all Rational coefficients, no floats)
# ---------------------------------------------------------------------------
T2: dict[tuple, sp.Expr] = {}

# Odd-weight doubles (all reduce to depth <= 1 by parity theorem)
# w3
T2[2,1] = sp.Rational(1,8)*pi2*ln2 - sp.Rational(7,16)*z3
# w5
T2[4,1] = (sp.Rational(1,96)*pi4*ln2 - sp.Rational(1,64)*pi2*z3
           - sp.Rational(31,64)*z5)
T2[2,3] = sp.Rational(1,16)*pi2*z3 - sp.Rational(31,64)*z5
# w7
T2[4,3] = (sp.Rational(1,128)*pi4*z3 - sp.Rational(5,128)*pi2*z5
           - sp.Rational(127,256)*z7)
T2[6,1] = (- sp.Rational(127,256)*z7 - sp.Rational(1,256)*pi2*z5
           - sp.Rational(1,768)*pi4*z3 + sp.Rational(1,960)*pi6*ln2)
T2[2,5] = (- sp.Rational(127,256)*z7 + sp.Rational(13,128)*pi2*z5
           - sp.Rational(1,384)*pi4*z3)
# w9
T2[6,3] = (- sp.Rational(511,1024)*z9 - sp.Rational(21,1024)*pi2*z7
           - sp.Rational(1,512)*pi4*z5 + sp.Rational(1,1280)*pi6*z3)
T2[8,1] = (- sp.Rational(511,1024)*z9 - sp.Rational(1,1024)*pi2*z7
           - sp.Rational(1,3072)*pi4*z5 - sp.Rational(1,7680)*pi6*z3
           + sp.Rational(17,161280)*pi8*ln2)
T2[2,7] = (- sp.Rational(511,1024)*z9 + sp.Rational(15,128)*pi2*z7
           - sp.Rational(1,768)*pi4*z5 - sp.Rational(1,3840)*pi6*z3)
T2[4,5] = (- sp.Rational(511,1024)*z9 - sp.Rational(35,1024)*pi2*z7
           + sp.Rational(13,1536)*pi4*z5)
# w11
T2[8,3] = (- sp.Rational(2047,4096)*z11 - sp.Rational(9,1024)*pi2*z9
           - sp.Rational(5,4096)*pi4*z7 - sp.Rational(1,5120)*pi6*z5
           + sp.Rational(17,215040)*pi8*z3)
T2[4,7] = (- sp.Rational(2047,4096)*z11 - sp.Rational(21,1024)*pi2*z9
           + sp.Rational(53,6144)*pi4*z7 - sp.Rational(1,7680)*pi6*z5)

# w10 even doubles: t2(4,6), t2(7,3), t2(8,2) — unidentified, enter W10 remainder
# These are marked as None and handled in the assembly loop as W10
T2_W10 = {(4,6), (7,3), (8,2)}

# Even-weight doubles via stuffle: t2(b,b) = (lam(b)^2 - lam(2b))/2
T2[2,2] = (lam_sym(2)**2 - lam_sym(4)) / 2
T2[4,4] = (lam_sym(4)**2 - lam_sym(8)) / 2
# antisymmetric a6 = t2(4,2)-t2(2,4) = -(1/1792)pi^6 + (7/32)z3^2
a6_sym = - sp.Rational(1,1792)*pi6 + sp.Rational(7,32)*z3**2
# t2(b,c)+t2(c,b) = lam(b)*lam(c) - lam(b+c)
T2[4,2] = (lam_sym(4)*lam_sym(2) - lam_sym(6) + a6_sym) / 2
T2[2,4] = (lam_sym(2)*lam_sym(4) - lam_sym(6) - a6_sym) / 2
# w8 even: from stageB
T2[6,2] = (sp.Rational(2,13)*t53 + sp.Rational(217,1664)*z3*z5
           - sp.Rational(11,645120)*pi8)
T2[2,6] = (- sp.Rational(2,13)*t53 - sp.Rational(217,1664)*z3*z5
           + sp.Rational(3,71680)*pi8)
# w4
T2[3,1] = (- Li4 + sp.Rational(23,5760)*pi4
           + sp.Rational(1,24)*pi2*ln2**2 - sp.Rational(1,24)*ln2**4)
# w6 diagonal
T2[3,3] = (lam_sym(3)**2 - lam_sym(6)) / 2
# w7 odd — also in G-table: t2(7,3)
T2[7,3] = (- sp.Rational(2047,4096)*z11 - sp.Rational(9,1024)*pi2*z9
           - sp.Rational(5,4096)*pi4*z7 - sp.Rational(1,5120)*pi6*z5
           + sp.Rational(17,215040)*pi8*z3)
# w7 odd — t2(2,7), t2(4,7) already set above; t2(7,1) is a depth-2 generator => t71
# t2(8,1) already set above; t2(8,2) is a w10 object
# Let's also fill depth-2 generators (these appear symbolically as t51, t53, t71)
T2[5,1] = t51
T2[5,3] = t53
T2[7,1] = t71

# ---------------------------------------------------------------------------
# 6. t3 substitution expressions
# ---------------------------------------------------------------------------
# From closure trailing-1 JSON (all verified >= 1e-216 residual):
T3_211 = (sp.Rational(1,2)*Li4 + sp.Rational(1,48)*ln2**4
          + sp.Rational(1,24)*pi2*ln2**2 - sp.Rational(19,5760)*pi4)
T3_411 = (sp.Rational(1,192)*pi4*ln2**2 - sp.Rational(1,64)*pi2*ln2*z3
          - sp.Rational(11,322560)*pi6 - sp.Rational(1,2)*t51
          - sp.Rational(7,128)*z3**2)
T3_321 = (sp.Rational(3,64)*pi2*ln2*z3 - sp.Rational(1,9216)*pi6
          - sp.Rational(1,2)*t51 - sp.Rational(49,256)*z3**2)
T3_341 = (sp.Rational(1,768)*pi4*ln2*z3 + sp.Rational(5,128)*pi2*ln2*z5
          + sp.Rational(61,2903040)*pi8 - sp.Rational(1,256)*pi2*z3**2
          + sp.Rational(1,2)*t53 - sp.Rational(1,2)*t71
          - sp.Rational(217,512)*z3*z5)

# From stageB:
T3_422 = (- sp.Rational(1,13)*t53 - sp.Rational(217,3328)*z3*z5
          + sp.Rational(3,512)*pi2*z3**2 - sp.Rational(1,5806080)*pi8)
T3_323 = (- sp.Rational(217,512)*z3*z5 + sp.Rational(3,128)*pi2*z3**2
          + sp.Rational(1,48384)*pi8)

# t3 values derived via stuffle relations:
# R1: lam(2)*t2(2,2) = 3*t3(2,2,2) + t2(4,2) + t2(2,4)
T3_222 = (lam_sym(2)*T2[2,2] - T2[4,2] - T2[2,4]) / 3

# R4: lam(2)*t2(4,2) = t3(2,4,2) + 2*t3(4,2,2) + t2(6,2) + t2(4,4)
T3_242 = lam_sym(2)*T2[4,2] - 2*T3_422 - T2[6,2] - T2[4,4]

# R5: lam(3)*t2(2,3) = t3(3,2,3) + 2*t3(2,3,3) + t2(5,3) + t2(2,6)
T3_233 = (lam_sym(3)*T2[2,3] - T3_323 - T2[5,3] - T2[2,6]) / 2

# PAIR substitutions (the sum is identified; individual values cancel):
# In linear table: -2048*t3(2,1,3) + -2048*t3(2,3,1) = -2048*(t3(2,1,3)+t3(2,3,1))
# R2: t3(2,1,3)+t3(2,3,1) = lam(3)*t2(2,1) - t2(5,1) - t2(2,4) - t3(3,2,1)
PAIR_213_231_SUM = (lam_sym(3)*T2[2,1] - T2[5,1] - T2[2,4] - T3_321)

# In linear table: +2048*t3(4,1,3) + +2048*t3(4,3,1) = 2048*(t3(4,1,3)+t3(4,3,1))
# R6: t3(4,1,3)+t3(4,3,1) = lam(3)*t2(4,1) - t2(7,1) - t2(4,4) - t3(3,4,1)
PAIR_413_431_SUM = (lam_sym(3)*T2[4,1] - T2[7,1] - T2[4,4] - T3_341)

# ---------------------------------------------------------------------------
# 7. S3_linear substitution (from the S3_linear_table)
# ---------------------------------------------------------------------------
tab = json.loads((HERE / "data" / "s3_decomp_table.json").read_text())
lin_table = {ast.literal_eval(k): sp.Rational(int(v[0]), int(v[1]))
             for k, v in tab["S3_linear_table"].items()}

def sym_key(k: tuple) -> sp.Expr:
    """Return symbolic expression for table key k."""
    if k == ('one',):
        return sp.Integer(1)
    if k[0] == 'lam':
        return lam_sym(k[1])
    if k[0] == 'lam2':
        return lam_sym(k[1]) * lam_sym(k[2])
    if k[0] == 't2':
        if (k[1], k[2]) in T2_W10:
            return None  # W10 object
        return T2[k[1], k[2]]
    # t3 cases handled specially below
    raise ValueError(f"unknown key: {k}")


# Build S3_linear symbolically. Handle t3 pairs specially.
# The t3 pairs that enter as SUMS:
# {t3(2,1,3), t3(2,3,1)} — both with coeff -2048 => -2048*(sum) = -2048*PAIR_213_231_SUM
# {t3(4,1,3), t3(4,3,1)} — both with coeff +2048 => 2048*(sum) = 2048*PAIR_413_431_SUM

# Individual t3 entries (not in pairs):
T3_map = {
    ('t3',2,1,1): T3_211,
    ('t3',4,1,1): T3_411,
    ('t3',2,2,2): T3_222,
    ('t3',2,3,3): T3_233,
    ('t3',2,4,2): T3_242,
    ('t3',4,2,2): T3_422,
    # Pair sums handled below; W10 marked as None
    ('t3',4,3,3): None,   # W10
    ('t3',4,4,2): None,   # W10
}

S3_linear_sym = sp.Integer(0)
t3_pair_213_231_done = False
t3_pair_413_431_done = False

for k, coeff in lin_table.items():
    if k[0] != 't3':
        sk = sym_key(k)
        if sk is None:
            pass  # W10 t2 object, skip in symbolic part
        else:
            S3_linear_sym += coeff * sk
    elif k == ('t3',4,3,3) or k == ('t3',4,4,2):
        pass  # W10: skip in symbolic part, add numerically later
    elif k == ('t3',2,1,3) or k == ('t3',2,3,1):
        if not t3_pair_213_231_done:
            # Both have coeff -2048; add once as -2048*(pair_sum)
            S3_linear_sym += sp.Rational(-2048,1) * PAIR_213_231_SUM
            t3_pair_213_231_done = True
    elif k == ('t3',4,1,3) or k == ('t3',4,3,1):
        if not t3_pair_413_431_done:
            S3_linear_sym += sp.Rational(2048,1) * PAIR_413_431_SUM
            t3_pair_413_431_done = True
    else:
        expr = T3_map.get(k)
        if expr is None:
            raise ValueError(f"unhandled t3 key: {k}")
        S3_linear_sym += coeff * expr

S3_linear_sym = sp.expand(S3_linear_sym)
print("S3_linear symbolic assembled. Collecting...")

# ---------------------------------------------------------------------------
# 8. S_min product terms
# ---------------------------------------------------------------------------
product_terms_sym = (- 16*S_min_sym*P_sym
                     - 8*P_sym**3
                     + 8*P_sym*Q_sym
                     + R3_sym)
product_terms_sym = sp.expand(product_terms_sym)
print("Product terms assembled.")

# ---------------------------------------------------------------------------
# 9. Full S^(3) symbolic (modulo W10)
# ---------------------------------------------------------------------------
S3_sym_modW10 = sp.expand(S3_linear_sym + product_terms_sym)
print(f"S^(3) symbolic (mod W10) assembled. Number of terms: "
      f"{len(S3_sym_modW10.as_ordered_terms())}")

# Collect by symbol
print("\nCollecting by symbol...")
# Display coefficient of each fundamental constant
all_symbols = [pi, ln2, z3, z5, z7, z9, z11, Li4, t51, t53, t71]
collected = sp.collect(S3_sym_modW10, all_symbols, evaluate=False)

# ---------------------------------------------------------------------------
# 10. Verify individual substitutions at 60 dps
# ---------------------------------------------------------------------------
print("\n=== Individual substitution verifications (60 dps) ===")
mpmath.mp.dps = DPS_INDIV

def mp_eval(expr: sp.Expr) -> mpmath.mpf:
    """Evaluate sympy expression numerically."""
    f = sp.lambdify(
        (pi, ln2, z3, z5, z7, z9, z11, Li4, t51, t53, t71),
        expr, "mpmath")
    return f(
        mpmath.pi, mpmath.log(2), mpmath.zeta(3), mpmath.zeta(5),
        mpmath.zeta(7), mpmath.zeta(9), mpmath.zeta(11),
        mpmath.polylog(4, mpmath.mpf(1)/2),
        cval_mp("t2(5,1)"), cval_mp("t2(5,3)"), cval_mp("t2(7,1)"))

INDIV_RESULTS = {}
INDIV_GATE_PASS = True

def check(name: str, expr: sp.Expr, cache_key: str, tol: float = 1e-50):
    global INDIV_GATE_PASS
    val_sym = mp_eval(expr)
    val_cache = cval_mp(cache_key)
    res = float(abs(val_sym - val_cache))
    ok = res <= tol
    INDIV_GATE_PASS = INDIV_GATE_PASS and ok
    status = "PASS" if ok else "FAIL"
    print(f"  {name}: residual {res:.3e} -> {status}")
    INDIV_RESULTS[name] = {"residual": f"{res:.3e}", "pass": ok}
    return ok

# Check all identified t2 values
check("t2(2,1)", T2[2,1], "t2(2,1)", tol=1e-50)
check("t2(4,1)", T2[4,1], "t2(4,1)", tol=1e-50)
check("t2(2,3)", T2[2,3], "t2(2,3)", tol=1e-50)
check("t2(4,3)", T2[4,3], "t2(4,3)", tol=1e-50)
check("t2(6,1)", T2[6,1], "t2(6,1)", tol=1e-50)
check("t2(2,5)", T2[2,5], "t2(2,5)", tol=1e-50)
check("t2(6,3)", T2[6,3], "t2(6,3)", tol=1e-50)
check("t2(8,1)", T2[8,1], "t2(8,1)", tol=1e-50)
check("t2(2,7)", T2[2,7], "t2(2,7)", tol=1e-50)
check("t2(4,5)", T2[4,5], "t2(4,5)", tol=1e-50)
check("t2(8,3)", T2[8,3], "t2(8,3)", tol=1e-50)
check("t2(4,7)", T2[4,7], "t2(4,7)", tol=1e-50)
check("t2(2,2)", T2[2,2], "t2(2,2)", tol=1e-50)
check("t2(4,4)", T2[4,4], "t2(4,4)", tol=1e-50)
check("t2(4,2)", T2[4,2], "t2(4,2)", tol=1e-50)
check("t2(2,4)", T2[2,4], "t2(2,4)", tol=1e-50)
check("t2(6,2)", T2[6,2], "t2(6,2)", tol=1e-50)
check("t2(2,6)", T2[2,6], "t2(2,6)", tol=1e-50)
check("t2(3,1)", T2[3,1], "t2(3,1)", tol=1e-50)
check("t2(3,3)", T2[3,3], "t2(3,3)", tol=1e-50)

# Check identified t3 values
check("t3(2,1,1)", T3_211, "t3(2,1,1)", tol=1e-50)
check("t3(4,1,1)", T3_411, "t3(4,1,1)", tol=1e-50)
check("t3(3,2,1)", T3_321, "t3(3,2,1)", tol=1e-50)
check("t3(3,4,1)", T3_341, "t3(3,4,1)", tol=1e-50)
check("t3(4,2,2)", T3_422, "t3(4,2,2)", tol=1e-50)
check("t3(3,2,3)", T3_323, "t3(3,2,3)", tol=1e-50)
check("t3(2,2,2)", T3_222, "t3(2,2,2)", tol=1e-50)
check("t3(2,4,2)", T3_242, "t3(2,4,2)", tol=1e-50)
check("t3(2,3,3)", T3_233, "t3(2,3,3)", tol=1e-50)

# Check the pair sums
val_pair1 = mp_eval(PAIR_213_231_SUM)
cache_pair1 = cval_mp("t3(2,1,3)") + cval_mp("t3(2,3,1)")
res1 = float(abs(val_pair1 - cache_pair1))
ok1 = res1 <= 1e-50
INDIV_GATE_PASS = INDIV_GATE_PASS and ok1
print(f"  PAIR t3(2,1,3)+t3(2,3,1): residual {res1:.3e} -> {'PASS' if ok1 else 'FAIL'}")
INDIV_RESULTS["PAIR_t3(2,1,3)+t3(2,3,1)"] = {"residual": f"{res1:.3e}", "pass": ok1}

val_pair2 = mp_eval(PAIR_413_431_SUM)
cache_pair2 = cval_mp("t3(4,1,3)") + cval_mp("t3(4,3,1)")
res2 = float(abs(val_pair2 - cache_pair2))
ok2 = res2 <= 1e-50
INDIV_GATE_PASS = INDIV_GATE_PASS and ok2
print(f"  PAIR t3(4,1,3)+t3(4,3,1): residual {res2:.3e} -> {'PASS' if ok2 else 'FAIL'}")
INDIV_RESULTS["PAIR_t3(4,1,3)+t3(4,3,1)"] = {"residual": f"{res2:.3e}", "pass": ok2}

# Check S_min
val_smin = mp_eval(S_min_sym)
cache_smin = mpmath.mpf("2.479936938034222554413579500829382144687925786617288458378798726559552777818374912328558931430046996")
res_smin = float(abs(val_smin - cache_smin))
ok_smin = res_smin <= 1e-50
INDIV_GATE_PASS = INDIV_GATE_PASS and ok_smin
print(f"  S_min: residual {res_smin:.3e} -> {'PASS' if ok_smin else 'FAIL'}")
INDIV_RESULTS["S_min"] = {"residual": f"{res_smin:.3e}", "pass": ok_smin}

print(f"\n  Individual gate: {'ALL PASS' if INDIV_GATE_PASS else 'SOME FAIL'}")

# ---------------------------------------------------------------------------
# 11. Main gate: symbolic_part + W10_numerical vs canonical at 200 dps
# ---------------------------------------------------------------------------
print("\n=== Main gate (200 dps) ===")
mpmath.mp.dps = DPS_GATE

sym_val = mp_eval(S3_sym_modW10)

# W10 numerical values from cache
# t3 W10 objects:
W10_coeff_433 = mpmath.mpf("-2048")
W10_coeff_442 = mpmath.mpf("-1024")
W10_433 = cval_mp("t3(4,3,3)")
W10_442 = cval_mp("t3(4,4,2)")
# t2 W10 objects (with their table coefficients):
W10_t2_46_coeff = mpmath.mpf(lin_table[('t2',4,6)])   # = -1536
W10_t2_73_coeff = mpmath.mpf(lin_table[('t2',7,3)])   # = -1024
W10_t2_82_coeff = mpmath.mpf(lin_table[('t2',8,2)])   # = -512
W10_t2_46 = cval_mp("t2(4,6)")
W10_t2_73 = cval_mp("t2(7,3)")
W10_t2_82 = cval_mp("t2(8,2)")
W10_num = (W10_coeff_433 * W10_433 + W10_coeff_442 * W10_442
           + W10_t2_46_coeff * W10_t2_46
           + W10_t2_73_coeff * W10_t2_73
           + W10_t2_82_coeff * W10_t2_82)
print(f"  W10 numerical contribution: {mpmath.nstr(W10_num, 20)}")

S3_total = sym_val + W10_num
residual = abs(S3_total - CANONICAL)
gate_pass = residual <= mpmath.mpf("1e-180")
print(f"  S3 symbolic + W10_num = {mpmath.nstr(S3_total, 40)}")
print(f"  canonical              = {mpmath.nstr(CANONICAL, 40)}")
print(f"  residual               = {mpmath.nstr(residual, 5)}")
print(f"  GATE: {'PASS' if gate_pass else 'FAIL'}")

# ---------------------------------------------------------------------------
# 12. Analyze which depth-2 generators survive in the full assembly
# ---------------------------------------------------------------------------
print("\n=== Coefficient analysis of depth-2 generators ===")
S3_exp = sp.expand(S3_sym_modW10)

def coeff_of(sym: sp.Symbol) -> sp.Expr:
    """Get the polynomial coefficient of sym in S3_exp (treating sym as a basis element)."""
    return sp.Poly(S3_exp, sym).nth(1)

# Check coefficients of t51, t53, t71
c51 = S3_exp.coeff(t51)
c53 = S3_exp.coeff(t53)
c71 = S3_exp.coeff(t71)
c_Li4 = S3_exp.coeff(Li4)
c_ln2 = S3_exp.coeff(ln2)
c_z7  = S3_exp.coeff(z7)
c_z9  = S3_exp.coeff(z9)
c_z11 = S3_exp.coeff(z11)

print(f"  coefficient of t51 (t2(5,1)): {sp.simplify(c51)}")
print(f"  coefficient of t53 (t2(5,3)): {sp.simplify(c53)}")
print(f"  coefficient of t71 (t2(7,1)): {sp.simplify(c71)}")
print(f"  coefficient of Li4:           {sp.simplify(c_Li4)}")
print(f"  coefficient of z7 (linear):   {sp.simplify(c_z7)}")
print(f"  coefficient of z9 (linear):   {sp.simplify(c_z9)}")
print(f"  coefficient of z11 (linear):  {sp.simplify(c_z11)}")

# Evaluate coefficients numerically to check for cancellation
mpmath.mp.dps = 60
def num_coeff(sym: sp.Symbol) -> str:
    c = S3_exp.coeff(sym)
    if c == 0:
        return "0"
    f = sp.lambdify((pi, ln2, z3, z5, z7, z9, z11, Li4, t51, t53, t71), c, "mpmath")
    v = f(mpmath.pi, mpmath.log(2), mpmath.zeta(3), mpmath.zeta(5),
          mpmath.zeta(7), mpmath.zeta(9), mpmath.zeta(11),
          mpmath.polylog(4, mpmath.mpf(1)/2),
          cval_mp("t2(5,1)"), cval_mp("t2(5,3)"), cval_mp("t2(7,1)"))
    return mpmath.nstr(v, 15)

print("\n  Numerical evaluation of depth-2 coefficients:")
print(f"    coeff(t51) * t51_val = {num_coeff(t51)} * {mpmath.nstr(cval_mp('t2(5,1)'), 10)}")
print(f"    coeff(t53) * t53_val = {num_coeff(t53)} * {mpmath.nstr(cval_mp('t2(5,3)'), 10)}")
print(f"    coeff(t71) * t71_val = {num_coeff(t71)} * {mpmath.nstr(cval_mp('t2(7,1)'), 10)}")

# Check if any depth-2 generator cancels
print(f"\n  t51 survives: {sp.simplify(c51) != 0}")
print(f"  t53 survives: {sp.simplify(c53) != 0}")
print(f"  t71 survives: {sp.simplify(c71) != 0}")

# ---------------------------------------------------------------------------
# 13. Display the closed form compactly
# ---------------------------------------------------------------------------
print("\n=== Closed form S^(3) mod W10 ===")
print("S^(3) = [Symbolic_Part] + W10")
print()
print("W10 = -2048 * t3(4,3,3) - 1024 * t3(4,4,2)  [unidentified weight-10 depth-3 objects]")
print()
# Show the collected form by generator type
print("Generators present in symbolic part:")
for sym, name in [(pi, "pi (pure)"), (ln2, "ln2"), (z3, "z3"), (z5, "z5"),
                  (z7, "z7"), (z9, "z9"), (z11, "z11"), (Li4, "Li4(1/2)"),
                  (t51, "t2(5,1)"), (t53, "t2(5,3)"), (t71, "t2(7,1)")]:
    c = sp.expand(S3_exp.coeff(sym))
    if c != 0:
        print(f"  {name}: coefficient is nonzero")
    else:
        print(f"  {name}: CANCELS (coefficient = 0)")

# Pure constant term
pure = S3_exp.subs([(t51,0),(t53,0),(t71,0),(Li4,0),(ln2,0),
                    (z3,0),(z5,0),(z7,0),(z9,0),(z11,0)])
print(f"\n  Pure pi-polynomial term: {sp.simplify(pure)}")

# ---------------------------------------------------------------------------
# 14. Save output
# ---------------------------------------------------------------------------
mpmath.mp.dps = DPS_GATE

OUT["gate_pass"] = bool(gate_pass)
OUT["gate_residual"] = mpmath.nstr(residual, 5)
OUT["canonical_200dps"] = CANONICAL_STR
OUT["S3_total_200dps"] = mpmath.nstr(S3_total, 60)
OUT["sym_part_200dps"] = mpmath.nstr(sym_val, 60)
OUT["W10_numerical"] = mpmath.nstr(W10_num, 60)
OUT["W10_coefficients"] = {
    "t3(4,3,3)": -2048, "t3(4,4,2)": -1024,
    "t2(4,6)": int(lin_table[('t2',4,6)]),
    "t2(7,3)": int(lin_table[('t2',7,3)]),
    "t2(8,2)": int(lin_table[('t2',8,2)]),
}
OUT["indiv_gate_pass"] = INDIV_GATE_PASS
OUT["indiv_results"] = INDIV_RESULTS

# Depth-2 generator survival
OUT["depth2_generators"] = {
    "t2(5,1)_survives": bool(sp.simplify(c51) != 0),
    "t2(5,3)_survives": bool(sp.simplify(c53) != 0),
    "t2(7,1)_survives": bool(sp.simplify(c71) != 0),
    "Li4_survives": bool(sp.simplify(c_Li4) != 0),
}

# Store symbolic coefficient formulas as strings
OUT["t51_coefficient"] = str(sp.simplify(c51))
OUT["t53_coefficient"] = str(sp.simplify(c53))
OUT["t71_coefficient"] = str(sp.simplify(c71))
OUT["Li4_coefficient"] = str(sp.simplify(c_Li4))

# Full S^(3) mod W10 as sympy string
OUT["S3_symbolic_mod_W10"] = str(S3_sym_modW10)

out_path = HERE / "data" / "s3_symbolic_assembly.json"
out_path.write_text(json.dumps(OUT, indent=2))
print(f"\nSaved: {out_path}")
print(f"\nFinal gate: {'PASS' if gate_pass else 'FAIL'}")
print(f"Individual gates: {'ALL PASS' if INDIV_GATE_PASS else 'SOME FAIL'}")
