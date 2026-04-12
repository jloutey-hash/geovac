"""
Track alpha-J: Arithmetic Origin of F = pi^2/6

Phase 4F alpha sprint. Tests whether F in K = pi*(B + F - Delta) arises
from a Dirichlet series of per-shell Casimir b(n) = n^2(n^2 - 1)/2 or
related lattice weights on the (n, l) shell lattice.

Background: Phases 4B-4E ruled out sphere-spectral mechanisms. The
structural reason is Jacobi inversion: round-sphere Laplace-Beltrami
spectra produce pi-linear (not pi^2) asymptotics and rational
Seeley-DeWitt data. F must come from an ARITHMETIC source.

Key observation:
  D_B(s) = Sum_{n=1}^inf b(n) * n^{-s}
         = (1/2) * [zeta(s-4) - zeta(s-2)]
At s=6: D_B(6) = (1/2)[zeta(2) - zeta(4)] = pi^2/12 - pi^4/180.
The leading term is F/2.

Even simpler: D_{n^2}(s) = Sum_n n^2 * n^{-s} = zeta(s-2).
At s=4: D_{n^2}(4) = zeta(2) = F EXACTLY.
"""

import json
from pathlib import Path
import sympy as sp
from sympy import Rational, zeta, pi, Symbol, simplify, nsimplify, oo, Sum, symbols
import mpmath
mpmath.mp.dps = 50


# ----------------------------------------------------------------------
# Setup: targets and (n,l) shell lattice at n_max = 3
# ----------------------------------------------------------------------

TARGETS_SYM = {
    "F = pi^2/6": pi**2 / 6,
    "2F = pi^2/3": pi**2 / 3,
    "F/2 = pi^2/12": pi**2 / 12,
    "F - Delta": pi**2 / 6 - Rational(1, 40),
    "B + F": 42 + pi**2 / 6,
    "K/pi": 42 + pi**2 / 6 - Rational(1, 40),
    "Delta = 1/40": Rational(1, 40),
    "B = 42": Rational(42),
    "B - N_init = 40": Rational(40),
}

TARGETS_NUM = {
    k: float(mpmath.mpf(str(mpmath.nstr(mpmath.mpf(sp.N(v, 50)), 40))))
    for k, v in TARGETS_SYM.items()
}

# The 6 (n, l) cells at n_max = 3
CELLS = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]

# Per-shell Casimir contribution (Paper 2): b(n) = n^2 * (n^2 - 1) / 2
def b_shell(n):
    return Rational(n * n * (n * n - 1), 2)

# Per-shell basis count: N(n) = n(n+1)(2n+1)/6 (sum of l*(l+1) up to n? -
# actually Paper 2 uses N(m) = m(m+1)(2m+1)/6 = sum of k^2 for k=1..m)
def N_shell(m):
    return Rational(m * (m + 1) * (2 * m + 1), 6)


# ----------------------------------------------------------------------
# SUBTASK 1: Finite shell-lattice Epstein zetas (sanity baseline)
# ----------------------------------------------------------------------

def subtask_1():
    results = {}

    def Q1(n, l):  # n^2 - 1
        return n * n - 1

    def Q2(n, l):  # l(l+1)
        return l * (l + 1)

    def Q3(n, l):  # n^2
        return n * n

    def w_1(n, l):
        return 1

    def w_2l1(n, l):
        return 2 * l + 1

    def w_n2(n, l):
        return n * n

    def w_mixed(n, l):
        return (2 * l + 1) * l * (l + 1)

    quadforms = [("Q1=n^2-1", Q1), ("Q2=l(l+1)", Q2), ("Q3=n^2", Q3)]
    weights = [
        ("w=1", w_1),
        ("w=2l+1", w_2l1),
        ("w=n^2", w_n2),
        ("w=(2l+1)l(l+1)", w_mixed),
    ]

    table = []
    for qname, Q in quadforms:
        for wname, w in weights:
            for s in [1, 2]:
                total = Rational(0)
                for (n, l) in CELLS:
                    qv = Q(n, l)
                    wv = w(n, l)
                    if qv == 0 or wv == 0:
                        continue
                    total += Rational(wv) * Rational(qv) ** (-s)
                entry = {
                    "Q": qname,
                    "weight": wname,
                    "s": s,
                    "symbolic": str(total),
                    "numeric": float(total),
                }
                # Check against rational targets B=42, 40
                for tname, tval in [("B=42", 42), ("B-N_init=40", 40),
                                    ("N_init=2", 2), ("3", 3)]:
                    if total == Rational(tval):
                        entry["HIT"] = tname
                table.append(entry)
    results["epstein_table"] = table

    # Find any hits
    hits = [e for e in table if "HIT" in e]
    results["rational_hits"] = hits
    return results


# ----------------------------------------------------------------------
# SUBTASK 2: Dirichlet series of per-shell Casimir b(n) = n^2(n^2-1)/2
# ----------------------------------------------------------------------

def subtask_2():
    results = {}
    s = Symbol('s')

    # Symbolic: D_B(s) = (1/2)[zeta(s-4) - zeta(s-2)]
    D_B_sym = Rational(1, 2) * (zeta(s - 4) - zeta(s - 2))
    results["D_B_formula"] = "(1/2)*(zeta(s-4) - zeta(s-2))"

    # Tabulate for integer s = 4..12 (s=5 has pole at zeta(1))
    table = []
    for s_val in range(4, 13):
        try:
            val_sym = D_B_sym.subs(s, s_val)
            val_simp = sp.simplify(val_sym)
            val_num = complex(val_sym.evalf(30))
        except Exception as e:
            val_simp = f"error: {e}"
            val_num = None

        entry = {
            "s": s_val,
            "symbolic": str(val_simp),
            "numeric_real": val_num.real if val_num is not None else None,
        }

        # Also tabulate 2*D_B(s), 4*D_B(s)
        entry["2*D_B"] = 2 * val_num.real if val_num is not None else None
        entry["4*D_B"] = 4 * val_num.real if val_num is not None else None

        # Comparison with targets
        if val_num is not None:
            matches = []
            for tname, tval in TARGETS_NUM.items():
                for mult, label in [(1, ""), (2, "2*"), (4, "4*"), (0.5, "(1/2)*")]:
                    v = mult * val_num.real
                    if abs(v - tval) < 1e-10:
                        matches.append(f"EXACT {label}D_B({s_val}) = {tname}")
                    elif abs(v - tval) < 1e-6:
                        matches.append(f"NEAR {label}D_B({s_val}) ~ {tname} (err {abs(v-tval):.2e})")
            entry["matches"] = matches
        table.append(entry)

    # Explicit s=6 check
    val6 = D_B_sym.subs(s, 6)
    val6_simp = sp.simplify(val6)
    results["D_B_at_6_symbolic"] = str(val6_simp)
    results["D_B_at_6_numeric"] = float(val6.evalf(30))
    results["2_D_B_at_6_symbolic"] = str(sp.simplify(2 * val6))
    results["2_D_B_at_6_numeric"] = float((2 * val6).evalf(30))
    results["F_numeric"] = float((pi**2 / 6).evalf(30))
    results["2_D_B_6_minus_F_sym"] = str(sp.simplify(2 * val6 - pi**2 / 6))
    results["2_D_B_6_minus_F_num"] = float((2 * val6 - pi**2 / 6).evalf(30))

    # s=5 pole: regularized
    # zeta(s-4) at s=5 is zeta(1) -> pole
    # Residue is 1/2 (from the zeta(1) residue = 1, times 1/2 prefactor)
    # Regularized (constant term) = (gamma - zeta(3))/2
    gamma = sp.EulerGamma
    D_B_reg_5 = (gamma - zeta(3)) / 2
    results["D_B_reg_at_5_symbolic"] = str(D_B_reg_5)
    results["D_B_reg_at_5_numeric"] = float(D_B_reg_5.evalf(30))
    results["D_B_residue_at_5"] = "1/2 (from zeta(1) pole)"

    results["table"] = table
    return results


# ----------------------------------------------------------------------
# SUBTASK 3: Variant Dirichlet series
# ----------------------------------------------------------------------

def subtask_3():
    results = {}
    s = Symbol('s')

    # (a) D_N(s) = Sum N(n) * n^{-s}
    # N(n) = n(n+1)(2n+1)/6 = (2n^3 + 3n^2 + n)/6
    # D_N(s) = (1/6)[2*zeta(s-3) + 3*zeta(s-2) + zeta(s-1)]
    D_N_sym = Rational(1, 6) * (2 * zeta(s - 3) + 3 * zeta(s - 2) + zeta(s - 1))
    results["D_N_formula"] = "(1/6)*(2*zeta(s-3) + 3*zeta(s-2) + zeta(s-1))"
    table_N = []
    for s_val in [4, 5, 6, 7, 8, 10]:
        try:
            val = D_N_sym.subs(s, s_val)
            try:
                val_num = float(val.evalf(30))
            except Exception:
                val_num = None
            table_N.append({"s": s_val, "symbolic": str(sp.simplify(val)), "numeric": val_num})
        except Exception as e:
            table_N.append({"s": s_val, "error": str(e)})
    results["D_N_table"] = table_N

    # (b) D_{n^2}(s) = Sum n^2 * n^{-s} = zeta(s-2)
    # KEY: at s=4, this equals zeta(2) = pi^2/6 = F EXACTLY
    results["D_n2_formula"] = "Sum_n n^2 * n^{-s} = zeta(s-2)"
    table_n2 = []
    for s_val in [3, 4, 5, 6, 7, 8, 10]:
        val_sym = zeta(s_val - 2)
        try:
            val_simp = sp.simplify(val_sym)
            val_num = float(val_sym.evalf(30))
        except Exception:
            val_simp = "pole"
            val_num = None
        matches = []
        if val_num is not None:
            for tname, tval in TARGETS_NUM.items():
                if abs(val_num - tval) < 1e-10:
                    matches.append(f"EXACT = {tname}")
        table_n2.append({
            "s": s_val,
            "symbolic": str(val_simp),
            "numeric": val_num,
            "matches": matches,
        })
    results["D_n2_table"] = table_n2

    # Emphasize s=4
    results["D_n2_at_4_symbolic"] = str(sp.simplify(zeta(2)))
    results["D_n2_at_4_exact"] = "pi^2/6 = F"

    # Naturalness of s = 4: possible interpretations
    naturalness = {
        "s = 4 = d_max": "Paper 0: d_max = 4 (packing axiom)",
        "s = 4 = dim(S^3) + 1": "dim(S^3) = 3, +1 = 4",
        "s = 4 = 2*N_init": "N_init = 2 (Paper 0), 2*2 = 4",
        "s = 4 = dim(R^4) = ambient": "S^3 embeds in R^4",
        "s = 4 = (n_max = 3) + 1": "one past shell cutoff",
        "s = 4 = Casimir power": "(s-2) = 2 = eigenvalue index of zeta(2)",
    }
    results["s_equals_4_interpretations"] = naturalness

    # (d) Sum |lambda_n| * n^{-s} = Sum (n^2-1) * n^{-s} = zeta(s-2) - zeta(s)
    D_lam_sym = zeta(s - 2) - zeta(s)
    results["D_lambda_formula"] = "Sum (n^2-1) * n^{-s} = zeta(s-2) - zeta(s)"
    table_lam = []
    for s_val in [3, 4, 5, 6, 7, 8]:
        val = D_lam_sym.subs(s, s_val)
        try:
            val_simp = sp.simplify(val)
            val_num = float(val.evalf(30))
        except Exception:
            val_simp = "pole"
            val_num = None
        table_lam.append({
            "s": s_val,
            "symbolic": str(val_simp),
            "numeric": val_num,
        })
    results["D_lambda_table"] = table_lam

    # (f) Sum l_max(n) * n^{-s} = Sum (n-1) * n^{-s} = zeta(s-1) - zeta(s)
    D_lmax_sym = zeta(s - 1) - zeta(s)
    results["D_lmax_formula"] = "Sum (n-1) * n^{-s} = zeta(s-1) - zeta(s)"
    table_lmax = []
    for s_val in [3, 4, 5, 6, 7]:
        val = D_lmax_sym.subs(s, s_val)
        try:
            val_simp = sp.simplify(val)
            val_num = float(val.evalf(30))
        except Exception:
            val_simp = "pole"
            val_num = None
        table_lmax.append({
            "s": s_val,
            "symbolic": str(val_simp),
            "numeric": val_num,
        })
    results["D_lmax_table"] = table_lmax

    # Search all variant Dirichlet series at s = 3..10 for target hits
    all_hits = []
    for name, formula in [
        ("D_B", Rational(1, 2) * (zeta(s - 4) - zeta(s - 2))),
        ("D_N", D_N_sym),
        ("D_n2", zeta(s - 2)),
        ("D_lambda", zeta(s - 2) - zeta(s)),
        ("D_lmax", zeta(s - 1) - zeta(s)),
    ]:
        for s_val in range(3, 13):
            try:
                val = formula.subs(s, s_val)
                try:
                    val_num = float(val.evalf(30))
                except Exception:
                    continue
                if not (val_num == val_num):  # nan
                    continue
                if abs(val_num) > 1e15:
                    continue
                for tname, tval in TARGETS_NUM.items():
                    for mult, mlabel in [(1, ""), (2, "2*"), (Rational(1, 2), "(1/2)*")]:
                        v = float(mult) * val_num
                        if abs(v - tval) < 1e-10:
                            all_hits.append({
                                "series": name,
                                "s": s_val,
                                "mult": str(mlabel),
                                "target": tname,
                                "symbolic": str(sp.simplify(mult * val)),
                            })
            except Exception:
                pass
    results["all_exact_hits"] = all_hits
    return results


# ----------------------------------------------------------------------
# SUBTASK 4: Selection principle in Dirichlet language
# ----------------------------------------------------------------------

def subtask_4():
    results = {}
    s = Symbol('s')

    # D_B(s)/D_N(s)
    # For the selection principle test, use the SIMPLER version:
    # the Paper 2 selection principle is B/N = 3 at m = 3 (a TRUNCATED
    # sum ratio). Test whether a Dirichlet analog exists.
    #
    # Interpretation 1: D_B(s) / D_{n^2}(s) = D_B(s) / zeta(s-2)
    #   = (1/2)*(zeta(s-4)/zeta(s-2) - 1)
    #   Set = 3 => zeta(s-4)/zeta(s-2) = 7
    #
    # Interpretation 2: D_B(s)/D_N(s) where D_N is the full Dirichlet
    # series of N(n). Also compute.

    # Interpretation 1: ratio vs n^2 weight
    ratio_1 = Rational(1, 2) * (zeta(s - 4) / zeta(s - 2) - 1)
    results["ratio_1_formula"] = "D_B(s)/D_{n^2}(s) = (1/2)*(zeta(s-4)/zeta(s-2) - 1)"

    def ratio_1_num(sv):
        return 0.5 * (float(mpmath.zeta(sv - 4)) / float(mpmath.zeta(sv - 2)) - 1)

    table_ratio = []
    for s_val in [4, 5, 6, 7, 8, 10, 12]:
        try:
            r = ratio_1_num(s_val)
            table_ratio.append({"s": s_val, "ratio": r, "eq_to_3": abs(r - 3) < 1e-6})
        except Exception as e:
            table_ratio.append({"s": s_val, "error": str(e)})
    results["ratio_1_table"] = table_ratio

    # Root-find: ratio_1 = 3 => zeta(s-4)/zeta(s-2) = 7
    # Using mpmath at 50 dps
    mpmath.mp.dps = 50
    def f(sv):
        sv_m = mpmath.mpf(sv)
        return mpmath.zeta(sv_m - 4) / mpmath.zeta(sv_m - 2) - 7

    # Scan first
    scan = []
    for sv in [4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 9, 10, 15, 20, 30, 50, 100]:
        try:
            val = float(f(sv))
            scan.append({"s": sv, "f(s) = zeta(s-4)/zeta(s-2) - 7": val})
        except Exception as e:
            scan.append({"s": sv, "error": str(e)})
    results["ratio_1_root_scan"] = scan

    # Try to root-find
    root_1 = None
    try:
        # Find sign changes
        test_points = [4.5, 5, 5.5, 6, 6.5, 7, 8, 10, 15, 20, 50, 100]
        sign_changes = []
        prev_val = None
        prev_s = None
        for sv in test_points:
            try:
                val = float(f(sv))
                if prev_val is not None and prev_val * val < 0:
                    sign_changes.append((prev_s, sv))
                prev_val = val
                prev_s = sv
            except Exception:
                prev_val = None
                prev_s = None
        results["ratio_1_sign_changes"] = sign_changes

        if sign_changes:
            a, b = sign_changes[0]
            root_1 = mpmath.findroot(f, (mpmath.mpf(a) + mpmath.mpf(b)) / 2)
            results["ratio_1_root_s_star"] = str(root_1)
            results["ratio_1_root_s_star_float"] = float(root_1)
        else:
            results["ratio_1_root_s_star"] = "no sign change in scan"
    except Exception as e:
        results["ratio_1_root_error"] = str(e)

    # Interpretation 2: D_B(s)/D_N(s)
    D_B_sym = Rational(1, 2) * (zeta(s - 4) - zeta(s - 2))
    D_N_sym = Rational(1, 6) * (2 * zeta(s - 3) + 3 * zeta(s - 2) + zeta(s - 1))

    def ratio_2_num(sv):
        svm = mpmath.mpf(sv)
        num = 0.5 * (float(mpmath.zeta(svm - 4)) - float(mpmath.zeta(svm - 2)))
        den = (2 * float(mpmath.zeta(svm - 3)) + 3 * float(mpmath.zeta(svm - 2)) + float(mpmath.zeta(svm - 1))) / 6
        return num / den

    table_r2 = []
    for s_val in [4.5, 5, 5.5, 6, 6.5, 7, 8, 10, 12]:
        try:
            r = ratio_2_num(s_val)
            table_r2.append({"s": s_val, "ratio": r, "eq_to_3": abs(r - 3) < 1e-6})
        except Exception as e:
            table_r2.append({"s": s_val, "error": str(e)})
    results["ratio_2_table"] = table_r2

    # Root find ratio_2 = 3
    def f2(sv):
        return ratio_2_num(float(sv)) - 3

    try:
        prev_val = None
        prev_s = None
        for sv in [4.5, 5, 5.5, 6, 6.5, 7, 8, 10, 12, 15, 20, 50]:
            try:
                val = f2(sv)
                if prev_val is not None and prev_val * val < 0:
                    root_2 = mpmath.findroot(lambda x: ratio_2_num(float(x)) - 3, (mpmath.mpf(prev_s) + mpmath.mpf(sv)) / 2)
                    results["ratio_2_root_s_star"] = str(root_2)
                    results["ratio_2_root_s_star_float"] = float(root_2)
                    break
                prev_val = val
                prev_s = sv
            except Exception:
                prev_val = None
                prev_s = None
    except Exception as e:
        results["ratio_2_root_error"] = str(e)

    return results


# ----------------------------------------------------------------------
# SUBTASK 5: Truncation correction at n_max = 3
# ----------------------------------------------------------------------

def subtask_5():
    results = {}

    # For D_{n^2}(s) = zeta(s-2), at s = 4: Sum n^{-2} = zeta(2) = F
    # Truncated: sum_{n=1}^3 n^{-2} = 1 + 1/4 + 1/9 = 49/36
    trunc_n2_4 = Rational(1) + Rational(1, 4) + Rational(1, 9)
    results["D_n2_4_trunc_nmax3_sym"] = str(trunc_n2_4)
    results["D_n2_4_trunc_nmax3_num"] = float(trunc_n2_4)
    tail_n2 = pi**2 / 6 - trunc_n2_4
    tail_n2_sym = sp.simplify(tail_n2)
    tail_n2_num = float(tail_n2.evalf(30))
    results["D_n2_4_tail_sym"] = str(tail_n2_sym)
    results["D_n2_4_tail_num"] = tail_n2_num
    results["D_n2_4_tail_vs_Delta"] = {
        "tail": tail_n2_num,
        "Delta": 1/40,
        "ratio_tail_over_Delta": tail_n2_num * 40,
        "Delta_over_tail": (1/40) / tail_n2_num,
    }

    # D_B(6) = (1/2)[zeta(2) - zeta(4)]
    # Truncated: sum_{n=1}^3 b(n) * n^{-6}
    # b(1)=0, b(2)=2^2*(2^2-1)/2 = 4*3/2 = 6, b(3)=9*8/2 = 36
    trunc_B_6 = Rational(0) + Rational(6, 2**6) + Rational(36, 3**6)
    trunc_B_6_simp = sp.simplify(trunc_B_6)
    results["D_B_6_trunc_nmax3_sym"] = str(trunc_B_6_simp)
    results["D_B_6_trunc_nmax3_num"] = float(trunc_B_6)
    D_B_6_full = (pi**2 / 12) - (pi**4 / 180)
    tail_B = D_B_6_full - trunc_B_6
    tail_B_num = float(tail_B.evalf(30))
    results["D_B_6_tail_sym"] = str(sp.simplify(tail_B))
    results["D_B_6_tail_num"] = tail_B_num
    results["D_B_6_tail_vs_Delta"] = {
        "tail": tail_B_num,
        "Delta": 1/40,
        "ratio_tail_over_Delta": tail_B_num * 40,
        "Delta_over_tail": (1/40) / tail_B_num,
    }

    # Also: does truncated D_n2_4 match B/something?
    # 49/36 ~ 1.361. 42/something?  42/something=1.361 => something = 30.857
    # 40/something=1.361 => something = 29.39
    results["misc_ratios"] = {
        "49/36_times_pi^2": float((Rational(49, 36) * pi**2).evalf(30)),
        "F_minus_49/36": tail_n2_num,
        "F_div_49/36": float((pi**2 / 6 / Rational(49, 36)).evalf(30)),
    }

    # Test whether (F - trunc) * 42 = 1, or similar
    # tail ~ 0.2838. 0.2838 * 42 = 11.92. Nope.
    # But: F - trunc = sum_{n=4}^inf n^{-2}
    # = zeta(2) - (1 + 1/4 + 1/9) = zeta(2) - 49/36
    # Does this equal anything clean? Let's see: check various formulas
    clean_checks = {
        "tail_times_1": tail_n2_num,
        "tail_times_2": 2 * tail_n2_num,
        "tail_times_42": 42 * tail_n2_num,
        "tail_times_40": 40 * tail_n2_num,
        "Delta/tail": (1/40) / tail_n2_num,
        "1/tail": 1 / tail_n2_num,
    }
    results["D_n2_4_tail_clean_checks"] = clean_checks

    # What about a DIFFERENT cutoff? Try n_max = 2 (only n=1,2):
    # trunc = 1 + 1/4 = 5/4. Tail = pi^2/6 - 5/4 ~ 0.3949. Not clean either.
    trunc_n2_4_nmax2 = Rational(1) + Rational(1, 4)
    tail_nmax2 = pi**2/6 - trunc_n2_4_nmax2
    results["D_n2_4_trunc_nmax2_sym"] = str(trunc_n2_4_nmax2)
    results["D_n2_4_tail_nmax2_num"] = float(tail_nmax2.evalf(30))

    return results


# ----------------------------------------------------------------------
# Run all subtasks
# ----------------------------------------------------------------------

def main():
    print("Track alpha-J: Arithmetic Origin of F")
    print("=" * 60)

    all_results = {
        "meta": {
            "phase": "4F",
            "track": "alpha-J",
            "date": "2026-04-10",
            "targets_symbolic": {k: str(v) for k, v in TARGETS_SYM.items()},
            "targets_numeric": TARGETS_NUM,
        }
    }

    print("\nSubtask 1: Finite Epstein zetas...")
    r1 = subtask_1()
    all_results["subtask_1"] = r1
    print(f"  rational hits: {len(r1['rational_hits'])}")

    print("\nSubtask 2: Dirichlet series D_B(s)...")
    r2 = subtask_2()
    all_results["subtask_2"] = r2
    print(f"  D_B(6) = {r2['D_B_at_6_symbolic']}")
    print(f"         = {r2['D_B_at_6_numeric']}")
    print(f"  2*D_B(6) - F = {r2['2_D_B_6_minus_F_sym']}")

    print("\nSubtask 3: Variant Dirichlet series...")
    r3 = subtask_3()
    all_results["subtask_3"] = r3
    print(f"  D_n2(4) exact hit: {r3['D_n2_at_4_exact']}")
    print(f"  Total exact hits: {len(r3['all_exact_hits'])}")
    for h in r3["all_exact_hits"]:
        print(f"    {h['series']}(s={h['s']}) mult={h['mult']} = {h['target']}")

    print("\nSubtask 4: Selection principle root...")
    r4 = subtask_4()
    all_results["subtask_4"] = r4
    if "ratio_1_root_s_star_float" in r4:
        print(f"  Ratio 1 root s* = {r4['ratio_1_root_s_star_float']}")
    else:
        print(f"  Ratio 1: {r4.get('ratio_1_root_s_star', 'no root')}")
    if "ratio_2_root_s_star_float" in r4:
        print(f"  Ratio 2 root s* = {r4['ratio_2_root_s_star_float']}")

    print("\nSubtask 5: Truncation correction at n_max=3...")
    r5 = subtask_5()
    all_results["subtask_5"] = r5
    print(f"  D_n2(4) tail = {r5['D_n2_4_tail_num']}")
    print(f"  Delta = 1/40 = 0.025")
    print(f"  tail / Delta = {r5['D_n2_4_tail_vs_Delta']['ratio_tail_over_Delta']}")
    print(f"  D_B(6) tail = {r5['D_B_6_tail_num']}")

    # Save
    outdir = Path("c:/Users/jlout/OneDrive/Desktop/Project_Geometric/debug/data/track_alpha_phase4f")
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "track_j_arithmetic.json"

    # Convert any non-JSON types
    def serialize(obj):
        if isinstance(obj, dict):
            return {k: serialize(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [serialize(x) for x in obj]
        elif isinstance(obj, tuple):
            return [serialize(x) for x in obj]
        elif isinstance(obj, (int, float, str, bool)) or obj is None:
            return obj
        elif isinstance(obj, complex):
            return {"real": obj.real, "imag": obj.imag}
        else:
            return str(obj)

    with open(outpath, "w") as f:
        json.dump(serialize(all_results), f, indent=2)

    print(f"\nSaved to {outpath}")
    return all_results


if __name__ == "__main__":
    main()
