"""Door 5: the zeta_{D^2} / chi_{-4} integer-s ladder forward-run.

Forcing-catalogue forward-run (docs/forcing_catalogue.md Door 5 / F-7).

Sweeps the integer-s ladder of every Dirac-spectral Dirichlet object in the
GeoVac QED apparatus and asks, for each value:

    (a) confirm the exact closed form (skeleton -> exact arithmetic;
        quarter-integer Hurwitz is exact, PSLQ used ONLY as Layer-2
        cross-check / identification);
    (b) does the value map to a known QED / spectral observable
        (anomalous moment, Lamb-shift term, Sommerfeld sum D_5/D_6,
        vacuum-polarization beta coefficient)?
    (c) is there a value the framework DETERMINES that has NOT yet been
        connected to physics?

Objects swept (all on the unit S^3 Camporesi-Higuchi Dirac spectrum
|lambda_n| = n + 3/2, g_n = 2(n+1)(n+2), n >= 0):

  Z2(s)      = zeta_{D^2}(s) = sum_n g_n / |lambda_n|^{2s}   (T9: pi^even only)
  D(s)       = sum_n g_n / |lambda_n|^s   (Dirac Dirichlet series, s-power)
  Deven(s)   = even-n part of D(s)
  Dodd(s)    = odd-n part of D(s)
  Ddiff(s)   = Deven(s) - Dodd(s) = 2^{s-1}(beta(s) - beta(s-2))   (chi_{-4})

Convergence: D(s) needs s > 3 (g_n ~ 2n^2, |lambda_n|^s ~ n^s);
Z2(s) = D(2s) needs 2s > 3, i.e. s >= 2.

All values are skeleton (Layer 1) quantities: forced, exact. Every
transcendental that appears in a closed form is tagged with its Paper 18
tier / master-Mellin sub-engine (M1 Hopf-base measure / M2 Seeley-DeWitt /
M3 vertex-parity Hurwitz) and its Paper 34 projection.

Run: python debug/door5_zeta_d2_ladder.py
"""

from __future__ import annotations

import json
from fractions import Fraction
from typing import Dict, List, Optional

import mpmath
import sympy as sp

mpmath.mp.dps = 120


# ---------------------------------------------------------------------------
# Exact (sympy) closed forms for the Dirac spectral objects
# ---------------------------------------------------------------------------

def Z2_exact(s: int) -> sp.Expr:
    """zeta_{D^2}(s) exact via T9 theorem.

    Z2(s) = sum_n g_n / |lambda_n|^{2s}
          = 2^{2s-1} [ lambda_Dir(2s-2) - lambda_Dir(2s) ]
    with lambda_Dir(2k) = (1 - 2^{-2k}) zeta(2k).

    For s >= 2 this is a two-term polynomial in pi^2 with rational
    coefficients (T9). Returns a closed-form sympy expression.
    """
    def lam(two_k: int) -> sp.Expr:
        return (1 - sp.Rational(1, 2 ** two_k)) * sp.zeta(two_k)

    # T9: two-term polynomial in pi^2.  sp.zeta(even) auto-evaluates to
    # rational*pi^even, so sp.expand gives the clean two-term pi-form.
    val = sp.Integer(2) ** (2 * s - 1) * (lam(2 * s - 2) - lam(2 * s))
    return sp.expand(val)


def D_exact(s: int) -> sp.Expr:
    """Dirac Dirichlet series D(s) = sum_n g_n / |lambda_n|^s, s-power.

    g_n = 2(n+1)(n+2), |lambda_n| = n + 3/2.  Writing u = n + 3/2 (so
    u runs over 3/2, 5/2, ...), g_n = 2u^2 - 1/2.  Hence

        D(s) = sum_{u} (2u^2 - 1/2) u^{-s}
             = 2 sum u^{2-s} - 1/2 sum u^{-s}
             = 2 zeta(s-2, 3/2) - 1/2 zeta(s, 3/2).

    Hurwitz at half-integer 3/2: zeta(s, 3/2) = (2^s - 1) zeta(s) - 2^s.
    The constant -2^s pieces cancel between the two terms, leaving

        D(s) = 2(2^{s-2} - 1) zeta(s-2) - (1/2)(2^s - 1) zeta(s).

    Needs s > 3 for convergence (here s >= 4 for clean even-shift even-zeta;
    s=4,5,... are the physical ladder).  Returns sympy closed form.
    """
    z_sm2 = sp.zeta(s - 2)
    z_s = sp.zeta(s)
    val = 2 * (sp.Integer(2) ** (s - 2) - 1) * z_sm2 \
        - sp.Rational(1, 2) * (sp.Integer(2) ** s - 1) * z_s
    return sp.nsimplify(sp.simplify(val), [sp.pi])


def Deven_exact(s: int) -> sp.Expr:
    """Even-n part of D(s) via quarter-integer Hurwitz at 3/4.

    Even n = 2k: lambda_{2k} = 2k + 3/2 = 2(k + 3/4) = 2u, u = k + 3/4.
    g_{2k} = 2(2k+1)(2k+2) = 8u^2 - 1/2 (substituting k = u - 3/4).

        Deven(s) = sum_k (8u^2 - 1/2) (2u)^{-s}
                 = 2^{-s} [ 8 zeta(s-2, 3/4) - 1/2 zeta(s, 3/4) ].
    """
    a = sp.Rational(3, 4)
    val = sp.Integer(2) ** (-s) * (
        8 * sp.zeta(s - 2, a) - sp.Rational(1, 2) * sp.zeta(s, a)
    )
    return val


def Dodd_exact(s: int) -> sp.Expr:
    """Odd-n part of D(s) via quarter-integer Hurwitz at 5/4.

    Odd n = 2k+1: lambda_{2k+1} = 2k + 5/2 = 2(k + 5/4) = 2u, u = k + 5/4.
    g_{2k+1} = 2(2k+2)(2k+3) = 8u^2 - 1/2.

        Dodd(s) = 2^{-s} [ 8 zeta(s-2, 5/4) - 1/2 zeta(s, 5/4) ].
    """
    a = sp.Rational(5, 4)
    val = sp.Integer(2) ** (-s) * (
        8 * sp.zeta(s - 2, a) - sp.Rational(1, 2) * sp.zeta(s, a)
    )
    return val


def Ddiff_chi4_exact(s: int) -> sp.Expr:
    """chi_{-4} closed form: Deven(s) - Dodd(s) = 2^{s-1}(beta(s) - beta(s-2)).

    beta = Dirichlet beta (L(.,chi_{-4})).  Returns the RHS closed form
    (the proven Paper 28 Theorem 'chi_4 identity').
    """
    return sp.Integer(2) ** (s - 1) * (
        _sympy_beta(s) - _sympy_beta(s - 2)
    )


def _sympy_beta(s: int) -> sp.Expr:
    """Dirichlet beta(s) as a sympy expression (closed form where known)."""
    # Known closed forms
    if s == 1:
        return sp.pi / 4
    if s == 3:
        return sp.pi ** 3 / 32
    if s == 5:
        return 5 * sp.pi ** 5 / 1536
    if s == 7:
        return 61 * sp.pi ** 7 / 184320
    if s <= 0:
        # beta(0) = 1/2, beta(-1) = 1/2, beta(-2k) = E_{2k}/2, beta(1-2k) via Euler
        # Use Hurwitz form for negative/zero generic
        return (sp.zeta(s, sp.Rational(1, 4)) - sp.zeta(s, sp.Rational(3, 4))) \
            / sp.Integer(4) ** s
    # Even s >= 2 (Catalan G = beta(2), beta(4), beta(6) -- no elementary closed form)
    return sp.Function('beta')(s)


# ---------------------------------------------------------------------------
# Numerical evaluation (mpmath, exact-precision) for cross-check + PSLQ
# ---------------------------------------------------------------------------

def _g_n(n: int) -> mpmath.mpf:
    return mpmath.mpf(2) * (n + 1) * (n + 2)


def _lam(n: int) -> mpmath.mpf:
    return mpmath.mpf(n) + mpmath.mpf(3) / 2


def Z2_num(s: int) -> mpmath.mpf:
    """zeta_{D^2}(s) by direct Hurwitz (exact-precision)."""
    # D(2s) with the s-power formula evaluated at 2s
    return D_num(2 * s)


def D_num(s: int) -> mpmath.mpf:
    """D(s) by Hurwitz at 3/2 (exact-precision)."""
    return (2 * mpmath.hurwitz(s - 2, mpmath.mpf(3) / 2)
            - mpmath.mpf(1) / 2 * mpmath.hurwitz(s, mpmath.mpf(3) / 2))


def Deven_num(s: int) -> mpmath.mpf:
    return mpmath.power(2, -s) * (
        8 * mpmath.hurwitz(s - 2, mpmath.mpf(3) / 4)
        - mpmath.mpf(1) / 2 * mpmath.hurwitz(s, mpmath.mpf(3) / 4)
    )


def Dodd_num(s: int) -> mpmath.mpf:
    return mpmath.power(2, -s) * (
        8 * mpmath.hurwitz(s - 2, mpmath.mpf(5) / 4)
        - mpmath.mpf(1) / 2 * mpmath.hurwitz(s, mpmath.mpf(5) / 4)
    )


def beta_num(s: int) -> mpmath.mpf:
    return (mpmath.hurwitz(s, mpmath.mpf(1) / 4)
            - mpmath.hurwitz(s, mpmath.mpf(3) / 4)) / mpmath.power(4, s)


# ---------------------------------------------------------------------------
# PSLQ identification (Layer-2 cross-check only)
# ---------------------------------------------------------------------------

def pslq_identify(value: mpmath.mpf, basis_terms: List, basis_labels: List[str],
                  tol: float = 1e-90, maxcoeff: int = 10 ** 8) -> Optional[Dict]:
    """Identify value as a Q-linear combo of basis_terms. Returns components."""
    basis = [value] + list(basis_terms)
    try:
        rel = mpmath.pslq(basis, tol=tol, maxcoeff=maxcoeff, maxsteps=10 ** 5)
    except Exception:
        rel = None
    if rel is None or rel[0] == 0:
        return None
    comps = {}
    recon = mpmath.mpf(0)
    for i in range(1, len(rel)):
        if rel[i] != 0:
            c = mpmath.mpf(-rel[i]) / mpmath.mpf(rel[0])
            comps[basis_labels[i - 1]] = str(Fraction(-rel[i], rel[0]))
            recon += c * basis_terms[i - 1]
    return {
        "components": comps,
        "residual": float(abs(value - recon)),
    }


# ---------------------------------------------------------------------------
# Sommerfeld D_p (hydrogen) -- already published through D_6.  Re-state the
# closed forms for the ladder table (these are the "physical-energy" Dirichlet
# sums; D(s) above is the "free spectrum" one).
# ---------------------------------------------------------------------------

SOMMERFELD_DP = {
    2: "-5/4 zeta(2) + zeta(3)",
    3: "19/8 zeta(4) - 11/4 zeta(5)",
    4: "-205/64 zeta(6) + 71/8 zeta(7) - 9/2 zeta(3)zeta(4)",
    5: "497/128 zeta(8) - 467/16 zeta(9) + 385/32 zeta(3)zeta(6) + 75/8 zeta(4)zeta(5)",
    6: ("-2289/512 zeta(10) + 1589/16 zeta(11) - 1617/64 zeta(3)zeta(8) "
        "- 1785/64 zeta(4)zeta(7) - 2065/64 zeta(5)zeta(6)"),
}


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def run_sweep(s_lo: int = 2, s_hi: int = 12) -> Dict:
    pi2 = mpmath.pi ** 2
    pi4 = mpmath.pi ** 4
    pi6 = mpmath.pi ** 6
    pi8 = mpmath.pi ** 8
    pi10 = mpmath.pi ** 10
    pi12 = mpmath.pi ** 12
    pi14 = mpmath.pi ** 14
    pi16 = mpmath.pi ** 16
    pi18 = mpmath.pi ** 18
    pi20 = mpmath.pi ** 20
    pi22 = mpmath.pi ** 22
    pi24 = mpmath.pi ** 24
    pi_even = [pi2, pi4, pi6, pi8, pi10, pi12, pi14, pi16, pi18, pi20, pi22, pi24]
    pi_even_lab = [f"pi^{2*k}" for k in range(1, 13)]

    G = mpmath.catalan

    rows: List[Dict] = []

    # ---------------- Z2(s) = zeta_{D^2}(s) ladder ----------------
    z2_rows = []
    for s in range(2, s_hi + 1):
        exact = Z2_exact(s)
        num = Z2_num(s)
        # cross-check exact vs num
        exact_num = mpmath.mpf(str(sp.N(exact, 100)))
        relerr = float(abs(exact_num - num) / abs(num)) if num != 0 else 0.0
        # PSLQ against pi^even basis (should land in span)
        ident = pslq_identify(num, pi_even, pi_even_lab, tol=1e-80)
        z2_rows.append({
            "s": s,
            "closed_form": str(exact),
            "value": mpmath.nstr(num, 30),
            "exact_vs_num_relerr": relerr,
            "pslq_pi_even": ident,
            "T9_pi_even_only": ident is not None,
        })

    # ---------------- D(s) ladder (free Dirac Dirichlet) ----------------
    d_rows = []
    for s in range(4, s_hi + 1):
        exact = D_exact(s)
        num = D_num(s)
        exact_num = mpmath.mpf(str(sp.N(exact, 100)))
        relerr = float(abs(exact_num - num) / abs(num)) if num != 0 else 0.0
        d_rows.append({
            "s": s,
            "closed_form": str(exact),
            "value": mpmath.nstr(num, 30),
            "exact_vs_num_relerr": relerr,
        })

    # ---------------- Deven / Dodd / Ddiff ladder ----------------
    split_rows = []
    for s in range(4, s_hi + 1):
        de = Deven_num(s)
        do = Dodd_num(s)
        ddiff = de - do
        # chi_4 RHS
        chi4_rhs = mpmath.power(2, s - 1) * (beta_num(s) - beta_num(s - 2))
        chi4_resid = float(abs(ddiff - chi4_rhs))
        # PSLQ Deven against {pi^even, beta(s), beta(s-2)}
        beta_s = beta_num(s)
        beta_sm2 = beta_num(s - 2)
        split_basis = pi_even[: (s)] + [beta_s, beta_sm2]
        split_lab = pi_even_lab[: (s)] + [f"beta({s})", f"beta({s-2})"]
        ident_even = pslq_identify(de, split_basis, split_lab, tol=1e-70)
        split_rows.append({
            "s": s,
            "Deven": mpmath.nstr(de, 30),
            "Dodd": mpmath.nstr(do, 30),
            "Ddiff": mpmath.nstr(ddiff, 30),
            "chi4_closed_form": str(Ddiff_chi4_exact(s)),
            "chi4_residual": chi4_resid,
            "chi4_verified": chi4_resid < 1e-60,
            "pslq_Deven": ident_even,
        })

    return {
        "Z2_ladder": z2_rows,
        "D_ladder": d_rows,
        "split_ladder": split_rows,
        "sommerfeld_Dp": SOMMERFELD_DP,
    }


# ---------------------------------------------------------------------------
# Observable map: which ladder values are connected vs unconnected
# ---------------------------------------------------------------------------

def observable_map() -> Dict:
    """Map each ladder value to a known QED/spectral observable, if any.

    'connected' = already pinned to a named observable in the corpus.
    'unconnected' = framework determines the value but no physical neighbor
                    is identified.
    """
    return {
        # --- Z2(s) ---
        "Z2(2) = pi^2 - pi^4/12": {
            "status": "connected",
            "observable": (
                "one-loop QED functional determinant / vacuum-polarization "
                "spectral input on S^3; the s=2 value enters zeta'_{D^2}(0) "
                "regularization. Paper 28 Theorem T9, Paper 50 F-theorem "
                "(scalar/Dirac spectral-zeta derivative)."
            ),
            "paper34_projection": "spectral action (proj #6)",
            "transcendental_tier": "calibration_pi (M2 Seeley-DeWitt)",
        },
        "Z2(3) = pi^4/3 - pi^6/30": {
            "status": "partially-connected",
            "observable": (
                "next spectral-action heat-kernel moment (a_4-order on S^3); "
                "enters the s=3 Mellin moment Tr(D^{-6}). Not individually "
                "pinned to a measured observable -- it is an internal "
                "spectral-action coefficient, same status as a_2."
            ),
            "paper34_projection": "spectral action (proj #6)",
            "transcendental_tier": "calibration_pi (M2)",
        },
        "Z2(s>=4)": {
            "status": "unconnected-internal",
            "observable": (
                "higher Mellin moments Tr(D^{-2s}); each is an exact two-term "
                "pi^even closed form (T9) but maps to UV-suppressed "
                "spectral-action coefficients with no isolated physical "
                "measurement. Internal to the heat-kernel expansion."
            ),
            "paper34_projection": "spectral action (proj #6)",
            "transcendental_tier": "calibration_pi (M2)",
        },
        # --- D(s) free Dirichlet (s-power) ---
        "D(4) = pi^2 - pi^4/12 = Z2(2)": {
            "status": "connected",
            "observable": (
                "D(4) coincides with Z2(2): the d_max=4 Fock evaluation. "
                "This is the F=zeta(2)-adjacent value; Paper 2 alpha apparatus "
                "uses the n^2 Fock degeneracy Dirichlet at s=4."
            ),
            "paper34_projection": "Hopf/spectral action",
            "transcendental_tier": "calibration_pi (M2)",
        },
        "D(5)": {
            "status": "connected (free) / distinct from Sommerfeld D_5",
            "observable": (
                "free CH-spectrum Dirichlet at s=5: D(5) = 14 zeta(3) - "
                "31 zeta(5)/2 (Paper 28 flat-space section, R-independent). "
                "First ODD-s free-Dirichlet value: introduces odd-zeta "
                "zeta(3), zeta(5) at the FREE spectrum level (mechanism = "
                "half-integer Hurwitz shift, M3 vertex-parity)."
            ),
            "paper34_projection": "vertex parity / half-integer Hurwitz (M3)",
            "transcendental_tier": "spinor-intrinsic odd-zeta (M3)",
        },
        # --- chi_4 difference ---
        "Ddiff(s) = 2^{s-1}(beta(s)-beta(s-2))": {
            "status": "connected (identity)",
            "observable": (
                "spectral-zeta -> Dirichlet-L bridge (Paper 28 Theorem chi_4). "
                "Even-s: Catalan G, beta(4), beta(6); odd-s: pi-powers "
                "(beta(odd) closed form). NOT pinned to a measured observable "
                "-- it is the L(s,chi_{-4}) content exposed by vertex parity."
            ),
            "paper34_projection": "vertex parity (M3)",
            "transcendental_tier": "Dirichlet-L / M3 vertex-parity Hurwitz",
        },
        # --- Sommerfeld D_p (hydrogen physical-energy) ---
        "D_2 = -5/4 zeta(2) + zeta(3)": {
            "status": "connected",
            "observable": (
                "(Z alpha)^4 degeneracy-weighted Dirac fine-structure shell "
                "sum; zeta(3) is the quantum-defect injection. Paper 28."
            ),
            "paper34_projection": "Dirac fine structure",
            "transcendental_tier": "mixed even/odd zeta (M3)",
        },
        "D_5, D_6": {
            "status": "connected",
            "observable": (
                "(Z alpha)^{10}, (Z alpha)^{12} hydrogen fine-structure "
                "Sommerfeld sums; product survival rule (Paper 28 Obs). "
                "Physical = high-order relativistic energy expansion."
            ),
            "paper34_projection": "Dirac fine structure",
            "transcendental_tier": "Euler-Zagier MZV (M3)",
        },
    }


def audit_D_odd_ladder() -> Dict:
    """Audit the D(s) odd-s ladder coefficient pattern.

    Observed:  D(2m+1) = a_m * zeta(2m-1) - b_m * zeta(2m+1)
      D(5)  = 14 z(3)  - (31/2) z(5)
      D(7)  = 62 z(5)  - (127/2) z(7)
      D(9)  = 254 z(7) - (511/2) z(9)
      D(11) = 1022 z(9)- (2047/2) z(11)

    Predicted from the closed form D(s) = 2(2^{s-2}-1)zeta(s-2)
                                          - (1/2)(2^s-1)zeta(s):
      a = 2(2^{s-2}-1) = 2^{s-1}-2  (Mersenne-shifted)
      b = (1/2)(2^s-1)             (Mersenne/2)
    so at s=2m+1: a = 2^{2m}-2, b = (2^{2m+1}-1)/2.

    Audit: confirm the coefficient pattern is FORCED by the closed form
    (not curve-fit), across all odd s.
    """
    checks = []
    for s in range(5, 13, 2):
        a_pred = Fraction(2 * (2 ** (s - 2) - 1))
        b_pred = Fraction(2 ** s - 1, 2)
        # extract coeffs from sympy closed form
        exact = D_exact(s)
        c_sm2 = exact.coeff(sp.zeta(s - 2))
        c_s = exact.coeff(sp.zeta(s))
        a_act = Fraction(int(sp.numer(c_sm2)), int(sp.denom(c_sm2)))
        b_act = -Fraction(int(sp.numer(c_s)), int(sp.denom(c_s)))
        checks.append({
            "s": s,
            "a_pred": str(a_pred), "a_actual": str(a_act),
            "b_pred": str(b_pred), "b_actual": str(b_act),
            "match": (a_pred == a_act) and (b_pred == b_act),
        })
    return {
        "claim": "D(odd s) = (2^{s-1}-2) zeta(s-2) - ((2^s-1)/2) zeta(s)",
        "free_params": 0,
        "mechanism": ("Forced by Hurwitz-at-3/2 reduction; the -2^s constant "
                      "pieces cancel between the two terms (NOT fitted)."),
        "checks": checks,
        "all_match": all(c["match"] for c in checks),
    }


if __name__ == "__main__":
    print("=" * 70)
    print("Door 5: zeta_{D^2} / chi_{-4} integer-s ladder forward-run")
    print("=" * 70)

    result = run_sweep(s_lo=2, s_hi=12)

    print("\n--- Z2(s) = zeta_{D^2}(s) ladder (T9: pi^even only) ---")
    for r in result["Z2_ladder"]:
        ok = "PI^EVEN" if r["T9_pi_even_only"] else "?"
        print(f"  s={r['s']:2d}: {r['closed_form']:<28s} "
              f"(relerr {r['exact_vs_num_relerr']:.1e}) [{ok}]")

    print("\n--- D(s) free Dirac Dirichlet (s-power) ladder ---")
    for r in result["D_ladder"]:
        print(f"  s={r['s']:2d}: {r['closed_form']:<40s} "
              f"(relerr {r['exact_vs_num_relerr']:.1e})")

    print("\n--- Deven/Dodd/Ddiff chi_4 ladder ---")
    for r in result["split_ladder"]:
        v = "VERIFIED" if r["chi4_verified"] else "FAIL"
        print(f"  s={r['s']:2d}: Ddiff = {r['chi4_closed_form']:<35s} "
              f"[chi4 {v}, resid {r['chi4_residual']:.1e}]")

    print("\n--- Sommerfeld D_p (hydrogen physical energy) ---")
    for p, form in result["sommerfeld_Dp"].items():
        print(f"  D_{p} = {form}")

    audit = audit_D_odd_ladder()
    result["D_odd_ladder_audit"] = audit
    print("\n--- AUDIT: D(odd s) coefficient pattern (free-param=0) ---")
    print(f"  claim: {audit['claim']}")
    for c in audit["checks"]:
        m = "OK" if c["match"] else "MISMATCH"
        print(f"  s={c['s']:2d}: a={c['a_actual']:>8s} (pred {c['a_pred']:>8s}), "
              f"b={c['b_actual']:>8s} (pred {c['b_pred']:>8s}) [{m}]")
    print(f"  all_match: {audit['all_match']}")

    obs = observable_map()
    result["observable_map"] = obs

    print("\n--- Observable map (connected vs unconnected) ---")
    for k, v in obs.items():
        print(f"  [{v['status']:<35s}] {k}")

    # Write JSON
    out_path = "debug/data/door5_zeta_d2_ladder.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)
    print(f"\nWrote {out_path}")
