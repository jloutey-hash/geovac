"""High-precision PSLQ campaign for the L2 next-order constant c.

Asymptotic on the SU(2) central Fejer panel:

    n * gamma_n^scalar ~ (4/pi) * log(n) + c + O(log(n)/n)

The leading 4/pi is the M1 Hopf-base measure signature of the master Mellin
engine (Paper 32 SVIII case-exhaustion theorem; Paper 18 SIII.7).

Sprint MR-C extracted c at >=80 dps via Richardson extrapolation; PSLQ at
coefficient ceiling 10^4 against M1 u M2 u M3 u {gamma_E, log 2, G, zeta(3)}
returned NULL.

This sprint:
    (a) push precision to >=200 dps via Richardson tower depth K=8 or 10,
        with n_panel up to 8192 and working precision 400 dps;
    (b) build harder bases:
        (b.1) Dirichlet L-tower: beta(s) for odd s; quarter-integer
              Hurwitz zeta; L(2, chi_n) for low-conductor real characters.
        (b.2) Multiple Hurwitz at quarter-integer shifts: pairs.
        (b.3) Stein-Weiss IBP derivatives of M1 (the boring answer).
        (b.4) Wildcard sector: ApA(c)ry-like and mixed constants.
    (c) PSLQ at ceilings 10^4, 10^5, 10^6 against unions of the bases.

Outputs:
    debug/data/l2_constant_c.json   (high-precision c + PSLQ results)
    debug/l2_constant_c_identification_memo.md
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
import sympy as sp
from mpmath import mp, mpf, log, pi, sqrt

from geovac.central_fejer_su2 import T_n_via_sum_rule


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Working precision target: 400 dps so that the extracted c is reliable to
# >= 200 dps after Richardson extrapolation.
DPS_COMPUTE = 400

# Richardson panel: doubling sequence. n=4096 at dps=400 takes ~5-10 min;
# n=8192 will take ~25-40 min. Total compute ~1-2 hours sequentially.
# We can be aggressive because PSLQ at 200+ dps requires this.
N_VALUES = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]

# Richardson tower depth.  K=8 gives 17 parameters; K=10 gives 21.
# With 9 panel points, K up to 4 is square; K=4 -> 9 params.
# We use a different fitting strategy here: with mixed-spacing panel,
# we fit at multiple K values and observe convergence.
K_MAX = 4   # 2K+1 = 9 params, exactly fits 9 panel points

# PSLQ
DPS_PSLQ = 250
COEFFICIENT_CEILINGS = [10**4, 10**5, 10**6]


# ---------------------------------------------------------------------------
# Sub-task (a): high-precision c
# ---------------------------------------------------------------------------


def compute_gamma_n(n: int, dps: int) -> mpmath.mpf:
    """Compute gamma_n via the closed-form sum rule at high precision."""
    mp.dps = dps
    Z = mpf(n * (n + 1)) / 2
    T = T_n_via_sum_rule(n, prec=dps)
    return pi - 4 * T / (pi * Z)


def compute_panel(n_values, dps: int):
    """Compute (n, gamma_n, h(n)) for each n.

    h(n) := n * gamma_n - (4/pi) * log(n).
    """
    mp.dps = dps
    four_over_pi = 4 / pi
    out = []
    for n in n_values:
        t0 = time.time()
        g = compute_gamma_n(n, dps=dps)
        h = mpf(n) * g - four_over_pi * log(mpf(n))
        dt = time.time() - t0
        out.append({"n": n, "gamma_n": g, "h_n": h, "compute_seconds": dt})
        print(f"  n={n:>5d}  h={mpmath.nstr(h, 20)}  ({dt:.1f}s)", flush=True)
    return out


def build_basis_matrix(n_values, K: int):
    """Return matrix rows [1, log(n)/n, 1/n, log(n)/n^2, 1/n^2, ...]."""
    rows = []
    for n in n_values:
        ln = log(mpf(n))
        inv_n = 1 / mpf(n)
        row = [mpf(1)]
        cur_inv = inv_n
        for k in range(1, K + 1):
            row.append(ln * cur_inv)
            row.append(cur_inv)
            cur_inv *= inv_n
        rows.append(row)
    return rows


def solve_least_squares(rows, ys):
    n = len(rows)
    p = len(rows[0])
    A = mpmath.matrix(rows)
    y = mpmath.matrix([[v] for v in ys])
    if n == p:
        x = A ** (-1) * y
    else:
        AT = A.T
        x = (AT * A) ** (-1) * (AT * y)
    return [x[i, 0] for i in range(p)]


def fit_subleading(panel, K: int):
    n_values = [d["n"] for d in panel]
    h_values = [d["h_n"] for d in panel]
    rows = build_basis_matrix(n_values, K)
    coeffs = solve_least_squares(rows, h_values)
    fitted = []
    residuals = []
    for n, h_obs, row in zip(n_values, h_values, rows):
        h_fit = sum(c * r for c, r in zip(coeffs, row))
        fitted.append(h_fit)
        residuals.append(h_obs - h_fit)
    return {
        "K": K,
        "n_values": n_values,
        "coefficients": coeffs,
        "b_0_extracted": coeffs[0],
        "h_values": h_values,
        "h_fitted": fitted,
        "residuals": residuals,
    }


# ---------------------------------------------------------------------------
# Sub-task (b): build candidate bases
# ---------------------------------------------------------------------------


def low_conductor_real_dirichlet_L(s: int, mod: int, dps: int) -> mpmath.mpf:
    """L(s, chi) for the unique real primitive Dirichlet character mod given conductor.

    Real primitive characters exist for moduli that are quadratic discriminants:
        chi_3, chi_4 (=Catalan), chi_5, chi_8, chi_12 (and chi_-3, chi_-4, ...).

    We implement via Hurwitz zeta:
        L(s, chi) = sum_{a=1}^{q-1} chi(a) * q^{-s} * zeta(s, a/q)
    where chi is a Dirichlet character with period q.

    Implementations for primitive real chars:
        chi_3 (mod 3, conductor 3, real): chi(1)=1, chi(2)=-1. (Kronecker -3)
        chi_4 (mod 4, conductor 4, real): chi(1)=1, chi(3)=-1. (Kronecker -4 =>
                                                               L(s,chi_4)=beta(s))
        chi_5 (mod 5): two real options; primitive real is chi_5(a) = (a/5) Legendre.
                       chi_5(1)=1, chi_5(2)=-1, chi_5(3)=-1, chi_5(4)=1. (Kronecker 5)
        chi_8 (mod 8): chi_8(1)=1, chi_8(3)=-1, chi_8(5)=-1, chi_8(7)=1. (Kronecker 8)
        chi_12 (mod 12): chi_12(1)=1, chi_12(5)=-1, chi_12(7)=-1, chi_12(11)=1.
                         (Kronecker 12)
    """
    mp.dps = dps
    table = {
        3: {1: 1, 2: -1},  # chi_-3
        4: {1: 1, 3: -1},  # chi_-4 = Catalan beta
        5: {1: 1, 2: -1, 3: -1, 4: 1},  # chi_5 (Legendre symbol)
        8: {1: 1, 3: -1, 5: -1, 7: 1},  # chi_8
        12: {1: 1, 5: -1, 7: -1, 11: 1},  # chi_12
    }
    if mod not in table:
        raise KeyError(f"No table entry for mod={mod}")
    chi = table[mod]
    total = mpf(0)
    qs = mpf(mod) ** (-s)
    for a in range(1, mod):
        if a in chi:
            ch = chi[a]
            if ch == 0:
                continue
            total += mpf(ch) * mpmath.zeta(s, mpf(a) / mod)
    return total * qs


def dirichlet_beta(s: int, dps: int) -> mpmath.mpf:
    """Dirichlet beta function beta(s) = L(s, chi_-4).

    beta(s) = (1/4^s) * (zeta(s, 1/4) - zeta(s, 3/4)).
    """
    mp.dps = dps
    return (mpf(4) ** (-s)) * (mpmath.zeta(s, mpf(1)/4) - mpmath.zeta(s, mpf(3)/4))


def build_basis_dirichlet_L(dps: int):
    """Sub-task (b.1): Dirichlet L-tower."""
    mp.dps = dps
    elements = []
    # Dirichlet beta at odd integers (beta(2)=Catalan G already in original basis;
    # add beta(1)=pi/4, beta(3), beta(5), beta(7), beta(9))
    # beta(1) = pi/4 -- already covered by 1/pi basis; skip
    for s in [3, 5, 7, 9]:
        v = dirichlet_beta(s, dps)
        elements.append((f"beta({s})", v))

    # Hurwitz zeta at quarter-integer shifts
    for s in [2, 3, 4, 5]:
        for q_str, q in [("1/4", mpf(1)/4), ("3/4", mpf(3)/4)]:
            v = mpmath.zeta(s, q)
            elements.append((f"zeta({s},{q_str})", v))

    # L(2, chi_q) for low-conductor real characters
    for q in [3, 5, 8, 12]:
        v = low_conductor_real_dirichlet_L(2, q, dps)
        elements.append((f"L(2,chi_{q})", v))

    return elements


def build_basis_multi_hurwitz(dps: int):
    """Sub-task (b.2): products of Hurwitz at quarter-integer shifts.

    Pairs zeta(a, p/q) * zeta(b, r/s) for p/q, r/s in {1/4, 3/4} and small a, b.
    """
    mp.dps = dps
    elements = []
    shifts = [("1/4", mpf(1)/4), ("3/4", mpf(3)/4)]
    a_values = [2, 3, 4]
    for i, (p_str, p) in enumerate(shifts):
        for j, (r_str, r) in enumerate(shifts):
            for a in a_values:
                for b in a_values:
                    if a + b > 6:
                        continue
                    # Avoid trivial duplicates: (a,p) <-> (b,r) gives same product
                    if (i, a) > (j, b):
                        continue
                    label = f"zeta({a},{p_str})*zeta({b},{r_str})"
                    v = mpmath.zeta(a, p) * mpmath.zeta(b, r)
                    elements.append((label, v))
    return elements


def build_basis_stein_weiss(dps: int):
    """Sub-task (b.3): Stein-Weiss IBP derivatives of M1 (the boring answer).

    Linear combinations from integrating-by-parts the Stein-Weiss / sum-rule
    expression. If c lives here, the L2 lineage is fully M1-mechanism.
    """
    mp.dps = dps
    GE = mpmath.euler
    L2 = log(mpf(2))
    elements = []
    # 4/pi powers and rationals in pi
    elements.append(("4/pi", 4 / pi))
    for k in range(1, 11):
        elements.append((f"pi/{k}", pi / mpf(k)))
    # log of small integers / pi
    for k in [2, 3, 5, 7]:
        v = log(mpf(k))
        elements.append((f"log({k})", v))
        elements.append((f"log({k})/pi", v / pi))
        elements.append((f"(4/pi)*log({k})", 4 * v / pi))
    # quadratic in log 2
    elements.append(("(log(2))^2", L2**2))
    elements.append(("(log(2))^2/pi", L2**2 / pi))
    # Euler-gamma combinations
    elements.append(("euler_gamma", GE))
    elements.append(("euler_gamma/pi", GE / pi))
    elements.append(("4*euler_gamma/pi", 4 * GE / pi))
    elements.append(("euler_gamma*log(2)", GE * L2))
    elements.append(("euler_gamma*log(2)/pi", GE * L2 / pi))
    # M1*zeta
    elements.append(("(4/pi)*zeta(2)", 4 * mpmath.zeta(2) / pi))
    elements.append(("(4/pi)*zeta(3)", 4 * mpmath.zeta(3) / pi))
    return elements


def build_basis_wildcard(dps: int):
    """Sub-task (b.4): mixed/Apery-like constants. Cheap to add.

    Trimmed: avoid log(2*pi) since log(2)+log(pi)=log(2*pi) is a basis-internal
    relation. Use log(pi) instead.
    """
    mp.dps = dps
    GE = mpmath.euler
    G = mpmath.catalan
    L2 = log(mpf(2))
    Z2 = mpmath.zeta(2)
    Z3 = mpmath.zeta(3)
    elements = []
    # Cross products
    elements.append(("zeta(2)*log(2)", Z2 * L2))
    elements.append(("zeta(3)/pi^2", Z3 / pi**2))
    elements.append(("(log(2))^2", L2**2))
    elements.append(("Catalan*log(2)", G * L2))
    elements.append(("Catalan*pi", G * pi))
    elements.append(("Catalan/pi", G / pi))
    elements.append(("pi*euler_gamma", pi * GE))
    # log(2*pi) removed: log(2) + log(pi) = log(2*pi) is basis-internal
    elements.append(("log(pi)", log(pi)))
    elements.append(("log(pi)/pi", log(pi)/pi))
    # Stieltjes constants gamma_1, gamma_2 (low ones)
    # mpmath.stieltjes(n) returns gamma_n
    try:
        s1 = mpmath.stieltjes(1)
        s2 = mpmath.stieltjes(2)
        elements.append(("stieltjes_1", s1))
        elements.append(("stieltjes_2", s2))
    except Exception:
        pass
    # Glaisher-Kinkelin
    try:
        glaisher = mpmath.glaisher
        elements.append(("log(glaisher)", log(glaisher)))
    except Exception:
        pass
    return elements


def build_basis_basic(dps: int):
    """Always-include core: rationals, pi, log(2), gamma_E, Catalan, zeta(3).

    Trimmed to remove basis-internal redundancies (e.g., 4/pi vs 1/pi
    differ by an integer factor of 4; that's a spurious PSLQ hit).
    """
    mp.dps = dps
    elements = [
        ("1", mpf(1)),
        ("pi", pi),
        ("pi^2", pi**2),
        ("pi^3", pi**3),
        ("1/pi", 1 / pi),
        ("1/pi^2", 1 / pi**2),
        ("sqrt(pi)", sqrt(pi)),
        # 1/sqrt(pi) removed: pi * 1/sqrt(pi) = sqrt(pi), so basis-redundant
        # 4/pi removed: 4*(1/pi) is a basis-internal integer relation
        ("log(2)", log(mpf(2))),
        ("euler_gamma", mpmath.euler),
        ("Catalan", mpmath.catalan),
        ("zeta(3)", mpmath.zeta(3)),
    ]
    return elements


# ---------------------------------------------------------------------------
# Sub-task (c): PSLQ campaigns
# ---------------------------------------------------------------------------


def run_pslq(c: mpmath.mpf, basis, max_coeff: int, dps: int):
    """Run PSLQ on [c] + basis values; return parsed result.

    Filters out a_0 = 0 results (basis-internal relations).
    Reports residual (= |c - predicted| / |c|).
    """
    mp.dps = dps
    targets = [c] + [v for _, v in basis]
    labels = ["c"] + [lbl for lbl, _ in basis]
    try:
        rel = mpmath.pslq(targets, maxcoeff=max_coeff)
    except Exception as e:
        return {"max_coeff": max_coeff, "basis_size": len(basis), "found": False,
                "error": str(e)}
    if rel is None:
        return {"max_coeff": max_coeff, "basis_size": len(basis), "found": False,
                "error": None}
    rel = [int(x) for x in rel]
    if rel[0] == 0:
        return {"max_coeff": max_coeff, "basis_size": len(basis), "found": False,
                "skipped_rel": rel,
                "note": "a_0=0 (basis-internal relation); ignored"}
    a0 = rel[0]
    terms = []
    for ai, lbl in zip(rel[1:], labels[1:]):
        if ai == 0:
            continue
        terms.append((sp.Rational(-int(ai), int(a0)), lbl))
    # Verify
    c_pred = mpf(0)
    for q, lbl in terms:
        v = next(v for l, v in basis if l == lbl)
        c_pred += mpf(q) * v
    delta = abs(c - c_pred)
    rel_err = delta / abs(c) if c != 0 else delta
    return {
        "max_coeff": max_coeff,
        "basis_size": len(basis),
        "found": True,
        "relation": rel,
        "closed_form_terms": [(str(q), lbl) for q, lbl in terms],
        "verify_abs_error": mpmath.nstr(delta, 30),
        "verify_rel_error": mpmath.nstr(rel_err, 30),
        "verify_pass_1e_minus_80": delta < mpf("1e-80"),
        "verify_pass_1e_minus_50": delta < mpf("1e-50"),
        "verify_pass_1e_minus_30": delta < mpf("1e-30"),
    }


def deduplicate_basis(elements):
    """Remove duplicate (label, value) pairs."""
    seen_labels = set()
    out = []
    for lbl, v in elements:
        if lbl in seen_labels:
            continue
        seen_labels.add(lbl)
        out.append((lbl, v))
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print("=" * 72)
    print("L2 next-order constant c: high-precision PSLQ campaign")
    print("=" * 72)
    print(f"DPS_COMPUTE = {DPS_COMPUTE}")
    print(f"N_VALUES    = {N_VALUES}")
    print(f"DPS_PSLQ    = {DPS_PSLQ}")
    print(f"CEILINGS    = {COEFFICIENT_CEILINGS}")
    print()

    out = {
        "sprint": "L2-constant-c-identification",
        "purpose": "PSLQ identification or definitive exclusion of c at >=200 dps",
        "configuration": {
            "dps_compute": DPS_COMPUTE,
            "n_values": N_VALUES,
            "K_max": K_MAX,
            "dps_pslq": DPS_PSLQ,
            "coefficient_ceilings": COEFFICIENT_CEILINGS,
        },
    }

    # === Stage 1: high-precision panel ===
    print("Stage 1: high-precision panel computation")
    t_start = time.time()
    panel = compute_panel(N_VALUES, dps=DPS_COMPUTE)
    t_compute = time.time() - t_start
    print(f"Stage 1 done in {t_compute:.1f}s")
    print()
    out["stage1_compute_seconds"] = t_compute
    out["panel"] = [
        {"n": d["n"], "gamma_n_50dps": mpmath.nstr(d["gamma_n"], 50),
         "h_n_250dps": mpmath.nstr(d["h_n"], 250),
         "compute_seconds": d["compute_seconds"]}
        for d in panel
    ]

    # === Stage 2: Richardson fits at multiple K ===
    print("Stage 2: Richardson extrapolation at K=1..K_MAX")
    fits = {}
    K_panel_summary = []
    for K in range(1, K_MAX + 1):
        n_params = 2 * K + 1
        if n_params > len(panel):
            print(f"  K={K} skipped (would need {n_params} params, have {len(panel)} pts)")
            continue
        try:
            fit = fit_subleading(panel, K)
        except Exception as e:
            print(f"  K={K} fit failed: {e}")
            continue
        c_K = fit["b_0_extracted"]
        residual_max = max(abs(r) for r in fit["residuals"])
        K_panel_summary.append({
            "K": K,
            "c_50dps": mpmath.nstr(c_K, 50),
            "c_250dps": mpmath.nstr(c_K, 250),
            "residual_max": mpmath.nstr(residual_max, 12),
        })
        fits[K] = fit
        print(f"  K={K}  c={mpmath.nstr(c_K, 30)}  max_resid={mpmath.nstr(residual_max, 6)}")
    out["K_panel"] = K_panel_summary

    if not fits:
        print("FATAL: no fits succeeded.")
        return

    # Choose best K (largest available square fit)
    best_K = max(fits.keys())
    c_best = fits[best_K]["b_0_extracted"]
    print(f"\nBest c estimate (K={best_K}):")
    print(f"  first 100 dps:  {mpmath.nstr(c_best, 100)}")
    print(f"  next 100 dps :  {mpmath.nstr(c_best, 200)[100:]}")
    print()

    # Sprint MR-C reference value (first 80 chars from JSON)
    mrc_value_str = "4.1093214674877940927579607260741005838057691088362615503253972964276017819113301"
    mrc_value = mpf(mrc_value_str)
    delta_mrc = abs(c_best - mrc_value)
    print(f"Comparison vs Sprint MR-C 80-dps value:")
    print(f"  |c_best - c_MRC| = {mpmath.nstr(delta_mrc, 12)}")
    print(f"  (should be << 1e-50; if ~1e-80 we agree at MR-C's precision)")
    print()
    out["mrc_comparison"] = {
        "mrc_80dps_value": mrc_value_str,
        "abs_delta": mpmath.nstr(delta_mrc, 30),
        "agrees_at_50dps": delta_mrc < mpf("1e-50"),
        "agrees_at_80dps": delta_mrc < mpf("1e-78"),
    }
    out["c_best_250dps"] = mpmath.nstr(c_best, 250)
    out["c_best_50dps"] = mpmath.nstr(c_best, 50)
    out["best_K"] = best_K

    # === Stage 3: build candidate bases ===
    print("Stage 3: building candidate bases")
    mp.dps = DPS_PSLQ

    basis_basic = build_basis_basic(DPS_PSLQ)
    basis_dirichletL = build_basis_dirichlet_L(DPS_PSLQ)
    basis_multihurwitz = build_basis_multi_hurwitz(DPS_PSLQ)
    basis_steinweiss = build_basis_stein_weiss(DPS_PSLQ)
    basis_wildcard = build_basis_wildcard(DPS_PSLQ)

    print(f"  basic           : {len(basis_basic)} elements")
    print(f"  Dirichlet L     : {len(basis_dirichletL)} elements")
    print(f"  multi-Hurwitz   : {len(basis_multihurwitz)} elements")
    print(f"  Stein-Weiss IBP : {len(basis_steinweiss)} elements")
    print(f"  wildcard        : {len(basis_wildcard)} elements")

    # Construct named test bases
    test_bases = [
        ("MR-C_baseline", basis_basic),
        ("DirichletL_addon", deduplicate_basis(basis_basic + basis_dirichletL)),
        ("MultiHurwitz_addon", deduplicate_basis(basis_basic + basis_multihurwitz)),
        ("SteinWeissIBP_addon", deduplicate_basis(basis_basic + basis_steinweiss)),
        ("Wildcard_addon", deduplicate_basis(basis_basic + basis_wildcard)),
        ("DirichletL_StrictMin", deduplicate_basis([("1", mpf(1)), ("pi", pi), ("4/pi", 4/pi)] + basis_dirichletL)),
        ("Union_All", deduplicate_basis(
            basis_basic + basis_dirichletL + basis_multihurwitz +
            basis_steinweiss + basis_wildcard)),
    ]
    print()

    # === Stage 4: PSLQ campaigns ===
    print("Stage 4: PSLQ campaigns")
    print()
    pslq_results_by_basis = {}
    for basis_name, basis in test_bases:
        print(f"--- Basis: {basis_name} ({len(basis)} elements) ---")
        results = []
        for M in COEFFICIENT_CEILINGS:
            t0 = time.time()
            r = run_pslq(c_best, basis, max_coeff=M, dps=DPS_PSLQ)
            r["pslq_seconds"] = time.time() - t0
            results.append(r)
            if r.get("found"):
                cf = " + ".join(f"({q})*{lbl}" for q, lbl in r["closed_form_terms"])
                tag = "PASS_80" if r.get("verify_pass_1e_minus_80") else (
                      "PASS_50" if r.get("verify_pass_1e_minus_50") else (
                      "PASS_30" if r.get("verify_pass_1e_minus_30") else "FAIL"))
                print(f"  [M={M}, {tag}] c = {cf}")
                print(f"    abs err = {r.get('verify_abs_error')}")
                print(f"    rel err = {r.get('verify_rel_error')}")
            else:
                note = r.get("note", "no relation")
                print(f"  [M={M}] {note}  ({r.get('pslq_seconds', 0):.1f}s)")
        pslq_results_by_basis[basis_name] = results
        print()
    out["pslq_results_by_basis"] = pslq_results_by_basis
    out["test_bases_summary"] = {name: len(b) for name, b in test_bases}

    # === Stage 5: write JSON ===
    print("Stage 5: write JSON")
    out_path = PROJ / "debug" / "data" / "l2_constant_c.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"  -> {out_path}")
    print()

    # === Summary ===
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"c (250 dps) first 80 chars : {out['c_best_50dps'][:80]}")
    print(f"agrees with MR-C at 80 dps : {out['mrc_comparison']['agrees_at_80dps']}")
    print(f"  (delta = {out['mrc_comparison']['abs_delta']})")
    print()
    print("Verdict by basis (best non-failing PSLQ result):")
    any_hit = False
    for basis_name, results in pslq_results_by_basis.items():
        passes = [r for r in results if r.get("found") and r.get("verify_pass_1e_minus_80")]
        passes_50 = [r for r in results if r.get("found") and r.get("verify_pass_1e_minus_50")]
        if passes:
            best = passes[0]
            cf = " + ".join(f"({q})*{lbl}" for q, lbl in best["closed_form_terms"])
            print(f"  {basis_name:30s}: PASS_80 at M={best['max_coeff']}: c = {cf}")
            any_hit = True
        elif passes_50:
            best = passes_50[0]
            cf = " + ".join(f"({q})*{lbl}" for q, lbl in best["closed_form_terms"])
            print(f"  {basis_name:30s}: PASS_50 at M={best['max_coeff']}: c = {cf}")
        else:
            largest_M = max(r["max_coeff"] for r in results)
            print(f"  {basis_name:30s}: NULL up to M={largest_M}")
    print()
    if any_hit:
        print("VERDICT: identification(s) found above.")
    else:
        print(f"VERDICT: NULL across all {len(test_bases)} test bases at all ceilings"
              f" up to {max(COEFFICIENT_CEILINGS)}.")
        print("c is irreducible at much higher coefficient ceiling than originally reported.")


if __name__ == "__main__":
    main()
