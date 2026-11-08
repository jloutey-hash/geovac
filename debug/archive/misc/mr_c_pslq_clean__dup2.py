"""Clean PSLQ for MR-C c, with non-redundant bases.

Loads c from the JSON output of mr_c_l2_subleading.py and runs PSLQ
against several non-redundant bases at multiple max_coeff levels.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
import sympy as sp
from mpmath import mp, mpf, log, pi, sqrt


def load_c():
    p = PROJ / "debug" / "data" / "mr_c_l2_subleading.json"
    with open(p) as f:
        d = json.load(f)
    return mpmath.mpf(d["c_best"])


def make_basis(name: str, dps: int):
    mp.dps = dps
    L2 = log(mpf(2))
    GE = mpmath.euler
    G = mpmath.catalan
    Z3 = mpmath.zeta(3)
    SQ_PI = sqrt(pi)

    bases = {
        # ---- M2-strict (no internal redundancies) ----
        "M2_strict_min": [
            ("1", mpf(1)),
            ("pi", pi),
            ("pi^2", pi**2),
            ("1/pi", 1/pi),
            ("1/pi^2", 1/pi**2),
            ("sqrt(pi)", SQ_PI),
            ("1/sqrt(pi)", 1/SQ_PI),
        ],
        # ---- M2 + log 2 (Stein-Weiss) ----
        "M2_plus_log2": [
            ("1", mpf(1)),
            ("pi", pi),
            ("pi^2", pi**2),
            ("1/pi", 1/pi),
            ("1/pi^2", 1/pi**2),
            ("sqrt(pi)", SQ_PI),
            ("log(2)", L2),
            ("log(2)/pi", L2/pi),
            ("log(2)*pi", L2*pi),
            ("log(2)^2", L2**2),
        ],
        # ---- M2 + log 2 + Euler-gamma (typical Euler-Maclaurin output) ----
        "M2_plus_log2_euler": [
            ("1", mpf(1)),
            ("pi", pi),
            ("pi^2", pi**2),
            ("1/pi", 1/pi),
            ("1/pi^2", 1/pi**2),
            ("sqrt(pi)", SQ_PI),
            ("log(2)", L2),
            ("log(2)/pi", L2/pi),
            ("log(2)*pi", L2*pi),
            ("euler_gamma", GE),
            ("euler_gamma/pi", GE/pi),
            ("euler_gamma*pi", GE*pi),
        ],
        # ---- Extended with M3 (Catalan, beta) and odd zeta ----
        "M2_M3_extended": [
            ("1", mpf(1)),
            ("pi", pi),
            ("pi^2", pi**2),
            ("1/pi", 1/pi),
            ("sqrt(pi)", SQ_PI),
            ("log(2)", L2),
            ("log(2)/pi", L2/pi),
            ("log(2)*pi", L2*pi),
            ("euler_gamma", GE),
            ("euler_gamma/pi", GE/pi),
            ("Catalan", G),
            ("Catalan/pi", G/pi),
            ("zeta(3)", Z3),
            ("zeta(3)/pi", Z3/pi),
        ],
        # ---- Stirling-type basis (1/2 log(2pi), Glaisher-Kinkelin, ...) ----
        "stirling_glaisher": [
            ("1", mpf(1)),
            ("pi", pi),
            ("pi^2", pi**2),
            ("log(pi)", log(pi)),
            ("log(2pi)", log(2*pi)),
            ("log(2)", L2),
            ("log(2)/pi", L2/pi),
            ("euler_gamma", GE),
            ("euler_gamma/pi", GE/pi),
            ("log_glaisher", mpmath.log(mpmath.glaisher)),
        ],
    }
    if name not in bases:
        raise KeyError(f"Unknown basis {name}; have {list(bases.keys())}")
    return bases[name]


def pslq_search(c, basis, max_coeffs=(50, 200, 1000, 10000, 100000), dps=120):
    mp.dps = dps
    targets = [c] + [v for _, v in basis]
    labels = ["c"] + [lbl for lbl, _ in basis]

    results = []
    for M in max_coeffs:
        try:
            rel = mpmath.pslq(targets, maxcoeff=M)
        except Exception as e:
            results.append({"max_coeff": M, "found": False, "error": str(e)})
            continue
        if rel is None:
            results.append({"max_coeff": M, "found": False})
            continue
        rel = [int(x) for x in rel]
        if rel[0] == 0:
            # Spurious basis-internal relation; drop and continue
            results.append({"max_coeff": M, "found": False, "skipped_rel": rel,
                            "note": "a_0=0 (basis-internal relation); ignored"})
            continue

        # c = -1/a0 * sum_{i>=1} a_i * b_i
        a0 = rel[0]
        terms = []
        for ai, lbl in zip(rel[1:], labels[1:]):
            if ai == 0:
                continue
            terms.append((sp.Rational(-int(ai), int(a0)), lbl))

        # verify numerically
        c_pred = mpf(0)
        for q, lbl in terms:
            v = next(v for l, v in basis if l == lbl)
            c_pred += mpf(q) * v
        delta = abs(c - c_pred)
        rel_err = delta / abs(c) if c != 0 else delta

        results.append({
            "max_coeff": M,
            "found": True,
            "relation": rel,
            "closed_form_terms": [(str(q), lbl) for q, lbl in terms],
            "verify_abs_error": mpmath.nstr(delta, 30),
            "verify_rel_error": mpmath.nstr(rel_err, 30),
            "verify_pass": delta < mpmath.mpf("1e-50"),
        })
        # If we got a hit at this M and it's small enough, also keep trying smaller M
    return results


def main():
    DPS = 120
    mp.dps = DPS
    c = load_c()
    print(f"c = {mpmath.nstr(c, 60)}")
    print(f"c - pi = {mpmath.nstr(c - pi, 30)}")
    print(f"c / pi = {mpmath.nstr(c / pi, 30)}")
    print(f"c - 4 = {mpmath.nstr(c - 4, 30)}")
    print(f"c * pi = {mpmath.nstr(c * pi, 30)}")
    print(f"c * pi^2 = {mpmath.nstr(c * pi**2, 30)}")
    print()

    out = {"c": mpmath.nstr(c, 80), "dps": DPS, "results_by_basis": {}}
    for basis_name in [
        "M2_strict_min",
        "M2_plus_log2",
        "M2_plus_log2_euler",
        "M2_M3_extended",
        "stirling_glaisher",
    ]:
        basis = make_basis(basis_name, DPS)
        print(f"\n=== Basis: {basis_name} ({len(basis)} elements) ===")
        for lbl, v in basis:
            print(f"  {lbl:24s} = {mpmath.nstr(v, 12)}")
        res = pslq_search(c, basis, dps=DPS)
        for r in res:
            if r.get("found"):
                cf = " + ".join(f"({co})*{lbl}" for co, lbl in r["closed_form_terms"])
                pas = "PASS" if r.get("verify_pass") else "FAIL_VERIFY"
                print(f"  [maxcoeff={r['max_coeff']}, {pas}] c = {cf}")
                print(f"    abs err = {r.get('verify_abs_error')}, rel err = {r.get('verify_rel_error')}")
            else:
                note = r.get("note", "")
                print(f"  [maxcoeff={r['max_coeff']}] no relation  {note}")
        out["results_by_basis"][basis_name] = [
            {k: v for k, v in r.items() if k != "rel"} for r in res
        ]

    # save
    out_path = PROJ / "debug" / "data" / "mr_c_pslq_clean.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\n-> {out_path}")


if __name__ == "__main__":
    main()
