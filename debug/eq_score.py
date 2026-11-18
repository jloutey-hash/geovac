"""Score record EQ (docs/predictions/2026-11-17_exchange_d_pp.md) from the eq_run outputs (suffix EQS).
Reads only the fields the record names; writes debug/data/eq_score.json. Usage: python debug/eq_score.py"""
import json
from fractions import Fraction
from pathlib import Path

D = Path(__file__).resolve().parent / "data"
P1, P2 = 67108859, 33554393
CASES = [("pppi0", 1, 1), ("3ddel", 1, 1), ("ppsig", 1, 1), ("3dsig", 1, 1),
         ("pppi0", 1, 2), ("3ddel", 1, 2), ("ppsig", 1, 2)]


def load(name):
    f = D / name
    return json.loads(f.read_text()) if f.exists() else None


def allowed_rates(a, b):
    return {Fraction(x) for x in (a + b, -(a + b), 2 * a, -2 * a, 2 * b, -2 * b, 0, a - b, b - a)}


def score_case(case, a, b):
    tag = f"{case}_{a}_{b}"
    lad, ho, nu = load(f"eq_ladder_{tag}_EQS.json"), load(f"eq_heldout_{tag}_EQS.json"), load(f"eq_numeric_{tag}_EQS.json")
    runs = {p: load(f"eq_run_{tag}_{p}_EQS.json") for p in (P1, P2)}
    out = {"case": tag, "complete": all(x is not None for x in (lad, ho, nu, runs[P1], runs[P2]))}
    if not out["complete"]:
        return out
    errors = [p for p, r in runs.items() if "error" in r]
    out["error_on_scored_prime"] = errors
    mo = {p: r.get("minimal_operator") for p, r in runs.items()}
    shapes = {p: (m["order"], m["degree"], tuple(m["coeff_degrees"])) if m else None for p, m in mo.items()}
    out["order_degree"] = [shapes[P1][0], shapes[P1][1]] if shapes[P1] else None
    out["shapes_agree"] = shapes[P1] == shapes[P2] and shapes[P1] is not None
    out["first_hit_null_both"] = all(r.get("first_hit") is None for r in runs.values())
    tp = ho["third_prime"]
    determined = tp["fit_equations"] >= tp["unknowns"] + 50
    out["third_prime_determined"] = determined
    # third-prime terms count only when the fit is determined (record: determinacy rule; underdetermined -> BORDERLINE)
    third_ok = (tp["fit_nullity"] == 1 and tp["fit_matches_crt"] and tp["crt_kills_all_rows_incl_heldout"]) if determined else True
    eq1 = (out["shapes_agree"] and lad.get("scored_primes_agree") is True and lad.get("k0") is not None and lad["k0"] <= 92
           and lad.get("k0_plus_4_agrees") is True and ho["exact_Q"]["residual_zero"] is True
           and third_ok and nu["numeric_pass"] is True)
    out["EQ1"] = bool(eq1 and determined)
    out["EQ1_instrument_fail"] = bool(out["shapes_agree"] and not eq1)
    # scored primes disagreeing in shape (both operators found) is VOID per the gate
    out["shape_mismatch"] = bool(shapes[P1] is not None and shapes[P2] is not None and shapes[P1] != shapes[P2])
    out["k0"], out["accepted"], out["crt_bits"] = lad.get("k0"), len(lad.get("accepted", [])), ho.get("crt_bits")
    order = shapes[P1][0] if shapes[P1] else None
    if (a, b) == (1, 1):
        out["EQ2a"] = order == 4
    else:
        out["EQ2b"] = order is not None and 9 <= order <= 13
    verdicts = []
    for p, r in runs.items():
        for pt in r.get("finite_nonzero_points", []):
            verdicts.append((p, str(pt.get("verdict"))))
    out["finite_point_verdicts"] = sorted({v for _, v in verdicts})
    genuine = {p: any(v.startswith("genuine") for q, v in verdicts if q == p) for p in (P1, P2)}
    out["EQ3"] = all(v.startswith("apparent") or v.startswith("exponents_only: apparent-type") for _, v in verdicts)
    out["genuine_both"] = genuine[P1] and genuine[P2]
    out["genuine_any"] = genuine[P1] or genuine[P2]
    r0 = {p: r["R0"] for p, r in runs.items()}
    out["R0_exponents"] = r0[P1].get("exponents")
    out["EQ4"] = all(z.get("fuchsian") is True and z.get("exponents_not_in_Fp", 1) == 0
                     and z.get("log_free_local_solutions", 99) < order and z.get("f_has_lnR") is True
                     for z in r0.values())
    out["log_free"] = r0[P1].get("log_free_local_solutions")
    inf = {p: r["infinity"] for p, r in runs.items()}
    rates_ok = lad.get("scored_rates_all_rational")
    rates_ok = all(rates_ok.values()) if isinstance(rates_ok, dict) else bool(rates_ok)
    sr = lad.get("scored_rates", {})
    rate_lists = [sorted(Fraction(x) for x, mult in v for _ in range(mult)) for v in sr.values()] if isinstance(sr, dict) else []
    same = len(rate_lists) == 2 and rate_lists[0] == rate_lists[1]
    inside = all(x in allowed_rates(a, b) for lst in rate_lists for x in lst)
    eq5 = all(i["rank_le_1"] for i in inf.values()) and rates_ok and same and inside
    if (a, b) == (1, 1) and rate_lists:
        eq5 = eq5 and set(rate_lists[0]) <= {Fraction(0), Fraction(2), Fraction(-2)}
    out["rates"] = [str(x) for x in rate_lists[0]] if rate_lists else None
    out["EQ5"] = bool(eq5)
    lc = ho.get("f_log_components_exact", {})
    out["lnR_onset"] = {k: v for k, v in lc.items()}
    out["numeric_worst_dps80"] = max(float(nu[R]["abs_diff_frob_dps80"]) for R in ("1.4", "3.015"))
    return out


def main():
    gate = load("eq_gate_EQS.json")
    res = {"EQ0_gate_pass": bool(gate and gate.get("gate_pass")), "cases": [score_case(*c) for c in CASES]}
    cs = res["cases"]
    if all(c["complete"] for c in cs):
        void = (not res["EQ0_gate_pass"]) or any(c["error_on_scored_prime"] or c["EQ1_instrument_fail"] or c.get("shape_mismatch")
                                                 or (c["case"].endswith("_1_1") and (c["order_degree"] or [0])[0] >= 6)
                                                 or (c["case"].endswith("_1_1") and c["genuine_any"]) for c in cs)
        stop = any(c["first_hit_null_both"] or (c["case"].endswith("_1_2") and c["genuine_both"]) for c in cs)
        go = all(c["EQ1"] and c.get("EQ2a", c.get("EQ2b")) and c["EQ3"] and c["EQ4"] and c["EQ5"] for c in cs)
        res["verdict"] = "VOID" if void else "STOP" if stop else "GO" if go else "BORDERLINE"
    else:
        res["verdict"] = "INCOMPLETE"
    (D / "eq_score.json").write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
