"""MG1-MG3: assemble debug/data/mg_pieces.json (everything exact, as strings) from the per-gamma stage outputs,
for the PM recheck.  Per gamma: the kite operator in z and its singular points, Y_log's annihilator L, the
composition chain and factors, the Legendre map, the K3 map, local data, residuals on coefficients, on phi and
on Y, the 1/200 controls, and the mod-p checks.  Singular-point classification: a rational root of the leading
coefficient is APPARENT when its exponents are distinct non-negative integers, all log-free; otherwise GENUINE.
A non-linear irreducible factor is apparent when the mod-p check finds r distinct non-negative integer
log-free exponents at its roots.
usage: python mg_assemble.py"""
from __future__ import annotations
import json
from fractions import Fraction
from pathlib import Path

DATA = Path(__file__).parent / "data"
GAMMAS = {"0.5": {"gamma": "1/2", "gstr": "1_2", "series": "mg_series_g0.5_M1600_d600_rho0.5.json", "rho": "1/2"},
          "3": {"gamma": "3", "gstr": "3", "series": "mg_series_g3_M1600_d600_rho2.json", "rho": "2"}}


def load(name):
    f = DATA / name
    return json.loads(f.read_text()) if f.exists() else None


def classify_point(exps, classes):
    """exps: list of exponent strings; classes: list of lists of (e, mult, logfree)"""
    es = [Fraction(e) for e in exps]
    classes = [cl["exponents"] if isinstance(cl, dict) else cl for cl in classes]   # survey: dicts; kite local: tuple lists
    logfree = all(ok for cl in classes for e, m, ok in cl)
    distinct = len(set(es)) == len(es)
    nonneg_int = all(e.denominator == 1 and e >= 0 for e in es)
    return "apparent" if (logfree and distinct and nonneg_int) else "genuine"


def classify_lead(local_points, lead, apparent_check, order):
    out = {"genuine": [], "apparent": [], "nonlinear_factors": []}
    for pt, v in local_points.items():
        out[classify_point(v["exps"], v["mono"] if "mono" in v else v["classes"])].append(pt)
    for key, ac in (apparent_check or {}).items():
        checks = ac["check"] if isinstance(ac, dict) and "check" in ac else ac
        deg = ac.get("degree") if isinstance(ac, dict) else None
        appar = all(len(set(c["int_exps"])) == order and all(e >= 0 for e in c["int_exps"]) and all(c["logfree"]) for c in checks)
        out["nonlinear_factors"].append({"factor_head": key, "degree": deg, "apparent": appar, "check": checks})
    return out


def kite_block(tag, g):
    op = load(f"mg_5_op_kite_g{tag}.json")
    loc = load(f"mg_kite_local_g{tag}.json")
    onphi = load(f"mg_9k_onphi_g{tag}.json")
    srch = {k: v for k, v in (load(f"mg_5_search_phi_{Path(g['series']).stem}_N800.json") or {}).items()}
    modp = [v for v in [load(p.name) for p in DATA.glob(f"mg_5e_modp_check_kite_g{tag}*.json")] if v]
    b = {"operator_Qform_theta_z": op, "search_orders": srch, "local_data": loc, "on_phi_quadrature": onphi, "modp_checks": modp}
    if loc:
        pts = {k: {"exps": v["exps"], "mono": v["classes"]} for k, v in loc["A"]["points"].items()}
        # z = 0 and infinity are singular by the kite's own exponents (non-integer / negative at 0; MUM)
        cl = classify_lead(pts, loc["leading_factors"], loc.get("apparent_check"), loc["order"])
        z0_exps = [Fraction(e) for e in loc["A"]["z0_exponents"]]
        cl["z0"] = "genuine" if (any(e.denominator != 1 or e < 0 for e in z0_exps) or
                                 not all(ok for c in loc["A"]["z0_classes"] for e, m, ok in c)) else "apparent"
        cl["predicted_points"] = loc["predicted_points"]
        cl["genuine_finite_nonzero"] = sorted(cl["genuine"], key=lambda s: Fraction(s))
        b["singular_points"] = cl
    return b


def y_block(tag, g):
    op = load(f"mg_5_op_L9_g{tag}.json")
    srch = load(f"mg_5_search_ylog_{Path(g['series']).stem}_N800.json")
    surv = load(f"mg_6_survey_g{tag}.json")
    chain = load(f"mg_7_pieces_g{tag}.json")
    pieces = load(f"mg_8_pieces_g{tag}.json")
    ord1 = load(f"mg_8b_order1_g{tag}.json")
    onY = load(f"mg_9_onY_g{tag}.json")
    k3 = load(f"mg_10_k3map_g{tag}.json")
    lmap = load(f"mg_map_factor_C2b_g{g['gstr']}.json")
    modp = [v for v in [load(p.name) for p in DATA.glob(f"mg_5e_modp_check_L9_g{tag}*.json")] if v]
    b = {"annihilator_Qform": op, "search_orders": srch, "survey": surv, "chain_and_factors": chain,
         "factor_local_data": pieces, "order_one_solutions": ord1, "on_Y_quadrature": onY, "k3_map": k3,
         "legendre_map": lmap, "modp_checks": modp}
    if surv:
        Lrep = surv["L"]
        pts = {k: {"exps": v["exps"], "mono": v["mono"]} for k, v in Lrep["points"].items()}
        cl = classify_lead(pts, Lrep["lead"], surv.get("apparent"), Lrep["order"])
        cl["genuine_x"] = sorted(cl["genuine"], key=lambda s: Fraction(s))
        b["singular_points"] = cl
    if chain:
        b["composition_series"] = {"chain_orders": chain["summary"]["chain_orders"],
                                   "factor_orders": chain["summary"]["factor_orders"],
                                   "divisions_exact": chain["summary"]["chain_divisions_exact"]}
    if lmap and "map" in lmap:
        m = lmap["map"]
        b["legendre_map_summary"] = {"twist": m["twist"], "E": m["E"], "deg_N0": m["dN"], "nullity": m["nullity"],
                                     "exact": m["exact"], "perturbed_1e-3_breaks": m["perturbed_1e-3_breaks"],
                                     "z": f"(t^2 + 1 - gamma^2)/(2t) at gamma = {g['gamma']}"}
    return b


def sqrt_frac(s):
    q = Fraction(s)
    from math import isqrt
    n, d = isqrt(q.numerator), isqrt(q.denominator)
    return str(Fraction(n, d)) if Fraction(n, d) ** 2 == q else None


out = {"record": "docs/predictions/2026-10-23_gamma1_mechanism.md", "pipeline": "debug/mg_*.py (copies of cf_g2_*/ga1_* stages)",
       "conventions": {"kite": "sum_i z^i Q_i(theta_z), theta_z = z d/dz, acting on phi(z) = Phi~(sqrt z)",
                       "Y": "sum_i t^(2i) Q_i(theta_t) on Y_log = t^3 sum binom(-3/2,n) phi_n t^(2n); theta-form lists are "
                            "sum_k c_k(x) theta_t^k with x = t^2, flint polynomial strings in x",
                       "maps": "Legendre: A = (N0 + N1 d/dx)/E from L_Leg(gamma) (degree -1/2 in z = (t^2+1-gamma^2)/(2t), twisted "
                               "by x^(1/2)) into C2b; K3: A = t^s (x-4)^-b sum_j N_j(x) theta^j from F_3's operator into F"},
       "controls": {"gamma2_legendre_map": load("mg_map_factor_C2b_g2.json"),
                    "gamma2_chain_reproduction": (load("mg_7_pieces_ctrl_g2.json") or {}).get("summary"),
                    "gamma2_L9_rational": load("mg_5_op_ctrl_g2_L9.json"),
                    "gamma2_kite": {"op": load("mg_5_op_kite_ctrl_g2.json"), "local": load("mg_kite_local_ctrl_g2.json"),
                                    "search": load("mg_5_search_phi_cf_g2_series_g2_M1600_d600_rho2_N800.json")},
                    "gamma1.5_kite_local": load("mg_kite_local_ctrl_g1.5.json"),
                    "gamma1.5_kite_on_phi": load("mg_9k_onphi_ctrl_g1.5.json"),
                    "gamma1_kite": {"op": load("mg_5_op_kite_g1.json"), "local": load("mg_kite_local_g1.json"),
                                    "search": load("mg_5_search_phi_mg_series_g1_M256_d120_rho1_N128.json")}},
       "gammas": {}}
for tag, g in GAMMAS.items():
    ser = load(g["series"])
    blk = {"gamma": g["gamma"], "series": {k: v for k, v in (ser or {}).items() if k != "phi_n"},
           "kite": kite_block(tag, g), "Y": y_block(tag, g)}
    if "singular_points" in blk["Y"]:
        blk["Y"]["singular_points"]["genuine_t_positive"] = [sqrt_frac(s) for s in blk["Y"]["singular_points"]["genuine_x"]]
    out["gammas"][g["gamma"]] = blk
(DATA / "mg_pieces.json").write_text(json.dumps(out, indent=1, default=str))
for gname, blk in out["gammas"].items():
    print("gamma =", gname)
    k = blk["kite"]
    if k.get("operator_Qform_theta_z"):
        print("  kite: order", k["operator_Qform_theta_z"]["r"], "z-degree", k["operator_Qform_theta_z"]["d"],
              "| genuine finite points", k.get("singular_points", {}).get("genuine_finite_nonzero"),
              "| apparent", k.get("singular_points", {}).get("apparent"),
              "| nonlinear", [(f["degree"], f["apparent"]) for f in k.get("singular_points", {}).get("nonlinear_factors", [])])
    y = blk["Y"]
    if y.get("annihilator_Qform"):
        print("  Y: order", y["annihilator_Qform"]["r"], "x-degree", y["annihilator_Qform"]["d"],
              "| genuine t>0", y.get("singular_points", {}).get("genuine_t_positive"),
              "| apparent x", y.get("singular_points", {}).get("apparent"),
              "| nonlinear", [(f["degree"], f["apparent"]) for f in y.get("singular_points", {}).get("nonlinear_factors", [])])
    if y.get("composition_series"):
        print("  series:", y["composition_series"])
    if y.get("legendre_map_summary"):
        print("  Legendre map:", y["legendre_map_summary"])
    if y.get("k3_map"):
        print("  K3 map:", {kk: y["k3_map"].get(kk) for kk in ("t_power", "x_minus_4_power", "numerator_degree", "identity_exact")})
print("written", DATA / "mg_pieces.json")
