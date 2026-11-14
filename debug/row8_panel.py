"""Record R8 (docs/predictions/2026-11-14_row8_two_centre_l2.md): scoring driver.

Two-centre integrals over hydrogenic orbitals with l <= 2 (rows 2, 3a, 3b) run through the
unchanged R-operator analysis of debug/two_centre_r_operators.py.  The panel is the record's,
fixed at registration.  Writes debug/data/row8_panel.json.

Usage:  python debug/row8_panel.py [--only CASEID ...]
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "debug"))

import row8_one_electron as E  # noqa: E402
import two_centre_r_operators as T  # noqa: E402
from geovac import qfd_core as Q  # noqa: E402
from geovac.two_center_eri import aabb_closed_form, hybrid_closed_form  # noqa: E402

R = T.R
F = Fraction
RV = sp.Rational(7, 5)


def o(s):
    """'3d-1' -> (3, 2, -1)."""
    n = int(s[0])
    l = "spd".index(s[1])
    m = int(s[2:]) if len(s) > 2 else 0
    return (n, l, m)


def A(s, Z):
    return ("A", F(Z)) + o(s)


def B(s, Z):
    return ("B", F(Z)) + o(s)


ROW2_PAIRS = {  # id: (a on A, Z_A, b on B, Z_B)
    "E1": ("2p0", 3, "1s", 1), "E2": ("2p1", 2, "2p1", 2), "E3": ("3d0", 3, "2p0", 2),
    "E4": ("3d2", 3, "3d2", F(3, 2)), "E5": ("3p1", 3, "3d1", 2), "E6": ("2s", 2, "3d0", 3),
}
ROW2_DENS = {  # id: (a, a', Z) both on A, potential at B
    "F1": ("2p0", "2p0", 3), "F2": ("2p1", "2p1", 2), "F3": ("3d0", "3s", 3),
    "F4": ("3d1", "2p1", 3), "F5": ("3d2", "3d2", 3), "F6": ("2p0", "2s", 2),
}
ROW3A = {  # id: (a, b, Z_A, c, d, Z_B)
    "C1": ("3d0", "3d0", 3, "1s", "1s", 1), "C2": ("2p1", "2p1", 2, "2p0", "2p0", 2),
    "C3": ("3d1", "3d1", 3, "3d2", "3d2", F(3, 2)), "C4": ("2p0", "3d0", 3, "1s", "2p0", 2),
    "C5": ("3d2", "3d0", 3, "3d0", "3d2", 3), "C6": ("2p1", "2p0", 2, "3d0", "3d1", 3),
    "C7": ("2s", "3p0", 3, "2p0", "3s", 2), "C8": ("3d2", "3d2", 3, "3d2", "3d2", 3),
}
ROW3B = {  # id: (a, b, c, Z_A, d, Z_B, relation)
    "H1": ("3d0", "3d0", "3d0", 3, "1s", 1, "eq"), "H2": ("2p1", "2p1", "2p1", 2, "2p1", 2, "eq"),
    "H3": ("3d0", "3d0", "1s", 3, "1s", 1, "lt"), "H4": ("2p0", "2p0", "2p0", 2, "3d0", 3, "eq"),
    "H5": ("3d1", "3d1", "3d1", 3, "2p1", 3, "gt"), "H6": ("3d2", "3d2", "3s", 3, "2s", 1, "lt"),
    "H7": ("2p0", "3d0", "2s", 2, "1s", 1, "eq"), "H8": ("3p1", "3d1", "3d0", 3, "2p0", 4, "gt"),
    "H9": ("2s", "2s", "3d0", 2, "3d0", 3, "gt"), "H10": ("3d2", "3d0", "3d-1", 3, "3d1", 3, "eq"),
}


def rate(s, Z):
    return F(Z) / int(s[0])


def panel():
    cases = []
    for cid, (a, za, b, zb) in ROW2_PAIRS.items():
        oa, ob = A(a, za), B(b, zb)
        for kind, fn in (("S", E.S), ("VA", E.VA), ("VB", E.VB), ("T", E.T_lap)):
            cases.append((f"{cid}{kind}", "2", f"<{a}[{za}]|{kind}|{b}_B[{zb}]>",
                          lambda fn=fn, oa=oa, ob=ob: fn(oa, ob), (oa, ob, kind)))
    for cid, (a, a2, z) in ROW2_DENS.items():
        oa, ob = A(a, z), A(a2, z)
        cases.append((cid, "2", f"<{a}|1/r_B|{a2}>[{z}]", lambda oa=oa, ob=ob: E.VB(oa, ob),
                      (oa, ob, "VB")))
    for cid, (a, b, za, c, d, zb) in ROW3A.items():
        cases.append((cid, "3a", f"({a} {b}|{c} {d})[{za};{zb}]",
                      lambda a=a, b=b, za=za, c=c, d=d, zb=zb: aabb_closed_form(
                          F(za), o(a), o(b), F(zb), o(c), o(d)), None))
    for cid, (a, b, c, za, d, zb, rel) in ROW3B.items():
        ac, ad = rate(c, za), rate(d, zb)
        got = "eq" if ac == ad else ("gt" if ad > ac else "lt")
        assert got == rel, (cid, ac, ad)
        cases.append((cid, "3b", f"({a} {b}|{c} {d}_B)[{za};{zb}] {rel}",
                      lambda a=a, b=b, c=c, za=za, d=d, zb=zb: hybrid_closed_form(
                          F(za), o(a), o(b), o(c), F(zb), o(d)), None))
    return cases


def controls():
    out = {}
    for name, fn in (("K1_dilog", lambda: sp.exp(-R) * (1 + 1 / R) + sp.polylog(2, sp.exp(-R))),
                     ("K1_E1sq", lambda: sp.exp(-R) + sp.expint(1, R) ** 2)):
        try:
            T.decompose(fn())
            out[name] = "ACCEPTED (control FAILED)"
        except AssertionError:
            out[name] = "REFUSED"
    out["K2_pole"] = T.finite_singular([R - 2, R - 3])
    out["K2_apparent"] = T.finite_singular([R - 5, R - 3])
    rec, _ = T.analyse("K3", "ctrl", "", lambda: Q.exchange_closed_form(
        F(1), (1, 0, 0), (1, 0, 0), F(1), (1, 0, 0), (1, 0, 0), R, tau_max=6))
    out["K3_exchange_lnR0"] = rec["f_has_lnR_at_0"]
    rec, _ = T.analyse("K4", "ctrl", "", lambda: sp.exp(-R) * sp.log(R) + sp.exp(-3 * R) / R)
    out["K4_planted_lnR0"] = rec["f_has_lnR_at_0"]
    out["pass"] = (out["K1_dilog"] == out["K1_E1sq"] == "REFUSED"
                   and out["K2_pole"][0]["apparent"] is False
                   and out["K2_apparent"][0]["apparent"] is True
                   and out["K3_exchange_lnR0"] and out["K4_planted_lnR0"])
    return out


def all_apparent(fin):
    for f in fin:
        if "apparent" in f:
            if not f["apparent"]:
                return False
        else:
            for exps in f["exponents_numeric"]:
                vals = [complex(sp.N(sp.sympify(x))) for x in exps]
                ints = [round(v.real) for v in vals]
                if any(abs(v - k) > 1e-8 for v, k in zip(vals, ints)):
                    return False
                if len(set(ints)) != len(ints) or min(ints) < 0:
                    return False
    return True


def main():
    only = set(sys.argv[sys.argv.index("--only") + 1:]) if "--only" in sys.argv else None
    out_path = ROOT / "debug" / "data" / "row8_panel.json"
    out = {"controls": None, "cases": []}
    if only is None:
        out["controls"] = controls()
        print("controls", out["controls"]["pass"], flush=True)
    for cid, row, label, fn, num in panel():
        if only and cid not in only:
            continue
        t0 = time.time()
        expr = fn()
        if sp.expand(expr) == 0:
            out["cases"].append({"id": cid, "row": row, "label": label, "zero": True})
            print(cid, "ZERO", flush=True)
            continue
        try:
            rec, expr = T.analyse(label, row, cid, lambda e=expr: e)
        except AssertionError as exc:
            rec = {"refused": str(exc)[:300]}
        rec.update({"id": cid, "row": row, "label": label})
        if "refused" not in rec:
            kinds = {a.split("|")[1] for a in rec["atoms"]}
            rec["kinds"] = sorted(kinds)
            rec["has_gamma"] = any("EulerGamma" in k for k in rec["constants"].values()) or \
                any("EulerGamma" in str(v) for v in rec["constants"])
            rec["all_apparent"] = all_apparent(rec["finite_singular"])
            rec["n_rates"] = len(set(rec["atom_rates"]))
        if num is not None and num[2] != "T":
            oa, ob, kind = num
            ex = float(sp.re(sp.N(expr.subs(R, RV), 30)))
            nu = E.numeric(oa, ob, float(RV), kind)
            rec["numeric_check_R7_5"] = {"exact": ex, "quad": nu, "absdiff": abs(ex - nu)}
        if num is not None and num[2] == "T":
            rec["T_routes_equal"] = bool(sp.simplify(E.T_lap(num[0], num[1]) - E.T_eig(num[0], num[1])) == 0)
        rec["wall_s"] = round(time.time() - t0, 1)
        out["cases"].append(rec)
        print(cid, row, rec.get("kinds"), "ord", rec.get("order_QR"), "rates", rec.get("n_rates"),
              "lnR0", rec.get("f_has_lnR_at_0"), "fuchs", rec.get("R0_fuchsian"),
              "app", rec.get("all_apparent"), "sym", rec.get("verified_symbolic"),
              rec.get("numeric_check_R7_5", {}).get("absdiff"), rec.get("T_routes_equal", ""),
              rec.get("refused", ""), f"{rec['wall_s']}s", flush=True)
        if only is None:
            out_path.write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    if only is None:
        out_path.write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
        print("wrote", out_path)


if __name__ == "__main__":
    main()
