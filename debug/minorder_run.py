"""Exact minimal-order lower bounds with true constants for the Paper 58 two-centre panels.

Usage:  python debug/minorder_run.py GROUP [CASEID ...]
GROUP in: p1s (the v6.7.26 1s panel, 30 cases), r8 (record R8, 48 cases), xeq (equal-exponent exchange, 7 cases).
Writes debug/data/minorder_<GROUP>.json (rewritten after every case).
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

import minorder_core as MC  # noqa: E402
from geovac import qfd_core as Q  # noqa: E402
from geovac import two_center_one_electron as E  # noqa: E402
from geovac.two_center_eri import aabb_closed_form, hybrid_closed_form  # noqa: E402

R = MC.R
F = Fraction
S1 = (1, 0, 0)


def o(s):
    return (int(s[0]), "spd".index(s[1]), int(s[2:]) if len(s) > 2 else 0)


def A(s, Z):
    return ("A", F(Z)) + o(s)


def B(s, Z):
    return ("B", F(Z)) + o(s)


def panel_1s():
    """Same list, same order, as debug/two_centre_r_operators.py build_panel()."""
    qa = lambda z: ("A", F(z), 1)   # noqa: E731
    qb = lambda z: ("B", F(z), 1)   # noqa: E731
    cases = []
    for za, zb in [(1, 1), (1, 2), (2, 1), (F(3, 2), 1)]:
        tag = f"zA={za},zB={zb}"
        cases.append((f"S {tag}", lambda za=za, zb=zb: Q.overlap(qa(za), qb(zb), R)))
        cases.append((f"T {tag}", lambda za=za, zb=zb: Q.kinetic(qa(za), qb(zb), R)))
        cases.append((f"VA(ab) {tag}", lambda za=za, zb=zb: Q._inv_r(qa(za), qb(zb), "A", R)))
        cases.append((f"VB(aa) {tag}", lambda za=za, zb=zb: Q._inv_r(qa(za), qa(za), "B", R)))
        cases.append((f"(AA|BB)1s {tag}", lambda za=za, zb=zb: aabb_closed_form(F(za), S1, S1, F(zb), S1, S1)))
        cases.append((f"(AA|AB)1s {tag}", lambda za=za, zb=zb: hybrid_closed_form(F(za), S1, S1, S1, F(zb), S1)))
    p0 = (2, 1, 0)
    for za, zb, rel in ((3, 3, "eq"), (3, 4, "gt"), (3, 2, "lt")):
        cases.append((f"(2p0 2p0|1s 1s_B) {rel} ZA={za},ZB={zb}",
                      lambda za=za, zb=zb: hybrid_closed_form(F(za), p0, p0, S1, F(zb), S1)))
    for z in (1, F(3, 2), 2):
        cases.append((f"(AB|AB)1s z={z}",
                      lambda z=z: Q.exchange_closed_form(F(z), S1, S1, F(z), S1, S1, R, tau_max=6)))
    stored = json.loads((ROOT / "debug" / "data" / "two_centre_r_operators.json").read_text())["cases"]
    assert len(stored) == len(cases)
    return [(nm, fn, int(st["order_QR"])) for (nm, fn), st in zip(cases, stored)]


ROW2_PAIRS = {"E1": ("2p0", 3, "1s", 1), "E2": ("2p1", 2, "2p1", 2), "E3": ("3d0", 3, "2p0", 2),
              "E4": ("3d2", 3, "3d2", F(3, 2)), "E5": ("3p1", 3, "3d1", 2), "E6": ("2s", 2, "3d0", 3)}
ROW2_DENS = {"F1": ("2p0", "2p0", 3), "F2": ("2p1", "2p1", 2), "F3": ("3d0", "3s", 3),
             "F4": ("3d1", "2p1", 3), "F5": ("3d2", "3d2", 3), "F6": ("2p0", "2s", 2)}
ROW3A = {"C1": ("3d0", "3d0", 3, "1s", "1s", 1), "C2": ("2p1", "2p1", 2, "2p0", "2p0", 2),
         "C3": ("3d1", "3d1", 3, "3d2", "3d2", F(3, 2)), "C4": ("2p0", "3d0", 3, "1s", "2p0", 2),
         "C5": ("3d2", "3d0", 3, "3d0", "3d2", 3), "C6": ("2p1", "2p0", 2, "3d0", "3d1", 3),
         "C7": ("2s", "3p0", 3, "2p0", "3s", 2), "C8": ("3d2", "3d2", 3, "3d2", "3d2", 3)}
ROW3B = {"H1": ("3d0", "3d0", "3d0", 3, "1s", 1), "H2": ("2p1", "2p1", "2p1", 2, "2p1", 2),
         "H3": ("3d0", "3d0", "1s", 3, "1s", 1), "H4": ("2p0", "2p0", "2p0", 2, "3d0", 3),
         "H5": ("3d1", "3d1", "3d1", 3, "2p1", 3), "H6": ("3d2", "3d2", "3s", 3, "2s", 1),
         "H7": ("2p0", "3d0", "2s", 2, "1s", 1), "H8": ("3p1", "3d1", "3d0", 3, "2p0", 4),
         "H9": ("2s", "2s", "3d0", 2, "3d0", 3), "H10": ("3d2", "3d0", "3d-1", 3, "3d1", 3)}


def panel_r8():
    """Same list, same order, as debug/row8_panel.py panel()."""
    cases = []
    for cid, (a, za, b, zb) in ROW2_PAIRS.items():
        oa, ob = A(a, za), B(b, zb)
        for kind, fn in (("S", E.S), ("VA", E.VA), ("VB", E.VB), ("T", E.T_lap)):
            cases.append((f"{cid}{kind}", lambda fn=fn, oa=oa, ob=ob: fn(oa, ob)))
    for cid, (a, a2, z) in ROW2_DENS.items():
        oa, ob = A(a, z), A(a2, z)
        cases.append((cid, lambda oa=oa, ob=ob: E.VB(oa, ob)))
    for cid, (a, b, za, c, d, zb) in ROW3A.items():
        cases.append((cid, lambda a=a, b=b, za=za, c=c, d=d, zb=zb: aabb_closed_form(
            F(za), o(a), o(b), F(zb), o(c), o(d))))
    for cid, (a, b, c, za, d, zb) in ROW3B.items():
        cases.append((cid, lambda a=a, b=b, c=c, za=za, d=d, zb=zb: hybrid_closed_form(
            F(za), o(a), o(b), o(c), F(zb), o(d))))
    stored = json.loads((ROOT / "debug" / "data" / "row8_panel.json").read_text())["cases"]
    assert len(stored) == len(cases)
    for (nm, _), st in zip(cases, stored):
        assert st["params"] == nm, (st["params"], nm)
    return [(nm, fn, int(st["order_QR"])) for (nm, fn), st in zip(cases, stored)]


XEQ = {  # equal exponents a = b = 1; hydrogenic nodeless orbitals, rate Z/n
    "1s": (1, S1, S1, 1, S1, S1),
    "2psig": (2, (2, 1, 0), S1, 1, (2, 1, 0), S1),
    "2ppi": (2, (2, 1, 1), S1, 1, (2, 1, -1), S1),
    "3dsig": (3, (3, 2, 0), S1, 1, (3, 2, 0), S1),
    "3ddel": (3, (3, 2, 2), S1, 1, (3, 2, -2), S1),
    "ppsig": (2, (2, 1, 0), (2, 1, 0), 2, (2, 1, 0), (2, 1, 0)),
    "pppi0": (2, (2, 1, 1), (2, 1, 1), 2, (2, 1, 1), (2, 1, 1)),
}


def panel_xeq():
    def make(spec):
        ZA, oa, ob, ZB, oc, od = spec

        def fn():
            terms = Q.exchange_closed_form(F(ZA), oa, ob, F(ZB), oc, od, R, tau_max=6, return_terms=True)
            assert all(t == 0 for t in terms[-2:]), "tau sum did not stop"
            return sp.Add(*terms)
        return fn
    return [(nm, make(spec), 4) for nm, spec in XEQ.items()]


GROUPS = {"p1s": panel_1s, "r8": panel_r8, "xeq": panel_xeq}


def main():
    group = sys.argv[1]
    only = set(sys.argv[2:])
    out_path = ROOT / "debug" / "data" / f"minorder_{group}.json"
    out = {"group": group, "cases": {}}
    if only and out_path.exists():
        out = json.loads(out_path.read_text())
    for nm, fn, r in GROUPS[group]():
        if only and nm not in only:
            continue
        t0 = time.time()
        f = fn()
        t1 = time.time()
        rec = MC.analyse(f, r)
        # one order below the claim must also hold trivially; one order ABOVE must fail when r = closure size
        rec["build_s"] = round(t1 - t0, 1)
        rec["analyse_s"] = round(time.time() - t1, 1)
        out["cases"][nm] = rec
        print(f"{nm:38s} r={r} closure={rec['n_closure']} rates={rec['n_rates_closure']} "
              f"Etype={rec['n_E_type_atoms']} lower_proved={rec['proved_lower']} "
              f"upper={rec.get('upper_proved')} minor={rec.get('minor_interval')} "
              f"consts={list(rec['constants'].values())} [{rec['build_s']}+{rec['analyse_s']} s]", flush=True)
        out_path.write_text(json.dumps(out, indent=1))
    ok = all(c["proved_lower"] for c in out["cases"].values())
    ex = all(c.get("upper_proved") for c in out["cases"].values())
    print(f"ALL LOWER BOUNDS PROVED: {ok};  ALL UPPER BOUNDS PROVED: {ex}", flush=True)


if __name__ == "__main__":
    main()
