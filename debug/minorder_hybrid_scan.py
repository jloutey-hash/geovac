"""Hybrid (AA|AB) with an L >= 1 same-centre pair: order in R against the exponent ratio, looking for
exceptional ratios (rate coincidences) where the order drops.  Two families, Z_B scanned:

  fam1  (2p0 2p0 | 1s 1s_B),  Z_A = 3:  pair rate a_p = 3, a_c = 3, a_d = Z_B
  fam2  (2p1 2p1 | 2p1 2p1_B), Z_A = 2:  pair rate a_p = 2, a_c = 1, a_d = Z_B / 2

For each Z_B: order r = smallest r with generic rank over Q(R) of the split vectors of f..f^(r) equal to r (an
order-r operator exists, exact), then the TRUE-constant lower bound at r (minorder_core.analyse).
Writes debug/data/minorder_hybrid_scan.json.   Usage: python debug/minorder_hybrid_scan.py [fam1|fam2]
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "debug"))

import minorder_core as MC  # noqa: E402
from geovac.two_center_eri import hybrid_closed_form  # noqa: E402

F = Fraction
RO, FR = MC.RO, MC.FR
S1, P0, P1 = (1, 0, 0), (2, 1, 0), (2, 1, 1)

FAMS = {
    "fam1": (lambda zb: hybrid_closed_form(F(3), P0, P0, S1, F(zb), S1),
             [F(1), F(3, 2), F(2), F(5, 2), F(3), F(4), F(5), F(6), F(7), F(9), F(12)],
             lambda zb: (F(3), F(3), F(zb))),
    "fam2": (lambda zb: hybrid_closed_form(F(2), P1, P1, P1, F(zb), P1),
             [F(1), F(2), F(3), F(4), F(6), F(8), F(10)],
             lambda zb: (F(2), F(1), F(zb) / 2)),
}


def generic_order(f, rmax=8):
    atoms, reps = RO.decompose(f)
    ds = [RO.to_field(atoms, reps, "split")]
    for r in range(1, rmax + 1):
        ds.append(RO.d_atoms(ds[-1]))
        keys = sorted(set().union(*[set(D) for D in ds]), key=str)
        M = DomainMatrix([[D.get(kk, FR.zero) for D in ds] for kk in keys], (len(keys), r + 1), FR)
        if M.rank() == r:
            return r
    return None


def main():
    which = sys.argv[1:] or list(FAMS)
    path = ROOT / "debug" / "data" / "minorder_hybrid_scan.json"
    out = json.loads(path.read_text()) if path.exists() else {}
    for fam in which:
        fn, zbs, rates = FAMS[fam]
        out.setdefault(fam, {})
        for zb in zbs:
            t0 = time.time()
            f = fn(zb)
            r = generic_order(f)
            rec = MC.analyse(f, r, upper=False)
            ap, ac, ad = rates(zb)
            nr = len(RO.atom_rates(f))
            row = {"a_p": str(ap), "a_c": str(ac), "a_d": str(ad), "order_generic_constants": r,
                   "n_rates": nr, "rates": sorted(str(x) for x in RO.atom_rates(f)),
                   "order_minus_rates": r - nr, "true_lower_proved": rec["proved_lower"],
                   "f_atoms": rec["f_atoms"], "n_closure": rec["n_closure"],
                   "constants": list(rec["constants"].values()), "seconds": round(time.time() - t0, 1)}
            out[fam][str(zb)] = row
            print(f"{fam} ZB={str(zb):5s} a_d={str(ad):4s} order={r} rates={nr} {row['rates']} "
                  f"order-rates={r - nr} true_lower={rec['proved_lower']} closure={rec['n_closure']} "
                  f"kinds={sorted({a.split(',')[1].strip() for a in rec['f_atoms']})} [{row['seconds']} s]",
                  flush=True)
            path.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
