#!/usr/bin/env python
"""Re-measure the balanced lambda column at each molecule's OWN geometry.

The published column (Papers 14 tab:second_row, 20 tab:molecules) was computed
at R = 3.015 -- LiH's bond length -- for every molecule, because
`build_balanced_hamiltonian` defaulted R to that value and ignored the spec.
See debug/qa/balanced_lambda_geometry_finding.md.

This prints the corrected column together with the geometry each value belongs
to, so the free-side input is named rather than inherited.

Usage:  python debug/qa/remeasure_balanced_lambda.py
"""
from __future__ import annotations

import sys
import time
import warnings

from geovac import molecular_spec as MS
from geovac.balanced_coupled import build_balanced_hamiltonian

# molecule -> (spec factory, published lambda at the WRONG R = 3.015)
ROWS = [
    ("NaH",  MS.nah_spec,   20.6),
    ("MgH2", MS.mgh2_spec, 110.5),
    ("HCl",  MS.hcl_spec,  869.7),
    ("H2S",  MS.h2s_spec,  879.6),
    ("PH3",  MS.ph3_spec,  895.3),
    ("SiH4", MS.sih4_spec, 914.1),
    ("KH",   MS.kh_spec,    32.2),
    ("CaH2", MS.cah2_spec, 128.0),
    ("HBr",  MS.hbr_spec,  877.0),
    ("H2Se", MS.h2se_spec, 883.4),
    ("AsH3", MS.ash3_spec, 885.2),
    ("GeH4", MS.geh4_spec, 877.3),
]


def lambda_ni(res) -> float:
    return sum(abs(complex(c)) for t, c in res["qubit_op"].terms.items() if t)


def main() -> int:
    warnings.simplefilter("ignore")
    print(f"{'mol':<6}{'R (bohr)':>10}{'published':>11}{'corrected':>11}"
          f"{'delta':>9}{'delta %':>9}   {'Pauli':>7}{'QWC':>7}")
    print("-" * 72)
    rows = []
    for name, factory, published in ROWS:
        spec = factory(max_n=2)
        R = getattr(spec, "R", None)
        t0 = time.perf_counter()
        res = build_balanced_hamiltonian(factory(max_n=2))   # R from the spec
        dt = time.perf_counter() - t0
        lam = lambda_ni(res)
        d = lam - published
        print(f"{name:<6}{R:>10.3f}{published:>11.1f}{lam:>11.4f}"
              f"{d:>9.3f}{100 * d / published:>8.2f}%   "
              f"{res['N_pauli']:>7}{res['n_qwc_balanced']:>7}   [{dt:.0f}s]")
        rows.append((name, R, published, lam, res["N_pauli"],
                     res["n_qwc_balanced"]))

    print("\n--- LaTeX (lambda column, corrected; R named) ---")
    for name, R, _pub, lam, n_p, n_q in rows:
        print(f"  {name:<6} & {R:.3f} & {n_p:,} & {lam:.1f} & {n_q:,} \\\\"
              .replace(",", "{,}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
