#!/usr/bin/env python
"""Re-measure Paper 19 tab:balanced_census from the live builders.

Every cell is pair-diagonal vintage (pre-2026-08-29), and the balanced cells
are additionally at whatever geometry the old R=3.015 default supplied.  The
table is re-measured rather than mapped, because Paper 19's BeH2/H2O balanced
counts (2,652 / 5,798) do not correspond to Paper 14 tab:balanced's
(8,867 / 19,741) by any convention -- so a substitution would be a guess.

Conventions reported for every cell so the table's own can be matched:
  Pauli  -- non-identity and identity-inclusive
  lambda -- non-identity and identity-inclusive, composed both with PK
            in-Hamiltonian and electronic-only (PK partitioned out)

Usage:  python debug/qa/remeasure_paper19_census.py
"""
from __future__ import annotations

import sys
import time
import warnings

from geovac import molecular_spec as MS
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.composed_qubit import build_composed_hamiltonian

ROWS = [
    ("LiH",  MS.lih_spec,  2, 30),
    ("BeH2", MS.beh2_spec, 3, 50),
    ("H2O",  MS.h2o_spec,  5, 70),
]


def norms(op):
    tot = sum(abs(complex(c)) for c in op.terms.values())
    ni = sum(abs(complex(c)) for t, c in op.terms.items() if t)
    return tot, ni


def main() -> int:
    warnings.simplefilter("ignore")
    print("Paper 19 tab:balanced_census -- re-measured "
          "(each molecule at its OWN geometry)\n")
    for name, factory, blocks, Q in ROWS:
        t0 = time.perf_counter()
        bal = build_balanced_hamiltonian(factory(max_n=2))
        comp_pk = build_composed_hamiltonian(
            factory(max_n=2), pk_in_hamiltonian=True, verbose=False)
        comp_el = build_composed_hamiltonian(
            factory(max_n=2), pk_in_hamiltonian=False, verbose=False)

        bal_tot, bal_ni = norms(bal["qubit_op"])
        pk_tot, pk_ni = norms(comp_pk["qubit_op"])
        el_tot, el_ni = norms(comp_el["qubit_op"])

        n_comp = comp_pk["N_pauli"]
        n_bal = bal["N_pauli"]
        dt = time.perf_counter() - t0

        print(f"{name}  blocks={blocks}  Q={Q}   R={getattr(factory(max_n=2), 'R', None)}"
              f"   [{dt:.0f}s]")
        print(f"   composed Pauli : {n_comp} non-id / {n_comp + 1} with id")
        print(f"   balanced Pauli : {n_bal} (builder returns identity-inclusive)")
        print(f"   Pauli ratio    : {n_bal / (n_comp + 1):.2f}x  "
              f"(identity-inclusive both)")
        print(f"   lambda comp PK : {pk_ni:.1f} non-id / {pk_tot:.1f} with id")
        print(f"   lambda comp el : {el_ni:.1f} non-id / {el_tot:.1f} with id")
        print(f"   lambda balanced: {bal_ni:.1f} non-id / {bal_tot:.1f} with id")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
