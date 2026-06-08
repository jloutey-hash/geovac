"""Mechanism diagnostic: which ERI terms break kappa-parity?

We have established that P_kappa does NOT commute with H_rel (residuals
10^-2 to 10^-1 on LiH/BeH/CaH at max_n=2). This script identifies which
ERI 4-tuples carry odd kappa-sign-flip and quantifies their weight,
producing a clean structural reason for the failure.

A two-body term a^dag b^dag d c shifts N_{kappa<0} by
    Delta = [kappa_a<0] + [kappa_b<0] - [kappa_c<0] - [kappa_d<0]
in {-2,-1,0,+1,+2}. For [H, P_kappa] = 0 we need Delta even mod 2 for
every nonzero ERI element (PARITY-EVEN terms commute with the Z-string).
Delta odd (= +-1) terms anticommute with P_kappa and produce the
nonzero commutator.

We also verify that h1 (one-body) is kappa-sign diagonal (no
odd-parity terms).
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Tuple

from geovac.composed_qubit_relativistic import build_composed_hamiltonian_relativistic
from geovac.molecular_spec import lih_spec_relativistic


def main() -> None:
    spec = lih_spec_relativistic(max_n=2, include_pk=True)
    result = build_composed_hamiltonian_relativistic(
        spec, alpha_num=7.2973525693e-3, pk_in_hamiltonian=True,
        verbose=False, include_breit=False,
    )

    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sprint_kappa_parity_diagnostic import (
        enumerate_relativistic_orbital_table,
    )
    table, Q = enumerate_relativistic_orbital_table(spec)

    def kappa_of(q: int) -> int:
        return table[q][2]

    eri = result["eri_sparse"]
    print(f"LiH_rel max_n=2: Q={Q}, |eri_sparse|={len(eri)}")

    # Classify ERI terms by Delta parity
    counts = defaultdict(int)
    weight_sums = defaultdict(float)
    for (a, b, c, d), val in eri.items():
        delta = (int(kappa_of(a) < 0) + int(kappa_of(b) < 0)
                 - int(kappa_of(c) < 0) - int(kappa_of(d) < 0))
        counts[delta] += 1
        weight_sums[delta] += abs(val)

    print("\nERI tuples by kappa-sign-flip count Delta:")
    print(f"  {'Delta':>6} {'count':>10} {'sum |val|':>14}")
    for d in sorted(counts):
        print(f"  {d:>6} {counts[d]:>10} {weight_sums[d]:>14.6e}")

    odd_count = sum(counts[d] for d in counts if d % 2 != 0)
    odd_weight = sum(weight_sums[d] for d in counts if d % 2 != 0)
    print(f"\n  Odd-Delta tuples (anticommute with P_kappa):  "
          f"{odd_count} / {len(eri)}  weight={odd_weight:.4e}")

    # Sample a few odd-Delta examples
    print("\n  Sample odd-Delta ERI terms:")
    shown = 0
    for (a, b, c, d), val in eri.items():
        delta = (int(kappa_of(a) < 0) + int(kappa_of(b) < 0)
                 - int(kappa_of(c) < 0) - int(kappa_of(d) < 0))
        if delta % 2 != 0:
            la = table[a]
            lb = table[b]
            lc = table[c]
            ld = table[d]
            print(f"    (a={a},b={b},c={c},d={d}) Delta={delta:+d} val={val:+.4e}")
            print(f"      a: n={la[1]} kappa={la[2]} 2mj={la[3]}")
            print(f"      b: n={lb[1]} kappa={lb[2]} 2mj={lb[3]}")
            print(f"      c: n={lc[1]} kappa={lc[2]} 2mj={lc[3]}")
            print(f"      d: n={ld[1]} kappa={ld[2]} 2mj={ld[3]}")
            shown += 1
            if shown >= 3:
                break

    # Also verify h1 has no off-diagonal odd-Delta entries
    h1 = result["h1_diag"]  # diagonal only here; check via fermion_op terms
    fop = result["fermion_op"]
    h1_terms = [(t, c) for t, c in fop.terms.items()
                if len(t) == 2 and t[0][1] == 1 and t[1][1] == 0]
    h1_odd_count = 0
    h1_odd_weight = 0.0
    for term, coef in h1_terms:
        p, q = term[0][0], term[1][0]
        delta = int(kappa_of(p) < 0) - int(kappa_of(q) < 0)
        if delta % 2 != 0:
            h1_odd_count += 1
            h1_odd_weight += abs(coef)
    print(f"\n  h1 odd-Delta terms:  {h1_odd_count}  weight={h1_odd_weight:.4e}")

    print("\nVERDICT: nonzero odd-Delta ERI weight = structural reason P_kappa")
    print("does not commute. The jj-coupled X_k Gaunt selection allows")
    print("processes p_{3/2} -> p_{1/2} at k=2 (parity-allowed, kappa-flip).")


if __name__ == "__main__":
    main()
