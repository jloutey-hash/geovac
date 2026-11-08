"""Step 0 of the #3 wiring sprint: phase-by-phase cost of `recondition_energy`.

WHY this exists rather than trusting the memo.  `debug/sprint_direct_build_memo.md`
records the split as one_body_mp 9.6 s / _factored_cob 7.4 s / vee_mp 3.7 s -- but
that is measured at (3,3,1) ONLY, and the phases scale differently in N (vee_mp's
X-table grows with p_max/l_neumann, the cob with N^2*(Nr+Na), the normalized solve
with N^2 mpf divisions).  Before claiming the direct one-body engine takes "the
whole calculation to seconds" at the (5,5)+delta headline, measure which phases
actually dominate at more than one truncation.

Two costs the memo does not name are instrumented here: the O(N^2) mpf
`H = H1 + V + Sf*S` assembly loop, and `_normalized_solve`'s 2N^2 mpf divisions
(7.6M of them at N = 1944).

The projection printed at the end is what the wiring would leave behind:
  removed   = one_body_mp + H_assembly_mpf + ONE of the two cob calls
              + (normalized_solve, which becomes float64)
  remaining = vee_mp + ONE cob call (on V alone) + transforms + a float64 solve
because the change of basis is LINEAR, so H_o = H1_o + cob(V) + Sf*S_o.

Usage:  python debug/direct_wire_baseline.py [j_max l_max mu_max [basis]]
        python debug/direct_wire_baseline.py            # (3,3,1) then (4,4,2)
"""
from __future__ import annotations

import sys
import time
from typing import Dict, Tuple

import mpmath as mp
import numpy as np

from geovac import neumann_vee_general_m as ngm
from geovac import prolate_recondition as pr


def baseline(j_max: int, l_max: int, mu_max: int, alpha: float = 1.0,
             basis: str = "laguerre_legendre", R: float = pr.R_DEFAULT,
             l_neumann: int = 0, dps: int = pr.DEFAULT_DPS
             ) -> Tuple[int, float, float, Dict[str, float]]:
    """Replicates recondition_energy's body with a timer on every phase."""
    t: Dict[str, float] = {}
    with mp.workdps(dps):
        idx = pr._product_index(j_max, l_max, mu_max)
        fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        N = len(fns)
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1)
                  if (l + m) % 2 == 0])
        Nmu = mu_max + 1
        if l_neumann <= 0:
            l_neumann = 2 * l_max + 4 * mu_max + 10
        n_mom = 4 * j_max + 4 * (mu_max + 1) + 4 * l_max + 16
        A = ngm._mono_moments(2.0 * alpha, n_mom)

        t0 = time.time()
        S, H1 = pr.one_body_mp(fns, alpha, R, A)
        t["one_body_mp"] = time.time() - t0
        print(f"    one_body_mp      {t['one_body_mp']:8.1f}s", flush=True)

        t0 = time.time()
        V = pr.vee_mp(fns, alpha, R, l_neumann)
        t["vee_mp"] = time.time() - t0
        print(f"    vee_mp           {t['vee_mp']:8.1f}s", flush=True)

        t0 = time.time()
        Sf = 1.0 / mp.mpf(R)
        H = np.empty((N, N), object)
        for i in range(N):
            for j in range(N):
                H[i, j] = H1[i, j] + V[i, j] + Sf * S[i, j]
        t["H_assembly_mpf"] = time.time() - t0
        print(f"    H_assembly_mpf   {t['H_assembly_mpf']:8.1f}s", flush=True)

        t0 = time.time()
        Tr, Ta = pr._transforms_per_mu(basis, j_max, l_max, mu_max, alpha)
        t["transforms"] = time.time() - t0
        print(f"    transforms       {t['transforms']:8.1f}s", flush=True)

        t0 = time.time()
        S_o = pr._factored_cob(S, Nmu, Nr, Na, Tr, Ta)
        t["cob_S"] = time.time() - t0
        print(f"    cob_S            {t['cob_S']:8.1f}s", flush=True)

        t0 = time.time()
        H_o = pr._factored_cob(H, Nmu, Nr, Na, Tr, Ta)
        t["cob_H"] = time.time() - t0
        print(f"    cob_H            {t['cob_H']:8.1f}s", flush=True)

        t0 = time.time()
        E, cond_norm, nk, _sweep = pr._normalized_solve(S_o, H_o)
        t["normalized_solve"] = time.time() - t0
        print(f"    normalized_solve {t['normalized_solve']:8.1f}s", flush=True)

    de = 100.0 * (-1.0 - E) / pr.DE_EXACT
    t["_cond_norm"] = cond_norm
    t["_n_kept"] = float(nk)
    return N, float(E), de, t


def report(j_max: int, l_max: int, mu_max: int, basis: str) -> None:
    print(f"\n=== ({j_max},{l_max}) mu<={mu_max}  basis={basis} ===", flush=True)
    t_all = time.time()
    N, E, de, t = baseline(j_max, l_max, mu_max, basis=basis)
    wall = time.time() - t_all
    phases = [k for k in t if not k.startswith("_")]
    total = sum(t[k] for k in phases)

    print(f"\n  N={N}  E={E:.7f}  D_e%={de:.3f}  "
          f"cond(norm)={t['_cond_norm']:.1e}  kept={int(t['_n_kept'])}")
    print(f"  wall {wall:.1f}s (phases sum {total:.1f}s)\n")
    print("  phase                 seconds    % of total")
    for k in phases:
        print(f"    {k:18s} {t[k]:8.1f}     {100 * t[k] / total:5.1f}")

    # What the wiring removes: the mpf one-body, the mpf H assembly, ONE of the
    # two cob calls (H_o = H1_o + cob(V) + Sf*S_o needs cob on V alone), and the
    # mpf normalized solve (becomes float64).
    removed = t["one_body_mp"] + t["H_assembly_mpf"] + t["cob_S"] + t["normalized_solve"]
    remaining = t["vee_mp"] + t["cob_H"] + t["transforms"]
    print(f"\n  PROJECTION for the direct-engine wiring:")
    print(f"    removed   {removed:8.1f}s  ({100 * removed / total:.1f}%)"
          f"   [one_body_mp + H_assembly + 1 cob + mpf solve]")
    print(f"    remaining {remaining:8.1f}s  ({100 * remaining / total:.1f}%)"
          f"   [vee_mp + 1 cob(V) + transforms] + float64 build/solve")
    if remaining > 0:
        print(f"    implied speedup  ~{total / remaining:.1f}x"
              f"   (NOT counting the float64 build/solve cost, ~0.2s + ~1s)")


if __name__ == "__main__":
    if len(sys.argv) >= 4:
        b = sys.argv[4] if len(sys.argv) >= 5 else "laguerre_legendre"
        report(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), b)
    else:
        for (j, l, mu) in [(3, 3, 1), (4, 4, 2)]:
            report(j, l, mu, "laguerre_legendre")
