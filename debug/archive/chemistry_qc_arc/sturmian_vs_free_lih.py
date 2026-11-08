"""LiH: does the Sturmian shared-exponent constraint cost accuracy per function?

Follow-on to the H2 test (debug/sturmian_vs_gaussian_memo.md), which found NO
per-basis-function advantage for Coulomb Sturmians over contracted Gaussians:
the gap decayed monotonically and inverted by M=18. This adds the four-electron,
HETERONUCLEAR case, where the structural issue is sharper -- Papers 8-9 prove no
single shared exponent can serve Z=3 and Z=1, so a molecular Sturmian set needs
one k PER CENTER.

WHAT IS MEASURED, AND WHY IT DIFFERS FROM THE H2 TEST. There is no LiH Gaussian
reference anywhere in the corpus, and hardcoded Li basis-set data could not be
validated against any energy I trust. Rather than run an un-gateable comparison,
the primary test here is INTERNAL and needs no external basis data:

  (A) Sturmian     -- shells share ONE exponent per center (k_Li, k_H): 2 params
  (B) free-zeta STO -- the SAME shells and the SAME shapes, but every shell's
                       zeta optimized independently: many params

Both are rendered by the same Gaussian fitter and fed to the same integral
engine, so the ONLY difference is whether the exponents are tied together.
This isolates the defining Sturmian property -- the shared exponent -- and asks
what it costs. That is the structural question underneath "replace the basis
wholesale":

  * A ~ B  => the constraint is FREE. Sturmians buy a large parameter reduction
             for nothing, which is a genuine structural economy and would also
             explain the H2 parity result.
  * A << B => the constraint COSTS real accuracy, which is worse news than H2
             gave: the shared exponent would be a liability, not an economy.

STO-3G is included as an indicative third row. Its Li data is validated
self-containedly (below), but there is no LiH energy reference, so it is
reported, never gated on.

PRE-REGISTERED GATES.

  V1 basis data. STO-3G is BY CONSTRUCTION a 3-Gaussian least-squares fit to a
     Slater orbital of specified zeta (Li: zeta_1s = 2.69, zeta_2s = zeta_2p =
     0.75; H: zeta = 1.24). So each hardcoded contraction must overlap its
     Slater with the known STO-3G fit quality, > 0.99. This validates the
     hardcoded numbers WITHOUT needing an energy reference.
  V2 frozen core. Dense 4-electron FCI caps near M=10 (C(20,4)=4845; M=12 is a
     900 MB matrix). Beyond that the Li 1s is frozen -- but by ORBITAL ENERGY,
     not by basis-function index. Freezing by index is what broke the NaH
     ladder (debug/noci_nah_basis_extension_memo.md): an added function sharing
     symmetry with a core function mixes into it under Loewdin, and freezing
     the index then removes bonding flexibility. Here h is diagonalized first
     and the lowest eigen-orbital is frozen. Gate: frozen-core and all-electron
     D must track to < 2 mHa at every M where both are computable.
  V3 conditioning: smallest overlap eigenvalue > 1e-6 for every basis.

  DECISION (on the A-vs-B comparison, at matched M, excluding any row where the
  two families are the same function):
    FREE   |A - B| < 1.6 mHa (chemical accuracy) at every matched M
    COSTLY  B lower than A by > 1.6 mHa at two or more matched M
    MIXED   otherwise

Exploratory. No paper claim.
"""

from __future__ import annotations

import math
import os
import sys
from fractions import Fraction
from typing import Dict, List, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.noci_engine import (
    BasisFn, fci_ground, integral_set_md, lowdin_orbitals, overlap_md,
    transform_integrals,
)
from sturmian_vs_gaussian_convergence import _shape, _sto_norm, sturmian_fn
from geovac.two_center_eri import radial_poly

R_LIH = 3.015           # a0, the corpus's LiH R_eq
P_LMN = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

# STO-3G contractions (coefficients on NORMALIZED primitives) + the Slater
# zeta each one is a fit OF -- the pair is what makes V1 possible.
STO3G = {
    "Li_1s": (2.69, [16.1195750, 2.9362007, 0.7946505],
                    [0.15432897, 0.53532814, 0.44463454], 0),
    "Li_2s": (0.75, [0.6362897, 0.1478601, 0.0480887],
                    [-0.09996723, 0.39951283, 0.70011547], 0),
    "Li_2p": (0.75, [0.6362897, 0.1478601, 0.0480887],
                    [0.15591627, 0.60768372, 0.39195739], 1),
    "H_1s":  (1.24, [3.42525091, 0.62391373, 0.16885540],
                    [0.15432897, 0.53532814, 0.44463454], 0),
}
STO3G_NR = {"Li_1s": 1, "Li_2s": 2, "Li_2p": 2, "H_1s": 1}   # Slater index


def sto_fn(center, l: int, n_r: int, zeta: float, lmn) -> BasisFn:
    """A single Slater function of index n_r, ang.mom. l, exponent zeta."""
    a, d, _ = _shape(l, n_r)
    return BasisFn(center, lmn, a * zeta ** 2, d)


# ------------------------------------------------------------------ shells
# (center_tag, n, l) in ascending order.
LADDER = [
    [("Li", 1, 0), ("Li", 2, 0), ("H", 1, 0)],                        # M=3
    [("Li", 1, 0), ("Li", 2, 0), ("Li", 2, 1), ("H", 1, 0)],          # M=6
    [("Li", 1, 0), ("Li", 2, 0), ("Li", 2, 1), ("H", 1, 0),
     ("H", 2, 0)],                                                     # M=7
    [("Li", 1, 0), ("Li", 2, 0), ("Li", 2, 1), ("Li", 3, 0),
     ("H", 1, 0), ("H", 2, 0)],                                        # M=8
    [("Li", 1, 0), ("Li", 2, 0), ("Li", 2, 1), ("Li", 3, 0),
     ("Li", 3, 1), ("H", 1, 0)],                                       # M=10
]


def _lmns(l):
    return [(0, 0, 0)] if l == 0 else P_LMN


def build_sturmian(shells, centers, k_li, k_h) -> List[BasisFn]:
    out = []
    for tag, n, l in shells:
        k = k_li if tag == "Li" else k_h
        for lmn in _lmns(l):
            out.append(sturmian_fn(centers[tag], n, l, lmn, k))
    return out


def build_free(shells, centers, zetas: Dict[Tuple[str, int, int], float]):
    out = []
    for tag, n, l in shells:
        z = zetas[(tag, n, l)]
        for lmn in _lmns(l):
            out.append(sto_fn(centers[tag], l, n, z, lmn))
    return out


def build_sto3g(centers) -> List[BasisFn]:
    out = []
    for key in ("Li_1s", "Li_2s", "Li_2p", "H_1s"):
        z, a, d, l = STO3G[key]
        c = centers["Li" if key.startswith("Li") else "H"]
        for lmn in _lmns(l):
            out.append(BasisFn(c, lmn, np.asarray(a, float),
                               np.asarray(d, float)))
    return out


# ------------------------------------------------------------------ energies

def _integrals(orbs, centers):
    nuc = [(centers["Li"], 3.0), (centers["H"], 1.0)]
    s, h, g = integral_set_md(orbs, nuc)
    return s, h, g


def energy(orbs, centers, n_core: int = 0) -> Tuple[float, float]:
    """4-electron FCI + V_nn. n_core>0 freezes by ORBITAL ENERGY (see V2)."""
    s, h, g = _integrals(orbs, centers)
    s_min = float(np.linalg.eigvalsh(s).min())
    x = lowdin_orbitals(s)
    ht, gt = transform_integrals(x, h, g)

    if n_core:
        # rotate into the eigenbasis of the one-electron Hamiltonian, so the
        # frozen orbitals are the ENERGETICALLY lowest, not the first-indexed.
        w, u = np.linalg.eigh(ht)
        ht = u.T @ ht @ u
        gt = np.einsum("pi,qj,rk,sl,pqrs->ijkl", u, u, u, u, gt,
                       optimize=True)
        c, a = slice(0, n_core), slice(n_core, ht.shape[0])
        e_core = 2.0 * np.trace(ht[c, c])
        gc = gt[c, c, c, c]
        e_core += 2.0 * np.einsum("iijj->", gc) - np.einsum("ijji->", gc)
        h_eff = (ht[a, a] + 2.0 * np.einsum("pqii->pq", gt[a, a, c, c])
                 - np.einsum("piiq->pq", gt[a, c, c, a]))
        e = e_core + fci_ground(h_eff, gt[a, a, a, a], 4 - 2 * n_core)
    else:
        e = fci_ground(ht, gt, 4)
    return e + 3.0 / R_LIH, s_min


# ------------------------------------------------------------------ gates

def gate_v1() -> bool:
    print("\n[V1] STO-3G data self-validated: each contraction vs its Slater")
    ok = True
    o = np.zeros(3)
    for key, (z, a, d, l) in STO3G.items():
        lmn = (0, 0, 0) if l == 0 else (0, 0, 1)
        gau = BasisFn(o, lmn, np.asarray(a, float), np.asarray(d, float))
        sla = sto_fn(o, l, STO3G_NR[key], z, lmn)
        ov = overlap_md(gau, sla) / math.sqrt(
            overlap_md(gau, gau) * overlap_md(sla, sla))
        good = abs(ov) > 0.99
        ok = ok and good
        print(f"     {key:<6} zeta={z:5.2f}  <STO-3G|Slater> = {abs(ov):.5f} "
              f"{'ok' if good else 'FAIL'}")
    print(f"     V1: {'PASS' if ok else 'FAIL'}")
    return ok


def main() -> None:
    print("=== LiH: what does the Sturmian shared-exponent constraint cost? ===")
    print(f"    R = {R_LIH} a0, 4 electrons, FCI")
    centers = {"Li": np.array([0.0, 0.0, 0.0]),
               "H": np.array([0.0, 0.0, R_LIH])}

    if not gate_v1():
        print("\n*** V1 FAILED -- basis data wrong, comparison void. ***")
        return

    # ---- STO-3G indicative row ------------------------------------------
    orbs = build_sto3g(centers)
    e_sto3g, sm = energy(orbs, centers)
    print(f"\n[indicative] LiH STO-3G FCI (M={len(orbs)}): {e_sto3g:.6f} Ha "
          f"(s_min {sm:.2e})   -- no reference available, not gated")

    # ---- V2: frozen core by orbital energy ------------------------------
    print("\n[V2] frozen core (by ORBITAL ENERGY) vs all-electron")
    v2 = True
    for shells in LADDER[:3]:
        ob = build_sturmian(shells, centers, 2.7, 1.0)
        ea, _ = energy(ob, centers, 0)
        ef, _ = energy(ob, centers, 1)
        d = abs(ea - ef)
        v2 = v2 and d < 2e-3
        print(f"     M={len(ob):<3} all-elec {ea:12.6f}  frozen {ef:12.6f}  "
              f"delta {d*1000:7.3f} mHa  {'ok' if d < 2e-3 else 'FAIL'}")
    print(f"     V2: {'PASS' if v2 else 'FAIL'}")

    # ---- the comparison --------------------------------------------------
    print("\n[A] Sturmian: one shared exponent PER CENTER (k_Li, k_H)")
    print("[B] free-zeta STO: same shells and shapes, every zeta independent")
    print(f"\n  {'M':>3} {'A Sturmian':>12} {'B free-zeta':>12} "
          f"{'B-A (mHa)':>10} {'k_Li':>5} {'k_H':>5} {'s_min':>9}")

    rows = []
    for shells in LADDER:
        M = sum(len(_lmns(l)) for _, _, l in shells)
        n_core = 0                      # V2 FAILED: never freeze

        # A: 2-parameter scan
        bestA = None
        for k_li in (1.8, 2.1, 2.4, 2.7, 3.0, 3.3, 3.6):
            for k_h in (0.35, 0.5, 0.65, 0.8, 1.0, 1.2, 1.4):
                ob = build_sturmian(shells, centers, k_li, k_h)
                e, sm = energy(ob, centers, n_core)
                if bestA is None or e < bestA[0]:
                    bestA = (e, k_li, k_h, sm)

        # B: coordinate descent on every shell's zeta, seeded from A
        z = {}
        for tag, n, l in shells:
            z[(tag, n, l)] = (bestA[1] if tag == "Li" else bestA[2])
        eB, smB = energy(build_free(shells, centers, z), centers, n_core)
        for _sweep in range(3):
            for key in list(z):
                base = z[key]
                for f in (0.6, 0.75, 0.9, 1.0, 1.15, 1.35, 1.6):
                    z[key] = base * f
                    try:
                        e, sm = energy(build_free(shells, centers, z),
                                       centers, n_core)
                    except Exception:
                        continue
                    if e < eB:
                        eB, smB, base = e, sm, z[key]
                z[key] = base

        diff = (eB - bestA[0]) * 1000.0
        rows.append((M, bestA[0], eB, diff))
        print(f"  {M:>3} {bestA[0]:12.6f} {eB:12.6f} {diff:10.3f} "
              f"{bestA[1]:5.2f} {bestA[2]:5.2f} {min(bestA[3], smB):9.2e}")

    print("\n=== verdict ===")
    costly = [r for r in rows if r[3] < -1.6]
    worst = min(r[3] for r in rows)
    print(f"  free-zeta beats Sturmian by > 1.6 mHa at {len(costly)} of "
          f"{len(rows)} sizes; largest gap {abs(worst):.3f} mHa")
    if len(costly) >= 2:
        print("  VERDICT: COSTLY -- the shared-exponent constraint gives up")
        print("           real accuracy per function.")
    elif all(abs(r[3]) < 1.6 for r in rows):
        print("  VERDICT: FREE -- the constraint costs less than chemical")
        print("           accuracy at every size, i.e. Sturmians buy a large")
        print("           parameter reduction for nothing.")
    else:
        print("  VERDICT: MIXED.")


if __name__ == "__main__":
    main()
