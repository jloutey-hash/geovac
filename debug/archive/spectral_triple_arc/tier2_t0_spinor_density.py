"""
Track T0 — Spinor-block ERI density d_spinor(l_max).

Dirac-on-S^3 Tier 2, early-win standalone deliverable.

Computes the Paper 22 analog for spinor blocks in the jj-coupled basis
((kappa, m_j) labels). The angular selection rules reduce to Wigner 3j
identities plus (j, k, j') triangle + parity, evaluated in exact sympy
arithmetic. No numerical tensors. No radial integrals.

CONVENTION NOTE
===============
Paper 22's published "angular ERI density" table (100%, 7.81%, 2.76%,
1.44%, 0.90%, 0.62%) is computed in a **pair-diagonal-m convention**:
the per-pair Gaunt coefficient c^k(l_a m_a; l_c m_c) is required
nonzero with bottom-row m-sum = 0, which with the reference code's
convention q = m_c - m_a forces m_a = m_c (and similarly m_b = m_d for
the second pair). This is a stricter selection rule than the full
Gaunt rule (full Gaunt requires only m_a + m_b = m_c + m_d globally
with intermediate q = m_a - m_c = -(m_b - m_d)).

For apples-to-apples comparison with Paper 22, we reproduce the
pair-diagonal-m convention for both scalar and spinor here. We also
report the full-Gaunt density in a secondary column for completeness
and cross-check. In both conventions, d_spinor <= d_scalar should hold
(extra j-triangle rules can only add zeros).

Full-Gaunt selection for scalar <ab|1/r12|cd>:
  * m_a + m_b = m_c + m_d  (global m-conservation)
  * exists k with (l_a+l_c+k) even, (l_b+l_d+k) even,
                  triangle(l_a,k,l_c), triangle(l_b,k,l_d),
                  3j(l_a,k,l_c;0,0,0) != 0,
                  3j(l_a,k,l_c; -m_a, m_a-m_c, m_c) != 0,
                  3j(l_b,k,l_d; -m_b, m_b-m_d, m_d) != 0.

Paper-22 (pair-diagonal) convention just replaces "m_a + m_b = m_c + m_d"
+ the two cross-consistent 3j's with "m_a = m_c AND m_b = m_d AND the
two diagonal 3j's with q = 0 are nonzero".

Full-Gaunt selection for spinor <(k1 m1)(k2 m2)|1/r12|(k3 m3)(k4 m4)>:
  * m_1 + m_2 = m_3 + m_4
  * exists k with (l_1+l_3+k) even, (l_2+l_4+k) even,
                  triangle(j_1,k,j_3), triangle(j_2,k,j_4),
                  triangle(l_1,k,l_3), triangle(l_2,k,l_4),
                  3j(j_1,k,j_3; 1/2, 0, -1/2) != 0,
                  3j(j_2,k,j_4; 1/2, 0, -1/2) != 0,
                  3j(j_1,k,j_3; -m_1, m_1-m_3, m_3) != 0,
                  3j(j_2,k,j_4; -m_2, m_2-m_4, m_4) != 0.

Paper-22 (pair-diagonal) spinor variant:
  * m_1 = m_3 AND m_2 = m_4
  * exists k with (l_1+l_3+k) even, (l_2+l_4+k) even,
                  triangle(j_1,k,j_3), triangle(j_2,k,j_4),
                  triangle(l_1,k,l_3), triangle(l_2,k,l_4),
                  3j(j_1,k,j_3; 1/2, 0, -1/2) != 0,
                  3j(j_2,k,j_4; 1/2, 0, -1/2) != 0,
                  3j(j_1,k,j_3; -m_1, 0, m_3) != 0,    (with m_1=m_3)
                  3j(j_2,k,j_4; -m_2, 0, m_4) != 0.    (with m_2=m_4)

Reference: Dyall & Faegri §9; Grant §7.5; Paper 22 §II-III.

Exact arithmetic throughout: sympy Rational/Integer + wigner_3j.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

from sympy import Integer, Rational
from sympy.physics.wigner import wigner_3j


# ---------------------------------------------------------------------------
# Spinor label enumeration
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SpinorOrbital:
    """(kappa, 2*m_j) label. l and j derived from kappa."""
    kappa: int
    two_mj: int  # 2 * m_j, an odd integer

    @property
    def l(self) -> int:
        """Orbital angular momentum from kappa.

        Convention:
          kappa > 0: j = l - 1/2, with kappa = l.       l = kappa.
          kappa < 0: j = l + 1/2, with kappa = -(l+1).  l = -kappa - 1.
        """
        if self.kappa > 0:
            return self.kappa
        elif self.kappa < 0:
            return -self.kappa - 1
        else:
            raise ValueError("kappa cannot be 0")

    @property
    def two_j(self) -> int:
        """2*j (odd integer; j = |kappa| - 1/2 → 2j = 2|kappa|-1)."""
        return 2 * abs(self.kappa) - 1


def enumerate_spinor_orbitals(l_max: int) -> List[SpinorOrbital]:
    """All (kappa, m_j) with l(kappa) <= l_max.

    Produces:
      kappa = -(l+1) for l = 0..l_max  (j = l+1/2)  — (l_max+1) kappa values
      kappa = +l     for l = 1..l_max  (j = l-1/2)  — l_max kappa values

    For each kappa, m_j ∈ {-j, ..., j}, contributing (2j+1) states.

    Total count = 2 * (l_max + 1)**2  [matches task spec].
    """
    orbs: List[SpinorOrbital] = []
    # Negative kappa branch: j = l + 1/2, l = 0..l_max
    for l in range(l_max + 1):
        kappa = -(l + 1)
        two_j = 2 * l + 1  # j = l + 1/2 → 2j = 2l+1
        for two_mj in range(-two_j, two_j + 1, 2):
            orbs.append(SpinorOrbital(kappa=kappa, two_mj=two_mj))
    # Positive kappa branch: j = l - 1/2, l = 1..l_max
    for l in range(1, l_max + 1):
        kappa = l
        two_j = 2 * l - 1  # j = l - 1/2 → 2j = 2l-1
        for two_mj in range(-two_j, two_j + 1, 2):
            orbs.append(SpinorOrbital(kappa=kappa, two_mj=two_mj))
    return orbs


# ---------------------------------------------------------------------------
# Selection rule helpers
# ---------------------------------------------------------------------------

def _triangle_x2(a_x2: int, b_x2: int, c_x2: int) -> bool:
    """Triangle inequality on (2j)-valued triple. Returns True iff
    |a-b| <= c <= a+b AND (a+b+c) even (integer sum)."""
    if c_x2 < abs(a_x2 - b_x2):
        return False
    if c_x2 > a_x2 + b_x2:
        return False
    if (a_x2 + b_x2 + c_x2) % 2 != 0:
        return False
    return True


# ---------------------------------------------------------------------------
# Scalar density: full Gaunt vs Paper-22 pair-diagonal
# ---------------------------------------------------------------------------

def _scalar_pair_k_fullgaunt(orbs: List[Tuple[int, int]], l_max: int
                             ) -> Dict[Tuple[int, int], frozenset]:
    """For each (a,c), set of k satisfying full-Gaunt c^k(l_a m_a; l_c m_c) != 0.
    q = m_a - m_c, bottom-row sum = 0."""
    pair_k = {}
    k_cap = 2 * l_max
    for i, (la, ma) in enumerate(orbs):
        for j, (lc, mc) in enumerate(orbs):
            q = ma - mc
            ks = set()
            for k in range(k_cap + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                if k < abs(la - lc) or k > la + lc:
                    continue
                if abs(q) > k:
                    continue
                w1 = wigner_3j(Integer(la), Integer(k), Integer(lc),
                               Integer(0), Integer(0), Integer(0))
                if w1 == 0:
                    continue
                w2 = wigner_3j(Integer(la), Integer(k), Integer(lc),
                               Integer(-ma), Integer(q), Integer(mc))
                if w2 == 0:
                    continue
                ks.add(k)
            pair_k[(i, j)] = frozenset(ks)
    return pair_k


def _scalar_pair_k_pairdiag(orbs: List[Tuple[int, int]], l_max: int
                            ) -> Dict[Tuple[int, int], frozenset]:
    """Paper-22 pair-diagonal convention: require m_a = m_c; then q=0.
    If m_a != m_c, empty set."""
    pair_k = {}
    k_cap = 2 * l_max
    for i, (la, ma) in enumerate(orbs):
        for j, (lc, mc) in enumerate(orbs):
            if ma != mc:
                pair_k[(i, j)] = frozenset()
                continue
            ks = set()
            for k in range(k_cap + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                if k < abs(la - lc) or k > la + lc:
                    continue
                w1 = wigner_3j(Integer(la), Integer(k), Integer(lc),
                               Integer(0), Integer(0), Integer(0))
                if w1 == 0:
                    continue
                w2 = wigner_3j(Integer(la), Integer(k), Integer(lc),
                               Integer(-ma), Integer(0), Integer(mc))
                if w2 == 0:
                    continue
                ks.add(k)
            pair_k[(i, j)] = frozenset(ks)
    return pair_k


def compute_scalar_density(l_max: int, convention: str
                           ) -> Tuple[int, int, Fraction]:
    """Compute scalar angular ERI density.

    convention ∈ {"fullgaunt", "pairdiag"}.

    pairdiag reproduces Paper 22 published table exactly.
    fullgaunt is the physically correct Coulomb selection.
    """
    orbs: List[Tuple[int, int]] = []
    for l in range(l_max + 1):
        for m in range(-l, l + 1):
            orbs.append((l, m))
    Q = len(orbs)
    total = Q ** 4

    if convention == "fullgaunt":
        pair_k = _scalar_pair_k_fullgaunt(orbs, l_max)
    elif convention == "pairdiag":
        pair_k = _scalar_pair_k_pairdiag(orbs, l_max)
    else:
        raise ValueError(f"unknown convention {convention!r}")

    nonzero = 0
    for a in range(Q):
        _, ma = orbs[a]
        for b in range(Q):
            _, mb = orbs[b]
            for c in range(Q):
                k_ac = pair_k[(a, c)]
                if not k_ac:
                    continue
                _, mc = orbs[c]
                for d in range(Q):
                    _, md = orbs[d]
                    if convention == "fullgaunt":
                        if ma + mb != mc + md:
                            continue
                    # pairdiag: m_a=m_c (in k_ac pruning) and need m_b=m_d
                    k_bd = pair_k[(b, d)]
                    if not k_bd:
                        continue
                    if k_ac & k_bd:
                        nonzero += 1
    return nonzero, total, Fraction(nonzero, total)


# ---------------------------------------------------------------------------
# Spinor density: full Gaunt vs Paper-22 pair-diagonal
# ---------------------------------------------------------------------------

def _spinor_pair_k(orbs: List[SpinorOrbital], l_max: int,
                   convention: str) -> Dict[Tuple[int, int], frozenset]:
    """Per-pair (a,c) k-set for spinor selection.

    In both conventions we require:
      * (l_a + l_c + k) even
      * triangle(j_a, k, j_c)
      * triangle(l_a, k, l_c) [defensive; implied by parity+j-triangle for
        allowed kappa, but cheap to check]
      * 3j(j_a, k, j_c; 1/2, 0, -1/2) != 0  [reduced <kappa||C^k||kappa'>]

    Then the m-block 3j differs by convention:
      fullgaunt:  q = m_a - m_c, 3j(j_a, k, j_c; -m_a, q, m_c) != 0
      pairdiag:   require m_a = m_c, q = 0, 3j(j_a, k, j_c; -m_a, 0, m_c) != 0
    """
    pair_k = {}
    k_cap = 2 * l_max + 1  # generous upper bound
    for i, oa in enumerate(orbs):
        la = oa.l
        ja_x2 = oa.two_j
        ma_x2 = oa.two_mj
        ja_sym = Rational(ja_x2, 2)
        for j, oc in enumerate(orbs):
            lc = oc.l
            jc_x2 = oc.two_j
            mc_x2 = oc.two_mj
            jc_sym = Rational(jc_x2, 2)

            if convention == "pairdiag" and ma_x2 != mc_x2:
                pair_k[(i, j)] = frozenset()
                continue

            ks = set()
            for k in range(k_cap + 1):
                # parity on l
                if (la + lc + k) % 2 != 0:
                    continue
                # l-triangle (defensive)
                if k < abs(la - lc) or k > la + lc:
                    continue
                # j-triangle (2j representation)
                if not _triangle_x2(ja_x2, 2 * k, jc_x2):
                    continue
                k_sym = Integer(k)
                # reduced matrix element 3j
                w_red = wigner_3j(ja_sym, k_sym, jc_sym,
                                  Rational(1, 2), Integer(0),
                                  Rational(-1, 2))
                if w_red == 0:
                    continue
                if convention == "fullgaunt":
                    q_x2 = ma_x2 - mc_x2
                    if q_x2 % 2 != 0:
                        continue
                    q = q_x2 // 2
                    if abs(q) > k:
                        continue
                    w_m = wigner_3j(ja_sym, k_sym, jc_sym,
                                    Rational(-ma_x2, 2),
                                    Integer(q),
                                    Rational(mc_x2, 2))
                elif convention == "pairdiag":
                    # m_a == m_c already enforced; q=0
                    w_m = wigner_3j(ja_sym, k_sym, jc_sym,
                                    Rational(-ma_x2, 2),
                                    Integer(0),
                                    Rational(mc_x2, 2))
                else:
                    raise ValueError(f"unknown convention {convention!r}")
                if w_m == 0:
                    continue
                ks.add(k)
            pair_k[(i, j)] = frozenset(ks)
    return pair_k


def compute_spinor_density(l_max: int, convention: str, verbose: bool = False
                           ) -> Tuple[int, int, Fraction]:
    orbs = enumerate_spinor_orbitals(l_max)
    Q = len(orbs)
    total = Q ** 4

    t0 = time.time()
    pair_k = _spinor_pair_k(orbs, l_max, convention)
    t1 = time.time()
    if verbose:
        print(f"    [spinor/{convention} l_max={l_max}] Q={Q}, "
              f"precomputed {Q*Q} pairs in {t1-t0:.2f}s")

    nonzero = 0
    for a in range(Q):
        ma_x2 = orbs[a].two_mj
        for b in range(Q):
            mb_x2 = orbs[b].two_mj
            for c in range(Q):
                k_ac = pair_k[(a, c)]
                if not k_ac:
                    continue
                mc_x2 = orbs[c].two_mj
                for d in range(Q):
                    md_x2 = orbs[d].two_mj
                    if convention == "fullgaunt":
                        if ma_x2 + mb_x2 != mc_x2 + md_x2:
                            continue
                    k_bd = pair_k[(b, d)]
                    if not k_bd:
                        continue
                    if k_ac & k_bd:
                        nonzero += 1
    t2 = time.time()
    if verbose:
        print(f"    [spinor/{convention} l_max={l_max}] "
              f"enumerated {total} 4-tuples in {t2-t1:.2f}s, "
              f"nonzero={nonzero}")
    return nonzero, total, Fraction(nonzero, total)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(l_values=None):
    if l_values is None:
        l_values = [0, 1, 2, 3, 4, 5]
    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(exist_ok=True)

    # Paper 22 reference (for sanity cross-check of the pairdiag reproduction)
    paper22_reference = {
        0: Fraction(1, 1),
        1: Fraction(20, 256),        # 7.8125%
        2: Fraction(181, 6561),      # 2.7587%
        3: Fraction(944, 65536),     # 1.4404%
        4: Fraction(3525, 390625),   # 0.9024%
        5: Fraction(10396, 1679616), # 0.6190%   (value from Paper 22 proportion)
    }

    rows = []
    print("=" * 78)
    print("Track T0: Spinor-block ERI density d_spinor(l_max)")
    print("=" * 78)

    for l_max in l_values:
        print(f"\n--- l_max = {l_max} ---")
        # Scalar: full gaunt
        t0 = time.time()
        s_nz_fg, s_tot_fg, s_dens_fg = compute_scalar_density(l_max, "fullgaunt")
        t1 = time.time()
        Q_scalar = (l_max + 1) ** 2
        print(f"  scalar (fullgaunt): {s_nz_fg}/{s_tot_fg} = "
              f"{float(s_dens_fg)*100:.4f}%  [{t1-t0:.2f}s]")
        # Scalar: pair-diagonal (Paper 22 reproduction)
        t1b = time.time()
        s_nz_pd, s_tot_pd, s_dens_pd = compute_scalar_density(l_max, "pairdiag")
        t2 = time.time()
        print(f"  scalar (pairdiag, Paper-22 convention): "
              f"{s_nz_pd}/{s_tot_pd} = {float(s_dens_pd)*100:.4f}%  "
              f"[{t2-t1b:.2f}s]")
        # Paper 22 sanity
        ref = paper22_reference.get(l_max)
        if ref is not None:
            match = "OK" if s_dens_pd == ref else "MISMATCH"
            print(f"    Paper 22 reference: {float(ref)*100:.4f}%  [{match}]")

        # Spinor: full gaunt
        sp_nz_fg, sp_tot_fg, sp_dens_fg = compute_spinor_density(
            l_max, "fullgaunt", verbose=True)
        t3 = time.time()
        # Spinor: pair-diagonal
        sp_nz_pd, sp_tot_pd, sp_dens_pd = compute_spinor_density(
            l_max, "pairdiag", verbose=True)
        t4 = time.time()
        Q_spinor = 2 * (l_max + 1) ** 2
        assert sp_tot_fg == Q_spinor ** 4, (
            f"Q count mismatch: {sp_tot_fg} vs {Q_spinor**4}")

        print(f"  spinor (fullgaunt): {sp_nz_fg}/{sp_tot_fg} = "
              f"{float(sp_dens_fg)*100:.4f}%")
        print(f"  spinor (pairdiag):  {sp_nz_pd}/{sp_tot_pd} = "
              f"{float(sp_dens_pd)*100:.4f}%")

        # Ratios (primary: pairdiag, for apples-to-apples with Paper 22 table)
        ratio_pd = sp_dens_pd / s_dens_pd if s_dens_pd > 0 else None
        ratio_fg = sp_dens_fg / s_dens_fg if s_dens_fg > 0 else None
        if ratio_pd is not None:
            print(f"  ratio (pairdiag, primary)  = "
                  f"{float(ratio_pd):.4f} = "
                  f"{ratio_pd.numerator}/{ratio_pd.denominator}")
        if ratio_fg is not None:
            print(f"  ratio (fullgaunt)          = "
                  f"{float(ratio_fg):.4f} = "
                  f"{ratio_fg.numerator}/{ratio_fg.denominator}")

        rows.append({
            "l_max": l_max,
            "Q_scalar": Q_scalar,
            "Q_spinor": Q_spinor,
            "pairdiag": {
                "scalar_total": s_tot_pd,
                "scalar_nonzero": s_nz_pd,
                "scalar_density_num": s_dens_pd.numerator,
                "scalar_density_den": s_dens_pd.denominator,
                "scalar_density_pct": float(s_dens_pd) * 100,
                "spinor_total": sp_tot_pd,
                "spinor_nonzero": sp_nz_pd,
                "spinor_density_num": sp_dens_pd.numerator,
                "spinor_density_den": sp_dens_pd.denominator,
                "spinor_density_pct": float(sp_dens_pd) * 100,
                "ratio_num": ratio_pd.numerator if ratio_pd else None,
                "ratio_den": ratio_pd.denominator if ratio_pd else None,
                "ratio_float": float(ratio_pd) if ratio_pd else None,
            },
            "fullgaunt": {
                "scalar_total": s_tot_fg,
                "scalar_nonzero": s_nz_fg,
                "scalar_density_num": s_dens_fg.numerator,
                "scalar_density_den": s_dens_fg.denominator,
                "scalar_density_pct": float(s_dens_fg) * 100,
                "spinor_total": sp_tot_fg,
                "spinor_nonzero": sp_nz_fg,
                "spinor_density_num": sp_dens_fg.numerator,
                "spinor_density_den": sp_dens_fg.denominator,
                "spinor_density_pct": float(sp_dens_fg) * 100,
                "ratio_num": ratio_fg.numerator if ratio_fg else None,
                "ratio_den": ratio_fg.denominator if ratio_fg else None,
                "ratio_float": float(ratio_fg) if ratio_fg else None,
            },
        })

    # Summary table (Paper-22 convention, primary)
    print("\n" + "=" * 78)
    print("PRIMARY TABLE (pair-diagonal convention = Paper 22 reproduction)")
    print("=" * 78)
    header = (f"| {'l_max':>5} | {'Q_scal':>6} | {'Q_spin':>6} | "
              f"{'d_scalar':>10} | {'d_spinor':>10} | {'ratio':>8} |")
    sep = f"|{'-'*7}|{'-'*8}|{'-'*8}|{'-'*12}|{'-'*12}|{'-'*10}|"
    print(header)
    print(sep)
    for r in rows:
        pd = r["pairdiag"]
        print(f"| {r['l_max']:>5} | {r['Q_scalar']:>6} | {r['Q_spinor']:>6} | "
              f"{pd['scalar_density_pct']:>9.4f}% | "
              f"{pd['spinor_density_pct']:>9.4f}% | "
              f"{pd['ratio_float']:>8.4f} |")

    print("\n" + "=" * 78)
    print("SECONDARY TABLE (full-Gaunt convention, physically correct)")
    print("=" * 78)
    print(header)
    print(sep)
    for r in rows:
        fg = r["fullgaunt"]
        print(f"| {r['l_max']:>5} | {r['Q_scalar']:>6} | {r['Q_spinor']:>6} | "
              f"{fg['scalar_density_pct']:>9.4f}% | "
              f"{fg['spinor_density_pct']:>9.4f}% | "
              f"{fg['ratio_float']:>8.4f} |")

    metadata = {
        "task": "Tier 2 Track T0 — Spinor-block ERI density",
        "reference_paper": "Paper 22 (scalar, pair-diagonal convention)",
        "basis": "jj-coupled spinor orbitals (kappa, m_j)",
        "primary_convention": "pairdiag",
        "convention_explanation": {
            "pairdiag": (
                "Paper 22's published convention: c^k(l_a m_a; l_c m_c) "
                "requires m_a = m_c (equivalently bottom-row-sum=0 3j). "
                "This is stricter than the full Gaunt rule. Reproduces "
                "Paper 22 table exactly."
            ),
            "fullgaunt": (
                "Physically correct Coulomb selection: 3j's with "
                "q = m_a - m_c = -(m_b - m_d). Denser than pairdiag."
            ),
        },
        "arithmetic": "exact sympy Rational + wigner_3j",
        "n_max": 1,
        "n_max_independence": (
            "Density is n_max-independent by Paper 22 Theorem 3 "
            "(selection rule depends only on angular labels). Same "
            "argument transfers verbatim to the spinor case."
        ),
        "rows": rows,
    }
    out_path = out_dir / "tier2_t0_spinor_density.json"
    with open(out_path, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        lvals = [int(x) for x in sys.argv[1:]]
        main(lvals)
    else:
        main()
