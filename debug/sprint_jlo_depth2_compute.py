"""
Sprint JLO-Depth2 driver — Reading A vs Reading B disambiguation via the
JLO entire cyclic cocycle at degree n = 2.

Context
-------
This is the named follow-on from v3.83.0 Reading C-strong
(`debug/sprint_na1_offdiag_substrate_memo.md`, 2026-06-06).  The trace-
functional Mellin probe `Tr(f_1(D) gamma f_2(D))` collapsed bit-exactly
to depth-1 on BOTH the CH-diagonal AND the off-diagonal CH substrates,
ruling the trace-of-single-gamma-insertion probe as structurally
insensitive to the primitive (Reading A) vs deconcatenation
(Reading B) distinction.

Paper 55 Remark `rem:reading_AB_jlo_cocycle` names the JLO depth-2
cocycle

    chi_D(a_0, a_1, a_2)
      = int_{Delta_2} Tr(a_0 e^{-s_0 D^2} [D, a_1] e^{-s_1 D^2}
                                          [D, a_2] e^{-s_2 D^2}) ds

as the structurally appropriate probe because TWO non-D operators are
inserted between heat-kernel factors, not one.  The collapse argument
(which says Tr(A gamma B) for diagonal A, B picks up only diag(gamma))
does NOT apply at depth 2:

    Tr(a_0 e^{-s_0 D^2} [D, a_1] e^{-s_1 D^2} [D, a_2] e^{-s_2 D^2})

is a product of FOUR matrices in the eigenbasis of D (three heat
kernels + a_0 are diagonal; [D, a_1] and [D, a_2] are NOT diagonal),
so the trace can carry off-diagonal entries of the inserted operators.
The cocycle structure carries operator ordering algebraically rather
than projecting onto the trace diagonal.

The Reading A / B / INCONCLUSIVE decision gate
-----------------------------------------------
Three diagnostics:

(I) Swap symmetry test.  Compute chi_D(1, e_s, e_t) and chi_D(1, e_t,
    e_s) on bit-exact sympy at n_max = 2 across the 25 ordered pairs
    (e_s, e_t) of sector idempotents.
    - Reading A (abelian / primitive / cocommutative): the swap leaves
      the cocycle invariant up to a (purely combinatorial) sign.
      Mechanism: cocommutative Hopf algebra coproducts are symmetric,
      so primitive products factor symmetrically.
    - Reading B (shuffle / deconcatenation): the swap produces
      asymmetric content not derivable from a symmetric primitive
      product structure.
    - Reading C-strong already verified at the trace level — but the
      JLO cocycle is NOT a trace of a single-gamma insertion at
      depth-2; the collapse argument fails by counting the inserted
      operators.

(II) Primitive-product factorisation test.  Reading A predicts the
     depth-2 cocycle factorises through depth-1 cocycle data:
       chi_D(1, e_s, e_t) =? combination of {chi_D(1, e_s),
                                              chi_D(1, e_t),
                                              chi_D(1, 1)}
     under a fixed factorisation rule (primitive product /
     unshuffled).  Operationally: PSLQ at multiple precisions
     between chi_D^{depth-2}(1, e_s, e_t) and the depth-1 panel.
     - A hit: chi_D^{(2)} is in the Q-span of depth-1 data times
       rational coefficients (Reading A).
     - No hit + algebraic cocycle content not derivable from
       depth-1: Reading B (genuinely new depth-2 content).
     - Hits with combinations that exactly match the shuffle
       coproduct (chi_D(1, e_s) * chi_D(1, e_t) +
        chi_D(1, e_t) * chi_D(1, e_s)): also Reading A
       (primitive-product shape).

(III) Cyclic-invariance audit.  The JLO cocycle is entire cyclic
      (modulo coboundary).  Check chi_D(1, e_s, e_t) vs
      chi_D(e_t, 1, e_s) and chi_D(e_s, e_t, 1) — these should agree
      bit-exactly at the cocycle level (a cocommutativity vs cyclic
      diagnostic).

A POSITIVE Reading A verdict requires:
  - (I) swap-symmetric (sign-corrected match) for at least 80% of
    panel pairs
  - (II) primitive-product factorisation PSLQ hit at cross-precision
  - (III) cyclic identity bit-exact

A POSITIVE Reading B verdict requires:
  - (I) swap-asymmetric: at least one panel pair gives genuinely
    distinct content (not just a sign flip)
  - (II) no primitive-product PSLQ hit + algebraic structure visible
  - (III) cyclic identity still holds (Reading B is cyclic but not
    cocommutative)

An INCONCLUSIVE verdict reports what additional structure is
needed.

Implementation
--------------
We use the existing JLO computation infrastructure
`debug/compute_jlo_nmax2.py` via `jlo_cochain_coeffs` (rederived
inline here so this sprint is self-contained) and adapt it for the
A-vs-B-specific panel:
  - Leading t^0 coefficient: (1/2!) Tr_s(a_0 [D, a_1] [D, a_2])
  - Next-order t^1: (-1/3!) sum over partitions of m_total = 1
    of Tr_s(a_0 D^{2m_0} [D, a_1] D^{2m_1} [D, a_2] D^{2m_2})

Both EVEN (with gamma) and ODD (without gamma) flavours computed.
The EVEN JLO is the relevant cocycle for the chirality-graded
spectral triple; the ODD JLO is also computed as a control.

Substrate options:
  - CH-diagonal D = Lambda + kappa*A_graph with kappa = -1/16
    (current Q5'-Stage1 default at n_max = 2)
  - For the depth-2 probe, the off-diagonal entries of [D, e_s] are
    automatically present because [Lambda, e_s] = 0 only if e_s is
    sector-diagonal in the (n, l) basis (which it is for sector
    idempotents), but [kappa*A_graph, e_s] is generically non-zero
    when A_graph connects different sectors (which it does — the
    Fock-projected S^3 adjacency carries off-diagonal in the
    (n, l) sector structure).

So with the default substrate D = Lambda + kappa * A_graph the
commutator [D, e_s] = kappa * [A_graph, e_s] is purely off-diagonal
in the sector structure; the depth-2 cocycle naturally inherits
this off-diagonal content.

PSLQ defaults
-------------
Cross-precision (50, 100 dps) on the primitive-product test, ceiling
10^6, maxsteps 2000.

Output
------
- `debug/data/jlo_depth2_results.json`: full panel of cocycle values
  (sympy Rational), swap-symmetry table, cyclic-invariance audit,
  PSLQ panel for primitive-product test.
- Verdict logged to stdout and JSON.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import sympy as sp
from sympy import Integer, Matrix, Rational, eye as sp_eye, factorial as sp_fact, zeros as sp_zeros

from geovac.spectral_triple import FockSpectralTriple


# =====================================================================
# Section 1: Substrate construction (sympy Rational)
# =====================================================================


def sector_idempotent(st: FockSpectralTriple, sector_idx: int) -> Matrix:
    """Return the sector idempotent e_s as a dim_H x dim_H diagonal
    rational matrix:  e_s[i, i] = 1 iff state i belongs to sector s."""
    N = st.dim_H
    M = sp_zeros(N, N)
    for i in range(N):
        if st._state_to_sector[i] == sector_idx:
            M[i, i] = Integer(1)
    return M


def commutator(D: Matrix, a: Matrix) -> Matrix:
    """Compute [D, a] = D a - a D."""
    return D * a - a * D


def partitions(total: int, n_parts: int) -> List[Tuple[int, ...]]:
    """All (m_0, ..., m_{n_parts-1}) with sum = total, m_i >= 0."""
    if n_parts == 1:
        return [(total,)]
    out = []
    for m0 in range(total + 1):
        for sub in partitions(total - m0, n_parts - 1):
            out.append((m0,) + sub)
    return out


# =====================================================================
# Section 2: JLO depth-2 cocycle
# =====================================================================


def jlo_depth2_coeffs(
    D: Matrix,
    D2: Matrix,
    gamma: Matrix,
    a0: Matrix,
    a1: Matrix,
    a2: Matrix,
    M_max: int = 2,
    use_gamma: bool = True,
) -> List[sp.Expr]:
    """Compute the t-power-series coefficients of the JLO depth-2
    cochain

      chi_D(a_0, a_1, a_2; t) = sum_{m = 0}^{M_max} t^m * c_m

    where

      c_m = (-1)^m / (m + 2)! *
            sum_{(m_0, m_1, m_2) | sum = m}
              Tr_s(a_0 D^{2 m_0} [D, a_1] D^{2 m_1} [D, a_2] D^{2 m_2})

    and Tr_s = Tr(gamma . *) for the EVEN flavour (use_gamma=True) and
    Tr_s = Tr otherwise.

    Returns list of length M_max + 1, coefficients c_m, m = 0, ..., M_max.
    """
    comm_a1 = commutator(D, a1)
    comm_a2 = commutator(D, a2)

    # Precompute D^{2k}
    D2_powers: List[Matrix] = [sp_eye(D.shape[0])]
    for _ in range(M_max):
        D2_powers.append(D2_powers[-1] * D2)

    sign_gamma = gamma if use_gamma else sp_eye(D.shape[0])

    coeffs: List[sp.Expr] = []
    for m_total in range(M_max + 1):
        s = Integer(0)
        for (m0, m1, m2) in partitions(m_total, 3):
            prod = a0 * D2_powers[m0] * comm_a1 * D2_powers[m1] * comm_a2 * D2_powers[m2]
            tr = (sign_gamma * prod).trace()
            s = s + tr
        prefac = Rational((-1) ** m_total, sp_fact(m_total + 2))
        coeffs.append(sp.simplify(prefac * s))
    return coeffs


def jlo_depth1_coeffs(
    D: Matrix,
    D2: Matrix,
    gamma: Matrix,
    a0: Matrix,
    a1: Matrix,
    M_max: int = 2,
    use_gamma: bool = True,
) -> List[sp.Expr]:
    """Compute the t-power-series coefficients of the JLO depth-1
    cochain chi_D(a_0, a_1; t).  Used for the primitive-product test."""
    comm_a1 = commutator(D, a1)
    D2_powers: List[Matrix] = [sp_eye(D.shape[0])]
    for _ in range(M_max):
        D2_powers.append(D2_powers[-1] * D2)
    sign_gamma = gamma if use_gamma else sp_eye(D.shape[0])
    coeffs: List[sp.Expr] = []
    for m_total in range(M_max + 1):
        s = Integer(0)
        for (m0, m1) in partitions(m_total, 2):
            prod = a0 * D2_powers[m0] * comm_a1 * D2_powers[m1]
            tr = (sign_gamma * prod).trace()
            s = s + tr
        prefac = Rational((-1) ** m_total, sp_fact(m_total + 1))
        coeffs.append(sp.simplify(prefac * s))
    return coeffs


def jlo_depth0_coeffs(
    D: Matrix,
    D2: Matrix,
    gamma: Matrix,
    a0: Matrix,
    M_max: int = 2,
    use_gamma: bool = True,
) -> List[sp.Expr]:
    """Compute the t-power-series coefficients of the JLO depth-0
    cochain chi_D(a_0; t).  Used as control."""
    D2_powers: List[Matrix] = [sp_eye(D.shape[0])]
    for _ in range(M_max):
        D2_powers.append(D2_powers[-1] * D2)
    sign_gamma = gamma if use_gamma else sp_eye(D.shape[0])
    coeffs: List[sp.Expr] = []
    for m_total in range(M_max + 1):
        prod = a0 * D2_powers[m_total]
        tr = (sign_gamma * prod).trace()
        prefac = Rational((-1) ** m_total, sp_fact(m_total))
        coeffs.append(sp.simplify(prefac * tr))
    return coeffs


# =====================================================================
# Section 3: Panel construction
# =====================================================================


def build_full_depth2_panel(
    st: FockSpectralTriple,
    D: Matrix,
    D2: Matrix,
    gamma: Matrix,
    M_max: int = 1,
) -> Dict:
    """For each ordered triple (a_0, a_1, a_2) with
    a_0 in {1, e_0, ..., e_4} and (a_1, a_2) in {(e_s, e_t) : s, t in 0..4},
    compute the depth-2 JLO cocycle as exact sympy rationals.

    M_max = 1 captures the leading and next-to-leading t-coefficients;
    full sweep at all 5^3 = 125 triples (sympy-exact) takes ~3-5 minutes
    at n_max = 2.

    To keep wall time bounded, we run two reduced panels:
      Panel A: a_0 = 1, (a_1, a_2) sweeping all 25 (e_s, e_t) ordered
               pairs.  Primary Reading A/B test.
      Panel B: a_0 = e_s, (a_1, a_2) = (e_t, e_u) for a small targeted
               sub-panel (10 representative triples) — sanity / extension
               check.
    """
    N = st.dim_H
    unit = sp_eye(N)
    sectors = list(range(st.n_sectors))
    e_list = [sector_idempotent(st, s) for s in sectors]

    # Panel A: a_0 = 1, all 25 ordered (s, t) pairs
    panel_A: Dict[str, Dict] = {}
    print(f"    [Panel A] {len(sectors) ** 2} ordered (s, t) pairs at a_0 = 1...")
    for s in sectors:
        for t in sectors:
            key = f"(1, e_{s}, e_{t})"
            coeffs_even = jlo_depth2_coeffs(D, D2, gamma, unit, e_list[s], e_list[t], M_max=M_max, use_gamma=True)
            coeffs_odd = jlo_depth2_coeffs(D, D2, gamma, unit, e_list[s], e_list[t], M_max=M_max, use_gamma=False)
            panel_A[key] = {
                "s": s,
                "t": t,
                "even": [str(c) for c in coeffs_even],
                "odd": [str(c) for c in coeffs_odd],
                "even_rational": [c for c in coeffs_even],  # kept as sympy for downstream analysis
                "odd_rational": [c for c in coeffs_odd],
            }

    # Panel B: a_0 = e_s, (a_1, a_2) = (e_t, e_u) — smaller targeted panel
    panel_B: Dict[str, Dict] = {}
    panel_B_triples = [
        (0, 1, 2), (0, 2, 1),  # swap test with a_0 = e_0
        (1, 0, 2), (1, 2, 0),
        (2, 0, 1), (2, 1, 0),
        (0, 3, 4), (0, 4, 3),  # swap test with non-adjacent sectors
        (1, 1, 1),             # all-same diagnostic
        (2, 2, 2),
    ]
    print(f"    [Panel B] {len(panel_B_triples)} targeted (s, t, u) triples...")
    for (s, t, u) in panel_B_triples:
        key = f"(e_{s}, e_{t}, e_{u})"
        coeffs_even = jlo_depth2_coeffs(D, D2, gamma, e_list[s], e_list[t], e_list[u], M_max=M_max, use_gamma=True)
        coeffs_odd = jlo_depth2_coeffs(D, D2, gamma, e_list[s], e_list[t], e_list[u], M_max=M_max, use_gamma=False)
        panel_B[key] = {
            "s": s, "t": t, "u": u,
            "even": [str(c) for c in coeffs_even],
            "odd": [str(c) for c in coeffs_odd],
            "even_rational": [c for c in coeffs_even],
            "odd_rational": [c for c in coeffs_odd],
        }
    return {"panel_A": panel_A, "panel_B": panel_B}


def build_depth1_panel(
    st: FockSpectralTriple,
    D: Matrix,
    D2: Matrix,
    gamma: Matrix,
    M_max: int = 2,
) -> Dict:
    """Compute depth-1 JLO cocycle coefficients chi_D(a_0, a_1; t) for
    (a_0, a_1) = (1, e_s) and (e_s, 1) for each sector s.  Used as
    factorisation candidates in Reading A test."""
    N = st.dim_H
    unit = sp_eye(N)
    panel: Dict[str, Dict] = {}
    for s in range(st.n_sectors):
        e_s = sector_idempotent(st, s)
        coeffs_even_us = jlo_depth1_coeffs(D, D2, gamma, unit, e_s, M_max=M_max, use_gamma=True)
        coeffs_odd_us = jlo_depth1_coeffs(D, D2, gamma, unit, e_s, M_max=M_max, use_gamma=False)
        coeffs_even_su = jlo_depth1_coeffs(D, D2, gamma, e_s, unit, M_max=M_max, use_gamma=True)
        coeffs_odd_su = jlo_depth1_coeffs(D, D2, gamma, e_s, unit, M_max=M_max, use_gamma=False)
        panel[f"(1, e_{s})"] = {
            "even": [str(c) for c in coeffs_even_us],
            "odd": [str(c) for c in coeffs_odd_us],
            "even_rational": [c for c in coeffs_even_us],
            "odd_rational": [c for c in coeffs_odd_us],
        }
        panel[f"(e_{s}, 1)"] = {
            "even": [str(c) for c in coeffs_even_su],
            "odd": [str(c) for c in coeffs_odd_su],
            "even_rational": [c for c in coeffs_even_su],
            "odd_rational": [c for c in coeffs_odd_su],
        }
    return panel


# =====================================================================
# Section 4: Reading A vs Reading B diagnostics
# =====================================================================


def diagnostic_swap_symmetry(panel_A: Dict) -> Dict:
    """Compare chi_D(1, e_s, e_t) and chi_D(1, e_t, e_s) for each
    ordered pair.

    For each (s, t) with s != t:
      sym_pair := chi_D(1, e_s, e_t) + chi_D(1, e_t, e_s)   (cocommutative)
      asym_pair := chi_D(1, e_s, e_t) - chi_D(1, e_t, e_s)  (deconcat)

    Reading A: asym_pair = 0 (cocommutative product structure).
    Reading B: asym_pair != 0 at SOME pairs (genuine deconcatenation).
    """
    n_sectors = max(d["s"] for d in panel_A.values()) + 1
    M_max = len(next(iter(panel_A.values()))["even_rational"]) - 1
    swap_table: Dict[str, Dict] = {}
    # For each coefficient order m = 0, ..., M_max
    for m in range(M_max + 1):
        order_table: Dict[str, Dict] = {}
        n_zero_asym_even = 0
        n_zero_asym_odd = 0
        n_total = 0
        max_abs_asym_even = Integer(0)
        max_abs_asym_odd = Integer(0)
        for s in range(n_sectors):
            for t in range(s + 1, n_sectors):  # strict ordering, only test s < t
                k_st = f"(1, e_{s}, e_{t})"
                k_ts = f"(1, e_{t}, e_{s})"
                if k_st not in panel_A or k_ts not in panel_A:
                    continue
                v_st_even = panel_A[k_st]["even_rational"][m]
                v_ts_even = panel_A[k_ts]["even_rational"][m]
                v_st_odd = panel_A[k_st]["odd_rational"][m]
                v_ts_odd = panel_A[k_ts]["odd_rational"][m]
                asym_even = sp.simplify(v_st_even - v_ts_even)
                asym_odd = sp.simplify(v_st_odd - v_ts_odd)
                sym_even = sp.simplify(v_st_even + v_ts_even)
                sym_odd = sp.simplify(v_st_odd + v_ts_odd)
                order_table[f"({s}, {t})"] = {
                    "v(s,t)_even": str(v_st_even),
                    "v(t,s)_even": str(v_ts_even),
                    "asym_even": str(asym_even),
                    "sym_even": str(sym_even),
                    "v(s,t)_odd": str(v_st_odd),
                    "v(t,s)_odd": str(v_ts_odd),
                    "asym_odd": str(asym_odd),
                    "sym_odd": str(sym_odd),
                    "asym_even_is_zero": (asym_even == 0),
                    "asym_odd_is_zero": (asym_odd == 0),
                }
                n_total += 1
                if asym_even == 0:
                    n_zero_asym_even += 1
                if asym_odd == 0:
                    n_zero_asym_odd += 1
                if abs(asym_even) > max_abs_asym_even:
                    max_abs_asym_even = abs(asym_even)
                if abs(asym_odd) > max_abs_asym_odd:
                    max_abs_asym_odd = abs(asym_odd)
        swap_table[f"t^{m}"] = {
            "n_pairs_total": n_total,
            "n_asym_zero_even": n_zero_asym_even,
            "n_asym_zero_odd": n_zero_asym_odd,
            "frac_symmetric_even": float(n_zero_asym_even) / max(n_total, 1),
            "frac_symmetric_odd": float(n_zero_asym_odd) / max(n_total, 1),
            "max_abs_asym_even": str(max_abs_asym_even),
            "max_abs_asym_odd": str(max_abs_asym_odd),
            "pair_detail": order_table,
        }
    return swap_table


def diagnostic_cocommutativity_vs_shuffle(panel_A: Dict, depth1_panel: Dict, st: FockSpectralTriple) -> Dict:
    """Test whether chi_D(1, e_s, e_t) at t^0 admits a primitive-product
    factorisation through depth-1 cocycle data.

    The two natural candidate "Reading A" reductions are:

    (a) Pure primitive product:
          chi_D(1, e_s, e_t)|_{t^0} ?= chi_D(1, e_s)|_{t^?} * chi_D(1, e_t)|_{t^?}
        — but this is a product of scalars (cochains are scalars at
        fixed t), so we test PSLQ identifications between the depth-2
        value and rational multiples / quadratic combinations of depth-1
        values.

    (b) Symmetric primitive product (shuffle-2-of-depth-1):
          chi_D(1, e_s, e_t) + chi_D(1, e_t, e_s)
        ?= rational * [chi_D(1, e_s)|_{t^0} + chi_D(1, e_t)|_{t^0}]^2
        or similar shape.

    Here we run exact sympy Q-linear-relation tests between the depth-2
    panel at t^0 and t^1 against the depth-1 panel at t^0 and t^1.
    A genuine Reading A would show a CLEAN Q-linear identity (no large
    denominators).  Reading B predicts no such reduction.

    Output: for each (s, t) sector pair, the chi_D(1, e_s, e_t) at t^0
    and the depth-1 panel values; ratios computed where finite.
    """
    n_sectors = st.n_sectors

    # Extract depth-1 t^0 and t^1 values (even flavour)
    d1_even_t0: Dict[int, sp.Expr] = {}
    d1_even_t1: Dict[int, sp.Expr] = {}
    for s in range(n_sectors):
        d1_even_t0[s] = depth1_panel[f"(1, e_{s})"]["even_rational"][0]
        d1_even_t1[s] = depth1_panel[f"(1, e_{s})"]["even_rational"][1]

    # Depth-1 t^0 should be (1/1!) Tr(gamma [D, e_s]) = 0 (commutator
    # of bounded ops with diagonal gamma has zero diagonal in the
    # eigenbasis where gamma is constant on the relevant sub-block).
    # But the M_3 source Tr(gamma D e^{-tD^2}) at t^0 is Tr(gamma D)
    # which is NOT what chi_D(1, e_s) gives.

    factor_table: Dict = {}
    for s in range(n_sectors):
        for t in range(n_sectors):
            key = f"(1, e_{s}, e_{t})"
            d2_t0_even = panel_A[key]["even_rational"][0]
            d2_t0_odd = panel_A[key]["odd_rational"][0]
            entry: Dict = {
                "d2_t0_even": str(d2_t0_even),
                "d2_t0_odd": str(d2_t0_odd),
                "d1_t0_even_s": str(d1_even_t0[s]),
                "d1_t0_even_t": str(d1_even_t0[t]),
                "d1_t1_even_s": str(d1_even_t1[s]),
                "d1_t1_even_t": str(d1_even_t1[t]),
            }
            # Reading A candidates at t^0:
            # (i) Q-multiple of d1_t1_s * d1_t1_t (primitive product candidate)
            #     since d1_t0 = (1/1!)*Tr(gamma [D, e_s]) is typically 0
            # (ii) Q-multiple of (d1_t1_s + d1_t1_t)
            cand_i = sp.simplify(d1_even_t1[s] * d1_even_t1[t])
            cand_ii = sp.simplify(d1_even_t1[s] + d1_even_t1[t])
            if cand_i != 0:
                ratio_i = sp.simplify(d2_t0_even / cand_i)
            else:
                ratio_i = sp.nan
            if cand_ii != 0:
                ratio_ii = sp.simplify(d2_t0_even / cand_ii)
            else:
                ratio_ii = sp.nan
            entry["ratio_d2t0_over_d1t1_s_times_d1t1_t"] = str(ratio_i)
            entry["ratio_d2t0_over_d1t1_s_plus_d1t1_t"] = str(ratio_ii)
            factor_table[key] = entry

    # Aggregated test: is the t^0 depth-2 value a Q-rational reduction of
    # the depth-1 data?  We check whether the family {chi_D(1, e_s, e_t)|_{t^0}}
    # over all (s, t) lies in the Q-span of {chi_D(1, e_s)|_{t^1} * chi_D(1, e_t)|_{t^1}}.
    # Construct the linear-algebra question over Q.
    rhs_vector = []
    lhs_vectors = []  # each row: vector of candidates
    key_list = []
    for s in range(n_sectors):
        for t in range(n_sectors):
            key = f"(1, e_{s}, e_{t})"
            rhs_vector.append(panel_A[key]["even_rational"][0])
            lhs_vectors.append([
                d1_even_t1[s] * d1_even_t1[t],
                d1_even_t1[s] + d1_even_t1[t],
                d1_even_t0[s] * d1_even_t0[t],
                d1_even_t0[s] + d1_even_t0[t],
                Integer(1),  # constant
            ])
            key_list.append(key)

    M = sp.Matrix(lhs_vectors)
    rhs = sp.Matrix(rhs_vector)
    # Solve M * x = rhs in Q (least-squares unique if rank consistent)
    # Use the augmented-matrix rank test.
    M_aug = M.row_join(rhs)
    rank_M = M.rank()
    rank_M_aug = M_aug.rank()
    solvable = (rank_M == rank_M_aug)
    if solvable and rank_M == 5:
        # over-determined and consistent: unique solution
        try:
            sol = M.solve(rhs)
            sol_str = [str(sp.simplify(x)) for x in sol]
        except Exception as e:
            sol_str = [f"solve failed: {e}"]
    elif solvable:
        sol_str = [f"system consistent with rank {rank_M} < 5 (under-determined)"]
    else:
        sol_str = ["INCONSISTENT — no Q-linear reduction of depth-2 t^0 to depth-1 panel"]

    return {
        "factor_table": factor_table,
        "primitive_reduction_test": {
            "rank_M": int(rank_M),
            "rank_M_aug": int(rank_M_aug),
            "system_solvable": bool(solvable),
            "solution_or_diag": sol_str,
            "n_rows": int(len(lhs_vectors)),
            "candidate_basis": [
                "d1_t1_s * d1_t1_t (primitive product candidate)",
                "d1_t1_s + d1_t1_t (sum candidate)",
                "d1_t0_s * d1_t0_t",
                "d1_t0_s + d1_t0_t",
                "constant 1",
            ],
            "key_list": key_list,
        },
    }


def diagnostic_diagonal_subblock(panel_A: Dict, st: FockSpectralTriple) -> Dict:
    """Reading C-strong already established trace-of-single-gamma-insertion
    collapses.  The JLO depth-2 has TWO commutators, so the collapse
    argument doesn't apply termwise.  However, on the COMMUTATIVE algebra
    A = C^5 of sector idempotents on the diagonal CH substrate
    (kappa = -1/16, Lambda + kappa * A_graph), we expect:
      [Lambda, e_s] = 0 (Lambda diagonal, e_s sector-diagonal in eigenbasis)
      [kappa*A, e_s] is purely OFF-DIAGONAL between sectors
    so the depth-2 cocycle at t^0 reduces to
      (1/2) Tr(gamma * a_0 * kappa^2 * [A, e_s] [A, e_t])
    which is a fully off-diagonal-times-off-diagonal product — the
    trace is non-zero only when [A, e_s] [A, e_t] picks up diagonal
    structure (which happens only when s = t or via specific
    sector-coupling patterns).

    This is the structural content we extract here.

    For (s, t) with s != t and both e_s, e_t connected by A_graph,
    chi_D(1, e_s, e_t)|_{t^0} = (1/2) * kappa^2 * Tr(gamma [A, e_s] [A, e_t])
    """
    diagonal_results = {}
    for key, entry in panel_A.items():
        s_idx = entry["s"]
        t_idx = entry["t"]
        t0_even = entry["even_rational"][0]
        t0_odd = entry["odd_rational"][0]
        diagonal_results[key] = {
            "s": s_idx,
            "t": t_idx,
            "same_sector": (s_idx == t_idx),
            "t0_even": str(t0_even),
            "t0_odd": str(t0_odd),
            "t0_even_iszero": (t0_even == 0),
            "t0_odd_iszero": (t0_odd == 0),
        }
    n_total = len(diagonal_results)
    n_even_zero = sum(1 for v in diagonal_results.values() if v["t0_even_iszero"])
    n_odd_zero = sum(1 for v in diagonal_results.values() if v["t0_odd_iszero"])
    return {
        "n_panel_pairs": n_total,
        "n_t0_even_zero": n_even_zero,
        "n_t0_odd_zero": n_odd_zero,
        "frac_t0_even_zero": n_even_zero / n_total,
        "frac_t0_odd_zero": n_odd_zero / n_total,
        "panel_details": diagonal_results,
    }


def diagnostic_cyclic_invariance(
    st: FockSpectralTriple,
    D: Matrix,
    D2: Matrix,
    gamma: Matrix,
    M_max: int = 1,
) -> Dict:
    """Test cyclic invariance of chi_D modulo a coboundary:
      chi_D(a_0, a_1, a_2) ?= +/- chi_D(a_2, a_0, a_1) ?= +/- chi_D(a_1, a_2, a_0)

    For commutative algebra A, JLO is cyclic at the cochain level
    (not just modulo coboundary).  This is a load-bearing sanity check
    — if it fails bit-exact, the computation has a bug.

    Take a small panel of triples.
    """
    N = st.dim_H
    unit = sp_eye(N)
    e_list = [sector_idempotent(st, s) for s in range(st.n_sectors)]
    triples_to_test = [
        (unit, e_list[0], e_list[1]),
        (unit, e_list[0], e_list[2]),
        (unit, e_list[1], e_list[2]),
        (e_list[0], e_list[1], e_list[2]),
        (e_list[1], e_list[0], e_list[2]),
    ]
    cyc_table = {}
    for idx, (a, b, c) in enumerate(triples_to_test):
        # Original
        v_abc_even = jlo_depth2_coeffs(D, D2, gamma, a, b, c, M_max=M_max, use_gamma=True)
        # Cyclic shift 1: (c, a, b)
        v_cab_even = jlo_depth2_coeffs(D, D2, gamma, c, a, b, M_max=M_max, use_gamma=True)
        # Cyclic shift 2: (b, c, a)
        v_bca_even = jlo_depth2_coeffs(D, D2, gamma, b, c, a, M_max=M_max, use_gamma=True)
        # On commutative A, even JLO at degree 2 has parity (-1)^2 = +1 under cyclic shift
        # (cyclic with degree-n sign in entire-cyclic-cohomology conventions:
        #  (b'phi)(a_0, ..., a_n) := phi(a_0 a_1, a_2, ..., a_n) - phi(a_0, a_1 a_2, ...) +
        #                             ...  + (-1)^n phi(a_n a_0, a_1, ..., a_{n-1}))
        # In our normalisation chi_D itself satisfies
        #   chi_D(a_n, a_0, ..., a_{n-1}) = (-1)^n chi_D(a_0, ..., a_n)
        # so at n = 2 we expect chi_D(a_2, a_0, a_1) = + chi_D(a_0, a_1, a_2)
        # (sign +1 from (-1)^2)
        diffs_abc_cab = [sp.simplify(v_abc_even[m] - v_cab_even[m]) for m in range(M_max + 1)]
        diffs_abc_bca = [sp.simplify(v_abc_even[m] - v_bca_even[m]) for m in range(M_max + 1)]
        cyc_table[f"triple_{idx}"] = {
            "v_abc": [str(c_) for c_ in v_abc_even],
            "v_cab": [str(c_) for c_ in v_cab_even],
            "v_bca": [str(c_) for c_ in v_bca_even],
            "diff_abc_cab": [str(d) for d in diffs_abc_cab],
            "diff_abc_bca": [str(d) for d in diffs_abc_bca],
            "cyclic_holds_abc_cab": all(d == 0 for d in diffs_abc_cab),
            "cyclic_holds_abc_bca": all(d == 0 for d in diffs_abc_bca),
        }
    return cyc_table


# =====================================================================
# Section 5: PSLQ cross-check (depth-2 vs depth-1 panel + π and ζ)
# =====================================================================

def pslq_test_depth2_against_depth1(
    panel_A: Dict,
    depth1_panel: Dict,
    st: FockSpectralTriple,
    dps_list: Tuple[int, ...] = (50, 100),
) -> Dict:
    """For each depth-2 panel entry, run PSLQ at each dps in dps_list
    against a basis of depth-1 panel values, looking for low-coefficient
    Q-linear relations.

    Reading A is supported if a CLEAN relation (coeff ceiling <= 100,
    cross-precision agreement) emerges.  Reading B is supported if no
    such relation emerges across all panel entries.

    Implementation note: on the FINITE-DIMENSIONAL substrate at n_max = 2,
    every depth-2 value is already in Q (rational), so the "PSLQ at high
    dps" question is really "is the rational value a low-coefficient
    Q-linear combination of depth-1 rationals?".  We answer this by direct
    sympy linear-algebra (which is exact), not by mpmath PSLQ — the
    rational-arithmetic check is stronger.
    """
    import mpmath as mp
    n_sectors = st.n_sectors
    # Build basis of depth-1 panel rational values at t^0, t^1, t^2 (even flavour)
    basis_q: List[sp.Expr] = []
    basis_labels: List[str] = []
    for s in range(n_sectors):
        for m in range(3):
            basis_q.append(depth1_panel[f"(1, e_{s})"]["even_rational"][m])
            basis_labels.append(f"d1_t{m}_(1,e_{s})")
            basis_q.append(depth1_panel[f"(e_{s}, 1)"]["even_rational"][m])
            basis_labels.append(f"d1_t{m}_(e_{s},1)")
    # Add constant 1 to the basis
    basis_q.append(Integer(1))
    basis_labels.append("constant_1")

    pslq_table = {}
    for key, entry in panel_A.items():
        target_q_t0 = entry["even_rational"][0]
        if target_q_t0 == 0:
            pslq_table[key] = {
                "target_t0_even_q": "0",
                "trivial_zero_target": True,
            }
            continue
        # Direct Q-linear relation: does target_q_t0 sit in Q-span of basis_q?
        # Convert to mp.mpf at each dps for PSLQ; also check exact sympy
        # rational-linear combinations.
        target_q_t0 = sp.simplify(target_q_t0)
        # Exact Q-divisibility check: is target = sum c_i * basis_i for small c_i?
        # We seek minimal-norm rational solution via Smith normal form / nullspace.
        # Build vector v_target = target_q_t0, matrix M of basis vectors (here all scalars).
        # In Q-span on scalars, the question is simply whether
        # target / gcd(basis) = integer combination of basis / gcd values.
        # We use a more direct approach: sympy linsolve.
        coeffs_var = sp.symbols(f"c_0:{len(basis_q)}", rational=True)
        eq = sum(coeffs_var[i] * basis_q[i] for i in range(len(basis_q))) - target_q_t0
        # Find any rational solution (under-determined system)
        try:
            sol = sp.solve(eq, coeffs_var[0])  # solve for c_0 in terms of others
            sol_str = str(sol)
        except Exception as e:
            sol_str = f"sol failed: {e}"
        # Now run mpmath PSLQ at each dps for the integer-relation hunt
        pslq_results = {}
        for dps in dps_list:
            with mp.workdps(dps):
                # Convert sympy rationals to mp.mpf at this precision
                vec = [mp.mpf(sp.Rational(b).p) / mp.mpf(sp.Rational(b).q) if b != 0 else mp.mpf(0) for b in basis_q]
                target_mp = mp.mpf(sp.Rational(target_q_t0).p) / mp.mpf(sp.Rational(target_q_t0).q)
                vec_for_pslq = [target_mp] + vec
                try:
                    rel = mp.pslq(vec_for_pslq, maxcoeff=10**6, maxsteps=2000)
                    pslq_results[f"dps_{dps}"] = {
                        "relation_found": rel is not None,
                        "relation": [int(r) for r in rel] if rel is not None else None,
                    }
                except Exception as e:
                    pslq_results[f"dps_{dps}"] = {"error": str(e)}
        pslq_table[key] = {
            "target_t0_even_q": str(target_q_t0),
            "pslq_per_dps": pslq_results,
        }
    return {
        "basis_labels": basis_labels,
        "basis_values_q": [str(b) for b in basis_q],
        "pslq_table": pslq_table,
    }


# =====================================================================
# Section 6: Main
# =====================================================================


def main() -> None:
    print("=" * 70)
    print("Sprint JLO-Depth2 — Reading A vs Reading B disambiguation")
    print("via the JLO entire cyclic cocycle at depth n = 2")
    print("=" * 70)

    t_global = time.time()

    print("\n[1] Building Camporesi-Higuchi spectral triple at n_max = 2...")
    st = FockSpectralTriple(n_max=2)
    D = st.dirac_operator
    Lam = st.diagonal_part
    gamma = st.grading
    N = st.dim_H
    print(f"    dim_H = {N}, n_sectors = {st.n_sectors}, sectors = {st.sectors}")
    print(f"    kappa = {st._kappa}")
    print(f"    |D - Lambda| = {sp.simplify((D - Lam).norm())}")

    D2 = D * D
    print(f"    Tr(D^2) = {(D2).trace()}, Tr(Lambda^2) = {(Lam * Lam).trace()}")

    print("\n[2] Building depth-2 JLO panel (full 25 ordered pairs at a_0 = 1 + 10 targeted triples)...")
    t0 = time.time()
    full_panel = build_full_depth2_panel(st, D, D2, gamma, M_max=1)
    print(f"    done in {time.time() - t0:.1f}s")

    print("\n[3] Building depth-1 JLO panel (for primitive-product test)...")
    t0 = time.time()
    depth1_panel = build_depth1_panel(st, D, D2, gamma, M_max=2)
    print(f"    done in {time.time() - t0:.1f}s")

    print("\n[4] Swap-symmetry diagnostic (Reading A vs B test (I))...")
    swap_table = diagnostic_swap_symmetry(full_panel["panel_A"])
    for order_label, order_data in swap_table.items():
        print(f"    {order_label}: even sym frac = {order_data['frac_symmetric_even']:.3f}, "
              f"odd sym frac = {order_data['frac_symmetric_odd']:.3f}, "
              f"max|asym_even| = {order_data['max_abs_asym_even']}, "
              f"max|asym_odd| = {order_data['max_abs_asym_odd']}")

    print("\n[5] Cyclic invariance audit...")
    cyc_table = diagnostic_cyclic_invariance(st, D, D2, gamma, M_max=1)
    n_cyc_hold = sum(1 for v in cyc_table.values() if v["cyclic_holds_abc_cab"] and v["cyclic_holds_abc_bca"])
    print(f"    {n_cyc_hold} / {len(cyc_table)} cyclic-identity triples hold bit-exact")

    print("\n[6] Diagonal-subblock structural extraction...")
    diag_table = diagnostic_diagonal_subblock(full_panel["panel_A"], st)
    print(f"    n_panel_pairs = {diag_table['n_panel_pairs']}, "
          f"t^0 even zero count = {diag_table['n_t0_even_zero']} "
          f"({diag_table['frac_t0_even_zero']*100:.1f}%), "
          f"t^0 odd zero count = {diag_table['n_t0_odd_zero']} "
          f"({diag_table['frac_t0_odd_zero']*100:.1f}%)")

    print("\n[7] Primitive-product factorisation test (Reading A vs B test (II))...")
    primprod = diagnostic_cocommutativity_vs_shuffle(full_panel["panel_A"], depth1_panel, st)
    pr = primprod["primitive_reduction_test"]
    print(f"    Q-linear-reduction rank: M = {pr['rank_M']}, M_aug = {pr['rank_M_aug']}")
    print(f"    System solvable in Q-span: {pr['system_solvable']}")
    print(f"    Diagnostic: {pr['solution_or_diag']}")

    print("\n[8] PSLQ cross-check at depth-2 t^0 vs depth-1 panel + constant...")
    pslq_data = pslq_test_depth2_against_depth1(full_panel["panel_A"], depth1_panel, st, dps_list=(50, 100))
    n_pslq_hits = 0
    n_pslq_targets = 0
    for key, entry in pslq_data["pslq_table"].items():
        if entry.get("trivial_zero_target"):
            continue
        n_pslq_targets += 1
        hits = entry.get("pslq_per_dps", {})
        if any(h.get("relation_found") for h in hits.values()):
            n_pslq_hits += 1
    print(f"    PSLQ hits: {n_pslq_hits} / {n_pslq_targets} non-zero targets")

    # Decision gate
    print("\n[9] Reading A vs Reading B verdict...")
    # Use t^1 (the leading non-trivial order) for the verdict — t^0 is structurally zero
    t1_data = swap_table["t^1"]
    swap_frac_even_t1 = t1_data["frac_symmetric_even"]
    max_asym_even_t1 = t1_data["max_abs_asym_even"]
    # If max asymmetric at t^1 is non-zero AND fraction of asymmetric pairs > 0, Reading B leaning
    # If frac == 1.0 (all swap-symmetric at t^1), Reading A leaning
    reading_A_signal = (swap_frac_even_t1 == 1.0)
    reading_B_signal = (swap_frac_even_t1 < 1.0) and max_asym_even_t1 != "0"
    if reading_A_signal:
        verdict = "READING A leaning (all swap-symmetric at t^1)"
    elif reading_B_signal:
        verdict = "READING B leaning (genuine swap-asymmetric content at t^1)"
    else:
        verdict = "INCONCLUSIVE"
    print(f"    Verdict: {verdict}")
    print(f"    swap_frac_even at t^1 = {swap_frac_even_t1:.3f}")
    print(f"    max|asym_even| at t^1 = {max_asym_even_t1}")

    # Strip sympy-native fields for JSON
    def strip_rational(panel):
        out = {}
        for k, v in panel.items():
            new = {kk: vv for kk, vv in v.items() if kk not in ("even_rational", "odd_rational")}
            out[k] = new
        return out
    panel_A_clean = strip_rational(full_panel["panel_A"])
    panel_B_clean = strip_rational(full_panel["panel_B"])
    depth1_clean = strip_rational(depth1_panel)

    out = {
        "sprint": "Q5'-JLO-Depth2",
        "date": "2026-06-06",
        "purpose": "Reading A vs Reading B disambiguation via JLO entire cyclic cocycle at depth 2",
        "context": {
            "v3.83.0_reading_C_strong": True,
            "trace_functional_collapse_already_established": True,
            "JLO_cocycle_carries_operator_ordering_algebraically": True,
        },
        "substrate": {
            "n_max": 2,
            "dim_H": N,
            "n_sectors": st.n_sectors,
            "sectors": st.sectors,
            "kappa": str(st._kappa),
            "Tr_D_squared": str(D2.trace()),
            "Tr_Lambda_squared": str((Lam * Lam).trace()),
            "norm_D_minus_Lambda": str(sp.simplify((D - Lam).norm())),
        },
        "panel_A_depth2_a0_unit": panel_A_clean,
        "panel_B_depth2_targeted": panel_B_clean,
        "depth1_panel": depth1_clean,
        "swap_symmetry_diagnostic": swap_table,
        "cyclic_invariance_audit": cyc_table,
        "diagonal_subblock_extraction": diag_table,
        "primitive_product_test": primprod,
        "pslq_test": pslq_data,
        "verdict": verdict,
        "verdict_summary": {
            "swap_frac_even_t1": swap_frac_even_t1,
            "max_abs_asym_even_t1": max_asym_even_t1,
            "n_cyclic_identities_hold": n_cyc_hold,
            "n_cyclic_triples_tested": len(cyc_table),
            "primitive_product_system_solvable": pr["system_solvable"],
            "pslq_hits": n_pslq_hits,
            "pslq_targets": n_pslq_targets,
        },
        "wall_seconds": time.time() - t_global,
    }
    out_path = Path("debug/data/jlo_depth2_results.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\n[OUT] Data written to {out_path}")
    print(f"[TIME] Total wall: {time.time() - t_global:.1f}s")


if __name__ == "__main__":
    main()
