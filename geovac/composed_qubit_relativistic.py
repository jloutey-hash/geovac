"""Spin-ful composed qubit Hamiltonian (Track T3, Dirac-on-S^3 Tier 2).

Wires T1 closed-form spinor matrix elements (``geovac.dirac_matrix_elements``)
and T2 Breit-Pauli spin-orbit (``geovac.spin_orbit``) into the composed-geometry
pipeline.  Produces a Jordan-Wigner-mapped qubit Hamiltonian in the
(κ, m_j)-labeled basis — one spin-orbital per DiracLabel, no separate spin
doubling in the FermionOperator (the spin is already in κ via j = |κ|−1/2).

Algebraic-first design
----------------------
* Radial integrals R^k(n₁l₁, n₂l₂, n₃l₃, n₄l₄): reuse the scalar
  ``_compute_rk_integrals_block`` path from composed_qubit (exact algebraic
  ``hypergeometric_slater.get_rk_float`` for n_max ≤ 3 + Fraction fallback).
* Angular jj-coupled Gaunt coefficient X_k(κ_a m_a; κ_c m_c) is full-Gaunt
  (Dyall §9, Grant §7.5):

      X_k(a,c) = (−1)^{j_a−m_j,a} · √((2j_a+1)(2j_c+1))
                 · 3j(j_a k j_c; ½ 0 −½)
                 · 3j(j_a k j_c; −m_j,a q m_j,c)
                 · π(l_a+k+l_c even)

  where q = m_j,a − m_j,c.  Implemented symbolically in
  ``sympy.physics.wigner.wigner_3j`` then cast to float once.
* h₁: diagonal hydrogenic ``-Z²/(2n²)`` (α → 0 limit) + PK Gaussian barrier
  reused from ``composed_qubit._compute_pk_matrix_elements`` (diagonal in
  (l,m), and κ determines l so PK splits cleanly between κ-branches) +
  T2 spin-orbit ``so_diagonal_matrix_element(n, κ, Z, α)``.

Block diagonality: the loop over blocks never couples different blocks in
the ERI — same composed architecture invariant as the scalar builder.

Qubit count: Q = 2·(Σ_b M_b) where M_b is the scalar spatial-orbital count
of block b (each spatial orbital expands into its two j-branches, each with
(2j+1) m_j values — the total is exactly 2·M_b for non-spin-restricted
scalar comparison, because 2(l+1)² = 2·#{(n,l,m) with l ≤ n−1}).

Author: GeoVac Development Team
Track:  T3, Dirac-on-S^3 Tier 2 sprint
Date:   2026-04-15
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import sympy as sp
from sympy import Integer, Rational, sqrt
from sympy.physics.wigner import wigner_3j

from openfermion import FermionOperator, jordan_wigner

from geovac.composed_qubit import (
    _compute_rk_integrals_block,
    _compute_pk_matrix_elements,
    _enumerate_states,
)
from geovac.dirac_matrix_elements import DiracLabel, kappa_to_l, kappa_to_j
from geovac.spin_orbit import so_diagonal_matrix_element


# ---------------------------------------------------------------------------
# Spinor-label enumeration (DiracLabel list per block)
# ---------------------------------------------------------------------------

def enumerate_dirac_labels(max_n: int, l_min: int = 0) -> List[DiracLabel]:
    """All DiracLabel's with 1 ≤ n_fock ≤ max_n and l ≥ l_min.

    For each n, l (l_min ≤ l < n) and each κ branch (κ = −(l+1) for
    j = l+1/2; κ = +l for j = l−1/2, l ≥ 1), emit every m_j ∈ {−j..j}.
    """
    labels: List[DiracLabel] = []
    for n in range(1, max_n + 1):
        for l in range(max(l_min, 0), n):
            kappa_options = [-(l + 1)]
            if l >= 1:
                kappa_options.append(l)
            for kappa in kappa_options:
                two_j = 2 * abs(kappa) - 1
                for two_m_j in range(-two_j, two_j + 1, 2):
                    labels.append(DiracLabel(n_fock=n, kappa=kappa,
                                             two_m_j=two_m_j))
    return labels


# ---------------------------------------------------------------------------
# jj-coupled Gaunt angular coefficient X_k(κ_a m_a; κ_c m_c)
# ---------------------------------------------------------------------------

_X_CACHE: Dict[Tuple[int, int, int, int, int], float] = {}


def jj_angular_Xk(kappa_a: int, two_m_a: int,
                  kappa_c: int, two_m_c: int, k: int) -> float:
    """Full-Gaunt jj-coupled angular coefficient X_k(a,c).

    X_k = Π(l_a+k+l_c even)
          · (−1)^{(2j_a − two_m_a)/2}
          · √((2j_a+1)(2j_c+1))
          · 3j(j_a k j_c; ½ 0 −½)
          · 3j(j_a k j_c; −m_j,a  q  m_j,c)
    with q = m_j,a − m_j,c = (two_m_a − two_m_c)/2.

    Returns 0 if any selection rule (parity, triangle, q-integer) fails.
    Caches float value (the exact sympy expression is a combination of
    √rationals, coerced via ``float`` at the last step).
    """
    key = (kappa_a, two_m_a, kappa_c, two_m_c, k)
    if key in _X_CACHE:
        return _X_CACHE[key]

    la = kappa_to_l(kappa_a)
    lc = kappa_to_l(kappa_c)
    if (la + lc + k) % 2 != 0:
        _X_CACHE[key] = 0.0
        return 0.0
    two_j_a = 2 * abs(kappa_a) - 1
    two_j_c = 2 * abs(kappa_c) - 1
    # j-triangle (in 2j representation)
    if 2 * k < abs(two_j_a - two_j_c) or 2 * k > (two_j_a + two_j_c):
        _X_CACHE[key] = 0.0
        return 0.0
    q_two = two_m_a - two_m_c
    if q_two % 2 != 0:
        _X_CACHE[key] = 0.0
        return 0.0
    q = q_two // 2
    if abs(q) > k:
        _X_CACHE[key] = 0.0
        return 0.0

    j_a = Rational(two_j_a, 2)
    j_c = Rational(two_j_c, 2)
    k_sp = Integer(k)
    # Reduced 3j (m-independent part)
    w_red = wigner_3j(j_a, k_sp, j_c,
                      Rational(1, 2), Integer(0), Rational(-1, 2))
    if w_red == 0:
        _X_CACHE[key] = 0.0
        return 0.0
    # m-dependent 3j
    w_m = wigner_3j(j_a, k_sp, j_c,
                    Rational(-two_m_a, 2), Integer(q), Rational(two_m_c, 2))
    if w_m == 0:
        _X_CACHE[key] = 0.0
        return 0.0

    # Phase (−1)^{j_a − m_j,a}, always an integer because 2(j_a − m_j,a) is
    # an even integer (j_a half-int, m_j,a half-int, diff integer).
    phase_exp = (two_j_a - two_m_a) // 2
    phase = (-1) ** phase_exp

    prefac = sqrt(Integer(two_j_a + 1) * Integer(two_j_c + 1))
    expr = Integer(phase) * prefac * w_red * w_m
    val = float(expr)
    _X_CACHE[key] = val
    return val


# ---------------------------------------------------------------------------
# Block ERI build in spinor basis
# ---------------------------------------------------------------------------

def _build_spinor_eri_block(
    Z: float,
    spinor_labels: List[DiracLabel],
    rk_cache: Dict[Tuple[int, ...], float],
) -> Dict[Tuple[int, int, int, int], float]:
    """Build physicist-notation ERI ``<ab|1/r12|cd>`` on spinor basis.

    Selection rules: Dyall §9 jj-coupled full-Gaunt (see jj_angular_Xk).
    Returns a sparse dict keyed by (a, b, c, d) block-local spin-orbital
    indices.
    """
    Q = len(spinor_labels)
    eri: Dict[Tuple[int, int, int, int], float] = {}

    # Precompute X_k tables and per-(a,c) nonzero-k sets
    from collections import defaultdict
    ac_Xk: Dict[Tuple[int, int], List[Tuple[int, float]]] = defaultdict(list)

    # Determine k_max
    l_max = max(lab.l for lab in spinor_labels) if spinor_labels else 0
    k_max = 2 * l_max

    for a in range(Q):
        la = spinor_labels[a]
        for c in range(Q):
            lc = spinor_labels[c]
            for k in range(k_max + 1):
                val = jj_angular_Xk(la.kappa, la.two_m_j,
                                    lc.kappa, lc.two_m_j, k)
                if abs(val) > 1e-15:
                    ac_Xk[(a, c)].append((k, val))

    # Iterate 4-tuples using ac/bd grouping
    for (a, c), acs in ac_Xk.items():
        la = spinor_labels[a]
        lc = spinor_labels[c]
        for (b, d), bds in ac_Xk.items():
            lb = spinor_labels[b]
            ld = spinor_labels[d]
            # Global m-conservation
            if la.two_m_j + lb.two_m_j != lc.two_m_j + ld.two_m_j:
                continue
            val = 0.0
            # Only contributions with k_ac == k_bd add
            for k_ac, x_ac in acs:
                for k_bd, x_bd in bds:
                    if k_ac != k_bd:
                        continue
                    k = k_ac
                    rk_key = (la.n_fock, la.l, lb.n_fock, lb.l,
                              lc.n_fock, lc.l, ld.n_fock, ld.l, k)
                    rk = rk_cache.get(rk_key)
                    if rk is None:
                        continue
                    val += x_ac * x_bd * rk
            if abs(val) > 1e-14:
                eri[(a, b, c, d)] = val

    return eri


# ---------------------------------------------------------------------------
# Main relativistic builder
# ---------------------------------------------------------------------------

def build_composed_hamiltonian_relativistic(
    spec,
    alpha_num: float = 7.2973525693e-3,  # CODATA α
    pk_in_hamiltonian: bool = True,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Build the spin-ful composed qubit Hamiltonian for ``spec``.

    Each block in ``spec.blocks`` is promoted from (n, l, m) scalar orbitals
    to (n, κ, m_j) Dirac spinors.  Block-diagonal ERI structure is preserved:
    no cross-block two-body terms are generated.  h₁ is:

      h₁(lab) = −Z²/(2 n²)           (scalar kinetic, α→0 limit)
                + h_PK(lab, lab')    (diagonal in (l, m_j); PK barrier)
                + H_SO(n, κ, Z, α)   (T2 Breit-Pauli, diagonal)

    The one-body h₁ has off-diagonal elements only from PK (same (l, m_j),
    different n).  PK's angular selection is δ(l,l')·δ(m_j,m_j') in the
    spinor basis — we enforce this by grouping labels by (κ, m_j) and
    reusing the scalar PK matrix per (l, m_ℓ) sector for each (κ, m_j).

    Returns dict with keys: M (spinor orbital count per block — here the
    block-level spinor count), Q (total qubits), N_pauli, qubit_op,
    fermion_op, h1_diag (float array), eri_sparse (dict), 1-norm metrics,
    wall_time_s.

    Notes
    -----
    * ``alpha_num`` is the numerical fine-structure constant.  Passing
      α = 0 makes the SO term vanish (non-rel limit); useful for regression
      tests that check scalar consistency.
    * The builder returns Pauli counts EXCLUDING the identity (to match
      the scalar ``build_composed_hamiltonian`` ``N_pauli`` convention).
    """
    t0 = time.perf_counter()

    # ---- Phase 1: enumerate spinor labels per sub-block ---------------
    sub_blocks: List[Dict[str, Any]] = []
    q_offset = 0
    for blk in spec.blocks:
        l_min = getattr(blk, 'l_min', 0)
        center_labels = enumerate_dirac_labels(blk.max_n, l_min=l_min)
        sub_blocks.append({
            'label': blk.label + '_center',
            'Z': blk.Z_center,
            'dirac_labels': center_labels,
            'offset': q_offset,
            'n_orbitals': len(center_labels),
            'block': blk,
            'side': 'center',
        })
        q_offset += len(center_labels)

        if blk.has_h_partner:
            partner_max_n = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_labels = enumerate_dirac_labels(partner_max_n, l_min=0)
            sub_blocks.append({
                'label': blk.label + '_partner',
                'Z': blk.Z_partner,
                'dirac_labels': partner_labels,
                'offset': q_offset,
                'n_orbitals': len(partner_labels),
                'block': blk,
                'side': 'partner',
            })
            q_offset += len(partner_labels)

    Q = q_offset  # total qubits = total spinor spin-orbitals (JW 1:1)
    if verbose:
        print(f"[build_composed_rel] {spec.name}: Q={Q} spinor orbitals")
        for sb in sub_blocks:
            print(f"  {sb['label']}: {sb['n_orbitals']} orbs (Z={sb['Z']:.1f})")

    # ---- Phase 2: diagonal one-body h1 --------------------------------
    h1 = np.zeros((Q, Q))
    h1_so = np.zeros(Q)  # track SO separately for 1-norm reporting

    for sb in sub_blocks:
        Z_sb = float(sb['Z'])
        off = sb['offset']
        for i, lab in enumerate(sb['dirac_labels']):
            # Non-rel kinetic-plus-Coulomb: -Z²/(2n²) (same across κ branches)
            h1[off + i, off + i] = -Z_sb ** 2 / (2.0 * lab.n_fock ** 2)
            # Spin-orbit (T2), pure diagonal.  Evaluate α symbolically then
            # substitute numerical.
            so_expr = so_diagonal_matrix_element(
                lab.n_fock, lab.kappa, Z=Integer(int(round(Z_sb))),
                alpha=sp.Float(alpha_num))
            so_val = float(so_expr)
            h1[off + i, off + i] += so_val
            h1_so[off + i] = so_val

    # ---- Phase 2b: PK barrier on center orbitals ----------------------
    # Diagonal in (l, m_j), per-(κ,m_j) group shares the scalar (l, m_ℓ)
    # PK radial matrix.  We build the scalar PK once per (sb, l) group
    # using enumerated (n, l, m) scalar states and index into it by n.
    for sb in sub_blocks:
        blk = sb['block']
        if sb['side'] != 'center':
            continue
        if not (blk.pk_A > 0 and blk.pk_B > 0):
            continue
        # Scalar PK matrix for the full (n, l, m) enumeration at this block
        scalar_states = _enumerate_states(blk.max_n, l_min=getattr(blk, 'l_min', 0))
        pk_scalar = _compute_pk_matrix_elements(
            float(sb['Z']), scalar_states, blk.pk_A, blk.pk_B)
        # Index scalar states by (n, l, m) -> row
        scalar_idx = {(n, l, m): k for k, (n, l, m) in enumerate(scalar_states)}
        # Group Dirac labels by (l, two_m_j) — PK is diagonal in (l, m_j),
        # and within that group couples different n's.
        from collections import defaultdict
        lm_group: Dict[Tuple[int, int], List[Tuple[int, DiracLabel]]] = defaultdict(list)
        for i, lab in enumerate(sb['dirac_labels']):
            lm_group[(lab.l, lab.two_m_j)].append((i, lab))
        for (l, two_mj), entries in lm_group.items():
            # Pick any scalar m with |m| ≤ l to recover PK(n, l) values
            # (PK is m-independent because V_PK is spherically symmetric).
            scalar_m = 0 if l == 0 else 0
            for i_a, lab_a in entries:
                na = lab_a.n_fock
                key_a = (na, l, scalar_m)
                if key_a not in scalar_idx:
                    continue
                for i_b, lab_b in entries:
                    nb = lab_b.n_fock
                    key_b = (nb, l, scalar_m)
                    if key_b not in scalar_idx:
                        continue
                    val = pk_scalar[scalar_idx[key_a], scalar_idx[key_b]]
                    if abs(val) > 1e-14 and pk_in_hamiltonian:
                        off = sb['offset']
                        h1[off + i_a, off + i_b] += val

    # ---- Phase 3: ERI per block (block-diagonal) ----------------------
    eri_sparse: Dict[Tuple[int, int, int, int], float] = {}
    for sb in sub_blocks:
        Z_sb = float(sb['Z'])
        labs = sb['dirac_labels']
        off = sb['offset']
        if len(labs) == 0:
            continue
        # Scalar (n, l, m) states needed for radial integrals (m=0 proxy).
        unique_nl = sorted({(lab.n_fock, lab.l) for lab in labs})
        scalar_proxy = [(n, l, 0) for (n, l) in unique_nl]
        rk_cache = _compute_rk_integrals_block(Z_sb, scalar_proxy)
        if verbose:
            print(f"  ERI {sb['label']}: Q_sb={len(labs)}, R^k entries={len(rk_cache)}")
        phys = _build_spinor_eri_block(Z_sb, labs, rk_cache)
        for (a, b, c, d), val in phys.items():
            eri_sparse[(off + a, off + b, off + c, off + d)] = val

    # ---- Phase 4: FermionOperator assembly ----------------------------
    # Hermitize ERI: for a real Hermitian 1/r12, ⟨ab|V|cd⟩ = ⟨cd|V|ab⟩ and
    # our selection-rule traversal covers all 4-tuples, but individual
    # symbolic X_k evaluations can absorb sign phases asymmetrically.
    # Symmetrize explicitly so the resulting FermionOperator is Hermitian.
    sym_eri: Dict[Tuple[int, int, int, int], float] = {}
    for (a, b, c, d), val in eri_sparse.items():
        key = (a, b, c, d)
        alt = (c, d, a, b)
        sym_eri[key] = sym_eri.get(key, 0.0) + 0.5 * val
        sym_eri[alt] = sym_eri.get(alt, 0.0) + 0.5 * val

    fermion_op = FermionOperator((), float(spec.nuclear_repulsion_constant))
    # One-body (spinor -> single spin-orbital, no spin doubling).
    for p in range(Q):
        for q in range(Q):
            hv = h1[p, q]
            if abs(hv) < 1e-12:
                continue
            fermion_op += FermionOperator(((p, 1), (q, 0)), hv)
    # Two-body: physicist ⟨ab|cd⟩ → 0.5 a† b† d c (standard convention).
    for (a, b, c, d), val in sym_eri.items():
        # Pauli exclusion: a == b or c == d gives zero.
        if a == b or c == d:
            continue
        if abs(val) < 1e-14:
            continue
        fermion_op += FermionOperator(
            ((a, 1), (b, 1), (d, 0), (c, 0)),
            0.5 * val,
        )

    qubit_op = jordan_wigner(fermion_op)

    # ---- Metrics ------------------------------------------------------
    all_terms = qubit_op.terms
    N_pauli = len(all_terms) - (1 if () in all_terms else 0)
    # 1-norm
    lam_total = sum(abs(c) for c in all_terms.values())
    lam_ni = sum(abs(c) for t, c in all_terms.items() if t != ())
    # QWC groups
    from geovac.measurement_grouping import qwc_groups
    qwc = qwc_groups(qubit_op)

    # Block-diagonality sanity (no cross-block (a,b,c,d) should appear)
    block_of = [None] * Q
    for sb in sub_blocks:
        for i in range(sb['n_orbitals']):
            block_of[sb['offset'] + i] = sb['label']
    cross_block_count = sum(
        1 for (a, b, c, d) in eri_sparse.keys()
        if not (block_of[a] == block_of[b] == block_of[c] == block_of[d])
    )

    wall = time.perf_counter() - t0
    if verbose:
        print(f"  Q={Q}, N_pauli={N_pauli}, lambda_ni={lam_ni:.4f}, "
              f"QWC={len(qwc)}, wall={wall:.2f}s")
        print(f"  cross-block ERI entries (should be 0): {cross_block_count}")

    return {
        'Q': Q,
        'N_pauli': N_pauli,
        'lambda_total': float(lam_total),
        'lambda_ni': float(lam_ni),
        'qwc_groups': len(qwc),
        'qubit_op': qubit_op,
        'fermion_op': fermion_op,
        'h1_diag': np.diag(h1).copy(),
        'h1_so_diag': h1_so,
        'eri_sparse': eri_sparse,
        'cross_block_eri_count': cross_block_count,
        'nuclear_repulsion': float(spec.nuclear_repulsion_constant),
        'wall_time_s': wall,
        'spec_name': spec.name,
        'blocks': [
            {'label': sb['label'], 'n_orbitals': sb['n_orbitals'],
             'Z': float(sb['Z'])}
            for sb in sub_blocks
        ],
        'alpha_num': alpha_num,
    }
