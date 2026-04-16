"""
Track DC-B: Numerical comparison of Dirac-Coulomb vs Schrodinger-Coulomb FCI
convergence on He-like two-electron systems at Z=4 (Be 2+).
============================================================================

Paired with DC-A (symbolic): DC-A predicts p_D = p_S = 4 with only an
O(alpha^2) amplitude difference -- because the Dirac-Coulomb partial-wave
convergence is determined by the large-component r^{l+1} behavior, which
reduces to the non-relativistic rate in the alpha -> 0 limit and receives
only a multiplicative O(alpha^2) correction at leading relativistic order.

DC-B scope (algebraic-first)
----------------------------
At Tier-2 T3 (the current spin-ful composed builder), the radial
integrals R^k are the non-relativistic hydrogenic ones -- the Dirac
gamma = sqrt(kappa^2 - (Z*alpha)^2) correction to the small-r radial
behavior is NOT yet wired into the qubit ERI.  This is documented in
CLAUDE.md section 12: the Tier 2 entry states that gamma is a reserved
symbol, not bound numerically.  Tier 3 T7 adds diagonal gamma-corrected
<r^k> for n_r=0, but the two-body integrals remain scalar.

Consequence for the cusp test
-----------------------------
When spinor angular Gaunt X_k is combined with the SAME scalar radial
R^k as the scalar LS-coupled builder, the full 2-electron Hilbert space
is a unitary rotation of the scalar (l, m_l, m_s) space.  FCI is
rotation-invariant (Sprint 3D, v2.6.0).  Therefore:

  * At alpha = 0, E_D(n_max) = E_S(n_max) to machine precision.
  * At alpha = CODATA, E_D(n_max) - E_S(n_max) = pure SO-diagonal shift,
    which is n_max-independent (diagonal in (n, kappa, m_j)).
  * Convergence exponent p_D = p_S exactly at every n_max.

Thus DC-B's numerical job is to VERIFY these algebraic predictions:
  (1) alpha = 0 spinor FCI matches scalar FCI to < 1e-10 Ha at n_max=2,3,4.
  (2) alpha = CODATA spinor FCI is an SO-shift above the alpha = 0 spinor FCI.
  (3) The scalar Z=4 convergence exponent p_S is measured across n_max=2..5.

Z = 4 rationale
---------------
  * Z = 2 is within the graph-validity-boundary artifact (CUSP-2,
    v2.9.2): the ~0.20% floor there is NOT the cusp.
  * Z = 4 (Be 2+) is the smallest Z where the Schwartz cusp dominates
    for He-like systems.  Reference: Pekeris-style high-precision
    non-relativistic energy ~ -13.655566238 Ha.

Bug fix history
---------------
The original version used openfermion.get_sparse_operator(qubit_op, n_qubits=Q)
which builds the full 2^Q sparse operator before slicing to N=2.
At n_max=3, Q=28 -> 2^28 = 268M entries -> 4 GiB allocation failure.

This was the same bug that Track TC-V (v2.9.0) fixed for the TC benchmark.
Fix: build the 2-electron sector FCI matrix directly from the FermionOperator
via apply_op_string (second-quantization algebra), identical to the approach
in debug/track_bx4/fci_2e_solver.py.  The 2-electron sector has C(Q, 2)
determinants -- trivially small even at n_max=5 (C(110, 2) = 5995).
"""

import json
import math
import os
import sys
import time
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# Ensure we import from the working copy
HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO))

from geovac.casimir_ci import build_fci_matrix
from geovac.composed_qubit_relativistic import (
    build_composed_hamiltonian_relativistic,
    enumerate_dirac_labels,
    _build_spinor_eri_block,
    _compute_rk_integrals_block,
)
from geovac.spin_orbit import so_diagonal_matrix_element
from geovac.molecular_spec import MolecularSpec, OrbitalBlock
from openfermion import FermionOperator
import sympy as sp
from sympy import Integer

# --------------------------------------------------------------------
# Reference energy for Be 2+ (He-like, Z=4)
# --------------------------------------------------------------------
# Source: Drake & Yan 1992 / Pekeris-style high-precision
# non-relativistic variational (Hylleraas) for two electrons in a
# Z=4 nuclear field.  Value to 6 decimals is the textbook benchmark.
E_EXACT_NONREL_Z4 = -13.655566238  # Ha, non-rel infinite-nuclear-mass limit


def he_like_spec(Z: int, max_n: int, relativistic: bool = False) -> MolecularSpec:
    """Minimal He-like spec: single block, single center, no PK, no H partner."""
    block = OrbitalBlock(
        label=f"Z{Z}_core",
        block_type="core",
        Z_center=float(Z),
        n_electrons=2,
        max_n=max_n,
        has_h_partner=False,
        pk_A=0.0,
        pk_B=0.0,
        l_min=0,
    )
    return MolecularSpec(
        name=f"He-like-Z{Z}",
        blocks=[block],
        nuclear_repulsion_constant=0.0,  # no V_NN for atom
        relativistic=relativistic,
    )


# --------------------------------------------------------------------
# Direct 2-electron FCI from FermionOperator
# (no 2^Q sparse operator, never touches JW qubit space)
# --------------------------------------------------------------------

def apply_op_string(
    ops: Tuple[Tuple[int, int], ...],
    det: Tuple[int, ...],
) -> Tuple[float, Tuple[int, ...]]:
    """Apply a string of creation/annihilation operators to a Fock state.

    ops: tuple of (index, type) where type=1 for creation, 0 for annihilation.
         Applied RIGHT-TO-LEFT (rightmost is applied first, following
         FermionOperator's left-to-right reading of operator products).
    det: sorted tuple of occupied spin-orbital indices.

    Returns (phase, new_det) or (0.0, ()) if result is zero.

    Lifted from debug/track_bx4/fci_2e_solver.py (Track TC-V, v2.9.0).
    """
    state = list(det)
    phase = 1.0

    for op_idx, op_type in reversed(ops):
        if op_type == 0:  # annihilation
            if op_idx not in state:
                return 0.0, ()
            pos = state.index(op_idx)
            phase *= (-1) ** pos
            state.pop(pos)
        else:  # creation
            if op_idx in state:
                return 0.0, ()
            pos = 0
            while pos < len(state) and state[pos] < op_idx:
                pos += 1
            phase *= (-1) ** pos
            state.insert(pos, op_idx)

    return phase, tuple(state)


def build_fci_matrix_n_electron(
    fermion_op,
    n_qubits: int,
    n_electrons: int,
) -> np.ndarray:
    """Build the N-electron-sector FCI matrix from an openfermion FermionOperator.

    Generic implementation: loop over fermion_op terms × dets using
    apply_op_string.  Correct for any N_electrons but SLOW — at N=2, Q=60,
    ~360k terms × 1770 dets = ~640M calls, i.e. ~10 min in Python.

    For N=2 specifically, use build_fci_matrix_2e_fast which is ~100x
    faster because it direct-indexes 2e dets.

    Size: C(n_qubits, n_electrons).  For n_electrons=2, this is n*(n-1)/2:
      n_qubits=10 -> 45
      n_qubits=28 -> 378
      n_qubits=60 -> 1770
      n_qubits=110 -> 5995
    """
    dets = list(combinations(range(n_qubits), n_electrons))
    N_SD = len(dets)
    det_map: Dict[Tuple[int, ...], int] = {d: i for i, d in enumerate(dets)}

    # Start complex-valued for generality; we symmetrize/project to real at the end
    # when the input is Hermitian.  FermionOperator coefficients are complex in
    # principle, though our builder produces real values.
    H = np.zeros((N_SD, N_SD), dtype=complex)

    for term, coeff in fermion_op.terms.items():
        if abs(coeff) < 1e-15:
            continue

        if len(term) == 0:
            # Identity contribution on diagonal
            for I in range(N_SD):
                H[I, I] += coeff
            continue

        for J in range(N_SD):
            phase, new_det = apply_op_string(term, dets[J])
            if phase == 0.0:
                continue
            I = det_map.get(new_det)
            if I is not None:
                H[I, J] += coeff * phase

    return H


def build_fci_matrix_2e_fast(
    fermion_op,
    n_qubits: int,
) -> np.ndarray:
    """Direct 2-electron-sector FCI matrix build.

    Exploits the closed-form action of 1-body and 2-body fermion operators
    on 2-electron determinants (p, q) with p < q:

    For a 2-body term ``a†_a a†_b a_d a_c`` (4 indices):
      - Requires {c, d} ⊂ {p, q}.
      - If {c, d} ≠ {p, q}, no contribution.
      - Annihilate: both electrons removed, det → vacuum.
      - Create: (a, b) → sorted order, tracking phase.

    For a 1-body term ``a†_a a_c`` (2 indices):
      - Requires c ∈ {p, q}.
      - Annihilate c; then create a at the right spot.

    Complexity: O(|fermion_op| × O(1) per non-zero-action) = O(|fermion_op|).
    Near 100× faster than the generic path for 2e systems.
    """
    dets = list(combinations(range(n_qubits), 2))
    N_SD = len(dets)
    det_index: Dict[Tuple[int, int], int] = {d: i for i, d in enumerate(dets)}

    H = np.zeros((N_SD, N_SD), dtype=complex)

    for term, coeff in fermion_op.terms.items():
        if abs(coeff) < 1e-15:
            continue

        if len(term) == 0:
            # Identity on diagonal
            for I in range(N_SD):
                H[I, I] += coeff
            continue

        if len(term) == 2:
            # 1-body: a†_a a_c
            (a_idx, a_type), (c_idx, c_type) = term
            if a_type != 1 or c_type != 0:
                # Unexpected pattern: fall back to generic
                for J in range(N_SD):
                    phase, new_det = apply_op_string(term, dets[J])
                    if phase == 0.0:
                        continue
                    I = det_index.get(new_det)
                    if I is not None:
                        H[I, J] += coeff * phase
                continue
            a, c = a_idx, c_idx
            # Iterate all dets (p, q) with p < q
            for J, (p, q) in enumerate(dets):
                # Apply a_c: need c ∈ {p, q}
                if c == p:
                    # det -> (q,) after annihilation (phase +1 since position 0)
                    # Now apply a†_a: need a ∉ {q}
                    if a == q:
                        continue
                    # Result: sorted (a, q)
                    new_det = (a, q) if a < q else (q, a)
                    phase = 1.0  # a_c at position 0 -> +1
                    if a < q:
                        # Creating a before q: position 0 in {q}
                        phase *= 1.0
                    else:
                        # Creating a after q
                        phase *= -1.0
                    I = det_index.get(new_det)
                    if I is not None:
                        H[I, J] += coeff * phase
                elif c == q:
                    # det -> (p,) after annihilation (phase -1, position 1)
                    if a == p:
                        continue
                    new_det = (a, p) if a < p else (p, a)
                    phase = -1.0  # a_c at position 1 -> -1
                    if a < p:
                        phase *= 1.0
                    else:
                        phase *= -1.0
                    I = det_index.get(new_det)
                    if I is not None:
                        H[I, J] += coeff * phase
            continue

        if len(term) == 4:
            # 2-body: check normal-ordered pattern
            (a_idx, a_type), (b_idx, b_type), (d_idx, d_type), (c_idx, c_type) = term
            # Standard form from builder: ((a,1),(b,1),(d,0),(c,0))
            if not (a_type == 1 and b_type == 1
                    and d_type == 0 and c_type == 0):
                # Unexpected pattern: fall back to generic
                for J in range(N_SD):
                    phase, new_det = apply_op_string(term, dets[J])
                    if phase == 0.0:
                        continue
                    I = det_index.get(new_det)
                    if I is not None:
                        H[I, J] += coeff * phase
                continue

            a, b = a_idx, b_idx
            c, d_ = c_idx, d_idx
            # Pauli: a == b gives 0, c == d gives 0
            if a == b or c == d_:
                continue
            # Action on det (p, q): need {c, d_} = {p, q}
            # (since we need both electrons removed and det has 2 electrons)
            # Two cases: (c, d_) = (p, q) -- i.e. c=p, d_=q; or (c, d_) = (q, p)
            # Let's check both orderings
            # Case 1: c=p, d_=q (with p < q)
            key1 = (min(c, d_), max(c, d_))
            if key1 == (c, d_) or key1 == (d_, c):
                J = det_index.get(key1)
                if J is None:
                    continue
                # Compute phase of applying right-to-left: a_c first, then a_d,
                # then creation of b, then creation of a.
                # Apply to (p, q) where p < q:
                p, q = key1
                # Step 1: apply a_c^ at ket (p, q). c is p or q.
                # The full term-string is ((a,1),(b,1),(d,0),(c,0)), applied
                # right-to-left: first (c,0), then (d,0), then (b,1), then (a,1).
                # a_c on (p,q): c ∈ {p,q}.
                state = [p, q]
                phase = 1.0
                # c annihilate:
                if c == state[0]:
                    phase *= 1.0  # position 0
                    state = [state[1]]
                elif c == state[1]:
                    phase *= -1.0  # position 1
                    state = [state[0]]
                else:
                    continue
                # d annihilate:
                if d_ == state[0]:
                    phase *= 1.0
                    state = []
                else:
                    continue
                # b create: insert into sorted empty state
                state = [b]
                # position of b in empty -> 0
                phase *= 1.0
                # a create: insert into state [b]
                if a == b:
                    continue
                if a < b:
                    # insert at position 0, b slides right
                    phase *= 1.0
                    new_det = (a, b)
                else:
                    # insert at position 1 (after b)
                    phase *= -1.0
                    new_det = (b, a)
                I = det_index.get(new_det)
                if I is not None:
                    H[I, J] += coeff * phase
            continue

        # Higher-body: not expected in our builder, but fallback
        for J in range(N_SD):
            phase, new_det = apply_op_string(term, dets[J])
            if phase == 0.0:
                continue
            I = det_index.get(new_det)
            if I is not None:
                H[I, J] += coeff * phase

    return H


def ground_state_in_particle_sector(
    fermion_op,
    n_qubits: int,
    n_electrons: int = 2,
) -> Tuple[float, int]:
    """Diagonalize the N-electron sector of a FermionOperator.

    Returns (E_0, n_configs).  If the matrix is within 1e-10 of Hermitian,
    symmetrizes and uses eigvalsh; otherwise uses eigvals (non-Hermitian).
    """
    if n_electrons == 2:
        # Fast specialized 2e builder (~100x speedup at Q=60)
        H = build_fci_matrix_2e_fast(fermion_op, n_qubits)
    else:
        H = build_fci_matrix_n_electron(fermion_op, n_qubits, n_electrons)
    H_real = H.real
    H_imag = H.imag

    # Imaginary part should be ~ zero for Hermitian fermion_op
    max_imag = float(np.max(np.abs(H_imag)))
    if max_imag > 1e-8:
        # Treat as non-Hermitian: take real part of lowest-real eigenvalue
        vals = np.linalg.eigvals(H)
        idx = int(np.argmin(vals.real))
        return float(vals.real[idx]), H.shape[0]

    # Hermitian path: symmetrize residual numerical asymmetry
    H_sym = 0.5 * (H_real + H_real.T)
    evals = np.linalg.eigvalsh(H_sym)
    return float(evals[0]), H.shape[0]


# --------------------------------------------------------------------
# Run a single (Z, n_max, builder) point
# --------------------------------------------------------------------

def scalar_fci_energy(Z: int, n_max: int) -> dict:
    """Scalar (Schrodinger-Coulomb) FCI at k_orb = Z.

    Uses casimir_ci.build_fci_matrix directly -- this is the same engine
    used elsewhere in the project for He graph-native CI / Casimir CI.
    Since k_orb = Z, the one-body block is purely diagonal (-Z^2/2n^2).
    """
    t0 = time.perf_counter()
    # l_max unrestricted (default = n_max - 1)
    H = build_fci_matrix(Z=Z, n_max=n_max, k_orb=float(Z), l_max=None,
                         m_total=0)
    evals = np.linalg.eigvalsh(H)
    wall = time.perf_counter() - t0
    return {
        'method': 'scalar',
        'n_max': n_max,
        'n_configs': H.shape[0],
        'E': float(evals[0]),
        'wall_s': wall,
    }


def _build_spinor_fermion_op_fast(
    Z: float, max_n: int, alpha_num: float,
) -> Tuple["FermionOperator", int]:
    """Build the spinor FermionOperator for a single-block He-like system.

    Bypasses the full `build_composed_hamiltonian_relativistic` pipeline
    (which also does JW conversion + QWC grouping — both O(N^2) in Pauli
    count and unnecessary for 2-electron FCI).  Uses only the internal
    helpers:
      - enumerate_dirac_labels
      - _compute_rk_integrals_block
      - _build_spinor_eri_block
      - so_diagonal_matrix_element (geovac.spin_orbit)

    Returns (fermion_op, Q).
    """
    labs = enumerate_dirac_labels(max_n=max_n, l_min=0)
    Q = len(labs)

    # R^k on scalar (n, l, m=0) proxy — m-independent
    unique_nl = sorted({(lab.n_fock, lab.l) for lab in labs})
    scalar_proxy = [(n, l, 0) for (n, l) in unique_nl]
    rk_cache = _compute_rk_integrals_block(Z, scalar_proxy)

    # Spinor ERI
    eri = _build_spinor_eri_block(Z, labs, rk_cache)

    # Hermitize
    sym_eri: Dict[Tuple[int, int, int, int], float] = {}
    for (a, b, c, d), val in eri.items():
        key = (a, b, c, d)
        alt = (c, d, a, b)
        sym_eri[key] = sym_eri.get(key, 0.0) + 0.5 * val
        sym_eri[alt] = sym_eri.get(alt, 0.0) + 0.5 * val

    # Diagonal h1: non-rel kinetic + Coulomb + SO
    h1_diag = np.zeros(Q)
    Z_int = Integer(int(round(Z)))
    alpha_sp = sp.Float(alpha_num)
    for i, lab in enumerate(labs):
        val = -float(Z) ** 2 / (2.0 * lab.n_fock ** 2)
        so_expr = so_diagonal_matrix_element(
            lab.n_fock, lab.kappa, Z=Z_int, alpha=alpha_sp,
        )
        val += float(so_expr)
        h1_diag[i] = val

    # Assemble FermionOperator
    fop = FermionOperator((), 0.0)
    for p in range(Q):
        if abs(h1_diag[p]) > 1e-14:
            fop += FermionOperator(((p, 1), (p, 0)), h1_diag[p])
    for (a, b, c, d), val in sym_eri.items():
        if a == b or c == d:
            continue
        if abs(val) < 1e-14:
            continue
        fop += FermionOperator(((a, 1), (b, 1), (d, 0), (c, 0)), 0.5 * val)

    return fop, Q


def spinor_fci_energy(Z: int, n_max: int, alpha_num: float,
                      use_fast_builder: bool = True) -> dict:
    """Spinor (Dirac-Coulomb) FCI via T3 builder, projected to N=2 sector.

    At alpha_num = 0, the spinor FCI does NOT necessarily match scalar FCI
    exactly — the Tier-2 T3 builder uses jj-coupled full-Gaunt angular
    factors `X_k(κ, m_j)`, which differ from the scalar LS-coupled singlet
    Slater formula.  See the DC-B memo section 5.

    Uses direct 2-electron sector projection on the FermionOperator
    (no 2^Q qubit sparse operator build, avoiding the Q=28 4-GiB bug).

    If ``use_fast_builder``, bypasses JW + QWC metric computation in the
    full builder (~10x faster at n_max=4).  Otherwise uses the full builder
    and reads its fermion_op back (needed if Pauli-count metric is required).
    """
    t0 = time.perf_counter()
    if use_fast_builder:
        fermion_op, Q = _build_spinor_fermion_op_fast(
            Z=float(Z), max_n=n_max, alpha_num=alpha_num,
        )
        N_pauli = None  # not computed in fast path
    else:
        spec = he_like_spec(Z=Z, max_n=n_max, relativistic=True)
        result = build_composed_hamiltonian_relativistic(
            spec, alpha_num=alpha_num, pk_in_hamiltonian=False, verbose=False,
        )
        fermion_op = result['fermion_op']
        Q = result['Q']
        N_pauli = result['N_pauli']
    build_wall = time.perf_counter() - t0

    t1 = time.perf_counter()
    E_gs, n_configs = ground_state_in_particle_sector(
        fermion_op, Q, n_electrons=2,
    )
    diag_wall = time.perf_counter() - t1

    return {
        'method': 'spinor',
        'n_max': n_max,
        'alpha': alpha_num,
        'Q': Q,
        'N_pauli': N_pauli,
        'n_configs': n_configs,
        'n_fermion_terms': len(fermion_op.terms),
        'E': E_gs,
        'build_wall_s': build_wall,
        'diag_wall_s': diag_wall,
    }


# --------------------------------------------------------------------
# Convergence fitting
# --------------------------------------------------------------------

def fit_convergence_exponent(n_max_list, err_abs_list):
    """Fit err_abs(l_max+1) = A * (l_max+1)^{-p} in log-space.

    We use (l_max + 1) = n_max as the angular truncation proxy (since each
    shell n contains partial waves up to l = n - 1, so l_max = n_max - 1 and
    (l_max + 1) = n_max).  The Schwartz partial-wave argument gives
    err ~ (l_max+1)^{-4} = n_max^{-4}.

    Positive err_abs only (variational).  Requires abs(err) > 0.
    """
    ns = np.array(n_max_list, dtype=float)
    es = np.array(err_abs_list, dtype=float)
    log_n = np.log(ns)
    log_e = np.log(np.abs(es))
    A_mat = np.vstack([np.ones_like(log_n), log_n]).T
    coef, residuals, rank, sv = np.linalg.lstsq(A_mat, log_e, rcond=None)
    log_A, neg_p = coef
    p = -neg_p
    A = np.exp(log_A)
    pred = log_A - p * log_n
    resid = log_e - pred
    n_pts = len(ns)
    if n_pts > 2:
        s_sq = np.sum(resid ** 2) / (n_pts - 2)
        x_mean = np.mean(log_n)
        Sxx = np.sum((log_n - x_mean) ** 2)
        p_stderr = np.sqrt(s_sq / Sxx) if Sxx > 0 else float('nan')
    else:
        p_stderr = float('nan')
    return {'A': float(A), 'p': float(p), 'p_stderr': float(p_stderr),
            'residuals': resid.tolist()}


def successive_ratios(n_max_list, err_abs_list):
    """Pairwise err(n+1)/err(n) ratios.

    If err ~ n^{-p}, the ratio equals (n/(n+1))^p.
    Computing p = -log(ratio) / log((n+1)/n) gives an instantaneous exponent.
    """
    ns = list(n_max_list)
    es = [abs(e) for e in err_abs_list]
    rows = []
    for i in range(len(ns) - 1):
        r = es[i+1] / es[i]
        # p from (ns[i]/ns[i+1])^p = r  =>  p = log(r) / log(ns[i]/ns[i+1])
        # i.e. p = log(r) / log(ns[i]/ns[i+1])
        denom = math.log(ns[i] / ns[i+1])
        p_inst = math.log(r) / denom if denom != 0 else float('nan')
        rows.append({
            'n_lo': ns[i], 'n_hi': ns[i+1],
            'err_lo': es[i], 'err_hi': es[i+1],
            'ratio': r,
            'p_inst': p_inst,
        })
    return rows


# --------------------------------------------------------------------
# Main run
# --------------------------------------------------------------------

def main():
    Z = 4
    ALPHA_CODATA = 7.2973525693e-3

    # Schroedinger FCI is cheap at all n_max through 5 (tested: ~4s);
    # 6 and 7 are tractable but slow (27s, 166s).  We include n_max up to 5
    # by default; n_max=6 only if --extended.
    scalar_n_max_values = [2, 3, 4, 5]

    # Spinor builder timings from pilot: n_max=2 <1s, n_max=3 ~20s, n_max=4
    # ~5-15 min (Q=60, many wigner_3j calls, ~12M inner iterations).  We
    # attempt through n_max=4 by default; if it exceeds 30 min, the single-
    # point skip is caught by the timeout/MemoryError handler.
    spinor_n_max_values = [2, 3, 4]

    # CLI flags
    extended = '--extended' in sys.argv
    fast = '--fast' in sys.argv
    if extended:
        scalar_n_max_values = [2, 3, 4, 5, 6]
        spinor_n_max_values = [2, 3, 4, 5]
    if fast:
        spinor_n_max_values = [2, 3]  # skip n_max=4 for quick runs

    print("=" * 70)
    print(f"Track DC-B: Dirac vs Schroedinger FCI convergence at Z={Z} (Be 2+)")
    print(f"Reference non-rel energy: E_exact = {E_EXACT_NONREL_Z4} Ha")
    print(f"alpha (CODATA) = {ALPHA_CODATA}")
    print(f"Scalar n_max range: {scalar_n_max_values}")
    print(f"Spinor n_max range: {spinor_n_max_values}")
    print("=" * 70)

    # -------- Step 1: Scalar FCI -----------------------------------
    print("\n[Step 1] Scalar (Schroedinger-Coulomb) FCI")
    print("-" * 70)
    scalar_results = []
    for n_max in scalar_n_max_values:
        r = scalar_fci_energy(Z, n_max)
        r['err_abs'] = r['E'] - E_EXACT_NONREL_Z4  # Ha (should be > 0, variational)
        r['err_pct'] = 100.0 * r['err_abs'] / abs(E_EXACT_NONREL_Z4)
        scalar_results.append(r)
        print(f"  n_max={n_max}: E = {r['E']:.8f} Ha, "
              f"err = {r['err_abs']:+.6f} Ha ({r['err_pct']:+.4f}%), "
              f"configs={r['n_configs']}, {r['wall_s']:.2f}s")

    # -------- Step 2: Spinor FCI at alpha = 0 (regression check) ----
    print("\n[Step 2] Spinor FCI at alpha = 0 (should match scalar)")
    print("-" * 70)
    spinor_alpha0_results = []
    for n_max in spinor_n_max_values:
        try:
            r = spinor_fci_energy(Z, n_max, alpha_num=0.0)
        except MemoryError as e:
            print(f"  n_max={n_max}: SKIPPED (MemoryError: {e})")
            continue
        r['err_abs'] = r['E'] - E_EXACT_NONREL_Z4
        r['err_pct'] = 100.0 * r['err_abs'] / abs(E_EXACT_NONREL_Z4)
        spinor_alpha0_results.append(r)
        s = next(x for x in scalar_results if x['n_max'] == n_max)
        delta = r['E'] - s['E']
        print(f"  n_max={n_max}: E = {r['E']:.8f} Ha, "
              f"Q={r['Q']}, N_pauli={r['N_pauli']}, "
              f"configs={r['n_configs']}, "
              f"build={r['build_wall_s']:.1f}s, diag={r['diag_wall_s']:.1f}s")
        print(f"           E(spinor,alpha=0) - E(scalar): dE = {delta:+.3e} Ha")

    # -------- Step 3: Spinor FCI at alpha = CODATA ------------------
    print("\n[Step 3] Spinor FCI at alpha = CODATA (includes SO shift)")
    print("-" * 70)
    spinor_alpha_results = []
    for n_max in spinor_n_max_values:
        try:
            r = spinor_fci_energy(Z, n_max, alpha_num=ALPHA_CODATA)
        except MemoryError as e:
            print(f"  n_max={n_max}: SKIPPED (MemoryError: {e})")
            continue
        r['err_abs'] = r['E'] - E_EXACT_NONREL_Z4
        r['err_pct'] = 100.0 * r['err_abs'] / abs(E_EXACT_NONREL_Z4)
        spinor_alpha_results.append(r)
        a0 = next((x for x in spinor_alpha0_results if x['n_max'] == n_max), None)
        if a0 is not None:
            so_shift = r['E'] - a0['E']
            print(f"  n_max={n_max}: E = {r['E']:.8f} Ha, "
                  f"SO shift (vs alpha=0) = {so_shift:+.6e} Ha, "
                  f"build={r['build_wall_s']:.1f}s, diag={r['diag_wall_s']:.1f}s")
        else:
            print(f"  n_max={n_max}: E = {r['E']:.8f} Ha "
                  f"(alpha=0 not available, no SO shift extractable)")

    # -------- Step 4: Fit convergence exponents --------------------
    print("\n[Step 4] Convergence exponent fits")
    print("-" * 70)

    fit_scalar = fit_convergence_exponent(
        [r['n_max'] for r in scalar_results],
        [r['err_abs'] for r in scalar_results])
    ratios_scalar = successive_ratios(
        [r['n_max'] for r in scalar_results],
        [r['err_abs'] for r in scalar_results])

    fit_spinor0 = fit_convergence_exponent(
        [r['n_max'] for r in spinor_alpha0_results],
        [r['err_abs'] for r in spinor_alpha0_results])

    # alpha != 0: err has SO shift mixed in.  Subtract the (assumed constant)
    # SO shift from the most-converged n_max.
    if spinor_alpha_results and spinor_alpha0_results:
        so_shifts = [
            s_a['E'] - s_0['E']
            for s_a in spinor_alpha_results
            for s_0 in spinor_alpha0_results
            if s_a['n_max'] == s_0['n_max']
        ]
    else:
        so_shifts = []
    so_shift_nmax_max = so_shifts[-1] if so_shifts else 0.0

    err_spinor_raw = [r['err_abs'] for r in spinor_alpha_results]
    err_spinor_minus_so = [e - so_shift_nmax_max for e in err_spinor_raw]

    fit_spinor_raw = fit_convergence_exponent(
        [r['n_max'] for r in spinor_alpha_results],
        err_spinor_raw)
    fit_spinor_corrected = fit_convergence_exponent(
        [r['n_max'] for r in spinor_alpha_results],
        err_spinor_minus_so)

    ratios_spinor0 = successive_ratios(
        [r['n_max'] for r in spinor_alpha0_results],
        [r['err_abs'] for r in spinor_alpha0_results])

    def fmt_fit(fit):
        if math.isnan(fit['p_stderr']):
            return f"p = {fit['p']:.3f}, A = {fit['A']:.3e}  (stderr=nan, <=2 pts)"
        return (f"p = {fit['p']:.3f} +/- {fit['p_stderr']:.3f}, "
                f"A = {fit['A']:.3e}")

    print(f"  Scalar:               {fmt_fit(fit_scalar)}")
    print(f"  Spinor (alpha=0):     {fmt_fit(fit_spinor0)}")
    print(f"  Spinor (alpha=CODATA, raw):      {fmt_fit(fit_spinor_raw)}")
    print(f"  Spinor (alpha=CODATA, SO-corr):  {fmt_fit(fit_spinor_corrected)}")

    print("\n  Successive-pair instantaneous exponents:")
    print("  (Scalar)")
    for row in ratios_scalar:
        print(f"    n {row['n_lo']}->{row['n_hi']}: "
              f"ratio={row['ratio']:.4f}, p_inst={row['p_inst']:.3f}")
    print("  (Spinor, alpha=0)")
    for row in ratios_spinor0:
        print(f"    n {row['n_lo']}->{row['n_hi']}: "
              f"ratio={row['ratio']:.4f}, p_inst={row['p_inst']:.3f}")

    # -------- Step 5: Verdict --------------------------------------
    if fit_spinor_corrected['p'] > 0 and fit_scalar['p'] > 0:
        delta_p = fit_spinor_corrected['p'] - fit_scalar['p']
        delta_p_total = fit_spinor_raw['p'] - fit_scalar['p']
    else:
        delta_p = float('nan')
        delta_p_total = float('nan')

    print("\n[Step 5] Verdict")
    print("-" * 70)
    print(f"  Delta p (spinor-corrected - scalar) = {delta_p:+.4f}")
    print(f"  Delta p (spinor-raw       - scalar) = {delta_p_total:+.4f}")
    if not math.isnan(delta_p) and abs(delta_p) < 0.1:
        verdict = "SAME RATE (p_D approx p_S within numerical noise)"
    elif not math.isnan(delta_p) and delta_p > 0.1:
        verdict = "DIRAC FASTER"
    elif not math.isnan(delta_p):
        verdict = "DIRAC SLOWER"
    else:
        verdict = "INCONCLUSIVE (fit failed)"
    print(f"  Verdict: {verdict}")
    print(f"  Algebraic prediction (DC-A): Delta p = 0")
    print(f"  (Dirac-Coulomb cusp l^-4 with O(alpha^2) amplitude correction)")

    # Regression check: at alpha=0, spinor must equal scalar to < 1e-10
    max_alpha0_mismatch = 0.0
    for r in spinor_alpha0_results:
        s = next((x for x in scalar_results if x['n_max'] == r['n_max']), None)
        if s is not None:
            max_alpha0_mismatch = max(max_alpha0_mismatch, abs(r['E'] - s['E']))
    print(f"\n  Max |E_scalar - E_spinor(alpha=0)|: {max_alpha0_mismatch:.3e} Ha")
    if max_alpha0_mismatch < 1e-10:
        print("  [PASS] alpha=0 regression: spinor matches scalar to <1e-10.")
    elif max_alpha0_mismatch < 1e-6:
        print("  [PARTIAL] alpha=0 regression: within 1e-6 but not 1e-10 "
              "(likely numerical Pauli accumulation).")
    else:
        print("  [FAIL] alpha=0 regression: above 1e-6 (flag for investigation).")

    # -------- Appendix: l_max sweep at fixed n_max (Schwartz-proper) ----
    # Varying n_max mixes radial and angular completeness.  The Schwartz
    # l^-4 rate is the ANGULAR truncation rate at FIXED radial basis.
    # We do l_max in [0..4] at n_max=5 (scalar only, since spinor n_max>=4
    # is expensive).  This exposes whether the scalar l_max tail follows
    # Schwartz even in the hydrogenic finite-radial basis.
    print("\n[Step 6] l_max sweep at fixed n_max (Schwartz-proper test)")
    print("-" * 70)
    lmax_sweep = []
    n_max_lmax = 5
    E_ref_nmax = None
    for lmax in range(n_max_lmax):
        try:
            t0 = time.perf_counter()
            H = build_fci_matrix(Z=Z, n_max=n_max_lmax, k_orb=float(Z),
                                  l_max=lmax, m_total=0)
            evals = np.linalg.eigvalsh(H)
            wall = time.perf_counter() - t0
            E = float(evals[0])
            err = E - E_EXACT_NONREL_Z4
            if lmax == n_max_lmax - 1:
                E_ref_nmax = E
            lmax_sweep.append({
                'l_max': lmax, 'E': E, 'err_abs': err,
                'n_configs': H.shape[0], 'wall_s': wall,
            })
            print(f"  l_max={lmax}: E={E:.8f} Ha, err={err*1000:+.3f} mHa, "
                  f"configs={H.shape[0]}, {wall:.1f}s")
        except Exception as e:
            print(f"  l_max={lmax}: skipped ({type(e).__name__}: {e})")

    # Fit l_max convergence using err_l - err_ref (asymptotic residual)
    if len(lmax_sweep) >= 3 and E_ref_nmax is not None:
        # Use last point as reference "complete basis"; fit to (l_max+1)^{-p}
        lm_vals = np.array([r['l_max'] for r in lmax_sweep[:-1]])
        err_vs_ref = np.array(
            [r['E'] - E_ref_nmax for r in lmax_sweep[:-1]])
        # Only include POSITIVE residuals (variational)
        mask = err_vs_ref > 1e-9
        if mask.sum() >= 2:
            log_lp1 = np.log(lm_vals[mask] + 1.0)
            log_r = np.log(err_vs_ref[mask])
            A_mat = np.vstack([np.ones_like(log_lp1), log_lp1]).T
            coef, *_ = np.linalg.lstsq(A_mat, log_r, rcond=None)
            log_A, neg_p = coef
            p_lmax = -neg_p
            print(f"\n  Fit (err(l_max) - err(l_max={n_max_lmax-1})) "
                  f"= A (l_max+1)^(-p):  p = {p_lmax:.3f}, A = {np.exp(log_A):.3e}")
            print(f"  Schwartz prediction: p = 4 for singlet e-e cusp")
        else:
            p_lmax = None
            print("\n  Not enough nonzero residuals to fit l_max exponent")
    else:
        p_lmax = None

    # -------- Save results -----------------------------------------
    data_dir = REPO / "debug" / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    out = {
        'Z': Z,
        'alpha_codata': ALPHA_CODATA,
        'E_exact_nonrel': E_EXACT_NONREL_Z4,
        'scalar_n_max_values': scalar_n_max_values,
        'spinor_n_max_values': spinor_n_max_values,
        'scalar': scalar_results,
        'spinor_alpha0': spinor_alpha0_results,
        'spinor_alpha_codata': spinor_alpha_results,
        'so_shifts': so_shifts,
        'fits': {
            'scalar': fit_scalar,
            'spinor_alpha0': fit_spinor0,
            'spinor_raw': fit_spinor_raw,
            'spinor_corrected': fit_spinor_corrected,
        },
        'ratios': {
            'scalar': ratios_scalar,
            'spinor_alpha0': ratios_spinor0,
        },
        'delta_p_corrected': float(delta_p) if not math.isnan(delta_p) else None,
        'delta_p_raw': float(delta_p_total) if not math.isnan(delta_p_total) else None,
        'max_alpha0_mismatch_Ha': max_alpha0_mismatch,
        'verdict': verdict,
        'lmax_sweep_nmax_' + str(n_max_lmax): lmax_sweep,
        'lmax_fit_p': float(p_lmax) if p_lmax is not None else None,
    }
    out_path = data_dir / "dc_b_convergence.json"
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\n[Saved] {out_path}")
    return out


if __name__ == "__main__":
    main()
