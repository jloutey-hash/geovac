"""Sprint Z2 tapering — apply Hopf-U(1) m -> -m symmetry tapering to all
28 composed-library molecules.

Mechanism
---------
1. **Particle-number / spin-z Z2's.**  Jordan-Wigner of a number-conserving,
   S_z-conserving fermionic Hamiltonian has two free Z-string stabilizers,
   the parity of N_alpha and N_beta.  These reduce qubit count by 2
   (the standard Bravyi-Kitaev tapering result).

2. **Hopf-U(1) m -> -m Z2.**  Paper 29 §5.3 identified the m_l -> -m_l
   reflection as the only sub-action of the Paper 25 Hopf U(1) that
   commutes with a real-integer adjacency.  Under m-conservation
   (m1+m2 = m3+m4) the Hamiltonian commutes with P : m -> -m.

   In the native (n, l, m) basis P is a permutation of spatial orbitals,
   NOT a Pauli Z-string.  To use it for tapering we first rotate to the
   real symmetric / antisymmetric combinations

       phi_+(n,l,|m|)  = [phi_{n,l,m} + phi_{n,l,-m}] / sqrt(2)   (P = +1)
       phi_-(n,l,|m|)  = [phi_{n,l,m} - phi_{n,l,-m}] / sqrt(2)   (P = -1)

   for m > 0.  Orbitals with m = 0 are P-fixed (P = +1).

   In this basis P is a number operator, P = (-1)^{N_-} where N_- is
   the electron count on antisymmetric spatial orbitals.  Under JW
   this becomes the Z-string

       P = product_{q in antisym spin-orbitals} Z_q .

This driver:

  - rebuilds h1/eri/Hcore in the symmetric / antisymmetric basis,
  - constructs the three Z-string stabilizers,
  - calls openfermion.transforms.taper_off_qubits,
  - compares ground-state energies before/after for systems small enough
    to diagonalise directly (Q_tapered <= 22),
  - writes a JSON panel and a memo-ready table for all 28 systems.

Usage::

    python debug/sprint_z2_tapering.py [--gs-cutoff Q]

Author: GeoVac PM (sub-agent)
Date: 2026-06-04
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse.linalg as spla

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from openfermion import QubitOperator, jordan_wigner  # noqa: E402
from openfermion.linalg import get_sparse_operator  # noqa: E402
from openfermion.transforms import taper_off_qubits  # noqa: E402
from openfermion.utils import commutator  # noqa: E402

from geovac.qubit_encoding import build_fermion_op_from_integrals  # noqa: E402
from geovac.composed_qubit import _enumerate_states, build_composed_hamiltonian  # noqa: E402

# ----------------------------------------------------------------------------
# Build / rebuild Hamiltonian via molecular_spec
# ----------------------------------------------------------------------------

# Map: ecosystem-export canonical name -> spec factory key (callable on
# geovac.molecular_spec, or a special builder).  Mirrors
# geovac.ecosystem_export._SYSTEM_REGISTRY.

HYDRIDE_Z = {
    'LiH': 3, 'BeH2': 4, 'CH4': 6, 'NH3': 7, 'H2O': 8, 'HF': 9,
    'NaH': 11, 'MgH2': 12, 'SiH4': 14, 'PH3': 15, 'H2S': 16, 'HCl': 17,
    'KH': 19, 'CaH2': 20,
    'GeH4': 32, 'AsH3': 33, 'H2Se': 34, 'HBr': 35,
}
MULTI_CENTER = {
    'LiF': 'lif_spec', 'CO': 'co_spec', 'N2': 'n2_spec',
    'F2': 'f2_spec', 'NaCl': 'nacl_spec',
    # CH2O, C2H2, C2H6 — listed in CLAUDE.md §1.5 ecosystem inventory but
    # not present in `geovac.molecular_spec` factories at this commit;
    # `ecosystem_export._SYSTEM_REGISTRY` does NOT include them either,
    # so they are excluded from the panel.  When the factories ship, add
    # them back here.
}
TM_HYDRIDE_Z = {
    'ScH': 21, 'TiH': 22, 'VH': 23, 'CrH': 24, 'MnH': 25,
    'FeH': 26, 'CoH': 27, 'NiH': 28, 'CuH': 29, 'ZnH': 30,
}

# Full 28-molecule library (Paper 14 inventory: He + H2 + 18 main-group
# hydrides + 8 multi-center + 10 TM hydrides = 38 candidates, but the
# Paper 14 "Pauli scaling N_Pauli = 11.10 * Q" census draws on 28 of these
# for the composed table.  We include every molecule available through
# ecosystem_export so the table is comprehensive.)

ALL_MOLECULES = (
    list(HYDRIDE_Z.keys()) + list(MULTI_CENTER.keys())
    + list(TM_HYDRIDE_Z.keys()) + ['SrH', 'BaH', 'H2']
)
# 'He' is single-center and has only s orbitals at max_n=2 -> no m>0
# pairs, so P is trivial; include it for completeness but mark tapering=0.
ALL_MOLECULES.append('He')


def build_orbital_table(
    name: str, max_n: int = 2, R: Optional[float] = None,
) -> Tuple[List[Tuple[Any, int, int, int]], np.ndarray, np.ndarray, float]:
    """Build the per-molecule integrals (h1, eri) in the NATIVE (n,l,m)
    basis, along with an orbital description list.

    Returns
    -------
    orbital_table : list of (sub_block_label, n, l, m)
    h1 : (M, M) ndarray (real symmetric)
    eri : (M, M, M, M) ndarray (chemist (pq|rs))
    nuc_rep : float
    """
    from geovac import molecular_spec as ms

    # Build the spec
    if name in HYDRIDE_Z:
        spec = ms.hydride_spec(HYDRIDE_Z[name], R=R, max_n=max_n)
    elif name in MULTI_CENTER:
        factory = getattr(ms, MULTI_CENTER[name])
        kwargs = {'max_n': max_n}
        if R is not None:
            kwargs['R'] = R
        spec = factory(**kwargs)
    elif name in TM_HYDRIDE_Z:
        spec = ms.transition_metal_hydride_spec(TM_HYDRIDE_Z[name], R=R)
    elif name in ('SrH', 'BaH'):
        if name == 'SrH':
            spec = ms.srh_spec(R=R, max_n=max_n)
        else:
            spec = ms.bah_spec(R=R, max_n=max_n)
    elif name == 'H2':
        from geovac.molecular_spec import MolecularSpec, OrbitalBlock
        Rh = R if R is not None else 1.4
        spec = MolecularSpec(
            name='H2',
            blocks=[OrbitalBlock(
                label='H2_bond', block_type='bond_pair', Z_center=1.0,
                n_electrons=2, max_n=max_n,
            )],
            nuclear_repulsion_constant=1.0 / Rh,
        )
    elif name == 'He':
        from geovac.molecular_spec import MolecularSpec, OrbitalBlock
        spec = MolecularSpec(
            name='He',
            blocks=[OrbitalBlock(
                label='He_core', block_type='atomic', Z_center=2.0,
                n_electrons=2, max_n=max_n,
            )],
            nuclear_repulsion_constant=0.0,
        )
    else:
        raise ValueError(f'unknown molecule {name}')

    # The build_composed_hamiltonian path returns h1, eri, nuc_rep in the
    # NATIVE (n,l,m) basis.  We re-use it to harvest those tensors and
    # the sub-block ordering.
    result = build_composed_hamiltonian(spec, pk_in_hamiltonian=True)

    # Re-enumerate the orbital_table the same way the builder does
    # (mirrors composed_qubit.build_composed_hamiltonian phase 1).
    # IMPORTANT: include block index so duplicate labels across atoms
    # (e.g. N2's two 'N_core' blocks) are distinguished.  Multi-center
    # specs may reuse block labels across nuclei; the sub-block key must
    # be (blk_idx, side) not (label, side).
    orbital_table: List[Tuple[Any, int, int, int]] = []
    for blk_idx, blk in enumerate(spec.blocks):
        l_min = getattr(blk, 'l_min', 0)
        sb_key_center = (blk_idx, blk.label, 'center')
        for (n, l, m) in _enumerate_states(blk.max_n, l_min=l_min):
            orbital_table.append((sb_key_center, n, l, m))
        if blk.has_h_partner:
            partner_max_n = (
                blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            )
            sb_key_partner = (blk_idx, blk.label, 'partner')
            for (n, l, m) in _enumerate_states(partner_max_n):
                orbital_table.append((sb_key_partner, n, l, m))

    h1 = result['h1']
    if 'h1_cross_block' in result and result['h1_cross_block'] is not None:
        # cross_block is added inside build_composed_hamiltonian already
        pass
    if result.get('h1_pk') is not None:
        # PK is added to h1 inside build when pk_in_hamiltonian=True
        pass
    eri = result['eri']
    nuc_rep = result['nuclear_repulsion']

    assert h1.shape == (len(orbital_table), len(orbital_table)), (
        f'{name}: M={len(orbital_table)} != h1.shape={h1.shape}'
    )

    return orbital_table, h1, eri, nuc_rep


# ----------------------------------------------------------------------------
# Symmetric / antisymmetric basis rotation
# ----------------------------------------------------------------------------

def build_pm_rotation(
    orbital_table: List[Tuple[Any, int, int, int]],
) -> Tuple[np.ndarray, np.ndarray]:
    """Construct the orthogonal rotation U that maps native orbitals
    to symmetric (P = +1) and antisymmetric (P = -1) combinations under
    the m -> -m reflection P.

    The action of the reflection is *per sub-block / per (n, l)*: it pairs
    (n, l, +m) with (n, l, -m) for m > 0 within the same sub-block.

    Returns
    -------
    U : (M, M) orthogonal matrix
        Rows = new basis (sym then antisym), Cols = native.
    parity : (M,) int array
        +1 for symmetric, -1 for antisymmetric, 0 if the orbital was not
        pairable (should never happen given the (n,l,m) packing).
    """
    M = len(orbital_table)
    U = np.zeros((M, M))
    parity = np.zeros(M, dtype=int)

    # Index orbitals so we can locate the (sub_block, n, l, -m) partner
    idx_of = {orb: i for i, orb in enumerate(orbital_table)}

    handled = [False] * M
    for i, (sb, n, l, m) in enumerate(orbital_table):
        if handled[i]:
            continue
        if m == 0:
            U[i, i] = 1.0
            parity[i] = +1
            handled[i] = True
            continue

        # Find the partner (sb, n, l, -m)
        partner_key = (sb, n, l, -m)
        if partner_key not in idx_of:
            # Should not happen for our (n,l,m) packings
            raise RuntimeError(
                f'orbital {orbital_table[i]} has no m -> -m partner'
            )
        j = idx_of[partner_key]
        if j == i:
            U[i, i] = 1.0
            parity[i] = +1
            handled[i] = True
            continue

        # Assign sym -> the row of the +m orbital, antisym -> -m row.
        # Concretely: if i is the +m orbital, place
        #     U[i, i] = U[i, j] = 1/sqrt(2)   (symmetric)
        #     U[j, i] = 1/sqrt(2),  U[j, j] = -1/sqrt(2)   (antisymmetric)
        # We use the convention that the row index with +m is +,
        # the one with -m is -.
        if m > 0:
            i_plus, i_minus = i, j
        else:
            i_plus, i_minus = j, i
        inv_sqrt2 = 1.0 / np.sqrt(2.0)
        U[i_plus, i_plus] = inv_sqrt2
        U[i_plus, i_minus] = inv_sqrt2
        U[i_minus, i_plus] = inv_sqrt2
        U[i_minus, i_minus] = -inv_sqrt2
        parity[i_plus] = +1
        parity[i_minus] = -1
        handled[i_plus] = True
        handled[i_minus] = True

    # Sanity: U U^T == I
    err = np.max(np.abs(U @ U.T - np.eye(M)))
    if err > 1e-12:
        raise RuntimeError(f'rotation not orthogonal, err={err:.2e}')

    return U, parity


def rotate_h1_eri(
    h1: np.ndarray, eri: np.ndarray, U: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply the basis rotation to h1 and eri.

    For h1 (one-body): h1' = U @ h1 @ U.T
    For eri (two-body, chemist (pq|rs)):
        eri'[p,q,r,s] = sum_{abcd} U[p,a] U[q,b] U[r,c] U[s,d] eri[a,b,c,d]
    """
    h1p = U @ h1 @ U.T
    # einsum for the 4-index tensor.  All real.
    erip = np.einsum('pa,qb,rc,sd,abcd->pqrs', U, U, U, U, eri,
                     optimize='optimal')
    return h1p, erip


# ----------------------------------------------------------------------------
# Stabilizer construction in Jordan-Wigner basis
# ----------------------------------------------------------------------------

def build_stabilizers(
    parity: np.ndarray,
) -> Tuple[QubitOperator, QubitOperator, QubitOperator]:
    """Build the three Z-string stabilizers:

      Z_alpha = product over alpha spin-orbitals of Z_q   (parity of N_alpha)
      Z_beta  = product over beta  spin-orbitals of Z_q   (parity of N_beta)
      P       = product over antisymmetric spin-orbitals (both spins) of Z_q

    JW convention is sp = 2 * p + sigma (sigma=0 alpha, sigma=1 beta).
    """
    M = len(parity)

    # alpha parity: qubits 0, 2, 4, ... (sigma = 0)
    Z_alpha = QubitOperator(())
    for p in range(M):
        Z_alpha *= QubitOperator(((2 * p, 'Z'),))
    # beta parity: qubits 1, 3, 5, ...
    Z_beta = QubitOperator(())
    for p in range(M):
        Z_beta *= QubitOperator(((2 * p + 1, 'Z'),))
    # m-reflection parity: antisymmetric spatial orbitals, both spins
    P_op = QubitOperator(())
    n_anti = int(np.sum(parity == -1))
    if n_anti == 0:
        # No antisymmetric orbital exists (e.g. He at max_n=2 only has s).
        # P is trivially +I; do NOT add it as a stabilizer.
        return Z_alpha, Z_beta, QubitOperator(())  # caller checks for empty
    for p in range(M):
        if parity[p] == -1:
            P_op *= QubitOperator(((2 * p, 'Z'),))
            P_op *= QubitOperator(((2 * p + 1, 'Z'),))

    return Z_alpha, Z_beta, P_op


def commutes(op_a: QubitOperator, op_b: QubitOperator,
             atol: float = 1e-10) -> bool:
    c = commutator(op_a, op_b)
    for coef in c.terms.values():
        if abs(coef) > atol:
            return False
    return True


# ----------------------------------------------------------------------------
# Spectrum (small systems only)
# ----------------------------------------------------------------------------

def ground_state_energy(
    qubit_op: QubitOperator, n_qubits: int, max_q: int = 18,
) -> Optional[float]:
    """Compute the ground-state energy of qubit_op as a sparse Hermitian
    matrix.  Returns None if n_qubits is too large.

    For n_qubits > max_q (default 18 -> 2^18 = 262k dim) we skip:
    the construction is a similarity transform on the rotated H_rot,
    so the spectrum is preserved by mathematics regardless of size.
    """
    if n_qubits <= 0:
        # constant operator
        return float(np.real(qubit_op.terms.get((), 0.0)))
    if n_qubits > max_q:
        return None
    H_sparse = get_sparse_operator(qubit_op, n_qubits=n_qubits)
    if H_sparse.shape[0] <= 1024:
        # dense diag for small ones
        E = np.linalg.eigvalsh(H_sparse.toarray())
        return float(E[0])
    else:
        # sparse Lanczos
        vals = spla.eigsh(H_sparse, k=1, which='SA', tol=1e-12)[0]
        return float(vals[0])


# ----------------------------------------------------------------------------
# Per-molecule tapering driver
# ----------------------------------------------------------------------------

def n_qubits_of(qubit_op: QubitOperator) -> int:
    """Largest qubit index + 1 in qubit_op (identity term returns 0)."""
    max_q = -1
    for term in qubit_op.terms:
        for q, _ in term:
            if q > max_q:
                max_q = q
    return max_q + 1


def n_pauli_nonid(qubit_op: QubitOperator) -> int:
    """Number of non-identity Pauli terms."""
    return len([t for t in qubit_op.terms if len(t) > 0])


def one_norm_nonid(qubit_op: QubitOperator) -> float:
    return float(sum(
        abs(c) for t, c in qubit_op.terms.items() if len(t) > 0
    ))


def taper_molecule(
    name: str,
    max_n: int = 2,
    R: Optional[float] = None,
    gs_cutoff: int = 22,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Run the full tapering pipeline for one molecule and return a result dict."""
    t0 = time.perf_counter()
    out: Dict[str, Any] = {'name': name, 'max_n': max_n, 'R': R}

    try:
        # 1. Build native (n,l,m) integrals
        orbital_table, h1_n, eri_n, nuc_rep = build_orbital_table(
            name, max_n=max_n, R=R,
        )
        M = len(orbital_table)
        Q_naive = 2 * M
        out['M'] = M
        out['Q_naive'] = Q_naive

        # 2. Naive JW Hamiltonian
        fop_native = build_fermion_op_from_integrals(h1_n, eri_n, nuc_rep)
        qop_native = jordan_wigner(fop_native)
        out['N_pauli_naive'] = n_pauli_nonid(qop_native)
        out['lambda_naive'] = one_norm_nonid(qop_native)

        # 3. Build P-eigenbasis rotation
        U, parity = build_pm_rotation(orbital_table)
        n_sym = int(np.sum(parity == +1))
        n_anti = int(np.sum(parity == -1))
        out['n_sym'] = n_sym
        out['n_anti'] = n_anti
        h1_r, eri_r = rotate_h1_eri(h1_n, eri_n, U)
        fop_rot = build_fermion_op_from_integrals(h1_r, eri_r, nuc_rep)
        qop_rot = jordan_wigner(fop_rot)
        out['N_pauli_rotated'] = n_pauli_nonid(qop_rot)
        out['lambda_rotated'] = one_norm_nonid(qop_rot)

        # 4. Stabilizers in JW basis
        Z_alpha, Z_beta, P_op = build_stabilizers(parity)

        # Always-present alpha/beta parity
        stabilizers = [Z_alpha, Z_beta]
        n_stab_z_only = 2
        # m-reflection P is a stabilizer only if (i) there is at least one
        # antisymmetric orbital, and (ii) it commutes with the rotated H.
        # In the (P-eigenbasis) Hamiltonian (ii) is true by construction,
        # but we test numerically for safety.
        has_P = (n_anti > 0) and (len(P_op.terms) > 0)
        if has_P:
            comm_with_P = commutes(qop_rot, P_op)
            out['P_commutes'] = bool(comm_with_P)
            if comm_with_P:
                stabilizers.append(P_op)
        else:
            out['P_commutes'] = False  # trivially: no antisym orbital exists

        out['n_stabilizers'] = len(stabilizers)

        # 5. Taper — try both signs for each stabilizer.  taper_off_qubits
        # always projects onto the +1 eigenspace; the physical sector may
        # differ.  We sweep all 2^n_stab sign choices, taper each, find the
        # minimum ground-state energy across sectors, and take that as the
        # tapered ground state.
        n_stab = len(stabilizers)
        sector_results: List[Dict[str, Any]] = []
        best = None
        for sec in range(2 ** n_stab):
            signs = [(+1 if (sec >> k) & 1 == 0 else -1)
                     for k in range(n_stab)]
            signed_stabs = [s * sgn for s, sgn in zip(stabilizers, signs)]
            qop_sec = taper_off_qubits(qop_rot, signed_stabs)
            Q_sec = n_qubits_of(qop_sec)
            sector_results.append({
                'signs': signs, 'Q': Q_sec,
                'N_pauli': n_pauli_nonid(qop_sec),
                'lambda': one_norm_nonid(qop_sec),
            })
            if Q_sec <= gs_cutoff:
                E_sec = ground_state_energy(qop_sec, Q_sec, max_q=gs_cutoff)
                sector_results[-1]['E'] = E_sec
                if E_sec is not None and (best is None or E_sec < best['E']):
                    best = {
                        'E': E_sec, 'Q': Q_sec, 'signs': signs,
                        'qop': qop_sec,
                    }
            else:
                sector_results[-1]['E'] = None

        # Tapered Hamiltonian Q / Pauli / 1-norm metrics: these are the
        # SAME across all sign sectors (the sign only changes individual
        # coefficient signs, not the term set or the qubit count).  Take
        # sector 0 as representative.
        rep = sector_results[0]
        Q_tapered = rep['Q']
        out['Q_tapered'] = Q_tapered
        out['delta_Q'] = Q_naive - Q_tapered
        out['N_pauli_tapered'] = rep['N_pauli']
        out['lambda_tapered'] = rep['lambda']
        out['sector_table'] = sector_results

        # 6. Spectrum check (only for small enough).  Pass gs_cutoff
        # through so the function does not silently fall back to its
        # internal default cap.
        if Q_naive <= gs_cutoff:
            E_native = ground_state_energy(qop_native, Q_naive, max_q=gs_cutoff)
            out['E_naive'] = E_native
        else:
            out['E_naive'] = None

        if best is not None:
            out['E_tapered'] = best['E']
            out['best_sector_signs'] = best['signs']
            if out['E_naive'] is not None:
                out['delta_E'] = abs(best['E'] - out['E_naive'])
                out['spectrum_preserved'] = bool(out['delta_E'] < 1e-8)
            else:
                out['delta_E'] = None
                out['spectrum_preserved'] = None
        else:
            out['E_tapered'] = None
            out['delta_E'] = None
            out['spectrum_preserved'] = None

        out['status'] = 'ok'
    except Exception as e:
        out['status'] = 'error'
        out['error'] = str(e)
        out['traceback'] = traceback.format_exc()
        if verbose:
            print(f'  ERROR on {name}: {e}')
    finally:
        out['wall_s'] = time.perf_counter() - t0
    return out


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gs-cutoff', type=int, default=22,
                    help='Max Q_naive for which to attempt full diagonalization '
                         '(default 22).')
    ap.add_argument('--max-n', type=int, default=2,
                    help='max_n for spec builders (default 2).')
    ap.add_argument('--molecules', nargs='+', default=None,
                    help='Subset of molecules; defaults to all 28+ in library.')
    ap.add_argument('--out', default='debug/data/sprint_z2_tapering.json')
    args = ap.parse_args()

    molecules = args.molecules if args.molecules else ALL_MOLECULES
    print(f'Running on {len(molecules)} systems with max_n={args.max_n}, '
          f'gs-cutoff={args.gs_cutoff}.')

    results: List[Dict[str, Any]] = []
    for i, name in enumerate(molecules):
        print(f'[{i+1:2d}/{len(molecules)}] {name} ...')
        res = taper_molecule(
            name, max_n=args.max_n, gs_cutoff=args.gs_cutoff, verbose=True,
        )
        if res['status'] == 'ok':
            de = res.get('delta_E')
            de_str = (f'|dE|={de:.2e}' if de is not None else
                      'dE=skip')
            print(
                f'  Q {res["Q_naive"]:3d}->{res["Q_tapered"]:3d} '
                f'(-{res["delta_Q"]}), '
                f'Pauli {res["N_pauli_naive"]:5d}->{res["N_pauli_tapered"]:5d}, '
                f'L {res["lambda_naive"]:.3f}->{res["lambda_tapered"]:.3f}, '
                f'P_comm={res["P_commutes"]}, {de_str}, '
                f'wall={res["wall_s"]:.1f}s'
            )
        else:
            print(f'  [ERROR] {res.get("error")}')
        results.append(res)

    # Write JSON
    out_path = PROJECT_ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f'\nWrote {out_path}')

    # Print summary table
    print('\nSummary (markdown):')
    print('| Mol | Q_naive | Q_tap | dQ | Pauli_naive | Pauli_tap | L_naive | L_tap | P_comm | |dE| |')
    print('|:----|:-------:|:-----:|:--:|------------:|----------:|--------:|------:|:------:|-----:|')
    for r in results:
        if r['status'] != 'ok':
            print(f"| {r['name']} | ERROR | | | | | | | | |")
            continue
        de = r.get('delta_E')
        if de is None:
            de_s = 'n/a'
        else:
            de_s = f'{de:.1e}'
        print(
            f"| {r['name']} | {r['Q_naive']} | {r['Q_tapered']} | "
            f"{r['delta_Q']} | {r['N_pauli_naive']} | {r['N_pauli_tapered']} | "
            f"{r['lambda_naive']:.3f} | {r['lambda_tapered']:.3f} | "
            f"{r['P_commutes']} | {de_s} |"
        )

    n_ok = sum(1 for r in results if r['status'] == 'ok')
    n_delta3 = sum(
        1 for r in results
        if r['status'] == 'ok' and r.get('delta_Q', 0) >= 3
    )
    n_delta2 = sum(
        1 for r in results
        if r['status'] == 'ok' and r.get('delta_Q', 0) >= 2
    )
    n_preserved = sum(
        1 for r in results
        if r['status'] == 'ok' and r.get('spectrum_preserved') is True
    )
    n_checked = sum(
        1 for r in results
        if r['status'] == 'ok' and r.get('spectrum_preserved') is not None
    )
    print(f'\n{n_ok}/{len(results)} ran without error.')
    print(f'{n_delta2} systems tapered >= 2 qubits.')
    print(f'{n_delta3} systems tapered >= 3 qubits (m-reflection applied).')
    print(f'{n_preserved}/{n_checked} systems passed spectrum-preservation gate.')

    return 0


if __name__ == '__main__':
    sys.exit(main())
