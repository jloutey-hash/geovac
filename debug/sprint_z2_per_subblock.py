"""Sprint Z2 per-sub-block tapering -- extend the global Hopf-U(1)
m -> -m tapering to one INDEPENDENT P_i per sub-block.

Background
----------
The global tapering sprint (debug/sprint_z2_tapering.py, 2026-06-04)
implements the global product P_global = prod_i P_i and tapers
3 qubits (alpha-parity, beta-parity, P_global).

This sprint replaces P_global with the per-sub-block reflections
P_i, one for each sub-block. By construction:

  * sub-blocks have disjoint qubit supports under JW,
    so all P_i commute pairwise: [P_i, P_j] = 0 trivially.
  * each P_i commutes with H_rot because the standard composed
    builder has ZERO cross-block ERIs (ERI tensor is filled
    only via eri[off+a, off+r, off+b, off+s] += val with all
    indices in the same sub-block) and the only inter-block
    coupling is via PK Gaussian barrier on h1 which is
    diagonal in (n, l).
  * Z_alpha, Z_beta and all P_i are Z-strings, all mutually
    commuting.

Total stabilizers = 2 + n_sub_blocks_with_antisym
qubit reduction = 2 + n_sub_blocks_with_antisym (vs. 3 in the
global sprint).

The global P_global = prod_i P_i is now LINEARLY DEPENDENT on
the per-sub-block P_i's, so we DO NOT include it -- we replace
it with the finer-grained set.

This driver:
  - reuses the rotation + JW pipeline of sprint_z2_tapering.py
  - constructs one P_i per sub-block (each is a Z-string on
    antisymmetric orbitals of that sub-block only)
  - verifies [P_i, H_rot] = 0 numerically for each i (audit
    discipline -- the block-index bug in the global sprint
    surfaced only on multi-center diatomics; per-sub-block may
    have analogous bugs that only surface on specific topologies)
  - sweeps 2^(2+n_sub_blocks) sign sectors and picks min energy
  - compares to E_naive at Q_naive <= gs_cutoff

Author: GeoVac per-sub-block tapering sprint (2026-06-04)
"""

from __future__ import annotations

import argparse
import json
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

# Reuse the existing global-sprint utilities -- exact same builder
# pipeline, same rotation U, so the per-sub-block extension only
# changes which Z-strings are passed as stabilizers.
import importlib.util  # noqa: E402
_spec = importlib.util.spec_from_file_location(
    "sprint_z2_tapering",
    str(Path(__file__).resolve().parent / "sprint_z2_tapering.py"),
)
_st = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_st)

build_orbital_table = _st.build_orbital_table
build_pm_rotation = _st.build_pm_rotation
rotate_h1_eri = _st.rotate_h1_eri
commutes = _st.commutes
ground_state_energy = _st.ground_state_energy
n_qubits_of = _st.n_qubits_of
n_pauli_nonid = _st.n_pauli_nonid
one_norm_nonid = _st.one_norm_nonid
ALL_MOLECULES = _st.ALL_MOLECULES


# ----------------------------------------------------------------------------
# Per-sub-block stabilizer construction
# ----------------------------------------------------------------------------

def build_per_subblock_stabilizers(
    orbital_table: List[Tuple[Any, int, int, int]],
    parity: np.ndarray,
) -> Tuple[QubitOperator, QubitOperator, List[QubitOperator], List[Any]]:
    """Build Z_alpha, Z_beta, and one P_i per sub-block.

    Each P_i is a Z-string over the antisymmetric spatial orbitals
    of sub-block i (both spins). Sub-blocks with zero antisymmetric
    orbitals contribute NO P_i (the corresponding reflection is
    identity, trivial stabilizer, do not include).

    Returns
    -------
    Z_alpha, Z_beta : QubitOperator
        Standard alpha/beta parity Z-strings.
    P_per_block : list of QubitOperator
        One per sub-block that has >= 1 antisymmetric orbital.
        Empty list if no sub-block has antisymmetric orbitals.
    sb_keys : list
        The sub-block key (typically (blk_idx, label, side)) for
        each entry of P_per_block -- for diagnostics only.
    """
    M = len(orbital_table)

    # alpha / beta parity
    Z_alpha = QubitOperator(())
    for p in range(M):
        Z_alpha *= QubitOperator(((2 * p, 'Z'),))
    Z_beta = QubitOperator(())
    for p in range(M):
        Z_beta *= QubitOperator(((2 * p + 1, 'Z'),))

    # Group orbital indices by sub-block, then collect antisymmetric ones
    from collections import OrderedDict
    sb_to_antisym: "OrderedDict[Any, List[int]]" = OrderedDict()
    for p, (sb, n, l, m) in enumerate(orbital_table):
        if sb not in sb_to_antisym:
            sb_to_antisym[sb] = []
        if parity[p] == -1:
            sb_to_antisym[sb].append(p)

    P_per_block: List[QubitOperator] = []
    sb_keys: List[Any] = []
    for sb, anti_indices in sb_to_antisym.items():
        if len(anti_indices) == 0:
            continue  # nothing to stabilize
        P_i = QubitOperator(())
        for p in anti_indices:
            P_i *= QubitOperator(((2 * p, 'Z'),))
            P_i *= QubitOperator(((2 * p + 1, 'Z'),))
        P_per_block.append(P_i)
        sb_keys.append(sb)

    return Z_alpha, Z_beta, P_per_block, sb_keys


def verify_pairwise_commute(ops: List[QubitOperator],
                            atol: float = 1e-10) -> Tuple[bool, List[Tuple[int, int]]]:
    """Test [op_i, op_j] = 0 for all i < j. Returns (all_commute, fail_pairs)."""
    fails: List[Tuple[int, int]] = []
    for i in range(len(ops)):
        for j in range(i + 1, len(ops)):
            if not commutes(ops[i], ops[j], atol=atol):
                fails.append((i, j))
    return (len(fails) == 0), fails


# ----------------------------------------------------------------------------
# Per-molecule tapering pipeline
# ----------------------------------------------------------------------------

def taper_molecule_per_subblock(
    name: str,
    max_n: int = 2,
    R: Optional[float] = None,
    gs_cutoff: int = 22,
    verbose: bool = True,
    # Hard cap on # sectors to sweep -- 2^14 = 16384 already prohibitive
    sector_cap_log2: int = 12,
) -> Dict[str, Any]:
    """Run the per-sub-block tapering pipeline for one molecule."""
    t0 = time.perf_counter()
    out: Dict[str, Any] = {'name': name, 'max_n': max_n, 'R': R}

    try:
        # 1. Native integrals (n,l,m) basis -- identical to global sprint
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

        # 3. P-eigenbasis rotation (same as global sprint)
        U, parity = build_pm_rotation(orbital_table)
        out['n_sym'] = int(np.sum(parity == +1))
        out['n_anti'] = int(np.sum(parity == -1))
        h1_r, eri_r = rotate_h1_eri(h1_n, eri_n, U)
        fop_rot = build_fermion_op_from_integrals(h1_r, eri_r, nuc_rep)
        qop_rot = jordan_wigner(fop_rot)
        out['N_pauli_rotated'] = n_pauli_nonid(qop_rot)
        out['lambda_rotated'] = one_norm_nonid(qop_rot)

        # 4. Per-sub-block stabilizers
        Z_alpha, Z_beta, P_per_block, sb_keys = build_per_subblock_stabilizers(
            orbital_table, parity,
        )
        out['n_sub_blocks'] = len(set(orb[0] for orb in orbital_table))
        out['n_sub_blocks_with_antisym'] = len(P_per_block)
        out['sub_block_keys'] = [str(k) for k in sb_keys]

        # AUDIT: every P_i must commute with H_rot.
        # By construction this is true (cross-block ERIs are zero in the
        # standard composed builder), but the audit catches structural
        # bugs analogous to the block-index bug in the global sprint.
        P_commutes_with_H: List[bool] = []
        P_commutator_norms: List[float] = []
        for P_i in P_per_block:
            c = commutator(qop_rot, P_i)
            max_coef = max(
                (abs(v) for v in c.terms.values()), default=0.0
            )
            P_commutes_with_H.append(max_coef < 1e-10)
            P_commutator_norms.append(float(max_coef))
        out['P_per_block_commutes_with_H'] = P_commutes_with_H
        out['P_per_block_commutator_max_coef'] = P_commutator_norms
        out['all_P_commute_with_H'] = all(P_commutes_with_H)

        # AUDIT: all stabilizers must mutually commute. Z_alpha, Z_beta
        # are global Z-strings; P_i are sub-block-restricted Z-strings.
        # All Z-strings on the same qubit set commute, so this is trivial,
        # but verify.
        all_stabs = [Z_alpha, Z_beta] + P_per_block
        all_commute, fail_pairs = verify_pairwise_commute(all_stabs)
        out['all_stabs_pairwise_commute'] = bool(all_commute)
        out['failing_stabilizer_pairs'] = fail_pairs

        # Filter out any P_i that fails the commute-with-H test (defensive
        # -- should never happen in standard composed builder).
        good_P = [P for P, ok in zip(P_per_block, P_commutes_with_H) if ok]
        bad_count = len(P_per_block) - len(good_P)
        out['n_P_dropped'] = bad_count

        stabilizers = [Z_alpha, Z_beta] + good_P
        n_stab = len(stabilizers)
        out['n_stabilizers'] = n_stab

        if n_stab > sector_cap_log2:
            out['status'] = 'skip_too_many_sectors'
            out['reason'] = f'n_stab={n_stab} > sector_cap_log2={sector_cap_log2}'
            out['wall_s'] = time.perf_counter() - t0
            return out

        # 5. Sweep all 2^n_stab sign sectors, find min-energy ground state.
        sector_results: List[Dict[str, Any]] = []
        best = None
        n_sectors = 2 ** n_stab

        # Compute representative sector (signs all +) to get Q/N_pauli/lambda
        # (sector-invariant: signs only flip individual coefficients)
        rep_signed = [+1 * s for s in stabilizers]
        rep_qop = taper_off_qubits(qop_rot, rep_signed)
        Q_rep = n_qubits_of(rep_qop)
        N_pauli_rep = n_pauli_nonid(rep_qop)
        lambda_rep = one_norm_nonid(rep_qop)

        if Q_rep <= gs_cutoff:
            # Sweep sectors only if small enough
            for sec in range(n_sectors):
                signs = [(+1 if (sec >> k) & 1 == 0 else -1)
                         for k in range(n_stab)]
                signed_stabs = [s * sgn for s, sgn in zip(stabilizers, signs)]
                qop_sec = taper_off_qubits(qop_rot, signed_stabs)
                Q_sec = n_qubits_of(qop_sec)
                E_sec = ground_state_energy(qop_sec, Q_sec, max_q=gs_cutoff)
                sector_results.append({
                    'signs': signs, 'Q': Q_sec, 'E': E_sec,
                })
                if E_sec is not None and (best is None or E_sec < best['E']):
                    best = {'E': E_sec, 'Q': Q_sec, 'signs': signs}

        # Tapered Q / N_pauli / lambda from the representative sector
        out['Q_tapered'] = Q_rep
        out['delta_Q'] = Q_naive - Q_rep
        out['N_pauli_tapered'] = N_pauli_rep
        out['lambda_tapered'] = lambda_rep

        # 6. Spectrum check vs naive
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
                out['relative_delta_E'] = (
                    out['delta_E'] / max(abs(out['E_naive']), 1e-30)
                )
                out['spectrum_preserved'] = bool(out['delta_E'] < 1e-8)
            else:
                out['delta_E'] = None
                out['relative_delta_E'] = None
                out['spectrum_preserved'] = None
        else:
            out['E_tapered'] = None
            out['delta_E'] = None
            out['relative_delta_E'] = None
            out['spectrum_preserved'] = None

        # Keep sector_table compact for JSON
        out['sector_table_len'] = len(sector_results)
        # store full table only if not too big
        if len(sector_results) <= 64:
            out['sector_table'] = sector_results

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
    ap.add_argument('--gs-cutoff', type=int, default=20)
    ap.add_argument('--max-n', type=int, default=2)
    ap.add_argument('--molecules', nargs='+', default=None)
    ap.add_argument('--out', default='debug/data/sprint_z2_per_subblock.json')
    ap.add_argument('--sector-cap-log2', type=int, default=14,
                    help='Skip molecules with > 2^N sign sectors '
                         '(N=14 -> 16384 sectors).')
    args = ap.parse_args()

    molecules = args.molecules if args.molecules else ALL_MOLECULES
    print(
        f'Running per-sub-block tapering on {len(molecules)} systems, '
        f'max_n={args.max_n}, gs-cutoff={args.gs_cutoff}.'
    )

    results: List[Dict[str, Any]] = []
    for i, name in enumerate(molecules):
        print(f'[{i+1:2d}/{len(molecules)}] {name} ...')
        res = taper_molecule_per_subblock(
            name, max_n=args.max_n, gs_cutoff=args.gs_cutoff, verbose=True,
            sector_cap_log2=args.sector_cap_log2,
        )
        if res['status'] == 'ok':
            de = res.get('delta_E')
            de_str = (f'|dE|={de:.2e}' if de is not None else 'dE=skip')
            print(
                f'  Q {res["Q_naive"]:3d}->{res["Q_tapered"]:3d} '
                f'(-{res["delta_Q"]}), n_sb={res["n_sub_blocks"]}, '
                f'n_sb_anti={res["n_sub_blocks_with_antisym"]}, '
                f'n_stab={res["n_stabilizers"]}, '
                f'all_P_comm={res["all_P_commute_with_H"]}, '
                f'{de_str}, wall={res["wall_s"]:.1f}s'
            )
        elif res['status'].startswith('skip'):
            print(f'  SKIP: {res.get("reason")}')
        else:
            print(f'  [ERROR] {res.get("error")}')
        results.append(res)

    out_path = PROJECT_ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f'\nWrote {out_path}')

    # Summary table
    print('\nSummary (markdown):')
    print(
        '| Mol | n_sb | Q_naive | Q_tap | dQ_per_sb | Pauli_per_sb | '
        'lambda_per_sb | all_P_comm | |dE|/|E| |'
    )
    print(
        '|:----|:----:|:-------:|:-----:|:---------:|:------------:|'
        ':-------------:|:----------:|:-------:|'
    )
    for r in results:
        if r['status'] == 'ok':
            de = r.get('relative_delta_E')
            de_s = 'n/a' if de is None else f'{de:.1e}'
            print(
                f"| {r['name']} | {r['n_sub_blocks']} | "
                f"{r['Q_naive']} | {r['Q_tapered']} | {r['delta_Q']} | "
                f"{r['N_pauli_tapered']} | {r['lambda_tapered']:.3f} | "
                f"{r['all_P_commute_with_H']} | {de_s} |"
            )
        else:
            status = r['status']
            n_sb = r.get('n_sub_blocks', '?')
            Q = r.get('Q_naive', '?')
            print(f"| {r['name']} | {n_sb} | {Q} | {status} | | | | | |")

    n_ok = sum(1 for r in results if r['status'] == 'ok')
    n_dq_ge_5 = sum(
        1 for r in results
        if r['status'] == 'ok' and r.get('delta_Q', 0) >= 5
    )
    n_dq_ge_4 = sum(
        1 for r in results
        if r['status'] == 'ok' and r.get('delta_Q', 0) >= 4
    )
    n_dq_ge_3 = sum(
        1 for r in results
        if r['status'] == 'ok' and r.get('delta_Q', 0) >= 3
    )
    n_checked = sum(
        1 for r in results
        if r['status'] == 'ok' and r.get('spectrum_preserved') is not None
    )
    n_preserved = sum(
        1 for r in results
        if r['status'] == 'ok' and r.get('spectrum_preserved') is True
    )
    print(f'\n{n_ok}/{len(results)} ran without error.')
    print(f'{n_dq_ge_3} systems achieved dQ >= 3.')
    print(f'{n_dq_ge_4} systems achieved dQ >= 4.')
    print(f'{n_dq_ge_5} systems achieved dQ >= 5.')
    print(f'{n_preserved}/{n_checked} systems passed spectrum-preservation gate.')

    return 0


if __name__ == '__main__':
    sys.exit(main())
