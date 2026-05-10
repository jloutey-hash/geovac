"""
Screened-Schrödinger valence basis correction (Track 3, W1c-residual)
=====================================================================

Replaces the hydrogenic Z_orb=1 valence orbital h1 diagonal entries on
frozen-core centers (Na, Mg, ..., Sr, Ba, ...) with the actual radial
Schrödinger eigenvalues of the FrozenCore Z_eff(r) potential.

Background
----------

The composed-qubit balanced builder uses ``Z_orb = Z_eff_valence`` (=1
for Na, =2 for Mg, etc., from atomic_classifier) as a generic
hydrogenic radial parameter on the heavy-atom center side of bond
blocks and lone-pair blocks. The h1 diagonal is then

    h1[i, i] = -Z_orb^2 / (2 * block_n^2)

This is **not** the actual valence eigenvalue. For Na 3s, the
hydrogenic Z=1 1s gives -0.5 Ha; the actual screened Na 3s eigenvalue
is -0.170 Ha (a 0.33 Ha overbinding bias).

The diagonal h1 corrections computed here are the dominant
single-electron energy correction; cross-center V_ne and within-block
ERIs are not corrected by this module (full-basis substitution would
require recomputing Slater radial integrals against the screened
wavefunction; deferred to a follow-on if this minimal correction shows
binding signature).

Convention
----------

Following the relativistic-composed convention
(geovac/composed_qubit_relativistic.py), the physical n maps as

    physical_n = block_n + n_val_offset
    n_val_offset = period - 1  (set in molecular_spec.hydride_spec)

so for Na ([Ne] core, period 3): block_n=1 -> physical_n=3 (Na 3s),
block_n=2 -> physical_n=4 (Na 4s/4p). For [Kr] core (period 5):
block_n=1 -> physical_n=5, block_n=2 -> physical_n=6.

The block_n=2 with l=1 entries map to "hydrogenic 4p"-class virtuals,
NOT to "lowest 3p above core". This matches existing framework
convention so backward-compat tests remain bit-exact when the
correction flag is OFF.

When the engineering closure flag is ON, this module replaces
the h1 diagonal entries on the heavy-atom-side sub-blocks of frozen-
core hydrides with the eigenvalue of ``_solve_screened_radial_log``
at ``physical_n = block_n + n_val_offset``.

For l=0 we use the log-grid solver (``_solve_screened_radial_log``)
with ``allow_l0=True`` per the Sprint Cs-HFS infrastructure; for l>=1
we use the standard FD solver (``_solve_screened_radial``).

Author: GeoVac Development Team (Track 3, post-PK W1c-residual sprint)
Date: 2026-05-09
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from geovac.molecular_spec import OrbitalBlock, MolecularSpec
from geovac.neon_core import (
    FrozenCore,
    _solve_screened_radial,
    _solve_screened_radial_log,
)


# Module-level cache for screened eigenvalues (avoid solving the same
# (Z, l, n_phys) problem repeatedly for every R-point of a PES sweep).
# Key: (Z_nuc_int, l, n_phys, n_grid). Value: eigenvalue (Ha).
_screened_eigenvalue_cache: Dict[Tuple[int, int, int, int], float] = {}


def _detect_frozen_core_Z(Z_nuc: float) -> bool:
    """True if Z_nuc has a frozen core (auto-detected)."""
    Z_int = int(round(Z_nuc))
    if Z_int < 11:
        return False
    try:
        FrozenCore(Z=Z_int)
        return True
    except ValueError:
        return False


def screened_valence_eigenvalue(
    Z_nuc: int,
    block_n: int,
    l: int,
    n_val_offset: int,
    n_grid: int = 4000,
    r_max: float = 80.0,
) -> float:
    """Eigenvalue of the FrozenCore Z_eff(r) Schrödinger equation.

    Parameters
    ----------
    Z_nuc : int
        Bare nuclear charge of the heavy atom.
    block_n : int
        Block-index principal quantum number (1, 2, 3, ...).
    l : int
        Orbital angular momentum.
    n_val_offset : int
        Period - 1 for the heavy atom (set by hydride_spec).
    n_grid : int
        Radial grid resolution.
    r_max : float
        Radial extent (bohr).

    Returns
    -------
    float
        Schrödinger eigenvalue in Ha.

    Raises
    ------
    ValueError
        If Z_nuc has no registered FrozenCore, or if (block_n, l) is
        invalid (l >= block_n + n_val_offset).
    """
    Z_int = int(round(Z_nuc))
    if not _detect_frozen_core_Z(Z_int):
        raise ValueError(
            f"Z_nuc={Z_nuc} has no FrozenCore; the screened-valence "
            f"correction only applies to frozen-core centers (Z>=11)."
        )

    n_phys = block_n + n_val_offset
    if n_phys < l + 1:
        raise ValueError(
            f"Invalid (block_n={block_n}, l={l}, n_val_offset={n_val_offset}): "
            f"physical_n={n_phys} < l+1={l+1}."
        )

    cache_key = (Z_int, l, n_phys, n_grid)
    if cache_key in _screened_eigenvalue_cache:
        return _screened_eigenvalue_cache[cache_key]

    if l == 0:
        # s-states: log-grid solver (handles 1/r origin singularity)
        energy, _u, _r, _R0 = _solve_screened_radial_log(
            Z=Z_int, l=0, n_target=n_phys,
            n_grid=n_grid, r_max=r_max,
        )
    else:
        # l >= 1: standard FD solver (centrifugal barrier regularizes origin)
        energy, _u, _r = _solve_screened_radial(
            Z=Z_int, l=l, n_target=n_phys,
            n_grid=n_grid, r_max=r_max,
        )

    energy = float(energy)
    _screened_eigenvalue_cache[cache_key] = energy
    return energy


def screened_valence_h1_diagonal(
    Z_orb: float,
    Z_nuc: float,
    block_n: int,
    l: int,
    n_val_offset: int,
    n_grid: int = 4000,
    r_max: float = 80.0,
) -> Tuple[float, float, float]:
    """Compute the h1 diagonal correction for a frozen-core valence orbital.

    Returns the screened (true) eigenvalue, the hydrogenic baseline, and
    the correction (screened - hydrogenic) so callers can apply Δh1 directly.

    Parameters
    ----------
    Z_orb : float
        The orbital's nuclear charge as used in
        ``-Z_orb^2 / (2 * block_n^2)`` (typically Z_eff_valence).
    Z_nuc : float
        Bare nuclear charge (used for FrozenCore lookup).
    block_n, l, n_val_offset : int
        See ``screened_valence_eigenvalue``.

    Returns
    -------
    e_screened : float
        Actual eigenvalue (Ha).
    e_hydrogenic : float
        Hydrogenic baseline -Z_orb^2 / (2 * block_n^2) (Ha).
    delta_h1 : float
        e_screened - e_hydrogenic (Ha). Add this to existing h1 diagonal
        to convert from hydrogenic to screened.

    Raises
    ------
    ValueError
        If Z_nuc has no registered FrozenCore.
    """
    e_screened = screened_valence_eigenvalue(
        int(round(Z_nuc)), block_n, l, n_val_offset,
        n_grid=n_grid, r_max=r_max,
    )
    e_hydrogenic = -float(Z_orb) ** 2 / (2.0 * block_n ** 2)
    return e_screened, e_hydrogenic, e_screened - e_hydrogenic


def compute_screened_h1_correction_block(
    block: OrbitalBlock,
    states: List[Tuple[int, int, int]],
    n_grid: int = 4000,
    r_max: float = 80.0,
) -> Optional[np.ndarray]:
    """Compute the diagonal h1 correction for one center-side sub-block.

    Returns a (n_states, n_states) diagonal correction matrix, or None
    if the block is not a frozen-core valence block.

    Frozen-core detection: ``block.Z_nuc_center >= 11`` and
    ``block.n_val_offset > 0`` (set by ``hydride_spec`` for second-row+
    hydrides). First-row blocks (Z=1-10) and explicit-core blocks have
    n_val_offset=0 and are skipped.

    The correction is applied to all (n, l, m) states with
    ``physical_n = block_n + n_val_offset >= l + 1`` (always true for
    standard hydride enumeration).

    Parameters
    ----------
    block : OrbitalBlock
        The center-side sub-block.
    states : list of (n, l, m)
        Orbital quantum numbers in this sub-block.
    n_grid, r_max : int, float
        Radial solver parameters.

    Returns
    -------
    delta_h1 : ndarray of shape (n_states, n_states), or None
        Diagonal correction matrix (only diagonal entries nonzero).
        None if the block is not a frozen-core valence block.
    """
    Z_nuc = float(block.Z_nuc_center)
    if not _detect_frozen_core_Z(Z_nuc) or block.n_val_offset == 0:
        return None

    Z_orb = float(block.Z_center)
    n_val_offset = int(block.n_val_offset)

    n_states = len(states)
    delta = np.zeros((n_states, n_states))

    for i, (n_block, l, m) in enumerate(states):
        try:
            _, _, dh = screened_valence_h1_diagonal(
                Z_orb, Z_nuc, n_block, l, n_val_offset,
                n_grid=n_grid, r_max=r_max,
            )
            delta[i, i] = dh
        except (ValueError, RuntimeError):
            # Unbound state at this (n_phys, l) under the screened
            # potential -- leave the hydrogenic baseline unchanged for
            # this orbital and warn-on-debug.
            delta[i, i] = 0.0

    return delta


def apply_screened_valence_correction(
    spec: MolecularSpec,
    h1: np.ndarray,
    sub_blocks: List[Dict],
    n_grid: int = 4000,
    r_max: float = 80.0,
    verbose: bool = False,
) -> Tuple[np.ndarray, Dict]:
    """Apply the screened-eigenvalue h1 correction to a balanced Hamiltonian.

    Parameters
    ----------
    spec : MolecularSpec
        Molecular specification (used to identify frozen-core blocks).
    h1 : ndarray of shape (M, M)
        Existing one-body Hamiltonian (will not be modified; a corrected
        copy is returned).
    sub_blocks : list of dict
        Sub-block geometry from ``balanced_coupled._get_block_geometry``.
    n_grid, r_max : int, float
        Radial solver parameters.
    verbose : bool
        Print per-orbital correction trace.

    Returns
    -------
    h1_corrected : ndarray
        New h1 with the diagonal correction applied.
    info : dict
        Diagnostic info: per-block corrections, total trace shift.
    """
    h1_corrected = h1.copy()
    info: Dict = {
        'block_corrections': [],
        'total_trace_shift': 0.0,
        'n_orbitals_corrected': 0,
    }

    for sb in sub_blocks:
        if sb['side'] != 'center':
            continue
        parent_block = spec.blocks[sb['parent_block']]
        states = sb['states']
        off = sb['offset']

        delta = compute_screened_h1_correction_block(
            parent_block, states,
            n_grid=n_grid, r_max=r_max,
        )
        if delta is None:
            continue

        n_st = len(states)
        for i in range(n_st):
            h1_corrected[off + i, off + i] += delta[i, i]

        trace_shift = float(np.trace(delta))
        n_corrected = int(np.sum(np.abs(np.diag(delta)) > 1e-15))
        info['total_trace_shift'] += trace_shift
        info['n_orbitals_corrected'] += n_corrected
        info['block_corrections'].append({
            'sub_block': sb['label'],
            'parent_label': parent_block.label,
            'Z_nuc': float(parent_block.Z_nuc_center),
            'Z_orb': float(parent_block.Z_center),
            'n_val_offset': int(parent_block.n_val_offset),
            'n_states': n_st,
            'trace_shift_ha': trace_shift,
            'per_orbital': [
                {
                    'state': (int(n), int(l), int(m)),
                    'delta_h1': float(delta[i, i]),
                }
                for i, (n, l, m) in enumerate(states)
            ],
        })

        if verbose:
            print(
                f"  [screened-valence] {sb['label']}: "
                f"Z_nuc={parent_block.Z_nuc_center} Z_orb={parent_block.Z_center} "
                f"n_val_offset={parent_block.n_val_offset} "
                f"trace shift = {trace_shift:+.4f} Ha"
            )

    return h1_corrected, info


def clear_eigenvalue_cache() -> None:
    """Clear the module-level eigenvalue cache (use between tests)."""
    _screened_eigenvalue_cache.clear()
