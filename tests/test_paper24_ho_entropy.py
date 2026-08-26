"""
Verification test for Paper 24 (Bargmann-Segal lattice), Section
sec:entanglement-rigidity -- which is a RETRACTION.

This file previously certified the two-fermion HO entanglement-rigidity
corollary.  That corollary is RETRACTED (2026-08-22) and this file now
pins the corrected physics plus a tripwire on the fix.

WHAT WAS CLAIMED (withdrawn):

    The spatial 1-RDM von Neumann entropy of the closed-shell
    two-fermion ground state on the Bargmann-Segal graph is
    identically zero for any central two-body interaction V(r_12),
    and the kinetic HO Hamiltonian commutes with the interaction:

        S_full = 0,   || [H_HO, V] ||_F / || H_HO ||_F < 1e-15.

WHY THE MECHANISM IS FALSE.  The claim rested on "any central V(r_12)
preserves the total HO quantum number N_tot".  The Moshinsky-Talmi
BRACKET conserves N -- that is a property of the HO coordinate
transformation.  The potential MATRIX ELEMENT does not.  A central
V(r_rel) is diagonal in the CM quantum numbers and in l_rel, but it
COUPLES different relative-n: the bra bracket forces
2n + l + 2N_CM + Lam = N_bra and the ket bracket forces
2n' + l + 2N_CM + Lam = N_ket, so n - n' = (N_bra - N_ket)/2.  The
N-changing elements are exactly the n != n' ones, and they are large --
<0,0|V|1,0> = +17.2 MeV against a -0.55 MeV diagonal for the Minnesota
singlet at b=1.

Paper 24's own text should have exposed this: it argued the Coulomb
case is nontrivial because 1/r_12 "does not preserve any HO-like
total-quanta quantum number" -- but 1/r_12 IS a central potential.
Both statements cannot hold.

WHERE THE ZERO CAME FROM.  An undocumented guard in
geovac/nuclear/moshinsky.py returned 0.0 whenever the bra and ket total
quantum numbers differed, discarding every coupling above.  Removed
2026-08-22; reproducible on demand via `conserve_N=True`, which
test_paper24_n_tot_guard_stays_removed uses as a tripwire.

METHODOLOGICAL NOTE.  The corroborating evidence was itself the
symptom.  The result read as robust partly because E0 came out
independent of N_max -- but that independence was the same guard,
freezing the ground state inside one N_tot block.  Corrected, E0
decreases with basis size, as a variational calculation must.  A
quantity that is EXACTLY zero, or EXACTLY invariant, where the physics
predicts only "small" or "weakly dependent", is a reason to ask what
has been restricted.

Corrected reference values (frozen, hbar*omega = 10 MeV, Minnesota
singlet):
    N_max = 2:  n_spatial=10, n_configs=15, E0=21.6538 MeV, S=0.0671
    N_max = 3:  n_spatial=20, n_configs=42, E0=21.6279 MeV, S=0.0716
    N_max = 4:                              E0=21.5442 MeV, S=0.0833
    (spatial trace-1 1-RDM, natural log, nats -- the same
    convention and the bit-identical builder used by Paper 27
    tab:ep2b; occupations (1.9800, 0.00666, ...) at N_max=2)
    rel commutator norm ~ 0.74 / 0.63 / 0.67 at N_max = 2 / 3 / 4

The pi-freeness results of Paper 24 are unaffected -- they concern the
graph's eigenvalues and adjacency weights, not the interacting
two-fermion state.  The asserted structural duality with Paper 27's
nonzero Coulomb scaling is withdrawn with the corollary.

The ground-state machinery is the production module
geovac.nuclear.ho_two_fermion; the 1-RDM and von Neumann entropy are
computed by the small self-contained helpers below (no debug imports).
"""

from __future__ import annotations

import os
import sys
from typing import List, Tuple

import numpy as np
import pytest

# Ensure project root on path (mirrors conftest).
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from geovac.nuclear.ho_two_fermion import (  # noqa: E402
    build_decomposed_ho_hamiltonians,
)


# ---------------------------------------------------------------------------
# Self-contained 1-RDM / entropy helpers (no archived-debug dependency).
# These reproduce debug/archive/misc/entanglement_geometry.py for the
# spatial-singlet two-fermion CI, so the Paper 24 corollary is backed
# entirely by in-suite production code.
# ---------------------------------------------------------------------------

def _build_singlet_1rdm(
    ci_coeffs: np.ndarray,
    configs: List[Tuple[int, int]],
    n_spatial: int,
) -> np.ndarray:
    """Spatial 1-RDM (Tr = 2) from spin-adapted singlet CI coefficients.

    Each config I = (i, j) is the spatial singlet
    [phi_i(1) phi_j(2) + phi_j(1) phi_i(2)] / N_IJ, with N_IJ = sqrt(2)
    for i != j and 1 for i == j.  Tracing out electron 2 gives

        rho_{a,c} = 2 * sum_{I,J} c_I c_J / (N_I N_J)
                      * sum_{(a,b) in I, (c,d) in J} delta(b, d).
    """
    rho = np.zeros((n_spatial, n_spatial))
    for I, (i, j) in enumerate(configs):
        for J, (p, q) in enumerate(configs):
            coeff = ci_coeffs[I] * ci_coeffs[J]
            if abs(coeff) < 1e-16:
                continue
            N_I = np.sqrt(2.0) if i != j else 1.0
            N_J = np.sqrt(2.0) if p != q else 1.0
            I_perms = [(i, j)] + ([(j, i)] if i != j else [])
            J_perms = [(p, q)] + ([(q, p)] if p != q else [])
            for a, b in I_perms:
                for c, d in J_perms:
                    if b == d:
                        rho[a, c] += coeff / (N_I * N_J)
    rho *= 2.0  # spin trace: closed-shell singlet, Tr = N_electrons = 2
    return rho


def _von_neumann_entropy(rho: np.ndarray) -> Tuple[float, np.ndarray]:
    """Return (S, occupations) for a 2-electron spatial 1-RDM (Tr = 2).

    Natural-orbital occupations n_i sum to 2; the normalized 1-RDM has
    eigenvalues lambda_i = n_i / 2 and S = -sum lambda_i log lambda_i.
    """
    occ = np.linalg.eigvalsh(rho)
    occ = np.sort(occ)[::-1]
    occ = np.maximum(occ, 0.0)
    lam = occ / 2.0
    S = 0.0
    for x in lam:
        if x > 1e-15:
            S -= x * np.log(x)
    return float(S), occ


# Frozen reference values (hbar*omega = 10 MeV, Minnesota NN singlet).
_REFERENCE = {
    # E0 CORRECTED 2026-08-22: the published 22.185 (identical at both
    # N_max) came from the N_tot truncation in moshinsky.py. See the
    # Paper 24 retraction.
    2: {'n_spatial': 10, 'n_configs': 15, 'E0_MeV': 21.6538},
    3: {'n_spatial': 20, 'n_configs': 42, 'E0_MeV': 21.6279},
}


def _build(N_max: int):
    return build_decomposed_ho_hamiltonians(N_max=N_max, hw=10.0)


@pytest.mark.parametrize('N_max', [2, 3])
def test_paper24_ho_entanglement_is_nonzero(N_max):
    """S_full is NONZERO: the GS is not a single Slater determinant.

    Was test_paper24_ho_zero_entanglement_entropy, asserting S = 0 as
    the entanglement-rigidity leg of the Paper 24 corollary. That
    corollary is RETRACTED (2026-08-22): the zero was produced by an
    undocumented N_tot guard in geovac/nuclear/moshinsky.py, not by
    physics. A central V(r_rel) is diagonal in the CM quantum numbers
    but couples different relative-n, and those couplings are large
    (+17.2 MeV off-diagonal against a -0.55 MeV diagonal).
    """
    data = _build(N_max)
    ref = _REFERENCE[N_max]

    # Build shape sanity.
    assert data['n_spatial'] == ref['n_spatial']
    assert data['n_configs'] == ref['n_configs']

    H_full = data['H_full']
    configs = data['configs']
    n_spatial = data['n_spatial']

    eigs, vecs = np.linalg.eigh(H_full)
    ci = vecs[:, 0]

    rho = _build_singlet_1rdm(ci, configs, n_spatial)
    # 1-RDM is a valid 2-electron density: trace = 2.
    assert np.trace(rho) == pytest.approx(2.0, abs=1e-10), \
        f'Tr(rho) = {np.trace(rho):.6f} at N_max={N_max}'

    S, occ = _von_neumann_entropy(rho)

    # The ground state is no longer a single Slater determinant, so the
    # leading occupation is close to but not exactly 2.
    assert occ[0] == pytest.approx(2.0, abs=5e-2), \
        f'top occupation = {occ[0]:.6f} at N_max={N_max}'
    for k in range(1, min(4, len(occ))):
        assert abs(occ[k]) < 5e-2, \
            f'occupation[{k}] = {occ[k]:.3e} nonzero at N_max={N_max}'

    # CORRECTED 2026-08-22. The published S = 0 was an artifact of an
    # undocumented N_tot truncation in moshinsky.py (Paper 24 retraction).
    # Measured with the correct two-body element: 0.0671 (N_max=2),
    # 0.0716 (3), 0.0833 (4) -- nonzero and slowly INCREASING with
    # basis. Band pinned tight enough to fail on the retracted zero
    # AND on a drift of more than ~40%.
    assert 0.04 < S < 0.11, \
        f'S_full = {S:.4f} outside the corrected band at N_max={N_max}'

    # Ground-state energy against the corrected reference.
    assert eigs[0] == pytest.approx(ref['E0_MeV'], abs=5e-2), \
        f'GS energy = {eigs[0]:.4f} MeV (expected {ref["E0_MeV"]}) ' \
        f'at N_max={N_max}'


@pytest.mark.parametrize('N_max', [2, 3])
def test_paper24_ho_kinetic_interaction_do_not_commute(N_max):
    """[H_HO, V] is O(1), not zero.

    Was test_paper24_ho_kinetic_interaction_commute. The claimed
    block-diagonality of V in N_tot is false: the Moshinsky-Talmi
    BRACKET conserves N, the potential MATRIX ELEMENT does not.
    RETRACTED 2026-08-22 with the rest of the corollary.
    """
    data = _build(N_max)
    H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
    H_vee = data['H_vee_full']

    # The one-body off-diagonal block is identically zero in the
    # Bargmann eigenbasis (pure HO is diagonal).
    assert np.linalg.norm(data['H_h1_offdiag']) < 1e-12, \
        'H_h1_offdiag should be identically zero for the pure HO'

    C = H_kin @ H_vee - H_vee @ H_kin
    rel_norm = np.linalg.norm(C) / np.linalg.norm(H_kin)
    # CORRECTED 2026-08-22: H_HO and V do NOT commute. Measured ratio
    # 0.74 / 0.63 / 0.67 at N_max = 2 / 3 / 4.
    assert 0.3 < rel_norm < 1.2, \
        f'|| [H_HO, V] ||_F / || H_HO ||_F = {rel_norm:.3e} ' \
        f'exceeds noise floor at N_max={N_max}'


def test_paper24_ho_gs_energy_is_variational():
    """E0 DECREASES with basis size, as a variational calculation must.

    Was test_paper24_ho_gs_energy_basis_independent, which asserted
    exact equality at N_max=2 and 3 and was cited as corroboration of
    the entanglement corollary. It was in fact the corollary's own
    symptom: the N_tot guard froze the ground state inside a single
    block, so enlarging the basis could not lower it. RETRACTED
    2026-08-22.
    """
    e2 = np.linalg.eigvalsh(_build(2)['H_full'])[0]
    e3 = np.linalg.eigvalsh(_build(3)['H_full'])[0]
    # CORRECTED 2026-08-22. E0 basis-INdependence was the truncation's
    # own signature: the ground state was frozen inside one N_tot block.
    # Corrected, E0 DECREASES with basis size, as a variational
    # calculation must (21.6538 -> 21.6279 -> 21.5442 MeV).
    assert e3 < e2, \
        f'GS energy basis-dependent: N_max=2 -> {e2:.6f}, ' \
        f'N_max=3 -> {e3:.6f} MeV'


def test_paper24_n_tot_guard_stays_removed():
    """Tripwire on the 2026-08-22 moshinsky fix.

    The retracted corollary was an artifact of `if N_bra != N_ket:
    return 0.0` in geovac/nuclear/moshinsky.py. This test reproduces the
    artifact on demand via the explicit `conserve_N` flag and confirms
    it is NOT the default -- so a silent reinstatement of the guard (or a
    default flip) fails here rather than quietly restoring a
    published-but-false zero.
    """
    import inspect
    from geovac.nuclear.moshinsky import lab_to_relative_matrix_element
    from geovac.nuclear.minnesota import minnesota_matrix_element_relative

    sig = inspect.signature(lab_to_relative_matrix_element)
    assert sig.parameters['conserve_N'].default is False, \
        'conserve_N must default to False; the N_tot guard is the artifact'

    # An N-changing lab-frame element: <0s,0s| V |0s,1s>, i.e. N_bra = 0
    # against N_ket = 2. The guard returned exactly 0.0 for every element
    # of this kind, discarding couplings an order of magnitude above the
    # diagonal terms it kept.
    kw = dict(n1=0, l1=0, n2=0, l2=0, n3=0, l3=0, n4=1, l4=0,
              L=0, S=0, V_rel_func=minnesota_matrix_element_relative, b=1.0)
    free = lab_to_relative_matrix_element(**kw, conserve_N=False)
    guarded = lab_to_relative_matrix_element(**kw, conserve_N=True)

    assert guarded == 0.0, 'conserve_N=True should reproduce the old guard'
    assert abs(free) > 1.0, \
        f'N-changing element should be O(1-10) MeV, got {free:.4f}'
