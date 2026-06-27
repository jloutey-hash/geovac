"""
Tests for channel-resolved 1-RDM extraction from Level 4 wavefunctions.

Validates:
1. Diagonal consistency with existing channel weights
2. Hermiticity (real symmetry) of gamma matrix
3. Trace = 1.0 (normalized eigenvector, not electron count)
4. Off-diagonal elements are non-negligible
"""

import numpy as np
import pytest

from geovac.inter_fiber_coupling import (
    extract_channel_data,
    extract_channel_1rdm,
    exchange_inter_fiber_energy,
    full_exchange_inter_fiber_energy,
    kinetic_orthogonalization_energy,
)
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK


@pytest.fixture(scope='module')
def beh2_setup():
    """Set up core + PK for BeH2 1-RDM tests."""
    core = CoreScreening(Z=4, l_max=2, n_alpha=200)
    core.solve(verbose=False)
    pk = AbInitioPK(core, n_core=2)
    pk_potentials = [pk.pk_dict(atom='A')]
    return {
        'pk_potentials': pk_potentials,
        'Z_eff': 2.0,
        'Z_ligand': 1.0,
        'l_max': 2,
        'n_alpha': 100,
    }


@pytest.fixture(scope='module')
def rdm_at_R3(beh2_setup):
    """Compute 1-RDM at R=3.0 (near equilibrium).

    Note: ``solve_level4_h2_multichannel`` no longer accepts a
    ``pk_potentials`` keyword (PK is now partitioned classically in the
    composed_qubit pipeline; see CLAUDE.md §3 "Polyatomic coupling"
    entries).  The underlying Level 4 wavefunction is computed without
    PK; PK is then applied at the channel-data extraction step where
    the inter_fiber_coupling API still supports it.  This exercises
    the 1-RDM and exchange machinery (the still-live infrastructure)
    on a Level 4 solution; the specific PK-augmented Level 4
    architecture for BeH2 binding has been superseded.
    """
    s = beh2_setup
    R = 3.0
    l4 = solve_level4_h2_multichannel(
        R=R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'], n_Re=300,
        verbose=False,
    )
    # Get channel data (shared between 1-RDM and existing functions).
    # The inter_fiber_coupling functions still accept pk_potentials.
    ch_data = extract_channel_data(
        l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'],
        pk_potentials=s['pk_potentials'], n_sample_Re=12,
    )
    rdm = extract_channel_1rdm(
        l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'],
        pk_potentials=s['pk_potentials'],
        channel_data=ch_data,
    )
    return {'rdm': rdm, 'ch_data': ch_data, 'R': R}


# ---------------------------------------------------------------------------
# 1. Diagonal consistency: gamma_{l1,l1} matches sum of channel weights
# ---------------------------------------------------------------------------

def test_1rdm_diagonal_consistency(rdm_at_R3):
    """Diagonal of gamma must match channel weights summed over l2.

    For each R_e sample, gamma_{l1,l1} = sum_{l2} ch_weight[(l1,l2)].
    """
    rdm = rdm_at_R3['rdm']
    ch_data = rdm_at_R3['ch_data']

    gamma_matrix = rdm['gamma_matrix']  # (n_l, n_l, n_sample)
    l_values = rdm['l_values']
    channels = rdm['channels']
    ch_weights = rdm['ch_weights']  # (n_sample, n_ch)
    n_sample = gamma_matrix.shape[2]

    l_to_idx = {l: i for i, l in enumerate(l_values)}

    for k in range(n_sample):
        for l1 in l_values:
            i = l_to_idx[l1]
            # Sum channel weights for all channels with this l1
            expected = 0.0
            for ic, ch in enumerate(channels):
                if ch[0] == l1:
                    expected += ch_weights[k, ic]
            actual = gamma_matrix[i, i, k]
            assert abs(actual - expected) < 1e-10, (
                f"R_e sample {k}, l1={l1}: gamma_diag={actual:.10f}, "
                f"sum_ch_weights={expected:.10f}, diff={abs(actual-expected):.2e}"
            )


# ---------------------------------------------------------------------------
# 2. Hermiticity: gamma is real symmetric
# ---------------------------------------------------------------------------

def test_1rdm_hermiticity(rdm_at_R3):
    """gamma_{l1,l1'} = gamma_{l1',l1} (real symmetric) at each R_e."""
    gamma_matrix = rdm_at_R3['rdm']['gamma_matrix']
    n_sample = gamma_matrix.shape[2]

    for k in range(n_sample):
        gamma_k = gamma_matrix[:, :, k]
        assert np.allclose(gamma_k, gamma_k.T, atol=1e-12), (
            f"Gamma not symmetric at sample {k}: "
            f"max diff = {np.max(np.abs(gamma_k - gamma_k.T)):.2e}"
        )

    # Also check the weighted average
    gamma_avg = rdm_at_R3['rdm']['gamma_avg']
    assert np.allclose(gamma_avg, gamma_avg.T, atol=1e-12), (
        f"Average gamma not symmetric: "
        f"max diff = {np.max(np.abs(gamma_avg - gamma_avg.T)):.2e}"
    )


# ---------------------------------------------------------------------------
# 3. Trace: Tr(gamma) = 1.0 (eigenvector normalization)
# ---------------------------------------------------------------------------

def test_1rdm_trace(rdm_at_R3):
    """Tr(gamma) should equal 1.0 at each R_e sample (eigenvector is
    normalized so sum of all channel weights = 1)."""
    gamma_matrix = rdm_at_R3['rdm']['gamma_matrix']
    n_sample = gamma_matrix.shape[2]

    for k in range(n_sample):
        tr = np.trace(gamma_matrix[:, :, k])
        assert abs(tr - 1.0) < 1e-8, (
            f"Tr(gamma) at sample {k} = {tr:.10f}, expected 1.0"
        )

    # Also check weighted average
    tr_avg = np.trace(rdm_at_R3['rdm']['gamma_avg'])
    assert abs(tr_avg - 1.0) < 1e-8, (
        f"Tr(gamma_avg) = {tr_avg:.10f}, expected 1.0"
    )


# ---------------------------------------------------------------------------
# 4. Off-diagonal elements are non-negligible
# ---------------------------------------------------------------------------

def test_1rdm_offdiag_nonzero(rdm_at_R3):
    """Off-diagonal elements should be > 1% of diagonal norm.

    If off-diagonal is negligible, the 1-RDM is already diagonal in the
    channel basis and there's nothing to gain from including off-diagonal
    exchange terms.
    """
    offdiag_ratio = rdm_at_R3['rdm']['offdiag_ratio']
    assert offdiag_ratio > 0.01, (
        f"Off-diagonal ratio = {offdiag_ratio:.4f} (<1%), "
        f"suggesting 1-RDM is nearly diagonal"
    )


# ---------------------------------------------------------------------------
# Full exchange fixtures and tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def full_exchange_at_R3(beh2_setup):
    """Compute full exchange energy at R=3.0."""
    s = beh2_setup
    R = 3.0
    l4 = solve_level4_h2_multichannel(
        R=R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'], n_Re=300,
        verbose=False, pk_potentials=s['pk_potentials'],
    )
    result = full_exchange_inter_fiber_energy(
        l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'],
        pk_potentials=s['pk_potentials'],
        n_sample_Re=10,
    )
    return {'result': result, 'l4': l4, 'R': R, 'setup': s}


# ---------------------------------------------------------------------------
# 4b. PK reaches the inter-fiber angular solve (regression guard)
# ---------------------------------------------------------------------------

def test_extract_channel_data_propagates_pk(beh2_setup):
    """``extract_channel_data`` must pass PK into the angular eigenproblem.

    Regression guard for the v2.7.0 shim (documented in
    ``docs/molecular_refactor_handoff.md`` §4): ``extract_channel_data``
    silently dropped ``pk_potentials`` when calling
    ``solve_angular_multichannel``.  After v3.56.0 restored PK to the radial
    wavefunction F (via ``solve_level4_h2_multichannel``) but left the shim,
    the inter-fiber angular eigenvectors became PK-blind while F was PK-aware.
    That mismatch fed PK-blind eigenvectors into the off-diagonal 1-RDM
    exchange (``full_exchange``), drifting BeH2 R_eq from 2.80 (11.7%) to
    3.00 (19.7%) — collapsing ``full_exchange`` onto the ``exchange`` minimum.

    If PK is propagated, the channel eigenvectors MUST change when PK is
    supplied vs. withheld.  A bit-identical result means the shim regressed.
    """
    s = beh2_setup
    R = 3.0
    l4 = solve_level4_h2_multichannel(
        R=R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'], n_Re=300,
        verbose=False, pk_potentials=s['pk_potentials'],
    )
    cd_pk = extract_channel_data(
        l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'],
        pk_potentials=s['pk_potentials'], n_sample_Re=8,
    )
    cd_nopk = extract_channel_data(
        l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'],
        pk_potentials=None, n_sample_Re=8,
    )
    # Compare the mid-range angular eigenvector (sign-agnostic — eigh phase).
    k = len(cd_pk['Re_samples']) // 2
    v_pk = np.abs(cd_pk['vec_2d_list'][k])
    v_nopk = np.abs(cd_nopk['vec_2d_list'][k])
    max_diff = float(np.max(np.abs(v_pk - v_nopk)))
    assert max_diff > 1e-4, (
        f"extract_channel_data appears to ignore pk_potentials "
        f"(max eigenvector diff = {max_diff:.2e}); the PK-drop shim has "
        f"regressed — inter-fiber angular eigenvectors are PK-blind while the "
        f"radial wavefunction F is PK-aware, which silently breaks the "
        f"off-diagonal 1-RDM full_exchange R_eq."
    )


# ---------------------------------------------------------------------------
# 5. Full exchange is negative (attractive)
# ---------------------------------------------------------------------------

def test_full_exchange_negative(full_exchange_at_R3):
    """K_AB < 0: exchange between identical fibers is attractive."""
    E_exch = full_exchange_at_R3['result']['E_exchange']
    assert E_exch < 0, (
        f"Full exchange energy = {E_exch:.6f} Ha, expected negative (attractive)"
    )


# ---------------------------------------------------------------------------
# 6. Full exchange exceeds diagonal-only
# ---------------------------------------------------------------------------

def test_full_exchange_offdiag_nonzero(full_exchange_at_R3):
    """Off-diagonal exchange contribution is non-negligible.

    The off-diagonal 1-RDM elements couple channels with different l1,
    contributing exchange terms weighted by (-1)^{l1+l1'} and cross-channel
    F^0. These may add to or subtract from the diagonal exchange depending
    on parity, but should be > 1% of the diagonal magnitude.

    NOTE: this is a sanity check, NOT the PK-propagation regression guard.
    The off-diagonal stays nonzero even in the PK-blind regression (the bug
    was wrong-valued PK-blind eigenvectors, not a zeroed off-diagonal). The
    actual regression guards are ``test_extract_channel_data_propagates_pk``
    (fast) and ``test_full_exchange_reduces_req_error`` (R_eq inward of
    exchange).
    """
    result = full_exchange_at_R3['result']
    E_offdiag = abs(result['E_exchange_offdiag'])
    E_diag = abs(result['E_exchange_diag'])
    ratio = E_offdiag / E_diag if E_diag > 1e-15 else 0.0
    assert ratio > 0.01, (
        f"|K_offdiag| = {E_offdiag:.6f}, |K_diag| = {E_diag:.6f}, "
        f"ratio = {ratio:.4f}: off-diagonal should be > 1% of diagonal"
    )


# ---------------------------------------------------------------------------
# 7. Full exchange R-dependence: less negative at large R
# ---------------------------------------------------------------------------

def test_full_exchange_r_dependence(beh2_setup):
    """K_full(R) becomes less negative at large R (exchange decays)."""
    s = beh2_setup
    energies = []
    for R in [2.5, 4.5]:
        l4 = solve_level4_h2_multichannel(
            R=R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
            l_max=s['l_max'], n_alpha=s['n_alpha'], n_Re=300,
            verbose=False, pk_potentials=s['pk_potentials'],
        )
        result = full_exchange_inter_fiber_energy(
            l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
            l_max=s['l_max'], n_alpha=s['n_alpha'],
            pk_potentials=s['pk_potentials'],
            n_sample_Re=10,
        )
        energies.append(result['E_exchange'])

    E_short, E_long = energies
    # Both should be negative, but E_long should be less negative
    assert E_short < E_long, (
        f"E_exch(R=2.5) = {E_short:.6f}, E_exch(R=4.5) = {E_long:.6f}: "
        f"exchange should decay at large R"
    )


# ---------------------------------------------------------------------------
# 8. Full exchange reduces R_eq error (improvement over S*F^0)
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_full_exchange_reduces_req_error(beh2_setup):
    """Full 1-RDM exchange pins BeH2 R_eq at ~11.7% AND beats S*F^0 exchange.

    Two assertions, both regression guards for the PK-propagation shim
    (see ``test_extract_channel_data_propagates_pk``):

    1. ``full_exchange`` R_eq error < 15% (Paper 17 Table: 2.80 bohr / 11.7%).
       The PK-blind regression drifted this to 3.00 bohr / 19.7%.
    2. ``full_exchange`` R_eq is STRICTLY inward of the diagonal ``exchange``
       (S*F^0) R_eq.  When the off-diagonal 1-RDM path is PK-blind it collapses
       onto the ``exchange`` minimum (both 3.00 bohr) and this strict
       inequality fails — the sharpest catch for the "full == exchange"
       regression class.
    """
    from geovac.composed_triatomic import ComposedTriatomicSolver

    R_grid = np.arange(2.2, 3.7, 0.2)  # brackets both minima (2.80 / 3.00)

    solver_full = ComposedTriatomicSolver.BeH2(
        l_max=2, interbond_mode='full_exchange', verbose=False)
    solver_full.solve_core()
    pes_full = solver_full.scan_pes(R_grid=R_grid, n_Re=300)
    R_eq_full = pes_full['R_eq']

    solver_exch = ComposedTriatomicSolver.BeH2(
        l_max=2, interbond_mode='exchange', verbose=False)
    solver_exch.solve_core()
    pes_exch = solver_exch.scan_pes(R_grid=R_grid, n_Re=300)
    R_eq_exch = pes_exch['R_eq']

    R_ref = 2.507  # bohr
    error_pct = abs(R_eq_full - R_ref) / R_ref * 100.0

    assert error_pct < 15.0, (
        f"full_exchange R_eq = {R_eq_full:.3f} bohr, error = {error_pct:.1f}% "
        f"(Paper 17: 2.80 bohr / 11.7%; PK-blind regression gives 3.00/19.7%)"
    )
    assert R_eq_full < R_eq_exch - 1e-9, (
        f"full_exchange R_eq ({R_eq_full:.3f}) must be strictly inward of "
        f"exchange R_eq ({R_eq_exch:.3f}); equality means the off-diagonal "
        f"1-RDM exchange has collapsed onto the diagonal S*F^0 mode "
        f"(PK-propagation regression)."
    )


# ---------------------------------------------------------------------------
# 9. Parity symmetry: gamma^B constructed correctly
# ---------------------------------------------------------------------------

def test_full_exchange_parity_symmetry(full_exchange_at_R3):
    """gamma^B_{l1',l1} = (-1)^{l1+l1'} * gamma^A_{l1',l1}.

    Verify the parity matrix is correct and that the exchange
    contraction uses the right symmetry.
    """
    result = full_exchange_at_R3['result']
    gamma_avg = result['gamma_avg']
    parity_matrix = result['parity_matrix']
    l_values = result['l_values']

    # Check parity matrix values
    n_l = len(l_values)
    for i, l1 in enumerate(l_values):
        for j, l1p in enumerate(l_values):
            expected = (-1) ** (l1 + l1p)
            assert parity_matrix[i, j] == expected, (
                f"Parity[{l1},{l1p}] = {parity_matrix[i,j]}, "
                f"expected {expected}"
            )

    # Construct gamma^B explicitly and verify
    gamma_B = parity_matrix * gamma_avg.T  # (-1)^{l1+l1'} * gamma^A_{l1',l1}

    # gamma^B should also be symmetric (since gamma^A is symmetric
    # and parity_matrix is symmetric)
    assert np.allclose(gamma_B, gamma_B.T, atol=1e-12), (
        f"gamma^B not symmetric: max diff = "
        f"{np.max(np.abs(gamma_B - gamma_B.T)):.2e}"
    )

    # Tr(gamma^A * gamma^B) should equal the exchange integrand sum
    # (before multiplying by F^0)
    trace_product = np.sum(gamma_avg * gamma_B)
    # This should be a well-defined finite number
    assert np.isfinite(trace_product), "Tr(gamma^A * gamma^B) is not finite"


# ---------------------------------------------------------------------------
# 10. Kinetic orthogonalization: positive (repulsive)
# ---------------------------------------------------------------------------

def test_kinetic_orthog_positive(beh2_setup):
    """DeltaT > 0: Löwdin orthogonalization raises kinetic energy."""
    s = beh2_setup
    R = 3.0
    l4 = solve_level4_h2_multichannel(
        R=R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'], n_Re=300,
        verbose=False, pk_potentials=s['pk_potentials'],
    )
    kin = kinetic_orthogonalization_energy(
        l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'],
        pk_potentials=s['pk_potentials'],
        n_sample_Re=10,
    )
    assert kin['delta_T'] > 0, (
        f"delta_T = {kin['delta_T']:.6f}, expected positive (repulsive)"
    )


# ---------------------------------------------------------------------------
# 11. Kinetic orthogonalization: R-dependence (varies with R)
# ---------------------------------------------------------------------------

def test_kinetic_orthog_r_dependence(beh2_setup):
    """DeltaT(R) is R-dependent (not constant).

    The kinetic correction depends on both the inter-fiber overlap
    (which decreases with R) and the l>=1 channel mixing (which can
    increase with R). The net R-dependence is non-trivial but must
    be present — a constant correction would indicate a bug.
    """
    s = beh2_setup
    delta_Ts = []
    for R in [2.5, 4.5]:
        l4 = solve_level4_h2_multichannel(
            R=R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
            l_max=s['l_max'], n_alpha=s['n_alpha'], n_Re=300,
            verbose=False, pk_potentials=s['pk_potentials'],
        )
        kin = kinetic_orthogonalization_energy(
            l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
            l_max=s['l_max'], n_alpha=s['n_alpha'],
            pk_potentials=s['pk_potentials'],
            n_sample_Re=10,
        )
        delta_Ts.append(kin['delta_T'])

    dT_short, dT_long = delta_Ts
    # Both positive
    assert dT_short > 0 and dT_long > 0, (
        f"DeltaT should be positive at both R: "
        f"DeltaT(2.5)={dT_short:.6f}, DeltaT(4.5)={dT_long:.6f}"
    )
    # R-dependent (differ by > 1%)
    rel_diff = abs(dT_short - dT_long) / max(dT_short, dT_long)
    assert rel_diff > 0.01, (
        f"DeltaT(R=2.5) = {dT_short:.6f}, DeltaT(R=4.5) = {dT_long:.6f}: "
        f"relative difference {rel_diff:.4f} < 1%, correction appears constant"
    )


# ---------------------------------------------------------------------------
# 12. Kinetic orthogonalization: magnitude is reasonable fraction of exchange
# ---------------------------------------------------------------------------

def test_kinetic_orthog_magnitude(full_exchange_at_R3, beh2_setup):
    """0.1% < |DeltaT / E_exchange| < 200%: physically reasonable."""
    s = beh2_setup
    R = full_exchange_at_R3['R']
    l4 = full_exchange_at_R3['l4']
    E_exch = abs(full_exchange_at_R3['result']['E_exchange'])

    kin = kinetic_orthogonalization_energy(
        l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'],
        pk_potentials=s['pk_potentials'],
        n_sample_Re=10,
    )
    ratio = kin['delta_T'] / E_exch if E_exch > 1e-15 else 0.0
    assert 0.001 < ratio < 2.0, (
        f"|DeltaT/E_exch| = {ratio:.4f}, expected between 0.1% and 200%"
    )


# ---------------------------------------------------------------------------
# 13. Full exchange + kinetic mode runs and gives bound state
# ---------------------------------------------------------------------------

def test_full_exchange_kinetic_runs(beh2_setup):
    """PES scan with 'full_exchange_kinetic' completes and gives bound state."""
    from geovac.composed_triatomic import ComposedTriatomicSolver

    solver = ComposedTriatomicSolver.BeH2(
        l_max=2, interbond_mode='full_exchange_kinetic', verbose=False)
    solver.solve_core()
    R_grid = np.arange(2.0, 4.5, 0.3)
    pes = solver.scan_pes(R_grid=R_grid, n_Re=300)

    assert pes['D_e'] > 0, (
        f"D_e = {pes['D_e']:.6f}: molecule not bound with "
        f"full_exchange_kinetic mode"
    )


# ---------------------------------------------------------------------------
# Run standalone diagnostic
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import json

    print("Setting up BeH2 core + PK...")
    core = CoreScreening(Z=4, l_max=2, n_alpha=200)
    core.solve(verbose=False)
    pk = AbInitioPK(core, n_core=2)
    pk_potentials = [pk.pk_dict(atom='A')]

    Z_eff, Z_B, l_max, n_alpha = 2.0, 1.0, 2, 100

    diag_data = {}
    for R in [2.0, 2.5, 3.0, 4.0, 5.0]:
        print(f"\n--- R = {R:.1f} bohr ---")
        l4 = solve_level4_h2_multichannel(
            R=R, Z_A=Z_eff, Z_B=Z_B,
            l_max=l_max, n_alpha=n_alpha, n_Re=300,
            verbose=False, pk_potentials=pk_potentials,
        )
        rdm = extract_channel_1rdm(
            l4, R, Z_A=Z_eff, Z_B=Z_B,
            l_max=l_max, n_alpha=n_alpha,
            pk_potentials=pk_potentials, n_sample_Re=15,
        )

        gamma_avg = rdm['gamma_avg']
        l_vals = rdm['l_values']
        tr = np.trace(gamma_avg)

        print(f"  l_values: {l_vals.tolist()}")
        print(f"  Tr(gamma) = {tr:.6f}")
        print(f"  ||offdiag|| / ||diag|| = {rdm['offdiag_ratio']:.4f}")
        print(f"  gamma_avg:")
        for i, l1 in enumerate(l_vals):
            row = "    "
            for j, l1p in enumerate(l_vals):
                row += f"{gamma_avg[i,j]:+.6f}  "
            print(row)

        diag_data[str(R)] = {
            'l_values': l_vals.tolist(),
            'gamma_avg': gamma_avg.tolist(),
            'trace': float(tr),
            'offdiag_ratio': rdm['offdiag_ratio'],
        }

    # Save diagnostic
    import os
    os.makedirs('debug/data', exist_ok=True)
    with open('debug/data/beh2_1rdm_diagnostic.json', 'w') as f:
        json.dump(diag_data, f, indent=2)
    print("\nDiagnostic saved to debug/data/beh2_1rdm_diagnostic.json")
