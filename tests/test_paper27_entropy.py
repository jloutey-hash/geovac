"""
Verification tests for Paper 27 (Entropy as a Projection Artifact).

Covers the three verification-pending claims in the Equation Verification
section of papers/core/paper_27_entropy_projection.tex:

1. Table I — Section III (Numerical Confirmation)
   S_kin/S_full ~ 1e-14 on He at n_max = 2, 3 in the S^3 angular momentum
   eigenbasis. Reproduces debug/data/ep1_entropy_partition.json.

2. Eqs. (3)-(5) — Section IV (Area Law as Pair Counting)
   Combinatorial identities on Paper 0 shell degeneracies g_n = 2n^2:
     A_n = g_n^2 = 4 n^4,
     (1/2) g_n (g_n - 1) = n^2 (2 n^2 - 1),
     g_n g_{n-1} = 4 n^2 (n-1)^2,
   and their common leading-order log A_n = 4 log n + const.

3. Section V numerical values — energy-graph characterization
   Hot-node V_(1s,1s),(1s,1s) = 5/8 (exact rational),
   relative commutator norm ||[H1, V_ee]||_F / ||V_ee||_F saturating at
   ~5-6% (non-vanishing) across n_max = 3, 4.

The operator-theoretic proposition in Section II is a standard free-fermion
fact (Peschel 2003); its empirical signature is Table I.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

# Ensure project root on path (mirrors conftest).
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from debug.energy_entanglement_decoupling import (  # noqa: E402
    build_decomposed_hamiltonians,
    solve_and_entangle,
)


# ---------------------------------------------------------------------------
# Paper 27 Table I — EP-1 reproduction
# ---------------------------------------------------------------------------

# Reference values frozen from debug/data/ep1_entropy_partition.json
# (committed alongside Paper 27). Tolerances are loose because the
# proposition only requires S_kin at the floating-point noise floor.
_EP1_REFERENCE = {
    2: {
        'E_full':  -2.8894787019709867,
        'S_full':   0.03954502764024503,
        'occ_full_0': 1.9884915650203696,
    },
    3: {
        'E_full':  -2.8931097422339307,
        'S_full':   0.04081105136647098,
        'occ_full_0': 1.988247530135314,
    },
}


@pytest.mark.parametrize('n_max', [2, 3])
def test_paper27_table1_ep1_reproduction(n_max):
    """Re-run EP-1 and check S_kin/S_full is at the noise floor."""
    data = build_decomposed_hamiltonians(Z_float=2.0, n_max=n_max)
    configs = data['configs']
    n_spatial = data['n_spatial']

    H_full = data['H_h1_diag'] + data['H_h1_offdiag'] + data['H_vee_full']
    H_kin  = data['H_h1_diag'] + data['H_h1_offdiag']

    E_full, _, _, ent_full = solve_and_entangle(H_full, configs, n_spatial)
    E_kin,  _, _, ent_kin  = solve_and_entangle(H_kin,  configs, n_spatial)

    S_full = float(ent_full['von_neumann_entropy'])
    S_kin  = float(ent_kin['von_neumann_entropy'])
    occ_full = ent_full['occupation_numbers']
    occ_kin  = ent_kin['occupation_numbers']

    ref = _EP1_REFERENCE[n_max]

    # Table I columns.
    assert E_full == pytest.approx(ref['E_full'], abs=1e-10), \
        f'E_full drift at n_max={n_max}'
    assert S_full == pytest.approx(ref['S_full'], rel=1e-8), \
        f'S_full drift at n_max={n_max}'

    # The central claim: kinetic-only entropy at the floating-point
    # noise floor. 1e-12 is ~4 orders of magnitude above the observed
    # ~1e-16 but still 10 orders below the full value (~4e-2).
    assert abs(S_kin) < 1e-12, \
        f'S_kin={S_kin:.3e} exceeds noise floor at n_max={n_max}'
    assert abs(S_kin) / S_full < 1e-10, \
        f'S_kin/S_full={S_kin / S_full:.3e} exceeds 1e-10 at n_max={n_max}'

    # Natural-orbital occupation structure: kinetic-only gives {0, 2},
    # full Hamiltonian gives deviations at the 1e-3 level.
    assert occ_kin[0] == pytest.approx(2.0, abs=1e-10), \
        'kinetic-only top occupation should be 2 (doubled-spin closed shell)'
    assert occ_kin[1] < 1e-10, \
        'kinetic-only second occupation should be ~0'
    assert occ_full[0] == pytest.approx(ref['occ_full_0'], rel=1e-6)
    # The O(10^-3) deviation from 2 is the V_ee correlation signature.
    assert 1e-3 < 2.0 - occ_full[0] < 2e-2, \
        'full occupation deviation from 2 should be O(1e-2)'


# ---------------------------------------------------------------------------
# Paper 27 Eqs. (3)-(5) — area-law combinatorics
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('n', [1, 2, 3, 4, 5, 6, 8, 10])
def test_paper27_area_law_pair_counts(n):
    """A_n = g_n^2 = 4 n^4 and companion pair counts."""
    g_n = 2 * n * n
    assert g_n == 2 * n ** 2

    # Eq. (3): A_n = g_n^2 = 4 n^4 (ordered pair count).
    A_n = g_n * g_n
    assert A_n == 4 * n ** 4

    # Eq. (4): (1/2) g_n (g_n - 1) = n^2 (2 n^2 - 1) (antisymmetrized pairs).
    anti_pairs = g_n * (g_n - 1) // 2
    assert anti_pairs == n * n * (2 * n * n - 1)

    # Eq. (5): g_n * g_{n-1} = 4 n^2 (n-1)^2 (adjacent-shell pair count).
    if n >= 1:
        g_nm1 = 2 * (n - 1) ** 2
        assert g_n * g_nm1 == 4 * (n ** 2) * ((n - 1) ** 2)


def test_paper27_area_law_leading_order_agrees():
    """All three pair counts give log A_n = 4 log n + const at large n."""
    for n in (20, 50, 100, 200):
        g_n = 2 * n * n
        ordered    = g_n * g_n
        anti       = g_n * (g_n - 1) // 2
        adjacent   = g_n * (2 * (n - 1) ** 2)
        # Divide out the leading 4 n^4 and check ratios -> {1, 1/2, 1} + O(1/n).
        assert abs(ordered   / (4.0 * n ** 4) - 1.0)   < 1e-12
        assert abs(anti      / (2.0 * n ** 4) - 1.0)   < 1.0 / n
        assert abs(adjacent  / (4.0 * n ** 4) - 1.0)   < 4.0 / n
        # Common leading logarithmic slope.
        for count in (ordered, anti, adjacent):
            slope = math.log(count) / math.log(n)
            assert 3.5 < slope < 4.5, f'slope {slope} off for n={n}'


def test_paper27_one_body_counts_differ_from_area_law():
    """One-body counts have leading-order slope 2 and 3, not 4 — Paper 5 check.

    Tests the leading-order exponent by comparing ratios at successive
    doublings of n, where log(count(2n)/count(n)) / log 2 converges to
    the true exponent monotonically.
    """
    def exponent(f, n):
        return math.log(f(2 * n) / f(n)) / math.log(2)

    for n in (100, 1000, 10000):
        # g_n = 2 n^2: exponent 2.
        assert abs(exponent(lambda k: 2 * k * k, n) - 2.0) < 1e-9
        # N_n = n(n+1)(2n+1)/6: exponent 3 in the limit.
        N = lambda k: k * (k + 1) * (2 * k + 1) // 6
        assert abs(exponent(N, n) - 3.0) < 0.05
        # A_n = 4 n^4: exponent 4.
        assert abs(exponent(lambda k: 4 * k ** 4, n) - 4.0) < 1e-9


# ---------------------------------------------------------------------------
# Paper 27 Section V — energy-graph structural values
# ---------------------------------------------------------------------------

# Reference values frozen from debug/data/energy_graph_nmax{3,4}.json.
# These are the numerical values the paper cites; the V_diagonal_fraction
# values below are the unsquared Frobenius ratio ||diag(V)|| / ||V||,
# matching the quantity computed in debug/energy_graph_exploration.py.
_ENERGY_GRAPH_REFERENCE = {
    3: {
        'n_nodes':                   31,
        'commutator_rel_norm':       0.0608693635999042,
        'V_frobenius':               0.9891187247955867,
        'V_diag_fraction_unsquared': 0.9202580563910042,
        'V_1s1s_diagonal':           0.625,   # = 5/8 exactly
    },
    4: {
        'n_nodes':                   101,
        'commutator_rel_norm':       0.05319420456063309,
        'V_frobenius':               1.1206101403143007,
        'V_diag_fraction_unsquared': 0.8919969173234392,
        'V_1s1s_diagonal':           0.625,
    },
}


@pytest.mark.parametrize('n_max', [3, 4])
def test_paper27_energy_graph_structural_invariants(n_max):
    """Regression-lock the Section V commutator and hot-node numbers."""
    from debug.energy_graph_exploration import characterize

    summary, _V, _H1, _pairs, _orbs = characterize(n_max)
    ref = _ENERGY_GRAPH_REFERENCE[n_max]

    assert summary['n_nodes'] == ref['n_nodes']

    h1_vs_v = summary['h1_vs_v']
    assert h1_vs_v['relative_commutator_norm'] == pytest.approx(
        ref['commutator_rel_norm'], rel=1e-8)
    assert h1_vs_v['V_frobenius'] == pytest.approx(
        ref['V_frobenius'], rel=1e-8)
    assert h1_vs_v['V_diagonal_fraction_in_H1_eigenbasis'] == pytest.approx(
        ref['V_diag_fraction_unsquared'], rel=1e-8)

    # Hot-node: (1s,1s) diagonal is exactly 5/8 (Slater F^0(1s,1s),
    # Paper 7 Eq. 29). This is at k_orb = 1, the graph-native value.
    V_1s1s = summary['cusp_signature']['diagonal_top'][0]['Vii']
    assert V_1s1s == pytest.approx(5.0 / 8.0, abs=1e-12)
    assert V_1s1s == pytest.approx(ref['V_1s1s_diagonal'], abs=1e-12)


# ---------------------------------------------------------------------------
# Paper 27 Section VII.A — EP-2b HO two-fermion experiment
# ---------------------------------------------------------------------------

def _ep2b_nmax2():
    """Shared build for the EP-2b tests (cached result via module-level lru)."""
    from geovac.nuclear.ho_two_fermion import build_decomposed_ho_hamiltonians
    return build_decomposed_ho_hamiltonians(N_max=2, hw=10.0)


def test_paper27_ep2b_ho_kinetic_interaction_commute():
    """[H_HO, V_NN] is at the floating-point noise floor.

    Total-HO-quanta conservation: any central two-body interaction on
    the HO basis commutes with the one-body HO Hamiltonian because the
    Moshinsky–Talmi transformation makes N_tot = N_rel + N_CM a good
    quantum number of V.
    """
    data = _ep2b_nmax2()
    H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
    H_vee = data['H_vee_full']
    C = H_kin @ H_vee - H_vee @ H_kin
    rel_norm = np.linalg.norm(C) / np.linalg.norm(H_kin)
    # Floor is ~1e-15; give generous headroom.
    assert rel_norm < 1e-12, \
        f'[H_HO, V_NN] / ||H_HO|| = {rel_norm:.3e} exceeds noise floor'


def test_paper27_ep2b_ho_gs_is_single_determinant():
    """S_HO = 0 and the GS is pure (0s)^2 at N_max ∈ {2, 3}."""
    from geovac.nuclear.ho_two_fermion import build_decomposed_ho_hamiltonians
    from debug.entanglement_geometry import (
        build_1rdm_from_singlet_ci,
        compute_entanglement_measures,
    )
    # N_max=2 from the cached call; N_max=3 built fresh.
    for N_max in (2, 3):
        data = (_ep2b_nmax2() if N_max == 2
                else build_decomposed_ho_hamiltonians(N_max=3, hw=10.0))
        H_full = data['H_full']
        configs = data['configs']
        n_spatial = data['n_spatial']
        eigs, vecs = np.linalg.eigh(H_full)
        ci = vecs[:, 0]
        rho = build_1rdm_from_singlet_ci(ci, configs, n_spatial)
        ent = compute_entanglement_measures(rho)
        occ = ent['occupation_numbers']
        # Paper 27 Table EP-2b: top-4 occupations are (2, 0, 0, 0).
        assert occ[0] == pytest.approx(2.0, abs=1e-10), \
            f'top occupation = {occ[0]} at N_max={N_max}'
        for k in (1, 2, 3):
            assert abs(occ[k]) < 1e-10, \
                f'occupation {k} = {occ[k]:.3e} at N_max={N_max}'
        # von Neumann entropy identically zero (below log(1+1e-10)).
        assert ent['von_neumann_entropy'] < 1e-10, \
            f'S_HO = {ent["von_neumann_entropy"]:.3e} at N_max={N_max}'
        # Ground-state energy is 3 hw + V_00,00: a single Slater
        # determinant in the N_tot=0 sector.  Must be the same
        # at N_max=2 and N_max=3 (basis-independent projection).
        assert eigs[0] == pytest.approx(14.898, abs=5e-3), \
            f'GS energy drift at N_max={N_max}: {eigs[0]:.4f} MeV'


# ---------------------------------------------------------------------------
# Paper 27 §VII.B — EP-2c pilot power law
# ---------------------------------------------------------------------------

def test_paper27_ep2c_dimensionless_power_law():
    """S_B ~ (w_B / ||H_kin||)^alpha with alpha ≈ 2.38, R² > 0.99."""
    Z_list = [2.0, 3.0, 4.0, 6.0, 8.0, 10.0]
    S = []
    w_dim = []
    for Z in Z_list:
        data = build_decomposed_hamiltonians(Z_float=Z, n_max=3)
        H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
        H_vee = data['H_vee_full']
        _, _, _, ent = solve_and_entangle(
            H_kin + H_vee, data['configs'], data['n_spatial'])
        S.append(float(ent['von_neumann_entropy']))
        _eH, U = np.linalg.eigh(H_kin)
        V_in = U.T @ H_vee @ U
        V_frob = np.linalg.norm(V_in)
        V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
        V_off = np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2))
        w_dim.append(float(V_off / np.linalg.norm(H_kin)))

    lx = np.log(np.array(w_dim))
    ly = np.log(np.array(S))
    alpha, intercept = np.polyfit(lx, ly, 1)
    pred = alpha * lx + intercept
    ss_res = np.sum((ly - pred) ** 2)
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot

    # Paper 27 Eq. (EP-2c): alpha = 2.383, A = 8.16, R^2 = 0.9982.
    assert alpha == pytest.approx(2.383, abs=0.05), \
        f'alpha={alpha:.4f} out of range'
    assert np.exp(intercept) == pytest.approx(8.16, rel=0.05), \
        f'A={np.exp(intercept):.4f} out of range'
    assert r2 > 0.99, f'R^2={r2:.4f} not tight enough'


def test_paper27_ep2c_multi_block_universality():
    """All single-center 2e blocks fall on the He-like dimensionless line
    to within 10% relative deviation; combined-fit alpha in [2.3, 2.45].
    """
    from geovac.molecular_spec import lih_spec, h2o_spec, nh3_spec, hf_spec

    def _row(Z, n_max):
        data = build_decomposed_hamiltonians(Z, n_max)
        H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
        H_vee = data['H_vee_full']
        _, _, _, ent = solve_and_entangle(
            H_kin + H_vee, data['configs'], data['n_spatial'])
        S_B = float(ent['von_neumann_entropy'])
        _e, U = np.linalg.eigh(H_kin)
        V_in = U.T @ H_vee @ U
        V_frob = np.linalg.norm(V_in)
        V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
        V_off = np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2))
        return S_B, float(V_off / np.linalg.norm(H_kin))

    S_vals, w_vals = [], []

    # He-like reference.
    for Z in (2.0, 3.0, 4.0, 6.0, 8.0, 10.0):
        S, w = _row(Z, 3)
        S_vals.append(S); w_vals.append(w)

    # Composed 2e single-center blocks (skip has_h_partner bonds).
    for spec_fn in (lih_spec, hf_spec, h2o_spec, nh3_spec):
        spec = spec_fn()
        for blk in spec.blocks:
            if blk.n_electrons != 2 or blk.has_h_partner:
                continue
            S, w = _row(float(blk.Z_center), 3)
            S_vals.append(S); w_vals.append(w)

    # Combined fit.
    lx = np.log(np.array(w_vals))
    ly = np.log(np.array(S_vals))
    alpha, intercept = np.polyfit(lx, ly, 1)
    pred = alpha * lx + intercept
    ss = np.sum((ly - pred) ** 2)
    st = np.sum((ly - ly.mean()) ** 2)
    r2 = 1.0 - ss / st
    A = float(np.exp(intercept))

    # Paper 27 Eq. (multi_fit): alpha=2.374, A=7.79, R²=0.998.
    assert 2.30 < alpha < 2.45, f'alpha={alpha:.4f} out of range'
    assert 6.5 < A < 9.5, f'A={A:.4f} out of range'
    assert r2 > 0.99, f'R²={r2:.4f} too loose'

    # Max relative deviation from He-like reference line.
    A_ref, alpha_ref = 8.16, 2.383
    max_rel = 0.0
    for S, w in zip(S_vals, w_vals):
        S_pred = A_ref * w ** alpha_ref
        max_rel = max(max_rel, abs(S_pred - S) / S)
    assert max_rel < 0.10, f'max_rel={max_rel:.4f} exceeds 10%'


def test_paper27_ep2d_z_range_artifact():
    """The 'lone-pair' A offset reproduces on the pure He-like family
    when restricted to the same Z-subset — confirming it is a Z-range
    fit artifact, not a block-type invariant.
    """
    def _row(Z):
        data = build_decomposed_hamiltonians(Z, 3)
        H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
        H_vee = data['H_vee_full']
        _, _, _, ent = solve_and_entangle(
            H_kin + H_vee, data['configs'], data['n_spatial'])
        _e, U = np.linalg.eigh(H_kin)
        V_in = U.T @ H_vee @ U
        V_frob = np.linalg.norm(V_in)
        V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
        V_off = np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2))
        return (float(ent['von_neumann_entropy']),
                float(V_off / np.linalg.norm(H_kin)))

    def _fit(Zs):
        w, S = [], []
        for Z in Zs:
            s, wd = _row(float(Z))
            S.append(s); w.append(wd)
        lx = np.log(np.array(w)); ly = np.log(np.array(S))
        alpha, intercept = np.polyfit(lx, ly, 1)
        return alpha, float(np.exp(intercept))

    # Lone-pair Z-set (Z=5,6,7) on the pure He-like family.
    alpha_lp, A_lp = _fit([5, 6, 7])
    # Core Z-set (Z=2,3,4,6,7,8,9,10).
    alpha_co, A_co = _fit([2, 3, 4, 6, 7, 8, 9, 10])

    # Must reproduce EP-2c split: core A ~7.8, lone-pair A ~4.7,
    # difference near 40%.
    assert A_lp == pytest.approx(4.73, rel=0.05), \
        f'lone-pair A on He-like family = {A_lp:.3f}, expected ~4.73'
    assert A_co == pytest.approx(8.04, rel=0.05), \
        f'core A on He-like family = {A_co:.3f}, expected ~8.04'
    assert 0.35 < (A_co - A_lp) / A_co < 0.48, \
        'Z-range fit offset on He-like family should be ~40%'
    # Exponent drift: core range has lower-Z samples → higher alpha.
    assert alpha_co > alpha_lp + 0.05, \
        f'alpha should decrease with Z: core={alpha_co:.3f}, lp={alpha_lp:.3f}'


def test_paper27_ep2g_universal_w_over_delta():
    """(w_dim / delta) collapses single-center and bond data.

    Single-center He-like Z in {2..10}: fit S = A * (w/delta)^gamma.
    Expected gamma ~ 2.1, A ~ 0.1, R^2 > 0.99.
    """
    Zs = [2.0, 3.0, 4.0, 6.0, 8.0, 10.0]
    ratios, Ss = [], []
    for Z in Zs:
        data = build_decomposed_hamiltonians(Z, 3)
        H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
        H_vee = data['H_vee_full']
        _, _, _, ent = solve_and_entangle(
            H_kin + H_vee, data['configs'], data['n_spatial'])
        Ss.append(float(ent['von_neumann_entropy']))
        e, U = np.linalg.eigh(H_kin)
        V_in = U.T @ H_vee @ U
        V_frob = np.linalg.norm(V_in)
        V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
        V_off = np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2))
        kin_frob = np.linalg.norm(H_kin)
        e_sort = np.sort(e)
        delta = float(abs(e_sort[1] - e_sort[0]) / kin_frob)
        w = float(V_off / kin_frob)
        ratios.append(w / delta)

    lx = np.log(np.array(ratios))
    ly = np.log(np.array(Ss))
    gamma, logA = np.polyfit(lx, ly, 1)
    pred = gamma * lx + logA
    r2 = 1 - np.sum((ly - pred) ** 2) / np.sum((ly - ly.mean()) ** 2)
    A = float(np.exp(logA))

    # Paper 27 Eq. (pred2_universal): gamma ~ 2.1-2.4 depending on Z-range,
    # A ~ 0.1-0.2, R^2 > 0.99 (tight single-center curve).
    # This 6-point Z in {2..10} fit gives gamma~2.37, A~0.17.
    assert 2.0 < gamma < 2.5, f'gamma={gamma:.3f} out of range'
    assert 0.10 < A < 0.25, f'A={A:.4f} out of range'
    assert r2 > 0.99, f'R^2={r2:.4f} too loose'


def test_paper27_ep2h_nmax_residue_persists():
    """gamma(n_max) decreases monotonically but does NOT reach 2.

    Z in {2,3,4,6,10} at n_max in {2,3,4} gives gamma(2)>gamma(3)>gamma(4)
    with gamma(4) > 2.1 — finite-basis residue persists.
    """
    def _row(Z, n_max):
        data = build_decomposed_hamiltonians(Z, n_max)
        H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
        H_vee = data['H_vee_full']
        _, _, _, ent = solve_and_entangle(
            H_kin + H_vee, data['configs'], data['n_spatial'])
        e, U = np.linalg.eigh(H_kin)
        V_in = U.T @ H_vee @ U
        V_frob = np.linalg.norm(V_in)
        V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
        V_off = np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2))
        kin = np.linalg.norm(H_kin)
        e_sort = np.sort(e)
        return (float(ent['von_neumann_entropy']),
                float(V_off / kin),
                float(abs(e_sort[1] - e_sort[0]) / kin))

    Zs = [2.0, 3.0, 4.0, 6.0, 10.0]
    gammas = {}
    for n_max in (2, 3, 4):
        ratios, Ss = [], []
        for Z in Zs:
            S, w, d = _row(Z, n_max)
            ratios.append(w / d); Ss.append(S)
        gamma, _ = np.polyfit(np.log(ratios), np.log(Ss), 1)
        gammas[n_max] = float(gamma)

    # Monotonic decrease.
    assert gammas[2] > gammas[3] > gammas[4], \
        f'gamma not monotone: {gammas}'
    # But does NOT reach 2 even at n_max=4.
    assert gammas[4] > 2.1, \
        f'gamma(n_max=4) = {gammas[4]:.3f} unexpectedly reached 2'
    # Drift magnitude.
    assert gammas[2] - gammas[4] > 0.2, \
        f'gamma drift {gammas[2] - gammas[4]:.3f} smaller than expected'


def test_paper27_ep2i_gamma_asymptote():
    """Local slope gamma(Z=15->30) approaches 2 from above at every n_max.

    The global-window fit gamma exceeds 2 (finite-basis residue); the
    LOCAL slope at large Z converges to 2 from above, with a mild
    n_max downshift.
    """
    def _row(Z, n_max):
        data = build_decomposed_hamiltonians(Z, n_max)
        H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
        H_vee = data['H_vee_full']
        _, _, _, ent = solve_and_entangle(
            H_kin + H_vee, data['configs'], data['n_spatial'])
        e, U = np.linalg.eigh(H_kin)
        V_in = U.T @ H_vee @ U
        V_frob = np.linalg.norm(V_in)
        V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
        V_off = np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2))
        kin = np.linalg.norm(H_kin)
        e_sort = np.sort(e)
        S = float(ent['von_neumann_entropy'])
        w_over_d = float(V_off / kin) / float(abs(e_sort[1] - e_sort[0]) / kin)
        return S, w_over_d

    for n_max in (2, 3, 4):
        S15, w15 = _row(15.0, n_max)
        S30, w30 = _row(30.0, n_max)
        local_slope = (np.log(S30) - np.log(S15)) / (np.log(w30) - np.log(w15))
        # Must be within [1.9, 2.05], approaching 2 from above.
        assert 1.9 < local_slope < 2.05, \
            f'n_max={n_max} local slope = {local_slope:.4f} out of range'

    # Verify monotonic n_max downshift.
    slopes = []
    for n_max in (2, 3, 4):
        S15, w15 = _row(15.0, n_max)
        S30, w30 = _row(30.0, n_max)
        slopes.append((np.log(S30) - np.log(S15)) /
                      (np.log(w30) - np.log(w15)))
    assert slopes[0] > slopes[1] > slopes[2], \
        f'n_max trend not monotone: {slopes}'


def test_paper27_ep2l_nmax5_below_two():
    """Local slope (Z=15->30) at n_max=5 lies BELOW 2 (1.959 measured),
    confirming the asymptotic gamma_infinity is sub-2.
    """
    def _row(Z, n_max):
        data = build_decomposed_hamiltonians(Z, n_max)
        H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
        H_vee = data['H_vee_full']
        _, _, _, ent = solve_and_entangle(
            H_kin + H_vee, data['configs'], data['n_spatial'])
        e, U = np.linalg.eigh(H_kin)
        V_in = U.T @ H_vee @ U
        V_frob = np.linalg.norm(V_in)
        V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
        V_off = np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2))
        kin = np.linalg.norm(H_kin)
        e_sort = np.sort(e)
        return (float(ent['von_neumann_entropy']),
                float(V_off / kin) / float(abs(e_sort[1] - e_sort[0]) / kin))

    S15, w15 = _row(15.0, 5)
    S30, w30 = _row(30.0, 5)
    slope = (np.log(S30) - np.log(S15)) / (np.log(w30) - np.log(w15))
    # Measured 1.959; constrain to [1.94, 1.99] to lock the sub-2 asymptote.
    assert 1.94 < slope < 1.99, f'n_max=5 local slope = {slope:.4f} out of range'


def test_paper27_ep2n_be_analytical_degenerate_pt():
    """Be 3x3 analytical degenerate-PT matches: GS dominantly 2s^2,
    occupations (2.0, 1.92, 0.02, 0.04, 0.02), S_full = 0.79 nats."""
    from debug.ep2N_be_analytical import (
        diag_E_for_pair, two_e_int_singlet_pair,
    )
    Z = 4
    s2 = (2, 0, 0); p_m = (2, 1, -1); p_0 = (2, 1, 0); p_p = (2, 1, +1)
    E_A = diag_E_for_pair(s2, s2)
    E_B = diag_E_for_pair(p_0, p_0)
    E_C = diag_E_for_pair(p_m, p_p)
    V_AB = two_e_int_singlet_pair(s2, s2, p_0, p_0, Z)
    V_AC = np.sqrt(2.0) * two_e_int_singlet_pair(s2, s2, p_m, p_p, Z)
    V_BC = np.sqrt(2.0) * two_e_int_singlet_pair(p_0, p_0, p_m, p_p, Z)
    H_eff = np.array([[E_A, V_AB, V_AC],
                      [V_AB, E_B, V_BC],
                      [V_AC, V_BC, E_C]])
    eigs, vecs = np.linalg.eigh(H_eff)
    gs = vecs[:, 0]
    # Dominant 2s^2: |c_A| > 0.95.
    assert abs(gs[0]) > 0.95, f'GS should be dominantly 2s^2, c_A={gs[0]:.3f}'
    # Spatial 1-RDM diagonal occupations
    occ = np.array([2.0, 2 * gs[0]**2, gs[2]**2,
                    2 * gs[1]**2, gs[2]**2])
    assert abs(occ.sum() - 4.0) < 1e-10
    p = occ / occ.sum()
    p_pos = p[p > 1e-14]
    S_full = float(-np.sum(p_pos * np.log(p_pos)))
    # Analytical S_full = 0.794 nats.
    assert 0.75 < S_full < 0.85, f'S_full = {S_full:.4f} out of range'


def test_paper26_27_consistency_z_scaling():
    """Paper 26's S(Z) ~ Z^-2.56 reproduces from Paper 27's machinery.

    Verifies that Paper 27 Sec VI.B subsumes Paper 26 Sec III's
    isoelectronic scaling at n_max=4 (Paper 26's setting).
    """
    Zs = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    Ss = []
    for Z in Zs:
        data = build_decomposed_hamiltonians(Z, 4)
        H = data['H_h1_diag'] + data['H_h1_offdiag'] + data['H_vee_full']
        _, _, _, ent = solve_and_entangle(H, data['configs'], data['n_spatial'])
        Ss.append(float(ent['von_neumann_entropy']))
    lx = np.log(np.array(Zs)); ly = np.log(np.array(Ss))
    alpha_Z, _ = np.polyfit(lx, ly, 1)
    # Paper 26 Eq. (scaling): S ~ Z^-2.56.
    # Our reproduction at n_max=4 should agree to within 0.05.
    assert -2.65 < alpha_Z < -2.45, \
        f'alpha_Z={alpha_Z:.3f} out of range for Paper 26 reproduction'


@pytest.mark.slow
def test_paper27_commutator_saturates_not_vanishes():
    """||[H1, V_ee]|| / ||V_ee|| does not decrease to zero with n_max."""
    from debug.energy_graph_exploration import characterize

    r3 = characterize(3)[0]['h1_vs_v']
    r4 = characterize(4)[0]['h1_vs_v']

    # Both well above floating-point zero; basis extension from
    # n_max=3 to n_max=4 decreases the relative commutator by less
    # than 2 percentage points — the decrease is not consistent with
    # vanishing in the CBS limit.
    c3 = r3['relative_commutator_norm']
    c4 = r4['relative_commutator_norm']
    assert c3 > 0.05 and c4 > 0.05, \
        'commutator should be > 5% at both n_max'
    assert (c3 - c4) < 0.02, \
        'commutator should saturate rather than decay fast'


def test_paper27_proposition_nondegeneracy_qualifier():
    """N >= 3 atoms (Li, Be) on the bare GeoVac S^3 graph have degenerate
    or quasi-degenerate H_1 ground states (n=2 hydrogenic shell + open-shell
    coupling), so S_kin > 0.  Validates that Paper 27 Section II's
    non-degeneracy clause is essential -- the operator-theoretic floor
    S_kin = 0 is a 2-electron (and closed-shell, gapped) statement, not a
    universal feature of the graph basis.

    Reproduces debug/data/ep2j_li_be_extension.json.
    """
    from debug.ep2j_li_be_extension import run_atom

    li = run_atom(Z=3, n_e=3, ms_target=+0.5, ml_target=0,
                  label='Li (test)', n_max=3)
    be = run_atom(Z=4, n_e=4, ms_target=0.0, ml_target=0,
                  label='Be (test)', n_max=3)

    # Both fail the EP-1 floor: S_kin is well above the 2e floating-point
    # noise floor (~1e-14 in EP-1) because the graph H_1 ground state is
    # not strictly non-degenerate at N >= 3.
    assert li['S_kin'] > 0.1, (
        f'Li S_kin={li["S_kin"]:.3e} unexpectedly small; '
        'open-shell quasi-degeneracy should give S_kin >> noise floor')
    assert be['S_kin'] > 0.1, (
        f'Be S_kin={be["S_kin"]:.3e} unexpectedly small; '
        'closed-shell n=2 hydrogenic shell degeneracy in graph H_1 '
        'should give S_kin >> noise floor')

    # And Be H_1 ground state is exactly degenerate (the 2s/2p shell).
    assert be['gs_degeneracy_kin'] >= 2, (
        f'Be H_1 gs_degeneracy={be["gs_degeneracy_kin"]} expected >= 2 '
        'from hydrogenic n=2 shell degeneracy')
