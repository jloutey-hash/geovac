r"""Fast regression backing for the LiH VMC-over-FCI build (debug/lih_vmc.py, Paper 19
"Fixed-geometry energy versus well shape").

The rigorous variational LiH energy -8.047 Ha is a STOCHASTIC VMC number (reproduce with
`python debug/lih_vmc.py gate6` / `final_run`; chronicled in CHANGELOG v5.15.22).  What is
DETERMINISTICALLY testable -- and what guarantees the VMC samples the correct Psi_CI * J --
is the machinery below:
  1. the analytic real-space orbital equals the engine's own grid evaluator (so the VMC
     orbital IS the function the ERIs were built from -> gate-6 must hold);
  2. the analytic orbital gradient and Laplacian match finite differences;
  3. the determinant sign convention sigma_I = (-1)^{inv(grouped spin-orbital order)}.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "debug"))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import lih_vmc as V  # noqa: E402


def _orbs():
    import mpmath as mp
    from prolate_mixed_eri import sto_orbital, valence_prolate_orbital, ZC_LI
    from prolate_allelectron_analytic_fci import sto_orbital_B
    from prolate_allelectron_c4 import from_sigma, valence_pi_orbital
    R = 3.015
    return R, {
        "Li1s_core": from_sigma(sto_orbital(ZC_LI, R, is_core=True)),
        "core2_4.5": from_sigma(sto_orbital(mp.mpf('4.5'), R, is_core=True)),
        "H1s": from_sigma(sto_orbital_B(mp.mpf('1.0'), R)),
        "bond21": from_sigma(valence_prolate_orbital(2, 1, mp.mpf('1.0'))),
        "pi_p00": valence_pi_orbital(0, 0, mp.mpf('1.0'), +1),
        "pi_m10": valence_pi_orbital(1, 0, mp.mpf('1.0'), -1),
    }


def test_orbital_gradient_matches_fd():
    R, orbs = _orbs()
    rng = np.random.default_rng(1)
    pts = rng.uniform(-2, 2, size=(6, 3)); pts[:, 2] += rng.uniform(-1, 1, 6)
    h = 1e-6
    for name, orb in orbs.items():
        prim = V.extract_prim(orb)
        _, grad = V.orbital_value_grad(prim, pts, R)
        gfd = np.empty_like(grad)
        for d in range(3):
            dp = np.zeros(3); dp[d] = h
            vp, _ = V.orbital_value_grad(prim, pts + dp, R)
            vm, _ = V.orbital_value_grad(prim, pts - dp, R)
            gfd[:, d] = (vp - vm) / (2 * h)
        err = np.max(np.abs(grad - gfd)) / max(np.max(np.abs(grad)), 1e-300)
        assert err < 1e-6, f"{name}: grad vs FD {err:.1e}"


def test_orbital_laplacian_matches_fd():
    R, orbs = _orbs()
    rng = np.random.default_rng(3)
    pts = rng.uniform(-2, 2, size=(6, 3)); pts[:, 2] += rng.uniform(-1, 1, 6)
    h = 1e-5
    for name, orb in orbs.items():
        prim = V.extract_prim(orb)
        v0, _, lap = V.orbital_vgl(prim, pts, R)
        lfd = np.zeros_like(lap)
        for d in range(3):
            dp = np.zeros(3); dp[d] = h
            vp, _, _ = V.orbital_vgl(prim, pts + dp, R)
            vm, _, _ = V.orbital_vgl(prim, pts - dp, R)
            lfd += (vp + vm - 2 * v0) / h ** 2
        err = np.max(np.abs(lap - lfd)) / max(np.max(np.abs(lap)), 1e-300)
        assert err < 1e-4, f"{name}: lap vs FD {err:.1e}"


def test_orbital_value_equals_engine_grid():
    """The real-space orbital must equal the engine's own _orb_on_grid (phi_azi=0)."""
    from prolate_allelectron_c4 import _orb_on_grid, _grid_template
    R, orbs = _orbs()
    tmpl = _grid_template(R, N_grid=24, xi_max=10.0)
    xi, eta = tmpl['xi'], tmpl['eta']
    ii, jj = len(xi) // 2, len(eta) // 2
    xi0, eta0 = xi[ii], eta[jj]
    rho = (R / 2) * np.sqrt(max((xi0 ** 2 - 1) * (1 - eta0 ** 2), 0.0))
    zc = (R / 2) * xi0 * eta0
    rc = np.array([[rho, 0.0, zc]])
    for name, orb in orbs.items():
        prim = V.extract_prim(orb)
        val, _ = V.orbital_value_grad(prim, rc, R)
        grid = _orb_on_grid(orb, tmpl)['psi'][ii, jj]   # (xi,eta) part, phi stripped
        err = abs(val[0] - grid) / max(abs(grid), 1e-300)
        assert err < 1e-9, f"{name}: real-space vs grid {err:.1e}"


def test_determinant_sign_convention():
    """sigma_I = (-1)^{inversions of the grouped spin-orbital order}; check known cases."""
    # alpha-orbs {0}, beta-orbs {0} -> [0,1] sorted -> +1
    assert V._perm_sign_inv([0, 1]) == 1.0
    # alpha {1}, beta {0} -> grouped [2,1] -> 1 inversion -> -1
    assert V._perm_sign_inv([2, 1]) == -1.0
    # alpha {0,1}, beta {0,1} -> grouped [0,2,1,3] -> inversions: (2,1) -> 1 -> -1
    assert V._perm_sign_inv([0, 2, 1, 3]) == -1.0
    # already sorted -> +1
    assert V._perm_sign_inv([0, 1, 2, 3]) == 1.0


def test_jastrow_cusp_finite_local_energy():
    """The Pade Jastrow keeps E_L finite as two opposite-spin electrons coalesce
    (the -1/2 lap(lnJ) ~ -A/r cancels V's +1/r); a divergence here = a cusp bug.
    Uses a synthetic 1-determinant wavefunction (no ERI build needed)."""
    R = 3.015
    # minimal wf: 2 orbitals, closed-shell-like, na=nb=2 needs 2 orbs; build by hand
    import mpmath as mp
    from prolate_mixed_eri import sto_orbital, valence_prolate_orbital, ZC_LI
    from prolate_allelectron_analytic_fci import sto_orbital_B
    from prolate_allelectron_c4 import from_sigma
    orbs = [from_sigma(sto_orbital(ZC_LI, R, is_core=True)),
            from_sigma(sto_orbital_B(mp.mpf('1.0'), R))]
    prims = [V.extract_prim(o) for o in orbs]
    T = np.eye(2)
    apairs = [(0, 1)]; bpairs = [(0, 1)]
    CS = np.array([[1.0]])
    wf = V.LiHWavefunction(prims, T, np.array([1.0]), apairs, bpairs, CS,
                           R, 3.0, 1.0, 3.0 / R, -8.0, 2, 2, 2)
    jas = V.PadeJastrow(1.5, 2, 2)
    rng = np.random.default_rng(0)
    r = V._init_walkers(wf, 1, rng)
    for dr in (0.5, 0.1, 0.02, 0.005):
        r[:, 2] = r[:, 0] + np.array([dr, 0.0, 0.0])   # e2(down) near e0(up)
        EL, _ = V.local_energy_analytic(wf, r, jas)
        assert np.all(np.isfinite(EL)), f"E_L not finite at r_02={dr}"
        assert abs(EL[0]) < 1e4, f"E_L blew up at r_02={dr}: {EL[0]}"


if __name__ == "__main__":
    test_orbital_gradient_matches_fd()
    test_orbital_laplacian_matches_fd()
    test_orbital_value_equals_engine_grid()
    test_determinant_sign_convention()
    test_jastrow_cusp_finite_local_energy()
    print("all lih_vmc deterministic gates PASS")
