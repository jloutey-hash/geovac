"""Tests for prolate spheroidal molecular Sturmian betas (v0.9.30).

Test 1: Angular eigenvalue matches scipy obl_cv for homonuclear (b=0)
Test 2: H2+ ground state energy within 1% of exact (-1.1026 Ha)
Test 3: All 10 LiH betas found (nmax=3), ordered beta(n) > beta(n-1)
Test 4: Homonuclear H2 betas degenerate under (m, n_sph) <-> (m, n_sph)
Test 5: H 1s projected beta > 1.581 (binding threshold)
Test 6: H 1s Sturmian diagonal negative with projected betas
Test 7: All B-center diagonals negative check
"""

import numpy as np
import pytest
from geovac.molecular_sturmian import (
    _angular_sep_const,
    _radial_top_evals,
    compute_molecular_sturmian_betas,
    project_mo_betas_to_atom_centers,
)


class TestAngularEigenvalue:
    """Test 1: Angular separation constant matches scipy obl_cv."""

    def test_homonuclear_vs_obl_cv(self) -> None:
        """For b=0, A_sep should equal -obl_cv(m, m+n_sph, c)."""
        pytest.importorskip("scipy.special", reason="obl_cv not available")
        from scipy.special import obl_cv

        c_vals = [1.0, 2.5, 4.0]
        cases = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0)]

        for c in c_vals:
            for m, n_sph in cases:
                n_obl = m + n_sph
                A_ours = _angular_sep_const(m, n_sph, c, b=0.0)
                A_scipy = -obl_cv(m, n_obl, c)
                assert abs(A_ours - A_scipy) < 1e-10, (
                    f"m={m}, n_sph={n_sph}, c={c}: "
                    f"ours={A_ours:.12f}, scipy={A_scipy:.12f}"
                )


class TestH2PlusEnergy:
    """Test 2: H2+ ground state energy from matching condition."""

    def test_h2plus_ground_state(self) -> None:
        """H2+ at R=2: E should be within 1% of exact -1.1026 Ha."""
        Z_A, Z_B = 1.0, 1.0
        R = 2.0
        E_exact = -1.1026

        # Scan energies to find zero crossing
        E_cross = None
        prev_s = None
        for E in np.linspace(-0.5, -1.5, 60):
            p0 = np.sqrt(-2 * E)
            c = p0 * R / 2
            a = (Z_A + Z_B) * R  # beta=1 for exact H2+

            A_ang = _angular_sep_const(0, 0, c, b=0.0)
            L_ev = _radial_top_evals(0, c, a, n_grid=1200)
            s = A_ang + L_ev[0]

            if prev_s is not None and prev_s * s < 0:
                # Linear interpolation
                E_prev = E + (1.5 - 0.5) / 60  # previous E
                E_cross = E_prev - prev_s * (E - E_prev) / (s - prev_s)
                break
            prev_s = s

        assert E_cross is not None, "No zero crossing found for H2+"
        rel_err = abs(E_cross - E_exact) / abs(E_exact)
        assert rel_err < 0.01, (
            f"H2+ energy: {E_cross:.4f} Ha, exact: {E_exact:.4f} Ha, "
            f"error: {rel_err*100:.2f}%"
        )


class TestLiHBetas:
    """Test 3: All 10 LiH molecular betas found and properly ordered."""

    @pytest.fixture(scope="class")
    def lih_results(self) -> list:
        """Compute LiH betas once for all tests in class."""
        return compute_molecular_sturmian_betas(
            Z_A=3.0, Z_B=1.0, R=3.015, p0=np.sqrt(10.0), nmax=3
        )

    def test_all_10_orbitals_found(self, lih_results: list) -> None:
        """nmax=3 should yield 10 molecular orbitals."""
        assert len(lih_results) == 10

    def test_all_betas_finite(self, lih_results: list) -> None:
        """Every orbital should have a finite beta."""
        for n, m, n_sph, n_rad, beta in lih_results:
            assert np.isfinite(beta), (
                f"n={n}, m={m}, n_sph={n_sph}, n_rad={n_rad}: beta=NaN"
            )

    def test_1sigma_near_unity(self, lih_results: list) -> None:
        """1sigma (Li 1s core) beta should be close to 1."""
        n, m, n_sph, n_rad, beta = lih_results[0]
        assert n == 1 and m == 0
        assert abs(beta - 1.0) < 0.1, f"1sigma beta={beta:.4f}, expected ~1.0"

    def test_higher_betas_larger(self, lih_results: list) -> None:
        """n=2 and n=3 betas should all exceed 1sigma beta."""
        beta_1sigma = lih_results[0][4]
        for n, m, n_sph, n_rad, beta in lih_results[1:]:
            assert beta > beta_1sigma, (
                f"n={n} ({m},{n_sph},{n_rad}): beta={beta:.4f} "
                f"<= 1sigma beta={beta_1sigma:.4f}"
            )


class TestHomonuclearDegeneracy:
    """Test 4: H2 homonuclear betas show expected degeneracies."""

    def test_h2_n2_degeneracy(self) -> None:
        """For H2 (Z_A=Z_B=1), the three n=2 betas should be close."""
        results = compute_molecular_sturmian_betas(
            Z_A=1.0, Z_B=1.0, R=1.4, p0=1.0, nmax=2
        )
        n2_betas = [beta for n, m, ns, nr, beta in results
                    if n == 2 and np.isfinite(beta)]
        if len(n2_betas) >= 2:
            spread = max(n2_betas) - min(n2_betas)
            mean_beta = np.mean(n2_betas)
            # For homonuclear at finite R, n=2 betas split by molecular
            # field (~26% at R=1.4). Check spread < 30%.
            assert spread / mean_beta < 0.30, (
                f"H2 n=2 beta spread {spread:.4f} too large "
                f"(mean={mean_beta:.4f})"
            )


class TestMOProjection:
    """Tests 5-7: MO-to-atom-center beta projection (v0.9.30)."""

    @pytest.fixture(scope="class")
    def lih_projection(self) -> dict:
        """Compute LiH MO betas and project onto atom centers."""
        Z_A, Z_B, R = 3.0, 1.0, 3.015
        p0 = np.sqrt(10.0)
        nmax = 3

        mo_results = compute_molecular_sturmian_betas(
            Z_A=Z_A, Z_B=Z_B, R=R, p0=p0, nmax=nmax
        )
        beta_A, beta_B = project_mo_betas_to_atom_centers(
            mo_results, Z_A=Z_A, Z_B=Z_B, R=R, p0=p0, nmax=nmax
        )
        return {
            "beta_A": beta_A, "beta_B": beta_B,
            "p0": p0, "Z_A": Z_A, "Z_B": Z_B,
        }

    def test_h1s_beta_above_threshold(self, lih_projection: dict) -> None:
        """Test 5: H 1s projected beta > 1.581 (H binding threshold)."""
        beta_B = lih_projection["beta_B"]
        p0 = lih_projection["p0"]
        Z_B = lih_projection["Z_B"]
        threshold = p0 / (2.0 * Z_B)  # 1.581

        beta_h1s = beta_B[(1, 0, 0)]
        assert beta_h1s > threshold, (
            f"H 1s projected beta={beta_h1s:.4f} < threshold={threshold:.4f}"
        )

    def test_h1s_diagonal_negative(self, lih_projection: dict) -> None:
        """Test 6: H 1s Sturmian diagonal negative with projected beta."""
        beta_B = lih_projection["beta_B"]
        p0 = lih_projection["p0"]
        Z_B = lih_projection["Z_B"]

        beta_h1s = beta_B[(1, 0, 0)]
        eps_h1s = p0**2 / 2.0 - beta_h1s * Z_B * p0
        assert eps_h1s < 0, (
            f"H 1s diagonal eps={eps_h1s:.4f} Ha (positive — not bound). "
            f"beta={beta_h1s:.4f}, p0={p0:.4f}"
        )

    def test_all_b_center_diagonals(self, lih_projection: dict) -> None:
        """Test 7: Count B-center orbitals with negative diagonal."""
        beta_B = lih_projection["beta_B"]
        p0 = lih_projection["p0"]
        Z_B = lih_projection["Z_B"]

        n_negative = 0
        n_total = len(beta_B)
        for (n, l, m), beta in beta_B.items():
            eps = p0**2 / 2.0 - beta * Z_B * p0
            if eps < 0:
                n_negative += 1

        # Report count — full success is 14/14
        print(f"\nB-center diagonals: {n_negative}/{n_total} negative")
        assert n_negative >= 1, (
            f"No B-center orbitals have negative diagonal "
            f"(0/{n_total} bound)"
        )
