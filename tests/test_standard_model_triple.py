"""Tests for the full Connes Standard Model on the GeoVac spectral triple (Sprint G4a).

Tests verify:
1. Finite triple (A_F = C (+) H (+) M_3(C)) axioms at KO-dim 6
2. Combined almost-commutative triple axioms at KO-dim 1
3. Inner fluctuation decomposition into U(1) x SU(2) x SU(3) + Higgs
4. Falsifier: Higgs sector non-trivial when Yukawa imposed
5. Consistency with H1 electroweak slice (lepton sector should match)
"""

import numpy as np
import pytest

from geovac.standard_model_triple import (
    GELL_MANN,
    StandardModelACTriple,
    StandardModelFiniteTriple,
    gell_mann_matrices,
    sm_gauge_only,
    sm_one_generation,
    standard_model_triple,
)


# ======================================================================
# Gell-Mann matrices
# ======================================================================

class TestGellMann:
    def test_count(self):
        assert len(GELL_MANN) == 8

    def test_hermitian(self):
        for i, g in enumerate(GELL_MANN):
            assert np.allclose(g, g.conj().T), f"lambda_{i+1} not Hermitian"

    def test_traceless(self):
        for i, g in enumerate(GELL_MANN):
            assert abs(np.trace(g)) < 1e-14, f"lambda_{i+1} not traceless"

    def test_orthogonality(self):
        for i in range(8):
            for j in range(8):
                tr = np.real(np.trace(GELL_MANN[i] @ GELL_MANN[j]))
                expected = 2.0 if i == j else 0.0
                assert abs(tr - expected) < 1e-13, f"Tr(lam_{i+1} lam_{j+1}) = {tr}, expected {expected}"

    def test_shape(self):
        for g in GELL_MANN:
            assert g.shape == (3, 3)


# ======================================================================
# StandardModelFiniteTriple
# ======================================================================

class TestFiniteTriple:
    @pytest.fixture
    def ft(self):
        return StandardModelFiniteTriple(
            yukawa_nu=0.0, yukawa_e=0.1, yukawa_u=0.2, yukawa_d=0.05
        )

    @pytest.fixture
    def ft_zero(self):
        return StandardModelFiniteTriple()

    def test_dimensions(self, ft):
        assert ft.dim_matter == 16
        assert ft.dim_H_F == 32

    def test_D_F_hermitian(self, ft):
        D = ft.dirac_F()
        assert D.shape == (32, 32)
        assert np.allclose(D, D.conj().T)

    def test_D_F_eigenvalues(self, ft):
        D = ft.dirac_F()
        eigs = np.sort(np.abs(np.linalg.eigvalsh(D)))
        # Eigenvalues should include ±y_e, ±y_u, ±y_d (each doubled from matter/anti)
        # plus zeros (from y_nu = 0) and tripled quarks
        assert eigs[0] < 1e-14  # zero from y_nu = 0

    def test_D_F_zero_yukawa(self, ft_zero):
        D = ft_zero.dirac_F()
        assert np.allclose(D, 0)

    def test_J_F_squared(self, ft):
        J = ft.real_structure_F()
        assert np.allclose(J @ J, np.eye(32))

    def test_gamma_F_squared(self, ft):
        g = ft.chirality_F()
        assert np.allclose(g @ g, np.eye(32))

    def test_gamma_F_hermitian(self, ft):
        g = ft.chirality_F()
        assert np.allclose(g, g.conj().T)

    def test_J_F_gamma_F_anticommute(self, ft):
        J = ft.real_structure_F()
        g = ft.chirality_F()
        assert np.allclose(J @ g + g @ J, 0), "{J_F, gamma_F} should be 0 (KO-dim 6)"

    def test_gamma_F_D_F_anticommute(self, ft):
        g = ft.chirality_F()
        D = ft.dirac_F()
        assert np.allclose(g @ D + D @ g, 0), "{gamma_F, D_F} should be 0"

    def test_J_F_D_F_commute(self, ft):
        J = ft.real_structure_F()
        D = ft.dirac_F()
        JDJinv = J @ np.conj(D) @ J.T
        assert np.allclose(JDJinv, D), "J_F D_F J_F^{-1} should equal D_F (KO-dim 6)"

    def test_algebra_action_on_matter_only(self, ft):
        I3 = np.eye(3, dtype=np.complex128)
        pi = ft.algebra_action(1.0 + 0.5j, (1.0, 0.2, 0.3, 0.1), I3)
        assert pi.shape == (32, 32)
        assert np.allclose(pi[16:32, :], 0), "pi_F should act on matter only"

    def test_algebra_action_block_structure(self, ft):
        I3 = np.eye(3, dtype=np.complex128)
        pi = ft.algebra_action(1.0, (1.0, 0.0, 0.0, 0.0), I3)
        matter = pi[0:16, 0:16]
        # Lepton block (0:4, 0:4) should be block_diag(q, diag(lam, conj_lam))
        lepton = matter[0:4, 0:4]
        assert abs(lepton[0, 0] - 1.0) < 1e-14  # q = 1 -> I_2
        assert abs(lepton[2, 2] - 1.0) < 1e-14  # lam = 1

    def test_lepton_quark_decoupled(self, ft):
        m = 0.5 * np.eye(3, dtype=np.complex128)
        pi = ft.algebra_action(1.0, (1.0, 0.0, 0.0, 0.0), m)
        matter = pi[0:16, 0:16]
        assert np.allclose(matter[0:4, 4:16], 0), "lepton-quark off-diagonal should be 0"
        assert np.allclose(matter[4:16, 0:4], 0)

    def test_color_action(self, ft):
        m = np.diag([1.0, 2.0, 3.0]).astype(np.complex128)
        pi = ft.algebra_action(1.0, (1.0, 0.0, 0.0, 0.0), m)
        matter = pi[0:16, 0:16]
        quark = matter[4:16, 4:16]
        # u_L colors at (0,3), (1,4), (2,5) within quark block
        # With flavor=I, quark = I_4 ⊗ m, so quark[0,0] = m[0,0] = 1
        assert abs(quark[0, 0] - 1.0) < 1e-14
        assert abs(quark[1, 1] - 2.0) < 1e-14
        assert abs(quark[2, 2] - 3.0) < 1e-14


# ======================================================================
# StandardModelACTriple — axioms
# ======================================================================

class TestSMACTripleAxioms:
    @pytest.fixture
    def T1(self):
        return sm_one_generation(1)

    @pytest.fixture
    def T2(self):
        return sm_one_generation(2)

    def test_dimensions_nmax1(self, T1):
        assert T1.dim_GV == 4
        assert T1.dim_F == 32
        assert T1.dim_H == 128

    def test_dimensions_nmax2(self, T2):
        assert T2.dim_GV == 16
        assert T2.dim_F == 32
        assert T2.dim_H == 512

    def test_D_hermitian(self, T2):
        D = T2.dirac_combined()
        assert np.allclose(D, D.conj().T)

    def test_J_squared_nmax1(self, T1):
        axioms = T1.verify_axioms()
        assert axioms["J_squared_residual"] < 1e-12

    def test_J_squared_nmax2(self, T2):
        axioms = T2.verify_axioms()
        assert axioms["J_squared_residual"] < 1e-12

    def test_JD_relation_nmax1(self, T1):
        axioms = T1.verify_axioms()
        assert axioms["JD_relation_residual"] < 1e-12

    def test_JD_relation_nmax2(self, T2):
        axioms = T2.verify_axioms()
        assert axioms["JD_relation_residual"] < 1e-12

    def test_order_zero_nmax1(self, T1):
        axioms = T1.verify_axioms()
        assert axioms["order_zero_max_residual"] < 1e-10

    def test_order_zero_nmax2(self, T2):
        axioms = T2.verify_axioms()
        assert axioms["order_zero_max_residual"] < 1e-10

    def test_order_one_nmax1(self, T1):
        axioms = T1.verify_axioms()
        assert axioms["order_one_max_residual"] < 1e-10

    def test_order_one_nmax2(self, T2):
        axioms = T2.verify_axioms()
        assert axioms["order_one_max_residual"] < 1e-10


# ======================================================================
# Gauge group identification
# ======================================================================

class TestGaugeGroup:
    @pytest.fixture
    def T2(self):
        return sm_one_generation(2)

    def test_gauge_group_full_sm(self, T2):
        census = T2.gauge_group_census(n_samples=30)
        assert census["SU2_lepton"], "SU(2) should be present in lepton sector"
        assert census["SU2_quark"], "SU(2) should be present in quark sector"
        assert census["SU3"], "SU(3) should be present"
        assert census["U1"], "U(1) should be present"
        assert census["gauge_group"] == "U(1) x SU(2) x SU(3)"

    def test_no_lepton_quark_mixing(self, T2):
        census = T2.gauge_group_census(n_samples=30)
        assert not census["lepton_quark_mixing"], "leptons and quarks should not mix"

    def test_gauge_trivial_nmax1(self):
        T1 = sm_one_generation(1)
        census = T1.gauge_group_census(n_samples=20)
        assert census["gauge_group"] == "trivial", \
            "At n_max=1, [D_GV, M_GV]=0 so gauge should be trivial"


# ======================================================================
# Falsifier (Higgs sector)
# ======================================================================

class TestFalsifier:
    def test_higgs_nonzero_with_yukawa(self):
        T = sm_one_generation(2)
        higgs_zero, reason, data = T.check_natural_negative(n_random_generators=20)
        assert not higgs_zero, f"Higgs should be non-trivial with imposed Yukawa: {reason}"

    def test_higgs_zero_without_yukawa(self):
        T = sm_gauge_only(2)
        higgs_zero, reason, data = T.check_natural_negative(n_random_generators=20)
        assert higgs_zero, f"Higgs should be zero without Yukawa: {reason}"

    def test_gauge_nonzero_nmax2(self):
        T = sm_one_generation(2)
        _, _, data = T.check_natural_negative(n_random_generators=20)
        assert data["gauge_mean"] > 1e-6, "gauge should be non-trivial at n_max >= 2"


# ======================================================================
# Decomposition structure
# ======================================================================

class TestDecomposition:
    @pytest.fixture
    def omega_nmax2(self):
        T = sm_one_generation(2)
        rng = np.random.default_rng(42)
        I3 = np.eye(3, dtype=np.complex128)
        generators = [(
            T.gv_multiplier(1), 0.5 + 0.3j, (0.1, 0.2, 0.3, 0.4),
            I3 + 0.1 * rng.standard_normal((3, 3)),
            T.gv_multiplier(2), 0.2 - 0.1j, (0.4, 0.3, 0.2, 0.1),
            I3 + 0.1 * rng.standard_normal((3, 3)),
        )]
        omega = T.inner_fluctuation_one_form(generators)
        return T, omega

    def test_matter_antimatter_decoupled(self, omega_nmax2):
        T, omega = omega_nmax2
        dc = T.decompose_fluctuation(omega)
        assert np.linalg.norm(dc["mat_anti_off"]) < 1e-10
        assert np.linalg.norm(dc["anti_mat_off"]) < 1e-10

    def test_lepton_quark_decoupled(self, omega_nmax2):
        T, omega = omega_nmax2
        dc = T.decompose_fluctuation(omega)
        assert np.linalg.norm(dc["lepton_quark_off"]) < 1e-10
        assert np.linalg.norm(dc["quark_lepton_off"]) < 1e-10

    def test_higgs_in_LR_blocks(self, omega_nmax2):
        T, omega = omega_nmax2
        dc = T.decompose_fluctuation(omega)
        h_lepton = np.linalg.norm(dc["lepton_higgs_LR"]) + np.linalg.norm(dc["lepton_higgs_RL"])
        h_quark = np.linalg.norm(dc["quark_higgs_LR"]) + np.linalg.norm(dc["quark_higgs_RL"])
        assert h_lepton > 0 or h_quark > 0, "Higgs should be in off-diagonal L-R blocks"


# ======================================================================
# Color content extraction
# ======================================================================

class TestColorContent:
    def test_su3_coefficients_nonzero(self):
        T = sm_one_generation(2)
        rng = np.random.default_rng(77)
        I3 = np.eye(3, dtype=np.complex128)
        m = I3 + 0.3 * (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))
        generators = [(
            T.gv_multiplier(1), 0.5, (0.1, 0.2, 0.3, 0.4), m,
            T.gv_multiplier(3), 0.3, (0.4, 0.1, 0.2, 0.3), I3,
        )]
        omega = T.inner_fluctuation_one_form(generators)
        color = T.extract_color_content(omega)
        su3_norm = np.linalg.norm(color["su3_coeffs_L"])
        assert su3_norm > 1e-6, "SU(3) content should be non-zero with non-trivial m"


# ======================================================================
# Consistency with H1 electroweak slice
# ======================================================================

class TestH1Consistency:
    def test_lepton_yukawa_matches_h1(self):
        from geovac.almost_commutative import ElectroweakFiniteTriple
        ew = ElectroweakFiniteTriple(yukawa_nu=0.0, yukawa_e=0.1)
        sm = StandardModelFiniteTriple(yukawa_nu=0.0, yukawa_e=0.1, yukawa_u=0.0, yukawa_d=0.0)

        D_ew = ew.dirac_F()  # 8x8
        D_sm = sm.dirac_F()  # 32x32

        # SM lepton block (0:4, 0:4 in matter, 16:20, 16:20 in antimatter)
        D_sm_lepton_matter = D_sm[0:4, 0:4]
        D_ew_lepton_matter = D_ew[0:4, 0:4]
        assert np.allclose(D_sm_lepton_matter, D_ew_lepton_matter), \
            "SM lepton Yukawa should match H1 electroweak"

    def test_chirality_lepton_matches_h1(self):
        from geovac.almost_commutative import ElectroweakFiniteTriple
        ew = ElectroweakFiniteTriple()
        sm = StandardModelFiniteTriple()

        g_ew = ew.chirality_F()  # 8x8
        g_sm = sm.chirality_F()  # 32x32

        g_ew_matter = g_ew[0:4, 0:4]
        g_sm_lepton_matter = g_sm[0:4, 0:4]
        assert np.allclose(g_sm_lepton_matter, g_ew_matter), \
            "SM lepton chirality should match H1"


# ======================================================================
# Convenience constructors
# ======================================================================

class TestConstructors:
    def test_sm_gauge_only(self):
        T = sm_gauge_only(1)
        assert T.dim_H == 128
        D_F = T.D_F()
        assert np.allclose(D_F, 0)

    def test_sm_one_generation(self):
        T = sm_one_generation(1)
        D_F = T.D_F()
        assert not np.allclose(D_F, 0)

    def test_standard_model_triple(self):
        T = standard_model_triple(1, yukawa_e=0.1)
        assert T.dim_H == 128
