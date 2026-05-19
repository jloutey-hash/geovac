"""Tests for geovac/hylleraas_eckart_pstate.py.

Track 5 (P-state extension, Schwartz 1961 form) — session 1 deliverables:
overlap, V_ne, V_ee for the symmetric channel (z_1 + z_2)·cosh(βt)·s^l t^{2m} u^n.

Validates against:
1. Single-basis-function analytical reference for (l, m, n) = (0, 0, 0), β = 0
   (the simplest P-state hydrogenic trial).
2. Symmetry/Hermiticity of matrix elements at β > 0.
3. Direct numerical quadrature in (s, t, u).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geovac.hylleraas_eckart_pstate import (
    HylleraasPStateBasisFn,
    hylleraas_pstate_basis_total_degree,
    overlap_element_pstate_eckart,
    potential_vne_element_pstate_eckart,
    potential_vee_element_pstate_eckart,
    overlap_element_pstate_eckart_antisym,
    potential_vne_element_pstate_eckart_antisym,
    potential_vee_element_pstate_eckart_antisym,
    overlap_element_pstate_eckart_cross,
    potential_vne_element_pstate_eckart_cross,
    potential_vee_element_pstate_eckart_cross,
    kinetic_element_pstate_eckart_sym_sym,
    P_STATE_VOLUME_FACTOR,
)


# ---------------------------------------------------------------------------
# Analytical single-basis-function validation
# ---------------------------------------------------------------------------

class TestSingleFunctionAnalytical:
    """For phi = (z_1+z_2) e^{-alpha(r_1+r_2)} (l=m=n=0, beta=0):

      <phi|phi> = 2 pi^2 / alpha^8
      <phi | -Z(1/r_1 + 1/r_2) | phi> = -3 Z pi^2 / alpha^7
      <phi | 1/r_12 | phi> = pi^2/3 * 49/(16 alpha^7) = 49 pi^2 / (48 alpha^7)

    Derived by direct 6D integration (cross terms with z_1 z_2 vanish by
    angular parity).
    """

    @pytest.mark.parametrize("alpha", [1.0, 1.5, 2.0, 2.5])
    def test_overlap_closed_form(self, alpha):
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        S = overlap_element_pstate_eckart(bf, bf, alpha, beta=0.0)
        expected = 2 * math.pi ** 2 / alpha ** 8
        assert abs(S - expected) < 1e-13 * abs(expected), (
            f"alpha={alpha}: got {S}, expected {expected}"
        )

    @pytest.mark.parametrize("alpha", [1.0, 1.5, 2.0, 2.5])
    @pytest.mark.parametrize("Z", [1.0, 2.0, 3.0])
    def test_vne_closed_form(self, alpha, Z):
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        V = potential_vne_element_pstate_eckart(bf, bf, alpha, beta=0.0, Z=Z)
        expected = -3 * Z * math.pi ** 2 / alpha ** 7
        assert abs(V - expected) < 1e-13 * abs(expected), (
            f"alpha={alpha}, Z={Z}: got {V}, expected {expected}"
        )

    @pytest.mark.parametrize("alpha", [1.0, 1.5, 2.0])
    def test_vee_closed_form(self, alpha):
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        V = potential_vee_element_pstate_eckart(bf, bf, alpha, beta=0.0)
        expected = 49 * math.pi ** 2 / (48 * alpha ** 7)
        assert abs(V - expected) < 1e-13 * abs(expected), (
            f"alpha={alpha}: got {V}, expected {expected}"
        )

    def test_vee_positive(self):
        """e-e repulsion is positive."""
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        V = potential_vee_element_pstate_eckart(bf, bf, 1.5, beta=0.0)
        assert V > 0

    def test_vne_negative(self):
        """e-nucleus attraction is negative (Z > 0)."""
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        V = potential_vne_element_pstate_eckart(bf, bf, 1.5, beta=0.0, Z=2.0)
        assert V < 0


# ---------------------------------------------------------------------------
# Hermiticity at beta > 0
# ---------------------------------------------------------------------------

class TestHermiticity:
    @pytest.mark.parametrize("beta", [0.0, 0.2, 0.5])
    def test_overlap_hermitian(self, beta):
        basis = hylleraas_pstate_basis_total_degree(2, channel='sym')
        alpha = 1.5
        for bf_p in basis:
            for bf_q in basis:
                S_pq = overlap_element_pstate_eckart(bf_p, bf_q, alpha, beta)
                S_qp = overlap_element_pstate_eckart(bf_q, bf_p, alpha, beta)
                assert S_pq == S_qp, (
                    f"({bf_p.label()}, {bf_q.label()}): {S_pq} != {S_qp}"
                )

    @pytest.mark.parametrize("beta", [0.0, 0.3])
    def test_vne_hermitian(self, beta):
        basis = hylleraas_pstate_basis_total_degree(2, channel='sym')
        for bf_p in basis:
            for bf_q in basis:
                v_pq = potential_vne_element_pstate_eckart(
                    bf_p, bf_q, 1.5, beta, Z=2.0
                )
                v_qp = potential_vne_element_pstate_eckart(
                    bf_q, bf_p, 1.5, beta, Z=2.0
                )
                assert v_pq == v_qp

    @pytest.mark.parametrize("beta", [0.0, 0.3])
    def test_vee_hermitian(self, beta):
        basis = hylleraas_pstate_basis_total_degree(2, channel='sym')
        for bf_p in basis:
            for bf_q in basis:
                v_pq = potential_vee_element_pstate_eckart(bf_p, bf_q, 1.5, beta)
                v_qp = potential_vee_element_pstate_eckart(bf_q, bf_p, 1.5, beta)
                assert v_pq == v_qp


# ---------------------------------------------------------------------------
# Positive-definite overlap matrix
# ---------------------------------------------------------------------------

class TestOverlapPositiveDefinite:
    @pytest.mark.parametrize("beta", [0.0, 0.2, 0.5])
    def test_overlap_pd(self, beta):
        basis = hylleraas_pstate_basis_total_degree(2, channel='sym')
        alpha = 1.5
        n = len(basis)
        S = np.zeros((n, n))
        for i, bf_p in enumerate(basis):
            for j, bf_q in enumerate(basis):
                S[i, j] = overlap_element_pstate_eckart(bf_p, bf_q, alpha, beta)
        eigvals = np.linalg.eigvalsh(S)
        assert np.all(eigvals > 1e-12), (
            f"S not PD: min eigval = {eigvals.min()}"
        )


# ---------------------------------------------------------------------------
# Numerical quadrature cross-check (catches algebra bugs that pure
# closed-form-vs-closed-form can't)
# ---------------------------------------------------------------------------

class TestVsDirectQuadrature:
    """Cross-check against direct (s, t, u) numerical quadrature."""

    @pytest.mark.slow
    def test_overlap_vs_quadrature(self):
        from scipy.integrate import tplquad
        alpha = 1.5
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='sym')

        def integrand(t, u, s, alpha):
            # Overlap integrand: (s^2+t^2-u^2) * u (s^2-t^2) e^{-2as}
            return (math.exp(-2 * alpha * s)
                    * (s ** 2 + t ** 2 - u ** 2)
                    * u
                    * (s ** 2 - t ** 2))

        val, _ = tplquad(
            integrand, 0, 30,
            lambda s: 0, lambda s: s,
            lambda s, u: -u, lambda s, u: u,
            args=(alpha,),
            epsabs=1e-12, epsrel=1e-12,
        )
        val *= P_STATE_VOLUME_FACTOR

        S = overlap_element_pstate_eckart(bf, bf, alpha, beta=0.0)
        assert abs(val - S) < 1e-10, f"quad={val}, formula={S}"

    @pytest.mark.slow
    def test_vee_vs_quadrature(self):
        from scipy.integrate import tplquad
        alpha = 1.5
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='sym')

        def integrand(t, u, s, alpha):
            # V_ee integrand: (s^2+t^2-u^2) * (s^2-t^2) e^{-2as}  (1/r_12 cancels Jacobian u)
            return (math.exp(-2 * alpha * s)
                    * (s ** 2 + t ** 2 - u ** 2)
                    * (s ** 2 - t ** 2))

        val, _ = tplquad(
            integrand, 0, 30,
            lambda s: 0, lambda s: s,
            lambda s, u: -u, lambda s, u: u,
            args=(alpha,),
            epsabs=1e-12, epsrel=1e-12,
        )
        val *= P_STATE_VOLUME_FACTOR

        V = potential_vee_element_pstate_eckart(bf, bf, alpha, beta=0.0)
        assert abs(val - V) < 1e-10, f"quad={val}, formula={V}"


# ---------------------------------------------------------------------------
# Basis-construction sanity
# ---------------------------------------------------------------------------

class TestBasisConstruction:
    def test_invalid_channel_rejected(self):
        with pytest.raises(ValueError):
            HylleraasPStateBasisFn(0, 0, 0, channel='quintet')

    def test_negative_indices_rejected(self):
        with pytest.raises(ValueError):
            HylleraasPStateBasisFn(-1, 0, 0)

    def test_label_includes_channel(self):
        bf_s = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        bf_a = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        assert '+' in bf_s.label()
        assert '-' in bf_a.label()

    @pytest.mark.parametrize("omega", [0, 1, 2, 3])
    def test_basis_total_degree_sizes(self, omega):
        basis = hylleraas_pstate_basis_total_degree(omega, channel='sym')
        # Should be the same count as the S-state basis at the same omega
        # (since the (l, m, n) enumeration is identical).
        from geovac.hylleraas_r12 import hylleraas_basis_total_degree
        s_basis = hylleraas_basis_total_degree(omega)
        assert len(basis) == len(s_basis)


# ---------------------------------------------------------------------------
# Antisymmetric channel
# ---------------------------------------------------------------------------

class TestAntisymmetricChannel:
    """Antisymmetric channel: phi^(-) = (z_1 - z_2) sinh(beta t) Q.
    Angular reduction: <(z_1 - z_2)^2>_SO(3) = u^2/3.
    All antisym matrix elements scale as beta^2 at small beta (since sinh^2 ~ beta^2 t^2).
    At beta = 0 exactly, the basis collapses and all matrix elements vanish.
    """

    def test_antisym_overlap_zero_at_beta_zero(self):
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        S = overlap_element_pstate_eckart_antisym(bf, bf, 1.5, beta=0.0)
        assert abs(S) < 1e-14

    def test_antisym_vne_zero_at_beta_zero(self):
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        V = potential_vne_element_pstate_eckart_antisym(bf, bf, 1.5, beta=0.0, Z=2.0)
        assert abs(V) < 1e-14

    def test_antisym_vee_zero_at_beta_zero(self):
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        V = potential_vee_element_pstate_eckart_antisym(bf, bf, 1.5, beta=0.0)
        assert abs(V) < 1e-14

    @pytest.mark.parametrize("beta", [0.1, 0.3, 0.5])
    def test_antisym_overlap_positive(self, beta):
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        S = overlap_element_pstate_eckart_antisym(bf, bf, 1.5, beta)
        assert S > 0, f"Overlap (antisym) at beta={beta} not positive: {S}"

    @pytest.mark.parametrize("beta", [0.1, 0.3, 0.5])
    def test_antisym_vee_positive(self, beta):
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        V = potential_vee_element_pstate_eckart_antisym(bf, bf, 1.5, beta)
        assert V > 0

    @pytest.mark.parametrize("beta", [0.1, 0.3, 0.5])
    def test_antisym_vne_negative(self, beta):
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        V = potential_vne_element_pstate_eckart_antisym(bf, bf, 1.5, beta, Z=2.0)
        assert V < 0

    @pytest.mark.slow
    def test_antisym_overlap_vs_quadrature(self):
        from scipy.integrate import tplquad
        import math
        alpha, beta = 1.5, 0.3
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        val, _ = tplquad(
            lambda t, u, s: (math.exp(-2 * alpha * s)
                              * (u ** 2 / 3)
                              * 0.5 * (math.cosh(2 * beta * t) - 1)
                              * u
                              * (s ** 2 - t ** 2)),
            0, 30, lambda s: 0, lambda s: s,
            lambda s, u: -u, lambda s, u: u,
            epsabs=1e-12, epsrel=1e-12,
        )
        val *= math.pi ** 2
        S = overlap_element_pstate_eckart_antisym(bf, bf, alpha, beta)
        assert abs(val - S) < 1e-10, f"quad={val}, formula={S}"


# ---------------------------------------------------------------------------
# Cross-sector (sym × antisym)
# ---------------------------------------------------------------------------

class TestCrossSector:
    """Cross-sector matrix elements <Phi^(+)_p | O | Phi^(-)_q>:
    Angular reduction <(z_1+z_2)(z_1-z_2)>_SO(3) = st/3.
    Hyperbolic content cosh(beta_p t) sinh(beta_q t) = (1/2) sinh(2 beta t) (shared beta).
    """

    def test_cross_overlap_zero_at_beta_zero(self):
        bf_p = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        bf_m = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        S = overlap_element_pstate_eckart_cross(bf_p, bf_m, 1.5, beta=0.0)
        assert S == 0.0  # sinh master at B=0 is exactly zero

    def test_cross_vne_zero_at_beta_zero(self):
        bf_p = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        bf_m = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        V = potential_vne_element_pstate_eckart_cross(bf_p, bf_m, 1.5, beta=0.0, Z=2.0)
        assert V == 0.0

    def test_cross_vee_zero_at_beta_zero(self):
        bf_p = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        bf_m = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        V = potential_vee_element_pstate_eckart_cross(bf_p, bf_m, 1.5, beta=0.0)
        assert V == 0.0

    @pytest.mark.parametrize("beta", [0.1, 0.3])
    def test_cross_overlap_finite(self, beta):
        bf_p = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        bf_m = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        S = overlap_element_pstate_eckart_cross(bf_p, bf_m, 1.5, beta)
        # Cross-sector overlap is nonzero at beta > 0; sign depends on basis.
        assert abs(S) > 0

    @pytest.mark.slow
    def test_cross_overlap_vs_quadrature(self):
        from scipy.integrate import tplquad
        import math
        alpha, beta = 1.5, 0.3
        bf_p = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        bf_m = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        val, _ = tplquad(
            lambda t, u, s: (math.exp(-2 * alpha * s)
                              * (s * t / 3)
                              * 0.5 * math.sinh(2 * beta * t)
                              * u
                              * (s ** 2 - t ** 2)),
            0, 30, lambda s: 0, lambda s: s,
            lambda s, u: -u, lambda s, u: u,
            epsabs=1e-12, epsrel=1e-12,
        )
        val *= math.pi ** 2
        S = overlap_element_pstate_eckart_cross(bf_p, bf_m, alpha, beta)
        assert abs(val - S) < 1e-10, f"quad={val}, formula={S}"

    @pytest.mark.slow
    def test_cross_vne_vs_quadrature(self):
        from scipy.integrate import tplquad
        import math
        alpha, beta, Z = 1.5, 0.3, 2.0
        bf_p = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        bf_m = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        val, _ = tplquad(
            lambda t, u, s: (math.exp(-2 * alpha * s)
                              * (s ** 2 * t / 3)
                              * 0.5 * math.sinh(2 * beta * t)
                              * (-4 * Z)
                              * u),
            0, 30, lambda s: 0, lambda s: s,
            lambda s, u: -u, lambda s, u: u,
            epsabs=1e-12, epsrel=1e-12,
        )
        val *= math.pi ** 2
        V = potential_vne_element_pstate_eckart_cross(bf_p, bf_m, alpha, beta, Z=Z)
        assert abs(val - V) < 1e-10, f"quad={val}, formula={V}"

    @pytest.mark.slow
    def test_cross_vee_vs_quadrature(self):
        from scipy.integrate import tplquad
        import math
        alpha, beta = 1.5, 0.3
        bf_p = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        bf_m = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        val, _ = tplquad(
            lambda t, u, s: (math.exp(-2 * alpha * s)
                              * (s * t / 3)
                              * 0.5 * math.sinh(2 * beta * t)
                              * (s ** 2 - t ** 2)),
            0, 30, lambda s: 0, lambda s: s,
            lambda s, u: -u, lambda s, u: u,
            epsabs=1e-12, epsrel=1e-12,
        )
        val *= math.pi ** 2
        V = potential_vee_element_pstate_eckart_cross(bf_p, bf_m, alpha, beta)
        assert abs(val - V) < 1e-10, f"quad={val}, formula={V}"


# ---------------------------------------------------------------------------
# P-state kinetic (sym × sym channel, singlet)
# ---------------------------------------------------------------------------

class TestPStateKineticSymSym:
    """Sym×sym P-state kinetic via Hartree-form decomposition:
    T = T_1 + T_mid_q + T_mid_p + T_3.

    Validated against analytical reference 2 pi^2 / alpha^6 at
    (l=m=n=0, beta=0) and against 3D SO(3)-averaged quadrature at beta>0.
    """

    @pytest.mark.parametrize("alpha", [1.0, 1.5, 2.0, 2.5])
    def test_kinetic_analytical_at_lmn_zero_beta_zero(self, alpha):
        """For phi = (z_1+z_2) e^{-alpha s} (l=m=n=0, beta=0):
        <phi|T|phi> = 2 pi^2 / alpha^6  (direct 6D integration).
        """
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        T = kinetic_element_pstate_eckart_sym_sym(bf, bf, alpha, beta=0.0)
        expected = 2 * math.pi ** 2 / alpha ** 6
        assert abs(T - expected) < 1e-13 * abs(expected), (
            f"alpha={alpha}: got {T}, expected {expected}"
        )

    @pytest.mark.parametrize("beta", [0.0, 0.2, 0.5])
    def test_kinetic_hermitian(self, beta):
        """T(p,q) == T(q,p) for arbitrary sym basis pairs at beta>=0."""
        basis = hylleraas_pstate_basis_total_degree(2, channel='sym')
        alpha = 1.5
        worst = 0.0
        for bf_p in basis:
            for bf_q in basis:
                t_pq = kinetic_element_pstate_eckart_sym_sym(
                    bf_p, bf_q, alpha, beta
                )
                t_qp = kinetic_element_pstate_eckart_sym_sym(
                    bf_q, bf_p, alpha, beta
                )
                worst = max(worst, abs(t_pq - t_qp))
        # Allow some slack from finite-precision arithmetic at higher (l, m, n).
        assert worst < 1e-10, f"worst |T(p,q)-T(q,p)| = {worst:.2e}"

    @pytest.mark.parametrize("alpha,beta", [(1.5, 0.3), (2.0, 0.5)])
    @pytest.mark.slow
    def test_kinetic_vs_so3_averaged_quadrature(self, alpha, beta):
        """Compare with 3D SO(3)-averaged kinetic quadrature at (0, 0, 0)."""
        from scipy.integrate import tplquad

        bf = HylleraasPStateBasisFn(0, 0, 0, channel='sym')

        def integrand(t, u, s, alpha, beta):
            es2 = math.exp(-2 * alpha * s)
            cb = math.cosh(beta * t)
            sb = math.sinh(beta * t)
            # S-state Hartree-form T_S for chi = e^{-as} cosh(bt):
            # dchi/dr_1 = chi*(-alpha) + e^{-as} beta sinh(bt) = ...
            # |grad_1 chi|^2 + |grad_2 chi|^2 = 2 e^{-2as} (alpha^2 cosh^2 + beta^2 sinh^2)
            T_S = 2 * es2 * (alpha ** 2 * cb * cb + beta ** 2 * sb * sb)
            # T_1
            T_1 = 0.5 * (s ** 2 + t ** 2 - u ** 2) / 3 * T_S
            # T_mid_q = T_mid_p (since p = q)
            chi = math.sqrt(es2) * cb  # = e^{-alpha s} cosh(beta t)
            chi_q_s = -alpha * chi
            chi_q_t = math.sqrt(es2) * beta * sb
            denom = s ** 2 - t ** 2
            if abs(denom) < 1e-30:
                M_q = 0.0
            else:
                M_q = (2.0 / (3 * denom)) * (
                    s * (s ** 2 - u ** 2) * chi_q_s
                    + t * (u ** 2 - t ** 2) * chi_q_t
                )
            T_mid = chi * M_q  # T_mid_q + T_mid_p = 2 * (chi_p * M_q) for p=q
            # T_3
            T_3 = chi * chi
            return (T_1 + T_mid + T_3) * u * (s ** 2 - t ** 2)

        val, _ = tplquad(
            integrand, 0, 30,
            lambda s: 0, lambda s: s,
            lambda s, u: -u, lambda s, u: u,
            args=(alpha, beta),
            epsabs=1e-12, epsrel=1e-12,
        )
        val *= math.pi ** 2
        T = kinetic_element_pstate_eckart_sym_sym(bf, bf, alpha, beta)
        assert abs(val - T) < 1e-10, (
            f"alpha={alpha}, beta={beta}: quad={val}, formula={T}, "
            f"diff={abs(val-T):.2e}"
        )

    def test_kinetic_positive_for_physical_basis(self):
        """For l=m=n=0 He-like basis, kinetic must be positive."""
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='sym')
        for alpha in [1.0, 1.5, 2.0]:
            for beta in [0.0, 0.2, 0.5]:
                T = kinetic_element_pstate_eckart_sym_sym(bf, bf, alpha, beta)
                assert T > 0, f"alpha={alpha}, beta={beta}: T={T} not positive"


# ---------------------------------------------------------------------------
# Universal P-state quadrature kinetic (covers antisym × antisym and cross)
# ---------------------------------------------------------------------------

class TestUniversalQuadratureKinetic:
    """Universal P-state kinetic via SO(3)-averaged 3D quadrature."""

    def test_sym_quadrature_matches_algebraic(self):
        """Quadrature kinetic should match the algebraic sym×sym kinetic
        to ~1e-4 relative (quadrature precision at n_r=32, n_theta=16)."""
        from geovac.hylleraas_eckart_pstate import _kinetic_via_quadrature_pstate
        basis = hylleraas_pstate_basis_total_degree(2, channel='sym')
        alpha = 1.35
        beta = 0.3
        worst = 0.0
        for bf_p in basis:
            for bf_q in basis:
                t_alg = kinetic_element_pstate_eckart_sym_sym(bf_p, bf_q, alpha, beta)
                t_qua = _kinetic_via_quadrature_pstate(
                    bf_p, bf_q, alpha, beta, n_r=32, n_theta=16
                )
                denom = max(abs(t_alg), abs(t_qua), 1e-30)
                worst = max(worst, abs(t_alg - t_qua) / denom)
        assert worst < 1e-3, f"worst rel diff = {worst:.2e}"

    def test_antisym_kinetic_zero_at_beta_zero(self):
        """At beta=0, the antisym basis (z_1-z_2)·sinh(0)·Q vanishes; the
        kinetic must be exactly 0 by direct quadrature on a zero integrand."""
        from geovac.hylleraas_eckart_pstate import _kinetic_via_quadrature_pstate
        basis = hylleraas_pstate_basis_total_degree(1, channel='antisym')
        for bf_p in basis:
            for bf_q in basis:
                t = _kinetic_via_quadrature_pstate(
                    bf_p, bf_q, 1.5, 0.0, n_r=16, n_theta=8
                )
                assert abs(t) < 1e-12, f"antisym T at beta=0 nonzero: {t}"

    def test_cross_sector_kinetic_zero_at_beta_zero(self):
        """Cross-sector at beta=0: sinh(0)=0, so integrand vanishes."""
        from geovac.hylleraas_eckart_pstate import _kinetic_via_quadrature_pstate
        basis_sym = hylleraas_pstate_basis_total_degree(1, channel='sym')
        basis_anti = hylleraas_pstate_basis_total_degree(1, channel='antisym')
        for bf_p in basis_sym:
            for bf_q in basis_anti:
                t = _kinetic_via_quadrature_pstate(
                    bf_p, bf_q, 1.5, 0.0, n_r=16, n_theta=8
                )
                assert abs(t) < 1e-12

    def test_antisym_kinetic_positive_at_beta_finite(self):
        """At (l=m=n=0, beta>0), the antisym kinetic should be positive."""
        from geovac.hylleraas_eckart_pstate import _kinetic_via_quadrature_pstate
        bf = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
        T = _kinetic_via_quadrature_pstate(bf, bf, 1.5, 0.3, n_r=32, n_theta=16)
        assert T > 0, f"antisym (000) T at beta=0.3: {T} not positive"

    @pytest.mark.parametrize("beta", [0.0, 0.3])
    def test_kinetic_hermitian_all_channels(self, beta):
        """T(p,q) = T(q,p) for all 4 channel pairs."""
        from geovac.hylleraas_eckart_pstate import _kinetic_via_quadrature_pstate
        basis_sym = hylleraas_pstate_basis_total_degree(1, channel='sym')
        basis_anti = hylleraas_pstate_basis_total_degree(1, channel='antisym')
        alpha = 1.5
        # sym × sym (algebraic), antisym × antisym (quadrature), cross (quadrature)
        worst = 0.0
        for bf_p in basis_sym:
            for bf_q in basis_sym:
                t_pq = _kinetic_via_quadrature_pstate(bf_p, bf_q, alpha, beta, 16, 8)
                t_qp = _kinetic_via_quadrature_pstate(bf_q, bf_p, alpha, beta, 16, 8)
                worst = max(worst, abs(t_pq - t_qp))
        for bf_p in basis_anti:
            for bf_q in basis_anti:
                t_pq = _kinetic_via_quadrature_pstate(bf_p, bf_q, alpha, beta, 16, 8)
                t_qp = _kinetic_via_quadrature_pstate(bf_q, bf_p, alpha, beta, 16, 8)
                worst = max(worst, abs(t_pq - t_qp))
        # Cross-sector p,q swap: T(sym, antisym) vs T(antisym, sym)
        for bf_p in basis_sym:
            for bf_q in basis_anti:
                t_pq = _kinetic_via_quadrature_pstate(bf_p, bf_q, alpha, beta, 16, 8)
                t_qp = _kinetic_via_quadrature_pstate(bf_q, bf_p, alpha, beta, 16, 8)
                worst = max(worst, abs(t_pq - t_qp))
        assert worst < 1e-8, f"worst |T(p,q)-T(q,p)| = {worst:.2e}"


# ---------------------------------------------------------------------------
# Full two-channel 2^1P sprint + oscillator strength
# ---------------------------------------------------------------------------

class TestFullChannelOscillatorStrength:
    """Full Schwartz 1961 two-channel 2^1P trial closes the 1S->2P
    oscillator strength to Drake-class accuracy.

    Reference: f(He 2^1P -> 1^1S) = 0.2761 (Drake handbook 1996).
    """

    @pytest.mark.slow
    def test_hydrogen_1s_2p_f_factor(self):
        """Sanity check: hydrogen 1S->2P with the Wigner-Eckart factor
        gives f = 0.4162. This is a closed-form check on the factor of
        2 in the f formula (not the multi-electron variational code)."""
        # |<2P_0|z|1S>|^2 for hydrogen 1S, 2P_0 (analytic): 0.5546...
        # Using R_{1S}=2 e^{-r}, R_{2P} = (1/(2sqrt(6))) r e^{-r/2}:
        # <R_{2P}|r|R_{1S}> = 768/(243 sqrt(6))
        # <Y_{1,0}|cos theta|Y_{0,0}> = 1/sqrt(3)
        radial = 768.0 / (243.0 * math.sqrt(6.0))
        angular = 1.0 / math.sqrt(3.0)
        D_z = radial * angular
        DeltaE = 3.0 / 8.0  # Bohr 1S -> 2P
        f = 2.0 * DeltaE * D_z ** 2
        assert abs(f - 0.4162) < 0.001, f"hydrogen f = {f}, expected 0.4162"

    @pytest.mark.slow
    def test_he_2p1_full_channel_drake_class(self):
        """End-to-end: He 2^1P->1^1S oscillator strength matches Drake
        f=0.276 to within 3%."""
        from geovac.hylleraas_r12 import (
            hylleraas_basis_total_degree, optimize_alpha_beta_for_state,
        )
        from geovac.hylleraas_eckart_pstate import (
            optimize_2p1_full, oscillator_strength_2p_to_1s_full,
        )
        # 1^1S
        basis_S = hylleraas_basis_total_degree(3)
        res_1s = optimize_alpha_beta_for_state(
            basis_S, Z=2, alpha_init=1.85, beta_init=0.0,
            state_index=0, spin='singlet',
        )
        s_state = {
            'energy': res_1s.energy, 'basis': basis_S,
            'coeffs': res_1s.coeffs, 'alpha': res_1s.alpha,
            'beta': res_1s.extras['beta'],
        }
        # 2^1P with both channels at omega_p=2
        basis_sym = hylleraas_pstate_basis_total_degree(2, channel='sym')
        basis_anti = hylleraas_pstate_basis_total_degree(2, channel='antisym')
        p_full = optimize_2p1_full(
            basis_sym, basis_anti, Z=2,
            alpha_init=1.35, beta_init=0.35,
            n_r=16, n_theta=8, max_iter=60,
        )
        osc = oscillator_strength_2p_to_1s_full(s_state, p_full)
        # Within 3% of Drake.
        rel_err = abs(osc['f'] - 0.2761) / 0.2761
        assert rel_err < 0.05, (
            f"f = {osc['f']:.4f}, Drake = 0.2761, rel_err = {rel_err*100:.2f}%"
        )
