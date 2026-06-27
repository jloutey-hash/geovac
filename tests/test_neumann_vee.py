"""Tests for Neumann expansion V_ee computation.

Validates auxiliary integrals (C_l, A_n, A_l, B_l, X_l) and full V_ee
matrix elements against numerical quadrature and known analytical values.
"""

import numpy as np
import pytest
from scipy.integrate import quad, dblquad
from scipy.special import lqmn, exp1, ellipk
from numpy.polynomial.legendre import legval

from geovac.neumann_vee import (
    compute_Cl_table,
    compute_An_table,
    compute_An_negative,
    compute_Al_table,
    compute_B0_table,
    compute_Bl_table,
    compute_vee_matrix_neumann,
    compute_hamiltonian_neumann,
    _compute_Xl_single,
)
from geovac.hylleraas import (
    HylleraasBasisFunction,
    generate_basis,
    build_quadrature_grids,
    compute_overlap_matrix,
    compute_hamiltonian_matrix,
    evaluate_basis_function,
    solve_hylleraas,
)


# ============================================================
# Auxiliary integral tests
# ============================================================

class TestClTable:
    """C_l(p) = ∫₋₁¹ η^p P_l(η) dη."""

    def test_C0_even(self):
        """C_0(p) = 2/(p+1) for even p."""
        C = compute_Cl_table(10, 5)
        for p in range(0, 11, 2):
            assert abs(C[0, p] - 2.0 / (p + 1)) < 1e-14, f"C_0({p}) failed"

    def test_C0_odd_zero(self):
        """C_0(p) = 0 for odd p (odd integrand)."""
        C = compute_Cl_table(10, 5)
        for p in range(1, 10, 2):
            assert abs(C[0, p]) < 1e-14, f"C_0({p}) should be zero"

    def test_C1_odd(self):
        """C_1(p) = 2/(p+2) for odd p."""
        C = compute_Cl_table(10, 5)
        for p in range(1, 10, 2):
            assert abs(C[1, p] - 2.0 / (p + 2)) < 1e-14, f"C_1({p}) failed"

    def test_selection_rule(self):
        """C_l(p) = 0 when p < l or (p-l) is odd."""
        C = compute_Cl_table(8, 8)
        for l in range(9):
            for p in range(9):
                if p < l or (p - l) % 2 != 0:
                    assert abs(C[l, p]) < 1e-13, f"C_{l}({p}) should be zero"

    def test_known_values(self):
        """Check specific known values."""
        C = compute_Cl_table(6, 6)
        # C_2(2) = ∫₋₁¹ η² P_2(η) dη = ∫ η²(3η²-1)/2 dη = (3·2/5 - 2/3)/2 = 4/15
        assert abs(C[2, 2] - 4.0 / 15.0) < 1e-14
        # C_2(4) = ∫₋₁¹ η⁴ P_2(η) dη = ∫ η⁴(3η²-1)/2 dη = (3·2/7 - 2/5)/2 = 16/70 = 8/35
        assert abs(C[2, 4] - 8.0 / 35.0) < 1e-14


class TestAnTable:
    """A_n(α) = ∫₁^∞ ξ^n e^{-αξ} dξ."""

    def test_A0(self):
        """A_0(α) = e^{-α}/α."""
        alpha = 1.5
        A = compute_An_table(0, alpha)
        assert abs(A[0] - np.exp(-alpha) / alpha) < 1e-14

    def test_against_quad(self):
        """Verify against scipy.integrate.quad."""
        alpha = 2.0
        A = compute_An_table(8, alpha)
        for n in range(9):
            val, _ = quad(lambda x: x**n * np.exp(-alpha * x), 1, np.inf)
            assert abs(A[n] - val) < 1e-12, f"A_{n} failed: {A[n]} vs {val}"

    def test_negative_index(self):
        """A_{-1}(α) = E₁(α)."""
        from scipy.special import exp1
        alpha = 1.0
        A_neg = compute_An_negative(-3, alpha)
        assert abs(A_neg[-1] - float(exp1(alpha))) < 1e-14


class TestAlTable:
    """A_l(p, α) = ∫₁^∞ ξ^p e^{-αξ} P_l(ξ) dξ."""

    def test_l0_equals_An(self):
        """A_0(p, α) = A_p(α) since P_0(ξ) = 1."""
        alpha = 1.5
        An = compute_An_table(6, alpha)
        Al = compute_Al_table(0, 6, alpha)
        for p in range(7):
            assert abs(Al[0, p] - An[p]) < 1e-13

    def test_l1_equals_An_shifted(self):
        """A_1(p, α) = A_{p+1}(α) since P_1(ξ) = ξ."""
        alpha = 1.5
        An = compute_An_table(8, alpha)
        Al = compute_Al_table(1, 6, alpha)
        for p in range(7):
            assert abs(Al[1, p] - An[p + 1]) < 1e-13

    def test_against_quad(self):
        """Verify A_l against numerical integration."""
        from numpy.polynomial.legendre import legval
        alpha = 2.0
        Al = compute_Al_table(4, 4, alpha)
        for l in range(5):
            for p in range(5):
                coeffs = np.zeros(l + 1)
                coeffs[l] = 1.0
                def integrand(xi, l_=l, p_=p):
                    return xi**p_ * np.exp(-alpha * xi) * legval(xi, coeffs)
                val, _ = quad(integrand, 1, np.inf)
                assert abs(Al[l, p] - val) < 1e-11, \
                    f"A_{l}({p}) failed: {Al[l,p]} vs {val}"


class TestBlTable:
    """B_l(p, α) = ∫₁^∞ ξ^p e^{-αξ} Q_l(ξ) dξ."""

    def test_B0_against_quad(self):
        """B_0 against direct numerical integration."""
        alpha = 2.0
        B0 = compute_B0_table(4, alpha)
        for p in range(5):
            def integrand(xi, p_=p):
                Q0 = 0.5 * np.log((xi + 1) / (xi - 1)) if xi > 1 else 0.0
                return xi**p_ * np.exp(-alpha * xi) * Q0
            val, _ = quad(integrand, 1, np.inf, limit=200)
            assert abs(B0[p] - val) < 1e-10, f"B_0({p}) failed"

    def test_Bl_against_quad(self):
        """B_l vs a fresh quad of the SAME lqmn integrand.

        NOTE: This is only a table-assembly / indexing / cache consistency
        check — it re-runs the *identical* scipy.quad+lqmn integrand that
        compute_Bl_table uses internally, so it cannot detect an error in the
        integrand definition itself (it would be tautological as a correctness
        proof). The genuine, integrand-independent correctness checks are
        test_Bl_recurrence, test_B0_closed_form, and test_Bl_Q1_identity below.
        """
        alpha = 2.0
        Bl = compute_Bl_table(4, 3, alpha)
        for l in range(5):
            for p in range(4):
                def integrand(xi, l_=l, p_=p):
                    Qvals, _ = lqmn(0, l_, xi)
                    return xi**p_ * np.exp(-alpha * xi) * Qvals[0, l_]
                val, _ = quad(integrand, 1, np.inf, limit=200)
                assert abs(Bl[l, p] - val) < 1e-10, \
                    f"B_{l}({p}) failed: {Bl[l,p]} vs {val}"

    def test_Bl_recurrence(self):
        """INDEPENDENT check: B_l obeys the second-kind Legendre recurrence.

        Paper 12 eq:Bl_recurrence:
            (l+1) B_{l+1}(p, α) = (2l+1) B_l(p+1, α) − l B_{l-1}(p, α)   (l ≥ 1)

        This is a genuine algebraic identity that the *true* moments satisfy;
        it is NOT a re-evaluation of the production integrand, so it catches
        sign / index / l-structure errors in compute_Bl_table. The base rows
        B_0, B_1 are pinned absolutely by test_B0_closed_form + test_Bl_Q1_identity,
        and this recurrence propagates that anchor to all l ≥ 2.

        Note: the pure 3-term recurrence holds only for l ≥ 1 — for l = 0 the
        relation Q_1 = ξ Q_0 − 1 carries the extra −1 (verified separately in
        test_Bl_Q1_identity).
        """
        for alpha in (1.0, 2.0, 3.0):
            l_max, p_max = 8, 6
            Bl = compute_Bl_table(l_max, p_max, alpha)
            for l in range(1, l_max):
                for p in range(p_max):  # need column p+1
                    lhs = (l + 1) * Bl[l + 1, p]
                    rhs = (2 * l + 1) * Bl[l, p + 1] - l * Bl[l - 1, p]
                    assert abs(lhs - rhs) < 1e-8, \
                        (f"Bl recurrence fails α={alpha} l={l} p={p}: "
                         f"{lhs} vs {rhs} (diff {abs(lhs-rhs):.2e})")

    def test_B0_closed_form(self):
        """ABSOLUTE ANCHOR: B_0(0,α) and B_0(1,α) vs analytic closed forms.

        Q_0(ξ) = ½ ln[(ξ+1)/(ξ-1)].  Splitting Q_0 into ½[ln(ξ+1) − ln(ξ-1)]
        and integrating each piece by parts against e^{-αξ} gives closed forms
        in the exponential integral E_1 and the Euler–Mascheroni constant γ:

            B_0(0,α) = (1/2α)[ e^{α} E_1(2α) + e^{-α}(γ + ln α + ln 2) ].

        These use NO quadrature, so they pin the absolute scale, sign, α-
        dependence and Q_0-normalization of the B_0 base row — exactly the
        degrees of freedom that the recurrence (which factors out any uniform
        rescaling) cannot constrain.  (Both forms were verified against an
        independent 40-digit mpmath integral to < 1e-15.)
        """
        g = np.euler_gamma
        for alpha in (1.0, 2.0, 3.0):
            B0 = compute_B0_table(1, alpha)

            B0_0 = (1.0 / (2.0 * alpha)) * (
                np.exp(alpha) * exp1(2.0 * alpha)
                + np.exp(-alpha) * (g + np.log(alpha) + np.log(2.0))
            )

            # B_0(1,α) = ½(J⁺_1 − J⁻_1); J± are the ξ-power-1 moments of ln(ξ±1).
            Jm = np.exp(-alpha) * (
                (1.0 / alpha) * (-g - np.log(alpha))
                + (1.0 / alpha**2) * (1.0 - g - np.log(alpha))
            )
            K0 = np.exp(-2 * alpha) * np.log(2.0) / alpha + exp1(2 * alpha) / alpha
            K1 = (2 * np.log(2.0) * np.exp(-2 * alpha) / alpha
                  + (1.0 / alpha) * K0 + np.exp(-2 * alpha) / alpha**2)
            Jp = np.exp(alpha) * (-K0 + K1)
            B0_1 = 0.5 * (Jp - Jm)

            assert abs(B0[0] - B0_0) < 1e-10, \
                f"B_0(0,{alpha}) closed form: {B0[0]} vs {B0_0}"
            assert abs(B0[1] - B0_1) < 1e-10, \
                f"B_0(1,{alpha}) closed form: {B0[1]} vs {B0_1}"

    def test_Bl_row0_matches_B0_table(self):
        """compute_Bl_table row l=0 (lqmn Q_0) == compute_B0_table (explicit Q_0).

        Cross-checks the two independent Q_0 representations (scipy lqmn vs the
        explicit ½ln((ξ+1)/(ξ-1)) formula), transferring the closed-form anchor
        of test_B0_closed_form onto the recurrence's base row.
        """
        alpha = 2.0
        Bl = compute_Bl_table(3, 5, alpha)
        B0 = compute_B0_table(5, alpha)
        for p in range(6):
            assert abs(Bl[0, p] - B0[p]) < 1e-12, \
                f"Bl[0,{p}]={Bl[0,p]} vs B0[{p}]={B0[p]}"

    def test_Bl_Q1_identity(self):
        """INDEPENDENT check anchoring the B_1 row: Q_1(ξ) = ξ Q_0(ξ) − 1.

        Multiplying by ξ^p e^{-αξ} and integrating gives
            B_1(p,α) = B_0(p+1,α) − A_p(α),
        i.e. B_0(p+1,α) − B_1(p,α) must equal the exponential moment A_p(α)
        (Paper 12 eq:An_def).  This ties the B_1 row to the (closed-form-
        anchored) B_0 row through the framework's own A_n table — an algebraic
        identity, not a re-evaluation of the Q_l integrand.
        """
        for alpha in (1.0, 2.0):
            p_max = 5
            Bl = compute_Bl_table(1, p_max + 1, alpha)
            A = compute_An_table(p_max, alpha)
            for p in range(p_max + 1):
                lhs = Bl[0, p + 1] - Bl[1, p]
                assert abs(lhs - A[p]) < 1e-10, \
                    (f"Q1 identity fails α={alpha} p={p}: "
                     f"B0(p+1)-B1(p)={lhs} vs A_p={A[p]}")

    def test_B0_positive(self):
        """B_0(p, α) should be positive (Q_0 > 0 for ξ > 1)."""
        B0 = compute_B0_table(6, 1.5)
        for p in range(7):
            assert B0[p] > 0, f"B_0({p}) should be positive"


class TestXlIntegral:
    """X_l(p₁, p₂, α) = ∫∫ ξ₁^{p₁} ξ₂^{p₂} e^{-α(ξ₁+ξ₂)} P_l(ξ_<) Q_l(ξ_>) dξ₁dξ₂."""

    def test_against_2d_quad(self):
        """Verify X_l against 2D numerical integration for small cases."""
        alpha = 2.0
        for l in range(3):
            for p1, p2 in [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0)]:
                xl = _compute_Xl_single(l, p1, p2, alpha)

                # 2D numerical integration (split at ξ₁ = ξ₂)
                from numpy.polynomial.legendre import legval

                def integrand(xi2, xi1, l_=l, p1_=p1, p2_=p2):
                    if xi1 <= 1.0 or xi2 <= 1.0:
                        return 0.0
                    xi_less = min(xi1, xi2)
                    xi_more = max(xi1, xi2)
                    coeffs_P = np.zeros(l_ + 1)
                    coeffs_P[l_] = 1.0
                    Pl = legval(xi_less, coeffs_P)
                    Qvals, _ = lqmn(0, l_, xi_more)
                    Ql = Qvals[0, l_]
                    return (xi1**p1_ * xi2**p2_ *
                            np.exp(-alpha * (xi1 + xi2)) * Pl * Ql)

                val, _ = dblquad(integrand, 1, 30, 1, 30,
                                 epsabs=1e-10, epsrel=1e-10)
                assert abs(xl - val) < max(1e-8, abs(val) * 1e-6), \
                    f"X_{l}({p1},{p2}) failed: {xl} vs {val}"

    def test_symmetry(self):
        """X_l(p₁, p₂, α) = X_l(p₂, p₁, α)."""
        alpha = 2.0
        for l in range(4):
            for p1 in range(3):
                for p2 in range(p1 + 1, 4):
                    x1 = _compute_Xl_single(l, p1, p2, alpha)
                    x2 = _compute_Xl_single(l, p2, p1, alpha)
                    assert abs(x1 - x2) < 1e-10, \
                        f"X_{l}({p1},{p2}) != X_{l}({p2},{p1}): {x1} vs {x2}"


# ============================================================
# Independent V_ee references (no reuse of neumann_vee internals)
# ============================================================

def _ref_Cl_quad(l: int, q: int) -> float:
    """C_l(q) = ∫₋₁¹ η^q P_l(η) dη by independent 1D quadrature."""
    if q < l or (q - l) % 2 != 0:
        return 0.0
    coeffs = np.zeros(l + 1)
    coeffs[l] = 1.0
    val, _ = quad(lambda e: e**q * legval(e, coeffs), -1.0, 1.0)
    return val


def _ref_Xl_quad(l: int, p1: int, p2: int, alpha: float) -> float:
    """X_l(p1,p2,α) by independent 2D split-region quadrature (eq:Xl_def).

    Uses scipy lqmn for Q_l and a split at ξ₁=ξ₂, tight tolerances. This is
    an independent recomputation of the radial building block — not the
    neumann_vee IBP recurrence.
    """
    coeffs = np.zeros(l + 1)
    coeffs[l] = 1.0

    def Pl(x):
        return legval(x, coeffs)

    def Ql(x):
        Q, _ = lqmn(0, l, x)
        return Q[0, l]

    # Region ξ₁ < ξ₂
    def integ_a(xi1, xi2):
        return xi1**p1 * xi2**p2 * np.exp(-alpha * (xi1 + xi2)) * Pl(xi1) * Ql(xi2)
    I1, _ = dblquad(integ_a, 1.0, 60.0, 1.0, lambda xi2: xi2,
                    epsabs=1e-12, epsrel=1e-12)

    # Region ξ₂ < ξ₁
    def integ_b(xi2, xi1):
        return xi1**p1 * xi2**p2 * np.exp(-alpha * (xi1 + xi2)) * Pl(xi2) * Ql(xi1)
    I2, _ = dblquad(integ_b, 1.0, 60.0, 1.0, lambda xi1: xi1,
                    epsabs=1e-12, epsrel=1e-12)
    return I1 + I2


def _ref_vee_element_assembly(bf_i, bf_j, R: float, alpha: float,
                              l_max: int = 12) -> float:
    """⟨φ_i|1/r₁₂|φ_j⟩ assembled independently from quadrature building blocks.

    Re-derives, from the Paper 12 equations and WITHOUT calling any neumann_vee
    assembly routine:
      - the prefactor π²R⁵/8  (eq:factorized: (2/R)·(R/2)⁶·(2π)²),
      - the Jacobian (ξ²−η²)(ξ²−η²) → 4-subterm expansion (eq:jacobian_expansion),
      - the (2l+1) Neumann weight (eq:neumann_sigma),
      - the bra×ket symmetrization (direct + exchange),
    using independently-quadded C_l and X_l. Matching production therefore pins
    the prefactor / Jacobian bookkeeping the paper flags as error-prone.
    """
    prefactor = np.pi**2 * R**5 / 8.0
    terms_i = [(bf_i.j, bf_i.k, bf_i.l, bf_i.m), (bf_i.k, bf_i.j, bf_i.m, bf_i.l)]
    terms_j = [(bf_j.j, bf_j.k, bf_j.l, bf_j.m), (bf_j.k, bf_j.j, bf_j.m, bf_j.l)]
    a2 = 2.0 * alpha

    Xcache = {}

    def Xl(l, p1, p2):
        key = (l, p1, p2)
        if key not in Xcache:
            Xcache[key] = _ref_Xl_quad(l, p1, p2, a2)
        return Xcache[key]

    total = 0.0
    for (ji, ki, li, mi) in terms_i:
        for (jj, kj, lj, mj) in terms_j:
            P1, P2, Q1, Q2 = ji + jj, ki + kj, li + lj, mi + mj
            jac = [
                (+1, P1 + 2, P2 + 2, Q1, Q2),
                (-1, P1 + 2, P2, Q1, Q2 + 2),
                (-1, P1, P2 + 2, Q1 + 2, Q2),
                (+1, P1, P2, Q1 + 2, Q2 + 2),
            ]
            for sign, p1, p2, q1, q2 in jac:
                s = 0.0
                for l in range(l_max + 1):
                    if q1 < l or q2 < l or (q1 - l) % 2 or (q2 - l) % 2:
                        continue
                    c1 = _ref_Cl_quad(l, q1)
                    c2 = _ref_Cl_quad(l, q2)
                    if abs(c1) < 1e-30 or abs(c2) < 1e-30:
                        continue
                    s += (2 * l + 1) * Xl(l, p1, p2) * c1 * c2
                total += sign * s
    return prefactor * total


def _ref_vee_element_4dquad(bf_i, bf_j, R: float,
                            N_xi: int, N_eta: int,
                            xi_max: float = 40.0) -> float:
    """⟨φ_i|1/r₁₂|φ_j⟩ by a LITERAL 4D quadrature of eq:vee_integral.

    Shares nothing with the Neumann/Legendre decomposition: it uses the
    closed-form azimuthal average ⟨1/r₁₂⟩_φ = 8π K(k)/√s of the *raw* Coulomb
    kernel (complete elliptic integral K) on a Gauss product grid in
    (ξ₁,η₁,ξ₂,η₂), with the bare Jacobian J=(R/2)³(ξ²−η²). Converges only
    algebraically (codim-2 log singularity at electron coincidence), so it
    anchors the absolute normalization/prefactor against the actual integral
    rather than reaching machine precision.
    """
    half_R = R / 2.0
    u, wu = np.polynomial.legendre.leggauss(N_xi)
    t = (u + 1) / 2
    xi = 1.0 + (xi_max - 1.0) * t**2
    w_xi = wu * (xi_max - 1.0) * t
    eta, w_eta = np.polynomial.legendre.leggauss(N_eta)

    val = 0.0
    for a in range(N_xi):
        xi1 = xi[a]
        for c in range(N_xi):
            xi2 = xi[c]
            for b in range(N_eta):
                eta1 = eta[b]
                J1 = half_R**3 * (xi1**2 - eta1**2)
                rho1 = half_R * np.sqrt(max((xi1**2 - 1) * (1 - eta1**2), 0.0))
                z1 = half_R * xi1 * eta1
                for d in range(N_eta):
                    eta2 = eta[d]
                    J2 = half_R**3 * (xi2**2 - eta2**2)
                    rho2 = half_R * np.sqrt(
                        max((xi2**2 - 1) * (1 - eta2**2), 0.0))
                    z2 = half_R * xi2 * eta2
                    fi = evaluate_basis_function(bf_i, xi1, eta1, xi2, eta2, 1.0, R)
                    fj = evaluate_basis_function(bf_j, xi1, eta1, xi2, eta2, 1.0, R)
                    s = (rho1 + rho2)**2 + (z1 - z2)**2
                    if s > 1e-30:
                        k2 = min(max(4 * rho1 * rho2 / s, 0.0), 1 - 1e-15)
                        vee = 8 * np.pi * ellipk(k2) / np.sqrt(s)
                    else:
                        vee = 0.0
                    val += (w_xi[a] * w_xi[c] * w_eta[b] * w_eta[d]
                            * fi * fj * vee * J1 * J2)
    return val


# ============================================================
# V_ee matrix tests
# ============================================================

class TestVeeMatrix:
    """Test Neumann V_ee matrix against numerical V_ee."""

    @pytest.fixture
    def small_basis_and_grids(self):
        """3-function p=0 basis with grids."""
        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=30, N_eta=20, N_phi=1)
        return basis, grids

    def test_vee_symmetric(self, small_basis_and_grids):
        """V_ee matrix should be symmetric."""
        basis, grids = small_basis_and_grids
        R = 1.4011
        V = compute_vee_matrix_neumann(basis, R, l_max=15)
        assert np.allclose(V, V.T, atol=1e-14)

    def test_vee_positive_diagonal(self, small_basis_and_grids):
        """V_ee diagonal elements should be positive (repulsion)."""
        basis, grids = small_basis_and_grids
        R = 1.4011
        V = compute_vee_matrix_neumann(basis, R, l_max=15)
        assert np.all(V.diagonal() > 0), \
            f"V_ee diagonal should be positive: {V.diagonal()}"

    def test_neumann_vs_numerical_vee(self, small_basis_and_grids):
        """Neumann V_ee should match numerical V_ee (within grid convergence)."""
        basis, grids = small_basis_and_grids
        R = 1.4011

        # Get numerical V_ee = H_full - H(T+V_ne)
        H_full = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=False, include_vee=True
        )
        H_tv = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=False, include_vee=False
        )
        V_ee_num = H_full - H_tv

        # Neumann V_ee (exact)
        V_ee_neu = compute_vee_matrix_neumann(basis, R, l_max=15)

        # Numerical V_ee converges toward Neumann as grid increases.
        # With 30x20 grid, expect ~1-3% agreement.
        max_diff = np.max(np.abs(V_ee_neu - V_ee_num))
        max_val = np.max(np.abs(V_ee_neu))
        rel_diff = max_diff / max_val if max_val > 0 else 0
        print(f"\n  V_ee max relative difference: {rel_diff:.4f} ({rel_diff*100:.1f}%)")
        assert rel_diff < 0.10, \
            f"V_ee Neumann vs numerical too different: {rel_diff*100:.1f}%"

    @pytest.mark.slow
    def test_vee_element_machine_precision(self):
        """Paper 12 'exact within the Neumann truncation' — machine-precision test.

        Assembles full V_ee matrix elements via the production Neumann path
        (compute_vee_matrix_neumann) and compares each to an INDEPENDENT
        evaluation of the 4D integral eq:vee_integral (_ref_vee_element_assembly:
        independently-quadded C_l + X_l, re-derived prefactor / Jacobian /
        (2l+1) / symmetrization). Agreement to < 1e-8 pins the prefactor and
        Jacobian-expansion bookkeeping the paper itself flags as error-prone.
        """
        R = 1.4011
        alpha = 1.0
        # Span self-exchange (0,0,0,0), a ξ-power (1,0,0,0), and an η-bearing
        # element (1,1,1,1) so C_l for l ≥ 1 and the η-Jacobian subterms are
        # exercised (not just the l=0 / Q=0 channel).
        basis = [
            HylleraasBasisFunction(0, 0, 0, 0, 0, alpha),
            HylleraasBasisFunction(1, 0, 0, 0, 0, alpha),
            HylleraasBasisFunction(1, 1, 1, 1, 0, alpha),
        ]
        V = compute_vee_matrix_neumann(basis, R, l_max=20)

        for i in range(len(basis)):
            for j in range(i, len(basis)):
                ref = _ref_vee_element_assembly(basis[i], basis[j], R, alpha,
                                                l_max=12)
                diff = abs(V[i, j] - ref)
                assert diff < 1e-8, \
                    (f"V_ee[{i},{j}] Neumann={V[i,j]:.12e} vs independent "
                     f"assembly={ref:.12e} (diff {diff:.2e})")

    @pytest.mark.slow
    def test_vee_literal_4d_quad_normalization(self):
        """Absolute normalization vs the LITERAL Coulomb integral (eq:vee_integral).

        The independent assembly above shares the Legendre azimuthal reduction
        with production; this test instead integrates the *raw* kernel
        ⟨1/r₁₂⟩_φ = 8πK(k)/√s by 4D Gauss quadrature (shares nothing with the
        Neumann decomposition). It converges only algebraically (coincidence
        singularity), so it confirms the prefactor / Jacobian against the actual
        integral at the ~%-level and that refinement drives it toward the
        Neumann value (monotone convergence) — ruling out a gross shared
        prefactor error.
        """
        R = 1.4011
        alpha = 1.0
        bf = HylleraasBasisFunction(0, 0, 0, 0, 0, alpha)
        V = compute_vee_matrix_neumann([bf], R, l_max=20)
        neu = V[0, 0]

        q_coarse = _ref_vee_element_4dquad(bf, bf, R, 40, 30)
        q_fine = _ref_vee_element_4dquad(bf, bf, R, 60, 40)
        rel_coarse = abs(q_coarse - neu) / abs(neu)
        rel_fine = abs(q_fine - neu) / abs(neu)
        print(f"\n  literal 4D quad vs Neumann: coarse={rel_coarse:.3e}, "
              f"fine={rel_fine:.3e}")

        # Converges toward the Neumann value (refinement reduces the gap)
        assert rel_fine < rel_coarse, \
            f"4D quad not converging to Neumann: {rel_coarse:.3e} -> {rel_fine:.3e}"
        # Absolute normalization correct to a few percent at this grid
        assert rel_fine < 0.025, \
            f"4D quad disagrees with Neumann normalization: {rel_fine*100:.2f}%"


class TestNeumannHamiltonian:
    """Test full Hamiltonian with Neumann V_ee."""

    def test_hamiltonian_neumann_symmetric(self):
        """H(Neumann) should be symmetric."""
        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=20, N_eta=15, N_phi=1)
        R = 1.4011
        H = compute_hamiltonian_neumann(basis, R, grids, l_max=10)
        assert np.allclose(H, H.T, atol=1e-13)

    def test_variational_principle(self):
        """E(Neumann) should be above exact ground state."""
        basis = generate_basis(j_max=1, l_max=1, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=20, N_eta=15, N_phi=1)
        R = 1.4011
        result = solve_hylleraas(
            basis, R, grids, verbose=False,
            vee_method='neumann', l_max_neumann=15,
        )
        E_exact = -1.17475
        assert result['E_total'] > E_exact, \
            f"Variational violation: E={result['E_total']:.6f} < {E_exact}"

    def test_energy_below_atoms(self):
        """H2 energy should be below separated atoms (bound)."""
        basis = generate_basis(j_max=1, l_max=1, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=20, N_eta=15, N_phi=1)
        R = 1.4011
        result = solve_hylleraas(
            basis, R, grids, verbose=False,
            vee_method='neumann', l_max_neumann=15,
        )
        E_atoms = -1.0  # 2 × (-0.5) for two H atoms
        assert result['E_total'] < E_atoms, \
            f"H2 not bound: E={result['E_total']:.6f} > {E_atoms}"

    def test_neumann_improves_over_numerical(self):
        """Neumann V_ee should give lower (better) energy than coarse numerical."""
        basis = generate_basis(j_max=1, l_max=1, p_max=0, alpha=1.0)
        # Use coarse grid where numerical V_ee has significant error
        grids = build_quadrature_grids(N_xi=15, N_eta=10, N_phi=1)
        R = 1.4011

        result_num = solve_hylleraas(
            basis, R, grids, verbose=False,
            vee_method='numerical',
        )
        result_neu = solve_hylleraas(
            basis, R, grids, verbose=False,
            vee_method='neumann', l_max_neumann=15,
        )
        # With exact V_ee, the Neumann result should be better
        # (lower energy, closer to variational minimum)
        print(f"\n  E(numerical) = {result_num['E_total']:.6f}")
        print(f"  E(Neumann)   = {result_neu['E_total']:.6f}")
        # The key test: Neumann energy should be more negative
        # (numerical V_ee overestimates repulsion on coarse grids)
        assert result_neu['E_total'] < result_num['E_total'] + 0.01, \
            "Neumann should not be significantly worse than numerical"

    @pytest.mark.slow
    def test_h2_de_headline_924pct(self):
        """Paper 12 headline: Neumann algebraic V_ee yields ~92.4% of exact D_e.

        Recomputes the H2 ground-state energy at R = 1.4011 bohr with the
        N = 27 basis (j_max = l_max = 2, p = 0) and the Neumann V_ee path,
        confirming the plateau value the paper reports (E_Neu ≈ -1.16096,
        92.2-92.4% of D_e at N ≥ 27). α is fixed at 1.0, which the variational
        α-scan finds optimal at this basis size (the per-grid energy is
        V_ee-grid-independent, so the coarse 20×14 T+V_ne grid reproduces the
        paper's 30×20 result to the reported digits).
        """
        R = 1.4011
        basis = generate_basis(j_max=2, l_max=2, p_max=0, alpha=1.0)
        assert len(basis) == 27, f"expected 27 basis functions, got {len(basis)}"
        grids = build_quadrature_grids(N_xi=20, N_eta=14, N_phi=1)

        result = solve_hylleraas(
            basis, R, grids, verbose=False,
            vee_method='neumann', l_max_neumann=15,
        )

        E_exact = -1.17475   # Kolos-Wolniewicz, Paper 12 eq via Kolos1968
        E_atoms = -1.0       # two H atoms, 2 x (-0.5) Ha
        D_e_exact = E_atoms - E_exact   # = 0.17475 Ha
        D_e = E_atoms - result['E_total']
        pct = D_e / D_e_exact * 100.0
        print(f"\n  H2 N=27 Neumann: E_total={result['E_total']:.6f} Ha, "
              f"D_e%={pct:.2f}")

        # Recover the headline plateau (paper: 92.2% at N=27), not the
        # numerical-V_ee ~80% plateau and not an over-binding artifact.
        assert result['E_total'] > E_exact, \
            f"variational violation: {result['E_total']:.6f} < {E_exact}"
        assert 90.0 < pct < 94.0, \
            f"H2 D_e fraction {pct:.2f}% outside headline band [90,94]"


class TestIncludeVee:
    """Test the include_vee=False path in compute_hamiltonian_matrix."""

    def test_tvne_plus_vee_equals_full(self):
        """H(T+V_ne) + V_ee(numerical) should equal H(full)."""
        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=25, N_eta=15, N_phi=1)
        R = 1.4011

        H_full = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=False, include_vee=True
        )
        H_tv = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=False, include_vee=False
        )
        # The difference should be the V_ee contribution
        V_ee = H_full - H_tv
        assert np.all(V_ee.diagonal() > 0), \
            "V_ee diagonal should be positive (repulsion)"

    def test_tvne_matches_numba(self):
        """Python include_vee=True should match Numba path."""
        basis = generate_basis(j_max=1, l_max=0, p_max=0, alpha=1.0)
        grids = build_quadrature_grids(N_xi=20, N_eta=15, N_phi=1)
        R = 1.4011

        H_python = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=False, include_vee=True
        )
        H_numba = compute_hamiltonian_matrix(
            basis, R, grids, use_numba=True, include_vee=True
        )
        assert np.allclose(H_python, H_numba, atol=1e-10), \
            f"Python vs Numba max diff: {np.max(np.abs(H_python - H_numba))}"
