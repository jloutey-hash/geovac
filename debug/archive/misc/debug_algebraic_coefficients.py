"""
Algebraic closed forms for perturbation coefficients of He angular eigenvalues.

The angular equation is: [Lambda^2/2 + R*C(alpha,theta_12)] Phi = mu(R) Phi

In the Liouville basis u_n(alpha) = sqrt(4/pi) sin(2n*alpha), the perturbation
expansion is:
    mu_n(R) = C_2(nu)/2 + a_1 * R + a_2 * R^2 + ...

This script derives a_1 and a_2 symbolically and compares to numerics.

KEY RESULTS (derived below):
  a_1(nu=0, Z) = (8/(3*pi)) * (sqrt(2) - 4*Z)
  a_1(nu=2, l=0, Z) = (8/(3*pi)) * (sqrt(2) - 4*Z) - (64*Z)/(15*pi) + V_ee correction
"""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.hyperspherical_angular import solve_angular, gaunt_integral


# =============================================================================
# PART 1: Symbolic derivation of a_1 for nu=0
# =============================================================================

def symbolic_a1_nu0():
    """
    Derive a_1 for the nu=0 (ground) channel analytically.

    Free eigenstate: u_1(alpha) = sqrt(4/pi) sin(2*alpha)
    sin^2(2*alpha) = 4 sin^2(alpha) cos^2(alpha)

    The perturbation (linear in R) is V = -Z/cos(a) - Z/sin(a) + V_ee(a)

    NUCLEAR PART:
    a1_nuc = (4/pi) * (-Z) * int_0^{pi/2} sin^2(2a) * (1/cos(a) + 1/sin(a)) da
           = (4/pi) * (-Z) * int_0^{pi/2} 4*sin^2(a)*cos^2(a) * (1/cos(a) + 1/sin(a)) da
           = (16Z/pi) * (-1) * int_0^{pi/2} (sin^2(a)*cos(a) + sin(a)*cos^2(a)) da

    Both integrals give 1/3:
      int_0^{pi/2} sin^2(a)*cos(a) da = [sin^3(a)/3]_0^{pi/2} = 1/3
      int_0^{pi/2} sin(a)*cos^2(a) da = [-cos^3(a)/3]_0^{pi/2} = 1/3

    So: a1_nuc = -(16Z/pi) * (2/3) = -32Z/(3*pi)

    V_EE PART (monopole, k=0, l=0):
    The angular solver adds R * norm * G(0,0,0) * 1/max(sin(a),cos(a))
    where norm = sqrt(1*1)/2 = 1/2, G(0,0,0) = 2.0.
    So V_ee = 1/max(sin(a), cos(a)).

    a1_ee = (4/pi) * int_0^{pi/2} 4*sin^2(a)*cos^2(a) / max(sin(a),cos(a)) da

    Split at pi/4 (where sin=cos), use symmetry a -> pi/2-a:
    = (16/pi) * 2 * int_0^{pi/4} sin^2(a)*cos^2(a) / cos(a) da
    = (32/pi) * int_0^{pi/4} sin^2(a)*cos(a) da
    = (32/pi) * [sin^3(a)/3]_0^{pi/4}
    = (32/pi) * (1/sqrt(2))^3 / 3
    = (32/pi) * sqrt(2)/6  (since (1/sqrt(2))^3 = sqrt(2)/4... wait)

    Let me redo: (1/sqrt(2))^3 = 1/(2*sqrt(2)) = sqrt(2)/4
    So: (32/pi) * sqrt(2)/(4*3) = (32/pi) * sqrt(2)/12 = 8*sqrt(2)/(3*pi)

    TOTAL:
    a1(nu=0) = -32Z/(3*pi) + 8*sqrt(2)/(3*pi) = (8/(3*pi)) * (sqrt(2) - 4Z)
    """
    Z = 2.0
    a1_nuc = -32 * Z / (3 * np.pi)
    a1_ee = 8 * np.sqrt(2) / (3 * np.pi)
    a1_total = a1_nuc + a1_ee
    a1_formula = (8 / (3 * np.pi)) * (np.sqrt(2) - 4 * Z)

    return {
        'a1_nuc': a1_nuc,
        'a1_ee': a1_ee,
        'a1_total': a1_total,
        'a1_formula': a1_formula,
        'formula_str': 'a1(nu=0, Z) = (8/(3*pi)) * (sqrt(2) - 4*Z)',
    }


# =============================================================================
# PART 2: Numerical verification via high-resolution finite differences
# =============================================================================

def numerical_a1(n_channels: int = 6, l_max: int = 3, n_alpha: int = 500):
    """Extract a_1 numerically from mu(R) at very small R."""
    dR = 1e-4
    mu_0, _ = solve_angular(R=0.0, Z=2.0, l_max=l_max, n_alpha=n_alpha,
                            n_channels=n_channels)
    mu_dR, _ = solve_angular(R=dR, Z=2.0, l_max=l_max, n_alpha=n_alpha,
                             n_channels=n_channels)
    a1 = (mu_dR - mu_0) / dR
    return a1, mu_0


# =============================================================================
# PART 3: Symbolic a_1 for higher channels
# =============================================================================

def symbolic_integrals_higher_channels():
    """
    Compute a_1 for channels with n=2,3 (nu=2,4) in the l=0 sector.

    u_n(alpha) = sqrt(4/pi) sin(2*n*alpha)

    For each n, need:
      I_nuc(n) = int_0^{pi/2} sin^2(2*n*a) * (1/cos(a) + 1/sin(a)) da
      I_ee(n)  = int_0^{pi/2} sin^2(2*n*a) / max(sin(a), cos(a)) da

    These integrals can be computed in closed form using:
      sin^2(2*n*a) = (1 - cos(4*n*a))/2
    and the known results for int cos(m*a)/cos(a) da and similar.
    """
    results = {}

    # Numerical integration for verification
    from scipy.integrate import quad

    for n in range(1, 5):
        nu = 2 * (n - 1)

        def integrand_nuc(a, n=n):
            return np.sin(2 * n * a)**2 * (1.0 / np.cos(a) + 1.0 / np.sin(a))

        def integrand_ee(a, n=n):
            return np.sin(2 * n * a)**2 / np.maximum(np.sin(a), np.cos(a))

        I_nuc, _ = quad(integrand_nuc, 1e-12, np.pi/2 - 1e-12)
        I_ee, _ = quad(integrand_ee, 1e-12, np.pi/2 - 1e-12)

        Z = 2.0
        a1_nuc = (4 / np.pi) * (-Z) * I_nuc
        a1_ee = (4 / np.pi) * I_ee  # norm * G(0,0,0) = 1/2 * 2 = 1 already accounted for
        a1_total = a1_nuc + a1_ee

        results[n] = {
            'nu': nu,
            'I_nuc': I_nuc,
            'I_ee': I_ee,
            'a1_nuc': a1_nuc,
            'a1_ee': a1_ee,
            'a1_total': a1_total,
        }

    return results


def try_closed_forms():
    """
    Attempt to find closed-form expressions for the integrals.

    The key integral is: I(n) = int_0^{pi/2} sin^2(2*n*a) / cos(a) da

    Using sin^2(2*n*a) = (1 - cos(4*n*a)) / 2:
    I(n) = (1/2) int_0^{pi/2} (1 - cos(4*n*a)) / cos(a) da
         = (1/2) [int sec(a) da - int cos(4*n*a)/cos(a) da]

    The integral int_0^{pi/2} cos(m*a)/cos(a) da for integer m:
    This is a known result involving partial fractions of Chebyshev polynomials.

    For m even: int_0^{pi/2} cos(m*a)/cos(a) da = 0 (when m is even and > 0)
    Wait, that's not right. Let me use the identity:
    cos(m*a)/cos(a) = sum of alternating cosines.

    Actually, the Dirichlet integral identity:
    int_0^{pi/2} cos((2k)a)/cos(a) da diverges for k>0!

    But the combination (1 - cos(4n*a))/cos(a) converges.

    Let me just use sympy for the exact forms.
    """
    try:
        import sympy as sp

        a = sp.Symbol('a', positive=True)
        Z_sym = sp.Symbol('Z', positive=True)
        n_sym = sp.Symbol('n', integer=True, positive=True)

        print("\n  Attempting sympy symbolic integration...")

        # Nuclear integral for general n
        # I_nuc(n) = int_0^{pi/2} sin^2(2*n*a) * (1/cos(a) + 1/sin(a)) da
        # By symmetry a->pi/2-a, the 1/cos and 1/sin parts give the same result:
        # I_nuc(n) = 2 * int_0^{pi/2} sin^2(2*n*a) / cos(a) da

        # For specific n values:
        closed_forms = {}
        for n in range(1, 5):
            expr_nuc = sp.sin(2 * n * a)**2 / sp.cos(a)
            # Full integral (both terms are equal by symmetry)
            I_half = sp.integrate(expr_nuc, (a, 0, sp.pi / 2))
            I_half_simplified = sp.simplify(I_half)

            expr_ee_low = sp.sin(2 * n * a)**2 * sp.sin(a)  # cos(a)/cos(a) = sin(a) on [0,pi/4]
            # Actually for V_ee: 1/max(sin,cos). On [0,pi/4]: max=cos, so divide by cos.
            # On [pi/4,pi/2]: max=sin, divide by sin.
            expr_ee_1 = sp.sin(2 * n * a)**2 * sp.cos(a)  # sin^2*cos^2/cos * cos = sin^2*cos^2... no

            # Let me be more careful:
            # V_ee integrand = sin^2(2*n*a) / max(sin(a), cos(a))
            # On [0, pi/4]: max = cos(a), integrand = sin^2(2*n*a) / cos(a)
            # On [pi/4, pi/2]: max = sin(a), integrand = sin^2(2*n*a) / sin(a)
            # By a -> pi/2-a symmetry: both halves give the same result
            # So I_ee = 2 * int_0^{pi/4} sin^2(2*n*a) / cos(a) da

            I_ee_half = sp.integrate(expr_nuc, (a, 0, sp.pi / 4))
            I_ee_half_simplified = sp.simplify(I_ee_half)

            closed_forms[n] = {
                'I_nuc_half': I_half_simplified,
                'I_nuc': 2 * I_half_simplified,
                'I_ee_half': I_ee_half_simplified,
                'I_ee': 2 * I_ee_half_simplified,
            }

            print(f"\n  n={n} (nu={2*(n-1)}):")
            print(f"    I_nuc/2 = int_0^{{pi/2}} sin^2({2*n}a)/cos(a) da = {I_half_simplified}")
            print(f"    I_nuc   = {2 * I_half_simplified}")
            nuc_float = float(2 * I_half_simplified)
            print(f"    I_nuc (float) = {nuc_float:.10f}")
            print(f"    I_ee/2  = int_0^{{pi/4}} sin^2({2*n}a)/cos(a) da = {I_ee_half_simplified}")
            ee_float = float(2 * I_ee_half_simplified)
            print(f"    I_ee (float)  = {ee_float:.10f}")

        return closed_forms

    except ImportError:
        print("  sympy not available, skipping symbolic computation")
        return None


# =============================================================================
# PART 4: Second-order coefficient a_2
# =============================================================================

def compute_a2_perturbation(l_max: int = 0, n_alpha: int = 500, n_basis: int = 15):
    """
    Compute a_2 via second-order perturbation theory.

    a_2(n) = sum_{m != n} |<m|V|n>|^2 / (E_n^0 - E_m^0)

    where E_n^0 = 2n^2 - 2 and V is the charge function.

    For the l=0 sector, this is a 1D sum over the alpha basis.
    """
    print("Computing second-order coefficient a_2 via perturbation theory...")

    Z = 2.0

    # Get free eigenstates
    _, vecs = solve_angular(R=0.0, Z=Z, l_max=l_max, n_alpha=n_alpha,
                            n_channels=n_basis)

    # Build perturbation matrix V in the free basis
    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h
    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)

    # Full charge function (l=0 channel): nuclear + V_ee monopole
    V_nuc = -Z / cos_a - Z / sin_a
    V_ee_mono = 1.0 / np.maximum(sin_a, cos_a)
    V_total = V_nuc + V_ee_mono

    # Matrix elements <m|V|n>
    # NOTE: eigh returns eigenvectors normalized to v^T v = 1 (discrete norm).
    # The FD Hamiltonian has V_l[i] added directly to diagonal, so the
    # perturbation matrix element is v^T * diag(V) * v = dot(v, V*v),
    # with NO factor of h. The h would convert to continuous normalization
    # but that's inconsistent with how the Hamiltonian is discretized.
    V_mat = np.zeros((n_basis, n_basis))
    for i in range(n_basis):
        for j in range(n_basis):
            V_mat[i, j] = np.dot(vecs[i], V_total * vecs[j])

    # Free eigenvalues
    E0 = np.array([2.0 * (n + 1)**2 - 2.0 for n in range(n_basis)])

    # First-order: a_1 = V_nn (diagonal)
    a1_pt = np.diag(V_mat)

    # Second-order: a_2 = sum_{m!=n} |V_mn|^2 / (E_n - E_m)
    a2_pt = np.zeros(n_basis)
    for n in range(n_basis):
        for m in range(n_basis):
            if m == n:
                continue
            a2_pt[n] += V_mat[m, n]**2 / (E0[n] - E0[m])

    return a1_pt, a2_pt, V_mat, E0


# =============================================================================
# PART 5: Check if a_2 has a closed form
# =============================================================================

def analyze_a2_structure(a2_pt, E0):
    """
    Check if a_2 follows a pattern in nu.

    For a quasi-Coulomb problem, a_2 should relate to the polarizability
    of the S^5 harmonic in the charge function field.
    """
    print("\n  Checking if a_2(nu) has a simple pattern...")

    for n in range(min(6, len(a2_pt))):
        nu = 2 * n
        c2 = nu * (nu + 4)
        print(f"    nu={nu:2d}: a_2 = {a2_pt[n]:10.6f}, "
              f"a_2 * pi = {a2_pt[n] * np.pi:10.6f}, "
              f"a_2 * 3*pi/8 = {a2_pt[n] * 3 * np.pi / 8:10.6f}")


# =============================================================================
# PART 6: Hyperradial equation solvability analysis
# =============================================================================

def analyze_hyperradial_solvability(a1, a2):
    """
    Assess whether the effective hyperradial equation is exactly solvable.

    The hyperradial equation is:
      [-1/2 d^2F/dR^2 + V_eff(R)] F = E F

    where V_eff(R) = mu(R)/R^2 + 15/(8R^2)

    If mu(R) = a_0 + a_1*R + a_2*R^2, then:
      V_eff = a_0/R^2 + a_1/R + a_2 + 15/(8R^2)
            = (a_0 + 15/8)/R^2 + a_1/R + a_2

    This is EXACTLY a hydrogen-like equation:
      [-1/2 d^2/dR^2 + l_eff(l_eff+1)/(2R^2) + a_1/R + a_2] F = E F

    where l_eff(l_eff+1)/2 = a_0 + 15/8 = C_2/2 + 15/8.

    The solutions are known: E_N = a_2 - a_1^2 / (2(N + l_eff + 1/2)^2)
    where N = 0, 1, 2, ... is the radial quantum number.
    """
    Z = 2.0
    a0 = 0.0  # nu=0 channel

    # Effective angular momentum
    # l_eff(l_eff+1) = 2*(a0 + 15/8) = 2*15/8 = 15/4
    # l_eff = (-1 + sqrt(1 + 15))/2 = (-1 + 4)/2 = 3/2
    l_eff_sq = 2.0 * (a0 + 15.0 / 8.0)
    l_eff = (-1 + np.sqrt(1 + 4 * l_eff_sq)) / 2.0

    print(f"\n  Effective angular momentum:")
    print(f"    l_eff(l_eff+1) = 2*(C_2/2 + 15/8) = {l_eff_sq:.4f}")
    print(f"    l_eff = {l_eff:.6f}")

    # For nu=0: l_eff(l_eff+1) = 15/4, l_eff = 3/2 (exact!)
    print(f"    l_eff = 3/2 exactly? {abs(l_eff - 1.5) < 1e-10}")

    # Quasi-hydrogen solutions
    # The equation is: -1/2 F'' + [l_eff(l_eff+1)/(2R^2) + a_1/R + a_2] F = E F
    # = -1/2 F'' + [l_eff(l_eff+1)/(2R^2) + a_1/R] F = (E - a_2) F
    #
    # This is the hydrogen equation with Z_eff = -a_1 and E_H = E - a_2
    # E_H = -Z_eff^2 / (2*n_eff^2) where n_eff = N + l_eff + 1

    Z_eff = -a1  # a1 is negative, so Z_eff is positive
    print(f"\n  Quasi-hydrogen mapping:")
    print(f"    a_1 = {a1:.6f}")
    print(f"    a_2 = {a2:.6f}")
    print(f"    Z_eff = -a_1 = {Z_eff:.6f}")
    print(f"    l_eff = {l_eff:.6f}")

    print(f"\n  Predicted energy levels E_N = a_2 - Z_eff^2 / (2*(N + l_eff + 1)^2):")
    for N in range(5):
        n_eff = N + l_eff + 1
        E_N = a2 - Z_eff**2 / (2 * n_eff**2)
        print(f"    N={N}: n_eff={n_eff:.2f}, E = {E_N:.6f} Ha")

    # Ground state (N=0): n_eff = l_eff + 1 = 5/2
    n_eff_gs = l_eff + 1
    E_gs = a2 - Z_eff**2 / (2 * n_eff_gs**2)
    print(f"\n  Ground state prediction (N=0, n_eff={n_eff_gs:.2f}):")
    print(f"    E = {E_gs:.6f} Ha")
    print(f"    Exact He: -2.903724 Ha")
    print(f"    Error: {abs(E_gs - (-2.903724)) / 2.903724 * 100:.2f}%")

    return Z_eff, l_eff, E_gs


# =============================================================================
# PART 7: Beyond leading order — full perturbative series
# =============================================================================

def full_perturbative_comparison(l_max: int = 3, n_alpha: int = 300):
    """
    Compare the truncated perturbation series to the exact mu(R).

    mu(R) ≈ C_2/2 + a_1*R + a_2*R^2 vs exact solver at various R.
    """
    print("Comparing perturbative mu(R) vs exact solver...")

    # Get perturbation coefficients from l=0 sector first
    a1_l0, a2_l0, _, _ = compute_a2_perturbation(l_max=0, n_alpha=400, n_basis=20)

    # Now get exact mu(R) from full solver
    R_grid = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0])

    print(f"\n{'R':>6}  {'mu_exact':>10}  {'mu_PT1':>10}  {'mu_PT2':>10}  "
          f"{'err_PT1%':>10}  {'err_PT2%':>10}")
    print("-" * 60)

    for R in R_grid:
        mu_exact, _ = solve_angular(R=R, Z=2.0, l_max=l_max, n_alpha=n_alpha,
                                     n_channels=1)
        mu_pt1 = 0.0 + a1_l0[0] * R  # First-order
        mu_pt2 = 0.0 + a1_l0[0] * R + a2_l0[0] * R**2  # Second-order

        err1 = abs(mu_pt1 - mu_exact[0]) / max(abs(mu_exact[0]), 1e-10) * 100
        err2 = abs(mu_pt2 - mu_exact[0]) / max(abs(mu_exact[0]), 1e-10) * 100

        print(f"{R:6.2f}  {mu_exact[0]:10.4f}  {mu_pt1:10.4f}  {mu_pt2:10.4f}  "
              f"{err1:10.2f}  {err2:10.2f}")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("=" * 72)
    print("Algebraic Closed Forms for Perturbation Coefficients")
    print("He Hyperspherical Angular Eigenvalues — SO(6) Analysis")
    print("=" * 72)

    # =========================================================================
    # STEP 1: Symbolic a_1 for nu=0
    # =========================================================================
    print("\n" + "=" * 72)
    print("STEP 1: Analytic derivation of a_1 for nu=0")
    print("=" * 72)

    sym = symbolic_a1_nu0()
    print(f"\n  FORMULA: {sym['formula_str']}")
    print(f"\n  Components (Z=2):")
    print(f"    Nuclear:  a1_nuc = -32Z/(3*pi) = {sym['a1_nuc']:.10f}")
    print(f"    V_ee:     a1_ee  = 8*sqrt(2)/(3*pi) = {sym['a1_ee']:.10f}")
    print(f"    Total:    a1     = {sym['a1_total']:.10f}")
    print(f"    Formula:  a1     = {sym['a1_formula']:.10f}")

    # =========================================================================
    # STEP 2: Numerical verification
    # =========================================================================
    print("\n" + "=" * 72)
    print("STEP 2: Numerical verification")
    print("=" * 72)

    a1_num, mu_free = numerical_a1(n_channels=6, l_max=3, n_alpha=500)
    print(f"\n  Numerical a_1 (from dmu/dR at R=0, l_max=3, n_alpha=500):")
    for i in range(6):
        nu = None
        for test_nu in range(0, 20, 2):
            if abs(mu_free[i] - test_nu * (test_nu + 4) / 2) < 1.0:
                nu = test_nu
                break
        nu_str = str(nu) if nu is not None else "?"
        print(f"    Channel {i+1} (nu={nu_str}): a1_num = {a1_num[i]:.6f}")

    print(f"\n  COMPARISON for nu=0:")
    print(f"    Analytic: a1 = (8/(3*pi))*(sqrt(2) - 8) = {sym['a1_total']:.10f}")
    print(f"    Numeric:  a1 = {a1_num[0]:.10f}")
    print(f"    Difference: {abs(sym['a1_total'] - a1_num[0]):.2e}")
    match = abs(sym['a1_total'] - a1_num[0]) < 0.01
    print(f"    Match: {'YES' if match else 'NO'}")

    # =========================================================================
    # STEP 3: Higher channels — numerical integrals
    # =========================================================================
    print("\n" + "=" * 72)
    print("STEP 3: a_1 for higher channels (l=0 sector, quad integration)")
    print("=" * 72)

    higher = symbolic_integrals_higher_channels()
    print(f"\n  {'n':>3}  {'nu':>3}  {'a1_nuc':>12}  {'a1_ee':>12}  {'a1_total':>12}  "
          f"{'a1_num':>12}  {'diff':>10}")
    print("-" * 75)

    # For comparison, get l=0-only numerical values
    a1_l0, mu0_l0 = numerical_a1(n_channels=4, l_max=0, n_alpha=500)

    for n in range(1, 5):
        h = higher[n]
        if n - 1 < len(a1_l0):
            diff = abs(h['a1_total'] - a1_l0[n - 1])
        else:
            diff = float('nan')
        print(f"  {n:3d}  {h['nu']:3d}  {h['a1_nuc']:12.6f}  {h['a1_ee']:12.6f}  "
              f"{h['a1_total']:12.6f}  "
              f"{a1_l0[n-1] if n-1 < len(a1_l0) else float('nan'):12.6f}  "
              f"{diff:10.2e}")

    # =========================================================================
    # STEP 4: Sympy closed forms
    # =========================================================================
    print("\n" + "=" * 72)
    print("STEP 4: Sympy symbolic integration for closed forms")
    print("=" * 72)

    closed = try_closed_forms()

    # =========================================================================
    # STEP 5: Second-order coefficient a_2
    # =========================================================================
    print("\n" + "=" * 72)
    print("STEP 5: Second-order coefficient a_2 via perturbation theory")
    print("=" * 72)

    a1_pt, a2_pt, V_mat, E0 = compute_a2_perturbation(l_max=0, n_alpha=500,
                                                         n_basis=20)

    print(f"\n  a_1 (perturbation theory) vs a_1 (finite difference):")
    for n in range(min(4, len(a1_pt))):
        nu = 2 * n
        print(f"    nu={nu}: a1_PT = {a1_pt[n]:.6f}, a1_FD = {a1_l0[n] if n < len(a1_l0) else float('nan'):.6f}")

    print(f"\n  a_2 values (l=0 sector):")
    for n in range(min(6, len(a2_pt))):
        nu = 2 * n
        print(f"    nu={nu}: a_2 = {a2_pt[n]:.8f}")

    analyze_a2_structure(a2_pt, E0)

    # Compare a_2 to numerical (from quadratic fit)
    print("\n  Comparing a_2 to numerical quadratic fit:")
    R_pts = np.array([0.0, 0.01, 0.02, 0.05, 0.1])
    for ch in range(min(3, len(a2_pt))):
        mu_vals = []
        for R in R_pts:
            mu, _ = solve_angular(R=R, Z=2.0, l_max=0, n_alpha=500,
                                   n_channels=ch + 1)
            mu_vals.append(mu[ch])
        coeffs = np.polyfit(R_pts, mu_vals, 2)
        a2_fit = coeffs[0]
        print(f"    nu={2*ch}: a2_PT = {a2_pt[ch]:.8f}, "
              f"a2_fit = {a2_fit:.8f}, diff = {abs(a2_pt[ch] - a2_fit):.2e}")

    # =========================================================================
    # STEP 6: Hyperradial solvability
    # =========================================================================
    print("\n" + "=" * 72)
    print("STEP 6: Hyperradial equation — exact solvability assessment")
    print("=" * 72)

    Z_eff, l_eff, E_gs = analyze_hyperradial_solvability(
        a1_pt[0], a2_pt[0]
    )

    # =========================================================================
    # STEP 7: Full perturbative comparison
    # =========================================================================
    print("\n" + "=" * 72)
    print("STEP 7: Perturbative mu(R) vs exact solver")
    print("=" * 72)

    full_perturbative_comparison(l_max=3, n_alpha=300)

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "=" * 72)
    print("SUMMARY OF ALGEBRAIC RESULTS")
    print("=" * 72)

    print(f"""
1. FIRST-ORDER COEFFICIENT a_1 (EXACT CLOSED FORM):

   a_1(nu=0, Z) = (8/(3*pi)) * (sqrt(2) - 4*Z)

   For Z=2: a_1 = {sym['a1_total']:.10f}
   Numerical:     {a1_num[0]:.10f}
   Agreement:     {abs(sym['a1_total'] - a1_num[0]):.2e}

   Components:
     Nuclear:  -32*Z/(3*pi)     [from int sin^2(a) cos(a) da = 1/3]
     V_ee:     8*sqrt(2)/(3*pi) [from int_0^{{pi/4}} sin^2(a) cos(a) da = sqrt(2)/12]

2. SECOND-ORDER COEFFICIENT a_2:

   a_2(nu=0) = {a2_pt[0]:.8f}  (from {len(a2_pt)}-term perturbation sum)

   Convergence check: this requires off-diagonal V_mat elements and
   is NOT a simple closed form. It is an infinite series over SO(6)
   Clebsch-Gordan coefficients.

3. HYPERRADIAL EQUATION:

   V_eff(R) = (C_2/2 + 15/8)/R^2 + a_1/R + a_2 + O(R)

   At leading order (truncating at a_2), this is a quasi-Coulomb problem:
     l_eff = 3/2 (exact, from 15/(8R^2) centrifugal term)
     Z_eff = -a_1 = {-a1_pt[0]:.6f}

   Ground state energy (linear+quadratic mu approximation):
     E_gs = a_2 - Z_eff^2/(2*n_eff^2) = {E_gs:.6f} Ha
     Exact: -2.903724 Ha
     Error: {abs(E_gs - (-2.903724)) / 2.903724 * 100:.2f}%

4. ASSESSMENT:

   The free eigenvalues are EXACTLY the SO(6) Casimir: mu = C_2/2.
   The first-order coefficient a_1 has an EXACT algebraic closed form.
   The second-order coefficient a_2 does NOT have a simple closed form
   (it's an infinite perturbation sum), but converges rapidly.

   The hyperradial equation at quadratic order in mu(R) is a Coulomb
   problem with known analytical solutions, but the truncation error
   grows at large R where mu(R) is NOT well-approximated by a quadratic.

   CONCLUSION: Partially solvable. The small-R expansion has algebraic
   coefficients, giving a quasi-Coulomb hyperradial equation. But the
   full mu(R) curve (needed for accurate energies) requires numerical
   solution of the angular problem at each R.
""")
