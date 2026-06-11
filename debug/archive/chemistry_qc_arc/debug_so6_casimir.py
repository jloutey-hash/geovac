"""
SO(6) Casimir analysis of He hyperspherical angular eigenvalues.

Investigates whether angular eigenvalues mu(R) have a closed-form
expression in terms of SO(6) Casimir operators.

Key insight: The hyperangular coordinates (alpha, theta_1, phi_1, theta_2, phi_2)
parametrize S^5. The isometry group is SO(6). The free (R=0) angular eigenvalues
are EXACTLY Lambda^2/2 = nu(nu+4)/2, i.e., half the SO(6) quadratic Casimir.

At finite R, the charge function C(alpha, theta_12) breaks SO(6) symmetry,
mixing different nu sectors. This analysis quantifies the deviation.
"""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.hyperspherical_angular import solve_angular


def so6_casimir_c2(nu: int) -> float:
    """Quadratic Casimir of SO(6) for symmetric representation (nu, 0, 0)."""
    return float(nu * (nu + 4))


def free_eigenvalue(n: int) -> float:
    """Free (R=0) angular eigenvalue: mu = 2n^2 - 2 = C2(2(n-1)) / 2."""
    return 2.0 * n**2 - 2.0


def analyze_free_spectrum(l_max: int = 5, n_alpha: int = 400, n_channels: int = 12):
    """
    Step 1: Verify free (R=0) eigenvalues match SO(6) Casimir formula.

    At R=0, the angular equation is just Lambda^2/2 on S^5 (restricted to L=0).
    Eigenvalues should be nu(nu+4)/2 with nu = 0, 2, 4, ...
    """
    print("=" * 72)
    print("STEP 1: Free spectrum (R=0) vs SO(6) Casimir")
    print("=" * 72)

    mu, _ = solve_angular(R=0.0, Z=2.0, l_max=l_max, n_alpha=n_alpha,
                          n_channels=n_channels)

    print(f"\n{'n':>3}  {'nu':>3}  {'C2=nu(nu+4)':>12}  {'C2/2':>8}  "
          f"{'mu_num':>10}  {'delta':>10}  {'match?':>7}")
    print("-" * 72)

    # The free eigenvalues for L=0, 1S sector
    # With l_max channels, we get states with nu = 0, 2, 4, ...
    # But also higher-l states can contribute. Need to match carefully.

    # Expected free eigenvalues: 2n^2-2 for l=0 channel, plus
    # contributions from l>0 channels at higher energies
    # For L=0, the angular quantum numbers are (nu, l) with nu >= 2l
    # and nu even. The eigenvalue is nu(nu+4)/2.
    # Degeneracy: for each nu, l ranges from 0 to nu/2.

    # Build expected spectrum for L=0, 1S
    expected = []
    for nu in range(0, 2 * (l_max + 1) + 10, 2):  # even nu only for 1S
        for l in range(min(nu // 2 + 1, l_max + 1)):
            # For L=0, l1=l2=l, and the radial quantum number n_alpha
            # gives additional states within each (nu, l) sector
            # Actually, for fixed l, the alpha quantum number gives
            # n_r = 0, 1, 2, ... with nu = 2*n_r + 2*l
            if nu >= 2 * l:
                n_r = (nu - 2 * l) // 2
                expected.append({
                    'nu': nu,
                    'l': l,
                    'n_r': n_r,
                    'C2': so6_casimir_c2(nu),
                    'mu_exact': nu * (nu + 4) / 2.0,
                })

    expected.sort(key=lambda x: x['mu_exact'])

    # Match numerical to expected
    used = set()
    for i, mu_i in enumerate(mu):
        # Find closest expected
        best_j = None
        best_dist = float('inf')
        for j, e in enumerate(expected):
            if j in used:
                continue
            dist = abs(mu_i - e['mu_exact'])
            if dist < best_dist:
                best_dist = dist
                best_j = j

        if best_j is not None and best_dist < 2.0:
            e = expected[best_j]
            used.add(best_j)
            match = "YES" if best_dist < 0.5 else "~" if best_dist < 2.0 else "NO"
            print(f"{i+1:3d}  {e['nu']:3d}  {e['C2']:12.1f}  {e['C2']/2:8.1f}  "
                  f"{mu_i:10.4f}  {mu_i - e['mu_exact']:10.4f}  {match:>7}")
        else:
            print(f"{i+1:3d}  {'?':>3s}  {'':>12s}  {'':>8s}  "
                  f"{mu_i:10.4f}  {'':>10s}  {'???':>7}")

    return mu


def analyze_finite_R(R_values=None, l_max: int = 3, n_alpha: int = 200,
                     n_channels: int = 6):
    """
    Step 2: Track eigenvalues at finite R — how do they deviate from C2/2?

    At finite R, the charge function C(alpha, theta_12) mixes SO(6) sectors.
    We track mu(R) and compare to the free values.
    """
    print("\n" + "=" * 72)
    print("STEP 2: Angular eigenvalues at finite R")
    print("=" * 72)

    if R_values is None:
        R_values = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]

    results = {}
    for R in R_values:
        mu, _ = solve_angular(R=R, Z=2.0, l_max=l_max, n_alpha=n_alpha,
                              n_channels=n_channels)
        results[R] = mu

    # Print table
    header = f"{'Channel':>8}"
    for R in R_values:
        header += f"  {'R='+str(R):>10}"
    print(f"\n{header}")
    print("-" * len(header))

    for ch in range(n_channels):
        row = f"{ch+1:8d}"
        for R in R_values:
            row += f"  {results[R][ch]:10.4f}"
        print(row)

    # Compute deviation from free values
    print(f"\n{'Channel':>8}  {'mu_free':>10}  {'nu':>4}  ", end="")
    for R in R_values[1:]:
        print(f"  {'dmu(R='+str(R)+')':>12}", end="")
    print()
    print("-" * 100)

    free_mu = results[0.0]
    for ch in range(n_channels):
        nu = None
        # Match to nearest C2/2
        for test_nu in range(0, 40, 2):
            if abs(free_mu[ch] - test_nu * (test_nu + 4) / 2) < 1.0:
                nu = test_nu
                break
        nu_str = str(nu) if nu is not None else "?"
        print(f"{ch+1:8d}  {free_mu[ch]:10.4f}  {nu_str:>4s}  ", end="")
        for R in R_values[1:]:
            dmu = results[R][ch] - free_mu[ch]
            print(f"  {dmu:12.4f}", end="")
        print()

    return results


def analyze_charge_function_decomposition(l_max: int = 3, n_alpha: int = 200):
    """
    Step 3: Decompose the charge function into SO(6) tensor components.

    The charge function C = -Z/r1 - Z/r2 + 1/r12 in hyperspherical coords:
    - Nuclear: C_nuc = -Z(1/cos(alpha) + 1/sin(alpha)) / R
    - Electron-electron: C_ee = 1/(R * sqrt(1 - sin(2alpha)*cos(theta_12)))

    In the Liouville basis, the nuclear part is diagonal in l (depends only on alpha).
    The e-e part couples different l via Gaunt integrals (Legendre expansion in theta_12).

    The nuclear part transforms as a scalar under angular rotations but is
    non-trivial in the hyperangle alpha. It belongs to specific SO(6)
    representations involving the alpha coordinate.
    """
    print("\n" + "=" * 72)
    print("STEP 3: Charge function SO(6) decomposition")
    print("=" * 72)

    # The nuclear potential in the Liouville-transformed equation contributes
    # V_nuc(alpha) = R * (-Z/cos(alpha) - Z/sin(alpha))
    # This is diagonal in l and only depends on alpha.
    #
    # The key question: what SO(6) representation does this belong to?
    #
    # On S^5, the nuclear potential -Z/r_i = -Z/(R*cos(alpha)) and -Z/(R*sin(alpha))
    # The functions 1/cos(alpha) and 1/sin(alpha) can be expanded in
    # hyperspherical harmonics Y_{nu,l}(alpha).
    #
    # For the L=0 sector, 1/cos(alpha) in the alpha basis {sin((n+1)*pi*alpha/(pi/2))}
    # has matrix elements that mix all n values.
    # This means the nuclear potential couples ALL SO(6) sectors.

    # Compute the nuclear potential matrix in the free eigenstate basis
    print("\nNuclear potential matrix elements <nu|V_nuc|nu'> in free basis:")
    print("(V_nuc = -Z/cos(a) - Z/sin(a), evaluated at l=0)")

    n_free = 6
    Z = 2.0

    # Get free eigenstates at R=0
    _, vecs_free = solve_angular(R=0.0, Z=Z, l_max=0, n_alpha=n_alpha,
                                 n_channels=n_free)

    # Build nuclear potential on alpha grid
    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    # V_nuc in the Liouville equation (for R*C, the R factors out)
    V_nuc = -Z / np.cos(alpha) - Z / np.sin(alpha)

    # Matrix elements <i|V_nuc|j>
    V_mat = np.zeros((n_free, n_free))
    for i in range(n_free):
        for j in range(n_free):
            V_mat[i, j] = h * np.dot(vecs_free[i], V_nuc * vecs_free[j])

    # Print the matrix
    header = "        "
    for j in range(n_free):
        nu_j = 2 * j
        header += f"  {'nu='+str(nu_j):>8}"
    print(header)

    for i in range(n_free):
        nu_i = 2 * i
        row = f"nu={nu_i:2d}  "
        for j in range(n_free):
            row += f"  {V_mat[i, j]:8.4f}"
        print(row)

    # Analyze: is V_nuc predominantly in a specific SO(6) irrep?
    print("\n  Diagonal dominance: |V_diag|/|V_total| = "
          f"{np.sum(np.abs(np.diag(V_mat))):.4f} / "
          f"{np.sum(np.abs(V_mat)):.4f} = "
          f"{np.sum(np.abs(np.diag(V_mat)))/np.sum(np.abs(V_mat)):.4f}")

    # Check if off-diagonal follows a pattern: V_{nu, nu'} ~ <nu|V|nu'>
    # For a rank-k SO(6) tensor, only nu' in [|nu-k|, nu+k] are non-zero
    print("\n  Selection rules check (which Dnu = nu'-nu transitions are nonzero):")
    for dnu in range(0, 2 * n_free, 2):
        vals = []
        for i in range(n_free):
            j = i + dnu // 2
            if j < n_free:
                vals.append(V_mat[i, j])
        if vals:
            max_val = max(abs(v) for v in vals)
            print(f"    Dnu={dnu:2d}: max |V| = {max_val:.6f}"
                  f"  {'<-- significant' if max_val > 0.01 else ''}")

    # Now do the same for the e-e potential
    print("\n\nElectron-electron potential matrix elements (l=0 channel, k=0 term):")
    print("V_ee(alpha) = 1/max(sin(a),cos(a)) for k=0 Legendre term")

    min_sc = np.minimum(np.sin(alpha), np.cos(alpha))
    max_sc = np.maximum(np.sin(alpha), np.cos(alpha))
    V_ee_k0 = 1.0 / max_sc  # k=0 monopole term

    V_ee_mat = np.zeros((n_free, n_free))
    for i in range(n_free):
        for j in range(n_free):
            V_ee_mat[i, j] = h * np.dot(vecs_free[i], V_ee_k0 * vecs_free[j])

    header = "        "
    for j in range(n_free):
        nu_j = 2 * j
        header += f"  {'nu='+str(nu_j):>8}"
    print(header)
    for i in range(n_free):
        nu_i = 2 * i
        row = f"nu={nu_i:2d}  "
        for j in range(n_free):
            row += f"  {V_ee_mat[i, j]:8.4f}"
        print(row)

    return V_mat, V_ee_mat


def analyze_perturbative_structure(l_max: int = 3, n_alpha: int = 300,
                                   n_channels: int = 6):
    """
    Step 4: Test if mu(R) follows a perturbative Casimir expansion.

    Hypothesis: mu_n(R) = C2_n/2 + a_n * R + b_n * R^2 + ...
    where a_n = <n|V_charge|n> (first-order correction).

    If a_n is determined by SO(6) Clebsch-Gordan coefficients, then
    there IS a closed-form algebraic relationship.
    """
    print("\n" + "=" * 72)
    print("STEP 4: Perturbative structure mu(R) = C2/2 + a*R + b*R^2 + ...")
    print("=" * 72)

    # Compute mu at many small R values to extract perturbation coefficients
    R_small = np.array([0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5])
    mu_data = np.zeros((n_channels, len(R_small)))

    for i, R in enumerate(R_small):
        mu, _ = solve_angular(R=R, Z=2.0, l_max=l_max, n_alpha=n_alpha,
                              n_channels=n_channels)
        mu_data[:, i] = mu

    # For each channel, fit mu(R) = a0 + a1*R + a2*R^2
    print(f"\n{'Ch':>3}  {'nu':>3}  {'a0 (=C2/2)':>12}  {'a1 (slope)':>12}  "
          f"{'a2 (curv.)':>12}  {'C2/2 exact':>12}  {'a0 err':>10}")
    print("-" * 80)

    for ch in range(n_channels):
        # Polynomial fit
        coeffs = np.polyfit(R_small, mu_data[ch], 2)
        a2, a1, a0 = coeffs

        # Expected nu
        nu = None
        for test_nu in range(0, 40, 2):
            if abs(a0 - test_nu * (test_nu + 4) / 2) < 1.0:
                nu = test_nu
                break

        c2_exact = nu * (nu + 4) / 2.0 if nu is not None else float('nan')
        nu_str = str(nu) if nu is not None else "?"

        print(f"{ch+1:3d}  {nu_str:>3s}  {a0:12.4f}  {a1:12.4f}  "
              f"{a2:12.4f}  {c2_exact:12.4f}  {a0 - c2_exact:10.6f}")

    # Check if a1 coefficients have a pattern
    print("\n  The a1 coefficients (linear in R) represent the diagonal matrix")
    print("  element of the charge function in the free eigenstate basis.")
    print("  If a1_n = f(nu_n) for some simple function f, the algebraic")
    print("  structure extends to first order in perturbation theory.")

    return mu_data, R_small


def analyze_large_R_asymptotics(l_max: int = 3, n_alpha: int = 300, n_channels: int = 4):
    """
    Step 5: Check large-R asymptotics against hydrogen quantum numbers.

    As R -> infinity, the system separates into He+(1s) + e-(nl).
    The angular eigenvalue should approach: mu(R) -> -Z^2 R^2 / (2n^2)
    for channel asymptotic to the n-th hydrogen state.
    """
    print("\n" + "=" * 72)
    print("STEP 5: Large-R asymptotics and hydrogen quantum numbers")
    print("=" * 72)

    R_large = np.array([5.0, 10.0, 20.0, 50.0, 100.0])
    Z = 2.0

    print(f"\n{'Ch':>3}", end="")
    for R in R_large:
        print(f"  {'V_eff(R='+str(int(R))+')':>14}", end="")
    print(f"  {'asymptotic n':>14}")
    print("-" * 90)

    for R in R_large:
        mu, _ = solve_angular(R=R, Z=Z, l_max=l_max, n_alpha=n_alpha,
                              n_channels=n_channels)
        if R == R_large[0]:
            mu_all = {R: mu}
        else:
            mu_all[R] = mu

    for ch in range(n_channels):
        row = f"{ch+1:3d}"
        for R in R_large:
            V_eff = mu_all[R][ch] / R**2 + 15.0 / (8.0 * R**2)
            row += f"  {V_eff:14.6f}"

        # Extract asymptotic n from largest R
        R_last = R_large[-1]
        V_eff_last = mu_all[R_last][ch] / R_last**2 + 15.0 / (8.0 * R_last**2)
        # V_eff -> -Z^2/(2n^2) for He+ threshold
        if V_eff_last < 0:
            n_eff = Z / np.sqrt(-2 * V_eff_last)
            row += f"  {n_eff:14.4f}"
        else:
            row += f"  {'unbound':>14}"
        print(row)

    print(f"\n  He+ thresholds: -Z^2/(2n^2) = ", end="")
    for n in range(1, 5):
        print(f"-{Z**2/(2*n**2):.4f} (n={n})  ", end="")
    print()


def summarize():
    """Print the theoretical summary."""
    print("\n" + "=" * 72)
    print("THEORETICAL SUMMARY")
    print("=" * 72)
    print("""
The hyperspherical angular equation for He at fixed hyperradius R is:

    [Lambda^2/2 + R * C(alpha, theta_12)] Phi = mu(R) Phi

where Lambda^2 is the grand angular momentum on S^5 and C is the charge
function containing nuclear attraction and electron repulsion.

SO(6) STRUCTURE:
  - S^5 = SO(6)/SO(5), so angular eigenstates carry SO(6) labels
  - For 1S (L=0, singlet), the allowed grand angular momentum is nu = 0, 2, 4, ...
  - The free (R=0) eigenvalues are EXACTLY mu = nu(nu+4)/2 = C_2[SO(6)]/2
  - This is the quadratic Casimir of the symmetric representation (nu, 0, 0)

CHARGE FUNCTION SYMMETRY BREAKING:
  - At finite R, the charge function C breaks SO(6) -> SO(3) x SO(3) x Z_2
    (individual angular momenta l1, l2 and particle exchange)
  - The L=0 restriction further reduces to: alpha dependence + Legendre coupling
  - Nuclear part (-Z/r1 - Z/r2) couples ALL even-Dnu sectors
  - Electron-electron part (1/r12) couples via Gaunt integrals with Delta_l selection

ANSWER: The free eigenvalues have a CLOSED-FORM Casimir expression.
At finite R, there is NO simple closed-form in terms of SO(6) Casimirs alone,
because the charge function is NOT an SO(6) scalar — it transforms as a
specific (non-trivial) SO(6) tensor that mixes all even-nu representations.

However, the PERTURBATIVE expansion mu(R) = C_2/2 + a_1*R + a_2*R^2 + ...
has coefficients determined by SO(6) Clebsch-Gordan coefficients of the
charge function tensor. This provides an algebraic (though not closed-form)
relationship.
""")


if __name__ == "__main__":
    print("SO(6) Casimir Analysis of He Angular Eigenvalues")
    print("=" * 72)

    # Step 1: Free spectrum
    mu_free = analyze_free_spectrum(l_max=3, n_alpha=400, n_channels=8)

    # Step 2: Finite R tracking
    results = analyze_finite_R(l_max=3, n_alpha=200, n_channels=6)

    # Step 3: Charge function decomposition
    V_nuc, V_ee = analyze_charge_function_decomposition(l_max=3, n_alpha=200)

    # Step 4: Perturbative structure
    mu_data, R_small = analyze_perturbative_structure(l_max=3, n_alpha=300,
                                                       n_channels=6)

    # Step 5: Large-R asymptotics
    analyze_large_R_asymptotics(l_max=3, n_alpha=300, n_channels=4)

    # Summary
    summarize()
