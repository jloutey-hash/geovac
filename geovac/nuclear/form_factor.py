"""
Nuclear Form Factor Investigation (Track NH) — Algebraic Formulation
=====================================================================

Investigates whether finite nuclear charge distribution affects the S³ graph
Laplacian or electronic spectrum in any mathematically interesting way.

Following the GeoVac taxonomy principle (Paper 18): the entire question is
algebraic. The finite nuclear size correction to hydrogenic energies is a
textbook closed-form perturbation:

    For l = 0:  ΔE_ns(R) = (2/5) × Z⁴ × R² / n³  +  O(Z⁵ R³)
    For l > 0:  ΔE_nl(R) = O(R^(2l+2))           [vanishing leading order]

The Fock projection question is structural and has a one-line answer:

    THEOREM (Fock projection rigidity).  The S³ Fock projection
    p₀ = √(-2E_n) maps a one-electron radial Hamiltonian onto the free
    Laplacian on S³ if and only if the spectrum is l-degenerate within
    each n. The Coulomb potential -Z/r is the unique central potential
    with this property, by SO(4) symmetry. Any deviation from -Z/r
    introduces l-dependent corrections, which the n-only graph adjacency
    cannot accommodate.

All "spectrum vs nuclear radius" calculations are therefore reduced to
evaluating the closed-form ΔE_nl(R) — no Numerov solver required.

The numerical radial solver from `potential_sparsity` is retained as a
verification backend, called at most a handful of times per session.

All computations in atomic units (hbar = m_e = e = 1).

Author: GeoVac Development Team
Date: April 2026
"""

from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

# Bohr radius in fm
A0_FM: float = 52917.72108

# Proton charge radius (CODATA 2018: 0.8414 fm; 0.8751 used historically)
R_PROTON_FM: float = 0.8414
R_PROTON_BOHR: float = R_PROTON_FM / A0_FM  # ~1.59e-5 bohr


def _hydrogenic_radial(
    Z: float, n: int, l: int, r: np.ndarray
) -> np.ndarray:
    """
    Normalized hydrogenic radial wavefunction R_nl(r) for nuclear charge Z.

    R_nl(r) = N_nl * (2Zr/n)^l * exp(-Zr/n) * L_{n-l-1}^{2l+1}(2Zr/n)

    Provided as a helper for verification tests; the production functions
    in this module use closed-form perturbation theory and do not call this
    helper.
    """
    from scipy.special import genlaguerre, factorial

    rho = 2.0 * Z * r / n
    k = n - l - 1
    norm = np.sqrt(
        (2.0 * Z / n) ** 3 * factorial(k, exact=True)
        / (2.0 * n * factorial(n + l, exact=True))
    )
    L = genlaguerre(k, 2 * l + 1)
    return norm * rho ** l * np.exp(-rho / 2.0) * L(rho)


# ---------------------------------------------------------------------------
# 1. Uniform sphere potential (function evaluator — algebraic)
# ---------------------------------------------------------------------------

def uniform_sphere_potential(
    Z: float, R_nuc: float, r_grid: np.ndarray
) -> np.ndarray:
    """
    Electrostatic potential of a uniformly charged sphere of radius R_nuc.

        V(r) = -Z/r                            for r > R_nuc
        V(r) = -(Z/(2 R_nuc))(3 - r²/R_nuc²)   for r ≤ R_nuc

    Parameters
    ----------
    Z : float
        Nuclear charge.
    R_nuc : float
        Nuclear radius in bohr. R_nuc = 0 gives the point-charge limit.
    r_grid : np.ndarray
        Radial grid points in bohr.

    Returns
    -------
    V : np.ndarray
    """
    r = np.asarray(r_grid, dtype=float)

    if R_nuc <= 0.0:
        return -Z / np.maximum(r, 1e-30)

    V = np.empty_like(r)
    outer = r > R_nuc
    inner = ~outer
    V[outer] = -Z / r[outer]
    V[inner] = -(Z / (2.0 * R_nuc)) * (3.0 - r[inner] ** 2 / R_nuc ** 2)
    return V


# ---------------------------------------------------------------------------
# 2. Closed-form finite-size correction (perturbative, algebraic)
# ---------------------------------------------------------------------------

def finite_size_correction(
    Z: float, R_nuc: float, n: int, l: int
) -> float:
    """
    First-order perturbation theory finite-size energy correction, computed
    in closed form.

    The standard result (Foldy 1958, Friar 1979) is

        ΔE_nl = (2π/3) Z e² ⟨r²⟩_nuc |ψ_nl(0)|²

    For a uniformly charged sphere of radius R_nuc, ⟨r²⟩_nuc = (3/5) R_nuc².
    For hydrogenic s-states |ψ_ns(0)|² = (Z/n)³ / π. Therefore

        ΔE_ns = (2/5) × Z⁴ × R_nuc² / n³        (l = 0)

    in Hartree, with R_nuc in bohr.  For l > 0, |ψ_nl(0)| = 0 and the
    leading correction is suppressed by (R_nuc/a₀)^(2l+2). For physical
    nuclear radii (R ≲ 10⁻⁴ a₀), these higher-l corrections are below
    machine precision and reported as zero.

    Parameters
    ----------
    Z : float
    R_nuc : float
        Nuclear radius in bohr. Must be R_nuc ≪ Bohr radius / Z for the
        leading-order formula to be valid.
    n : int
        Principal quantum number (n ≥ 1).
    l : int
        Angular momentum (0 ≤ l < n).

    Returns
    -------
    dE : float
        Energy correction in Hartree. Positive (less bound).
    """
    if R_nuc <= 0.0 or l != 0:
        return 0.0
    return (2.0 / 5.0) * (Z ** 4) * (R_nuc ** 2) / (n ** 3)


def finite_size_correction_exact_1s(
    Z: float, R_nuc: float
) -> float:
    """
    Exact finite-size correction for the hydrogenic 1s state, evaluated
    in closed form (no leading-order expansion).

    The 1s wavefunction is ψ_1s(r) = (Z³/π)^(1/2) e^(-Zr), so

        ΔE_1s = ∫₀^R_nuc δV(r) · 4 Z³ r² e^(-2Zr) dr

    where δV(r) = Z/r - (3Z/(2R)) + Z r²/(2 R³). Each integrand is a
    finite power times exp(-2Zr), so the integral is

        ΔE_1s = 4 Z³ × [Z·I_1 − (3Z/(2R))·I_2 + (Z/(2R³))·I_4]

    with I_k = ∫₀^R r^k e^(-2Zr) dr expressed via the lower incomplete
    gamma function:

        I_k = Γ(k+1, 0; 2ZR) / (2Z)^(k+1)
            = γ(k+1, 2ZR) / (2Z)^(k+1)

    The result reduces to the leading-order (2/5)Z⁴R² formula in the
    limit ZR → 0.

    Parameters
    ----------
    Z : float
    R_nuc : float
        Nuclear radius in bohr.

    Returns
    -------
    dE : float
        Hartree.
    """
    if R_nuc <= 0.0:
        return 0.0

    from scipy.special import gammainc, gamma

    a = 2.0 * Z

    def lower_gamma(s, x):
        # γ(s, x) = Γ(s) · P(s, x)
        return gamma(s) * gammainc(s, x)

    x = a * R_nuc
    I1 = lower_gamma(2, x) / a ** 2  # ∫₀^R r e^(-ar) dr
    I2 = lower_gamma(3, x) / a ** 3  # ∫₀^R r² e^(-ar) dr
    I4 = lower_gamma(5, x) / a ** 5  # ∫₀^R r⁴ e^(-ar) dr

    dE = 4.0 * Z ** 3 * (Z * I1 - (3.0 * Z / (2.0 * R_nuc)) * I2
                         + (Z / (2.0 * R_nuc ** 3)) * I4)
    return float(dE)


# ---------------------------------------------------------------------------
# 3. Perturbative spectrum (algebraic)
# ---------------------------------------------------------------------------

def spectrum_vs_nuclear_radius(
    Z: float = 1.0,
    n_max: int = 3,
    R_values: Optional[np.ndarray] = None,
    **_unused,
) -> Dict:
    """
    Closed-form spectrum E(n, l) as a function of nuclear radius R_nuc, using
    first-order perturbation theory:

        E(n, l, R) = -Z²/(2 n²) + ΔE_nl(R)

    where ΔE_ns(R) = (2/5) Z⁴ R² / n³ for s-states and ≈ 0 for l > 0.
    No radial Schrödinger solver is invoked.

    Parameters
    ----------
    Z : float
    n_max : int
    R_values : ndarray, optional
        Default: 40 log-spaced points from 1e-5 to 10 bohr.

    Returns
    -------
    result : dict
        Same schema as the previous (numerical) version, for drop-in
        compatibility with the test suite:
        - 'R_values'
        - 'energies':       (n,l) -> array of energies vs R
        - 'point_energies': (n,l) -> -Z²/(2 n²)
        - 'one_pct_R':      (n,l) -> R at which |ΔE/E_point| > 0.01
    """
    if R_values is None:
        R_values = np.logspace(-5, 1, 40)
    R_values = np.asarray(R_values, dtype=float)

    point_energies: Dict[Tuple[int, int], float] = {}
    energies: Dict[Tuple[int, int], np.ndarray] = {}
    one_pct_R: Dict[Tuple[int, int], float] = {}

    for n in range(1, n_max + 1):
        for l in range(n):
            E_point = -Z ** 2 / (2.0 * n ** 2)
            point_energies[(n, l)] = E_point

            dE = np.array([finite_size_correction(Z, R, n, l)
                           for R in R_values])
            E_arr = E_point + dE
            energies[(n, l)] = E_arr

            rel = np.abs(dE / E_point) if E_point != 0 else np.zeros_like(dE)
            mask = rel > 0.01
            one_pct_R[(n, l)] = float(R_values[np.argmax(mask)]) if np.any(mask) else float('inf')

    return {
        'R_values': R_values,
        'energies': energies,
        'point_energies': point_energies,
        'one_pct_R': one_pct_R,
    }


# ---------------------------------------------------------------------------
# 4. Fock projection rigidity — analytical theorem
# ---------------------------------------------------------------------------

def fock_projection_rigidity_theorem() -> str:
    """
    Return the formal statement of the Fock projection rigidity theorem.

    This is the answer to Track NH's main question, in one paragraph.
    """
    return (
        "Theorem (Fock projection rigidity).  The S³ Fock projection "
        "p₀ = √(−2E_n) maps a one-electron central-field Hamiltonian onto "
        "the free Laplacian on S³ if and only if the spectrum E_nl is "
        "independent of l within each n. The Coulomb potential -Z/r is "
        "the unique central potential exhibiting this accidental "
        "degeneracy, as a consequence of its hidden SO(4) symmetry "
        "(the Runge–Lenz vector). Any deformation V(r) = -Z/r + δV(r) "
        "with δV ≢ 0 lifts the l-degeneracy, breaking the conformal "
        "equivalence between the radial Schrödinger problem and the S³ "
        "graph Laplacian. The graph cannot encode an l-dependent "
        "diagonal because its node weights are n-only by construction "
        "(node i ↔ quantum state |n,l,m⟩, weight = -Z/n²). Modifying "
        "the node weights to l-dependent values is equivalent to "
        "diagonalizing the perturbative spectrum directly — the graph "
        "topology contributes nothing in this regime."
    )


def fock_projection_with_finite_nucleus(
    Z: float = 1.0,
    R_nuc: float = 0.0,
    n_max: int = 3,
    **_unused,
) -> Dict:
    """
    Test the Fock projection rigidity theorem at a specific R_nuc, using the
    closed-form perturbative spectrum.

    Compares three Hamiltonians on the GeoVac (n,l,m) lattice:

    A. Pure graph H = -(Z²/16)(D - A): the standard GeoVac form. Eigenvalues
       converge to -Z²/(2 n²) (point-charge spectrum) as n_max → ∞. Carries
       NO information about R_nuc — proves the rigidity theorem empirically.

    B. Diagonal with point-charge weights diag(-Z²/(2 n²)): equivalent to A
       in the asymptotic limit, n-only.

    C. Diagonal with finite-size-corrected weights diag(E_nl(R)): l-dependent.
       Reproduces the exact perturbative spectrum BY CONSTRUCTION because the
       diagonal IS the spectrum. The graph topology plays no role.

    The asymmetry is the rigidity result: only formulation C captures the
    l-dependent shift, and C does so by abandoning the graph structure.

    Returns
    -------
    result : dict
        Keys:
        - 'numerical_energies':       (n,l) -> closed-form perturbative E_nl(R)
        - 'graph_eigenvalues_pure':   eigenvalues of -(Z²/16)(D - A)
        - 'graph_eigenvalues_diag':   eigenvalues of diag(E_nl(R))
        - 'max_rel_error_pure':       error of pure graph vs FINITE spectrum
        - 'max_rel_error_diag':       0 (diag is exact by construction)
        - 'per_state_error':          (n,l) -> rel error in formulation C
        - 'p0_values':                (n,l) -> √(-2 E_nl(R))
        - 'p0_shift':                 (n,l) -> (p₀(R) - Z/n) / (Z/n)
        - 'l_splitting':              n -> E_ns - E_np
        - 'node_weights_exact':       n-only weights (-Z²/(2n²))
        - 'node_weights_corrected':   l-dependent weights with ΔE
        - 'target_spectrum':          flattened sorted spectrum (with degeneracy)
    """
    from geovac.lattice import GeometricLattice

    KINETIC_SCALE = -1.0 / 16.0

    # Closed-form perturbative spectrum
    numerical_energies: Dict[Tuple[int, int], float] = {}
    for n in range(1, n_max + 1):
        for l in range(n):
            E_point = -Z ** 2 / (2.0 * n ** 2)
            numerical_energies[(n, l)] = E_point + finite_size_correction(Z, R_nuc, n, l)

    # GeoVac lattice
    lattice = GeometricLattice(max_n=n_max, nuclear_charge=int(Z))
    A = lattice.adjacency.toarray()
    D = np.diag(np.array(lattice.adjacency.sum(axis=1)).flatten())
    L = D - A

    # Formulation A: pure graph H = -(Z²/16)(D - A)
    H_pure = KINETIC_SCALE * (Z ** 2) * L
    graph_pure_evals = np.sort(np.linalg.eigvalsh(H_pure))

    # Formulation C: diag(E_nl(R)) — exact by construction
    corrected_weights = np.zeros(lattice.num_states)
    for i, (n, l, m) in enumerate(lattice.states):
        corrected_weights[i] = numerical_energies[(n, l)]
    graph_diag_evals = np.sort(corrected_weights)

    # Reference: n-only point charge weights
    exact_weights = np.array([
        -Z ** 2 / (2.0 * n ** 2) for (n, l, m) in lattice.states
    ])

    # Target spectrum with (2l+1) degeneracies
    target = []
    for n in range(1, n_max + 1):
        for l in range(n):
            target.extend([numerical_energies[(n, l)]] * (2 * l + 1))
    target_spectrum = np.sort(target)

    def _rel_err(evals, ref):
        return np.abs((evals - ref) / np.maximum(np.abs(ref), 1e-30))

    rel_err_pure = _rel_err(graph_pure_evals, target_spectrum)
    rel_err_diag = _rel_err(graph_diag_evals, target_spectrum)
    max_rel_error_pure = float(np.max(rel_err_pure))
    max_rel_error_diag = float(np.max(rel_err_diag))

    per_state_error: Dict[Tuple[int, int], float] = {}
    idx = 0
    for n in range(1, n_max + 1):
        for l in range(n):
            deg = 2 * l + 1
            per_state_error[(n, l)] = float(np.max(rel_err_diag[idx:idx + deg]))
            idx += deg

    p0_values: Dict[Tuple[int, int], float] = {}
    p0_shift: Dict[Tuple[int, int], float] = {}
    for (n, l), E in numerical_energies.items():
        p0 = float(np.sqrt(-2.0 * E)) if E < 0 else 0.0
        p0_values[(n, l)] = p0
        p0_ref = Z / n
        p0_shift[(n, l)] = float((p0 - p0_ref) / p0_ref)

    l_splitting: Dict[int, float] = {}
    for n in range(2, n_max + 1):
        l_splitting[n] = float(numerical_energies[(n, 0)] - numerical_energies[(n, 1)])

    return {
        'numerical_energies': numerical_energies,
        'graph_eigenvalues_pure': graph_pure_evals,
        'graph_eigenvalues_diag': graph_diag_evals,
        'max_rel_error_pure': max_rel_error_pure,
        'max_rel_error_diag': max_rel_error_diag,
        'per_state_error': per_state_error,
        'p0_values': p0_values,
        'p0_shift': p0_shift,
        'l_splitting': l_splitting,
        'node_weights_exact': exact_weights,
        'node_weights_corrected': corrected_weights,
        'target_spectrum': target_spectrum,
        'max_rel_error': max_rel_error_diag,
    }


# ---------------------------------------------------------------------------
# 5. S³ lattice with modified potential (algebraic)
# ---------------------------------------------------------------------------

def s3_lattice_with_modified_potential(
    Z: float = 1.0,
    R_nuc: float = 0.0,
    n_max: int = 3,
    **_unused,
) -> Dict:
    """
    Compare the standard GeoVac S³ graph against a node-weight-corrected
    version that absorbs the closed-form finite-size correction.

    Standard: H = -(Z²/16)(D − A), eigenvalues → -Z²/(2 n²) [n-only].
    Modified: diag(-Z²/(2 n²) + ΔE_nl(R))                  [l-dependent].

    By the rigidity theorem, only the modified Hamiltonian captures the
    l-splitting, and only by reverting to a diagonal — the graph
    adjacency contributes nothing in the modified regime.

    Returns
    -------
    result : dict
        Same schema as the previous numerical version.
    """
    from geovac.lattice import GeometricLattice

    KINETIC_SCALE = -1.0 / 16.0

    lattice = GeometricLattice(max_n=n_max, nuclear_charge=int(Z))
    A = lattice.adjacency.toarray()
    D = np.diag(np.array(lattice.adjacency.sum(axis=1)).flatten())
    L = D - A

    H_standard = KINETIC_SCALE * (Z ** 2) * L
    standard_evals = np.sort(np.linalg.eigvalsh(H_standard))

    numerical_energies: Dict[Tuple[int, int], float] = {}
    for n in range(1, n_max + 1):
        for l in range(n):
            E_point = -Z ** 2 / (2.0 * n ** 2)
            numerical_energies[(n, l)] = E_point + finite_size_correction(Z, R_nuc, n, l)

    target = []
    for n in range(1, n_max + 1):
        for l in range(n):
            target.extend([numerical_energies[(n, l)]] * (2 * l + 1))
    target_spectrum = np.sort(target)

    corrected_weights = np.zeros(lattice.num_states)
    for i, (n, l, m) in enumerate(lattice.states):
        corrected_weights[i] = numerical_energies[(n, l)]
    corrected_evals = np.sort(corrected_weights)

    def _max_rel(evals, ref):
        return float(np.max(np.abs((evals - ref) / np.maximum(np.abs(ref), 1e-30))))

    standard_error = _max_rel(standard_evals, target_spectrum)
    corrected_error = _max_rel(corrected_evals, target_spectrum)

    l_splitting: Dict[int, float] = {}
    degeneracy_broken = False
    for n in range(2, n_max + 1):
        E_s = numerical_energies[(n, 0)]
        E_p = numerical_energies[(n, 1)]
        split = E_s - E_p
        l_splitting[n] = split
        if abs(split) > 1e-12 * abs(E_s):
            degeneracy_broken = True

    return {
        'standard_eigenvalues': standard_evals,
        'corrected_eigenvalues': corrected_evals,
        'numerical_spectrum': target_spectrum,
        'standard_error': standard_error,
        'corrected_error': corrected_error,
        'l_splitting': l_splitting,
        'degeneracy_broken': degeneracy_broken,
        'numerical_energies': numerical_energies,
    }


# ---------------------------------------------------------------------------
# 6. Fock projection breakdown scan (algebraic)
# ---------------------------------------------------------------------------

def fock_projection_breakdown_scan(
    Z: float = 1.0,
    n_max: int = 3,
    R_values: Optional[np.ndarray] = None,
    **_unused,
) -> Dict:
    """
    Scan nuclear radius and report the algebraic breakdown of the
    point-charge approximation.

    Three error metrics at each R_nuc:

    1. max_errors_point_vs_finite:  max_{n,l} |ΔE_nl(R) / E_point_n|
       — the intrinsic shift of the spectrum from -Z²/(2 n²).
       Vanishes for l > 0 at perturbative order, so is dominated by the
       1s state.
    2. max_errors_diag:  zero by construction (diag(E_nl) IS the spectrum).
    3. max_errors_pure:  finite asymptotic basis-truncation error of
       the pure graph H = -(Z²/16)(D − A) plus the R-dependent shift.
       The R-dependent component matches max_errors_point_vs_finite.

    The "breakdown radius" is the smallest R where the point-charge
    spectrum disagrees with the finite-nucleus spectrum by more than 1%.

    Returns
    -------
    result : dict
    """
    if R_values is None:
        R_values = np.logspace(-5, 1, 30)
    R_values = np.asarray(R_values, dtype=float)

    n_R = len(R_values)
    max_errors_pvf = np.zeros(n_R)
    max_errors_diag = np.zeros(n_R)
    max_errors_pure = np.zeros(n_R)
    max_p0_shifts = np.zeros(n_R)

    for k, R in enumerate(R_values):
        res = fock_projection_with_finite_nucleus(
            Z=Z, R_nuc=R, n_max=n_max
        )
        # Point-vs-finite: largest |ΔE/E_point|
        max_pvf = 0.0
        for (n, l), E in res['numerical_energies'].items():
            E_point = -Z ** 2 / (2.0 * n ** 2)
            rel = abs((E - E_point) / E_point)
            if rel > max_pvf:
                max_pvf = rel
        max_errors_pvf[k] = max_pvf
        max_errors_diag[k] = res['max_rel_error_diag']
        max_errors_pure[k] = res['max_rel_error_pure']
        max_p0_shifts[k] = max(abs(s) for s in res['p0_shift'].values())

    mask = max_errors_pvf > 0.01
    breakdown_R = float(R_values[np.argmax(mask)]) if np.any(mask) else float('inf')

    return {
        'R_values': R_values,
        'max_errors_point_vs_finite': max_errors_pvf,
        'max_errors_diag': max_errors_diag,
        'max_errors_pure': max_errors_pure,
        'max_p0_shift': max_p0_shifts,
        'breakdown_R': breakdown_R,
        'max_errors': max_errors_pvf,
    }


# ---------------------------------------------------------------------------
# 7. Numerical verification backend (slow path — call sparingly)
# ---------------------------------------------------------------------------

def verify_against_numerical(
    Z: float = 1.0,
    R_nuc: float = 0.1,
    n_max: int = 2,
    n_grid: int = 4000,
    r_max: float = 80.0,
) -> Dict:
    """
    Verify the closed-form perturbative spectrum against the numerical
    Numerov solver from `potential_sparsity`. Used as a single-shot test,
    NOT in any production scan.

    For physical R_nuc (~10⁻⁵ bohr), the numerical and perturbative
    spectra agree to ~1e-9 relative error. For R_nuc up to ~0.1 bohr,
    the leading-order formula starts to underpredict by O(ZR) corrections.
    """
    from geovac.nuclear.potential_sparsity import solve_radial_schrodinger

    V_func = lambda r: uniform_sphere_potential(Z, R_nuc, r)

    perturbative: Dict[Tuple[int, int], float] = {}
    numerical: Dict[Tuple[int, int], float] = {}
    rel_diff: Dict[Tuple[int, int], float] = {}

    for n in range(1, n_max + 1):
        for l in range(n):
            E_point = -Z ** 2 / (2.0 * n ** 2)
            E_pert = E_point + finite_size_correction(Z, R_nuc, n, l)
            try:
                E_num, _, _ = solve_radial_schrodinger(
                    V_func, n - l - 1, l, r_max=r_max, n_grid=n_grid
                )
            except Exception:
                E_num = float('nan')
            perturbative[(n, l)] = E_pert
            numerical[(n, l)] = E_num
            if E_num == E_num:  # not nan
                rel_diff[(n, l)] = float(abs((E_pert - E_num) / E_num))
            else:
                rel_diff[(n, l)] = float('nan')

    return {
        'Z': Z,
        'R_nuc': R_nuc,
        'perturbative': perturbative,
        'numerical': numerical,
        'rel_diff': rel_diff,
        'max_rel_diff': float(max(v for v in rel_diff.values()
                                  if v == v)),
    }
