"""
Prolate Spheroidal Lattice for H2+ and one-electron diatomics.

The molecular analog of Fock's S3 projection: the two-center Coulomb
problem separates in prolate spheroidal coordinates (xi, eta), and the
eigenvalues are obtained by matching separated 1D equations — no free
parameters beyond grid resolution.

Coordinate system (nuclei A, B at foci, distance R apart):
  xi  in [1, inf)  — confocal ellipsoidal (r_A + r_B = R*xi)
  eta in [-1, +1]  — confocal hyperboloidal (r_A - r_B = R*eta)
  phi in [0, 2*pi) — azimuthal

Separated equations for m=0 (sigma states):
  xi:  d/dxi[(xi^2-1)dF/dxi] + (A + a*xi - c^2*xi^2) F = 0
  eta: d/deta[(1-eta^2)dG/deta] + (-A + c^2*eta^2 + b*eta) G = 0

where a = R*(Z_A+Z_B), b = R*(Z_B-Z_A), c^2 = -R^2*E_elec/2.
Self-consistency: both equations share the separation constant A.

The angular equation is solved spectrally (Legendre basis, exact).
The radial equation is solved with self-adjoint FD (Neumann at xi=1).
Root-finding in c^2 matches the two equations.

Reference: Bates, Ledsham & Stewart, Phil. Trans. A 246, 215 (1953).
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from typing import Dict, Optional, Tuple
import time

from geovac.molecular_sturmian import _angular_sep_const


class ProlateSpheroidalLattice:
    """
    Separated prolate spheroidal solver for one-electron diatomics.

    Solves H2+ (or HeH2+, etc.) by matching the separated xi and eta
    equations in prolate spheroidal coordinates. No free parameters.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    Z_A : int
        Charge of nucleus A.
    Z_B : int
        Charge of nucleus B.
    N_xi : int
        Number of radial grid points.
    xi_max : float
        Maximum xi value.
    m : int
        Azimuthal quantum number (0 for sigma states).
    """

    def __init__(
        self,
        R: float,
        Z_A: int = 1,
        Z_B: int = 1,
        N_xi: int = 5000,
        xi_max: float = 25.0,
        m: int = 0,
    ):
        self.R = R
        self.Z_A = Z_A
        self.Z_B = Z_B
        self.N_xi = N_xi
        self.xi_max = xi_max
        self.m = m

        # Derived
        self._a = R * (Z_A + Z_B)
        self._b = R * (Z_B - Z_A)

    def _radial_top_eigenvalue(
        self, c2: float, A: float
    ) -> float:
        """Top eigenvalue of the radial operator.

        L_xi = d/dxi[(xi^2-1)dF/dxi] + (A + a*xi - c^2*xi^2 - m^2/(xi^2-1))

        Returns the largest eigenvalue; should be 0 at the correct c^2.
        """
        N = self.N_xi
        xi_min = 1.0 + 5e-4
        h = (self.xi_max - xi_min) / (N + 1)
        xi = xi_min + (np.arange(N) + 1) * h

        xi2_1 = xi**2 - 1
        q = A + self._a * xi - c2 * xi**2
        if self.m != 0:
            q -= self.m**2 / xi2_1

        p_plus = (xi + h / 2)**2 - 1
        p_minus = (xi - h / 2)**2 - 1

        diag = -(p_plus + p_minus) / h**2 + q
        off = p_plus[:-1] / h**2

        # Neumann BC at xi=1 for m=0
        if self.m == 0:
            diag[0] = -p_plus[0] / h**2 + q[0]

        evals = eigh_tridiagonal(diag, off, eigvals_only=True)
        return np.max(evals)

    def _angular_separation_constant(self, c2: float) -> float:
        """Angular separation constant A for given c^2.

        Uses the spectral (Legendre basis) solver from molecular_sturmian.
        """
        c = np.sqrt(max(c2, 1e-15))
        return _angular_sep_const(
            self.m, 0, c, b=self._b, n_basis=50
        )

    def _residual(self, c2: float) -> float:
        """Residual for root-finding: should be 0 at the correct c^2."""
        A = self._angular_separation_constant(c2)
        return self._radial_top_eigenvalue(c2, A)

    def solve(self) -> Tuple[float, float, float]:
        """Solve for the ground state electronic energy.

        Returns
        -------
        E_elec : float
            Electronic energy (Ha), not including V_NN.
        c2 : float
            Separation parameter c^2 = -R^2*E_elec/2.
        A : float
            Angular separation constant.
        """
        # Adaptive search range for c^2
        # c^2 = -R^2*E/2 > 0 for bound states
        # Upper bound: E can't be more negative than -2*(Z_A+Z_B)^2/R^2 (united atom)
        c2_max = self.R**2 * (self.Z_A + self.Z_B)**2
        c2_min = 0.01

        # Verify bracket
        f_lo = self._residual(c2_min)
        if f_lo < 0:
            raise ValueError(
                f"Residual negative at c2_min={c2_min}. "
                "No bound state found."
            )

        # Find upper bracket
        f_hi = self._residual(c2_max)
        if f_hi > 0:
            # Extend search
            for c2_try in np.linspace(c2_max, c2_max * 2, 20):
                if self._residual(c2_try) < 0:
                    c2_max = c2_try
                    break
            else:
                raise ValueError(
                    f"Could not bracket root up to c2={c2_max*2}"
                )

        c2_sol = brentq(self._residual, c2_min, c2_max, xtol=1e-12)
        E_elec = -2.0 * c2_sol / self.R**2
        A = self._angular_separation_constant(c2_sol)

        return E_elec, c2_sol, A

    def total_energy(self) -> float:
        """Electronic energy + nuclear repulsion."""
        E_elec, _, _ = self.solve()
        return E_elec + self.Z_A * self.Z_B / self.R


def scan_h2plus_pes(
    R_values: np.ndarray,
    Z_A: int = 1,
    Z_B: int = 1,
    N_xi: int = 5000,
    xi_max: float = 25.0,
    m: int = 0,
    verbose: bool = True,
) -> Dict[str, np.ndarray]:
    """Scan H2+ PES over a range of R values."""
    n = len(R_values)
    E_elec = np.zeros(n)
    E_total = np.zeros(n)

    t0 = time.time()
    for i, R in enumerate(R_values):
        lattice = ProlateSpheroidalLattice(
            R=R, Z_A=Z_A, Z_B=Z_B,
            N_xi=N_xi, xi_max=max(xi_max, R * 3), m=m
        )
        try:
            E_el, _, _ = lattice.solve()
            E_elec[i] = E_el
            E_total[i] = E_el + Z_A * Z_B / R
        except ValueError:
            E_elec[i] = np.nan
            E_total[i] = np.nan
        if verbose:
            print(f"  R={R:6.3f}  E_elec={E_elec[i]:10.6f}  "
                  f"E_total={E_total[i]:10.6f}")

    dt = time.time() - t0
    if verbose:
        print(f"  PES scan: {n} points in {dt:.1f}s")

    return {
        'R': R_values,
        'E_elec': E_elec,
        'E_total': E_total,
        'V_NN': Z_A * Z_B / R_values,
    }


def fit_spectroscopic_constants(
    R_values: np.ndarray, E_values: np.ndarray
) -> Dict[str, float]:
    """Fit PES near minimum to extract R_eq, E_min, D_e, k."""
    valid = ~np.isnan(E_values)
    R_v = R_values[valid]
    E_v = E_values[valid]

    idx_min = np.argmin(E_v)

    if idx_min == 0 or idx_min == len(R_v) - 1:
        return {
            'R_eq': R_v[idx_min],
            'E_min': E_v[idx_min],
            'D_e': -0.5 - E_v[idx_min],
            'k': 0.0,
            'boundary': True,
        }

    lo = max(0, idx_min - 3)
    hi = min(len(R_v), idx_min + 4)
    coeffs = np.polyfit(R_v[lo:hi], E_v[lo:hi], 2)
    R_eq = -coeffs[1] / (2 * coeffs[0])
    E_min = np.polyval(coeffs, R_eq)

    return {
        'R_eq': R_eq,
        'E_min': E_min,
        'D_e': -0.5 - E_min,
        'k': 2 * coeffs[0],
        'boundary': False,
    }
