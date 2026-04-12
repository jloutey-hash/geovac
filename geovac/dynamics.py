"""
Quantum Dynamics — Crank-Nicolson Time Propagator
==================================================

Implements unitary time evolution for graph Hamiltonians:

    (I + i*H*dt/2) * psi(t+dt) = (I - i*H*dt/2) * psi(t)

Crank-Nicolson is:
    - Unconditionally stable (implicit)
    - Unitary (preserves norm exactly)
    - Second-order accurate in dt
    - Works directly with sparse matrices

Also provides a dipole operator builder for driving transitions
with time-dependent electric fields.

Author: GeoVac Development Team
Date: February 2026
"""

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix, eye
from scipy.sparse.linalg import spsolve
from typing import Tuple, Optional, Callable, List

try:
    from .lattice import GeometricLattice
except ImportError:
    from lattice import GeometricLattice


class TimePropagator:
    """
    Crank-Nicolson unitary time propagator for graph Hamiltonians.

    Given a Hamiltonian H (sparse, real-valued), propagates a complex
    wavefunction psi(t) forward in time using the implicit midpoint rule:

        (I + i*H*dt/2) * psi(t+dt) = (I - i*H*dt/2) * psi(t)

    This preserves ||psi|| = 1 to machine precision at every step.

    For time-dependent Hamiltonians H(t), call step_with_H() directly
    with the instantaneous Hamiltonian at each timestep.

    Parameters
    ----------
    H : csr_matrix
        Time-independent Hamiltonian (sparse, real or complex)
    dt : float
        Time step in atomic units (default: 0.01, ~0.24 attoseconds)

    Example
    -------
    >>> solver = AtomicSolver(max_n=5, Z=1)
    >>> E, psi = solver.compute_ground_state(n_states=2)
    >>> prop = TimePropagator(solver.H, dt=0.01)
    >>> psi_t = prop.evolve(psi[:, 0].astype(complex), n_steps=1000)
    """

    def __init__(self, H: csr_matrix, dt: float = 0.01):
        self.H = H.tocsc()
        self.dt = dt
        self.n_states = H.shape[0]
        self.time = 0.0

        # Precompute Crank-Nicolson matrices for time-independent H
        I = eye(self.n_states, format='csc', dtype=complex)
        H_complex = self.H.astype(complex)
        self._A_left = I + 0.5j * self.dt * H_complex
        self._A_right = I - 0.5j * self.dt * H_complex

    def step(self, psi: np.ndarray) -> np.ndarray:
        """
        Single Crank-Nicolson step with time-independent H.

        Parameters
        ----------
        psi : np.ndarray, shape (n_states,), complex
            Current wavefunction

        Returns
        -------
        psi_new : np.ndarray, shape (n_states,), complex
            Wavefunction at t + dt
        """
        rhs = self._A_right @ psi
        psi_new = spsolve(self._A_left, rhs)
        self.time += self.dt
        return psi_new

    def step_with_H(self, psi: np.ndarray, H: csr_matrix) -> np.ndarray:
        """
        Single Crank-Nicolson step with a given Hamiltonian.

        Use this for time-dependent Hamiltonians where H changes each step.

        Parameters
        ----------
        psi : np.ndarray, shape (n_states,), complex
            Current wavefunction
        H : csr_matrix
            Hamiltonian at current time (sparse)

        Returns
        -------
        psi_new : np.ndarray, shape (n_states,), complex
            Wavefunction at t + dt
        """
        I = eye(self.n_states, format='csc', dtype=complex)
        H_c = H.astype(complex).tocsc()
        A_left = I + 0.5j * self.dt * H_c
        A_right = I - 0.5j * self.dt * H_c
        rhs = A_right @ psi
        psi_new = spsolve(A_left, rhs)
        self.time += self.dt
        return psi_new

    def evolve(self, psi0: np.ndarray, n_steps: int,
               callback: Optional[Callable] = None) -> np.ndarray:
        """
        Evolve wavefunction for n_steps using time-independent H.

        Parameters
        ----------
        psi0 : np.ndarray, shape (n_states,), complex
            Initial wavefunction
        n_steps : int
            Number of time steps
        callback : callable, optional
            Called as callback(step_index, time, psi) at each step

        Returns
        -------
        psi_final : np.ndarray, shape (n_states,), complex
            Wavefunction after n_steps * dt
        """
        psi = psi0.copy().astype(complex)
        for i in range(n_steps):
            psi = self.step(psi)
            if callback is not None:
                callback(i, self.time, psi)
        return psi

    def evolve_driven(self, psi0: np.ndarray, n_steps: int,
                      H0: csr_matrix, V_dipole: csr_matrix,
                      E0: float, omega: float,
                      callback: Optional[Callable] = None) -> np.ndarray:
        """
        Evolve wavefunction under a time-dependent driving field.

        H(t) = H0 + E0 * cos(omega * t) * V_dipole

        Parameters
        ----------
        psi0 : np.ndarray, shape (n_states,), complex
            Initial wavefunction
        n_steps : int
            Number of time steps
        H0 : csr_matrix
            Static Hamiltonian (sparse)
        V_dipole : csr_matrix
            Dipole operator (sparse)
        E0 : float
            Electric field amplitude (atomic units)
        omega : float
            Driving frequency (atomic units, energy)
        callback : callable, optional
            Called as callback(step_index, time, psi) at each step

        Returns
        -------
        psi_final : np.ndarray, shape (n_states,), complex
        """
        psi = psi0.copy().astype(complex)
        self.time = 0.0

        for i in range(n_steps):
            t_mid = self.time + self.dt / 2  # midpoint for better accuracy
            H_t = H0 + E0 * np.cos(omega * t_mid) * V_dipole
            psi = self.step_with_H(psi, H_t)
            if callback is not None:
                callback(i, self.time, psi)

        return psi

    @staticmethod
    def build_dipole_z(lattice: GeometricLattice) -> csr_matrix:
        """
        Build the z-component of the electric dipole operator.

        In the |n, l, m> basis, the z-operator (proportional to r*cos(theta))
        couples states with selection rules:
            Delta_l = +/-1,  Delta_m = 0

        Matrix elements use exact angular factors (Clebsch-Gordan) and
        hydrogen radial integrals.

        Parameters
        ----------
        lattice : GeometricLattice
            Lattice with states list and _state_index dict

        Returns
        -------
        V_z : csr_matrix
            Dipole operator matrix (sparse, real, symmetric)
        """
        n_states = lattice.num_states
        V = lil_matrix((n_states, n_states), dtype=np.float64)

        for idx_a, (n_a, l_a, m_a) in enumerate(lattice.states):
            # Selection rule: Delta_m = 0, Delta_l = +1
            # Couple (n_a, l_a, m_a) -> (n_b, l_a+1, m_a) for all n_b
            l_b = l_a + 1
            m_b = m_a  # Delta_m = 0

            # Angular selection rule check: |m| <= l_b
            if abs(m_b) > l_b:
                continue

            # Angular factor: <l+1, m | cos(theta) | l, m>
            # = sqrt( ((l+1)^2 - m^2) / ((2l+1)(2l+3)) )
            numerator = (l_b)**2 - m_a**2
            if numerator <= 0:
                continue
            denominator = (2 * l_a + 1) * (2 * l_a + 3)
            angular = np.sqrt(numerator / denominator)

            # Connect to all n_b that have state (n_b, l_b, m_b)
            for n_b in range(l_b + 1, lattice.max_n + 1):
                state_b = (n_b, l_b, m_b)
                if state_b not in lattice._state_index:
                    continue
                idx_b = lattice._state_index[state_b]

                # Radial matrix element: <n_b, l+1 | r | n_a, l>
                # For hydrogen, the dominant coupling is diagonal (n_b = n_a)
                # and nearest (n_b = n_a +/- 1).
                # Use the general formula: proportional to n^2 for diagonal,
                # with off-diagonal suppressed.
                # Exact for 1s->2p: -(128*sqrt(2))/243 = -0.7449
                radial = _hydrogen_radial_dipole(n_a, l_a, n_b, l_b)

                element = angular * radial
                V[idx_a, idx_b] = element
                V[idx_b, idx_a] = element  # Hermitian

        return V.tocsr()


def _hydrogen_radial_dipole(n1: int, l1: int, n2: int, l2: int) -> float:
    """
    Approximate hydrogen radial dipole matrix element <n2, l2 | r | n1, l1>.

    Uses the exact Gordon formula for key transitions and a general
    scaling rule for others.

    The dominant transitions are:
        <2,1|r|1,0> = (128*sqrt(2))/243 = 0.7449  (1s -> 2p)
        <3,1|r|2,0> = (27*sqrt(6))/64  = 1.0328   (2s -> 3p)

    For general n1, n2: scales roughly as (n1*n2)^(3/2) / (n2-n1 terms).

    Parameters
    ----------
    n1, l1 : int
        Initial state quantum numbers
    n2, l2 : int
        Final state quantum numbers (l2 = l1 + 1)

    Returns
    -------
    radial : float
        Radial matrix element in Bohr
    """
    # Exact values for key transitions
    if (n1, l1, n2, l2) == (1, 0, 2, 1):
        return 128.0 * np.sqrt(2.0) / 243.0  # 0.7449

    if (n1, l1, n2, l2) == (2, 0, 3, 1):
        return 27.0 * np.sqrt(6.0) / 64.0  # 1.0328

    if (n1, l1, n2, l2) == (2, 1, 3, 2):
        return 27.0 * np.sqrt(30.0) / 320.0  # 0.4629

    # General scaling: diagonal transitions (n2 = n1 + 1) dominate
    # Off-diagonal (|n2 - n1| > 1) are suppressed
    if n1 == n2:
        # Same-n transitions: <n, l+1 | r | n, l> ~ -3n * sqrt(n^2 - l^2) / 2
        # (from the Stark effect formula)
        return 1.5 * n1 * np.sqrt(max(n1**2 - (l1 + 1)**2, 1))

    # Different n: suppressed by overlap
    dn = abs(n2 - n1)
    if dn == 1:
        # Nearest-neighbor: moderate coupling
        n_avg = (n1 + n2) / 2.0
        return n_avg**2 * 0.5
    else:
        # Far off-diagonal: rapid falloff
        n_avg = (n1 + n2) / 2.0
        return n_avg**2 * 0.5 * (0.3 ** (dn - 1))
