"""
Hyperfine Geometric Impedance

Computes the impedance mismatch Δκ between electron and nuclear lattices
for hyperfine splitting calculations.

Theory:
-------
The hyperfine splitting arises from magnetic coupling between electron
and nuclear spins. The geometric impedance captures the phase space mismatch:

  Δκ = S_electron / S_nuclear

where:
- S_electron = symplectic capacity of electron lattice
- S_nuclear = symplectic capacity of nuclear lattice (mass-scaled)

The mass scaling accounts for the nuclear recoil:
  S_nuclear_scaled = S_base × (m_proton / m_electron)

This impedance appears in the hyperfine energy formula:
  ΔE_HFS = E₀ × α² × Δκ × g_p × C

Reference:
----------
old_research_archive/src/hyperfine_impedance.py
Papers/Paper_4_Universality.tex (Section 5)

Author: GeoVac Development Team
Date: February 14, 2026
Version: 0.1.0-alpha
Status: Experimental
"""

import numpy as np
from typing import Dict

from ADSCFT.bulk.paraboloid_lattice import ParaboloidLattice
from ADSCFT.bulk.symplectic import compute_shell_capacity, compute_impedance_mismatch


# Physical constants
MASS_RATIO_PROTON_ELECTRON = 1836.15267343
MASS_RATIO_MUON_ELECTRON = 206.768


class HyperfineImpedanceCalculator:
    """
    Calculate geometric impedance for hyperfine splitting.

    Computes the symplectic capacity ratio that determines the
    strength of electron-nuclear magnetic coupling.

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number
    n_shell : int
        Principal quantum number for capacity calculation
    lepton_mass_ratio : float
        Lepton mass / electron mass

    Attributes
    ----------
    lattice : ParaboloidLattice
        Geometric lattice
    S_electron : float
        Electron symplectic capacity
    S_nuclear : float
        Nuclear symplectic capacity (mass-scaled)
    impedance_mismatch : float
        Δκ = S_electron / S_nuclear

    Examples
    --------
    >>> calc = HyperfineImpedanceCalculator(max_n=10, n_shell=1)
    >>> result = calc.compute_impedance()
    >>> print(f"Δκ = {result['impedance_mismatch']:.6e}")
    Δκ = 5.448e-04
    """

    def __init__(
        self,
        max_n: int = 10,
        n_shell: int = 1,
        lepton_mass_ratio: float = 1.0
    ):
        """
        Initialize hyperfine impedance calculator.

        Parameters
        ----------
        max_n : int
            Maximum principal quantum number
        n_shell : int
            Principal quantum number for calculation
        lepton_mass_ratio : float
            Lepton mass / electron mass (1.0 or 206.768)
        """
        self.max_n = max_n
        self.n_shell = n_shell
        self.lepton_mass_ratio = lepton_mass_ratio

        # Create lattice
        self.lattice = ParaboloidLattice(max_n=max_n)

        # Results
        self.S_electron: float = 0.0
        self.S_nuclear: float = 0.0
        self.impedance_mismatch: float = 0.0

    def compute_electron_capacity(self) -> float:
        """
        Compute electron symplectic capacity.

        Returns
        -------
        S_electron : float
            Electron lattice capacity
        """
        # For n=1, extrapolate from n=2 (pole singularity)
        if self.n_shell == 1:
            S_2 = compute_shell_capacity(self.lattice, n=2)
            self.S_electron = S_2 * (1.0 / 2.0)**4
        else:
            self.S_electron = compute_shell_capacity(self.lattice, n=self.n_shell)

        return self.S_electron

    def compute_nuclear_capacity(self) -> float:
        """
        Compute nuclear symplectic capacity (mass-scaled).

        The nuclear lattice has the same geometric structure but
        is mass-scaled by m_p / m_lepton.

        Returns
        -------
        S_nuclear : float
            Nuclear lattice capacity (mass-scaled)
        """
        # Base capacity (same as electron)
        S_base = self.S_electron if self.S_electron > 0 else self.compute_electron_capacity()

        # Mass scaling
        mass_ratio = MASS_RATIO_PROTON_ELECTRON / self.lepton_mass_ratio
        self.S_nuclear = S_base * mass_ratio

        return self.S_nuclear

    def compute_impedance(self) -> Dict[str, float]:
        """
        Compute complete impedance analysis.

        Returns
        -------
        results : dict
            All computed quantities
        """
        # Compute capacities
        S_e = self.compute_electron_capacity()
        S_n = self.compute_nuclear_capacity()

        # Compute impedance
        self.impedance_mismatch = compute_impedance_mismatch(
            S_electron=S_e,
            S_nuclear=S_e,  # Same base geometry
            mass_ratio=MASS_RATIO_PROTON_ELECTRON / self.lepton_mass_ratio
        )

        results = {
            'S_electron': S_e,
            'S_nuclear': S_n,
            'impedance_mismatch': self.impedance_mismatch,
            'n_shell': self.n_shell,
            'lepton_mass_ratio': self.lepton_mass_ratio,
            'mass_scaling': MASS_RATIO_PROTON_ELECTRON / self.lepton_mass_ratio
        }

        return results

    def print_results(self, verbose: bool = True):
        """Print formatted results"""
        results = self.compute_impedance()

        if verbose:
            print("="*70)
            print("HYPERFINE GEOMETRIC IMPEDANCE")
            print("="*70)
            print(f"\nLattice Configuration:")
            print(f"  max_n = {self.max_n}")
            print(f"  n_shell = {self.n_shell}")
            print(f"  Lepton mass ratio = {self.lepton_mass_ratio:.3f}")

            print(f"\nSymplectic Capacities:")
            print(f"  S_electron = {results['S_electron']:.10f}")
            print(f"  S_nuclear  = {results['S_nuclear']:.10f}")
            print(f"  Mass scaling = {results['mass_scaling']:.6f}")

            print(f"\nImpedance Mismatch:")
            print(f"  Δκ = {results['impedance_mismatch']:.6e}")
            print(f"  1/Δκ = {1.0/results['impedance_mismatch']:.6e}")

            print("="*70)

        return results


if __name__ == "__main__":
    import sys
    import io

    # Set UTF-8 encoding for Windows
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

    # Test for electronic hydrogen
    print("\nElectronic Hydrogen (n=1)")
    calc_e = HyperfineImpedanceCalculator(max_n=15, n_shell=1, lepton_mass_ratio=1.0)
    result_e = calc_e.print_results(verbose=True)

    # Test for muonic hydrogen
    print("\n\nMuonic Hydrogen (n=1)")
    calc_mu = HyperfineImpedanceCalculator(
        max_n=15,
        n_shell=1,
        lepton_mass_ratio=MASS_RATIO_MUON_ELECTRON
    )
    result_mu = calc_mu.print_results(verbose=True)

    # Compare
    print("\n" + "="*70)
    print("COMPARISON")
    print("="*70)
    ratio = result_mu['impedance_mismatch'] / result_e['impedance_mismatch']
    print(f"Δκ(μ) / Δκ(e) = {ratio:.6f}")
    print(f"Mass ratio effect: {MASS_RATIO_MUON_ELECTRON:.3f}")
    print("="*70)
