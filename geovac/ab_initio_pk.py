"""
Ab initio Phillips-Kleinman pseudopotential from core wavefunctions.

Derives the Gaussian PK barrier V_PK(r) = A * exp(-B*r^2) / r^2 entirely
from quantities available in the solved CoreScreening object:

    B = 1 / (2 * r_inv2^2)       -- width from <1/r^2>-weighted core radius
    A = |E_core/N_core - E_val| * N_core  -- height from core-valence energy gap

The effective radius r_inv2 = sqrt(<n(r)>/<n(r)/r^2>) weights the inner core
-- where the PK barrier must act -- rather than the diffuse tail that dominates
the RMS radius.  This gives B values 5-15x larger than the RMS formula and
produces LiH R_eq within ~6% of experiment.

No experimental molecular data is used.  Only atomic properties: Z, the core
wavefunction (via CoreScreening), and the first ionization energy of the
neutral atom (for the valence energy estimate).

References:
  - Phillips & Kleinman, Phys. Rev. 116, 287 (1959)
  - Paper 15, Sec. V (pseudopotential in Level 4 solver)
"""

import numpy as np
from typing import Optional, Dict

from geovac.core_screening import CoreScreening


# First ionization energies of neutral atoms (Hartree), used as valence
# energy estimates.  These are ATOMIC properties, not molecular ones.
# Source: NIST Atomic Spectra Database.
_IONIZATION_ENERGIES: Dict[int, float] = {
    1: 0.5,        # H  -> H+
    2: 0.9036,     # He -> He+
    3: 0.1981,     # Li -> Li+
    4: 0.3426,     # Be -> Be+
    5: 0.3049,     # B  -> B+
    6: 0.4138,     # C  -> C+
    7: 0.5341,     # N  -> N+
    8: 0.5003,     # O  -> O+
    9: 0.6404,     # F  -> F+
    10: 0.7925,    # Ne -> Ne+
    11: 0.1889,    # Na -> Na+
    12: 0.2810,    # Mg -> Mg+
}


class AbInitioPK:
    """
    Derive Phillips-Kleinman pseudopotential parameters from core solution.

    Zero fitted parameters -- everything from the core wavefunction and
    atomic ionization energies.

    The PK pseudopotential enforces Pauli exclusion between valence and core
    electrons by creating a repulsive barrier where the core density is large.
    The Gaussian/r^2 form V = A * exp(-B*r^2) / r^2 is a standard choice that:
    - Diverges as 1/r^2 near the nucleus (matching the centrifugal barrier)
    - Decays exponentially beyond the core radius

    Parameters
    ----------
    core : CoreScreening
        Solved core screening object (must have .solve() called).
    n_core : int
        Number of core electrons (default 2).
    E_val_override : float or None
        Override the valence energy estimate. If None, uses the first
        ionization energy of the neutral atom with Z_neutral = Z - n_core + n_val
        where n_val is the number of valence electrons in the neutral atom.

    Attributes
    ----------
    A : float
        PK barrier height (Ha * bohr^2).
    B : float
        PK barrier width exponent (1/bohr^2).
    r_core : float
        RMS radius of the core electron density (bohr).
    E_core_per_electron : float
        Core energy per electron (Ha).
    E_val_est : float
        Estimated valence orbital energy (Ha, negative).
    """

    def __init__(
        self,
        core: CoreScreening,
        n_core: int = 2,
        E_val_override: Optional[float] = None,
        B_method: str = 'inv2',
    ) -> None:
        core._check_solved()
        self._core = core
        self._n_core = n_core
        self._Z = core.Z
        self._B_method = B_method

        # --- Compute all B candidate radii ---
        self._B_candidates = self._compute_B_candidates()

        # --- Select B from the chosen method ---
        if B_method not in self._B_candidates:
            raise ValueError(
                f"B_method must be one of {list(self._B_candidates.keys())},"
                f" got '{B_method}'"
            )
        chosen = self._B_candidates[B_method]
        self._r_core = chosen['r']
        self._B = chosen['B']

        # --- Valence energy estimate ---
        E_core_total = core.energy
        self._E_core_per_electron = E_core_total / n_core

        if E_val_override is not None:
            self._E_val_est = E_val_override
        else:
            # Use ionization energy of the neutral atom
            # For Li (Z=3, n_core=2): neutral atom is Li, IE = 0.198 Ha
            # For Be (Z=4, n_core=2): neutral atom is Be, IE = 0.343 Ha
            # The valence orbital energy is approximately -IE
            Z_int = int(round(self._Z))
            if Z_int in _IONIZATION_ENERGIES:
                self._E_val_est = -_IONIZATION_ENERGIES[Z_int]
            else:
                # Fallback: hydrogenic estimate for valence n=2
                Z_eff_val = self._Z - n_core
                self._E_val_est = -Z_eff_val**2 / 8.0  # -Z_eff^2/(2*n^2), n=2

        # --- A from core-valence energy gap ---
        # The PK barrier height ~ |E_core/N_core - E_val| * N_core
        # This is the kinetic energy cost of orthogonalizing the valence
        # wavefunction to the core.
        self._A = abs(self._E_core_per_electron - self._E_val_est) * self._n_core

    def _compute_B_candidates(self) -> Dict[str, Dict[str, float]]:
        """
        Compute B from several core radius definitions.

        Returns a dict mapping method name -> {'r': radius, 'B': 1/(2*r^2)}.
        """
        from scipy.interpolate import CubicSpline
        from scipy.optimize import brentq

        r = self._core.r_grid
        n_r = self._core.n_r
        N_core = self._core.N_core

        candidates: Dict[str, Dict[str, float]] = {}

        # --- RMS radius (original formula) ---
        r2_avg = np.trapezoid(r**2 * n_r, r) / np.trapezoid(n_r, r)
        r_rms = float(np.sqrt(r2_avg))
        candidates['rms'] = {'r': r_rms, 'B': 1.0 / (2.0 * r_rms**2)}

        # --- Peak radius: where r^2 * n(r) is maximum ---
        # This is the most probable radial distance for a core electron.
        radial_density = r**2 * n_r
        r_peak = float(r[np.argmax(radial_density)])
        if r_peak > 1e-10:
            candidates['peak'] = {'r': r_peak, 'B': 1.0 / (2.0 * r_peak**2)}

        # --- Median radius: where N_core(r) = N_total/2 ---
        # For a 2-electron core, this is where N_core = 1.0.
        N_half = self._n_core / 2.0
        N_max = N_core[-1]
        if N_max > N_half:
            N_spline = CubicSpline(r, N_core)
            # Find bracketing interval
            idx_cross = np.searchsorted(N_core, N_half)
            r_lo = r[max(idx_cross - 1, 1)]
            r_hi = r[min(idx_cross + 1, len(r) - 1)]
            r_median = float(brentq(lambda x: N_spline(x) - N_half, r_lo, r_hi))
            candidates['median'] = {
                'r': r_median, 'B': 1.0 / (2.0 * r_median**2),
            }

        # --- Expectation of 1/r^2 weighted by density ---
        # r_eff = sqrt( integral n(r) dr / integral n(r)/r^2 dr )
        # Weights the inner core more heavily than RMS.
        mask = r > 0.01  # avoid 1/r^2 singularity
        num = np.trapezoid(n_r[mask], r[mask])
        den = np.trapezoid(n_r[mask] / r[mask]**2, r[mask])
        if den > 0:
            r_inv2 = float(np.sqrt(num / den))
            candidates['inv2'] = {
                'r': r_inv2, 'B': 1.0 / (2.0 * r_inv2**2),
            }

        return candidates

    @property
    def B_candidates(self) -> Dict[str, Dict[str, float]]:
        """All computed B candidates with their radii and B values."""
        return self._B_candidates

    @property
    def B_method(self) -> str:
        """The method used for the B parameter."""
        return self._B_method

    @property
    def A(self) -> float:
        """PK barrier height (Ha * bohr^2)."""
        return self._A

    @property
    def B(self) -> float:
        """PK barrier width exponent (1/bohr^2)."""
        return self._B

    @property
    def r_core(self) -> float:
        """Effective core radius used for B (bohr)."""
        return self._r_core

    @property
    def E_core_per_electron(self) -> float:
        """Core energy per electron (Ha)."""
        return self._E_core_per_electron

    @property
    def E_val_est(self) -> float:
        """Estimated valence orbital energy (Ha, negative)."""
        return self._E_val_est

    @property
    def Z(self) -> float:
        """Nuclear charge."""
        return self._Z

    @property
    def n_core(self) -> int:
        """Number of core electrons."""
        return self._n_core

    def V_pk(self, r: np.ndarray) -> np.ndarray:
        """
        Evaluate the PK pseudopotential at distance r.

        V_PK(r) = A * exp(-B * r^2) / r^2

        Parameters
        ----------
        r : ndarray
            Radial distances (bohr).  Values <= 0 return 0.

        Returns
        -------
        ndarray
            PK potential values (Ha).
        """
        r = np.asarray(r, dtype=float)
        result = np.zeros_like(r)
        safe = r > 1e-15
        result[safe] = self._A * np.exp(-self._B * r[safe]**2) / r[safe]**2
        return result

    def pk_dict(self, atom: str = 'A') -> dict:
        """
        Return PK parameters as a dict compatible with the Level 4 solver.

        Parameters
        ----------
        atom : str
            Which nucleus the core belongs to ('A' or 'B').

        Returns
        -------
        dict
            {'C_core': A, 'beta_core': B, 'atom': atom}
        """
        return {
            'C_core': self._A,
            'beta_core': self._B,
            'atom': atom,
        }

    def summary(self) -> str:
        """Return a human-readable summary of the ab initio PK parameters."""
        method_desc = {
            'rms': 'RMS core radius',
            'peak': 'peak of r^2*n(r)',
            'median': 'median core radius (N=1)',
            'inv2': '<1/r^2>-weighted radius',
        }
        desc = method_desc.get(self._B_method, self._B_method)
        lines = [
            f"Ab initio PK for Z={self._Z:.0f} ({self._n_core} core electrons):",
            f"  r_core  = {self._r_core:.4f} bohr  ({desc})",
            f"  B       = {self._B:.4f} bohr^-2  (= 1/(2*r_core^2),"
            f" method='{self._B_method}')",
            f"  E_core  = {self._core.energy:.6f} Ha  (total core energy)",
            f"  E_core/e= {self._E_core_per_electron:.4f} Ha  (per electron)",
            f"  E_val   = {self._E_val_est:.4f} Ha  (valence estimate)",
            f"  A       = {self._A:.4f} Ha*bohr^2"
            f"  (= |{self._E_core_per_electron:.3f} - ({self._E_val_est:.3f})|"
            f" * {self._n_core})",
        ]
        # Show all candidates
        lines.append("  B candidates:")
        for name, info in self._B_candidates.items():
            marker = " <-- selected" if name == self._B_method else ""
            lines.append(
                f"    {name:8s}: r={info['r']:.4f}, B={info['B']:.4f}{marker}"
            )
        return "\n".join(lines)
