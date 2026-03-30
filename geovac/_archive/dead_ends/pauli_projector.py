"""
Pauli projector for core-valence separation in the Level 4 solver.

Removes core orbital contamination from the valence angular basis by
projecting out the component of each angular channel that overlaps
with the 1s² core.  The core 1s orbital is approximated as a hydrogenic
1s with effective Z (Clementi-Raimondi Z_eff ≈ 2.69 for Li).

The projector operates in the angular channel space at fixed hyperradius
R_e and internuclear distance R:

    P = I - |v⟩⟨v| / ⟨v|v⟩

where v_ν = ⟨φ_core(r₁)φ_core(r₂) | Φ_ν(α, θ₁, θ₂; ρ)⟩ is the overlap
of the core 1s² state with angular channel ν.

References:
  - Paper 15, Sec. V (multichannel expansion)
  - core_screening.py (CoreScreening class)
"""

import numpy as np
from typing import Optional, Tuple, List
from math import sqrt, factorial

from geovac.level4_multichannel import _channel_list, _channel_list_extended


class CoreValenceProjector:
    """
    Project out core 1s² contamination from Level 4 angular Hamiltonian.

    The core 1s² wavefunction in hyperspherical coordinates is:

        ψ_core(r₁, r₂) = φ_{1s}(r₁) φ_{1s}(r₂)

    where φ_{1s}(r) = 2 Z_eff^{3/2} exp(-Z_eff r) / √(4π).

    In the (α, l₁, l₂) basis of the Level 4 solver, only the
    (l₁=0, l₂=0) channel has nonzero overlap with the s-wave core.

    Parameters
    ----------
    Z_eff_core : float
        Effective nuclear charge for the core 1s orbital.
        Default 2.69 (Clementi-Raimondi for Li).
    l_max : int
        Maximum angular momentum in the Level 4 channel expansion.
    n_alpha : int
        Number of alpha grid points.
    homonuclear : bool
        Whether the molecule is homonuclear (affects channel list).
    m_max : int
        Maximum |m| per electron.
    l_max_per_m : dict or None
        Per-|m| angular momentum limit.
    """

    def __init__(
        self,
        Z_eff_core: float = 2.69,
        l_max: int = 2,
        n_alpha: int = 100,
        homonuclear: bool = False,
        m_max: int = 0,
        l_max_per_m: Optional[dict] = None,
    ) -> None:
        self.Z_eff_core = Z_eff_core
        self.l_max = l_max
        self.n_alpha = n_alpha
        self.homonuclear = homonuclear
        self.m_max = m_max
        self.l_max_per_m = l_max_per_m

        # Build channel list
        if m_max == 0:
            channels_2 = _channel_list(l_max, homonuclear=homonuclear)
            self.channels_4 = [(l1, 0, l2, 0) for l1, l2 in channels_2]
        else:
            self.channels_4 = _channel_list_extended(
                l_max, m_max, l_max_per_m, homonuclear=homonuclear,
            )
        self.n_ch = len(self.channels_4)

        # Build alpha grid (same convention as Level 4 solver)
        h = (np.pi / 2) / (n_alpha + 1)
        self.alpha_grid = (np.arange(n_alpha) + 1) * h
        self.h_alpha = h

        self.N = self.n_ch * n_alpha

    def _core_overlap_vector(self, R_e: float) -> np.ndarray:
        """
        Compute the overlap of the 1s² core with each angular basis function.

        The hydrogenic 1s radial function is:
            R_{1s}(r) = 2 Z^{3/2} exp(-Z r)

        In hyperspherical coords: r₁ = R_e cos α, r₂ = R_e sin α.
        The overlap with the (l₁=0, l₂=0, m₁=0, m₂=0) channel at each
        alpha point is:

            v(α) = R_{1s}(R_e cos α) × R_{1s}(R_e sin α)
                 = 4 Z³ exp(-Z R_e (cos α + sin α))

        For channels with l₁ ≠ 0 or l₂ ≠ 0, the overlap is zero because
        the 1s orbital is purely s-wave.

        Parameters
        ----------
        R_e : float
            Electronic hyperradius (bohr).

        Returns
        -------
        v : ndarray of shape (N,)
            Overlap vector in the full angular basis.
        """
        Z = self.Z_eff_core
        v = np.zeros(self.N)

        cos_a = np.cos(self.alpha_grid)
        sin_a = np.sin(self.alpha_grid)

        # Core 1s² profile on the alpha grid
        # R_{1s}(r) = 2 Z^{3/2} exp(-Z r), unnormalized in the angular sense
        core_profile = 4.0 * Z**3 * np.exp(
            -Z * R_e * (cos_a + sin_a)
        )

        # Only the (0, 0, 0, 0) channel has overlap
        for ic, (l1, m1, l2, m2) in enumerate(self.channels_4):
            if l1 == 0 and l2 == 0 and m1 == 0 and m2 == 0:
                start = ic * self.n_alpha
                v[start:start + self.n_alpha] = core_profile
                break

        return v

    def build_projector(self, R_e: float) -> np.ndarray:
        """
        Build the Pauli projection matrix P = I - |v⟩⟨v| / ⟨v|v⟩.

        P projects OUT the core component, leaving only the valence space.

        Parameters
        ----------
        R_e : float
            Electronic hyperradius (bohr).

        Returns
        -------
        P : ndarray of shape (N, N)
            Projection matrix with properties: P² = P, P = Pᵀ,
            rank = N - 1, eigenvalues ∈ {0, 1}.
        """
        v = self._core_overlap_vector(R_e)
        v_norm_sq = np.dot(v, v)

        if v_norm_sq < 1e-30:
            return np.eye(self.N)

        # P = I - |v⟩⟨v| / ⟨v|v⟩
        P = np.eye(self.N) - np.outer(v, v) / v_norm_sq
        return P

    def project(self, H_angular: np.ndarray, R_e: float) -> np.ndarray:
        """
        Apply the Pauli projector to the angular Hamiltonian.

        Returns P @ H @ P where P removes the core overlap.

        Parameters
        ----------
        H_angular : ndarray of shape (N, N)
            Angular Hamiltonian matrix.
        R_e : float
            Electronic hyperradius (bohr).

        Returns
        -------
        H_projected : ndarray of shape (N, N)
            Projected Hamiltonian with core contamination removed.
        """
        P = self.build_projector(R_e)
        return P @ H_angular @ P

    @classmethod
    def from_core_screening(
        cls,
        core_screening,
        l_max: int = 2,
        n_alpha: int = 100,
        homonuclear: bool = False,
        m_max: int = 0,
        l_max_per_m: Optional[dict] = None,
    ) -> "CoreValenceProjector":
        """
        Construct from a CoreScreening object.

        Uses the nuclear charge Z from the CoreScreening to set
        Z_eff_core ≈ Z - 0.31 (Clementi-Raimondi 1s screening for Li-like).

        Parameters
        ----------
        core_screening : CoreScreening
            Solved core screening object.
        """
        # Clementi-Raimondi: for 1s in Li, screening constant ≈ 0.31
        # Z_eff = Z - 0.31 for Z >= 3
        Z = core_screening.Z
        if Z >= 3:
            Z_eff_core = Z - 0.31
        else:
            Z_eff_core = Z

        return cls(
            Z_eff_core=Z_eff_core, l_max=l_max, n_alpha=n_alpha,
            homonuclear=homonuclear, m_max=m_max,
            l_max_per_m=l_max_per_m,
        )
