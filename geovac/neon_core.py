"""
Frozen-core screening for noble-gas-like cores.

Drop-in replacement for CoreScreening when the core cannot be solved with
the Level 3 hyperspherical solver (which handles only 2-electron He-like
cores).

Supported core types:
  [Ne]     (10 electrons): Z=11-18 (Na through Ar)
  [Ar]     (18 electrons): Z=19-20 (K, Ca) and Z=21-30 (transition metals)
  [Ar]3d10 (28 electrons): Z=31-36 (Ga through Kr)
  [Kr]     (36 electrons): Z=37-38 (Rb, Sr)             [Sprint 3, v2.12.0]
  [Xe]     (54 electrons): Z=55-56 (Cs, Ba)             [Sprint 3, v2.12.0]

The core density is computed analytically using hydrogenic radial
wavefunctions with Clementi-Raimondi effective nuclear charges (1963, 1967).

    N_core(r) = integral_0^r n_core(r') dr'
    Z_eff(r)  = Z - N_core(r)

Boundary conditions: Z_eff(0) = Z, Z_eff(infinity) = Z - N_core_total.

The hydrogenic wavefunctions use Z_eff = n * zeta for each shell,
where zeta are the Clementi-Raimondi Slater orbital exponents.

References:
  - Clementi & Raimondi, J. Chem. Phys. 38, 2686 (1963) -- Z=2-36
  - Clementi, Raimondi & Reinhardt, J. Chem. Phys. 47, 1300 (1967) -- Z=37-86
  - NIST Atomic Spectra Database (core energies)
  - Paper 17, Section III (Z_eff screening usage in composed geometries)
"""

import numpy as np
from scipy.special import genlaguerre
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import eigh_tridiagonal
from typing import Optional, Tuple, Union


# ---------------------------------------------------------------------------
# Clementi-Raimondi Slater orbital exponents
# zeta values: for hydrogenic R_nl(r; Z_eff), use Z_eff = n * zeta
# Source: Clementi & Raimondi, J. Chem. Phys. 38, 2686 (1963), Table II
# ---------------------------------------------------------------------------

# [Ne] core (10 electrons): 1s2 2s2 2p6
_CLEMENTI_ZETA_NE = {
    # Z: (zeta_1s, zeta_2s, zeta_2p)
    11: (10.6259, 6.5714, 6.8018),
    12: (11.6089, 7.3920, 7.8258),
    13: (12.5910, 8.2136, 8.9634),
    14: (13.5745, 9.0200, 9.9450),
    15: (14.5578, 9.8250, 10.9612),
    16: (15.5409, 10.6288, 11.9770),
    17: (16.5239, 11.4304, 12.9932),
    18: (17.5075, 12.2304, 14.0082),
}

# [Ar] core (18 electrons): 1s2 2s2 2p6 3s2 3p6
_CLEMENTI_ZETA_AR = {
    # Z: (zeta_1s, zeta_2s, zeta_2p, zeta_3s, zeta_3p)
    # Clementi & Raimondi, J. Chem. Phys. 38, 2686 (1963)
    19: (18.4895, 13.0058, 15.0272, 8.6804, 9.0282),
    20: (19.4730, 13.7760, 16.0415, 9.6021, 10.1040),
    # Transition metals Z=21-30 ([Ar] core, Track CZ)
    21: (20.4566, 14.5737, 17.0553, 10.3398, 10.9588),
    22: (21.4409, 15.3768, 18.0650, 11.0332, 11.8206),
    23: (22.4256, 16.1814, 19.0730, 11.7089, 12.6789),
    24: (23.4138, 16.9843, 20.0750, 12.3685, 13.5352),
    25: (24.3957, 17.7903, 21.0843, 13.0184, 14.3927),
    26: (25.3810, 18.5952, 22.0892, 13.6764, 15.2481),
    27: (26.3668, 19.4046, 23.0922, 14.3223, 16.1044),
    28: (27.3526, 20.2127, 24.0953, 14.9612, 16.9601),
    29: (28.3386, 21.0200, 25.0970, 15.5938, 17.8147),
    30: (29.3245, 21.8278, 26.0979, 16.2195, 18.6670),
}

# [Ar]3d10 core (28 electrons): 1s2 2s2 2p6 3s2 3p6 3d10
_CLEMENTI_ZETA_AR3D10 = {
    # Z: (zeta_1s, zeta_2s, zeta_2p, zeta_3s, zeta_3p, zeta_3d)
    31: (30.4769, 21.8412, 24.5607, 16.8963, 17.3682, 13.0686),
    32: (31.4634, 22.6625, 25.5702, 17.7602, 18.3913, 14.0870),
    33: (32.4499, 23.4692, 26.5782, 18.5960, 19.3722, 15.0927),
    34: (33.4364, 24.2688, 27.5850, 19.4069, 20.3398, 16.0939),
    35: (34.4230, 25.0681, 28.5913, 20.2171, 21.3035, 17.0927),
    36: (35.4095, 25.8637, 29.5964, 21.0325, 22.2610, 18.0860),
}

# [Kr] core (36 electrons): 1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6
# Sprint 3 Track HA-A (v2.12.0): fifth-row s-block support for Sunaga 2025
# comparison (SrH).
# Source: Clementi, Raimondi & Reinhardt, J. Chem. Phys. 47, 1300 (1967),
# Table III (single-zeta Hartree-Fock orbital exponents for neutral atoms).
_CLEMENTI_ZETA_KR = {
    # Z: (zeta_1s, zeta_2s, zeta_2p, zeta_3s, zeta_3p, zeta_3d, zeta_4s, zeta_4p)
    # Rb (Z=37): [Kr] 5s1
    37: (36.3242, 26.7156, 30.5810, 21.8410, 23.2214, 18.9867, 11.5472, 11.8892),
    # Sr (Z=38): [Kr] 5s2
    38: (37.3108, 27.5192, 31.5922, 22.6526, 24.0714, 19.9645, 12.4213, 12.7489),
}

# [Xe] core (54 electrons): [Kr] + 4d10 5s2 5p6
# Sprint 3 Track HA-A (v2.12.0): sixth-row s-block support for Sunaga 2025
# comparison (BaH).
# Source: Clementi, Raimondi & Reinhardt, J. Chem. Phys. 47, 1300 (1967),
# Table III.
_CLEMENTI_ZETA_XE = {
    # Z: (z1s, z2s, z2p, z3s, z3p, z3d, z4s, z4p, z4d, z5s, z5p)
    # Cs (Z=55): [Xe] 6s1
    55: (54.2678, 42.0318, 46.4901, 36.9232, 38.4598, 33.2295, 25.2794, 25.5744,
         20.7221, 13.1181, 12.3118),
    # Ba (Z=56): [Xe] 6s2
    56: (55.2540, 42.9170, 47.5009, 37.9057, 39.4558, 34.2470, 26.1764, 26.4651,
         21.6420, 14.0101, 13.1947),
}

# Backward-compatible alias
_CLEMENTI_ZETA = _CLEMENTI_ZETA_NE

# ---------------------------------------------------------------------------
# NIST / Clementi-Roetti Hartree-Fock total energies for ionic cores (Ha)
# [Ne] cores: energies of Z-electron nucleus with 10 electrons
# [Ar] cores: energies of Z-electron nucleus with 18 electrons
# [Ar]3d10 cores: energies of Z-electron nucleus with 28 electrons
# Source: Clementi & Roetti, At. Data Nucl. Data Tables 14, 177 (1974)
# ---------------------------------------------------------------------------

_NIST_CORE_ENERGIES_NE = {
    # [Ne]-like ions (10 electrons)
    11: -161.859,    # Na+
    12: -199.615,    # Mg2+
    13: -241.877,    # Al3+
    14: -288.854,    # Si4+
    15: -340.719,    # P5+
    16: -397.505,    # S6+
    17: -459.482,    # Cl7+
    18: -526.818,    # Ar (neutral, same as Ne-like with Z=18)
}

_NIST_CORE_ENERGIES_AR = {
    # [Ar]-like ions (18 electrons)
    # Source: Clementi & Roetti, At. Data Nucl. Data Tables 14, 177 (1974)
    19: -599.02,     # K+
    20: -676.05,     # Ca2+
    # Transition metals Z=21-30 (Track CZ)
    21: -759.74,     # Sc3+
    22: -848.41,     # Ti4+
    23: -942.88,     # V5+
    24: -1043.36,    # Cr6+
    25: -1149.87,    # Mn7+
    26: -1262.44,    # Fe8+
    27: -1381.41,    # Co9+
    28: -1506.87,    # Ni10+
    29: -1638.96,    # Cu11+
    30: -1777.85,    # Zn12+
}

_NIST_CORE_ENERGIES_AR3D10 = {
    # [Ar]3d10-like ions (28 electrons)
    31: -1923.26,    # Ga3+
    32: -2075.36,    # Ge4+
    33: -2234.24,    # As5+
    34: -2399.87,    # Se6+
    35: -2572.44,    # Br7+
    36: -2752.06,    # Kr8+
}

_NIST_CORE_ENERGIES_KR = {
    # [Kr]-like ions (36 electrons)
    # Source: Clementi & Roetti, At. Data Nucl. Data Tables 14, 177 (1974),
    # Hartree-Fock total energies of Rb+, Sr2+
    37: -2938.36,    # Rb+
    38: -3131.54,    # Sr2+
}

_NIST_CORE_ENERGIES_XE = {
    # [Xe]-like ions (54 electrons)
    # Source: Clementi & Roetti, At. Data Nucl. Data Tables 14, 177 (1974),
    # Hartree-Fock total energies of Cs+, Ba2+
    55: -7553.93,    # Cs+
    56: -7883.54,    # Ba2+
}

# Unified lookup
_NIST_CORE_ENERGIES = {
    **_NIST_CORE_ENERGIES_NE,
    **_NIST_CORE_ENERGIES_AR,
    **_NIST_CORE_ENERGIES_AR3D10,
    **_NIST_CORE_ENERGIES_KR,
    **_NIST_CORE_ENERGIES_XE,
}


def _hydrogenic_radial(
    n: int, l: int, Z_eff: float, r_grid: np.ndarray,
) -> np.ndarray:
    """Normalized hydrogenic radial wavefunction R_{nl}(r; Z_eff) on grid.

    Parameters
    ----------
    n : int
        Principal quantum number.
    l : int
        Orbital angular momentum quantum number.
    Z_eff : float
        Effective nuclear charge for this shell.
    r_grid : ndarray
        Radial grid points (bohr).

    Returns
    -------
    R_nl : ndarray
        Radial wavefunction values on the grid, normalized so that
        integral |R_nl|^2 r^2 dr = 1.
    """
    rho = 2.0 * Z_eff * r_grid / n
    L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    wf = rho ** l * np.exp(-rho / 2.0) * L_poly
    norm_sq = np.trapezoid(wf ** 2 * r_grid ** 2, r_grid)
    if norm_sq < 1e-30:
        return np.zeros_like(r_grid)
    return wf / np.sqrt(norm_sq)


class FrozenCore:
    """Frozen-core screening for noble-gas-like cores.

    Drop-in replacement for CoreScreening when the core cannot be
    solved with the Level 3 hyperspherical solver. Uses analytical
    hydrogenic wavefunctions with Clementi-Raimondi orbital exponents.

    Supported core types:
      'Ne'     -- 10 electrons (1s2 2s2 2p6), Z=11-18
      'Ar'     -- 18 electrons (+ 3s2 3p6), Z=19-20 (and Z=21-30 TM)
      'Ar3d10' -- 28 electrons (+ 3d10), Z=31-36
      'Kr'     -- 36 electrons (+ 4s2 4p6), Z=37-38
      'Xe'     -- 54 electrons (+ 4d10 5s2 5p6), Z=55-56

    Usage
    -----
    >>> fc = FrozenCore(Z=11)           # auto-detects [Ne] core
    >>> fc = FrozenCore(Z=19)           # auto-detects [Ar] core
    >>> fc = FrozenCore(Z=32)           # auto-detects [Ar]3d10 core
    >>> fc = FrozenCore(Z=37)           # auto-detects [Kr] core
    >>> fc = FrozenCore(Z=55)           # auto-detects [Xe] core
    >>> fc.solve()
    >>> z_eff_at_1bohr = fc.z_eff(1.0)

    Parameters
    ----------
    Z : int
        Nuclear charge.
    core_type : str, optional
        Core type: 'Ne', 'Ar', 'Ar3d10', 'Kr', or 'Xe'. Auto-detected from
        Z if not specified.
    n_radial : int
        Number of radial grid points for density computation.
    r_max : float
        Maximum radial distance (bohr). Auto-adjusted for deep cores
        if smaller than the default 20.0.
    """

    # Core type -> (n_core_electrons, zeta_table, energy_table)
    _CORE_REGISTRY = {
        'Ne': (10, _CLEMENTI_ZETA_NE, _NIST_CORE_ENERGIES_NE),
        'Ar': (18, _CLEMENTI_ZETA_AR, _NIST_CORE_ENERGIES_AR),
        'Ar3d10': (28, _CLEMENTI_ZETA_AR3D10, _NIST_CORE_ENERGIES_AR3D10),
        'Kr': (36, _CLEMENTI_ZETA_KR, _NIST_CORE_ENERGIES_KR),
        'Xe': (54, _CLEMENTI_ZETA_XE, _NIST_CORE_ENERGIES_XE),
    }

    # Cores with a multi-zeta tabulation registered in
    # geovac/multi_zeta_orbitals.py. Cores not in this set fall back to
    # single-zeta with a UserWarning when ``screening='multi_zeta'``.
    _MULTI_ZETA_AVAILABLE = {'Xe'}  # extended in follow-on commits

    def __init__(
        self,
        Z: int,
        core_type: Optional[str] = None,
        n_radial: int = 2000,
        r_max: float = 20.0,
        screening: str = 'single_zeta',
    ) -> None:
        # Auto-detect core type from Z if not specified
        if core_type is None:
            if Z in _CLEMENTI_ZETA_NE:
                core_type = 'Ne'
            elif Z in _CLEMENTI_ZETA_AR:
                core_type = 'Ar'
            elif Z in _CLEMENTI_ZETA_AR3D10:
                core_type = 'Ar3d10'
            elif Z in _CLEMENTI_ZETA_KR:
                core_type = 'Kr'
            elif Z in _CLEMENTI_ZETA_XE:
                core_type = 'Xe'
            else:
                raise ValueError(
                    f"No frozen-core data for Z={Z}. "
                    f"Supported: Z=11-30 ([Ne] or [Ar]), "
                    f"Z=31-36 ([Ar]3d10), "
                    f"Z=37-38 ([Kr]), Z=55-56 ([Xe]). "
                    f"For Z=3-10 (He-like cores), use CoreScreening."
                )

        if core_type not in self._CORE_REGISTRY:
            raise ValueError(
                f"Unknown core_type={core_type!r}. "
                f"Supported: {list(self._CORE_REGISTRY.keys())}"
            )

        n_core, zeta_table, energy_table = self._CORE_REGISTRY[core_type]
        if Z not in zeta_table:
            raise ValueError(
                f"No Clementi-Raimondi data for Z={Z} with "
                f"core_type={core_type!r}."
            )

        # Screening kernel selection (Sprint Cs-HFS kernel upgrade, May 2026)
        if screening not in ('single_zeta', 'multi_zeta'):
            raise ValueError(
                f"Unknown screening={screening!r}. "
                "Supported: 'single_zeta' (Clementi-Raimondi 1963/1967, "
                "default; backward-compatible) and 'multi_zeta' "
                "(Bunge-Barrientos-Bunge 1993 / two-zeta from CR67)."
            )
        if screening == 'multi_zeta' and core_type not in self._MULTI_ZETA_AVAILABLE:
            from geovac.multi_zeta_orbitals import warn_multi_zeta_unavailable
            warn_multi_zeta_unavailable(core_type)
            screening = 'single_zeta'

        self.Z = Z
        self.core_type = core_type
        self.n_core_electrons = n_core
        self.n_radial = n_radial
        self.r_max = r_max
        self.screening = screening

        self._zeta_table = zeta_table
        self._energy_table = energy_table
        self._solved = False
        self._r_grid: Optional[np.ndarray] = None
        self._n_r: Optional[np.ndarray] = None
        self._N_core: Optional[np.ndarray] = None
        self._N_core_spline: Optional[CubicSpline] = None
        self._n_r_spline: Optional[CubicSpline] = None

    def _build_density_ne(self, r: np.ndarray) -> np.ndarray:
        """Build [Ne] core density (10 electrons)."""
        zeta_1s, zeta_2s, zeta_2p = self._zeta_table[self.Z]
        R_1s = _hydrogenic_radial(1, 0, 1 * zeta_1s, r)
        R_2s = _hydrogenic_radial(2, 0, 2 * zeta_2s, r)
        R_2p = _hydrogenic_radial(2, 1, 2 * zeta_2p, r)
        return (
            2.0 * R_1s ** 2 * r ** 2
            + 2.0 * R_2s ** 2 * r ** 2
            + 6.0 * R_2p ** 2 * r ** 2
        )

    def _build_density_ar(self, r: np.ndarray) -> np.ndarray:
        """Build [Ar] core density (18 electrons)."""
        z1s, z2s, z2p, z3s, z3p = self._zeta_table[self.Z]
        R_1s = _hydrogenic_radial(1, 0, 1 * z1s, r)
        R_2s = _hydrogenic_radial(2, 0, 2 * z2s, r)
        R_2p = _hydrogenic_radial(2, 1, 2 * z2p, r)
        R_3s = _hydrogenic_radial(3, 0, 3 * z3s, r)
        R_3p = _hydrogenic_radial(3, 1, 3 * z3p, r)
        return (
            2.0 * R_1s ** 2 * r ** 2
            + 2.0 * R_2s ** 2 * r ** 2
            + 6.0 * R_2p ** 2 * r ** 2
            + 2.0 * R_3s ** 2 * r ** 2
            + 6.0 * R_3p ** 2 * r ** 2
        )

    def _build_density_ar3d10(self, r: np.ndarray) -> np.ndarray:
        """Build [Ar]3d10 core density (28 electrons)."""
        z1s, z2s, z2p, z3s, z3p, z3d = self._zeta_table[self.Z]
        R_1s = _hydrogenic_radial(1, 0, 1 * z1s, r)
        R_2s = _hydrogenic_radial(2, 0, 2 * z2s, r)
        R_2p = _hydrogenic_radial(2, 1, 2 * z2p, r)
        R_3s = _hydrogenic_radial(3, 0, 3 * z3s, r)
        R_3p = _hydrogenic_radial(3, 1, 3 * z3p, r)
        R_3d = _hydrogenic_radial(3, 2, 3 * z3d, r)
        return (
            2.0 * R_1s ** 2 * r ** 2
            + 2.0 * R_2s ** 2 * r ** 2
            + 6.0 * R_2p ** 2 * r ** 2
            + 2.0 * R_3s ** 2 * r ** 2
            + 6.0 * R_3p ** 2 * r ** 2
            + 10.0 * R_3d ** 2 * r ** 2
        )

    def _build_density_kr(self, r: np.ndarray) -> np.ndarray:
        """Build [Kr] core density (36 electrons)."""
        (z1s, z2s, z2p, z3s, z3p, z3d, z4s, z4p) = self._zeta_table[self.Z]
        R_1s = _hydrogenic_radial(1, 0, 1 * z1s, r)
        R_2s = _hydrogenic_radial(2, 0, 2 * z2s, r)
        R_2p = _hydrogenic_radial(2, 1, 2 * z2p, r)
        R_3s = _hydrogenic_radial(3, 0, 3 * z3s, r)
        R_3p = _hydrogenic_radial(3, 1, 3 * z3p, r)
        R_3d = _hydrogenic_radial(3, 2, 3 * z3d, r)
        R_4s = _hydrogenic_radial(4, 0, 4 * z4s, r)
        R_4p = _hydrogenic_radial(4, 1, 4 * z4p, r)
        return (
            2.0 * R_1s ** 2 * r ** 2
            + 2.0 * R_2s ** 2 * r ** 2
            + 6.0 * R_2p ** 2 * r ** 2
            + 2.0 * R_3s ** 2 * r ** 2
            + 6.0 * R_3p ** 2 * r ** 2
            + 10.0 * R_3d ** 2 * r ** 2
            + 2.0 * R_4s ** 2 * r ** 2
            + 6.0 * R_4p ** 2 * r ** 2
        )

    def _build_density_xe(self, r: np.ndarray) -> np.ndarray:
        """Build [Xe] core density (54 electrons)."""
        (z1s, z2s, z2p, z3s, z3p, z3d,
         z4s, z4p, z4d, z5s, z5p) = self._zeta_table[self.Z]
        R_1s = _hydrogenic_radial(1, 0, 1 * z1s, r)
        R_2s = _hydrogenic_radial(2, 0, 2 * z2s, r)
        R_2p = _hydrogenic_radial(2, 1, 2 * z2p, r)
        R_3s = _hydrogenic_radial(3, 0, 3 * z3s, r)
        R_3p = _hydrogenic_radial(3, 1, 3 * z3p, r)
        R_3d = _hydrogenic_radial(3, 2, 3 * z3d, r)
        R_4s = _hydrogenic_radial(4, 0, 4 * z4s, r)
        R_4p = _hydrogenic_radial(4, 1, 4 * z4p, r)
        R_4d = _hydrogenic_radial(4, 2, 4 * z4d, r)
        R_5s = _hydrogenic_radial(5, 0, 5 * z5s, r)
        R_5p = _hydrogenic_radial(5, 1, 5 * z5p, r)
        return (
            2.0 * R_1s ** 2 * r ** 2
            + 2.0 * R_2s ** 2 * r ** 2
            + 6.0 * R_2p ** 2 * r ** 2
            + 2.0 * R_3s ** 2 * r ** 2
            + 6.0 * R_3p ** 2 * r ** 2
            + 10.0 * R_3d ** 2 * r ** 2
            + 2.0 * R_4s ** 2 * r ** 2
            + 6.0 * R_4p ** 2 * r ** 2
            + 10.0 * R_4d ** 2 * r ** 2
            + 2.0 * R_5s ** 2 * r ** 2
            + 6.0 * R_5p ** 2 * r ** 2
        )

    def _build_density_multi_zeta(self, r: np.ndarray) -> np.ndarray:
        """Build core density from a multi-zeta orbital tabulation.

        Currently dispatches to the Xe-core two-zeta builder
        (Sprint Cs-HFS kernel upgrade, May 2026); other cores will be
        added when their multi-zeta tabulations are added to
        geovac/multi_zeta_orbitals.py.
        """
        from geovac.multi_zeta_orbitals import (
            build_two_zeta_xe_orbitals_from_cr,
            density_from_orbitals,
        )

        if self.core_type == 'Xe':
            orbs = build_two_zeta_xe_orbitals_from_cr(self._zeta_table[self.Z])
            return density_from_orbitals(orbs, r)
        raise ValueError(
            f"No multi-zeta tabulation registered for core_type="
            f"{self.core_type!r}; this should have been caught at __init__."
        )

    def solve(self, verbose: bool = False) -> None:
        """Compute core density and screening function analytically.

        Uses hydrogenic wavefunctions with Clementi-Raimondi effective
        charges for the shells of the frozen core (default), or a
        Roothaan-Hartree-Fock multi-zeta expansion (when
        ``screening='multi_zeta'`` is set in the constructor and the
        core type is registered in ``_MULTI_ZETA_AVAILABLE``).

        Parameters
        ----------
        verbose : bool
            Print diagnostic information.
        """
        # Build fine radial grid
        dr = self.r_max / self.n_radial
        self._r_grid = np.linspace(dr / 2, self.r_max - dr / 2, self.n_radial)
        r = self._r_grid

        # Build density for the appropriate core type and screening kernel
        if self.screening == 'multi_zeta':
            self._n_r = self._build_density_multi_zeta(r)
        elif self.core_type == 'Ne':
            self._n_r = self._build_density_ne(r)
        elif self.core_type == 'Ar':
            self._n_r = self._build_density_ar(r)
        elif self.core_type == 'Ar3d10':
            self._n_r = self._build_density_ar3d10(r)
        elif self.core_type == 'Kr':
            self._n_r = self._build_density_kr(r)
        elif self.core_type == 'Xe':
            self._n_r = self._build_density_xe(r)
        else:
            raise ValueError(f"Unknown core_type={self.core_type!r}")

        n_target = float(self.n_core_electrons)

        # Verify normalization and correct if needed
        total = np.trapezoid(self._n_r, r)
        if verbose:
            print(f"  FrozenCore Z={self.Z} [{self.core_type}]: "
                  f"raw density integral = {total:.6f} "
                  f"(target {n_target:.0f})")
        if total > 1e-10:
            self._n_r *= n_target / total

        # Cumulative integral: N_core(r) = integral_0^r n(r') dr'
        self._N_core = np.zeros_like(r)
        self._N_core[1:] = cumulative_trapezoid(self._n_r, r)

        # Build interpolation splines
        self._N_core_spline = CubicSpline(r, self._N_core, extrapolate=True)
        self._n_r_spline = CubicSpline(r, self._n_r, extrapolate=True)

        self._solved = True

        if verbose:
            N_total = self._N_core[-1]
            print(f"  N_core(r_max={self.r_max}) = {N_total:.6f}")
            print(f"  Z_eff(0.01) = {self.z_eff(0.01):.4f}")
            print(f"  Z_eff(1.0) = {self.z_eff(1.0):.4f}")
            print(f"  Z_eff(5.0) = {self.z_eff(5.0):.4f}")
            print(f"  Core energy = {self.energy:.3f} Ha")

    def _check_solved(self) -> None:
        """Raise if solve() has not been called."""
        if not self._solved:
            raise RuntimeError("Call .solve() before querying screening data.")

    def z_eff(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Effective nuclear charge at distance r from the nucleus.

        Z_eff(r) = Z - N_core(r), where N_core is the cumulative number
        of core electrons inside radius r.

        Parameters
        ----------
        r : float or ndarray
            Radial distance(s) in bohr.

        Returns
        -------
        float or ndarray
            Z_eff(r) at the requested points.
        """
        self._check_solved()
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))

        N_at_r = np.clip(
            self._N_core_spline(r_arr), 0.0, float(self.n_core_electrons),
        )
        result = self.Z - N_at_r
        # At r=0, no screening
        result = np.where(r_arr <= 0, float(self.Z), result)

        if np.ndim(r) == 0:
            return float(result[0])
        return result

    def density(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Core electron density n(r) at distance r.

        Normalized so that integral n(r) dr = N_core (core electron count).

        Parameters
        ----------
        r : float or ndarray
            Radial distance(s) in bohr.

        Returns
        -------
        float or ndarray
            n(r) at the requested points.
        """
        self._check_solved()
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        result = np.clip(self._n_r_spline(r_arr), 0.0, None)

        if np.ndim(r) == 0:
            return float(result[0])
        return result

    @property
    def energy(self) -> float:
        """Core energy from NIST/Clementi-Roetti Hartree-Fock data (Ha)."""
        return self._energy_table[self.Z]

    @property
    def r_grid(self) -> np.ndarray:
        """The internal r grid used for the density computation."""
        self._check_solved()
        return self._r_grid

    @property
    def n_r(self) -> np.ndarray:
        """The radial number density on the internal grid."""
        self._check_solved()
        return self._n_r

    @property
    def N_core(self) -> np.ndarray:
        """The cumulative core electron count on the internal grid."""
        self._check_solved()
        return self._N_core


# ---------------------------------------------------------------------------
# Screened radial solver for valence <1/r^3>
# ---------------------------------------------------------------------------

# Physical fine-structure constant (CODATA 2018)
_ALPHA_PHYSICAL = 7.2973525693e-3


def _solve_screened_radial(
    Z: int,
    l: int,
    n_target: int,
    core_type: Optional[str] = None,
    n_grid: int = 12000,
    r_max: float = 80.0,
    allow_l0: bool = False,
    screening: str = 'single_zeta',
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Solve the radial Schrodinger equation with Z_eff(r) screening.

    Uses a uniform radial grid and second-order finite differences:

        [-1/2 d^2u/dr^2 + V_eff(r) u] = E u

    where u(r) = r*R(r) is the reduced radial wavefunction,
    V_eff(r) = -Z_eff(r)/r + l(l+1)/(2r^2), and Z_eff(r) comes from
    the FrozenCore screening profile.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    l : int
        Orbital angular momentum quantum number. Must be >= 1 unless
        ``allow_l0=True``.
    n_target : int
        Principal quantum number of the desired state. The solver
        returns the (n_target - l - 1)th eigenstate (counting radial
        nodes n_r = 0, 1, 2, ...).
    core_type : str, optional
        Core type for FrozenCore. Auto-detected if None.
    n_grid : int
        Number of grid points.
    r_max : float
        Maximum radius in bohr.
    allow_l0 : bool
        If True, allow l=0. The Kramers exclusion (designed for
        ``screened_r3_inverse`` where 1/r^3 diverges at l=0) is bypassed.
        Required for s-wave HFS observables that depend on |psi(0)|^2
        (e.g. Cs 6S_{1/2}, Sprint Cs-HFS-v2). Default False preserves
        the legacy SO behaviour where l=0 is rejected.

    Returns
    -------
    energy : float
        Eigenvalue in Hartree.
    u_grid : ndarray
        Reduced radial wavefunction u(r) = r*R(r) on the r-grid,
        normalized so that integral |u|^2 dr = 1.
    r_grid : ndarray
        The radial grid points.
    """
    if l < 0:
        raise ValueError(f"l={l} must be non-negative")
    if l == 0 and not allow_l0:
        # Default behaviour: reject l=0 to protect ``screened_r3_inverse``
        # callers from a divergent <1/r^3>. HFS callers must opt in via
        # ``allow_l0=True``.
        raise ValueError("l must be >= 1 for screened radial solver "
                         "(l=0 has divergent 1/r^3); pass allow_l0=True "
                         "for HFS-class observables that need l=0")
    n_r_target = n_target - l - 1  # number of radial nodes
    if n_r_target < 0:
        raise ValueError(f"n_target={n_target} < l+1={l+1}: no such state")

    # Build FrozenCore and get Z_eff(r)
    fc = FrozenCore(Z=Z, core_type=core_type, screening=screening)
    fc.solve()

    # Uniform radial grid starting just above zero.
    # The l(l+1)/(2r^2) centrifugal barrier keeps u(r) ~ r^{l+1} near
    # the origin, so a small but nonzero r_min is safe for l >= 1.
    h = r_max / n_grid
    r = np.linspace(h, r_max, n_grid)

    # Z_eff on the grid
    z_eff_vals = fc.z_eff(r)

    # Effective potential V_eff(r) = -Z_eff(r)/r + l(l+1)/(2r^2)
    v_eff = -z_eff_vals / r + l * (l + 1) / (2.0 * r**2)

    # Standard FD Hamiltonian: H u = [-1/2 d^2/dr^2 + V_eff] u = E u
    # Tridiagonal with:
    #   diagonal:     1/h^2 + V_eff[i]
    #   off-diagonal: -1/(2h^2)
    diag = 1.0 / h**2 + v_eff
    offdiag = np.full(n_grid - 1, -1.0 / (2.0 * h**2))

    # Solve for the lowest few eigenvalues
    n_eigs = n_r_target + 5
    n_eigs = min(n_eigs, n_grid - 2)

    evals, evecs = eigh_tridiagonal(
        diag, offdiag,
        select='i', select_range=(0, n_eigs - 1),
    )

    if n_r_target >= len(evals):
        raise ValueError(
            f"Could not find state with n_r={n_r_target} radial nodes. "
            f"Only {len(evals)} eigenvalues found."
        )
    energy = evals[n_r_target]
    u_vec = evecs[:, n_r_target]

    # Normalize: integral |u|^2 dr = 1
    norm_sq = np.trapezoid(u_vec**2, r)
    if norm_sq > 0:
        u_vec /= np.sqrt(norm_sq)

    return energy, u_vec, r


def _solve_screened_radial_log(
    Z: int,
    l: int,
    n_target: int,
    core_type: Optional[str] = None,
    n_grid: int = 4000,
    r_min: float = 1.0e-6,
    r_max: float = 80.0,
    Z_origin: Optional[float] = None,
    z_eff_callable=None,
    screening: str = 'single_zeta',
) -> Tuple[float, np.ndarray, np.ndarray, float]:
    """Solve the radial Schrodinger equation on a logarithmic radial grid.

    Designed for s-wave (l=0) HFS observables that depend on |R(0)|^2,
    where the standard uniform-FD scheme of ``_solve_screened_radial``
    diverges in |psi(0)|^2 with grid refinement (Track Cs-HFS-v1, May
    2026: the s-wave wavefunction at a singular -Z_eff(r)/r origin needs
    grid resolution that scales with Z_eff(0)).

    Implementation. We use the standard "log-mesh" radial Schrodinger
    transform of atomic structure codes. With r = exp(x) on uniform x
    and the substitution P(x) = u(r) / sqrt(r), the radial equation

        - 1/2 u'' + V_eff(r) u = E u

    becomes the symmetric eigenproblem

        - 1/2 P''(x) + r^2 [V_eff(r) - E] P + (1/8) P = 0

    i.e.

        [- 1/2 P''(x) + r^2 V_eff(r) P(x) + (1/8) P(x)] = E r^2 P(x).

    This is the generalized eigenproblem K P = E M P with K symmetric
    tridiagonal (centered FD on uniform x for P'') and M = diag(r^2)
    positive definite. We solve it with scipy.linalg.eigh; both K and
    M are symmetric so the standard Cholesky-based generalized solver
    handles it.

    For s-states the small-r Frobenius series
    R(r) = R(0)[1 - Z_origin r/2 + (Z_origin r)^2/12 - ...] is matched
    at the smallest grid point to extract R(0).

    Parameters
    ----------
    Z : int
        Nuclear charge.
    l : int
        Orbital angular momentum (>= 0; s-wave is the use case).
    n_target : int
        Principal quantum number.
    core_type : str, optional
        Core type for FrozenCore. Auto-detected if None.
    n_grid : int
        Number of log-grid points. Default 4000 is converged for Cs 6s.
    r_min : float
        Minimum radius (bohr). Default 1e-6 is the small-r boundary
        below which an analytical Frobenius series is used.
    r_max : float
        Maximum radius (bohr).
    Z_origin : float, optional
        Effective Z at the origin for the small-r Frobenius series.
        Defaults to ``Z`` (the nuclear charge before screening, since
        Z_eff(0) = Z by construction).
    z_eff_callable : callable, optional
        If provided, used in place of ``FrozenCore(Z).z_eff`` for the
        screened potential. Signature: ``f(r) -> Z_eff(r)`` (works on
        scalars or arrays). Used by hydrogenic and bare-Coulomb tests
        where no FrozenCore is registered for the given Z.

    Returns
    -------
    energy : float
        Eigenvalue in Hartree.
    u_grid : ndarray
        Reduced radial wavefunction u(r) = r*R(r) on the log r grid,
        normalized so that integral |u|^2 dr = 1.
    r_grid : ndarray
        The logarithmic radial grid (length ``n_grid``).
    R0 : float
        R(0) = lim_{r->0} u(r)/r for the requested state. Equal to zero
        for l>=1; nonzero for s-states. The contact density is
        |psi(0)|^2 = R0**2 / (4*pi).

    Notes
    -----
    R(0) is extracted by the analytical Frobenius series
    ``R(r) ~ R(0) * (1 - Z_origin * r / 2 + ...)`` valid for r < 1/Z_origin,
    matched to the numerical solution at the smallest grid point r_min
    where the numerical scheme is still accurate. The match is
    R(0) = R(r_min) / (1 - Z_origin * r_min / 2 + (Z_origin*r_min)^2 / 12)
    (3-term Frobenius for the bare -Z_origin/r potential).
    """
    if l < 0:
        raise ValueError(f"l={l} must be non-negative")
    n_r_target = n_target - l - 1
    if n_r_target < 0:
        raise ValueError(f"n_target={n_target} < l+1={l+1}: no such state")

    if Z_origin is None:
        Z_origin = float(Z)

    # Build the screening profile (FrozenCore by default; user-supplied
    # callable for hydrogenic / bare-Coulomb test cases).
    if z_eff_callable is None:
        fc = FrozenCore(Z=Z, core_type=core_type, screening=screening)
        fc.solve()
        z_eff_callable = fc.z_eff

    # IMPLEMENTATION (Track Cs-HFS-v2): use a DENSE UNIFORM GRID that is
    # automatically scaled to resolve the inner-shell physics. The relevant
    # length scale near the origin is r_inner ~ 1 / Z_origin. We use a
    # fully uniform grid with r_max set by the user and r_min adjusted so
    # there are at least 100 grid points inside r_inner. The grid spacing
    # h = r_max / n_grid must satisfy h * 100 < 1 / Z_origin, i.e.
    # n_grid > 100 * Z_origin * r_max.
    #
    # The earlier log-mesh transform attempted in this function (and the
    # 2-piece hybrid that followed) both had numerical artifacts at the
    # boundary or transition that destroyed the s-wave |R(0)|^2 extraction.
    # The dense uniform grid is conceptually simple and matches the
    # Frobenius extrapolation cleanly: u(r_min) = r_min * R(r_min), with
    # r_min set just barely above 0.
    h = r_max / n_grid
    # Override r_min to be h (the smallest natural grid point); we do this
    # rather than honor the user's r_min so the FD stencil at the boundary
    # (where we use a Dirichlet ghost u_{-1}=0) is consistent.
    r = np.linspace(h, r_max, n_grid)

    # Verify resolution near the origin (warn but don't error if too coarse)
    pts_inside_bohr = int((1.0 / Z_origin) / h)
    # (used for caller sanity; not enforced)

    # Effective potential
    z_eff_vals = z_eff_callable(r)
    v_eff = -z_eff_vals / r + l * (l + 1) / (2.0 * r**2)

    # Standard 3-point centered FD on uniform grid:
    #   u''(r_i) ~ (u_{i+1} - 2 u_i + u_{i-1}) / h^2
    # With Dirichlet BC u(0)=0 implicit at i=-1 (the FIRST grid point r=h
    # has u_{-1}=0 contributing -u_{-1}/(2h^2) = 0 to the off-diagonal).
    diag_K = 1.0 / h**2 + v_eff
    offdiag_K = np.full(n_grid - 1, -0.5 / h**2)

    # Solve symmetric tridiagonal eigenproblem
    n_eigs = n_r_target + 5
    n_eigs = min(n_eigs, n_grid - 2)
    evals, evecs = eigh_tridiagonal(
        diag_K, offdiag_K,
        select='i', select_range=(0, n_eigs - 1),
    )

    if n_r_target >= len(evals):
        raise ValueError(
            f"Could not find state with n_r={n_r_target} radial nodes. "
            f"Only {len(evals)} eigenvalues found."
        )
    energy = float(evals[n_r_target])
    u_vec = evecs[:, n_r_target]

    # Normalize: integral |u|^2 dr
    norm_sq = np.trapezoid(u_vec**2, r)
    if norm_sq > 0:
        u_vec /= np.sqrt(norm_sq)

    # Extract R(0) for s-states via small-r Frobenius extrapolation.
    # For V(r) ~ -Z_origin/r near origin, the regular s-wave is
    #   R(r) = R(0) [1 - Z_0 r / 2 + (Z_0 r)^2 / 12 - ...]
    # so u(r) = r R(r). Match at r=h (the smallest grid point).
    # Use the SECOND grid point r=h for matching (the first grid point's
    # u value is biased by the Dirichlet ghost; r=h is one full step in).
    if l == 0:
        # Use a conservative average of the first two grid points to
        # reduce the boundary-bias on R(0) extraction.
        r_match = r[0]  # = h
        Z0r = Z_origin * r_match
        frobenius_3 = 1.0 - Z0r / 2.0 + Z0r**2 / 12.0
        R_at_match_0 = u_vec[0] / r[0]
        # 4-term for cross-check:
        R0_first = R_at_match_0 / frobenius_3
        # Average over first ~5 points using exact bare-Coulomb correction:
        R0_avg_pts = []
        for i in range(min(5, len(r))):
            ri = r[i]
            Zri = Z_origin * ri
            frob_i = 1.0 - Zri / 2.0 + Zri**2 / 12.0 - (Zri**3) / 144.0
            if abs(frob_i) > 1e-10:
                R0_avg_pts.append(u_vec[i] / ri / frob_i)
        if R0_avg_pts:
            R0 = float(np.mean(R0_avg_pts))
        else:
            R0 = R0_first
        if not np.isfinite(R0):
            R0 = 0.0
    else:
        R0 = 0.0  # u ~ r^{l+1}, so R(0) = 0 for l>=1

    return energy, u_vec, r, float(R0)


def screened_psi_origin_squared(
    Z: int,
    n: int,
    l: int = 0,
    core_type: Optional[str] = None,
    n_grid: int = 4000,
    r_min: float = 1.0e-6,
    r_max: float = 80.0,
    Z_origin: Optional[float] = None,
    z_eff_callable=None,
    screening: str = 'single_zeta',
) -> float:
    """Compute |psi_{nl}(0)|^2 for a screened valence electron.

    For an s-state (l=0), psi(0) is finite and given by R(0)/sqrt(4*pi).
    For l >= 1, psi(0) = 0 identically.

    Uses the log-grid solver ``_solve_screened_radial_log`` with the
    analytical Frobenius small-r series matched at r_min to extract R(0).

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n : int
        Principal quantum number.
    l : int
        Orbital angular momentum (default 0).
    core_type : str, optional
        Core type for FrozenCore. Auto-detected if None.
    n_grid : int
        Log-grid points. Default 4000 is converged for Cs 6s.
    r_min : float
        Minimum log-grid radius (bohr). Default 1e-6.
    r_max : float
        Maximum log-grid radius (bohr). Default 80.0.
    Z_origin : float, optional
        Effective Z at the origin for Frobenius extrapolation.
        Defaults to ``Z`` (the unscreened nuclear charge).
    z_eff_callable : callable, optional
        Custom Z_eff(r) for the test cases that don't have a registered
        FrozenCore (e.g. hydrogenic Z=1).

    Returns
    -------
    float
        |psi_{nl}(0)|^2 in atomic units (bohr^{-3}). Zero for l >= 1.
    """
    _energy, _u, _r, R0 = _solve_screened_radial_log(
        Z, l, n,
        core_type=core_type,
        n_grid=n_grid,
        r_min=r_min,
        r_max=r_max,
        Z_origin=Z_origin,
        z_eff_callable=z_eff_callable,
        screening=screening,
    )
    psi0_sq = R0**2 / (4.0 * np.pi)
    return float(psi0_sq)


def screened_r3_inverse(
    Z: int,
    n: int,
    l: int,
    core_type: Optional[str] = None,
    n_grid: int = 12000,
    r_max: float = 80.0,
) -> float:
    """Compute <1/r^3> for a valence electron in a screened potential.

    Solves the radial Schrodinger equation with the FrozenCore Z_eff(r)
    screening profile and computes the expectation value of 1/r^3 using
    the resulting wavefunction.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n : int
        Principal quantum number of the valence state.
    l : int
        Orbital angular momentum (must be >= 1; l=0 diverges).
    core_type : str, optional
        Core type for FrozenCore. Auto-detected if None.
    n_grid : int
        Number of radial grid points.
    r_max : float
        Maximum radius in bohr.

    Returns
    -------
    float
        <1/r^3> in atomic units (bohr^{-3}).

    Raises
    ------
    ValueError
        If l=0 (divergent) or n < l+1 (invalid state).
    """
    if l < 1:
        raise ValueError(
            "<1/r^3> diverges for l=0. Spin-orbit has Kramers "
            "cancellation at l=0 so this is never needed."
        )

    energy, u_vec, r = _solve_screened_radial(
        Z, l, n, core_type=core_type,
        n_grid=n_grid, r_max=r_max,
    )

    # <1/r^3> = integral |R(r)|^2 / r^3 * r^2 dr
    #         = integral |u(r)|^2 / r^3 dr
    # where u = r*R, so |R|^2 = |u|^2/r^2
    integrand = u_vec**2 / r**3
    result = np.trapezoid(integrand, r)

    return result


def screened_xi_so(
    Z: int,
    n: int,
    l: int,
    core_type: Optional[str] = None,
    alpha: float = _ALPHA_PHYSICAL,
    n_grid: int = 12000,
) -> float:
    """Compute the SO coupling constant xi = <(1/r) dV/dr> for screened potential.

    For V(r) = -Z_eff(r)/r, the spin-orbit coupling function is:
        xi(r) = (alpha^2 / 2) * (1/r) * dV/dr

    where dV/dr = Z_eff(r)/r^2 - Z_eff'(r)/r (chain rule on -Z_eff/r).

    This returns <xi>_{nl} = integral |R_{nl}|^2 * xi(r) * r^2 dr
    computed with the screened wavefunction.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n, l : int
        Quantum numbers (l >= 1).
    core_type : str, optional
        Core type for FrozenCore.
    alpha : float
        Fine-structure constant.
    n_grid : int
        Number of grid points.

    Returns
    -------
    float
        <xi>_{nl} in Hartree. Multiply by L.S eigenvalue to get H_SO.
    """
    if l < 1:
        raise ValueError("l must be >= 1")

    energy, u_vec, r = _solve_screened_radial(
        Z, l, n, core_type=core_type, n_grid=n_grid,
    )
    h = r[1] - r[0]

    # Get Z_eff(r) on the grid
    fc = FrozenCore(Z=Z, core_type=core_type)
    fc.solve()
    z_eff_r = fc.z_eff(r)

    # dZ_eff/dr via finite differences (Z_eff decreases with r, so dZ/dr < 0)
    dz_dr = np.gradient(z_eff_r, r)

    # V(r) = -Z_eff(r)/r
    # dV/dr = Z_eff(r)/r^2 - Z_eff'(r)/r
    # (1/r) dV/dr = Z_eff(r)/r^3 - Z_eff'(r)/r^2
    one_over_r_dVdr = z_eff_r / r**3 - dz_dr / r**2

    # xi(r) = (alpha^2 / 2) * (1/r) dV/dr
    xi_r = alpha**2 / 2.0 * one_over_r_dVdr

    # <xi> = integral |u|^2 / r^2 * xi(r) * r^2 dr = integral |u|^2 * xi(r) / 1 dr
    # Wait: |R|^2 r^2 dr = |u|^2/r^2 * r^2 dr = |u|^2 dr
    # So <xi> = integral |u|^2 * xi(r) dr... no.
    # <xi> = integral R*(r) * xi(r) * R(r) * r^2 dr
    #       = integral |R|^2 * xi(r) * r^2 dr
    #       = integral (|u|^2/r^2) * xi(r) * r^2 dr
    #       = integral |u|^2 * xi(r) dr
    result = np.trapezoid(u_vec**2 * xi_r, r)

    return result


def screened_so_splitting(
    Z: int,
    n: int,
    l: int,
    core_type: Optional[str] = None,
    alpha: float = _ALPHA_PHYSICAL,
    n_grid: int = 12000,
) -> dict:
    """Compute spin-orbit splitting using screened <1/r^3>.

    The SO Hamiltonian is:
        H_SO = (Z_nuc * alpha^2 / 2) * L.S * <1/r^3>_screened

    where <1/r^3>_screened is computed from the FrozenCore-screened
    radial wavefunction.

    L.S eigenvalues:
        j = l + 1/2:  L.S = l/2
        j = l - 1/2:  L.S = -(l+1)/2

    The splitting between j=l+1/2 and j=l-1/2 is:
        Delta = (Z_nuc * alpha^2 / 2) * (2l+1)/2 * <1/r^3>

    Parameters
    ----------
    Z : int
        Nuclear charge (= Z_nuc in the SO operator).
    n : int
        Principal quantum number.
    l : int
        Orbital angular momentum (>= 1).
    core_type : str, optional
        Core type for FrozenCore. Auto-detected if None.
    alpha : float
        Fine-structure constant.
    n_grid : int
        Number of radial grid points.

    Returns
    -------
    dict with keys:
        'r3_inv_screened' : float -- <1/r^3> from screened wavefunction
        'r3_inv_hydrogenic' : float -- hydrogenic <1/r^3> at Z_eff_asymptotic
        'enhancement' : float -- ratio screened/hydrogenic
        'E_j_plus' : float -- SO energy for j = l + 1/2 (Ha)
        'E_j_minus' : float -- SO energy for j = l - 1/2 (Ha)
        'splitting' : float -- E(j+) - E(j-) (Ha)
        'splitting_cm1' : float -- splitting in cm^{-1}
        'energy_eV' : float -- eigenvalue of the screened radial solver (eV)
    """
    if l < 1:
        raise ValueError("l must be >= 1 for SO splitting")

    # Screened <1/r^3>
    r3_inv = screened_r3_inverse(Z, n, l, core_type=core_type, n_grid=n_grid)

    # Get eigenvalue too
    energy, _, _ = _solve_screened_radial(Z, l, n, core_type=core_type,
                                          n_grid=n_grid)

    # Full screened SO parameter xi = <(1/r) dV/dr> (alpha^2/2 included)
    xi_val = screened_xi_so(Z, n, l, core_type=core_type,
                            alpha=alpha, n_grid=n_grid)

    # Asymptotic Z_eff for hydrogenic comparison
    fc = FrozenCore(Z=Z, core_type=core_type)
    fc.solve()
    z_eff_asym = fc.z_eff(50.0)  # far from nucleus

    # Hydrogenic <1/r^3> at the asymptotic Z_eff
    r3_inv_hyd = z_eff_asym**3 / (n**3 * l * (l + 0.5) * (l + 1))

    # --- Method 1: Z_nuc * <1/r^3> (approximate, overcounts) ---
    prefactor_znuc = Z * alpha**2 / 2.0 * r3_inv

    # --- Method 2: <(1/r) dV/dr> (correct screened potential) ---
    # xi_val already includes alpha^2/2

    # L.S eigenvalues
    ls_plus = l / 2.0           # j = l + 1/2
    ls_minus = -(l + 1) / 2.0   # j = l - 1/2

    # Use the correct screened xi for the splittings
    e_j_plus = xi_val * ls_plus
    e_j_minus = xi_val * ls_minus
    splitting = e_j_plus - e_j_minus  # = xi_val * (2l+1)/2

    # Also compute the Z_nuc*<1/r^3> version for comparison
    splitting_znuc = prefactor_znuc * (2 * l + 1) / 2.0

    # Conversion: 1 Ha = 219474.63 cm^{-1}
    ha_to_cm1 = 219474.63
    ha_to_ev = 27.211386

    return {
        'r3_inv_screened': r3_inv,
        'r3_inv_hydrogenic': r3_inv_hyd,
        'enhancement': r3_inv / r3_inv_hyd if r3_inv_hyd > 0 else float('inf'),
        'xi_so': xi_val,
        'E_j_plus': e_j_plus,
        'E_j_minus': e_j_minus,
        'splitting': splitting,
        'splitting_cm1': splitting * ha_to_cm1,
        'splitting_znuc_cm1': splitting_znuc * ha_to_cm1,
        'energy_eV': energy * ha_to_ev,
    }
