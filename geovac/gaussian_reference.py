"""
Gaussian-Basis Reference Hamiltonians
======================================

Provides genuine Gaussian-basis molecular Hamiltonians from published
integral values for benchmarking against GeoVac lattice encodings.

Each system function returns a dict with:
    h1              : np.ndarray (M, M)     — one-electron integrals (spatial MO basis)
    eri             : np.ndarray (M,M,M,M)  — two-electron integrals, CHEMIST notation (pq|rs)
    nuclear_repulsion : float               — V_NN
    n_electrons     : int
    n_spatial       : int                   — number of spatial orbitals M
    description     : str
    literature_energy : float               — reference FCI energy (Ha)
    source          : str                   — literature citation

The ERIs are in chemist notation:
    (pq|rs) = integral phi_p(1) phi_q(1) (1/r12) phi_r(2) phi_s(2) d1 d2

Compatible with qubit_encoding.build_fermion_op_from_integrals().

Author: GeoVac Development Team
Date: March 2026
"""

from pathlib import Path
from typing import Dict, Any, Tuple

import numpy as np

from openfermion import jordan_wigner

from geovac.qubit_encoding import build_fermion_op_from_integrals


# ---------------------------------------------------------------------------
# H2 STO-3G at R = 1.4 bohr
# ---------------------------------------------------------------------------

def h2_sto3g(R: float = 1.4) -> Dict[str, Any]:
    """
    H2 in STO-3G minimal basis at bond length R (bohr).

    Returns MO-basis integrals from Szabo & Ostlund, *Modern Quantum
    Chemistry* (Dover, 1996), Chapter 3 worked example.  The two spatial
    MOs are sigma_g (bonding) and sigma_u (antibonding).

    At R = 1.4 bohr the FCI energy (exact in this basis) is -1.1373 Ha.

    Parameters
    ----------
    R : float
        Bond length in bohr.  Published integrals are for R = 1.4;
        other values are NOT supported (integrals are R-dependent).

    Returns
    -------
    dict
        Keys: h1, eri, nuclear_repulsion, n_electrons, n_spatial,
              description, literature_energy, source.

    Raises
    ------
    ValueError
        If R != 1.4 (integrals only published at this geometry).

    References
    ----------
    Szabo & Ostlund, *Modern Quantum Chemistry*, Dover 1996, Table 3.15.
    MO integrals after restricted Hartree-Fock in minimal basis.
    """
    if abs(R - 1.4) > 1e-6:
        raise ValueError(
            f"STO-3G H2 integrals are only available at R=1.4 bohr, got R={R}. "
            "Other geometries require an integral engine (e.g. PySCF)."
        )

    n_spatial = 2
    nuclear_repulsion = 1.0 / R  # Z_A * Z_B / R = 1/1.4

    # One-electron integrals in canonical MO basis (Szabo & Ostlund Table 3.15)
    h1 = np.array([
        [-1.2528, 0.0],
        [0.0, -0.4756],
    ])

    # Two-electron integrals in chemist notation (pq|rs), MO basis
    # Only unique nonzero elements listed; full 8-fold permutational
    # symmetry enforced.
    eri = np.zeros((n_spatial, n_spatial, n_spatial, n_spatial))

    # Coulomb integrals
    eri[0, 0, 0, 0] = 0.6746   # (σg σg | σg σg)
    eri[1, 1, 1, 1] = 0.6975   # (σu σu | σu σu)
    eri[0, 0, 1, 1] = 0.6632   # (σg σg | σu σu) = J12
    eri[1, 1, 0, 0] = 0.6632   # symmetry

    # Exchange integrals
    eri[0, 1, 0, 1] = 0.1813   # (σg σu | σg σu) = K12
    eri[1, 0, 1, 0] = 0.1813   # symmetry
    eri[0, 1, 1, 0] = 0.1813   # (σg σu | σu σg)
    eri[1, 0, 0, 1] = 0.1813   # symmetry

    return {
        'h1': h1,
        'eri': eri,
        'nuclear_repulsion': nuclear_repulsion,
        'n_electrons': 2,
        'n_spatial': n_spatial,
        'description': 'H2 STO-3G at R=1.4 bohr (Szabo & Ostlund)',
        'literature_energy': -1.1373,
        'source': (
            'Szabo & Ostlund, Modern Quantum Chemistry, Dover 1996, '
            'Table 3.15 (MO integrals after RHF)'
        ),
    }


# ---------------------------------------------------------------------------
# He STO-3G
# ---------------------------------------------------------------------------

def he_sto3g() -> Dict[str, Any]:
    """
    Helium atom in STO-3G basis (1 spatial orbital, 2 spin-orbitals).

    Uses analytical Slater-type orbital integrals with the STO-3G
    exponent zeta = 1.6875 (Hehre, Stewart & Pople, JCP 51, 2657, 1969).
    For a single-zeta basis, the STO integrals are known exactly:

        h11 = T + V_ne = zeta^2/2 - Z*zeta
        (11|11) = 5*zeta/8

    These are the exact STO values; the 3-Gaussian fit introduces
    errors < 0.001 Ha that are negligible for benchmarking purposes.

    With only 1 spatial orbital the FCI is identical to Hartree-Fock
    (no correlation — only one possible 2-electron determinant).

    Returns
    -------
    dict
        Keys: h1, eri, nuclear_repulsion, n_electrons, n_spatial,
              description, literature_energy, source.

    References
    ----------
    Hehre, Stewart & Pople, JCP 51, 2657 (1969) — STO-3G basis.
    Slater integral formulas: Szabo & Ostlund, Appendix A.
    """
    Z = 2
    zeta = 1.6875  # STO-3G exponent for He

    # Analytical STO integrals
    h11 = zeta**2 / 2.0 - Z * zeta   # T + V_ne = -1.951171875
    j11 = 5.0 * zeta / 8.0           # (11|11) = 1.054687500

    n_spatial = 1

    h1 = np.array([[h11]])

    eri = np.zeros((n_spatial, n_spatial, n_spatial, n_spatial))
    eri[0, 0, 0, 0] = j11

    # E = 2*h11 + (11|11) = -2.847656250 Ha
    literature_energy = 2.0 * h11 + j11

    return {
        'h1': h1,
        'eri': eri,
        'nuclear_repulsion': 0.0,
        'n_electrons': 2,
        'n_spatial': n_spatial,
        'description': 'He STO-3G (zeta=1.6875, analytical STO integrals)',
        'literature_energy': literature_energy,
        'source': (
            'Analytical STO integrals with zeta=1.6875 from '
            'Hehre, Stewart & Pople, JCP 51, 2657 (1969). '
            'Formulas: Szabo & Ostlund, Appendix A.'
        ),
    }


# ---------------------------------------------------------------------------
# H2 6-31G at R = 1.4 bohr — STUB
# ---------------------------------------------------------------------------

def h2_631g(R: float = 1.4) -> Dict[str, Any]:
    """
    H2 in 6-31G basis at bond length R (bohr).  4 spatial orbitals,
    8 qubits after Jordan-Wigner.

    .. warning::
        NOT YET IMPLEMENTED.  Requires MO integrals from an integral
        engine (PySCF, Psi4, or Gaussian).  PySCF is unavailable on
        the current Windows build environment.

    TODO
    ----
    Fill in h1 (4x4) and eri (4x4x4x4) from PySCF output::

        from pyscf import gto, scf, ao2mo
        mol = gto.M(atom='H 0 0 0; H 0 0 1.4', basis='6-31g', unit='Bohr')
        mf = mol.RHF().run()
        h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
        eri = ao2mo.full(mol, mf.mo_coeff)  # chemist notation

    The 6-31G H2 FCI energy at R=1.4 bohr is approximately -1.1515 Ha.
    (Sherrill & Schaefer, Adv Quantum Chem 34, 1999.)

    Parameters
    ----------
    R : float
        Bond length in bohr.

    Raises
    ------
    NotImplementedError
        Always — integrals not yet available.
    """
    raise NotImplementedError(
        "H2 6-31G integrals require PySCF or equivalent integral engine. "
        "Run the PySCF snippet in this function's docstring to generate "
        "the h1 (4x4) and eri (4x4x4x4) arrays, then hardcode them here."
    )


# ---------------------------------------------------------------------------
# He cc-pVDZ — computed MO integrals
# ---------------------------------------------------------------------------
# Integrals computed from cc-pVDZ basis set [2s1p] (Woon & Dunning, 1994)
# using a single-center Gaussian integral engine with real spherical
# harmonics, Slater R^k radial integrals, and real Gaunt coefficients.
# See debug/compute_he_gaussian_integrals.py for the integral engine.
#
# FCI energy: -2.8875948311 Ha (exact within basis for 2 electrons)
# Published reference: -2.8877 Ha (Woon & Dunning, JCP 100, 2975, 1994)
# Difference: 0.0001 Ha — within numerical integration tolerance.

def he_cc_pvdz() -> Dict[str, Any]:
    """
    Helium atom in cc-pVDZ basis (5 spatial orbitals, 10 spin-orbitals).

    Basis: cc-pVDZ [2s1p] from Woon & Dunning, JCP 100, 2975 (1994).
    MO basis: eigenstates of the core Hamiltonian (T + V_ne).
    Orbital ordering: 1s, 2s, 2p_x, 2p_y, 2p_z.

    FCI energy (exact within basis): -2.8876 Ha.
    156 Pauli terms after Jordan-Wigner encoding.

    Returns
    -------
    dict
        Keys: h1, eri, nuclear_repulsion, n_electrons, n_spatial,
              description, literature_energy, source.
    """
    n_spatial = 5

    # One-electron integrals in MO basis (diagonal — eigenstates of h_core)
    h1 = np.diag([-1.99362334, -0.03762721, 0.78499729, 0.78499729, 0.78499729])

    # Two-electron integrals in chemist notation (pq|rs), MO basis
    # Computed from cc-pVDZ Gaussians via Slater R^k expansion with
    # real Gaunt coefficients. Only unique nonzero elements listed.
    eri = np.zeros((n_spatial, n_spatial, n_spatial, n_spatial))

    # s-s block
    eri[0, 0, 0, 0] = 1.24534987
    eri[0, 0, 1, 1] = 0.81704098;  eri[1, 1, 0, 0] = 0.81704098
    eri[1, 1, 1, 1] = 0.62811238
    eri[0, 0, 0, 1] = 0.33726197;  eri[0, 1, 0, 0] = 0.33726197
    eri[0, 0, 1, 0] = 0.33726197;  eri[1, 0, 0, 0] = 0.33726197
    eri[0, 1, 0, 1] = 0.18695371;  eri[1, 0, 1, 0] = 0.18695371
    eri[0, 1, 1, 0] = 0.18695371;  eri[1, 0, 0, 1] = 0.18695371
    eri[1, 1, 0, 1] = 0.16134384;  eri[0, 1, 1, 1] = 0.16134384
    eri[1, 1, 1, 0] = 0.16134384;  eri[1, 0, 1, 1] = 0.16134384

    # s-p Coulomb (same for all p components by spherical symmetry)
    for p in range(2, 5):
        eri[0, 0, p, p] = 1.03082782;  eri[p, p, 0, 0] = 1.03082782
        eri[1, 1, p, p] = 0.70538742;  eri[p, p, 1, 1] = 0.70538742
        eri[0, 1, p, p] = 0.18794441;  eri[p, p, 0, 1] = 0.18794441
        eri[1, 0, p, p] = 0.18794441;  eri[p, p, 1, 0] = 0.18794441

    # s-p exchange
    for p in range(2, 5):
        eri[0, p, 0, p] = 0.18720434;  eri[p, 0, p, 0] = 0.18720434
        eri[0, p, p, 0] = 0.18720434;  eri[p, 0, 0, p] = 0.18720434
        eri[1, p, 1, p] = 0.02447355;  eri[p, 1, p, 1] = 0.02447355
        eri[1, p, p, 1] = 0.02447355;  eri[p, 1, 1, p] = 0.02447355
        # Cross s-p exchange
        eri[0, p, 1, p] = -0.00099483;  eri[p, 0, p, 1] = -0.00099483
        eri[0, p, p, 1] = -0.00099483;  eri[p, 0, 1, p] = -0.00099483
        eri[1, p, 0, p] = -0.00099483;  eri[p, 1, p, 0] = -0.00099483
        eri[1, p, p, 0] = -0.00099483;  eri[p, 1, 0, p] = -0.00099483

    # p-p block: Coulomb and exchange (spherical symmetry)
    for p in range(2, 5):
        eri[p, p, p, p] = 1.04053090
        for q in range(2, 5):
            if q != p:
                eri[p, p, q, q] = 0.91311895
                eri[p, q, p, q] = 0.06370597
                eri[p, q, q, p] = 0.06370597

    literature_energy = -2.8875948311

    return {
        'h1': h1,
        'eri': eri,
        'nuclear_repulsion': 0.0,
        'n_electrons': 2,
        'n_spatial': n_spatial,
        'description': 'He cc-pVDZ [2s1p], 5 spatial orbitals (computed)',
        'literature_energy': literature_energy,
        'source': (
            'MO integrals computed from cc-pVDZ basis (Woon & Dunning, '
            'JCP 100, 2975, 1994) using single-center Gaussian integral '
            'engine with real Gaunt coefficients. '
            'See debug/compute_he_gaussian_integrals.py.'
        ),
    }


# ---------------------------------------------------------------------------
# He cc-pVTZ — cached MO integrals
# ---------------------------------------------------------------------------

def he_cc_pvtz() -> Dict[str, Any]:
    """
    Helium atom in cc-pVTZ basis (14 spatial orbitals, 28 spin-orbitals).

    Basis: cc-pVTZ [3s2p1d] from Woon & Dunning, JCP 100, 2975 (1994).
    MO integrals loaded from cache (computed by the single-center
    Gaussian integral engine).

    FCI energy (exact within basis): -2.9002 Ha.
    21,607 Pauli terms after Jordan-Wigner encoding.

    Returns
    -------
    dict
        Keys: h1, eri, nuclear_repulsion, n_electrons, n_spatial,
              description, literature_energy, source.

    Raises
    ------
    FileNotFoundError
        If the cached integral file is not found. Run
        debug/compute_he_gaussian_integrals.py to regenerate.
    """
    cache_path = Path(__file__).parent / 'cache' / 'he_cc_pvtz_mo_integrals.npz'
    if not cache_path.exists():
        raise FileNotFoundError(
            f"cc-pVTZ integral cache not found at {cache_path}. "
            "Run debug/compute_he_gaussian_integrals.py to generate it."
        )

    data = np.load(cache_path)
    h1 = data['h1_mo']
    eri = data['eri_mo']
    literature_energy = float(data['fci_energy'])

    return {
        'h1': h1,
        'eri': eri,
        'nuclear_repulsion': 0.0,
        'n_electrons': 2,
        'n_spatial': 14,
        'description': 'He cc-pVTZ [3s2p1d], 14 spatial orbitals (computed)',
        'literature_energy': literature_energy,
        'source': (
            'MO integrals computed from cc-pVTZ basis (Woon & Dunning, '
            'JCP 100, 2975, 1994) using single-center Gaussian integral '
            'engine with real Gaunt coefficients. '
            'See debug/compute_he_gaussian_integrals.py.'
        ),
    }


# Gaussian scaling fit (H2 data from Paper 14) — kept for reference
_GAUSS_ALPHA = 4.591
_GAUSS_C = 0.02087


def _estimate_gaussian_pauli_terms(n_qubits: int) -> int:
    """Estimate Pauli term count for a dense-ERI Gaussian Hamiltonian."""
    return int(round(_GAUSS_C * n_qubits ** _GAUSS_ALPHA))


# ---------------------------------------------------------------------------
# Convenience: build OpenFermion operators from any reference system
# ---------------------------------------------------------------------------

def build_qubit_hamiltonian(
    system: Dict[str, Any],
) -> Tuple[Any, Any, int]:
    """
    Build FermionOperator and QubitOperator from a reference system dict.

    Parameters
    ----------
    system : dict
        Output of h2_sto3g(), he_sto3g(), etc.

    Returns
    -------
    fermion_op : FermionOperator
    qubit_op : QubitOperator
    n_pauli_terms : int
    """
    fermion_op = build_fermion_op_from_integrals(
        system['h1'],
        system['eri'],
        system['nuclear_repulsion'],
    )
    qubit_op = jordan_wigner(fermion_op)
    n_pauli_terms = len(qubit_op.terms)
    return fermion_op, qubit_op, n_pauli_terms
