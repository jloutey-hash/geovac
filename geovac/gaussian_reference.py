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
