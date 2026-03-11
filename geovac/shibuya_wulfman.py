"""
Shibuya-Wulfman nuclear attraction integrals on S3.

Cross-atom nuclear attraction matrix elements for a hydrogen-like orbital
|n,l,m> on center A in the field of nucleus B at distance R:

    V^SW_{(n'l'm'),(nlm)} = -(Z_B / p0) * D^(n)_{(l'm'),(lm)}(gamma) * f(n, gamma)

where:
    p0      = momentum-shell parameter (energy-shell radius on S3)
    gamma   = bond_angle(R, p0) = arccos((p0^2 - 1/R^2)/(p0^2 + 1/R^2))
    D^(n)   = SO(4) Wigner D-matrix element (opposite-sign rotation)
    f(n, gamma) = shell-dependent form factor encoding R-dependence

The D-matrix is block-diagonal in n (same-shell only). Cross-n coupling
would require the full SO(4,2) conformal group and is not implemented.

Form factor:
    For n=1, the exact Shibuya-Wulfman result gives f(1, gamma) = 1 when
    p0 varies with R self-consistently. Since we use a fixed p0, we need
    the form factor to encode the R-dependence that p0 would otherwise
    carry. We use f(n, gamma) = sin(gamma), which gives:
        - f -> 0 as R -> infinity (gamma -> 0): correct dissociation limit
        - f maximal at intermediate R: correct bonding regime
        - f -> 0 as R -> 0 (gamma -> pi): correct united-atom limit

References:
    Shibuya & Wulfman, Proc. R. Soc. A 286, 376 (1965)
    Aquilanti et al., Chem. Phys. Lett. 318, 619 (2000)
    Paper 8: Bond Sphere Geometry (GeoVac, 2026), Sections IX-XII
"""

from __future__ import annotations

import numpy as np
from typing import Optional, List, Tuple

from .wigner_so4 import wigner_D_so4, bond_angle, d_matrix_block


# ---------------------------------------------------------------------------
#  Form factor
# ---------------------------------------------------------------------------

def sw_form_factor(n: int, gamma: float, sturmian: bool = False) -> float:
    """Shell-dependent form factor f(n, gamma) for SW nuclear attraction.

    Encodes the R-dependence of the cross-atom nuclear attraction when
    using a fixed (non-self-consistent) momentum scale p0.

    f(n, gamma) = sin(gamma)

    This gives the correct limiting behavior:
        - gamma -> 0 (R -> inf): f -> 0  (atoms separate, no attraction)
        - gamma -> pi (R -> 0):  f -> 0  (united-atom limit)
        - gamma ~ pi/2:          f ~ 1   (maximum bonding)

    Note: For a self-consistent p0(R) determination, f(1, gamma) = 1
    would be exact (Shibuya-Wulfman 1965). The sin(gamma) form factor
    compensates for the fixed-p0 approximation.

    In the Sturmian basis (sturmian=True), all functions share a single
    p0 and live on the same S3. The SW cross-center integrals are exact
    D-matrix elements with no form factor (Paper 9, Sec. V).

    Parameters
    ----------
    n : int
        Principal quantum number of the shell.
    gamma : float
        Bond angle in radians (from bond_angle(R, p0)).
    sturmian : bool
        If True, return 1.0 (exact in Sturmian basis). Default False.

    Returns
    -------
    float
        Form factor value in [0, 1].
    """
    if sturmian:
        return 1.0  # exact in Sturmian basis, Paper 9 Sec. V
    return np.sin(gamma)


# ---------------------------------------------------------------------------
#  Single matrix element
# ---------------------------------------------------------------------------

def sw_nuclear_attraction(
    np_: int, lp: int, mp: int,
    n: int, l: int, m: int,
    Z_B: float, p0: float, gamma: float
) -> float:
    """Compute a single Shibuya-Wulfman nuclear attraction matrix element.

    V^SW = -(Z_B / p0) * D^(n)_{(l'm'),(lm)}(gamma) * f(n, gamma)

    Only same-n shells couple (D-matrix is block-diagonal in n).
    Returns 0 for n' != n.

    Parameters
    ----------
    np_ : int
        Principal quantum number of the bra state (n').
    lp : int
        Orbital angular momentum of the bra state (l').
    mp : int
        Magnetic quantum number of the bra state (m').
    n : int
        Principal quantum number of the ket state.
    l : int
        Orbital angular momentum of the ket state.
    m : int
        Magnetic quantum number of the ket state.
    Z_B : float
        Nuclear charge of the attracting center.
    p0 : float
        Momentum-shell parameter.
    gamma : float
        Bond angle in radians.

    Returns
    -------
    float
        Matrix element in Hartree.
    """
    if np_ != n:
        return 0.0

    D_elem = wigner_D_so4(n, lp, mp, l, m, gamma)
    f = sw_form_factor(n, gamma)

    return -(Z_B / p0) * D_elem * f


# ---------------------------------------------------------------------------
#  Full coupling matrix (one direction: orbitals on A, nucleus B)
# ---------------------------------------------------------------------------

def sw_coupling_matrix(
    n_max: int, Z_B: float, p0: float, gamma: float
) -> np.ndarray:
    """Build the full SW nuclear attraction matrix for orbitals on center A
    in the field of nucleus B.

    V^SW_{ij} = -(Z_B / p0) * D^(n)_{(l'_i, m'_i),(l_j, m_j)}(gamma) * f(n, gamma)

    The matrix is n_spatial x n_spatial where n_spatial = sum_{n=1}^{n_max} n^2.
    States are ordered by (n, l, m) with l = 0..n-1, m = -l..l.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.
    Z_B : float
        Nuclear charge of the attracting center.
    p0 : float
        Momentum-shell parameter.
    gamma : float
        Bond angle in radians.

    Returns
    -------
    np.ndarray
        Shape (n_spatial, n_spatial). SW nuclear attraction matrix in Hartree.
    """
    # Build state list: (n, l, m) ordered as in GeometricLattice
    states: List[Tuple[int, int, int]] = []
    for n_val in range(1, n_max + 1):
        for l_val in range(n_val):
            for m_val in range(-l_val, l_val + 1):
                states.append((n_val, l_val, m_val))

    n_spatial = len(states)
    V = np.zeros((n_spatial, n_spatial))

    # Process shell by shell using d_matrix_block for efficiency.
    # The D-matrix is orthogonal (D^T = D^{-1}), not symmetric. But the
    # nuclear attraction V_B = -Z_B/r_B is Hermitian, so V must be symmetric.
    # We symmetrize each block: V_block = scale * (D + D^T) / 2.
    offset = 0
    for n_val in range(1, n_max + 1):
        n_sq = n_val * n_val
        D_block = d_matrix_block(n_val, gamma)  # n^2 x n^2
        D_sym = (D_block + D_block.T) / 2.0     # Hermitian symmetrization
        f = sw_form_factor(n_val, gamma)
        scale = -(Z_B / p0) * f

        V[offset:offset + n_sq, offset:offset + n_sq] = scale * D_sym
        offset += n_sq

    return V


# ---------------------------------------------------------------------------
#  Symmetrized two-center coupling matrix
# ---------------------------------------------------------------------------

def sw_coupling_matrix_AB(
    n_max: int, Z_A: float, Z_B: float,
    p0: float, gamma: float
) -> np.ndarray:
    """Build symmetrized SW coupling matrix for a diatomic A-B.

    The combined spatial basis has 2 * n_spatial orbitals:
        indices 0..n_spatial-1          : orbitals on atom A
        indices n_spatial..2*n_spatial-1: orbitals on atom B

    Cross-atom blocks are symmetrized for heteronuclear molecules:
        V_cross = (V_AB + V_BA) / 2
    where V_AB uses Z_B (A orbitals feel nucleus B) and V_BA uses Z_A
    (B orbitals feel nucleus A). This ensures the combined matrix is
    symmetric (Hermitian) even when Z_A != Z_B.

    Diagonal blocks are zero (same-atom nuclear attraction is already
    in the atomic Hamiltonian).

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.
    Z_A : float
        Nuclear charge of center A.
    Z_B : float
        Nuclear charge of center B.
    p0 : float
        Momentum-shell parameter.
    gamma : float
        Bond angle in radians.

    Returns
    -------
    np.ndarray
        Shape (2*n_spatial, 2*n_spatial). Combined cross-atom coupling
        matrix in Hartree.
    """
    # Cross-nuclear: orbitals on A feel Z_B, orbitals on B feel Z_A
    V_AB = sw_coupling_matrix(n_max, Z_B, p0, gamma)  # A orbs, B nucleus
    V_BA = sw_coupling_matrix(n_max, Z_A, p0, gamma)  # B orbs, A nucleus

    # Symmetrize: average the two contributions so that
    # V[i_A, j_B] = V[j_B, i_A] even when Z_A != Z_B.
    # Both V_AB and V_BA are individually symmetric (Hermitian D-blocks),
    # so V_cross = (V_AB + V_BA)/2 is also symmetric.
    V_cross = (V_AB + V_BA) / 2.0

    n_spatial = V_AB.shape[0]
    n_total = 2 * n_spatial
    V = np.zeros((n_total, n_total))

    V[:n_spatial, n_spatial:] = V_cross
    V[n_spatial:, :n_spatial] = V_cross  # symmetric

    return V
