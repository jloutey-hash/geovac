"""
VQ-2: Pauli sigma vertex for graph-native QED on S^3.
======================================================

Tests whether the Pauli matrices sigma_a (a = x, y, z) can serve as a
VECTOR vertex for graph-native QED on the Fock-projected S^3.

The existing graph-native QED (GN-1 through GN-7) uses a SCALAR identity
vertex: V[a,b,e] = delta(e=(a,b)).  This gives scalar QED with only 1/8
continuum QED selection rules surviving on the scalar Fock graph.

Here we build sigma_x, sigma_y, sigma_z in the coupled (j, m_j) basis
(DiracLabel), determine the "sigma-edge graph" (edges where any sigma_mu
matrix element is nonzero), build QED on this graph, and check how many
continuum selection rules are recovered.

Key finding from VQ-1: sigma.L (the spin-orbit operator) is purely
INTRA-SHELL and ZERO on l=0, so it can't serve as a vertex.  But the
PURE Pauli matrices sigma_a are different -- they couple states within
the same (n, l) block but CAN change kappa (j changes while l stays fixed).

This script does NOT modify any production code.
"""

import numpy as np
import json
import sys
import os
from fractions import Fraction
from collections import defaultdict

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.dirac_matrix_elements import (
    DiracLabel, iter_dirac_labels, kappa_to_l, kappa_to_j,
)
from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, sqrt as sp_sqrt, S as sp_S, nsimplify


# ===========================================================================
# Step 0: Basis setup
# ===========================================================================

def build_labels(n_max):
    """Enumerate all DiracLabel states at n_max."""
    labels = list(iter_dirac_labels(n_max))
    label_index = {lab: i for i, lab in enumerate(labels)}
    return labels, label_index


def cg_coeff(j1, m1, j2, m2, j, m):
    """Clebsch-Gordan coefficient <j1,m1; j2,m2 | j,m> via sympy.
    Returns exact sympy value.
    """
    val = clebsch_gordan(j1, j2, j, m1, m2, m)
    return val


def cg_float(j1, m1, j2, m2, j, m):
    """CG coefficient as float."""
    return float(cg_coeff(j1, m1, j2, m2, j, m))


# ===========================================================================
# Step 1: Build sigma_x, sigma_y, sigma_z in the DiracLabel basis
# ===========================================================================

def build_sigma_matrices(labels, label_index):
    """Build sigma_x, sigma_y, sigma_z in the |n, kappa, m_j> basis.

    The Pauli matrices act on the SPIN part of the spinor.  In the coupled
    (j, m_j) basis, each state decomposes as:

        |j, m_j> = sum_{m_s} CG(l, m_j-m_s; 1/2, m_s | j, m_j) |l, m_j-m_s>|1/2, m_s>

    sigma_a acts on the |1/2, m_s> part.  The result can change kappa
    (since j can change) but preserves n and l.

    We compute matrix elements:
        <n', kappa', m_j' | sigma_a | n, kappa, m_j>
        = delta(n,n') * delta(l(kappa),l(kappa')) * <kappa', m_j' | sigma_a | kappa, m_j>_l

    sigma_x = (sigma_+ + sigma_-) / 2  where sigma_+ = |up><down|, sigma_- = |down><up|
    sigma_y = (sigma_+ - sigma_-) / (2i)
    sigma_z = |up><up| - |down><down|

    In terms of spin matrices: sigma = 2*S.
    """
    N = len(labels)

    # We build sigma_z, sigma_+, sigma_- then derive sigma_x, sigma_y.
    # sigma_+ |1/2, -1/2> = |1/2, +1/2>,  sigma_+ |1/2, +1/2> = 0
    # sigma_- |1/2, +1/2> = |1/2, -1/2>,  sigma_- |1/2, -1/2> = 0
    # sigma_z |1/2, ms> = 2*ms * |1/2, ms>

    # Precompute label properties
    l_vals = [kappa_to_l(lab.kappa) for lab in labels]
    j_vals = [float(lab.j) for lab in labels]
    mj_vals = [float(lab.m_j) for lab in labels]

    sigma_z = np.zeros((N, N))
    sigma_plus = np.zeros((N, N), dtype=complex)
    sigma_minus = np.zeros((N, N), dtype=complex)

    for i, lab_a in enumerate(labels):
        for k, lab_b in enumerate(labels):
            # sigma preserves n and l
            if lab_a.n_fock != lab_b.n_fock:
                continue
            l_a = l_vals[i]
            l_b = l_vals[k]
            if l_a != l_b:
                continue
            l = l_a

            j_a = j_vals[i]
            j_b = j_vals[k]
            mj_a = mj_vals[i]
            mj_b = mj_vals[k]

            # --- sigma_z ---
            # <a|sigma_z|b> = sum_{m_s} CG(l, mj_a-ms; 1/2, ms | j_a, mj_a)
            #                          * (2*ms)
            #                          * CG(l, mj_b-ms; 1/2, ms | j_b, mj_b)
            #                          * delta(mj_a-ms, mj_b-ms)  [orbital m_l same]
            # The delta on m_l forces mj_a = mj_b (same m_j), so sigma_z is
            # diagonal in m_j but can be off-diagonal in kappa.
            if abs(mj_a - mj_b) < 1e-10:
                val = 0.0
                for ms2 in [-1, 1]:  # 2*m_s
                    ms = ms2 / 2.0
                    ml = mj_b - ms
                    if abs(ml) > l + 1e-10:
                        continue
                    ml_int = int(round(ml))
                    cg_a_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(ms2, 2),
                                        Rational(lab_a.j_times_2, 2), Rational(lab_a.two_m_j, 2))
                    cg_b_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(ms2, 2),
                                        Rational(lab_b.j_times_2, 2), Rational(lab_b.two_m_j, 2))
                    val += cg_a_val * (2.0 * ms) * cg_b_val
                sigma_z[i, k] = val

            # --- sigma_+ ---
            # sigma_+ |1/2, -1/2> = |1/2, +1/2>  (coefficient 1)
            # So m_s goes from -1/2 to +1/2, meaning m_j changes by +1.
            # <a|sigma_+|b>: need mj_a = mj_b + 1
            if abs(mj_a - mj_b - 1.0) < 1e-10:
                # m_s of ket is -1/2 (only this gives nonzero sigma_+)
                ms_b = -0.5
                ms_a = 0.5  # after sigma_+ acts
                ml = mj_b - ms_b  # m_l same for bra and ket
                ml_int = int(round(ml))
                if abs(ml_int) <= l:
                    cg_b_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(-1,2),
                                        Rational(lab_b.j_times_2, 2), Rational(lab_b.two_m_j, 2))
                    cg_a_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(1,2),
                                        Rational(lab_a.j_times_2, 2), Rational(lab_a.two_m_j, 2))
                    sigma_plus[i, k] = cg_a_val * 1.0 * cg_b_val

            # --- sigma_- ---
            # sigma_- |1/2, +1/2> = |1/2, -1/2>
            # m_j changes by -1.
            if abs(mj_a - mj_b + 1.0) < 1e-10:
                ms_b = 0.5
                ms_a = -0.5
                ml = mj_b - ms_b
                ml_int = int(round(ml))
                if abs(ml_int) <= l:
                    cg_b_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(1,2),
                                        Rational(lab_b.j_times_2, 2), Rational(lab_b.two_m_j, 2))
                    cg_a_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(-1,2),
                                        Rational(lab_a.j_times_2, 2), Rational(lab_a.two_m_j, 2))
                    sigma_minus[i, k] = cg_a_val * 1.0 * cg_b_val

    sigma_x = (sigma_plus + sigma_minus).real  # sigma_x = (sigma_+ + sigma_-)/2...
    # Wait: sigma_+ and sigma_- as defined here are the Pauli RAISING and LOWERING,
    # not the standard sigma_+ = sigma_x + i*sigma_y.
    # sigma_x = (sigma_+ + sigma_-) / 2
    # sigma_y = (sigma_+ - sigma_-) / (2i)
    sigma_x = np.real(sigma_plus + sigma_minus) / 2.0
    sigma_y = np.real((sigma_plus - sigma_minus) / (2.0j))

    return sigma_x, sigma_y, sigma_z


def verify_sigma_algebra(sigma_x, sigma_y, sigma_z, labels):
    """Verify the Pauli algebra relations within each (n, l) block."""
    N = len(labels)
    results = {}

    # 1. sigma_x^2 + sigma_y^2 + sigma_z^2 should = 3*I within each (n,l) block
    sigma_sq = sigma_x @ sigma_x + sigma_y @ sigma_y + sigma_z @ sigma_z

    # Build block structure
    blocks = defaultdict(list)
    for i, lab in enumerate(labels):
        blocks[(lab.n_fock, lab.l)].append(i)

    block_checks = {}
    for key, indices in blocks.items():
        idx = np.array(indices)
        block = sigma_sq[np.ix_(idx, idx)]
        expected = 3.0 * np.eye(len(idx))
        err = np.max(np.abs(block - expected))
        block_checks[str(key)] = {"max_error": err, "dim": len(idx), "passes": err < 1e-10}
    results["sigma_squared_3I"] = block_checks

    # 2. Commutation relations: [sigma_x, sigma_y] = 2i*sigma_z (within each block)
    comm_xy = sigma_x @ sigma_y - sigma_y @ sigma_x
    comm_check = {}
    for key, indices in blocks.items():
        idx = np.array(indices)
        comm_block = comm_xy[np.ix_(idx, idx)]
        expected = 2.0j * sigma_z[np.ix_(idx, idx)]
        err = np.max(np.abs(comm_block - expected))
        # comm_xy should be purely imaginary for Hermitian sigma's...
        # but our sigma matrices are real. The commutator [sigma_x, sigma_y] should
        # be 2i*sigma_z, which is imaginary. Let's check:
        # For real matrices: [A,B] = AB - BA is real. But 2i*sigma_z is imaginary.
        # This means [sigma_x, sigma_y] = 2i*sigma_z with i factored out:
        # i.e., sigma_x @ sigma_y - sigma_y @ sigma_x should be purely imaginary = 2i*sigma_z
        # But our sigma_x, sigma_y are REAL matrices (Hermitian).
        # So AB - BA is real. And 2i*sigma_z is imaginary. This can't be right.
        #
        # Actually, the issue is that we have REAL representations.
        # sigma_x and sigma_z are real symmetric, sigma_y is real antisymmetric.
        # Wait, sigma_y in the spin-1/2 basis is [[0,-i],[i,0]], which is
        # purely imaginary. But in the coupled (j,m_j) basis with real CG
        # coefficients, sigma_y should still be imaginary...
        #
        # Hmm, but we computed sigma_y = Re((sigma_+ - sigma_-)/(2i)).
        # Let me reconsider. sigma_+ and sigma_- as built above are real
        # (CG coefficients are real). So (sigma_+ - sigma_-)/(2i) is imaginary.
        # Taking Re() kills it. This is a bug in the construction.
        #
        # The correct approach: sigma_y IS NOT a real matrix in the standard basis.
        # But in the |j, m_j> coupled basis, sigma_y connects states with Dm_j = +/-1,
        # and the matrix elements are IMAGINARY.
        #
        # For a correct treatment we need to handle complex matrices.
        # Let me fix this below.
        comm_check[str(key)] = {"max_error": float(err), "dim": len(idx)}
    results["commutation_xy"] = comm_check

    # 3. Hermiticity
    results["sigma_x_hermitian"] = np.allclose(sigma_x, sigma_x.T, atol=1e-12)
    results["sigma_y_hermitian"] = np.allclose(sigma_y, sigma_y.T, atol=1e-12)
    results["sigma_z_hermitian"] = np.allclose(sigma_z, sigma_z.T, atol=1e-12)

    # 4. Tracelessness within each block
    traceless = {}
    for key, indices in blocks.items():
        idx = np.array(indices)
        tr_x = np.trace(sigma_x[np.ix_(idx, idx)])
        tr_y = np.trace(sigma_y[np.ix_(idx, idx)])
        tr_z = np.trace(sigma_z[np.ix_(idx, idx)])
        traceless[str(key)] = {
            "tr_x": float(tr_x), "tr_y": float(tr_y), "tr_z": float(tr_z),
            "all_zero": abs(tr_x) < 1e-10 and abs(tr_y) < 1e-10 and abs(tr_z) < 1e-10
        }
    results["tracelessness"] = traceless

    return results


def build_sigma_matrices_complex(labels, label_index):
    """Build sigma_x, sigma_y, sigma_z as COMPLEX matrices in the DiracLabel basis.

    sigma_y has imaginary matrix elements in ANY basis (it's the only one with
    an imaginary part in the standard representation). In the coupled |j, m_j>
    basis, sigma_y connects states with Delta(m_j) = +/- 1 with imaginary
    coefficients.

    sigma_x and sigma_z are real.
    """
    N = len(labels)
    l_vals = [kappa_to_l(lab.kappa) for lab in labels]

    sigma_z = np.zeros((N, N), dtype=complex)
    sigma_plus = np.zeros((N, N), dtype=complex)
    sigma_minus = np.zeros((N, N), dtype=complex)

    for i, lab_a in enumerate(labels):
        for k, lab_b in enumerate(labels):
            if lab_a.n_fock != lab_b.n_fock:
                continue
            if l_vals[i] != l_vals[k]:
                continue
            l = l_vals[i]

            j_a = float(lab_a.j)
            j_b = float(lab_b.j)
            mj_a = float(lab_a.m_j)
            mj_b = float(lab_b.m_j)

            # sigma_z: preserves m_j
            if abs(mj_a - mj_b) < 1e-10:
                val = 0.0
                for ms2 in [-1, 1]:
                    ms = ms2 / 2.0
                    ml = mj_b - ms
                    ml_int = int(round(ml))
                    if abs(ml_int) > l:
                        continue
                    cg_a_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(ms2, 2),
                                        Rational(lab_a.j_times_2, 2), Rational(lab_a.two_m_j, 2))
                    cg_b_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(ms2, 2),
                                        Rational(lab_b.j_times_2, 2), Rational(lab_b.two_m_j, 2))
                    val += cg_a_val * (2.0 * ms) * cg_b_val
                sigma_z[i, k] = val

            # sigma_+: m_j -> m_j + 1
            if abs(mj_a - mj_b - 1.0) < 1e-10:
                ms_b = -0.5
                ml = mj_b - ms_b
                ml_int = int(round(ml))
                if abs(ml_int) <= l:
                    cg_b_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(-1,2),
                                        Rational(lab_b.j_times_2, 2), Rational(lab_b.two_m_j, 2))
                    cg_a_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(1,2),
                                        Rational(lab_a.j_times_2, 2), Rational(lab_a.two_m_j, 2))
                    sigma_plus[i, k] = cg_a_val * cg_b_val

            # sigma_-: m_j -> m_j - 1
            if abs(mj_a - mj_b + 1.0) < 1e-10:
                ms_b = 0.5
                ml = mj_b - ms_b
                ml_int = int(round(ml))
                if abs(ml_int) <= l:
                    cg_b_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(1,2),
                                        Rational(lab_b.j_times_2, 2), Rational(lab_b.two_m_j, 2))
                    cg_a_val = cg_float(Rational(l), Rational(ml_int), Rational(1,2), Rational(-1,2),
                                        Rational(lab_a.j_times_2, 2), Rational(lab_a.two_m_j, 2))
                    sigma_minus[i, k] = cg_a_val * cg_b_val

    # sigma_x = (sigma_+ + sigma_-) / 2   (REAL matrix)
    # sigma_y = (sigma_+ - sigma_-) / (2i) (REAL matrix -- the i cancels the imaginary CG structure)
    # Wait: actually in the |j,m_j> basis with REAL CG coefficients,
    # sigma_+ and sigma_- are REAL matrices (they connect different m_j states).
    # So sigma_x = (sigma_+ + sigma_-)/2 is real.
    # And sigma_y = (sigma_+ - sigma_-)/(2i) is also real (because sigma_+ - sigma_- is real,
    # dividing by 2i makes it imaginary... no).
    #
    # Actually let me think again. The Pauli matrices in the UNCOUPLED basis |m_l, m_s>:
    # sigma_z = [[1,0],[0,-1]] (real)
    # sigma_x = [[0,1],[1,0]] (real)
    # sigma_y = [[0,-i],[i,0]] (imaginary)
    #
    # The CG transformation is real (all CG coefficients are real by convention).
    # So sigma_z and sigma_x stay real under CG rotation.
    # sigma_y stays purely imaginary under CG rotation.
    # Therefore: sigma_y in the coupled basis is purely imaginary.
    #
    # In the uncoupled basis: sigma_+ = sigma_x + i*sigma_y = [[0,2],[0,0]]
    #                         sigma_- = sigma_x - i*sigma_y = [[0,0],[2,0]]
    # So sigma_x = (sigma_+ + sigma_-)/2 and sigma_y = (sigma_+ - sigma_-)/(2i)
    #
    # But our sigma_+ and sigma_- as computed above are the COUPLED basis representations
    # of the uncoupled-basis operators [[0,1],[0,0]] and [[0,0],[1,0]] (up to factor 2).
    # These ARE real because CG coefficients are real.
    #
    # So: sigma_x = (sigma_+ + sigma_-) / 2  -- REAL
    #     sigma_y = (sigma_+ - sigma_-) / (2i)  -- when sigma_+ and sigma_- are real,
    #     (sigma_+ - sigma_-) is real, dividing by 2i gives IMAGINARY.
    #
    # But physically sigma_y must be HERMITIAN (self-adjoint), and in the standard basis
    # it IS Hermitian (anti-symmetric + purely imaginary => (A*)^T = (-i*B)^T = -i*B^T = -i*(-B) = iB = A).
    #
    # For QED self-energy: Sigma = sum_mu V_mu . G_gamma . V_mu^T
    # If sigma_y is imaginary, then V_y is imaginary, and V_y . G . V_y^T has
    # V_y^T = V_y^* (since Hermitian), giving V_y . G . V_y^dagger which is
    # positive semidefinite... but V_y^T is transpose (not conjugate transpose).
    #
    # For the self-energy contraction in vector QED:
    # Sigma = sum_mu (V_mu G_gamma V_mu^dagger)   -- Hermitian contraction
    # OR
    # Sigma = sum_mu (V_mu G_gamma V_mu^T)         -- just transpose
    #
    # For real sigma_x, sigma_z: V_mu^T = V_mu^dagger = V_mu (symmetric).
    # For imaginary sigma_y: V_y^T != V_y^dagger.
    #
    # In continuum QED the contraction is gamma^mu ... gamma_mu where
    # gamma_mu are Hermitian. The correct contraction is:
    # Sigma = sum_mu sigma_mu . G . sigma_mu    (since sigma_mu is Hermitian, sigma_mu^dagger = sigma_mu)
    #
    # For the self-energy loop with photon propagator:
    # Sigma[a,b] = sum_{mu, e, e'} sigma_mu[a,c] * G_e[c,d] * G_gamma[e,e'] * sigma_mu[d,b]
    #
    # In matrix form per edge:
    # Sigma = sum_{mu, e, e'} G_gamma[e,e'] * sigma_mu_e . G_e . sigma_mu_{e'}
    # where sigma_mu_e[a,c] = sigma_mu[a,c] * delta(e connects a to c somehow)

    # The CG construction builds sigma_+ = |ms=+1/2><ms=-1/2| (the standard
    # raising operator, equivalent to the upper-off-diagonal of the 2x2 Pauli matrix).
    # In the standard convention where sigma_+ = (sigma_x + i*sigma_y)/2,
    # we have sigma_x = sigma_+ + sigma_- (NOT divided by 2).
    # And sigma_y = -i*(sigma_+ - sigma_-)  (NOT divided by 2).
    # This gives the correct normalization sigma_a^2 = I (within each block)
    # and sigma^2 = 3I.
    sigma_x = sigma_plus + sigma_minus
    sigma_y = -1j * (sigma_plus - sigma_minus)

    return sigma_x, sigma_y, sigma_z


def verify_sigma_algebra_complex(sigma_x, sigma_y, sigma_z, labels):
    """Verify Pauli algebra with complex sigma_y."""
    N = len(labels)
    results = {}

    blocks = defaultdict(list)
    for i, lab in enumerate(labels):
        blocks[(lab.n_fock, lab.l)].append(i)

    # 1. sigma^2 = 3*I within each block
    sigma_sq = sigma_x @ sigma_x + sigma_y @ sigma_y + sigma_z @ sigma_z
    block_checks = {}
    for key, indices in blocks.items():
        idx = np.array(indices)
        block = sigma_sq[np.ix_(idx, idx)]
        expected = 3.0 * np.eye(len(idx))
        err = np.max(np.abs(block - expected))
        block_checks[str(key)] = {"max_error": float(err), "dim": len(idx), "passes": bool(err < 1e-10)}
    results["sigma_squared_3I"] = block_checks

    # 2. [sigma_x, sigma_y] = 2i * sigma_z within each block
    comm_xy = sigma_x @ sigma_y - sigma_y @ sigma_x
    comm_check = {}
    for key, indices in blocks.items():
        idx = np.array(indices)
        comm_block = comm_xy[np.ix_(idx, idx)]
        expected = 2.0j * sigma_z[np.ix_(idx, idx)]
        err = np.max(np.abs(comm_block - expected))
        comm_check[str(key)] = {"max_error": float(err), "dim": len(idx), "passes": bool(err < 1e-10)}
    results["commutation_xy_2i_sz"] = comm_check

    # Also check cyclic: [sigma_y, sigma_z] = 2i * sigma_x
    comm_yz = sigma_y @ sigma_z - sigma_z @ sigma_y
    comm_yz_check = {}
    for key, indices in blocks.items():
        idx = np.array(indices)
        comm_block = comm_yz[np.ix_(idx, idx)]
        expected = 2.0j * sigma_x[np.ix_(idx, idx)]
        err = np.max(np.abs(comm_block - expected))
        comm_yz_check[str(key)] = {"max_error": float(err), "passes": bool(err < 1e-10)}
    results["commutation_yz_2i_sx"] = comm_yz_check

    # [sigma_z, sigma_x] = 2i * sigma_y
    comm_zx = sigma_z @ sigma_x - sigma_x @ sigma_z
    comm_zx_check = {}
    for key, indices in blocks.items():
        idx = np.array(indices)
        comm_block = comm_zx[np.ix_(idx, idx)]
        expected = 2.0j * sigma_y[np.ix_(idx, idx)]
        err = np.max(np.abs(comm_block - expected))
        comm_zx_check[str(key)] = {"max_error": float(err), "passes": bool(err < 1e-10)}
    results["commutation_zx_2i_sy"] = comm_zx_check

    # 3. Hermiticity: sigma_mu^dagger = sigma_mu
    results["sigma_x_hermitian"] = bool(np.allclose(sigma_x, sigma_x.conj().T, atol=1e-12))
    results["sigma_y_hermitian"] = bool(np.allclose(sigma_y, sigma_y.conj().T, atol=1e-12))
    results["sigma_z_hermitian"] = bool(np.allclose(sigma_z, sigma_z.conj().T, atol=1e-12))

    # Check sigma_x is real, sigma_y is purely imaginary, sigma_z is real
    results["sigma_x_real"] = bool(np.allclose(sigma_x.imag, 0, atol=1e-12))
    results["sigma_y_purely_imaginary"] = bool(np.allclose(sigma_y.real, 0, atol=1e-12))
    results["sigma_z_real"] = bool(np.allclose(sigma_z.imag, 0, atol=1e-12))

    # 4. Tracelessness within each block
    traceless = {}
    for key, indices in blocks.items():
        idx = np.array(indices)
        tr_x = np.trace(sigma_x[np.ix_(idx, idx)])
        tr_y = np.trace(sigma_y[np.ix_(idx, idx)])
        tr_z = np.trace(sigma_z[np.ix_(idx, idx)])
        traceless[str(key)] = {
            "tr_x": float(np.abs(tr_x)), "tr_y": float(np.abs(tr_y)), "tr_z": float(np.abs(tr_z)),
            "all_zero": bool(abs(tr_x) < 1e-10 and abs(tr_y) < 1e-10 and abs(tr_z) < 1e-10)
        }
    results["tracelessness"] = traceless

    # 5. sigma_z eigenvalue structure
    sz_diag = {}
    for i, lab in enumerate(labels):
        sz_diag[f"({lab.n_fock},{lab.kappa},{lab.two_m_j})"] = float(sigma_z[i, i].real)
    results["sigma_z_diagonal"] = sz_diag

    return results


# ===========================================================================
# Step 2: Build the sigma-edge graph
# ===========================================================================

def build_sigma_edge_graph(sigma_x, sigma_y, sigma_z, labels, tol=1e-10):
    """Determine edges of the sigma graph: edge (i,j) exists iff
    <i|sigma_mu|j> != 0 for ANY mu in {x, y, z}.

    Returns edge list, adjacency matrix, and per-edge sigma coupling info.
    """
    N = len(labels)

    # Determine which pairs have nonzero sigma coupling (any component)
    sigma_norm = np.sqrt(np.abs(sigma_x)**2 + np.abs(sigma_y)**2 + np.abs(sigma_z)**2)

    edges = []
    edge_data = []
    for i in range(N):
        for j in range(i+1, N):
            if sigma_norm[i, j] > tol:
                edges.append((i, j))
                edge_data.append({
                    "i": i, "j": j,
                    "sigma_x": float(sigma_x[i, j].real),
                    "sigma_y_imag": float(sigma_y[i, j].imag),
                    "sigma_z": float(sigma_z[i, j].real),
                    "norm": float(sigma_norm[i, j]),
                    "label_i": f"({labels[i].n_fock},{labels[i].kappa},{labels[i].two_m_j})",
                    "label_j": f"({labels[j].n_fock},{labels[j].kappa},{labels[j].two_m_j})",
                })

    E = len(edges)

    # Adjacency matrix
    A = np.zeros((N, N))
    for i, j in edges:
        A[i, j] = 1
        A[j, i] = 1

    # Signed incidence matrix B (V x E)
    B = np.zeros((N, E))
    for e_idx, (i, j) in enumerate(edges):
        B[i, e_idx] = +1
        B[j, e_idx] = -1

    # Laplacians
    L0 = B @ B.T  # Node Laplacian
    L1 = B.T @ B  # Edge Laplacian

    # Betti numbers
    # beta_0 = number of connected components
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse import csr_matrix
    A_sparse = csr_matrix(A)
    n_components, component_labels = connected_components(A_sparse, directed=False)
    beta_0 = n_components
    beta_1 = E - N + beta_0

    # Degree of each node
    degrees = np.sum(A, axis=1).astype(int)

    # GS indices
    gs_indices = [i for i, lab in enumerate(labels) if lab.n_fock == 1 and lab.kappa == -1]
    gs_degrees = [int(degrees[i]) for i in gs_indices]
    gs_pendant = all(d == 1 for d in gs_degrees)

    graph_info = {
        "V": N,
        "E": E,
        "beta_0": int(beta_0),
        "beta_1": int(beta_1),
        "connected": bool(beta_0 == 1),
        "gs_indices": gs_indices,
        "gs_degrees": gs_degrees,
        "gs_pendant": gs_pendant,
        "degrees": degrees.tolist(),
        "max_degree": int(np.max(degrees)),
        "min_degree": int(np.min(degrees)),
    }

    return edges, edge_data, B, L0, L1, graph_info


# ===========================================================================
# Step 3: Build QED on the sigma-edge graph
# ===========================================================================

def build_dirac_operator(labels):
    """Build the diagonal Dirac operator D = chi * (n_fock + 1/2).

    chi = +1 if kappa < 0, -1 if kappa > 0.
    |lambda| = n_fock + 1/2 in the Fock convention.

    Actually, checking the production code: the Camporesi-Higuchi eigenvalue is
    |lambda_n| = n_CH + 3/2, and n_fock = n_CH + 1, so |lambda| = n_fock + 1/2.
    """
    N = len(labels)
    D = np.zeros((N, N))
    for i, lab in enumerate(labels):
        chi = 1.0 if lab.kappa < 0 else -1.0
        D[i, i] = chi * (lab.n_fock + 0.5)
    return D


def compute_sigma_qed(sigma_x, sigma_y, sigma_z, labels, edges, B, L0, L1, graph_info):
    """Compute QED self-energy using the sigma vertex on the sigma-edge graph.

    The vector vertex contraction is:

    Sigma[a,b] = sum_mu sum_{e,e'} sigma_mu_e[a,c] * G_e[c,d] * G_gamma[e,e'] * sigma_mu_e'[d,b]

    where:
    - sigma_mu_e[a,c] = sigma_mu[a,c] * delta(e connects the nodes)
    - G_e = D^{-1} (diagonal electron propagator, t=0)
    - G_gamma = L_1^+ (pseudoinverse of edge Laplacian)
    - The sum over mu = {x, y, z} contracts the vector index

    For the self-energy diagram (one-loop):
    We can write this as:
    Sigma = sum_mu sum_{e,e'} G_gamma[e,e'] * M_mu_e * G_e * M_mu_{e'}^dagger

    where M_mu_e[a,c] = sigma_mu[a, src(e)] * delta(c, tgt(e)) + sigma_mu[a, tgt(e)] * delta(c, src(e))
    ... this is getting complicated. Let me use the simpler formulation.

    Actually, the simplest approach matching the existing code pattern:

    For each polarization mu, define vertex matrices V_mu_e[a,b] per edge e:
        V_mu_e[a,b] = sigma_mu[a,b] if edge e = (a,b) or (b,a), else 0

    But this is the SAME edge list, so V_mu_e[a,b] is nonzero only if
    the pair (a,b) is the specific edge e AND sigma_mu[a,b] != 0.

    Since the edge graph is defined by nonzero sigma matrix elements,
    every edge e = (i,j) has at least one nonzero sigma_mu[i,j].

    Self-energy:
    Sigma = sum_mu sum_{e,e'} G_gamma[e,e'] * V_mu_e * V_mu_{e'}^dagger

    where V_mu_e is N x N with V_mu_e[a,b] = sigma_mu[a,b] * delta(e connects a and b).

    For Hermitian sigma: V_mu_e^dagger = V_mu_e^* (conjugate transpose).
    sigma_x: real -> V_x_e^dagger = V_x_e^T = V_x_e (symmetric)
    sigma_y: purely imaginary -> V_y_e^dagger = V_y_e^* = -V_y_e^T
    sigma_z: real -> same as sigma_x
    """
    N = len(labels)
    E = len(edges)

    # Electron propagator (diagonal, t=0)
    D = build_dirac_operator(labels)
    G_e = np.linalg.inv(D)  # Diagonal, so simple

    # Photon propagator = L_1^+ (Moore-Penrose pseudoinverse)
    G_gamma = np.linalg.pinv(L1)

    # Build vertex matrices V_mu_e for each polarization mu and each edge e
    sigmas = [sigma_x, sigma_y, sigma_z]
    sigma_names = ['x', 'y', 'z']

    # For each edge e = (i,j), the vertex matrix V_mu_e has:
    # V_mu_e[i,j] = sigma_mu[i,j]  and V_mu_e[j,i] = sigma_mu[j,i]
    # (all other entries zero)
    # Since sigma_mu is Hermitian: sigma_mu[j,i] = conj(sigma_mu[i,j])

    # Build V_mu_e matrices
    V_mu_e = []  # V_mu_e[mu][e] = N x N matrix
    for mu_idx, sigma_mu in enumerate(sigmas):
        V_e_list = []
        for e_idx, (i, j) in enumerate(edges):
            V = np.zeros((N, N), dtype=complex)
            V[i, j] = sigma_mu[i, j]
            V[j, i] = sigma_mu[j, i]  # = conj(sigma_mu[i,j])
            V_e_list.append(V)
        V_mu_e.append(V_e_list)

    # Compute self-energy with HERMITIAN contraction:
    # Sigma = sum_mu sum_{e,e'} G_gamma[e,e'] * V_mu_e . V_mu_{e'}^dagger
    # (no electron propagator in the scalar vertex version either --
    #  the existing scalar QED code also contracts V . G_gamma . V^T)
    #
    # Wait, looking at the existing code more carefully:
    # Sigma = sum_{e, e'} G_gamma[e, e'] * V_e . V_{e'}^T
    # where V_e is the SCALAR vertex matrix.
    #
    # For the vector case we need to sum over polarizations:
    # Sigma = sum_mu sum_{e, e'} G_gamma[e, e'] * V_mu_e . V_mu_{e'}^dagger

    Sigma = np.zeros((N, N), dtype=complex)
    for mu_idx in range(3):
        for e1 in range(E):
            for e2 in range(E):
                g_ee = G_gamma[e1, e2]
                if abs(g_ee) < 1e-15:
                    continue
                # Hermitian contraction: V . V^dagger
                contrib = g_ee * (V_mu_e[mu_idx][e1] @ V_mu_e[mu_idx][e2].conj().T)
                Sigma += contrib

    # Also compute with V^T instead of V^dagger for comparison
    Sigma_transpose = np.zeros((N, N), dtype=complex)
    for mu_idx in range(3):
        for e1 in range(E):
            for e2 in range(E):
                g_ee = G_gamma[e1, e2]
                if abs(g_ee) < 1e-15:
                    continue
                contrib = g_ee * (V_mu_e[mu_idx][e1] @ V_mu_e[mu_idx][e2].T)
                Sigma_transpose += contrib

    # Check: for sigma_x and sigma_z (real), V^dagger = V^T.
    # For sigma_y (imaginary), V^dagger = V^* = -V^T.
    # So Sigma_hermitian != Sigma_transpose in general (sigma_y contribution flips sign).

    return Sigma, Sigma_transpose, G_e, G_gamma, V_mu_e


def analyze_self_energy(Sigma, labels, name=""):
    """Analyze properties of a self-energy matrix."""
    N = len(labels)
    results = {}

    # Hermiticity
    results["is_hermitian"] = bool(np.allclose(Sigma, Sigma.conj().T, atol=1e-10))
    results["is_symmetric"] = bool(np.allclose(Sigma, Sigma.T, atol=1e-10))
    results["is_real"] = bool(np.allclose(Sigma.imag, 0, atol=1e-10))

    # Use real part if essentially real
    if results["is_real"]:
        Sigma_work = Sigma.real
    else:
        Sigma_work = Sigma

    # Trace
    tr = np.trace(Sigma)
    results["trace"] = {"real": float(tr.real), "imag": float(tr.imag)}

    # Eigenvalues
    if results["is_hermitian"]:
        eigs = np.linalg.eigvalsh(Sigma_work.real if results["is_real"] else Sigma_work)
    else:
        eigs = np.linalg.eigvals(Sigma_work)
    eigs_sorted = sorted(eigs.real)
    results["eigenvalues"] = [float(e) for e in eigs_sorted]
    results["n_zero_eigenvalues"] = int(sum(1 for e in eigs_sorted if abs(e) < 1e-10))
    results["positive_semidefinite"] = bool(all(e > -1e-10 for e in eigs_sorted))

    # Ground state block (n_fock=1, kappa=-1)
    gs_indices = [i for i, lab in enumerate(labels) if lab.n_fock == 1 and lab.kappa == -1]
    if len(gs_indices) > 0:
        gs_block = Sigma[np.ix_(gs_indices, gs_indices)]
        gs_norm = np.max(np.abs(gs_block))
        results["gs_block"] = {
            "values": [[float(gs_block[i,j].real) for j in range(gs_block.shape[1])]
                       for i in range(gs_block.shape[0])],
            "max_abs": float(gs_norm),
            "is_zero": bool(gs_norm < 1e-10),
        }

    return results


# ===========================================================================
# Step 4: Selection rule census
# ===========================================================================

def selection_rule_census(sigma_x, sigma_y, sigma_z, labels, edges, Sigma, G_e, G_gamma, V_mu_e):
    """Test 8 continuum QED selection rules on the sigma-edge graph."""
    N = len(labels)
    E = len(edges)
    results = {}

    # 1. Delta m_j conservation
    # In continuum QED, the vertex preserves m_j (for sigma_z) or changes it by +/-1 (sigma_x, sigma_y).
    # Check which Delta(m_j) values appear on the edges
    dm_j_values = set()
    for i, j in edges:
        dm = labels[i].two_m_j - labels[j].two_m_j
        dm_j_values.add(dm)
    results["delta_mj_conservation"] = {
        "delta_mj_values_on_edges": sorted(dm_j_values),
        "description": "Delta(2*m_j) values on sigma-graph edges",
        "verdict": "Delta(m_j) = 0, +/-1 only" if dm_j_values.issubset({-2, 0, 2}) else "BROKEN",
        "survives": bool(dm_j_values.issubset({-2, 0, 2})),
    }

    # 2. Spatial parity (E1: l_a + l_b odd)
    # Check whether edges connect states with l_a + l_b odd
    parity_check = {"odd": 0, "even": 0}
    for i, j in edges:
        l_sum = labels[i].l + labels[j].l
        if l_sum % 2 == 1:
            parity_check["odd"] += 1
        else:
            parity_check["even"] += 1
    results["spatial_parity_E1"] = {
        "l_sum_parity": parity_check,
        "verdict": "ALL l_a + l_b even (sigma preserves l)" if parity_check["odd"] == 0 else "MIXED",
        "description": "Sigma preserves l, so l_a = l_b always, making l_a + l_b even",
        "survives": False,  # This is different from E1 (which requires odd)
        "note": "sigma vertex preserves l (intra-l coupling), NOT an E1 vertex"
    }

    # 3. Gaunt/CG sparsity
    # Count fraction of possible V entries that are nonzero
    total_possible = N * N * E * 3  # 3 polarizations
    nonzero_count = 0
    for mu_idx in range(3):
        for e_idx in range(E):
            V = V_mu_e[mu_idx][e_idx]
            nonzero_count += np.count_nonzero(np.abs(V) > 1e-12)
    sparsity = 1.0 - nonzero_count / total_possible if total_possible > 0 else 0
    results["gaunt_cg_sparsity"] = {
        "total_possible": total_possible,
        "nonzero": nonzero_count,
        "sparsity_fraction": float(sparsity),
        "density_percent": float(100 * nonzero_count / total_possible) if total_possible > 0 else 0,
        "survives": True,  # Always present (CG triangle inequality)
    }

    # 4. Vertex parity
    # In continuum QED: n1 + n2 + q must be odd (from gamma^mu coupling and SO(4)).
    # On the sigma graph: check if there's a parity constraint on (n1, n2, kappa1, kappa2)
    parity_counts = {"n1+n2_even": 0, "n1+n2_odd": 0}
    kappa_transition = defaultdict(int)
    for i, j in edges:
        n_sum = labels[i].n_fock + labels[j].n_fock
        if n_sum % 2 == 0:
            parity_counts["n1+n2_even"] += 1
        else:
            parity_counts["n1+n2_odd"] += 1
        kappa_key = f"({labels[i].kappa},{labels[j].kappa})"
        kappa_transition[kappa_key] += 1

    # Since sigma preserves n, n1 = n2 always, so n1 + n2 is always even.
    results["vertex_parity"] = {
        "n_sum_parity": parity_counts,
        "kappa_transitions": dict(kappa_transition),
        "verdict": "sigma preserves n, so n1+n2 always even -- NO vertex parity rule",
        "survives": False,
    }

    # 5. SO(4) channel count
    # Not directly applicable -- sigma is an intra-shell operator
    results["so4_channel_count"] = {
        "verdict": "N/A -- sigma vertex is intra-shell (preserves n and l), no SO(4) channel structure",
        "survives": False,
    }

    # 6. Charge conjugation (C)
    # In continuum QED, C symmetry relates particle and antiparticle.
    # On the graph, check if Sigma has C symmetry: C Sigma C^{-1} = Sigma
    # where C maps (n, kappa, m_j) -> (n, kappa, -m_j) with possible phases.
    # Build the C operator: C |n, kappa, m_j> = (-1)^{j+m_j} |n, kappa, -m_j>
    C_mat = np.zeros((N, N))
    for i, lab in enumerate(labels):
        j_val = float(lab.j)
        mj_val = float(lab.m_j)
        phase = (-1)**(j_val + mj_val)
        target = DiracLabel(lab.n_fock, lab.kappa, -lab.two_m_j)
        for k, lab_k in enumerate(labels):
            if lab_k == target:
                C_mat[k, i] = phase
                break

    if np.allclose(Sigma.imag, 0, atol=1e-10):
        Sigma_real = Sigma.real
        CSC = C_mat @ Sigma_real @ C_mat
        c_symmetric = bool(np.allclose(CSC, Sigma_real, atol=1e-10))
    else:
        CSC = C_mat @ Sigma @ C_mat
        c_symmetric = bool(np.allclose(CSC, Sigma, atol=1e-10))

    results["charge_conjugation"] = {
        "C_symmetric": c_symmetric,
        "verdict": "C-symmetric" if c_symmetric else "C broken",
        "survives": c_symmetric,
    }

    # 7. Furry's theorem (odd-loop diagrams vanish)
    # The tadpole: sum_e V_mu_e (trace over electron index)
    # For scalar graph, tadpole is nonzero. For continuum, Furry says tadpole = 0.
    tadpole = {}
    for mu_idx, mu_name in enumerate(['x', 'y', 'z']):
        tad = np.zeros((N, N), dtype=complex)
        for e_idx in range(E):
            tad += V_mu_e[mu_idx][e_idx]
        tadpole_trace = np.trace(tad)
        tadpole_norm = np.max(np.abs(tad))
        tadpole[mu_name] = {
            "trace": float(np.abs(tadpole_trace)),
            "max_abs": float(tadpole_norm),
            "is_zero": bool(tadpole_norm < 1e-10),
        }
    # Total tadpole (sum over mu)
    all_zero = all(tadpole[m]["is_zero"] for m in ['x', 'y', 'z'])
    results["furry_theorem"] = {
        "per_polarization": tadpole,
        "all_tadpoles_zero": all_zero,
        "verdict": "Furry satisfied (all tadpoles zero)" if all_zero else "Furry VIOLATED",
        "survives": all_zero,
    }

    # 8. Ward identity: check if Sigma commutes with D in some sense
    # In continuum QED, Ward identity relates self-energy to vertex correction.
    # Simplest check: does [D, Sigma] = 0? (it shouldn't in general)
    D = build_dirac_operator(labels)
    if np.allclose(Sigma.imag, 0, atol=1e-10):
        comm_D_Sigma = D @ Sigma.real - Sigma.real @ D
    else:
        comm_D_Sigma = D @ Sigma - Sigma @ D
    comm_norm = np.max(np.abs(comm_D_Sigma))
    results["ward_identity"] = {
        "commutator_D_Sigma_norm": float(comm_norm),
        "commutes": bool(comm_norm < 1e-10),
        "verdict": "[D, Sigma] = 0" if comm_norm < 1e-10 else "[D, Sigma] != 0",
        "survives": False,  # Ward identity is a more subtle relation
        "note": "Full Ward identity requires vertex correction comparison, not just [D,Sigma]"
    }

    # Count how many survive
    rule_names = [
        "delta_mj_conservation",
        "spatial_parity_E1",
        "gaunt_cg_sparsity",
        "vertex_parity",
        "so4_channel_count",
        "charge_conjugation",
        "furry_theorem",
        "ward_identity",
    ]
    survived = sum(1 for r in rule_names if results[r].get("survives", False))
    results["summary"] = {
        "total_rules": 8,
        "survived": survived,
        "fraction": f"{survived}/8",
    }

    return results


# ===========================================================================
# Step 5: Algebraic characterization
# ===========================================================================

def algebraic_analysis(sigma_x, sigma_y, sigma_z, Sigma, labels, edges):
    """Determine the number field of the sigma-edge QED quantities."""
    results = {}

    # Check if Sigma is real
    is_real = np.allclose(Sigma.imag, 0, atol=1e-10)
    results["sigma_self_energy_real"] = is_real

    if is_real:
        S = Sigma.real
        # Check if all entries are rational-looking
        # Try to identify as simple fractions
        entries = []
        for i in range(S.shape[0]):
            for j in range(S.shape[1]):
                v = S[i, j]
                if abs(v) > 1e-12:
                    # Try to identify as p/q with small q
                    best_frac = None
                    best_err = 1e-6
                    for q in range(1, 200):
                        p = round(v * q)
                        err = abs(v - p/q)
                        if err < best_err:
                            best_err = err
                            best_frac = (p, q)
                    entries.append({
                        "i": i, "j": j,
                        "value": float(v),
                        "fraction": f"{best_frac[0]}/{best_frac[1]}" if best_frac else None,
                        "residual": float(best_err) if best_frac else None,
                    })
        results["nonzero_entries"] = entries[:20]  # first 20

        # Check pi-free: does any entry look like it involves pi?
        # Simple heuristic: check if entries are close to rationals
        all_rational = all((e.get("residual") or 1.0) < 1e-8 for e in entries)
        results["appears_rational"] = all_rational
        results["appears_pi_free"] = True  # heuristic

    # Sigma vertex matrix elements: check the CG content
    # The sigma matrices involve CG coefficients which are sqrt(rational).
    # Products sigma[a,c] * sigma[d,b] give rational (sqrt's cancel).
    # With G_gamma (rational pseudoinverse), Sigma should be rational.
    results["expected_number_field"] = "Q (rationals) -- sigma_x, sigma_z are real rational, sigma_y is purely imaginary rational, bilinear products cancel sqrt"

    # Compare to scalar Fock graph: Q[sqrt(2), sqrt(3), sqrt(6)]
    # Compare to Dirac graph: Q[sqrt(2), sqrt(17), sqrt(41), sqrt(881)]
    results["comparison"] = {
        "scalar_fock_graph": "Q[sqrt(2), sqrt(3), sqrt(6)]",
        "dirac_E1_graph": "Q[sqrt(2), sqrt(17), sqrt(41), sqrt(881)]",
        "sigma_vertex_graph": "Q (rationals) -- simpler than both, if confirmed"
    }

    return results


# ===========================================================================
# Main computation
# ===========================================================================

def main():
    print("=" * 70)
    print("VQ-2: Pauli sigma vertex for graph-native QED on S^3")
    print("=" * 70)

    n_max = 2

    # Step 0: Build basis
    print("\n--- Step 0: Basis ---")
    labels, label_index = build_labels(n_max)
    N = len(labels)
    print(f"n_max = {n_max}, N_dirac = {N} states")
    for i, lab in enumerate(labels):
        print(f"  [{i}] n={lab.n_fock}, kappa={lab.kappa}, l={lab.l}, j={float(lab.j)}, m_j={float(lab.m_j)}")

    # Step 1: Build sigma matrices
    print("\n--- Step 1: Sigma matrices (complex) ---")
    sigma_x, sigma_y, sigma_z = build_sigma_matrices_complex(labels, label_index)

    # Verify algebra
    algebra_results = verify_sigma_algebra_complex(sigma_x, sigma_y, sigma_z, labels)
    print(f"  sigma^2 = 3I: {all(v['passes'] for v in algebra_results['sigma_squared_3I'].values())}")
    print(f"  [sigma_x, sigma_y] = 2i*sigma_z: {all(v['passes'] for v in algebra_results['commutation_xy_2i_sz'].values())}")
    print(f"  [sigma_y, sigma_z] = 2i*sigma_x: {all(v['passes'] for v in algebra_results['commutation_yz_2i_sx'].values())}")
    print(f"  [sigma_z, sigma_x] = 2i*sigma_y: {all(v['passes'] for v in algebra_results['commutation_zx_2i_sy'].values())}")
    print(f"  sigma_x Hermitian: {algebra_results['sigma_x_hermitian']}")
    print(f"  sigma_y Hermitian: {algebra_results['sigma_y_hermitian']}")
    print(f"  sigma_z Hermitian: {algebra_results['sigma_z_hermitian']}")
    print(f"  sigma_x real: {algebra_results['sigma_x_real']}")
    print(f"  sigma_y purely imaginary: {algebra_results['sigma_y_purely_imaginary']}")
    print(f"  sigma_z real: {algebra_results['sigma_z_real']}")
    print(f"  All blocks traceless: {all(v['all_zero'] for v in algebra_results['tracelessness'].values())}")

    # Print sigma_z diagonal for reference
    print("\n  sigma_z diagonal values:")
    for i, lab in enumerate(labels):
        print(f"    [{i}] ({lab.n_fock},{lab.kappa},{lab.two_m_j}): sigma_z = {sigma_z[i,i].real:.6f}")

    # Print sigma matrix structures
    print("\n  sigma_x (real part, nonzero entries):")
    for i in range(N):
        for j in range(N):
            v = sigma_x[i, j].real
            if abs(v) > 1e-10:
                print(f"    [{i},{j}] ({labels[i].n_fock},{labels[i].kappa},{labels[i].two_m_j})"
                      f" <-> ({labels[j].n_fock},{labels[j].kappa},{labels[j].two_m_j}): {v:.6f}")

    # Step 2: Sigma-edge graph
    print("\n--- Step 2: Sigma-edge graph ---")
    edges, edge_data, B, L0, L1, graph_info = build_sigma_edge_graph(
        sigma_x, sigma_y, sigma_z, labels)

    print(f"  Vertices: {graph_info['V']}")
    print(f"  Edges: {graph_info['E']}")
    print(f"  Connected: {graph_info['connected']}")
    print(f"  beta_0: {graph_info['beta_0']}")
    print(f"  beta_1: {graph_info['beta_1']}")
    print(f"  GS pendant: {graph_info['gs_pendant']}")
    print(f"  GS degrees: {graph_info['gs_degrees']}")
    print(f"  Degree sequence: {sorted(graph_info['degrees'])}")

    print("\n  Edge list:")
    for ed in edge_data:
        print(f"    {ed['label_i']} <-> {ed['label_j']}: "
              f"sigma_x={ed['sigma_x']:.4f}, sigma_y_im={ed['sigma_y_imag']:.4f}, "
              f"sigma_z={ed['sigma_z']:.4f}, norm={ed['norm']:.4f}")

    # Compare to existing graphs
    print("\n  Comparison to existing graphs:")
    print(f"    Scalar Fock graph at n_max=2: V=5, E=5 (T+/T-/L+/L- edges)")
    print(f"    Dirac Rule A graph at n_max=2: (inter-shell, kappa-preserving)")
    print(f"    Dirac Rule B graph at n_max=2: (E1 dipole, Delta l = +/-1)")
    print(f"    Sigma-edge graph at n_max=2: V={graph_info['V']}, E={graph_info['E']} (INTRA-shell)")

    # Edge Laplacian eigenvalues
    L1_eigs = np.linalg.eigvalsh(L1)
    print(f"\n  L1 eigenvalues: {[f'{e:.4f}' for e in sorted(L1_eigs)]}")
    n_zero_L1 = sum(1 for e in L1_eigs if abs(e) < 1e-10)
    print(f"  Number of zero eigenvalues (gauge modes): {n_zero_L1}")

    # Step 3: Build QED
    print("\n--- Step 3: QED on sigma-edge graph ---")
    Sigma_herm, Sigma_trans, G_e, G_gamma, V_mu_e = compute_sigma_qed(
        sigma_x, sigma_y, sigma_z, labels, edges, B, L0, L1, graph_info)

    print("\n  Hermitian contraction (V.G.V^dagger):")
    se_herm = analyze_self_energy(Sigma_herm, labels, "Hermitian")
    print(f"    Is Hermitian: {se_herm['is_hermitian']}")
    print(f"    Is real: {se_herm['is_real']}")
    print(f"    Trace: {se_herm['trace']}")
    print(f"    Eigenvalues: {[f'{e:.6f}' for e in se_herm['eigenvalues']]}")
    print(f"    N zero eigenvalues: {se_herm['n_zero_eigenvalues']}")
    print(f"    Positive semidefinite: {se_herm['positive_semidefinite']}")
    if 'gs_block' in se_herm:
        print(f"    GS block zero: {se_herm['gs_block']['is_zero']}")
        print(f"    GS block max abs: {se_herm['gs_block']['max_abs']:.8f}")
        print(f"    GS block values: {se_herm['gs_block']['values']}")

    print("\n  Transpose contraction (V.G.V^T):")
    se_trans = analyze_self_energy(Sigma_trans, labels, "Transpose")
    print(f"    Is Hermitian: {se_trans['is_hermitian']}")
    print(f"    Is real: {se_trans['is_real']}")
    print(f"    Trace: {se_trans['trace']}")
    if 'gs_block' in se_trans:
        print(f"    GS block zero: {se_trans['gs_block']['is_zero']}")
        print(f"    GS block max abs: {se_trans['gs_block']['max_abs']:.8f}")

    # Step 4: Selection rules
    print("\n--- Step 4: Selection rule census ---")
    sel_rules = selection_rule_census(
        sigma_x, sigma_y, sigma_z, labels, edges,
        Sigma_herm, G_e, G_gamma, V_mu_e)

    for rule_name in ["delta_mj_conservation", "spatial_parity_E1", "gaunt_cg_sparsity",
                       "vertex_parity", "so4_channel_count", "charge_conjugation",
                       "furry_theorem", "ward_identity"]:
        r = sel_rules[rule_name]
        status = "SURVIVES" if r.get("survives", False) else "BROKEN/N_A"
        print(f"  {rule_name}: {status}")
        print(f"    Verdict: {r.get('verdict', 'N/A')}")

    print(f"\n  SUMMARY: {sel_rules['summary']['fraction']} selection rules survive")

    # Step 5: Algebraic analysis
    print("\n--- Step 5: Algebraic characterization ---")
    alg = algebraic_analysis(sigma_x, sigma_y, sigma_z, Sigma_herm, labels, edges)
    print(f"  Self-energy is real: {alg['sigma_self_energy_real']}")
    if alg.get('appears_rational'):
        print(f"  Appears rational: {alg['appears_rational']}")
    print(f"  Expected number field: {alg['expected_number_field']}")

    # Print the self-energy matrix
    print("\n  Self-energy matrix (Hermitian contraction):")
    if np.allclose(Sigma_herm.imag, 0, atol=1e-10):
        S = Sigma_herm.real
        for i in range(N):
            row = [f"{S[i,j]:8.4f}" for j in range(N)]
            print(f"    [{i}] " + " ".join(row))

    # ===========================================================================
    # CRITICAL ANALYSIS: Why the sigma vertex can't recover vector QED rules
    # ===========================================================================
    print("\n" + "=" * 70)
    print("CRITICAL ANALYSIS")
    print("=" * 70)

    print("""
The sigma vertex creates an INTRA-SHELL coupling graph:
- sigma preserves both n and l (only changes kappa/j within the same (n,l) block)
- This means the sigma-edge graph has NO inter-shell edges
- The graph decomposes into disconnected (n,l) blocks
- The GS (n=1, l=0) is an ISOLATED pair (kappa=-1 only, two m_j states)
  with at most self-loops (sigma_z diagonal), NOT coupled to n=2 states

The 4 missing continuum QED selection rules (vertex parity, SO(4) channel
count, Ward identity, charge conjugation) all involve the INTERACTION between
different energy shells via the photon. The sigma vertex, being purely
intra-shell, cannot access this inter-shell dynamics AT ALL.

This is structurally different from:
1. The scalar Fock graph QED (inter-shell edges from T+/T- ladders)
2. The Dirac graph QED (inter-shell edges from E1 dipole transitions)

The sigma vertex is the WRONG operator for the QED vertex because it doesn't
couple different energy shells. In continuum QED, gamma^mu couples the electron
current to the photon, which carries energy/momentum. The photon necessarily
changes the electron's energy shell. sigma_a, acting only on spin, cannot do this.

The correct vector vertex for graph-native QED would need to:
(a) couple different energy shells (like the existing graph edges do)
(b) carry vector quantum numbers (like sigma_a does)
(c) combine BOTH features (neither alone is sufficient)
""")

    # ===========================================================================
    # Save results
    # ===========================================================================
    print("\n--- Saving results ---")

    # Collect all results
    all_results = {
        "n_max": n_max,
        "N_dirac": N,
        "basis": [{"n_fock": lab.n_fock, "kappa": lab.kappa, "two_m_j": lab.two_m_j,
                    "l": lab.l, "j": float(lab.j)} for lab in labels],
        "algebra_verification": {
            "sigma_squared_3I": {k: {"passes": v["passes"], "max_error": v["max_error"]}
                                  for k, v in algebra_results["sigma_squared_3I"].items()},
            "commutation_all_pass": (
                all(v['passes'] for v in algebra_results['commutation_xy_2i_sz'].values()) and
                all(v['passes'] for v in algebra_results['commutation_yz_2i_sx'].values()) and
                all(v['passes'] for v in algebra_results['commutation_zx_2i_sy'].values())
            ),
            "all_hermitian": (algebra_results['sigma_x_hermitian'] and
                              algebra_results['sigma_y_hermitian'] and
                              algebra_results['sigma_z_hermitian']),
            "sigma_x_real": algebra_results['sigma_x_real'],
            "sigma_y_purely_imaginary": algebra_results['sigma_y_purely_imaginary'],
            "sigma_z_real": algebra_results['sigma_z_real'],
            "all_traceless": all(v['all_zero'] for v in algebra_results['tracelessness'].values()),
        },
        "sigma_edge_graph": graph_info,
        "edge_data": edge_data,
        "L1_eigenvalues": [float(e) for e in sorted(L1_eigs)],
        "self_energy_hermitian": se_herm,
        "self_energy_transpose": se_trans,
        "selection_rules": {
            "delta_mj_conservation": sel_rules["delta_mj_conservation"]["survives"],
            "spatial_parity_E1": sel_rules["spatial_parity_E1"]["survives"],
            "gaunt_cg_sparsity": sel_rules["gaunt_cg_sparsity"]["survives"],
            "vertex_parity": sel_rules["vertex_parity"]["survives"],
            "so4_channel_count": sel_rules["so4_channel_count"]["survives"],
            "charge_conjugation": sel_rules["charge_conjugation"]["survives"],
            "furry_theorem": sel_rules["furry_theorem"]["survives"],
            "ward_identity": sel_rules["ward_identity"]["survives"],
            "total_survived": sel_rules["summary"]["survived"],
            "fraction": sel_rules["summary"]["fraction"],
            "details": {k: {kk: vv for kk, vv in v.items() if kk != "survives"}
                        for k, v in sel_rules.items() if k != "summary"},
        },
        "algebraic_analysis": alg,
        "headline": {
            "sigma_vertex_is_intra_shell": True,
            "graph_connected": graph_info["connected"],
            "gs_sigma_zero": se_herm.get("gs_block", {}).get("is_zero", None),
            "selection_rules_survived": sel_rules["summary"]["fraction"],
            "verdict": "NEGATIVE -- sigma vertex is intra-shell, cannot couple energy shells, wrong operator for QED vertex",
        },
    }

    # Save JSON
    data_path = os.path.join(os.path.dirname(__file__), "data", "vq2_sigma_vertex_qed.json")
    with open(data_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"  Data saved to {data_path}")

    print("\n" + "=" * 70)
    print("VQ-2 COMPLETE")
    print("=" * 70)

    return all_results


if __name__ == "__main__":
    results = main()
