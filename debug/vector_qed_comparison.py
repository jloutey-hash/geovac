"""
Vector QED vs Scalar Graph-Native QED Comparison
=================================================

Builds explicit vector QED on the finite Fock graph electron basis
and compares to the existing scalar graph-native QED (geovac/graph_qed_*.py).

The vector QED has:
  - Photon modes (q, L=1, M_L) with q=1..q_max
  - Vertex parity rule: n_a + n_b + q = odd
  - Angular momentum conservation via Gaunt integral / 3j symbols
  - Parity selection: l_a + l_b + 1 must be odd => |l_a - l_b| = odd
  - Triangle rule: |l_a - l_b| <= 1 <= l_a + l_b

The scalar graph QED has:
  - Photon on edges of the Fock scalar graph (nearest-neighbor only)
  - Vertex via CG projection (no parity, no channel count)
  - No vertex parity rule

Key question: Does the scalar graph correctly capture low-q (IR) photon
contributions, with discrepancy growing at high-q (UV)?
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy import sqrt as sp_sqrt, Rational, N as sp_N

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.graph_qed_self_energy import compute_self_energy
from geovac.graph_qed_photon import build_fock_graph
from geovac.lattice import GeometricLattice


# ===========================================================================
# Step 1: Electron basis (scalar Fock graph nodes)
# ===========================================================================

def build_electron_basis(n_max: int) -> List[Tuple[int, int, int]]:
    """Build the scalar Fock electron basis: (n, l, m) states.

    States: n=1..n_max, l=0..n-1, m=-l..l
    Total: n_max^2 states.
    """
    states = []
    for n in range(1, n_max + 1):
        for l in range(0, n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
    return states


# ===========================================================================
# Step 2: Vector photon modes
# ===========================================================================

def build_photon_modes(q_max: int) -> List[Tuple[int, int]]:
    """Build vector photon modes (q, M_L) with L=1 fixed.

    q = 1, 2, ..., q_max
    M_L = -1, 0, +1

    Photon propagator: G_photon(q) = 1 / omega_q where omega_q = q(q+2)
    (transverse vector Laplacian eigenvalue on unit S^3, co-exact 1-forms).
    """
    modes = []
    for q in range(1, q_max + 1):
        for M_L in [-1, 0, 1]:
            modes.append((q, M_L))
    return modes


def photon_propagator(q: int) -> float:
    """Photon propagator 1/omega_q where omega_q = q(q+2)."""
    omega = q * (q + 2)
    return 1.0 / omega


# ===========================================================================
# Step 3: Vector vertex tensor
# ===========================================================================

def vertex_parity_allowed(n_a: int, n_b: int, q: int) -> bool:
    """Check vertex parity: n_a + n_b + q must be odd.

    This uses the Fock convention where n starts at 1.
    In the CH convention (n starts at 0), this would be (n_a-1) + (n_b-1) + q odd.

    Actually, the continuum vertex parity rule is stated in CH convention:
      n_ext + n_int + q is odd, where n_ext, n_int are CH indices (starting at 0).

    With Fock convention (n starts at 1), we need:
      (n_a - 1) + (n_b - 1) + q = odd
    => n_a + n_b + q = odd + 2 = odd  (parity preserved)

    So in Fock convention: n_a + n_b + q is odd.
    """
    return (n_a + n_b + q) % 2 == 1


def triangle_allowed(n_a: int, n_b: int, q: int) -> bool:
    """Check SO(4) triangle inequality for photon mode q.

    In CH convention: |n_ext - n_int| <= q <= n_ext + n_int.
    In Fock convention (shift by 1): |(n_a-1) - (n_b-1)| <= q <= (n_a-1) + (n_b-1)
    => |n_a - n_b| <= q <= n_a + n_b - 2
    """
    n_a_ch = n_a - 1
    n_b_ch = n_b - 1
    return abs(n_a_ch - n_b_ch) <= q <= n_a_ch + n_b_ch


def angular_coupling(l_a: int, m_a: int, l_b: int, m_b: int, M_L: int) -> float:
    """Gaunt-type angular coupling for vector (E1) photon.

    A(l_a, m_a; l_b, m_b; M_L) =
        (-1)^{m_a} * sqrt((2l_a+1)*(2l_b+1)*3/(4*pi))
        * wigner_3j(l_a, 1, l_b, -m_a, M_L, m_b)
        * wigner_3j(l_a, 1, l_b, 0, 0, 0)

    Selection rules from the 3j symbols:
      - m_b = m_a + M_L (from first 3j: -m_a + M_L + m_b = 0)
        Wait, the 3j symbol (l_a, 1, l_b; -m_a, M_L, m_b) requires
        -m_a + M_L + m_b = 0 => m_b = m_a - M_L
      - l_a + 1 + l_b must be even for the second 3j to be nonzero
        (all m's are zero => sum of j's must be even)
        => l_a + l_b + 1 even => l_a + l_b is odd => |l_a - l_b| is odd
      - Triangle: |l_a - l_b| <= 1 <= l_a + l_b

    Note: The m-conservation from first 3j is:
      3j(l_a, 1, l_b; -m_a, M_L, m_b) requires (-m_a) + M_L + m_b = 0
      => m_b = m_a - M_L
    """
    # Check m-conservation: m_b = m_a - M_L
    if m_b != m_a - M_L:
        return 0.0

    # Check parity: l_a + l_b + 1 must be even (for second 3j with all-zero m)
    if (l_a + l_b + 1) % 2 != 0:
        return 0.0

    # Check triangle for L=1 photon
    if abs(l_a - l_b) > 1 or l_a + l_b < 1:
        return 0.0

    # Compute the 3j symbols
    threej_1 = float(wigner_3j(l_a, 1, l_b, -m_a, M_L, m_b))
    threej_2 = float(wigner_3j(l_a, 1, l_b, 0, 0, 0))

    if abs(threej_1) < 1e-15 or abs(threej_2) < 1e-15:
        return 0.0

    # Prefactor
    prefactor = ((-1)**m_a) * np.sqrt((2*l_a + 1) * (2*l_b + 1) * 3 / (4 * np.pi))

    return prefactor * threej_1 * threej_2


def compute_vector_vertex(
    states: List[Tuple[int, int, int]],
    q_max: int,
    use_radial_weight: bool = True,
) -> np.ndarray:
    """Compute the vector vertex tensor V_vector[a, b, (q, M_L)].

    V_vector[a, b, (q,M)] = R(n_a, l_a; n_b, l_b; q) * A(l_a, m_a; l_b, m_b; M_L)

    Parameters
    ----------
    states : list of (n, l, m)
    q_max : int
    use_radial_weight : bool
        If True, use R = 1/sqrt(omega_q) as radial weight.
        If False, use R = 1 (structure-only mode).

    Returns
    -------
    V : ndarray of shape (N_states, N_states, 3*q_max)
        V[a, b, photon_idx] where photon_idx = (q-1)*3 + (M_L+1)
    """
    N = len(states)
    N_photon = 3 * q_max
    V = np.zeros((N, N, N_photon))

    for a_idx, (n_a, l_a, m_a) in enumerate(states):
        for b_idx, (n_b, l_b, m_b) in enumerate(states):
            for q in range(1, q_max + 1):
                # Check vertex parity
                if not vertex_parity_allowed(n_a, n_b, q):
                    continue

                # Check SO(4) triangle
                if not triangle_allowed(n_a, n_b, q):
                    continue

                # Radial coupling weight
                if use_radial_weight:
                    R = 1.0 / np.sqrt(q * (q + 2))
                else:
                    R = 1.0

                for M_L in [-1, 0, 1]:
                    A = angular_coupling(l_a, m_a, l_b, m_b, M_L)
                    if abs(A) < 1e-15:
                        continue

                    photon_idx = (q - 1) * 3 + (M_L + 1)
                    V[a_idx, b_idx, photon_idx] = R * A

    return V


# ===========================================================================
# Step 4: Compute vector self-energy
# ===========================================================================

def compute_vector_self_energy(
    states: List[Tuple[int, int, int]],
    q_max: int,
    use_radial_weight: bool = True,
) -> Tuple[np.ndarray, Dict]:
    """Compute the vector QED one-loop self-energy.

    Sigma_vector[a, b] = sum_{c, q, M_L} V[a, c, (q,M)] * G_photon(q) * V[c, b, (q,M)]

    Returns Sigma matrix and per-q decomposition.
    """
    N = len(states)
    V = compute_vector_vertex(states, q_max, use_radial_weight)

    Sigma = np.zeros((N, N))
    Sigma_per_q = {}  # Sigma decomposed by photon mode q

    for q in range(1, q_max + 1):
        G_q = photon_propagator(q)
        Sigma_q = np.zeros((N, N))

        for M_L in [-1, 0, 1]:
            photon_idx = (q - 1) * 3 + (M_L + 1)
            # V[:, :, photon_idx] is the vertex matrix for this photon mode
            V_mode = V[:, :, photon_idx]
            # Sigma_q += G_q * V_mode @ V_mode.T
            # Note: Sum over intermediate state c:
            # Sigma[a,b] = sum_c V[a,c,mode] * G_q * V[c,b,mode]
            #            = G_q * (V_mode @ V_mode.T)[a,b]
            # Wait -- V[c,b,mode] means we need V transposed in second index:
            # Actually Sigma[a,b] = sum_c V[a,c,mode] * V[b,c,mode] * G_q
            #                     = G_q * sum_c V[a,c,mode] * V[b,c,mode]
            #                     = G_q * (V_mode @ V_mode.T)[a,b]
            # But the second vertex should be V[c,b] not V[b,c].
            # Let's be careful:
            # Sigma[a,b] = sum_{c,mode} V[a,c,mode] * G_photon * V[c,b,mode]
            # = G_q * sum_c V[a,c,mode] * V[c,b,mode]
            # = G_q * (V_mode @ V_mode)[a,b]   -- matrix multiplication!
            #
            # Wait, that's V_mode[a,c] * V_mode[c,b] = (V_mode @ V_mode)[a,b]
            # which is V_mode squared, not V @ V^T.
            #
            # Actually, the self-energy has the vertex structure:
            # electron a -> (emit photon mode) -> intermediate c -> (absorb photon mode) -> electron b
            # First vertex: V[a, c, mode]  (electron a emits, becomes c)
            # Second vertex: V[c, b, mode] (electron c absorbs, becomes b)
            #
            # So: Sigma[a,b] = G_q * sum_c V[a,c,mode] * V[c,b,mode]
            #                = G_q * (V @ V)[a,b] where V = V_mode here
            #
            # This IS matrix multiplication V_mode @ V_mode (V squared, not V @ V^T).
            # UNLESS the second vertex is the TRANSPOSE (absorption = time-reversed emission).
            #
            # In the scalar graph code (graph_qed_self_energy.py line 366):
            #   Sigma += G_gamma[e1, e2] * (V_mats[e1] * V_mats[e2].T)
            # There V_mats[e2].T is used because the vertex tensor is defined as
            # V[a, b, e] = coupling when electron goes a->b with photon e,
            # and the SECOND vertex has electron going c->b, so we need V[c, b, e]
            # which for the second vertex is V_e2[c, b] = V_e2.T[b, c] transposed.
            #
            # Actually let's look more carefully at the scalar code.
            # V_mats[e] is N x N where V_mats[e][a, b] = V[a, b, e].
            # Sigma = sum_{e1,e2} G[e1,e2] * V_e1 * V_e2^T
            # => Sigma[a, b] = sum_{e1,e2,c} G[e1,e2] * V[a,c,e1] * V[b,c,e2]
            #                = sum_{e1,e2,c} G[e1,e2] * V[a,c,e1] * V^T[c,b,e2]
            # Hmm, this uses V_e2^T which means contracting over the FIRST index of V_e2.
            #
            # For vector QED with diagonal photon propagator G[mode,mode] = G_q:
            # Sigma[a,b] = sum_{c, mode} V[a,c,mode] * G_q * V[b,c,mode]
            #            = G_q * sum_c V[a,c,mode] * V[b,c,mode]
            #            = G_q * (V_mode @ V_mode^T)[a,b]
            #
            # So it IS V @ V^T when the coupling is symmetric: the second vertex
            # couples (b -> c) in the time-reversed sense, meaning we read it as
            # V[b, c, mode] which gives V_mode[b, c].
            # Then sum_c V_mode[a,c] * V_mode[b,c] = (V_mode @ V_mode^T)[a,b].

            Sigma_q += G_q * (V_mode @ V_mode.T)

        Sigma_per_q[q] = Sigma_q
        Sigma += Sigma_q

    info = {
        'q_max': q_max,
        'n_states': N,
        'use_radial_weight': use_radial_weight,
    }

    return Sigma, Sigma_per_q


# ===========================================================================
# Step 5: Load/compute scalar self-energy
# ===========================================================================

def compute_scalar_self_energy_numpy(n_max: int) -> np.ndarray:
    """Compute the scalar graph-native self-energy at n_max (numpy path).

    This replicates the scalar self-energy computation from
    geovac/graph_qed_self_energy.py but returns in the SCALAR FOCK basis
    (n, l, m) rather than the Dirac (n, kappa, m_j) basis.

    The scalar graph self-energy is:
      Sigma[v1, v2] = sum_{e1, e2} G_gamma[e1, e2] * delta(v1 in e1) * delta(v2 in e2)

    where the photon propagator G_gamma = L_1^+ is the pseudoinverse of the
    edge Laplacian.

    Actually, the scalar self-energy in graph_qed_self_energy.py works in the
    DIRAC basis. For a fair comparison, we need to project back to the scalar
    Fock basis OR compute both in the same basis.

    Let's compute the scalar self-energy directly on the Fock graph nodes
    (bypassing the Dirac projection) for a clean structural comparison.
    """
    # Build the Fock graph
    fock_data = build_fock_graph(n_max)
    V = fock_data.V
    E = fock_data.E
    edges = fock_data.edges
    states = fock_data.states

    # Build incidence matrix B (V x E)
    B = np.zeros((V, E), dtype=float)
    for k, (i, j) in enumerate(edges):
        B[i, k] = 1.0
        B[j, k] = -1.0

    # Edge Laplacian L1 = B^T @ B (E x E)
    L1 = B.T @ B

    # Photon propagator = pseudoinverse of L1
    # Use eigendecomposition for the pseudoinverse
    eigvals, eigvecs = np.linalg.eigh(L1)
    G_gamma = np.zeros((E, E))
    for i in range(E):
        if eigvals[i] > 1e-10:
            G_gamma += (1.0 / eigvals[i]) * np.outer(eigvecs[:, i], eigvecs[:, i])

    # Scalar vertex on the Fock graph:
    # The vertex couples two NODES through an EDGE.
    # V_scalar[v1, v2, e] = 1 if edge e connects v1 and v2, else 0
    # (This is the simplest graph vertex: identity coupling on each edge)

    # Build vertex matrices V_e (V x V) for each edge e
    # V_e[v1, v2] = 1 if e = (v1, v2) or e = (v2, v1), else 0
    # i.e., V_e = |v1><v2| + |v2><v1| for edge e = (v1, v2)

    # Self-energy: Sigma[a, b] = sum_{e1, e2} G_gamma[e1, e2] * V_e1[a, :] @ V_e2[:, b]
    # But V_e[a, c] is nonzero only when a is an endpoint of e and c is the other endpoint.
    # So V_e1[a, c] = delta((a,c) is edge e1 or (c,a) is edge e1)

    # More precisely:
    # V_e[a, b] = 1 if edge e connects a to b (undirected)
    # So V_e is a symmetric matrix with exactly two nonzero entries: V_e[i,j] = V_e[j,i] = 1
    # for edge e = (i, j).

    # Sigma[a, b] = sum_{e1, e2, c} V_e1[a, c] * G_gamma[e1, e2] * V_e2[c, b]
    #             = sum_{e1, e2} G_gamma[e1, e2] * (V_e1 @ V_e2)[a, b]

    Sigma_scalar = np.zeros((V, V))

    # Pre-build the V_e matrices (sparse but just store endpoint pairs)
    for e1 in range(E):
        i1, j1 = edges[e1]
        for e2 in range(E):
            g_ee = G_gamma[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            i2, j2 = edges[e2]
            # (V_e1 @ V_e2)[a, b] = sum_c V_e1[a,c] * V_e2[c,b]
            # V_e1[a,c] nonzero only for (a,c) = (i1,j1) or (j1,i1)
            # V_e2[c,b] nonzero only for (c,b) = (i2,j2) or (j2,i2)
            # So: (V_e1 @ V_e2)[a,b] = V_e1[a, i2]*V_e2[i2, b] + V_e1[a, j2]*V_e2[j2, b]
            # V_e2[i2, b]: nonzero when b = j2 (value 1)
            # V_e2[j2, b]: nonzero when b = i2 (value 1)
            # V_e1[a, i2]: nonzero when (a, i2) is an endpoint pair of e1
            # V_e1[a, j2]: nonzero when (a, j2) is an endpoint pair of e1

            # Enumerate the 4 possible contributions:
            # 1. V_e1[i1, j1]=1, so if j1==i2: V_e1[i1, i2]=1, V_e2[i2, j2]=1 => Sigma[i1, j2] += g
            # 2. V_e1[i1, j1]=1, so if j1==j2: V_e1[i1, j2]=1, V_e2[j2, i2]=1 => Sigma[i1, i2] += g
            # 3. V_e1[j1, i1]=1, so if i1==i2: V_e1[j1, i2]=1, V_e2[i2, j2]=1 => Sigma[j1, j2] += g
            # 4. V_e1[j1, i1]=1, so if i1==j2: V_e1[j1, j2]=1, V_e2[j2, i2]=1 => Sigma[j1, i2] += g

            # Wait, let me be more careful.
            # V_e1[a, c] = 1 iff {a, c} = {i1, j1} (edge e1 connects a and c)
            # V_e2[c, b] = 1 iff {c, b} = {i2, j2} (edge e2 connects c and b)
            #
            # So (V_e1 @ V_e2)[a, b] = sum_c V_e1[a,c] * V_e2[c,b]
            # Nonzero contributions need c such that {a,c}={i1,j1} AND {c,b}={i2,j2}
            #
            # Case c = j1 (requires a = i1):
            #   Need {j1, b} = {i2, j2} => b = i2 or b = j2 (and j1 is the other)
            #   If j1 == i2: b = j2, contribute to Sigma[i1, j2]
            #   If j1 == j2: b = i2, contribute to Sigma[i1, i2]
            #
            # Case c = i1 (requires a = j1):
            #   Need {i1, b} = {i2, j2} =>
            #   If i1 == i2: b = j2, contribute to Sigma[j1, j2]
            #   If i1 == j2: b = i2, contribute to Sigma[j1, i2]

            if j1 == i2:
                Sigma_scalar[i1, j2] += g_ee
            if j1 == j2:
                Sigma_scalar[i1, i2] += g_ee
            if i1 == i2:
                Sigma_scalar[j1, j2] += g_ee
            if i1 == j2:
                Sigma_scalar[j1, i2] += g_ee

    return Sigma_scalar


# ===========================================================================
# Step 6: Comparison analysis
# ===========================================================================

def selection_rule_census(
    states: List[Tuple[int, int, int]],
    Sigma_scalar: np.ndarray,
    Sigma_vector: np.ndarray,
    threshold: float = 1e-12,
) -> Dict:
    """Count selection rule violations.

    For each pair (a, b), check which vector QED selection rules
    are violated by the scalar self-energy.
    """
    N = len(states)

    # For each element, classify why vector says zero or nonzero
    violations = {
        'vertex_parity': 0,        # scalar nonzero, vector zero due to parity
        'm_conservation': 0,       # scalar nonzero, vector zero due to delta_m
        'l_parity': 0,             # scalar nonzero, vector zero due to l_a+l_b odd
        'triangle_l': 0,           # scalar nonzero, vector zero due to triangle on l
        'triangle_n': 0,           # scalar nonzero, vector zero due to SO(4) triangle
        'scalar_only_nonzero': 0,  # scalar nonzero, vector zero (any reason)
        'vector_only_nonzero': 0,  # vector nonzero, scalar zero
        'both_nonzero': 0,         # both nonzero
        'both_zero': 0,            # both zero
    }

    # Detailed violation tracking
    violation_details = []

    for a in range(N):
        for b in range(N):
            s_val = abs(Sigma_scalar[a, b])
            v_val = abs(Sigma_vector[a, b])

            s_nz = s_val > threshold
            v_nz = v_val > threshold

            if s_nz and v_nz:
                violations['both_nonzero'] += 1
            elif not s_nz and not v_nz:
                violations['both_zero'] += 1
            elif s_nz and not v_nz:
                violations['scalar_only_nonzero'] += 1
                # Diagnose WHY vector is zero
                n_a, l_a, m_a = states[a]
                n_b, l_b, m_b = states[b]

                # Check each selection rule
                # The self-energy Sigma[a,b] involves a sum over intermediate c
                # and photon modes. Vector Sigma[a,b] = 0 means NO intermediate
                # state c connects a to b through any allowed photon.
                #
                # For the diagonal (a=b), vertex parity is the main constraint.
                # For off-diagonal, it's more complex.

                # Simple diagnostic: check if ANY q allows the transition
                any_q_allowed = False
                for q in range(1, max(n_a, n_b) * 2 + 2):
                    if vertex_parity_allowed(n_a, n_b, q) and triangle_allowed(n_a, n_b, q):
                        any_q_allowed = True
                        break

                detail = {
                    'a': states[a], 'b': states[b],
                    'scalar_val': float(s_val),
                    'reason': []
                }

                if not any_q_allowed:
                    # No direct vertex connects a to b (but self-energy is second-order!)
                    pass

                # Check m-conservation constraint: for Sigma[a,b] to be nonzero,
                # need an intermediate c such that V[a,c,mode] and V[c,b,mode] nonzero.
                # Angular coupling needs: m_c = m_a - M_L and m_b = m_c - M_L' (summed)
                # For same-mode: m_b = m_a - 2*M_L (not the right constraint)
                # Actually the sum over M_L relaxes this.

                # For self-energy (same photon mode absorbed and emitted):
                # Sigma involves V[a,c,(q,M)] * V[c,b,(q,M)] summed over c and (q,M)
                # But V[c,b,(q,M)] means electron c -> b emitting photon (q,M)
                # Wait no: V[c,b,mode] = coupling for (c -> b + photon_mode)?
                # Or is it V[b,c,mode] = coupling for (b -> c + photon)?
                # In our convention V[a,b,mode] = coupling strength for
                # "electron a, electron b, photon mode" interaction.
                # The self-energy has: emit at first vertex, absorb at second.
                # First vertex: electron a -> electron c + photon
                # Second vertex: electron c + photon -> electron b
                # So: V[a,c,mode] * V[b,c,mode] (contracting over c)
                # which gives (V @ V^T)[a,b]

                violation_details.append(detail)

            elif v_nz and not s_nz:
                violations['vector_only_nonzero'] += 1

    return violations


def diagonal_comparison(
    states: List[Tuple[int, int, int]],
    Sigma_scalar: np.ndarray,
    Sigma_vector: np.ndarray,
) -> Dict:
    """Compare diagonal elements (self-energy shifts) between scalar and vector."""
    N = len(states)
    diag_scalar = np.diag(Sigma_scalar)
    diag_vector = np.diag(Sigma_vector)

    entries = []
    for i in range(N):
        entries.append({
            'state': states[i],
            'n': states[i][0],
            'l': states[i][1],
            'm': states[i][2],
            'Sigma_scalar': float(diag_scalar[i]),
            'Sigma_vector': float(diag_vector[i]),
            'ratio': float(diag_vector[i] / diag_scalar[i]) if abs(diag_scalar[i]) > 1e-15 else None,
        })

    # Check if proportional
    if np.all(np.abs(diag_scalar) > 1e-15):
        ratios = diag_vector / diag_scalar
        ratio_cv = np.std(ratios) / np.mean(ratios) if np.mean(ratios) != 0 else float('inf')
    else:
        ratios = None
        ratio_cv = None

    return {
        'entries': entries,
        'mean_ratio': float(np.mean(ratios)) if ratios is not None else None,
        'ratio_cv': float(ratio_cv) if ratio_cv is not None else None,
        'is_proportional': ratio_cv is not None and ratio_cv < 0.01,
        'frobenius_scalar': float(np.linalg.norm(diag_scalar)),
        'frobenius_vector': float(np.linalg.norm(diag_vector)),
    }


def per_q_decomposition(
    states: List[Tuple[int, int, int]],
    Sigma_scalar: np.ndarray,
    Sigma_per_q: Dict[int, np.ndarray],
) -> Dict:
    """Decompose vector self-energy by photon mode q and compare to scalar."""
    results = {}

    cumulative = np.zeros_like(Sigma_scalar)
    frob_scalar = np.linalg.norm(Sigma_scalar)

    for q in sorted(Sigma_per_q.keys()):
        Sigma_q = Sigma_per_q[q]
        cumulative += Sigma_q

        frob_q = np.linalg.norm(Sigma_q)
        frob_cumulative = np.linalg.norm(cumulative)

        # Compare cumulative to scalar
        diff = cumulative - Sigma_scalar
        frob_diff = np.linalg.norm(diff)

        # Correlation
        if frob_scalar > 1e-15 and frob_cumulative > 1e-15:
            corr = np.sum(cumulative * Sigma_scalar) / (frob_cumulative * frob_scalar)
        else:
            corr = 0.0

        results[q] = {
            'frobenius_this_q': float(frob_q),
            'frobenius_cumulative': float(frob_cumulative),
            'frobenius_diff_from_scalar': float(frob_diff),
            'fraction_of_scalar': float(frob_cumulative / frob_scalar) if frob_scalar > 0 else 0,
            'correlation_with_scalar': float(corr),
        }

    return results


# ===========================================================================
# Main comparison driver
# ===========================================================================

def run_comparison(n_max: int, q_max: int = None, use_radial_weight: bool = True) -> Dict:
    """Run the full vector vs scalar QED comparison.

    Parameters
    ----------
    n_max : int
        Electron basis truncation (n_max^2 states).
    q_max : int, optional
        Maximum photon mode. Default: 2*n_max.
    use_radial_weight : bool
        Whether to weight radial coupling by 1/sqrt(omega_q).
    """
    if q_max is None:
        q_max = 2 * n_max

    print(f"\n{'='*70}")
    print(f"Vector vs Scalar QED Comparison: n_max={n_max}, q_max={q_max}")
    print(f"{'='*70}")

    # Build electron basis
    states = build_electron_basis(n_max)
    N = len(states)
    print(f"\nElectron basis: {N} states (n_max={n_max})")
    for i, s in enumerate(states):
        print(f"  [{i}] n={s[0]}, l={s[1]}, m={s[2]}")

    # Compute vector self-energy
    print(f"\nComputing vector self-energy (q_max={q_max}, radial_weight={use_radial_weight})...")
    Sigma_vector, Sigma_per_q = compute_vector_self_energy(states, q_max, use_radial_weight)

    # Compute scalar self-energy on same basis
    print("Computing scalar graph self-energy on Fock basis...")
    Sigma_scalar = compute_scalar_self_energy_numpy(n_max)

    # Ensure dimensions match
    assert Sigma_scalar.shape == (N, N), f"Shape mismatch: scalar {Sigma_scalar.shape} vs expected ({N},{N})"
    assert Sigma_vector.shape == (N, N), f"Shape mismatch: vector {Sigma_vector.shape}"

    # Basic properties
    print(f"\n--- Matrix Properties ---")
    print(f"Sigma_scalar: Frobenius norm = {np.linalg.norm(Sigma_scalar):.6f}")
    print(f"Sigma_vector: Frobenius norm = {np.linalg.norm(Sigma_vector):.6f}")
    print(f"Sigma_scalar symmetric: {np.allclose(Sigma_scalar, Sigma_scalar.T, atol=1e-12)}")
    print(f"Sigma_vector symmetric: {np.allclose(Sigma_vector, Sigma_vector.T, atol=1e-12)}")

    # Sparsity pattern
    threshold = 1e-12
    nnz_scalar = np.sum(np.abs(Sigma_scalar) > threshold)
    nnz_vector = np.sum(np.abs(Sigma_vector) > threshold)
    print(f"\nNonzero entries: scalar={nnz_scalar}/{N*N}, vector={nnz_vector}/{N*N}")
    print(f"Sparsity: scalar={1-nnz_scalar/(N*N):.3f}, vector={1-nnz_vector/(N*N):.3f}")

    # Selection rule census
    print(f"\n--- Selection Rule Census ---")
    census = selection_rule_census(states, Sigma_scalar, Sigma_vector, threshold)
    for key, count in census.items():
        if count > 0:
            print(f"  {key}: {count}")

    # Fraction of scalar Frobenius norm that is "allowed" by vector rules
    # Elements where both are nonzero
    mask_both = (np.abs(Sigma_scalar) > threshold) & (np.abs(Sigma_vector) > threshold)
    frob_allowed = np.linalg.norm(Sigma_scalar[mask_both])
    frob_total = np.linalg.norm(Sigma_scalar)
    print(f"\nFraction of scalar ||.||_F in 'allowed' positions: "
          f"{frob_allowed/frob_total:.4f}" if frob_total > 0 else "N/A")

    # Diagonal comparison
    print(f"\n--- Diagonal (Self-Energy Shifts) ---")
    diag_comp = diagonal_comparison(states, Sigma_scalar, Sigma_vector)
    for entry in diag_comp['entries']:
        s = entry['state']
        r_str = f"{entry['ratio']:.4f}" if entry['ratio'] is not None else "N/A"
        print(f"  ({s[0]},{s[1]},{s[2]}): scalar={entry['Sigma_scalar']:.6f}, "
              f"vector={entry['Sigma_vector']:.6f}, ratio={r_str}")

    if diag_comp['mean_ratio'] is not None:
        print(f"\n  Mean ratio (vector/scalar): {diag_comp['mean_ratio']:.6f}")
        print(f"  Ratio CV: {diag_comp['ratio_cv']:.6f}")
        print(f"  Proportional: {diag_comp['is_proportional']}")

    # Per-q decomposition
    print(f"\n--- Per-q Decomposition ---")
    per_q = per_q_decomposition(states, Sigma_scalar, Sigma_per_q)
    print(f"  {'q':<4} {'||Sigma_q||':<14} {'||cumul||':<14} {'||diff||':<14} {'frac_scalar':<14} {'corr':<10}")
    for q in sorted(per_q.keys()):
        d = per_q[q]
        print(f"  {q:<4} {d['frobenius_this_q']:<14.6f} {d['frobenius_cumulative']:<14.6f} "
              f"{d['frobenius_diff_from_scalar']:<14.6f} {d['fraction_of_scalar']:<14.4f} "
              f"{d['correlation_with_scalar']:<10.6f}")

    # Print full matrices
    print(f"\n--- Sigma_scalar matrix ---")
    print(np.array2string(Sigma_scalar, precision=4, suppress_small=True))

    print(f"\n--- Sigma_vector matrix ---")
    print(np.array2string(Sigma_vector, precision=4, suppress_small=True))

    # Eigenvalues
    eig_scalar = np.sort(np.linalg.eigvalsh(Sigma_scalar))
    eig_vector = np.sort(np.linalg.eigvalsh(Sigma_vector))
    print(f"\n--- Eigenvalues ---")
    print(f"  Scalar: {eig_scalar}")
    print(f"  Vector: {eig_vector}")

    # Ground state analysis
    # The ground state is (1, 0, 0) = index 0
    print(f"\n--- Ground State Analysis ---")
    gs_idx = 0
    print(f"  Sigma_scalar[GS, GS] = {Sigma_scalar[gs_idx, gs_idx]:.8f}")
    print(f"  Sigma_vector[GS, GS] = {Sigma_vector[gs_idx, gs_idx]:.8f}")
    print(f"  Scalar GS row nonzeros: {np.sum(np.abs(Sigma_scalar[gs_idx, :]) > threshold)}")
    print(f"  Vector GS row nonzeros: {np.sum(np.abs(Sigma_vector[gs_idx, :]) > threshold)}")

    # Vertex parity check for ground state
    # For GS (n=1), self-energy needs V[GS, c, mode] * V[GS, c, mode]
    # Angular coupling from GS (l=0, m=0): need l_c=1 (dipole), m_c = 0 - M_L
    # Vertex parity for (n_GS=1, n_c, q): 1 + n_c + q must be odd
    # => n_c + q must be even
    # Triangle: |0 - (n_c-1)| <= q <= 0 + (n_c-1) => q = n_c - 1
    # So: n_c + (n_c - 1) = 2*n_c - 1 = always odd! Never even!
    # This means vertex parity FORBIDS any coupling from GS for the only
    # triangle-allowed q value!
    print(f"\n  VERTEX PARITY ANALYSIS for GS (n=1):")
    print(f"  For GS (n_Fock=1, n_CH=0): triangle requires q = n_c - 1")
    print(f"  Parity requires 1 + n_c + q = odd => n_c + q even")
    print(f"  With q = n_c - 1: n_c + (n_c-1) = 2*n_c - 1 = ALWAYS ODD")
    print(f"  => Vertex parity FORBIDS all couplings from GS!")
    print(f"  => Sigma_vector[GS, :] should be identically zero.")
    print(f"  Actual max |Sigma_vector[GS, :]| = {np.max(np.abs(Sigma_vector[gs_idx, :])):.2e}")

    # Assemble results
    results = {
        'n_max': n_max,
        'q_max': q_max,
        'use_radial_weight': use_radial_weight,
        'n_states': N,
        'states': [list(s) for s in states],
        'matrix_properties': {
            'frobenius_scalar': float(np.linalg.norm(Sigma_scalar)),
            'frobenius_vector': float(np.linalg.norm(Sigma_vector)),
            'nnz_scalar': int(nnz_scalar),
            'nnz_vector': int(nnz_vector),
            'symmetric_scalar': bool(np.allclose(Sigma_scalar, Sigma_scalar.T, atol=1e-12)),
            'symmetric_vector': bool(np.allclose(Sigma_vector, Sigma_vector.T, atol=1e-12)),
        },
        'selection_rule_census': census,
        'fraction_scalar_in_allowed': float(frob_allowed/frob_total) if frob_total > 0 else 0,
        'diagonal_comparison': diag_comp,
        'per_q_decomposition': per_q,
        'eigenvalues_scalar': [float(x) for x in eig_scalar],
        'eigenvalues_vector': [float(x) for x in eig_vector],
        'ground_state': {
            'scalar_sigma_gs': float(Sigma_scalar[gs_idx, gs_idx]),
            'vector_sigma_gs': float(Sigma_vector[gs_idx, gs_idx]),
            'vector_gs_row_max': float(np.max(np.abs(Sigma_vector[gs_idx, :]))),
            'vertex_parity_forbids_gs': True,
        },
        'Sigma_scalar_matrix': Sigma_scalar.tolist(),
        'Sigma_vector_matrix': Sigma_vector.tolist(),
    }

    return results


# ===========================================================================
# Main
# ===========================================================================

if __name__ == "__main__":
    output_dir = Path(__file__).parent / "data"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = {}

    # Run at n_max=2 (4 states)
    results_2 = run_comparison(n_max=2, q_max=4, use_radial_weight=True)
    all_results['n_max_2'] = results_2

    print("\n\n")

    # Run at n_max=3 (9 states)
    results_3 = run_comparison(n_max=3, q_max=6, use_radial_weight=True)
    all_results['n_max_3'] = results_3

    # Also run with R=1 (no radial weight) for structural comparison
    print("\n\n")
    print("="*70)
    print("STRUCTURAL COMPARISON (R=1, no radial weight)")
    print("="*70)
    results_struct = run_comparison(n_max=2, q_max=4, use_radial_weight=False)
    all_results['n_max_2_no_radial'] = results_struct

    # Save
    output_path = output_dir / "vector_qed_comparison.json"
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n\nResults saved to: {output_path}")

    # Final summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\n1. Ground state protection:")
    print(f"   Scalar: Sigma[GS,GS] = {results_2['ground_state']['scalar_sigma_gs']:.6f} (NONZERO)")
    print(f"   Vector: Sigma[GS,GS] = {results_2['ground_state']['vector_sigma_gs']:.6f} (should be ZERO)")
    print(f"   => Vertex parity PROTECTS the ground state in vector QED")
    print(f"   => The scalar graph BREAKS this protection (pendant-edge theorem)")

    print(f"\n2. Selection rule violations (n_max=2):")
    c = results_2['selection_rule_census']
    print(f"   Scalar-only nonzero: {c['scalar_only_nonzero']} elements")
    print(f"   Vector-only nonzero: {c['vector_only_nonzero']} elements")
    print(f"   Both nonzero: {c['both_nonzero']} elements")

    print(f"\n3. Frobenius norm fraction in allowed positions: "
          f"{results_2['fraction_scalar_in_allowed']:.4f}")

    print(f"\n4. Per-q convergence: Does scalar match low-q vector?")
    pq = results_2['per_q_decomposition']
    for q in sorted(pq.keys()):
        d = pq[q]
        print(f"   q={q}: cumul correlation with scalar = {d['correlation_with_scalar']:.4f}")
