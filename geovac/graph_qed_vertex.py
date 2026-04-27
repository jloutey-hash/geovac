"""
Electron-photon vertex coupling and one-loop vacuum polarization on the
GeoVac graph.
==========================================================================

The electron lives on the Dirac graph (nodes = |n, kappa, m_j>).
The photon lives on edges of the scalar Fock graph.
The vertex coupling bridges these two spaces.

Vertex construction
-------------------
Each Dirac state |n, kappa, m_j> has orbital content (n, l, m_l)
through the Clebsch-Gordan decomposition

    |j, m_j> = sum_{m_l, m_s} <l, m_l, 1/2, m_s | j, m_j> |l, m_l> |m_s>

This defines a projection matrix P[a, v] from Dirac node a to scalar
Fock node v = (n, l, m), where the projection coefficient is the CG
coefficient for m_l = m(v) and m_s = m_j(a) - m(v), evaluated at the
correct (l, j) of the Dirac state.

The vertex tensor V[a, b, e] (Dirac state a, Dirac state b, photon
edge e) is:

    V[a, b, e] = sum_{v1, v2} P[a, v1] * P[b, v2] * delta(e = (v1, v2))

This is a purely graph-native construction:
- P uses CG coefficients (algebraic, pi-free)
- The edge connectivity comes from the Fock graph structure
- No continuum formulas are used

Physics: the vertex couples the electron CURRENT (bilinear in Dirac
states, projected to orbital content) to the photon field (on Fock
edges). This is the graph analog of the gamma^mu coupling in
continuum QED, where gamma^mu projects the spinor bilinear onto
the vector (gauge field) representation.

One-loop vacuum polarization
-----------------------------
The one-loop vacuum polarization tensor Pi[e1, e2] is:

    Pi[e1, e2] = sum_{a,b,c,d} V[a,b,e1] * G_e[b,c] * V[c,d,e2] * G_e[d,a]

This is a finite matrix trace:

    Pi = V^T . G_e . V . G_e    (contracted over electron indices)

where V is reshaped as (N_e, N_dirac^2) and G_e is the electron
propagator. The result is an E x E matrix (photon edge space).

At finite n_max this is an EXACT FINITE TRACE -- no infinite sums,
no regularization, no renormalization. All entries are algebraic
(pi-free by the algebraic-integer theorem for rational matrix
operations).

Transcendental taxonomy (Paper 18)
-----------------------------------
- Projection matrix P: algebraic (CG coefficients are sqrt(rational))
- Vertex tensor V: algebraic (products of CG coefficients)
- Vacuum polarization Pi: rational (trace over products of rational
  propagator entries and algebraic vertex entries; the sqrt's cancel
  in the bilinear pairing)

References
----------
- GeoVac graph_qed_propagator.py (electron propagator)
- GeoVac graph_qed_photon.py (photon propagator)
- GeoVac fock_graph_hodge.py (Hodge infrastructure)
- GeoVac dirac_matrix_elements.py (DiracLabel, kappa_to_l)
- GeoVac qed_vertex.py (continuum vertex for comparison)
- GeoVac qed_vacuum_polarization.py (continuum VP for comparison)
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import sympy as sp
from sympy import Integer, Matrix, Rational, sqrt, zeros as sp_zeros
from sympy.physics.quantum.cg import CG

from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels, kappa_to_l
from geovac.graph_qed_propagator import (
    DiracGraphOperator,
    electron_propagator,
)
from geovac.graph_qed_photon import build_fock_graph, FockGraphData
from geovac.lattice import GeometricLattice

__all__ = [
    "build_projection_matrix",
    "build_vertex_tensor",
    "compute_vacuum_polarization",
    "vertex_pi_free_certificate",
    "vertex_summary",
    "analyze_vertex_and_vp",
]


# ---------------------------------------------------------------------------
# Step 1: Projection matrix P[a, v] from Dirac to Fock
# ---------------------------------------------------------------------------

def _cg_coefficient(l: int, m_l: int, m_s: Rational,
                    j: Rational, m_j: Rational) -> sp.Expr:
    """Clebsch-Gordan coefficient <l, m_l, 1/2, m_s | j, m_j>.

    Returns the exact sympy expression (typically sqrt(rational)).
    """
    s_half = Rational(1, 2)
    result = CG(l, m_l, s_half, m_s, j, m_j).doit()
    return sp.nsimplify(result, rational=False)


def build_projection_matrix(
    n_max: int,
) -> Tuple[Matrix, List[DiracLabel], List[Tuple[int, int, int]]]:
    """Build the projection matrix P from Dirac states to scalar Fock states.

    P[a, v] = CG coefficient projecting Dirac state a onto scalar Fock node v.

    Each Dirac state |n, kappa, m_j> decomposes into orbital content
    |n, l, m_l> x |m_s> via Clebsch-Gordan coupling:

        |j, m_j> = sum_{m_s} <l, m_l, 1/2, m_s | j, m_j> |l, m_l> |m_s>

    where m_l = m_j - m_s.

    The projection P[a, v] is nonzero when:
    - Dirac state a has n_fock = n(v) and l(kappa) = l(v)
    - m(v) = m_j(a) - m_s for some valid m_s in {-1/2, +1/2}

    In that case, P[a, v] = <l(v), m(v), 1/2, m_s | j(a), m_j(a)>.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    P : sympy.Matrix (N_dirac x V_fock)
        Projection matrix with algebraic (sqrt-rational) entries.
    dirac_labels : list of DiracLabel
        Ordered Dirac state labels (rows of P).
    fock_states : list of (n, l, m) tuples
        Ordered scalar Fock state labels (columns of P).
    """
    dirac_labels = list(iter_dirac_labels(n_max))
    N_dirac = len(dirac_labels)

    lat = GeometricLattice(max_n=n_max, topological_weights=False)
    fock_states = lat.states
    V_fock = len(fock_states)
    fock_index = {s: i for i, s in enumerate(fock_states)}

    P = sp_zeros(N_dirac, V_fock)

    for a_idx, dlab in enumerate(dirac_labels):
        n_a = dlab.n_fock
        l_a = dlab.l
        j_a = dlab.j  # Rational
        m_j_a = dlab.m_j  # Rational

        # Two possible spin projections
        for two_ms in [-1, 1]:
            m_s = Rational(two_ms, 2)
            m_l = m_j_a - m_s  # Rational; must be integer for scalar state

            # Check m_l is integer
            if m_l.denominator != 1:
                continue
            m_l_int = int(m_l)

            # Check m_l in valid range for l_a
            if abs(m_l_int) > l_a:
                continue

            # Scalar Fock node
            fock_node = (n_a, l_a, m_l_int)
            if fock_node not in fock_index:
                continue

            v_idx = fock_index[fock_node]

            # CG coefficient
            cg = _cg_coefficient(l_a, m_l_int, m_s, j_a, m_j_a)
            P[a_idx, v_idx] = P[a_idx, v_idx] + cg

    return P, dirac_labels, fock_states


# ---------------------------------------------------------------------------
# Step 2: Vertex tensor V[a, b, e]
# ---------------------------------------------------------------------------

def build_vertex_tensor(
    n_max: int,
    P: Optional[Matrix] = None,
    dirac_labels: Optional[List[DiracLabel]] = None,
    fock_data: Optional[FockGraphData] = None,
) -> Tuple[List[Tuple[int, int, int, sp.Expr]], int, int, int]:
    """Build the vertex coupling tensor V[a, b, e].

    V[a, b, e] = sum_{v1, v2} P[a, v1] * P[b, v2] * delta(e connects v1, v2)

    where P is the Dirac-to-Fock projection matrix.

    Rather than storing a dense 3D array, returns a sparse list of
    nonzero entries (a, b, e, value).

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.
    P : sympy.Matrix, optional
        Pre-computed projection matrix. If None, computed internally.
    dirac_labels : list, optional
        Dirac labels for P rows. If None, computed with P.
    fock_data : FockGraphData, optional
        Pre-computed Fock graph data. If None, computed internally.

    Returns
    -------
    entries : list of (a, b, e, value) tuples
        Nonzero entries of V[a, b, e].
    N_dirac : int
        Number of Dirac states.
    V_fock : int
        Number of Fock nodes.
    E_fock : int
        Number of Fock edges (photon modes).
    """
    if P is None:
        P, dirac_labels, fock_states = build_projection_matrix(n_max)
    else:
        lat = GeometricLattice(max_n=n_max, topological_weights=False)
        fock_states = lat.states

    if fock_data is None:
        fock_data = build_fock_graph(n_max)

    N_dirac = P.rows
    V_fock = P.cols
    E_fock = fock_data.E
    edges = fock_data.edges

    entries = []

    for e_idx, (v1, v2) in enumerate(edges):
        # For each edge e connecting Fock nodes v1 and v2,
        # V[a, b, e] = P[a, v1] * P[b, v2] + P[a, v2] * P[b, v1]
        # (symmetrized because edge is undirected)
        for a in range(N_dirac):
            p_a_v1 = P[a, v1]
            p_a_v2 = P[a, v2]
            if p_a_v1 == 0 and p_a_v2 == 0:
                continue
            for b in range(N_dirac):
                p_b_v1 = P[b, v1]
                p_b_v2 = P[b, v2]

                val = p_a_v1 * p_b_v2 + p_a_v2 * p_b_v1
                val = sp.nsimplify(sp.expand(val), rational=False)
                if val != 0:
                    entries.append((a, b, e_idx, val))

    return entries, N_dirac, V_fock, E_fock


def vertex_tensor_to_matrices(
    entries: List[Tuple[int, int, int, sp.Expr]],
    N_dirac: int,
    E_fock: int,
) -> List[Matrix]:
    """Convert sparse vertex entries to a list of N_dirac x N_dirac matrices,
    one per photon edge.

    V_e[a, b] = V[a, b, e]

    Returns
    -------
    list of sympy.Matrix, length E_fock
        Each matrix is V_e (N_dirac x N_dirac).
    """
    V_mats = [sp_zeros(N_dirac, N_dirac) for _ in range(E_fock)]
    for a, b, e, val in entries:
        V_mats[e][a, b] = V_mats[e][a, b] + val
    return V_mats


# ---------------------------------------------------------------------------
# Step 3: One-loop vacuum polarization
# ---------------------------------------------------------------------------

def compute_vacuum_polarization(
    n_max: int,
    t: sp.Expr = Rational(0),
    exact: bool = True,
) -> Dict:
    """Compute the one-loop vacuum polarization tensor Pi[e1, e2].

    Pi[e1, e2] = Tr_{electron}[ V_{e1}^T . G_e . V_{e2} . G_e ]
               = sum_{a,b,c,d} V[a,b,e1] * G_e[b,c] * V[c,d,e2] * G_e[d,a]

    This is a FINITE MATRIX TRACE: no infinite sums, no regularization.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.
    t : sympy.Rational
        Dirac graph hopping parameter (default 0 = free diagonal propagator).
    exact : bool
        Use exact sympy arithmetic (suitable for n_max <= 2).

    Returns
    -------
    dict with keys:
        'Pi' : sympy.Matrix or np.ndarray
            Vacuum polarization tensor (E x E).
        'Pi_entries' : list of lists of str
            String representation of Pi entries.
        'Pi_trace' : str
            Trace of Pi.
        'Pi_eigenvalues' : list
            Eigenvalues of Pi.
        'pi_free' : bool
            True if all Pi entries are rational (no pi).
        'vertex_nonzero_count' : int
            Number of nonzero vertex entries.
        'vertex_sparsity' : float
            Fraction of zero vertex entries.
        'N_dirac' : int
        'E_fock' : int
        'n_max' : int
        'hopping_t' : str
    """
    # Build projection matrix
    P, dirac_labels, fock_states = build_projection_matrix(n_max)

    # Build Fock graph
    fock_data = build_fock_graph(n_max)
    E = fock_data.E

    # Build vertex tensor (sparse)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )

    # Build vertex matrices V_e (one per edge)
    V_mats = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    # Build electron propagator
    op = DiracGraphOperator(n_max=n_max, t=t)
    G_e, is_rational = electron_propagator(op, exact=exact)

    # Compute Pi[e1, e2] = Tr(V_{e1}^T . G_e . V_{e2} . G_e)
    if exact:
        Pi = sp_zeros(E, E)
        for e1 in range(E):
            for e2 in range(E):
                # Pi[e1, e2] = Tr(V_{e1}^T . G_e . V_{e2} . G_e)
                # = sum_{a,b,c,d} V_{e1}[b,a] * G_e[b,c] * V_{e2}[c,d] * G_e[d,a]
                # But this is: Tr(V1^T @ G @ V2 @ G)
                prod = V_mats[e1].T * G_e * V_mats[e2] * G_e
                Pi[e1, e2] = sp.nsimplify(prod.trace(), rational=False)
    else:
        # numpy path
        V_mats_np = [np.array(vm.tolist(), dtype=float) for vm in V_mats]
        G_np = np.array(G_e.tolist(), dtype=float) if isinstance(G_e, Matrix) else G_e
        Pi_np = np.zeros((E, E))
        for e1 in range(E):
            for e2 in range(E):
                prod = V_mats_np[e1].T @ G_np @ V_mats_np[e2] @ G_np
                Pi_np[e1, e2] = np.trace(prod)
        Pi = Pi_np

    # Analyze Pi
    if exact:
        Pi_entries = [[str(Pi[i, j]) for j in range(E)] for i in range(E)]
        Pi_trace = str(Pi.trace())
        Pi_eigs = list(Pi.eigenvals().keys())
        Pi_eigenvalues = [str(sp.nsimplify(ev, rational=False)) for ev in Pi_eigs]
        pi_free = _check_pi_free_matrix(Pi)
    else:
        Pi_entries = [[f"{Pi[i, j]:.10g}" for j in range(E)] for i in range(E)]
        Pi_trace = f"{np.trace(Pi):.10g}"
        Pi_eigenvalues = [f"{ev:.10g}" for ev in np.linalg.eigvalsh(Pi)]
        pi_free = None  # Cannot certify from floats

    total_possible = N_dirac * N_dirac * E_fock
    vertex_nnz = len(entries)
    vertex_sparsity = 1.0 - vertex_nnz / total_possible if total_possible > 0 else 1.0

    return {
        'Pi': Pi,
        'Pi_entries': Pi_entries,
        'Pi_trace': Pi_trace,
        'Pi_eigenvalues': Pi_eigenvalues,
        'pi_free': pi_free,
        'vertex_nonzero_count': vertex_nnz,
        'vertex_total_possible': total_possible,
        'vertex_sparsity': vertex_sparsity,
        'N_dirac': N_dirac,
        'V_fock': V_fock,
        'E_fock': E_fock,
        'n_max': n_max,
        'hopping_t': str(t),
    }


# ---------------------------------------------------------------------------
# Step 4: Pi-free certificate
# ---------------------------------------------------------------------------

def _check_pi_free_matrix(M: Matrix) -> bool:
    """Check that all entries of a sympy Matrix are algebraic (pi-free).

    Returns True if no entry contains sp.pi, sp.E, sp.EulerGamma,
    or other transcendental constants.
    """
    transcendental_atoms = (sp.pi, sp.E, sp.EulerGamma, sp.Catalan)
    for entry in M:
        expr = sp.sympify(entry)
        for atom in sp.preorder_traversal(expr):
            if atom in transcendental_atoms:
                return False
            if isinstance(atom, (sp.zeta, sp.log, sp.exp)):
                return False
    return True


def _check_rational_matrix(M: Matrix) -> bool:
    """Check that all entries are rational (not just algebraic)."""
    for entry in M:
        s = sp.nsimplify(entry, rational=True)
        if not isinstance(s, (sp.Rational, sp.Integer,
                              sp.core.numbers.Zero,
                              sp.core.numbers.One,
                              sp.core.numbers.NegativeOne)):
            return False
    return True


def vertex_pi_free_certificate(n_max: int) -> Dict:
    """Verify that vertex tensor and vacuum polarization are pi-free.

    Checks:
    1. Projection matrix P has algebraic (sqrt-rational) entries.
    2. Vertex tensor V has algebraic entries.
    3. Vacuum polarization Pi has rational entries (sqrt's cancel).
    4. No transcendental (pi, e, zeta) appears anywhere.

    Parameters
    ----------
    n_max : int

    Returns
    -------
    dict with boolean and diagnostic fields.
    """
    # Projection matrix
    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    P_pi_free = _check_pi_free_matrix(P)

    # Check if P entries are algebraic (sqrt of rational is OK)
    P_entries_algebraic = True
    for entry in P:
        expr = sp.sympify(entry)
        if expr.free_symbols:
            P_entries_algebraic = False
            break

    # Vertex tensor
    fock_data = build_fock_graph(n_max)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_vals = [val for _, _, _, val in entries]
    V_pi_free = all(
        not any(atom in (sp.pi, sp.E, sp.EulerGamma, sp.Catalan)
                for atom in sp.preorder_traversal(sp.sympify(v)))
        for v in V_vals
    ) if V_vals else True

    # Vacuum polarization
    vp = compute_vacuum_polarization(n_max, t=Rational(0), exact=True)

    return {
        'n_max': n_max,
        'P_pi_free': P_pi_free,
        'P_entries_algebraic': P_entries_algebraic,
        'P_shape': (P.rows, P.cols),
        'V_pi_free': V_pi_free,
        'V_nonzero_count': len(entries),
        'Pi_pi_free': vp['pi_free'],
        'Pi_rational': _check_rational_matrix(vp['Pi']) if isinstance(vp['Pi'], Matrix) else None,
        'Pi_trace': vp['Pi_trace'],
        'all_pass': P_pi_free and V_pi_free and (vp['pi_free'] is True),
    }


# ---------------------------------------------------------------------------
# Step 5: Vertex structure summary
# ---------------------------------------------------------------------------

def vertex_summary(n_max: int) -> Dict:
    """Comprehensive summary of vertex coupling structure at n_max.

    Returns
    -------
    dict with:
        - Projection matrix P dimensions and sparsity
        - Vertex tensor nonzero count and selection rules
        - Vertex entries grouped by edge
        - Selection rule analysis
    """
    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    fock_data = build_fock_graph(n_max)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )

    # P analysis
    P_nnz = sum(1 for entry in P if entry != 0)
    P_total = P.rows * P.cols

    # Group vertex entries by edge
    edge_groups: Dict[int, List[Dict]] = {}
    for a, b, e, val in entries:
        if e not in edge_groups:
            edge_groups[e] = []
        d_a = dirac_labels[a]
        d_b = dirac_labels[b]
        edge_groups[e].append({
            'a': a,
            'b': b,
            'a_label': f"(n={d_a.n_fock}, k={d_a.kappa}, 2mj={d_a.two_m_j})",
            'b_label': f"(n={d_b.n_fock}, k={d_b.kappa}, 2mj={d_b.two_m_j})",
            'value': str(val),
        })

    # Classify which electron transitions each photon edge mediates
    edge_analysis = []
    for e_idx in range(E_fock):
        v1, v2 = fock_data.edges[e_idx]
        s1 = fock_data.states[v1]
        s2 = fock_data.states[v2]
        dn = s2[0] - s1[0]
        dm = s2[2] - s1[2]

        edge_type = "T" if dn != 0 else "L"

        n_couplings = len(edge_groups.get(e_idx, []))

        # Analyze quantum number changes in the electron transitions
        delta_n_set = set()
        delta_l_set = set()
        delta_mj_set = set()
        for coupling in edge_groups.get(e_idx, []):
            a_idx = coupling['a']
            b_idx = coupling['b']
            da = dirac_labels[a_idx]
            db = dirac_labels[b_idx]
            delta_n_set.add(db.n_fock - da.n_fock)
            delta_l_set.add(db.l - da.l)
            delta_mj_set.add(db.two_m_j - da.two_m_j)

        edge_analysis.append({
            'edge_index': e_idx,
            'fock_nodes': (s1, s2),
            'edge_type': edge_type,
            'dn_photon': dn,
            'dm_photon': dm,
            'n_couplings': n_couplings,
            'electron_dn_values': sorted(delta_n_set),
            'electron_dl_values': sorted(delta_l_set),
            'electron_dmj_values': sorted(delta_mj_set),
        })

    return {
        'n_max': n_max,
        'N_dirac': N_dirac,
        'V_fock': V_fock,
        'E_fock': E_fock,
        'P_shape': (P.rows, P.cols),
        'P_nonzero': P_nnz,
        'P_total': P_total,
        'P_sparsity': 1.0 - P_nnz / P_total if P_total > 0 else 1.0,
        'vertex_nonzero': len(entries),
        'vertex_total': N_dirac * N_dirac * E_fock,
        'vertex_sparsity': 1.0 - len(entries) / (N_dirac * N_dirac * E_fock)
        if N_dirac * N_dirac * E_fock > 0 else 1.0,
        'edge_analysis': edge_analysis,
    }


# ---------------------------------------------------------------------------
# Step 6: Full analysis driver
# ---------------------------------------------------------------------------

def analyze_vertex_and_vp(
    n_max: int = 2,
    t: sp.Expr = Rational(0),
) -> Dict:
    """Run the full vertex coupling and vacuum polarization analysis.

    Parameters
    ----------
    n_max : int
        Graph truncation level (2 recommended for exact sympy).
    t : sympy.Rational
        Dirac graph hopping parameter.

    Returns
    -------
    dict with comprehensive analysis, suitable for JSON serialization.
    """
    summary = vertex_summary(n_max)
    vp = compute_vacuum_polarization(n_max, t=t, exact=True)
    cert = vertex_pi_free_certificate(n_max)

    # Remove the non-serializable sympy Matrix from vp
    vp_serial = {k: v for k, v in vp.items() if k != 'Pi'}

    return {
        'module': 'geovac.graph_qed_vertex',
        'n_max': n_max,
        'hopping_t': str(t),
        'vertex_structure': summary,
        'vacuum_polarization': vp_serial,
        'pi_free_certificate': cert,
    }


# ---------------------------------------------------------------------------
# CLI / script entry point
# ---------------------------------------------------------------------------

def _run_analysis_and_save(output_path: Optional[Path] = None) -> Dict:
    """Run full analysis at n_max=2 and save JSON."""

    if output_path is None:
        output_path = Path(__file__).parent.parent / "debug" / "data" / "gn4_graph_qed_vertex.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("Analyzing vertex coupling and vacuum polarization at n_max=2 ...")
    results = analyze_vertex_and_vp(n_max=2, t=Rational(0))

    # Serialize edge_analysis fock_nodes tuples
    for ea in results['vertex_structure']['edge_analysis']:
        ea['fock_nodes'] = [list(n) for n in ea['fock_nodes']]

    with output_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"Saved: {output_path}")
    return results


if __name__ == "__main__":
    _run_analysis_and_save()
