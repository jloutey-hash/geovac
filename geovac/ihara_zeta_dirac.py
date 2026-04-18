"""
Ihara zeta of the Dirac-S^3 graph (Track RH-C).
================================================

This module builds a **spin-ful** adjacency graph on the Dirac-orbital
node set (n_fock, kappa, m_j) and delegates all Ihara-zeta / Ramanujan
machinery to :mod:`geovac.ihara_zeta` (Track RH-A).

Motivation
----------
Track RH-A (``debug/ihara_zeta_memo.md``) established that the scalar
S^3 Coulomb graph and the S^5 Bargmann-Segal graph are both
Ramanujan in the weak Kotani-Sunada sense. Paper 29 promotes that
observation.  Track RH-C asks: **does the spin-ful graph inherit
the Ramanujan property, and does the κ-structure visibly split the
zero set?**

Node set
--------
Dirac labels (n_fock, κ, m_j) as built by
:func:`geovac.dirac_matrix_elements.iter_dirac_labels`.  With
n_fock ≤ n_max one gets:

    n_max = 1 : only (1, κ=-1, m_j=±1/2) → 2 nodes     (one s_{1/2} orbital)
    n_max = 2 : adds l=0 at n=2 (2 nodes) and p_{1/2}, p_{3/2}
                (2 + 4 nodes) → 10 total nodes
    n_max = 3 : adds n=3, l=0,1,2 → 28 total nodes

This is the physical Dirac-orbital set at level n_max: each scalar
(n, l, m) state lifts to a full (j, m_j) multiplet.  |κ| = j + 1/2.
The count at n_max=2 is 2 + 8 = 10 (l=0:2, p1/2:2, p3/2:4).

Adjacency rules (EXPLICIT — the key design choice)
--------------------------------------------------
This module computes Ihara zetas for TWO physically reasonable
spin-ful adjacency rules and compares them, because there is no
canonical scalar-S^3-equivalent Dirac adjacency.  Both rules lift
the scalar (n ↔ n±1, m ↔ m±1) Fock Coulomb adjacency to spinors
in a different way:

- **Rule A (spinor lift, "scalar-analog"):**
  Edge iff either (i) Δn_fock=±1, Δκ=0, Δm_j=0, or
  (ii) Δn_fock=0, Δκ=0, Δ(2 m_j)=±2 (i.e. Δm_j = ±1).
  This is the straightforward lift of the scalar ladders
  n ↔ n±1 and m ↔ m±1 to the spinor basis at fixed κ; κ itself
  is never touched.

- **Rule B (dipole, E1 photon-emission):**
  Edge iff Δn_fock ∈ {-1, 0, +1} and the transition is
  dipole-allowed:

    (i) parity flip (Δl = ±1, where l = kappa_to_l(κ));
    (ii) |Δj| ≤ 1, not (j = j' = 0) [automatic here since j ≥ 1/2];
    (iii) |Δ m_j| ≤ 1.

  This matches the standard atomic E1 selection rule and is the
  relativistic analog of the m ↔ m±1 / n ↔ n±1 ladders once one
  passes to (κ, m_j).  Self-loops (Δn=Δκ=Δm_j=0) excluded.  No
  spin-flip-only edges at fixed n,l.

Rule A preserves κ, so the graph **splits into |κ|-sectors** the
same way the scalar Coulomb graph splits into l-shells.  Rule B
mixes κ via l = l±1 parity flips and induces a single connected
component beyond n_max=1.  Neither rule is obviously canonical;
reporting both is itself a scientific finding.

Transcendental taxonomy (Paper 18 / Paper 24)
---------------------------------------------
Every quantity built here is a combinatorial graph invariant.
The Ihara zeta depends only on 0/1 connectivity, so the choice of
spinor convention (pair-diagonal vs full-Gaunt, Paper 22) does not
enter.  The output is a rational-coefficient polynomial in the
symbolic variable s, exactly as for the scalar graphs in RH-A.

References
----------
- R. Camporesi and A. Higuchi, "On the eigenfunctions of the Dirac
  operator on spheres and real hyperbolic spaces", J. Geom. Phys.
  20 (1996) 1-18.  (Spectrum |λ_n| = n + 3/2, degeneracies.)
- R. Szmytkowski, "Recurrence and differential relations for
  spherical spinors", J. Math. Chem. 42 (2007) 397-445.  (κ-basis
  angular matrix elements; dipole selection rules in (κ, m_j).)
- GeoVac Paper 14 §V (Tier 2 spin-ful composed qubit encoding).
- GeoVac Paper 23 (Dirac graph construction context).
- GeoVac Paper 29 (scalar Ramanujan observation, April 2026).
- ``geovac.ihara_zeta`` — Track RH-A machinery (reused verbatim).
- ``geovac.dirac_s3`` — Camporesi-Higuchi spectrum, π-free certificate.
- ``geovac.dirac_matrix_elements`` — DiracLabel, iter_dirac_labels.

Author: GeoVac Track RH-C (Sprint "hey-buddy-we-need-crystalline-sprout"),
        April 2026.
"""

from __future__ import annotations

from typing import Dict, List, Literal, Optional, Sequence, Tuple

import numpy as np
import sympy as sp

from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels, kappa_to_l
from geovac.ihara_zeta import (
    functional_equation_report,
    hashimoto_matrix,
    ihara_zeta_bass,
    is_ramanujan,
    zeta_zeros_from_hashimoto,
)


AdjacencyRule = Literal["A", "B"]


__all__ = [
    "build_dirac_s3_graph",
    "dirac_s3_ihara_zeta",
    "describe_adjacency_rule",
    "per_kappa_components",
]


# ---------------------------------------------------------------------------
# Adjacency rule descriptions
# ---------------------------------------------------------------------------

def describe_adjacency_rule(rule: AdjacencyRule) -> str:
    """Return a human-readable description of the adjacency rule."""
    if rule == "A":
        return (
            "Rule A (spinor lift / scalar-analog): "
            "Δn_fock=±1 at fixed (κ, m_j), OR Δm_j=±1 at fixed (n_fock, κ). "
            "κ is preserved. Graph splits into κ-sectors, analogous to the "
            "l-shell splitting of the scalar S^3 Coulomb graph."
        )
    if rule == "B":
        return (
            "Rule B (dipole / E1 selection rule): "
            "Δn_fock ∈ {-1, 0, +1}, parity flip (Δl = ±1 where l = kappa_to_l(κ)), "
            "|Δj| ≤ 1, |Δm_j| ≤ 1, excluding the trivial (Δn=Δκ=Δm_j=0) "
            "self-loop. Mixes κ via l ↔ l±1, producing a connected graph "
            "beyond n_max=1."
        )
    raise ValueError(f"unknown adjacency rule {rule!r} (expected 'A' or 'B')")


# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------

def _dirac_labels(n_max: int) -> List[DiracLabel]:
    """All DiracLabel's with n_fock ≤ n_max, in canonical order."""
    if n_max < 1:
        raise ValueError("n_max must be ≥ 1")
    return list(iter_dirac_labels(n_max))


def _edge_rule_A(a: DiracLabel, b: DiracLabel) -> bool:
    """Rule A: scalar-analog ladders at fixed κ."""
    if a.kappa != b.kappa:
        return False
    dn = a.n_fock - b.n_fock
    dmj2 = a.two_m_j - b.two_m_j
    # n-ladder: Δn=±1, Δm_j = 0
    if abs(dn) == 1 and dmj2 == 0:
        return True
    # m-ladder: Δn = 0, Δ(2 m_j) = ±2 (i.e. Δm_j = ±1)
    if dn == 0 and abs(dmj2) == 2:
        return True
    return False


def _edge_rule_B(a: DiracLabel, b: DiracLabel) -> bool:
    """Rule B: E1 dipole selection.

    Δn_fock ∈ {-1, 0, +1}, Δl = ±1, |Δj| ≤ 1, |Δm_j| ≤ 1, not a self-loop.
    """
    dn = a.n_fock - b.n_fock
    if abs(dn) > 1:
        return False
    # parity flip: Δl = ±1
    la = kappa_to_l(a.kappa)
    lb = kappa_to_l(b.kappa)
    if abs(la - lb) != 1:
        return False
    # |Δj| ≤ 1.  2j = 2|κ|-1 so Δ(2j) ∈ {-2, 0, +2}.
    dtwo_j = (2 * abs(a.kappa) - 1) - (2 * abs(b.kappa) - 1)
    if abs(dtwo_j) > 2:
        return False
    # |Δm_j| ≤ 1 (in 2·m_j units, ≤ 2)
    dmj2 = a.two_m_j - b.two_m_j
    if abs(dmj2) > 2:
        return False
    # Exclude self-loop (should not happen if labels are distinct, but guard).
    if a == b:
        return False
    return True


_RULE_DISPATCH = {"A": _edge_rule_A, "B": _edge_rule_B}


def build_dirac_s3_graph(
    n_max: int,
    adjacency_rule: AdjacencyRule,
) -> Tuple[np.ndarray, List[DiracLabel], np.ndarray, str]:
    """Build the Dirac-S^3 spin-ful adjacency graph.

    Parameters
    ----------
    n_max : int
        Maximum Fock principal quantum number to include (n_fock = 1..n_max).
    adjacency_rule : {"A", "B"}
        Which edge rule to apply.  See ``describe_adjacency_rule``.

    Returns
    -------
    adjacency : numpy.ndarray of shape (V, V), dtype=int, symmetric 0/1
        Adjacency matrix of the spin-ful graph.  Zero diagonal.
    labels : list of DiracLabel
        Canonical-ordered node labels, len = V.
    degree_sequence : numpy.ndarray of shape (V,), dtype=int
        Degree of each node (sum of each row of adjacency).
    description : str
        Human-readable sentence summarising the rule and size.
    """
    if adjacency_rule not in _RULE_DISPATCH:
        raise ValueError(
            f"adjacency_rule must be 'A' or 'B', got {adjacency_rule!r}")
    labels = _dirac_labels(n_max)
    V = len(labels)
    edge_test = _RULE_DISPATCH[adjacency_rule]
    A = np.zeros((V, V), dtype=int)
    # pair enumeration
    for i in range(V):
        for j in range(i + 1, V):
            if edge_test(labels[i], labels[j]):
                A[i, j] = 1
                A[j, i] = 1
    deg = A.sum(axis=1).astype(int)
    desc = (
        f"Dirac-S^3 graph, rule {adjacency_rule}, n_max={n_max}: "
        f"V={V}, E={int(A.sum()) // 2}, "
        f"max deg={deg.max() if V else 0}, min deg={deg.min() if V else 0}."
    )
    return A, labels, deg, desc


# ---------------------------------------------------------------------------
# Per-κ decomposition
# ---------------------------------------------------------------------------

def per_kappa_components(
    labels: Sequence[DiracLabel],
    A: np.ndarray,
) -> Dict[int, Tuple[np.ndarray, List[int]]]:
    """Extract per-κ sub-blocks of the adjacency.

    For each distinct κ value, returns the induced subgraph on the
    nodes with that κ (using only edges between those nodes).  For
    Rule A this is the full connected decomposition of the graph.
    For Rule B the per-κ induced subgraphs are disconnected from
    each other (κ is mixed by the dipole edges), but per-κ blocks
    are still well-defined combinatorially.

    Returns
    -------
    dict mapping κ → (A_sub, indices)
        A_sub : sub-adjacency of that κ sector.
        indices : list of original indices (into ``labels``) for this sector.
    """
    kappa_values = sorted(set(lab.kappa for lab in labels))
    out: Dict[int, Tuple[np.ndarray, List[int]]] = {}
    for k in kappa_values:
        idx = [i for i, lab in enumerate(labels) if lab.kappa == k]
        if len(idx) == 0:
            continue
        sub = A[np.ix_(idx, idx)].copy()
        out[k] = (sub, idx)
    return out


# ---------------------------------------------------------------------------
# Bit-exact π-free certificate: adjacency spectrum is rational
# ---------------------------------------------------------------------------

def _adjacency_spectrum_is_rational(A: np.ndarray) -> Tuple[bool, List[sp.Rational]]:
    """Verify that adjacency eigenvalues are exact rationals/algebraic.

    For an integer matrix, eigenvalues are roots of an integer-coefficient
    polynomial (the characteristic polynomial).  They are "π-free" in
    the Paper 24 sense iff they are all rational.  This is only true
    for very restricted topologies (e.g. complete multipartite); for
    most graphs the eigenvalues are algebraic but NOT rational.  For
    the purposes of a π-free certificate on an integer-adjacency
    graph, we report:

      - strict rationality: all eigenvalues are in ℚ (very strong);
      - weak π-freeness: the characteristic polynomial has integer
        coefficients (always true for integer A), so eigenvalues are
        algebraic over ℚ (no π).  This is weaker but ALWAYS holds.

    Returns (strict_rational, sympy_eigenvalues_if_rational).
    """
    M = sp.Matrix(A.tolist())
    charpoly = M.charpoly()  # Poly over s
    # Roots in sympy exact form
    roots = sp.roots(charpoly, multiple=False)
    if roots is None or len(roots) == 0:
        return False, []
    is_all_rational = all(sp.ask(sp.Q.rational(r)) is True for r in roots)
    if is_all_rational:
        return True, [sp.Rational(r) for r in roots for _ in range(int(roots[r]))]
    return False, []


def _charpoly_is_integer_coefficient(A: np.ndarray) -> bool:
    """Weaker π-free certificate: characteristic polynomial has integer coeffs.

    For any integer matrix this is trivially true; we verify it explicitly
    via sympy.
    """
    if A.size == 0:
        return True
    M = sp.Matrix(A.tolist())
    # M.charpoly(gen) returns a PurePoly in the requested generator directly,
    # avoiding the "lambda symbol leaking via as_expr" bug.
    s_local = sp.symbols("s")
    cp = M.charpoly(s_local)
    return all(getattr(c, "is_Integer", False) for c in cp.all_coeffs())


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def dirac_s3_ihara_zeta(
    n_max: int,
    adjacency_rule: AdjacencyRule,
    *,
    factor: bool = True,
) -> Dict:
    """Compute the full Ihara zeta analysis for a Dirac-S^3 spin-ful graph.

    Reuses the RH-A machinery verbatim; does not re-implement Ihara-Bass
    or the Hashimoto eigenvalue solve.

    Parameters
    ----------
    n_max : int
        Maximum Fock principal quantum number.
    adjacency_rule : {"A", "B"}
        Which spin-ful adjacency to use.
    factor : bool, default True
        If True, attempt to factor the Ihara-Bass polynomial over ℚ
        and return the factored form.  Factoring may be slow for
        large n_max; set False to skip.

    Returns
    -------
    dict with keys:
        n_max, adjacency_rule, description,
        V, E, c, r_betti1, degree_sequence,
        zeta_inverse_expanded (sympy),
        zeta_inverse_factored (sympy, if factor=True),
        hashimoto_spectrum (numpy, sorted by |μ|),
        spectral_radius, max_abs_nontrivial,
        q_max, sqrt_q_max,
        ramanujan_verdict (bool), ramanujan_deviation (float),
        ramanujan_explanation,
        zeros (numpy complex, sorted by |s|),
        per_kappa_sizes (dict κ → (V_κ, E_κ, components)),
        charpoly_integer_coefficient (bool),
        adjacency_spectrum_all_rational (bool),
        functional_equation_report (dict).

    Notes
    -----
    - Bit-exact rational arithmetic is used for the Bass polynomial;
      the Hashimoto eigenvalue solve is numpy-numeric (as in RH-A).
    - All keys are JSON-serialisable where needed (see `dirac_s3_data_dict`).
    """
    # Build graph
    A, labels, deg, desc = build_dirac_s3_graph(n_max, adjacency_rule)
    V = int(A.shape[0])
    E = int(A.sum()) // 2

    # Reuse RH-A machinery
    s = sp.symbols("s")
    zeta_sym = ihara_zeta_bass(A)  # sympy rational function
    # sp.cancel puts 1/zeta in reduced (numerator/denominator) form and, when
    # the result is a polynomial, collapses it cleanly — so sp.Poly(...)
    # downstream works without tripping over un-simplified rational forms.
    zeta_inv = sp.cancel(1 / zeta_sym)
    zeta_inv_factored = None
    if factor:
        try:
            zeta_inv_factored = sp.factor(zeta_inv)
        except Exception:  # pragma: no cover
            zeta_inv_factored = None

    # Hashimoto spectrum and derived quantities
    T = hashimoto_matrix(A).astype(float)
    if T.size == 0:
        hashi_spec = np.array([], dtype=complex)
    else:
        hashi_spec = np.linalg.eigvals(T)
    # Sort by absolute value descending for readability
    order = np.argsort(-np.abs(hashi_spec)) if hashi_spec.size else np.array([], dtype=int)
    hashi_spec_sorted = hashi_spec[order] if order.size else hashi_spec
    if hashi_spec.size:
        mags = np.abs(hashi_spec_sorted)
        rho = float(mags[0])
        non_trivial = mags < rho - 1e-9
        mu_nt_max = float(mags[non_trivial].max()) if non_trivial.any() else 0.0
    else:
        rho = 0.0
        mu_nt_max = 0.0

    if A.sum() == 0:
        # Zero-edge graph: trivially Ramanujan (no non-trivial eigenvalues
        # to bound). Upstream is_ramanujan cannot handle an empty Hashimoto
        # spectrum, so we short-circuit here.
        is_ram, dev, explanation = True, 0.0, "trivial: empty graph (no edges)"
    else:
        is_ram, dev, explanation = is_ramanujan(A)
    q_max = int(deg.max()) - 1 if V else 0
    sqrt_q_max = float(np.sqrt(max(q_max, 0)))

    # Zeros of zeta = 1/μ for non-zero μ
    if hashi_spec.size:
        zeros = zeta_zeros_from_hashimoto(A)
        # sort by |s|
        zeros = zeros[np.argsort(np.abs(zeros))]
    else:
        zeros = np.array([], dtype=complex)

    # Per-κ structure
    per_k = per_kappa_components(labels, A)
    per_kappa_sizes: Dict[int, Dict] = {}
    for k, (sub, idx) in per_k.items():
        Vk = sub.shape[0]
        Ek = int(sub.sum()) // 2
        # Count components within this κ-induced subgraph
        from geovac.ihara_zeta import _count_components
        ck = _count_components(sub) if Vk else 0
        per_kappa_sizes[int(k)] = {
            "V": Vk,
            "E": Ek,
            "components": ck,
        }

    # π-free certificates
    all_rational, _rational_evals = _adjacency_spectrum_is_rational(A)
    charpoly_int = _charpoly_is_integer_coefficient(A)

    # Functional equation report
    fe = functional_equation_report(A)

    return {
        "n_max": n_max,
        "adjacency_rule": adjacency_rule,
        "description": desc,
        "V": V,
        "E": E,
        "c": int(fe["c"]),
        "r_betti1": int(fe["r_betti1"]),
        "degree_sequence": deg.tolist(),
        "zeta_inverse_expanded": zeta_inv,
        "zeta_inverse_factored": zeta_inv_factored,
        "hashimoto_spectrum": hashi_spec_sorted,
        "spectral_radius": rho,
        "max_abs_nontrivial": mu_nt_max,
        "q_max": q_max,
        "sqrt_q_max": sqrt_q_max,
        "ramanujan_verdict": bool(is_ram),
        "ramanujan_deviation": float(dev),
        "ramanujan_explanation": explanation,
        "zeros": zeros,
        "per_kappa_sizes": per_kappa_sizes,
        "charpoly_integer_coefficient": bool(charpoly_int),
        "adjacency_spectrum_all_rational": bool(all_rational),
        "functional_equation_report": fe,
        # Forward the node label triples for reproducibility
        "labels_triple": [(lab.n_fock, lab.kappa, lab.two_m_j) for lab in labels],
    }
