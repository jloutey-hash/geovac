"""Guards for claims re-tiered or newly printed in the /qa trunk FULL #11
remediation (2026-09-30).  Each test names the wrong answer it rejects.

  1. Paper 7: the variationally optimised single-exponent helium energy is
     k* = 27/16, E* = -729/256 exactly; the retired -11/8 is half the TOTAL
     energy at k = Z (not a halved one-body sum, which gives -3/4), and Paper 7
     must say so.  Rejects any non-optimised k and the wrong error story.
  2. Paper 1: on the binary production lattice the node degree is NOT
     monotone in l -- it depends on m and on which shell neighbours exist.
     Rejects the retired "higher l => more connections" claim.
  3. geovac.lattice_spectrum.lambda_max_from_operator refuses an operator
     with l-changing (off-block) weight, where its blockwise top would be
     wrong.  Rejects silently returning the blockwise value.
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sps
import sympy as sp

from geovac.lattice import GeometricLattice
from geovac.lattice_spectrum import lambda_max_from_operator


# ---------------------------------------------------------------- 1. E*
def _he_single_exponent_energy(k: sp.Expr, Z: int) -> sp.Expr:
    """Two 1s electrons with exponent k in nuclear charge Z:
    kinetic 2*(k^2/2), attraction 2*(-Z k), repulsion F0(1s,1s) = 5k/8."""
    return k**2 - 2 * Z * k + sp.Rational(5, 8) * k


def test_paper7_he_optimised_single_exponent_is_minus_729_over_256():
    k = sp.symbols("k", positive=True)
    E = _he_single_exponent_energy(k, 2)
    kstar = sp.solve(sp.diff(E, k), k)
    assert kstar == [sp.Rational(27, 16)]
    Estar = sp.simplify(E.subs(k, kstar[0]))
    assert Estar == sp.Rational(-729, 256)
    # n_max = 1 fixed-k point (k = Z) is -11/4, above E*: variational order
    assert E.subs(k, 2) == sp.Rational(-11, 4)
    assert E.subs(k, 2) > Estar
    # The retired printing's -11/8 is HALF THE TOTAL energy at k = Z (corrected
    # DELTA #24: an earlier version of this test, and Paper 7, attributed it to
    # halving the one-body sum, which gives -3/4 instead).
    assert E.subs(k, 2) / 2 == sp.Rational(-11, 8)
    one_body_halved = (k**2 - 2 * 2 * k) / 2 + sp.Rational(5, 8) * k
    assert one_body_halved.subs(k, 2) == sp.Rational(-3, 4)


def test_paper7_prints_the_exact_he_values_and_the_right_error_story():
    """Paper 7 must carry -729/256, 27/16 and -11/4, and must not print the
    retired -11/8 (half the k = Z total) or the wrong 'halved the one-body sum'
    story (the edit bar removed the retraction clause; history is in CHANGELOG)."""
    import os
    p = os.path.join(os.path.dirname(__file__), "..", "papers", "group3_foundations",
                     "Paper_7_Dimensionless_Vacuum.tex")
    t = open(p, encoding="utf-8").read()
    assert "-729/256" in t and "27/16" in t and "-11/4" in t
    assert "-11/8" not in t
    assert "halved the one-body sum" not in t


# ------------------------------------------------------- 2. degree vs l
def _degrees(n_max: int) -> dict:
    lat = GeometricLattice(max_n=n_max)
    A = lat.adjacency
    A = abs(A.toarray()) if hasattr(A, "toarray") else abs(np.asarray(A))
    d = A.sum(axis=1)
    return {tuple(s): float(d[i]) for i, s in enumerate(lat.states)}


def test_paper1_degree_not_monotone_in_l():
    deg = _degrees(4)
    # 2p m=0 has degree 3 > 2s = 2 = 2p m=+-1
    assert deg[(2, 0, 0)] == 2 and deg[(2, 1, 0)] == 3
    assert deg[(2, 1, 1)] == deg[(2, 1, -1)] == 2
    # interior shell n=3 at n_max=4: max degree per l is 2, 4, 3 for l = 0, 1, 2
    mx = {l: max(v for (n, ll, m), v in deg.items() if n == 3 and ll == l) for l in (0, 1, 2)}
    assert mx == {0: 2.0, 1: 4.0, 2: 3.0}
    mn = {l: min(v for (n, ll, m), v in deg.items() if n == 3 and ll == l) for l in (0, 1, 2)}
    assert mn == {0: 2.0, 1: 3.0, 2: 2.0}   # the printed ranges 2, 3-4, 2-3
    # the top shell (n = n_max = 4) is lowered by truncation: no T+ neighbour
    assert max(v for (n, l, m), v in deg.items() if n == 4 and l == 1) < mx[1]
    # the retired claim (degree increasing in l) is false here
    assert not (mx[0] <= mx[1] <= mx[2])


def test_paper1_degree_plateau_in_larger_interior_shell():
    """Paper 1: in an interior shell the max degree rises from l=0, plateaus
    from l=1 to l=n-2, and drops at l=n-1 (n=5 at n_max=6: 2, 4, 4, 4, 3).
    Rejects the retired 'peaks at l=1' and 'increasing in l' readings."""
    deg = _degrees(6)
    mx = [max(v for (n, l, m), v in deg.items() if n == 5 and l == ll) for ll in range(5)]
    assert mx == [2.0, 4.0, 4.0, 4.0, 3.0]


def test_paper1_cg_degree_not_monotone_in_l_once_m_counted():
    """Paper 1 Sec. III: on the CG construction the m=0 nodes satisfy
    D_2p > D_2s and D_3p > D_3s, but degree is not monotone in l once m != 0
    nodes are counted (3d m=+-2 = 3.225 < 3p m=0 = 5.410)."""
    import importlib.util
    import os
    spec = importlib.util.spec_from_file_location(
        "deg_table", os.path.join(os.path.dirname(__file__), "test_trunk_qa_degree_table.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    states, idx, A = mod.build_A_cg(5)
    d = A.sum(axis=1)
    D = {s: float(d[idx[s]]) for s in states}
    assert D[(2, 1, 0)] > D[(2, 0, 0)] and D[(3, 1, 0)] > D[(3, 0, 0)]
    assert abs(D[(3, 2, 2)] - 3.225) < 1e-3 and abs(D[(3, 1, 0)] - 5.410) < 1e-3
    assert D[(3, 2, 2)] < D[(3, 1, 0)]


# --------------------------------------------- 3. off-block refusal
def test_lambda_max_from_operator_refuses_off_block_weight():
    lat = GeometricLattice(max_n=4)
    A = lat.adjacency
    A = sps.csr_matrix(A) if not sps.issparse(A) else A.tocsr()
    D = sps.diags(np.asarray(abs(A).sum(axis=1)).ravel())
    L = (D - A).tocsr()
    states = lat.states
    good = lambda_max_from_operator(L, states)
    dense_top = float(np.max(np.linalg.eigvalsh(L.toarray())))
    assert abs(good - dense_top) < 1e-10
    # couple an l=0 node to an l=1 node: an l-changing edge
    i = [k for k, s in enumerate(states) if tuple(s) == (2, 0, 0)][0]
    j = [k for k, s in enumerate(states) if tuple(s) == (2, 1, 0)][0]
    Lb = L.tolil()
    Lb[i, j] = Lb[j, i] = -1.0
    Lb[i, i] += 1.0
    Lb[j, j] += 1.0
    with pytest.raises(ValueError, match="off-block"):
        lambda_max_from_operator(Lb.tocsr(), states)
