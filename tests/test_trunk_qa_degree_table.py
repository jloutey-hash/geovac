"""
TRUNK QA — Claim 6 (STRETCH): Paper 1 sec:III degree table (D_2s=0.854,
D_2p=3.416) from A = |T+|+|T-|+|L+|+|L-| (CG-magnitude weights).

This construction is NOT in production. We implement it here from Paper 1's
documented Biedenharn-Louck matrix elements:
  <l,m+1|L+|l,m> = sqrt((l-m)(l+m+1))
  <n+1,l|T+|n,l> = sqrt((n-l)(n+l+1)/4)
with A built by element-wise absolute value (undirected edges), and the
weighted degree D_ii = sum_j |A_ij|, exactly as Paper 1 sec:IX.B describes.
We then check whether it reproduces the paper's Table degrees.

FINDINGS (reported honestly):
  Faithful CG-magnitude construction at n_max=10:
    D_(1,0,0) = 0.707   D_(2,0,0) = 1.932   D_(2,1,0) = 3.828   D_(3,0,0)=2.957
  Paper 1 Table:
    D_(1,0,0) = 0.854   D_(2,0,0) = 0.854   D_(2,1,0) = 3.416   D_(3,0,0)=0.854
  Two structural mismatches:
   (1) The paper claims ALL s-states share D = 0.854 (n-independent). The
       faithful build gives s-state degree GROWING with n (0.71, 1.93, 2.96).
   (2) D_2p/D_2s = 4.0 exactly in the paper; faithful build gives ~1.98.
  No standard normalization (symmetric D^-1/2 A D^-1/2, binary) reproduces
  0.854 either. => The sec:III degree numbers do NOT follow from the paper's
  own stated operators. DOWNGRADE the sec:III degree table.

  What DOES reproduce from the construction (and is honestly backable):
   - p-states have strictly higher weighted degree than s-states at the same n
     (D_2p > D_2s), i.e. higher-l = higher connectivity. The QUALITATIVE
     differential-connectivity claim survives; the specific NUMBERS do not.
"""

from __future__ import annotations

import math

import numpy as np

import pytest

PAPER_D_2S = 0.854
PAPER_D_2P = 3.416


def build_A_cg(max_n: int):
    """A = |T+|+|T-|+|L+|+|L-| from Paper 1's documented matrix elements."""
    states = [(n, l, m) for n in range(1, max_n + 1)
              for l in range(n) for m in range(-l, l + 1)]
    idx = {s: i for i, s in enumerate(states)}
    N = len(states)
    A = np.zeros((N, N))
    for (n, l, m) in states:
        i = idx[(n, l, m)]
        if m < l:                                   # L+/L- (abs, both directions)
            j = idx[(n, l, m + 1)]
            a = math.sqrt((l - m) * (l + m + 1))
            A[i, j] += a
            A[j, i] += a
        if (n + 1, l, m) in idx:                    # T+/T- (abs, both directions)
            j = idx[(n + 1, l, m)]
            a = math.sqrt((n - l) * (n + l + 1) / 4.0)
            A[i, j] += a
            A[j, i] += a
    return states, idx, A


def _degrees(max_n: int):
    states, idx, A = build_A_cg(max_n)
    D = A.sum(axis=1)
    return idx, D


def test_qualitative_differential_connectivity_holds():
    """The QUALITATIVE claim: p-states have higher weighted degree than
    s-states at the same n. This reproduces and is BACKED."""
    idx, D = _degrees(10)
    assert D[idx[(2, 1, 0)]] > D[idx[(2, 0, 0)]], "expect D_2p > D_2s"
    assert D[idx[(3, 1, 0)]] > D[idx[(3, 0, 0)]], "expect D_3p > D_3s"


def test_paper_specific_degree_NUMBERS_do_not_reproduce():
    """The SPECIFIC numbers 0.854 / 3.416 are NOT produced by the paper's
    own documented construction. This is the DOWNGRADE finding."""
    idx, D = _degrees(10)
    d2s = D[idx[(2, 0, 0)]]
    d2p = D[idx[(2, 1, 0)]]
    # Faithful build: ~1.93 / ~3.83, far from 0.854 / 3.416.
    assert abs(d2s - PAPER_D_2S) > 0.5, f"D_2s={d2s:.3f} unexpectedly near 0.854"
    assert abs(d2p - PAPER_D_2P) > 0.3, f"D_2p={d2p:.3f} unexpectedly near 3.416"


def test_paper_s_state_constancy_claim_is_false():
    """Paper 1's Table lists D_1s = D_2s = D_3s = 0.854 (n-independent).
    The faithful construction gives s-state degree GROWING with n.
    The constancy claim is falsified by the paper's own operators."""
    idx, D = _degrees(10)
    d1s = D[idx[(1, 0, 0)]]
    d2s = D[idx[(2, 0, 0)]]
    d3s = D[idx[(3, 0, 0)]]
    # If the paper were right these would be equal. They are not.
    assert not (abs(d1s - d2s) < 0.05 and abs(d2s - d3s) < 0.05), (
        f"s-state degrees claimed constant but are {d1s:.3f},{d2s:.3f},{d3s:.3f}"
    )
    assert d1s < d2s < d3s, (
        f"s-state degree should grow with n: {d1s:.3f},{d2s:.3f},{d3s:.3f}"
    )


def test_ratio_is_not_exactly_four():
    """Paper's D_2p/D_2s = 4.0 exactly; faithful build gives ~1.98."""
    idx, D = _degrees(10)
    ratio = D[idx[(2, 1, 0)]] / D[idx[(2, 0, 0)]]
    assert abs(ratio - 4.0) > 1.5, f"ratio {ratio:.3f} unexpectedly near 4.0"


@pytest.mark.xfail(reason="Paper 1 sec:III degree numbers 0.854/3.416 are not "
                          "reproduced by the paper's own stated operators "
                          "(faithful build gives 1.93/3.83). DOWNGRADE.")
def test_paper_degree_table_reproduces_DOWNGRADE_MARKER():
    """xfail marker: if someone supplies the missing normalization that yields
    0.854/3.416, this XPASSes and flags the construction for review."""
    idx, D = _degrees(10)
    assert abs(D[idx[(2, 0, 0)]] - PAPER_D_2S) < 0.05
    assert abs(D[idx[(2, 1, 0)]] - PAPER_D_2P) < 0.05
