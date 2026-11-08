"""Archived 2026-11-08: the ST-SU3 Bargmann graph / plaquette census.

These pinned the census quoted by the archived Paper 25 Sec. VII.A (7,765
plaquettes at N_max = 3).  No live paper cites the census since Paper 25 was
retired into Paper 30 (2026-11-07); the N_max = 5 graph counts that Paper 24
still quotes stay live in tests/test_su3_wilson_s5.py
(TestPaper24S5GraphCounts).  Reference data: debug/data/st_su3_plaquettes.json.
"""
from collections import Counter

import numpy as np
import pytest

from geovac.su3_wilson_s5 import bargmann_adjacency_dense, enumerate_plaquettes


def _graph_counts(N_max):
    """(V, E, c, beta_1) of the Bargmann graph at N_max."""
    A = bargmann_adjacency_dense(N_max)
    V = A.shape[0]
    E = int(A.sum()) // 2
    seen = set()
    c = 0
    for s in range(V):
        if s in seen:
            continue
        c += 1
        stack = [s]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.extend(int(y) for y in np.nonzero(A[x])[0] if y not in seen)
    return V, E, c, E - V + c


def test_graph_counts():
    # (V, E, c, beta_1) per debug/data/st_su3_plaquettes.json
    assert _graph_counts(2) == (10, 15, 1, 6)
    assert _graph_counts(3) == (20, 42, 1, 23)
    assert _graph_counts(4) == (35, 90, 1, 56)


def test_plaquette_census_nmax2():
    A = bargmann_adjacency_dense(2)
    pl = enumerate_plaquettes(A, max_length=8, both_orientations=False)
    cnt = Counter(len(P) for P in pl)
    assert dict(cnt) == {4: 15, 6: 21, 8: 146}
    assert len(pl) == 182


def test_plaquette_census_nmax3():
    """N_max = 3, lengths <= 8: 88 + 526 + 7151 = 7,765 plaquettes."""
    A = bargmann_adjacency_dense(3)
    pl = enumerate_plaquettes(A, max_length=8, both_orientations=False)
    cnt = Counter(len(P) for P in pl)
    assert dict(cnt) == {4: 88, 6: 526, 8: 7151}
    assert len(pl) == 7765


def test_plaquette_census_nmax4_short():
    A = bargmann_adjacency_dense(4)
    pl = enumerate_plaquettes(A, max_length=6, both_orientations=False)
    cnt = Counter(len(P) for P in pl)
    assert dict(cnt) == {4: 272, 6: 2934}


@pytest.mark.slow
def test_plaquette_census_nmax4_full():
    """N_max = 4 up to L = 8 (about 40 s)."""
    A = bargmann_adjacency_dense(4)
    pl = enumerate_plaquettes(A, max_length=8, both_orientations=False)
    cnt = Counter(len(P) for P in pl)
    assert dict(cnt) == {4: 272, 6: 2934, 8: 64021}
    assert len(pl) == 67227
