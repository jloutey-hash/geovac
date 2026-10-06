r"""Light operator algebra for en_factors.py (no grid / kernel imports; unit-tested fast).

Backs (Paper 19, LiH e-n extension of the CI-relaxed marriage): the symbolic bookkeeping of the
electron-nucleus columns -- sums of products of one-electron weight tokens (per electron) times at
most two pair-kernel edges -- and the closed-form Laplacian of a w-weighted geminal (g_lap), whose
grad w . grad f cross term is decomposed into one-electron tokens (u, dzw, z, and the circulation
+-1 tokens drwP/drwM, rhoP/rhoM).  Token -> real-space function (exact, phi = azimuth):
    w, lapw = the switch and its Laplacian;  u = rho d_rho w + z d_z w;  dzw = d_z w;  z;
    drwP/drwM = d_rho w e^{+-i phi};  rhoP/rhoM = rho e^{+-i phi}.
Edge keys: 'f' = e^{-g r}, 'Y' = e^{-g r}/r, 'coul' = 1/r.
"""
from __future__ import annotations

from itertools import combinations
from typing import Dict, Optional, Sequence, Tuple

PAIRS6 = [(1, 2), (3, 4), (1, 3), (1, 4), (2, 3), (2, 4)]   # == enumerators.PAIRS6
ELEC = (1, 2, 3, 4)


# =========================================================================== #
# operator algebra: sums of products of one-electron weights (tokens per electron)
# =========================================================================== #
# An operator is {term: coef}, term = 4-tuple (electrons 1..4) of SORTED token tuples, token in
# {'w','lapw'}; () = 1.  E.g. W1 = {(('w',),(),(),()): 1, ((),('w',),(),()): 1, ...}.
Op = Dict[Tuple[Tuple[str, ...], ...], float]
ONE: Op = {((), (), (), ()): 1.0}


def _single(e: int, toks: Tuple[str, ...]) -> Tuple[Tuple[str, ...], ...]:
    t = [(), (), (), ()]
    t[e - 1] = tuple(sorted(toks))
    return tuple(t)


def op_add(*ops: Op, coefs: Optional[Sequence[float]] = None) -> Op:
    out: Op = {}
    for k, o in enumerate(ops):
        c = 1.0 if coefs is None else coefs[k]
        for t, v in o.items():
            out[t] = out.get(t, 0.0) + c * v
    return {t: v for t, v in out.items() if v != 0.0}


def op_mul(A: Op, B: Op) -> Op:
    out: Op = {}
    for ta, va in A.items():
        for tb, vb in B.items():
            t = tuple(tuple(sorted(ta[e] + tb[e])) for e in range(4))
            out[t] = out.get(t, 0.0) + va * vb
    return {t: v for t, v in out.items() if v != 0.0}


def W1_op() -> Op:
    return {_single(e, ('w',)): 1.0 for e in ELEC}


def W2_op() -> Op:
    return op_mul_pairs([('w',), ('w',)])


def op_mul_pairs(toks2) -> Op:
    """sum_{i<j} tok_a(i) tok_b(j) + tok_b(i) tok_a(j) (symmetrised; = sum_{i<j} w_i w_j for w,w)."""
    out: Op = {}
    for (i, j) in combinations(ELEC, 2):
        for (a, b) in ((toks2[0], toks2[1]), (toks2[1], toks2[0])):
            t = [(), (), (), ()]
            t[i - 1] = tuple(sorted(a)); t[j - 1] = tuple(sorted(b))
            t = tuple(t)
            out[t] = out.get(t, 0.0) + 0.5
    return out


def lap_op(A: Op) -> Op:
    """sum_i lap_i A for A a product of w-weights: lap_i of a single w_i factor -> 'lapw';
    w_i^2 -> 2 w lapw + 2 |grad w|^2 is NOT needed by phase 1 (asserted)."""
    out: Op = {}
    for t, v in A.items():
        for e in range(4):
            toks = t[e]
            if not toks:
                continue
            assert toks == ('w',), f"lap of {toks} not supported"
            nt = list(t); nt[e] = ('lapw',)
            nt = tuple(nt)
            out[nt] = out.get(nt, 0.0) + v
    return out


def base_name(kind: str, toks: Tuple[str, ...]) -> str:
    return kind if not toks else kind + "*" + "*".join(toks)



TOKCIRC = {'w': 0, 'lapw': 0, 'u': 0, 'dzw': 0, 'z': 0, 'drwP': 1, 'drwM': -1, 'rhoP': 1, 'rhoM': -1}
GOp = Dict[tuple, float]


def _edge(i: int, j: int, key: str) -> Tuple[int, int, str]:
    return (min(i, j), max(i, j), key)


def g_from(op: Op) -> GOp:
    return {(t, ()): c for t, c in op.items()}


def g_geminal(kind: str) -> GOp:
    """'F' (plain), 'cc' (w_i w_j f_ij), 'c1' ((w_i + w_j) f_ij)."""
    out: GOp = {}
    for (i, j) in PAIRS6:
        e = (_edge(i, j, 'f'),)
        if kind == 'F':
            terms = [((), ())]
        elif kind == 'cc':
            terms = [(('w',), ('w',))]
        elif kind == 'c1':
            terms = [(('w',), ()), ((), ('w',))]
        else:
            raise ValueError(kind)
        for (ti, tj) in terms:
            t = [(), (), (), ()]; t[i - 1] = ti; t[j - 1] = tj
            k = (tuple(t), e)
            out[k] = out.get(k, 0.0) + 1.0
    return out


def g_mul(A: GOp, B: GOp) -> GOp:
    out: GOp = {}
    for (ta, ea), va in A.items():
        for (tb, eb), vb in B.items():
            t = tuple(tuple(sorted(ta[e] + tb[e])) for e in range(4))
            k = (t, tuple(sorted(ea + eb)))
            out[k] = out.get(k, 0.0) + va * vb
    return {k: v for k, v in out.items() if v != 0.0}


def g_lap(B: GOp, g: float) -> GOp:
    """sum_i lap_i B for B a sum of (w-weights) x (at most one f edge); exp geminal."""
    out: GOp = {}

    def add(t, e, c):
        k = (tuple(tuple(sorted(x)) for x in t), tuple(sorted(e)))
        out[k] = out.get(k, 0.0) + c

    for (t, edges), c in B.items():
        assert len(edges) <= 1 and all(k == 'f' for (_i, _j, k) in edges), edges
        for e in range(4):                                     # (lap a) f
            if t[e]:
                assert t[e] == ('w',), t[e]
                nt = list(t); nt[e] = ('lapw',)
                add(nt, edges, c)
        if not edges:
            continue
        (i, j, _k) = edges[0]
        add(t, ((i, j, 'f'),), 2.0 * g * g * c)                # a (lap_i + lap_j) f
        add(t, ((i, j, 'Y'),), -4.0 * g * c)
        for (k, l) in ((i, j), (j, i)):                         # 2 (da/dw_k) grad w_k . grad_k f
            if t[k - 1] != ('w',):
                continue
            base = list(t); base[k - 1] = ()
            for (tk, tl, cf) in (('u', None, -2.0 * g), ('dzw', 'z', 2.0 * g),
                                 ('drwP', 'rhoM', g), ('drwM', 'rhoP', g)):
                nt = list(base); nt[k - 1] = (tk,)
                if tl is not None:
                    nt[l - 1] = tuple(t[l - 1]) + (tl,)
                add(nt, ((i, j, 'Y'),), cf * c)
    return {k: v for k, v in out.items() if v != 0.0}


def _route(term) -> str:
    """'D' (dense diagram factors, m = 0 weights, f / Y edges) or 'G' (gvee_sparse on VariantRef2)."""
    toks, edges = term
    if any(k == 'coul' for (_i, _j, k) in edges):
        return 'G'
    if any(x in ('u', 'dzw', 'z', 'drwP', 'drwM', 'rhoP', 'rhoM') for tt in toks for x in tt):
        return 'G'
    return 'D'
