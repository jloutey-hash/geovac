"""Guard: FCI-atoms Table III (``tab:convergence_detail``) NNZ column.

Paper: papers/group2_quantum_chemistry/paper_fci_atoms.tex, Table III and the
sparsity paragraph ("for He at n_max = 5 it is 2.1% (759,739 nonzeros out of
5,995^2 entries)").  The NNZ is that of the exact-h1 (physical, diagonal
-Z^2/2n^2) slater_full Hamiltonian, entries with |H_IJ| >= 1e-8 kept.

Two routes are compared with the printed numbers:

  (1) the production assembly: ``geovac.lattice_index.LatticeIndex`` +
      ``geovac.direct_ci.DirectCISolver.assemble_hamiltonian``;
  (2) an INDEPENDENT count built here from scratch: angular factors from
      ``sympy.physics.wigner.gaunt`` (exact), radial R^k from a fine-grid
      cumulative-Simpson Hartree-potential integration written in this file,
      and an excitation-driven Slater-Condon count.  Nothing in route (2)
      calls the code under test.

Fire-test: a wrong selection rule in production, e.g. the historical m-transfer
sign bug in ``LatticeIndex._ck_coefficient`` (``q = ma - mc`` -> ``q = mc - ma``,
which silently drops every m-changing multipole), must make route (1) miss the
printed NNZ while route (2) still hits it.
"""
from __future__ import annotations

import contextlib
import io
import itertools
import math
import warnings
from functools import lru_cache
from typing import Dict, List, Tuple

import numpy as np
import pytest
from scipy.integrate import cumulative_simpson
from scipy.special import genlaguerre

THRESH = 1e-8

# (Z, n_electrons, n_max) -> (N_SD, NNZ) as printed in Table III.
TABLE_III: Dict[Tuple[int, int, int], Tuple[int, int]] = {
    (2, 2, 4): (1770, 82910),
    (2, 2, 5): (5995, 759739),
    (3, 3, 3): (3276, 111192),
}


# ----------------------------------------------------------------- route (2)
def _states(nmax: int) -> List[Tuple[int, int, int]]:
    return [(n, l, m) for n in range(1, nmax + 1)
            for l in range(n) for m in range(-l, l + 1)]


@lru_cache(maxsize=None)
def _gaunt(l1: int, l2: int, l3: int, m1: int, m2: int, m3: int) -> float:
    from sympy.physics.wigner import gaunt
    return float(gaunt(l1, l2, l3, m1, m2, m3))


def _radial_rk(nmax: int) -> Tuple[Dict[Tuple[int, int], int], np.ndarray]:
    """R^k[k, P(ac), P(bd)] at Z = 1 on a fine grid (scales as Z)."""
    nl = [(n, l) for n in range(1, nmax + 1) for l in range(n)]
    idx = {p: i for i, p in enumerate(nl)}
    npts = 100001                       # odd: composite Simpson weights
    r = np.linspace(0.0, 250.0, npts)
    hstep = r[1] - r[0]
    wts = np.full(npts, 2.0)
    wts[1::2] = 4.0
    wts[0] = wts[-1] = 1.0
    wts *= hstep / 3.0
    R = np.empty((len(nl), npts))
    for i, (n, l) in enumerate(nl):
        x = 2.0 * r / n
        R[i] = np.exp(-x / 2) * x ** l * genlaguerre(n - l - 1, 2 * l + 1)(x)
        R[i] /= math.sqrt(np.dot(wts, R[i] ** 2 * r ** 2))
    npair = len(nl) ** 2
    rho = np.einsum("ar,cr->acr", R, R).reshape(npair, npts)
    kmax = 2 * (nmax - 1)
    out = np.zeros((kmax + 1, npair, npair))
    rs = r[1:]
    w = rho * r ** 2 * wts               # outer-integral measure, incl. r = 0
    for k in range(kmax + 1):
        # integrands vanish at r = 0 for every Gaunt-allowed (pair, k)
        f_in = np.zeros_like(rho)
        f_out = np.zeros_like(rho)
        f_in[:, 1:] = rho[:, 1:] * rs ** (k + 2)
        f_out[:, 1:] = rho[:, 1:] * rs ** (1 - k)
        inner = cumulative_simpson(f_in, x=r, axis=1, initial=0.0)
        outer = (np.dot(f_out, wts)[:, None]
                 - cumulative_simpson(f_out, x=r, axis=1, initial=0.0))
        Y = np.zeros_like(rho)
        Y[:, 1:] = inner[:, 1:] * rs ** (-k - 1) + outer[:, 1:] * rs ** k
        out[k] = Y @ w.T                  # out[k][ac, bd]
        del f_in, f_out, inner, outer, Y
    return idx, out


def _spatial_eri(Z: int, nmax: int) -> np.ndarray:
    """eri[a,b,c,d] = <ab|1/r12|cd>, electron 1 in (a,c), electron 2 in (b,d)."""
    st = _states(nmax)
    n = len(st)
    idx, rk = _radial_rk(nmax)
    nl_n = int(round(math.sqrt(rk.shape[1])))
    P = np.array([idx[(s[0], s[1])] for s in st])
    m = np.array([s[2] for s in st])
    # mask[a,b,c,d] = (m_c - m_a) == (m_b - m_d)
    mask = (m[None, None, :, None] - m[:, None, None, None]) == \
           (m[None, :, None, None] - m[None, None, None, :])
    eri = np.zeros((n, n, n, n))
    for k in range(rk.shape[0]):
        F = np.zeros((n, n))   # F[x,y] = int Y_x^* Y_{k,mx-my} Y_y
        for x, (_, lx, mx) in enumerate(st):
            for y, (_, ly, my) in enumerate(st):
                q = mx - my
                if abs(q) <= k and (lx + k + ly) % 2 == 0:
                    F[x, y] = (-1) ** mx * _gaunt(lx, k, ly, -mx, q, my)
        if not F.any():
            continue
        Rk = rk[k].reshape(nl_n, nl_n, nl_n, nl_n)
        Rfull = Rk[P[:, None, None, None], P[None, None, :, None],
                   P[None, :, None, None], P[None, None, None, :]]
        # A1 = int Y_a^* Y_kq^* Y_c = F[c,a]; A2 = F[b,d]
        eri += (4 * math.pi / (2 * k + 1)) * np.einsum(
            "ca,bd->abcd", F, F) * mask * Rfull
    return Z * eri


def independent_nnz(Z: int, nel: int, nmax: int) -> Tuple[int, int]:
    st = _states(nmax)
    nsp = len(st)
    nso = 2 * nsp
    eri = _spatial_eri(Z, nmax)
    h = np.array([-Z * Z / (2.0 * s[0] ** 2) for s in st])
    sp = np.arange(nso) // 2
    sg = np.arange(nso) % 2

    def g_slice(p: int, q: int) -> np.ndarray:
        """<pq||rs> over all spin-orbitals r, s."""
        d = (eri[sp[p], sp[q]][sp[:, None], sp[None, :]]
             * ((sg[:, None] == sg[p]) & (sg[None, :] == sg[q])))
        x = (eri[sp[p], sp[q]][sp[None, :], sp[:, None]]
             * ((sg[None, :] == sg[p]) & (sg[:, None] == sg[q])))
        return d - x

    dets = list(itertools.combinations(range(nso), nel))
    nnz = 0
    iu = np.triu_indices(nso, 1)
    for D in dets:
        occ = np.zeros(nso, bool)
        occ[list(D)] = True
        virt = ~occ
        diag = h[sp[list(D)]].sum()
        for p, q in itertools.combinations(D, 2):
            diag += g_slice(p, q)[p, q]
        nnz += abs(diag) >= THRESH
        G = {(p, q): g_slice(p, q) for p, q in itertools.combinations(D, 2)}
        for p in D:   # singles p -> r: sum_q <pq||rq>
            me = np.zeros(nso)
            for q in D:
                if q == p:
                    continue
                gpq = G[(p, q)] if p < q else -G[(q, p)]
                me += gpq[:, q]
            nnz += int(np.count_nonzero((np.abs(me) >= THRESH) & virt))
        vv = virt[iu[0]] & virt[iu[1]]
        for (p, q), g in G.items():   # doubles pq -> rs
            nnz += int(np.count_nonzero((np.abs(g[iu]) >= THRESH) & vv))
    return len(dets), nnz


# ----------------------------------------------------------------- route (1)
def production_nnz(Z: int, nel: int, nmax: int) -> Tuple[int, int]:
    from geovac.lattice_index import LatticeIndex
    from geovac.direct_ci import DirectCISolver
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        li = LatticeIndex(n_electrons=nel, max_n=nmax, nuclear_charge=Z,
                          vee_method="slater_full", h1_method="exact")
        H = DirectCISolver(li).assemble_hamiltonian()
    return li.n_sd, int(H.nnz)


FAST = [(2, 2, 4), (3, 3, 3)]


@pytest.mark.parametrize("key", FAST, ids=lambda k: f"Z{k[0]}_n{k[2]}")
def test_production_nnz_matches_table_iii(key):
    assert production_nnz(*key) == TABLE_III[key]


@pytest.mark.parametrize("key", FAST, ids=lambda k: f"Z{k[0]}_n{k[2]}")
def test_independent_nnz_matches_table_iii(key):
    assert independent_nnz(*key) == TABLE_III[key]


@pytest.mark.slow
def test_he_nmax5_nnz_both_routes():
    key = (2, 2, 5)
    assert production_nnz(*key) == TABLE_III[key]
    assert independent_nnz(*key) == TABLE_III[key]
