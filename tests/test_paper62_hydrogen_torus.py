"""Paper 62, Sec. hring, torus paragraph (record TR, docs/predictions/2026-11-02_hydrogen_torus.md,
GO): regression guard.  Written 2026-11-02 as a separate step after the paper (CLAUDE.md Sec. 9,
guard-writing rule).

Rebuilds small tori of M = 3 rings (N = 6, 7 fast, plus N = 8 slow; centre spacing 1.2 A, in-ring
spacing 1.0 A, the record's control geometry; one Slater 1s per H fitted by 3 Gaussians, Lowdin
site orbitals).  Nothing is imported from debug/; the basis helpers are the already-guarded ones in
tests/test_water_ccsd.py, and the Z2 counter and real DFT are written out here.

  * the toroidal winding k is exact: the k-breaking share of h and of (pq|rs) is below 1e-13;
  * Z-type Z2 counts site/tor/full are 2/4/5 at even N and 2/3/4 at odd N (hand-written from the
    record's derivation: two spin-number parities; plus the k-parity pieces, one more at even N;
    plus one from m);
  * slow (about 17 s): the poloidal winding m leaks with a |delta m| = 1 local slope against r/R
    that rises toward thinner tori and stays below 1.05, for h and for g.  The lower bound 0.50 is
    measured on these fatter tori (slopes h 0.575, 0.628; g 0.596, 0.672, measured 2026-11-02), not
    the record's TR3a bound 0.60, which was set on thinner ones (its controls start at N = 8);
  * slow (about 1 min): the record's first scored step itself, M = 3, N = 10 -> 11 at s = 1.4 A.
    The four leak shares match the stored scored rows to 1e-9, and both slopes lie in TR3a's
    [0.60, 1.05] (h 0.765, g 0.834).  Added 2026-11-02 (PI-approved) so that the record's own numbers
    reproduce, not only the mechanism.

Fire-tested 2026-11-02 with debug/qa/fire_test.py; every plant FIRED:
  1. an inner-outer asymmetry in the k direction (ring j = 0 pushed out by 1%);
  2. the tor labelling built as a DFT over the tube index a instead of the ring index j;
  3. the hand-written even-N tor count given as 3 (the odd-N value);
  4. the hand-written slope band centred on 2 (the |delta m| = 2 order) instead of 1;
  5. the (pq|rs) momentum change mislabelled (-1, 1, 1, -1 in place of -1, 1, -1, 1), slope test.
The scored-step test was fire-tested the same day, 4/4 FIRED: the control spacing 1.2 in place of 1.4;
the in-ring spacing 1.05 in place of 1.0; the |delta m| = 1 class counted one-sided (dm == 1); the
Slater exponent 1.30 in place of 1.24.
"""
import importlib.util
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62_tr', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)
GE = w.GE

ANG = 1.0 / 0.52917721067
M_TUBE = 3
S_RING, D_TUBE = 1.2, 1.0
NS = (6, 7)
ZERO = 1e-13
Z2_EXPECTED = {'even': {'site': 2, 'tor': 4, 'full': 5}, 'odd': {'site': 2, 'tor': 3, 'full': 4}}
SLOPE_BAND = (0.50, 1.05)   # lower bound measured here, see docstring
NS_SLOPE = (6, 7, 8)
S_SCORED = 1.4
TR3A_BAND = (0.60, 1.05)
# Scored rows of record TR (debug/data/torus_study.jsonl, sha256 c2de2d26...f3d1), M = 3, s = 1.4, d = 1.0:
# r/R and the Delta m = 1 shares of h and of (pq|rs).
SCORED = {10: (0.2548729212664928, 0.0679133026281965, 0.08782289376115512),
          11: (0.23236909647411566, 0.06327654672594189, 0.08130481997577564)}


def torus(N, M, s=S_RING):
    R = s * ANG / (2 * np.sin(np.pi / N))
    r = D_TUBE * ANG / (2 * np.sin(np.pi / M))
    nuc = []
    for j in range(N):
        ph = 2 * np.pi * j / N
        for a in range(M):
            th = 2 * np.pi * a / M
            rho = R + r * np.cos(th)
            nuc.append((np.array([rho * np.cos(ph), rho * np.sin(ph), r * np.sin(th)]), 1.0))
    return nuc, r / R


def site_integrals(nuc):
    fns, blocks = [], []
    for C, _ in nuc:
        w.sto_shell(C, 1, 0, 1.24, 3, fns, blocks)
    T = w.spherical_map(fns, blocks)
    S, h, g = w.ao_integrals(fns, nuc, T)
    return GE.transform_integrals(GE.lowdin_orbitals(S), h, g)


def real_dft(N):
    j = np.arange(N)
    cols = [np.full(N, 1 / np.sqrt(N))]
    for k in range(1, (N - 1) // 2 + 1):
        cols.append(np.sqrt(2.0 / N) * np.cos(2 * np.pi * k * j / N))
        cols.append(np.sqrt(2.0 / N) * np.sin(2 * np.pi * k * j / N))
    if N % 2 == 0:
        cols.append((-1.0) ** j / np.sqrt(N))
    return np.array(cols).T


def z2_count(h, g, tol=1e-10):
    """Number of independent Z-type Z2 symmetries of the JW qubit Hamiltonian: 2M minus the GF(2)
    rank of the X-parts of its terms (spin orbital 2p + s)."""
    M = h.shape[0]
    rows = set()
    for p, q in np.argwhere(np.abs(h) > tol).tolist():
        if p != q:
            for s in (0, 1):
                rows.add((1 << (2 * p + s)) ^ (1 << (2 * q + s)))
    for p, q, r, s in np.argwhere(np.abs(g) > tol).tolist():
        for a in (0, 1):
            for b in (0, 1):
                v = (1 << (2 * p + a)) ^ (1 << (2 * q + a)) ^ (1 << (2 * r + b)) ^ (1 << (2 * s + b))
                if v:
                    rows.add(v)
    pivots = {}
    for v in rows:
        while v:
            top = v.bit_length() - 1
            if top in pivots:
                v ^= pivots[top]
            else:
                pivots[top] = v
                break
    return 2 * M - len(pivots)


def km_shares(N, M, h, g):
    """k-breaking share, and the k-conserving |delta m| = 1 share, of h and of (pq|rs)."""
    U = np.kron(np.exp(2j * np.pi * np.outer(np.arange(N), np.arange(N)) / N) / np.sqrt(N),
                np.exp(2j * np.pi * np.outer(np.arange(M), np.arange(M)) / M) / np.sqrt(M))
    hc = U.conj().T @ h @ U
    gc = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, U.conj(), U, U.conj(), U, optimize=True)
    k, m = np.repeat(np.arange(N), M), np.tile(np.arange(M), N)
    out = {}
    for name, X, sg in (('h', hc, (-1, 1)), ('g', gc, (-1, 1, -1, 1))):
        idx = np.meshgrid(*([np.arange(N * M)] * X.ndim), indexing='ij')
        dk = sum(s * k[i] for s, i in zip(sg, idx)) % N
        dm = sum(s * m[i] for s, i in zip(sg, idx)) % M
        p2 = np.abs(X) ** 2
        tot = p2.sum()
        out[name + '_k'] = float(np.sqrt(p2[dk != 0].sum() / tot))
        out[name + '_m1'] = float(np.sqrt(p2[(dk == 0) & (np.minimum(dm, M - dm) == 1)].sum() / tot))
    return out


_CACHE = {}


def build(N):
    if N not in _CACHE:
        nuc, rR = torus(N, M_TUBE)
        h, g = site_integrals(nuc)
        labs = {'site': np.eye(N * M_TUBE), 'tor': np.kron(real_dft(N), np.eye(M_TUBE)),
                'full': np.kron(real_dft(N), real_dft(M_TUBE))}
        z2 = {}
        for name, V in labs.items():
            z2[name] = z2_count(*GE.transform_integrals(V, h, g))
        _CACHE[N] = {'rR': rR, 'shares': km_shares(N, M_TUBE, h, g), 'z2': z2}
    return _CACHE[N]


@pytest.fixture(scope='module')
def tori():
    return {N: build(N) for N in NS}


def test_toroidal_winding_is_exact(tori):
    """Wrong answer rejected: a torus whose small rings are not all alike (an inner-outer or
    radius asymmetry that varies with j), so k is no longer a symmetry."""
    for N, t in tori.items():
        for name in ('h', 'g'):
            assert t['shares'][name + '_k'] < ZERO, (N, name, t['shares'][name + '_k'])


def test_z2_counts_site_tor_full(tori):
    """Wrong answers rejected: a tor labelling that is not the k transform over the ring index; a
    Z2 derivation that forgets the extra real k = N/2 piece at even N; the m label credited with
    no extra parity."""
    for N, t in tori.items():
        assert t['z2'] == Z2_EXPECTED['even' if N % 2 == 0 else 'odd'], (N, t['z2'])


@pytest.mark.slow
def test_poloidal_leak_slope_rises_toward_one():
    """Wrong answers rejected: a poloidal leak that does not shrink with r/R; a leak of order
    (r/R)^2 for |delta m| = 1; slopes that fall toward thinner tori."""
    for name in ('h', 'g'):
        slopes = []
        for N1, N2 in zip(NS_SLOPE, NS_SLOPE[1:]):
            a, b = build(N1), build(N2)
            slopes.append(np.log(b['shares'][name + '_m1'] / a['shares'][name + '_m1'])
                          / np.log(b['rR'] / a['rR']))
        assert all(SLOPE_BAND[0] <= s <= SLOPE_BAND[1] for s in slopes), (name, slopes)
        assert all(s2 > s1 for s1, s2 in zip(slopes, slopes[1:])), (name, slopes)


@pytest.mark.slow
def test_scored_step_reproduces_the_record():
    """Wrong answers rejected: a guard that passes on its own geometry while the record's scored
    numbers no longer reproduce (a changed basis, spacing or share definition); a scored slope
    outside the record's own TR3a band."""
    got = {}
    for N, (rR, hm1, gm1) in SCORED.items():
        nuc, rR_now = torus(N, M_TUBE, S_SCORED)
        sh = km_shares(N, M_TUBE, *site_integrals(nuc))
        assert abs(rR_now - rR) < 1e-12, (N, rR_now)
        assert abs(sh['h_m1'] - hm1) < 1e-9 and abs(sh['g_m1'] - gm1) < 1e-9, (N, sh)
        assert sh['h_k'] < ZERO and sh['g_k'] < ZERO, (N, sh)
        got[N] = (rR_now, sh)
    for name in ('h', 'g'):
        s = (np.log(got[11][1][name + '_m1'] / got[10][1][name + '_m1'])
             / np.log(got[11][0] / got[10][0]))
        assert TR3A_BAND[0] <= s <= TR3A_BAND[1], (name, s)
