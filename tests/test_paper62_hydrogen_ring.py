"""Paper 62, Sec. hring, the H6 ring (record HR, docs/predictions/2026-11-01_hydrogen_ring.md, STOP):
regression guard.  Written 2026-11-02 as a separate step after the paper (CLAUDE.md Sec. 9,
guard-writing rule).

Inputs are a stored copy (tests/data/paper62_hring_2026_11_02.json, built by
debug/build_hring_guard_data_2026_11_02.py from debug/data/hring_study.jsonl and
hring_exponents.json, sha256 inside); nothing is imported from debug/.

  * fast: Table tab:hring's "vs DZP" column, the 4.9 mHa mean-field part and the 1.5 mHa given
    back, the STOP against the registered 6.2 mHa bar, the matched-accuracy turnover (both larger
    Gaussian bases below the Slater DZP, cc-pVTZ below), the 6.3x / 3.6x cost of the cheapest
    winner, the up-to-1.4% spread of canonical Pauli counts, and the label results (Z2 5 vs 3;
    4.37x at 1e-10; site fewer at 1e-2 and at 1e-3, ring fewer at 1e-4; identical fixed-label
    counts in the Slater and cc-pVDZ spaces), all derived from the stored rows, never copied.
  * slow (about 1 min): the re-optimized cc-pVDZ row recomputed end to end (RHF, CCSD, CCSD(T)
    to 1e-6) with the original spin-orbital solver of tests/test_water_ccsd.py, an independent
    route from the memory-lean copy that produced the row.
  * slow (about 1 min): the published cc-pVDZ space relabelled from scratch (site and ring
    labels): Z2 counts 3 and 5, Pauli counts at 1e-10 and 1e-2 reproduce exactly.
  The Slater DZP row (45 min) and the larger bases are READ.

Fire-tested 2026-11-02 with debug/qa/fire_test.py; every plant FIRED:
  1. fast: the Slater DZP CCSD(T) taken as the CCSD (gap measured on the wrong energy);
  2. fast: the gap scored against the published cc-pVDZ instead of the re-optimized one;
  3. fast: the cheapest winner taken as aug-cc-pVDZ instead of cc-pVDZ+TZpol;
  4. fast: the ring and site label rows swapped;
  5. slow energy: the published p exponent 0.727 in place of the re-optimized one;
  6. slow labels: the ring transform taken over the orbital type instead of the site index;
  7. slow labels: the local p frame left Cartesian (no rotation to radial / tangential).
"""
import importlib.util
import json
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62_hr', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)
_spec2 = importlib.util.spec_from_file_location('p62_torus_for_hr', os.path.join(HERE, 'test_paper62_hydrogen_torus.py'))
tor = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(tor)          # real_dft and the Z2 counter, already guarded
GE = w.GE

D = json.load(open(os.path.join(HERE, 'data', 'paper62_hring_2026_11_02.json')))
F = D['final']
L = D['labels']
DZ_BASIS = os.path.join(HERE, 'data', 'paper62_cc_pvdz_H_O_Ne_2026_10_31.json')
ANG = 1.0 / 0.52917721067
CC = 'CCSD(T)'
BAR = 6.2                                  # record HR, HR-A1 (mHa)
DZP, DZ_OPT = 'sto_dzp', 'ccpvdz_polopt'


def mha(a, b, key=CC):
    return 1e3 * (F[a][key] - F[b][key])


# ------------------------------------------------------------------------------- fast
def test_equal_size_gap_misses_the_bar():
    """Table vs DZP: cc-pVDZ +4.8, cc-pVDZ* +3.5; the Slater DZP wins but by less than 6.2 mHa.

    Wrong answers rejected: the gap measured on the CCSD or RHF energy; the gap scored against the
    published cc-pVDZ; a Slater DZP above the Gaussian (wrong sign)."""
    assert round(mha('ccpvdz', DZP), 1) == 4.8
    gap = mha(DZ_OPT, DZP)
    assert round(gap, 1) == 3.5 and 0 < gap < BAR


def test_mean_field_part_exceeds_the_gap():
    """The RHF gap is 4.9 mHa, larger than the whole gap; correlation gives 1.5 mHa back.

    Wrong answer rejected: the correlation part taken with the wrong sign (adding to the gap)."""
    rhf = mha(DZ_OPT, DZP, 'RHF')
    assert round(rhf, 1) == 4.9
    assert round(rhf - mha(DZ_OPT, DZP), 1) == 1.5


def test_matched_accuracy_turns_over():
    """aug-cc-pVDZ -1.1, cc-pVDZ+TZpol -3.6, cc-pVTZ -17.5 mHa; the cheapest winner (by Pauli
    count) is cc-pVDZ+TZpol at 6.3x the Pauli terms and 3.6x the 1-norm of the Slater DZP.

    Wrong answers rejected: the cheapest winner taken as aug-cc-pVDZ (more functions, more terms);
    the cost ratios taken against the re-optimized cc-pVDZ."""
    assert [round(mha(b, DZP), 1) for b in ('augccpvdz', 'ccpvdz_tzpol', 'ccpvtz')] == [-1.1, -3.6, -17.5]
    winners = [b for b in F if mha(b, DZP) < 0]
    cheapest = min(winners, key=lambda b: F[b]['qubit']['n_pauli'])
    assert cheapest == 'ccpvdz_tzpol' and F[cheapest]['nbf'] == 48
    assert round(F[cheapest]['qubit']['n_pauli'] / F[DZP]['qubit']['n_pauli'], 1) == 6.3
    assert round(F[cheapest]['qubit']['one_norm'] / F[DZP]['qubit']['one_norm'], 1) == 3.6


def test_canonical_counts_spread_fixed_labels_agree():
    """Canonical counts of the three 30-function rows differ by up to 1.4%; in fixed labellings the
    Slater and cc-pVDZ counts at 1e-10 are identical (149,591 ring, 653,971 site).

    Wrong answer rejected: a radial-shape effect (fixed-label counts differing between spaces)."""
    c = [F[b]['qubit']['n_pauli'] for b in ('ccpvdz', DZ_OPT, DZP)]
    assert len(set(c)) == 3 and round(100 * (max(c) - min(c)) / min(c), 1) == 1.4
    for lab, n in (('ring', 149591), ('site', 653971)):
        assert L[DZP][lab]['pauli_1e-10'] == L['ccpvdz'][lab]['pauli_1e-10'] == n


def test_ring_labels_exact_site_labels_local():
    """Both spaces: Z2 5 (ring) vs 3 (site); site/ring at 1e-10 = 4.37; 87% of ring integrals
    exactly zero; site fewer at 1e-2 (Slater 4,267 vs 7,987) and 1e-3, ring fewer at 1e-4.

    Wrong answer rejected: the ring and site rows swapped (every inequality reverses)."""
    for sp in (DZP, 'ccpvdz'):
        s, r = L[sp]['site'], L[sp]['ring']
        assert (s['z2'], r['z2']) == (3, 5)
        assert round(s['pauli_1e-10'] / r['pauli_1e-10'], 2) == 4.37
        assert round(r['exact_zero_eri_frac'], 2) == 0.87
        assert s['pauli_0.01'] < r['pauli_0.01'] and s['pauli_0.001'] < r['pauli_0.001']
        assert s['pauli_0.0001'] > r['pauli_0.0001']
    assert (L[DZP]['site']['pauli_0.01'], L[DZP]['ring']['pauli_0.01']) == (4267, 7987)


# ------------------------------------------------------------------------------- slow
def ring_nuclei():
    Rc = D['geometry']['spacing_ang'] * ANG / (2 * np.sin(np.pi / 6))
    return [(np.array([Rc * np.cos(2 * np.pi * k / 6), Rc * np.sin(2 * np.pi * k / 6), 0.0]), 1.0)
            for k in range(6)]


def h_basis(p=None):
    data = json.load(open(DZ_BASIS))
    shells = json.loads(json.dumps(data['elements']['1']['electron_shells']))
    if p is not None:
        for sh in shells:
            if sh['angular_momentum'][0] == 1:
                sh['exponents'] = [repr(p)]
    path = os.path.join(HERE, 'data', '_p62_hring_tmp_%d.json' % os.getpid())
    json.dump({'elements': {'1': {'electron_shells': shells}}}, open(path, 'w'))
    try:
        return w.gaussian_basis(path, ring_nuclei())
    finally:
        os.remove(path)


@pytest.mark.slow
def test_reoptimized_ccpvdz_row_recomputes():
    """RHF, CCSD and CCSD(T) of cc-pVDZ* reproduce the stored row to 1e-6 Ha."""
    nuc = ring_nuclei()
    fns, T = h_basis(D['p_exponents']['gto_p'])
    S, h, g = w.ao_integrals(fns, nuc, T)
    Ehf, eps, C = w.rhf(S, h, g, 6, w.vnn(nuc))
    fs, anti = w.spin_orbital(h, g, C, eps)
    Ec, Et = w.ccsd(fs, anti, 6)
    got = {'RHF': Ehf, 'CCSD': Ehf + Ec, CC: Ehf + Ec + Et}
    for k, v in got.items():
        assert abs(v - F[DZ_OPT][k]) < 1e-6, (k, v, F[DZ_OPT][k])


@pytest.mark.slow
def test_published_ccpvdz_labels_recompute():
    """Site and ring labels of the published cc-pVDZ space, built here: Z2 3 and 5; Pauli counts
    at 1e-10 and 1e-2 equal the stored ones exactly."""
    nuc = ring_nuclei()
    fns, T = h_basis()
    assert np.allclose(T, np.eye(T.shape[0]))             # s and p only: Cartesian = spherical
    per = T.shape[0] // 6
    assert per == 5
    for j in range(6):                                    # s, s', px, py, pz on each site, in order
        assert all(np.allclose(fns[j * per + t].center, nuc[j][0]) for t in range(per))
        assert [fns[j * per + t].lmn for t in range(2, 5)] == [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    for j in range(6):                                    # p triple -> (radial, tangential, z)
        th = 2 * np.pi * j / 6
        c, s = np.cos(th), np.sin(th)
        T[j * per + 2: j * per + 4, j * per + 2: j * per + 4] = [[c, -s], [s, c]]
    S, h, g = w.ao_integrals(fns, nuc, T)
    hs, gs = GE.transform_integrals(GE.lowdin_orbitals(S), h, g)
    U = np.zeros((6 * per, 6 * per))
    Fr = tor.real_dft(6)
    for t in range(per):
        U[t::per, t * 6:(t + 1) * 6] = Fr
    hr, gr = GE.transform_integrals(U, hs, gs)
    enuc = w.vnn(nuc)
    for name, (hh, gg) in (('site', (hs, gs)), ('ring', (hr, gr))):
        ref = L['ccpvdz'][name]
        assert tor.z2_count(hh, gg) == ref['z2'], name
        for tol in (1e-10, 1e-2):
            n = int(w.jw_pauli_stats(hh, gg, enuc, tol=tol)['n_pauli'])
            assert n == ref['pauli_%g' % tol], (name, tol, n, ref['pauli_%g' % tol])
