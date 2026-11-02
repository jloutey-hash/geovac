"""Paper 62, Sec. r3, the 3d set on the atom's correlated energy (record SL-3d,
docs/predictions/2026-11-01_scale_ladder_3d.md, GO): regression guard.  Written 2026-11-01
as a separate step after the paper (CLAUDE.md Sec. 9, guard-writing rule).

Inputs are a stored copy (tests/data/paper62_sl3d_2026_11_01.json, built by
debug/build_sl3d_guard_data_2026_11_01.py from debug/data/scale_ladder*.json(l)); nothing
is imported from debug/.

  * fast: the atom-set 3d values and the kits the water rows used; the four water/atom
    ratios, the band (re-derived from the two controls), SL3d-M, the RHF-part signs and the
    19.0 / 21.7 mHa positions of R5Hd, all derived from the stored energies by the record's
    own formulas, never copied.
  * slow (atom, 225 s measured at one thread): the free O atom at rung R3, UHF + fc-CCSD(T) in the
    Gaussian-fitted kit, recomputed at the tied 3d, the atom-set 3d and the atom-set 3d
    times 0.97 and 1.03: the stored energies reproduce and the atom-set 3d is the minimum.
  * slow (water, 854 s measured at 4 threads): rung R3d end to end, RHF and fc-CCSD(T) to 1e-6.

Fire-tested 2026-11-01 with debug/qa/fire_test.py; every plant FIRED:
  1. fast kits: the R3d row given the tied 3d (2.2174722116380785) in place of the atom-set one;
  2. fast ratios: R5Hd scored against the R5 water row instead of R5H;
  3. fast ratios: the atom UHF change used as A in place of the fc-CCSD(T) change;
  4. fast ratios: the band's lower bound mis-copied as 0.852;
  5. fast step: SL3d-M with the sign of the step flipped;
  6. fast step: SL3d-M referenced to the tied-3d R2 row instead of R2d;
  7. slow atom: the minimum judged on the UHF energy instead of fc-CCSD(T);
  8. slow atom: the tied 3d stored as the atom-set one;
  9. slow water: R3d built with the tied 3d.
"""
import importlib.util
import json
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62_sl3d', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)
_spec2 = importlib.util.spec_from_file_location('p62_ladder_for_sl3d', os.path.join(HERE, 'test_paper62_scale_ladder.py'))
lad = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(lad)        # the R3 probe's one-centre O atom (UHF + fc-CCSD(T)), already guarded

D = json.load(open(os.path.join(HERE, 'data', 'paper62_sl3d_2026_11_01.json')))
O_PRIMS = [tuple(p) for p in D['o_prims']]
H_PRIMS = [tuple(p) for p in D['h_prims']]
KEY = 'CCSD(T)_fc'
BAND = (0.844, 1.156)                 # record SL-3d, SL3d-T (3 x 0.052, the larger control |rho - 1|)
# scored case -> (tied-3d water row it is compared with, atom rung that set its 3d)
SCORED = {'R3d': ('R3', 'R3'), 'R4d': ('R4', 'R4'), 'R5d': ('R5', 'R5'), 'R5Hd': ('R5H', 'R5')}
# Paper 62 sec:r3 and the record's Disposition
PAPER_K3D = {'R2': 2.80, 'R3': 2.80, 'R4': 2.79, 'R5': 2.74}
PAPER_RATIO = {'R3d': 1.039, 'R4d': 1.052, 'R5d': 0.905, 'R5Hd': 0.848}


def expand(R):
    sp = D['registered_sp'][R]
    return [float(sp['ks'][g]) for g in sp['groups']]


def atom_change(R, key=KEY):
    a = D['atom'][R]
    return 1e3 * (a['new'][key] - a['old'][key])


def water_change(case, e=KEY):
    o, _ = SCORED[case]
    return 1e3 * (D['new_water'][case][e] - D['old_water'][o][e])


# ------------------------------------------------------------------------------- fast
def test_atom_set_3d_values_and_the_kits_each_water_row_used():
    """The atom puts the 3d at 2.80 / 2.80 / 2.79 / 2.74 (R2..R5), the old 3d was the 3p
    (tied), the atom-set 3d lowers the atom's fc-CCSD(T) and is the lowest of all stored R3
    evaluations, and every new water row uses the registered s/p exponents, its rung's
    atom-set 3d and the registered H set (R5Hd: R5's 3d, H2-set H).

    Wrong answers rejected: a water row built with the tied 3d in place of the atom-set one;
    R5Hd given a 3d other than R5's; s/p exponents re-optimized with the new 3d; an atom-set
    3d that is not the atom's correlated minimum among the points evaluated.
    """
    for R, k in PAPER_K3D.items():
        a = D['atom'][R]
        assert abs(a['k3d_new'] - k) < 0.005, (R, a['k3d_new'])
        assert a['s_p_exponents'] == expand(R)
        assert a['k3d_old'] == expand(R)[4]                       # tied to the 3p
        assert a['new'][KEY] < a['old'][KEY]
    evals = D['atom_evals_R3']
    assert min(e[KEY] for e in evals) > D['atom']['R3']['new'][KEY] - 1e-7
    for case, (o, R) in SCORED.items():
        r = D['new_water'][case]
        assert r['o5'] == expand(R), case
        assert r['k3d'] == D['atom'][R]['k3d_new'], case
        h = D['h_h2set'] if case == 'R5Hd' else D['h_k1']
        assert r['h_ks'] == h['h_ks'] and r['h_groups'] == h['h_groups'], case
    assert D['new_water']['R2d']['k3d'] == D['atom']['R2']['k3d_new']


def test_transfer_ratios_lie_in_the_band_the_controls_set():
    """SL3d-T: water change over atom fc-CCSD(T) change is 1.039, 1.052, 0.905, 0.848 on
    R3d, R4d, R5d, R5Hd, all inside 0.844-1.156; the half-width is three times the larger
    control |rho - 1| (C1 0.963, C2 1.052, recomputed from stored energies); R5Hd sits
    within 0.005 of the lower edge.  Also the reported parts: the water RHF part is negative
    on R3d, R4d and positive on R5d, R5Hd; the water correlation part is below the atom's.

    Wrong answers rejected: R5Hd scored against R5 instead of R5H (ratio 1.124, inside the
    band but not the paper's); the atom UHF change in place of the correlated one; a band
    bound mis-copied (0.852 puts R5Hd outside); a control ratio that does not reproduce.
    """
    c1 = D['control_c1']
    a1 = 1e3 * (c1['atom']['R3dB'][KEY] - c1['atom']['R3'][KEY])
    w1 = 1e3 * (c1['water_R3dB'][KEY] - D['old_water']['R3'][KEY])
    a2 = atom_change('R2')
    w2 = 1e3 * (D['new_water']['R2d'][KEY] - D['old_water']['R2'][KEY])
    rho1, rho2 = w1 / a1, w2 / a2
    assert abs(rho1 - 0.963) < 5e-4 and abs(rho2 - 1.052) < 5e-4, (rho1, rho2)
    half = 3 * max(abs(rho1 - 1), abs(rho2 - 1))
    assert abs(half - (BAND[1] - 1)) < 1e-3 and abs((1 - BAND[0]) - (BAND[1] - 1)) < 1e-12, half
    for case, (o, R) in SCORED.items():
        A = atom_change(R)
        rho = water_change(case) / A
        assert abs(rho - PAPER_RATIO[case]) < 5e-4, (case, rho)
        assert BAND[0] <= rho <= BAND[1], (case, rho)
        rhf = water_change(case, 'RHF')
        assert (rhf < 0) == (case in ('R3d', 'R4d')), (case, rhf)
        A_corr = A - atom_change(R, 'UHF')
        assert (water_change(case) - rhf) / A_corr < 1, case
    r5h = water_change('R5Hd') / atom_change('R5')
    assert 0 < r5h - BAND[0] < 0.005, r5h


def test_r3_stays_above_r2_and_r5hd_sits_between_ccpvdz_and_the_dzp():
    """SL3d-M: W(R3d) - W(R2d) = +0.781 mHa, so the R2 -> R3 step stays non-monotone; the
    best atom-set rung is R5Hd at -76.260011 Ha, 19.0 mHa below cc-pVDZ and 21.7 above the
    Slater DZP (fc-CCSD(T)).

    Wrong answers rejected: the step with its sign flipped; the step taken against the
    tied-3d R2 row (-7.78 mHa, monotone); a best rung other than R5Hd; the two distances
    swapped or taken against the all-electron energy.
    """
    nw = D['new_water']
    step = 1e3 * (nw['R3d'][KEY] - nw['R2d'][KEY])
    assert step > 0 and abs(step - 0.781) < 1e-3, step
    best = min(nw, key=lambda k: nw[k][KEY])
    assert best == 'R5Hd' and abs(nw[best][KEY] - (-76.260011)) < 5e-7
    below = 1e3 * (D['reference_fc']['ccpvdz'] - nw[best][KEY])
    above = 1e3 * (nw[best][KEY] - D['reference_fc']['slater_dzp'])
    assert abs(below - 19.0) < 0.05 and abs(above - 21.7) < 0.05, (below, above)


# ------------------------------------------------------------------------------- slow
@pytest.mark.slow
def test_r3_atom_set_3d_reproduces_and_is_the_correlated_minimum():
    """Rung R3, free O atom (M_S = 1, UHF + fc-CCSD(T), 12 G per s/p primitive, 8 per 3d),
    recomputed at the tied 3d, the atom-set 3d and that 3d times 0.97 and 1.03 (four
    evaluations at OMP_NUM_THREADS=1, 225 s in all measured 2026-11-01).  The stored UHF
    and fc-CCSD(T) energies reproduce to 1e-6, A(R3) = -10.240 mHa, and both neighbours lie
    above the atom-set point by more than 0.05 mHa at fc-CCSD(T).

    Wrong answers rejected: the minimum judged on the atom's UHF energy (whose minimum lies
    below 2.71, so the 0.97 point is lower); a tied 3d stored as the atom-set one; an atom
    change that is not the recomputed correlated one.
    """
    a = D['atom']['R3']
    k = a['k3d_new']
    sp = a['s_p_exponents']
    E = {}
    for name, kk in (('old', a['k3d_old']), ('new', k), ('lo', 0.97 * k), ('hi', 1.03 * k)):
        uhf, corr = lad._o_atom(tuple(sp) + (kk,))
        E[name] = {'UHF': uhf, KEY: uhf + corr}
    for name in ('old', 'new'):
        for e in ('UHF', KEY):
            assert abs(E[name][e] - a[name][e]) < 1e-6, (name, e, E[name][e], a[name][e])
    assert abs(1e3 * (E['new'][KEY] - E['old'][KEY]) - (-10.240)) < 1e-3
    crit = KEY
    for side in ('lo', 'hi'):
        assert E[side][crit] - E['new'][crit] > 5e-5, (side, 1e3 * (E[side][crit] - E['new'][crit]))


def _water_r3d():
    r = D['new_water']['R3d']
    nuc = w.water_nuclei()
    fns, cols_of = [], []
    for center, Z in nuc:
        if Z == 8.0:
            for (n, l), kk in zip(O_PRIMS, r['o5']):
                w._shell_kit(center, [(l, {n: 1.0}, kk)], 12, 8, fns, [], cols_of)
            w._shell_kit(center, [(2, {3: 1.0}, r['k3d'])], 12, 8, fns, [], cols_of)
        else:
            for (n, l), g in zip(H_PRIMS, r['h_groups']):
                w._shell_kit(center, [(l, {n: 1.0}, r['h_ks'][g])], 12, 8, fns, [], cols_of)
    T = w._kit_columns(fns, cols_of)
    assert T.shape[1] == 24
    S, h, g = w.ao_integrals(fns, nuc, T)
    res = w.run(S, h, g, 10, w.vnn(nuc), frozen_core=1)
    return float(res['RHF']), float(res['CCSD(T)'])


@pytest.mark.slow
def test_r3d_water_recomputes_end_to_end():
    """Rung R3d (the case SL3d-M is about), water at 24 functions, recomputed from the stored
    exponents with the builder and solver of tests/test_water_ccsd.py: RHF and fc-CCSD(T) to
    1e-6 of the stored row, and the ratio against the stored tied-3d R3 row and the stored
    atom change is 1.039 (854 s at OMP_NUM_THREADS=4 measured 2026-11-01, integrals dominate).

    Wrong answers rejected: a stored row computed with the tied 3d (the record's old kit);
    an all-electron energy stored as frozen core; a row whose exponents do not produce it.
    """
    rhf, fc = _water_r3d()
    r = D['new_water']['R3d']
    assert abs(rhf - r['RHF']) < 1e-6 and abs(fc - r[KEY]) < 1e-6, (rhf, fc)
    rho = 1e3 * (fc - D['old_water']['R3'][KEY]) / atom_change('R3')
    assert abs(rho - 1.039) < 5e-4, rho
