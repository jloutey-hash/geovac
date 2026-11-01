"""Paper 62, Sec. dzp (Table tab:dzp): regression guard for the headline, the Slater DZP basis against cc-pVDZ on
water at equal size (24 functions).  Written 2026-11-01 as a separate step (CLAUDE.md Sec. 9, guard-writing rule),
closing the gap that every other Paper 62 guard READ the DZP energy as a registered literal.

Recomputed here, end to end, from the exponents stored beside this test (tests/data/paper62_water_dzp_2026_11_01.json,
a copy of debug/data/water_sto_exponents.json) and the cc-pVDZ basis file stored for the matched-accuracy guard:
  * fast: the free O atom's radial average-of-configuration energy at the stored primitives (the rule that set them);
  * slow: Slater DZP at 12 Gaussians per primitive and 8 per polarization function, RHF and fc-CCSD(T); cc-pVDZ the
    same; the equal-size gap 40.7 mHa and its mean-field share 28.8 mHa (registry p62_dzp_vs_ccpvdz_mha,
    p62_dzp_vs_ccpvdz_scf_mha);
  * slow: the convergence claim, 12 -> 16 Gaussians per primitive moves RHF and fc-CCSD(T) by about 0.07 mHa.
Fire-tested 2026-11-01 (debug/qa/fire_test.py), 3/3 fired: O s exponents scaled by 1.01 (fast test); the gap
taken as 39.5 mHa, the re-optimized-polarization value; 8 Gaussians on the primitives.  The 16 G test was not
fire-tested separately (it shares energies() with the headline test).
"""
import importlib.util
import json
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62_dzp', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)

EX = json.load(open(os.path.join(HERE, 'data', 'paper62_water_dzp_2026_11_01.json')))
DZ_BASIS = os.path.join(HERE, 'data', 'paper62_cc_pvdz_H_O_Ne_2026_10_31.json')
# Paper 62 Table tab:dzp
DZP_RHF, DZP_FC = -76.055573, -76.281739            # p62_dzp_rhf, p62_dzp_fc (12 G)
DZP16_RHF, DZP16_FC = -76.055644, -76.281811        # the 16 G row
CCPVDZ_RHF, CCPVDZ_FC = -76.026795, -76.241025      # p62_ccpvdz_fc
GAP_MHA, GAP_SCF_MHA = 40.7, 28.8
_CACHE = {}


def contracted_ex(d):
    sp = [(n, z) for n, z in zip(w.O_S_N, d['Os'])]
    _, Cs, cp = w.radial_atom(sp, [(2, z) for z in d['Op']])
    return dict(d, cs=Cs, cp=cp, O3d=d['O3d_fc'], H2p=d['H2p_fc'])


def energies(name, ng=12, ng_pol=8):
    key = (name, ng, ng_pol)
    if key not in _CACHE:
        nuc = w.water_nuclei()
        if name == 'dzp':
            fns, T = w.slater_dzp_contracted(nuc, contracted_ex(EX), ng, ng_pol=ng_pol)
        else:
            fns, T = w.gaussian_basis(DZ_BASIS, nuc)
        assert T.shape[1] == 24
        S, h, g = w.ao_integrals(fns, nuc, T)
        r = w.run(S, h, g, 10, w.vnn(nuc), frozen_core=1)
        _CACHE[key] = (float(r['RHF']), float(r['CCSD(T)']))
    return _CACHE[key]


def test_oxygen_primitives_reproduce_the_free_atom_energy():
    """The stored O primitives give the stored free-atom average energy (the rule that set them, nothing tuned on
    water), and the contraction is the atomic 1s, 2s, 2p.

    Wrong answers rejected: exponents copied with a slip (any one moved by 1%), and the p set used as s.
    """
    sp = [(n, z) for n, z in zip(w.O_S_N, EX['Os'])]
    E, Cs, cp = w.radial_atom(sp, [(2, z) for z in EX['Op']])
    assert abs(E - EX['E_O_atom_avg_radial']) < 1e-7, E
    assert np.shape(Cs) == (4, 2) and np.shape(cp) == (3,)


@pytest.mark.slow
def test_dzp_and_ccpvdz_recompute_and_the_headline_gap_follows():
    """Slater DZP (12 G, polarization 8 G) and cc-pVDZ on water, recomputed: RHF and fc-CCSD(T) to 1e-6 of Table
    tab:dzp, so the equal-size gap is 40.7 mHa and its mean-field share 28.8 mHa (about 25 min at 4 threads).

    Wrong answers rejected: the RHF-optimized polarization exponents in place of the fc ones, 8 Gaussians on the
    primitives (the HF control's shape), an all-electron correlation energy in place of frozen core, and the gap
    taken against cc-pVDZ with re-optimized polarization (39.5 mHa).
    """
    d_rhf, d_fc = energies('dzp')
    c_rhf, c_fc = energies('ccpvdz')
    assert abs(d_rhf - DZP_RHF) < 1e-6 and abs(d_fc - DZP_FC) < 1e-6, (d_rhf, d_fc)
    assert abs(c_rhf - CCPVDZ_RHF) < 1e-6 and abs(c_fc - CCPVDZ_FC) < 1e-6, (c_rhf, c_fc)
    assert abs(1000 * (c_fc - d_fc) - GAP_MHA) < 0.05
    assert abs(1000 * (c_rhf - d_rhf) - GAP_SCF_MHA) < 0.05


@pytest.mark.slow
def test_twelve_gaussians_are_converged_to_a_tenth_of_a_millihartree():
    """12 -> 16 Gaussians per primitive (polarization at 8) moves RHF and fc-CCSD(T) by 0.07 mHa, and the 16 G row
    reproduces Table tab:dzp to 1e-6 (about 35 min more at 4 threads).

    Wrong answers rejected: a fit that has not converged at 12 (a move above 0.1 mHa), and a 16 G row that does not
    match the table.
    """
    r12, f12 = energies('dzp')
    r16, f16 = energies('dzp', ng=16)
    assert abs(r16 - DZP16_RHF) < 1e-6 and abs(f16 - DZP16_FC) < 1e-6, (r16, f16)
    assert abs(r16 - r12) < 1e-4 and abs(f16 - f12) < 1e-4, (r16 - r12, f16 - f12)
    assert abs(f16 - f12) > 1e-5
