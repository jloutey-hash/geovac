"""Paper 62, Sec. native (Table tab:native, row "Sturmian, one k per atom"): regression guard for the
abstract's "a single Coulomb-Sturmian scale per atom by more than one hartree".  Written 2026-11-02 as a
separate step (CLAUDE.md Sec. 9); until then the row was READ (docs/qa/paper_62.done.md, C8 row 4).

The kit: Coulomb Sturmians with one shared exponent per centre (tests/test_water_ccsd.py::sturmian_kit),
k_H = 1 and k_O = 3.7809502390032312, the scale that minimizes the radial oxygen energy (the value the
original run used, debug/water_qubit_resources.py, commit bef395e0; nothing is imported from debug/).

  * fast: the radial O atom in that one-scale space is -71.40617 Ha (paper: -71.406, 3.4 Ha above the
    Slater primitives' -74.768), and k_O is its minimum (1% either side is higher);
  * slow (water end to end, 24 functions): RHF -74.2879 and fc-CCSD(T) -74.560280 to 1e-6, so the kit is
    more than one hartree above cc-pVDZ (READ, registry p62_ccpvdz_fc; recomputed by
    tests/test_paper62_matched_accuracy.py) and above STO-3G's -75.0123.

Fire-tested 2026-11-02 with debug/qa/fire_test.py; every plant FIRED:
  1. fast: the scale set 1% off its minimum (3.743);
  2. fast: the p shells dropped from the shared scale's space (s only);
  3. slow: the two-scale kit's valence scale 2.281 used as the one shared O scale.
"""
import importlib.util
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62_ss', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)

K_O = 3.7809502390032312
ATOM_E = -71.40617                  # Paper 62 Sec. native: -71.406
SLATER_PRIM_ATOM = -74.768          # Paper 62 Sec. native, the Slater primitives
RHF, FC = -74.2879260204534, -74.56028035153616
CCPVDZ_FC = -76.241025              # p62_ccpvdz_fc
STO3G_FC = -75.0123                 # Table tab:native


def atom(k):
    return w.radial_atom([(1, k), (2, k), (3, k)], [(2, k), (3, k)])[0]


def test_one_scale_oxygen_atom_and_its_minimum():
    """Wrong answers rejected: a scale that is not the radial-atom minimum; a shared-scale space missing
    shells (fewer than s n <= 3, p n <= 3); an atom energy within 3 Ha of the Slater primitives."""
    e = atom(K_O)
    assert abs(e - ATOM_E) < 1e-5, e
    assert atom(0.99 * K_O) > e and atom(1.01 * K_O) > e
    assert round(e - SLATER_PRIM_ATOM, 1) == 3.4


@pytest.mark.slow
def test_one_scale_water_is_more_than_a_hartree_above_cc_pvdz():
    """Water, 24 functions, RHF and fc-CCSD(T) to 1e-6 of the stored run; more than 1 Ha above cc-pVDZ and
    above STO-3G.  Wrong answers rejected: the two-scale kit's valence scale as the shared scale; an
    all-electron energy in place of frozen core."""
    nuc = w.water_nuclei()
    fns, T = w.sturmian_kit(nuc, {8.0: K_O, 1.0: 1.0})
    assert T.shape[1] == 24
    S, h, g = w.ao_integrals(fns, nuc, T)
    r = w.run(S, h, g, 10, w.vnn(nuc), frozen_core=1)
    assert abs(float(r['RHF']) - RHF) < 1e-6, r['RHF']
    fc = float(r['CCSD(T)'])
    assert abs(fc - FC) < 1e-6, fc
    assert fc - CCPVDZ_FC > 1.0 and fc > STO3G_FC
