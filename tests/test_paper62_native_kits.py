"""Paper 62, Sec. native (Table of exponential kits): regression guard for the hydrogenic kit's water energy, the
abstract's "the hydrogenic Z/n rule by 196 mHa" (registry p62_hyd_fc, p62_hyd_above_ccpvdz_mha).  Written 2026-11-01
as a separate step (CLAUDE.md Sec. 9); before it, only the kit's one-electron spectrum was guarded
(tests/test_water_ccsd.py::test_hydrogenic_kit_is_the_hydrogen_spectrum).

Recomputed here: the bare-Z hydrogenic kit (decay Z/n, 24 functions, 12 Gaussians per function, 8 per cusp-free
polarization function, zero parameters), RHF and fc-CCSD(T).  cc-pVDZ is READ (registry p62_ccpvdz_fc; recomputed by
tests/test_paper62_matched_accuracy.py and tests/test_paper62_dzp_headline.py).
Not guarded here: the Z_eff = 6 variant, the single-scale Sturmian kit (its scale is a scan, "more than one hartree")
and the two-scale kit, which is ladder rung R2 (record docs/predictions/2026-10-30_scale_ladder.md).
Fire-tested 2026-11-01 (debug/qa/fire_test.py), 1/1 fired: the Z_eff = 6 valence kit in place of bare Z.
"""
import importlib.util
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62_kits', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)

HYD_FC = -76.045244                 # p62_hyd_fc
CCPVDZ_FC = -76.241025              # p62_ccpvdz_fc
GAP_MHA = 196


@pytest.mark.slow
def test_hydrogenic_kit_water_energy_and_its_loss_to_cc_pvdz():
    """Bare-Z hydrogenic kit on water: 24 functions, fc-CCSD(T) to 1e-6 of the registered value, so it sits 196 mHa
    above cc-pVDZ (47 min at 4 threads, measured 2026-11-01).

    Wrong answers rejected: the Z_eff = 6 valence variant (-75.947, 294 mHa above), a kit with the H 2p dropped
    (fewer than 24 functions), and an all-electron energy in place of frozen core.
    """
    nuc = w.water_nuclei()
    fns, T = w.hydrogenic_kit(nuc)
    assert T.shape[1] == 24
    S, h, g = w.ao_integrals(fns, nuc, T)
    r = w.run(S, h, g, 10, w.vnn(nuc), frozen_core=1)
    fc = float(r['CCSD(T)'])
    assert abs(fc - HYD_FC) < 1e-6, fc
    assert abs(round(1000 * (fc - CCPVDZ_FC)) - GAP_MHA) == 0
