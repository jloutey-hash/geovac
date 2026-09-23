"""Guards for the Route C LiH ENERGY story (CHANGELOG v5.15.18-.20): the productionization
(sparse FCI + float64 ERI) and the additive-F12 near-chemical estimate (-8.062).

The geometry headline (R_eq +0.2%) and the C4 pi-ERI machinery are guarded by
test_prolate_allelectron_fci.py + the CHANGELOG chronicle; this file backs the ENERGY claims
Paper 19 rests on.

Fire-tested claims (what wrong answer each rejects):
  (1) sparse FCI == dense FCI bit-exact  -> rejects a wrong connected-pair filter / sign error
      in fci_fast (would shift the ground-state energy).
  (2) He-like cusp machinery reproduces KNOWN He exact  -> rejects a broken Hylleraas volume
      element / kinetic form (the machinery the additive estimate rests on).
  (3) Li core cusp ~30 mHa  -> rejects a wrong core-cusp number (the load-bearing 30.1 mHa that
      takes -8.032 -> -8.062).
  (4) [slow] float64 ERI + sparse FCI reproduces the banked mpf LiH energy to <50 uHa
      -> rejects a float64-assembly regression in the productionized engine.
"""
import os
import sys

import numpy as np
import pytest

_DEBUG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "debug")
if _DEBUG not in sys.path:
    sys.path.insert(0, _DEBUG)
sys.argv = [sys.argv[0]]   # debug drivers parse argv at import

import fci_fast                                              # noqa: E402
from prolate_allelectron_fci import fci_energy               # noqa: E402
from lih_r12_ceiling_probe import energy as helike_energy    # noqa: E402

RADIAL = ['t2', 's']
RICH = ['u', 't2', 's', 'u2', 'ut2']


def test_sparse_fci_bitexact_vs_dense():
    """(1) fci_fast.fci_energy_fast == dense fci_energy (M=8, 4e).  Rejects a connected-pair
    filter that drops real determinant couplings or a Slater-Condon sign error."""
    rng = np.random.default_rng(1)
    M = 8
    h1 = rng.standard_normal((M, M)); h1 = h1 + h1.T
    e = rng.standard_normal((M, M, M, M)) * 0.1
    e = e + e.transpose(1, 0, 3, 2)
    e = e + e.transpose(2, 3, 0, 1)                          # chemist (pq|rs)=(rs|pq)
    Ed, nd = fci_energy(h1, e, M, 4)
    Ef, ndf = fci_fast.fci_energy_fast(h1, e, M, 4)
    assert nd == ndf
    assert abs(Ed - Ef) < 1e-8, f"sparse {Ef} vs dense {Ed} (dE={Ed-Ef:.2e})"


def test_additive_f12_he_control():
    """(2) He (Z=2) rich Hylleraas basis reproduces KNOWN exact -2.90372 to <2 mHa.
    Validates the 2e cusp machinery the additive-F12 estimate uses.  Rejects a broken
    volume element / gradient-kinetic form (He would miss by >2 mHa)."""
    E_rich = helike_energy(2.0, 27 / 16, RICH)[1]
    assert abs(E_rich - (-2.90372)) < 2e-3, f"He rich {E_rich} vs exact -2.90372"
    # and it must be variational (above exact)
    assert E_rich > -2.90372 - 1e-6


def test_additive_f12_core_cusp():
    """(3) Li core (Z=3) r12 cusp = E(radial) - E(rich) is ~30 mHa (the load-bearing number
    that carries -8.032 -> -8.062).  Rejects a wrong core cusp outside [25, 35] mHa."""
    E_rad = helike_energy(3.0, 2.6875, RADIAL)[1]
    E_rich = helike_energy(3.0, 2.6875, RICH)[1]
    cusp = E_rad - E_rich
    assert 0.025 < cusp < 0.035, f"Li core cusp {cusp*1e3:.1f} mHa outside [25,35]"
    # rich must stay above the exact core (-7.27991): variational
    assert E_rich > -7.27991 - 1e-4


@pytest.mark.slow
def test_float64_engine_reproduces_mpf_lih():
    """(4) The productionized fast engine (float64 ERI + sparse FCI) reproduces the banked
    mpf LiH energy (-7.99468, minimal-valence M=6, R=3.015) to <50 uHa.  Rejects a float64
    X-table assembly regression.  ~60s."""
    import prolate_energy_ladder as L
    import prolate_float_eri as F
    import lih_core2exp_probe as P
    L.build_eri_tensor_m = F.build_eri_tensor_m_f
    L.fci_energy = fci_fast.fci_energy_fast
    Et, M, Mk, nd, cond, dt = P.run(Jb=1, Lb=0, npi=0, Jpi=0, Lpi=0, alpha=1.0, core2=None)
    assert Et > -8.070, f"variational violation: {Et} below exact"
    assert abs(Et - (-7.99468)) < 5e-5, f"float64 engine {Et} vs banked mpf -7.99468"
