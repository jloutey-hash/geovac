r"""Route C LiH driver: the core-enriched orbital sets and the energy call.

Promoted from ``debug/lih_core2exp_probe.py`` (LiH productionization step 2a, 2026-10-04):
ONLY ``R``, ``E_EXACT``, ``lih_orbs`` and ``run`` -- the functions tests, the golden oracle and
the downstream drivers (lih_marriage_phase0, lih_pi_m16, lih_pi_m16_primcache, lih_vmc,
float_eri.validate/push) call.  Code verbatim; numerics unchanged (bit-identical against
tests/data/lih_golden/golden.json).  The probe's scan/pairs/bank/confirm experiment drivers
stay in debug/ and import from here.

``lih_orbs`` builds the Route C LiH set at R = 3.015: the analytic Li 1s core (Slater
2.6875) + optional extra core STO exponents (``core2``), Li 2s / H 1s / H^- anchors, a
bond-centred sigma ladder xi^j eta^l e^{-alpha xi} (j<=Jbond, l<=Lbond) and optional pi
shells.  ``run`` assembles it through ``energy_ladder.assemble_rebased`` (Laguerre x
Legendre re-based FCI) at Z_A=3, Z_B=1, 4 electrons, V_nn = 3/R.

Paper 19 results this backs (with the fast engine, ``use_fast_engine``):
  run(2,1,0,0,0,1.0,core2=(4.5,1.6))  -> -8.02335  (sigma core-enriched reference, M=12)
  run(2,1,1,1,0,1.0,core2=(4.5,1.6))  -> -8.02905  (M=16 orbital ceiling; pi reference)
  run(3,1,1,2,0,1.0,core2=(4.5,1.6))  -> -8.0317   (M=20 pure-orbital push, "-8.032")
  run(1,0,0,0,0,1.0)                   -> -7.99468  (float64-vs-mpf guard, M=6)

IMPORT SIDE EFFECT: sets the global ``mp.mp.dps = 60`` before importing ``energy_ladder``
(exactly as the probe did).
"""
from __future__ import annotations

import time
from typing import List, Optional, Sequence, Tuple

import mpmath as mp
mp.mp.dps = 60

from geovac.prolate_allelectron import energy_ladder as L            # noqa: E402
from geovac.prolate_allelectron.c4 import (                          # noqa: E402
    OrbitalM, from_sigma, valence_prolate_orbital, valence_pi_orbital, ZC_LI)
from geovac.prolate_allelectron.mixed_eri import sto_orbital         # noqa: E402
from geovac.prolate_allelectron.analytic_fci import sto_orbital_B    # noqa: E402

R = 3.015
E_EXACT = -8.070


def lih_orbs(Jbond: int, Lbond: int, npi: int, Jpi: int, Lpi: int, alpha: float,
             core2: Optional[Sequence[float]] = None) -> Tuple[List[OrbitalM], List[Tuple]]:
    """Route C LiH ladder + optional 2nd-core STO exponents (list of zetas)."""
    orbs, tags = [], []

    def add(o, key):
        orbs.append(o); tags.append((key, 0, 0))

    add(from_sigma(sto_orbital(ZC_LI, R, is_core=True)), ('coreLi', 0, 0))
    for k, z2 in enumerate(core2 or []):
        add(from_sigma(sto_orbital(mp.mpf(z2), R, is_core=True)), (f'coreLi{k+2}', 0, 0))
    add(from_sigma(sto_orbital(mp.mpf('0.65'), R)), ('Li2s', 0, 0))
    add(from_sigma(sto_orbital_B(mp.mpf('1.0'), R)), ('H1s', 0, 0))
    add(from_sigma(sto_orbital_B(mp.mpf('0.70'), R)), ('Hm', 0, 0))
    for j in range(Jbond + 1):
        for l in range(Lbond + 1):
            orbs.append(from_sigma(valence_prolate_orbital(j, l, mp.mpf(alpha))))
            tags.append((('bond', 0, Jbond, Lbond), j, l))
    if npi > 0:
        for ms in (+1, -1):
            for j in range(Jpi + 1):
                for l in range(Lpi + 1):
                    orbs.append(valence_pi_orbital(j, l, mp.mpf(alpha), ms))
                    tags.append((('bondpi', ms, Jpi, Lpi), j, l))
    return orbs, tags


def run(Jb: int, Lb: int, npi: int, Jpi: int, Lpi: int, alpha: float,
        core2: Optional[Sequence[float]] = None, tol: float = 1e-11
        ) -> Tuple[float, int, int, int, float, float]:
    """Return (E_tot, M, M_kept, n_det, cond(S), wall seconds) for one Route C LiH basis."""
    orbs, tags = lih_orbs(Jb, Lb, npi, Jpi, Lpi, alpha, core2)
    t = time.time()
    Et, Ee, Mk, nd, cond = L.assemble_rebased(orbs, tags, R, 3.0, 1.0, 4, 3.0 / R,
                                              alpha, tol=tol)
    return Et, len(orbs), Mk, nd, cond, time.time() - t


def use_fast_engine() -> None:
    """Select the productionized engine for ``run``: float64 ERI assembly + sparse FCI.

    Rebinds the module globals ``energy_ladder.build_eri_tensor_m`` and
    ``energy_ladder.fci_energy`` -- exactly the patch the slow Route C test, the golden
    oracle and ``float_eri.validate/push`` apply (it is process-global and not undone).
    """
    from geovac.prolate_allelectron import float_eri as F
    from geovac.prolate_allelectron import fci_fast
    L.build_eri_tensor_m = F.build_eri_tensor_m_f
    L.fci_energy = fci_fast.fci_energy_fast
