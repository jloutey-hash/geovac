"""Prolate all-electron diatomic stack (Paper 19's Route C LiH and its energy engine).

Promoted from debug/ in the LiH productionization, step 2a (2026-10-04); the promoted code is
bit-identical to the debug/ originals against ``tests/data/lih_golden/golden.json``.

Modules (import them explicitly -- this ``__init__`` imports nothing, because several modules
set the GLOBAL mpmath precision at import time and the import order is part of what the
golden values depend on):

==========================  =========================================================
``fci_grid``                grid MO integrals, generalized-m azimuthal kernels, dense
                            Slater-Condon determinant FCI (``fci_energy``)
``fci_fast``                sparse connected-pair FCI + Lanczos (``fci_energy_fast``)
``mixed_eri``               C2: analytic two-rate sigma core-valence ERIs, ``Orbital``
                            (import sets ``mp.mp.dps = 60``)
``analytic_fci``            C3: unified sigma one-body + sigma all-electron FCI (dps 60)
``c4``                      C4: general-m (pi) ERIs/one-body, ``OrbitalM``, the sigma+pi
                            R_eq scan (dps 60)
``float_eri``               float64 X-table assembly of the C4 ERIs (fast engine)
``energy_ladder``           Laguerre x Legendre re-based assembly ``assemble_rebased``
                            (dps 60)
``routec``                  the Route C LiH orbital sets (``lih_orbs``) and ``run``;
                            ``use_fast_engine`` (dps 60)
``helike_cusp``             He-like Hylleraas cusp probe (the 30.1 mHa Li core cusp)
``atomcentered_core``       C1 gate (import sets ``mp.dps = 50``!)
``core_hartree``            frozen-core Hartree screening one-body (``ZC_LI``)
``eri_core_grid_diagnostic``  C2 sizing diagnostic (grid ERIs vs exact)
``pure_product_core_eri``   Paper 19 Table tab:pure_product_core
``_paths``                  ``DEBUG_DATA`` (logs/data stay in debug/data)
==========================  =========================================================
"""
