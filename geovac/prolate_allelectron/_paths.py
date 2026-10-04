"""Single source of truth for the on-disk locations the prolate all-electron stack touches.

Decision (LiH productionization step 2a, 2026-10-04): the stack's data files stay in
``debug/data/`` and are resolved RELATIVE TO THE REPO ROOT through ``DEBUG_DATA`` -- nothing
is moved into, or ever written into, ``geovac/``.  What uses it:

* the ``__main__`` drivers of ``analytic_fci``, ``c4``, ``energy_ladder``,
  ``eri_core_grid_diagnostic``, ``pure_product_core_eri`` WRITE their logs there (as the debug/
  originals did);
* ``pure_product_core_eri.gate`` / ``run_cross`` READ
  ``prolate_separation_constant_stageC_ints.npz`` (a diagnostic cross-check only).

No test-backed or golden-backed computation of this package reads a data file.
"""
from __future__ import annotations

import os

REPO_ROOT: str = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEBUG_DATA: str = os.path.join(REPO_ROOT, "debug", "data")
