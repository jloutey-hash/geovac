"""Second-route check of debug/he_spectrum_phys_remeasure.py energies via the
production builder build_graph_native_fci(..., h1_method='exact', subblock=...).
Appends to debug/data/residual_sweep.log."""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from geovac.casimir_ci import build_graph_native_fci  # noqa: E402

HA_CM = 219474.6313632
lines = ["--- cross-check (build_graph_native_fci, h1_method='exact') ---"]
for n in (3, 5):
    s = np.sort(np.linalg.eigvalsh(build_graph_native_fci(
        2, n, spin="singlet", subblock=(0, 0), h1_method="exact")))
    t = np.sort(np.linalg.eigvalsh(build_graph_native_fci(
        2, n, spin="triplet", subblock=(0, 0), h1_method="exact")))
    p = np.sort(np.linalg.eigvalsh(build_graph_native_fci(
        2, n, spin="singlet", subblock=(0, 1), h1_method="exact")))
    lines.append(f"  n={n}: E11S={s[0]:.6f} E21S={s[1]:.6f} E23S={t[0]:.6f} "
                 f"dE={(s[1]-t[0])*HA_CM:.1f} cm-1  E(2^1P, sp M_L=0)={p[0]:.6f}")
with open(ROOT / "debug" / "data" / "residual_sweep.log", "a", encoding="utf-8") as fp:
    fp.write("\n".join(lines) + "\n")
print("\n".join(lines))
