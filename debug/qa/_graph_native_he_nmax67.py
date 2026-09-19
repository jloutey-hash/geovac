r"""Settle the validation_benchmarks.md monotone inconsistency: compute the
graph-native He CI error at n_max = 6 and 7, the unmeasured rungs.

SELF-VALIDATING, because the danger here is stamping a NEW wrong number by
using the wrong pipeline.  There are several He pipelines (self-consistent
variational-k `convergence_table`; `solve_variational_fast`; the fixed-k=Z
graph-native `build_graph_native_fci`).  Only the LAST is the 0.19% headline
("zero-parameter, exact algebraic integrals, k=Z").  So:

  1. compute n_max=5 first and require it to reproduce the KNOWN ladder value
     0.24963% (the parent /qa run's reproduced n_max=5) to within 0.01% abs.
     If it does not, the pipeline or convention is wrong -> ABORT, report
     nothing for 6/7.  A number produced by an unvalidated path is worse than
     the honest "not measured" flag already in validation_benchmarks.md.
  2. only on a passing anchor, compute n_max=6 and n_max=7.

The monotone question it settles: validation_benchmarks.md records n_max=7
"< 0.20%", n_max=8 "0.207%", n_max=9 "0.201%" -- impossible on a
monotone-decreasing variational ladder.  A measured n_max=7 says which row is
wrong (or whether the whole ladder needs re-measuring).

E_exact = -2.903724377034 Ha (Pekeris/Drake), the value casimir_ci and the
certified artifact both use.

Run (background):  python debug/qa/_graph_native_he_nmax67.py
Writes: debug/data/graph_native_he_nmax67.txt
"""
from __future__ import annotations

import time
import numpy as np

from geovac.casimir_ci import build_graph_native_fci

E_EXACT = -2.903724377034
ANCHOR_NMAX = 5
ANCHOR_PCT = 0.24963      # parent /qa run's reproduced value
ANCHOR_TOL = 0.01         # abs % tolerance on the anchor
OUT = "debug/data/graph_native_he_nmax67.txt"


def err_pct(n_max: int) -> tuple[float, float, int]:
    t0 = time.time()
    H = build_graph_native_fci(Z=2, n_max=n_max, m_total=0, spin="singlet")
    dim = H.shape[0]
    E0 = float(np.linalg.eigvalsh(H)[0])
    pct = abs((E0 - E_EXACT) / E_EXACT) * 100.0
    return E0, pct, dim, time.time() - t0


def main() -> None:
    lines = []

    def log(s: str) -> None:
        print(s, flush=True)
        lines.append(s)

    log("=== graph-native He CI: settle the n_max=6/7 monotone question ===")
    log(f"E_exact = {E_EXACT} Ha; anchor: n_max={ANCHOR_NMAX} must reproduce "
        f"{ANCHOR_PCT}% (+-{ANCHOR_TOL})")

    E5, p5, dim5, t5 = err_pct(ANCHOR_NMAX)
    log(f"n_max=5: dim={dim5}  E0={E5:.8f}  err={p5:.5f}%  [{t5:.0f}s]")
    if abs(p5 - ANCHOR_PCT) > ANCHOR_TOL:
        log(f"ABORT: n_max=5 gave {p5:.5f}%, not ~{ANCHOR_PCT}% -- this pipeline "
            f"is NOT the 0.19% graph-native headline route. Reporting nothing "
            f"for 6/7; the validation_benchmarks.md flag stands.")
        _write(lines)
        return

    log("anchor PASSED -- pipeline validated, computing 6 and 7")
    for n_max in (6, 7):
        E, p, dim, t = err_pct(n_max)
        log(f"n_max={n_max}: dim={dim}  E0={E:.8f}  err={p:.5f}%  [{t:.0f}s]")

    log("")
    log("READING (compare to validation_benchmarks.md rows n_max=7 '<0.20%', "
        "n_max=8 '0.207%', n_max=9 '0.201%'):")
    log("  a monotone-decreasing ladder is expected; the measured n_max=7 above "
        "says which of the three recorded rows is inconsistent. Do NOT edit the "
        "benchmark values from this alone -- report to the PM for adjudication.")
    _write(lines)


def _write(lines: list[str]) -> None:
    import io
    io.open(OUT, "w", encoding="utf-8", newline="").write("\n".join(lines) + "\n")
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
