r"""Sweep the retired graph-native He "0.19% @ n_max=7" headline to the measured
0.216% (E = -2.89746 Ha).  17 loci (16 from the C17 gate + synthesis.done.md's
"0.19% CI", which the gate's graph-native require_nearby cannot see).

Basis: build_graph_native_fci at n_max=7 = 0.21559% (E -2.89746, dim 1218),
validated NON-circularly against paper_fci_atoms's own n_max=6 = 0.23% (measured
0.22864%).  The papers' 0.19% (E -2.8983) is a stale, non-monotone outlier no
current path reproduces.  n_max=6 = 0.23% is CORRECT and left alone; the
adiabatic "0.19-0.20%" FLOOR is a different solver and is not touched.

All anchors are RAW strings (no escape interpretation); every replace is
count-checked (exactly one match) and fail-loud.

Run:  python debug/qa/_sweep_graphnative019.py
"""
from __future__ import annotations

import io
import sys

EDITS = [
    # --- Paper 7: the stale ENERGY + percent ---
    (r"papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
     r"energy is $-2.8983$~Ha (0.19\% error vs.\ exact $-2.903724$~Ha),",
     r"energy is $-2.89746$~Ha (0.216\% error vs.\ exact $-2.903724$~Ha),"),

    # --- Paper 13 ---
    (r"papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
     r"with exact rational Slater integrals achieves 0.19\% at $n_{\max} = 7$",
     r"with exact rational Slater integrals achieves 0.216\% at $n_{\max} = 7$"),
    (r"papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
     r"0.19\% at $n_{\max} = 7$ with zero free parameters (no orbital",
     r"0.216\% at $n_{\max} = 7$ with zero free parameters (no orbital"),

    # --- FCI-atoms (keep n_max=6 = 0.23%, correct only n_max=7) ---
    (r"papers/group2_quantum_chemistry/paper_fci_atoms.tex",
     r"0.23\% at $n_{\max} = 6$, and at $n_{\max} = 7$ reaches 0.19\%",
     r"0.23\% at $n_{\max} = 6$, and at $n_{\max} = 7$ reaches 0.216\%"),

    # --- Paper 18 ---
    (r"papers/group3_foundations/paper_18_exchange_constants.tex",
     r"Graph-native FCI achieves 0.19\% at $n_{\max} = 7$",
     r"Graph-native FCI achieves 0.216\% at $n_{\max} = 7$"),

    # --- group2 synthesis ---
    (r"papers/synthesis/group2_quantum_chemistry_synthesis.tex",
     r"non-variational); $0.19\%$ (graph-native CI) & \cite{loutey_paper13}",
     r"non-variational); $0.216\%$ (graph-native CI) & \cite{loutey_paper13}"),
    (r"papers/synthesis/group2_quantum_chemistry_synthesis.tex",
     r"integrals and reaches $0.19\%$ at $n_{\max}=7$ with no free",
     r"integrals and reaches $0.216\%$ at $n_{\max}=7$ with no free"),
    (r"papers/synthesis/group2_quantum_chemistry_synthesis.tex",
     r"cusp, $0.19\%$ graph-native)~\cite{loutey_paper13}",
     r"cusp, $0.216\%$ graph-native)~\cite{loutey_paper13}"),

    # --- INDEX ---
    (r"papers/INDEX.md",
     r"graph-native CI at 0.19% with zero parameters |",
     r"graph-native CI at 0.216% (n_max=7) with zero parameters |"),

    # --- claims_register ---
    (r"docs/claims_register.md",
     r"0.19% (zero-parameter graph-native CI, $n_{\max}{=}7$) | Paper 13 | MEASURED | vs. Pekeris exact $-2.9037$ Ha (code-confirmed 2026-06-26) |",
     r"0.216% (zero-parameter graph-native CI, $n_{\max}{=}7$; re-measured 2026-09-19, the prior 0.19% was stale) | Paper 13 | MEASURED | vs. Pekeris exact $-2.9037$ Ha (n_max=7 = 0.216%, dim 1218) |"),

    # --- claim_test_matrix: the NO-TEST coverage-gap row, now MEASURED ---
    (r"docs/claim_test_matrix.md",
     r"| 13 | He graph-native 0.19% @ n_max=7 | — (tests stop at n_max=3; n_max=5=0.25% confirmed) | **NO-TEST** | coverage gap (n_max=7 ~uncomputable in CI; note in-paper) |",
     r"| 13 | He graph-native **0.216%** @ n_max=7 (MEASURED 2026-09-19; the retired 0.19% was a stale, non-monotone outlier) | `debug/qa/_graph_native_he_nmax67.py` (n_max=5/6/7 = 0.2496/0.22864/0.21559%, self-anchored to the paper's own n_max=6=0.23%) | **DEBUG-DRIVER ONLY** | still no committed CI test (n_max=7 build ~31 min); the driver + validation_benchmarks flag are the record |"),
    (r"docs/claim_test_matrix.md",
     r"| FCI-A | He 0.26% (grid hybrid) / 0.19% (graph-native) |",
     r"| FCI-A | He 0.26% (grid hybrid) / 0.216% (graph-native, n_max=7, MEASURED 2026-09-19) |"),

    # --- group2.done.md (frozen DoD watch-notes; self-declaring) ---
    (r"docs/qa/group2.done.md",
     r"graph-native CI **0.19%** at $n_{\max}=7$, zero parameters. (Note:",
     r"graph-native CI **0.216%** at $n_{\max}=7$ (MEASURED 2026-09-19; the 0.19% this line carried was stale), zero parameters. (Note:"),
    (r"docs/qa/group2.done.md",
     r"graph-native CI He 0.19% / Be 0.71% / Li 1.03% (zero-parameter, exact",
     r"graph-native CI He 0.216% (2026-09-19; was 0.19%) / Be 0.71% / Li 1.03% (zero-parameter, exact"),

    # --- synthesis.done.md (the "0.19% CI" the gate can't see) ---
    (r"docs/qa/synthesis.done.md",
     r"He 0.004% cusp / 0.022% raw / 0.19% CI;",
     r"He 0.004% cusp / 0.022% raw / 0.216% CI;"),

    # --- CLAUDE.md best-results (S2) and hierarchy (S5 numeric cell) ---
    (r"CLAUDE.md",
     r"| He (graph-native CI) | 0.19% | Zero-parameter, exact algebraic integrals, n_max=7 | 13 |",
     r"| He (graph-native CI) | 0.216% (n_max=7) | Zero-parameter, exact algebraic integrals | 13 |"),
    (r"CLAUDE.md",
     r"0.19% (graph-native CI n_max=7, 0 params, exact algebraic integrals) | 13 |",
     r"0.216% (graph-native CI n_max=7, 0 params, exact algebraic integrals) | 13 |"),
]


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    by_file: dict = {}
    for path, old, new in EDITS:
        by_file.setdefault(path, []).append((old, new))

    failures = []
    for path, pairs in by_file.items():
        s = io.open(path, encoding="utf-8").read()
        n = 0
        for old, new in pairs:
            c = s.count(old)
            if c != 1:
                failures.append(f"{path}: anchor count {c} (need 1) -> {old[:70]!r}")
                continue
            s = s.replace(old, new, 1)
            n += 1
        if n:
            io.open(path, "w", encoding="utf-8", newline="").write(s)
        print(f"  {path}: {n}/{len(pairs)} applied")

    if failures:
        print("\nFAILURES (nothing written for the failed anchors):")
        for f in failures:
            print("   " + f)
        raise SystemExit(1)
    print(f"\nall {len(EDITS)} anchored replacements applied cleanly")


if __name__ == "__main__":
    main()
