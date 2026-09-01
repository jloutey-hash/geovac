#!/usr/bin/env python
"""Re-measure Paper 14 tab:second_row from the live builders.

WHY
---
The Stage-4 ledger carried the balanced QWC columns as an open item ("they are
not pending: the builder returns n_qwc_balanced directly ... being measured for
the full library").  Both paper captions now say the column was measured
2026-08-30 and the standing "pending" note closed.  A previous run's "fixed" is
a claim, not a fact (the /qa hard rule), so this recomputes every printed cell
instead of trusting the caption.

Verdict 2026-08-31: the columns are CORRECT and reproduce exactly.

THE CONVENTION TRAP (the reason this script exists as a permanent artifact)
--------------------------------------------------------------------------
`build_balanced_hamiltonian` returns a dict in which the identity term is
handled INCONSISTENTLY between the two architectures it reports:

    res["N_pauli_composed"]  EXCLUDES the identity   (NaH: 558)
    res["N_pauli"]           INCLUDES the identity   (NaH: 575)
    res["one_norm"]          INCLUDES the identity   (NaH: 193.44)

The papers print Pauli counts WITH the identity and 1-norms WITHOUT it, so a
faithful re-measurement is:

    Pauli composed  = res["N_pauli_composed"] + 1
    Pauli balanced  = res["N_pauli"]
    lambda_ni       = sum |c| over res["qubit_op"] terms, identity EXCLUDED
    QWC balanced    = res["n_qwc_balanced"]

Reading `one_norm` as lambda_ni gives 193.44 where the paper says 20.6, and
reading `N_pauli_composed` as the printed count is off by exactly one.  Neither
is a paper defect; both are the identity-in / identity-out convention this
corpus has been bitten by before (C21 check E).  Anyone re-measuring these
tables from the returned dict must apply the four lines above.

Usage:
    python debug/qa/verify_balanced_qwc.py            # all six rows
    python debug/qa/verify_balanced_qwc.py --quick    # NaH + MgH2 only
"""
from __future__ import annotations

import argparse
import sys
import time

from geovac import molecular_spec as MS
from geovac.balanced_coupled import build_balanced_hamiltonian

# Paper 14 tab:second_row, as printed (2026-08-30 measurement).
# molecule -> (Q, Pauli_composed, Pauli_balanced, lambda_ni_balanced, QWC_bal)
PRINTED = {
    "NaH":   (20,   559,    575,  20.6,    69),
    "MgH2":  (40,  1117,   4861, 110.5,   903),
    "HCl":   (50,  1396,   9824, 869.7,  1173),
    "H2S":   (60,  1675,  13863, 879.6,  1282),
    "PH3":   (70,  1954,  18854, 895.3,  1477),
    "SiH4":  (80,  2233,  24745, 914.1,  1551),
}

SPECS = {
    "NaH":  MS.nah_spec,
    "MgH2": MS.mgh2_spec,
    "HCl":  MS.hcl_spec,
    "H2S":  MS.h2s_spec,
    "PH3":  MS.ph3_spec,
    "SiH4": MS.sih4_spec,
}


def lambda_ni(qubit_op) -> float:
    """1-norm EXCLUDING the identity -- the convention the papers print."""
    return sum(abs(c) for k, c in qubit_op.terms.items() if k != ())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="only NaH and MgH2 (the two cheapest rows)")
    args = ap.parse_args()

    names = ["NaH", "MgH2"] if args.quick else list(PRINTED)

    print(f"{'mol':<6} {'Q':>3}  {'quantity':<14} {'printed':>10} "
          f"{'measured':>12}  verdict")
    print("-" * 66)

    bad = 0
    for name in names:
        t0 = time.perf_counter()
        try:
            res = build_balanced_hamiltonian(SPECS[name](), verbose=False)
        except Exception as exc:                       # noqa: BLE001
            print(f"{name:<6}  ERROR building: {type(exc).__name__}: {exc}")
            bad += 1
            continue
        dt = time.perf_counter() - t0

        Q, p_comp, p_bal, lam, qwc = PRINTED[name]
        checks = [
            ("Q",            Q,      res["Q"]),
            # +1: the composed count is reported non-identity, printed with it
            ("Pauli comp.",  p_comp, res["N_pauli_composed"] + 1),
            ("Pauli bal.",   p_bal,  res["N_pauli"]),
            ("QWC bal.",     qwc,    res["n_qwc_balanced"]),
        ]
        for label, printed, measured in checks:
            ok = printed == measured
            bad += (not ok)
            print(f"{name:<6} {Q:>3}  {label:<14} {printed:>10} "
                  f"{measured:>12}  {'ok' if ok else 'MISMATCH'}")

        got = lambda_ni(res["qubit_op"])
        ok = abs(round(got, 1) - lam) < 0.05
        bad += (not ok)
        print(f"{name:<6} {Q:>3}  {'lambda_ni bal.':<14} {lam:>10.1f} "
              f"{got:>12.4f}  {'ok' if ok else 'MISMATCH'}   [{dt:.1f}s]")
        print()

    if bad:
        print(f"RESULT: FAIL ({bad} cell(s) disagree with the live builder)")
        return 1
    print(f"RESULT: PASS (every printed cell in {len(names)} row(s) "
          f"reproduces from the builder)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
