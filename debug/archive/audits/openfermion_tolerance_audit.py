"""Audit: does OpenFermion's EQ_TOLERANCE (1e-8) accumulation drop change any library Pauli count?

SymbolicOperator.__iadd__ deletes a term whenever its running coefficient falls below 1e-8, so a term
assembled from small pieces can vanish (found on Paper 62 water, CHANGELOG v6.7.16). This script builds each
system twice -- as shipped, and with the drop threshold lowered to 1e-13 -- and compares the qubit operators.

Usage: python debug/openfermion_tolerance_audit.py [--systems LiH,H2O,...] [--out debug/data/...json]
"""
import argparse
import contextlib
import json
import time

import numpy as np
import openfermion.ops.operators.symbolic_operator as so

from geovac import ecosystem_export as ee

TIGHT = 1e-13


@contextlib.contextmanager
def tolerance(tol: float):
    saved = so.SymbolicOperator.__dict__["_issmall"]  # the staticmethod object itself
    orig = saved.__func__

    def _issmall(val, tol=tol):
        return orig(val, tol)

    so.SymbolicOperator._issmall = staticmethod(_issmall)
    try:
        yield
    finally:
        so.SymbolicOperator._issmall = saved


def terms(qop) -> dict:
    return {k: complex(v) for k, v in qop.terms.items()}


def compare(a: dict, b: dict) -> dict:
    """a = shipped, b = tight. Report terms missing from a, and coefficient drift."""
    missing = [k for k in b if k not in a and abs(b[k]) >= TIGHT * 10]
    extra = [k for k in a if k not in b]
    common = [k for k in a if k in b]
    drift = max((abs(a[k] - b[k]) for k in common), default=0.0)
    return {
        "n_shipped": len(a),
        "n_tight": len(b),
        "n_tight_above_1e-12": sum(1 for v in b.values() if abs(v) >= 1e-12),
        "missing": len(missing),
        "missing_abs_sum": float(sum(abs(b[k]) for k in missing)),
        "missing_max": float(max((abs(b[k]) for k in missing), default=0.0)),
        "extra": len(extra),
        "max_coef_drift": float(drift),
        "norm_bound": float(sum(abs(b[k]) for k in missing) + sum(abs(a[k] - b[k]) for k in common)
                            + sum(abs(a[k]) for k in extra)),
        "missing_examples": [[str(k), float(abs(b[k]))] for k in missing[:4]],
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", default="")
    ap.add_argument("--out", default="debug/data/openfermion_tolerance_audit.json")
    a = ap.parse_args()
    names = a.systems.split(",") if a.systems else sorted(set(ee._SYSTEM_REGISTRY.values()))
    out = {}
    for name in names:
        t0 = time.time()
        try:
            shipped = terms(ee.hamiltonian(name).to_openfermion())
            with tolerance(TIGHT):
                tight = terms(ee.hamiltonian(name).to_openfermion())
            r = compare(shipped, tight)
        except Exception as e:  # report, keep going
            r = {"error": repr(e)}
        r["seconds"] = round(time.time() - t0, 1)
        out[name] = r
        print(name, {k: v for k, v in r.items() if k != "missing_examples"}, flush=True)
        with open(a.out, "w") as f:
            json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
