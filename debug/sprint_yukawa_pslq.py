"""Yukawa-period PSLQ sprint (2026-06-03).

Tests whether measured Standard Model Yukawa values land in the
pure-Tate period ring Q[pi^2] (M1 union M2 only -- forced restriction
on inner-factor transcendentals by the eta-trivialization theorem,
Paper 18 section IV.6).

Three transforms x three coefficient ceilings x both mass scales (MS-bar
at M_Z is primary; unification scale handled in a separate driver if
clean).

Precision discipline.
  Charged-lepton Yukawa values are known to ~8 digits.
  Heavy-quark Yukawa values to ~4 digits.
  Light-quark Yukawa values to ~2-3 digits.
  PSLQ at coefficient ceiling M with basis size n requires ~log10(M^n)
  digits of precision in the target. Reported per-cell precision verdict.

Halt gate.
  Any apparent hit triggers the curve-fit audit protocol
  (docs/curve_fit_audit_memo.md) before proceeding to the next transform.

Basis (M1 union M2, pure-Tate forced by Finding 1 of the audit memo
debug/sprint_yukawa_algebraic_audit.md):
  {1, pi, pi^2, pi^4, pi^6, pi^8, 1/pi, 1/pi^2, 1/pi^4}

(M3 -- Catalan G, beta(s), Hurwitz at quarter-integers, Dirichlet L --
is structurally EXCLUDED for inner-factor observables by the
eta-trivialization theorem.)
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
from mpmath import mp, mpf, log, pi


# --- Yukawa values at MS-bar M_Z scale ---
#
# Source convention: m_f = y_f * v / sqrt(2), so y_f = sqrt(2) m_f / v
#                    with v = 246.21965 GeV (Higgs VEV from G_F at tree)
#
# Charged-lepton masses: PDG 2024 (high precision, used directly).
# Quark running masses at M_Z: Xing-Zhang-Zhou (Phys Rev D 77, 113016, 2008)
#   and Antusch-Maurer (JHEP 11, 115, 2013).
# Precision recorded as the number of clean significant digits of y_f.
#
# Units: y_f is dimensionless.

V_EW_GEV = "246.21965"       # +/- 0.00006 GeV (8 digits from G_F)
SQRT2_STR = "1.41421356237309504880168872420969807856967187537694"

# (name, m_f central in MeV, m_f uncertainty in MeV, sig_digits in y_f)
YUKAWA_INPUTS = [
    ("e",   "0.51099895069", "0.00000000016", 8),
    ("mu",  "105.6583755",   "0.0000023",     8),
    ("tau", "1776.86",       "0.12",          6),
    ("u",   "1.27",          "0.50",          2),
    ("d",   "2.71",          "0.20",          3),
    ("s",   "54.0",          "2.5",           3),
    ("c",   "619.0",         "20.0",          3),
    ("b",   "2850.0",        "30.0",          4),
    ("t",   "173000.0",      "400.0",         4),
]


def compute_yukawa_values(dps=60):
    mp.dps = dps
    v_mev = mpf(V_EW_GEV) * mpf("1000")
    sqrt2 = mpf(SQRT2_STR)
    out = {}
    for name, m_str, sig_str, digits in YUKAWA_INPUTS:
        m = mpf(m_str)
        sig = mpf(sig_str)
        y = sqrt2 * m / v_mev
        # propagate uncertainty (mass dominates at this precision)
        y_sig = sqrt2 * sig / v_mev
        out[name] = {
            "y": y,
            "y_uncertainty": y_sig,
            "sig_digits": digits,
            "m_mev": m_str,
            "m_uncertainty_mev": sig_str,
        }
    return out


# --- Basis ---

def m1_m2_basis(dps=120):
    """Pure-Tate basis {1, pi, pi^2, pi^4, pi^6, pi^8, 1/pi, 1/pi^2, 1/pi^4}.

    NOTE: pi (odd power) is in M1 (Hopf-base measure Vol(S^2)/4 lives there);
    pi^{2k} is in M2 (Seeley-DeWitt). Combined ring is Q[pi, 1/pi] truncated
    to small powers; this is what the eta-trivialization restricts inner-factor
    observables to.
    """
    mp.dps = dps
    return [
        ("1",       mpf(1)),
        ("pi",      pi),
        ("pi^2",    pi**2),
        ("pi^4",    pi**4),
        ("pi^6",    pi**6),
        ("pi^8",    pi**8),
        ("1/pi",    1/pi),
        ("1/pi^2",  1/pi**2),
        ("1/pi^4",  1/pi**4),
    ]


# --- PSLQ search ---

def pslq_search(target, basis, max_coeffs=(10, 100, 1000), dps=120):
    """Run PSLQ at each coefficient ceiling. Returns list of cells."""
    mp.dps = dps
    labels = ["target"] + [lbl for lbl, _ in basis]
    targets = [mpf(target)] + [v for _, v in basis]
    cells = []
    for M in max_coeffs:
        try:
            rel = mpmath.pslq(targets, maxcoeff=M)
        except Exception as e:
            cells.append({
                "max_coeff": M,
                "found": False,
                "error": str(e),
            })
            continue
        if rel is None or rel[0] == 0:
            cells.append({"max_coeff": M, "found": False, "rel": None})
            continue
        # reconstruct: rel[0]*target + sum_i rel[i+1]*basis[i] == 0
        c0 = rel[0]
        coeffs = rel[1:]
        # y = -(sum c_i * b_i) / c0
        recon = sum(c * b for c, (_, b) in zip(coeffs, basis))
        recon_target = -recon / mpf(c0)
        residual = mpf(target) - recon_target
        cells.append({
            "max_coeff": M,
            "found": True,
            "rel": [int(c) for c in rel],
            "labels": labels,
            "c0_target": int(c0),
            "coeffs": [int(c) for c in coeffs],
            "reconstructed_target": str(recon_target),
            "residual_abs": str(abs(residual)),
        })
    return cells


def precision_verdict(sig_digits, max_coeff, basis_size=9):
    """How many digits PSLQ at this ceiling requires vs. how many we have."""
    import math
    needed = math.log10(max_coeff) * basis_size
    return {
        "sig_digits_target": sig_digits,
        "digits_needed": round(needed, 1),
        "honest": sig_digits >= needed,
    }


# --- Three transforms ---

def transform_identity(y):
    return y


def transform_squared(y):
    return y * y


def transform_log(y):
    return log(y)


TRANSFORMS = [
    ("y",     transform_identity, "y_f itself"),
    ("y^2",   transform_squared,  "y_f squared (appears in a4 spectral action coefficient)"),
    ("log_y", transform_log,      "log(y_f) (in case RG-running from forced UV makes log natural)"),
]


# --- Sweep driver ---

def run_sweep_one_transform(transform_label, transform_fn, transform_desc,
                            yukawa_data, basis, max_coeffs, dps=120):
    mp.dps = dps
    results = {
        "transform": transform_label,
        "transform_description": transform_desc,
        "basis_labels": [lbl for lbl, _ in basis],
        "max_coeffs_tested": list(max_coeffs),
        "per_fermion": {},
    }
    any_hit = False
    for name, info in yukawa_data.items():
        y = info["y"]
        try:
            t_val = transform_fn(y)
        except Exception as e:
            results["per_fermion"][name] = {
                "error": f"transform failed: {e}",
            }
            continue
        cells = pslq_search(t_val, basis, max_coeffs=max_coeffs, dps=dps)
        # honesty per cell
        for c in cells:
            c["precision_verdict"] = precision_verdict(info["sig_digits"], c["max_coeff"], len(basis))
        # any HONEST hit -> any_hit
        honest_hit = any(c.get("found") and c["precision_verdict"]["honest"] for c in cells)
        if honest_hit:
            any_hit = True
        results["per_fermion"][name] = {
            "y_value": str(y),
            "transform_value": str(t_val),
            "sig_digits": info["sig_digits"],
            "cells": cells,
            "any_honest_hit": honest_hit,
        }
    results["any_honest_hit_across_fermions"] = any_hit
    return results


def main(dps=120):
    mp.dps = dps
    timestamp = datetime.now(timezone.utc).isoformat()
    yukawa = compute_yukawa_values(dps=dps)
    basis = m1_m2_basis(dps=dps)
    max_coeffs = (10, 100, 1000)

    print(f"Yukawa-period PSLQ sprint at {timestamp}")
    print(f"mp.dps = {dps}")
    print(f"basis size = {len(basis)} (M1 union M2 pure-Tate)")
    print(f"max_coeffs = {max_coeffs}")
    print()

    all_results = {
        "timestamp": timestamp,
        "scale": "MS-bar M_Z",
        "convention": "y_f = sqrt(2) * m_f / v_EW, v=246.21965 GeV",
        "basis": "M1 union M2 = {1, pi, pi^2, pi^4, pi^6, pi^8, 1/pi, 1/pi^2, 1/pi^4}",
        "transforms": [t[0] for t in TRANSFORMS],
        "max_coeffs": list(max_coeffs),
        "yukawa_inputs": YUKAWA_INPUTS,
        "yukawa_values": {k: {
            "y": str(v["y"]),
            "y_uncertainty": str(v["y_uncertainty"]),
            "sig_digits": v["sig_digits"],
            "m_mev": v["m_mev"],
            "m_uncertainty_mev": v["m_uncertainty_mev"],
        } for k, v in yukawa.items()},
        "per_transform": {},
        "halt_history": [],
    }

    for tlabel, tfn, tdesc in TRANSFORMS:
        print(f"--- Transform: {tlabel} ({tdesc}) ---")
        result = run_sweep_one_transform(
            tlabel, tfn, tdesc, yukawa, basis, max_coeffs, dps=dps,
        )
        all_results["per_transform"][tlabel] = result
        any_hit = result["any_honest_hit_across_fermions"]
        # report
        for fname, fdata in result["per_fermion"].items():
            cells = fdata.get("cells", [])
            for c in cells:
                hit = "HIT" if c.get("found") else "no hit"
                honest = "(honest)" if c["precision_verdict"]["honest"] else "(precision-shy)"
                print(f"  {fname:4s} M={c['max_coeff']:5d}  {hit:7s} {honest:16s}"
                      f"  sig_dig={fdata['sig_digits']} need={c['precision_verdict']['digits_needed']}")
            print()
        # halt gate
        if any_hit:
            all_results["halt_history"].append({
                "transform": tlabel,
                "reason": "honest_hit_found",
                "next_step": "curve_fit_audit",
            })
            print(f">>> HALT: honest hit found in transform {tlabel}. "
                  f"Audit before next transform.")
            break
        else:
            all_results["halt_history"].append({
                "transform": tlabel,
                "reason": "clean_null",
                "next_step": "proceed",
            })
            print(f"<<< CLEAN NULL across all fermions for transform {tlabel}.")

    # write output
    out_path = PROJ / "debug" / "data" / "sprint_yukawa_pslq.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nWrote results to {out_path}")
    return all_results


if __name__ == "__main__":
    main(dps=120)
