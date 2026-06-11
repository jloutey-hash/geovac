# -*- coding: utf-8 -*-
"""
B3 Phase-3 band-exhaustion probe -- TAG = legs
Characterizes j_max-dependence of D_max leg costs and baseline chain deficits
under band exhaustion (Paper-38 compression family, fixed continuum object).

Writes: debug/data/wh7_band_exh_legs.json
"""
import sys
import json
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh7_band_exhaustion_lib as lib

RNG_SEED = 42
np.random.seed(RNG_SEED)

THETA = 0.3
T_TOTAL = 1.0
BETA_DEFAULT = 1.0
BETA_PROBE = [0.5, 2.0]
REFS = ["null", "mixed"]

# ── helpers ──────────────────────────────────────────────────────────────────

def jkey(jmax):
    """Hashable string key for jmax (Sympy Rational)."""
    return str(jmax)


def flow_identity_check(cfg):
    """
    c12 should equal c23 by the flow-translation identity (KMS + flow commute
    at beta=1 under the wedge-KMS convention).  Returns |c12 - c23|.
    """
    c12 = lib.d_max(cfg["s1"], cfg["mid"])
    c23 = lib.d_max(cfg["mid"], cfg["s3"])
    return abs(c12 - c23)


def compute_row(jmax, ref_name, beta=BETA_DEFAULT):
    """Return dict of observables for one (jmax, ref, beta) triplet."""
    H = lib.reference_H(jmax, ref_name)
    h_norm = float(np.linalg.norm(H, 2))

    _, _, w = lib.wedge(jmax)
    wedge_dim = len(w)

    rho, Z = lib.kms_state(jmax, beta)
    # override: rebuild config at the requested beta
    import scipy.linalg as sl
    om0 = lib.conj(sl.expm(1j * THETA * H), rho)
    cfg = {
        "jmax": jmax,
        "Z": float(Z),
        "s1": om0,
        "mid": lib.conj(lib.flow_U(jmax, T_TOTAL / 2), om0),
        "s3": lib.conj(lib.flow_U(jmax, T_TOTAL), om0),
    }

    c12 = lib.d_max(cfg["s1"], cfg["mid"])
    c23 = lib.d_max(cfg["mid"], cfg["s3"])
    c13 = lib.d_max(cfg["s1"], cfg["s3"])
    deficit = c12 + c23 - c13
    id_err = abs(c12 - c23)

    return {
        "jmax": jkey(jmax),
        "ref": ref_name,
        "beta": beta,
        "wedge_dim": wedge_dim,
        "Z": float(Z),
        "H_norm": h_norm,
        "c12": c12,
        "c23": c23,
        "c13": c13,
        "deficit": deficit,
        "flow_id_err": id_err,
    }


# ── main sweep ────────────────────────────────────────────────────────────────

results = {}

print("=== B3 Phase-3 band-exhaustion legs probe ===")
print(f"{'jmax':>6}  {'ref':>6}  {'beta':>5}  {'wdim':>5}  "
      f"{'c12':>12}  {'c23':>12}  {'c13':>12}  "
      f"{'deficit':>12}  {'id_err':>10}  {'H_norm':>10}")

for ref in REFS:
    rows = []
    for jmax in lib.JMAX_LADDER:
        row = compute_row(jmax, ref, beta=BETA_DEFAULT)
        rows.append(row)
        print(f"{row['jmax']:>6}  {row['ref']:>6}  {row['beta']:>5.2f}  "
              f"{row['wedge_dim']:>5}  "
              f"{row['c12']:>12.6f}  {row['c23']:>12.6f}  {row['c13']:>12.6f}  "
              f"{row['deficit']:>12.6f}  {row['flow_id_err']:>10.2e}  "
              f"{row['H_norm']:>10.4f}")
    results[ref] = {"rows": rows}

# ── successive differences ────────────────────────────────────────────────────

print("\n=== Successive differences across ladder ===")
for ref in REFS:
    rows = results[ref]["rows"]
    diffs = {}
    keys = ["c12", "c23", "c13", "deficit", "H_norm"]
    for k in keys:
        vals = [r[k] for r in rows]
        d = [vals[i+1] - vals[i] for i in range(len(vals)-1)]
        diffs[k] = d
        print(f"  {ref:>6} {k:>10}: vals={[f'{v:.5f}' for v in vals]}  "
              f"diffs={[f'{x:+.5f}' for x in d]}")
    results[ref]["diffs"] = diffs

# ── convergence verdict ───────────────────────────────────────────────────────

def verdict(vals, diffs):
    """CONVERGENT / DIVERGENT / UNDECIDED."""
    if len(diffs) == 0:
        return "UNDECIDED"
    abs_diffs = [abs(d) for d in diffs]
    last_val = abs(vals[-1])
    # CONVERGENT: diffs decreasing AND last |diff| < 1e-2 * |last val|
    diffs_decreasing = all(abs_diffs[i] >= abs_diffs[i+1]
                           for i in range(len(abs_diffs)-1))
    last_small = (abs_diffs[-1] < 1e-2 * last_val) if last_val > 1e-14 else (abs_diffs[-1] < 1e-14)
    if diffs_decreasing and last_small:
        return "CONVERGENT"
    # DIVERGENT: values growing, diffs not shrinking
    vals_growing = all(vals[i] <= vals[i+1] for i in range(len(vals)-1))
    diffs_not_shrinking = not diffs_decreasing
    if vals_growing and diffs_not_shrinking:
        return "DIVERGENT"
    return "UNDECIDED"


print("\n=== Verdicts ===")
verdicts = {}
for ref in REFS:
    rows = results[ref]["rows"]
    diffs_dict = results[ref]["diffs"]
    verdicts[ref] = {}
    for k in ["c12", "c23", "c13", "deficit"]:
        vals = [r[k] for r in rows]
        diffs = diffs_dict[k]
        v = verdict(vals, diffs)
        verdicts[ref][k] = v
        print(f"  {ref:>6}  {k:>10}: {v}")
    results[ref]["verdicts"] = verdicts[ref]

# ── beta probe (null ref, two largest windows) ────────────────────────────────

print("\n=== Beta probe (ref=null, jmax=5/2 and 3) ===")
beta_probe_results = {}
for jmax in [lib.JMAX_LADDER[-2], lib.JMAX_LADDER[-1]]:  # 5/2 and 3
    beta_probe_results[jkey(jmax)] = {}
    for beta in BETA_PROBE:
        row = compute_row(jmax, "null", beta=beta)
        beta_probe_results[jkey(jmax)][beta] = row
        print(f"  jmax={row['jmax']:>4}  beta={beta:.2f}  "
              f"c12={row['c12']:.6f}  c23={row['c23']:.6f}  "
              f"c13={row['c13']:.6f}  deficit={row['deficit']:.6f}  "
              f"id_err={row['flow_id_err']:.2e}")

results["beta_probe"] = beta_probe_results

# ── flow-identity summary ─────────────────────────────────────────────────────

print("\n=== Flow-translation identity check (|c12 - c23|) ===")
for ref in REFS:
    rows = results[ref]["rows"]
    errs = [r["flow_id_err"] for r in rows]
    max_err = max(errs)
    print(f"  {ref}: max |c12-c23| across ladder = {max_err:.3e}  "
          f"(pass if < 1e-12: {'YES' if max_err < 1e-12 else 'NO -- inspect'})")
    results[ref]["flow_id_max_err"] = max_err

# ── save ──────────────────────────────────────────────────────────────────────

out_path = Path(__file__).resolve().parent / "data" / "wh7_band_exh_legs.json"
out_path.parent.mkdir(parents=True, exist_ok=True)

# Convert Sympy types for JSON
def to_json(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    raise TypeError(f"Not serializable: {type(obj)}")

with open(out_path, "w") as f:
    json.dump(results, f, default=to_json, indent=2)

print(f"\nSaved: {out_path}")
print("=== done ===")
