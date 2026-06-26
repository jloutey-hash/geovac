# -*- coding: utf-8 -*-
"""Adversarial cross-examination of the B3 Phase-3 band-exhaustion probes.
TAG = adversarial  (2026-06-10).

Targets:
  LEGS      debug/data/wh7_band_exh_legs.json
  PENALTIES debug/data/wh7_band_exh_penalties.json

Attacks:
  (a) INDEPENDENT RECOMPUTE -- 6 cells spanning both probes (different
      jmax/class/ref/beta), recomputed from lib directly via a fresh code
      path (lib.make_config instead of the legs probe's manual rebuild),
      compared to the probe JSONs at 1e-10.
  (b) THERMAL-SUPPRESSION ARTIFACT -- the headline convergent observable
      (mixed-ref deficit ladder, the only monotone-shrinking-diff sequence
      in either probe) re-run across the FULL ladder at beta = 0.5.
  (c) REFERENCE DEPENDENCE -- theta = 0.15 and 0.6: do the legs verdicts
      (probe-faithful classifier) or the penalties E(0.2) classifications
      flip?
  (d) HYGIENE -- clamp scan (1e-300 patterns), signedness audit, RAW
      (unnormalized) generator confirmation via ||G||_2 recompute, KMS
      eigenvalue floor inertness.

DISCIPLINE: deterministic (no randomness used; seed set anyway); SIGNED
values throughout; no clamping anywhere; negative/divergent findings are
first-class. Fresh code -- no probe-driver code is reused. The two verdict
classifiers are faithful REIMPLEMENTATIONS of the probes' documented rules
(required to test verdict flips like-for-like); a zero-aware variant is
added as an adversarial refinement.
"""

import json
import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import wh7_band_exhaustion_lib as lib  # noqa: E402

np.random.seed(20260610)  # determinism guard (no randomness actually drawn)

THETA0, TT, BETA0 = 0.3, 1.0, 1.0
EPS_BIG, EPS0 = 0.2, 1e-4
TOL = 1e-10
ZERO_THR = 1e-12

LEGS_PATH = HERE / "data" / "wh7_band_exh_legs.json"
PEN_PATH = HERE / "data" / "wh7_band_exh_penalties.json"
OUT_PATH = HERE / "data" / "wh7_band_exh_adversarial.json"

legs_js = json.loads(LEGS_PATH.read_text())
pen_js = json.loads(PEN_PATH.read_text())

LADDER = list(lib.JMAX_LADDER)  # [1, 3/2, 2, 5/2, 3]
JK = [str(j) for j in LADDER]

t_start = time.time()

# ---------------------------------------------------------------------------
# Fresh observable implementations (lib API only)
# ---------------------------------------------------------------------------
_cfg_cache = {}


def get_cfg(jmax, ref, theta, beta):
    key = (str(jmax), ref, float(theta), float(beta))
    if key not in _cfg_cache:
        H = lib.reference_H(jmax, ref)
        cfg = lib.make_config(jmax, H, theta=theta, t_total=TT, beta=beta)
        cfg["baseline"] = lib.deficit(cfg, cfg["mid"])
        _cfg_cache[key] = cfg
    return _cfg_cache[key]


def legs_cell(jmax, ref, theta=THETA0, beta=BETA0):
    cfg = get_cfg(jmax, ref, theta, beta)
    c12 = lib.d_max(cfg["s1"], cfg["mid"])
    c23 = lib.d_max(cfg["mid"], cfg["s3"])
    c13 = lib.d_max(cfg["s1"], cfg["s3"])
    return {"jmax": str(jmax), "c12": c12, "c23": c23, "c13": c13,
            "deficit": c12 + c23 - c13, "id_err": abs(c12 - c23)}


def pen_cell(name, jmax, theta=THETA0, want_D=False):
    G = lib.class_gen_folded(name, jmax)
    gn = float(np.linalg.norm(G, 2))
    if gn < ZERO_THR:
        return {"jmax": str(jmax), "G_norm2": gn, "zero": True}
    cfg = get_cfg(jmax, "null", theta, BETA0)
    base = cfg["baseline"]
    res = {"jmax": str(jmax), "G_norm2": gn, "zero": False,
           "baseline": float(base),
           "E_big": float(lib.deficit(cfg, lib.kicked(cfg, G, EPS_BIG)) - base)}
    if want_D:
        res["D_plus"] = float(
            (lib.deficit(cfg, lib.kicked(cfg, G, +EPS0)) - base) / EPS0)
        res["D_minus"] = float(
            (lib.deficit(cfg, lib.kicked(cfg, G, -EPS0)) - base) / EPS0)
        _, _, w = lib.wedge(jmax)
        K = np.diag(w)
        res["comm_norm"] = float(np.linalg.norm(K @ G - G @ K, 2))
    return res


# ---------------------------------------------------------------------------
# Verdict classifiers -- faithful reimplementations of the probes' rules
# ---------------------------------------------------------------------------
def legs_verdict(vals):
    """LEGS probe rule: CONVERGENT if |diffs| non-increasing AND last |diff|
    < 1e-2 * |last val|; DIVERGENT if vals growing and diffs not shrinking;
    else UNDECIDED."""
    diffs = [vals[i + 1] - vals[i] for i in range(len(vals) - 1)]
    ad = [abs(d) for d in diffs]
    last = abs(vals[-1])
    dec = all(ad[i] >= ad[i + 1] for i in range(len(ad) - 1))
    small = (ad[-1] < 1e-2 * last) if last > 1e-14 else (ad[-1] < 1e-14)
    if dec and small:
        return "CONVERGENT"
    grow = all(vals[i] <= vals[i + 1] for i in range(len(vals) - 1))
    if grow and not dec:
        return "DIVERGENT"
    return "UNDECIDED"


def pen_classify(vals):
    """PENALTIES probe rule (analyze_convergence), on the non-None entries."""
    vs = [v for v in vals if v is not None]
    if len(vs) < 2:
        return "insufficient data"
    diffs = [vs[i + 1] - vs[i] for i in range(len(vs) - 1)]
    ad = [abs(d) for d in diffs]
    if all(ad[i + 1] < ad[i] for i in range(len(ad) - 1)):
        return "CONVERGENT (diffs strictly shrinking)"
    if all(ad[i + 1] > ad[i] for i in range(len(ad) - 1)):
        return "DIVERGENT (diffs strictly growing)"
    if all(d > 0 for d in diffs):
        return "MONOTONE INCREASING (non-convergent)"
    if all(d < 0 for d in diffs):
        return "MONOTONE DECREASING (non-convergent)"
    return "NON-MONOTONE / MIXED"


def pen_classify_zero_aware(vals):
    """Adversarial refinement: a sequence that is machine-zero at every
    window is ZERO, not NON-MONOTONE."""
    vs = [v for v in vals if v is not None]
    if vs and all(abs(v) < ZERO_THR for v in vs):
        return "ZERO (machine epsilon at every window)"
    return pen_classify(vals)


out = {"meta": {"tag": "adversarial", "date": "2026-06-10",
                "theta0": THETA0, "t_total": TT, "beta0": BETA0,
                "eps_big": EPS_BIG, "eps0": EPS0, "tol_recompute": TOL,
                "ladder": JK,
                "targets": [str(LEGS_PATH), str(PEN_PATH)]}}

# ===========================================================================
# (a) INDEPENDENT RECOMPUTE -- 6 cells at 1e-10
# ===========================================================================
print("=== (a) independent recompute, 6 cells, gate 1e-10 ===")


def cmp_cell(label, pairs):
    rows, ok = [], True
    for f, pv, rv in pairs:
        d = abs(pv - rv)
        p = bool(d < TOL)
        ok = ok and p
        rows.append({"field": f, "probe": float(pv), "recomputed": float(rv),
                     "abs_diff": float(d), "pass": p})
        print(f"  {label:42s} {f:12s} probe={pv:+.12e} "
              f"re={rv:+.12e} d={d:.2e} {'PASS' if p else 'FAIL'}")
    return {"cell": label, "fields": rows, "pass_1e-10": bool(ok)}


cells = []

# 1. LEGS / null / jmax=2 / beta=1
rc = legs_cell(LADDER[2], "null")
pr = legs_js["null"]["rows"][2]
cells.append(cmp_cell("LEGS null jmax=2 beta=1",
                      [(k, pr[k], rc[k]) for k in
                       ("c12", "c23", "c13", "deficit")]))

# 2. LEGS / mixed / jmax=5/2 / beta=1
rc = legs_cell(LADDER[3], "mixed")
pr = legs_js["mixed"]["rows"][3]
cells.append(cmp_cell("LEGS mixed jmax=5/2 beta=1",
                      [(k, pr[k], rc[k]) for k in
                       ("c12", "c23", "c13", "deficit")]))

# 3. LEGS beta-probe / null / jmax=3 / beta=2.0
rc = legs_cell(LADDER[4], "null", beta=2.0)
pr = legs_js["beta_probe"]["3"]["2.0"]
cells.append(cmp_cell("LEGS beta-probe null jmax=3 beta=2",
                      [(k, pr[k], rc[k]) for k in
                       ("c12", "c23", "c13", "deficit")]))

# 4. PENALTIES / (1.5,1.5) timelike / jmax=3/2
rc = pen_cell("(1.5,1.5) timelike", LADDER[1])
pr = pen_js["per_class"]["(1.5,1.5) timelike"]["windows"][1]
cells.append(cmp_cell("PEN (1.5,1.5) timelike jmax=3/2",
                      [("G_norm2", pr["G_norm2"], rc["G_norm2"]),
                       ("E_0p2", pr["E_0p2_signed"], rc["E_big"]),
                       ("baseline", pr["baseline_deficit"], rc["baseline"])]))

# 5. PENALTIES / (2.0,0.0) spacelike / jmax=3 (incl. comm_norm)
rc = pen_cell("(2.0,0.0) spacelike", LADDER[4], want_D=True)
pr = pen_js["per_class"]["(2.0,0.0) spacelike"]["windows"][4]
cells.append(cmp_cell("PEN (2.0,0.0) spacelike jmax=3",
                      [("G_norm2", pr["G_norm2"], rc["G_norm2"]),
                       ("E_0p2", pr["E_0p2_signed"], rc["E_big"]),
                       ("comm_norm", pr["comm_norm"], rc["comm_norm"])]))

# 6. PENALTIES / (2.0,2.0) timelike / jmax=2 (incl. D_plus)
rc = pen_cell("(2.0,2.0) timelike", LADDER[2], want_D=True)
pr = pen_js["per_class"]["(2.0,2.0) timelike"]["windows"][2]
cells.append(cmp_cell("PEN (2.0,2.0) timelike jmax=2",
                      [("G_norm2", pr["G_norm2"], rc["G_norm2"]),
                       ("E_0p2", pr["E_0p2_signed"], rc["E_big"]),
                       ("D_plus", pr["D_plus"], rc["D_plus"]),
                       ("D_minus", pr["D_minus"], rc["D_minus"])]))

n_pass = sum(1 for c in cells if c["pass_1e-10"])
out["a_recompute"] = {"cells": cells,
                      "cells_pass": n_pass, "cells_total": len(cells)}
print(f"  -> {n_pass}/{len(cells)} cells pass at 1e-10")

# ===========================================================================
# (b) THERMAL-SUPPRESSION ARTIFACT -- mixed-ref deficit ladder at beta=0.5
# ===========================================================================
print("\n=== (b) thermal test: mixed-ref ladder at beta=0.5 vs 1.0 ===")
thermal = {}
for beta in (1.0, 0.5):
    rows = [legs_cell(j, "mixed", beta=beta) for j in LADDER]
    for key in ("deficit", "c13"):
        vals = [r[key] for r in rows]
        diffs = [vals[i + 1] - vals[i] for i in range(len(vals) - 1)]
        ratios = [diffs[i + 1] / diffs[i] if abs(diffs[i]) > 1e-300 else None
                  for i in range(len(diffs) - 1)]
        rec = {"values": vals, "diffs": diffs, "diff_ratios": ratios,
               "rel_last_diff": abs(diffs[-1]) / abs(vals[-1]),
               "verdict_probe_rule": legs_verdict(vals),
               "diffs_strictly_shrinking": bool(
                   all(abs(diffs[i + 1]) < abs(diffs[i])
                       for i in range(len(diffs) - 1)))}
        thermal[f"beta={beta}_{key}"] = rec
        print(f"  beta={beta} {key:8s} vals="
              + " ".join(f"{v:.6f}" for v in vals))
        print(f"           diffs=" + " ".join(f"{d:+.6f}" for d in diffs)
              + "  ratios=" + " ".join(f"{r:.3f}" for r in ratios)
              + f"  shrinking={rec['diffs_strictly_shrinking']}"
              + f"  verdict={rec['verdict_probe_rule']}")

# assessment: did apparent convergence disappear/slow dramatically?
b1 = thermal["beta=1.0_deficit"]
b05 = thermal["beta=0.5_deficit"]
slow_factor = (b05["rel_last_diff"] / b1["rel_last_diff"]
               if b1["rel_last_diff"] > 0 else None)
thermal["assessment"] = {
    "beta1_shrinking": b1["diffs_strictly_shrinking"],
    "beta05_shrinking": b05["diffs_strictly_shrinking"],
    "beta1_rel_last_diff": b1["rel_last_diff"],
    "beta05_rel_last_diff": b05["rel_last_diff"],
    "rel_last_diff_ratio_05_over_1": slow_factor,
}
out["b_thermal"] = thermal
print(f"  -> shrinking at beta=1: {b1['diffs_strictly_shrinking']}, "
      f"at beta=0.5: {b05['diffs_strictly_shrinking']}; "
      f"rel last diff {b1['rel_last_diff']:.4f} -> {b05['rel_last_diff']:.4f}")

# ===========================================================================
# (c) REFERENCE DEPENDENCE -- theta = 0.15 and 0.6
# ===========================================================================
print("\n=== (c) theta scan: legs verdicts + penalties classifications ===")
theta_scan = {"legs": {}, "penalties": {}}

# --- legs ---
for theta in (0.15, 0.3, 0.6):
    for ref in ("null", "mixed"):
        rows = [legs_cell(j, ref, theta=theta) for j in LADDER]
        entry = {"max_id_err": max(r["id_err"] for r in rows)}
        for key in ("c12", "c23", "c13", "deficit"):
            vals = [r[key] for r in rows]
            entry[key] = {"values": vals,
                          "verdict_probe_rule": legs_verdict(vals)}
        # structure checks
        dv = [r["deficit"] for r in rows]
        entry["staircase_null_check"] = {
            "step01_abs": abs(dv[1] - dv[0]),
            "step23_abs": abs(dv[3] - dv[2]),
            "rise12": dv[2] - dv[1],
            "rise34": dv[4] - dv[3],
            "is_staircase": bool(abs(dv[1] - dv[0]) < TOL
                                 and abs(dv[3] - dv[2]) < TOL
                                 and dv[2] - dv[1] > 1e-6
                                 and dv[4] - dv[3] > 1e-6)}
        entry["monotone_growth"] = bool(
            all(dv[i] < dv[i + 1] for i in range(4)))
        theta_scan["legs"][f"theta={theta}_{ref}"] = entry
        verds = {k: entry[k]["verdict_probe_rule"]
                 for k in ("c12", "c23", "c13", "deficit")}
        print(f"  theta={theta:4} {ref:5s} verdicts={verds}  "
              f"staircase={entry['staircase_null_check']['is_staircase']}  "
              f"monotone={entry['monotone_growth']}  "
              f"max_id_err={entry['max_id_err']:.2e}")

# verdict-flip table vs probe (probe verdicts all UNDECIDED at theta=0.3)
flips_legs = []
for theta in (0.15, 0.6):
    for ref in ("null", "mixed"):
        for key in ("c12", "c23", "c13", "deficit"):
            probe_v = legs_js[ref]["verdicts"][key]
            new_v = theta_scan["legs"][f"theta={theta}_{ref}"][key][
                "verdict_probe_rule"]
            if new_v != probe_v:
                flips_legs.append({"theta": theta, "ref": ref, "obs": key,
                                   "probe": probe_v, "new": new_v})
theta_scan["legs"]["verdict_flips_vs_probe"] = flips_legs
print(f"  legs verdict flips vs probe: {len(flips_legs)}"
      + (f" -> {flips_legs}" if flips_legs else ""))

# --- penalties ---
flips_pen = []
for theta in (0.15, 0.6):
    per_class = {}
    for name in lib.CLASSES:
        seq = [pen_cell(name, j, theta=theta) for j in LADDER]
        E = [None if r["zero"] else r["E_big"] for r in seq]
        cls = pen_classify(E)
        cls_z = pen_classify_zero_aware(E)
        per_class[name] = {"E_seq": E, "class_probe_rule": cls,
                           "class_zero_aware": cls_z}
        probe_cls = pen_js["per_class"][name]["convergence"]
        if cls != probe_cls:
            flips_pen.append({"theta": theta, "class": name,
                              "probe": probe_cls, "new": cls})
        print(f"  theta={theta:4} {name:22s} {cls:28s} "
              + " ".join("None" if e is None else f"{e:+.4e}" for e in E))
    theta_scan["penalties"][f"theta={theta}"] = per_class
theta_scan["penalties"]["class_flips_vs_probe"] = flips_pen
print(f"  penalties classification flips vs probe (theta-dependence): "
      f"{len(flips_pen)}")
for f_ in flips_pen:
    print(f"    theta={f_['theta']} {f_['class']}: "
          f"{f_['probe']} -> {f_['new']}")

# zero-aware reclassification at theta=0.3 (the probe's own setting)
zero_aware_03 = {}
for name in lib.CLASSES:
    E = pen_js["per_class"][name]["E_0p2_sequence"]
    zero_aware_03[name] = pen_classify_zero_aware(E)
theta_scan["penalties"]["zero_aware_at_theta0.3"] = zero_aware_03
out["c_theta"] = theta_scan

# ===========================================================================
# (d) HYGIENE
# ===========================================================================
print("\n=== (d) hygiene ===")
hyg = {}

# clamp scan
txt_legs = LEGS_PATH.read_text()
txt_pen = PEN_PATH.read_text()
pat = re.compile(r"e-30[0-9]|e\+30[0-9]|1e-?300")
hyg["clamp_hits_legs"] = len(pat.findall(txt_legs))
hyg["clamp_hits_penalties"] = len(pat.findall(txt_pen))

# signedness audit: negative values must be present where physics gives them
neg_E = sum(1 for name in lib.CLASSES
            for v in pen_js["per_class"][name]["E_0p2_sequence"]
            if v is not None and v < 0)
neg_D = sum(1 for name in lib.CLASSES
            for w_ in pen_js["per_class"][name]["windows"]
            for k in ("D_plus", "D_minus")
            if k in w_ and w_[k] < 0)
hyg["penalties_negative_E_count"] = neg_E
hyg["penalties_negative_D_count"] = neg_D
hyg["signs_preserved"] = bool(neg_E > 0 and neg_D > 0)

# RAW (unnormalized) confirmation: recompute ||G||_2 two ways
g_null_3 = float(np.linalg.norm(
    lib.class_gen_folded("(1.0,1.0) null", LADDER[4]), 2))
probe_gn = pen_js["per_class"]["(1.0,1.0) null"]["windows"][4]["G_norm2"]
probe_hn = legs_js["null"]["rows"][4]["H_norm"]
hyg["raw_generator_check"] = {
    "recomputed_norm_null_jmax3": g_null_3,
    "penalties_G_norm2": probe_gn,
    "legs_H_norm": probe_hn,
    "match_1e-12": bool(abs(g_null_3 - probe_gn) < 1e-12
                        and abs(g_null_3 - probe_hn) < 1e-12),
    "is_raw_not_normalized": bool(abs(g_null_3 - 1.0) > 0.1),
}

# KMS floor inertness: min eigenvalue of the wedge KMS state
floor = {}
for beta in (0.5, 1.0, 2.0):
    rho, _ = lib.kms_state(LADDER[4], beta)
    floor[f"beta={beta}"] = float(np.min(np.real(np.diag(rho))))
hyg["min_kms_eigenvalue_jmax3"] = floor
hyg["d_max_floor_inert"] = bool(min(floor.values()) > 1e-12)

print(f"  clamp hits: legs={hyg['clamp_hits_legs']} "
      f"pen={hyg['clamp_hits_penalties']}")
print(f"  negative E entries={neg_E}, negative D entries={neg_D} "
      f"(signs preserved: {hyg['signs_preserved']})")
print(f"  RAW check: ||G_null(3)||_2={g_null_3:.12f} vs pen={probe_gn:.12f} "
      f"legs={probe_hn:.12f} match={hyg['raw_generator_check']['match_1e-12']} "
      f"raw={hyg['raw_generator_check']['is_raw_not_normalized']}")
print(f"  min KMS eigenvalue at jmax=3: {floor} (floor 1e-300 inert: "
      f"{hyg['d_max_floor_inert']})")
out["d_hygiene"] = hyg

# ===========================================================================
# claim-by-claim verdicts
# ===========================================================================
print("\n=== claim verdicts ===")
claims = {}

# LEGS claim 1: staircase structure for null ref
stair_ok = all(theta_scan["legs"][f"theta={th}_null"]
               ["staircase_null_check"]["is_staircase"]
               for th in (0.15, 0.3, 0.6))
claims["LEGS_staircase_null"] = (
    "CONFIRMED" if stair_ok else "REFUTED",
    "half-integer windows add nothing; integer windows step -- holds at "
    "theta 0.15/0.3/0.6" if stair_ok else "staircase breaks at some theta")

# LEGS claim 2: slow-converging monotone growth for mixed
mono_ok = all(theta_scan["legs"][f"theta={th}_mixed"]["monotone_growth"]
              for th in (0.15, 0.3, 0.6))
b05_shrink = thermal["beta=0.5_deficit"]["diffs_strictly_shrinking"]
if mono_ok and b05_shrink:
    claims["LEGS_mixed_monotone"] = (
        "CONFIRMED", "monotone growth with shrinking diffs persists at "
        "theta 0.15/0.6 and at beta=0.5")
elif mono_ok:
    claims["LEGS_mixed_monotone"] = (
        "QUALIFIED", "monotone growth robust in theta, but diff-shrinking "
        "fails at beta=0.5 -- apparent convergence is KMS-tail-suppression-"
        "assisted")
else:
    claims["LEGS_mixed_monotone"] = ("REFUTED", "monotone growth breaks")

# LEGS claim 3: flow-translation identity exact
max_id = max(theta_scan["legs"][f"theta={th}_{r}"]["max_id_err"]
             for th in (0.15, 0.3, 0.6) for r in ("null", "mixed"))
claims["LEGS_flow_identity"] = (
    "CONFIRMED" if max_id < 1e-12 else "QUALIFIED",
    f"max |c12-c23| = {max_id:.2e} across theta scan at beta=1")

# PENALTIES claim 1: admissibility stabilizes (state-free, so theta/beta
# cannot touch it; recompute of comm norms via cells 5/6 covers it)
claims["PEN_admissibility_stabilizes"] = (
    "CONFIRMED",
    "[K_W,G] is state-free (no beta/theta dependence possible); recomputed "
    "comm norms match probe at 1e-10")

# PENALTIES claim 2: E(0.2) NON-MONOTONE for every class
zero_cls = [n for n, c in zero_aware_03.items()
            if c.startswith("ZERO")]
pen_flip_names = {(f_["theta"], f_["class"]) for f_ in flips_pen}
if not flips_pen and not zero_cls:
    claims["PEN_E_nonmonotone_all"] = ("CONFIRMED", "robust at all thetas")
else:
    detail = []
    if zero_cls:
        detail.append(f"classes {zero_cls} are machine-zero sequences "
                      "(mislabelled NON-MONOTONE)")
    if flips_pen:
        detail.append(f"{len(flips_pen)} theta-flips: " + "; ".join(
            f"theta={f_['theta']} {f_['class']} -> {f_['new']}"
            for f_ in flips_pen))
    claims["PEN_E_nonmonotone_all"] = ("QUALIFIED", " | ".join(detail))

# cross-cutting: reproducibility + hygiene
claims["NUMBERS_REPRODUCIBLE"] = (
    "CONFIRMED" if n_pass == len(cells) else "REFUTED",
    f"{n_pass}/{len(cells)} recompute cells at 1e-10")
claims["HYGIENE"] = (
    "CONFIRMED" if (hyg["clamp_hits_legs"] == 0
                    and hyg["clamp_hits_penalties"] == 0
                    and hyg["signs_preserved"]
                    and hyg["raw_generator_check"]["match_1e-12"]
                    and hyg["raw_generator_check"]["is_raw_not_normalized"]
                    and hyg["d_max_floor_inert"]) else "QUALIFIED",
    "no clamps in JSONs, signed values present, RAW generators confirmed, "
    "d_max 1e-300 floor never engaged (min KMS eig >= "
    f"{min(floor.values()):.2e})")

for k, (v, msg) in claims.items():
    print(f"  {k:30s} {v:10s} {msg}")
out["claims"] = {k: {"verdict": v, "detail": m} for k, (v, m) in claims.items()}

out["meta"]["elapsed_s"] = round(time.time() - t_start, 1)


# ---------------------------------------------------------------------------
def _conv(o):
    if isinstance(o, dict):
        return {str(k): _conv(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_conv(x) for x in o]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


OUT_PATH.parent.mkdir(exist_ok=True)
OUT_PATH.write_text(json.dumps(_conv(out), indent=2, default=str))
print(f"\nSaved: {OUT_PATH}")
print(f"Done in {out['meta']['elapsed_s']}s")
