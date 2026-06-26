# -*- coding: utf-8 -*-
"""Band-exhaustion driver: per-class kick PENALTIES and ADMISSIBILITY.
TAG = penalties  (2026-06-10, B3 Phase-3).

Protocol:
  - Reference: lib.reference_H(jmax, 'null')   theta=0.3  t_total=1.0  beta=1.0
  - For each class in lib.CLASSES x each jmax in lib.JMAX_LADDER:
      G = lib.class_gen_folded(name, jmax)   [RAW, unnormalized]
      ||G||_2
      baseline = deficit(cfg, mid)           [should be ~0 by construction]
      E(0.2) = deficit(cfg, kicked(cfg,G,0.2)) - baseline   [SIGNED]
      D_plus  = (deficit(cfg, kicked(cfg,G,+e0)) - baseline) / e0   e0=1e-4
      D_minus = (deficit(cfg, kicked(cfg,G,-e0)) - baseline) / e0
      evenness residual: max(|c12(+e0)-c23(-e0)|, |c12(-e0)-c23(+e0)|)
          c12(e) = d_max(s1, kicked(cfg,G,e))
          c23(e) = d_max(kicked(cfg,G,e), s3)
      [K_W, G] operator norm   (K_W = diag(w) in wedge basis)
  - Convergence check: signed successive differences in E(0.2) shrinking?
  - Admissibility stabilization: commuting vs not, window-by-window.
  - Decision gate: record convergent/divergent/non-monotone/stable/unstable.

Hard rules: deterministic (no randomness here), signed values, never clamp.
"""

import sys, json, time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh7_band_exhaustion_lib as lib

# ── constants ──────────────────────────────────────────────────────────────────
THETA   = 0.3
T_TOTAL = 1.0
BETA    = 1.0
EPS_BIG = 0.2          # for E(eps)
EPS0    = 1e-4         # one-sided quotients
COMM_THR = 1e-10       # commutes if [K,G] < this
ZERO_THR = 1e-12       # generator folds to zero

# ── helpers ────────────────────────────────────────────────────────────────────

def commutator_norm(G, w):
    K = np.diag(w)
    return float(np.linalg.norm(K @ G - G @ K, 2))


def run_class(name, jmax):
    """Return a dict of all observables for one (class, jmax) pair."""
    _, _, w = lib.wedge(jmax)
    G = lib.class_gen_folded(name, jmax)          # RAW folded
    gnorm = float(np.linalg.norm(G, 2))

    result = {"jmax": str(jmax), "wedge_dim": len(w), "G_norm2": gnorm}

    if gnorm < ZERO_THR:
        result["folded_to_zero"] = True
        result["comm_norm"] = 0.0
        result["admissible"] = "commutes (zero generator)"
        return result

    result["folded_to_zero"] = False

    # reference
    H_ref = lib.reference_H(jmax, "null")
    cfg   = lib.make_config(jmax, H_ref, theta=THETA, t_total=T_TOTAL, beta=BETA)

    # baseline: deficit at mid (should be ~0 for a geodesic midpoint, but measure it)
    baseline = lib.deficit(cfg, cfg["mid"])
    result["baseline_deficit"] = float(baseline)

    # E(0.2) signed
    s2_big = lib.kicked(cfg, G, EPS_BIG)
    E_big  = lib.deficit(cfg, s2_big) - baseline
    result["E_0p2_signed"] = float(E_big)

    # one-sided quotients at eps0
    s2_pp = lib.kicked(cfg, G, +EPS0)
    s2_pm = lib.kicked(cfg, G, -EPS0)
    d_plus  = (lib.deficit(cfg, s2_pp) - baseline) / EPS0
    d_minus = (lib.deficit(cfg, s2_pm) - baseline) / EPS0
    result["D_plus"]  = float(d_plus)
    result["D_minus"] = float(d_minus)

    # evenness residual
    c12_pp = lib.d_max(cfg["s1"], s2_pp)
    c23_pp = lib.d_max(s2_pp,     cfg["s3"])
    c12_pm = lib.d_max(cfg["s1"], s2_pm)
    c23_pm = lib.d_max(s2_pm,     cfg["s3"])
    ev_res = max(abs(c12_pp - c23_pm), abs(c12_pm - c23_pp))
    result["evenness_residual"] = float(ev_res)

    # commutator [K_W, G]
    cn = commutator_norm(G, w)
    result["comm_norm"] = float(cn)
    result["admissible"] = "commutes" if cn < COMM_THR else "non-commuting"

    return result


def analyze_convergence(vals):
    """Given ordered list of E(0.2) across the jmax ladder, classify convergence."""
    if len(vals) < 2:
        return "insufficient data"
    diffs = [vals[i+1] - vals[i] for i in range(len(vals)-1)]
    abs_diffs = [abs(d) for d in diffs]
    if all(abs_diffs[i+1] < abs_diffs[i] for i in range(len(abs_diffs)-1)):
        return "CONVERGENT (diffs strictly shrinking)"
    if all(abs_diffs[i+1] > abs_diffs[i] for i in range(len(abs_diffs)-1)):
        return "DIVERGENT (diffs strictly growing)"
    if all(d > 0 for d in diffs):
        return "MONOTONE INCREASING (non-convergent)"
    if all(d < 0 for d in diffs):
        return "MONOTONE DECREASING (non-convergent)"
    return "NON-MONOTONE / MIXED"


def analyze_admissibility(records):
    """Given ordered records (by jmax ladder), check if admissibility classification stabilizes."""
    states = [r["admissible"] for r in records if not r.get("folded_to_zero", False)]
    if not states:
        return "all zero"
    last = states[-1]
    # stabilizes if from some index onward they're all the same
    for start in range(len(states)):
        if all(s == states[start] for s in states[start:]):
            return f"STABLE from window {start} onward: {last}"
    return f"UNSTABLE (alternating): {states}"


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    t0 = time.time()
    output = {
        "meta": {
            "tag": "penalties",
            "date": "2026-06-10",
            "theta": THETA,
            "t_total": T_TOTAL,
            "beta": BETA,
            "eps_big": EPS_BIG,
            "eps0": EPS0,
            "comm_threshold": COMM_THR,
            "zero_threshold": ZERO_THR,
            "jmax_ladder": [str(j) for j in lib.JMAX_LADDER],
            "classes": list(lib.CLASSES.keys()),
        },
        "selftest": lib.selftest(),   # bake in validation check
        "per_class": {},
        "summary": {},
        "caveats": [],
    }

    # Sprint-3b admissibility cross-check table (expected from task spec)
    sb_expected = {
        # mu'=0 classes commute at EVERY window
        "(1.0,0.0) spacelike": "commutes_all",
        "(2.0,0.0) spacelike": "commutes_all",
        # (2.0,1.0) spacelike zero ONLY at jmax=1
        "(2.0,1.0) spacelike": "zero_at_1_only",
        # (2.0,2.0) timelike commutes ONLY at jmax=1
        "(2.0,2.0) timelike": "commutes_at_1_only",
        # others: no special expectation
        "(0.5,0.5) spacelike": "non-commuting_all",
        "(1.0,1.0) null":      "non-commuting_all",
        "(1.5,1.5) timelike":  "non-commuting_all",
    }

    for name in lib.CLASSES:
        print(f"\n  === {name} ===")
        class_rows = []
        for jmax in lib.JMAX_LADDER:
            r = run_class(name, jmax)
            class_rows.append(r)
            if r.get("folded_to_zero"):
                print(f"    jmax={jmax}  FOLDS TO ZERO")
            else:
                print(f"    jmax={jmax}  ||G||={r['G_norm2']:.4f}  "
                      f"E(0.2)={r['E_0p2_signed']:+.6f}  "
                      f"D+={r['D_plus']:+.4f}  D-={r['D_minus']:+.4f}  "
                      f"ev={r['evenness_residual']:.3e}  "
                      f"[K,G]={r['comm_norm']:.3e}  adm={r['admissible']}")

        # per-class analysis
        e_vals = [r["E_0p2_signed"] for r in class_rows if not r.get("folded_to_zero")]
        conv   = analyze_convergence(e_vals)
        adm    = analyze_admissibility(class_rows)

        # Sprint-3b cross-check
        comm_seq = [r["admissible"] for r in class_rows]
        zero_seq = [r.get("folded_to_zero", False) for r in class_rows]

        sb_check = "PASS"
        sb_detail = ""
        exp = sb_expected.get(name, "")
        if exp == "commutes_all":
            if not all("commutes" in s for s in comm_seq):
                sb_check = "FAIL"
                sb_detail = f"expected all commute, got {comm_seq}"
        elif exp == "zero_at_1_only":
            if not (zero_seq[0] and not any(zero_seq[1:])):
                sb_check = "FAIL"
                sb_detail = f"expected zero only at jmax=1, got {zero_seq}"
        elif exp == "commutes_at_1_only":
            # commutes at jmax=1 (or zero); non-commuting for rest
            ok1 = ("commutes" in comm_seq[0]) or zero_seq[0]
            ok_rest = all("non-commuting" in s for s in comm_seq[1:])
            if not (ok1 and ok_rest):
                sb_check = "FAIL"
                sb_detail = f"expected commute at 1 only, got {comm_seq}"
        elif exp == "non-commuting_all":
            if any("commutes" in s and "commutes (zero" not in s for s in comm_seq):
                sb_check = "FAIL"
                sb_detail = f"expected all non-commuting, got {comm_seq}"

        output["per_class"][name] = {
            "windows": class_rows,
            "E_0p2_sequence": [r["E_0p2_signed"] if not r.get("folded_to_zero") else None
                               for r in class_rows],
            "convergence": conv,
            "admissibility": adm,
            "sprint3b_check": sb_check,
            "sprint3b_detail": sb_detail,
        }
        print(f"    -> convergence: {conv}")
        print(f"    -> admissibility: {adm}")
        print(f"    -> Sprint-3b: {sb_check} {sb_detail}")

    # Global caveats
    output["caveats"] = [
        ("(1.0,1.0) null coincides with the reference direction (reference_H uses "
         "class_gen_folded('(1.0,1.0) null', jmax)). The kick G == H_ref at every "
         "window, so kicked(cfg, G, eps) rotates the KMS orbit in the SAME direction "
         "as the reference perturbation. E(0.2) measures deficit of a second rotation "
         "on top of the first; it is structurally different from a transverse kick. "
         "D_plus and D_minus are directional derivatives along the reference orbit "
         "itself and thus measure orbit curvature, not cross-direction penalty."),
        ("Baseline deficit is computed as deficit(cfg, cfg['mid']). For a true geodesic "
         "midpoint this should be zero; small nonzero values (|baseline| < 1e-6) are "
         "numerical. Large baseline would indicate make_config does not produce a "
         "geodesic triple. Values reported in per-window rows."),
        ("E(0.2) signed: positive = kick INCREASES the sum of d_max legs (penalty); "
         "negative = kick DECREASES it (improvement / shortcut). Negative E is a "
         "first-class physical result indicating the kicked state is a better triangle "
         "midpoint than the mid-orbit reference."),
        ("D_plus, D_minus are one-sided finite-difference derivatives at eps0=1e-4 "
         "relative to baseline. They approximate the directional derivative of the "
         "deficit functional in the +G and -G directions respectively."),
        ("Admissibility threshold: [K_W, G] < 1e-10 -> commutes. This is a strict "
         "criterion; values in the 1e-13 to 1e-15 range are machine-epsilon zeros."),
        ("Band exhaustion convention: RAW (unnormalized) generators. Growing jmax "
         "with FIXED (b, mu', mu) is an honest window-compression of a fixed SU(2) "
         "object. Norms grow; penalties scale accordingly. Normalized comparisons "
         "would require dividing by G_norm2 at each window."),
    ]

    # Summary table
    output["summary"] = {
        name: {
            "convergence": output["per_class"][name]["convergence"],
            "admissibility": output["per_class"][name]["admissibility"],
            "sprint3b_check": output["per_class"][name]["sprint3b_check"],
            "E_0p2_at_jmax3": (output["per_class"][name]["E_0p2_sequence"][-1]
                                if output["per_class"][name]["E_0p2_sequence"][-1] is not None
                                else "zero"),
        }
        for name in lib.CLASSES
    }

    elapsed = time.time() - t0
    output["meta"]["elapsed_s"] = round(elapsed, 1)
    print(f"\nDone in {elapsed:.1f}s")

    # Serialize (convert any sympy objects to strings in selftest dims)
    def convert(o):
        if isinstance(o, dict):
            return {str(k): convert(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [convert(x) for x in o]
        if isinstance(o, (np.floating, np.integer)):
            return float(o)
        if isinstance(o, float):
            return o
        if isinstance(o, bool):
            return o
        return o

    out_path = Path(__file__).resolve().parent / "data" / "wh7_band_exh_penalties.json"
    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(convert(output), f, indent=2, default=str)
    print(f"Results written to {out_path}")
    return output


if __name__ == "__main__":
    main()
