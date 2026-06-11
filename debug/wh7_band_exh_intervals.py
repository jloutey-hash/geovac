# -*- coding: utf-8 -*-
"""Band-exhaustion Phase-3: Interval functional and orbit geometry.
TAG = intervals
Writes debug/data/wh7_band_exh_intervals.json.

Reference: lib.reference_H(jmax, 'null'), theta=0.3, beta=1.0.

Structural finding (pre-run diagnostic):
  The 'null' class generator C^1_{1,0} has mu'=1, mu=0 (integer spin-1).
  In the folded wedge, this maps 2|m'| -> 2|m'+1|, so parity(weight) is
  PRESERVED (even stays even, odd stays odd).  The KMS state rho is diagonal
  on weight eigenstates, so om0 = e^{i theta H_null} rho e^{-...} inherits
  this even-weight parity: all odd-even cross blocks are exactly zero.
  Consequence: the orbit under K_W (weights all even) has effective period PI
  (not 2pi) -- confirmed numerically to < 5e-17.
  The mixed reference C^{1/2}_{1/2,1/2} crosses parity -> period = 2pi.
  We report both T_pi and T_2pi measurements for documentation; use the
  MEASURED period (pi for null) as the working interval for recovery/additivity.

Discipline (hard rules):
 - Deterministic: fixed offsets, no random state.
 - Report SIGNED values; never clamp negatives.
 - Negative or divergent findings are first-class results -- report them.
 - All final numbers come from the JSON (written at the end).
"""
import sys
import json
import numpy as np
from pathlib import Path
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh7_band_exhaustion_lib as lib

# ── helpers ─────────────────────────────────────────────────────────────────

def flow_state(jmax, rho0, t):
    """Evolve rho0 forward by time t under the wedge boost K_W."""
    U = lib.flow_U(jmax, t)
    return lib.conj(U, rho0)


def grid_parabolic_argmin(vals, ts):
    """Grid argmin + parabolic refinement. Returns (t_min, val_min)."""
    k = int(np.argmin(vals))
    if k == 0 or k == len(vals) - 1:
        return float(ts[k]), float(vals[k])
    ya, yb, yc = vals[k-1], vals[k], vals[k+1]
    denom = ya - 2*yb + yc
    if abs(denom) < 1e-15:
        return float(ts[k]), float(yb)
    frac = 0.5 * (ya - yc) / denom
    dt = float(ts[1]) - float(ts[0])
    t_min = float(ts[k]) + frac * dt
    val_min = float(yb) - (ya - yc)**2 / (8 * denom)
    return t_min, val_min


def interval_recovery(om_a, om_b, jmax, period, n_grid=600):
    """Recover Delta-tau between om_a and om_b via grid argmin + parabolic refine.
    Uses BOTH d_max-matching and trace-matching.
    Returns dict with recovered_dmax, recovered_trace, cost_dmax, cost_trace.
    """
    ts = np.linspace(0.0, period, n_grid, endpoint=False)
    dmax_vals = np.zeros(n_grid)
    td_vals   = np.zeros(n_grid)
    for i, t in enumerate(ts):
        om_at = flow_state(jmax, om_a, t)
        dmax_vals[i] = lib.d_max(om_at, om_b)
        td_vals[i]   = lib.trace_dist(om_at, om_b)

    t_dm, cost_dm = grid_parabolic_argmin(dmax_vals, ts)
    t_td, cost_td = grid_parabolic_argmin(td_vals, ts)
    return {
        "recovered_dmax":  t_dm,
        "recovered_trace": t_td,
        "cost_dmax":  cost_dm,
        "cost_trace": cost_td,
    }


def interval_length_dmax(om_a, om_b, jmax, period, n_grid=600):
    """Estimate ell(a, b) = recovered Delta-tau via d_max-matching."""
    r = interval_recovery(om_a, om_b, jmax, period, n_grid)
    return r["recovered_dmax"]


def measure_period(jmax, om0, max_T=4*np.pi, n_pts=800):
    """Measure first return period via trace_dist sweep.
    Returns (detected_period, td_at_pi, td_at_2pi, full_ts, full_vals).
    """
    T_pi  = float(np.pi)
    T_2pi = 2.0 * float(np.pi)
    ts = np.linspace(0.0, max_T, n_pts, endpoint=False)
    td_vals = np.array([lib.trace_dist(flow_state(jmax, om0, t), om0) for t in ts])
    td_pi  = float(lib.trace_dist(flow_state(jmax, om0, T_pi),  om0))
    td_2pi = float(lib.trace_dist(flow_state(jmax, om0, T_2pi), om0))

    # First return: skip t < 0.05 to avoid t=0, find first minimum below threshold
    skip = n_pts // 20
    thresh = max(1e-8, 1e-3 * td_vals[skip:].max())
    close_to_zero = np.where(td_vals[skip:] < thresh)[0]
    if len(close_to_zero) > 0:
        detected = float(ts[skip + close_to_zero[0]])
    else:
        detected = T_2pi  # fallback
    return detected, td_pi, td_2pi, ts.tolist(), td_vals.tolist()


# ── deterministic ordered triples for additivity ────────────────────────────
# Five triples (ta, Delta1, Delta2); all offsets < 0.25*period (no wraparound).
# Using fractions of a canonical scale=1.0; actual t = scale * period.
TRIPLE_SCALES = [
    (0.05, 0.10, 0.08),
    (0.08, 0.12, 0.10),
    (0.10, 0.15, 0.12),
    (0.12, 0.08, 0.15),
    (0.06, 0.11, 0.09),
]


def additivity_check(om0, jmax, period, triple_scales, n_grid=600):
    """For each (sa, sd1, sd2) in fraction-of-period units:
    build om_a=flow(om0, sa*period), om_b=flow(om_a, sd1*period),
    om_c=flow(om_b, sd2*period) and check ell(a,b)+ell(b,c)=ell(a,c).
    Returns (list of detail dicts, max_abs_error).
    """
    errors = []
    for (sa, sd1, sd2) in triple_scales:
        ta   = sa  * period
        d1   = sd1 * period
        d2   = sd2 * period
        om_a = flow_state(jmax, om0, ta)
        om_b = flow_state(jmax, om_a, d1)
        om_c = flow_state(jmax, om_b, d2)
        ell_ab = interval_length_dmax(om_a, om_b, jmax, period, n_grid)
        ell_bc = interval_length_dmax(om_b, om_c, jmax, period, n_grid)
        ell_ac = interval_length_dmax(om_a, om_c, jmax, period, n_grid)
        err = (ell_ab + ell_bc) - ell_ac   # SIGNED
        errors.append({
            "ta": float(ta), "d1": float(d1), "d2": float(d2),
            "true_d1":  float(d1),
            "true_d2":  float(d2),
            "ell_ab": float(ell_ab),
            "ell_bc": float(ell_bc),
            "ell_ac": float(ell_ac),
            "additivity_error_signed": float(err),
        })
    max_err = max(abs(e["additivity_error_signed"]) for e in errors)
    return errors, float(max_err)


# ── main loop ────────────────────────────────────────────────────────────────

def run():
    results = {}

    for jmax in lib.JMAX_LADDER:
        jkey = str(float(jmax))
        wlabs, V, w = lib.wedge(jmax)
        dim = len(w)
        print(f"\n{'='*65}")
        print(f"  jmax = {float(jmax)}  wedge_dim = {dim}")
        print(f"{'='*65}")

        rho0, Z = lib.kms_state(jmax, beta=1.0)
        H_ref = lib.reference_H(jmax, "null")
        om0 = lib.conj(expm(1j * 0.3 * H_ref), rho0)
        print(f"  Z = {Z:.6f}")

        # ── Check parity structure of H_ref and om0 ─────────────────────────
        odd_idx  = [i for i, ww in enumerate(w) if int(round(ww)) % 2 != 0]
        even_idx = [i for i, ww in enumerate(w) if int(round(ww)) % 2 == 0]
        h_parity_break = 0.0
        if len(odd_idx) > 0 and len(even_idx) > 0:
            h_parity_break = float(np.linalg.norm(H_ref[np.ix_(odd_idx, even_idx)]))
        om0_parity_break = 0.0
        if len(odd_idx) > 0 and len(even_idx) > 0:
            om0_parity_break = float(np.linalg.norm(om0[np.ix_(odd_idx, even_idx)]))
        print(f"  H_null odd-even block norm (parity check) = {h_parity_break:.4e}")
        print(f"  om0    odd-even block norm (parity check) = {om0_parity_break:.4e}")
        parity_preserved = (h_parity_break < 1e-10 and om0_parity_break < 1e-10)
        print(f"  Even-weight parity preserved: {parity_preserved}")

        # ── Orbit period measurement ─────────────────────────────────────────
        T_pi  = float(np.pi)
        T_2pi = 2.0 * float(np.pi)

        detected_period, td_pi, td_2pi, ts_sweep, td_sweep = \
            measure_period(jmax, om0, max_T=4*np.pi, n_pts=800)

        # Working period: use detected (should be pi due to parity)
        period = detected_period
        print(f"  td@pi  = {td_pi:.4e}  (threshold: machine eps)")
        print(f"  td@2pi = {td_2pi:.4e}")
        print(f"  Detected working period = {period:.6f}  (pi={T_pi:.6f}, 2pi={T_2pi:.6f})")

        # ── Orbit injectivity ────────────────────────────────────────────────
        ts_inj = np.linspace(0.1, period - 0.1, 40)
        td_inj   = np.array([lib.trace_dist(flow_state(jmax, om0, t), om0) for t in ts_inj])
        dmax_inj = np.array([lib.d_max(flow_state(jmax, om0, t), om0)      for t in ts_inj])

        min_td   = float(np.min(td_inj))
        min_dmax = float(np.min(dmax_inj))
        max_td   = float(np.max(td_inj))
        print(f"  Injectivity: min td={min_td:.4e}  min dmax={min_dmax:.4e}  max td={max_td:.6f}")

        # ── Trace-scale law: max(td_profile) * Z ────────────────────────────
        trace_scale = max_td * Z
        print(f"  Trace-scale law: max(td) * Z = {trace_scale:.6f}")

        # ── Interval recovery ────────────────────────────────────────────────
        # Ensure true_delta < period/2 so recovery is unambiguous
        true_delta = min(0.7 * period, period * 0.35)
        ta_start   = period * 0.2
        om_a = flow_state(jmax, om0, ta_start)
        om_b = flow_state(jmax, om_a, true_delta)

        rec = interval_recovery(om_a, om_b, jmax, period, n_grid=600)
        err_dmax  = rec["recovered_dmax"]  - true_delta
        err_trace = rec["recovered_trace"] - true_delta
        rel_err_dmax  = err_dmax  / true_delta
        rel_err_trace = err_trace / true_delta
        print(f"  Recovery (true_delta={true_delta:.4f}):")
        print(f"    d_max:  recovered={rec['recovered_dmax']:.6f},  err={err_dmax:+.4e}  "
              f"({rel_err_dmax*100:+.3f}%)")
        print(f"    trace:  recovered={rec['recovered_trace']:.6f},  err={err_trace:+.4e}  "
              f"({rel_err_trace*100:+.3f}%)")

        # ── Additivity ───────────────────────────────────────────────────────
        add_details, max_add_err = additivity_check(
            om0, jmax, period, TRIPLE_SCALES, n_grid=600
        )
        print(f"  Additivity max |error| = {max_add_err:.4e}")
        for i, d in enumerate(add_details):
            print(f"    triple {i}: ell(a,b)={d['ell_ab']:.4f}, ell(b,c)={d['ell_bc']:.4f}, "
                  f"ell(a,c)={d['ell_ac']:.4f}, err={d['additivity_error_signed']:+.4e}")

        results[jkey] = {
            "jmax": float(jmax),
            "wedge_dim": dim,
            "Z": float(Z),
            "parity_structure": {
                "H_null_odd_even_norm": h_parity_break,
                "om0_odd_even_norm":    om0_parity_break,
                "even_weight_parity_preserved": bool(parity_preserved),
                "consequence": "effective_period_is_pi_not_2pi" if parity_preserved
                               else "effective_period_is_2pi",
            },
            "orbit_period": {
                "T_pi":   T_pi,
                "T_2pi":  T_2pi,
                "td_at_pi":           td_pi,
                "td_at_2pi":          td_2pi,
                "detected_period":    detected_period,
                "working_period":     float(period),
            },
            "injectivity": {
                "min_trace_dist": min_td,
                "min_dmax":       min_dmax,
                "max_trace_dist": max_td,
                "inj_ts":         ts_inj.tolist(),
                "inj_td_profile": td_inj.tolist(),
                "inj_dmax_profile": dmax_inj.tolist(),
            },
            "trace_scale_law": {
                "max_td_profile": max_td,
                "Z":              float(Z),
                "max_td_times_Z": float(trace_scale),
            },
            "interval_recovery": {
                "true_delta":       float(true_delta),
                "ta_start":         float(ta_start),
                "recovered_dmax":   rec["recovered_dmax"],
                "recovered_trace":  rec["recovered_trace"],
                "cost_dmax":        rec["cost_dmax"],
                "cost_trace":       rec["cost_trace"],
                "error_dmax":        float(err_dmax),
                "error_trace":       float(err_trace),
                "rel_error_dmax":    float(rel_err_dmax),
                "rel_error_trace":   float(rel_err_trace),
            },
            "additivity": {
                "triples": add_details,
                "max_abs_error": float(max_add_err),
            },
        }

    return results


if __name__ == "__main__":
    import io
    out_path = Path(__file__).resolve().parent / "data" / "wh7_band_exh_intervals.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    res = run()

    # Summary table
    print()
    print("=" * 80)
    hdr = f"{'jmax':>5}  {'dim':>5}  {'td@pi':>10}  {'td@2pi':>10}  "
    hdr += f"{'period':>8}  {'rec_err_dm%':>12}  {'rec_err_td%':>12}  {'add_err':>10}"
    print(hdr)
    print("-" * 80)
    for jkey, d in res.items():
        print(
            f"{d['jmax']:>5.1f}  "
            f"{d['wedge_dim']:>5}  "
            f"{d['orbit_period']['td_at_pi']:>10.3e}  "
            f"{d['orbit_period']['td_at_2pi']:>10.3e}  "
            f"{d['orbit_period']['working_period']:>8.5f}  "
            f"{d['interval_recovery']['rel_error_dmax']*100:>+12.4f}  "
            f"{d['interval_recovery']['rel_error_trace']*100:>+12.4f}  "
            f"{d['additivity']['max_abs_error']:>10.4e}"
        )

    print()
    print("Trace-scale law (max_td * Z):")
    for jkey, d in res.items():
        print(f"  jmax={d['jmax']:.1f}: {d['trace_scale_law']['max_td_times_Z']:.6f}")

    with open(out_path, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nJSON written to {out_path}")
