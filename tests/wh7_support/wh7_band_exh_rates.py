#!/usr/bin/env python3
"""
WH7 B3 Phase-3 band-exhaustion -- RATE IDENTIFICATION (TAG=rates).

Reads the two probe JSONs:
    debug/data/wh7_band_exh_legs.json       (legs + beta probe)
    debug/data/wh7_band_exh_penalties.json  (per-class penalties)

and classifies the decay of successive differences d_k = x(j_{k+1}) - x(j_k)
against four candidate rate classes:
    (i)   gamma-type (Paper 38 action-seminorm):  d ~ log(2j)/(2j)
    (ii)  thermal:                                d ~ exp(-2*beta*j) = exp(-beta*(2j))
    (iii) power law:                              d ~ (2j)^{-p}      (p fitted)
    (iv)  1/Z (Gibbs flattening):                 d ~ 1/Z(window)
plus two diagnostics outside the requested list (clearly labelled):
    (v)   geometric / free-rate exponential:      d ~ exp(-c*(2j))   (c fitted)
    (vi)  rational form 1/((2j+1)(2j+2))          (amplitude only; motivated by the
          EXACT closed form discovered for H_norm(null), see below)

GATE AUDIT: the probes flagged NO observable formally CONVERGENT.  The mixed-
reference observables (c12,c23,c13,deficit,H_norm) have strictly positive,
monotone-decreasing diffs (converging-but-ungated); the null-reference
observables are staircases whose integer-step subsequence has only 2 positive
diffs.  We classify the former, treat the latter separately (single-ratio
comparison only), and never clamp: machine-eps diffs are reported signed as-is
and treated as exact freezes, not as small positives.

Deterministic: no randomness anywhere.  numpy + scipy.special.zeta only.
"""
import json
import math
import os
from fractions import Fraction

import numpy as np
from scipy.special import zeta as hurwitz_zeta

HERE = os.path.dirname(os.path.abspath(__file__))            # .../debug
DATA = os.path.join(HERE, "data")
LEGS_PATH = os.path.join(DATA, "wh7_band_exh_legs.json")
PEN_PATH = os.path.join(DATA, "wh7_band_exh_penalties.json")
OUT_PATH = os.path.join(DATA, "wh7_band_exh_rates.json")

ZERO_TOL = 1e-12      # |diff| below this = exact freeze (machine epsilon)
FREEZE_TOL = 1e-10    # pair-freeze detection in penalty E-sequences
RMS_CONS = 0.05       # ratio-level rms log-residual below this -> CONSISTENT
RMS_EXCL = 0.20       # above this -> EXCLUDED
BETA = 1.0            # beta of the main ladder

JM = np.array([1.0, 1.5, 2.0, 2.5, 3.0])
TWOJ = 2.0 * JM                          # [2,3,4,5,6]
X_UP = TWOJ[1:]                          # diff k lives between windows k,k+1
X_LO = TWOJ[:-1]
X_MID = 0.5 * (X_UP + X_LO)


# ----------------------------------------------------------------- utilities
def gamma_model(x):
    return np.log(x) / x


def thermal_model(x):
    # d ~ exp(-2*beta*j) = exp(-beta * (2j))
    return np.exp(-BETA * x)


def amp_fit(d, m):
    """1-parameter (amplitude) fit in log space.  Returns (A, rms_log_resid)."""
    ln_d, ln_m = np.log(d), np.log(m)
    lnA = float(np.mean(ln_d - ln_m))
    res = ln_d - (lnA + ln_m)
    return math.exp(lnA), float(np.sqrt(np.mean(res ** 2)))


def power_fit(d, x):
    """2-parameter fit ln d = lnA - p ln x.  Returns (A, p, rms_log_resid)."""
    ln_d, ln_x = np.log(d), np.log(x)
    M = np.vstack([np.ones_like(ln_x), ln_x]).T
    coef, *_ = np.linalg.lstsq(M, ln_d, rcond=None)
    pred = M @ coef
    rms = float(np.sqrt(np.mean((ln_d - pred) ** 2)))
    return math.exp(coef[0]), float(-coef[1]), rms


def geo_fit(d, x):
    """2-parameter fit ln d = a - c x.  Returns (A, c, rms_log_resid)."""
    ln_d = np.log(d)
    M = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(M, ln_d, rcond=None)
    pred = M @ coef
    rms = float(np.sqrt(np.mean((ln_d - pred) ** 2)))
    return math.exp(coef[0]), float(-coef[1]), rms


def ratios(seq):
    a = np.asarray(seq, dtype=float)
    return (a[1:] / a[:-1]).tolist()


def rms_log_ratio(r_obs, r_pred):
    lr = np.log(np.asarray(r_obs) / np.asarray(r_pred))
    return float(np.sqrt(np.mean(lr ** 2)))


def monotone_flag(seq):
    a = np.asarray(seq, dtype=float)
    if np.all(np.diff(a) > 0):
        return "increasing"
    if np.all(np.diff(a) < 0):
        return "decreasing"
    return "non-monotone"


def nice_fraction(v, max_den=64, tol=1e-12):
    fr = Fraction(v).limit_denominator(max_den)
    if abs(float(fr) - v) < tol:
        return str(fr)
    return None


# ------------------------------------------------------- load the probe data
with open(LEGS_PATH, "r") as fh:
    legs = json.load(fh)
with open(PEN_PATH, "r") as fh:
    pen = json.load(fh)

Z_ladder = np.array([row["Z"] for row in legs["null"]["rows"]])  # same for both refs
assert np.allclose(Z_ladder, [row["Z"] for row in legs["mixed"]["rows"]])

out = {
    "meta": {
        "tag": "rates",
        "date": "2026-06-10",
        "inputs": ["wh7_band_exh_legs.json", "wh7_band_exh_penalties.json"],
        "jmax_ladder": ["1", "3/2", "2", "5/2", "3"],
        "two_j": TWOJ.tolist(),
        "beta_main": BETA,
        "zero_tol": ZERO_TOL,
        "rms_consistent": RMS_CONS,
        "rms_excluded": RMS_EXCL,
        "x_convention": "diffs evaluated at UPPER window 2j (sensitivity at lower/mid reported)",
        "free_parameter_count": {
            "gamma": 1, "thermal_beta_fixed": 1, "power": 2,
            "invZ_upper": 1, "invZ_lower": 1, "geometric_free_rate": 2,
            "rational_a1": 1,
        },
    },
    "gate_audit": {
        "formally_CONVERGENT_observables": [],
        "note": ("Legs probe verdicts: all UNDECIDED (null = staircase; mixed = "
                 "monotone but last diff/value 4-5% > 1% gate). Penalty probe: all 7 "
                 "classes NON-MONOTONE / MIXED. Rate classification below is therefore "
                 "performed on the CONVERGING-BUT-UNGATED set (mixed-reference c12/c23/"
                 "c13/deficit/H_norm: strictly positive monotone-decreasing diffs) and, "
                 "separately, on the null integer-step subsequence (2 diffs only). "
                 "This selection is post hoc -- see caveats."),
    },
}


# ------------------------------------------- H_norm(null) exact-form discovery
hn_null = np.array([row["H_norm"] for row in legs["null"]["rows"]])
d_hn_null = np.diff(hn_null)
prod_check = d_hn_null * (X_UP + 1.0) * (X_UP + 2.0)        # should be constant
const_dev = float(np.max(np.abs(prod_check - 2.0 * math.sqrt(6.0))))
closed_form = math.sqrt(6.0) * TWOJ / (TWOJ + 2.0)
closed_dev = float(np.max(np.abs(hn_null - closed_form)))
ratio_fracs = [nice_fraction(r) for r in ratios(d_hn_null)]
out["h_norm_null_exact"] = {
    "values": hn_null.tolist(),
    "diffs": d_hn_null.tolist(),
    "diff_ratios": ratios(d_hn_null),
    "diff_ratios_as_fractions": ratio_fracs,                 # expect 2/3, 5/7, 3/4
    "closed_form": "H_norm(jmax) = sqrt(6) * (2j)/(2j+2);  d_k = 2*sqrt(6)/((2j+1)(2j+2))",
    "d_times_(2j+1)(2j+2)": prod_check.tolist(),
    "max_dev_from_2sqrt6": const_dev,
    "max_dev_values_vs_closed_form": closed_dev,
    "continuum_norm_limit": math.sqrt(6.0),
    "limit_identification": "sqrt(6) = sqrt(b(b+1)(2b+1)) at b=1 (algebraic, pi-free)",
    "fitted_p_over_window": power_fit(d_hn_null, X_UP)[1],
    "asymptotic_p_of_closed_form": 2.0,
    "lesson": ("A bit-exact 1/((2j+1)(2j+2)) sequence (asymptotic p=2) FITS as "
               "p_eff~1.5 over the available window 2j in [3,6].  Any p~1.5 claim "
               "below is therefore a LOCAL effective exponent, window-degenerate "
               "with p=2 rational forms."),
}

# H_norm(mixed): rational scan + same machinery (reported, not load-bearing)
hn_mix = np.array([row["H_norm"] for row in legs["mixed"]["rows"]])
hn_mix_sq_fracs = [nice_fraction(v * v, max_den=1000, tol=1e-9) for v in hn_mix]
d_hn_mix = np.diff(hn_mix)
amp_a05 = d_hn_mix * (X_UP + 0.5) * (X_UP + 1.5)
out["h_norm_mixed_check"] = {
    "values": hn_mix.tolist(),
    "values_sq_low_denominator_fractions": hn_mix_sq_fracs,
    "diffs": d_hn_mix.tolist(),
    "diff_ratios": ratios(d_hn_mix),
    "amplitude_under_1/((2j+1/2)(2j+3/2))": amp_a05.tolist(),
    "amplitude_scatter_pct": float(100 * np.std(amp_a05) / np.mean(amp_a05)),
    "note": ("Amplitude under the half-shifted rational form is ~constant (~1% "
             "scatter) but NOT machine-exact (contrast null: 1e-15).  Suggestive "
             "only; continuum norm of the mixed multiplier not identified."),
}


# --------------------------------------------------- mixed-ref classification
def classify(d, label, Z):
    """Full 4-diff classification against all candidate models."""
    d = np.asarray(d, dtype=float)
    assert np.all(d > ZERO_TOL), f"{label}: non-positive diff fed to classifier"
    r_obs = ratios(d)
    rec = {"diffs": d.tolist(), "ratios_observed": r_obs,
           "ratio_trend": monotone_flag(r_obs),
           "span_decades": float(np.log10(d.max() / d.min()))}

    models = {}
    # (i) gamma at upper 2j
    m = gamma_model(X_UP)
    A, rms_v = amp_fit(d, m)
    models["gamma"] = {"k_free": 1, "A": A, "rms_value": rms_v,
                       "ratios_pred": ratios(m), "rms_ratio": rms_log_ratio(r_obs, ratios(m))}
    # (ii) thermal, beta fixed = 1
    m = thermal_model(X_UP)
    A, rms_v = amp_fit(d, m)
    models["thermal_beta1"] = {"k_free": 1, "A": A, "rms_value": rms_v,
                               "ratios_pred": ratios(m), "rms_ratio": rms_log_ratio(r_obs, ratios(m))}
    # (iii) power law (p free); also p under lower/mid conventions
    A, p, rms_v = power_fit(d, X_UP)
    _, p_lo, _ = power_fit(d, X_LO)
    _, p_mid, _ = power_fit(d, X_MID)
    m = X_UP ** (-p)
    models["power"] = {"k_free": 2, "A": A, "p_upper": p, "p_lower": p_lo, "p_mid": p_mid,
                       "rms_value": rms_v, "ratios_pred": ratios(m),
                       "rms_ratio": rms_log_ratio(r_obs, ratios(m))}
    # (iv) 1/Z, both alignments
    for tag, zz in (("invZ_upper", Z[1:]), ("invZ_lower", Z[:-1])):
        m = 1.0 / zz
        A, rms_v = amp_fit(d, m)
        models[tag] = {"k_free": 1, "A": A, "rms_value": rms_v,
                       "ratios_pred": ratios(m), "ratio_trend_pred": monotone_flag(ratios(m)),
                       "rms_ratio": rms_log_ratio(r_obs, ratios(m))}
    # (v) geometric / free-rate exponential (diagnostic, not in the 4-list)
    A, c, rms_v = geo_fit(d, X_UP)
    m = np.exp(-c * X_UP)
    models["geometric_free_rate"] = {"k_free": 2, "A": A, "c": c, "rms_value": rms_v,
                                     "ratios_pred": ratios(m),
                                     "rms_ratio": rms_log_ratio(r_obs, ratios(m))}
    # (vi) rational a=1 form (diagnostic; motivated by H_norm(null) exact form)
    m = 1.0 / ((X_UP + 1.0) * (X_UP + 2.0))
    A, rms_v = amp_fit(d, m)
    models["rational_a1"] = {"k_free": 1, "A": A, "rms_value": rms_v,
                             "ratios_pred": ratios(m),
                             "rms_ratio": rms_log_ratio(r_obs, ratios(m))}

    # verdicts: ratio test is primary; demand consistency with value fit
    best = min(models, key=lambda k: models[k]["rms_ratio"])
    best_rms = models[best]["rms_ratio"]
    for name, mm in models.items():
        rr = mm["rms_ratio"]
        shape_pred = monotone_flag(mm["ratios_pred"])
        shape_ok = (shape_pred == rec["ratio_trend"])
        if rr > RMS_EXCL or (not shape_ok and rr > 3 * max(best_rms, 1e-12)):
            v = "EXCLUDED"
        elif rr < RMS_CONS and shape_ok:
            v = "CONSISTENT"
        else:
            v = "MARGINAL"
        mm["shape_pred"] = shape_pred
        mm["verdict"] = v
    # consistency between ratio-test winner and value-fit winner
    best_value = min(models, key=lambda k: models[k]["rms_value"])
    rec["models"] = models
    rec["best_by_ratio_test"] = best
    rec["best_by_value_fit"] = best_value
    rec["ratio_value_consistent"] = bool(
        best == best_value
        or models[best_value]["rms_ratio"] < 1.5 * best_rms)
    return rec


mixed_diffs = legs["mixed"]["diffs"]
mixed_cls = {}
for obs in ("c12", "c23", "c13", "deficit", "H_norm"):
    mixed_cls[obs] = classify(mixed_diffs[obs], f"mixed.{obs}", Z_ladder)
out["mixed_classification"] = mixed_cls

# tail extrapolation for the headline observables (clearly caveated)
extrap = {}
for obs in ("c12", "c13", "deficit"):
    d = np.asarray(mixed_diffs[obs])
    last_val = legs["mixed"]["rows"][-1][obs]
    A_pow, p, _ = power_fit(d, X_UP)
    tail_pow = float(A_pow * hurwitz_zeta(p, 7.0)) if p > 1 else None
    m = 1.0 / ((X_UP + 1.0) * (X_UP + 2.0))
    A_rat, _ = amp_fit(d, m)
    tail_rat = float(A_rat / 8.0)        # telescoping sum_{x>=7} 1/((x+1)(x+2)) = 1/8
    extrap[obs] = {
        "value_at_jmax3": last_val,
        "power_fit": {"p": p, "tail": tail_pow,
                      "limit": (last_val + tail_pow) if tail_pow else None},
        "rational_a1": {"A": A_rat, "tail": tail_rat, "limit": last_val + tail_rat},
        "caveat": "tail is 35-55% of current value; limit uncertain at the +/-15% level",
    }
out["mixed_limit_extrapolation"] = extrap


# ------------------------------------------------ null-ref staircase handling
null_stairs = {}
for obs in ("c12", "c23", "c13", "deficit"):
    d = np.asarray(legs["null"]["diffs"][obs], dtype=float)
    frozen = [i for i, v in enumerate(d) if abs(v) < ZERO_TOL]
    pos = [(i, v) for i, v in enumerate(d) if v > ZERO_TOL]
    rec = {"diffs_signed": d.tolist(),
           "frozen_steps_indices": frozen,
           "frozen_steps_note": "reported signed, machine-epsilon; treated as exact freezes, NOT clamped",
           "n_positive_diffs": len(pos)}
    if len(pos) == 2:
        (i1, v1), (i2, v2) = pos
        x1, x2 = X_UP[i1], X_UP[i2]          # 2j of the two jumps (4 and 6)
        r = v2 / v1
        preds = {
            "gamma": gamma_model(x2) / gamma_model(x1),
            "thermal_beta1": math.exp(-BETA * (x2 - x1)),
            "invZ": Z_ladder[i1 + 1] / Z_ladder[i2 + 1],
            "rational_a1": ((x1 + 1) * (x1 + 2)) / ((x2 + 1) * (x2 + 2)),
            "power_p1.5": (x2 / x1) ** (-1.5),
        }
        rec.update({
            "jump_values": [v1, v2], "jump_2j": [float(x1), float(x2)],
            "observed_ratio": r,
            "implied_power_p_star": math.log(r) / math.log(x1 / x2),
            "model_predicted_ratios": {k: float(v) for k, v in preds.items()},
            "model_log_misfit": {k: float(abs(math.log(r / v))) for k, v in preds.items()},
            "verdict": ("UNDECIDED -- 2 diffs / 1 ratio: zero dof.  gamma and thermal "
                        "excluded by margin; power(p~1.5), 1/Z and rational_a1 are "
                        "numerically DEGENERATE on this single step."),
        })
    null_stairs[obs] = rec
# the degeneracy worth recording explicitly
null_stairs["degeneracy_note"] = {
    "Z(2j=6)/Z(2j=4)": float(Z_ladder[4] / Z_ladder[2]),
    "(6/4)^1.5": float((6 / 4) ** 1.5),
    "rational_a1 30/56": 30 / 56,
    "comment": ("Z itself grows ~ (2j)^{3/2} (smoothed) at beta=1, so candidates "
                "(iii) p~1.5 and (iv) 1/Z are partially the SAME hypothesis on the "
                "null subsequence; for the mixed ref they separate because Z "
                "zig-zags within integer/half-integer pairs while the mixed diffs "
                "are smooth -- there 1/Z is excluded by shape."),
}
out["null_staircase"] = null_stairs

# null H_norm is smooth -> full classification applies (and the form is EXACT)
out["null_H_norm_classification"] = classify(legs["null"]["diffs"]["H_norm"],
                                             "null.H_norm", Z_ladder)


# -------------------------------------------------------- beta discriminator
bp = legs["beta_probe"]
b_rows = {0.5: bp["3"]["0.5"], 2.0: bp["3"]["2.0"]}
b_rows_lo = {0.5: bp["5/2"]["0.5"], 2.0: bp["5/2"]["2.0"]}
r1 = {k: legs["null"]["rows"][i] for k, i in (("5/2", 3), ("3", 4))}
beta_disc = {}
for obs in ("c12", "deficit"):
    d_beta = {}
    d_beta["0.5"] = b_rows[0.5][obs] - b_rows_lo[0.5][obs]
    d_beta["1.0"] = r1["3"][obs] - r1["5/2"][obs]
    d_beta["2.0"] = b_rows[2.0][obs] - b_rows_lo[2.0][obs]
    obs_ratio_05 = d_beta["0.5"] / d_beta["1.0"]
    obs_ratio_20 = d_beta["2.0"] / d_beta["1.0"]
    # thermal prediction for the AMPLITUDE of the 5/2->3 diff: exp(-2*beta*3) at j=3
    pred_05 = math.exp(-2 * 0.5 * 3) / math.exp(-2 * 1.0 * 3)   # e^{+3}
    pred_20 = math.exp(-2 * 2.0 * 3) / math.exp(-2 * 1.0 * 3)   # e^{-6}
    beta_disc[obs] = {
        "diff_5/2_to_3_by_beta": d_beta,
        "observed_ratio_beta0.5_over_beta1": obs_ratio_05,
        "observed_ratio_beta2_over_beta1": obs_ratio_20,
        "thermal_pred_ratio_beta0.5": pred_05,
        "thermal_pred_ratio_beta2": pred_20,
        "exclusion_factor_beta0.5": pred_05 / obs_ratio_05,
        "exclusion_factor_beta2": obs_ratio_20 / pred_20,
        "verdict": ("THERMAL EXCLUDED: diffs GROW with beta (super-linear), "
                    "opposite in direction and orders of magnitude away from "
                    "exp(-2*beta*j) suppression."),
    }
beta_disc["scope_caveat"] = ("beta probe exists only for the NULL reference and only at "
                             "jmax in {5/2, 3} (one diff per beta) -- this excludes a "
                             "thermal MECHANISM for the staircase steps directly; for the "
                             "mixed ref thermal is excluded independently by the ratio "
                             "test (pred 0.368 vs obs 0.67-0.74).")
out["beta_discriminator"] = beta_disc


# --------------------------------------------------- penalty-layer rate audit
pen_audit = {}
b_parity = {"integer_b_pair_freeze": [], "half_integer_b_every_step": [], "other": []}
for cname, cdata in pen["per_class"].items():
    b_val = float(cname.split(",")[0].strip("("))
    E = [w.get("E_0p2_signed") for w in cdata["windows"]]
    E_seq = [e for e in E if e is not None]
    d = list(np.diff(E_seq))
    frozen = [i for i, v in enumerate(d) if abs(v) < FREEZE_TOL]
    pos = [(i, v) for i, v in enumerate(d) if v > FREEZE_TOL]
    neg = [(i, v) for i, v in enumerate(d) if v < -FREEZE_TOL]
    # pair-freeze staircase = freezes exactly at the within-pair steps
    all_eps = all(abs(e) < 1e-10 for e in E_seq)
    staircase = bool(len(E_seq) == 5 and abs(d[0]) < FREEZE_TOL
                     and abs(d[2]) < FREEZE_TOL and not all_eps)
    rec = {"b": b_val, "E_sequence": E, "diffs_signed": d,
           "frozen_steps": frozen, "n_pos": len(pos), "n_neg": len(neg),
           "monotone": monotone_flag(E_seq) if len(E_seq) > 1 else "n/a",
           "pair_freeze_staircase": staircase,
           "rate_classification": "NOT APPLICABLE (non-monotone / mixed-sign diffs)"}
    if all_eps:
        rec["rate_classification"] = "EXACT ZERO at every window (commuting class, frozen)"
    elif staircase and len(pos) == 2 and len(neg) == 0:
        v1, v2 = pos[0][1], pos[1][1]
        rec["integer_step_jumps"] = [v1, v2]
        rec["integer_step_ratio"] = v2 / v1
        rec["implied_power_p_star"] = math.log(v2 / v1) / math.log(4.0 / 6.0)
        rec["rate_classification"] = ("staircase: integer-step subsequence, 1 ratio, "
                                      "UNDECIDED (consistent with the legs-probe "
                                      "power/1/Z degenerate class)")
    elif staircase:
        rec["rate_classification"] = ("staircase BUT sign-mixed integer steps "
                                      "(growth then decrease or sign-crossing): "
                                      "NOT APPLICABLE")
    # b-parity bookkeeping
    if all_eps:
        b_parity["other"].append(cname + " (exact zero)")
    elif staircase:
        b_parity["integer_b_pair_freeze"].append(cname)
    elif b_val % 1 == 0.5:
        b_parity["half_integer_b_every_step"].append(cname)
    else:
        b_parity["other"].append(cname)
    pen_audit[cname] = rec

# special case: (2.0,1.0) folded to zero at window 1 -> staircase test on 4 windows
c21 = pen["per_class"]["(2.0,1.0) spacelike"]
E21 = [w.get("E_0p2_signed") for w in c21["windows"]][1:]
d21 = list(np.diff(E21))
pen_audit["(2.0,1.0) spacelike"]["post_zero_window_diffs"] = d21
pen_audit["(2.0,1.0) spacelike"]["note"] = (
    "first window folded to zero; over windows 2..5 the class pair-freezes at "
    "(2,5/2) and its integer-step jumps GROW (0.0553 -> 0.2805, ratio 5.07): "
    "divergent staircase, opposite of the null-ref legs staircase")

out["penalty_rate_audit"] = {
    "per_class": pen_audit,
    "b_parity_dichotomy": b_parity,
    "structural_note": (
        "Pair-freeze staircases (movement ONLY at half-integer->integer window "
        "steps) occur exactly for integer-b classes; half-integer-b classes "
        "((0.5,0.5), (1.5,1.5)) move at every step.  Same dichotomy as the legs "
        "probe (null ref b=1 staircases, mixed ref b=1/2 smooth).  Matches the "
        "(-1)^{2b} spin-statistics grading of B3 Phase-1: a (-1)^{2b}-even "
        "multiplier cannot couple the newly added half-integer shell to the "
        "integer shells.  OBSERVATION, not proven here."),
}


# ------------------------------------------------------------------ headline
p_c12 = mixed_cls["c12"]["models"]["power"]["p_upper"]
p_c13 = mixed_cls["c13"]["models"]["power"]["p_upper"]
p_def = mixed_cls["deficit"]["models"]["power"]["p_upper"]
out["headline"] = {
    "question": ("does the state-level layer converge at the metric-layer gamma "
                 "rate, the thermal rate, or something else?"),
    "answer": (
        "SOMETHING ELSE.  The state-level (D_max) diffs follow a power-law class "
        f"with local effective exponent p_eff ~ 1.5 (c12 {p_c12:.2f}, c13 {p_c13:.2f}, "
        f"deficit {p_def:.2f}).  Gamma-type log(2j)/(2j) is EXCLUDED (predicted "
        "per-step diff ratios 0.93-0.95 vs observed 0.67-0.74).  Thermal "
        "exp(-2*beta*j) is EXCLUDED twice over (fixed-beta ratio e^{-1}=0.368 vs "
        "observed 0.67-0.74; beta probe shows diffs GROW with beta, ~3 OoM away "
        "and wrong direction).  1/Z is EXCLUDED for the mixed ref by shape "
        "(zig-zag vs smooth) but DEGENERATE with power p~1.5 on the null "
        "integer-step subsequence (Z ~ (2j)^{3/2} smoothed at beta=1).  "
        "Constant identification is NOT claimed: the H_norm(null) exact closed "
        "form proves that an asymptotic-p=2 rational sequence fits as p_eff~1.49 "
        "over this window.  If the diffs are increments of a convergent sequence, "
        "the implied state-level ERROR rate is (2j)^{-(p_eff-1)} ~ (2j)^{-0.5}: "
        "SLOWER than the metric layer's gamma rate log(2j)/(2j)."),
}

out["caveats"] = [
    "GATE: no observable was formally CONVERGENT; classification performed on the "
    "converging-but-ungated mixed-reference set (post-hoc selection on apparent "
    "monotone convergence) -- selection-bias flagged per audit rule.",
    "Only 4 diffs per observable spanning <0.5 decades: this is a CLASS "
    "classification, not a constant identification.  p_eff ~ 1.5 is a local "
    "exponent; the bit-exact H_norm(null) counterexample (asymptotic p=2, fitted "
    "p_eff 1.49 on this window) shows the window cannot distinguish p=1.5 from "
    "rational p->2 forms.",
    "Free parameters: power and geometric carry 2 each vs 1 for gamma/thermal/1/Z; "
    "part of power's residual advantage is structural.  Its win rests on the "
    "ratio-TREND signature (observed ratios strictly increasing 0.67->0.74), which "
    "no 1-parameter candidate and no constant-ratio model reproduces.",
    "power vs geometric(free-rate): residual discrimination is weak (rms_ratio "
    "~0.03 vs ~0.04); the increasing-ratio trend favors power qualitatively.  "
    "A beta probe on the MIXED reference would discriminate further (geometric-"
    "thermal mechanisms must shift with beta).",
    "rational_a1 model was tried AFTER the H_norm(null) exact form was discovered "
    "(one effective look-elsewhere); it is amplitude-only (1 param) and fits the "
    "mixed ratios at rms ~0.015, better than the 2-param power fit -- reported as "
    "an observation, not promoted.",
    "Limit extrapolations carry 35-55% tails; quoted limits uncertain at +/-15%.",
    "beta discriminator data exist only for the null reference (one diff per "
    "beta); thermal exclusion for the mixed ref rests on the ratio test alone.",
    "Penalty layer: all 7 classes fail the monotonicity precondition; no rate "
    "class is assigned there.  The b-parity staircase dichotomy is an "
    "observation consistent with the (-1)^{2b} grading, not a theorem.",
]

with open(OUT_PATH, "w") as fh:
    json.dump(out, fh, indent=2)

# ------------------------------------------------------------- console digest
print("=== WH7 band-exhaustion RATES (TAG=rates) ===")
print(f"gate audit: formally CONVERGENT observables: NONE")
print(f"\nH_norm(null) EXACT: H = sqrt(6)*(2j)/(2j+2); "
      f"max dev {closed_dev:.2e}; d*(2j+1)(2j+2)-2sqrt6 max {const_dev:.2e}")
print(f"  fitted p over window: {out['h_norm_null_exact']['fitted_p_over_window']:.3f} "
      f"(true asymptotic p = 2)  <-- window-degeneracy lesson")
print("\nmixed-ref classification (ratio-test rms | verdict):")
for obs in ("c12", "c13", "deficit", "H_norm"):
    mm = mixed_cls[obs]["models"]
    line = f"  {obs:8s} p_eff={mm['power']['p_upper']:.3f}  "
    line += " ".join(f"{k}:{mm[k]['rms_ratio']:.3f}/{mm[k]['verdict'][:4]}"
                     for k in ("gamma", "thermal_beta1", "power", "invZ_upper",
                               "geometric_free_rate", "rational_a1"))
    print(line)
    print(f"           obs ratios {['%.4f' % r for r in mixed_cls[obs]['ratios_observed']]} "
          f"trend={mixed_cls[obs]['ratio_trend']} "
          f"consistent(ratio~value)={mixed_cls[obs]['ratio_value_consistent']}")
print("\nnull staircase (integer-step subsequence, single ratio):")
for obs in ("c12", "c13", "deficit"):
    r = null_stairs[obs]
    print(f"  {obs:8s} ratio={r['observed_ratio']:.4f} p*={r['implied_power_p_star']:.3f} "
          f"preds: " + " ".join(f"{k}={v:.4f}" for k, v in r["model_predicted_ratios"].items()))
print("\nbeta discriminator (null ref, 5/2->3 diff):")
for obs in ("c12", "deficit"):
    bd = beta_disc[obs]
    print(f"  {obs}: obs r(b=2)/r(b=1)={bd['observed_ratio_beta2_over_beta1']:.3f} "
          f"thermal pred {bd['thermal_pred_ratio_beta2']:.2e} "
          f"(off x{bd['exclusion_factor_beta2']:.0f}, wrong direction)")
print("\npenalty layer: all classes NON-MONOTONE -> no rate assigned")
print("  b-parity dichotomy:", json.dumps(b_parity, indent=2))
print(f"\nJSON written: {OUT_PATH}")
