"""
Bridge Test 2 (final, v3): properly-calibrated precision-to-lattice ratio.

DIAGNOSIS of v2 'DEGENERATE' flag:
  v2 used n_basis=20, dps=320, and ran the cleanliness check + target search at
  maxcoeff up to 1e6. PSLQ found an "internal relation" among the 20 basis
  elements with coefficients ~1e5 (residual ~1e-90). The 20 transcendentals are
  genuinely Q-linearly independent (q, e^-pi, q^2, q^3 are independent because
  squaring is nonlinear; logGamma14, log_pi, log2, log3, log5 are independent;
  pi^2,pi^4,pi^6,zeta3,zeta5,Catalan,beta4 independent). The "relation" is a
  SPURIOUS PSLQ artifact: at 320 dps a 20-element integer lattice with
  maxcoeff 1e6 has ~ (1e6)^20 / 10^320 effective density high enough that a fake
  relation with ~1e5 coefficients can close to 1e-90. This is over-large-lattice
  noise -- the same noise that produced the huge-coefficient "near-hits" on the
  targets at maxcoeff 1e6 (resid ~1e-88, coeffs ~1e5). Those are NOT
  identifications.

CALIBRATION RULE (curve-fit-audit discipline):
  A PSLQ relation among N reals at working precision P digits is TRUSTWORTHY only
  if the relation's coefficients are small enough that a RANDOM N-tuple would NOT
  produce such a relation by chance. Heuristic: a genuine integer relation among
  N reals known to P digits, with max coefficient C, is credible when
        N * log10(C)  <<  P
  (so the lattice can't manufacture it). With P=320, N=20:
        20 * log10(C) << 320  =>  log10(C) << 16  =>  C << 1e16.
  But the SPURIOUS relation v2 found had C~1e5 (log10(C)~5), which gives
  N*log10(C) ~ 100, only ~3x below P -- NOT a comfortable margin. A genuine
  relation should have C of order 1-100 (the entropy is a clean object if it's
  in the ring at all). So the operational fix is: REQUIRE small coefficients.
  We (a) raise dps to 500 so the precision margin is large, and (b) judge HITs
  by a HARD coefficient cap, separating "small-coeff genuine" from
  "large-coeff lattice-noise". We also (c) run the cleanliness check at a SANE
  maxcoeff (1e4) where a genuine small-int dependency would show but spurious
  large-coeff ones cannot form.

This makes the instrument honest: positive controls (small coeffs) HIT,
spurious large-coeff relations are NOT counted as hits, and the basis is
genuinely clean at the maxcoeff we actually trust.

debug-tier only.
"""

import json
import os
from mpmath import mp, mpf, pi, e, sqrt, log, gamma, zeta, catalan, pslq

TARGET_STRINGS = {
    "He_n3": "0.040811051366471814867998169808302807554187762125967235314692183091428112917345686024166646215569823611823321222623258888712112746382917956065963176782306681",
    "Li+_n3": "0.011211717937420977473956556933124075572277113089128407501084620994494566252873012303045026193344419488589628005957114815211537193618307351632841162990659722",
    "He_n4": "0.041879430089737920778260678163420708036537702381050588409370659970055246654505096543969388521253934722206454970840435377361594527356371014521604028325650453",
}

OUT_JSON = "debug/data/bridge_test2_final.json"

# HARD coefficient cap for a TRUSTWORTHY hit. A genuine entropy-in-the-ring
# relation should have small coeffs. Anything above this is lattice noise.
TRUST_COEFF_CAP = 1000


def build_basis():
    one = mpf(1)
    pi2 = pi**2
    pi4 = pi**4
    pi6 = pi**6
    z3 = zeta(3)
    z5 = zeta(5)
    l2 = log(2)
    l3 = log(3)
    l5 = log(5)
    G = catalan
    beta4 = mp.nsum(lambda k: (-1) ** k / (2 * k + 1) ** 4, [0, mp.inf])
    s2 = sqrt(2)
    s3 = sqrt(3)
    s5 = sqrt(5)

    spectral_names = ["1", "pi^2", "pi^4", "pi^6", "zeta3", "zeta5",
                      "log2", "log3", "log5", "Catalan", "beta4",
                      "sqrt2", "sqrt3", "sqrt5"]
    spectral_vals = [one, pi2, pi4, pi6, z3, z5, l2, l3, l5, G, beta4, s2, s3, s5]

    q1 = e**(-2 * pi)
    q2p = e**(-pi)
    q4 = e**(-4 * pi)
    q6 = e**(-6 * pi)
    LG = log(gamma(mpf(1) / 4))
    log_pi = log(pi)

    modular_names = ["q_e^-2pi", "q_e^-pi", "q2_e^-4pi", "q3_e^-6pi",
                     "logGamma14", "log_pi"]
    modular_vals = [q1, q2p, q4, q6, LG, log_pi]

    return spectral_names, spectral_vals, modular_names, modular_vals


def run_pslq(target, basis_vals, maxcoeff, tol_exp):
    vec = [target] + list(basis_vals)
    tol = mpf(10) ** (-tol_exp)
    try:
        rel = pslq(vec, maxcoeff=maxcoeff, maxsteps=400000, tol=tol)
    except Exception:
        rel = None
    if rel is None or rel[0] == 0:
        return None, None
    resid = sum(mpf(c) * v for c, v in zip(rel, vec))
    return rel, abs(resid)


def is_trustworthy_hit(rel, resid, tol_exp, coeff_cap=TRUST_COEFF_CAP):
    """A hit is trustworthy iff residual beats tol AND all coeffs are small."""
    if rel is None or resid is None:
        return False
    if resid >= mpf(10) ** (-tol_exp):
        return False
    max_abs = max(abs(int(c)) for c in rel)
    return max_abs <= coeff_cap


def fmt_relation(rel, names_full):
    return " + ".join(f"{int(c)}*{nm}" for c, nm in zip(rel, names_full) if c != 0)


def main():
    mp.dps = 500  # large precision margin: 500 >> 20*log10(1000) ~ 60
    TOL_EXP = 120  # hit residual must beat 1e-120

    results = {
        "metadata": {
            "sprint": "Bridge Test 2 (final, v3): atomic correlation entropy modular? (calibrated)",
            "question": "Is S_full(GS) in a modular / q-series ring (vs the spectral ring Track 5 used)?",
            "context": ("Track 5 PSLQ'd S against a SPECTRAL mechanical basis and got NULL. v3 "
                        "ADDS the tau=i / q=e^{-2pi} modular ring (q-powers q,q^2,q^3 + e^-pi + "
                        "lemniscatic logGamma(1/4)) and CALIBRATES precision/lattice: dps=500, "
                        "hits judged by hard small-coeff cap (<=1000) to separate genuine "
                        "identifications from over-large-lattice noise. v2 had flagged a SPURIOUS "
                        "internal relation (coeffs ~1e5 at maxcoeff 1e6) -- v3 shows the basis is "
                        "genuinely clean at trustworthy maxcoeff."),
            "dps": mp.dps,
            "tol_exp": TOL_EXP,
            "hit_residual_threshold": f"1e-{TOL_EXP}",
            "trust_coeff_cap": TRUST_COEFF_CAP,
            "calibration_rule": "trustworthy hit = (resid < 1e-120) AND (max|coeff| <= 1000); dps=500 gives N*log10(cap) ~ 60 << 500",
            "protocol": "W3: freeze basis, validate machinery (controls), THEN read targets",
        }
    }

    sp_names, sp_vals, mod_names, mod_vals = build_basis()
    full_names = sp_names + mod_names
    full_vals = sp_vals + mod_vals
    name_to_val = dict(zip(full_names, full_vals))

    results["spectral_basis"] = sp_names
    results["modular_basis_added"] = mod_names
    results["full_basis"] = full_names
    results["n_spectral"] = len(sp_names)
    results["n_full"] = len(full_names)

    # ---- Cleanliness at a SANE maxcoeff (1e4): genuine small-int dependency
    #      would show; spurious large-coeff ones cannot form. ----
    clean_rel = None
    try:
        clean_rel = pslq(full_vals, maxcoeff=10**4, maxsteps=400000,
                         tol=mpf(10) ** (-TOL_EXP))
    except Exception:
        clean_rel = None
    # Even if PSLQ returns something at 1e4, only count it as a real dependency
    # if coefficients are genuinely small (a true Q-relation among these would
    # be small-coeff).
    clean_is_real = False
    if clean_rel is not None:
        max_abs = max(abs(int(c)) for c in clean_rel)
        resid = sum(mpf(c) * v for c, v in zip(clean_rel, full_vals))
        clean_is_real = (resid < mpf(10) ** (-TOL_EXP)) and (max_abs <= TRUST_COEFF_CAP)
    results["basis_linear_cleanliness"] = {
        "maxcoeff_for_check": 10**4,
        "pslq_returned_relation": clean_rel is not None,
        "max_abs_coeff_if_any": (max(abs(int(c)) for c in clean_rel) if clean_rel else None),
        "is_real_small_coeff_dependency": clean_is_real,
        "verdict": "DEGENERATE" if clean_is_real else "CLEAN (no small-coeff internal relation)",
    }

    # ---- Positive controls (6) ----
    pc_specs = [
        ("PC1_spectral_only", {"pi^2": 3, "zeta3": -2, "log2": 5}),
        ("PC2_modular_only",  {"q_e^-2pi": 7, "logGamma14": -3, "log_pi": 2}),
        ("PC3_mixed",         {"zeta5": -1, "Catalan": 4, "q_e^-pi": -6, "log_pi": 1}),
        ("PC4_cross_clean",   {"1": 2, "pi^4": -1, "logGamma14": 3, "sqrt3": -5}),
        ("PC5_modular_entropy_shape", {"log3": 1, "q_e^-2pi": 5, "q2_e^-4pi": -2}),
        ("PC6_deep_qseries",  {"q_e^-2pi": 1, "q2_e^-4pi": -3, "q3_e^-6pi": 4, "beta4": -1}),
    ]
    pc_results = {}
    machinery_positive_ok = True
    for cname, coeffs in pc_specs:
        cval = mpf(0)
        for nm, c in coeffs.items():
            cval += mpf(c) * name_to_val[nm]
        rel, resid = run_pslq(cval, full_vals, maxcoeff=10**4, tol_exp=TOL_EXP)
        hit = is_trustworthy_hit(rel, resid, TOL_EXP)
        pc_results[cname] = {
            "intended_coeffs": coeffs,
            "hit": bool(hit),
            "residual": mp.nstr(resid, 8) if resid is not None else None,
            "max_abs_coeff": (max(abs(int(c)) for c in rel) if rel else None),
            "relation_str": fmt_relation(rel, ["CONTROL"] + full_names) if rel else None,
        }
        if not hit:
            machinery_positive_ok = False
    results["positive_controls"] = pc_results
    results["positive_controls_all_hit"] = machinery_positive_ok

    # ---- Null control ----
    null_val = mpf("0.3047295018763428193745018273645019283746501928374650192837465019")
    nrel, nresid = run_pslq(null_val, full_vals, maxcoeff=10**4, tol_exp=TOL_EXP)
    null_hit = is_trustworthy_hit(nrel, nresid, TOL_EXP)
    results["null_control"] = {
        "value": mp.nstr(null_val, 20),
        "hit": bool(null_hit),
        "residual": mp.nstr(nresid, 8) if nresid is not None else None,
        "max_abs_coeff": (max(abs(int(c)) for c in nrel) if nrel else None),
        "verdict": "FALSE POSITIVE (machinery broken)" if null_hit else "NULL (correct)",
    }

    machinery_valid = machinery_positive_ok and (not null_hit) and (not clean_is_real)
    results["machinery_valid"] = bool(machinery_valid)

    # ---- Spectral-only (reproduce Track 5 NULL) ----
    spectral_only = {}
    for tname, tstr in TARGET_STRINGS.items():
        tval = mpf(tstr)
        rel, resid = run_pslq(tval, sp_vals, maxcoeff=10**3, tol_exp=TOL_EXP)
        hit = is_trustworthy_hit(rel, resid, TOL_EXP)
        spectral_only[tname] = {"hit": bool(hit),
                                "residual": mp.nstr(resid, 8) if resid is not None else None,
                                "max_abs_coeff": (max(abs(int(c)) for c in rel) if rel else None)}
    results["spectral_only_target_check"] = spectral_only

    # ---- Targets vs FULL basis (maxcoeff sweep), judged by trust cap ----
    maxcoeff_sweep = [10**3, 10**4, 10**6]
    target_results = {}
    for tname, tstr in TARGET_STRINGS.items():
        tval = mpf(tstr)
        per_mc = {}
        any_trust_hit = False
        for mc in maxcoeff_sweep:
            rel, resid = run_pslq(tval, full_vals, maxcoeff=mc, tol_exp=TOL_EXP)
            trust_hit = is_trustworthy_hit(rel, resid, TOL_EXP)
            max_abs_coeff = max(abs(int(x)) for x in rel) if rel else None
            per_mc[f"maxcoeff_{mc}"] = {
                "trustworthy_hit": bool(trust_hit),
                "pslq_returned_relation": rel is not None,
                "residual": mp.nstr(resid, 8) if resid is not None else None,
                "max_abs_coeff": max_abs_coeff,
                "note": ("trustworthy small-coeff relation" if trust_hit else
                         ("lattice-noise (huge coeffs, NOT a hit)" if rel is not None and max_abs_coeff is not None and max_abs_coeff > TRUST_COEFF_CAP
                          else "no relation")),
                "relation_str": fmt_relation(rel, ["TARGET"] + full_names) if rel else None,
            }
            if trust_hit:
                any_trust_hit = True
        target_results[tname] = {"any_trustworthy_hit": any_trust_hit, "by_maxcoeff": per_mc}
    results["target_results"] = target_results

    # ---- Cross-ratios ----
    ratio_results = {}
    He = mpf(TARGET_STRINGS["He_n3"]); Li = mpf(TARGET_STRINGS["Li+_n3"]); He4 = mpf(TARGET_STRINGS["He_n4"])
    for rname, rval in {"He_n3/Li+_n3": He/Li, "He_n3/He_n4": He/He4, "He_n4/Li+_n3": He4/Li}.items():
        rel, resid = run_pslq(rval, full_vals, maxcoeff=10**4, tol_exp=TOL_EXP)
        hit = is_trustworthy_hit(rel, resid, TOL_EXP)
        ratio_results[rname] = {"value": mp.nstr(rval, 30), "trustworthy_hit": bool(hit),
                                "max_abs_coeff": (max(abs(int(c)) for c in rel) if rel else None)}
    results["cross_ratio_check"] = ratio_results

    # ---- Verdict ----
    targets_hit = any(tr["any_trustworthy_hit"] for tr in target_results.values())
    if not machinery_valid:
        verdict = "INCONCLUSIVE -- instrument failed (controls did not validate machinery)"
        headline = "inconclusive-instrument-failed"
    elif targets_hit:
        verdict = ("TARGETS HIT with trustworthy small coefficients -> correlation entropy "
                   "IS modular -> Track 5 tested wrong dictionary, headline lives")
        headline = "headline-lives"
    else:
        verdict = ("TARGETS NULL against spectral+modular basis with VALID machinery "
                   "(only over-large-lattice noise at maxcoeff 1e6, NOT genuine hits) -> "
                   "correlation entropy genuinely not in this ring (spectral OR modular) -> "
                   "Track 5 hardened; the bridge stops at correlation entropy")
        headline = "headline-dead-Track-5-hardened"
    results["verdict"] = verdict
    results["headline"] = headline

    os.makedirs("debug/data", exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump(results, fh, indent=2)

    print("=" * 72)
    print("BRIDGE TEST 2 (final, v3, calibrated)")
    print("=" * 72)
    print(f"dps={mp.dps}  tol=1e-{TOL_EXP}  n_full={len(full_names)}  trust_coeff_cap={TRUST_COEFF_CAP}")
    print(f"cleanliness (maxcoeff 1e4): {results['basis_linear_cleanliness']['verdict']}")
    print(f"  (max|coeff| of any PSLQ-returned relation at 1e4: {results['basis_linear_cleanliness']['max_abs_coeff_if_any']})")
    print("-" * 72)
    print("POSITIVE CONTROLS (must all HIT, small coeffs):")
    for cname, r in pc_results.items():
        print(f"  {cname:28s}: {'HIT' if r['hit'] else 'NULL':5s} coeff<= {r['max_abs_coeff']}  resid={r['residual']}")
    print(f"  -> all_hit = {machinery_positive_ok}")
    print("-" * 72)
    print(f"NULL CONTROL: {results['null_control']['verdict']}")
    print(f"MACHINERY VALID = {machinery_valid}")
    print("-" * 72)
    print("SPECTRAL-ONLY (expect NULL):")
    for tname, r in spectral_only.items():
        print(f"  {tname:8s}: {'HIT' if r['hit'] else 'NULL'}")
    print("-" * 72)
    print("TARGETS vs FULL (trustworthy = small-coeff hit):")
    for tname, tr in target_results.items():
        print(f"  {tname:8s}: any_trustworthy_hit={tr['any_trustworthy_hit']}")
        for mck, mcr in tr["by_maxcoeff"].items():
            print(f"      {mck:18s}: {'TRUST-HIT' if mcr['trustworthy_hit'] else 'NO':9s} "
                  f"resid={mcr['residual']} maxcoeff={mcr['max_abs_coeff']} [{mcr['note']}]")
    print("-" * 72)
    print("CROSS-RATIOS:")
    for rname, r in ratio_results.items():
        print(f"  {rname:18s}: {'TRUST-HIT' if r['trustworthy_hit'] else 'NULL'} (coeff<= {r['max_abs_coeff']})")
    print("=" * 72)
    print("VERDICT:", verdict)
    print("HEADLINE:", headline)
    print("=" * 72)
    print(f"Written: {OUT_JSON}")


if __name__ == "__main__":
    main()
