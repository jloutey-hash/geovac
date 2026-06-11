"""
Bridge Test 2 (final, v2): richer modular ring + per-tranche cleanliness control.

v1 (bridge_test2_final.py) gave: machinery VALID, spectral-only NULL (reproduces
Track 5), targets NULL against a MINIMAL modular extension {q=e^-2pi, e^-pi,
logGamma14, log_pi}.

v2 strengthens the MODULAR side to make the NULL conclusive:
  - Adds independent q-POWERS q^2=e^{-4pi}, q^3=e^{-6pi} (genuine q-series terms,
    Q-linearly independent of q, q2=e^-pi).
  - Adds a "modular-entropy-shaped" positive control PC5 = log(3) + 5*q - 2*q^2
    (this is EXACTLY the structural form bridge_test1 found for the wedge KMS
    entropy: S = log(integer) + q-series). If PSLQ recovers PC5, the basis can
    represent a genuine modular entropy -- so a NULL on the targets is a real
    statement that the atomic entropy is NOT of that form.
  - Keeps the basis Q-LINEARLY CLEAN: q, e^-pi, q^2=e^-4pi, q^3=e^-6pi are all
    independent (no q^2 = (q)^2 LINEAR relation -- squaring is nonlinear, so
    they ARE linearly independent over Q). logGamma14, log_pi, log2, log3
    independent. We verify cleanliness explicitly.

We do NOT add Gamma(1/4), varpi, eta(i) RAW alongside their logs (the v2
pathology). We do NOT add log_varpi or log_eta (Q-dependent on logGamma14,
log2, log_pi -- they would break cleanliness). The lemniscatic ring is fully
represented by the independent generator logGamma14 plus the spectral logs.

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

OUT_JSON = "debug/data/bridge_test2_final.json"  # overwrite v1 with the stronger run


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
    G = catalan        # = beta(2)  (M3)
    # beta(4) (M3, vertex-parity Dirichlet-L)
    beta4 = mp.nsum(lambda k: (-1) ** k / (2 * k + 1) ** 4, [0, mp.inf])
    s2 = sqrt(2)
    s3 = sqrt(3)
    s5 = sqrt(5)

    # SPECTRAL group (Track-5 family, cleaned & widened):
    spectral_names = ["1", "pi^2", "pi^4", "pi^6", "zeta3", "zeta5",
                      "log2", "log3", "log5", "Catalan", "beta4",
                      "sqrt2", "sqrt3", "sqrt5"]
    spectral_vals = [one, pi2, pi4, pi6, z3, z5, l2, l3, l5, G, beta4, s2, s3, s5]

    # MODULAR group: tau=i / q=e^{-2pi} ring.
    # q-powers (linearly independent over Q -- squaring is nonlinear):
    q1 = e**(-2 * pi)   # q
    q2p = e**(-pi)      # e^{-pi}  (= q^{1/2}, an independent transcendental)
    q4 = e**(-4 * pi)   # q^2
    q6 = e**(-6 * pi)   # q^3
    # lemniscatic / eta ring -- only the INDEPENDENT generator:
    LG = log(gamma(mpf(1) / 4))   # log Gamma(1/4)
    log_pi = log(pi)              # log pi

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


def fmt_relation(rel, names_full):
    terms = []
    for c, nm in zip(rel, names_full):
        if c != 0:
            terms.append(f"{int(c)}*{nm}")
    return " + ".join(terms)


def main():
    mp.dps = 320
    TOL_EXP = 90

    results = {
        "metadata": {
            "sprint": "Bridge Test 2 (final, v2): atomic correlation entropy modular?",
            "question": "Is S_full(GS) in a modular / q-series ring (vs the spectral ring Track 5 used)?",
            "context": ("Track 5 PSLQ'd S against a SPECTRAL mechanical basis (pi/zeta/log/Catalan) "
                        "and got NULL. This re-test ADDS the tau=i / q=e^{-2pi} modular ring "
                        "(q-powers + lemniscatic log-Gamma(1/4)) on a Q-linearly-clean basis, "
                        "validated by structurally-different positive controls including a "
                        "modular-entropy-shaped control PC5 = log3 + 5q - 2q^2 (the form "
                        "bridge_test1 found for the wedge KMS entropy)."),
            "dps": mp.dps,
            "tol_exp": TOL_EXP,
            "hit_residual_threshold": f"1e-{TOL_EXP}",
            "protocol": "W3: freeze basis, validate machinery (controls), THEN read targets",
            "maxsteps": 400000,
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

    # ---- Linear-cleanliness self-check (no internal small-int relation) --
    clean_rel = None
    try:
        clean_rel = pslq(full_vals, maxcoeff=10**6, maxsteps=400000,
                         tol=mpf(10) ** (-TOL_EXP))
    except Exception:
        clean_rel = None
    results["basis_linear_cleanliness"] = {
        "internal_relation_found": clean_rel is not None,
        "relation": [int(c) for c in clean_rel] if clean_rel else None,
        "relation_str": fmt_relation(clean_rel, full_names) if clean_rel else None,
        "verdict": "CLEAN" if clean_rel is None else "DEGENERATE -- internal dependency!",
    }

    # ---- Positive controls (>=5, structurally different) ----------------
    pc_specs = [
        ("PC1_spectral_only", {"pi^2": 3, "zeta3": -2, "log2": 5}),
        ("PC2_modular_only",  {"q_e^-2pi": 7, "logGamma14": -3, "log_pi": 2}),
        ("PC3_mixed",         {"zeta5": -1, "Catalan": 4, "q_e^-pi": -6, "log_pi": 1}),
        ("PC4_cross_clean",   {"1": 2, "pi^4": -1, "logGamma14": 3, "sqrt3": -5}),
        # PC5: modular-entropy-shaped (the bridge_test1 form: log(int) + q-series)
        ("PC5_modular_entropy_shape", {"log3": 1, "q_e^-2pi": 5, "q2_e^-4pi": -2}),
        # PC6: deep q-series (touches q^3)
        ("PC6_deep_qseries", {"q_e^-2pi": 1, "q2_e^-4pi": -3, "q3_e^-6pi": 4, "beta4": -1}),
    ]
    pc_results = {}
    machinery_positive_ok = True
    for cname, coeffs in pc_specs:
        cval = mpf(0)
        for nm, c in coeffs.items():
            cval += mpf(c) * name_to_val[nm]
        rel, resid = run_pslq(cval, full_vals, maxcoeff=10**4, tol_exp=TOL_EXP)
        hit = rel is not None and resid is not None and resid < mpf(10) ** (-TOL_EXP)
        pc_results[cname] = {
            "intended_coeffs": coeffs,
            "hit": bool(hit),
            "residual": mp.nstr(resid, 8) if resid is not None else None,
            "relation_str": fmt_relation(rel, ["CONTROL"] + full_names) if rel else None,
        }
        if not hit:
            machinery_positive_ok = False
    results["positive_controls"] = pc_results
    results["positive_controls_all_hit"] = machinery_positive_ok

    # ---- Null control ---------------------------------------------------
    null_val = mpf("0.3047295018763428193745018273645019283746501928374650192837465019")
    nrel, nresid = run_pslq(null_val, full_vals, maxcoeff=10**4, tol_exp=TOL_EXP)
    null_hit = nrel is not None and nresid is not None and nresid < mpf(10) ** (-TOL_EXP)
    results["null_control"] = {
        "value": mp.nstr(null_val, 20),
        "hit": bool(null_hit),
        "residual": mp.nstr(nresid, 8) if nresid is not None else None,
        "verdict": "FALSE POSITIVE (machinery broken)" if null_hit else "NULL (correct)",
    }

    machinery_valid = machinery_positive_ok and (not null_hit) and (clean_rel is None)
    results["machinery_valid"] = bool(machinery_valid)

    # ---- Spectral-only target check (reproduce Track 5 NULL) ------------
    spectral_only = {}
    for tname, tstr in TARGET_STRINGS.items():
        tval = mpf(tstr)
        rel, resid = run_pslq(tval, sp_vals, maxcoeff=10**3, tol_exp=TOL_EXP)
        hit = rel is not None and resid is not None and resid < mpf(10) ** (-TOL_EXP)
        spectral_only[tname] = {"hit": bool(hit),
                                "residual": mp.nstr(resid, 8) if resid is not None else None}
    results["spectral_only_target_check"] = spectral_only

    # ---- Targets vs FULL basis (maxcoeff sweep) -------------------------
    maxcoeff_sweep = [10**3, 10**4, 10**6]
    target_results = {}
    for tname, tstr in TARGET_STRINGS.items():
        tval = mpf(tstr)
        per_mc = {}
        any_hit = False
        for mc in maxcoeff_sweep:
            rel, resid = run_pslq(tval, full_vals, maxcoeff=mc, tol_exp=TOL_EXP)
            hit = rel is not None and resid is not None and resid < mpf(10) ** (-TOL_EXP)
            max_abs_coeff = max(abs(int(x)) for x in rel) if rel else None
            per_mc[f"maxcoeff_{mc}"] = {
                "hit": bool(hit),
                "residual": mp.nstr(resid, 8) if resid is not None else None,
                "relation_str": fmt_relation(rel, ["TARGET"] + full_names) if rel else None,
                "max_abs_coeff": max_abs_coeff,
            }
            if hit:
                any_hit = True
        target_results[tname] = {"any_hit": any_hit, "by_maxcoeff": per_mc}
    results["target_results"] = target_results

    # ---- HIT-survival (precision up + maxcoeff down) --------------------
    survival = {}
    for tname, tr in target_results.items():
        loose_hit = tr["by_maxcoeff"].get("maxcoeff_1000000", {}).get("hit", False) or \
                    tr["by_maxcoeff"].get("maxcoeff_10000", {}).get("hit", False)
        if loose_hit:
            mp.dps = 420
            tval = mpf(TARGET_STRINGS[tname])
            sp_n2, sp_v2, mod_n2, mod_v2 = build_basis()
            full_v2 = sp_v2 + mod_v2
            rel, resid = run_pslq(tval, full_v2, maxcoeff=10**3, tol_exp=130)
            hit = rel is not None and resid is not None and resid < mpf(10) ** (-130)
            survival[tname] = {"tested": True,
                               "survives_precision_up_maxcoeff_down": bool(hit),
                               "residual": mp.nstr(resid, 8) if resid is not None else None}
            mp.dps = 320
        else:
            survival[tname] = {"tested": False, "reason": "no loose hit to re-test"}
    results["hit_survival_check"] = survival

    # ---- Cross-ratios (extra: is any ratio modular?) --------------------
    ratio_results = {}
    He = mpf(TARGET_STRINGS["He_n3"])
    Li = mpf(TARGET_STRINGS["Li+_n3"])
    He4 = mpf(TARGET_STRINGS["He_n4"])
    ratios = {"He_n3/Li+_n3": He / Li, "He_n3/He_n4": He / He4, "He_n4/Li+_n3": He4 / Li}
    for rname, rval in ratios.items():
        rel, resid = run_pslq(rval, full_vals, maxcoeff=10**4, tol_exp=TOL_EXP)
        hit = rel is not None and resid is not None and resid < mpf(10) ** (-TOL_EXP)
        ratio_results[rname] = {"value": mp.nstr(rval, 30), "hit": bool(hit),
                                "residual": mp.nstr(resid, 8) if resid is not None else None,
                                "relation_str": fmt_relation(rel, ["RATIO"] + full_names) if rel else None}
    results["cross_ratio_check"] = ratio_results

    # ---- Verdict --------------------------------------------------------
    targets_any_hit = any(tr["any_hit"] for tr in target_results.values())
    targets_survive = any(s.get("survives_precision_up_maxcoeff_down", False)
                          for s in survival.values())
    if not machinery_valid:
        verdict = "INCONCLUSIVE -- instrument failed (controls did not validate machinery)"
        headline = "inconclusive-instrument-failed"
    elif targets_any_hit and targets_survive:
        verdict = ("TARGETS HIT and survive precision-up/maxcoeff-down -> correlation "
                   "entropy IS modular -> Track 5 tested wrong dictionary, headline lives")
        headline = "headline-lives"
    elif targets_any_hit and not targets_survive:
        verdict = ("TARGETS HIT at loose conditions but DID NOT survive tightening -> "
                   "selection-bias noise -> Track 5 hardened")
        headline = "headline-dead-Track-5-hardened (loose hits were noise)"
    else:
        verdict = ("TARGETS NULL against spectral+modular basis with VALID machinery -> "
                   "correlation entropy genuinely not in this ring (spectral OR modular) -> "
                   "Track 5 hardened, the bridge stops at correlation entropy")
        headline = "headline-dead-Track-5-hardened"
    results["verdict"] = verdict
    results["headline"] = headline

    os.makedirs("debug/data", exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump(results, fh, indent=2)

    print("=" * 72)
    print("BRIDGE TEST 2 (final, v2): atomic correlation entropy modular?")
    print("=" * 72)
    print(f"dps={results['metadata']['dps']}  tol=1e-{TOL_EXP}  n_full_basis={len(full_names)}")
    print(f"  spectral({len(sp_names)}): {sp_names}")
    print(f"  modular ({len(mod_names)}): {mod_names}")
    print(f"basis cleanliness: {results['basis_linear_cleanliness']['verdict']}")
    print("-" * 72)
    print("POSITIVE CONTROLS (must all HIT):")
    for cname, r in pc_results.items():
        print(f"  {cname:28s}: {'HIT' if r['hit'] else 'NULL':5s}  resid={r['residual']}")
    print(f"  -> all_hit = {machinery_positive_ok}")
    print("-" * 72)
    print(f"NULL CONTROL: {results['null_control']['verdict']}")
    print(f"MACHINERY VALID = {machinery_valid}")
    print("-" * 72)
    print("SPECTRAL-ONLY (expect NULL, reproduces Track 5):")
    for tname, r in spectral_only.items():
        print(f"  {tname:8s}: {'HIT' if r['hit'] else 'NULL'}")
    print("-" * 72)
    print("TARGETS vs FULL (spectral+modular):")
    for tname, tr in target_results.items():
        print(f"  {tname:8s}: any_hit={tr['any_hit']}")
        for mck, mcr in tr["by_maxcoeff"].items():
            print(f"      {mck:18s}: {'HIT' if mcr['hit'] else 'NULL':5s} "
                  f"resid={mcr['residual']} maxcoeff_in_rel={mcr['max_abs_coeff']}")
    print("-" * 72)
    print("CROSS-RATIOS vs FULL:")
    for rname, r in ratio_results.items():
        print(f"  {rname:18s}: {'HIT' if r['hit'] else 'NULL'}")
    print("-" * 72)
    print("HIT-survival:")
    for tname, s in survival.items():
        print(f"  {tname:8s}: {s}")
    print("=" * 72)
    print("VERDICT:", verdict)
    print("HEADLINE:", headline)
    print("=" * 72)
    print(f"Written: {OUT_JSON}")


if __name__ == "__main__":
    main()
