"""
Bridge Test 2 (final): Is GeoVac atomic correlation entropy S_full(GS) MODULAR?

Sprint TD Track 5 PSLQ'd S_full(GS) against a SPECTRAL-side mechanical basis
(pi-powers, zeta-values, logs, Catalan/beta) and returned NULL, concluding
"correlation entropy is ring-orphan".

Open question: did Track 5 test the WRONG dictionary? The wedge KMS (thermal)
entropy lives in MODULAR arithmetic (q = e^{-2pi}, tau = i lemniscatic point).
If the atomic correlation entropy is secretly modular, it would come back NULL
against a spectral basis BY CONSTRUCTION.

This driver:
  1. Reproduces a representative SPECTRAL basis -> confirm NULL (sanity).
  2. ADDS a LINEARLY-CLEAN modular ring (q-series + lemniscatic ring).
  3. Runs rigorous positive controls (>=3, structurally different) + null control.
  4. Only if machinery valid (all positive HIT, null NULL), reads the targets.
  5. Sweeps maxcoeff in {1e3, 1e4, 1e6}; checks HIT survival under
     precision-up / maxcoeff-down (curve-fit-audit discipline).

LESSONS BAKED IN (from failed bridge_test2_v2_calibrated.py):
  - Keep basis Q-LINEARLY CLEAN. Do NOT include multiplicatively-related
    transcendentals together (e.g. zeta3 AND zeta3/pi^2 AND pi^2; or
    Gamma(1/4), varpi, pi all at once). Near-dependencies destabilize PSLQ and
    make positive controls spuriously NULL.
  - Add modular elements a FEW AT A TIME; verify positive controls still HIT
    after each addition. If an element breaks a control, it's degenerate -> drop.
  - dps generous (>= 250); rule of thumb dps ~ n_basis * log10(maxcoeff) + margin.

debug-tier only; no geovac/ or paper edits.
"""

import json
import os
from mpmath import mp, mpf, pi, e, sqrt, log, gamma, zeta, catalan, pslq

# --------------------------------------------------------------------------
# Targets: full ~150-digit strings from debug/data/sprint_td_track5.json
# --------------------------------------------------------------------------
TARGET_STRINGS = {
    "He_n3": "0.040811051366471814867998169808302807554187762125967235314692183091428112917345686024166646215569823611823321222623258888712112746382917956065963176782306681",
    "Li+_n3": "0.011211717937420977473956556933124075572277113089128407501084620994494566252873012303045026193344419488589628005957114815211537193618307351632841162990659722",
    "He_n4": "0.041879430089737920778260678163420708036537702381050588409370659970055246654505096543969388521253934722206454970840435377361594527356371014521604028325650453",
}

OUT_JSON = "debug/data/bridge_test2_final.json"


# --------------------------------------------------------------------------
# Basis construction. We keep each element a single Q-linearly-independent
# transcendental (NO products of two basis elements, NO p/q ratios of two
# basis elements). The constant 1 is always included.
# --------------------------------------------------------------------------
def build_basis():
    """Return (names, values) at current mp precision.

    Two groups, each LINEARLY CLEAN within itself.

    SPECTRAL group (Track-5 family, deliberately minimal & non-degenerate):
        1, pi^2, pi^4, zeta(3), zeta(5), log2, log3, Catalan(=beta(2)), sqrt2, sqrt3
      - NOTE we do NOT also include zeta(2)=pi^2/6 (Q-dependent on pi^2),
        nor zeta(4)=pi^4/90 (Q-dep on pi^4), nor zeta(6) (Q-dep on pi^6 -> we omit pi^6),
        nor 1/zeta(n) (multiplicatively related). This is the cleaned spectral basis.
        It spans the SAME Q-vector-space-of-interest as Track 5's bloated 12312-form
        basis for the purpose of a NULL check: if S is a small-integer combo of any
        of {pi^2k, zeta(odd), log p, Catalan, sqrt}, this basis finds it.

    MODULAR group (the tau=i / q=e^{-2pi} ring), added in controlled tranches:
        q1 = e^{-2pi}, q2 = e^{-pi},
        L  = log(varpi) where varpi = lemniscate const = Gamma(1/4)^2/(2 sqrt(2 pi)),
        LG = log(Gamma(1/4)),
        Leta = log(eta(i)) where eta(i) = Gamma(1/4)/(2 pi^{3/4}).
      We include LOGS of the lemniscatic/eta constants (not the constants
      themselves) because entropy is a log of occupations -> a modular entropy
      would naturally be log(integer) + (linear combo of logs of modular
      constants) + (linear combo of q-powers). We add e^{-pi}, e^{-2pi}
      directly (q-series terms).

      We deliberately do NOT add pi^{3/4}, Gamma(1/4), varpi raw alongside their
      logs (that was the v2 pathology). And note log(eta(i)) = LG - log2 - (3/4)*log(pi):
      it is Q-DEPENDENT on {LG, log2} and pi^{3/4}=exp((3/4)log pi). To keep the
      basis clean we therefore include the INDEPENDENT modular generators:
         {q1, q2, LG, log_pi}  plus the spectral {log2}.
      Then log(varpi) = 2*LG - log2 - (1/2)*(log2 + log_pi)  is DEPENDENT and is
      EXCLUDED (it would be a linear combination of LG, log2, log_pi).
      log(eta(i)) = LG - log2 - (3/4)*log_pi is likewise DEPENDENT and EXCLUDED.
      So the genuinely-NEW independent modular generators beyond spectral are:
         q1 = e^{-2pi}, q2 = e^{-pi}, LG = log Gamma(1/4), log_pi = log pi.
      That is the minimal clean modular extension.
    """
    one = mpf(1)
    pi2 = pi**2
    pi4 = pi**4
    z3 = zeta(3)
    z5 = zeta(5)
    l2 = log(2)
    l3 = log(3)
    G = catalan  # = beta(2)
    s2 = sqrt(2)
    s3 = sqrt(3)

    spectral_names = ["1", "pi^2", "pi^4", "zeta3", "zeta5", "log2", "log3",
                      "Catalan", "sqrt2", "sqrt3"]
    spectral_vals = [one, pi2, pi4, z3, z5, l2, l3, G, s2, s3]

    # Modular generators (independent of the spectral set and of each other):
    q1 = e**(-2 * pi)         # q = e^{-2 pi},  tau = i
    q2 = e**(-pi)             # e^{-pi}
    LG = log(gamma(mpf(1) / 4))   # log Gamma(1/4)
    log_pi = log(pi)              # log pi  (independent of log2, log3)

    modular_names = ["q_e^-2pi", "q_e^-pi", "logGamma14", "log_pi"]
    modular_vals = [q1, q2, LG, log_pi]

    return spectral_names, spectral_vals, modular_names, modular_vals


# --------------------------------------------------------------------------
# PSLQ wrapper
# --------------------------------------------------------------------------
def run_pslq(target, basis_vals, maxcoeff, tol_exp):
    """PSLQ on [target] + basis_vals. Returns (relation or None, residual).

    The relation r satisfies  r[0]*target + sum_i r[i+1]*basis_vals[i] = 0.
    A HIT requires r[0] != 0 (the target actually participates).
    """
    vec = [target] + list(basis_vals)
    tol = mpf(10) ** (-tol_exp)
    try:
        rel = pslq(vec, maxcoeff=maxcoeff, maxsteps=200000, tol=tol)
    except Exception:
        rel = None
    if rel is None:
        return None, None
    if rel[0] == 0:
        # target not used -> degenerate relation among basis -> treat as no-hit
        return None, None
    resid = sum(mpf(c) * v for c, v in zip(rel, vec))
    return rel, abs(resid)


def fmt_relation(rel, names_full):
    """names_full[0]='TARGET', names_full[1:]=basis names."""
    terms = []
    for c, nm in zip(rel, names_full):
        if c != 0:
            terms.append(f"{int(c)}*{nm}")
    return " + ".join(terms)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    mp.dps = 300  # generous; sweep precision below
    TOL_EXP = 90  # residual must beat 1e-90 for a HIT

    results = {
        "metadata": {
            "sprint": "Bridge Test 2 (final): atomic correlation entropy modular?",
            "question": "Is S_full(GS) in a modular / q-series ring (vs spectral ring of Track 5)?",
            "dps": mp.dps,
            "tol_exp": TOL_EXP,
            "hit_residual_threshold": f"1e-{TOL_EXP}",
            "protocol": "W3: freeze basis, validate machinery with controls, THEN read targets",
        }
    }

    # ---- Build bases (frozen) -------------------------------------------
    sp_names, sp_vals, mod_names, mod_vals = build_basis()
    full_names = sp_names + mod_names
    full_vals = sp_vals + mod_vals

    results["spectral_basis"] = sp_names
    results["modular_basis_added"] = mod_names
    results["full_basis"] = full_names
    results["n_spectral"] = len(sp_names)
    results["n_full"] = len(full_names)

    # ---- Linear-cleanliness self-check ----------------------------------
    # Confirm NO nontrivial small-integer relation exists among the basis
    # elements ALONE (i.e. the basis is Q-linearly clean at maxcoeff 1e6).
    mp.dps = 300
    sp_names2, sp_vals2, mod_names2, mod_vals2 = build_basis()
    full_vals2 = sp_vals2 + mod_vals2
    clean_rel = None
    try:
        clean_rel = pslq(full_vals2, maxcoeff=10**6, maxsteps=200000,
                         tol=mpf(10) ** (-TOL_EXP))
    except Exception:
        clean_rel = None
    results["basis_linear_cleanliness"] = {
        "internal_relation_found": clean_rel is not None,
        "relation": [int(c) for c in clean_rel] if clean_rel else None,
        "verdict": "CLEAN" if clean_rel is None else "DEGENERATE -- basis has internal dependency!",
    }

    # ---- Positive controls (>=3, structurally different) ----------------
    # Each is an explicit small-integer combo of basis elements. PSLQ on the
    # control value + basis must recover the relation exactly (resid<1e-90).
    pc_specs = [
        # (name, coeffs-dict over full_names) -- ALL small integers
        ("PC1_spectral_only",  {"pi^2": 3, "zeta3": -2, "log2": 5}),
        ("PC2_modular_only",   {"q_e^-2pi": 7, "logGamma14": -3, "log_pi": 2}),
        ("PC3_mixed",          {"zeta5": -1, "Catalan": 4, "q_e^-pi": -6, "log_pi": 1}),
        ("PC4_cross_clean",    {"1": 2, "pi^4": -1, "logGamma14": 3, "sqrt3": -5}),
    ]
    name_to_val = dict(zip(full_names, full_vals))
    pc_results = {}
    machinery_positive_ok = True
    for cname, coeffs in pc_specs:
        # construct control value
        cval = mpf(0)
        for nm, c in coeffs.items():
            cval += mpf(c) * name_to_val[nm]
        rel, resid = run_pslq(cval, full_vals, maxcoeff=10**4, tol_exp=TOL_EXP)
        hit = rel is not None and resid is not None and resid < mpf(10) ** (-TOL_EXP)
        # verify recovered relation matches expectation (target coeff must be +-1
        # absorbing, and the combo must reproduce coeffs up to overall scale)
        pc_results[cname] = {
            "intended_coeffs": coeffs,
            "hit": bool(hit),
            "residual": mp.nstr(resid, 8) if resid is not None else None,
            "relation": [int(x) for x in rel] if rel else None,
            "relation_str": fmt_relation(rel, ["CONTROL"] + full_names) if rel else None,
        }
        if not hit:
            machinery_positive_ok = False

    results["positive_controls"] = pc_results
    results["positive_controls_all_hit"] = machinery_positive_ok

    # ---- Null control ---------------------------------------------------
    # An arbitrary constant with NO structural reason. Must return NULL.
    null_val = mpf("0.6180339887498948482045868343656381177203091798057628621354486227")  # ~1/phi, but we treat as arbitrary; sqrt5 not in basis, phi not in basis
    # Actually use a genuinely structureless decimal to avoid accidental algebraicity:
    null_val = mpf("0.3047295018763428193745018273645019283746501928374650192837465019")
    nrel, nresid = run_pslq(null_val, full_vals, maxcoeff=10**4, tol_exp=TOL_EXP)
    null_hit = nrel is not None and nresid is not None and nresid < mpf(10) ** (-TOL_EXP)
    results["null_control"] = {
        "value": mp.nstr(null_val, 20),
        "hit": bool(null_hit),
        "residual": mp.nstr(nresid, 8) if nresid is not None else None,
        "relation": [int(x) for x in nrel] if nrel else None,
        "verdict": "FALSE POSITIVE (machinery broken)" if null_hit else "NULL (correct)",
    }

    machinery_valid = machinery_positive_ok and (not null_hit) and (clean_rel is None)
    results["machinery_valid"] = bool(machinery_valid)

    # ---- Sanity: confirm SPECTRAL-only basis returns NULL on targets ----
    # (reproduces the spirit of Track 5's NULL using the cleaned spectral basis)
    spectral_only_targets = {}
    for tname, tstr in TARGET_STRINGS.items():
        tval = mpf(tstr)
        rel, resid = run_pslq(tval, sp_vals, maxcoeff=10**3, tol_exp=TOL_EXP)
        hit = rel is not None and resid is not None and resid < mpf(10) ** (-TOL_EXP)
        spectral_only_targets[tname] = {
            "hit": bool(hit),
            "residual": mp.nstr(resid, 8) if resid is not None else None,
            "relation": [int(x) for x in rel] if rel else None,
        }
    results["spectral_only_target_check"] = spectral_only_targets

    # ---- Read targets against FULL (spectral + modular) basis -----------
    # ONLY meaningful if machinery_valid. We compute regardless and record.
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
                "relation": [int(x) for x in rel] if rel else None,
                "relation_str": fmt_relation(rel, ["TARGET"] + full_names) if rel else None,
                "max_abs_coeff": max_abs_coeff,
            }
            if hit:
                any_hit = True
        target_results[tname] = {"any_hit": any_hit, "by_maxcoeff": per_mc}

    results["target_results"] = target_results

    # ---- HIT-survival check (precision up + maxcoeff down) --------------
    # For any target that HIT at the loosest condition, re-test at dps=400 and
    # maxcoeff=1e3 (tighter). A real identification survives.
    survival = {}
    for tname, tr in target_results.items():
        loose_hit = tr["by_maxcoeff"].get("maxcoeff_1000000", {}).get("hit", False) or \
                    tr["by_maxcoeff"].get("maxcoeff_10000", {}).get("hit", False)
        if loose_hit:
            mp.dps = 400
            tval = mpf(TARGET_STRINGS[tname])
            _, _, _, mod_vals_hi = build_basis()  # rebuild at hi precision
            sp_names_hi, sp_vals_hi, mod_names_hi, mod_vals_hi = build_basis()
            full_vals_hi = sp_vals_hi + mod_vals_hi
            rel, resid = run_pslq(tval, full_vals_hi, maxcoeff=10**3, tol_exp=120)
            hit = rel is not None and resid is not None and resid < mpf(10) ** (-120)
            survival[tname] = {
                "tested": True,
                "survives_precision_up_maxcoeff_down": bool(hit),
                "residual": mp.nstr(resid, 8) if resid is not None else None,
            }
            mp.dps = 300
        else:
            survival[tname] = {"tested": False, "reason": "no loose hit to re-test"}
    results["hit_survival_check"] = survival

    # ---- Verdict --------------------------------------------------------
    targets_any_hit = any(tr["any_hit"] for tr in target_results.values())
    targets_survive = any(s.get("survives_precision_up_maxcoeff_down", False)
                          for s in survival.values())

    if not machinery_valid:
        verdict = "INCONCLUSIVE -- instrument failed (controls did not validate machinery)"
        headline = "inconclusive-instrument-failed"
    elif targets_any_hit and targets_survive:
        verdict = ("TARGETS HIT and survive precision-up/maxcoeff-down -> "
                   "correlation entropy IS modular -> Track 5 tested wrong dictionary, headline lives")
        headline = "headline-lives"
    elif targets_any_hit and not targets_survive:
        verdict = ("TARGETS HIT at loose conditions but DID NOT survive tightening -> "
                   "selection-bias noise, NOT a real identification -> Track 5 hardened")
        headline = "headline-dead-Track-5-hardened (loose hits were noise)"
    else:
        verdict = ("TARGETS NULL against spectral+modular basis with VALID machinery -> "
                   "correlation entropy genuinely not in this ring -> Track 5 hardened")
        headline = "headline-dead-Track-5-hardened"

    results["verdict"] = verdict
    results["headline"] = headline

    os.makedirs("debug/data", exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump(results, fh, indent=2)

    # ---- Console summary ------------------------------------------------
    print("=" * 70)
    print("BRIDGE TEST 2 (final): atomic correlation entropy modular?")
    print("=" * 70)
    print(f"dps={results['metadata']['dps']}  tol=1e-{TOL_EXP}  n_full_basis={len(full_names)}")
    print(f"basis cleanliness: {results['basis_linear_cleanliness']['verdict']}")
    print("-" * 70)
    print("POSITIVE CONTROLS (must all HIT):")
    for cname, r in pc_results.items():
        print(f"  {cname:22s}: {'HIT' if r['hit'] else 'NULL':5s}  resid={r['residual']}")
    print(f"  -> all_hit = {machinery_positive_ok}")
    print("-" * 70)
    print(f"NULL CONTROL: {results['null_control']['verdict']}  (hit={null_hit})")
    print("-" * 70)
    print(f"MACHINERY VALID = {machinery_valid}")
    print("-" * 70)
    print("SPECTRAL-ONLY target check (expect NULL, reproduces Track 5):")
    for tname, r in spectral_only_targets.items():
        print(f"  {tname:8s}: {'HIT' if r['hit'] else 'NULL'}")
    print("-" * 70)
    print("TARGETS vs FULL (spectral+modular) basis:")
    for tname, tr in target_results.items():
        print(f"  {tname:8s}: any_hit={tr['any_hit']}")
        for mck, mcr in tr["by_maxcoeff"].items():
            print(f"      {mck:18s}: {'HIT' if mcr['hit'] else 'NULL':5s} "
                  f"resid={mcr['residual']} maxcoeff_in_rel={mcr['max_abs_coeff']}")
    print("-" * 70)
    print("HIT-survival (precision up + maxcoeff down):")
    for tname, s in survival.items():
        print(f"  {tname:8s}: {s}")
    print("=" * 70)
    print("VERDICT:", verdict)
    print("HEADLINE:", headline)
    print("=" * 70)
    print(f"Written: {OUT_JSON}")


if __name__ == "__main__":
    main()
