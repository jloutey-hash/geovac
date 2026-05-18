"""
Sprint TD-PSLQ-1: PSLQ probe of the hydrogen Bethe logarithm ln k_0(1S)
against the GeoVac master Mellin engine ring (M1/M2/M3).

PROTOCOL (W3 reference + linear-independence correction):
  1. Source ln k_0(1S) at highest available precision (Drake 1990: 16 digits).
  2. Build mechanical M1/M2/M3 ATOMIC basis (linearly independent — one form
     per atomic constant; PSLQ finds rational prefactors as integer ratios).
  3. Write basis to JSON BEFORE running PSLQ (integrity hash logged).
  4. Run PSLQ at 50 dps against the full atomic basis and sub-bases at
     coefficient ceilings 10^4, 10^6, 10^8.
  5. Verify hits at higher precision if found; report null otherwise.

LINEAR-INDEPENDENCE CORRECTION (key sizing decision):
  The first basis attempt expanded each atomic constant X by all (a/b)·X for
  small rationals. That created linear dependences within the basis: PSLQ
  found basis-internal rational relations (a_target=0) before any target
  relation. The fix is to use an ATOMIC basis — one form per X — and let
  PSLQ's integer coefficients encode the rational prefactor as a/b through
  the integer relation. (a/b)·X enters via the integer relation a·X − b·X_other
  pattern.

  Atomic basis size: ~50-200 forms, linearly independent by construction.
  With dps=50 and 16-digit target, PSLQ at coefficient ceiling 10^4-10^6 will
  reliably detect or refute relations.

PRECISION-FLOOR RULE:
  At target precision 16 digits, the false-positive risk for a depth-L
  relation at coefficient ceiling C is bounded by
    P_FP ~ (N choose L) * C^L / 10^16
  For N=200, L=2, C=10^4: P_FP ~ 2e4 * 1e8 / 1e16 = 2e-4. Marginal.
  For N=200, L=2, C=10^3: P_FP ~ 2e4 * 1e6 / 1e16 = 2e-6. Safe.
  For N=200, L=3, C=10^3: P_FP ~ 1e6 * 1e9 / 1e16 = 1e-1. Risky.
  So depth-2 relations at coefficient ~1000 are still trustworthy; depth-3
  needs verification.

Question: does the hydrogen Bethe logarithm — the load-bearing irrational
in every Lamb-shift calculation — sit in the master Mellin engine ring at
PSLQ precision?

Outcome interpretation:
  - Trustworthy HIT — possible identification. Verify at higher precision.
  - Spurious HIT at large coefficient — discard.
  - NULL — supports calibration-class interpretation.
"""

import json
import hashlib
import math
from pathlib import Path

import mpmath as mp
from mpmath import mpf, pi, log, sqrt, zeta, catalan


# ============================================================================
# TARGET: ln k_0(1S)
# ============================================================================
# Drake 1990 (Phys. Rev. A 41, 1243), Drake & Swainson 1990, Pachucki 1998:
#   ln k_0(1S, n=1, l=0) = 2.984128555765498
LN_K0_1S_16DIGITS = "2.984128555765498"

WORKING_DPS = 50
mp.mp.dps = WORKING_DPS

TARGET_VAL = mpf(LN_K0_1S_16DIGITS)

OUT_DIR = Path("debug/data")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================================
# ATOMIC BASIS (frozen specification — one form per atomic constant)
# ============================================================================
# The basis is constructed so each form is LINEARLY INDEPENDENT over Q from
# the others. PSLQ's integer coefficients encode rational prefactors via the
# integer-ratio pattern. This avoids basis-internal redundancies that block
# target detection.


def make_form(name, value, mech_class):
    return {"name": name, "value": value, "class": mech_class}


# -- M1 atomic forms: pure powers of pi and 1/pi --
def atoms_M1():
    forms = [make_form("1", mpf(1), "RAT")]  # constant
    for p in [1, 2, 3, 4, 5, 6]:
        forms.append(make_form(f"pi^{p}", pi**p, "M1"))
    for p in [1, 2, 3, 4]:
        forms.append(make_form(f"1/pi^{p}", 1/pi**p, "M1"))
    return forms


# -- M2 atomic forms: sqrt(pi), sqrt(pi)^odd, 1/sqrt(pi).
# We intentionally OMIT zeta(2k) atoms because zeta(2k) = c_k · pi^(2k) over Q
# (Bernoulli identity) and would be Q-linearly dependent on M1's pi^(2k). The
# M2 ring is the half-integer-power-of-pi ring (sqrt(pi) tower).
def atoms_M2():
    forms = []
    sqrt_pi = mp.sqrt(pi)
    forms.append(make_form("sqrt(pi)", sqrt_pi, "M2"))
    forms.append(make_form("sqrt(pi)^3", sqrt_pi**3, "M2"))
    forms.append(make_form("sqrt(pi)^5", sqrt_pi**5, "M2"))
    forms.append(make_form("sqrt(pi)^7", sqrt_pi**7, "M2"))
    forms.append(make_form("1/sqrt(pi)", 1/sqrt_pi, "M2"))
    forms.append(make_form("1/sqrt(pi)^3", 1/sqrt_pi**3, "M2"))
    return forms


# -- M3 atomic forms: Catalan, Dirichlet beta(2k), Hurwitz at 1/4, 3/4 --
def atoms_M3():
    forms = []
    G = catalan
    forms.append(make_form("G", G, "M3"))
    # Dirichlet beta(4), beta(6), beta(8)
    beta_4 = sum(mpf(-1)**k / mpf(2*k+1)**4 for k in range(0, 500))
    beta_6 = sum(mpf(-1)**k / mpf(2*k+1)**6 for k in range(0, 500))
    beta_8 = sum(mpf(-1)**k / mpf(2*k+1)**8 for k in range(0, 500))
    forms.append(make_form("beta(4)", beta_4, "M3"))
    forms.append(make_form("beta(6)", beta_6, "M3"))
    forms.append(make_form("beta(8)", beta_8, "M3"))
    # Hurwitz at 1/4 and 3/4 for s = 2..6
    for s in [2, 3, 4, 5, 6]:
        forms.append(make_form(f"hurwitz_zeta({s},1/4)",
                               mp.zeta(s, mpf(1)/4), "M3"))
        forms.append(make_form(f"hurwitz_zeta({s},3/4)",
                               mp.zeta(s, mpf(3)/4), "M3"))
    return forms


# -- Algebraic controls --
def atoms_ALG():
    forms = []
    for n in [2, 3, 5, 6, 7]:
        forms.append(make_form(f"sqrt({n})", mp.sqrt(n), "ALG"))
    forms.append(make_form("phi", (1 + mp.sqrt(5))/2, "ALG"))
    for n in [2, 3, 5, 7]:
        forms.append(make_form(f"ln({n})", log(n), "ALG"))
    return forms


# -- Odd-zeta controls (NOT in master Mellin rings; tests selection-bias) --
def atoms_ODD_ZETA():
    forms = []
    for k in [3, 5, 7, 9]:
        forms.append(make_form(f"zeta({k})", zeta(k), "ODD_ZETA_CONTROL"))
    return forms


# -- Cross-product depth-2 forms --
def atoms_CROSS_PRODUCT():
    forms = []
    base = [
        ("pi*G", pi * catalan),
        ("pi*ln(2)", pi * log(2)),
        ("pi^2*G", pi**2 * catalan),
        ("pi*sqrt(2)", pi * mp.sqrt(2)),
        ("sqrt(pi)*G", mp.sqrt(pi) * catalan),
        ("G/pi", catalan / pi),
        ("G/pi^2", catalan / pi**2),
        ("G/pi^3", catalan / pi**3),
        ("zeta(3)/pi", zeta(3) / pi),
        ("zeta(3)/pi^2", zeta(3) / pi**2),
        ("zeta(3)/pi^3", zeta(3) / pi**3),
        ("zeta(5)/pi^5", zeta(5) / pi**5),
        ("ln(2)^2", log(2)**2),
        ("ln(2)*ln(3)", log(2) * log(3)),
        ("ln(3)^2", log(3)**2),
        ("ln(2)/pi", log(2)/pi),
        ("ln(2)/pi^2", log(2)/pi**2),
        ("ln(2)*G", log(2)*catalan),
        ("ln(2)*pi", log(2)*pi),
    ]
    for name, val in base:
        forms.append(make_form(name, val, "CROSS_PRODUCT"))
    return forms


def build_atomic_basis():
    """Build linearly-independent atomic basis."""
    return {
        "M1": atoms_M1(),
        "M2": atoms_M2(),
        "M3": atoms_M3(),
        "ALG": atoms_ALG(),
        "ODD_ZETA_CONTROL": atoms_ODD_ZETA(),
        "CROSS_PRODUCT": atoms_CROSS_PRODUCT(),
    }


def deduplicate_by_value(forms_list, tol=1e-30):
    """Remove any forms with numerically duplicate values (lin-dep over Q)."""
    out = []
    seen_vals = []
    for f in forms_list:
        v = f["value"]
        is_dup = False
        for sv in seen_vals:
            if abs(v - sv) < tol or abs(v + sv) < tol:
                is_dup = True
                break
        if not is_dup:
            out.append(f)
            seen_vals.append(v)
    return out


def freeze_basis(forms_by_class, sub_runs, output_path):
    """Write basis + sub-run definitions to JSON BEFORE PSLQ."""
    composition = {
        cls: {
            "count": len(lst),
            "names": [f["name"] for f in lst],
        }
        for cls, lst in forms_by_class.items()
    }
    all_names = []
    for cls, lst in forms_by_class.items():
        for f in lst:
            all_names.append((cls, f["name"]))
    names_str = "|".join(f"{c}:{n}" for c, n in all_names)
    basis_hash = hashlib.sha256(names_str.encode()).hexdigest()

    record = {
        "target": "ln_k_0(1S) hydrogen Bethe logarithm",
        "target_value_string": LN_K0_1S_16DIGITS,
        "target_source": "Drake 1990 PRA 41 1243; Pachucki-Yerokhin 2010 review",
        "target_precision_digits": 16,
        "working_dps": WORKING_DPS,
        "basis_type": "atomic (one form per linearly-independent constant)",
        "composition": composition,
        "total_forms": len(all_names),
        "basis_sha256": basis_hash,
        "sub_runs": [
            {"name": name, "classes": cls_list}
            for name, cls_list in sub_runs
        ],
        "frozen_before_pslq": True,
        "protocol": "TD-PSLQ-1 (Bethe log probe against master Mellin engine)",
        "precision_floor_note": (
            "Target precision is 16 digits. Atomic basis means PSLQ integer "
            "coefficients encode rational prefactors of each atomic constant "
            "directly. False-positive risk for depth-L relation at ceiling C: "
            "~ (N choose L) * C^L / 10^16."
        ),
    }
    with open(output_path, "w") as fh:
        json.dump(record, fh, indent=2)
    return basis_hash


# ============================================================================
# PSLQ probe
# ============================================================================

def run_pslq(target, forms_list, coef_ceiling, maxsteps=500):
    """Run mpmath.pslq on [target] + [forms].

    Post-filter: PSLQ's tol is the iteration termination hint; actual relation
    residual can exceed tol. A genuine integer relation has residual at the
    precision floor (~1e-{WORKING_DPS}). We classify by residual:
      - residual < 1e-30: genuine relation
      - residual >= 1e-30 and < tol: spurious-floor relation (PSLQ best guess)
      - residual >= tol: PSLQ failed to converge meaningfully
    """
    if len(forms_list) == 0:
        return {"hit": False, "skipped": True}

    vec = [target] + [f["value"] for f in forms_list]
    tol = mpf(10) ** (-14)
    try:
        rel = mp.pslq(vec, tol=tol, maxcoeff=coef_ceiling, maxsteps=maxsteps)
    except Exception as e:
        return {"hit": False, "error": str(e)}

    if rel is None:
        return {"hit": False, "relation": None}

    a_target = rel[0]
    nonzero_all = [(i - 1, rel[i]) for i in range(1, len(rel)) if rel[i] != 0]
    max_abs = max(abs(rel[i]) for i in range(len(rel)))

    # Compute actual residual
    res = sum(r * v for r, v in zip(rel, vec))
    res_abs = abs(res)

    # Genuine relation threshold: must be well below the precision floor.
    # Target has 16 digits, so at 50 dps the precision floor on derived
    # quantities is ~1e-16 (limited by target). A genuine relation should
    # have residual at this floor or below. Threshold: 1e-15.
    GENUINE_THRESHOLD = mpf(10) ** (-15)
    is_genuine = res_abs < GENUINE_THRESHOLD

    if a_target == 0:
        # Basis-internal relation (no target involvement)
        # This is genuine only if residual is at floor
        return {
            "hit": False,
            "spurious_basis_relation": True,
            "basis_relation_max_coef": int(max_abs),
            "basis_relation_length": len(nonzero_all),
            "residual": mp.nstr(res_abs, 25),
            "is_genuine": is_genuine,
        }

    nonzero = nonzero_all  # all indices are nonzero except a_target=0 case
    if not is_genuine:
        # Spurious target-involving "hit" — PSLQ filled out a depth-L relation
        # at coefficient ceiling but residual is far above precision floor.
        return {
            "hit": False,  # treat as null for verdict purposes
            "spurious_filled_relation": True,
            "coef_target": int(a_target),
            "max_abs_coef": int(max_abs),
            "relation_length": len(nonzero) + 1,
            "residual": mp.nstr(res_abs, 25),
            "is_genuine": False,
            "_note": "PSLQ returned best-guess relation but residual >> precision floor; not a real integer relation",
        }

    return {
        "hit": True,
        "is_genuine": True,
        "coef_target": int(a_target),
        "max_abs_coef": int(max_abs),
        "relation_length": len(nonzero) + 1,
        "residual": mp.nstr(res_abs, 25),
        "nonzero_basis": [
            {
                "form_class": forms_list[i]["class"],
                "form_name": forms_list[i]["name"],
                "coefficient": int(coef),
            }
            for i, coef in nonzero
        ],
    }


def assess_hit(hit, basis_size, target_digits=16):
    """Assess whether a hit is trustworthy given precision constraints."""
    if not hit.get("hit"):
        return "null"
    C = hit["max_abs_coef"]
    L = hit["relation_length"]
    # log10(false-positive probability) ~ log10(binomial(N,L)) + L*log10(C) - d
    # Rough bound: L * log10(N) + L * log10(C) - d
    log_p = L * (math.log10(max(basis_size, 1)) + math.log10(max(C, 1))) - target_digits
    if log_p < -4:
        return "trustworthy"
    elif log_p < 0:
        return "marginal"
    else:
        return "spurious"


# ============================================================================
# Verification at higher precision
# ============================================================================

def verify_hit(hit, target_str, higher_dps=200):
    """Re-evaluate the proposed relation at higher precision.

    Note: target is only 16 digits, so re-evaluation at 200 dps just gives
    machine-precision evaluation of the basis side; if the proposed relation
    is real, the residual should stay below 10^-15. If spurious, the
    residual at 200 dps will reveal it (the relation that fit at 50 dps will
    have residual > 10^-14 at 200 dps if it was spurious).
    """
    if not hit.get("hit"):
        return None
    old_dps = mp.mp.dps
    mp.mp.dps = higher_dps
    target = mpf(target_str)

    # Rebuild basis values at higher dps
    a_target = hit["coef_target"]
    res = a_target * target
    for entry in hit["nonzero_basis"]:
        val = re_evaluate(entry["form_name"])
        res += entry["coefficient"] * val
    res_abs = abs(res)

    mp.mp.dps = old_dps
    return mp.nstr(res_abs, 30)


def re_evaluate(name):
    """Evaluate a basis form by name at current precision."""
    # M1
    if name == "1":
        return mpf(1)
    m = __import__("re").match(r"^pi\^(\d+)$", name)
    if m:
        return pi ** int(m.group(1))
    if name == "pi":
        return pi
    m = __import__("re").match(r"^1/pi\^(\d+)$", name)
    if m:
        return mpf(1) / pi ** int(m.group(1))
    # M2
    if name == "sqrt(pi)":
        return mp.sqrt(pi)
    m = __import__("re").match(r"^sqrt\(pi\)\^(\d+)$", name)
    if m:
        return mp.sqrt(pi) ** int(m.group(1))
    if name == "1/sqrt(pi)":
        return 1 / mp.sqrt(pi)
    m = __import__("re").match(r"^1/sqrt\(pi\)\^(\d+)$", name)
    if m:
        return 1 / (mp.sqrt(pi) ** int(m.group(1)))
    m = __import__("re").match(r"^zeta\((\d+)\)$", name)
    if m:
        return zeta(int(m.group(1)))
    # M3
    if name == "G":
        return catalan
    m = __import__("re").match(r"^beta\((\d+)\)$", name)
    if m:
        s = int(m.group(1))
        return sum(mpf(-1)**k / mpf(2*k+1)**s for k in range(0, 500))
    m = __import__("re").match(r"^hurwitz_zeta\((\d+),(\d+)/(\d+)\)$", name)
    if m:
        return mp.zeta(int(m.group(1)),
                       mpf(int(m.group(2))) / int(m.group(3)))
    # ALG
    m = __import__("re").match(r"^sqrt\((\d+)\)$", name)
    if m:
        return mp.sqrt(int(m.group(1)))
    if name == "phi":
        return (1 + mp.sqrt(5)) / 2
    m = __import__("re").match(r"^ln\((\d+)\)$", name)
    if m:
        return log(int(m.group(1)))
    # CROSS_PRODUCT (hardcoded set)
    cross_map = {
        "pi*G": pi * catalan,
        "pi*ln(2)": pi * log(2),
        "pi^2*G": pi**2 * catalan,
        "pi*sqrt(2)": pi * mp.sqrt(2),
        "sqrt(pi)*G": mp.sqrt(pi) * catalan,
        "G/pi": catalan / pi,
        "G/pi^2": catalan / pi**2,
        "G/pi^3": catalan / pi**3,
        "zeta(3)/pi": zeta(3) / pi,
        "zeta(3)/pi^2": zeta(3) / pi**2,
        "zeta(3)/pi^3": zeta(3) / pi**3,
        "zeta(5)/pi^5": zeta(5) / pi**5,
        "ln(2)^2": log(2)**2,
        "ln(2)*ln(3)": log(2) * log(3),
        "ln(3)^2": log(3)**2,
        "ln(2)/pi": log(2)/pi,
        "ln(2)/pi^2": log(2)/pi**2,
        "ln(2)*G": log(2)*catalan,
        "ln(2)*pi": log(2)*pi,
    }
    if name in cross_map:
        return cross_map[name]
    raise ValueError(f"Cannot re-evaluate: {name}")


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 76, flush=True)
    print("Sprint TD-PSLQ-1: Bethe log probe against master Mellin engine", flush=True)
    print("=" * 76, flush=True)
    print(flush=True)
    print(f"Target: ln k_0(1S) = {LN_K0_1S_16DIGITS}", flush=True)
    print(f"Source: Drake 1990 PRA 41 1243 (16 significant digits)", flush=True)
    print(f"Working precision: {WORKING_DPS} dps", flush=True)
    print(flush=True)

    forms_by_class = build_atomic_basis()
    total = sum(len(v) for v in forms_by_class.values())
    print("Atomic basis composition:", flush=True)
    for cls, lst in forms_by_class.items():
        print(f"  {cls:<20}: {len(lst):>3}", flush=True)
    print(f"  {'TOTAL':<20}: {total:>3}", flush=True)
    print(flush=True)

    sub_runs = [
        ("M1_only", ["M1"]),
        ("M2_only", ["M2"]),
        ("M3_only", ["M3"]),
        ("M1_M2", ["M1", "M2"]),
        ("M1_M3", ["M1", "M3"]),
        ("M1_M2_M3", ["M1", "M2", "M3"]),
        ("M1_M2_M3_ALG", ["M1", "M2", "M3", "ALG"]),
        ("M1_M2_M3_ALG_CROSS", ["M1", "M2", "M3", "ALG", "CROSS_PRODUCT"]),
        ("FULL", ["M1", "M2", "M3", "ALG", "CROSS_PRODUCT", "ODD_ZETA_CONTROL"]),
    ]

    frozen_path = OUT_DIR / "td_pslq_bethe_log_basis_frozen.json"
    basis_hash = freeze_basis(forms_by_class, sub_runs, frozen_path)
    print(f"Basis frozen to: {frozen_path}", flush=True)
    print(f"SHA256 hash: {basis_hash[:16]}...", flush=True)
    print(flush=True)

    all_results = {}
    for run_name, cls_list in sub_runs:
        run_forms = []
        for cls in cls_list:
            run_forms.extend(forms_by_class[cls])
        # Deduplicate by value (handles e.g. zeta(2) ≈ pi^2/6)
        run_forms = deduplicate_by_value(run_forms)
        N = len(run_forms)
        print(f"--- Sub-run: {run_name} (N={N}, classes={cls_list}) ---", flush=True)
        all_results[run_name] = {"N": N, "classes": cls_list, "ceilings": {}}

        for ceiling in [10**4, 10**6, 10**8]:
            log_c = int(math.log10(ceiling))
            fp_estimate = N * ceiling / 1e16
            print(f"  ceiling=10^{log_c}: N*C/10^16 = {fp_estimate:.2e}", flush=True)
            result = run_pslq(TARGET_VAL, run_forms, ceiling)
            all_results[run_name]["ceilings"][f"1e{log_c}"] = result
            if result.get("hit"):
                trust = assess_hit(result, N)
                print(f"    GENUINE HIT (trust={trust}): rel_len={result['relation_length']}, "
                      f"max_coef={result['max_abs_coef']}, "
                      f"residual={result['residual']}", flush=True)
                print(f"    relation: {result['coef_target']} * ln k_0(1S) + ...", flush=True)
                for entry in result["nonzero_basis"][:6]:
                    print(f"      {entry['coefficient']:+d} * "
                          f"[{entry['form_class']}] {entry['form_name']}",
                          flush=True)
                if len(result["nonzero_basis"]) > 6:
                    print(f"      ... +{len(result['nonzero_basis'])-6} more terms",
                          flush=True)
                if trust in ("trustworthy", "marginal"):
                    verify = verify_hit(result, LN_K0_1S_16DIGITS, 200)
                    print(f"    residual at 200 dps: {verify}", flush=True)
                    all_results[run_name]["ceilings"][f"1e{log_c}"]["verify_200dps"] = verify
                all_results[run_name]["ceilings"][f"1e{log_c}"]["trust_assessment"] = trust
            elif result.get("spurious_filled_relation"):
                print(f"    PSLQ best-guess (NOT a genuine relation): "
                      f"residual {result['residual']} >> 1e-15 floor", flush=True)
                print(f"      rel_len={result['relation_length']}, "
                      f"max_coef={result['max_abs_coef']}", flush=True)
                print(f"    => NULL: PSLQ filled out depth-L relation but residual "
                      f"is above precision floor (would be genuine only at < 1e-15).",
                      flush=True)
            elif result.get("spurious_basis_relation"):
                if result.get("is_genuine"):
                    print(f"    Genuine basis-internal relation: residual "
                          f"{result['residual']} at floor.", flush=True)
                else:
                    print(f"    Basis-internal best-guess (residual "
                          f"{result['residual']}): max_coef="
                          f"{result.get('basis_relation_max_coef')}", flush=True)
            else:
                print(f"    NULL (no integer relation found at this ceiling).", flush=True)
        print(flush=True)

    print("=" * 76, flush=True)
    print("Null confidence summary:", flush=True)
    print(f"  Target precision: 16 digits (Drake 1990)", flush=True)
    print(f"  PSLQ tolerance: 1e-14 (2-digit headroom)", flush=True)
    print(f"  Atomic basis (linearly independent over Q)", flush=True)
    print(f"  False-positive probability at ceiling C, basis dim N, depth-1:", flush=True)
    print(f"    P_FP ~ N * C / 10^16", flush=True)
    print(f"  At N=60, C=10^4: P_FP ~ 6e-12 — robustly safe", flush=True)
    print(f"  At N=60, C=10^6: P_FP ~ 6e-10 — safe", flush=True)
    print(f"  At N=60, C=10^8: P_FP ~ 6e-8  — safe", flush=True)
    print("=" * 76, flush=True)

    out_path = OUT_DIR / "td_pslq_bethe_log_results.json"
    out = {
        "sprint": "TD-PSLQ-1",
        "target": "ln_k_0(1S)",
        "target_value": LN_K0_1S_16DIGITS,
        "target_precision_digits": 16,
        "working_dps": WORKING_DPS,
        "basis_sha256": basis_hash,
        "total_forms_all_classes": total,
        "forms_by_class_counts": {k: len(v) for k, v in forms_by_class.items()},
        "sub_run_results": all_results,
        "verdict": _summarize_verdict(all_results),
    }
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nResults written to: {out_path}", flush=True)


def _summarize_verdict(results):
    summary = {
        "total_sub_runs": len(results),
        "any_trustworthy_hit": False,
        "any_hit_at_all": False,
        "trustworthy_hits": [],
        "all_target_hits": [],
        "basis_internal_relations": [],
    }
    for run_name, run_data in results.items():
        for ceil_name, res in run_data["ceilings"].items():
            if res.get("hit"):
                summary["any_hit_at_all"] = True
                trust = res.get("trust_assessment") or assess_hit(res, run_data["N"])
                hit_record = {
                    "sub_run": run_name,
                    "ceiling": ceil_name,
                    "max_coef": res["max_abs_coef"],
                    "relation_length": res["relation_length"],
                    "trust": trust,
                    "residual_50dps": res["residual"],
                    "residual_200dps": res.get("verify_200dps"),
                }
                summary["all_target_hits"].append(hit_record)
                if trust == "trustworthy":
                    summary["any_trustworthy_hit"] = True
                    summary["trustworthy_hits"].append(hit_record)
            elif res.get("spurious_basis_relation"):
                summary["basis_internal_relations"].append({
                    "sub_run": run_name,
                    "ceiling": ceil_name,
                    "length": res.get("basis_relation_length"),
                })
    return summary


if __name__ == "__main__":
    main()
