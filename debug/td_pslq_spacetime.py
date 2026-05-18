"""
Sprint TD-PSLQ-2: PSLQ probe of the hydrogen 1S self-energy coefficient
A_60(1S) — a spacetime/relativistic Layer-2 constant — against the GeoVac
master Mellin engine ring (M1/M2/M3).

TARGET RATIONALE:
  TD-PSLQ-1 tested ln k_0(1S) and returned clean null. But the Bethe log is
  bound-state QED virtual-state sum, NOT spacetime/relativistic. It was the
  wrong test for the PI's "spacetime corrections close Layer-2 residuals to
  bit-exactness" hypothesis.

  A_60(1S) = -30.92415(1) is the alpha(Zalpha)^6 self-energy coefficient.
  Mechanism: pure relativistic-correction-to-self-energy beyond Bethe-log.
  Spacetime origin (Dirac-Coulomb kinematics + radiative correction).

  Source: Jentschura, Mohr, Soff PRL 82, 53 (1999); also tabulated in
  Mohr-Jentschura-Pachucki 2000 (arXiv:physics/0001068).
  Precision: 7 significant digits (uncertainty in last digit).

  This is the CANONICAL spacetime Layer-2 constant in atomic QED.

PROTOCOL: identical to TD-PSLQ-1 — atomic basis, frozen before PSLQ, post-
filter for genuine relations at precision floor 1e-15 (consistent with 7-digit
target precision at 50 dps working).

PRECISION FLOOR:
  Target precision 7 digits is TIGHTER than Bethe log's 16 digits in terms of
  allowable coefficient ceiling. At 7-digit precision:
    P_FP ~ N * C / 10^7 for depth-1
    For N=70, C=10^4: P_FP ~ 7e-2 — BORDERLINE
    For N=70, C=10^3: P_FP ~ 7e-3 — safe-ish
    For N=70, C=10^2: P_FP ~ 7e-4 — safe
  So at this target's precision, we restrict to ceiling 10^4 with strong
  residual post-filter (must be < 1e-6 = 2 digit headroom on the target).
  Higher ceilings (10^6, 10^8) are run for completeness but their results
  are interpreted as suggestive rather than definitive.

  NOTE: this is why a clean-null verdict at this precision is necessarily
  weaker than the Bethe-log null at 16 digits. The structural-skeleton-scope
  interpretation must explicitly acknowledge this.

Question: does A_60(1S) — the canonical spacetime/relativistic Layer-2
constant — sit in the master Mellin engine ring at PSLQ precision?

Outcome interpretation:
  - HIT at coefficient < 10^4: partial confirmation of PI hypothesis;
    structural-skeleton picture revises to "spacetime corrections are in M1/M2/M3."
  - HIT at coefficient > 10^4: marginal; treat as suggestive only.
  - NULL: PI hypothesis NOT supported for spacetime self-energy channel;
    structural-skeleton-scope reading sharpens further.
"""

import json
import hashlib
import math
from pathlib import Path

import mpmath as mp
from mpmath import mpf, pi, log, sqrt, zeta, catalan


# ============================================================================
# TARGET: A_60(1S), the alpha(Zalpha)^6 self-energy coefficient
# ============================================================================
# Jentschura, Mohr, Soff PRL 82, 53 (1999); Mohr-Jentschura-Pachucki 2000:
#   A_60(1S) = -30.924 15(1)
A_60_1S_VALUE = "-30.92415"
A_60_1S_UNCERTAINTY = 1  # in units of last digit (i.e. ±0.00001)
A_60_1S_DIGITS = 7

WORKING_DPS = 50
mp.mp.dps = WORKING_DPS

TARGET_VAL = mpf(A_60_1S_VALUE)

OUT_DIR = Path("debug/data")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================================
# ATOMIC BASIS — extended for spacetime/relativistic content
# ============================================================================
# Same M1/M2/M3 atomic forms as TD-PSLQ-1, plus:
#   - Extra rationals (A_60 is naturally of magnitude ~30, so large-rational
#     identifications like 30 + small_correction are plausible)
#   - Powers of 2 and 3 (common in relativistic-correction closed forms;
#     e.g., 4/3, 11/24, 35/72 patterns from Dirac-Coulomb expansions)
#   - Bernoulli-like rationals (relativistic kinematic expansions)


def make_form(name, value, mech_class):
    return {"name": name, "value": value, "class": mech_class}


# -- M1 atomic forms: powers of pi and 1/pi --
def atoms_M1():
    forms = [make_form("1", mpf(1), "RAT")]
    for p in [1, 2, 3, 4, 5, 6]:
        forms.append(make_form(f"pi^{p}", pi**p, "M1"))
    for p in [1, 2, 3, 4]:
        forms.append(make_form(f"1/pi^{p}", 1/pi**p, "M1"))
    return forms


# -- M2 atomic forms: sqrt(pi) tower --
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
    beta_4 = sum(mpf(-1)**k / mpf(2*k+1)**4 for k in range(0, 500))
    beta_6 = sum(mpf(-1)**k / mpf(2*k+1)**6 for k in range(0, 500))
    beta_8 = sum(mpf(-1)**k / mpf(2*k+1)**8 for k in range(0, 500))
    forms.append(make_form("beta(4)", beta_4, "M3"))
    forms.append(make_form("beta(6)", beta_6, "M3"))
    forms.append(make_form("beta(8)", beta_8, "M3"))
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


# -- Odd-zeta controls (NOT in master Mellin rings) --
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
        # Additional spacetime-flavored cross products:
        ("pi^2*ln(2)", pi**2 * log(2)),
        ("pi^3*ln(2)", pi**3 * log(2)),
        ("pi^4/G", pi**4 / catalan),
        ("zeta(3)*pi", zeta(3) * pi),
        ("zeta(3)*pi^2", zeta(3) * pi**2),
        ("zeta(3)*ln(2)", zeta(3) * log(2)),
        ("ln(2)*sqrt(pi)", log(2) * mp.sqrt(pi)),
        ("G*sqrt(2)", catalan * mp.sqrt(2)),
    ]
    for name, val in base:
        forms.append(make_form(name, val, "CROSS_PRODUCT"))
    return forms


# -- Spacetime-flavored augmentations: large-magnitude pi and pi*ln(2)
#    combinations that arise naturally in relativistic-correction closed forms.
#    A_60 ~ -31, so we explicitly include forms with that magnitude scale
#    so PSLQ can identify a dominant single-term hit.
def atoms_SPACETIME_AUG():
    forms = []
    base = [
        # Large pi^k naturally near |A_60| ~ 30:
        ("pi^3", pi**3),                   # ~31.0
        ("10*pi^2", 10*pi**2),             # ~98.7
        ("pi^4/3", pi**4 / 3),             # ~32.5
        ("pi^4/4", pi**4 / 4),             # ~24.4
        ("10*log(2)*pi", 10*log(2)*pi),    # ~21.8
        # Bernoulli-style relativistic-kinematic constants:
        ("4/3", mpf(4)/3),
        ("11/24", mpf(11)/24),
        ("35/72", mpf(35)/72),
        ("139/144", mpf(139)/144),
        # Dirac-Coulomb closed-form patterns:
        ("pi^2*log(2)/6", pi**2*log(2)/6),
        ("pi^2*log(2)/3", pi**2*log(2)/3),
        ("pi^2/12", pi**2/12),
        ("7*zeta(3)/8", 7*zeta(3)/8),
        ("7*pi^4/360", 7*pi**4/360),
        ("31*pi^2/216", 31*pi**2/216),
        # Mixed:
        ("pi^2*ln(2)/2", pi**2*log(2)/2),
        ("3*zeta(3)/2", 3*zeta(3)/2),
        ("9/2", mpf(9)/2),
    ]
    for name, val in base:
        forms.append(make_form(name, val, "SPACETIME_AUG"))
    return forms


def build_atomic_basis():
    return {
        "M1": atoms_M1(),
        "M2": atoms_M2(),
        "M3": atoms_M3(),
        "ALG": atoms_ALG(),
        "ODD_ZETA_CONTROL": atoms_ODD_ZETA(),
        "CROSS_PRODUCT": atoms_CROSS_PRODUCT(),
        "SPACETIME_AUG": atoms_SPACETIME_AUG(),
    }


def deduplicate_by_value(forms_list, tol=1e-30):
    """Remove forms with numerically duplicate values."""
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
        "target": "A_60(1S) hydrogen alpha(Zalpha)^6 self-energy coefficient",
        "target_value_string": A_60_1S_VALUE,
        "target_source": (
            "Jentschura, Mohr, Soff PRL 82, 53 (1999); "
            "Mohr-Jentschura-Pachucki 2000 (arXiv:physics/0001068); "
            "Yerokhin-Pachucki-Patkos 2019 review (arXiv:1809.00462)"
        ),
        "target_precision_digits": A_60_1S_DIGITS,
        "target_uncertainty_last_digit": A_60_1S_UNCERTAINTY,
        "target_mechanism": "spacetime/relativistic — Dirac-Coulomb kinematics + one-loop radiative",
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
        "protocol": "TD-PSLQ-2 (spacetime probe of A_60(1S) against master Mellin engine)",
        "precision_floor_note": (
            "Target precision is 7 digits — tighter than TD-PSLQ-1's 16-digit "
            "Bethe log. At ceiling 10^4, depth-1 false-positive probability "
            "~ N * 10^4 / 10^7 = N * 1e-3. For N=80: ~8%. Borderline. "
            "Genuine relations must have residual < 1e-6 (1 digit headroom)."
        ),
    }
    with open(output_path, "w") as fh:
        json.dump(record, fh, indent=2)
    return basis_hash


# ============================================================================
# PSLQ probe
# ============================================================================

def run_pslq(target, forms_list, coef_ceiling, maxsteps=500):
    """Run mpmath.pslq on [target] + [forms]."""
    if len(forms_list) == 0:
        return {"hit": False, "skipped": True}

    vec = [target] + [f["value"] for f in forms_list]
    # PSLQ tolerance: 7-digit precision target → tol=1e-6 (1 digit headroom)
    # But want PSLQ to iterate to find any candidate, so set tol modestly tight
    tol = mpf(10) ** (-6)
    try:
        rel = mp.pslq(vec, tol=tol, maxcoeff=coef_ceiling, maxsteps=maxsteps)
    except Exception as e:
        return {"hit": False, "error": str(e)}

    if rel is None:
        return {"hit": False, "relation": None}

    a_target = rel[0]
    nonzero_all = [(i - 1, rel[i]) for i in range(1, len(rel)) if rel[i] != 0]
    max_abs = max(abs(rel[i]) for i in range(len(rel)))

    res = sum(r * v for r, v in zip(rel, vec))
    res_abs = abs(res)

    # GENUINE_THRESHOLD: tighter than tol. For 7-digit target, want residual
    # at least 1 digit below target precision floor: |residual| < 1e-6.
    # A truly real integer relation at this target precision can have
    # residual no smaller than ~ |max_coef| * 1e-7 (rounding limit of target).
    # So we use 1e-6 as our genuine cutoff.
    GENUINE_THRESHOLD = mpf(10) ** (-6)
    is_genuine = res_abs < GENUINE_THRESHOLD

    if a_target == 0:
        return {
            "hit": False,
            "spurious_basis_relation": True,
            "basis_relation_max_coef": int(max_abs),
            "basis_relation_length": len(nonzero_all),
            "residual": mp.nstr(res_abs, 25),
            "is_genuine": is_genuine,
        }

    nonzero = nonzero_all
    if not is_genuine:
        return {
            "hit": False,
            "spurious_filled_relation": True,
            "coef_target": int(a_target),
            "max_abs_coef": int(max_abs),
            "relation_length": len(nonzero) + 1,
            "residual": mp.nstr(res_abs, 25),
            "is_genuine": False,
            "_note": "PSLQ best-guess but residual >> 1e-6 floor; not a real relation",
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


def assess_hit(hit, basis_size, target_digits=7):
    """Assess trustworthiness given precision constraints."""
    if not hit.get("hit"):
        return "null"
    C = hit["max_abs_coef"]
    L = hit["relation_length"]
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

    Since target is only 7 digits, this re-evaluates the BASIS side at higher
    precision and checks whether the linear combination of high-precision basis
    constants gives a value within target's 7-digit precision floor.
    """
    if not hit.get("hit"):
        return None
    old_dps = mp.mp.dps
    mp.mp.dps = higher_dps
    target = mpf(target_str)

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
    import re as re_mod
    if name == "1":
        return mpf(1)
    m = re_mod.match(r"^pi\^(\d+)$", name)
    if m:
        return pi ** int(m.group(1))
    if name == "pi":
        return pi
    m = re_mod.match(r"^1/pi\^(\d+)$", name)
    if m:
        return mpf(1) / pi ** int(m.group(1))
    if name == "sqrt(pi)":
        return mp.sqrt(pi)
    m = re_mod.match(r"^sqrt\(pi\)\^(\d+)$", name)
    if m:
        return mp.sqrt(pi) ** int(m.group(1))
    if name == "1/sqrt(pi)":
        return 1 / mp.sqrt(pi)
    m = re_mod.match(r"^1/sqrt\(pi\)\^(\d+)$", name)
    if m:
        return 1 / (mp.sqrt(pi) ** int(m.group(1)))
    m = re_mod.match(r"^zeta\((\d+)\)$", name)
    if m:
        return zeta(int(m.group(1)))
    if name == "G":
        return catalan
    m = re_mod.match(r"^beta\((\d+)\)$", name)
    if m:
        s = int(m.group(1))
        return sum(mpf(-1)**k / mpf(2*k+1)**s for k in range(0, 500))
    m = re_mod.match(r"^hurwitz_zeta\((\d+),(\d+)/(\d+)\)$", name)
    if m:
        return mp.zeta(int(m.group(1)),
                       mpf(int(m.group(2))) / int(m.group(3)))
    m = re_mod.match(r"^sqrt\((\d+)\)$", name)
    if m:
        return mp.sqrt(int(m.group(1)))
    if name == "phi":
        return (1 + mp.sqrt(5)) / 2
    m = re_mod.match(r"^ln\((\d+)\)$", name)
    if m:
        return log(int(m.group(1)))
    cross_and_aug = {
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
        "pi^2*ln(2)": pi**2 * log(2),
        "pi^3*ln(2)": pi**3 * log(2),
        "pi^4/G": pi**4 / catalan,
        "zeta(3)*pi": zeta(3) * pi,
        "zeta(3)*pi^2": zeta(3) * pi**2,
        "zeta(3)*ln(2)": zeta(3) * log(2),
        "ln(2)*sqrt(pi)": log(2) * mp.sqrt(pi),
        "G*sqrt(2)": catalan * mp.sqrt(2),
        "10*pi^2": 10*pi**2,
        "pi^4/3": pi**4 / 3,
        "pi^4/4": pi**4 / 4,
        "10*log(2)*pi": 10*log(2)*pi,
        "4/3": mpf(4)/3,
        "11/24": mpf(11)/24,
        "35/72": mpf(35)/72,
        "139/144": mpf(139)/144,
        "pi^2*log(2)/6": pi**2*log(2)/6,
        "pi^2*log(2)/3": pi**2*log(2)/3,
        "pi^2/12": pi**2/12,
        "7*zeta(3)/8": 7*zeta(3)/8,
        "7*pi^4/360": 7*pi**4/360,
        "31*pi^2/216": 31*pi**2/216,
        "pi^2*ln(2)/2": pi**2*log(2)/2,
        "3*zeta(3)/2": 3*zeta(3)/2,
        "9/2": mpf(9)/2,
    }
    if name in cross_and_aug:
        return cross_and_aug[name]
    raise ValueError(f"Cannot re-evaluate: {name}")


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 76, flush=True)
    print("Sprint TD-PSLQ-2: spacetime probe of A_60(1S) against master Mellin engine", flush=True)
    print("=" * 76, flush=True)
    print(flush=True)
    print(f"Target: A_60(1S) = {A_60_1S_VALUE}", flush=True)
    print(f"Source: Jentschura-Mohr-Soff 1999; Mohr-Jentschura-Pachucki 2000", flush=True)
    print(f"Target precision: {A_60_1S_DIGITS} digits (uncertainty ±0.00001)", flush=True)
    print(f"Mechanism: alpha(Zalpha)^6 self-energy — spacetime/relativistic", flush=True)
    print(f"Working precision: {WORKING_DPS} dps", flush=True)
    print(f"Genuine-relation threshold: |residual| < 1e-6", flush=True)
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
        ("M1_M2_M3_ALG_CROSS_AUG", ["M1", "M2", "M3", "ALG", "CROSS_PRODUCT", "SPACETIME_AUG"]),
        ("FULL", ["M1", "M2", "M3", "ALG", "CROSS_PRODUCT", "SPACETIME_AUG", "ODD_ZETA_CONTROL"]),
    ]

    frozen_path = OUT_DIR / "td_pslq_spacetime_basis_frozen.json"
    basis_hash = freeze_basis(forms_by_class, sub_runs, frozen_path)
    print(f"Basis frozen to: {frozen_path}", flush=True)
    print(f"SHA256 hash: {basis_hash[:16]}...", flush=True)
    print(flush=True)

    all_results = {}
    for run_name, cls_list in sub_runs:
        run_forms = []
        for cls in cls_list:
            run_forms.extend(forms_by_class[cls])
        run_forms = deduplicate_by_value(run_forms)
        N = len(run_forms)
        print(f"--- Sub-run: {run_name} (N={N}, classes={cls_list}) ---", flush=True)
        all_results[run_name] = {"N": N, "classes": cls_list, "ceilings": {}}

        for ceiling in [10**4, 10**6, 10**8]:
            log_c = int(math.log10(ceiling))
            fp_estimate = N * ceiling / 10**A_60_1S_DIGITS
            print(f"  ceiling=10^{log_c}: N*C/10^{A_60_1S_DIGITS} = {fp_estimate:.2e}", flush=True)
            result = run_pslq(TARGET_VAL, run_forms, ceiling)
            all_results[run_name]["ceilings"][f"1e{log_c}"] = result
            if result.get("hit"):
                trust = assess_hit(result, N)
                print(f"    GENUINE HIT (trust={trust}): rel_len={result['relation_length']}, "
                      f"max_coef={result['max_abs_coef']}, "
                      f"residual={result['residual']}", flush=True)
                print(f"    relation: {result['coef_target']} * A_60(1S) + ...", flush=True)
                for entry in result["nonzero_basis"][:6]:
                    print(f"      {entry['coefficient']:+d} * "
                          f"[{entry['form_class']}] {entry['form_name']}",
                          flush=True)
                if len(result["nonzero_basis"]) > 6:
                    print(f"      ... +{len(result['nonzero_basis'])-6} more terms",
                          flush=True)
                if trust in ("trustworthy", "marginal"):
                    verify = verify_hit(result, A_60_1S_VALUE, 200)
                    print(f"    residual at 200 dps (basis-side): {verify}", flush=True)
                    all_results[run_name]["ceilings"][f"1e{log_c}"]["verify_200dps"] = verify
                all_results[run_name]["ceilings"][f"1e{log_c}"]["trust_assessment"] = trust
            elif result.get("spurious_filled_relation"):
                print(f"    PSLQ best-guess (NOT a genuine relation): "
                      f"residual {result['residual']} >> 1e-6 floor", flush=True)
                print(f"      rel_len={result['relation_length']}, "
                      f"max_coef={result['max_abs_coef']}", flush=True)
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
    print(f"  Target precision: 7 digits (Jentschura-Mohr 1999)", flush=True)
    print(f"  PSLQ tolerance: 1e-6 (1-digit headroom)", flush=True)
    print(f"  Genuine-relation threshold: |residual| < 1e-6", flush=True)
    print(f"  False-positive probability at ceiling C, basis dim N, depth-1:", flush=True)
    print(f"    P_FP ~ N * C / 10^7", flush=True)
    print(f"  At N=80, C=10^4: P_FP ~ 8e-2 — BORDERLINE", flush=True)
    print(f"  At N=80, C=10^6: P_FP ~ 8e+0 — UNRELIABLE", flush=True)
    print(f"  At N=80, C=10^8: P_FP ~ 8e+2 — UNRELIABLE", flush=True)
    print(f"  ==> verdict depends primarily on C=10^4 results", flush=True)
    print("=" * 76, flush=True)

    out_path = OUT_DIR / "td_pslq_spacetime_results.json"
    out = {
        "sprint": "TD-PSLQ-2",
        "target": "A_60(1S)",
        "target_value": A_60_1S_VALUE,
        "target_precision_digits": A_60_1S_DIGITS,
        "target_mechanism": "spacetime/relativistic self-energy alpha(Zalpha)^6",
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
