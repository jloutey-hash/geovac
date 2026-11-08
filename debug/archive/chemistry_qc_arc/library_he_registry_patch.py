"""Targeted, re-read-immediately-before patch of the He qubit keys in
debug/qa/numeric_registry.py (another agent edits the same file concurrently).
Each OLD block must match exactly once or the script aborts without writing."""
import os
import sys

P = os.path.join(os.path.dirname(__file__), "qa", "numeric_registry.py")
src = open(P, encoding="utf-8").read()

PROV = ("MEASURED 2026-09-25 debug/library_he_fix.py, hamiltonian('He') = "
        "LatticeIndex(slater_full, h1_method='exact') -- helium's PHYSICAL h1; "
        "2nd route composed builder on the He spec agrees")
EDITS = [
    # ---- Pauli ----
    ('''    # ---- He atomic series, exact global-M_L ERI rule --------------------
    "he_n2_pauli": dict(
        value=287, convention="non-identity Pauli terms", q=10,
        provenance="MEASURED 2026-08-30, JordanWignerEncoder on "
                   "LatticeIndex(He, max_n=2, slater_full)",
        aliases={288: "including identity"}),
    "he_n3_pauli": dict(
        value=14078, convention="non-identity Pauli terms", q=28,
        provenance="MEASURED 2026-08-30, same route",
        aliases={14079: "including identity"}),
    "he_n4_pauli": dict(
        value=250402, convention="non-identity Pauli terms", q=60,
        provenance="MEASURED 2026-08-30, same route",
        aliases={250403: "including identity"}),
    "he_n5_pauli": dict(
        value=2434441, convention="non-identity Pauli terms", q=110,
        provenance="MEASURED 2026-08-30, ecosystem hamiltonian('He', max_n=5); "
                   "28 min build",
        aliases={2434442: "including identity"}),
''',
     '''    # ---- He atomic series, exact global-M_L ERI rule, PHYSICAL h1 --------
    # Re-measured 2026-09-25: until then the library He (and these keys)
    # carried LatticeIndex(h1_method='hybrid'), which adds a non-physical
    # +1/16 hop on every GeometricLattice edge to the one-body operator.
    # In the k=Z hydrogenic basis T+V_ne is exactly diagonal, so the hop is
    # not helium's.  Retired hop-model values are in RETIRED below.
    "he_n2_pauli": dict(
        value=279, convention="non-identity Pauli terms", q=10,
        provenance=PROV_HE + " exactly (279)",
        aliases={280: "including identity"}),
    "he_n3_pauli": dict(
        value=14046, convention="non-identity Pauli terms", q=28,
        provenance=PROV_HE + " exactly (14046)",
        aliases={14047: "including identity"}),
    "he_n4_pauli": dict(
        value=250322, convention="non-identity Pauli terms", q=60,
        provenance="MEASURED 2026-09-25 debug/library_he_fix.py (physical h1)",
        aliases={250323: "including identity"}),
    "he_n5_pauli": dict(
        value=2434281, convention="non-identity Pauli terms", q=110,
        provenance="MEASURED 2026-09-25 debug/library_he_fix.py (physical h1; "
                   "~65 min build under load)",
        aliases={2434282: "including identity"}),
'''),
    # ---- lambda ----
    ('''    "he_n2_lambda": dict(
        value=11.175, convention="non-identity 1-norm (Ha)", q=10,
        provenance="MEASURED 2026-08-30"),
    "he_n3_lambda": dict(
        value=74.207, convention="non-identity 1-norm (Ha)", q=28,
        provenance="MEASURED 2026-08-30"),
    "he_n4_lambda": dict(
        value=275.718, convention="non-identity 1-norm (Ha)", q=60,
        provenance="MEASURED 2026-08-30"),
    "he_n5_lambda": dict(
        value=790.007, convention="non-identity 1-norm (Ha)", q=110,
        provenance="MEASURED 2026-08-30"),
''',
     '''    # lambda: each value = retired hop-model value - E/8 exactly (E = 3, 13,
    # 34, 70 lattice edges), an analytic cross-check the build reproduces.
    "he_n2_lambda": dict(
        value=10.800, convention="non-identity 1-norm (Ha)", q=10,
        provenance=PROV_HE + " to 5e-5 rel (10.80048)"),
    "he_n3_lambda": dict(
        value=72.582, convention="non-identity 1-norm (Ha)", q=28,
        provenance=PROV_HE + " to 1.3e-5 rel (72.58242)"),
    "he_n4_lambda": dict(
        value=271.468, convention="non-identity 1-norm (Ha)", q=60,
        provenance="MEASURED 2026-09-25 debug/library_he_fix.py (271.46811)"),
    "he_n5_lambda": dict(
        value=781.257, convention="non-identity 1-norm (Ha)", q=110,
        provenance="MEASURED 2026-09-25 debug/library_he_fix.py (781.25697)"),
'''),
    # ---- QWC ----
    ('''    "he_n2_qwc": dict(value=67, convention="greedy QWC groups", q=10,
                      provenance="MEASURED 2026-08-30"),
    "he_n3_qwc": dict(value=5569, convention="greedy QWC groups", q=28,
                      provenance="MEASURED 2026-08-30"),
    "he_n4_qwc": dict(value=86224, convention="greedy QWC groups", q=60,
                      provenance="MEASURED 2026-08-30, ~62 min"),
''',
     '''    "he_n2_qwc": dict(value=63, convention="greedy QWC groups", q=10,
                      provenance="MEASURED 2026-09-25 debug/library_he_fix.py "
                                 "(physical h1; composed route 64 -- greedy is "
                                 "order-sensitive)"),
    "he_n3_qwc": dict(value=5574, convention="greedy QWC groups", q=28,
                      provenance="MEASURED 2026-09-25 (physical h1)"),
    "he_n4_qwc": dict(value=86249, convention="greedy QWC groups", q=60,
                      provenance="MEASURED 2026-09-25 (physical h1), ~115 min"),
'''),
    # ---- exponents ----
    ('''    "exp_pauli_4pt": dict(
        value=3.773, convention="log-log exponent, 4 points Q=10..110",
        provenance="MEASURED 2026-08-30; R^2 1.0000, max|log resid| 0.0064"),
    "exp_lambda_4pt": dict(
        value=1.774, convention="log-log exponent, 4 points Q=10..110",
        provenance="MEASURED 2026-08-30; R^2 0.9997, max|log resid| 0.0421"),
    "exp_qwc_3pt": dict(
        value=4.013, convention="log-log exponent, 3 points Q=10..60",
        provenance="MEASURED 2026-08-30; R^2 0.9976.  Q=110 infeasible: "
                   "2.4M terms against O(N^2) greedy grouping"),
    "exp_pauli_3pt": dict(
        value=3.779, convention="log-log exponent, 3 points Q=10..60",
        provenance="MEASURED 2026-08-30; the subset the slow tests compute"),
    "exp_lambda_3pt": dict(
        value=1.792, convention="log-log exponent, 3 points Q=10..60",
        provenance="MEASURED 2026-08-30"),
''',
     '''    "exp_pauli_4pt": dict(
        value=3.785, convention="log-log exponent, 4 points Q=10..110",
        provenance="MEASURED 2026-09-25 debug/library_he_fits.py (physical h1); "
                   "R^2 1.0000, max|log resid| 0.0125"),
    "exp_lambda_4pt": dict(
        value=1.784, convention="log-log exponent, 4 points Q=10..110",
        provenance="MEASURED 2026-09-25 (physical h1); R^2 0.9997, "
                   "max|log resid| 0.0433"),
    "exp_qwc_3pt": dict(
        value=4.049, convention="log-log exponent, 3 points Q=10..60",
        provenance="MEASURED 2026-09-25 (physical h1); R^2 0.9972.  Q=110 "
                   "infeasible: 2.4M terms against O(N^2) greedy grouping"),
    "exp_pauli_3pt": dict(
        value=3.795, convention="log-log exponent, 3 points Q=10..60",
        provenance="MEASURED 2026-09-25 (physical h1); the subset the slow "
                   "tests compute"),
    "exp_lambda_3pt": dict(
        value=1.802, convention="log-log exponent, 3 points Q=10..60",
        provenance="MEASURED 2026-09-25 (physical h1)"),
'''),
]

for old, new in EDITS:
    n = src.count(old)
    if n != 1:
        sys.exit(f"ABORT: anchor matched {n}x: {old[:70]!r}")
    src = src.replace(old, new)

# PROV_HE constant, defined once just before MEASURED dict usage: put it at
# module top after the docstring imports -- insert before the first "MEASURED = {".
anchor = "\nMEASURED = {"
if src.count(anchor) != 1:
    sys.exit("ABORT: MEASURED anchor")
src = src.replace(anchor, "\nPROV_HE = " + repr(PROV) + "\n" + anchor, 1)

# RETIRED hop-model values.  Context-gated: 287/288 are common numerals.
ret_anchor = '''    11.29:   ("he_n2_lambda",       r"lambda|1-norm",          None),'''
if src.count(ret_anchor) != 1:
    sys.exit("ABORT: RETIRED anchor")
ret_new = ret_anchor + '''
    # Retired 2026-09-25: the +1/16 graph-hop ('hybrid') h1 He model.
    11.175:  ("he_n2_lambda",       r"lambda|1-norm",          r"retired|hop|hybrid"),
    11.18:   ("he_n2_lambda",       r"lambda|1-norm",          r"retired|hop|hybrid"),
    74.207:  ("he_n3_lambda",       r"lambda|1-norm",          r"retired|hop|hybrid"),
    74.21:   ("he_n3_lambda",       r"lambda|1-norm",          r"retired|hop|hybrid"),
    275.718: ("he_n4_lambda",       r"lambda|1-norm",          r"retired|hop|hybrid"),
    275.72:  ("he_n4_lambda",       r"lambda|1-norm",          r"retired|hop|hybrid"),
    790.007: ("he_n5_lambda",       r"lambda|1-norm",          r"retired|hop|hybrid"),
    790.01:  ("he_n5_lambda",       r"lambda|1-norm",          r"retired|hop|hybrid"),
    14078:   ("he_n3_pauli",        r"Pauli|terms",            r"retired|hop|hybrid"),
    14079:   ("he_n3_pauli",        r"Pauli|terms",            r"retired|hop|hybrid"),
    250402:  ("he_n4_pauli",        r"Pauli|terms",            r"retired|hop|hybrid"),
    250403:  ("he_n4_pauli",        r"Pauli|terms",            r"retired|hop|hybrid"),
    2434441: ("he_n5_pauli",        r"Pauli|terms",            r"retired|hop|hybrid"),
    2434442: ("he_n5_pauli",        r"Pauli|terms",            r"retired|hop|hybrid"),
    86224:   ("he_n4_qwc",          r"QWC|measurement group",  r"retired|hop|hybrid"),
    5569:    ("he_n3_qwc",          r"QWC|measurement group",  r"retired|hop|hybrid"),
    3.773:   ("exp_pauli_4pt",      r"exponent|alpha|scaling", r"retired|hop|hybrid"),
    1.774:   ("exp_lambda_4pt",     r"exponent|alpha|scaling", r"retired|hop|hybrid"),
    4.013:   ("exp_qwc_3pt",        r"exponent|alpha|scaling|QWC", r"retired|hop|hybrid"),
    3.779:   ("exp_pauli_3pt",      r"exponent|alpha|scaling", r"retired|hop|hybrid"),
    1.792:   ("exp_lambda_3pt",     r"exponent|alpha|scaling", r"retired|hop|hybrid"),'''
src = src.replace(ret_anchor, ret_new, 1)

with open(P, "w", encoding="utf-8") as fh:
    fh.write(src)
print("patched OK")
