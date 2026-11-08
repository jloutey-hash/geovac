"""FULL-run cert-blocker 1 (partial): register the declared-debt literals that
are TRACKED-reproducible, each with a value I MEASURED (not copied from prose).

Measured 2026-09-13 by debug/fullrun_measure_literals.py from tracked code:
  Q hydrogenic 1.1909 (paper 1.19)   Q sturmian 3.3328 (paper 3.33)   -- exact
  SW cond exponent 1.8518 on window N=12..16, R=2 (full range 1.80; asymptote N^2)
  L2 overlap exponent 1.7114 (paper printed 1.70 -> corrected to measured 1.71)
  water A_1 (raw, _water_A1) 1.9580 over N=12..192  (the "N^1.96 here" value)

NOT registered here, because the tracked helper does NOT reproduce them --
recorded as a finding, to promote-or-reclassify at the group2 review:
  * Gaussian ratio N^6 / N^1.4: tracked _gaussian_metric gives ratio-1.6 -> N^6.30
    (kappa ~1.5e6 at N=16, NOT the paper's ~1e5) and ratio-3 -> N^1.54 (NOT 1.4).
    The paper's own caveat calls these "illustrative rather than a fixed factor."
  * water A_1 N^1.97 (sec:molecular probe, 19.9->698 over N=6..36): a different
    three-center-SW construction from _water_A1; not in the tracked suite.
  * tab:resource d_inv/kappa row: [RESOURCE MODEL], O(1) cancels in ratios per
    the caption -- modelled, not a measured canonical value.

Write-first; LaTeX via raw strings. Idempotent.
"""
from __future__ import annotations

import sys

REG = "debug/qa/numeric_registry.py"
P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
DOD = "docs/qa/paper_60.done.md"

# ---- registry keys (inserted before the p60_onenorm_exponent entry) --------
REG_ANCHOR = '    "p60_onenorm_exponent": dict('
REG_NEW = '''    "p60_l2_inflation_hydrogenic": dict(
        value=1.19, convention="exponent: log-log fit of the naive-Loewdin "
                               "block-encoding 1-norm lambda vs Q=2N for the "
                               "HYDROGENIC (a=Z/n) atomic basis, N=1..5 (Q=2..10), "
                               "N=1 excluded (eq:blowup)",
        q=None,
        provenance="MEASURED 2026-09-13 via "
                   "geovac.sturmian_l2_encoding.fit_lambda_exponent(hydrogenic): "
                   "1.1909. Test-banded (1.0,1.4) in "
                   "tests/test_sturmian_l2_encoding.py::"
                   "test_paper60_l2_encoding_blowup_exponents."),
    "p60_l2_inflation_sturmian": dict(
        value=3.33, convention="exponent: same fit as p60_l2_inflation_hydrogenic "
                               "but the shared-scale STURMIAN basis; eq:blowup",
        q=None,
        provenance="MEASURED 2026-09-13 via fit_lambda_exponent(sturmian): 3.3328. "
                   "Test-banded (3.0,3.6), same test. The shared-scale set inflates "
                   "faster than hydrogenic by >1 in exponent -- the eq:blowup "
                   "mechanism (density of S^-1/2)."),
    "p60_sw_cond_exponent": dict(
        value=1.85, convention="exponent: log-log fit of cond of the momentum-space "
                               "Shibuya-Wulfman two-center metric vs N=2n, WINDOW "
                               "N=12..16 at R=2 bohr. A PRE-ASYMPTOTIC window reading, "
                               "NOT canonical: the full N=4..16 slope is 1.80 and the "
                               "DERIVED asymptote is N^2 (eq:sigma_law)",
        q=None,
        provenance="MEASURED 2026-09-13 via tests/test_paper60_sturmian.py::"
                   "_sw_metric_mom: 1.8518 on N=12..16 (full-range 1.8003). The SW "
                   "exponent is test-banded (1.5,2.2) as ~N^1.8 in "
                   "test_paper60_sw_better_conditioned_than_l2."),
    "p60_l2_overlap_exponent": dict(
        value=1.71, convention="exponent: log-log fit of cond of the shared-scale "
                               "Coulomb-Sturmian L2 overlap vs N=2,3,5,8 (the "
                               "3.0/5.83/13.93/32.16 sequence). Comparator to "
                               "p60_sw_cond_exponent",
        q=None,
        provenance="MEASURED 2026-09-13: 1.7114 on N=2,3,5,8. REGISTERED BECAUSE THE "
                   "PAPER PRINTED 1.70 (an approximate comparator); corrected to the "
                   "measured 1.71 at its one locus. Sequence pinned in "
                   "tests/test_paper60_sturmian.py::"
                   "test_paper60_l2_overlap_condition_number_grows."),
    "p60_water_a1_exponent": dict(
        value=1.96, convention="exponent: log-log fit of cond of the RAW water A_1 "
                               "block vs N=2n over N=12..192, via _water_A1. This is "
                               "the 'N^1.96 here' value of the sec:resource "
                               "preconditioner table; DISTINCT from the sec:molecular "
                               "three-center-SW probe (N^1.97, 19.9->698 over N=6..36, "
                               "a different construction)",
        q=None,
        provenance="MEASURED 2026-09-13 via tests/test_paper60_preconditioner.py::"
                   "_water_A1: 1.9580 over N=12..192 (cond 183/2696/41700). The raw "
                   "exponent is banded (1.85,2.05) in "
                   "test_uniform_banding_alone_is_not_the_lever."),
    "p60_onenorm_exponent": dict('''

# ---- paper annotations (unique-context anchors; \\gvq renders literal only) -
EDITS = [
    (P, "reg-q-hydro", r"\gvq{p60_l2_inflation_hydrogenic}",
     r"\lambda_{\text{hydrogenic}}\sim Q^{1.19},",
     r"\lambda_{\text{hydrogenic}}\sim Q^{\gvq{p60_l2_inflation_hydrogenic}{1.19}},"),

    (P, "reg-q-sturm", r"\gvq{p60_l2_inflation_sturmian}",
     r"\lambda_{\text{Sturmian}}\sim Q^{3.33},",
     r"\lambda_{\text{Sturmian}}\sim Q^{\gvq{p60_l2_inflation_sturmian}{3.33}},"),

    (P, "reg-sw", r"\gvq{p60_sw_cond_exponent}",
     r"$\mathrm{cond}(S)\sim N^{1.85}$, essentially the same rate as the $L^2$ overlap's",
     r"$\mathrm{cond}(S)\sim N^{\gvq{p60_sw_cond_exponent}{1.85}}$, essentially the same rate as the $L^2$ overlap's"),

    (P, "reg-l2", r"\gvq{p60_l2_overlap_exponent}",
     "overlap's\n$N^{1.70}$ (decisively power-law, not exponential)",
     "overlap's\n$N^{\\gvq{p60_l2_overlap_exponent}{1.71}}$ (decisively power-law, not exponential)"),

    (P, "reg-water", r"\gvq{p60_water_a1_exponent}",
     r"($N^{1.96}$ here);\ the preconditioned column is bounded",
     r"($N^{\gvq{p60_water_a1_exponent}{1.96}}$ here);\ the preconditioned column is bounded"),
]


def main() -> int:
    # registry
    with open(REG, encoding="utf-8") as fh:
        reg = fh.read()
    if '"p60_sw_cond_exponent"' in reg:
        print("  skip registry (already applied)")
    elif reg.count(REG_ANCHOR) != 1:
        print(f"  MISS registry anchor count={reg.count(REG_ANCHOR)}")
        return 3
    else:
        with open(REG, "w", encoding="utf-8") as fh:
            fh.write(reg.replace(REG_ANCHOR, REG_NEW, 1))
        print("  ok    registry: 5 keys added")

    # paper
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    applied, missed = [], []
    for _, name, marker, old, new in EDITS:
        if marker in t:
            print(f"  skip  {name} (already applied)")
            continue
        if t.count(old) != 1:
            missed.append((name, t.count(old)))
            continue
        t = t.replace(old, new)
        applied.append(name)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    for n in applied:
        print(f"  ok    {n}")
    for n, c in missed:
        print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)} annotations, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
