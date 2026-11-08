r"""Discrimination test for the eight C16 entries added in the round-3 DELTA.

REGISTRY DISCRIMINATION RULE (/qa, hard): whenever you ADD or EDIT a C16 entry
you must prove it discriminates BEFORE moving on -- run the pattern against
(a) the retired wording, where it must FIRE, and (b) the corrected wording,
where it must stay SILENT.  An entry that fires on nothing is worse than no
entry at all, because the gate then reports PASS and the class is believed
guarded.  This has bitten the corpus at least three times (two dead regexes in
C17, one in C16, and C11 which could not fail at all).

Escaping LaTeX inside a Python raw string inside a regex is three nested
escape layers.  Assume it is wrong until this script says otherwise.

The RETIRED strings below are copied from git (the pre-remediation text); the
CORRECTED strings are copied from the files as they now stand.  The scanner
joins newlines to spaces before matching, so these are tested both raw and
newline-joined.
"""
from __future__ import annotations

import io
import re
import sys

CHK = "debug/qa/check_retracted_terms.py"

# (entry id, retired wording -> must FIRE, corrected wording -> must stay SILENT)
CASES = [
    ("p15-sigma-pi-decoupled",
     r"""coupling also conserves the total $M = m_1 + m_2$, the $\sigma$ ($m = 0$)
and $\pi$ ($|m| = 1$) sectors are completely decoupled in the angular
eigenvalue problem.""",
     r"""$\pi$ channels contribute a near-constant additive offset to $D_e$
across $\sigma$ $l_{\max}$ (Table~\ref{tab:extended_convergence}):"""),

    ("p15-sigma-pi-decoupled",
     r"""independent of $\sigma$ $l_{\max}$, consistent with complete
$\sigma$--$\pi$ decoupling in the angular eigenvalue problem.""",
     r"""$\pi$ channels contribute a near-constant ${\sim}6.6$~pp offset
across $\sigma$ $l_{\max}$, rising to $7.2$~pp at $l_{\max} = 6$."""),

    ("p13-hyperspherical-cusp-advantage",
     r"""a manifold in the five-dimensional hyperangular space, not a point
singularity---the crucial structural advantage of these coordinates.""",
     r"""a manifold in the five-dimensional hyperangular space, not a point
singularity---a structural feature of these coordinates."""),

    ("p13-hyperspherical-cusp-advantage",
     r"""hyperradial equation.  This is the key advantage over single-electron
coordinate systems: the cusp becomes a \emph{boundary condition on""",
     r"""coordinate.  This is a difference from single-electron coordinate
systems, not a demonstrated advantage over them."""),

    ("p13-hyperspherical-cusp-advantage",
     r"""\theta_{12})$ place the electron-electron cusp at a boundary
    condition rather than a coordinate singularity.""",
     r"""\theta_{12})$ place the electron-electron cusp on a coalescence
    manifold, where it becomes a boundary condition on the angular
    problem at fixed $R$."""),

    ("p13-sparsity-causes-005",
     r"""This sparsity is ultimately responsible for the 0.05\% accuracy
of the single-channel ($l = 0$) adiabatic approximation: the""",
     r"""This sparsity is why the single-channel ($l = 0$) adiabatic
approximation works as well as it does:\ the ground channel is only"""),

    ("p13-graph-reproduces-hydrogenic",
     r"""loutey_paper7}.  The discrete graph Laplacian reproduces
    hydrogenic energies to $< 0.1\%$~\cite{loutey_paper0,""",
     r"""a property of the construction:\ $E_0 = \kappa \lambda_{\max}$
    holds by construction, so any such figure is a bound on the
    truncation deficit of the spectrum, not an accuracy against
    experiment"""),

    ("p18-mu-rho-r-algebraic",
     r"""The adiabatic eigenvalues $\mu(R)$ and $\mu(\rho,R)$ are
    algebraic functions whose coefficients inherit transcendentals""",
     r"""The Level-3 adiabatic eigenvalues $\mu(R)$ are
    algebraic functions whose coefficients inherit transcendentals"""),

    ("p18-mu-needed-for-sub-01",
     r"""algebraic framework, but the $\mu(R)$ parameterization is needed
to achieve sub-0.1\% accuracy.""",
     r"""algebraic framework, but the $\mu(R)$ parameterization was once
thought necessary for sub-$0.1\%$ accuracy."""),

    ("p12-envelope-insensitive",
     r"""$7.6\%$ gap---is insensitive to the choice:\ every variational
point on the two one-dimensional slices quoted above exceeds $98.4\%$,""",
     r"""$\sigma$ axis, are what close the gap---is robust across the grid.  The
\emph{size} of the closure is not."""),

    ("p12-grid-floor-955",
     r"""``essentially all'' of it.  The defensible statement is that the
qualitative effect is robust while the quantitative value ranges over
$95.5$--$99.1\%$ across the grid, the independent and well-conditioned
Gaussian route giving $99.10\%$.""",
     r"""$95.5\%$ is the cell at $\alpha = 1.10$, threshold $10^{-8}$.  It is the
weakest variational value only if $\alpha$ is restricted to $\ge 1.10$---which
excludes precisely the low-$\alpha$ half this paragraph has just identified as
pathological.  The qualitative conclusion is therefore \emph{not} robust
pointwise across the grid, and this paper does not claim that it is."""),

    ("p12-grid-floor-955",
     r"""slices above exceeds $98.4\%$, but the weakest variational value
anywhere on the full $(\alpha, \mathrm{threshold})$ grid is $95.5\%$,""",
     r"""Over the full declared range $\alpha \in [0.90, 1.30]$ the
weakest variational value is $76.11\%$, at $\alpha = 0.90$ and threshold
$10^{-8}$."""),

    ("p15-delta-cbs-above-97",
     r"""    add $+0.65$~pp at $l_{\max} = 4$ but at $8.7\times$ cost.
    CBS extrapolation including $\delta$ would push the limit above 97\%.""",
     r"""    add $+0.65$~pp at $l_{\max} = 4$ but at $8.7\times$ cost.
    No $\delta$ contribution to the $\sigma{+}\pi$ gap is claimed here."""),
]


def load_patterns() -> dict:
    """Pull the live patterns straight out of the gate, not a copy."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("c16", CHK)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    out = {}
    for entry in mod.REGISTRY:
        out[entry["id"]] = entry["pattern"]
    return out


def hits(pattern: str, text: str) -> bool:
    pat = re.compile(pattern, re.IGNORECASE)
    if pat.search(text):
        return True
    return bool(pat.search(text.replace("\n", " ")))


def main() -> int:
    pats = load_patterns()
    bad = 0
    print("C16 round-3 entry discrimination -- %d cases\n" % len(CASES))
    for eid, retired, corrected in CASES:
        if eid not in pats:
            print("  [MISSING]  %s -- not in REGISTRY" % eid)
            bad += 1
            continue
        p = pats[eid]
        fired = hits(p, retired)
        silent = not hits(p, corrected)
        ok = fired and silent
        bad += 0 if ok else 1
        print("  [%s]  %s" % ("OK  " if ok else "FAIL", eid))
        print("        fires on retired wording : %s" % ("YES" if fired else "NO  <-- DEAD PATTERN"))
        print("        silent on corrected text : %s" % ("YES" if silent else "NO  <-- FALSE POSITIVE"))
    print("")
    if bad:
        print("RESULT: FAIL -- %d case(s) did not discriminate" % bad)
        return 1
    print("RESULT: PASS -- every entry fires on the retired wording and stays")
    print("        silent on the corrected wording.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
