"""group2 Batch-4 remediation: Paper 58 F1 (claims reviewer) -- the polyatomic
scope paragraph (Sec I.A) mislabels FCI anchor figures as RHF and splices them
onto RHF-ladder endpoints as one basis-driven curve.

Verified against debug/sprint_poly3_genuine_polyatomic_memo.md:
  * +10.0% (BeH2) / +14.9% (H2O) are the L0 FCI ANCHORS (matched to composed
    FCI 11.7% / 19.4%), NOT RHF.
  * the RHF basis ladder runs 7.38 -> 1.40% (BeH2) and 8.77 -> 2.46% (H2O).
  * minimal-basis FCI is WORSE than minimal-basis RHF, and the memo (Sec 4)
    states part of the L3 residual is the method (RHF's neglected correlation),
    not the basis -- so "basis-limited rather than method-limited" is too strong.
Fix: separate the FCI-vs-FCI margin from the RHF basis ladder, label each
correctly, and soften the basis-limited claim. Takeaway (scope paragraph, angle
capability) preserved. Does not propagate (synthesis only says "open question").

Write-first; LaTeX raw strings; idempotent.
"""
from __future__ import annotations
import sys

P = "papers/group2_quantum_chemistry/paper_58_abelian_residue.tex"

OLD = r"""(BeH$_2$ and H$_2$O, basis ladder to $M = 19$): at matched minimal basis and
matched correlation treatment the present route improves the geometry only
modestly over generation~1 ($+10.0\%$ vs $11.7\%$ for BeH$_2$, $+14.9\%$ vs
$19.4\%$ for H$_2$O), but the error is basis-limited rather than
method-limited and falls to $+1.4\%$ and $+2.5\%$ respectively once the basis
is enlarged---and the H$_2$O bond angle, an observable generation~1 cannot
represent at all, comes out at $106.3^\circ$ against an experimental
$104.5^\circ$.  Those figures are RHF over Gaussian-fitted Slater shapes in
the same engine used for Sec.~\ref{sec:nah}; they are reported here as scope,
not as a result of this paper."""

NEW = r"""(BeH$_2$ and H$_2$O, basis ladder to $M = 19$).  At matched minimal basis and
matched correlation treatment (FCI on both sides) the present route improves
the geometry only modestly over generation~1---$+10.0\%$ vs $11.7\%$ for
BeH$_2$, $+14.9\%$ vs $19.4\%$ for H$_2$O.  Enlarging the basis under a common
RHF treatment then drives the $R_{\rm eq}$ error down from $+7.4\%$ to $+1.4\%$
(BeH$_2$) and $+8.8\%$ to $+2.5\%$ (H$_2$O)---so geometry is basis-limited over
this ladder, though RHF's neglect of correlation means the ladder is not at its
accuracy limit.  The H$_2$O bond angle, an observable generation~1 cannot
represent at all, comes out at $106.3^\circ$ against an experimental
$104.5^\circ$.  The ladder figures use RHF over Gaussian-fitted Slater shapes
in the same engine used for Sec.~\ref{sec:nah}; the matched-correlation anchors
are FCI.  They are reported here as scope, not as a result of this paper."""

MARKER = r"matched correlation treatment (FCI on both sides)"


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("  skip  paper58-f1 (already applied)"); return 0
    if t.count(OLD) != 1:
        print(f"  MISS  paper58-f1: anchor count={t.count(OLD)}"); return 3
    t = t.replace(OLD, NEW)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("  ok    paper58-f1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
