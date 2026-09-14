"""DELTA #2 remediation -- paper 60 + the group2 synthesis.

APPLIER BUG CLASS FIXED HERE.  `delta_fix_04_last_three.py` (and `03` before
it) accumulate edits in memory and `return 2` from inside the edit loop on a
single stale anchor -- BEFORE the write loop.  One stale anchor therefore
discards every edit in the script, silently, while the chronicle records the
remediation as applied.  That is exactly what happened: F4 and F7 were written
on 2026-09-13, never landed, and were recorded in CHANGELOG v5.11.7 as done.
Found by the DELTA #2 claim-impact reviewer.  This applier writes every edit
that DID match and reports the misses loudly with a nonzero exit.

FINDINGS APPLIED
----------------
M1 (LARGE) -- the water control sentence names a control that CANNOT support
   it.  `blockdiag(P,P)` is `I2 (x) tri(1,2,1)` and the rotation is `V (x) I`;
   they COMMUTE (verified to 3e-13 at N=192), so that control returns the
   identical spectrum in either frame and is blind to the rotation by
   construction.  The paper also mis-stated its exponent as 0.98 "halved from
   1.96" -- measured, the uniform column is N^1.950 against the raw N^1.967,
   i.e. banding alone does essentially NOTHING.  The discriminating control is
   the SELECTIVE preconditioner in the UNROTATED frame, and it is far stronger
   evidence than what was printed.  Measured (debug/p60_water_control_numbers.py,
   paper N = 2n):

       N      raw      uniform   sel/unrot   sel/rot
       24     698.8    729.2     1690        42.39
       48     2696     2766      2.03e4      43.62
       96     1.055e4  1.07e4    2.903e5     43.97
      192     4.17e4   4.201e4   4.437e6     44.06
      exponent  1.967    1.950     3.791       0.018

   So: wrong frame = 106x WORSE than untreated; right frame = flat at 44.

M2 (LARGE) -- the synthesis attaches the s-only K=105 posing-cost ladder to
   K=452, which is full spdf.  The paper disambiguates explicitly; the
   synthesis does not, and "105" appears nowhere in it.
M3 -- "Every rung is taken at the same n_max=10 family" is false for the
   leading 1s^2 rung, which is the 1x1 single configuration.
M4 -- "the largest computed basis K=452" is contradicted by the paper's own
   K=514 loci and by the registry's own alias.  K=452 is the largest basis at
   which BOTH roots were computed.
M5 -- fix #6 left "as measured below whose metric", a clause with no antecedent
   pointing the wrong way (the measurement is above, not below).
M6 -- the synthesis drops the "s-sector shared-scale" scope the paper carries
   at all three of its own loci.
M8 -- the water gain is stated without its control at the abstract, the
   conclusion and the synthesis.  Must move with M1.
N1 -- "derived" is the prohibited word for the conditioning law (PRIOR ART).

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
S = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
REG = "debug/qa/numeric_registry.py"

M1_OLD = (
    "A control confirms that it is the rotation that \\emph{bounds} the\n"
    "growth, not the preconditioning alone:\\ the naive $\\mathrm{blockdiag}(P,P)$\n"
    "without it leaves the growth \\emph{unbounded} ($2766\\to42008$ over the same\n"
    "range --- an exponent of $0.98$, halved from the raw column's $1.96$ but still\n"
    "growing, against the bounded $44$)."
)

M1_NEW = (
    "Two controls separate what is doing the work, and the first cannot do the\n"
    "job an earlier version of this paragraph gave it.  Uniform banding,\n"
    "$\\mathrm{blockdiag}(P,P)$, is $I_2\\otimes\\mathrm{tri}(1,2,1)$ while the\n"
    "rotation is $V\\otimes I$:\\ the two \\emph{commute} (to $3\\times10^{-13}$ at\n"
    "$N=192$), so that control returns the identical spectrum in the rotated and\n"
    "unrotated frames and is blind to the rotation by construction.  What it does\n"
    "establish is that banding \\emph{alone} is not the lever --- $2766\\to42008$\n"
    "over the same range, an exponent of $N^{1.95}$ against the raw column's\n"
    "$N^{1.97}$, i.e.\\ no material change.  The discriminating control is the\n"
    "\\emph{selective} preconditioner $\\mathrm{blockdiag}(P,I)$ applied in the\n"
    "\\emph{unrotated} frame, and it fails in the strong direction:\\ conditioning\n"
    "grows as $N^{3.79}$ and reaches $4.4\\times10^{6}$ at $N=192$, $106\\times$\n"
    "\\emph{worse} than leaving the metric untreated.  Alignment to the null\n"
    "direction is therefore load-bearing rather than incidental:\\ the same band is\n"
    "worse than nothing in the wrong frame and flat at $44$ in the right one."
)

EDITS = [
    (P, "M1-water-control", "blind to the rotation by construction", M1_OLD, M1_NEW),

    (P, "M3-chain-rung-scope", "beyond the leading $1s^2$ configuration",
     "Every rung is taken at the same $n_{\\max}=10$ family, so the sequence is a pure $\\ell_{\\max}$ ladder;",
     "Every rung beyond the leading $1s^2$ single configuration (a $1\\times1$ matrix) is taken at the same $n_{\\max}=10$ family, so that part of the sequence is a pure $\\ell_{\\max}$ ladder;"),

    (P, "M4-abstract-superlative", "the largest basis where both roots were computed, $K=452$, the\nground state",
     "and at the largest computed basis $K=452$ the\nground state sits",
     "and at the largest basis where both roots were computed, $K=452$, the\nground state sits"),

    (P, "M4-body-superlative", "At the largest basis where both roots were computed, $K=452$",
     "At the largest computed basis, $K=452$, the ground-state error is",
     "At the largest basis where both roots were computed, $K=452$, the ground-state error is"),

    (P, "M5-dangling-clause", "does not inherit:\\ the molecular metric is itself polynomially conditioned",
     "molecular problem does not inherit, as measured below whose metric is polynomially conditioned ($\\mathrm{cond}(S)\\sim",
     "molecular problem does not inherit:\\ the molecular metric is itself polynomially conditioned ($\\mathrm{cond}(S)\\sim"),

    (P, "N1-derived-word-a", "has an exact derivation rather than being a fit",
     "The conditioning law is in fact \\emph{derived},\nand the fitted exponent is a window reading.",
     "The conditioning law has an exact derivation rather than being a fit\n(the derivation is prior art --- see the \\textbf{[PRIOR ART]} paragraph below ---\nand what we claim is the identification), and the fitted exponent is a window\nreading."),

    (P, "N1-derived-word-b", "The water block obeys the same exact law",
     "The water block obeys the same derived law:",
     "The water block obeys the same exact law:"),

    (P, "M8-abstract-control", "provided the band is aligned to the symbol's null direction",
     "reaching water's $A_1$ block (bounded at $44$ out to $N=192$, against a raw\n$N^{1.97}$).",
     "reaching water's $A_1$ block (bounded at $44$ out to $N=192$, against a raw\n$N^{1.97}$), provided the band is aligned to the symbol's null direction:\\ in the\nunrotated frame the same band is $106\\times$ \\emph{worse} than no treatment at all."),

    (P, "M8-conclusion-control", "the same band applied in the unrotated frame is worse",
     "size without requiring equivalent centers, and it reaches water's $A_1$ block\n(Sec.~\\ref{sec:resource}).",
     "size without requiring equivalent centers, and it reaches water's $A_1$ block\n(Sec.~\\ref{sec:resource}), provided the band is aligned to the symbol's null\ndirection --- the same band applied in the unrotated frame is worse than leaving\nthe metric untreated."),

    (S, "M2-synthesis-basis-mix", "on the $s$-only ladder at $K=105$",
     "and at $K=452$ the ground state sits $4.28\\times$ above chemical\naccuracy while $2\\,^{1}S$ sits $1.08\\times$, the posing cost falling by\n$4.3\\times$, $3.0\\times$ and $2.6\\times$ across the first three rungs.",
     "and at $K=452$ --- the largest basis where both roots were computed, full\n$s\\!+\\!p\\!+\\!d\\!+\\!f$ --- the ground state sits $4.28\\times$ above chemical\naccuracy while $2\\,^{1}S$ sits $1.08\\times$;\\ and on the \\emph{$s$-only} ladder at\n$K=105$ the posing cost falls by $4.3\\times$, $3.0\\times$ and $2.6\\times$ across\nthe first three rungs.  (The two are different ladders and the paper says so;\\\nquoting one against the other is the defect this sentence was rewritten to\nremove.)"),

    (S, "M6-synthesis-scope", "Established on $s$-sector shared-scale bases at $M=2$",
     "penalty from $n^3$ to $n$.  Established for $M=2$ and non-collinear $M=3$;\\ the\ncollinear case is open.",
     "penalty from $n^3$ to $n$.  Established on $s$-sector shared-scale bases at\n$M=2$ and non-collinear $M=3$;\\ the collinear case is open."),

    (S, "M8-synthesis-control", "when the band is aligned to the symbol's null direction",
     "conditioning cost without requiring equivalent centers, reaching water's $A_1$\nblock,",
     "conditioning cost without requiring equivalent centers, reaching water's $A_1$\nblock when the band is aligned to the symbol's null direction (in the unrotated\nframe the same band is worse than no treatment at all),"),

    (REG, "M4-registry-superlative", "largest basis where both roots were computed",
     "at the largest computed basis K=452",
     "at the largest basis where both roots were computed, K=452 (the "
     "ground-state-only ladder reaches K=514)"),
]


def main() -> int:
    loaded: dict[str, str] = {}
    applied, skipped, missed = [], [], []
    for path, name, marker, old, new in EDITS:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
            skipped.append(name)
            continue
        if t.count(old) != 1:
            missed.append((name, t.count(old)))
            continue
        loaded[path] = t.replace(old, new)
        applied.append(name)
    # WRITE FIRST, report after -- a miss must never discard a match.
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    for n in applied:
        print(f"  ok    {n}")
    for n in skipped:
        print(f"  skip  {n} (already applied)")
    for n, c in missed:
        print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
