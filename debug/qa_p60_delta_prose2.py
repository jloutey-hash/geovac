"""STEP 4b: the floor, the margins, and the 'one fact' coupling.

M6 + M9 ARE ONE DEFECT, and the window data says which way to fix it.  I pulled
debug/data/p60_floor_windows_l3.json rather than trusting either report:

    ground   fit c across 20 windows   6.389 -> 6.473   (drift 1.31%, rising)
             Shanks across windows     6.818 -> 6.618   (descending)
    2^1S     fit c                     1.640 -> 1.647
             Shanks                    1.719 -> 1.676

So the bracket [6.47, 6.62] is built correctly:  lower endpoint = the SUPREMUM of
the windowed fits (the fit approaches from below), upper = the INFIMUM of Shanks
(descending from above).  The headline 6.44 is a different estimator -- the
single full-range K=74..514 fit -- dragged down by the small-K end where c is
lowest.  Since the fit approaches from below, the full-range value is strictly
the worse estimate, which is why 6.44 falls outside its own bracket.

And that is exactly why only the ground state is inconsistent:  2^1S quotes a
MEASURED ladder endpoint, the ground state quotes a FIT.  Which makes M9 -- "we
quote measured ladder endpoints rather than fitted values throughout" -- not a
stray wording slip but the paper's own stated policy, honoured for the excited
state and broken for the ground state.  One fix cures both.  The headline 6.44
is kept (it is a true statement about that fit, and the registry key means it)
but is no longer called "the floor" on its own.

M5.  The printed margins cannot be reproduced from the printed formula.  With
eps_p = eps_E/p_kappa ~ 8e-4 and the printed gaps you get 423, 52, 17, 8.  The
printed 254, 27, 8.6, 3.9 require TWO undeclared conventions:  a half-gap
separation criterion, and a per-root eps_p running 6.61e-4 -> 7.91e-4 rather than
the single ~8e-4.  Verified by direct arithmetic:  0.336/(2*1.5936e-3/2.40986)
= 254.0 exactly.  The conclusion is unchanged under either convention -- every
margin exceeds 1 through 4^1S -- so this is a reproducibility fix, not a
correction.

M3.  "Metric-free posing, sublinear 1-norm and accuracy floor are one fact, not
three" reinstates precisely the attribution the same section's own "substantive
finding" retracts:  the sublinearity belongs to the GROWTH RULE, so it cannot
also be one structural fact with the metric-free posing.  The abstract's
two-term version is correct and is the model.  What IS one fact with the posing
is the floor, and the 1-norm advantage *relative to freeing the scale*.
"""
import io

PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"


def patch(path, pairs, label):
    s = io.open(path, encoding="utf-8").read()
    for old, new in pairs:
        assert old in s, "%s: not found: %.65s" % (path, old.replace("\n", " "))
        s = s.replace(old, new, 1)
    io.open(path, "w", encoding="utf-8").write(s)
    print("  %-44s %s (%d)" % (label, path.split("/")[-1], len(pairs)))


print("M6 + M9 -- the floor sits outside its own bracket:")
patch(PAP, [
 ("$\\Delta E_\\infty=\\gvq{p60_energy_floor}{6.44}$~mHa (residual $<0.002$~mHa), with\n"
  "$6.446$~mHa still predicted at $K=10^{5}$.  Seven times the basis closes $16\\%$\n"
  "of the deficit.  The floor is $4.0\\times$ chemical accuracy\n($1.594$~mHa).",
  "$\\Delta E_\\infty=\\gvq{p60_energy_floor}{6.44}$~mHa (residual $<0.002$~mHa), with\n"
  "$6.446$~mHa still predicted at $K=10^{5}$.  Seven times the basis closes $16\\%$\n"
  "of the deficit.  That single full-range fit is dominated by its small-$K$ end,\n"
  "however, and because this fit family approaches the floor \\emph{from below} it\n"
  "is a lower estimate rather than a central one;\\ resolving it window by window\n"
  "and pairing it with a model-free extrapolation brackets the floor at\n"
  "$[6.47,6.62]$~mHa (below), and it is the bracket we quote.  Either way the\n"
  "floor is $\\approx\\!4\\times$ chemical accuracy ($1.594$~mHa)."),
 ("We nonetheless\nquote measured ladder endpoints rather than fitted values throughout, and record",
  "In the state-dependence\n"
  "comparisons below we therefore quote measured ladder endpoints rather than\n"
  "fitted values, and where a fitted floor is quoted we give its bracket rather\n"
  "than a single value.  We record"),
], "headline re-scoped + policy sentence made true")

print("\nM5 -- the margins need their two conventions declared:")
patch(PAP, [
 ("$\\epsilon_p=\\epsilon_E/p_\\kappa\\approx8\\times10^{-4}$ against gaps of $0.336$,\n"
  "$0.041$, $0.014$ and $0.006$ --- margins of $254$, $27$, $8.6$ and $3.9$ --- so",
  "$\\epsilon_p=\\epsilon_E/p_\\kappa(k)$, which runs $6.6$--$7.9\\times10^{-4}$ across\n"
  "the first four roots, against $p_\\kappa$ gaps of $0.336$, $0.041$, $0.014$ and\n"
  "$0.006$;\\ requiring half a gap of separation gives margins of $254$, $27$,\n"
  "$8.6$ and $3.9$ --- so"),
], "per-root eps_p + half-gap criterion")

print("\nM3 -- the 'one fact' coupling contradicts the section's own finding:")
patch(PAP, [
 ("\\emph{Metric-free posing, sublinear $1$-norm\nand accuracy floor are one fact, not three.}",
  "\\emph{Metric-free posing and accuracy floor are\n"
  "one fact, not two --- both are the scale lock, and the $1$-norm advantage\n"
  "\\emph{over the free-scale alternative} is its third consequence.  The\n"
  "sublinearity itself is not:\\ that belongs to the growth rule.}"),
], "paper")
patch(SYN, [
 ("Metric-free posing, sublinear $1$-norm and\naccuracy floor are one fact rather than three.",
  "Metric-free posing and accuracy floor are one\n"
  "fact rather than two, both being the scale lock;\\ the sublinearity is not part\n"
  "of that identity, since it belongs to the basis-growth rule."),
], "synthesis")
print("\ndone.")
