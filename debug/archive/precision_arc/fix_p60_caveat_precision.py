"""Three precision fixes, all found by verifying cited numbers against the data.

1. The bracket [6.47, 6.62] / [1.647, 1.676] is NOT a Shanks range -- it is the
   window-fit approaching from BELOW meeting Shanks descending from ABOVE. My
   prose attributed both endpoints to Shanks. Endpoints correct, attribution
   wrong.
2. The M4 caveat was quoted from the tightest tail window (0.818x, 342x) without
   saying so, and without the fact that materially weakens it: M4's own limit is
   WINDOW-UNSTABLE (0.692 -> 0.790 -> 0.818x) while M1's is stable to 3 dp
   (1.032 / 1.033 / 1.033x). That instability is itself evidence against M4, so
   the honest caveat is stronger for the paper, not weaker.
3. A retired value live in production code: geovac/sturmian_secular.py's
   docstring and docs/code_architecture.md both state "~K^0.84" flatly. That is
   the pre-2026-09-07 window fit on an unconverged radial box; the registry
   value is 0.82 and 0.84 is explicitly RETIRED.
"""
import io

PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
SEC = "geovac/sturmian_secular.py"
ARCH = "docs/code_architecture.md"

tex = io.open(PAP, encoding="utf-8").read()

OLD = r"""Across $21$ fit windows the free-floor value drifts by only
$1.3\%$ (ground) and $0.4\%$ ($2\,^{1}S$), and drifts \emph{upward} --- the fit
approaches from below --- while a model-free Shanks extrapolation descends from
above, bracketing $[6.47,6.62]$ and $[1.647,1.676]$~mHa.  Window drift is
therefore not the exposure;\ \emph{model family} is.  A two-parameter
$c+b/\ln K$ form, rejected at $340\times$ worse RMS, would place the $2\,^{1}S$
floor at $0.82\times$ chemical accuracy --- \emph{below} it.  We therefore quote
measured ladder endpoints rather than fitted asymptotes throughout, and record
that ``$2\,^{1}S$ saturates above chemical accuracy'' is a model-selection
conclusion, not a measurement."""

NEW = r"""Across $21$ fit windows the free-floor value drifts by only
$1.3\%$ (ground) and $0.4\%$ ($2\,^{1}S$), and drifts \emph{upward}:\ the fit
approaches from below.  A model-free Shanks extrapolation descends from above,
so the two together bracket the floor at $[6.47,6.62]$ and
$[1.647,1.676]$~mHa --- the lower endpoint from the fit, the upper from Shanks.
Window drift is therefore not the exposure;\ \emph{model family} is.  A
two-parameter $c+b/\ln K$ form would place the $2\,^{1}S$ floor \emph{below}
chemical accuracy.  It is rejected on two grounds rather than one:\ it fits
$96$--$342\times$ worse in RMS, and --- unlike the three-parameter form, whose
limit is stable to three decimals across every window ($1.032$, $1.033$,
$1.033\times$ chemical accuracy) --- \emph{its own limit is window-unstable},
running $0.69\to0.79\to0.82\times$ as the window tightens.  A model whose
asymptote moves with the window is not asserting an asymptote.  We nonetheless
quote measured ladder endpoints rather than fitted values throughout, and record
that ``$2\,^{1}S$ saturates above chemical accuracy'' remains a model-selection
conclusion rather than a measurement."""

assert OLD in tex, "caveat locus not found"
tex = tex.replace(OLD, NEW, 1)
io.open(PAP, "w", encoding="utf-8").write(tex)
print("paper: floor-caveat precision fixed (bracket attribution + M4 instability)")

# --- retired 0.84 live in production code and the architecture doc
sec = io.open(SEC, encoding="utf-8").read()
S_OLD = "approaching ``~K^0.84`` for the full s+p+d+f basis (Paper 60 ``eq:sublinear``);"
S_NEW = ("approaching ``~K^0.82`` over the fitted window K = 74..164 for the full\n"
         "        s+p+d+f basis (Paper 60 ``eq:sublinear``) -- a WINDOW fit, not an\n"
         "        asymptote:  the local slope falls monotonically to 0.766 by K = 514.\n"
         "        (``K^0.84`` is RETIRED -- it was measured on an unconverged 60-bohr\n"
         "        radial box whose relative error grew across the fit range.)")
assert S_OLD in sec, "sturmian_secular docstring locus not found"
sec = sec.replace(S_OLD, S_NEW, 1)
io.open(SEC, "w", encoding="utf-8").write(sec)
print("geovac/sturmian_secular.py: retired K^0.84 corrected in docstring")

arch = io.open(ARCH, encoding="utf-8").read()
A_OLD = "`gen_configs` (K^0.84, L²-divergence)"
A_NEW = "`gen_configs` (||M||_1 ~ K^0.82 on a window, L²-divergence; K^0.84 retired)"
assert A_OLD in arch, "code_architecture locus not found"
arch = arch.replace(A_OLD, A_NEW, 1)
io.open(ARCH, "w", encoding="utf-8").write(arch)
print("docs/code_architecture.md: retired K^0.84 corrected")
