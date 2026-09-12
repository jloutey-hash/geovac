"""The abstract advertises that BOTH documented cost risks of the genre (metric
conditioning; an outer energy search) are structurally absent.  That is true --
but only of the posing that carries the accuracy floor.  Eq.(scale_lock) says
the two facts are the same fact, so freeing the scale to reach ground-state
chemical accuracy hands BOTH risks back, not merely the 1-norm.  Say so in the
same breath, so the abstract does not contradict itself.
"""
import io

PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
tex = io.open(PAP, encoding="utf-8").read()

OLD = r"""Solving the \emph{same} span variationally with the scale
freed reaches chemical accuracy, at the price of the entire encoding advantage
($1$-norm growth from $K^{0.72}$ to $K^{2.75}$)."""

NEW = r"""Solving the \emph{same} span variationally with the scale
freed reaches chemical accuracy, but hands back \emph{both} of the cost risks
named above --- the $L^{2}$ metric returns, and the scale becomes an outer
search --- along with the encoding advantage itself ($1$-norm growth from
$K^{0.72}$ to $K^{2.75}$, $\mathrm{cond}(S)\sim K^{0.94}$).  The absence of
those two risks and the presence of the floor are therefore not independent
virtues and defects but one structural fact."""

assert OLD in tex, "abstract cost-risk locus not found"
tex = tex.replace(OLD, NEW, 1)
io.open(PAP, "w", encoding="utf-8").write(tex)
print("abstract: cost-risk consistency stated")
