"""C17 blocked the commit, correctly.

I wrote the free-scale comparison's locked exponent as `\\|M\\|_1 ~ K^{0.72}` --
the HEADLINE form the `paper60-atomic-sublinear-exponent` family guards, whose
canonical value is 0.82 (full s+p+d+f, window K=74-164, converged domain).  But
0.72 is a different quantity on a different ladder:  the s-ONLY free-scale
resource run over K=21-136.  Writing a non-headline quantity in the headline
notation was the actual defect;  the gate caught it.

Fix: register the free-scale comparison as a matched SET (all four exponents
come off one ladder, so they must move together), and restate every locus to
name the ladder explicitly instead of borrowing the headline form.
"""
import io

REG = "debug/qa/numeric_registry.py"
PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
CHG = "CHANGELOG.md"
DON = "docs/qa/paper_60.done.md"

# ----------------------------------------------------------------- registry
reg = io.open(REG, encoding="utf-8").read()
NEW = '''    "p60_freescale_set_sonly": dict(
        value=0.72, convention="exponent: log-log slope of the LOCKED ||M||_1 vs K "
                               "on the S-ONLY free-scale comparison ladder, "
                               "K = 21..136 (n_max 6..16). NOT the headline "
                               "eq:sublinear exponent -- that is 0.82, full "
                               "s+p+d+f, window K=74..164. Different sector, "
                               "different window, different ladder",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_freescale_resource.py. This is "
                   "the FIRST member of a matched SET measured on one ladder, and "
                   "the set must move together: locked ||M||_1 K^0.716, "
                   "||H(lambda*)||_1 K^1.952, whitened ||S^-1/2 H S^-1/2||_1 "
                   "K^2.749, cond(S) K^0.936. Quoting any one against a value from "
                   "another ladder is the defect C17 blocked on 2026-09-08, when "
                   "0.72 was written in the headline ||M||_1 ~ K^p form.",
        aliases={1.95: "||H(lambda*)||_1 exponent, same ladder",
                 2.75: "whitened exponent, same ladder",
                 0.94: "cond(S) exponent, same ladder"}),
'''
anchor = '    "p60_cond_S_converged": dict('
assert anchor in reg and "p60_freescale_set_sonly" not in reg
reg = reg.replace(anchor, NEW + anchor, 1)
io.open(REG, "w", encoding="utf-8").write(reg)
print("registry: p60_freescale_set_sonly added (matched set of four)")

# -------------------------------------------------------------------- paper
t = io.open(PAP, encoding="utf-8").read()
A_OLD = r"""$K^{0.72}$ to $K^{2.75}$, $\mathrm{cond}(S)\sim K^{0.94}$)."""
A_NEW = r"""$\gvq{p60_freescale_set_sonly}{K^{0.72}}$ to $K^{2.75}$ on the $s$-only
comparison ladder, $\mathrm{cond}(S)\sim K^{0.94}$)."""
assert A_OLD in t
t = t.replace(A_OLD, A_NEW, 1)

S_OLD = r"""encoding advantage entirely:\ $\|M\|_1\sim K^{0.72}$ becomes
$\|H(\lambda^{\ast})\|_1\sim K^{1.95}$ and
$\|S^{-1/2}HS^{-1/2}\|_1\sim K^{2.75}$, with $\mathrm{cond}(S)\sim K^{0.94}$ --- a
factor $5\times10^{3}$ at $K=136$."""
S_NEW = r"""encoding advantage entirely.  On the $s$-only comparison ladder
($K=21$--$136$;\ note this is a different sector and window from the headline
exponent of Eq.~\eqref{eq:sublinear}, so the four numbers below are quoted as one
matched set) the locked $1$-norm's
$\gvq{p60_freescale_set_sonly}{K^{0.72}}$ becomes $K^{1.95}$ for
$H(\lambda^{\ast})$ and $K^{2.75}$ for the whitened $S^{-1/2}HS^{-1/2}$, with
$\mathrm{cond}(S)\sim K^{0.94}$ --- a factor $5\times10^{3}$ at $K=136$."""
assert S_OLD in t
t = t.replace(S_OLD, S_NEW, 1)
io.open(PAP, "w", encoding="utf-8").write(t)
print("paper: 2 loci restated with the ladder named")

# ---------------------------------------------------------------- synthesis
s = io.open(SYN, encoding="utf-8").read()
Y_OLD = r"""from $K^{0.72}$ to $K^{2.75}$."""
Y_NEW = r"""from $K^{0.72}$ to $K^{2.75}$ on the $s$-only comparison ladder
($K=21$--$136$, not the headline window)."""
assert Y_OLD in s
s = s.replace(Y_OLD, Y_NEW, 1)
io.open(SYN, "w", encoding="utf-8").write(s)
print("synthesis: ladder named")

# ----------------------------------------------------------- changelog/done
c = io.open(CHG, encoding="utf-8").read()
C_OLD = "The price of freeing it is the whole encoding advantage: `‖M‖₁ ~ K^0.72` becomes `‖H(λ*)‖₁ ~ K^1.95` and `‖S^{-1/2}HS^{-1/2}‖₁ ~ K^2.75`, with `cond(S) ~ K^0.94` --- a factor 5e3 at K=136."
C_NEW = "The price of freeing it is the whole encoding advantage. On the **s-only comparison ladder** (K=21..136 — a different sector and window from the headline eq:sublinear exponent, so these four are one matched set) the locked `‖M‖₁ ~ K^0.72` becomes `‖H(λ*)‖₁ ~ K^1.95` and `‖S^{-1/2}HS^{-1/2}‖₁ ~ K^2.75`, with `cond(S) ~ K^0.94` --- a factor 5e3 at K=136."
assert C_OLD in c
c = c.replace(C_OLD, C_NEW, 1)
io.open(CHG, "w", encoding="utf-8").write(c)
print("CHANGELOG: ladder named")

d = io.open(DON, encoding="utf-8").read()
D_OLD = "> `‖·‖₁` from `K^0.72` to `K^2.75`. Backing:"
D_NEW = "> `‖·‖₁` from `K^0.72` to `K^2.75` (s-only ladder, K=21..136 — not the\n> headline window). Backing:"
assert D_OLD in d
d = d.replace(D_OLD, D_NEW, 1)
io.open(DON, "w", encoding="utf-8").write(d)
print("done-record: ladder named")
