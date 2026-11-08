"""Two claims I wrote this sprint are contradicted by the measurements.

(a) "the ratio widens with K" -- it NARROWS, monotonically: 5.17 (K=36),
    5.11 (55), 4.66 (78), 4.29 (105).  What widens is the ABSOLUTE gap
    (2.84 -> 3.23 mHa), and that flattens by K=105.
(b) "threefold to fourfold per rung" -- the measured per-rung reductions at
    K=105 are 4.29x, 3.04x, 2.58x.  The third rung is outside the range.

Both replaced with the measured numbers, at all five loci (four in the paper,
one in the registry provenance).  The synthesis carries (b) too and is fixed in
the same pass rather than left to a later sweep.
"""
import io

PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
REG = "debug/qa/numeric_registry.py"

tex = io.open(PAP, encoding="utf-8").read()

# --- (a) abstract
A_OLD = r"""cost falling threefold to fourfold per rung up the $^{1}S$
ladder."""
A_NEW = r"""cost falling by $4.3\times$, $3.0\times$ and $2.6\times$ across the
first three rungs of the $^{1}S$ ladder."""
assert A_OLD in tex, "abstract rung locus not found"
tex = tex.replace(A_OLD, A_NEW, 1)

# --- (b) Sec.4 "ratio widens"
S_OLD = r"""ratio widens with $K$.  The floor is a ground-state pathology rather than a"""
S_NEW = r"""ratio in fact \emph{narrows} slowly with basis size --- $5.17$ at $K=36$,
$4.66$ at $78$, $4.29$ at $105$ --- while the absolute separation widens
($2.84\to3.23$~mHa) and then flattens.  The floor is a ground-state pathology
rather than a"""
assert S_OLD in tex, "Sec.4 ratio locus not found"
tex = tex.replace(S_OLD, S_NEW, 1)

# --- (c) Sec.4 per-rung
R_OLD = r"""roots, a threefold-to-fourfold reduction per rung."""
R_NEW = r"""roots --- reductions of $4.3\times$, $3.0\times$ and $2.6\times$
per rung, themselves shrinking."""
assert R_OLD in tex, "Sec.4 rung locus not found"
tex = tex.replace(R_OLD, R_NEW, 1)

# --- (d) conclusion
C_OLD = r"""$2\,^{1}S$ at identical encoding cost, falling threefold to fourfold per rung up
the $^{1}S$ ladder."""
C_NEW = r"""$2\,^{1}S$ at identical encoding cost, and falling by $4.3\times$,
$3.0\times$ and $2.6\times$ across the first three rungs of the $^{1}S$ ladder."""
assert C_OLD in tex, "conclusion rung locus not found"
tex = tex.replace(C_OLD, C_NEW, 1)

io.open(PAP, "w", encoding="utf-8").write(tex)
print("paper: 4 loci corrected")

# --- (e) synthesis
syn = io.open(SYN, encoding="utf-8").read()
Y_OLD = r"""accuracy while $2\,^{1}S$ sits $1.12\times$, the posing cost falling threefold
to fourfold per rung up the ladder."""
Y_NEW = r"""accuracy while $2\,^{1}S$ sits $1.12\times$, the posing cost falling by
$4.3\times$, $3.0\times$ and $2.6\times$ across the first three rungs."""
assert Y_OLD in syn, "synthesis rung locus not found"
syn = syn.replace(Y_OLD, Y_NEW, 1)
io.open(SYN, "w", encoding="utf-8").write(syn)
print("synthesis: rung claim corrected")

# --- (f) registry provenance
reg = io.open(REG, encoding="utf-8").read()
P_OLD = '"4.3x SMALLER than the ground state and the ratio widens with K "'
P_NEW = ('"4.3x SMALLER than the ground state. The ratio NARROWS with K "\n'
         '                   "(5.17 at K=36, 5.11 at 55, 4.66 at 78, 4.29 at 105); "\n'
         '                   "what widens is the absolute separation, 2.84 -> 3.23 mHa, "\n'
         '                   "and that flattens by K=105. RETIRED: \'the ratio widens "\n'
         '                   "with K\', written 2026-09-08 and refuted by the backing "\n'
         '                   "test the same day. "')
assert P_OLD in reg, "registry provenance locus not found"
reg = reg.replace(P_OLD, P_NEW, 1)
io.open(REG, "w", encoding="utf-8").write(reg)
print("registry: p60_posing_cost_exc provenance corrected")
