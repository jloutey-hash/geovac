"""Fold the water A_1 transfer result into Paper 60, the matrix, and the register."""
from pathlib import Path

# ------------------------------------------------------------------- Paper 60
P = Path("papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex")
src = P.read_text(encoding="utf-8")

OLD = r"""\emph{Scope.} Measured on the homonuclear two-center $s$-sector symbol, for
which the parity blocks are $I\pm C$;\ whether the same construction reaches a
polyatomic block with symmetry-inequivalent centers (the
$\mathrm{H}_2\mathrm{O}$ $A_1$ case above) is untested here.  The end-to-end
resource claim additionally requires pricing the sine-transform circuit and a
block-encoding of $G$, neither of which is attempted in this paper;\ what is
established is the conditioning statement."""

NEW = r"""\textbf{[MEASURED]} \textbf{The lever reaches the polyatomic case the gerade
lever cannot}, and the structural reason is that the degeneracy's \emph{direction}
is geometry-independent.  At $\chi=\pi$ every block symbol tends to $j_0(0)=1$
whatever the separation, so for $M$ centers the matrix symbol degenerates to the
rank-one all-ones matrix and its null space has dimension $M-1$, spanned by a
fixed vector set that does not move with geometry or with basis size.  For
water's $A_1$ block, $A(\chi)=\bigl[\begin{smallmatrix}1&\sqrt2 a_{\rm OH}\\
\sqrt2 a_{\rm OH}&1+a_{\rm HH}\end{smallmatrix}\bigr]$ degenerates to
$\bigl[\begin{smallmatrix}1&\sqrt2\\\sqrt2&2\end{smallmatrix}\bigr]$, singular
with null direction $v\propto(\sqrt2,-1)$ and trace $3$.  Rotating the
two-dimensional block space by $v$ and applying $\mathrm{tri}(1,2,1)$ to that
component alone:
\begin{center}
\begin{tabular}{crr}
\toprule
$N$ & $\mathrm{cond}(A_1)$ & preconditioned\\
\midrule
$12$  & $183.0$    & $38.45$\\
$48$  & $2696.1$   & $43.62$\\
$192$ & $41699.7$  & $\mathbf{44.06}$\\
\bottomrule
\end{tabular}
\end{center}
The raw column reproduces the $N^{1.97}$ of the probe above independently
($N^{1.96}$ here);\ the preconditioned column is bounded, with increments
collapsing geometrically ($3.94,1.24,0.35,0.09$).  The constant is larger than
the diatomic $2.23$, but the growth---which is what makes the metric penalty
basis-dependent---is gone.  A control confirms the rotation is doing the work
rather than the preconditioning as such:\ the naive $\mathrm{blockdiag}(P,P)$
without it leaves the growth untouched ($2766\to42008$ over the same range).

\emph{Scope.} What is established is the \emph{conditioning} statement, on
$s$-sector shared-scale bases at $M=2$ and $M=3$.  The end-to-end resource claim
additionally requires pricing the sine-transform circuit and a block-encoding of
$G$, neither of which is attempted here;\ and by the block-structure result
above, none of this recovers $\ell$-selection, which no congruence can."""

assert src.count(OLD) == 1, "scope paragraph anchor not found"
P.write_text(src.replace(OLD, NEW), encoding="utf-8")
print("paper 60: scope paragraph replaced with the measured water transfer")

# ------------------------------------------------------------------- matrix row
M = Path("docs/claim_test_matrix.md")
s = M.read_text(encoding="utf-8")
ROW = ("| 60 | sec:resource (third lever, TRANSFER) — the breach reaches water's `A_1` block, the "
       "symmetry-inequivalent-center case where the gerade lever fails: `cond(A_1)` raw `~N^1.96` "
       "(independently reproducing the paper's `N^1.97`) vs preconditioned 38.45 -> 44.06 bounded "
       "over N=12..192. Mechanism: at `chi=pi` every block symbol -> `j0(0)=1`, so the matrix "
       "symbol is rank-one all-ones and the null space has dimension `M-1`, spanned by a FIXED "
       "geometry-independent direction | `tests/test_paper60_preconditioner.py`"
       "``::test_lever_transfers_to_water_A1_block`` + ``::test_water_needs_the_null_direction_rotation`` "
       "| tracked `geovac/sturmian_sigma_law.py` | **NEW 2026-09-12** | BACKED-SOUND. The CONTROL is "
       "the load-bearing half: naive `blockdiag(P,P)` without the rotation leaves the growth intact "
       "(2766 -> 42008), so the test excludes the reading that any preconditioner would do. "
       "Fire-tested: replacing the null-direction rotation with the identity FIRES. "
       "rests on: the chi=pi symbol limit (eq:sigma_law's own input) |")
A = "| 60 | sec:resource (third lever, LIMIT) —"
i = s.index(A)
M.write_text(s[:i] + ROW + "\n" + s[i:], encoding="utf-8")
print("claim_test_matrix: +1 transfer row")

# ---------------------------------------------------------------- walls register
W = Path("docs/walls/register.md")
s = W.read_text(encoding="utf-8")
OLDW = ("**Honest scope.** The breach is measured on the homonuclear two-center `s`-sector symbol, "
        "where the parity blocks are `I +- C`. Whether the construction reaches a polyatomic block "
        "with symmetry-inequivalent centers (water's `A_1`, `cond ~ N^1.97`) is **untested** -- and "
        "that is the case the gerade lever already fails, so it is the one that matters. Next probe.")
NEWW = ("**Scope, now measured.** The breach reaches **water's `A_1` block** -- the "
        "symmetry-inequivalent-center case where the gerade lever fails: raw `cond ~ N^1.96` "
        "(independently reproducing the paper's `N^1.97`) against a bounded `38.45 -> 44.06` over "
        "`N = 12..192`. It works because the degeneracy's DIRECTION is geometry-independent: at "
        "`chi = pi` every block symbol tends to `j0(0) = 1`, so for `M` centers the matrix symbol is "
        "the rank-one all-ones matrix and its null space has dimension `M-1`, fixed. Control: the "
        "naive `blockdiag(P,P)` without the null-direction rotation leaves the growth intact "
        "(`2766 -> 42008`), so the rotation is doing the work. Remaining scope: `s`-sector "
        "shared-scale bases at `M = 2, 3`; the end-to-end resource claim still needs the "
        "sine-transform circuit and a block-encoding of `G` priced.")
assert s.count(OLDW) == 1, "walls scope paragraph not found"
W.write_text(s.replace(OLDW, NEWW), encoding="utf-8")
print("walls register: scope updated from 'untested' to the measured transfer")
