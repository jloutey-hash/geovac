"""Replace the fabricated Serra bibitem with the verified one.

The previous attempt routed a backslash-bearing edit through a bash heredoc,
which halves backslashes (memory: feedback_no_heredoc_backslashes).  The guard
assertion failed, the correction did not apply, and the applier then ran with
the unverified reference.  Verified reference (AMS, vol 66 no 218, Apr 1997):

  S. Serra, "Optimal, quasi-optimal and superlinear band-Toeplitz
  preconditioners for asymptotically ill-conditioned positive definite Toeplitz
  systems," Math. Comp. 66, 651 (1997).
"""
from pathlib import Path

P = Path("papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex")
src = P.read_text(encoding="utf-8")

FABRICATED = (
    r"S.~Serra, ``Preconditioning strategies for asymptotically ill-conditioned block" "\n"
    r"Toeplitz systems,'' \textit{BIT}\ \textbf{34}, 579 (1994);\ and ``On the" "\n"
    r"extreme eigenvalues of Hermitian (block) Toeplitz matrices,''" "\n"
    r"\textit{Linear Algebra Appl.}\ \textbf{270}, 109 (1998)."
)
VERIFIED = (
    r"S.~Serra, ``Optimal, quasi-optimal and superlinear band-Toeplitz" "\n"
    r"preconditioners for asymptotically ill-conditioned positive definite Toeplitz" "\n"
    r"systems,'' \textit{Math.\ Comp.}\ \textbf{66}, 651 (1997)."
)

n = src.count(FABRICATED)
assert n == 1, f"fabricated bibitem matched {n} times -- inspect before proceeding"
P.write_text(src.replace(FABRICATED, VERIFIED), encoding="utf-8")

# prove the unverified strings are gone and the verified one is present
after = P.read_text(encoding="utf-8")
for bad in ("BIT}", "579 (1994)", "Linear Algebra Appl.}\\ \\textbf{270}"):
    assert bad not in after, f"residual unverified fragment: {bad}"
assert r"\textbf{66}, 651 (1997)" in after
print("Serra bibitem corrected; no unverified fragment remains.")
