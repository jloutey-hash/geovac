r"""RECONSTRUCTION step 4: Paper 26 SS V table + prose, provenance, molecular
values and the remaining retired-scaling mentions.

Final measured I_cv (post-order-fix, bipartite 2*S_core):
  He 0.751(*)  Li 0.2276  Be 3.654e-3  B 1.711e-3  C 6.648e-4
  N/O/F ~0 (exact core closure; 1s occupation exactly 2)
(*) He carried over -- 2 electrons, no core/valence split.
Molecular: S_core 0.008, S_bond 0.330, ratio ~40, R-independent.
"""
import io

P = "papers/group6_precision_observations/paper_26_entanglement.tex"
s = io.open(P, encoding="utf-8").read()
n = 0


def rep(old, new, label):
    global s, n
    assert s.count(old) == 1, f"{label}: {s.count(old)} matches"
    s = s.replace(old, new)
    n += 1
    print(f"  ok  {label}")


# ---- provenance paragraph --------------------------------------------
rep(r"""\emph{Provenance (register tiers).}  All four properties---the
energy--entanglement decoupling, the $S\sim Z^{-2.56}$ scaling, the
$42.4\%\!\to\!99.2\%$ basis-intrinsic sparsity step (with a
$Z$-independent nonzero ERI count), and the core--valence factorization
thresholds---are MEASURED (computed) results.  The reported entropy is
the von Neumann entropy of the normalized one-body reduced density
matrix throughout.""",
    r"""\emph{Provenance (register tiers).}  The energy--entanglement
decoupling, the $S\sim Z^{-2.56}$ scaling, the
$17.1\%\!\to\!100\%$ basis-intrinsic sparsity transition (with a
$Z$-independent nonzero ERI count), the hub-migration pattern, the
core--valence factorization thresholds, the $Z^{-0.85}$ decay of the
off-diagonal energy fraction (\S\ref{sec:threelayer}), and the
composed per-block $R$-independence with its $\approx 40\times$
core/bond ratio are MEASURED (computed) results; the $Z \geq 7$
core--valence vanishing at $n_{\max} = 2$ is stronger --- a derivation
from a machine-verified premise (exact core closure).  The reported
entropy is the von Neumann entropy of the normalized one-body reduced
density matrix throughout.""",
    "provenance paragraph")

# ---- SS V table + prose ----------------------------------------------
rep(r"""He & 2 & $0.751$ \\
Li & 3 & $0.213$ \\
Be & 4 & $0.002$ \\
B--Ne & 5--10 & $< 10^{-3}$ \\""",
    r"""He & 2 & $0.751$ \\
Li & 3 & $0.228$ \\
Be & 4 & $3.7\times10^{-3}$ \\
B  & 5 & $1.7\times10^{-3}$ \\
C  & 6 & $6.6\times10^{-4}$ \\
N--F & 7--9 & $< 10^{-14}$ (exact) \\""",
    "SS V table rows")

rep(r"""The core-valence mutual information drops sharply with $Z$.  By
$Z = 4$ (beryllium), $I_{\text{cv}} \approx 0.002$, and for $Z \geq 5$
it is below $10^{-3}$.  The composed factorization's assumption
of negligible core-valence entanglement is quantitatively justified
for all atoms from boron onward.""",
    r"""The core-valence mutual information drops sharply with $Z$, decaying
monotonically from $3.7\times10^{-3}$ at beryllium through
$1.7\times10^{-3}$ (B) to $6.6\times10^{-4}$ (C), and then reaching
the floating-point floor from nitrogen onward.  The vanishing at
$Z \geq 7$ is not a small number but an \emph{exact} property of this
basis:\ every determinant in the ground multiplet carries $1s^2$, so
$\rho_{1s}$ is pure and $I_{\text{cv}} = 0$ identically (measured
$\le 3\times10^{-15}$; the $1s$ occupation is exactly $2$).

This constrains the basis rather than justifying the composed
factorization.  Exact closure says the $1s$ orbital is
\emph{unentangled in this basis at this truncation}, which is a
statement about the angular momentum eigenbasis at $n_{\max} = 2$; it
is consistent with the composed fiber-bundle ansatz but does not
establish it, since the ansatz's content is about the molecular
setting where the core and valence blocks carry different effective
charges.""",
    "SS V prose")

# ---- molecular values -------------------------------------------------
rep(r"""  S_{\text{bond}}(R) = 0.303 \text{ nats}""",
    r"""  S_{\text{bond}}(R) = 0.330 \text{ nats}""",
    "eq:rindep value")
rep(r"""($Z_{\text{eff}} = 3$) has $S_{\text{core}} = 0.006$ nats, while""",
    r"""($Z_{\text{eff}} = 3$) has $S_{\text{core}} = 0.008$ nats, while""",
    "S_core")
rep(r"""$S_{\text{bond}} = 0.303$ nats.  The ratio is""",
    r"""$S_{\text{bond}} = 0.330$ nats.  The ratio is""",
    "S_bond")

# ---- conclusions ------------------------------------------------------
rep(r"""Fourth, core-valence mutual information falls to $\approx 2\times10^{-3}$
by $Z = 4$ and below $10^{-3}$ for $Z \geq 5$, quantitatively justifying
the composed fiber bundle factorization.

These results establish that the $O(Q^{2.5})$ Pauli scaling of
the GeoVac qubit Hamiltonians is a direct consequence of the""",
    r"""Fourth, core-valence mutual information falls to
$\approx 4\times10^{-3}$ by $Z = 4$, decays monotonically through
$Z = 6$, and vanishes exactly for $Z \geq 7$ --- a property of the
basis, consistent with but not establishing the composed fiber bundle
factorization.

These results establish that the Gaunt sparsity of
the GeoVac qubit Hamiltonians is a direct consequence of the""",
    "conclusions")

io.open(P, "w", encoding="utf-8").write(s)
print(f"\n{n} edits applied")
