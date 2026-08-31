r"""RECONSTRUCTION step 6: restore the last missing correction record --
the SS IV n_max provenance note (a prior-session correction preserved in the
pre-truncation PDF)."""
import io

P = "papers/group6_precision_observations/paper_26_entanglement.tex"
s = io.open(P, encoding="utf-8").read()

old = r"""We compute $I(i, j)$ from the FCI ground state for first-row atoms
in the graph eigenbasis at $n_{\max} = 2$."""
new = r"""We compute $I(i, j)$ from the FCI ground state for first-row atoms
in the graph eigenbasis at $n_{\max} = 2$.  (Corrected 2026-08-28:\
the values reported in this section previously came from
$n_{\max} = 3$ and $n_{\max} = 4$ runs while the text stated
$n_{\max} = 2$ --- one of them, $I(2s, 3s)$, referenced an orbital
absent from the $n_{\max} = 2$ basis.  Every number below is now the
$n_{\max} = 2$ value.  The migration pattern itself is unchanged, and
is cleaner in the corrected data.)"""
assert s.count(old) == 1, s.count(old)
io.open(P, "w", encoding="utf-8").write(s.replace(old, new))
print("  ok  SS IV n_max provenance note restored")
