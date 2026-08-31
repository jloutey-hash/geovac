r"""Paper 26 + tests, FINAL values (Condon-Shortley order fix, verified
against the independent casimir_ci route at matched exponent).

  I_cv:  Li 2.276e-1  Be 3.654e-3  B 1.711e-3  C 6.648e-4  N/O/F ~0 (exact)
  1s occ: Li 1.975985  Be 1.999712
  exact core closure from N (Z=7) onward; 0 support dets lack 1s^2 at N/O/F
  degeneracies: N 4, O 7, F 6
  hub MI (returned member): N 0.1475  O 1.0625  F 0.3105
  sampled max:              N 0.6925  O 1.3395  F 1.3178
  determinant residuals: N mixed (3.8e-14..0.817), O none exact
                         (0.033..0.956), F all exact (~1e-15)
"""
import io

n = 0


def rep(path, old, new, label):
    global n
    t = io.open(path, encoding="utf-8").read()
    assert t.count(old) == 1, f"{label}: {t.count(old)} matches"
    io.open(path, "w", encoding="utf-8").write(t.replace(old, new))
    n += 1
    print(f"  ok  {label}")


P = "papers/group6_precision_observations/paper_26_entanglement.tex"

# ---- hub values (third and final revision) ----------------------------
rep(P, r"""mutual information: $I(2p_{-1}, 2p_{+1}) = 0.566$ (N),
$0.014$ (O), $0.982$ (F) \emph{for the multiplet member the eigensolver
returns} --- these ground states are 4-, 3- and 4-fold degenerate, and
the value ranges from $0$ (single-determinant members) up to a
member-dependent maximum measured at $1.43$ (N), $1.386$ (O), $1.386$ (F)""",
    r"""mutual information: $I(2p_{-1}, 2p_{+1}) = 0.148$ (N),
$1.063$ (O), $0.311$ (F) \emph{for the multiplet member the eigensolver
returns} --- these ground states are 4-, 7- and 6-fold degenerate, and
the value ranges from $0$ (single-determinant members) up to a
member-dependent maximum measured at $0.69$ (N), $1.34$ (O), $1.32$ (F)""",
    "hub values (final)")

# ---- core-closure: the exact theorem is RESTORED ----------------------
rep(P, r"""What is degeneracy-robust is the core-valence statement, though it is
now a \emph{bound} rather than the exact theorem claimed previously.
For O and F every determinant in the ground multiplet carries $1s^2$, so
$\rho_{1s}$ is pure for any member and $I_{\text{cv}} = 0$ exactly.  For
N this premise fails under the exact Coulomb selection rule:\ sixteen
support determinants lack $1s^2$, because the restored $m$-changing
multipoles include the core$\,\to\,$valence double excitation
$\langle 1s\,1s | 2p_{-1}\,2p_{+1}\rangle$.  They carry
$5.8\times10^{-5}$ of the multiplet weight ($0.0014\%$), giving a
measured $I_{\text{cv}} \le 4.3\times10^{-4}$ over sampled members ---
so core closure at N is \emph{approximate at the $10^{-4}$ level}, not
exact.""",
    r"""What is degeneracy-robust is the core-valence statement, and it holds
as a \emph{theorem}:\ for N, O and F alike, every determinant in the
ground multiplet carries $1s^2$, so $\rho_{1s}$ is pure for any member
and $I_{\text{cv}} = 0$ exactly (measured $\le 3\times10^{-15}$).
(An intermediate 2026-08-29 revision reported this as failing at N, with
sixteen support determinants lacking $1s^2$; that was an artifact of an
incompletely corrected Gaunt assembly --- the second Condon--Shortley
factor was ordered $c^k(b,d)$ rather than $c^k(d,b)$, a sign error on the
$m$-changing terms.  With the assembly verified against an independent
evaluator the premise is restored at all three.)""",
    "core-closure theorem restored")

# ---- I_cv threshold (moves again: exact from N) -----------------------
rep(P, r"""$R$-independent, with core/bond entropy ratios of $40\times$ at""",
    r"""$R$-independent, with core/bond entropy ratios of $40\times$ at""",
    "abstract ratio (unchanged, checked)")

print(f"\n{n} paper edits applied")

# ---- tests -------------------------------------------------------------
T = "tests/test_paper26_entanglement.py"
t = io.open(T, encoding="utf-8").read()
pairs = [
    # I_cv row: B/C/N move; exact closure from N
    ("assert b == pytest.approx(1.51e-3, abs=2e-4), f'B I_cv = {b:.4e}'",
     "assert b == pytest.approx(1.71e-3, abs=2e-4), f'B I_cv = {b:.4e}'"),
    ("    assert abs(occ5 - 1.999887) < 1e-5, f'B 1s occ {occ5}'",
     "    assert abs(occ5 - 1.999876) < 1e-5, f'B 1s occ {occ5}'"),
    ("    row = [atom(Z, Z)['I_core_valence'] for Z in (4, 5, 6, 7)]",
     "    row = [atom(Z, Z)['I_core_valence'] for Z in (4, 5, 6)]"),
    ("    assert all(row[k] > row[k + 1] for k in range(3)), \\\n"
     "        f'core-valence MI is not monotone decreasing across Be..N: {row}'",
     "    assert all(row[k] > row[k + 1] for k in range(2)), \\\n"
     "        f'core-valence MI is not monotone decreasing across Be..C: {row}'"),
    ("    assert row[-1] > 1e-5, \\\n"
     "        'N I_cv collapsed to the floor -- m-changing multipoles dropped again'",
     "    assert row[-1] > 1e-5, \\\n"
     "        'C I_cv collapsed to the floor -- m-changing multipoles dropped again'"),
    ("    # Exact closure survives, but only from O onward.\n    for Z in (8, 9):",
     "    # Exact closure holds from N onward.\n    for Z in (7, 8, 9):"),
    ("    assert abs(li - 0.2271) < 5e-3, f'Li I_cv = {li:.4e} != 0.227 (bipartite)'",
     "    assert abs(li - 0.2276) < 5e-3, f'Li I_cv = {li:.4e} != 0.228 (bipartite)'"),
    ("    assert abs(occ4 - 1.99972) < 1e-4, f'Be 1s occ {occ4}'",
     "    assert abs(occ4 - 1.999712) < 1e-4, f'Be 1s occ {occ4}'"),
    ("    expected_deg = {7: 4, 8: 3, 9: 4}",
     "    expected_deg = {7: 4, 8: 7, 9: 6}"),
]
for old, new in pairs:
    assert t.count(old) == 1, f"test pin: {t.count(old)} for {old[:46]!r}"
    t = t.replace(old, new)
    n += 1
io.open(T, "w", encoding="utf-8").write(t)
print(f"{len(pairs)} test pins updated")
