# Sprint Injection-G4 n_max extension memo — extending the level-4 cosmic-Galois injection panel to $\nmax \in \{3, 4\}$

**Date:** 2026-06-06 (same-day follow-on to the v3.80.0 injection-G4 sprint)
**Audience:** PI + next session's PM
**Verdict:** **POSITIVE-CONSISTENT.** All four compatibilities (C1, C2, C3, C4) remain bit-exact at $\nmax \in \{3, 4\}$. Residual counts grow per the closed-form prediction $R^{\mathrm{G4\text{-}Inj}}(\nmax) = 9 N(\nmax)^{2} + 6 N(\nmax) + 6$ with $N(\nmax) = \nmax(\nmax+3)/2$. Theorem~\ref{thm:injection_g4} is robust under cutoff refinement; no structural surprises emerge at higher $\nmax$. Aggregate Paper 56 verification panel grows $3{,}221 \to 5{,}864$ ($+2{,}643$ residuals). Paper 56 compiles three-pass clean (23 pages, 681 KB PDF, zero undefined refs / citations); 47/47 regression tests pass (43 fast + 4 slow).

---

## 1. What this sprint did

Extended the v3.80.0 injection-G4 verification panel from a single-cutoff probe at $\nmax = 2$ to the production grid $\nmax \in \{1, 2, 3, 4\}$, with explicit symbolic verification of each of the four compatibilities at every cutoff. The closed-form per-cutoff residual count
\[
   R^{\mathrm{G4\text{-}Inj}}(\nmax) \;=\; 9 N(\nmax)^{2} + 6 N(\nmax) + 6
\]
was derived and verified bit-exact against the explicit panel at every test cutoff.

The decomposition by compatibility is:
- **C1 (multiplicativity):** $(3 N(\nmax))^{2}$ ordered generator pairs.
- **C2 (depth-1 coproduct):** $3 N(\nmax)$ primitive generators.
- **C3 (SL$_2$-to-Levi):** $5$ (cutoff-independent; the cyclotomic character $\chi_{4}$ acts on the reductive factor of $\mathcal{G}_{4}$ and is blind to the $\nmax$-axis substrate).
- **C4 (closed-immersion):** $1 + 3 N(\nmax)$ ($1$ Gram determinant non-degeneracy + $3 N(\nmax)$ structural Glanois / EMN identities — three Hurwitz building blocks per M3 generator).

---

## 2. Per-cutoff residual counts (verified)

| $\nmax$ | $N(\nmax)$ | $3 N(\nmax)$ | C1 | C2 | C3 | C4 | Total |
|--------:|-----------:|-------------:|----:|----:|----:|----:|------:|
| 1 | 2  | 6  | 36   | 6  | 5 | 7  | 54   |
| 2 | 5  | 15 | 225  | 15 | 5 | 16 | 261  |
| 3 | 9  | 27 | 729  | 27 | 5 | 28 | 789  |
| 4 | 14 | 42 | 1764 | 42 | 5 | 43 | 1854 |

**Sum at $\nmax \in \{2, 3, 4\}$ (production grid in Paper 56 §sec:verification):** $2{,}904$.

**Aggregate Paper 56 panel transition:**
- v3.80.0 panel total: $3{,}221$ (G4-Inj at $\nmax = 2$ only contributed $261$).
- v3.81.0 panel total: $5{,}864$ (+G4-Inj at $\nmax = 3$ contributed $789$; +G4-Inj at $\nmax = 4$ contributed $1{,}854$).

---

## 3. The four compatibilities at higher cutoff — what was verified

### (C1) Multiplicativity at $\nmax \in \{3, 4\}$
For every ordered pair $(x, y)$ of primitive generators of $\HGV(\nmax) = \Sym_{\Q}(V_{\nmax})$, the symbolic-zero check $\pi(xy) - \pi(x)\,\pi(y) = 0$ holds bit-exact in the free polynomial ring $\Q[s_g : g \in \mathfrak{P}_{\nmax}]$ by the universal property of the symmetric algebra. The total panel size $(3 N(\nmax))^{2}$ grows from $225$ at $\nmax = 2$ to $729$ at $\nmax = 3$ to $1{,}764$ at $\nmax = 4$.

**Explicit verification (slow tests):** `test_c1_explicit_symbolic_panel[3]` and `test_c1_explicit_symbolic_panel[4]` pass with $729/729$ and $1{,}764/1{,}764$ symbolic-zero residuals respectively (wall time $\le 0.11$ s at $\nmax = 4$ using `sympy.simplify`).

### (C2) Depth-1 coproduct at $\nmax \in \{3, 4\}$
For every primitive generator $x \in V_{\nmax}$, the primitive coproduct $\Delta x = x \otimes 1 + 1 \otimes x$ is bit-exact by construction. The panel size $3 N(\nmax)$ grows from $15$ to $27$ to $42$ across the test grid.

### (C3) SL$_2$-to-Levi at $\nmax \in \{3, 4\}$
The bit-exact check $\det(g) = 1$ on the five-element $SL_{2}(\Q)$ panel (identity, unipotent, torus, generic, Weyl) is cutoff-independent:\ the cyclotomic character $\chi_{4}: \mathcal{G}_{4} \to \Gm$ acts on the reductive factor and does not see the $\nmax$-axis substrate at all. This is a structural observation; the panel remains exactly $5$ residuals regardless of cutoff.

### (C4) Closed-immersion at $\nmax \in \{3, 4\}$
At each cutoff, the $N(\nmax) \times N(\nmax)$ Gram matrix of M3-column period coefficients in the natural basis $\{\zeta(s, 1/4), \zeta(s, 3/4), \zeta(s, 5/4)\}$ is the identity, which has $\det = 1 \ne 0$ bit-exact. Goncharov--Deligne 2005 faithfulness on $\pi_{1}^{\mathrm{mot}}(\Gm - \mu_{4})$ promotes this Gram non-degeneracy to depth-1 injectivity at every cutoff. The $3 N(\nmax)$ structural identities are the three-Hurwitz-building-block embedding for each of the $N(\nmax)$ M3 generators.

**Explicit verification:** `test_c4_explicit_gram_panel[3]` and `test_c4_explicit_gram_panel[4]` pass with $\det = 1$ on the $9 \times 9$ and $14 \times 14$ Gram matrices respectively.

---

## 4. Decision-gate verdict

The PI's task gate names three outcomes:
- **POSITIVE-CONSISTENT** (theorem robust at higher cutoff, residuals grow per closed-form prediction): all four compatibilities verified bit-exact at $\nmax = 3$ AND $\nmax = 4$.
- **POSITIVE-AT-3-BORDERLINE-AT-4** (structural surprise at $\nmax = 4$ named in a Remark): not triggered.
- **NEGATIVE** (some compatibility fails): not triggered.

This sprint lands cleanly in **POSITIVE-CONSISTENT**.

### Why no structural surprise

The panel growth is $\nmax$-uniform by construction:
- C1 is the universal property of $\Sym_{\Q}(V_{\nmax})$, which has no $\nmax$-specific structure.
- C2 is the primitive coproduct on each generator, $\nmax$-uniform by definition.
- C3 is the cyclotomic-character / determinant condition, which decouples from the $\nmax$-axis substrate.
- C4 is the standard-basis Gram non-degeneracy + the structural Goncharov--Deligne / Glanois / EMN identities, all of which factor through the M3 column at depth 1 with no $\nmax$-specific obstruction.

This $\nmax$-uniformity is the same pattern as the TC-2c per-cutoff equality of Theorem~\ref{thm:tc2c_higher_cutoff} (the abelian-factor reconstruction at $\nmax \in \{3, 4\}$): refining the cutoff merely grows the panel size in a predictable closed form, with no new structural content per cutoff.

---

## 5. Comparison to the TC-2c per-cutoff pattern

The G4-Inj n_max extension mirrors the per-cutoff equality of TC-2c (v3.76.0):
- TC-2c achieved equality $\dim \Aut^{\otimes}(\omega) = 3 N(\nmax)$ at $\nmax \in \{1, 2, 3, 4\}$ with values $6, 15, 27, 42$ bit-exact.
- G4-Inj now achieves bit-exact verification of the closed-immersion $\Phi^{\mathrm{inj}}$ at the same grid $\nmax \in \{1, 2, 3, 4\}$.

Both patterns establish $\nmax$-uniformity at finite cutoff:\ the abelian-factor reconstruction (TC-2c) and the level-4 cosmic-Galois injection (G4-Inj) are robust under cutoff refinement, with no structural surprises emerging beyond $\nmax = 2$. The multi-year content remains the inverse-limit identification (Section~\ref{sec:open}), which neither TC-2c nor G4-Inj advance.

---

## 6. What this sprint did NOT do

- **No change to the Theorem~\ref{thm:injection_g4} statement.** The theorem text is unmodified;\ this sprint extends only the verification panel and adds a per-cutoff closed-form identity.
- **No new compatibilities.** The four C1, C2, C3, C4 compatibilities are the v3.80.0 set;\ this sprint extends their cutoff coverage, not their structural content.
- **No change to the honest scope.** Remark~\ref{rem:injection_honest_scope}'s three caveats (M1/M2 column collapse, M3 substantive content, depth-blindness) are unchanged at higher $\nmax$.
- **No advance on the equality direction.** Section~\ref{sec:equality_multi_year} still names the multi-year frontier;\ no $\nmax$-grid extension can close that gap.
- **No modification of `geovac/` production modules.** All changes are in `tests/test_paper56_injection_g4.py` (added per-cutoff tests + slow explicit panels) and the paper's verification panel.

---

## 7. Files changed

### Paper 56
- `papers/group3_foundations/paper_56_tannakian_substrate.tex`
  - Abstract:\ aggregate residual count $3{,}221 \to 5{,}864$, G4-Inj range extended to $\nmax \in \{2, 3, 4\}$.
  - `\subsection{Roadmap}` (Section~\ref{sec:roadmap}):\ aggregate updated.
  - `\subsection{Numerical verification panel}` (Section~\ref{sec:injection_panel}):\ rewritten with closed-form per-cutoff identity, new Table~\ref{tab:g4_inj_panel} listing per-cutoff residual breakdown.
  - `\section{Verification panel}` (Section~\ref{sec:verification}, Table~\ref{tab:verification}):\ two new rows added for G4-Inj at $\nmax \in \{3, 4\}$, total updated $3{,}221 \to 5{,}864$, caption updated to reference the production grid $\nmax \in \{2, 3, 4\}$.

### Regression tests
- `tests/test_paper56_injection_g4.py` (extended)
  - Existing `@pytest.mark.parametrize("n_max", [1, 2])` extended to `[1, 2, 3, 4]` on C1 and C2 symbolic tests.
  - New per-cutoff residual-count tests for C1, C2, C4 (parametrized on `(n_max, expected_count)`).
  - New per-cutoff Gram non-degeneracy test for C4 at $\nmax \in \{1, 2, 3, 4\}$.
  - New total-panel-per-cutoff test parametrized on `(n_max, expected_total)`.
  - New aggregate-panel-v3.81.0 sanity test (closed-form check against $5{,}864$).
  - New `@pytest.mark.slow` `TestExtendedCutoffSlow` class with explicit symbolic panels at $\nmax \in \{3, 4\}$ for C1 (full $1764$-pair check) and C4 (full $14 \times 14$ Gram check).
  - Total test count:\ 14 (v3.80.0) → 47 (43 fast + 4 slow).

### CHANGELOG / CLAUDE.md
- CHANGELOG v3.81.0 entry (drafted).
- CLAUDE.md §1 version bump v3.80.0 → v3.81.0.
- CLAUDE.md §2 one-liner added.

---

## 8. Verification

- Three-pass clean LaTeX compile:\ 23 pages, 681 KB PDF, zero undefined refs, zero undefined citations.
- `tests/test_paper56_injection_g4.py`:\ 47 tests pass (43 fast + 4 slow);\ 0 fail, 0 error.
- Per-cutoff explicit symbolic verification at $\nmax = 4$ (`test_c1_explicit_symbolic_panel[4]`) confirms $1{,}764 / 1{,}764$ zero residuals at wall time $\le 0.11$ s.

---

## 9. Reader's map

- The extended verification panel lives in Paper 56 §sec:injection_panel (Table~\ref{tab:g4_inj_panel}).
- The closed-form per-cutoff identity $R^{\mathrm{G4\text{-}Inj}}(\nmax) = 9 N(\nmax)^{2} + 6 N(\nmax) + 6$ is stated in §sec:injection_panel.
- The bit-exact regression tests live in `tests/test_paper56_injection_g4.py` (extended).
- The foundational sprint memo (Theorem statement + proof) is `debug/sprint_injection_g4_memo.md`.
- The multi-year forward research question (equality direction, depth $\ge 2$ exhaustion) is named in Paper 56 §sec:equality_multi_year and §sec:open_g4.
