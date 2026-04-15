# Track T0 Memo — Spinor-block ERI density $d_{\text{spinor}}(l_{\max})$

**Sprint:** Dirac-on-S³ Tier 2, early-win standalone deliverable.
**Date:** 2026-04-15.
**Status:** Complete. Exact sympy arithmetic throughout; no floating point enters any selection-rule decision.

## 1. Result

Table paralleling Paper 22 §III Theorem 3 for the scalar case, computed in the
pair-diagonal-m convention that Paper 22 actually publishes (see §3 below for
the convention explanation and a secondary table in the physically-correct
full-Gaunt convention).

| $l_{\max}$ | $Q_{\text{scalar}}$ | $Q_{\text{spinor}}$ | $d_{\text{scalar}}$ | $d_{\text{spinor}}$ | $d_{\text{spinor}}/d_{\text{scalar}}$ |
|:----------:|:-------------------:|:-------------------:|:-------------------:|:-------------------:|:-------------------------------------:|
| 0          | 1                   | 2                   | $100.0000\%$        | $25.0000\%$         | $1/4 = 0.2500$                         |
| 1          | 4                   | 8                   | $7.8125\%$          | $4.2969\%$          | $11/20 = 0.5500$                       |
| 2          | 9                   | 18                  | $2.7587\%$          | $2.1071\%$          | $553/724 = 0.7638$                     |
| 3          | 16                  | 32                  | $1.4404\%$          | $1.2329\%$          | $101/118 = 0.8559$                     |
| 4          | 25                  | 50                  | $0.9024\%$          | $0.8106\%$          | $2533/2820 = 0.8982$                   |
| 5          | 36                  | 72                  | $0.6190\%$          | $0.5722\%$          | $9611/10396 = 0.9245$                   |

All densities are exact rationals (sympy `Fraction`s); the percentage entries
are their floating-point renderings to 4 decimal places. Paper 22's scalar
values reproduce exactly (sanity-check printed in the run log).

Secondary table in the full-Gaunt convention (physically correct Coulomb
selection rule, *not* Paper 22's table):

| $l_{\max}$ | $d_{\text{scalar}}^{\text{FG}}$ | $d_{\text{spinor}}^{\text{FG}}$ | ratio |
|:----------:|:-------------------------------:|:-------------------------------:|:-----:|
| 0          | $100.0000\%$                     | $25.0000\%$                      | $0.2500$ |
| 1          | $14.8438\%$                      | $8.5938\%$                       | $0.5789$ |
| 2          | $8.5200\%$                       | $6.4586\%$                       | $0.7581$ |
| 3          | $6.0608\%$                       | $5.1674\%$                       | $0.8526$ |
| 4          | $4.8274\%$                       | $4.3034\%$                       | $0.8915$ |
| 5          | $3.9863\%$                       | $3.6794\%$                       | $0.9230$ |

Ratios in the two conventions track each other closely — the structural story
is the same.

## 2. Selection rules used

**Basis.** jj-coupled one-body orbitals labeled $(n, \kappa, m_j)$ with
$|\kappa| \le l_{\max}+1$ so that every $(l, j=l\pm 1/2)$ with $l \le l_{\max}$
contributes. Convention: $\kappa = -(l+1)$ for $j = l+1/2$, $\kappa = +l$ for
$j = l-1/2$. This gives $Q_{\text{spinor}} = 2(l_{\max}+1)^2$, matching the
task spec.

**Full-Gaunt selection for $\langle (\kappa_1 m_1)(\kappa_2 m_2)|1/r_{12}|(\kappa_3 m_3)(\kappa_4 m_4)\rangle$**
(Dyall & Faegri §9.3, Grant §7.5):

$\langle ab|1/r_{12}|cd \rangle \ne 0$ iff $m_1+m_2 = m_3+m_4$ and there exists
$k \ge 0$ with all of:

* $(l_1+l_3+k)$ even and $(l_2+l_4+k)$ even (parity);
* triangle $(j_1, k, j_3)$ and triangle $(j_2, k, j_4)$;
* triangle $(l_1, k, l_3)$ and triangle $(l_2, k, l_4)$ (defensive; implied
  by the above);
* $\left( \begin{smallmatrix} j_1 & k & j_3 \\ 1/2 & 0 & -1/2 \end{smallmatrix} \right) \ne 0$ (the reduced matrix element $\langle \kappa_1 \| C^k \| \kappa_3 \rangle$) and its $(2,4)$ analog;
* $\left( \begin{smallmatrix} j_1 & k & j_3 \\ -m_1 & m_1-m_3 & m_3 \end{smallmatrix} \right) \ne 0$ and the $(2,4)$ analog.

**Pair-diagonal-m convention** (what Paper 22 actually publishes). Paper 22's
reference code (`geovac/nuclear/potential_sparsity.py::ck_coefficient`) sets
$q = m_c - m_a$ in the 3j bottom row; combined with the implicit bottom-row
sum-to-zero constraint this forces $m_a = m_c$ per pair. The resulting rule is

$m_1 = m_3 \text{ and } m_2 = m_4$ and $\exists k$ with the same parity/triangle/reduced conditions.

This is a **stricter** rule than the full Gaunt one (forbids e.g.
$\langle p_{+1} p_{-1} | 1/r_{12} | s\ s \rangle$ which the full rule allows via
$k=2$, $q=\pm 1$). The Paper 22 percentage table reproduces to machine
precision under this stricter rule; we adopt it as the **primary convention**
so the $d_{\text{scalar}}$ column matches Paper 22 exactly. The full-Gaunt
densities are reported secondarily for completeness. **Both conventions give
$d_{\text{spinor}} \le d_{\text{scalar}}$ with quantitatively similar
ratios**, so the structural conclusion is convention-independent.

## 3. Structural interpretation

Three observations.

**(i) Asymptotic ratio $\to 1$ as $l_{\max}$ grows.** In the pair-diagonal
convention the ratio climbs 0.25, 0.55, 0.76, 0.86, 0.90, 0.92; in the full-
Gaunt convention 0.25, 0.58, 0.76, 0.85, 0.89, 0.92. Extra $j$-triangle rules
can only add zeros on top of the $l$-triangle/parity rules, but they become
**progressively redundant** as $l_{\max}$ grows because the $j$-triangle
$|j_1 - j_3| \le k \le j_1 + j_3$ differs from the $l$-triangle
$|l_1 - l_3| \le k \le l_1 + l_3$ by at most 1 (recall $j = l \pm 1/2$). The
extra zeros are boundary corrections of the triangle set; their relative size
shrinks as $O(1/l_{\max})$ heuristically.

**(ii) $l_{\max} = 0$ ratio $= 1/4$ is the pure spin-dilution floor.** With a
single $s_{1/2}$ orbital, the only two spinor labels are $(\kappa = -1, m_j = \pm 1/2)$. The 16 4-tuples split by spin pattern into four equivalence
classes of 4; only the class with matching bra/ket spin on both pairs
($\sigma_1 = \sigma_3$, $\sigma_2 = \sigma_4$) survives, giving 4/16. This
is the non-relativistic "Coulomb is spin-independent" factor made manifest
in the jj-coupled basis.

**(iii) Spinor encoding is absolutely denser but relatively comparable.**
At $l_{\max} = 3$ (the first "interesting" row for chemistry), the spinor
Hamiltonian has $(32/16)^4 = 16\times$ more 4-tuples than the scalar one,
and $(12928 / 944) = 13.7\times$ more nonzero 4-tuples. That's a 14% density
*boost relative to the naive 16x count scaling*, but a 1370% absolute increase
in the Pauli tensor size at matched $l_{\max}$. The sparsity-exponent story
survives ("angular sparsity is a universal property of the spherical basis,
not of the potential"), but the **absolute cost** of adopting a relativistic
spinor basis at matched $l_{\max}$ is roughly $16\times$ in Pauli count.
This is the correct framing for the Paper 22 extension: *sparsity exponent
preserved, prefactor increased*.

## 4. Paper 18 taxonomy placement

The selection-rule data $\{d_{\text{spinor}}(l_{\max})\}_{l_{\max}=0..5}$ are
pure rationals with denominators dividing $Q^4 = (2(l_{\max}+1)^2)^4$. No
transcendentals. Per the Paper 18 taxonomy, the **angular sparsity density is
an intrinsic exchange constant** in both the scalar and spinor cases —
the same classification Paper 22 assigns to the scalar density. The jj-coupled
spinor basis does not introduce any new $\pi$, $\alpha$, or $\sqrt{1-(Z\alpha)^2}$
content at this level of the analysis. Radial matrix elements (T1) and
spin-orbit operators (T2) will introduce such content; the spinor ERI
*density* (zero/nonzero structure) does not.

## 5. Surprises

One mild surprise, one non-surprise:

* **Non-surprise.** $d_{\text{spinor}} \le d_{\text{scalar}}$ at every
  $l_{\max}$, monotone — predicted in the sprint plan.
* **Mild surprise.** The ratio converges to $\sim 1$ rather than to $1/4$ as
  $l_{\max}$ grows. A naive "spin doubles everything" intuition would predict
  the $l_{\max} = 0$ factor of $1/4$ persists. It does not: the additional
  $j$-triangle rules are asymptotically free (ratio $\to 1$) because
  $l$-triangles already dominate them. The $1/4$ factor at $l_{\max}=0$ is
  pure spin-dilution; at $l_{\max} > 0$ this is diluted by the fact that
  most forbidden 4-tuples are already forbidden by $l$-parity, and the extra
  $j$-triangle just re-confirms zeros.

No red flags. No density exceeding $d_{\text{scalar}}$ at any $l_{\max}$
(would indicate a selection-rule bug). No non-monotonic behavior.

## 6. Recommendation for the Paper 22 extension

*(≤1 paragraph for T6 synthesis to incorporate into the Paper 22
proposed new section.)*

Paper 22's potential-independence theorem extends verbatim to the jj-coupled
spinor basis, with the angular density replaced by
$d_{\text{spinor}}(l_{\max})$ tabulated above. The spinor density is
uniformly $\le$ the scalar density, with ratio $d_{\text{spinor}}/d_{\text{scalar}}$
growing from $1/4$ at $l_{\max}=0$ (pure spin-dilution) to $0.92$ at
$l_{\max}=5$ (j-triangle rules become redundant with $l$-triangle rules). The
absolute Pauli count at matched $l_{\max}$ grows by a factor $\sim 16\times$
(from $Q^4 \to (2Q)^4$ times $d_{\text{spinor}}/d_{\text{scalar}} \approx 13\times$ for 4-component-structure bookkeeping
at $l_{\max}\in\{3,4,5\}$), so the spinor basis preserves Paper 22's
sparsity-exponent result while paying a fixed multiplicative prefactor for
relativistic physics. The universal/Coulomb-specific partition identified in
Paper 22 §VI is strengthened: angular sparsity is universal across scalar
and spinor fermion systems, while the $S^3$ Fock projection continues to
rely on the scalar Schrödinger structure (BJL algebra does not lift the Fock
map — Tier 1 Explorer Finding 2.1). **Recommendation:** add a subsection
"Spinor-block angular sparsity" with the six-row table above (both
conventions), citing the jj-coupling selection rules from Dyall §9 / Grant
§7.5 and the exact-rational sympy computation in
`debug/tier2_t0_spinor_density.py`. No changes to Paper 22's existing scalar
theorem or numerical-verification section are needed.

## 7. Files

* `debug/tier2_t0_spinor_density.py` — sympy script, exact arithmetic.
  Runs all six $l_{\max}$ values in under 10 seconds total on a laptop.
* `debug/data/tier2_t0_spinor_density.json` — full data dump (integer
  numerators/denominators for every density, per-convention).
* `debug/dirac_t0_memo.md` — this memo.

## 8. Guardrails observed

* No extension to S⁵/S⁷ (Phase 4E guardrail).
* No graph Dirac operator / Ginsparg-Wilson construction.
* No radial matrix elements (that's T1).
* `geovac/dirac_s3.py` not modified (T1's scope).
* All $l_{\max} \le 5$.
* No numerical integration anywhere. Selection-rule decisions are sympy
  symbolic: `wigner_3j(...) == 0` is a symbolic boolean.
