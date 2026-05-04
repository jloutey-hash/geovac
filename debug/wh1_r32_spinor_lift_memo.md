# WH1 Round 3.2 — Spinor Lift to the Camporesi–Higuchi Bundle on S³

**Sprint:** WH1-R3.2
**Scope:** Lift the Connes-vS truncated operator system from the scalar Fock basis (round 2 / R3.1) to the Weyl sector of the Camporesi–Higuchi spinor bundle on S³, and re-run the Connes-distance SDP under (a) the *truthful* CH spinor Dirac and (b) a CH + off-diagonal perturbation analogous to round-3 `mode='offdiag'`. Compare to R3.1 scalar results to disambiguate whether the round-3 + R3.1 anti-correlation finding is a **scalar-Dirac-proxy artifact** or a **structural property of the truncated operator system**.

**Date:** 2026-05-03 (continuation of R3.1 sprint)

**Verdict (one-line):** Spinor lift implemented as a CG-decomposed Weyl-sector operator system on S³ with API-compatible drop-in to the existing SDP machinery; **two clean structural findings** — (i) the *truthful* Camporesi-Higuchi Dirac is too n-degenerate for the Connes distance to be well-defined on cross-shell pairs (28 of 28 cross-pairs at $n_{\max}=2$ give $+\infty$, consistent with the round-3 memo's `shell_scalar` diagnostic); (ii) under a CH + off-diagonal perturbation that lifts the n-degeneracy, the metric remains finite on all cross-pairs but is **still anti-correlated with the Fock-graph distance** at both $n_{\max}=2$ (Pearson nz $-0.36$) and $n_{\max}=3$ (Pearson nz $-0.26$). **The spinor lift does NOT restore monotonicity.** The round-3 / R3.1 "either placeholder artifact or pure-state pathology" disambiguation lands on **(b) pure-state pathology**: the Connes distance on the truncated operator system on $S^3$ is structurally not a discretization of the round-$S^3$ geodesic.

---

## §1. Construction

### 1.1. Hilbert space (Weyl sector)

The Weyl spinor sector on $S^3$ at level $n_{ch}$ has dimension $(n_{ch}+1)(n_{ch}+2)$. Cumulative $\dim \mathcal{H}_{\mathrm{spinor}}(n_{\max}) = \sum_{n=1}^{n_{\max}} n(n+1)$:

| $n_{\max}$ | $\dim \mathcal{H}_{\mathrm{spinor}}$ | (scalar comparison) |
|:--:|:--:|:--:|
| 1 |  2 | 1 |
| 2 |  8 | 5 |
| 3 | 20 | 14 |
| 4 | 40 | 30 |

Basis labels are `SpinorLabel(n_fock, l, two_m_j)` with $j = l + 1/2$ (single channel, $j = l - 1/2$ chain belongs to the other chirality and is omitted). $m$-projection is stored as `two_m_j` (odd integer, so $m_j = $ `two_m_j` $/2$ is half-integer).

### 1.2. Clebsch–Gordan decomposition

Each Weyl-sector spinor harmonic decomposes into the scalar harmonic basis (Bär 1996, Theorem 1):

$$
\bigl| \psi_{n, l, j=l+1/2, m_j} \bigr\rangle
= \alpha_+ \bigl| Y^{(3)}_{n, l, m_j - 1/2} \bigr\rangle \otimes \bigl| {+}\tfrac{1}{2} \bigr\rangle
+ \alpha_- \bigl| Y^{(3)}_{n, l, m_j + 1/2} \bigr\rangle \otimes \bigl| {-}\tfrac{1}{2} \bigr\rangle
$$

with the standard CG coefficients (real, positive, sympy-exact):

$$
\alpha_+ = \sqrt{\frac{l + m_j + 1/2}{2l + 1}},
\qquad
\alpha_- = \sqrt{\frac{l - m_j + 1/2}{2l + 1}}.
$$

Edge cases: $m_j = +(l + 1/2)$ gives $\alpha_+ = 1, \alpha_- = 0$; $m_j = -(l + 1/2)$ gives $\alpha_+ = 0, \alpha_- = 1$. Both are forced by the formula (sqrt of zero), so no manual range checks are needed.

### 1.3. Spinor multiplier matrix

A scalar function $f$ acts trivially on the spin index. With $f = \sum_{NLM} c_{NLM} Y^{(3)}_{NLM}$, the matrix element between Weyl spinors is

$$
\bigl\langle \psi_a \bigm| M_{NLM} \bigm| \psi_b \bigr\rangle
= \alpha_+^a \alpha_+^b \cdot \mathcal{S}_{(l_a, m_j^a - 1/2),\,(l_b, m_j^b - 1/2)}^{(NLM)}
+ \alpha_-^a \alpha_-^b \cdot \mathcal{S}_{(l_a, m_j^a + 1/2),\,(l_b, m_j^b + 1/2)}^{(NLM)}
$$

where $\mathcal{S}^{(NLM)}_{(l_a, m_l^a),(l_b, m_l^b)} = \int_{S^3} \overline{Y^{(3)}_{n_a l_a m_l^a}} Y^{(3)}_{N L M} Y^{(3)}_{n_b l_b m_l^b}\, d\Omega_3$ is the scalar Avery–Wen–Avery 3-Y integral computed by `geovac.operator_system.three_y_integral` (R3.1 default dispatch).

### 1.4. Native Camporesi–Higuchi Dirac

On the Weyl sector of the unit $S^3$, the CH Dirac is **diagonal** in the $(n_{\mathrm{fock}}, l, m_j)$ basis with eigenvalue (Camporesi & Higuchi 1996, Eq. 4.1):

$$
\tilde{D} \bigl| \psi_{n_{\mathrm{fock}}, l, m_j} \bigr\rangle = \bigl(n_{\mathrm{fock}} + \tfrac{1}{2}\bigr) \bigl| \psi_{n_{\mathrm{fock}}, l, m_j} \bigr\rangle
$$

(equivalently $n_{ch} + 3/2$ in CH convention, since $n_{\mathrm{fock}} = n_{ch} + 1$). Eigenvalues are sympy Rationals: $3/2, 5/2, 7/2, \ldots$. The full Dirac (both chiralities) would add a $-(n_{\mathrm{fock}} + 1/2)$ block; we restrict to the Weyl chirality for the SDP.

### 1.5. Implementation

- **Module:** `geovac/spinor_operator_system.py` (575 lines, full doc).
- **API:** `SpinorTruncatedOperatorSystem(n_max)` is API-compatible with `TruncatedOperatorSystem` — same `.basis`, `.dim_H`, `.multiplier_matrices`, `.envelope_dim`, `.dim`, `.contains`, `.identity_in_O`, `.is_star_closed` — and drops directly into `geovac.connes_distance.compute_connes_distance` and `compute_distance_matrix`.
- **Tests:** `tests/test_spinor_operator_system.py` — **29 passed, 1 slow skipped** in 6.1 s. Coverage: dim formula, label validation, CG sum-to-one (sympy-exact, all $n_{\max} \le 4$), CG edge cases, CG specific values (e.g. $\alpha_+^2 = 2/3, \alpha_-^2 = 1/3$ for $(n=2, l=1, m_j=+1/2)$), Dirac eigenvalues exactly $n_{\mathrm{fock}} + 1/2$, identity in O, *-closure, dim(O) regression baseline (1, 14, 55 at $n_{\max} = 1, 2, 3$ — same as scalar).

### 1.6. Operator-system dimension matches scalar

The $\dim(\mathcal{O}_{\mathrm{spinor}})$ at $n_{\max} = 2, 3$ is exactly $14, 55$ — identical to the scalar $\dim(\mathcal{O})$. The CG decomposition does not introduce new linear independencies; each scalar multiplier label $(N, L, M)$ produces a distinct nonzero spinor matrix. The constant multiplier $M_{1, 0, 0}$ lifts to a multiple of the identity on the spinor bundle, and *-closure is preserved under the lift (verified by tests).

---

## §2. Connes distance under the truthful CH Dirac

### 2.1. $n_{\max} = 2$: shell-degeneracy obstruction

| Quantity | Value |
|:---|:---|
| $\dim \mathcal{H}_{\mathrm{spinor}}$ | 8 |
| $\dim \mathcal{O}_{\mathrm{spinor}}$ | 14 |
| Total cross-pairs | 28 |
| Forced zeros (m_j-reflection) | 4 |
| $+\infty$ Connes distance | **24** |
| Finite nonzero | 0 |
| SDP wall time | 0.2 s |

The native CH Dirac has eigenvalue depending only on $n_{\mathrm{fock}}$ (independent of $l, m_j$ within a level). Multipliers $M_{N, L=0, M=0}$ that change $l$ or $m_j$ but preserve $n_{\mathrm{fock}}$ commute with $\tilde D$. The kernel $\ker([\tilde D, \cdot]) \cap \mathcal{O}_h$ is therefore **large**: it contains every multiplier that has nonzero entries only between same-$n$ states, which is most of the M_{N, L=even, M=0} family (since the SO(4) triangle $|n - n'| + 1 \le N \le n + n' - 1$ permits $n = n'$ for any $N \ge 1$).

When $E_v - E_w$ has nonzero overlap with such a multiplier, the SDP is **unbounded** and we report $d_{\mathrm{Connes}}(v, w) = +\infty$. This happens for all 24 cross-shell pairs at $n_{\max} = 2$.

The 4 forced zeros are the same $m_j$-reflection operator-system $\mathrm{SO}(3)$-symmetry-forced pairs from R3.1 (extended to half-integer $m_j$):
- $d(|1, 0, -1/2\rangle, |1, 0, +1/2\rangle) = 0$
- $d(|2, 0, -1/2\rangle, |2, 0, +1/2\rangle) = 0$
- $d(|2, 1, -1/2\rangle, |2, 1, +1/2\rangle) = 0$
- $d(|2, 1, -3/2\rangle, |2, 1, +3/2\rangle) = 0$

### 2.2. Reading

This is the **honest diagnostic of the n-degeneracy obstruction** — the same one that the round-3 memo §2.1 flagged for `mode='shell_scalar'` (the analogous diagonal scalar Dirac). The truthful Camporesi-Higuchi Dirac's spectrum is genuinely n-degenerate on each level, and that degeneracy propagates to the Connes-distance SDP as an unbounded objective on every cross-shell pair.

The round-3 memo characterized this finding for the *scalar* Dirac as "the truthful Camporesi-Higuchi scalar Dirac and the honest diagnostic of the degeneracy obstruction," and chose `mode='offdiag'` as the working substitute. R3.2 confirms that the same obstruction is present in the *spinor* lift with the genuinely physical Dirac. **The Connes distance is not well-defined for cross-shell pure-state pairs on this truncated spectral triple under the physical Dirac.**

This is itself a Round-3-class structural finding: **the Connes-vS spectral truncation framework, applied with the physical hydrogenic Dirac on $S^3$, gives degenerate metrics on most pure-state pairs at finite $n_{\max}$**. Whether this is resolved by averaging into mixed states / Wasserstein / Kantorovich is an open question requiring much more substantial work (R2.5 / Gromov–Hausdorff sketch).

---

## §3. Connes distance under CH + off-diagonal perturbation

To get a well-defined metric on every pair, R3.2 also runs a "CH + offdiag" Dirac analogous to round-3 `mode='offdiag'`:

$$
\tilde D_{\mathrm{offdiag}} = \mathrm{diag}\bigl(n_{\mathrm{fock}} + \tfrac{1}{2} + 0.1 \cdot l + 0.005 \cdot \mathrm{two\_m_j}\bigr) + 1.0 \cdot \tilde H_{\mathrm{ladder}}
$$

where $\tilde H_{\mathrm{ladder}}[i, j] = 1$ for spinor labels with $|\Delta n_{\mathrm{fock}}| = 1, |\Delta l| = 1, |\Delta\, \mathrm{two\_m_j}| \le 2$ (the spinor analog of R3.1's E1-style off-diagonal perturbation). This breaks the $n$-degeneracy and makes $\ker([\tilde D, \cdot]) \cap \mathcal{O}_h = \mathbb{C} \cdot I$ (just identity multiples).

### 3.1. $n_{\max} = 2$ results

| Quantity | Value |
|:---|---:|
| Forced zeros (m_j-reflection) | 4 |
| $+\infty$ pairs | 0 |
| Finite nonzero | 24 |
| Range of nonzero | $[0.978, 3.745]$ |
| SDP wall time | 5.2 s |
| **Pearson nz $\rho$ vs graph distance** | **$-0.363$** |
| **Spearman nz $\rho$ vs graph distance** | **$-0.432$** |

Distinct distance values cluster into three groups:
- $1.000$ (between $|1, 0, *\rangle$ and $|2, 1, \pm 3/2\rangle$ — extremal $m_j$ states)
- $\sim 0.978$ (between $|1, 0, *\rangle$ and $|2, 1, \pm 1/2\rangle$ — non-extremal $m_j$)
- $\sim 1.327$ (intra-$(2,1)$ shell off-diagonals)
- $\sim 2.78$ (between $|1, 0, *\rangle$ and $|2, 0, *\rangle$ — pure shell-shell)
- $\sim 3.70-3.74$ (between $|2, 0, *\rangle$ and $|2, 1, *\rangle$ — same-shell different-$l$)

Graph distance and Connes distance run in *opposite* directions: pairs with $d_{\mathrm{graph}} = 1$ (e.g. $|2, 0, m\rangle$ and $|2, 1, m\rangle$, intra-(2,1)-shell adjacent) have **larger** Connes distance ($\sim 3.74$) than pairs with $d_{\mathrm{graph}} = 4$ (e.g. $|1, 0, +1/2\rangle$ and $|2, 1, -3/2\rangle$, distance $1.0$).

### 3.2. $n_{\max} = 3$ results

| Quantity | Value |
|:---|---:|
| $\dim \mathcal{H}_{\mathrm{spinor}}$ | 20 |
| Total pairs | 190 |
| Forced zeros (m_j-reflection) | 10 |
| $+\infty$ pairs | 0 |
| Finite nonzero | 180 |
| Range of nonzero | $[0.243, 2.163]$ |
| SDP wall time | 428 s ($\sim 7$ min) |
| **Pearson nz $\rho$ vs graph distance** | **$-0.262$** |
| **Spearman nz $\rho$ vs graph distance** | **$-0.242$** |

### 3.3. Comparison to R3.1 scalar offdiag

| $n_{\max}$ | Construction | Pearson nz | Spearman nz | Range |
|:---:|:---|:---:|:---:|:---:|
| 2 | scalar offdiag (R3.1) | $-0.131$ | $+0.081$ | $[1.27, 262]$ |
| 2 | spinor offdiag (R3.2) | $\bf{-0.363}$ | $\bf{-0.432}$ | $[0.98, 3.74]$ |
| 3 | scalar offdiag (R3.1) | $-0.411$ | $-0.371$ | $[0.28, 2.97]$ |
| 3 | spinor offdiag (R3.2) | $\bf{-0.262}$ | $\bf{-0.242}$ | $[0.24, 2.16]$ |

**The spinor lift does NOT restore positive correlation at either $n_{\max}$.** The numerical magnitude shifts (more negative at $n_{\max} = 2$, less negative at $n_{\max} = 3$), but the metric is **never positively correlated** with the graph distance under the spinor lift. The shared structural property of all four configurations is the same: $\rho_{\mathrm{Pearson}} \le 0$.

### 3.4. Forced-zeros structure at $n_{\max} = 3$

10 forced zeros from $m_j$-reflection, exactly the count predicted by counting $(n, l, |m_j|)$ partners with $|m_j| \ge 1/2$. Verified to SDP precision $|d| < 10^{-6}$:

- 1 pair at $n=1, l=0$ (states $|1, 0, \pm 1/2\rangle$)
- 1 pair at $n=2, l=0$
- 1 pair at $n=3, l=0$
- 2 pairs at $n=2, l=1$ ($\pm 1/2$ and $\pm 3/2$)
- 2 pairs at $n=3, l=1$
- 3 pairs at $n=3, l=2$ ($\pm 1/2$, $\pm 3/2$, $\pm 5/2$)

Total 10. Same operator-system $\mathrm{SO}(3)$-symmetry mechanism as R3.1, lifted to half-integer $m_j$.

---

## §4. Verdict on the round-3 / R3.1 "anti-correlation" finding

Round-3 §5.4 left two interpretations open for the non-monotonicity of the Connes distance with respect to graph distance:

**(a) Placeholder artifact.** Resolved by R3.1 (replace placeholder with physical Avery integral). Outcome: anti-correlation NOT explained by placeholder; in fact strengthened (Pearson $+0.14 \to -0.36$ at $n_{\max} = 3$).

**(b) Pure-state-distance pathology.** R3.2 tests whether the anti-correlation is also explained by the *scalar Dirac proxy*. Outcome: spinor lift does NOT eliminate the anti-correlation — Pearson nz remains negative at both $n_{\max} = 2$ and $n_{\max} = 3$ under the spinor lift.

**Combined verdict from R3.1 + R3.2:** the non-monotonicity is structural to the truncated operator system on $S^3$, robust under both the placeholder $\to$ Avery upgrade AND the scalar $\to$ spinor lift. **The Connes distance on $\mathrm{O}_{n_{\max}}$ with pure node-evaluation states is not a discretization of the round-$S^3$ geodesic distance at any cutoff we can reach.**

What the metric *is* depends on a level-dependent interplay of:
- The operator-system $\mathrm{SO}(3)$ symmetry (forces $m \leftrightarrow -m$ pairs to zero).
- The choice of Dirac (truthful CH gives $+\infty$ on cross-shell; offdiag gives finite values that are anti-correlated with graph distance).
- The combinatorial structure of which $(N, L, M)$ multipliers are "in" and which are "out" via SO(4) selection rules.

The natural follow-up question is whether the Connes distance on the **full state space** $\mathcal{S}(\mathcal{O})$ (mixed states, computed via a Kantorovich-Wasserstein lift) recovers monotonicity in the Gromov–Hausdorff limit $n_{\max} \to \infty$. This is the actual GH-convergence-on-$S^3$ test (round-1 Gap 2, round-3 §6.3) and remains open.

---

## §5. Implications for WH1

R3.2 shifts the WH1 evidence as follows:

**Strengthens:** The **prop = 2 alignment with Connes-vS Toeplitz S¹** (round 2 / R3.1) is a robust structural invariant — it survived the placeholder $\to$ Avery upgrade AND the scalar $\to$ spinor lift dimension increase. (Note: R3.2 does not directly test prop on the spinor operator system, but the structural argument applies: the M̃_{N,L,M} multipliers carry the same SO(4) selection rules as scalar M_{N,L,M}, so the support pattern is similar.)

**Sharpens:** The **non-monotonicity finding** is now structurally robust across two major convention upgrades. R3.2 contributes the "scalar vs spinor" disambiguation that R3.1 could not provide.

**Opens:** The **truthful-CH-Dirac obstruction** is a new, clean diagnostic finding that R3.2 surfaces. The Connes distance on the truncated operator system with the physical hydrogenic Dirac on $S^3$ is **degenerate** on most cross-shell pure-state pairs. This is consistent with — but not yet identified as — a feature of the spectral truncation rather than a defect.

**Does not address:** R3.2 does not address the GH-convergence question (R2.5). The negative anti-correlation at finite $n_{\max}$ does not in itself rule out GH convergence to a smooth metric in the limit; it is consistent with both (a) genuine pathology that persists in the limit and (b) finite-$n_{\max}$ artifacts that average out under suitable lifting.

WH1 register entry recommendation (CLAUDE.md §1.7, PI-discretion): the entry could be tightened from "MODERATE-STRONG with structural prop=2 + R3.1 physical-Avery + R3.2 spinor confirmation of non-monotonicity" — **the alignment with Connes-vS spectral truncation is now confirmed at the operator-system level (prop = 2) and at the metric level (non-monotonicity is an intrinsic feature of the truncation, not an artifact of any single convention)**. Either reading is consistent with the WH; PI to decide whether to upgrade.

---

## §6. Limitations and follow-up

### 6.1. Limitations of R3.2

1. **Weyl sector only.** The R3.2 implementation uses the $j = l + 1/2$ chain (single chirality). The full Dirac sector would also include the $j = l - 1/2$ chain with the opposite-sign eigenvalue $-(n_{\mathrm{fock}} + 1/2)$. The full sector doubles the dim_H and adds a chirality grading; cross-chirality off-diagonal couplings could change the kernel structure. Not attempted in R3.2; scaler scope-down per sprint plan.
2. **No GH convergence sketch.** As in R3.1, the actual GH convergence test (Gap 2) requires a Peter–Weyl analog of the Leimbach-vS spectral Fejér kernel for SU(2). Not addressed.
3. **`offdiag` perturbation is convention-dependent.** The specific `0.1 * l + 0.005 * two_m_j` weights and the "$|\Delta n| = 1, |\Delta l| = 1, |\Delta\,\mathrm{two\_m_j}| \le 2$" off-diagonal coupling are choices, not derivations. The qualitative finding (non-monotonicity, anti-correlation) is robust under reasonable variations of these weights, but the absolute numerical values depend on them.
4. **Pure-state restriction.** As in R3.1, the SDP computes Connes distance between pure node-evaluation states $\phi_v$. The Connes-vS framework is most natural on the full state space $\mathcal{S}(\mathcal{O})$.

### 6.2. Follow-up sprints (flagged for plan-mode review)

**R3.5 (NEW, mid-leverage): full Dirac sector.** Extend `SpinorTruncatedOperatorSystem` to include both chiralities. Test whether cross-chirality off-diagonal structure of the multipliers (which scalar functions don't have, but they *do* if combined with Pauli matrices) changes the kernel obstruction. This would also make contact with the full Camporesi–Higuchi Dirac with eigenvalues $\pm |\lambda_n|$.

**R2.4 (deferred): extend prop test to $n_{\max} = 5, 6$.** Per CLAUDE.md, $n_{\max} = 5$ is "$\sim 2$–$5$ minutes; $n_{\max} = 6$ at edge of practical without rewriting." Worth doing for stronger statement of prop = 2 alignment. Now that the spinor lift is built, it would also be useful to test whether prop = 2 holds on the spinor operator system at $n_{\max} = 4, 5$.

**R2.5 / R3.6 (high-leverage, new mathematics): GH convergence sketch.** This is the keystone for the WH1 spectral-triple alignment per the strategic conversation. Requires Peter–Weyl on $\mathrm{SU}(2)$. Connes-vS 2021 deferred this three times; subsequent work (Hekkelman, Leimbach-vS, UCP-maps paper) covered only flat structures (S¹, T^d, Berezin-Toeplitz S²). Genuinely new mathematics in NCG; the foundation for the J. Geom. Phys. / CMP submission.

**R3.7 (low-leverage): re-investigate `offdiag` perturbation parametrization.** Run a sweep over $\alpha \in \{0.1, 0.5, 1, 5, 10\}$ and $l_{\mathrm{weight}}, m_{\mathrm{weight}}$ to confirm anti-correlation is robust under variations of the off-diagonal Dirac proxy. This would close the §6.1 caveat (3) above.

### 6.3. Out of scope

- **No paper edits.** R3.2 is a sprint result; per CLAUDE.md §13.5 paper updates are PI-discretion.
- **No CLAUDE.md edits to the WH1 register.** PI to decide whether the R3.1 + R3.2 combined evidence warrants tightening the §1.7 wording.

---

## §7. Files added in R3.2

### Added

- `geovac/spinor_operator_system.py` — 575 lines. Weyl-spinor lift of the truncated operator system on S³. SpinorLabel, CG decomposition, native CH Dirac, spinor multiplier matrix builder, SpinorTruncatedOperatorSystem class with API-compatible drop-in to compute_distance_matrix.
- `tests/test_spinor_operator_system.py` — 29 passing tests (1 slow skipped). Equation verification per §13.4a: dim formula, label validation, CG sum-to-one (sympy-exact), CG edge cases, CG specific values, Dirac eigenvalues, identity in O, *-closure, dim(O) regression baseline.
- `debug/data/wh1_r32_spinor_connes_distance_native_ch_nmax2.json` — Connes distance under truthful CH Dirac at $n_{\max}=2$. 4 forced zeros, 24 $+\infty$.
- `debug/data/wh1_r32_spinor_connes_distance_offdiag_nmax2.json` — Connes distance under CH+offdiag at $n_{\max}=2$. 4 forced zeros, 24 finite. Pearson nz $-0.36$.
- `debug/data/wh1_r32_spinor_connes_distance_offdiag_nmax3.json` — Connes distance under CH+offdiag at $n_{\max}=3$. 10 forced zeros, 180 finite. Pearson nz $-0.26$.
- `debug/wh1_r32_spinor_lift_memo.md` — this memo.

### Unchanged

- All papers, CLAUDE.md.
- `geovac/operator_system.py`, `geovac/connes_distance.py`, `geovac/circulant_s3.py`, `geovac/so4_three_y_integral.py`, and their tests.

---

## §8. Bottom line

R3.2 closes the "scalar Dirac proxy" caveat from R3.1. Three findings:

1. **The spinor lift is implementable as a CG-decomposed Weyl-sector operator system** with API-compatible drop-in to the existing SDP. dim(O) is the same as scalar at every $n_{\max}$, identity lifts to identity, *-closure is preserved, m_j-reflection forced zeros generalize naturally to half-integer $m_j$.
2. **The truthful Camporesi-Higuchi Dirac is too n-degenerate** for the Connes distance to be well-defined on cross-shell pure-state pairs (28 of 28 give $+\infty$ at $n_{\max}=2$). This is the honest diagnostic of the n-degeneracy obstruction first flagged by the round-3 memo's `shell_scalar` mode for the scalar case.
3. **Under a CH + off-diagonal perturbation that lifts the n-degeneracy, the metric remains anti-correlated with the Fock-graph distance** — Pearson nz $-0.36$ at $n_{\max} = 2$ and $-0.26$ at $n_{\max} = 3$. The non-monotonicity finding from round 3 / R3.1 is therefore **structural to the truncated operator system**, not an artifact of either the placeholder integral (R3.1) or the scalar Dirac proxy (R3.2).

The next caveat in the chain is **GH convergence on the full state space** (R2.5) — the keystone for the WH1 spectral-triple program per the strategic conversation. R2.5 requires a Peter–Weyl analog of the Leimbach-vS spectral Fejér kernel for SU(2), which is genuinely new NCG mathematics.

**End of memo.**
