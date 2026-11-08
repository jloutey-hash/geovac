# Sprint Q5'-ProSystem-Lockdown PS-2 — $U^*$-action compatibility with PS-1 transitions

**Date:** 2026-06-06 (single-thread sprint, second of four sub-tracks)
**Sprint:** PS-2 of Q5'-ProSystem-Lockdown (PS-1 closed v3.67.0; PS-3, PS-4 to follow)
**Driver:** `debug/compute_q5p_ps2_ustar_compatibility.py`
**Module:** `geovac/pro_system.py` (PS-2 additions ~280 lines; total file ~610 lines)
**Data:** `debug/data/sprint_q5p_ps2_ustar_compatibility.json`
**Wall time:** 0.21 s
**Discipline:** bit-exact `sympy.Rational` / `sympy.Integer` throughout. No floats. No PSLQ. No transcendentals introduced.

---

## 1. TL;DR

**Verdict: POSITIVE.** The Levi-decomposed cosmic-Galois group $U^* = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ (v3.63.0 L1) acts compatibly with PS-1's closed-form transition maps across the entire $n_{\max} \le 5$ panel.

**485 / 485 bit-exact zero residuals:**

- **10 Hopf-hom cofiltered axiom identities** ($\Phi_{m, k} = \Phi_{n, k} \cdot \Phi_{m, n}$ where $\Phi$ is the block-diagonal three-Mellin-slot lift of PS-1's $P_{m, k}$, verified for all triples in $\{1, \ldots, 5\}$).
- **435 $\mathbb{G}_a$ generator compatibility checks** (every translation generator of $\mathbb{G}_a^{3 N(m)}$ at each pair $(m, k) \le 5$ either survives cleanly to its $(m \to k)$ image when $n \le k$, or drops to zero when $n > k$; total generator count summed over the 10 pairs: $15 + 27 + 27 + 42 + 42 + 42 + 60 + 60 + 60 + 60 = 435$).
- **40 class-level $U^*$-action identities** (Interpretation C closure for $\chi$ and $\eta$ across all 10 pairs $\times$ 2 characters $\times$ 2 boolean conditions per pair:\ $P_{m, k}(U^* \cdot \psi^{(m)}) = U^* \cdot P_{m, k}(\psi^{(m)}) = \psi^{(k)}$).

**$SL_2$ commutativity is categorical, not per-cell:**\ the $SL_2$ factor acts on the Peter--Weyl decoration which lives on the $j_{\max}$ axis, independent of the $n_{\max}$ axis on which $P_{m, k}$ acts. Acting first and then truncating equals truncating first and then acting by independence of axes;\ no cell-level check is needed.

**$U^*$ lifts to the PS-1 pro-system by construction.** PS-2 closes the second sub-track of the Pro-System-Lockdown sprint;\ PS-3 (inverse limit $\mathcal{O}_\infty$ and extended $U^*$-action on the limit object) can proceed on this substrate.

---

## 2. Verdict against decision gate

| Gate | Verdict |
|:-----|:--------|
| POSITIVE | **selected** — bit-exact cofiltered axiom at the Hopf-hom level (10 triples), bit-exact $\mathbb{G}_a$ generator compatibility (435 generator checks), bit-exact class-level $U^*$-action compatibility on $\chi, \eta$ (40 identities), categorical $SL_2$ commutativity. Total:\ 485 zero residuals out of 485. |
| BORDERLINE | not selected — closure is bit-exact at the full panel. |
| STOP | rejected — every identity verified bit-exact;\ the falsifier extension to $n_{\max} = 5$ reproduces the compatibility without modification. |

---

## 3. What PS-2 adds beyond PS-1 and v3.66.0 FO3

PS-1 (v3.67.0) gave the closed-form algebra-level transitions $P_{m, k}: \mathcal{O}_m \to \mathcal{O}_k$ and the cofiltered axiom $P_{m, k} = P_{n, k} \cdot P_{m, n}$ bit-exact across the $n_{\max} \le 5$ panel. v3.63.0 L1 gave the Levi-decomposed $U^* = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ acting on the candidate Hopf algebra $\mathcal{H}_{\mathrm{GV}}(n_{\max})$. v3.66.0 FO3 closed Interpretation C of the $U^*$-action on $\chi, \eta, F(s)$ at the period-pairing level.

PS-2 connects these three results with three structural ingredients:

1. **Hopf-hom lift $\Phi_{m, k}$ of PS-1's $P_{m, k}$.** Defined as the algebra homomorphism $\mathcal{H}_{\mathrm{GV}}(m) \to \mathcal{H}_{\mathrm{GV}}(k)$ on the primitive generators $\{x_{(n, l), s}^{(m)}\}_{(n, l), s}$ (with $s \in \{0, 1, 2\}$ indexing the master Mellin engine slots M1, M3, M2) by

$$
\Phi_{m, k}(x_{(n, l), s}^{(m)}) = \begin{cases} x_{(n, l), s}^{(k)} & \text{if } n \le k, \\ 0 & \text{if } n > k. \end{cases}
$$

   Realised as a $3 N(k) \times 3 N(m)$ block-diagonal integer matrix $\mathrm{diag}(P_{m, k}, P_{m, k}, P_{m, k})$ — three independent copies of PS-1's transition, one per Mellin slot. Class:\ `HopfTransition(n_high, n_low)` in `geovac/pro_system.py`. Hopf-axiom compatibility (coproduct, counit, antipode) is structurally trivial on the abelian primitive substrate established in v3.61.0 Track A;\ this is the reason the substantive PS-2 content lives at the cofiltered-axiom and class-action levels rather than at the Hopf-axiom level.

2. **Cofiltered axiom at the Hopf-hom level.** For every triple $1 \le k < n < m \le 5$, the matrix identity $\Phi_{m, k} = \Phi_{n, k} \cdot \Phi_{m, n}$ holds bit-exact (10 zero residuals at matrix dimensions $3 N(k) \times 3 N(m)$). Reduces to PS-1's algebra-level cofiltered axiom by the block-diagonal structure, verified independently here.

3. **$U^*$-action commutes with $\Phi$ at the generator level.** The $\mathbb{G}_a^{3 N(m)}$ factor has one translation generator $e_{(n, l), s}^{(m)}$ per primitive generator $x_{(n, l), s}^{(m)}$ of $\mathcal{H}_{\mathrm{GV}}(m)$. For every pair $(m, k) \le 5$ and every generator $e_g^{(m)}$, the compatibility condition is

$$
\Phi_{m, k}\!\left(e_g^{(m)} \cdot \psi\right) = \begin{cases} e_g^{(k)} \cdot \Phi_{m, k}(\psi) & \text{if } n \le k, \\ 0 & \text{if } n > k, \end{cases}
$$

   bit-exact for every $\psi$ in $\mathcal{H}_{\mathrm{GV}}(m)$. Verified per-generator across all 10 pairs:\ $3 N(m)$ generators per cutoff $m$, summed across the 10 pairs gives 435 generator-level identities. At each pair $(m, k)$, exactly $3 N(k)$ generators survive and $3 (N(m) - N(k))$ generators are killed, with the survival pattern matching the PS-1 sector-locality structurally.

4. **$SL_2$ commutativity is a categorical statement.** The $SL_2$ factor of $U^* = \mathbb{G}_a \rtimes SL_2$ acts on the Peter--Weyl decoration of v3.63.0 L2, which is parameterized by the $j_{\max}$ axis. The transition $P_{m, k}$ acts on the $n_{\max}$ axis. The two axes are independent, so $SL_2 \circ P_{m, k} = P_{m, k} \circ SL_2$ by independence-of-axes;\ this is not something requiring a per-cell bit-exact check and is recorded as a categorical bit-exact statement.

5. **Class-level $U^*$-action on $\chi$ and $\eta$.** The Interpretation C closure (v3.66.0 FO3) established that $U^*$ acts trivially on $\chi$ and $\eta$ (depth-0 integer-valued cocycle classes). PS-2 verifies that this triviality survives the transition:\ both $P_{m, k}(U^* \cdot \psi^{(m)})$ and $U^* \cdot P_{m, k}(\psi^{(m)})$ reduce to $P_{m, k}(\psi^{(m)})$, which equals $\psi^{(k)}$ by PS-1's pull-back identity. Verified bit-exact for both $\chi$ and $\eta$ at all 10 pairs (20 identities $\times$ 2 boolean conditions per pair $=$ 40 zero residuals).

---

## 4. The block-diagonal Hopf transition $\Phi_{m, k}$

### 4.1 Substrate

The candidate Stage-2 Hopf algebra (v3.61.0 Track A) is

$$
\mathcal{H}_{\mathrm{GV}}(n_{\max}) = \mathrm{Sym}_{\mathbb{Q}}(V_{n_{\max}}),
\qquad
V_{n_{\max}} = \bigoplus_{(n, l), s} \mathbb{Q} \cdot x_{(n, l), s}^{(n_{\max})},
$$

with $(n, l)$ ranging over $\mathrm{sectors}(n_{\max})$ and $s \in \{0, 1, 2\}$ indexing the master Mellin engine slots. Dimension of $V_{n_{\max}}$:\ $3 N(n_{\max})$.

Per master Mellin engine (Paper 18 \S III.7):

| $s$ | Mechanism | Transcendental content |
|:---:|:----------|:------|
| 0 | M1 (Hopf-base measure) | $\pi$ at $s = d/2$ pole |
| 1 | M3 (vertex-parity Hurwitz) | $\zeta(3), \zeta(5), \ldots$ at quarter-integer shifts |
| 2 | M2 (Seeley--DeWitt) | $\pi^{2k}$ at integer $s$ |

### 4.2 Closed form

$$
\boxed{
\Phi_{m, k}(x_{(n, l), s}^{(m)}) = \begin{cases} x_{(n, l), s}^{(k)} & \text{if } n \le k, \\ 0 & \text{if } n > k. \end{cases}
}
$$

As a $3 N(k) \times 3 N(m)$ matrix:\ $\mathrm{diag}(P_{m, k}, P_{m, k}, P_{m, k})$. Hopf-axiom compatibility:

- **Coproduct:** $\Delta(\Phi(x)) = \Phi(x) \otimes 1 + 1 \otimes \Phi(x) = (\Phi \otimes \Phi)(\Delta(x))$ for every primitive $x$.
- **Counit:** $\varepsilon(\Phi(x)) = 0 = \Phi(\varepsilon(x))$ trivially on $V$.
- **Antipode:** $\Phi(S(x)) = \Phi(-x) = -\Phi(x) = S(\Phi(x))$.

Structural reason for the trivial axiom verification:\ $\mathcal{H}_{\mathrm{GV}}$ is abelian primitive (v3.61.0 Track A, forced by the sector-disjointness of Camporesi--Higuchi shell idempotents). The coproduct on $\mathrm{Sym}(V)$ is determined by $V$ and is uniformly primitive across $n_{\max}$.

### 4.3 Cofiltered axiom panel

| Triple $(m, n, k)$ | $\dim(\Phi)$ | Residual$^2$ | Verdict |
|:------------------:|:---------------------:|:------------:|:-------:|
| $(3, 2, 1)$ | $6 \times 27$  | 0 | OK |
| $(4, 2, 1)$ | $6 \times 42$  | 0 | OK |
| $(4, 3, 1)$ | $6 \times 42$  | 0 | OK |
| $(4, 3, 2)$ | $15 \times 42$ | 0 | OK |
| $(5, 2, 1)$ | $6 \times 60$  | 0 | OK |
| $(5, 3, 1)$ | $6 \times 60$  | 0 | OK |
| $(5, 3, 2)$ | $15 \times 60$ | 0 | OK |
| $(5, 4, 1)$ | $6 \times 60$  | 0 | OK |
| $(5, 4, 2)$ | $15 \times 60$ | 0 | OK |
| $(5, 4, 3)$ | $27 \times 60$ | 0 | OK |

10 / 10 bit-exact zero residuals. The Hopf-hom lift respects composition under PS-1's transitions on every triple.

---

## 5. $\mathbb{G}_a$ generator compatibility

### 5.1 Generators and survival pattern

At cutoff $m$, the pro-unipotent factor $\mathbb{G}_a^{3 N(m)}$ has $3 N(m)$ translation generators $e_{(n, l), s}^{(m)}$ for $(n, l) \in \mathrm{sectors}(m)$ and $s \in \{0, 1, 2\}$. Counts:

| $m$ | $3 N(m)$ |
|:---:|:--------:|
| 1 | 6   |
| 2 | 15  |
| 3 | 27  |
| 4 | 42  |
| 5 | 60  |

For each pair $(m, k)$ with $k < m$:
- Generators with $n \le k$ **survive**:\ image is $e_{(n, l), s}^{(k)}$. Count:\ $3 N(k)$.
- Generators with $n > k$ are **killed**:\ image is zero. Count:\ $3 (N(m) - N(k))$.

### 5.2 Per-pair bit-exact panel

| Pair $(m, k)$ | Total | Survived | Killed | Verdict |
|:-------------:|:-----:|:--------:|:------:|:-------:|
| $(2, 1)$ | 15 | 6  | 9  | OK |
| $(3, 1)$ | 27 | 6  | 21 | OK |
| $(3, 2)$ | 27 | 15 | 12 | OK |
| $(4, 1)$ | 42 | 6  | 36 | OK |
| $(4, 2)$ | 42 | 15 | 27 | OK |
| $(4, 3)$ | 42 | 27 | 15 | OK |
| $(5, 1)$ | 60 | 6  | 54 | OK |
| $(5, 2)$ | 60 | 15 | 45 | OK |
| $(5, 3)$ | 60 | 27 | 33 | OK |
| $(5, 4)$ | 60 | 42 | 18 | OK |

Total generator checks:\ $15 + 27 + 27 + 42 + 42 + 42 + 60 + 60 + 60 + 60 = 435$. All bit-exact OK.

### 5.3 Structural reason

The $\mathbb{G}_a^{3 N(m)}$ generator $e_{(n, l), s}^{(m)}$ translates in the $x_{(n, l), s}^{(m)}$ coordinate of $\mathcal{H}_{\mathrm{GV}}(m)$. Its image under $\Phi_{m, k}$ acts on $\mathcal{H}_{\mathrm{GV}}(k)$ by translating in $x_{(n, l), s}^{(k)}$ when $n \le k$, and acts as zero when $n > k$ (the generator's coordinate has been removed by truncation). This is the abelian translation analog of PS-1's sector-locality.

---

## 6. Class-level $U^*$-action on $\chi$ and $\eta$

### 6.1 Interpretation C trivially compatible

v3.66.0 FO3 closed Interpretation C of the $U^*$-action via period-pairing on the image $\pi(\mathcal{H}_{\mathrm{Levi}}) \subset \mathbb{C}$:\ $U^*$ acts on the cocycle-class periods by motivic Galois conjugation. For $\chi$ and $\eta$ (depth 0, integer-valued), $U^*$ acts as the identity:

$$
U^* \cdot \chi^{(m)} = \chi^{(m)}, \qquad U^* \cdot \eta^{(m)} = \eta^{(m)}.
$$

The pro-system-level compatibility condition is therefore

$$
\boxed{P_{m, k}\!\left(U^* \cdot \psi^{(m)}\right) = U^* \cdot P_{m, k}(\psi^{(m)}) = \psi^{(k)}}
$$

for $\psi \in \{\chi, \eta\}$ and all pairs $(m, k)$. Both sides reduce to $P_{m, k}(\psi^{(m)})$, which equals $\psi^{(k)}$ by PS-1's pull-back identity.

### 6.2 Bit-exact panel

| Pair $(m, k)$ | $\chi$ lhs $=$ rhs | $\chi$ lhs $=$ $\chi^{(k)}$ | $\eta$ lhs $=$ rhs | $\eta$ lhs $=$ $\eta^{(k)}$ |
|:-------------:|:------------------:|:---------------------------:|:------------------:|:---------------------------:|
| $(2, 1)$ | OK | OK | OK | OK |
| $(3, 1)$ | OK | OK | OK | OK |
| $(3, 2)$ | OK | OK | OK | OK |
| $(4, 1)$ | OK | OK | OK | OK |
| $(4, 2)$ | OK | OK | OK | OK |
| $(4, 3)$ | OK | OK | OK | OK |
| $(5, 1)$ | OK | OK | OK | OK |
| $(5, 2)$ | OK | OK | OK | OK |
| $(5, 3)$ | OK | OK | OK | OK |
| $(5, 4)$ | OK | OK | OK | OK |

Total class-level identities:\ 10 pairs $\times$ 2 characters $\times$ 2 conditions $=$ 40 bit-exact zero residuals.

---

## 7. $SL_2$ commutativity (categorical)

The Levi $U^* = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ has the semisimple factor $SL_2$ acting on the Peter--Weyl decoration of v3.63.0 L2 (FO1 extended at $j_{\max} = 3/2$). The Peter--Weyl decoration is parameterized by the $j_{\max}$ axis, **independent of the $n_{\max}$ axis** on which PS-1's $P_{m, k}$ acts. Acting first by $SL_2$ and then truncating in $n_{\max}$ produces the same result as truncating first and then acting; this is the independence-of-axes statement.

Concretely:\ for any $g \in SL_2$ and any $\psi \in \mathcal{H}_{\mathrm{GV}}^{\mathrm{Levi}}(m)$,

$$
\Phi_{m, k}^{\mathrm{Levi}}(g \cdot \psi) = g \cdot \Phi_{m, k}^{\mathrm{Levi}}(\psi),
$$

bit-exactly, by the tensor-product structure $\mathcal{H}_{\mathrm{GV}}^{\mathrm{Levi}}(n_{\max}) = \mathcal{H}_{\mathrm{GV}}(n_{\max}) \otimes \mathcal{H}^{J^*}$ where the first factor is $n_{\max}$-dependent and the second is not. $\Phi_{m, k}^{\mathrm{Levi}} = \Phi_{m, k} \otimes \mathrm{id}_{\mathcal{H}^{J^*}}$, which commutes with $g \otimes \mathrm{id}$ for any $g \in SL_2$ acting on the second factor.

This is recorded as a categorical bit-exact statement;\ no per-cell verification is required because the commutativity is structural, not contingent.

---

## 8. Honest scope

**Closed at theorem grade (bit-exact at finite cutoff $n_{\max} \le 5$):**

- Block-diagonal Hopf-hom lift $\Phi_{m, k}$ in closed form (10 triples bit-exact cofiltered).
- $\mathbb{G}_a^{3 N(m)}$ generator compatibility:\ 435 checks across the panel.
- Class-level $U^*$-action compatibility on $\chi$ and $\eta$:\ 40 identities bit-exact.
- $SL_2$ commutativity by categorical independence-of-axes (no per-cell check needed).
- All four cohere via the abelian primitive structure of $\mathcal{H}_{\mathrm{GV}}$ established in v3.61.0 Track A, the Levi product structure of v3.63.0 L1, and the Interpretation C closure of v3.66.0 FO3.

**PS-2 reduces, transparently:**

- The Hopf-axiom compatibility under $\Phi_{m, k}$ is structurally trivial on abelian primitive substrate; the substantive content lives at the cofiltered-axiom and class-action levels. **This is honest:**\ PS-2 is genuinely lighter than PS-1 because the Hopf algebra is abelian primitive (v3.61.0 Track A), which is itself a structural finding of the program.
- $SL_2$ commutativity is categorical;\ no cell-level falsification is possible without enriching the substrate to involve $j_{\max}$ truncations alongside $n_{\max}$ (PS-3 / PS-4 territory).
- The continuum-side $U^*$-action on $F(s)$ — the non-trivial part of Interpretation C (acting on M2 as the Tate subgroup, on M3 as the standard motivic Galois action on $\mathrm{MT}(\mathbb{Q}, 1)$ odd-zeta classes) — is **not** addressed at finite cutoff in PS-2;\ $F(s)$ is an inverse-limit object and lives in PS-3.

**Sprint-scale next steps (PS-3, PS-4):**

- **PS-3 — Inverse limit $\mathcal{O}_\infty$ and extended action.** Define $\mathcal{O}_\infty = \varprojlim \mathcal{O}_{n_{\max}}$ with appropriate topology;\ extend the closed-form transitions and the $U^*$-action to the limit;\ verify continuity and density;\ carry $F(s)$ to the limit and test the non-trivial $U^*$-action on M2 / M3 components.
- **PS-4 — Endomorphism rigidity / Tannakian readiness probe.** Characterise $\mathrm{End}_{\mathrm{compatible}}(\mathcal{O}_n)$ at each cutoff;\ map remaining categorical gaps for Tannakian closure proper (abelian setting, tensor structure, rigidity, fiber functor).

**Tannakian closure itself remains a multi-year frontier** (named at v3.66.0 sprint §5.2 follow-on register). PS-2 confirms the substrate from PS-1 is compatible with the Levi-decomposed cosmic-Galois action;\ PS-3 lifts to the continuum;\ PS-4 probes endomorphism rigidity; Tannakian closure is the formal lift on top.

**Numerical observation:**

- The 435 $\mathbb{G}_a$ generator checks split exactly $3 N(k)$ survivors $+$ $3 (N(m) - N(k))$ killed at every pair, with the survival pattern bit-exactly tracking PS-1's sector-locality. This is the abelian translation analog of PS-1's truncation rule:\ the same combinatorial pattern controls both the algebra-level transition $P_{m, k}$ and the Hopf-algebraic translation action of $\mathbb{G}_a^{3 N(m)}$.

**Curve-fit audit (`feedback_audit_numerical_claims`):**

- The block-diagonal Hopf-hom matrix is **derived** from the abelian primitive structure of $\mathcal{H}_{\mathrm{GV}}$ and PS-1's $P_{m, k}$ — three independent copies by tensor product over Mellin slots, not fitted. Zero free parameters.
- The 435 generator-compatibility checks are direct membership tests against the sector-locality rule, not curve-fit alignments.
- Interpretation C triviality on $\chi, \eta$ is the v3.66.0 FO3 published structural conclusion (depth-0 integer-valued $\Rightarrow$ $U^*$-fixed), not a parameter to fit.
- Selection bias:\ the verdict gate was articulated before running computations;\ the outcome matches the strongest gate option (POSITIVE) — the bit-exact compatibility is the prediction from PS-1 + v3.61.0 Track A + v3.63.0 L1 + v3.66.0 FO3.

**Discrete-for-skeleton compliance (`feedback_discrete_for_skeleton`):**

- Every matrix entry, every per-generator check, every class-level residual is bit-exact `sympy.Integer` / `sympy.Rational`. Zero floats. Zero PSLQ. Zero transcendentals introduced at finite cutoff.

**Tag transcendentals (`feedback_tag_transcendentals`):**

- Zero transcendentals appear in this sprint. PS-2 operates entirely on Layer 1 (the bit-exact skeleton);\ continuum-side $F(s)$ M2 / M3 components (Layer 2 transcendentals tagged in v3.66.0 FO2) are deferred to PS-3.

**WH1 PROVEN unaffected.**

**Hard prohibitions check (CLAUDE.md §13.5).** No changes to natural geometry hierarchy. No fitted/empirical parameters. No deletion of negative results from §3. No removal of "conjectural" label from Paper 2 combination rule.

---

## 9. Files

### Produced
- `geovac/pro_system.py` (PS-2 additions ~280 lines appended; total file 610 lines). Public API additions:\ `MELLIN_SLOTS`, `primitive_generators(n_max)`, `n_primitive_generators(n_max)`, `HopfTransition(n_high, n_low)`, `verify_hopf_cofiltered_axiom(m, n, k)`, `verify_Ga_generator_compatibility(m, k)`, `verify_class_action_compatibility(m, k, psi_high, psi_low)`.
- `debug/compute_q5p_ps2_ustar_compatibility.py` — driver (~280 lines, 0.21 s wall, bit-exact).
- `debug/data/sprint_q5p_ps2_ustar_compatibility.json` — bit-exact data dump:\ 10 cofiltered triples, 10 $\mathbb{G}_a$ pair results, 10 class-pair results $\times$ 2 characters.
- `debug/sprint_q5p_ps2_ustar_compatibility_memo.md` — this memo.
- `tests/test_pro_system.py` — 22 new tests appended (total 51, all pass in 0.93 s).

### Used (load-bearing inputs)
- `geovac/pro_system.py` PS-1 substrate (`TransitionMap`, `verify_cofiltered_axiom`).
- `geovac/spectral_triple.py` (`FockSpectralTriple` for $\chi, \eta$ class extraction).
- `debug/sprint_q5p_ps1_transitions_memo.md` (PS-1 closure).
- `debug/sprint_q5p_stage2_hopf_memo.md` (v3.61.0 Track A;\ abelian primitive substrate).
- `debug/sprint_q5p_levi_synthesis_memo.md` (v3.63.0 L1;\ Levi decomposition).
- `debug/sprint_q5p_fo2_fo3_mt_period_memo.md` (v3.66.0 FO3;\ Interpretation C closure for $\chi, \eta, F(s)$).

### Published references
- Connes, A.; Marcolli, M. ``Renormalization and motivic Galois theory.'' Int. Math. Res. Not. (2004), 76:\ 4073--4091. (arXiv:math/0409306.)
- Connes, A.; Kreimer, D. ``Hopf algebras, renormalization and noncommutative geometry.'' Comm. Math. Phys. 199 (1998), 203--242.
- Marcolli, M.; van Suijlekom, W. D. ``Gauge networks in noncommutative geometry.'' J. Geom. Phys. 75 (2014), 71--91 (= arXiv:1301.3480).

---

## 10. Paper-edit recommendations (PI to apply)

### 10.1 Paper 55 \S subsec:open_m2_m3 — ONE new paragraph after the PS-1 paragraph

Insert after the PS-1 paragraph from v3.67.0:

```latex
\emph{Pro-system lockdown PS-2:\ $U^*$-action compatibility with PS-1
transitions (Sprint Q5'-ProSystem-Lockdown, PS-2 sub-track,
2026-06-06;\ memo
\texttt{debug/sprint\_q5p\_ps2\_ustar\_compatibility\_memo.md};\
data
\texttt{debug/data/sprint\_q5p\_ps2\_ustar\_compatibility.json}).}
The block-diagonal Hopf-hom lift $\Phi_{m, k}: \mathcal{H}_{\mathrm{GV}}(m)
\to \mathcal{H}_{\mathrm{GV}}(k)$ of PS-1's $P_{m, k}$ acts on the
primitive generators by $\Phi(x_{(n, l), s}^{(m)}) = x_{(n, l), s}^{(k)}$
if $n \le k$ and $0$ otherwise, where $s \in \{0, 1, 2\}$ indexes the
master Mellin engine slots (M1, M3, M2 per Paper 18 \S III.7).  The
$3 N(k) \times 3 N(m)$ matrix of $\Phi_{m, k}$ is the block-diagonal
$\mathrm{diag}(P_{m, k}, P_{m, k}, P_{m, k})$.  Three bit-exact
verifications at $n_{\max} \le 5$:\ (i) the cofiltered axiom
$\Phi_{m, k} = \Phi_{n, k} \cdot \Phi_{m, n}$ holds bit-exact at all
10 triples (10 zero residuals);\ (ii) the 435 $\mathbb{G}_a^{3 N(m)}$
translation generators across the 10 pairs $(m, k) \le 5$ all
satisfy the sector-locality rule $\Phi_{m, k}(e_g^{(m)} \cdot x) =
e_g^{(k)} \cdot \Phi_{m, k}(x)$ if $n \le k$, $0$ otherwise, bit-exact;
(iii) the class-level $U^*$-action on $\chi, \eta$ (depth-0 trivial by
v3.66.0 FO3 Interpretation C) commutes with $\Phi_{m, k}$ bit-exact at
all 10 pairs (40 identities).  The $SL_2$ factor of $U^* =
\mathbb{G}_a \rtimes SL_2$ acts on the Peter--Weyl decoration
($j_{\max}$ axis), independent of the $n_{\max}$ axis on which
$\Phi_{m, k}$ acts;\ commutativity is categorical
(independence-of-axes), not requiring a per-cell check.  Total
bit-exact zero residuals across the panel:\ $10 + 435 + 40 = 485$.
The cosmic-Galois $U^*$ lifts to the PS-1 pro-system by construction;\
PS-3 carries the inverse limit $\mathcal{O}_\infty = \varprojlim
\mathcal{O}_{n_{\max}}$ and extends the $U^*$-action to it (including
the non-trivial part of Interpretation C on the continuum $F(s)$).
```

### 10.2 Paper 32 — no edit needed at PS-2

The Levi substrate (Paper 32 §VIII `rem:q5p_levi_synthesis_substrate` from v3.63.0) and the Interpretation C closure (`rem:q5p_interpretation_C_closure` from v3.66.0) already cover the framework's structural content;\ PS-2 is documentation of compatibility at the pro-system level, captured in Paper 55. If a forward pointer becomes useful at PS-3 (inverse limit), Paper 32 §VIII is the natural insertion point.

### 10.3 Paper 18 — no edit needed

PS-2 operates on the abelian primitive substrate at Layer 1;\ Paper 18's master Mellin engine §III.7 framing is upstream and unaffected.

---

## 11. One-line verdict

**POSITIVE.** The Levi-decomposed cosmic-Galois $U^* = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ acts compatibly with PS-1's closed-form transitions:\ block-diagonal Hopf-hom lift $\Phi_{m, k}$ satisfies the cofiltered axiom bit-exact at all 10 triples in $\{1, \ldots, 5\}$ (10 identities);\ all 435 $\mathbb{G}_a$ generators across the 10 pairs $(m, k) \le 5$ obey the sector-locality survival rule bit-exact;\ class-level $U^*$-action on $\chi, \eta$ (Interpretation C depth-0 trivial) commutes with $\Phi_{m, k}$ bit-exact (40 identities);\ $SL_2$ commutativity is categorical by independence of axes. 485 / 485 bit-exact zero residuals across the panel. The cosmic-Galois $U^*$ lifts to the PS-1 pro-system by construction;\ PS-3 (inverse limit and continuum $F(s)$ extension) can proceed on this substrate.
