# Sprint Q5'-ProSystem-Lockdown PS-1 — transition maps, cofiltered axiom, $n_{\max} = 5$ extension

**Date:** 2026-06-06 (single-thread sprint, first of four sub-tracks in the Pro-System-Lockdown plan)
**Sprint:** PS-1 of Q5'-ProSystem-Lockdown (PS-2, PS-3, PS-4 to follow, single-threaded per PI direction)
**Driver:** `debug/compute_q5p_ps1_transitions.py`
**Module:** `geovac/pro_system.py` (new, ~300 lines)
**Data:** `debug/data/sprint_q5p_ps1_transitions.json`
**Wall time:** 0.14 s
**Discipline:** bit-exact `sympy.Rational` / `sympy.Integer` throughout. No floats. No PSLQ. No transcendentals introduced.

---

## 1. TL;DR

**Verdict: POSITIVE.** The truncated Camporesi--Higuchi sector-idempotent pro-system $\{\mathcal{O}_{n_{\max}}\}_{n_{\max} \ge 1}$ has explicit closed-form transition maps $P_{m, k}: \mathcal{O}_m \to \mathcal{O}_k$ (algebra homomorphisms, realised as $N(k) \times N(m)$ 0/1 integer matrices in canonical lex order) that satisfy the cofiltered axiom $P_{m, k} = P_{n, k} \cdot P_{m, n}$ bit-exact for every triple $1 \le k < n < m \le 5$. Per-sector closed forms $\chi_{(n, l)}, \eta_{(n, l)}$ from v3.60.0 extend bit-exact to the new cell $n_{\max} = 5$ (6 new sectors, 20 total). Pull-back $P^*_{m, k}[\psi^{(m)}] = [\psi^{(k)}]$ holds bit-exact for ALL pairs $1 \le k < m \le 5$ (10 pairs $\times$ 2 characters), not only for consecutive pairs.

**130 / 130 bit-exact zero residuals:**

- 100 per-sector closed-form predictions ($50$ $\chi$ + $50$ $\eta$ across the five cutoffs).
- 10 cofiltered axiom triples ($\binom{5}{3}$ = 10 matrix-level identities $P_{m,k} = P_{n,k} \cdot P_{m,n}$).
- 20 all-pairs pull-back identities ($\binom{5}{2}$ = 10 pairs $\times$ 2 characters).

**PS-1 deliverable:** v3.60.0's "pull-back compatibility at consecutive cutoffs" is promoted to **a closed-form inverse system with explicit transition morphisms and bit-exact cofiltered structure**. The pro-system is now a closed-form object in the strict sense; PS-2 (U\* compatibility with transitions) can proceed on this substrate.

---

## 2. Verdict against decision gate

| Gate | Verdict |
|:-----|:--------|
| POSITIVE | **selected** — bit-exact cofiltered axiom at all 10 triples; bit-exact pull-back at all 10 pairs $\times$ 2 characters; closed-form extension to $n_{\max} = 5$ verified. Total: 130 zero residuals out of 130 identities tested. |
| BORDERLINE | not selected — closure is bit-exact at the full panel; no near-zero residuals, no partial closure. |
| STOP | rejected — every identity verified bit-exact; the falsifier extension to $n_{\max} = 5$ reproduces the closed-form structure without modification. |

---

## 3. What PS-1 adds beyond v3.60.0

The Sprint Q5'-Stage1-Prosystem result (v3.60.0) verified pull-back compatibility at the **value level** for the JLO HP$^{\mathrm{even}}$ class $\chi$ and the CM-$\eta$ residue class $\eta$ at consecutive cutoff pairs $(n+1, n) \in \{(1,2), (2,3), (3,4)\}$ via the bare equation $P^*_{n+1 \to n}[\psi^{(n+1)}] = [\psi^{(n)}]$ in each pair separately. The closed-form sector evolution rule was identified ($\chi$ and $\eta$ depend only on the sector label, not on $n_{\max}$).

PS-1 supplies three additional structural ingredients:

1. **The transition map as a closed-form object.** Each $P_{m, k}: \mathcal{O}_m \to \mathcal{O}_k$ is now realised as a $N(k) \times N(m)$ 0/1 integer matrix in canonical lex order on sectors, accessible via `TransitionMap(m, k).matrix` and `TransitionMap(m, k).apply_to_vector(v)` in `geovac/pro_system.py`. This is the explicit form of the algebra homomorphism, not merely a per-cell value check.

2. **The cofiltered axiom.** For every triple $1 \le k < n < m \le 5$, the matrix identity $P_{m, k} = P_{n, k} \cdot P_{m, n}$ holds bit-exact (10 zero residuals, computed by `verify_cofiltered_axiom(m, n, k)`). This is the load-bearing condition that promotes "compatible family at consecutive cutoffs" to "inverse system" in the categorical sense. Without it there is no well-defined pro-object to take a limit of; with it, PS-3 can construct $\mathcal{O}_\infty = \varprojlim \mathcal{O}_{n_{\max}}$ as a sequential limit and PS-2 can lift the cosmic-Galois $U^\ast$ to the pro-system.

3. **Falsifier extension to $n_{\max} = 5$.** A new cutoff cell (six new sectors $(5, 0), \ldots, (5, 5)$, total 20 sectors) is added to the pull-back and closed-form panels. All 6 new sector predictions for $\chi$ and 6 for $\eta$ match the closed forms bit-exact, and the four new pull-back pairs $P^*_{5, k}$ for $k \in \{1, 2, 3, 4\}$ all close bit-exact. The structural pattern from v3.60.0 ($\chi, \eta$ are sector-local, $n_{\max}$-independent per sector) is not a finite-cutoff artifact.

---

## 4. Closed-form objects (`geovac/pro_system.py`)

### 4.1 Sector enumeration

$$
\boxed{
\text{Sectors at cutoff } n_{\max}: \quad \{(n, l) : 1 \le n \le n_{\max}, \ 0 \le l \le n\}, \quad
N(n_{\max}) = \frac{n_{\max}(n_{\max} + 3)}{2}.
}
$$

| $n_{\max}$ | $N(n_{\max})$ | $\dim \mathcal{H}$ |
|:----------:|:-------------:|:------------------:|
| 1 | 2 | 4 |
| 2 | 5 | 16 |
| 3 | 9 | 40 |
| 4 | 14 | 80 |
| 5 | 20 | 140 |

The closed form $N(n_{\max}) = n_{\max}(n_{\max}+3)/2$ is verified bit-exact at every cutoff against the explicit `FockSpectralTriple.n_sectors`.

### 4.2 Transition map $P_{m, k}$

For sectors in canonical lex order $(1, 0), (1, 1), (2, 0), (2, 1), (2, 2), (3, 0), \ldots, (m, m)$, the matrix of $P_{m, k}: \mathcal{O}_m \to \mathcal{O}_k$ is the $N(k) \times N(m)$ projection

$$
[P_{m, k}]_{i, j} = \begin{cases} 1 & \text{if } 1 \le i \le N(k) \text{ and } i = j \\ 0 & \text{otherwise.} \end{cases}
$$

In matrix form $P_{m, k} = [I_{N(k)} \mid 0_{N(k) \times (N(m) - N(k))}]$ — the first $N(k)$ columns are the identity, the last $N(m) - N(k)$ columns are zero. As an algebra homomorphism it sends $e_{(n, l)} \mapsto e_{(n, l)}$ if $n \le k$, and $e_{(n, l)} \mapsto 0$ if $n > k$.

Examples at the new cutoff $n_{\max} = 5$:

| Pair $(m, k)$ | Shape | nnz | First-block trace |
|:-------------:|:-----:|:---:|:----:|
| $(5, 4)$ | $14 \times 20$ | 14 | 14 |
| $(5, 3)$ | $9 \times 20$  | 9  | 9  |
| $(5, 1)$ | $2 \times 20$  | 2  | 2  |

### 4.3 Cofiltered axiom

For every triple $1 \le k < n < m$, the closed-form transitions satisfy

$$
\boxed{P_{m, k} = P_{n, k} \cdot P_{m, n}}
$$

(matrix product, $N(k) \times N(m) = N(k) \times N(n) \cdot N(n) \times N(m)$). The structural reason is direct: $P_{m, k}$ keeps the first $N(k)$ coordinates of an $N(m)$-vector; composing "keep first $N(n)$" then "keep first $N(k)$" produces "keep first $N(k)$" identically.

Bit-exact verification at all 10 triples in $\{1, 2, 3, 4, 5\}$:

| Triple $(m, n, k)$ | $\dim(\text{matrix})$ | Residual$^2$ | Verdict |
|:------------------:|:---------------------:|:------------:|:-------:|
| $(3, 2, 1)$ | $2 \times 9$  | 0 | OK |
| $(4, 2, 1)$ | $2 \times 14$ | 0 | OK |
| $(4, 3, 1)$ | $2 \times 14$ | 0 | OK |
| $(4, 3, 2)$ | $5 \times 14$ | 0 | OK |
| $(5, 2, 1)$ | $2 \times 20$ | 0 | OK |
| $(5, 3, 1)$ | $2 \times 20$ | 0 | OK |
| $(5, 3, 2)$ | $5 \times 20$ | 0 | OK |
| $(5, 4, 1)$ | $2 \times 20$ | 0 | OK |
| $(5, 4, 2)$ | $5 \times 20$ | 0 | OK |
| $(5, 4, 3)$ | $9 \times 20$ | 0 | OK |

10 / 10 bit-exact zero residuals.

---

## 5. Per-sector closed forms at $n_{\max} = 5$

The v3.60.0 closed forms

$$
\chi_{(n, l)} = \begin{cases} +2 & l < n \\ -2n & l = n \end{cases}, \qquad
\eta_{(n, l)} = \begin{cases} (2l + 1)(2n + 1) & l < n \\ n(2n + 1) & l = n \end{cases}
$$

extend bit-exact to the new cell. At $n_{\max} = 5$ the six new sectors are $(5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5)$, and the closed forms give

| Sector | $\chi$ | $\eta$ |
|:------:|:------:|:------:|
| $(5, 0)$ | $+2$  | $11$ |
| $(5, 1)$ | $+2$  | $33$ |
| $(5, 2)$ | $+2$  | $55$ |
| $(5, 3)$ | $+2$  | $77$ |
| $(5, 4)$ | $+2$  | $99$ |
| $(5, 5)$ | $-10$ | $55$ |

all matched bit-exact by direct trace computation on `FockSpectralTriple(n_max=5)`. Totals at $n_{\max} = 5$:

- $\sum_s \chi_s = 5 \cdot 2 - 10 = 0$ (McKean--Singer index $\chi(S^3) = 0$ preserved).
- $\sum_s \eta_s = M_3(5) = 5 \cdot 6^2 \cdot 7 / 2 = 630$ (Sub-Sprint 2a polynomial bit-exact).

Cumulative across all five cutoffs:

| $n_{\max}$ | $\dim \mathcal{H}$ | $\sum \chi_s$ | $M_3 = \sum \eta_s$ |
|:----------:|:------------------:|:-------------:|:-------------------:|
| 1 | 4   | 0 | 6   |
| 2 | 16  | 0 | 36  |
| 3 | 40  | 0 | 120 |
| 4 | 80  | 0 | 300 |
| 5 | 140 | 0 | 630 |

Polynomial closed forms $M_1(n_{\max}) = 2 n_{\max}(n_{\max}+1)(n_{\max}+2) / 3$ and $M_3(n_{\max}) = n_{\max}(n_{\max}+1)^2(n_{\max}+2) / 2$ verified bit-exact at all five cutoffs.

---

## 6. All-pairs pull-back identity

For every pair $1 \le k < m \le 5$ and both characters $\psi \in \{\chi, \eta\}$,

$$
P^*_{m, k}[\psi^{(m)}] = [\psi^{(k)}]
$$

at the class-vector level (the dim-$N(m)$ vector restricts to the dim-$N(k)$ vector bit-exact). The 10 pairs are $(2,1), (3,1), (3,2), (4,1), (4,2), (4,3), (5,1), (5,2), (5,3), (5,4)$. All 20 identities (10 pairs $\times$ 2 characters) close bit-exact.

This is structurally a corollary of (cofiltered axiom + sector-locality), but verifying it independently across all pairs (not just consecutive ones) is the standard inverse-system falsifier and is part of the PS-1 deliverable.

---

## 7. Honest scope

**Closed at theorem grade (bit-exact at finite cutoff $n_{\max} \le 5$):**

- Transition map $P_{m, k}$ as a closed-form $N(k) \times N(m)$ 0/1 integer matrix.
- Cofiltered axiom $P_{m, k} = P_{n, k} \cdot P_{m, n}$ at all 10 triples in $\{1, \ldots, 5\}$.
- Per-sector closed forms $\chi_{(n, l)}, \eta_{(n, l)}$ at all 50 sectors across $n_{\max} \in \{1, \ldots, 5\}$.
- All-pairs pull-back $P^*_{m, k}[\psi^{(m)}] = [\psi^{(k)}]$ at all 10 pairs $\times$ 2 characters.
- Total sum invariants $M_1$ and $M_3$ at all five cutoffs.

**Sprint-scale next steps (PS-2, PS-3, PS-4 of this lockdown sprint):**

- **PS-2 — $U^\ast$ compatibility with transitions.** Verify $U^\ast$-action commutes with the closed-form transitions of PS-1, for each generator of $U^\ast = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ (v3.63.0 Levi decomposition) on each of the three Interpretation-C-closed characters $\chi$, $\eta$, $F(s)$. If bit-exact: $U^\ast$ lifts to the pro-system by construction.
- **PS-3 — Inverse limit $\mathcal{O}_\infty$ and extended action.** Define $\mathcal{O}_\infty = \varprojlim \mathcal{O}_{n_{\max}}$ with natural topology; extend the closed-form transitions and the $U^\ast$-action to the limit; verify continuity and density. Carries $F(s)$ (continuum Mellin lift, v3.66.0 FO2) to the limit object.
- **PS-4 — Endomorphism rigidity / Tannakian readiness.** Characterise $\mathrm{End}_{\mathrm{compatible}}(\mathcal{O}_n)$ at each finite cutoff (is it $U^\ast$ exactly, or $U^\ast$ + garbage?). Map the remaining categorical gaps for Tannakian closure proper: abelian setting, tensor structure, rigidity, fiber functor. Named-gap list becomes the next sprint after this one.

**Tannakian closure itself remains a multi-year frontier**, named in v3.66.0 sprint §5.2 follow-on register and in the v3.66.0 wrap-up conversation with the PI as "the canonical next step after the pro-system is locked down." PS-1 is the substrate layer of that lockdown; Tannakian closure is the formal lift on top of PS-3's inverse-limit object.

**Numerical observation:**

- The bit-exact cleanliness of the cofiltered axiom is the structurally cleanest possible inverse-system data:\ the transition matrices are 0/1 projections that compose under matrix multiplication to the direct projection, with no twists, no coboundary modifications, and no per-cutoff renormalization. This is sector-locality of the Camporesi--Higuchi shell decomposition made categorical:\ each shell's contribution to the algebra is independent of higher shells.

**Curve-fit audit (`feedback_audit_numerical_claims`):**

- The closed-form transition matrix is **derived** from the sector-locality of the Camporesi--Higuchi shell decomposition (sector idempotents at lower cutoff are sector idempotents at higher cutoff verbatim), not fitted. Zero free parameters.
- The cofiltered axiom $P_{m, k} = P_{n, k} \cdot P_{m, n}$ is a structural identity for sector-projection maps and is verified bit-exact at all 10 triples — direct matrix multiplication, not curve-fit alignment.
- All-pairs pull-back identity at the class-vector level is the same structural identity at the cocycle-data level; verified independently across 10 pairs $\times$ 2 characters.
- Selection bias: the verdict gate was articulated before running computations; the outcome matches the strongest gate option (POSITIVE), not a fall-back. The bit-exact cleanliness is the prediction from the sector-locality of the Camporesi--Higuchi shell decomposition; no other outcome was consistent with that structural input.

**Discrete-for-skeleton compliance (`feedback_discrete_for_skeleton`):**

- Every transition matrix entry, every per-sector class value, every cofiltered-axiom residual, and every pull-back identity is bit-exact `sympy.Integer` / `sympy.Rational`. Zero floats. Zero PSLQ. Zero transcendentals introduced at finite cutoff.

**Tag transcendentals (`feedback_tag_transcendentals`):**

- Zero transcendentals appear in this sprint. PS-1 operates entirely on Layer 1 (the bit-exact skeleton); the Mellin lift $F(s)$ (Layer 2) is deferred to PS-3.

**WH1 PROVEN unaffected.** This sprint extends v3.60.0 with the cofiltered-axiom dimension; it does not test or extend WH1's propinquity-side foundation.

**Hard prohibitions check (CLAUDE.md §13.5).** No changes to natural geometry hierarchy. No fitted/empirical parameters. No deletion of negative results from §3. No removal of "conjectural" label from Paper 2 combination rule.

---

## 8. Files

### Produced
- `geovac/pro_system.py` — new module (~300 lines). Public API: `sectors_at_cutoff(n_max)`, `N_sectors(n_max)`, `TransitionMap(n_high, n_low)`, `compose(P_outer, P_inner)`, `verify_cofiltered_axiom(m, n, k)`. Bit-exact sympy throughout.
- `debug/compute_q5p_ps1_transitions.py` — driver (~350 lines, 0.14 s wall, bit-exact).
- `debug/data/sprint_q5p_ps1_transitions.json` — bit-exact data dump: 5 cutoffs $\times$ sectors $\times$ class vectors; 10 cofiltered triples + residuals; 10 pull-back pairs $\times$ 2 characters; 3 transition matrix examples at $n_{\max} = 5$.
- `debug/sprint_q5p_ps1_transitions_memo.md` — this memo.

### Used (load-bearing inputs)
- `geovac/spectral_triple.py` (`FockSpectralTriple` for explicit $\Lambda, \gamma, A, D$ at each $n_{\max}$).
- `debug/sprint_q5p_prosystem_memo.md` (v3.60.0; closed-form sector evolution rule and consecutive-pair pull-back compatibility).
- `debug/sprint_q5p_hp_round3_followons_2026_06_06_memo.md` (v3.66.0; Interpretation C of $U^\ast$-action; $F(s) \subset \mathrm{MT}(\mathbb{Q}, 1)$; substrate locked for PS-2 forward).
- `debug/sprint_q5p_levi_arc_2026_06_06_memo.md` (v3.63.0; Levi decomposition $U^\ast = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$, substrate for PS-2 generators).

### Published references
- Connes, A.; van Suijlekom, W. D. ``Spectral truncations in noncommutative geometry and operator systems.'' Comm. Math. Phys. 383 (2021).
- Marcolli, M.; van Suijlekom, W. D. ``Gauge networks in noncommutative geometry.'' J. Geom. Phys. 75 (2014), 71--91 (= arXiv:1301.3480).
- Paper 38 \S VIII L4 (`papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex`; Berezin reconstruction substrate).
- Paper 55 \S\ref{subsec:open_m2_m3} (`papers/group3_foundations/paper_55_periods_of_geovac.tex`; Q5' Stage 1 construction context).

---

## 9. Paper-edit recommendations (PI to apply)

### 9.1 Paper 55 \S subsec:open_m2_m3 — ONE new paragraph after the v3.66.0 FO3 closure paragraph

Insert after the v3.66.0 Interpretation C paragraph:

```latex
\emph{Pro-system lockdown PS-1 (Sprint Q5'-ProSystem-Lockdown, PS-1
sub-track, 2026-06-06; memo
\texttt{debug/sprint\_q5p\_ps1\_transitions\_memo.md}; data
\texttt{debug/data/sprint\_q5p\_ps1\_transitions.json}; module
\texttt{geovac/pro\_system.py}).} The transition maps
$P_{m, k}: \mathcal{O}_m \to \mathcal{O}_k$ between truncated
Camporesi--Higuchi sector-idempotent algebras are realised as
closed-form $N(k) \times N(m)$ 0/1 integer matrices in canonical lex
order on sectors:\ $P_{m, k} = [I_{N(k)} \mid 0_{N(k) \times (N(m) -
N(k))}]$, sending $e_{(n, l)} \mapsto e_{(n, l)}$ if $n \le k$ and
$0$ otherwise. The cofiltered axiom $P_{m, k} = P_{n, k} \cdot
P_{m, n}$ is verified bit-exact for every triple
$1 \le k < n < m \le 5$ (10 matrix-level identities, 10 bit-exact
zero residuals), promoting the v3.60.0 result (pull-back
compatibility at consecutive cutoffs) to a closed-form inverse
system in the strict categorical sense.  The per-sector closed
forms $\chi_{(n, l)}$ and $\eta_{(n, l)}$ extend bit-exact to
$n_{\max} = 5$ (20 sectors total), with all-pairs pull-back
$P^*_{m, k}[\psi^{(m)}] = [\psi^{(k)}]$ holding bit-exact for the 10
pairs $\{(m, k) : 1 \le k < m \le 5\}$ and both characters
($\chi, \eta$):\ 20 bit-exact zero residuals at the class-vector
level.  Total bit-exact zero residuals across the panel:\
$100 + 10 + 20 = 130$.  This is the substrate layer of the
pro-system lockdown program;\ PS-2 verifies that the cosmic-Galois
$U^\ast = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ acts
compatibly with $P_{m, k}$, PS-3 constructs the inverse limit
$\mathcal{O}_\infty = \varprojlim \mathcal{O}_{n_{\max}}$ and
extends the $U^\ast$-action to it, PS-4 probes endomorphism
rigidity and Tannakian readiness.
```

### 9.2 Paper 32 — no edit needed

The transition-map closed form and cofiltered axiom are infrastructure for the Stage-1 pro-system; Paper 32 §VIII already references the v3.60.0 closure via `rem:q5p_prosystem_functoriality`. PS-1 is documented in Paper 55 (the construction paper for the Q5' arc). If a forward pointer in Paper 32 becomes useful at PS-3 (the inverse limit), that's the natural insertion point.

### 9.3 Paper 18 — no edit needed

PS-1 operates on Layer 1 (the bit-exact skeleton) before any Mellin lift; Paper 18's master Mellin engine §III.7 framing is upstream of this work and unaffected.

---

## 10. One-line verdict

**POSITIVE.** The truncated Camporesi--Higuchi sector-idempotent pro-system has closed-form transition maps $P_{m, k}$ as $N(k) \times N(m)$ 0/1 integer matrices in canonical lex order; the cofiltered axiom $P_{m, k} = P_{n, k} \cdot P_{m, n}$ holds bit-exact at all 10 triples in $\{1, \ldots, 5\}$; the per-sector closed forms $\chi_{(n, l)}, \eta_{(n, l)}$ and all-pairs pull-back identities extend bit-exact to $n_{\max} = 5$; 130 / 130 bit-exact zero residuals across the panel. PS-1 deliverable of the Pro-System-Lockdown sprint is closed; PS-2 ($U^\ast$ compatibility with transitions) can proceed on this closed-form substrate.
