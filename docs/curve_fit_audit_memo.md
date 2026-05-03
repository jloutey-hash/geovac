# Curve-Fit / Pattern-Match Audit Memo

**Date:** 2026-05-02
**Purpose:** Identify all PSLQ-shaped or "small-integer combination found
numerically" claims across the GeoVac papers, classify by epistemic
status, and surface recommended paper edits for the refactor sprint.

This memo is an audit, not an action plan. Specific edits are proposed
in §5 for PI review before execution.

---

## 1. Classification scheme

Three tiers:

| Tier | Label | Definition |
|:----:|:------|:-----------|
| **A** | **Derived** | Closed-form result with proof — symbolic identity, theorem, exact derivation. No PSLQ involved, or PSLQ used only to *confirm* a result independently derived. |
| **B** | **PSLQ-identified, structurally explained** | Numerical PSLQ match *plus* a structural argument showing why the form must hold (motivic constraint, exact Hurwitz reduction, Apéry-style irrationality, etc.). PSLQ finds the form; structure proves it. |
| **C** | **Pattern-matched** | Numerical match found by combinatorial or PSLQ search after the target was known, with no independent derivation. Form has at least as many free choices as the data constrains. The claim is "these known constants combine to hit this number to N digits." |

Tier C is where curve-fitting risk lives. A Tier C claim can be true,
but the evidence offered is consistent with chance.

The audit assesses **how each claim is presented**, not whether it
might eventually be derived. A Tier C claim that is honestly framed as
"a numerical pattern of unexplained origin" is fine; one framed as
"a concrete bridge between X and Y" is overstated.

---

## 2. Survey by paper

### 2.1 Paper 2 (α conjecture) — `papers/core/paper_2_alpha.tex` (move to `papers/observations/`)

| ID | Claim | Tier | Current framing | Notes |
|:---|:------|:----:|:----------------|:------|
| P2-1 | $K = \pi(B + F - \Delta)$ matches $1/\alpha$ at $8.8\times 10^{-8}$ | C | Explicitly conjectural; statistical validation section reports $p = 5.2 \times 10^{-9}$ from 1.92e9 formula search | Headline claim. The paper itself states "it does not, by itself, constitute a derivation." Three structural decompositions (B = Casimir trace, F = $D_{n^2}(d_{\max})$, Δ = $g_3^{\text{Dirac}}$) ARE Tier A — the *combination* is Tier C. |
| P2-2 | Spectral determinant near-miss: $4\pi^2 e^{\zeta(3)/(2\pi^2)} \approx 41.957 \approx B = 42$ (0.1%) | C | "Whether this near-miss reflects a deeper connection... or is merely a coincidence, remains an open question" | Honestly framed. Replacing $B = 42$ with this determinant makes CODATA agreement 11,000× worse — explicitly noted. |
| P2-3 | $R_{\text{predict}} = K - \alpha^{-1} - \alpha^2 \approx \pi^3 \alpha^3$ to 0.25% | C | "A numerical pattern of unexplained structural origin... a suggestive hint, not a derivation" | Honestly framed in Open Questions section. Structurally cross-checked against S⁵ Seeley-DeWitt and explicitly *not* found to derive. |
| P2-4 | Three ad-hoc bridge corrections to spectral determinant: $+\pi/72$, $+1/24$ | C | "Both remain *ad hoc*" | Honestly labeled. |
| P2-5 | $\Delta^{-1} = g_3^{\text{Dirac}} = 2(n+1)(n+2)\|_{n=3} = 40$ | A | Closed-form Camporesi-Higuchi degeneracy | Exact identity from spinor representation theory. Not curve-fit. |
| P2-6 | $F = D_{n^2}(d_{\max}) = \zeta_R(2) = \pi^2/6$ | A | Identified as Dirichlet series of Fock degeneracy at packing exponent | Exact symbolic identity (sympy). |
| P2-7 | $B = 42$ as Casimir truncation at $m = 3$ | A | Closed-form $\sum_{n,l} (2l+1) l(l+1)$ over $n \le 3$ | Exact integer from finite sum. |
| P2-8 | Statistical validation $p = 5.2 \times 10^{-9}$ | A* | Reproducible search documented | Search well-defined; *interpretation* depends on prior over "natural" formulas. The number is real; what it implies is debatable. |

**Paper 2 verdict:** Mostly honest about what is and isn't derived. The
headline conjecture (P2-1) is correctly labeled conjectural, with
structural decomposition of each ingredient and explicit acknowledgment
that the combination rule is not derived. P2-2 through P2-4 are
appropriately framed as near-misses or hints. **No major edits required
for Paper 2 itself**; the move from `core/` to `observations/` reflects
the conjectural nature of the headline more accurately than the current
Core placement.

---

### 2.2 Paper 28 (QED on S³) — `papers/observations/paper_28_qed_s3.tex`

| ID | Claim | Tier | Current framing | Notes |
|:---|:------|:----:|:----------------|:------|
| P28-1 | $c_2 = (2 - B\Delta - F\Delta - F/B)/5 = 19/100 - 41\pi^2/25200$ | C | "Concrete numerical bridge between Paper 2 and Paper 28" | **OVERSTATED.** PSLQ-identified at 8 digits on 13-digit data. T9 forces form to be rational + rational·π², so we are explaining 2 numbers. Form has 4 numerator terms + denominator 5 (≥ 5 free choices) — more than the 2 numbers it explains. **Smoking gun:** $c_3 = -5.946 \times 10^{-7}$ is nonzero at 200σ and the *same form does not predict it* (no extension found in same basis at $|\text{coeff}| \le 10^4$). A real structural identity would extend. This is the c₂ claim flagged in PI conversation. |
| P28-2 | $D_5 = (497/128)\zeta(8) - (467/16)\zeta(9) + (385/32)\zeta(3)\zeta(6) + (75/8)\zeta(4)\zeta(5)$ | B | PSLQ at 250 dps with subsequent symbolic verification of $c_5(n)$ for $n=1..15$ | Tier B because the per-shell $c_p(n)$ Sommerfeld coefficients have closed-form structure, and the PSLQ identification is verified at exact rationals across 15 shells. |
| P28-3 | $D_6$ analytical assembly from 7 weight-11 Euler sums | B | PSLQ at 400 dps with sub-basis decomposition; product survival rule extends to p=6 | Tier B. Each sub-basis decomposition individually PSLQ-identified, then assembled symbolically. The product survival rule (max(0, ⌊(2p−5)/2⌋) surviving products) extends across $D_2..D_6$ — *that* is a derived prediction, not a pattern. |
| P28-4 | $D_{\text{even}}(s) - D_{\text{odd}}(s) = 2^{s-1}[\beta(s) - \beta(s-2)]$ | A | Closed-form via two Hurwitz identities at quarter-integer shifts | Exact symbolic equality; PSLQ used only for confirmation. |
| P28-5 | $S_{\min} = 2.4799369380...$ irreducible to 200 digits | A | 15 PSLQ failures across 100+ basis | Tier A as a *negative* result: the irreducibility claim is independent of pattern matching. |
| P28-6 | Self-energy structural zero $\Sigma(n_{\text{ext}} = 0) = 0$ | A | Closed-form from vertex parity selection rule (n₁+n₂+q odd, n₂=0 → 2n₁ odd impossible) | Theorem. |
| P28-7 | T9 squared Dirac spectral zeta theorem | A | Symbolic proof via Bernoulli reduction | Theorem. |
| P28-8 | $F_2 \sim 3.507 \cdot n^{-0.573}$ power law | C | "Power-law fit, R² = 0.99990" | 5 data points, 2 parameters; near-perfect $R^2$ is expected. The interpretation $F_2 \to 0$ as $n \to \infty$ is a *fit-based extrapolation*, not derived. The Richardson extrapolation giving $F_\infty \approx -0.22$ (unphysical) is correctly noted as a sanity flag. |
| P28-9 | $C_{\text{VP}}/C_{\text{SE}} = 0.1035 \approx 3/29$ | C | "Closest simple rational" | Honestly framed. The ratio is constant to 0.83% across $n_{\max} = 3,4,5$ — the *constancy* is real, but the rational identification is just "nearest simple fraction." |
| P28-10 | Pendant-edge theorem $\Sigma(\text{GS}) = 2(n_{\max} - 1)/n_{\max}$ | A | Exact match across $n_{\max} = 2..6$; mechanism: pendant vertex degree-1, path-graph Laplacian inverse | Theorem. |
| P28-11 | Selection rule census 1/8 (scalar graph), 4/8 (Dirac graph), 7/8 (vector scalar), 8/8 (vector Dirac) | A | Each rule individually verified | Verified per-rule; partition is empirical census, not pattern. |

**Paper 28 verdict:** P28-1 (the c₂ formula) needs to be reframed.
The current "concrete bridge" language overstates; the data is
consistent with c₂ being a Tier C pattern-match. The proposed edit
(see §5) demotes this to "the same Paper 2 invariants appear in a
PSLQ-found form for c₂; absence of structural prediction for c₃
suggests this may be coincidental." P28-8 should add a caveat that
the "$F_2 \to 0$" extrapolation is fit-based and the Richardson
extrapolation gives an unphysical value, making the asymptotic
claim weak.

---

### 2.3 Paper 18 (exchange constants) — `papers/core/paper_18_exchange_constants.tex`

| ID | Claim | Tier | Current framing | Notes |
|:---|:------|:----:|:----------------|:------|
| P18-1 | Various near-misses listed in §VII (Hopf-twist $S^3$ vs $S^1\times S^2$ at $6.8\times 10^{-3}$, etc.) | C | "Cleanest near-miss" — labeled as such | Honest framing; documented as eliminated mechanisms. |
| P18-2 | Three-tier exchange constant taxonomy | A | Each transcendental classified by the operator/bundle/diagram axis it lives on | Taxonomic framework, not a numerical claim. |
| P18-3 | PSLQ-confirmed motivic-weight assignment for ζ(3) at Dirac-Dirichlet at $s = d_{\max}$ | B | Verified via Apéry Q-linear independence of ζ(2), ζ(3) | Tier B: PSLQ-found, structurally constrained by Apéry. |

**Paper 18 verdict:** Already honestly framed. No edits required from
this audit. (Separate Paper 18 v2.0 restructure is a different sprint
per WH2 / §1.7.)

---

### 2.4 Paper 14 (qubit encoding) — `papers/core/paper_14_qubit_encoding.tex`

Brief check: Paper 14's headline scaling claims (O(Q^2.5), 51-1712×
sparsity advantage) are **measured power-law fits to real data**, not
PSLQ identifications. The closest thing to pattern-matching is the
"composed N_Pauli = 11.11 × Q across 28 molecules" claim. This is an
empirical regularity stable to ±0.1, not a derived formula. Honestly
framed as "isostructural invariance." **No edits required.**

---

### 2.5 Other papers

Surveyed via grep for PSLQ / coincidence / near-miss patterns:

- **Papers 7, 22, 24, 27** (core theorems): no Tier C content found.
- **Paper 29** (Ramanujan/Ihara): all results are exact integer
  polynomial factorizations or proven Ramanujan-bound checks. Tier A.
- **Paper 30** (SU(2) Wilson): exact maximal-torus identity, exact
  weak-coupling kinetic identity, Monte Carlo data points reported as
  data not as PSLQ targets. Tier A/data.
- **Paper 25** (Hopf gauge structure): synthesis paper; reframings are
  observational, not pattern-matched.
- **Paper 5, 4, 3** (conjectures folder): speculative by construction;
  appropriately labeled.

No additional Tier C overclaims found.

---

## 3. The c₂ situation in detail

Paper 28 §`curvature_coefficients` claims:

> The same spectral objects that compose K = π(B + F − Δ) also
> control the curvature dependence of the one-loop anomalous magnetic
> moment on S³.

The actual evidence:

| Aspect | Status |
|:-------|:-------|
| T9 theorem forces $c_2 = a + b\pi^2$ for rationals $a, b$ | Real (theorem) |
| $a = 19/100$, $b = -41/25200$ extracted from numerical $c_2^{\text{apparent}}$ | Real (8-digit match on 13-digit data) |
| The combination $(2 - B\Delta - F\Delta - F/B)/5$ equals $a + b\pi^2$ | Real (sympy identity, given the four numerator terms and denominator 5) |
| The combination was *predicted* before extraction | **No** — found by PSLQ-style search after $a, b$ were known |
| The same combination predicts $c_3$ | **No** — $c_3$ measured at $200\sigma$ nonzero, no fit in same basis |
| The combination derives from the curvature heat kernel | **No** — Schwinger-DeWitt derivation flagged as future work in §`curvature_coefficients` |

A T9-compatible $c_2 = a + b\pi^2$ has 2 real degrees of freedom. The
proposed form $(c_0 - c_1 B\Delta - c_2 F\Delta - c_3 F/B)/c_4$ with
small-integer coefficients in $\mathbb{Z}$ has roughly $\binom{20}{4}
\sim 5000$ candidate forms when each coefficient ranges over $\{-5,
\ldots, 5\}$. Finding *one* form that matches $a$ to high precision
and *another* (or the same one) that matches $b$ to high precision is
expected by chance at this search size and precision.

The candidates listed in `g2_c2_verification.json`:

| Candidate | Match digits |
|:----------|:------------:|
| cross-invariant $(2 - B\Delta - F\Delta - F/B)/5$ | 7.0 |
| $(6\pi^2 - 58)/7$ | 4.6 |
| $25/144$ | 2.7 |
| $7/40 = 7\Delta$ | 2.2 |
| $19/100$ | 1.0 |

The cross-invariant beats the next-best by ~2.4 digits. That's
suggestive but not conclusive at the search size. The killer is the
$c_3$ extension failure: a real identity should keep going.

**Recommendation:** demote in Paper 28 from "concrete bridge" to
"numerical coincidence between PSLQ identification of $c_2$ and Paper 2
invariants; structural mechanism would require deriving $c_2$ from
Schwinger-DeWitt curvature expansion of the spinor heat kernel on S³,
which has not been done."

---

## 4. The α situation in detail

The Paper 2 headline claim $K = \pi(B + F - \Delta)$ is:

- **Honestly framed in the paper** as conjectural.
- **Structurally decomposed**: each of B, F, Δ has a derived spectral home.
- **Statistically validated** at $p = 5.2 \times 10^{-9}$ in a 1.92e9 formula search.
- **Combination rule still undisclosed** by any first-principles argument.

The combination is in Tier C. The decomposition is Tier A.

The audit does **not** recommend removing Paper 2 from the canon. It
recommends:
1. Moving from `core/` to `observations/` (the headline is observational).
2. Tightening cross-paper language that treats the K formula as more
   than conjectural (notably P28-1 above, which builds on it).
3. Continuing to flag Tier C claims as Tier C in synthesis prose.

---

## 5. Recommended paper edits (for PI review before execution)

### 5.1 High priority

**EDIT-1: Paper 28 §`curvature_coefficients` — demote c₂ "Paper 2 bridge"**

Current text (line 936-945):
> Since α cancels exactly in the ratio F₂/[α/(2π)], the curvature
> coefficients c₁, c₂ are pure geometric invariants of the Dirac
> spectral sum on S³. The appearance of all three Paper 2 invariants
> (B, F, Δ) in c₂ provides a concrete numerical bridge between the α
> construction (Paper 2) and the QED vertex structure (this paper):
> the same spectral objects that compose K = π(B + F − Δ) also
> control the curvature dependence of the one-loop anomalous
> magnetic moment on S³.

Proposed replacement:
> Since α cancels in the ratio F₂/[α/(2π)], c₂ is a pure geometric
> invariant of the Dirac spectral sum on S³, and the T9 theorem
> constrains it to the form $a + b\pi^2$ with rationals $a, b$. The
> PSLQ identification (Eq. \ref{eq:c2_identification}) writes $a$ and
> $b$ as a small-integer combination of B, F, Δ, but this should be
> read with caution: (i) the form has more free coefficients than the
> two real numbers it explains, so finding *some* small-integer
> combination is expected by chance at this search size; (ii) the
> *same* form does not predict $c_3$, which is measured nonzero at
> $200\sigma$ (Eq. \ref{eq:c3_numerical}) and admits no fit in the
> extended T9 basis. A genuine structural identity should extend. The
> coincidence may reflect a deeper Schwinger-DeWitt origin for both B,
> F, Δ and c₂ on S³, but this would require deriving c₂ from the
> spinor heat kernel curvature expansion — not done in this work and
> flagged as an open question.

### 5.2 Medium priority

**EDIT-2: Paper 28 §`f2_convergence` — caveat the F₂ → 0 extrapolation**

Add a sentence noting that the $F_2 \sim 3.507 \cdot n^{-0.573}$ fit
is a 5-point power-law, the Richardson extrapolation gives an
unphysical $F_\infty \approx -0.22$, and the asymptotic claim should
be read as suggestive of slow convergence, not as a derived limit.

**EDIT-3: Paper 28 abstract**

Remove the c₂ formula from the abstract. The current abstract
(lines 47-53) treats c₂ as a headline result; given §3 above, c₂ should
appear in the body but not as a top-line claim.

### 5.3 Low priority / housekeeping

**EDIT-4: CLAUDE.md §2 frontier**

Replace the c₂ entry treating it as a "first concrete bridge between
Paper 2 and Paper 28" with the demoted framing.

**EDIT-5: CLAUDE.md §6 paper inventory**

Update Paper 2's location (`core/` → `observations/`) and adjust the
priority/loading description.

**EDIT-6: CLAUDE.md §13.5 Hard Prohibitions**

The current text refers to Paper 2 as Core. Update to reflect the move.

---

## 6. What this audit does NOT change

- The 18 symbolic proofs of Paper 7 — bulletproof.
- The angular sparsity theorem (Paper 22) — proven.
- The Bargmann-Segal π-free certificate (Paper 24) — proven.
- The Pauli scaling and Gaussian baseline comparisons (Paper 14) — measured.
- The Camporesi-Higuchi spectrum, the SD two-term exactness, the T9 theorem,
  the self-energy structural zero, the pendant-edge theorem — all theorems.

The audit narrows the framing of pattern-matched claims. It does not
weaken the structural results, which are the framework's actual load-bearing
content.

---

## 7. Suggested next actions (for PI)

1. Approve EDIT-1, EDIT-2, EDIT-3 (Paper 28 c₂ and F₂ caveats).
2. Approve EDIT-4, EDIT-5, EDIT-6 (CLAUDE.md sync).
3. Approve Paper 2 move to `observations/`.
4. Decide whether to formally retire the active α-derivation program
   (currently "paused" per WH5; this audit supports retirement until
   a structural lead emerges).
5. Decide whether the next refactor sprint targets the universal/Coulomb
   partition write-up (per the conversation that prompted this audit).
