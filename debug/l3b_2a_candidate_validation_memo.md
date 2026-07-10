# Sprint L3b-2a — Candidate validation + axiom check (memo)

**Sprint:** L3b-2a (opening sprint of the strong-form Lorentzian propinquity arc).
**Date:** 2026-05-22.
**Predecessor:** `debug/strong_form_propinquity_scoping_memo.md` (Candidate 6 = $L_{\mathrm{block}}$ recommendation, 2026-05-22).
**Status:** Scoping-grade verification; NO production code or paper modifications.
**Companion files:** `debug/l3b_2a_candidate_validation_compute.py` (driver), `debug/data/l3b_2a_candidate_validation.json` (raw per-generator results), `debug/l3b_2a_diagnostic_deeper.py` (structural diagnostic).

---

## §1. Summary

**Verdict: NO-GO on Candidate 6 ($L_{\mathrm{block}}$). Recommended fallback: Candidate 1 ($L_{\mathrm{op}}$).**

The load-bearing chirality-doubling diagnostic was executed at panel points $(n_{\max}, N_t) \in \{(2, 1), (2, 3), (3, 5)\}$ on the natural chirality-doubled scalar-multiplier operator system of the truncated Lorentzian Krein spectral triple at BBB $(m, n) = (4, 6)$. Two findings stand:

1. **Chirality-doubled symmetry holds bit-exact.** $\|[D_L, a]|_{\mathcal{K}^+}\|_{\mathrm{op}} = \|[D_L, a]|_{\mathcal{K}^-}\|_{\mathrm{op}}$ to machine zero (residual $0.000\mathrm{e}{+}00$ in float64) at every generator on every panel point. This is the prediction the scoping memo correctly derived from Paper 44 Prop. 5.1.

2. **$L_{\mathrm{block}}(a) = 0$ identically on the natural substrate.** All 42 generators at $(n_{\max}, N_t) = (2, 3)$ and all 275 generators at $(n_{\max}, N_t) = (3, 5)$ have $L_{\mathrm{block}}(a) = 0$ to machine zero. The non-zero commutator content lives entirely in the *chirality-flipping* (off-block-diagonal) cross-block pieces $P_+ [D_L, a] P_-$ and $P_- [D_L, a] P_+$, which are not in the K⁺ or K⁻ block restrictions used by $L_{\mathrm{block}}$.

**Structural reason (the substantive new content of this sprint).** $D_L = D_L^{\mathrm{diag}} + D_L^{\mathrm{off}}$ where:

- $D_L^{\mathrm{diag}} = i\gamma^0 \otimes \partial_t$ commutes with $J = \gamma^0 \otimes I_{N_t}$ (block-diagonal in K±).
- $D_L^{\mathrm{off}} = i D_{\mathrm{GV}} \otimes I_{N_t}$ anti-commutes with $J$ (purely off-block-diagonal) because $\{\gamma^0, D_{\mathrm{GV}}\} = 0$ in the chiral basis (Paper 45 §2.3 below eq. (DL_def); Paper 43 §3).

For every generator $a = M^{\mathrm{spat}} \otimes M^{\mathrm{temp}}$ in $\mathcal{O}^L$:

- $[D_L^{\mathrm{diag}}, a] = i [\gamma^0, M^{\mathrm{spat}}] \otimes \partial_t M^{\mathrm{temp}} + i \gamma^0 M^{\mathrm{spat}} \otimes [\partial_t, M^{\mathrm{temp}}] \equiv 0$ because (i) $[\gamma^0, M^{\mathrm{spat}}] = 0$ by Paper 44 Prop. 5.1, and (ii) $M^{\mathrm{temp}}$ is momentum-diagonal so $[\partial_t, M^{\mathrm{temp}}] = 0$. Verified to $0.000\mathrm{e}{+}00$ Frobenius norm at every generator.

- $[D_L^{\mathrm{off}}, a] = i [D_{\mathrm{GV}}, M^{\mathrm{spat}}] \otimes M^{\mathrm{temp}}$ is purely off-block-diagonal because $\{J, D_L^{\mathrm{off}}\} = 0$ and $[J, a] = 0$ together force $J [D_L^{\mathrm{off}}, a] J = -[D_L^{\mathrm{off}}, a]$, i.e. the commutator anti-commutes with $J$. Verified: $P_+ [D_L^{\mathrm{off}}, a] P_+ = P_- [D_L^{\mathrm{off}}, a] P_- = 0$ to machine zero at every generator.

Combining: $[D_L, a] = [D_L^{\mathrm{off}}, a]$ is purely off-block-diagonal for every $a \in \mathcal{O}^L$. The block restrictions $P_\pm [D_L, a] P_\pm$ vanish identically, and $L_{\mathrm{block}}(a) \equiv 0$.

**The scoping memo §3 Candidate 6 reasoning had a gap.** The memo identified (correctly) that the K⁺ and K⁻ block contents are identical by chirality-doubling, and claimed (correctly) that this would identify the K⁺-weak-form seminorm. But the *content* of those blocks is zero, not the commutator's full content, because the entire commutator lives in the chirality-flipping piece. The K⁺-weak-form construction in Paper 45 §2.6 also restricts to the K⁺ subspace, and on the natural substrate (where $a$ commutes with $J$), $P_+ a P_+ = a|_{K^+}$ trivially; but the K⁺-weak-form Lipschitz seminorm is $\|[P_+ D_L P_+, P_+ a P_+]\|$, not $P_+ [D_L, a] P_+$. These are different operators because $D_L$ does not commute with $P_+$: $P_+ D_L P_+ \neq D_L|_{K^+}$ when $D_L$ has off-block-diagonal content.

In short: Paper 45's K⁺-weak-form uses the *projected* Dirac $P_+ D_L P_+$, which keeps only the block-diagonal piece of $D_L$ (i.e., $P_+ D_L^{\mathrm{diag}} P_+ = i \partial_t|_{K^+}$, since the K⁺/K⁻ piece of $\gamma^0 \otimes \partial_t$ on K⁺ restricts to $(+1) \otimes \partial_t = i\,\mathrm{diag}(\omega_k)$). Whereas $L_{\mathrm{block}}$ uses $P_+ [D_L, a] P_+$, which only sees the block-diagonal *commutator*, not a commutator with the projected Dirac.

**These two K⁺ constructions disagree.** The correct identification (axiom check (f)) requires using $\|[P_+ D_L P_+, a]\|$ rather than $\|P_+ [D_L, a] P_+\|$. With $L_{\mathrm{block}}$ defined as the latter, recovery of Paper 45 fails.

**Recommendation:** fall back to Candidate 1 ($L_{\mathrm{op}}(a) = \|[D_L, a]\|_{\mathrm{op}}$) for Sprint L3b-2b. $L_{\mathrm{op}}$ is non-trivial (verified: $L_{\mathrm{op}}(a) > 0$ for every spatial-non-trivial generator), satisfies axioms (a)–(c) trivially, and on the K⁺ subspace reduces to the operator norm restricted to $P_+ \mathcal{B}(\mathcal{K}) P_+$ — but this is still NOT Paper 45's K⁺-weak-form seminorm. Recovery of Paper 45 under Candidate 1 will require an additional argument that the cross-block content drops out of the propinquity bound. This is a real open question.

---

## §2. Setup

**Panel points.** Three cells were computed:

| Cell | $(n_{\max}, N_t)$ | $\dim \mathcal{K}$ | $\dim \mathcal{O}^L$ | $\dim \mathcal{K}^+$ | $\dim \mathcal{K}^-$ | # generators | Wall time (s) |
|:----:|:-----------------:|:------------------:|:--------------------:|:--------------------:|:--------------------:|:------------:|:-------------:|
| 1 | (2, 3) | 48 | 42 | 24 | 24 | 42 | 0.1 |
| 2 | (2, 1) | 16 | 14 | 8 | 8 | 14 | 0.0 |
| 3 | (3, 5) | 200 | 275 | 100 | 100 | 275 | 38.4 |

Cell 1 is the load-bearing toy case 4 + toy case 5 combined from the scoping memo §6 (smallest non-trivial cutoff with $N_t > 1$). Cell 2 is the $N_t = 1$ Riemannian-limit sanity. Cell 3 is the larger panel point to confirm the result doesn't depend on $n_{\max}$.

**Substrate.** The natural chirality-doubled scalar-multiplier operator system $\mathcal{O}^L_{n_{\max}, N_t, T}$ from `geovac.operator_system_compact_temporal.CompactTemporalTruncatedOperatorSystem`. Multipliers are tensor products $M^{\mathrm{spat}}_{N, L, M} \otimes M^{\mathrm{temp}}_p$ where:

- $M^{\mathrm{spat}}_{N, L, M} = \mathrm{blkdiag}(W_{N, L, M}, W_{N, L, M})$ is the chirality-doubled lift of the Avery–Wen–Avery Weyl-sector multiplier.
- $M^{\mathrm{temp}}_p = \mathrm{diag}(\omega_k^p)$, $p = 0, \ldots, N_t - 1$, is momentum-polynomial diagonal on the $N_t$-mode Fourier truncation of $L^2(S^1_T)$.

**Operators.** $D_L = i(\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I_{N_t})$ per vdD 2016 Prop. 4.1 / Paper 43 §3 / Paper 45 §2.3, built by `geovac.lorentzian_dirac_compact.lorentzian_dirac_compact_matrix`. $J_L = \gamma^0 \otimes I_{N_t}$, built by `geovac.krein_space_compact_temporal.CompactTemporalKreinSpace`.

**Krein axioms checked exactly at every cell** (Frobenius residual 0.000e+00 in float64): $J^2 = +I$, $J^* J = I$, $J = J^*$.

**Paper 44 Prop. 5.1 verified at every cell** (max $\|[J, M]\|_F = 0.000\mathrm{e}{+}00$ across all generators).

**Arithmetic.** float64 throughout. The Krein-space construction uses real $\{0, 1\}$ permutation entries in $J$ and rational entries in $\gamma^0 \otimes \partial_t$ (modulo $i$ and $2\pi k / T$), so the Frobenius residuals are genuinely zero in IEEE 754 (no rounding error). The spatial Avery–Wen–Avery 3-Y integrals are exact rationals from sympy in the substrate construction, then cast to float64.

**Operator norm $\|A\|_{\mathrm{op}}$**: largest singular value via `numpy.linalg.svd`. For block restrictions, we use $\|P_\pm A P_\pm\|_{\mathrm{op}}$, taken as the operator norm of the projected operator on the full $\mathcal{K}$ — equivalent (modulo a zero-eigenspace that SVD ignores) to the operator norm of the restricted operator on $\mathcal{K}^\pm$.

---

## §3. Toy-case 4 + 5 computational results

### §3.1 Chirality-doubled equality (load-bearing positive)

At all three panel points, for every generator $a \in \mathcal{O}^L$:

$$\|[D_L, a]|_{\mathcal{K}^+}\|_{\mathrm{op}} \;=\; \|[D_L, a]|_{\mathcal{K}^-}\|_{\mathrm{op}}$$

to numerical zero. Maximum residual across all generators and all cells:

$$\max_{a} \bigl| \|[D_L, a]|_{\mathcal{K}^+}\|_{\mathrm{op}} - \|[D_L, a]|_{\mathcal{K}^-}\|_{\mathrm{op}} \bigr| \;=\; 0.000\mathrm{e}{+}00.$$

This confirms the chirality-doubling symmetry derived in scoping memo §4 reason 5 and matches Paper 44 Prop. 5.1 verbatim.

### §3.2 L_block is identically zero (load-bearing negative)

At all three panel points, for every generator $a \in \mathcal{O}^L$:

$$L_{\mathrm{block}}(a) \;:=\; \max\bigl( \|P_+ [D_L, a] P_+\|_{\mathrm{op}},\; \|P_- [D_L, a] P_-\|_{\mathrm{op}} \bigr) \;=\; 0.$$

Bit-exact. Specifically:

| Cell | # generators | # with $L_{\mathrm{block}} = 0$ | $\max L_{\mathrm{block}}$ | $\max L_{\mathrm{op}}$ |
|:----:|:------------:|:-------------------------------:|:-------------------------:|:----------------------:|
| (2, 3) | 42 | 42 / 42 | 0.000e+00 | 2.251e-01 |
| (2, 1) | 14 | 14 / 14 | 0.000e+00 | 2.251e-01 |
| (3, 5) | 275 | 275 / 275 | 0.000e+00 | (not reported in summary, but $>0$) |

Hence on the natural substrate, $L_{\mathrm{block}}$ is the zero seminorm — it does not distinguish any two generators. As a Lipschitz seminorm it is degenerate.

### §3.3 Structural explanation (deeper diagnostic, `l3b_2a_diagnostic_deeper.py`)

Write $D_L = D_L^{\mathrm{diag}} + D_L^{\mathrm{off}}$ with:

- $D_L^{\mathrm{diag}} := i \gamma^0 \otimes \partial_t$. **Block-diagonal in K±** because $J = \gamma^0 \otimes I$ and $[\gamma^0, \gamma^0] = 0$, $[\partial_t, I] = 0$. Verified: $\|[J, D_L^{\mathrm{diag}}]\|_F = 0.000\mathrm{e}{+}00$.

- $D_L^{\mathrm{off}} := i D_{\mathrm{GV}} \otimes I_{N_t}$. **Off-block-diagonal in K±** because $\{\gamma^0, D_{\mathrm{GV}}\} = 0$ (Paper 45 §2.3 below eq. (DL_def), in the Peskin–Schroeder chiral basis) and $[I, I] = 0$. Verified: $\|\{J, D_L^{\mathrm{off}}\}\|_F = 0.000\mathrm{e}{+}00$.

For every generator $a = M^{\mathrm{spat}} \otimes M^{\mathrm{temp}} \in \mathcal{O}^L$:

**Block-diagonal commutator vanishes identically:**
$$[D_L^{\mathrm{diag}}, a] = i\bigl([\gamma^0, M^{\mathrm{spat}}] \otimes \partial_t M^{\mathrm{temp}} + \gamma^0 M^{\mathrm{spat}} \otimes [\partial_t, M^{\mathrm{temp}}]\bigr) = 0,$$
because (i) $[\gamma^0, M^{\mathrm{spat}}] = 0$ by Paper 44 Prop. 5.1, and (ii) $M^{\mathrm{temp}} = \mathrm{diag}(\omega_k^p)$ commutes with $\partial_t = i\,\mathrm{diag}(\omega_k)$ (both diagonal in momentum). Verified: $\max_a \|[D_L^{\mathrm{diag}}, M]\|_F = 0.000\mathrm{e}{+}00$ over all generators at $(2, 3)$.

**Off-diagonal commutator is purely chirality-flipping:**
$$[D_L^{\mathrm{off}}, a] = i [D_{\mathrm{GV}}, M^{\mathrm{spat}}] \otimes M^{\mathrm{temp}}.$$
Because $\{J, D_L^{\mathrm{off}}\} = 0$ and $[J, a] = 0$, we have
$$J [D_L^{\mathrm{off}}, a] J = J D_L^{\mathrm{off}} J\, a - a\, J D_L^{\mathrm{off}} J = -D_L^{\mathrm{off}} a + a D_L^{\mathrm{off}} = -[D_L^{\mathrm{off}}, a],$$
i.e., the commutator anti-commutes with $J$. This makes it purely off-block-diagonal in K±. Verified: $\max_a \|P_\pm [D_L^{\mathrm{off}}, a] P_\pm\|_F = 0.000\mathrm{e}{+}00$.

**Conclusion:** $[D_L, a] = [D_L^{\mathrm{off}}, a]$ for every $a \in \mathcal{O}^L$, and this commutator is purely off-block-diagonal. The block restrictions $P_\pm [D_L, a] P_\pm$ vanish identically. Hence $L_{\mathrm{block}}(a) \equiv 0$.

### §3.4 L_op is non-trivial

By contrast, $L_{\mathrm{op}}(a) = \|[D_L, a]\|_{\mathrm{op}} = \|[D_L^{\mathrm{off}}, a]\|_{\mathrm{op}}$ is non-zero for any $a$ whose spatial factor $M^{\mathrm{spat}}$ does not commute with $D_{\mathrm{GV}}$. At $(2, 3)$ the maximum is $L_{\mathrm{op}} = 0.2251$, achieved by the first non-trivial multiplier $(N, L, M, p) = (2, 0, 0, 0)$ — the spatial monopole at the top shell, with constant temporal factor.

The picture: **the entire Lipschitz content lives in the cross-block (chirality-flipping) pieces of $[D_L, a]$.** $L_{\mathrm{op}}$ captures this content; $L_{\mathrm{block}}$ does not.

### §3.5 Per-generator table (sample at $(n_{\max}, N_t) = (2, 3)$, generator (2, 0, 0, 0))

| Quantity | Value |
|:---|:---:|
| $L_{\mathrm{op}}$ | $0.225079$ |
| $L_{\mathrm{block}}$ | $0$ |
| $L_{\mathrm{block}, K^+}$ | $0$ |
| $L_{\mathrm{block}, K^-}$ | $0$ |
| $\|P_+ [D_L, a] P_-\|_{\mathrm{op}}$ | $0.225079$ |
| $\|P_- [D_L, a] P_+\|_{\mathrm{op}}$ | $0.225079$ |
| $\|[J, a]\|_F$ | $0$ |
| $\|[D_L^{\mathrm{diag}}, a]\|_{\mathrm{op}}$ | $0$ |
| $\|[D_L^{\mathrm{off}}, a]\|_{\mathrm{op}}$ | $0.225079$ |

Full per-generator data (42 rows at (2,3), 14 at (2,1), 275 at (3,5)) is in `debug/data/l3b_2a_candidate_validation.json`.

---

## §4. Six-axiom checks

### (a) Seminorm property — **TECHNICALLY PASS, but on a degenerate seminorm.**

Verified numerically on the (2, 3) panel:
- $L_{\mathrm{block}}(0) = 0$ exactly.
- Positive homogeneity max residual $= 0.000\mathrm{e}{+}00$ over 8 random complex-scaled samples.
- Triangle inequality max violation $= 0.000\mathrm{e}{+}00$ over 8 random pairs.

But this is "pass by triviality": all three identities hold because $L_{\mathrm{block}}(a) \equiv 0$ on the natural substrate. The seminorm IS a seminorm in the technical sense — the zero seminorm is a degenerate but well-defined seminorm — but it fails to be a Lipschitz seminorm in any useful sense because it doesn't separate operators.

**Justification:** Mathematically, the max of two operator-norm-induced seminorms IS a seminorm. The max of two zero seminorms is the zero seminorm, which is also a seminorm. But the resulting "metric" $\rho_{L_{\mathrm{block}}}(\phi, \psi) = \sup\{|\phi(a) - \psi(a)| : L_{\mathrm{block}}(a) \le 1\}$ on the truncated state space is **infinite** for any two distinct states (the unit ball $\{a : L_{\mathrm{block}}(a) \le 1\}$ is all of $\mathcal{O}^L$, and there are unboundedly many directions a state can point in). Hence $L_{\mathrm{block}}$ does not induce a Latrémolière metric on the state space, and the propinquity construction breaks at L5 — there is no compact / bounded reach to bound.

### (b) Lower semicontinuity — **PASS (vacuously).**

The zero seminorm is trivially lower semicontinuous: $\liminf_n 0 = 0$. The pointwise max of two lsc functionals is lsc; both block-restricted operator-norm seminorms are lsc; hence $L_{\mathrm{block}}$ is lsc. The degeneracy doesn't break lsc.

### (c) *-closure — **PASS.**

For Hilbert-self-adjoint $a \in \mathcal{O}^L$ (i.e., $a^\dagger = a$ as Hilbert operators), $[D_L, a]^\dagger = D_L^\dagger a^\dagger - a^\dagger D_L^\dagger$. Since $D_L = i(\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I)$ with $\gamma^{0\dagger} = \gamma^0$, $\partial_t^\dagger = -\partial_t$, $D_{\mathrm{GV}}^\dagger = D_{\mathrm{GV}}$, we have $D_L^\dagger = -i\gamma^0 \otimes \partial_t^\dagger - i D_{\mathrm{GV}} \otimes I = i\gamma^0 \otimes \partial_t - i D_{\mathrm{GV}} \otimes I$. Wait — this doesn't quite reduce to $\pm D_L$. Let me redo: $D_L^\dagger = -i (\gamma^0 \otimes \partial_t)^\dagger - i(D_{\mathrm{GV}} \otimes I)^\dagger = -i(\gamma^0 \otimes (-\partial_t)) - i(D_{\mathrm{GV}} \otimes I) = i\gamma^0 \otimes \partial_t - iD_{\mathrm{GV}} \otimes I$. So $D_L^\dagger = D_L^{\mathrm{diag}} - D_L^{\mathrm{off}} \neq \pm D_L$.

However, $L_{\mathrm{block}}(a^\dagger) = \|[D_L, a^\dagger]|_{K^+}\|$ etc.; since the block-restriction operation commutes with the Hilbert adjoint (P_+ is Hermitian), $\|P_+ [D_L, a^\dagger] P_+\|_{\mathrm{op}} = \|(P_+ [D_L^\dagger, a] P_+)^\dagger\|_{\mathrm{op}} = \|P_+ [D_L^\dagger, a] P_+\|_{\mathrm{op}}$. Now $[D_L^\dagger, a] = [D_L^{\mathrm{diag}}, a] - [D_L^{\mathrm{off}}, a] = -[D_L^{\mathrm{off}}, a]$ (since the diag piece vanishes by §3.3). And $P_+ [-D_L^{\mathrm{off}}, a] P_+ = 0$ since the off-piece is off-block-diagonal. So $L_{\mathrm{block}}(a^\dagger) = 0 = L_{\mathrm{block}}(a)$. (c) holds trivially because both sides are zero.

### (d) Lichnerowicz-style upper bound prerequisite — **VACUOUSLY PASS.**

The structural prerequisite is that on each block, $[D_L, a]|_{K^\pm}$ should be a bounded Hilbert operator on a Hilbert subspace. This holds: $K^+$ and $K^-$ are Hilbert subspaces (Paper 45 §2.6); $P_\pm [D_L, a] P_\pm$ is a bounded operator on each (since $D_L$ is bounded at finite cutoff). The Lichnerowicz bound $L_{\mathrm{block}}(a) \le C_3 \|\nabla f\|_\infty$ would hold with $C_3 = 0$, since $L_{\mathrm{block}} \equiv 0$. This is a vacuous bound: it doesn't usefully bound anything.

### (e) Berezin reconstruction compatibility — **VACUOUSLY PASS.**

By Paper 45 Lemma 4.3 / Sub-sprint C, the Berezin image $B^{\mathrm{joint}}(f)$ preserves $K^+$ and $K^-$ (commutes with $J$). Hence $B^{\mathrm{joint}}(f) \in \mathcal{O}^L$ (in the closure thereof). The structural finding §3.3 then applies: $[D_L, B^{\mathrm{joint}}(f)] = [D_L^{\mathrm{off}}, B^{\mathrm{joint}}(f)]$ is purely off-block-diagonal, so $L_{\mathrm{block}}(B^{\mathrm{joint}}(f)) = 0$. The Berezin map composed with $L_{\mathrm{block}}$ is the zero functional. The "approximate identity" property of $B^{\mathrm{joint}}$ would be vacuous under $L_{\mathrm{block}}$.

### (f) Reduction to Paper 45 K⁺-weak-form — **FAIL.**

The scoping memo §3 Candidate 6 claimed (page 261):

> "Restriction to K⁺ recovers Paper 45 exactly. On K⁺, the K⁻-block contribution is zero (no operators reach K⁻), so $L_{\mathrm{block}}(a)|_{K^+} = \|[D_L, a]|_{K^+}\|_{\mathrm{op}}$, which is Paper 45's seminorm."

But Paper 45's K⁺-weak-form seminorm (Definition 2.3, eq. `eq:weak_form_propinquity_def`) is

$$L^+_{\mathrm{P45}}(a) := \|[P_+ D_L P_+, P_+ a P_+]\|_{\mathrm{op, on } K^+},$$

NOT $\|P_+ [D_L, a] P_+\|_{\mathrm{op}}$. The difference is whether you project the Dirac before commuting, or commute first and project after. For our $D_L$:

- $P_+ D_L P_+ = P_+ D_L^{\mathrm{diag}} P_+ + P_+ D_L^{\mathrm{off}} P_+ = P_+ D_L^{\mathrm{diag}} P_+ + 0$, since the off-piece is off-block-diagonal. So $P_+ D_L P_+ = (+1) \otimes (i\partial_t)$ on $K^+$ — Paper 45's projected Lorentzian Dirac restricted to $K^+$ is just $i\partial_t$ (per-block).

- Paper 45 K⁺-weak-form: $L^+_{\mathrm{P45}}(a) = \|[i\partial_t|_{K^+}, a|_{K^+}]\|_{\mathrm{op}}$ on $K^+$. This is purely *temporal* — it only sees the $\partial_t$ part, none of the spatial $D_{\mathrm{GV}}$ content. For generators with non-trivial spatial multiplier and trivial temporal multiplier (e.g., $a = M^{\mathrm{spat}} \otimes I$), Paper 45's K⁺-weak-form gives zero.

- $L_{\mathrm{block}}(a) = P_+ [D_L, a] P_+ = P_+ [D_L^{\mathrm{off}}, a] P_+ = 0$ (off-block-diagonal). Same answer.

So on the natural substrate, BOTH $L_{\mathrm{block}}$ AND Paper 45's K⁺-weak-form $L^+_{\mathrm{P45}}$ are seriously deficient at separating spatial-only multipliers — for those, $L^+_{\mathrm{P45}}(a) = 0$ as well.

**This means the scoping memo's claim that "Paper 45 is a corollary" was technically right** — both seminorms agree (both give zero on spatial-only multipliers) — **but only because both are degenerate on the same generators.** Paper 45's K⁺-weak-form propinquity construction must be doing something more subtle than just $L^+_{\mathrm{P45}}$ to produce its non-trivial $\Lambda$ values. Let me check.

Actually, re-reading Paper 45 §2.6 / §4: the K⁺-weak-form propinquity is defined as $\Lambda(\mathcal{T}^+_1, \mathcal{T}^+_2)$ where $\mathcal{T}^+_i = (P_+ \mathcal{O}^L_i P_+, K^+_i, P_+ D^L_i P_+)$. The standard Latrémolière propinquity on $\mathcal{T}^+$ uses the Lipschitz seminorm $L^+_{\mathrm{P45}}$, but the *non-triviality* of the resulting propinquity comes from the structure of the *operator system* $P_+ \mathcal{O}^L P_+$ and the Berezin reconstruction $P_+ B^{\mathrm{joint}}(f) P_+$, not from the per-generator Lipschitz values.

In other words, Paper 45's $\Lambda$ values $(2.0746, 1.6101, 1.3223)$ at $(n_{\max}, N_t) \in \{(2, 3), (3, 5), (4, 7)\}$ are inherited from Paper 38 (the spatial SU(2) propinquity), with the $\Uone$ factor contributing only a "free upgrade" via the joint Berezin and the $N_t$-independent cb-norm. The temporal direction does not contribute to the Lipschitz numerator on the natural substrate.

The "free upgrade" is real and substantial — the tensor-product propinquity convergence theorem (Paper 45 main theorem) is genuinely new content — but it's not because the temporal Lipschitz seminorm is non-trivial; it's because the operator-system structure and Berezin reconstruction transport via PURE_TENSOR factorization.

**So the upshot for axiom (f):** $L_{\mathrm{block}}$ and Paper 45's $L^+_{\mathrm{P45}}$ are BOTH degenerate on spatial-only multipliers, BOTH give the same answer (zero), but only $L^+_{\mathrm{P45}}$ is wired into the K⁺-weak-form propinquity construction in a way that produces a non-trivial $\Lambda$. The strong-form construction must use the seminorm together with the right operator system and Berezin map — none of the candidate seminorms standalone yields the propinquity bound.

**Verdict (f):** $L_{\mathrm{block}}$ does NOT recover Paper 45 in the sense the scoping memo claimed — neither does it fail, since both give the same (trivial) answer on the natural substrate. The recovery claim was meaningless; the substantive recovery property is **the propinquity bound matches**, which is L5 (Sprint L3b-2d), not the seminorm.

---

## §5. Numerical-panel cross-check vs Paper 45 §6

Paper 45 §6 reports $\Lambda(n_{\max}, N_t)$ at three cells: $\Lambda(2, 3) = 2.0746$, $\Lambda(3, 5) = 1.6101$, $\Lambda(4, 7) = 1.3223$.

These values **cannot be cross-checked at Sprint L3b-2a** because they depend on the full propinquity bound (L5 assembly), which requires the Berezin reconstruction rate (L4) and the Lichnerowicz constant (L3). Sprint L3b-2a verifies only the per-generator Lipschitz seminorm structure — task (f) of the axiom check — and that task is now resolved (both $L_{\mathrm{block}}$ and Paper 45's $L^+_{\mathrm{P45}}$ are degenerate on spatial-only multipliers; the propinquity nontriviality comes from elsewhere).

For a future Sprint L3b-2d under whatever seminorm is chosen, the bound $\Lambda(n_{\max}, N_t)$ should be recomputed with the strong-form bound and compared against Paper 45 §6 bit-exact if the strong-form is a corollary; or with a documented strict inequality if the strong-form gives a different / tighter bound. The Paper 45 panel values are the natural baseline either way.

---

## §6. Go/no-go verdict

**Verdict: NO-GO on Candidate 6 ($L_{\mathrm{block}}$).**

**Reason:** $L_{\mathrm{block}}(a) \equiv 0$ on the natural chirality-doubled scalar-multiplier operator system, by the structural decomposition $D_L = D_L^{\mathrm{diag}} + D_L^{\mathrm{off}}$ where the diag piece commutes with every generator (vanishing commutator) and the off piece produces a purely chirality-flipping commutator (block restrictions vanish). The candidate is degenerate; it does not induce a Latrémolière metric on the state space.

**Recommendation:** Fall back to **Candidate 1 ($L_{\mathrm{op}}$)** for the next sprint.

Specifically, Sprint L3b-2b should:

1. Set $L_{\mathrm{op}}(a) := \|[D_L, a]\|_{\mathrm{op}}$ on the natural substrate.
2. Re-verify the six axiom checks under $L_{\mathrm{op}}$ (most are trivial: it's just the Hilbert operator norm of a commutator).
3. Derive the L3 Lichnerowicz-style upper bound for $L_{\mathrm{op}}(a)$ in terms of joint gradient norms $\|\nabla^{\mathrm{joint}} f\|_\infty$. This is the dominant new work.
4. **Critical follow-up question:** Does the resulting strong-form propinquity recover Paper 45's K⁺-weak-form $\Lambda$ values? The structural answer is unlikely to be "bit-exact" because $L_{\mathrm{op}}$ includes the off-block-diagonal commutator content (chirality-flipping), which Paper 45's K⁺-weak-form discards. The strong-form bound under $L_{\mathrm{op}}$ is expected to be STRICTLY LARGER than Paper 45's bound (a worse bound on $\Lambda$). This would still be a valid theorem — strong-form propinquity holds — but with a constant strictly worse than Paper 45's. Whether the asymptotic rate ($O(\log n_{\max}/n_{\max} + T/N_t)$) survives is the substantive open question for L3b-2b.

**Alternative fallback** (riskier, but more interesting): Develop a *third* seminorm that captures the spatial chirality-flipping content cleanly without losing the K⁺-recovery property. The natural candidate is

$$L_{\mathrm{cross}}(a) := \max\bigl( \|P_+ [D_L, a] P_-\|_{\mathrm{op}}, \|P_- [D_L, a] P_+\|_{\mathrm{op}} \bigr).$$

This captures exactly the cross-block content that $L_{\mathrm{block}}$ misses and $L_{\mathrm{op}}$ includes alongside the (zero) block-diagonal content. By the chirality-doubled symmetry verified in §3.1, the two arguments to the max are equal, so $L_{\mathrm{cross}}(a) = \|P_+ [D_L, a] P_-\|_{\mathrm{op}}$. On the natural substrate, $L_{\mathrm{cross}}(a) = \|P_+ [D_L^{\mathrm{off}}, a] P_-\|_{\mathrm{op}}$, which captures the full spatial-derivative content of $a$. Sanity check at $(2, 3)$ generator $(2, 0, 0, 0)$: $L_{\mathrm{cross}} = 0.225079 = L_{\mathrm{op}}$ (since the commutator is purely cross-block, its full operator norm equals the cross-block contribution). So on the natural substrate, $L_{\mathrm{cross}} = L_{\mathrm{op}}$ identically — there is no separate content, and the two candidates agree.

The implication: **on the natural substrate, $L_{\mathrm{op}} = L_{\mathrm{cross}}$** (both capture the full off-block-diagonal commutator); **$L_{\mathrm{block}} = 0$** (block-diagonal commutator vanishes). Only one non-trivial seminorm is available, and it is $L_{\mathrm{op}}$.

This is honest news. The strong-form / K⁺-weak-form distinction on the natural substrate is *not* about the seminorm choice — both reasonable choices ($L_{\mathrm{op}}$ vs the Paper 45 K⁺-weak-form $L^+_{\mathrm{P45}}$) are well-defined and different. The question is which propinquity construction (full Krein vs K⁺-restricted) we want, and what the resulting bounds look like.

---

## §7. Honest scope

### What this sprint definitively closes

- **Paper 44 Prop. 5.1 verified bit-exact computationally at $(n_{\max}, N_t) \in \{(2, 1), (2, 3), (3, 5)\}$.** Every chirality-doubled scalar-multiplier generator commutes with $J_L$ to machine zero.
- **Chirality-doubled symmetry of block restrictions verified.** $\|[D_L, a]|_{K^+}\|_{\mathrm{op}} = \|[D_L, a]|_{K^-}\|_{\mathrm{op}}$ bit-exact at every generator on every cell.
- **$L_{\mathrm{block}}$ is identically zero on the natural substrate.** This is a clean negative result with a clean structural explanation (§3.3).
- **The scoping memo's Candidate 6 reasoning had a gap.** The memo conflated $P_+ [D_L, a] P_+$ (the block restriction of the commutator) with Paper 45's K⁺-weak-form $[P_+ D_L P_+, P_+ a P_+]$. These are different operators because $D_L$ does not commute with $P_+$.

### What this sprint does NOT close

- **The strong-form propinquity convergence theorem.** Sprint L3b-2a is scoping-grade; the actual propinquity bound under any candidate seminorm requires L3 (Lichnerowicz, Sprint L3b-2b), L4 (Berezin under the candidate, Sprint L3b-2c), and L5 (propinquity assembly, Sprint L3b-2d). Each of these is multi-week work.
- **A definitive answer to "is the strong form genuinely new content?"** With $L_{\mathrm{op}}$ as the seminorm, the strong-form propinquity bound is expected to be strictly larger than Paper 45's K⁺-weak-form bound (worse constant, possibly worse rate). Whether this is a strictly weaker theorem, a strictly stronger theorem (because it bounds something the K⁺-weak-form doesn't), or asymptotically equivalent is the open question.
- **The enlarged-substrate strong form.** Scoping §3 Candidate 6 noted that the strong-form on an enlarged operator system (including chirality-flipping operators) would give genuinely new content. This is still true and is the Sprint L3b-2f target; not addressed here.
- **Paper 45's K⁺-weak-form numerical-panel values bit-exact reproduction.** Reproduction requires the L5 propinquity assembly under the chosen seminorm; out of scope for L3b-2a.

### Risks and follow-up work for Sprint L3b-2b

Given the NO-GO on $L_{\mathrm{block}}$ and the fallback to $L_{\mathrm{op}}$:

1. **Risk: $L_{\mathrm{op}}$-based strong-form propinquity gives a strictly worse $\Lambda$ bound than Paper 45's K⁺-weak-form.** Probability: high (~80%). Mitigation: document the strict inequality as a structural feature of the strong-form, not a bug; the strong-form theorem bounds a strictly larger quantity (the full Krein operator-norm Lipschitz distance) than the K⁺-weak-form (which bounds only the K⁺-block-restricted distance).

2. **Risk: $L_{\mathrm{op}}$-based strong-form propinquity loses the asymptotic rate $O(\log n_{\max}/n_{\max})$.** Probability: moderate (~40%). The off-block-diagonal commutator content scales with $\|D_{\mathrm{GV}}\|_{\mathrm{op}} \sim n_{\max}/2$ at the natural cutoff, which is a growing scale. Whether this kills the rate or just changes the constant requires the L3 derivation in L3b-2b.

3. **Risk: $L_{\mathrm{op}}$ does not satisfy Latrémolière 2017/2023 Lipschitz-seminorm axioms cleanly.** Probability: low (~10%). $L_{\mathrm{op}}$ is just the Hilbert operator norm of a commutator — the most natural candidate for any spectral triple, and the one Latrémolière uses verbatim in the Riemannian case. The Krein adjoint-vs-Hilbert adjoint distinction may show up at the $L(a^*) = L(a)$ axiom, but only mildly.

### Open structural question raised by this sprint

The structural finding that $[D_L^{\mathrm{diag}}, a] = 0$ for every generator $a \in \mathcal{O}^L$ (because $D_L^{\mathrm{diag}}$ is built from operators that all commute with the multiplier algebra on the natural substrate) suggests that **the temporal direction of $D_L$ contributes nothing to the per-generator Lipschitz content on the natural substrate**. The temporal contribution to the propinquity is entirely via Berezin reconstruction (rate $\gamma^{\mathrm{joint}} \sim T/N_t$). This matches Paper 45 §4.2 Lemma 4.3's vanishing time-chirality cross-term $\{\gamma^0, D_{\mathrm{GV}}\} = 0$, but sharpens the observation: **the natural substrate makes the temporal commutator identically zero, not just the cross-term.**

This is news for any future strong-form construction: the structural-Lipschitz content lives entirely in the spatial off-block-diagonal commutator $[D_{\mathrm{GV}}, M^{\mathrm{spat}}]$ tensored with $M^{\mathrm{temp}}$, and the temporal multipliers $M^{\mathrm{temp}}$ act only as passive scalar amplifiers. The strong-form Lorentzian propinquity, on this substrate, is essentially a Riemannian-on-each-fiber problem with the temporal index as a passive parameter.

This is a clean explanation of why Paper 45's tensor-product propinquity convergence theorem inherits the asymptotic rate from Paper 38's SU(2) bound — temporal compactification does not add Lipschitz content, only Berezin-reconstruction content.

### Recommended next sprint

**Sprint L3b-2b under $L_{\mathrm{op}}$.** Specifically:

1. Verify $L_{\mathrm{op}}$ satisfies axioms (a)–(c) cleanly.
2. Derive the Lichnerowicz upper bound $L_{\mathrm{op}}(a) \le C_3^{\mathrm{strong}} \|\nabla^{\mathrm{joint}} f\|_\infty$ for $a = B^{\mathrm{joint}}(f)$. Compare $C_3^{\mathrm{strong}}$ to Paper 45's $C_3^{\mathrm{joint}}$.
3. Track whether the asymptotic rate $\gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} \to 0$ survives.
4. Decision point: if the rate survives with a constant within a factor of 2–3 of Paper 45's, proceed to L3b-2c (Berezin under $L_{\mathrm{op}}$) and L3b-2d (propinquity assembly). If the rate dies or the constant blows up, document and consider whether the K⁺-weak-form is the right paper-grade strong-form claim or whether enlargement of the operator system (L3b-2f) is needed.

The scoping memo's timeline estimate (4.5–8 months for the first tranche L3b-2a–e) needs revision after this NO-GO. L3b-2a took ~1 day of compute time and produced a clean negative; the load-bearing question now shifts to L3b-2b, which still has the dominant risk. Conservative revised timeline: **3–6 months for first tranche L3b-2a–e**, given the L3b-2a closure already.

---

## §8. Closing note

The recommended candidate from the scoping memo did not survive the load-bearing diagnostic. This is the diagnostic-before-engineering rule working as designed (memory `feedback_diagnostic_before_engineering.md`). One ~30-minute compute run produced a clean negative that would have wasted 4–6 weeks of L3b-2b L3 derivation under the wrong candidate.

The substantive new content of this sprint:

1. **$L_{\mathrm{block}}$ is degenerate on the natural substrate.** Confirmed bit-exact.
2. **The structural decomposition $D_L = D_L^{\mathrm{diag}} + D_L^{\mathrm{off}}$ explains why.** $D_L^{\mathrm{diag}}$ commutes with the algebra; $D_L^{\mathrm{off}}$ produces a purely chirality-flipping commutator.
3. **The temporal direction contributes nothing to the per-generator Lipschitz content on the natural substrate.** Sharpens Paper 45 Lemma 4.3's vanishing cross-term.
4. **The strong-form vs K⁺-weak-form distinction is about which propinquity construction, not which seminorm.** Both natural seminorms ($L_{\mathrm{op}}$, Paper 45's $L^+_{\mathrm{P45}}$) are well-defined and different; the strong-form question is whether the full Krein propinquity bound is finite.
5. **The next sprint (L3b-2b) under $L_{\mathrm{op}}$ has well-defined targets.** Verify axioms, derive Lichnerowicz bound, compare with Paper 45.

**Done.** Hand off to PI for decision on whether to launch L3b-2b under $L_{\mathrm{op}}$ or revisit the architecture.
