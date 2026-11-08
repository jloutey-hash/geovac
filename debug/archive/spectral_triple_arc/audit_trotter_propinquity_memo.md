# Audit memo: sprint_trotter_propinquity

**Date:** 2026-06-04
**Auditor role:** adversarial; burden of proof on the sprint
**Sprint headline:** "STOP — propinquity-aware Trotter bound is uniformly LOOSER than naive Suzuki-Trotter at production parameters; sharpening is multi-step research, not sprint-scale."
**Sprint proposed Paper 20 footnote:** ~25-line LaTeX footnote framing the negative finding with a 4/π M1 master-Mellin-engine positive structural reading.

## Summary verdict

**CONFIRMED with one substantive scope tightening + two technical caveats** — the headline NEGATIVE is correct and robust at every L_H convention tested, but two of the headline numbers (the "19× worse than naive" comparison at n_max=5000, and the "$\gamma_{n_max} \sim (4/\pi)\log n/n$ leading constant" attribution as M1 master-Mellin signature) overclaim and should be tightened before the Paper 20 footnote lands.

## 1. Free-parameter count

The headline number "trunc_bound = 0.40 for LiH, 0.99 for BeH₂, 7.66 for H₂O at n_max=2" depends on:

| Parameter | Source | Free? |
|:----------|:-------|:------|
| γ_{n_max=2} = 2.0746 | `compute_propinquity_bound(2)`, Paper 38 L2 | **Not free** — independently computed at fixed precision |
| N_active ∈ {2, 4, 8} | Active-electron counts from spec | **Not free** — physics choice |
| L_H = mean(\|c_j\|) | **Sprint choice** of mean | **FREE — 4 candidates considered** |
| sim_time t = 1.0 a.u. | Sprint default | **FREE — 1 value tested** |
| ε = 10⁻³ (production target) | Sprint default | **FREE — 3 values tested (1e-1, 1e-3, 1e-6)** |
| Even split (ε_trunc = ε_trotter = ε/2) | Sprint choice | **FREE — alternative splits not tested** |

The most consequential free parameter is **L_H**. The sprint cycles through max/p90/mean/median in the driver, picks `mean`, and motivates it inline (lines 244-249 of driver) as "falsifier-safe representative." The other three are reported but not used.

This is a sub-percent literature-convention free parameter being absorbed by inline justification. **My re-run at all four conventions:** trunc_bound for H₂O ranges from 0.73 (median) to 3043 (max) — three orders of magnitude. The overshoot vs ε/2 = 5×10⁻⁴ budget ranges from 1459× (median) to 6×10⁶× (max). The verdict (NEGATIVE) **survives at every convention**, but the specific magnitude depends on L_H choice. The memo should report the L_H sensitivity explicitly rather than picking one.

## 2. Selection bias

Three alternatives tested at the bound level (max/p90/mean/median for L_H); sprint picks the middle value. This is honest. No cherry-picking on the headline verdict — every L_H convention gives the same NEGATIVE.

However, the **comparison frame** (propinquity-aware vs naive Suzuki-Trotter) has a structural bias: the sprint compares `r_prop` at the **feasibility-crossover n_max = 5000** to `r_naive` computed at the **production n_max = 2**. From the sprint's own driver comments (lines 86-87), λ and α_comm both grow with n_max. So:

- `r_prop(n_max=5000) = 13,934` uses production-cutoff α_comm (which is what the driver actually has access to)
- `r_naive(n_max=2) = 729` uses production-cutoff λ
- **Honest comparison:** `r_naive(n_max=5000)` would be much larger than 729 because λ scales with the basis. The "19× worse than naive" claim therefore **rigs the comparison in favor of naive** (gives naive the small-basis λ while the propinquity bound is being evaluated at the feasibility crossover).

This is a known limitation in the sprint's own `n_max_feasibility_sweep` function (driver lines 372-378): "We DO NOT rebuild the Hamiltonian at each n_max; instead we use the α_comm / L_H from the production n_max=2 point as the proxy." But the memo does not flag the same issue for r_naive. **Fix:** either (a) say "both bounds evaluated with PRODUCTION α/λ" and label the n_max=5000 row as illustrative, or (b) state that r_naive at n_max=5000 would be ~k× worse than 729 (some published growth rate), making the propinquity comparison even less favorable to propinquity. Either way the verdict (NEGATIVE) is unchanged and may strengthen.

## 3. Independent test

The sprint's bound is derived analytically; there is no held-out check beyond "FCI error at LiH n_max=2 is ~2% (Paper 20)" as a sanity floor for `||H_∞ - H_{n_max}||_op`. This is named in the memo §6 as "three orders of magnitude smaller than the propinquity-derived 0.40 upper bound" — and that *is* an independent test (FCI error is computed separately from the propinquity bound). It confirms looseness, which is consistent with the sprint's NEGATIVE.

**A genuinely independent test would be:** compute `||H_{n=3} - H_{n=2}||_op` directly (e.g., as an SVD or Frobenius norm difference on a small molecule) and compare to `(γ_{n=2} - γ_{n=3}) · N_active · L_H`. The sprint does not do this. This would distinguish "the bound is structurally loose by 1000×" from "the bound is structurally loose for a reason that goes away at moderate n_max." Recommended for any follow-on.

## 4. Robustness

I re-ran the bound calculation at the four L_H conventions (above). The NEGATIVE verdict is bit-stable at all four. The sweep across ε ∈ {1e-1, 1e-3, 1e-6} is in the JSON: at every ε and every L_H convention, every molecule fails feasibility at n_max=2.

**Not tested:** different sim_times t, different ε-split ratios (e.g., 1:99 instead of 50:50), different molecules outside the LiH/BeH₂/H₂O panel. At larger t (chemistry-relevant t ~ 10-100 to resolve energy gaps ~ 0.01-0.1 Ha), the bound only gets worse because it's linear in t. At larger t-asymmetric ε-splits favoring ε_trotter, the feasibility crossover only moves to larger n_max (because trunc_budget shrinks). Robustness is in the right direction for the NEGATIVE.

## 5. Honest scope

The sprint memo §7 is mostly honest about scope. Three sub-claims need tightening:

**Claim 1 (sprint §2/§9 footnote): "The 4/π constant identifies the M1 (Hopf-base measure) signature of the master Mellin engine of Paper 18."**

This is an over-attribution flagged by the **transcendental-tag memory rule** (`feedback_tag_transcendentals.md`). The L2 asymptotic γ_n ~ (4/π) log n/n is bit-correct (memory `l2_quantitative_rate_4_over_pi.md` is the load-bearing source). But invoking it inside the Trotter bound and claiming the Trotter cost "carries the M1 signature" overreaches. The Trotter cost depends on α_comm and λ (whose transcendental content is **not** M1 — they're Pauli-coefficient sums set by the JW transform of the composed Hamiltonian, with M2/M3 origin in the spectral content per Paper 28/Paper 34). The 4/π appears only in the **truncation bound**, which is **loose by 3 OoM**, so the M1 signature does not actually enter the binding constraint. The footnote should either drop the M1 framing or restrict it to "in the *truncation* bound, which is currently not load-bearing."

**Claim 2 (sprint §5.2 / verdict): "The propinquity-aware bound at feasibility is ~19× MORE expensive than naive Suzuki-Trotter."**

Per §2 above, the comparison frame is unfair to propinquity. The honest version: "evaluated at production-cutoff α_comm with feasibility-crossover γ, the propinquity-aware step count is 19× the naive count at production-cutoff λ. A like-for-like comparison at n_max=5000 would require evaluating λ at n_max=5000, which is computationally inaccessible and is the entire reason the propinquity-aware approach was tried in the first place."

**Claim 3 (verdict / §8): "the structural finding is genuinely new (single-particle propinquity → many-body Trotter via Duhamel + N_active linearity carries the 4/π M1 Hopf-base measure signature)."**

The Duhamel-plus-N_active step **is** standard (Childs-Su-Tran-Wiebe-Zhu 2021 §IV is the cited source). The novelty is plugging the **propinquity-derived γ** into that standard lift, not the lift itself. "Genuinely new" overclaims. The honest version: "applying the Paper 38 propinquity to the standard Duhamel+N_active lift gives a feasibility certificate; the certificate is too loose to be useful but is the first time the propinquity has been instantiated as a quantum-simulation cost bound." This is still publishable as a footnote; it just shouldn't be sold as a structural innovation.

## Track-specific check: propinquity-to-Trotter derivation

The chain is:

1. **L4(c) of Paper 38:** ‖M_f − P_{n_max} M_f P_{n_max}‖_op ≤ γ_{n_max} · ‖∇f‖_∞. **Rigorous, in Paper 38.**
2. **Sprint's lift:** "By the standard N_a-body norm bound (Childs-Su 2021 §IV), per-electron truncation propagates linearly to ‖H_∞ − H_{n_max}‖_op ≤ N_active · γ_{n_max} · L_H." **Citation is in the right published literature, but the bound assumes the composed Hamiltonian decomposes into commuting single-particle blocks — sprint §6.3 admits this is loose by ~2× due to cross-electron ERIs.** The L_H = max-Pauli-coefficient identification is also a separate inequality (Lipschitz-of-multiplier ≤ max-Pauli-coefficient) that is sloppy — the relationship is via Riesz representation on $S^3$, not a simple max. Sprint uses mean instead of max for "falsifier-safe" reasons, which makes the bound numerically smaller but no longer an upper bound in any rigorous sense. **The bound is NOT rigorous as stated.** Honest reframing: "a heuristic estimate, expected to be the right order of magnitude based on Paper 22 Gaunt sparsity."
3. **Duhamel:** ‖e^{-iAt} − e^{-iBt}‖ ≤ t‖A−B‖. **Rigorous.**
4. **CSTWZ Eq. 48 for ε_trotter:** **Rigorous, published.**

So the chain is **rigorous at steps 1, 3, 4** and **heuristic at step 2** (both the N_active linearity and the L_H = mean(\|c_j\|) identification). The sprint's memo §7.4 names the loose-L_H gap but the structural status (heuristic vs rigorous) is not flagged in the verdict.

**Recommendation:** the Paper 20 footnote should say "a heuristic propinquity-derived feasibility estimate" rather than "a propinquity-derived bound." Otherwise readers will assume it's a rigorous upper bound that has been numerically evaluated.

## Pass/fail items

- **PASS** verdict robustness across L_H ∈ {max, p90, mean, median}.
- **PASS** verdict robustness across ε ∈ {1e-1, 1e-3, 1e-6}.
- **PASS** named obstructions (3 sources of slack identified in §6).
- **PASS** honest scope statement in §7 (5 numbered caveats).
- **PASS** STOP decision based on quantitative budget overshoot, not interpretive judgment.
- **FAIL** rigorous-vs-heuristic status of the bound: stated as a "bound" throughout, actually heuristic at step 2 (N_active linearity + L_H = mean, neither rigorous on their own).
- **FAIL** r_naive vs r_prop comparison frame at n_max=5000 (rigged in favor of naive by using production-cutoff λ for naive but feasibility-crossover γ for propinquity).
- **FAIL** M1 signature attribution: the 4/π only enters the (loose, not-load-bearing) truncation bound, not the binding Trotter cost; the "structural finding carries M1 signature" claim is over-attribution.

## Recommended modifications to the sprint's verdict

1. Replace "propinquity-derived **bound**" with "propinquity-derived heuristic estimate" throughout, OR add an explicit "STATUS: heuristic, not rigorous (loose-L_H and N_active-linearity assumptions)" line.
2. In the §5.2 r_prop=13934 vs r_naive=729 comparison, add a sentence: "The 19× factor is evaluated with production-cutoff α_comm and feasibility-cutoff γ; an honest like-for-like comparison would require α_comm at n_max=5000, which is computationally inaccessible — this is the structural reason the propinquity-aware approach was attempted."
3. In the Paper 20 footnote, drop the "carries the M1 (Hopf-base measure) signature" framing OR rewrite as "the 4/π appears in the (currently-loose) truncation bound and may carry the M1 signature once the bound is tightened." The current phrasing is the **transcendental-tag memory rule** firing on anonymous-but-loose-tag of M1.
4. Add an explicit L_H sensitivity row to the §5.1 table — same result, three values per molecule (median/mean/p90) — so the reader sees the dependence rather than the inline justification.

The STOP verdict itself is **CORRECT and robust**. The four modifications above tighten the framing without changing the conclusion.

---

**Audit verdict line:** **CONFIRMED because** the NEGATIVE verdict (propinquity-aware bound looser than naive Suzuki-Trotter at production parameters) is bit-robust across L_H ∈ {max, p90, mean, median} and ε ∈ {1e-1, 1e-3, 1e-6}, while the named overclaims (rigorous-vs-heuristic, n_max=5000 comparison frame, M1 signature attribution) are framing issues that tighten but do not overturn the structural finding.
