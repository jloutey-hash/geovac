# Track TS-C: Three-Asymmetry Investigation

**Sprint:** TS-C (research memo, not a paper)
**Date:** 2026-05-04
**Goal:** Investigate the three highest-priority asymmetries surfaced by Track TS-B (`debug/track_ts_b_dictionary_memo.md`) between the Paper-34 projection taxonomy and the Paper-32 spectral-triple framing.

**Source files read:**
- `debug/tx_a_dimensional_axiomatization_memo.md` (TX-A claims)
- `debug/data/tx_a_composition_table.json` (TX-A data)
- `debug/tx_a_signature_arithmetic.py` (TX-A driver, re-run live)
- `debug/track_ts_b_dictionary_memo.md` (Track B asymmetry list)
- `papers/group5_qed_gauge/paper_2_alpha.tex` (B/F/Δ canonical decomposition)
- `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` (15-projection list, depth-prediction §VI, open question §VIII.3)
- `papers/group6_precision_observations/paper_35_time_as_projection.tex` (rest-mass + observation/temporal-window proposition)
- `papers/group3_foundations/paper_31_universal_coulomb_partition.tex` (universal/Coulomb A/D split)
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` (K-decomposition Sec. \ref{sec:K_decomposition})
- `debug/data/track_alpha_phase4b/track_d_analysis.md` (Hopf graph morphism, α-D)
- `debug/data/track_alpha_phase4f/track_j_analysis.md` (F = D_{n²}(d_max), α-J)

---

## 0. Executive Summary

| Sub-task | Verdict |
|:--|:--|
| **C-1 (TX-A memo vs data on `spectral_action ↔ observation_window`)** | **MEMO WRONG** (data is right). The composition is `commute` in both directions, full stop. The headline "UV/IR ordering anomaly captured by signature arithmetic" claim of the TX-A memo §4 is unsupported by the JSON and is also unsupported by the script logic that produced the JSON. The narrative in TX-A §4.4 should be retracted; the JSON is the ground truth. |
| **C-7 (shared implicit Hopf in K = π(B+F−Δ))** | **NOT CONFIRMED**. The shared-Hopf hypothesis fails three independent tests: (i) F = D_{n²}(d_max) is derived without invoking the Hopf base, (ii) Δ⁻¹ = g_3^Dirac is derived without invoking the Hopf base, (iii) the Hopf-quotient graph at n_max=3 (Phase 4B α-D) does not encode K through any independent spectral invariant. The 6-order-of-magnitude under-shoot of the depth-rule is a **real anomaly**; the resolution is not "shared Hopf," but most likely Paper-34 §VIII.3's existing reading: K is a coincidence below the framework's intrinsic resolution, with Hopf entering each chain explicitly through the **outer** π prefactor (not as a hidden additional projection). |
| **#14 (rest-mass as central deformation)** | **CONFIRMED** as a central deformation, with a sharper formulation than Track B's: rest-mass is an additive shift by a central element of the **enveloping** algebra (m²·𝟙 in the squared-Dirac sector for KG; a more substantive but still central γ⁰m for true Dirac). Of the candidate other "central" projections (hopf, sturmian, scale, drake_swainson), **only Drake-Swainson is a credible second resident** of a "central" tier; the others fail because they introduce non-central content. The proposal is that **Paper 31 acquires a third tier — the *parametric/central* sector — sitting between universal (A-only) and Coulomb-specific (D-substantive)**, with rest-mass and Drake-Swainson as charter members. |

**Most surprising finding (single paragraph, ~150 words):**

When Track B noticed that the TX-A memo's headline non-commutation (spectral_action ↔ observation_window) was missing from the JSON, the natural assumption was that the JSON had a stale layer-grading bug — fix the grading and the memo would be vindicated. The opposite turned out to be true: the script logic is *correct under any reasonable refinement*, and TX-A's narrative claim was a memo confabulation that retro-justifies the desired field-theoretic UV/IR ordering analogy without checking what the metadata produces. The JSON is honest: in TX-A's coarse-grained signature arithmetic, *every* output_layer-2 projection commutes with every other output_layer-2 projection, because none of them feed the other's input. Drake-Swainson is the unique non-commuting partner because it alone has input_layer = output_layer = 2. The framework's "UV/IR anomaly" is real physics that signature arithmetic *cannot see*. This is the cleanest illustration so far that the TX-A axiomatization captures less structure than the memo gave it credit for, and that Track B's call for a finer layer grading is the right next move.

---

## 1. Sub-task 1 — Verify C-1 (TX-A memo/data discrepancy)

### 1.1 The memo claim

`debug/tx_a_dimensional_axiomatization_memo.md` contains the following statement at §0 (Executive summary, last column of the T3 row):

> Most surprising entry: `spectral_action ; observation_window` and its reverse are both one-way, breaking commutation across the temperature/cutoff scale duality.

and at §4 ("Most surprising entry"):

> The `spectral_action ; observation_window` (and its reverse) are both flagged as one-way / asymmetric depending on which layer they target. ... Both target Layer 2, so neither can compose with the other in the "feed forward" sense (each consumes the bare graph independently). ... Signature arithmetic correctly flags this as a non-commutation, even though it cannot tell you *which* order is physically required.
>
> This is the cleanest illustration in the table that the dimensional axiomatization captures real structure: the failure to commute encodes the same content as the field-theoretic UV/IR ordering.

This claim is also reproduced verbatim in `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` §IX (Conclusion):

> the 15×15 composition table encodes the field-theoretic UV/IR ordering anomaly through the mutual one-way relation between spectral action and observation/temporal-window.

### 1.2 What the data say

`debug/data/tx_a_composition_table.json` returns the following:

```json
"spectral_action":      { ..., "observation_window": "commute" }
"observation_window":   { ..., "spectral_action":    "commute" }
```

Both directions are `commute`. There is no `one_way` flag for this pair. The aggregate counts (`summary_counts.one_way_AB = 11`, `one_way_BA = 11`) are entirely accounted for by drake_swainson's row (forward) and column (reverse): `drake_swainson` is the *unique* projection whose composition with each of the 11 non-idempotent non-Drake projections is asymmetric.

### 1.3 What the script does

Re-running `debug/tx_a_signature_arithmetic.py` (live, May 2026) reproduces the JSON bit-for-bit. The relevant function is `composition_relation(p1, p2)`:

```python
can_a_then_b = a["output_layer"] >= b["input_layer"]
can_b_then_a = b["output_layer"] >= a["input_layer"]
if not can_a_then_b and not can_b_then_a:
    return "not_composable"
# ... variable conflict check ...
if can_a_then_b and can_b_then_a:
    return "commute"
```

For `spectral_action` (input_layer=1, output_layer=2) and `observation_window` (input_layer=1, output_layer=2):
- `can_a_then_b`: 2 ≥ 1 ⇒ True
- `can_b_then_a`: 2 ≥ 1 ⇒ True
- Variable sets: `{Lambda, alpha}` ∩ `{beta_temp}` = ∅ ⇒ no conflict
- ⇒ returns `commute`

The script's logic is consistent. The JSON is what the script produces. The memo's narrative is what the agent *wished* the script produced.

### 1.4 Diagnosis

The TX-A memo claim is a **confabulation**: the agent that produced the memo had a strong physics prior (UV/IR ordering anomaly is famous in QFT, specifically the non-commutation of `Λ → ∞` and `β → ∞`) and retrofitted that intuition onto the metadata table, describing the spectral_action ↔ observation_window pair as "the most surprising entry" without verifying that this characterization matched the JSON it had just written.

The mechanism the agent invokes ("each consumes the bare graph independently, so neither can compose in the feed-forward sense") is a structurally sensible point in the physics, but it is **not what the script encodes**. The script encodes a layer-compatibility check (`output_layer ≥ input_layer`), and under that test, two different Layer-1→Layer-2 projections commute trivially: each can be applied to the bare graph, and then either the cutoff Λ or the inverse-temperature β is applied to the resulting Layer-2 spectrum. The signature arithmetic does not distinguish "Layer 2 consuming bare graph" from "Layer 2 acting on a Layer-2 output" because the layer index is a single integer, not a Layer-2-output-type tag.

The memo's intuition is right (the UV/IR limits do not commute in real QFT), but the way it would manifest at the metadata level requires a finer grading. Track B C-1 already proposed this: split output_layer = 2 into UV-cutoff Layer 2, thermal Layer 2, and flow-regulated Layer 2.

### 1.5 Recommended fixes

#### Recommended CLAUDE.md correction paragraph

(Currently CLAUDE.md §2 contains, verbatim, the TX-A memo paragraph including the spectacular UV/IR claim. The relevant text is the paragraph beginning "TX-A dimensional axiomatization (May 2026)"; recommend appending a single corrective sentence:)

> *Erratum (2026-05-04, Track TS-C):* TX-A's headline claim that `spectral_action` and `observation_window` are mutually one-way is a memo confabulation, not a JSON-supported result. The actual TX-A composition table (`debug/data/tx_a_composition_table.json`) flags the pair as `commute` in both directions, and the script logic correctly returns this under the coarse layer-grading TX-A used. The 22 one-way cells are exclusively drake_swainson's row and column. The QFT UV/IR ordering anomaly is real physics but is not captured by TX-A's metadata; Track B C-1 (`debug/track_ts_b_dictionary_memo.md`) flagged this as the highest-priority refinement of the layer index.

#### Recommended Paper 34 §IX (Conclusion) correction

(Currently lines 1238–1241 of `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` repeat the wrong claim. Recommend replacing:)

> *(replace)* "the 15$\times$15 composition table encodes the field-theoretic UV/IR ordering anomaly through the mutual one-way relation between spectral action and observation/temporal-window."

> *(with)* "the 15$\times$15 composition table flags Drake--Swainson as the unique projection with asymmetric composition relations, accounting for all 22 one-way cells. The field-theoretic UV/IR ordering anomaly between `spectral_action` and `observation_window` is real physics but is *not* encoded by the current coarse layer-grading; refining the output\_layer index to distinguish UV-cutoff, thermal, and flow Layer-2 outputs is the natural next step (Track TS-C 2026-05-04)."

#### Recommended TX-A memo correction

(In `debug/tx_a_dimensional_axiomatization_memo.md`, §4 "Most surprising entry": the entire paragraph claiming spectral_action ↔ observation_window is mutually one-way should be retracted in favor of:)

> *(replace)* The most surprising entry is the asymmetric position of drake_swainson. All 22 one-way cells in the JSON are drake_swainson's row (forward composition with 11 partners) plus drake_swainson's column (reverse composition with the same 11). drake_swainson is the unique projection with input_layer = output_layer = 2 (idempotent), making its layer behavior categorically different from every other projection. The pair `(spectral_action, observation_window)` is `commute` in the JSON, despite the suggestive QFT UV/IR ordering anomaly: both target output_layer=2 from input_layer=1, both have disjoint variable sets, and signature arithmetic at the current layer-grading granularity correctly returns `commute`. The non-commutation of `Λ → ∞` and `β → ∞` is real physics that TX-A cannot see; the resolution is a finer layer index, not a tweak to the existing one. *(See Track TS-C 2026-05-04 for the diagnosis.)*

### 1.6 Verdict

**MEMO WRONG, DATA RIGHT** — the JSON and the script that produced it are internally consistent; the TX-A memo's headline narrative is a confabulation that the PI should treat as retracted. The structural lesson (the framework has a real UV/IR ordering anomaly that signature arithmetic cannot see at this granularity) is preserved in Track B's C-1 seed.

---

## 2. Sub-task 2 — C-7 (Shared implicit Hopf in K = π(B+F−Δ))

### 2.1 The hypothesis

Track B §2.c.7 observed that under Paper-34's two-layer reading, the K = π(B+F−Δ) coincidence decomposes into three projection chains:

- **B = 42:** Fock alone (Phase 4B α-C: B is the truncated Casimir trace at m=3 on the (n,l) lattice).
- **F = π²/6:** Fock + spectral_action (or its Dirichlet sub-projection at s=4) (Phase 4F α-J: F = D_{n²}(d_max) = ζ_R(2)).
- **Δ = 1/40:** Fock + spinor_lift (Phase 4H SM-D: Δ⁻¹ = g_3^Dirac = 40, the third single-chirality Camporesi–Higuchi Dirac mode degeneracy).

Paper 34's Prediction 1 (§VI) states that error compounds with projection depth, with a 2–3% residual at depth 2–3. Observed K residual is **8.8 × 10⁻⁸**, six orders of magnitude tighter. Track B's hypothesis: the three chains share a fourth implicit projection (most likely Hopf), invoked but not counted in the depth.

### 2.2 Test 1 — F derivation

From `debug/data/track_alpha_phase4f/track_j_analysis.md`, the **clean derivation of F**:

> D_{n²}(s) = Σ_n n² · n^{−s} = ζ(s − 2)
>
> At s = 4: D_{n²}(4) = ζ(2) = π²/6 = F EXACTLY

The derivation invokes:
1. The Fock degeneracy g_n = n² (Paper 7 §VI).
2. The packing exponent s = d_max = 4 (Paper 0).
3. The Dirichlet zeta ζ_R(s−2).

**No reference to the Hopf bundle anywhere in the derivation.** The Hopf bundle structure S¹ → S³ → S² nowhere enters the Dirichlet computation: the (n,l) lattice and its degeneracy g_n = n² are spelled out in pure SO(4) representation theory, before the Hopf structure is invoked. The π in F = π²/6 is the *intrinsic* π of the Riemann zeta function ζ_R(2), arising from the Bernoulli number B_2 = 1/6 via Euler's identity ζ_R(2k) = (−1)^{k+1} (2π)^{2k} B_{2k} / (2(2k)!) — categorically different from the Hopf-base measure π = Vol(S²)/4.

**Test 1 verdict: NO HIDDEN HOPF in F.** The π in F is a Riemann-ζ π, not a Hopf-measure π.

### 2.3 Test 2 — Δ derivation

Δ⁻¹ = g_3^Dirac = 2(n+1)(n+2) at n=3 = 2·4·5 = 40 (Phase 4H SM-D). The derivation invokes:
1. The Camporesi–Higuchi Dirac spectrum on S³ (Paper 32 §3.3): |λ_n^D| = n + 3/2 with degeneracy g_n^D = 2(n+1)(n+2).
2. The cutoff n = 3 (the same triple-selected cutoff as B).

**No reference to the Hopf bundle.** The Camporesi–Higuchi spectrum is derived from SO(4) representation theory on the Dirac sector of the spectral triple, independent of any Hopf decomposition. The Hopf-equivariant decomposition (Phase 4I D4) gives a *cleaner labeling* (20 ⊕ 20 charge-parity split) but the multiplicity 40 is unchanged.

**Test 2 verdict: NO HIDDEN HOPF in Δ.** The 40 is a Dirac-mode multiplicity, not a Hopf base/fiber datum.

### 2.4 Test 3 — Hopf graph morphism (α-D)

`debug/data/track_alpha_phase4b/track_d_analysis.md` directly tested whether the Hopf-quotient graph at n_max = 3 carries any spectral invariant that produces K. The answer was *clean negative*:

> The raw graph Laplacian spectral invariants of L_{S³} and L_{S²} do NOT encode K. ... None of these hit K, K/π, B = 42, or B + F − Δ to better than 3% in any natural combination.
>
> The integer B = 42 IS recovered exactly, but as the degeneracy-weighted Casimir sum ... This is the quotient-space Casimir trace, not a graph Laplacian spectral invariant. It reproduces Paper 2's B = 42 construction exactly and therefore gives K = π·(B+F−Δ) to 4.8e-7 relative error. **But this is Paper 2's existing formula restated in graph language, not a new derivation.**

So the Hopf bundle is *not* providing an independent fourth chain that links B, F, Δ — it is at most a relabeling that recovers Paper 2's construction.

**Test 3 verdict: HOPF GRAPH MORPHISM IS TAUTOLOGICAL on K.** It does not give an independent derivation, hidden or otherwise.

### 2.5 What about the outer π prefactor?

Paper 2's combination rule is K = **π** · (B + F − Δ). The outer π *is* identified with the Hopf bundle measure (Sprint A α-PI, 2026-04-18; CLAUDE.md §1.7 WH5):

> α-PI POSITIVE PARTIAL identifies the external π in K = π(B+F−Δ) as the Hopf-bundle measure factor Vol(S²)/4 = Vol(S³)/Vol(S¹) = π, equivalently a₀² = π (squared Seeley-DeWitt zeroth coefficient)

This is the **single instance** of Hopf entering the K formula explicitly: as the outer prefactor. It is *already counted* as one of the projections (Hopf projection #2 in Paper 34's list). The K-chain in Paper 34's matches catalogue (Table \ref{tab:catalog}, "Conjectural / Observation status" section) is currently tagged:

> $K = \pi(B + F - \Delta) \approx 1/\alpha$ : Fock $\circ$ Hopf $\circ$ spinor

— i.e., depth 3, with Hopf as the explicit middle projection providing the outer π. There is no "hidden" extra Hopf projection beyond this.

### 2.6 Resolution of the depth-prediction anomaly

The Track B hypothesis "shared implicit Hopf across all three chains" is **wrong**. Hopf enters the K formula *once*, as the outer π prefactor, and is already counted in the depth-3 tag.

The 6-order-of-magnitude under-shoot of the depth-rule prediction (3% expected vs. 8.8 × 10⁻⁸ observed) is therefore a **real anomaly under Prediction 1**.

This is exactly what Paper 34 §VIII.3 says, and it corresponds verbatim to the three readings already enumerated there:

> Either the projection-depth model under-counts what's chained here, or something structural we have not identified is forcing the combination, or the match is a coincidence below the framework's intrinsic resolution.

The current sprint Track C-7 rules out the first reading (no hidden Hopf, no other obvious common projection — Phases 4B–4I exhausted nine attempts to find one). The K-CC sprint (May 2026, see CLAUDE.md §1.7 WH5(v)) ruled out a unified Connes–Chamseddine spectral-action mechanism. Twelve mechanisms total are now eliminated. The remaining live readings are:

- **(ii) Structural mechanism not yet identified**: still possible but increasingly thin given the eliminated mechanisms.
- **(iii) K is a coincidence below the framework's intrinsic resolution**: WH5's standing reading; Paper 2 in Observations status (2026-05-02 audit).

### 2.7 Implications for the depth-prediction

The depth-prediction is *not refuted* by K because Paper 34 already explicitly identifies K as the outlier (§VIII.3). What Track C-7 adds:

1. The "shared implicit Hopf" mechanism is closed (this sub-task).
2. The K-CC clean negative (CLAUDE.md §1.7 WH5(v)) closed the unified spectral-action mechanism.
3. Twelve mechanisms total now eliminated for "K = single-principle derivation."

Paper 34's depth-prediction therefore continues to apply to all *non-K* multi-projection matches in the catalogue, with K as the *named, explicitly-flagged* exception. This is the honest reading. The depth-rule does not have a known exception class beyond K itself.

### 2.8 Recommended Paper 34 §VIII.3 sharpening

Currently §VIII.3 (lines 1176–1191 of `paper_34_projection_taxonomy.tex`) reads:

> **$K = \pi(B + F - \Delta) \approx 1/\alpha$.** The Paper~2 conjecture is the most extreme multi-projection coincidence in the catalogue ... Under Prediction~\ref{pred:error_depth} a three-projection match should carry $\sim 3\%$ error, not $10^{-7}$. Either the projection-depth model under-counts what's chained here, or something structural we have not identified is forcing the combination, or the match is a coincidence below the framework's intrinsic resolution.

Recommend appending a Track-C-7 update sentence (PI to apply or redirect):

> **Track C-7 update (2026-05-04).** The "depth model under-counts" reading was tested by examining whether the three chains share a fourth implicit projection. They do not: F is derived from a Dirichlet zeta on the (n,l) lattice with no Hopf invocation (Phase 4F α-J); Δ⁻¹ is a Camporesi–Higuchi Dirac mode degeneracy with no Hopf invocation (Phase 4H SM-D); and the Hopf-quotient graph at n_max = 3 does not encode K through any independent spectral invariant beyond a Casimir relabeling (Phase 4B α-D). The outer π prefactor of K *is* identified with the Hopf-base measure (Sprint A α-PI), but this is already counted as the Hopf projection in the depth-3 chain. Under-counting is therefore not the mechanism. Combined with the K-CC clean negative on a unified Connes–Chamseddine spectral-action source (CLAUDE.md~\S 1.7 WH5(v)), twelve mechanisms total have been eliminated. The remaining live readings are (ii) and (iii) above; the standing interpretation is (iii), per WH5.

### 2.9 Verdict

**SHARED-HOPF NOT CONFIRMED.** The three chains for B, F, Δ do not share a hidden Hopf projection beyond the outer π prefactor that is already counted. The 6-order-of-magnitude under-shoot is a *real* and *named* exception to the depth-prediction; the resolution lives in the WH5 reading (K is a coincidence below the framework's intrinsic resolution), not in a fixable depth-counting error. Paper 34 §VIII.3 should be sharpened to record the Track C-7 closure.

---

## 3. Sub-task 3 — Rest-mass central deformation (Paper 34 projection #14)

### 3.1 Track B's framing

Track B §2.c.4 observed:

> Rest-mass #14 ... commutes with all 13 non-drake projections, introduces no transcendental ... and adds a single variable + dimension (m, mass). At the TX-A signature-arithmetic level it is "innocent." But at the spectral-triple level it modifies $D$: $D \mapsto D + m\cdot\mathbb{1}$ ... *Any* modification of $D$ is, by Paper 31's partition, structurally Coulomb-specific ... The asymmetry is: *the dictionary correctly tags rest-mass as universal at the V/D/T level, but the spectral-triple framing classifies it as `D`-touching, and Paper 31's partition would naively classify `D`-touching as Coulomb-specific.*
>
> The resolution is that *additive scalar deformations of $D$* are an exception to Paper 31's partition: they are $D$-modifications that commute with everything in $\mathcal{A}$ and therefore live in the centre $Z(\mathcal{A}) = \mathbb{C}\cdot\mathbb{1}$.

The Track C-3 question: is rest-mass the only projection that lives in this central tier, or are there others?

### 3.2 Sharpening "central deformation" precisely

The natural reading of "central deformation" in spectral-triple language is a deformation that takes the form

$$
D \mapsto D + c \cdot \mathbb{1}_{\mathcal{H}}
$$

where $c \in \mathbb{C}$ is a constant scalar and $\mathbb{1}_{\mathcal{H}}$ is the identity on the Hilbert space. Equivalently, the deformation is an element of the **centre of the algebra acting on $\mathcal{H}$**, which for the abelian $\mathcal{A}_{\mathrm{GV}} = \mathbb{C}^{V_{\mathrm{Fock}}}$ in Paper 32 is itself, and for the relevant operator algebra $\mathcal{B}(\mathcal{H})$ is $\mathbb{C} \cdot \mathbb{1}$.

For the **Klein-Gordon** sector, rest-mass is exactly this: the dispersion relation shifts as $\omega_n^2 \to \omega_n^2 + m^2$ (Paper 35 Eq. 3), which at the operator level is $D^2 \to D^2 + m^2 \mathbb{1}$. The centre of the squared-Dirac sector contains $m^2 \cdot \mathbb{1}$, and rest-mass is a deformation by this central element.

For the **true Dirac** sector, rest-mass is more subtle: the standard relativistic deformation is $D \to D + m \gamma^0$, where $\gamma^0$ is a chirality matrix that does *not* live in the centre of the spinor algebra. However, the *square* of this deformed operator is $D^2 + m^2 \mathbb{1} - m \{\gamma^0, D\}$, where the anticommutator $\{\gamma^0, D\}$ vanishes in the chiral basis when $D$ is the Camporesi–Higuchi spinor Dirac operator on S³. So rest-mass shifts $D^2$ centrally even in the Dirac case.

**The precise framing:** rest-mass is a **central deformation of the squared Dirac sector**, equivalently an additive shift of the Dirac eigenvalues by $m^2$ via the squared spectrum. Track B's "framing (i)" is correct at the squared-D level; "framing (ii)" (tensor product) is *not* what happens — rest-mass does not introduce a new Hilbert-space factor.

### 3.3 Why this exceptions Paper 31's partition

Paper 31's universal/Coulomb partition is the **A versus D split**: results that depend only on $\mathcal{A}$ are universal; results that depend on the specific form of $D$ are Coulomb-specific.

Rest-mass deforms $D$ but does so *centrally*: it adds $c \cdot \mathbb{1}$ to (the square of) $D$, which lives in $Z(\mathcal{B}(\mathcal{H})) = \mathbb{C} \cdot \mathbb{1}$. A central deformation has the property that

$$
[D + c\mathbb{1}, a] = [D, a] \quad \forall a \in \mathcal{A},
$$

i.e., the bounded-commutator structure (Connes axiom 2) is unchanged. The order-one axiom (axiom 4) is also unchanged because the central element commutes with everything. The spinor-bundle structure, real-structure J, KO-dimension — all are preserved. **Rest-mass changes nothing about the spectral triple's *axiomatic* structure**; it only re-scales the spectrum.

This is structurally why rest-mass is universal: any central potential admits an additive mass shift, and the shift does not engage any of the features that distinguish the Coulomb triple from the HO triple. Paper 31's universal sector is "things that depend on $\mathcal{A}$ only"; the natural extension is "things that depend on $\mathcal{A}$ and on the centre of the operator algebra." This is not yet articulated in Paper 31.

### 3.4 Candidate other central-deformation projections

Track B suggested testing hopf, sturmian, scale, drake_swainson. Adding the natural physics candidate global-phase / Wick rotation. Test each:

| Candidate | Is it $D \to D + c\mathbb{1}$? | Verdict |
|:--|:--|:--|
| **Hopf bundle (#2)** | NO — Hopf decomposes $\mathcal{H}$ into base+fibre and (in the almost-commutative extension to U(1) gauge) introduces *inner* derivations of $\mathcal{A} \otimes \mathbb{C}$. The U(1) action is by a non-central element (a $\mathrm{diag}(e^{i\theta_v})$ matrix). | NOT CENTRAL. Hopf is a substantive bundle decomposition; rejecting. |
| **Sturmian (#5)** | NO — Sturmian is a change of $\mathcal{H}$-basis at exponent $\lambda = Z/n$. It is a unitary on $\mathcal{H}$, and unitaries that are not multiples of $\mathbb{1}$ are not central. The Sturmian basis is rich (the basis indices re-organize via $\lambda$), not central. | NOT CENTRAL. Sturmian is a basis re-parameterization; rejecting. |
| **"Scale" (CLAUDE.md §1.4 focal length $p_0^2 = -2E$)** | This is *part of the Fock projection* (which conformally maps the energy-shell onto S³); it is not a separate projection in the dictionary. The energy-shell rescaling itself is a $\sqrt{-2E}$ scaling of the radial coordinate. This is a multiplicative rescaling, *not* a central shift; the Dirac operator transforms covariantly $D \to \lambda D$, not $D \to D + c\mathbb{1}$. | NOT CENTRAL (and not a projection on its own). |
| **Drake–Swainson (#13)** | The K subtraction-scale is *transient* and the final result is K-independent: this is structurally an additive shift of an intermediate spectral sum that *cancels* in the converged answer. The ledger shift is along $\mathbb{R} \cdot \mathbb{1}$ in the spectral-density domain — i.e., it modifies the integration measure by an additive constant that is then subtracted. **This is an additive, central, and ultimately self-cancelling shift in the flow tier — exactly the kind of structure rest-mass exhibits, with the *additional* property of cancellation.** | **CONFIRMED CENTRAL** (with cancellation). Drake-Swainson is a candidate second resident of the central tier. |
| **Wick rotation / global phase** (Track B speculative) | Wick rotation $t \to -i\tau$ is a coordinate change on the temporal direction; it is not a deformation of $D$. Global phase $\psi \to e^{i\phi}\psi$ is a $U(1)$ gauge in the abelian sector; it commutes with everything but acts by a unitary, not by $D \to D + c\mathbb{1}$. | NOT CENTRAL in the Track B sense; Wick rotation is structural (KO-dim transition; see Track B C-6) and global phase is a unitary not a $D$-shift. |

Two robust candidates emerge: **rest-mass** and **Drake–Swainson**. Both have the structural form "additive shift of $D$ (or its square, or its associated spectral measure) by a central element of the operator algebra." Both commute with every $a \in \mathcal{A}$. Both preserve all Connes axioms. Both are universal across central potentials (rest-mass: any system has a mass; Drake-Swainson: any flow-regularized observable can use this subtraction). Neither moves the operator system or the bundle structure.

### 3.5 Proposal: a third tier in Paper 31

Paper 31's current partition is binary: universal (depends only on $\mathcal{A}$) versus Coulomb-specific (depends on $D$). The Track B observation is that this is too coarse: there is a class of $D$-deformations that are not Coulomb-specific because they live in the centre.

**Proposal (PI to evaluate):** add a third tier — the **parametric** or **central** sector — sitting between universal and Coulomb-specific:

> **Universal sector** ($\mathcal{A}$-only): angular sparsity, composed scaling laws, selection rules. Independent of $D$.
>
> **Parametric / central sector** ($\mathcal{A}$ + central deformations of $D$): rest-mass shifts, flow-tier subtractions. Depend on the centre $Z$ of the operator algebra acting on $\mathcal{H}$, which is just $\mathbb{C} \cdot \mathbb{1}$ for an irreducible spectral triple. Universal across central potentials in the same way as the universal sector but additionally support the addition of scalar parameters.
>
> **Potential-specific sector** (substantive $D$): rigidity theorems, calibration $\pi$, transcendental content, the $\alpha$-conjecture. Depend on the specific non-central structure of $D$.

This proposal is *low cost* (descriptive, not predictive — it sorts existing results into three bins instead of two) and *high clarity* (it explains why rest-mass and Drake-Swainson behave like universal-sector projections at the V/D/T level despite being $D$-touching). It does not affect any of Paper 31's existing theorems.

### 3.6 Recommended Paper 31 amendment (PI to apply or redirect)

Two-sentence proposal for a new subsection under Paper 31 §IV (or §VI; PI to choose location):

> **Parametric / central sector.** A third tier of the framework's $\mathcal{A}/D$ partition consists of deformations of $D$ by central elements of the operator algebra acting on $\mathcal{H}$, equivalently additive scalar shifts of the spectrum that preserve all Connes axioms (bounded commutators, order-one, real structure, KO-dimension). The two charter members are the rest-mass projection (Paper 35; an additive shift $D^2 \to D^2 + m^2 \mathbb{1}$ that is universal across central potentials) and the Drake–Swainson asymptotic subtraction (Paper 34 §III.13; a transient, cancelling shift in the spectral-density domain). Both are formally $D$-touching but functionally universal; they support the introduction of physical scalar parameters into the framework without engaging any potential-specific structure of $D$.

### 3.7 Verdict

**CONFIRMED with sharpening.** Rest-mass is genuinely a central deformation of the squared-Dirac sector (Track B's framing (i) holds; framing (ii) does not). Of the candidate other "central" projections in the dictionary, only **Drake–Swainson** survives as a credible second resident; hopf, sturmian, scale, and Wick rotation all fail the central-shift test for distinct structural reasons. The proposal is to add a third tier to Paper 31's universal/Coulomb partition: **the parametric/central sector**, whose charter members are rest-mass and Drake-Swainson.

---

## 4. Cross-cutting observations

### 4.1 Two of three asymmetries point at the same structural blind spot

- **C-1** (TX-A memo vs data): the TX-A signature arithmetic at coarse layer-grading does not distinguish UV-cutoff Layer 2, thermal Layer 2, and flow Layer 2 outputs.
- **C-7** (K depth-rule violation): the depth-rule under-counts K's structural content if there is a hidden invariant the framework has not identified — but the present sprint shows there is *not* a hidden Hopf, leaving (ii) and (iii) as the live readings.

The two asymmetries differ in their relationship to the framework: C-1 says TX-A's grading is too coarse and proposes a fix; C-7 says K is anomalous *despite* the framework working correctly elsewhere, and the WH5 reading (K is a coincidence below intrinsic resolution) is the increasingly likely interpretation.

The shared observation between C-1 and C-7: **the TX-A axiomatization sees less than Track B / Paper 32 do**. C-1 is about the depth of the *layer index*; C-7 is about the depth of the *projection chain*. Both call for refinement, but in different directions.

### 4.2 The third tier proposal (Paper 31) is independently motivated by C-1's UV/IR observation

If the Paper 31 amendment in §3 is adopted (parametric/central sector), it has the side effect of giving a *language* for understanding why two outwardly-similar projections (`spectral_action` and `observation_window`) can have different properties: spectral_action is parametric (introduces $\Lambda$, an external scale; the spectral action *is* a central deformation in the test-function sense, $\mathrm{Tr}\, f(D^2/\Lambda^2)$), while observation_window is genuinely structural (compactifies the temporal direction, changes the topology of the ambient manifold, eventually introduces KO-dim transitions per Track B P-6). The TX-A signature arithmetic conflates these because both add (Λ, β) variables and (energy, time) dimensions; the third-tier framing distinguishes them at the level of *what kind of deformation is being applied to D*. This is a candidate refinement of C-1's coarse-grading objection that does not require enriching the layer index.

---

## 5. Report-back summary

- **Path to memo:** `debug/track_ts_c_asymmetry_investigation_memo.md`
- **Three verdicts as a triple:** (C-1: MEMO WRONG, DATA RIGHT; C-7: SHARED-HOPF NOT CONFIRMED; #14: CENTRAL DEFORMATION CONFIRMED with sharpening).
- **Recommended edits drafted:** (i) CLAUDE.md erratum on TX-A's headline UV/IR claim; (ii) Paper 34 §IX (Conclusion) correction replacing the same wrong claim; (iii) TX-A memo §4 retraction; (iv) Paper 34 §VIII.3 sharpening for K under Track-C-7 closure; (v) Paper 31 third-tier proposal (parametric/central sector with rest-mass + Drake-Swainson). PI to decide whether to apply each.
- **Most surprising finding (~150 words):** When Track B noticed that the TX-A memo's headline non-commutation (spectral_action ↔ observation_window) was missing from the JSON, the natural assumption was that the JSON had a stale layer-grading bug — fix the grading and the memo would be vindicated. The opposite turned out to be true: the script logic is *correct under any reasonable refinement*, and TX-A's narrative claim was a memo confabulation that retro-justifies the desired field-theoretic UV/IR ordering analogy without checking what the metadata produces. The JSON is honest: in TX-A's coarse-grained signature arithmetic, *every* output_layer-2 projection commutes with every other output_layer-2 projection, because none of them feed the other's input. Drake-Swainson is the unique non-commuting partner because it alone has input_layer = output_layer = 2. The framework's "UV/IR anomaly" is real physics that signature arithmetic *cannot see*. This is the cleanest illustration so far that the TX-A axiomatization captures less structure than the memo gave it credit for, and that Track B's call for a finer layer grading is the right next move.

---

*Memo word count: ~3300 words. Sprint TS-C complete. No production code or papers modified; recommended edits drafted as proposed text only, per sprint constraints.*
