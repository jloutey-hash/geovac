# Track SM-B Analysis: Does Σ_f N_c Q_f² = |λ_3| at n_max = 3?

**Date:** 2026-04-10
**Sprint:** α SM-B
**Verdict:** **FAIL (numerical coincidence; no structural map)**

---

## 1. The Observation

Per-generation charged-fermion electric-charge-squared trace (standard QED normalization):

    Σ_f N_c Q_f² (per generation)
      = (−1)² (charged lepton)
      + 3 · (2/3)² (up-type quark, 3 colors)
      + 3 · (−1/3)² (down-type quark, 3 colors)
      = 1 + 4/3 + 1/3
      = 8/3.

Three generations:

    Σ_f N_c Q_f² (SM) = 3 · (8/3) = 8.

And in Paper 2, the top-shell S³ Laplace eigenvalue at the selected cutoff n_max = 3 is

    |λ_3| = 3² − 1 = 8.

**Both quantities equal exactly 8.** The question is whether this identity has a structural origin or is a numerical coincidence.

---

## 2. Where 8 Appears in Paper 2

The value 8 appears in Paper 2 in exactly two places, both tied to |λ_3| = 8:

1. As the magnitude of the largest S³ Laplace eigenvalue included at n_max = 3.
2. Inside the boundary correction
   Δ = 1/(|λ_3| · N(2)) = 1/(8 · 5) = 1/40.

The value 8 does **not** appear additively in B = 42 (which decomposes as 0 + 6 + 36 from Casimir contributions of shells n = 1, 2, 3) nor in F = π²/6. In K = π(B + F − Δ), the 8 enters only as a denominator inside the tiny correction Δ.

Note that Paper 2's B closed form (Eq. 17) is
    B(m) = Σ_{n=1}^m n²(n²−1)/2 = Σ_{n=1}^m g_n |λ_n| / 2,
so B itself is a weighted sum of the |λ_n|. The specific value |λ_3| = 8 enters B with weight g_3/2 = 9/2 (contributing 36 to the total 42), not as an isolated 8.

---

## 3. Candidate Structural Mechanisms

### (a) Generation ↔ S³ shell identification

**Hypothesis:** Map generation index g = 1, 2, 3 to S³ shell index n = 1, 2, 3. Then Σ over generations becomes Σ over shells.

**Test:** The per-generation trace is 8/3 for *every* generation (QED is generation-blind in charge assignments — each generation has identical (−1, +2/3, −1/3) charges). But the natural per-shell quantities on S³ at shell n are:
- |λ_n| = n² − 1: 0, 3, 8 (varying with n)
- g_n = n²: 1, 4, 9
- (n² − 1)/n: 0, 3/2, 8/3 (varying with n)

The per-generation 8/3 is a **constant** across generations; every S³ per-shell quantity **varies** with n. Therefore there is no shell → generation map that reproduces the SM sum from a Paper 2 invariant. The sum 3 · (8/3) = 8 is recoverable only by adding three identical copies of a universal constant — which carries no S³ structure.

*Negative.*

### (b) Arithmetic identity 8/3 = |λ_3| / dim(S³)

**Hypothesis:** Per-generation 8/3 = |λ_3|/dim(S³), and three generations replicate dim(S³) = 3 to give |λ_3|.

**Test:** The identity (8/3) · 3 = 8 is trivial. For it to be structural, we would need:
1. An independent reason why each generation contributes exactly |λ_3|/dim(S³).
2. An independent reason why there are exactly dim(S³) = 3 generations.

Paper 2's selection principle B/N = 3 fixes dim(S³) as the Laplace-Beltrami base dimension, not as a multiplicity of replicated copies. The generation count in the SM is famously unexplained (the "family problem"); Paper 2 does not address it. Without independent derivations of both the per-gen value and the gen count, the "factorization" is a post-hoc rewrite.

*Negative.*

### (c) Anomaly cancellation / Casimir trace identification

**Hypothesis:** The QED β-function coefficient b_0 ∝ Σ N_c Q_f² is structurally a "trace of a charge Casimir" over matter representations. Paper 2's B = 42 is a "trace of an angular Casimir" over the Fock lattice. Maybe the two are the same trace in a unified gauge theory.

**Test:** The SM electric-charge operator is Q = T_3 + Y, with integer Q forced by the SU(2)_L × U(1)_Y embedding of U(1)_em. Fractional quark charges (±1/3, ±2/3) arise because SU(3)_color multiplets share U(1)_em charge with the integer lepton sector via hypercharge assignment. The Hopf bundle S¹ → S³ → S² realizes **U(1)**, not SU(3) × SU(2) × U(1); its first Chern class is integer-valued, so Dirac quantization produces **integer** electric charges from a Hopf U(1) alone. To get 1/3-fractional charges, the Hopf U(1) would need to be embedded into a larger non-abelian bundle (SU(5), SO(10), Pati–Salam). Paper 2 does not do this and the selection principle B/N = 3 provides no hook for non-abelian embedding.

Additionally, Σ N_c Q² = 8 for three SM gens is not the *only* charge-trace invariant. Related traces yield different integers:
- Σ N_c Y² (Weyl, per gen) = 10/3, three gens = 10 (not 8).
- Σ N_c Y² (Dirac, per gen) = 5/3, three gens = 5 (not 8).
- Σ N_c Q^4 (per gen) = 1 + 3·16/81 + 3·1/81 = 130/81, three gens = 130/27.

Only Σ N_c Q² happens to hit 8. The selection of which trace to compute is being done *after* seeing the target.

*Negative.*

### (d) Charge quantization from Hopf first Chern class

**Hypothesis:** The Hopf bundle's first Chern class c_1 = 1 forces charge quantization in integer units; the SM fractional charges are an artifact of the SU(3) color decomposition.

**Test:** This is true, but it cuts against the identification. Hopf c_1 gives integer charges only (Dirac quantization Qe·g_magnetic = n/2 ∈ (1/2)ℤ). The SM has Q ∈ {0, ±1/3, ±2/3, ±1}. A Hopf-only theory cannot produce the SM charge spectrum without an additional SU(3)_color structure that is not present in Paper 2. Furthermore, even if we accepted fractional charges, the arithmetic 8/3 per generation is specific to the (1, 3, 3) multiplicities — it reflects Dirac ± Gell-Mann charge assignments, not bundle topology.

*Negative.*

### (e) Δ = 1/(|λ_3| · N(2)) substitution

**Hypothesis:** Since Paper 2 explicitly uses 8 = |λ_3| in the denominator of Δ, maybe the SM equivalent would be Δ_SM = 1/(Σ N_c Q² · N(2)) = 1/(8 · 5) = 1/40 — same numerical value.

**Test:** This is tautological because the two 8s are equal. It does not provide structure — it renames |λ_3| as "charge-squared trace" without any independent identification. If tomorrow we found a different integer I and asked whether Δ = 1/(I · N(2)) worked, the answer would depend entirely on I, not on any SM property.

*Negative (rename, not derivation).*

---

## 4. Negative Test: Other Charge-Trace Integers

To sanity-check that "8 is special", I computed several charge-like traces over SM fermions:

| Trace | Per gen | 3 gens |
|:------|:--------|:-------|
| Σ N_c Q² | 8/3 | **8** |
| Σ N_c |Q| | 1 + 2 + 1 = 4 | 12 |
| Σ N_c Q⁴ | 130/81 | 130/27 |
| Σ N_c Y² (Weyl) | 10/3 | 10 |
| Σ N_c Y² (Dirac) | 5/3 | 5 |
| Σ N_c Q² including ν_R | 8/3 | 8 (ν_R has Q=0) |
| Σ all Weyl d.o.f. | 15 | 45 |

Only Σ N_c Q² (with or without ν_R) hits 8. This is consistent with "8 = 8" being a selection bias: we picked the one sum that matches. No other natural charge trace hits any of Paper 2's special integers (42, 14, 5, 40).

---

## 5. Deeper Obstruction: n_max = 3 Is Not the "Generation Number"

Paper 2's selection principle fixes n_max = 3 via B(m)/N(m) = 3 ⇔ (m + 2)(m − 1) = 10 ⇔ m = 3 (unique). The "3" here is **dim(S³)**, the manifold dimension of the S³ base. It is not a multiplicity count. The claim "n_max = 3 selects 3 generations" confuses two distinct appearances of the number 3:

1. dim(S³) = 3 (geometric dimension, intrinsic to S³)
2. N_gen = 3 (family index count, currently unexplained in the SM)

There is no principle in Paper 2 that connects these. And empirically, the SM works equally well as a QFT with N_gen = 1 or 2; the number 3 is an experimental input, not a consistency condition. If Paper 2 were to derive 3 generations, it would need an anomaly-cancellation or geometric-obstruction argument, which the selection principle is not.

---

## 6. What About the Per-Shell Weighting?

Paper 2's B = 42 can be written as B = Σ_{n=1}^3 g_n · |λ_n|/2 with per-shell contributions (0, 6, 36). If we imagine a "matter content" assigned per shell with weight g_n/2 = (1/2, 2, 9/2), and we want the weighted sum Σ g_n · W_n / 2 to equal Σ_f N_c Q_f² = 8, we get the equation

    (1/2)·W_1 + 2·W_2 + (9/2)·W_3 = 8.

Setting W_n = |λ_n| (the Paper 2 choice) gives B = 42, not 8. Setting W_n = 8/3 (constant, SM per-gen) gives (1/2 + 2 + 9/2)·(8/3) = 7·8/3 = 56/3 ≈ 18.67, not 8. No simple per-shell matter assignment reproduces both B = 42 and a side identity to 8.

---

## 7. Conclusion

**The equality Σ_f N_c Q_f² = 8 = |λ_3| is a numerical coincidence.** Specifically:

1. The SM 8 lives in a flavor/charge Hilbert space with SU(3)_color × SU(2)_L × U(1)_Y structure. The Paper 2 8 lives in a geometric Hilbert space (L²(S³)) with SO(4) structure. There is no natural map between the two.
2. Paper 2's Hopf bundle has integer first Chern class and cannot reproduce fractional 1/3 quark charges without an additional non-abelian embedding not present in the framework.
3. Per-generation contribution 8/3 is universal across generations; every per-shell S³ invariant varies with n. The "sum over 3 generations = sum over 3 shells" identification is inconsistent.
4. The SM has multiple natural charge traces; only Σ N_c Q² hits 8. Selection bias is the simpler explanation.
5. Replacing |λ_3| with Σ N_c Q² in Δ = 1/(|λ_3|·N(2)) is a rename, not a derivation; it does not change the numerical prediction.
6. Paper 2's n_max = 3 = dim(S³) is a manifold dimension, not a multiplicity count; conflating it with N_gen = 3 is not supported by any argument in Paper 2.

### Should Paper 2 mention this coincidence?

**No.** The reasoning:

- Paper 2 is already conjectural (clearly labeled). Adding a further speculative identification would dilute the focus of the paper, which is to present the spectral-geometric construction and document its current conjectural status.
- The coincidence Σ N_c Q² = 8 = |λ_3| is arithmetically exact but would take a non-trivial detour (hypercharge, anomaly cancellation, color embedding) to even state properly. A reader who is skeptical of Paper 2's main conjecture will find this aside even less rigorous, and a reader who is convinced will not need the additional hook.
- If a future derivation of Paper 2's combination rule produces a structural reason for 8/3 per generation (e.g., via a Kaluza–Klein reduction to SU(3) × SU(2) on S³ × internal manifold), then that future work can cite this coincidence. Until then, logging it in project memory is sufficient.
- Including numerical coincidences without structural backing has historically hurt reception of α-derivation papers (Eddington's 137, Wyler's Γ(3/2) formula); the guideline in the project's framing is to be parsimonious about such claims.

Recommendation: **log in `debug/data/track_alpha_sm/` and in CLAUDE.md Section 2 as an observed but shelved coincidence; do not mention in Paper 2.**

---

## 8. Status Summary Table

| Mechanism | Status | Reason |
|:---|:---|:---|
| (a) Gen ↔ shell map | NEGATIVE | Per-gen 8/3 constant, per-shell varies |
| (b) 8/3 = \|λ_3\|/dim(S³) | NEGATIVE | Trivial factorization, no independent derivation of gen count |
| (c) β-function = Casimir trace identification | NEGATIVE | Different representation spaces; Hopf bundle has no SU(3)×SU(2) |
| (d) Hopf Chern class → SM charges | NEGATIVE | First Chern class is integer; cannot produce ±1/3 |
| (e) Δ rename | NEGATIVE | Tautological; no new structure |

**Verdict: FAIL.** The equality 8 = 8 is a numerical coincidence arising from the arithmetic of (1 + 4/3 + 1/3) · 3 vs. 3² − 1. All five candidate mechanisms fail.

---

## 9. References

- Paper 2 (GeoVac): `papers/conjectures/paper_2_alpha.tex`, esp. Sec. III (Casimir B = 42, selection principle B/N = 3, Δ = 1/(|λ_3|·N(2))).
- Paper 7 (GeoVac): `papers/core/Paper_7_Dimensionless_Vacuum.tex` (S³ conformal geometry, Fock projection).
- Beta function (physics), Wikipedia: https://en.wikipedia.org/wiki/Beta_function_(physics)
- Southampton QFT lecture notes, Sec. 3 (One-loop QED counterterms): https://www.southampton.ac.uk/~doug/qft/aqft3.pdf
- Ramond, "Journeys Beyond the Standard Model", Ch. 6 (Standard Model one-loop structure): https://www.phys.ufl.edu/~ramond/JourneysChapter6_CUP.pdf
- Lohitsiri, "Anomalies and the Standard Model" (Cambridge thesis): https://www.repository.cam.ac.uk/bitstreams/bf4b85dc-2fc1-4e24-8b34-a9b9b4a86bf4/download
- Tong, "Standard Model" lecture notes, Sec. 5 (Electroweak Interactions and anomaly cancellation): https://www.damtp.cam.ac.uk/user/tong/sm/standardmodel5.pdf
- Costa & Fogli, "Charge quantization in the standard model with three generations of fermions", Phys. Rev. D 45, 1701 (1992).
- Urbantke, "The Hopf fibration—seven times in physics": https://www.fuw.edu.pl/~suszek/pdf/Urbantke2003.pdf
- Scaletta, "Fiber Bundles, The Hopf Map and Magnetic Monopoles": https://www.math.columbia.edu/~ums/pdf/Scaletta%20LecNotes%20UMS%202-3-10.pdf
- Color charge, Scholarpedia: http://www.scholarpedia.org/article/Color_charge
