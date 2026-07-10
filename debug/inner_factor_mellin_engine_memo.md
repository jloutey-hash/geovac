# Master Mellin Engine on Finite Spectral Triples (Inner-Factor Probe)

**Sprint:** Inner-factor Mellin engine investigation (May 2026)
**Driver:** `debug/finite_spectral_triple_engine.py`
**Data:** `debug/data/finite_spectral_triple_engine.json`
**Verdict:** **Mixed-positive: structural theorem on finite triples (η-class trivializes by construction); engine output reduces to Yukawa Dirichlet series; cuts the Krajewski candidate set but does NOT pick A_F uniquely.**

---

## 1. Question

Sprint TS-E1 (May 2026) closed the master Mellin engine on the *outer* (continuum or large-N truncation) GeoVac spectral triple as a case-exhaustion theorem: every π-source in any finite chain of Paper 34 projections engages one of three sub-mechanisms

- **M1** (k=0): Hopf-base measure / propinquity rates / 4/π = Vol(S²)/π² signature
- **M2** (k=2): Seeley–DeWitt / heat-kernel / spectral-action coefficients in √π·ℚ ⊕ π²·ℚ
- **M3** (k=1): vertex-parity Hurwitz / Dirichlet-L content / Catalan G, β(4)

each a sub-case of the master Mellin transform 𝓜[Tr(D^k · e^{-tD²})] at operator order k ∈ {0, 1, 2}.

**Question:** Does this engine extend to *finite* spectral triples — specifically, the inner factor A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ) of the Connes–Chamseddine Standard Model? If so, does it constrain or reverse-engineer A_F?

This is the natural follow-up to Sprint H1's "POSITIVE-THIN" verdict: H1 showed GeoVac's outer triple admits a Higgs structurally, but the Yukawa data must be supplied empirically. The question here is whether GeoVac's *taxonomic* machinery (not just the AC construction) says anything about which finite triples are eligible inputs.

---

## 2. Leg 1 — Engine on canonical finite triples

For a finite Hermitian D_F with spectrum {(λ_i, m_i)} (eigenvalue, multiplicity), the trace is a finite sum

  Tr(D_F^k · e^{-tD_F²}) = Σ_i m_i · λ_i^k · e^{-tλ_i²}

and its Mellin transform structure (with eigenvalues > 0) is

  𝓜[Σ_i c_i · e^{-tλ_i²}](s) = Γ(s) · Σ_i c_i · λ_i^{-2s}

i.e., a generalized Dirichlet series in the eigenvalues.

We tested four canonical finite spectral triples symbolically:

| Triple | dim H_F | Spectrum |
|---|---|---|
| SM lepton sector (ℂ ⊕ ℍ, one gen) | 8 | ±y_ν, ±y_e (each ×2) |
| Full SM (ℂ ⊕ ℍ ⊕ M_3(ℂ), one gen) | 32 | ±y_ν, ±y_e (×2); ±y_u, ±y_d (×6 color) |
| Comparator: ℂ ⊕ ℂ | 2 | ±m |
| Comparator: M_2(ℂ) ⊕ M_3(ℂ) (off-diag rank-2) | 5 | ±σ_1, ±σ_2, 0 |

**Result (verified symbolically in `debug/data/finite_spectral_triple_engine.json`):**

| Triple | k=0 (Tr e^{-tD²}) | k=1 (Tr D·e^{-tD²}) | k=2 (Tr D²·e^{-tD²}) |
|---|---|---|---|
| SM lepton | 4·e^{-ty_ν²} + 4·e^{-ty_e²} | **0** | 4y_ν²·e^{-ty_ν²} + 4y_e²·e^{-ty_e²} |
| Full SM | 4(e^{-ty_ν²} + e^{-ty_e²}) + 12(e^{-ty_u²} + e^{-ty_d²}) | **0** | 4(y_ν²·e^{-ty_ν²} + y_e²·e^{-ty_e²}) + 12(y_u²·e^{-ty_u²} + y_d²·e^{-ty_d²}) |
| ℂ ⊕ ℂ | 2·e^{-tm²} | **0** | 2m²·e^{-tm²} |
| M_2 ⊕ M_3 | 1 + 2e^{-tσ_1²} + 2e^{-tσ_2²} | **0** | 2σ_1²·e^{-tσ_1²} + 2σ_2²·e^{-tσ_2²} |

**The k=1 trace vanishes identically on all four cases.** This is not a numerical coincidence — it is structural and we prove it cleanly in Leg 2.

The k=0 and k=2 outputs are linear combinations of Gaussians at the Yukawa-eigenvalue scales with rational coefficients (multiplicities). The Mellin output is a Dirichlet sum

  Γ(s) · Σ_i (rational) · y_i^{-2s} (k=0)
  Γ(s) · Σ_i (rational) · y_i^{2-2s} (k=2)

which lives in the ring **ℚ[y_1^{-2s}, y_2^{-2s}, …, y_n^{-2s}]** — a *parameter-tied Dirichlet ring*, not in any of M1/M2/M3 from the outer-factor classification.

---

## 3. Leg 2 — η-vanishing as a structural theorem

**Theorem (η-class trivialization on Krajewski-class finite triples).** Let (A_F, H_F, D_F, J_F, γ_F) be a finite even spectral triple in the Krajewski–Paschke–Sitarz lineage (Krajewski 1998, Paschke–Sitarz 2000). The chirality grading γ_F satisfies {γ_F, D_F} = 0 (a defining axiom of an even spectral triple). Then for every k ∈ ℕ_odd,

  Tr(D_F^k · e^{-tD_F²}) ≡ 0  ∀ t > 0.

In particular, **𝓜[Tr(D_F · e^{-tD_F²})] ≡ 0** — the M3 sub-mechanism trivializes.

**Proof.** {γ_F, D_F} = 0 ⟹ γ_F · D_F · γ_F^{-1} = -D_F. Since γ_F² = 1, γ_F is a unitary involution. Conjugation by γ_F is a similarity transformation on H_F that sends D_F → -D_F while preserving traces. So

  Tr(D_F^k · e^{-tD_F²}) = Tr(γ_F · D_F^k · e^{-tD_F²} · γ_F^{-1})
                        = Tr((-D_F)^k · e^{-tD_F²})
                        = (-1)^k · Tr(D_F^k · e^{-tD_F²}).

For odd k, this forces Tr(D_F^k · e^{-tD_F²}) = -Tr(D_F^k · e^{-tD_F²}), hence zero. □

**Equivalent statement:** σ(D_F) is symmetric about 0 (each λ paired with -λ via γ_F). This was verified directly in Leg 2 of the JSON: Σ_i λ_i · m_i = 0 in all four cases.

**What this means for the Mellin engine:**

The case-exhaustion theorem of Sprint TS-E1 covered three operator orders k ∈ {0, 1, 2}. On the outer GeoVac triple, all three are non-trivially populated — k=0 by Hopf-base measure, k=1 by half-integer Hurwitz on the Camporesi–Higuchi spinor spectrum, k=2 by Seeley–DeWitt on the round-S³ Dirac. **On a finite spectral triple in the Krajewski lineage, only k ∈ {0, 2} are non-trivial**; k=1 is structurally zero by chirality grading.

This is a *real constraint*. The Connes–Chamseddine SM construction *requires* an even finite spectral triple (γ_F is the chirality matrix that distinguishes left from right and eventually becomes the weak-isospin γ_5 on the matter representation). So **any inner factor compatible with Connes-Chamseddine SM has trivial M3 content.**

The two-mechanism reduction is structural to the inner factor of any almost-commutative theory. The outer factor retains all three.

---

## 4. Leg 3 — Combined outer × inner factorization

For T_combined = T_GV ⊗ T_F with combined Dirac

  D = D_GV ⊗ 1_F + γ_GV ⊗ D_F

(the Connes–Chamseddine convention, also used in our Sprint H1 implementation),

  D² = (D_GV ⊗ 1)² + (γ_GV ⊗ D_F)² + {D_GV ⊗ 1, γ_GV ⊗ D_F}

The cross term is

  {D_GV ⊗ 1, γ_GV ⊗ D_F} = (D_GV γ_GV + γ_GV D_GV) ⊗ D_F = 0

because γ_GV anticommutes with D_GV (it is the chirality grading of the outer triple). So

  D² = D_GV² ⊗ 1_F + 1_GV ⊗ D_F².

**Verified symbolically** in `leg_3_combined_factorization` of the JSON output: for a 2×2 toy outer × 2×2 inner construction, D_total² produces exactly D_GV² ⊗ 1 + 1 ⊗ D_F² with zero residual.

**Consequence for the engine:**

  e^{-tD²} = e^{-tD_GV²} ⊗ e^{-tD_F²}
  Tr(e^{-tD²}) = Tr(e^{-tD_GV²}) · Tr(e^{-tD_F²})

The combined-engine Mellin output factorizes as

  𝓜[Tr(D^k · e^{-tD²})](s) = (engine on outer at order k) ★ (engine on inner at order k)

(where ★ is the Mellin convolution / Dirichlet product). For k=0 specifically, the outer factor lives in M1/M2 rings (depending on bundle) and the inner factor contributes Γ(s) · Σ_i c_i · y_i^{-2s}. The combined output is therefore

  outer-π-content × Yukawa-Dirichlet-content

a *product structure*, not a sum. The π-content of the SM spectral action comes entirely from the outer factor; the SM-distinguishing parameters (Yukawas) enter only through the inner Dirichlet sum. **This is exactly the Connes–Chamseddine bosonic spectral action structure**: gravity + gauge from outer Seeley–DeWitt, Yukawa-dependent terms from the inner trace.

The structural alignment with the GeoVac outer-factor taxonomy is clean: π's sit in their proper M_i ring, Yukawas are decoupled, and the combined output is a simple product.

---

## 5. Leg 4 — KO-dim restriction and Krajewski candidate cut

**Sprint H1 (verified in `geovac/almost_commutative.py` and `tests/test_almost_commutative.py`):** combined KO-dim 9 ≡ 1 (mod 8) gives ε(combined) = -1, ε'(combined) = +1, hence J² = -I and JD = +DJ exactly on truthful CH at n_max ∈ {1, 2, 3}. This requires inner KO-dim 6 (the standard SM convention, Connes–Marcolli Table 13.1).

**Krajewski classification at KO-dim 6 (mod 8):** the finite spectral triples of KO-dim 6 form a constrained sub-class. Up to small-dimensional cases:

- ℂ ⊕ ℂ — KO-dim 6 with appropriate J_F (matter–antimatter swap × conj).
- ℂ ⊕ ℍ — the "electroweak slice" of Sprint H1; KO-dim 6 with J_F = i σ_2 K on H.
- M_2(ℂ) ⊕ M_3(ℂ) — yes (with appropriate matter/antimatter doubling).
- ℂ ⊕ ℍ ⊕ M_3(ℂ) — Connes' SM choice; KO-dim 6.

Krajewski's classification is large but combinatorial: each candidate has a Krajewski diagram (bipartite graph encoding bimodule structure). The KO-6 cut is necessary but does not single out the SM choice. Within KO-6, the SM has the smallest non-trivial real-dimensional algebra producing color SU(3) plus electroweak SU(2) × U(1). It is "Connes–Chamseddine minimal" but not "Krajewski-unique".

**What the Mellin engine analysis adds beyond KO-dim:**

(a) **η-trivialization** (Leg 2) is structural to any Krajewski-class triple, not specific to SM. So the engine does *not* distinguish SM from other KO-6 candidates by k=1 content.

(b) **k=0 / k=2 distinguishing power** is bounded by the spectrum-multiplicity structure. SM A_F has spectrum {y_ν, y_e, y_u, y_d} with multiplicities {2, 2, 6, 6}. The **multiplicity ratio 6:2 = 3:1 is the SU(3) color factor.** This is the only place where "color = 3" appears in the trace structure. A KO-6 triple with simpler algebra (e.g., ℂ ⊕ ℍ alone) gives multiplicity 2:2 = 1:1, no color factor.

(c) **GeoVac-compatibility test:** when tensored with the GeoVac outer triple, the combined spectral action

  S = Tr f(D/Λ) = Tr f(√(D_GV² + D_F²)/Λ)

picks up Seeley–DeWitt coefficients on the outer side weighted by inner-trace factors. For SM A_F the inner trace at order Λ^4 gives 8 · (1+1+3+3) = 16 species (counting chirality), reproducing the standard "Λ^4 cosmological term" with the SM color factor 3 baked in. **This integer 3 is the same SU(3) color factor that Connes–Chamseddine derive from spectral action coefficient matching** (Chamseddine–Connes–Marcolli 2007, Adv. Theor. Math. Phys. 11). It is not new content from GeoVac, but it confirms the engine on GeoVac × A_F reproduces the SM normalization.

**Verdict on Leg 4:** KO-6 cut is a real combinatorial constraint inherited from H1. Engine analysis does NOT further reduce candidates beyond what Connes–Chamseddine already do (spectral action coefficient matching). Color factor 3 appears as the multiplicity-ratio in the engine output, consistent with but not predicted by GeoVac.

---

## 6. Three things the engine does say

After the four legs:

**(i) Two-mechanism reduction is structural.** Any inner factor compatible with Connes–Chamseddine has M3 ≡ 0. So the master Mellin engine on the *combined* outer × inner triple supports M1/M2 from the outer × Yukawa-Dirichlet from the inner, with M3 only on the outer. This is a clean partition.

**(ii) Yukawas are the inner factor's "calibration exchange constants."** The inner Mellin output is a generalized Dirichlet series in the Yukawa eigenvalues. Under Paper 18's exchange-constant taxonomy, this places Yukawa values in their own *parameter-tied tier*: not intrinsic (depends on free parameters), not calibration in the ζ_R(2k) sense (the ring is a free polynomial in y_i^{-2s}, not a fixed transcendental algebra), not embedding (no continuum integration). They constitute a new fourth tier — **"inner-factor input data"** — that the outer-factor taxonomy does not predict.

This sharpens H1/G3/G4a's "no autonomous SM-distinguishing data" finding: the autonomous data lives in the outer factor's algebraic-extension rings (ℚ[α²]/(γ²+(Zα)²-1) etc.), the SM-distinguishing data lives in the inner factor's Yukawa Dirichlet ring, and **the two rings have no overlap.** GeoVac's taxonomy controls the first; A_F controls the second; they tensor cleanly but neither determines the other.

**(iii) The spectral action factorization confirms the architecture.** D² = D_GV² ⊗ 1 + 1 ⊗ D_F² (cross term zero by chirality anticommutation) means the combined Mellin engine is multiplicatively separable. The π's of the SM spectral action come *entirely* from the outer factor — exactly what the Paper 35 / Paper 34 / Paper 18 taxonomy is built for. The Yukawas appear as "scaling shifts" in the spectral-action coefficients, not as new transcendentals.

---

## 7. What the engine does NOT do

- **Does not pick A_F uniquely.** KO-6 + η-vanishing + GeoVac-compatibility leaves a large Krajewski subclass. SM A_F is one minimal choice; other KO-6 triples (e.g., M_4(ℂ), various extensions) pass the same tests.
- **Does not predict Yukawa values.** The Dirichlet ring has free parameters; Connes' spectral action principle gives mass relations *at the unification scale* but Yukawas at low energy require running. Engine analysis is silent on this.
- **Does not provide a "second packing axiom."** GeoVac's packing axiom is a 2D→3D lift producing (n,l,m,s). No analogous axiom is here for color or generations. The empirical A_F choice remains the load-bearing input.
- **Does not derive generations.** Three generations is not a consequence of the engine; it is an additional empirical multiplicity factor in D_F.

---

## 8. Honest verdict

**MIXED-POSITIVE.**

Real structural findings, in order of strength:

1. **(POSITIVE)** η-vanishing is a structural theorem on any Krajewski-class finite spectral triple. The master Mellin engine *necessarily* has a 2/3 reduction on the inner factor (M1, M2 only; M3 ≡ 0).

2. **(POSITIVE)** Combined outer × inner Mellin engine factorizes multiplicatively. The π-content sits entirely in the outer factor (M1/M2/M3 rings); the SM-distinguishing data sits entirely in the inner factor (Yukawa Dirichlet ring); the two rings are orthogonal and the combined output is a clean product.

3. **(POSITIVE-THIN)** Sharpens H1/G3/G4a's "no autonomous SM-distinguishing data": the inner factor's Dirichlet ring is a *fourth tier* in Paper 18's exchange-constant taxonomy ("inner-factor input data") that the outer-factor taxonomy categorically does not predict.

4. **(NEGATIVE)** Engine does not pick A_F uniquely. KO-6 + η-vanishing + factorization compatibility leaves many candidates; SM A_F is one minimal choice but not *the* choice. No reverse-engineering of the inner geometry.

5. **(NEGATIVE)** No Yukawa prediction. The free parameters of the SM survive the engine analysis intact.

**Overall:** the master Mellin engine *does* extend to finite spectral triples and *does* give two real structural results on them (η-trivialization, multiplicative factorization). It does not reverse-engineer A_F. The honest position is that GeoVac's taxonomy provides:

- a **structural vocabulary** for placing SM parameters in a tier ("inner-factor input data") parallel to but disjoint from the outer-factor exchange constants;
- a **structural theorem** (η-trivialization) that constrains finite-triple candidates;
- a **factorization** that confirms the Connes–Chamseddine architecture is internally consistent with the master Mellin engine framework.

It does not provide what the user (PI) was hoping for: a derivation principle for A_F or its Yukawa data. That remains open. **The "second packing axiom" question is genuinely unsolved by this analysis.**

---

## 9. Recommended follow-ups

If the PI wants to push further, three tractable next moves, in order of payoff:

**(A) Paper 18 §IV taxonomic extension:** add a fourth tier to the exchange-constant grid for "inner-factor input data" — the parameter-tied Dirichlet ring of Yukawas. Cite this memo. Modest writing task, ~1 week. Sharpens the universal/Coulomb partition (Paper 31) into a three-way universal/Coulomb-specific/inner-factor partition.

**(B) Krajewski–Mellin compatibility audit:** systematically compute Tr(D_F^k · e^{-tD_F²}) for the small Krajewski classification (dim H_F ≤ 16 say, KO-6 only) and rank candidates by GeoVac-compatibility metrics (color factor reproduction, multiplicity ratio structure, spectral action normalization). ~2-3 weeks. Could land a "GeoVac selects a small subclass" result if the metrics are picky enough — or a clean "no further selection" negative.

**(C) Yukawa-ratio numerology audit:** compute the Mellin Dirichlet sum for SM A_F at integer s and test whether observed Yukawa ratios (m_t/m_b, m_τ/m_e, etc.) plug into algebraic identities. Koide's formula territory. **High risk of pareidolia, low expected payoff. Not recommended without strong prior.**

The PI's directive on (A) and (B) determines whether this thread continues. (C) is genuinely speculative.

---

## 10. Files

- `debug/finite_spectral_triple_engine.py` (computational driver, ~250 lines)
- `debug/data/finite_spectral_triple_engine.json` (symbolic results)
- `debug/inner_factor_mellin_engine_memo.md` (this file, ~3500 words)

No commits made. No `geovac/` modifications. No paper edits.

---

## 11. Out-of-scope flags

- **Multi-generation extension:** all calculations done at one generation. The 3-generation case multiplies multiplicities by 3 but does not change the structural findings (η still vanishes, ring is still parameter-tied Dirichlet).
- **Operator-order-3 (k=3) check:** by η-vanishing argument, ALL odd k vanish on Krajewski-class triples. So the engine on finite triples really has only two non-trivial operator orders, not three. The "two-mechanism reduction" is sharper than the "2/3 reduction" stated above; this could be folded into the §6(i) bullet if the PI agrees.
- **Connection to Paper 35 (time-as-projection):** the inner Yukawa Dirichlet ring is a *static* spectral object; it does not contain a temporal-window projection. So Paper 35's "π enters iff continuous integration over time/spectral parameter" still holds — the inner factor's lack of π is consistent with that prediction, since there is no continuous parameter to integrate over.
