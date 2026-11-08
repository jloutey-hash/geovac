# Sprint memo — the configuration operator is the composed basis's metric (2026-08-25, v5.1.2)

**Origin.** `/aha` pass run immediately after v5.1.0/v5.1.1. Phase A produced 13 candidates
(6 cross-corpus, 3 distant-field, 4 inversions); Phase B killed or re-aimed most and left one
whose falsifier was cheap enough to run inside the pass. This memo records the run, the two
corrections it forced, and the follow-on tracks — including one I proposed and then killed.

**Drivers.** `debug/aha_config_operator_is_metric.py` · data
`debug/data/aha_config_operator_is_metric.json` · test
`tests/test_paper32_config_operator_metric.py` (15 fast, 2.8 s).

---

## 1. The result

The v5.1.0 configuration operator `F = Σᵢ Pᵢ` — whose conical intersections and π Berry phase
were written up as belonging to "the abstract configuration operator, NOT the electronic
Hamiltonian" — **is the composed basis's overlap matrix**, up to similarity.

**Proof (three lines).** Let `G` be the union-basis Gram, blocks `G_ij = ⟨χᵢ|χⱼ⟩`, with each
*intra*-center block equal to `I`. Put `X = chol(G)ᵀ`, so `XᵀX = G` and the columns of `X` are
the basis vectors in an orthonormal ambient frame; let `W_k` be center `k`'s columns. Then

1. `W_kᵀW_k = G_kk = I` ⇒ `W_k` has orthonormal columns;
2. the pseudoinverse in `P_k = W_k W_k⁺` collapses to `W_kᵀ`;
3. `F = Σ_k W_k W_kᵀ = X Xᵀ = X (XᵀX) X⁻¹ = X G X⁻¹`.

So `F ≃ G` — **similar**, not merely cospectral. Every spectral statement about `F` (the CIs,
the band gaps, the Z₂ holonomy of its eigenvectors) is a statement about the **overlap metric
of the composed basis** over nuclear geometry.

**The hypothesis is real, not a technicality.** Intra-center orthonormality holds throughout
GeoVac: automatically for distinct-`l` sets on one center, and *exactly* for the
Shibuya–Wulfman metric (Paper 60 §325; measured `max|intra − I| = 7.07e-17`,
`aha_track1_findings.md`).

### Verification

| case | max·spec F − spec G· | ‖FX − XG‖_max |
|:--|:--|:--|
| linear BeH₂, 8 geometries incl. both CIs | ≤ 3.1e-15 | ≤ 5.6e-16 |
| bent BeH₂, 160°→90° incl. H₂O 104.5° | ≤ 2.7e-15 | ≤ 6.7e-16 |
| four centers (rank-3 per center) | 2.2e-15 | 8.9e-16 |

**Non-tautology control.** Perturbing one intra-center block off the identity by 0.3 breaks the
identity: `|Δspec|` 1.3e-15 → **6.7e-2**, six orders above tolerance. Step (1) is load-bearing.

**Four centers, blockwise.** `Tr(P₁P₂P₃P₄) = Tr(G₁₂G₂₃G₃₄G₄₁)` to 3.3e-16 — i.e.
`rem:four_center_modulus` had *already written* the holonomy on the overlap blocks `S_ij`
without noticing it was an instance of the general identity. Half-surfaced, now general.

---

## 2. Two corrections the run forced

**(a) The strong form of the lead is FALSE — and it is pinned as such.** I expected the CIs to
mark a conditioning blow-up, connecting to Paper 60's metric wall. They do not. The crossing
is **mid-spectrum** (eigenvalue ≈ 0.714 against `λ_min` ≈ 0.129) and `cond(G)` passes through
`d*` **monotonically**: 19.87 → 16.19 → 15.04 across (2.30, 2.445, 2.50). Conditioning and
holonomy are different functionals of the same metric. The paper says so explicitly so the
claim cannot regrow.

**(b) A measurement bug in my own promoted driver.** The scratchpad version reported the Löwdin
*return* drift; on promotion I rewrote it as a max over the loop, which measures **variation**,
not single-valuedness — and variation is nonzero for any non-constant loop (largest, in fact,
on the *CI-free* control, 1.08e-01). Both are now reported separately. The corrected statement:

| loop | eigvec return | Löwdin RETURN | (loop variation) |
|:--|--:|--:|--:|
| encircles central CI | **−1.0000** | 0.00e+00 | 4.69e-02 |
| encircles central CI (tight) | **−1.0000** | 0.00e+00 | 1.87e-02 |
| CI-free control | **+1.0000** | 0.00e+00 | 1.08e-01 |

**The orthogonalization criterion.** The eigenvectors of the *metric* are double-valued around
a CI. Any scheme that **selects** metric eigenvectors — canonical orthogonalization with a drop
threshold, the standard near-linear-dependence remedy — is path-dependent along a geometry scan
through such a point. Symmetric Löwdin `G^{-1/2}`, a function of the matrix rather than of its
eigenvectors, returns bit-exactly. This is a usable statement about composed-basis PES scans.

---

## 3. Tame ≠ composable (the scoping correction)

Phase A proposed that the multi-focal wall's threshold is 2→3 focal lengths, mirroring the
tame→wild threshold. **Refuted:** the multi-focal wall fires at **two** (Paper 57 §614: two
focal lengths give `c_ab ≠ c_a·c_b`; every instance in `multi_focal_wall_pattern` is 2-focal —
two electrons each with their own scale, proton magnetization × electron 1s, etc.).

What survives is a scoping correction the corpus needed anyway. The unified statement
"composition wall = non-commuting projections" now covers **two structurally different rungs**:

- **Rung 1 (two centers).** Projections do not commute — `‖[P_A,P_B]‖ = 1/2`, the maximal
  value — but are *completely classified* by principal angles (Halmos). The obstruction is
  quantitative and the cost of forcing commutation is known in closed form. **This is what the
  multi-focal wall is.**
- **Rung 2 (≥3 centers).** The classification itself fails (*-wild); the molecular triple sits
  at an irreducible `M₆` point. There is no normal form available to force.

Corollary worth recording: **no 3-focal-length precision observable has ever been tested in the
corpus.** The one ≥3-focal object that exists is the TC three-body operator (v5.0.2: non-abelian,
shared-vertex rank 4→9→16, does not collapse), which behaves as wildness would predict — but
that is a vertex angular-coupling object and this is a center-subspace object, so the
resemblance is a shape match pending an axis check, not a claim.

---

## 4. Killed, with reasons (banked for `/walls`)

- **"Commutant = the tapering-symmetry budget."** *Proposed and killed in this sprint.* The idea
  was that the commutant collapse (LiH dim 6 → BeH₂ dim 1) explains the §3 dead-end "non-abelian
  M-vS gauge reduces Pauli count = NEGATIVE". Two independent facts kill it: (i) that negative
  was **measured on LiH — two centers**, where the commutant is dim 6, so a three-center collapse
  cannot explain it; (ii) the Hopf Z₂ tapering **works on BeH₂** (verified: `hamiltonian('BeH2',
  tapered='global')` succeeds) where the center-projection commutant *is* trivial. So a
  nontrivial Hamiltonian symmetry coexists with a trivial center-projection commutant — the two
  commutants are different objects, and `commutant(A) ⊆ commutant(H)` is a bound in the useless
  direction.
- **"F's CIs mark a conditioning blow-up."** Killed by §2(a): `cond(G)` monotone through `d*`.
- **"Multi-focal threshold is 2→3."** Killed by §3: it fires at 2.
- **"Chemistry error tracks three-body irreducible content."** LiH has 5.3% `R_eq` error with
  *no* three-body term at all, and CHEM-ACCURACY localizes the defect at 100% `max_n`
  (angular-invariant). Wrong axis: angular-operator structure vs radial basis size.
- **"α's K = π(B+F−Δ) is a tame/wild statement."** Downgraded, not run: Paper 2
  §sec:rigidity_connection already states the 2+1 split (B, F on scalar S³; Δ on spinor), and
  it is an axis error regardless — B, F, Δ are numbers, not projections, so nothing generates
  an algebra and nothing can be wild.

**Distant-field tier yielded re-description, as the skill predicts:** Löwdin orthogonalization
is *equal temperament* (force the circle to close by smearing the comma uniformly; lose just
intonation exactly as you lose sparsity); the nested commutator is the three centers'
*Condorcet cycle*.

---

## 5. Applied

- **Paper 32** `rem:config_operator_is_metric` (new): the identity, its proof, the scope
  hypothesis, the non-tautology control, the blockwise four-center form, and both consequences
  — including the explicit negative so the conditioning claim cannot regrow. Cites Paper 60
  (new bibitem).
- **Paper 32** `rem:multicenter_composition`: new paragraph carrying the two-rung scoping
  ("tameness is not composability").
- **Tests** `tests/test_paper32_config_operator_metric.py` — 15 fast (2.8 s): molecular +
  abstract-random identity, the discriminating non-tautology control, the four-center blockwise
  holonomy (with a non-vacuity guard), the orthogonalization loops with a CI-free control, and
  the honest negative pinned as its own test.
- Paper 32 compiles clean (88 pp, exit 0, no undefined refs); C19 escape gate PASS; duration /
  retracted-terms / K-label / paper-test-refs gates PASS; 47 passed + 1 skipped across the four
  Paper 32 suites plus the 18 topological-integrity proofs.

## 6. The two named follow-ons — both RUN, both NEGATIVE, one with a finding underneath

### 6.1 Axis check: TC three-body vs ≥3-projection wildness — NOT the same obstruction

Driver `debug/aha_axis_check_tame_vs_wild.py` · data
`debug/data/aha_axis_check_tame_vs_wild.json`. Killed two independent ways.

**(a) The `/aha` reading misread the TC numbers.** The memo's "rank 4 → 9 → 16" runs over the
columns **l=1 | l=2 (L_corr=2) | l=2 (L_corr=3)** — growth in *angular momentum / basis*, not in
leg count. And the non-collapse is present at the smallest basis: **87.5% of external pairs
already have shared-vertex rank ≥ 2 at l=1**. The TC obstruction fires at **two** correlator legs
(`Y·Y` is already a multiplet). There is no 2→3 threshold on that axis at all.

**(b) A theorem separates them.** The TC obstruction is compact-group representation theory
(SO(3)/SU(2) coupling), and by **Peter–Weyl every compact group's representation algebra is
type I — tame**, however non-abelian. Three orthogonal projections are *not* a compact-group
representation; arbitrary subspace configurations escape Peter–Weyl, which is how they reach the
wild regime. Measured (commutant dims are predicted exactly by Schur, so the routine
self-validates):

| object | dim | max‖[A,B]‖ | commutant | verdict |
|:--|--:|--:|--:|:--|
| SU(2) on l = 0⊕1⊕2 | 9 | **2.000** | 3 | REDUCIBLE |
| SU(2) on l = 1⊕1 (multiplicity 2) | 6 | 1.000 | 4 | REDUCIBLE |
| SU(2) on l = 1⊕2⊕2 | 13 | 2.000 | 5 | REDUCIBLE |
| SU(2) on l = 1 alone | 3 | 1.000 | 1 | irreducible (single irrep) |
| BeH₂ **two** centers | 6 | 0.487 | 6 | REDUCIBLE (Halmos tame) |
| BeH₂ **three** centers | 6 | **0.487** | **1** | **IRREDUCIBLE** |

The SU(2) cases are *more* non-abelian than the projections and stay reducible; the three
projections are irreducible at *smaller* non-commutativity. **Non-abelianness and wildness are
independent properties.**

**The finding underneath (worth keeping).** This explains an existing corpus fact with a
mechanism: the framework's angular sparsity is graded, finite and **potential-independent**
(Paper 22 — depends only on `l_max`, never on `V(r)`) *because* it is compact-group
representation theory. And it is decision-relevant in the other direction: **no Gaunt-style
graded reduction should be expected to rescue multi-center composition**, because that axis is
not graded — it is wild. Written into Paper 32 `rem:multicenter_composition`; backed by
`tests/test_paper32_composition.py::test_compact_group_coupling_is_tame` (4 parametrized,
self-validating against Schur) and `::test_nonabelian_and_wild_are_independent`.

### 6.2 The WH7 Lorentzian-leg reading — KILLED by data already in the corpus

The proposal: WH7's compact modular circle Z₂ (spinor −1 vs scalar +1 at β/2, the μ₄ Hodge
circle) is the T-symmetric equator of the same monopole structure, with Wick rotation playing
the role of the U(1) lift. **No new computation was needed — the falsifier had already fired.**

The BeH₂ U(1) lift *does something*: a Peierls phase opens the empty σ_y axis with gap ∝ φ
(slope ≈ 0.55), converting Z₂ into Chern ±1. The analogy therefore predicts that a temporal Wick
rotation must open something on the modular circle. The measured result is the opposite:
`debug/sprint_lorentzian_toeplitz_kplus_probe_memo.md` records that a **genuine temporal Wick
involution `J_t = sign(D_t)` leaves the surviving seminorm bit-exactly `ω_q`** — unchanged from
the Euclidean value — and this is frozen in `tests/test_lorentzian_toeplitz_kplus.py`
(re-run this session: 11 passed). Signature-blind, bit-exactly. A lift that changes nothing is
not a lift.

**Mechanism of the error — a Z₂ pun**, and precisely the axis category error the axis map warns
about. The two Z₂'s have different sources: the BeH₂ one is Berry holonomy sourced by a
**codimension-2 degeneracy** in a ≥2-parameter family with a nonzero branching Jacobian; the
modular one is the **SU(2) double cover** (half-integer `m_j` ⇒ spinor picks up −1 at β/2). The
modular generator `K = diag(2m_j)` has a *fixed* integer spectrum at every cutoff — there is no
parameter family, no degeneracy, and hence nothing for a monopole to sit at. This is a fifth
member of the family the ℚ(i) seam audit already closed (four candidates, verdict CONVERGENT:
"each grabbed one descendant of the structure, not an ancestor").

## 7. Open

- Nothing from §6 — both follow-ons are closed. The residual open items are the pre-existing
  ones (Paper 32 Phase-4 re-review OWED; the `/walls` clusters still unworked).
