# Sprint C3 — $N_{\mathrm{gen}}$ non-selection theorem (2026-06-08)

## TL;DR

**Theorem-grade upgrade of the Direction 2 + Read 2 NO-GO scopings.** Paper 32 §VIII gains `thm:n_gen_non_selection` + `rem:n_gen_scope`. Statement: under the canonical Chamseddine--Connes SM representation, for every integer $N \ge 1$ the GeoVac spectral triple with inner Hilbert-space multiplicity $N$ satisfies every axiom in the framework's current set, hence no axiom selects a specific value of $N$. Conditional on canonical representation; unconditional version (any SM-phenomenologically-consistent rep correlates $N$ with algebra factor count) remains the open multi-year Read-2-sharpest-falsifier question. Paper 32 compiles three-pass clean at 76 pages.

Companion to H1 Yukawa non-selection (the only other structural non-selection theorem in the corpus prior to today) and the Forced-Count Theorem. Together these three theorems characterise the full free-side content of the inner-factor structural-skeleton boundary at the canonical-representation level.

## Sprint context

Sprint C3 was framed as "upgrade 'no handle' to 'structurally cannot.'" The empirical non-selection results in the corpus:

| Entry | Empirical status before C3 | Status after C3 |
|:------|:---------------------------|:----------------|
| F1 Yukawa values | THEOREM (H1 + Yukawa-PSLQ) | unchanged |
| **F2 $N_{\mathrm{gen}}$** | DIAGNOSTIC (Direction 2 + Read 2 NO-GO) | **THEOREM (this sprint)** |
| F3 Inner KO-dim | DIAGNOSTIC (Direction 2) | DIAGNOSTIC (named for future upgrade) |
| H6, H7 Chem inner data | DIAGNOSTIC (W1e period-class) | DIAGNOSTIC |
| E6 Combination rule $K$ | EMPIRICAL (12 mechanisms eliminated) | EMPIRICAL |
| E7, E8 LS-8a counterterms | EMPIRICAL | EMPIRICAL |
| G1–G6 Multi-focal walls | EMPIRICAL (6 instances) | EMPIRICAL |
| D5, D6 CC fine-tuning | EMPIRICAL ($\varphi(2)/\varphi(1)^2 \sim 10^{-124}$) | EMPIRICAL |
| F4–F7 Higgs VEV, CKM, PMNS, $\nu$ | EMPIRICAL (inherit F1/F4) | EMPIRICAL |

F2 was the most tractable upgrade target. The two NO-GO scopings (Direction 2 packing-reach + Read 2 representation-level) supplied the structural ingredients; formalisation into a Theorem block was the remaining work.

## Theorem statement (paraphrase; see Paper 32 `thm:n_gen_non_selection` for the formal version)

**Theorem.** Let $\mathcal{A}$ denote the framework's current axiom set:
1. Paper 0 packing axiom
2. Standard real-spectral-triple axioms (Hermiticity, $J$-reality, KO-dim chirality, order-one)
3. Hopf-tower truncation (Hurwitz, Doors 4b/4d/4e)
4. Bertrand's classical theorem (closed-orbit selection)
5. Upgrade B (sphere-Lie-group axiom)

Fix the canonical CCM SM representation of $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$ on the inner Hilbert space. For every integer $N \ge 1$, the GeoVac spectral triple $\mathcal{T}_N$ with $\dim \mathcal{H}_F = 32N$ satisfies every axiom in $\mathcal{A}$. Consequently $\mathcal{A}$ contains no morphism that selects a specific value of $N$.

**Proof structure.** Six bullets in Paper 32, one per axiom, each showing $N$-independence. The five framework axioms (Paper 0, standard RST, Hopf-tower, Bertrand, Upgrade B) operate on the outer manifold $S^3$ and the inner algebra $\mathcal{A}_F$, neither of which carries $N$ as a parameter. The canonical CCM rep distributes 32 fermion DOF per generation across the three algebra factors, with the multiplicity $N$ entering as an external Hilbert-space multiplicity rather than an emergent constraint.

The Read 2 scoping content enters as the structural closure of any candidate Hopf-tower-to-representation shortcut: three independent obstructions to correlating $N$ with algebra-factor count (DOF count, CKM/PMNS coupling structure, alternative-assignment phenomenology), each independently closing the shortcut.

## What this is and is not

**What it is.** A conditional structural non-selection theorem for $N_{\mathrm{gen}}$ under the canonical CCM SM representation choice. Theorem-grade in the same sense as H1 Yukawa non-selection: explicit moduli space ($\mathbb{N}_{\ge 1}$ instead of the H1 128-dim diagonal slice), structural argument that no axiom selects a point, three independent obstructions to candidate counter-arguments.

**What it isn't.** An unconditional impossibility theorem. The Read 2 sharpest-falsifier remains:\ exhibit an SM-phenomenologically-consistent representation of $\mathcal{A}_F$ on $\mathcal{H}_F$ where (a) the 32 fermion DOF per generation distribute across the three algebra factors in a way that makes the three generations the three rep blocks, and (b) the Yukawa-Higgs cross-generation coupling (CKM/PMNS) is naturally encoded as the algebra cross-factor coupling. No such construction is currently exhibited; constructing one (or proving it impossible) is multi-year structural-research grade.

The Theorem is therefore phrased as conditional on the canonical CCM rep, with the conditional made explicit in `rem:n_gen_scope`.

## Cross-references

**Updates Paper 32:**
- `thm:n_gen_non_selection` added to §VIII after `rem:forced_count_with_seam`
- `rem:n_gen_scope` added pairing the conditional theorem with the unconditional open question
- `\paragraph{Files.}` extended with cross-reference to this memo

**Updates needed in Paper 57:**
- §6.1 (open question on P5 tautology / non-triviality): add note that the F2 (N_gen) and F1 (Yukawa) catalogue entries are now both theorem-grade non-selection, strengthening the structural reading of P5's failure-mode decomposition for the inner-factor input-data family
- Catalogue entry F2 updated to note theorem-grade status

## What's left in the corpus on the C3 non-selection arc

Three theorem-grade non-selection results now exist:
1. **H1 Yukawa non-selection theorem** (Paper 32 §VIII.C, pre-C3): structurally forces 8 free Yukawa parameters per generation in the diagonal slice
2. **Forced-Count Theorem** (Paper 32 §VIII `thm:forced_count`, pre-C3): forces $\dim \mathcal{M}(D_F) = 128$ per generation
3. **$N_{\mathrm{gen}}$ non-selection theorem** (Paper 32 §VIII `thm:n_gen_non_selection`, C3, this sprint): formal structural verdict on the multiplicity factor

Remaining empirical / diagonostic non-selection entries that could be candidate future upgrade targets:
- **F3 inner KO-dimensional signature**: closely parallel to F2; same Direction 2 NO-GO structure. The conditional theorem would say "no axiom in $\mathcal{A}$ selects KO-dim signature". Sprint-scale, structurally analogous to F2's theorem.
- **E6 combination rule $K = \pi(B + F - \Delta)$**: 12 mechanisms eliminated for single-principle derivation; the structural theorem would be "no single morphism in $\mathcal{A}$ generates $K$ as a combination of the three independent spectral homes."
- **E7 / E8 LS-8a renormalization counterterms**: structural theorem would be "the spectral action axiom in $\mathcal{A}$ is single-cutoff; no morphism in $\mathcal{A}$ produces multi-cutoff renormalization counterterms."
- **G1–G6 multi-focal composition wall**: structural theorem would formalise the multi-focal-composition-wall pattern as a structural impossibility within $\mathcal{A}$. The most ambitious of the candidates; would unify six independent empirical instances under one theorem.
- **Chemistry-side $\eta$-trivialization analog**: multi-month NCG-research target, structural property of the FrozenCore $Z_{\rm eff}(r)$ pipeline mirroring $\{\gamma_F, D_F\} = 0$.

Of these, F3 (inner KO-dim) is the most tractable next target; sprint-scale, mirrors F2.

## Decision gate

The C3 arc has produced one theorem in this sprint. The user can:
1. **Pause** C3 here. Three theorem-grade non-selection results in the corpus; the empirical ones documented as such. Path #1 vein C is exhausted at sprint scale.
2. **Continue** to F3 inner KO-dim (sprint-scale mirror of F2).
3. **Continue** to E6 combination rule structural impossibility (medium-hard).
4. **Continue** to G1–G6 multi-focal composition unification theorem (hardest).
5. **Pivot** to a different vein within path #1, or to path #2 (Brown-Kleinschmidt outreach).

PI choice.

## Files

- **`papers/group1_operator_algebras/paper_32_spectral_triple.tex`** — `thm:n_gen_non_selection` + `rem:n_gen_scope` added to §VIII, three-pass clean compile at 76 pages.
- **`debug/sprint_c3_n_gen_non_selection_memo.md`** — this memo.
- (Pending) Paper 57 §6.1 + catalogue entry F2 update with theorem-grade reference.
- (Pending) CLAUDE.md + CHANGELOG updates.

## Honest scope

- **Theorem-grade conditional on canonical CCM rep.** The proof is rigorous within the conditional. The unconditional version remains multi-year.
- **No new mathematics.** The Direction 2 packing-reach NO-GO and the Read 2 representation-level NO-GO supplied all the structural content. This sprint crystallises them into a Theorem block in the same way Sprint Forced-Count Synthesis (2026-06-03) crystallised H1 into `thm:forced_count`.
- **Verified at theorem level by reference, not bit-exactly.** The proof sketch is six bullets; the structural content of each bullet is established in cited memos. No new computational verification was needed.
- **Theorem stands paired with empirical confirmation.** Direction 2 NO-GO + Read 2 NO-GO provide the empirical / diagnostic baseline; the Theorem provides the formal statement.
