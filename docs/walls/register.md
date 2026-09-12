# Walls Register — refined negative results

> The living home of the `/walls` pass (`.claude/commands/walls.md`). CLAUDE.md §3 is the append-only *ledger* of negatives (never deleted/modified, §13.5); **this register is the refinement layer on top of it** — re-verification status, hard/soft classification, and clustering into shared-mechanism findings. It points back to §3 rows and `debug/`/`memory/` memos; it never rewrites them.
>
> **Vocabulary.** Status: STANDING / SOFTENED / BREACHED / MIS-SCOPED / SUPERSEDED. Kind: HARD (proven structural impossibility) / SOFT (engineering·convention·precision·basis-limited, revisit-candidate) / OPEN-LEANING (tested-negative, a named better attempt exists). Cluster maturity: CRYSTALLIZED / FORMING / SINGLETON.

**Last run:** 2026-08-25 (+ `/aha` v5.1.2 kill-banking; bootstrap + accuracy-axis + composition-wall 3-body refinement — CHEM-ACCURACY worked fully, incl. the γ-determination arc and the polyatomic irreducible-three-body result; other clusters named, not yet worked).

---

## Cluster index

| Cluster | Maturity | Members | One-line shared mechanism |
|:--------|:---------|:-------:|:--------------------------|
| **CHEM-ACCURACY** | CRYSTALLIZED | 9 | Chemical accuracy and native sparsity are **two compounding walls** that co-locate favorably only at one center. |
| QED-SELECTION-RULES | *forming (unworked)* | ~4 | Graph photon is scalar; vector-photon selection rules need (L,M_L)-carrying photon (calibration, not skeleton). |
| ALPHA-COMBINATION | *forming (unworked)* | ~13 | B, F, Δ have independent spectral homes; their *combination* into α⁻¹ has no single-mechanism derivation (12+ eliminated). |
| GRAVITY-SUBSTRATE | *forming (unworked)* | ~12 | Finite-lattice substrate artifacts (Möbius α>1, A-coefficient, S_BH cutoff) vanish only in the a→0 continuum limit. |
| RELATIVISTIC-TAPERING | *forming (unworked)* | ~3 | jj-coupled Coulomb conserves M_J/κ-branch as *sums*, not parities → no original-basis Z-string Z₂ survives. |
| PERIODS-TRANSCENDENCE | *forming (unworked)* | ~5 | Two-electron/3-centre transcendence is set by *density topology* (elliptic, irregular), not the operator; no factorization closes it. |
| NUCLEAR-CONTINUUM-RESPONSE | *singleton-ish (unworked)* | ~3 | Sum-over-states / IR-extrapolation fragile for continuum response; the error is the interaction (Minnesota), not the solver. |
| MOLECULAR-BINDING-W1e | *forming (unworked)* | ~10 | Second-row over-binding lives at the integral-**specification** step, not correlation (DMRG=FCI); cheap engineering exhausted. |

### Banked kills awaiting clustering (`/aha` v5.1.2, 2026-08-25)

Candidates killed for a *decisive, stated* reason during the `/aha` pass on the v5.1.0
polyatomic arc. Filed here per the skill's B4 step; not yet clustered. Memo
`debug/sprint_config_operator_is_metric_memo.md` §4.

| Killed candidate | Decisive reason | Kind |
|:--|:--|:--|
| Commutant of the center-projections = the tapering-symmetry budget (would explain the M-vS gauge Pauli negative structurally) | The M-vS negative was measured on **LiH — two centers**, where the commutant is dim 6, not 1; and Hopf Z₂ tapering **works on BeH₂** (verified) where that commutant *is* trivial. `commutant(A) ⊆ commutant(H)` bounds in the useless direction — different objects. | HARD |
| The configuration operator's conical intersections mark a metric-**conditioning** blow-up (would unify with Paper 60's wall) | `cond(G)` is monotone through `d*` (19.87→16.19→15.04); the crossing is mid-spectrum (≈0.714 vs `λ_min`≈0.129). Conditioning and holonomy are different functionals of the same metric. Pinned as a test. | HARD |
| The multi-focal composition wall's threshold is 2→3 focal lengths | It fires at **two** (Paper 57 §614: `c_ab ≠ c_a·c_b`; every `multi_focal_wall_pattern` instance is 2-focal). Superseded by the two-rung scoping. | HARD |
| Chemistry error tracks three-body irreducible content (‖[[P,P],P]‖) rather than `max_n` | LiH has 5.3% `R_eq` error with **no three-body term at all** (2 centers), and CHEM-ACCURACY localizes the defect at 100% `max_n`, angular-invariant. Angular-operator axis vs radial-basis axis. | HARD |
| Paper 2's `K = π(B+F−Δ)` is a tame/wild (2-reducible, 3-wild) statement | Paper 2 §sec:rigidity_connection already states the 2+1 *sector* split (B, F on scalar S³; Δ on spinor); and it is an axis error regardless — B, F, Δ are numbers, not projections, so nothing generates an algebra. | HARD |
| The TC three-body non-collapse is the same obstruction as ≥3-projection wildness | Two independent kills. (i) The `/aha` reading misread the TC table: rank 4→9→16 runs over **l=1 / l=2(L_corr=2) / l=2(L_corr=3)** — growth in angular momentum, not leg count — and 87.5% of pairs already have rank ≥2 at l=1, so the TC obstruction fires at **two** legs. (ii) **Peter–Weyl**: SO(3)/SU(2) coupling is compact-group rep theory ⇒ type I ⇒ **tame** however non-abelian. Measured: SU(2) on l=0⊕1⊕2 has ‖[J,J]‖=2.0 yet commutant 3; three center-projections at ‖[P,P]‖=0.49 have commutant 1. Non-abelian and wild are independent. | HARD |
| WH7's compact modular-circle Z₂ is the T-symmetric equator of a Berry monopole, with Wick rotation as the U(1) lift | The falsifier had **already fired** in the corpus: a genuine temporal Wick involution `J_t = sign(D_t)` leaves the seminorm **bit-exactly `ω_q`**, unchanged from Euclidean (`tests/test_lorentzian_toeplitz_kplus.py`, 11 passed) — whereas the BeH₂ U(1) lift opens a gap ∝ φ. A lift that changes nothing is not a lift. **Z₂ pun**: Berry holonomy from a codimension-2 degeneracy vs the SU(2) double cover; `K = diag(2m_j)` has fixed integer spectrum, no parameter family, no degeneracy. Fifth member of the ℚ(i)-seam CONVERGENT family. | HARD |

> Only **CHEM-ACCURACY** is worked below. The rest are named from the §3 table so the next `/walls <cluster>` runs have targets; member counts are approximate first reads, not verified crystallizations.

---

## Cluster: CHEM-ACCURACY — why diatomics aren't chemically accurate

**Maturity: CRYSTALLIZED** · Owning papers: 13, 17, 19, 14 · Meta-findings: [[two_kinds_of_sparsity]], [[composition_wall_non_commuting_projections]], [[multi_focal_wall_pattern]]

### Framing diagnosis (not itself a wall): **exact ≠ accurate**
The v4.77–79 closed-form-ERI arc bought **decidability** (provable exact zeros → sparsity), **not accuracy**. Measured (CLAUDE.md orientation box, v4.73.0 A/B/C): the integral-exactness axis is worth ~0.003 Ha on H₂ and then flatlines; basis size (max_n) is worth ~0.050 Ha and rising; the chemistry defect is **100% max_n orbital-basis**, bit-invariant to angular and quadrature refinement. So the accuracy wall is a **Layer-2 (projection) basis limit**, and "make the integrals exact" polishes the wrong axis. This diagnosis reframes the whole cluster.

### The crystallized statement (the payoff)
GeoVac's chemical-accuracy limit is **two compounding walls, not one**:

- **Wall A — accuracy (the e-e cusp / Layer-2).** The electron–electron cusp is a sharp feature no product of one-electron functions can represent efficiently; smooth bases resolve it only slowly in max_n. GeoVac's Sturmians have the *nuclear* cusp built in but are **no better than Gaussians on the e-e cusp**. Wall A is **breachable in principle** by explicit-r₁₂ machinery (F12 / TC / geminals).
- **Wall B — sparsity/scaling (non-orthogonality).** GeoVac has **only symmetry sparsity**, which dies at two centers (SO(4)→axial→point group; N3b: cross-center zeros are m-rule-only, constant factor). Every attempt to keep sparsity while adding molecular coupling **relocates** a non-orthogonality cost — it never removes it (cost conservation / dual-basis theorem). Wall B is **HARD**.

**The compounding — the actual insight:** Wall A and Wall B co-locate *favorably* **only at one center**.
- **Atoms:** Wall A is cheaply breachable (xTC: 0 fill-in at every basis, v5.0.9) **while Wall B is absent** (single center → symmetry sparsity intact). → the atomic-xTC win.
- **Molecules:** breaching Wall A **re-triggers** Wall B — the explicit geminal is ~93% collinear with its reference pair, so κ(S)≈180–250, the *same* non-orthogonality wall (v5.0.7 R12-CI conditioning). Cheap-operator **or** enlarged-space accuracy, not both.

> **Therefore: chemical accuracy and native sparsity co-locate at one center and separate at two.** The framework's honest lanes are (i) **ATOMS**, where both walls are favorable, and (ii) **quantum-simulation STRUCTURE**, where sparsity/decidability is the *product* and accuracy is not the claim. This is why "cracking diatomic *accuracy*" is not the framework's game — and precisely where it *does* win.

### Operational consequence (dispatch rule)
- **(A) Atomic accuracy + sparsity** → **GO — now delivered across the whole monoatomic library (v5.0.11).** xTC sparsity swept over Z=1–56: fill-in is governed exactly by Unsöld's theorem (0 iff the reference density is spherical = closed / half-filled-high-spin / s-open, **29 of 56 atoms**); non-spherical p/d references give only a tiny structured fill-in (≤0.15% of blocks; open-d adds an L′=4 multipole that first bites at g). The atomic lane is mapped; the open frontier here is *accuracy* per-atom (the radial/1-norm side), not sparsity.
- **(B) Molecular chemical accuracy via a better basis/integral** → the **exact≠accurate trap**; STOP unless the proposal breaches Wall B — and the dual-basis theorem says it can't.
- **(C) Molecular quantum-simulation resource claims** → **GO.** Sparsity/decidability is the product; do not attach a chemical-accuracy claim.

### Falsifier for the crystallization
A molecular (≥2-center) method that reaches sub-mHa **and** keeps native integral sparsity **without** relocating the cost into a dense S⁻¹ᐟ²/S⁻¹ contraction or a device-measured many-body overlap. Would break both the two-walls compounding *and* [[two_kinds_of_sparsity]].

### Member walls

| Wall | Ledger ref | Status | Kind | Note |
|:-----|:-----------|:-------|:-----|:-----|
| e-e cusp = accuracy wall; chemistry defect 100% max_n | §2 v4.73.0 (A/B/C); §5 Level-3 note | **STANDING** | SOFT→HARD¹ | ¹SOFT-for-atoms (breached, xTC), HARD-for-molecules (re-triggers Wall B). |
| "cheap cusp treatment loses sparsity" | §2 v5.0.7 (TC tension) | **SOFTENED** | — | Breached for **atoms**: xTC handles the cusp AND keeps Gaunt sparsity (0 fill-in), Paper 14 §tc_atomic_sparsity, v5.0.9. Molecules still lose it at 2 centers. |
| Two kinds of sparsity (symmetry dies at 2 centers) | [[two_kinds_of_sparsity]]; §3 rows | **STANDING** | **HARD** | SO(4)→axial→point group; N3b m-rule-only, constant factor. The load-bearing Wall B. |
| Composition wall = non-commuting projections | [[composition_wall_non_commuting_projections]]; §2 v4.73.0, v5.1.0 | **STANDING** | **HARD** | ‖[P_A,P_B]‖=0.50; Löwdin=forced commutation=densification; dual-basis theorem (Artacho–del Bosch). **Refined + CLASSIFIED 2026-08-25:** LiH (2 centers) is *completely* pairwise (Halmos — angles 7.6°/44.7°/67.3°, REDUCIBLE, commutant dim 6). BeH₂ (3 centers): the three center-projections generate the FULL M₆ **irreducibly** (commutant dim **1**) — one indecomposable three-body block; nested commutator 0.32–0.42 at bonding, →0 far, survives charge asymmetry. **Qualitative reducible→irreducible jump** = the operator-algebra reason polyatomics are categorically harder (cannot be decomposed at all vs always ⊕≤2×2). **Geometric family (2026-08-25, sub-agent-verified):** F=ΣP_i has **exactly 3 conical intersections** in the linear region (central 2.445 + off-axis mirror pair, exact gaps ~1e-8) with a **π Berry phase**; bending enriches (curve d\*(θ), ≥3 branches, a structural Renner–Teller mechanism). **DISTINCT** from BeH₂'s physical electronic CI (bent Be+H₂ insertion) — the structural/generic cousin, NOT an identity (the shared π phase is generic). Algebra (M₆) + geometry (CI lattice) = same 3-center object two ways. Structural, not a chemistry lever. |
| Löwdin retrofit / non-orthogonal encoding | §3 (2026-07-07 ×2) | **STANDING** | **HARD** | Cost conservation: Löwdin→17.9× denser; biorthogonal→dual basis just as dense; NOCI→device-measured overlaps. Same wall in four currencies. |
| W1e — second-row over-binding | [[w1e_localized_hamiltonian_spec_level]]; §3 (F1–F6, DMRG, HF) | **STANDING** | **HARD** | Localized to the integral-**specification** step, not correlation (DMRG=FCI on NaH). Cheap engineering exhausted (F1–F6, Schmidt, core, kwargs, explicit-core HF). |
| TC cheap↔accurate don't co-locate | §2 v5.0.7; §3 (2026-08-23 ×3) | **STANDING** | **HARD** (conjunction) | The two-walls compounding *in the cusp currency*. |
| cheap-TC γ-selection criterion | §2 v5.0.8, v5.0.12; §3 (2026-08-23, 2026-08-24 ×2) | **STANDING** | **HARD** (cheap route)¹ | ¹Firmed 2026-08-24: no *determined* γ rescues the cheap non-Herm operator — global fixed/physics γ scatters (F12-std 5/2/59 mHa He/Li/Be), and a determined position-dependent γ(r)∝n^{1/3} makes it *worse* (slopes 80/134/348 vs 40/43/253). Mechanism: the non-Herm convective term feeds on the correlation factor's gradients. Determined γ works only variationally (R12-CI/F12 → κ(S)≈200). The fragility also *scales with electron count* (crossing slope 40→253, He→Be). |
| accurate per-atom accuracy (xTC recovery) | §2 v5.0.12 | **STANDING** | — | xTC at the *oracle* γ recovers 94–98% of the correlation gap incl. Be (4e) — the cusp recovery scales; only the *cheap self-selection* fails (above row). The determined γ(r)∝n^{1/3} exists (Giner μ(r) / Wagner–Gori-Giorgi avoidance radius) but lives on the variational side. |
| skeleton-native γ(r) from equal-area cells | §2 v5.0.12 (deferred) | **OPEN** | — | Does Paper 0's equal-area packing *give* γ(r)∝n^{1/3}? Not a scaling calc — Paper 0's "shells" are l-subshells (not real space), so it needs the full Fock momentum→position map. Genuine open skeleton-program derivation, not faked. |
| Elliptic/CM basis for chemical accuracy | §3 (2026-08-23 ×3) | **STANDING** | **HARD** | The bond's genus-1 elliptic period is *radial + inert*; the accuracy gap is the *e-e* cusp, which a radial basis principle can't reach. |

### This run's delta
Bootstrap. Net-new refinements over the raw §3 ledger:
1. **SOFTENED** the "cheap cusp loses sparsity" wall — the atomic-xTC breach (v5.0.9) was in the ledger as a *positive*, but the corresponding *negative* was still implicitly monolithic. Now scoped: breached for atoms, standing for molecules.
2. **CRYSTALLIZED** the two-walls-compounding statement, unifying the atomic-xTC win (v5.0.3–5.0.9) with the molecular TC tension (v5.0.7) via "Wall A is breachable exactly where Wall B is absent." Candidate for PI promotion (see below).


### Delta 2026-09-12 (PI-directed consolidation, v5.11.0)

**The composition wall was one row; it is three independent axes, and they now have different statuses.** Driver: the Paper 60 KMS/preconditioner arc (CHANGELOG v5.10.18-19), three literature scans, `tests/test_paper60_{kms_attribution,preconditioner}.py`.

The register carried `||[P_A,P_B]|| = 0.50` as a measurement. It is the *saturation value* of an exact formula, `max_k sigma_k sqrt(1 - sigma_k^2)`, over the **same** singular spectrum that carries `cond(S) = (1 + sigma_max)/(1 - sigma_max)`. So the commutator and the conditioning are one object, as Paper 60 already said. What was riding along with them, and should not have been, is the `l`-block-structure loss. Splitting the three:

| Axis | Status | Mechanism | Changed? |
|:--|:--|:--|:--|
| **Conditioning** (`cond(S) ~ n^2`) | **BREACHED** | The symbol's zero at `chi = pi` has known order and location, so a band-Toeplitz preconditioner in the sense of Serra (*Math. Comp.* **66**, 651 (1997)) removes it: `P = tri(1,2,1)` exactly, DST-I diagonalizable in closed form, `cond -> 2.23` **flat in n** against 19 127 at `n = 160`. Whitening-invariant, so the spectrum is unchanged. | **YES** -- this axis was thought hard and is not |
| **Locality** (`S^{-1/2}` dense) | **STANDING** (capped, not removable) | The *other* pole. Preconditioning cures `chi = pi` and cannot touch the `chi -> 0` chirp, which fixes the off-diagonal envelope at `j^{-5/4}`. Measured: `G^{-1/2}` bandwidth 11->17 at 1e-2 (vs a fixed 0.72`n` for `S^{-1/2}`) but 27->141 at 1e-3. Profile exponent n-STABLE for `G^{-1/2}` (-1.19), DRIFTING for `S^{-1/2}` (-0.90 -> -0.78). | scoped |
| **`l`-block structure** | **STANDING, HARD** | **Proposition D**: if `S` is not `l`-block diagonal, no block-diagonal congruence orthogonalizes it -- Loewdin, canonical or Cholesky. Holds at **every** `cond(S) > 1` and does not relax as `cond(S) -> 1+`. Not a functional of the sigma spectrum at all. | **STRENGTHENED** |

**Consequence for the cluster's dispatch rule.** Rule (B) ("molecular chemical accuracy via a better basis/integral -> STOP unless the proposal breaches Wall B") is unchanged in outcome but its *reason* is now sharper: a proposal that improves conditioning no longer counts as progress toward sparsity, because Proposition D makes those independent. Conversely a conditioning-only proposal should no longer be rejected on Wall-B grounds -- that axis is open, and the quantum-resource lane (rule C) is where it pays.

**Mechanism correction (2026-09-12, v5.11.2).** The conditioning axis's obstruction is `||sigma||_inf = 1` attained at `p = 0` (Ron-Shen), **not** completeness of the one-centre set. The Bessel deficit against the one-centre span plateaus at 0.380 / 0.696 / 0.907 for `kR = 1/2/4`, flat over `N = 16..256`, so that set is far from complete in the molecular metric. The corrected reading is stronger operationally: the near-dependence is **one direction**, which is why a fixed rank-`M-1` rotation removes it. Two escapes checked and closed: adding `l>0` cannot help (`sigma_max` is non-decreasing under submatrix extension), and excluding `l=0` does not either (an `l>=1` combination can still concentrate near `p=0`) -- the second is reasoning, not measurement. Backing `tests/test_paper60_one_direction.py`.

**Falsifier for the split.** A congruence that is simultaneously (i) `l`-block diagonal and (ii) orthogonalizing, on a metric with nonzero inter-center coupling -- which Proposition D forbids outright; or a locality repair that reaches the `chi -> 0` chirp, which would have to change the *symbol*, not the matrix.

**Scope, now measured.** The breach reaches **water's `A_1` block** -- the symmetry-inequivalent-center case where the gerade lever fails: raw `cond ~ N^1.96` (independently reproducing the paper's `N^1.97`) against a bounded `38.45 -> 44.06` over `N = 12..192`. It works because the degeneracy's DIRECTION is geometry-independent: at `chi = pi` every block symbol tends to `j0(0) = 1`, so for `M` centers the matrix symbol is the rank-one all-ones matrix and its null space has dimension `M-1`, fixed. Control: the naive `blockdiag(P,P)` without the null-direction rotation leaves the growth intact (`2766 -> 42008`), so the rotation is doing the work. Remaining scope: `s`-sector shared-scale bases at `M = 2, 3`; the end-to-end resource claim still needs the sine-transform circuit and a block-encoding of `G` priced.

---

## Promotion candidates (PI-gated)

**CHEM-ACCURACY two-walls compounding.** One-sentence mechanism: *chemical accuracy (Wall A, breachable) and native sparsity (Wall B, hard) co-locate favorably only at one center, so breaching the cusp for molecules re-triggers the non-orthogonality wall.* Support: ≥3 independent negatives (R12-CI conditioning v5.0.7; flexible-Jastrow cap v5.0.8; elliptic-basis null v5.0.1) + 1 breach (atomic xTC v5.0.9) + 2 crystallized meta-findings it subsumes. Falsifier: above. **Proposed home:** a sharpening of [[two_kinds_of_sparsity]] and/or a Paper 14 remark tying the atomic xTC win to the molecular tension. **PI decision required** to write it into a paper.

---

## Backlog — clusters to work next
Run `/walls <cluster>` on any of: QED-SELECTION-RULES · ALPHA-COMBINATION · GRAVITY-SUBSTRATE · RELATIVISTIC-TAPERING · PERIODS-TRANSCENDENCE · NUCLEAR-CONTINUUM-RESPONSE · MOLECULAR-BINDING-W1e. Each needs its members pulled from §3, re-verified against current state, tagged hard/soft, and a crystallization attempted.
