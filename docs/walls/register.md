# Walls Register — refined negative results

> The living home of the `/walls` pass (`.claude/commands/walls.md`). CLAUDE.md §3 is the append-only *ledger* of negatives (never deleted/modified, §13.5); **this register is the refinement layer on top of it** — re-verification status, hard/soft classification, and clustering into shared-mechanism findings. It points back to §3 rows and `debug/`/`memory/` memos; it never rewrites them.
>
> **Vocabulary.** Status: STANDING / SOFTENED / BREACHED / MIS-SCOPED / SUPERSEDED. Kind: HARD (proven structural impossibility) / SOFT (engineering·convention·precision·basis-limited, revisit-candidate) / OPEN-LEANING (tested-negative, a named better attempt exists). Cluster maturity: CRYSTALLIZED / FORMING / SINGLETON.

**Last run:** 2026-08-25 (+ `/aha` v5.1.2 kill-banking; bootstrap + accuracy-axis + composition-wall 3-body refinement — CHEM-ACCURACY worked fully, incl. the γ-determination arc and the polyatomic irreducible-three-body result; other clusters named, not yet worked). · **2026-09-19 (v5.14.6):** architecture-swap audit **composed → prolate-native diatomic**, PI-adjudicated — a new *kind* of pass (audit against an architecture change, not a shared-mechanism cluster); see the dedicated section below. · **2026-10-16 (v5.20.18):** `/walls PERIODS-TRANSCENDENCE`, PI-invoked — the cluster worked FORMING → CRYSTALLIZED (home: Paper 59 `sec:f12` + `sec:beyond`, seam J9); one member MIS-SCOPED (route-specific), one new OPEN-LEANING member filed; see its section below. · **2026-10-27 (v5.21.19):** `/walls PERIODS-TRANSCENDENCE`, PI-invoked — the crystallization's "never by the operator" clause MIS-SCOPED by four paper-owned results (Gaussian and Slater geminal kernels move the class); re-crystallized as "the integral's graph: density topology plus kernel-line masses"; P3e-2 kind settled to SOFT (missing R_V periods); see the cluster's 2026-10-27 delta.

---

## Cluster index

| Cluster | Maturity | Members | One-line shared mechanism |
|:--------|:---------|:-------:|:--------------------------|
| **CHEM-ACCURACY** | CRYSTALLIZED | 9 | Chemical accuracy and native sparsity are **two compounding walls** that co-locate favorably only at one center. |
| QED-SELECTION-RULES | *forming (unworked)* | ~4 | Graph photon is scalar; vector-photon selection rules need (L,M_L)-carrying photon (calibration, not skeleton). |
| ALPHA-COMBINATION | *forming (unworked)* | ~13 | B, F, Δ have independent spectral homes; their *combination* into α⁻¹ has no single-mechanism derivation (12+ eliminated). |
| GRAVITY-SUBSTRATE | *forming (unworked)* | ~12 | Finite-lattice substrate artifacts (Möbius α>1, A-coefficient, S_BH cutoff) vanish only in the a→0 continuum limit. |
| RELATIVISTIC-TAPERING | *forming (unworked)* | ~3 | jj-coupled Coulomb conserves M_J/κ-branch as *sums*, not parities → no original-basis Z-string Z₂ survives. |
| **PERIODS-TRANSCENDENCE** | CRYSTALLIZED (2026-10-16; re-scoped 2026-10-27) | 7 | For two-electron integrals over Slater densities the transcendence class is set by the **integral's graph: the density topology the kernels couple, plus the masses on the kernel lines**. A massless Coulomb kernel adds no threshold, so for Coulomb integrals the density topology alone sets the rung; a massive (geminal) kernel can move it. Never set by the basis, the number of nuclei, or a factorization; the three-centre radial core has no closed form through any factorization. A pattern across the tested instances. Home: Paper 59 `sec:f12` + `sec:beyond` + `sec:density_ring`; seam J9. |
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
| LiH R_eq drift = cross-V_ne R-slope = bond-orbital finite extent = **frozen-core rigidity** | §2 v5.14.8, v5.15.2, v5.15.6; `test_paper19_reqdrift_mechanism.py`, `test_prolate_lih_r12_coupled.py` | **STANDING** | **HARD** (parameter-free)¹ | Two-term V_NN↔cross-V_ne balance; the 8.8% drift is the *finite-extent screening* of the fixed Z_orb=1 bond orbital. Mechanism PINNED (v5.15.2): cross-attraction R-slope → +Z/R² (point-charge limit) as the orbital contracts, and that limit *is* the required force Z_A Z_B/R²; the fixed orbital delivers 93.9%. ¹Contraction is **inexpressible** in the parameter-free construction — h1 = −Z_orb²/2n² ties Z_orb to Z_nuc, so forcing Z_orb>1 drives E below exact (−9.6 vs −8.07). **LOCALIZED to the core (v5.15.6):** the prolate-native route with a *frozen* core (variational valence + r₁₂) was built and tested and STILL drifts outward, growing with angular basis (+3.4%→+5.5%, l=2→3), vs HeH⁺'s convergence (−4.5%→−0.5%) at the identical basis — so a variational *valence* is insufficient; the rigidity is the frozen core's static screening. Breach now requires a variational *core* = all-electron (N=4) prolate CI. The drift is absent in H2/He/HeH⁺ precisely because those have no frozen core. |

### This run's delta
Bootstrap. Net-new refinements over the raw §3 ledger:
1. **SOFTENED** the "cheap cusp loses sparsity" wall — the atomic-xTC breach (v5.0.9) was in the ledger as a *positive*, but the corresponding *negative* was still implicitly monolithic. Now scoped: breached for atoms, standing for molecules.
2. **CRYSTALLIZED** the two-walls-compounding statement, unifying the atomic-xTC win (v5.0.3–5.0.9) with the molecular TC tension (v5.0.7) via "Wall A is breachable exactly where Wall B is absent." Candidate for PI promotion (see below).


### Delta 2026-09-12 (PI-directed consolidation, v5.11.0)

**The composition wall was one row; it is three independent axes, and they now have different statuses.** Driver: the Paper 60 KMS/preconditioner arc (CHANGELOG v5.10.18-19), three literature scans, `tests/test_paper60_{kms_attribution,preconditioner}.py`.

The register carried `||[P_A,P_B]|| = 0.50` as a measurement. It is the *saturation value* of an exact formula, `max_k sigma_k sqrt(1 - sigma_k^2)`, over the **same** singular spectrum that carries `cond(S) = (1 + sigma_max)/(1 - sigma_max)`. So the commutator and the conditioning are one object, as Paper 60 already said. What was riding along with them, and should not have been, is the `l`-block-structure loss. Splitting the three:

| Axis | Status | Mechanism | Changed? |
|:--|:--|:--|:--|
| **Conditioning** (`cond(S) ~ n^2` on the `s`-sector) | **BREACHED on the `s`-sector; STANDING on a sigma basis for preconditioning** (relabelled twice 2026-10-09, PI-approved; second relabel: no preconditioner of the s-sector kind exists, canonical orthogonalisation is a partial remedy -- see the sigma-basis note below. First relabel: with `l > 0` cond is exponential in shells, gerade too, near-null directions multiply -- see the scope correction below) | On the `s`-sector: the symbol's zero at `chi = pi` has known order and location, so a band-Toeplitz preconditioner in the sense of Serra (*Math. Comp.* **66**, 651 (1997)) removes it: `P = tri(1,2,1)` exactly, DST-I diagonalizable in closed form, `cond -> 2.23` **flat in n** against 19 127 at `n = 160`. Whitening-invariant, so the spectrum is unchanged. | **YES** -- this axis was thought hard and is not |
| **Locality** (`S^{-1/2}` dense) | **STANDING** (capped, not removable) | The *other* pole. Preconditioning cures `chi = pi` and cannot touch the `chi -> 0` chirp, which fixes the off-diagonal envelope at `j^{-5/4}`. Measured: `G^{-1/2}` bandwidth 11->17 at 1e-2 (vs a fixed 0.72`n` for `S^{-1/2}`) but 27->141 at 1e-3. Profile exponent n-STABLE for `G^{-1/2}` (-1.19), DRIFTING for `S^{-1/2}` (-0.90 -> -0.78). | scoped |
| **`l`-block structure** | **STANDING, HARD** | **The block-diagonal congruence result** (re-attributed 2026-09-12: Loewdin/Slater-Koster 1954; the label "Proposition D" is retired [retracted 2026-09-12: p60-prop-d-as-new] and only the `l`-vs-`m` application is ours): if `S` is not `l`-block diagonal, no block-diagonal congruence orthogonalizes it -- Loewdin, canonical or Cholesky. Holds at **every** `cond(S) > 1` and does not relax as `cond(S) -> 1+`. Not a functional of the sigma spectrum at all. | **STRENGTHENED** |

**Consequence for the cluster's dispatch rule.** Rule (B) ("molecular chemical accuracy via a better basis/integral -> STOP unless the proposal breaches Wall B") is unchanged in outcome but its *reason* is now sharper: a proposal that improves conditioning no longer counts as progress toward sparsity, because that result makes those independent. Conversely a conditioning-only proposal should no longer be rejected on Wall-B grounds -- that axis is open, and the quantum-resource lane (rule C) is where it pays.

**Mechanism correction (2026-09-12, v5.11.2).** The conditioning axis's obstruction is `||sigma||_inf = 1` attained at `p = 0` (Ron-Shen), **not** completeness of the one-centre set. The Bessel deficit against the one-centre span plateaus at 0.380 / 0.696 / 0.907 for `kR = 1/2/4`, flat over `N = 16..256`, so that set is far from complete in the molecular metric [retracted 2026-10-09: p60-incomplete-in-molecular-metric] (s-sector only -- see the 2026-10-09 scope correction below). The corrected reading is stronger operationally: the near-dependence is **one direction**, which is why a fixed rank-`M-1` rotation removes it. Two escapes checked and closed: adding `l>0` cannot help (`sigma_max` is non-decreasing under submatrix extension), and excluding `l=0` does not either (an `l>=1` combination can still concentrate near `p=0`) -- the second is reasoning, not measurement. Backing `tests/test_paper60_one_direction.py`.

**Scope correction (2026-10-09; row relabelled the same day, PI-approved).** Every measurement behind the BREACHED row is on the **`s`-sector** sine basis. On a full sigma basis (`n <= N`, `l < n`, `m = 0`; `debug/sturmian_sigma_basis_sw.py`, position-space prolate quadrature, reproduces the `s`-sector plateau to 3 digits) the Bessel deficit does **not** plateau (0.158 -> 0.037 over `N = 2..6` at `kR = 1`), so the plateau above is the `l = 0` projection of an off-centre function, not incompleteness of the Coulomb-Sturmian set in the molecular metric. With `l > 0`: `cond(S)` is exponential in the number of shells (18 -> 1.1e11 over 3 -> 15 functions/centre at `kR = 3`), the gerade sector too (2.3 -> 2.1e7), and the near-null directions multiply (0,1,3,5,8 eigenvalues < 1e-2 for `N = 2..6`). So "one direction", the tri(1,2,1) preconditioner and the flat gerade constant are `s`-sector results, and the conditioning axis is **open** beyond the `s`-sector -- not breached (row status updated accordingly). The "adding `l>0` cannot help" sentence above was right about `sigma_max`; what it missed is that `l>0` is where the accuracy is (H2+ 5e-6 Ha at 15 functions/centre vs 1.3% `s`-only). Paper 60 sec:molecular carries the measurement.

**Sigma-basis verdict (2026-10-09, PI-approved: OPEN -> STANDING for preconditioning).** In the own-centre SW metric `<f|(p^2+k^2)/2k^2|g>` each centre's axial Sturmian set is complete and orthonormal and the other centre's is its translate, unitary in that metric (translation reading = Shibuya-Wulfman 1965), so the two-centre metric is a finite section of `[[I,U],[U*,I]]` (spectrum {0,2}): half the directions are redundancy, not a symbol zero. Any centre-block-diagonal SPD preconditioner gives `cond(P^-1 S) >= cond(S)`; per-l, per-shell and tri(1,2,1) levers gain constant factors; fixed-rank deflation cannot bound it. Partial remedy: canonical orthogonalisation (Avery's phase-out) at 1e-6 holds cond at 1e6 with < 1e-7 Ha H2+ loss through 36 functions/centre, but needs a classical eigendecomposition and a threshold that falls with the basis (thresholds 1e-2..1e-4 stall). Open follow-on: whether a shrinking-threshold schedule converges cleanly. Backing `tests/test_paper60_sigma_redundancy.py`; Paper 60 sec:molecular.

**Falsifier for the split.** A congruence that is simultaneously (i) `l`-block diagonal and (ii) orthogonalizing, on a metric with nonzero inter-center coupling -- which that result forbids outright; or a locality repair that reaches the `chi -> 0` chirp, which would have to change the *symbol*, not the matrix.

**Scope, now measured.** The breach reaches **water's `A_1` block** -- the symmetry-inequivalent-center case where the gerade lever fails: raw `cond ~ N^1.96` (independently reproducing the paper's `N^1.97`) against a bounded `38.45 -> 44.06` over `N = 12..192`. It works because the degeneracy's DIRECTION is geometry-independent: at `chi = pi` every block symbol tends to `j0(0) = 1`, so for `M` centers the matrix symbol is the rank-one all-ones matrix and its null space has dimension `M-1`, fixed. Control: the *discriminating* test is the SELECTIVE `blockdiag(P,I)` applied in the UNROTATED frame -- it grows as `N^3.79` and reaches `4.4e6` at `N=192`, **106x WORSE than untreated**, while the same band aligned to the null direction is flat at 44; so the ALIGNMENT, not the banding, is doing the work. (The uniform `blockdiag(P,P)` is NOT a valid control here: it commutes with the rotation to `3e-13` and so cannot discriminate it -- it runs `N^1.94`, essentially the raw `N^1.96`. Corrected 2026-09-13, DELTA #3.) Remaining scope: `s`-sector shared-scale bases at `M = 2` and **non-collinear** `M = 3`; the **collinear case is open and is NOT claimed** (corrected 2026-09-12: the unqualified "M = 2, 3" authorised dispatch into exactly the regime Paper 60 declines, and BeH2 and CO2 are collinear and in this corpus's own library). Mechanism: for collinear centres `P D2 P` is rank ONE, so only one of the `M-1` null directions opens at order `p^2` and the rest at 4, 6, ..., 2(M-1). **Resource claim closed 2026-09-12 (v5.11.3):** `G` has a DIRECT block-encoding — its symbol is a bounded ratio, `||ratio||_inf = 0.3716 = ||G||`, so a circulant-embedded Toeplitz-minus-Hankel encoding carries `alpha = O(1)` instead of the composed `O(n^2)`; the resulting whitening reaches the amplitude floor to 0.13% with residual conditioning 1.234. The metric penalty scales as `n` against the untreated `n^3`. What remains uncompiled is the circuit, not the construction.

---

### Delta 2026-09-20 (v5.15.2) — LiH R_eq-drift mechanism PINNED

v5.14.8 localized the balanced-LiH bond-length drift to one term (the cross-center
V_ne R-slope, 8.8% too weak) and *named* the mechanism. This delta *demonstrates*
it, upgrading a named mechanism to a proven WALL with the mechanism pinned (the
mission's WALL deliverable, §1.7). Memo `debug/sprint_balanced_reqdrift_mechanism_memo.md`.

- **The deficit is finite-extent screening, and its limit is the required force
  exactly.** Closed-form: the cross-attraction R-slope of a 1s orbital of exponent
  z is $dV/dR = d/dR\{-(Z/R)[1-(1+zR)e^{-2zR}]\}$, which rises monotonically to the
  point-charge limit $+Z/R^2$ as $z\to\infty$ — and $+Z/R^2$ **is** the required
  force $Z_A Z_B/R^2$ that cancels $dV_{NN}/dR$. The fixed Z_orb=1 orbital gives
  93.9% (6.1 pp of the measured 8.8 pp from one 1s term; the rest is 2s/2p
  admixture + the weaker Li→H cross term). Backing `test_paper19_reqdrift_mechanism.py`.
- **The wall is STRENGTHENED, not softened.** The obvious repair — relax the bond
  exponent — is *structurally inexpressible* in the parameter-free construction:
  h1 = −Z_orb²/2n² locks Z_orb to Z_nuc, so forcing Z_orb=2.3 on the Z=1 H centre
  deepens h1 (trace −11.0→−19.6 Ha) and drives E below the exact floor (−9.6 vs
  −8.07) — a variational catastrophe, not physics. Freeing the orbital needs h1
  from true kinetic+nuclear integrals = a non-parameter-free (variational) basis.
- **Corpus cross-check (the control group).** H2/He show no such drift *because*
  they use variational natural-geometry bases (prolate α_opt for H2; hyperspherical
  for He) where the orbital contracts freely. Confirms the drift is a recipe
  artifact, not fundamental.
- **Falsifier.** A parameter-free (zero-knob) balanced/composed diatomic that
  reaches sub-% R_eq without a variational, contraction-capable orbital would break
  the "finite-extent, inexpressible" reading.

### Delta 2026-09-20 (v5.15.6) — the prolate frozen-core route BUILT + TESTED: rigidity localizes to the core

The v5.15.2 falsifier flagged the "prolate-native ≥4e, frozen-core route" as UNBUILT.
The **frozen-core** half is now built and tested (`debug/lih_r12_coupled.py`,
`lih_r12_req_scan.py`; CHANGELOG v5.15.6), and it SHARPENS the wall rather than
breaching it:

- **The prolate route with a frozen core still drifts.** A variational, r₁₂-correlated
  prolate valence CI (frozen Li 1s² via closed-form Hartree V_H + a λ-plateau-verified
  Huzinaga projector, both now r₁₂-coupled) gives R_eq **+3.4% (l=2) → +5.5% (l=3)** —
  still outward, still growing with the angular basis. r₁₂ absorbs part of the angular
  over-polarization (+7.1%→+5.5% at l=3) but cannot reach a core d.o.f.
- **Controlled comparison (decisive).** The identical r₁₂-CI on HeH⁺ (two active
  electrons, NO frozen core) *converges*, −4.5%→−0.5%, across the same l=2→3. The one
  variable differing between drift and convergence is the frozen core.
- **Verdict:** the finite-extent rigidity lives in the **frozen core's static
  screening**, not only the valence bond orbital → a variational *valence* is
  insufficient. The breach path narrows: **prolate-native → prolate-native
  variational-CORE (all-electron N=4)**, which meets the N=4 explicit-r₁₂ wall.
- **Falsifier (updated).** An all-electron (variational-core) prolate ≥4e two-center CI
  that reaches sub-% R_eq would breach; if it *still* drifts for the same reason, the
  obstruction is deeper than the frozen core. This is the remaining unbuilt artifact.
- Honest baseline: prolate frozen-core (+5.5% @ l=3) improves on the *balanced*
  recipe (+8.8%, Paper 19) but does NOT beat *composed* l-dep-PK (5.3% @ l_max=2,
  Paper 17) — a mechanism diagnostic, not a new best LiH R_eq.

### Delta 2026-09-20 (v5.15.7) — all-electron variational-core prolate LiH: cure CONFIRMED (bracketed); clean R_eq is a NUMERICAL wall, NOT physics

The v5.15.6 falsifier called for an all-electron (variational-core) prolate ≥4e CI.
BUILT + validated (`debug/prolate_allelectron_fci.py`: multi-exponent prolate FCI +
generalized-m π ERIs; V_ee pairwise → avoids the N=4 explicit-r₁₂ 4-body wall; HF==eckart,
μ-kernels==numerical, H₂ converges). Two PI-directed review cycles then characterized the
outcome:
- **Cure CONFIRMED by bracketing.** The all-electron variational core removes the
  frozen-core rigidity — and R_eq changes SIGN, collapsing INWARD (mirror of the +5.5%
  outward frozen-core drift). The two errors bracket the truth (inward ← 3.015 → +5.5%
  frozen → +22% over-corrected), and standard all-electron QC gives good LiH R_eq. So
  "variational core = cure" STANDS. **HARD → the physics obstruction is retired**; what
  remains is numerical.
- **NEW numerical wall, mechanism pinned (two halves):** (i) the isolated Li²⁺ 1s² core
  (physically R-independent) swings **0.36 Ha** with R because the prolate grid's
  real-space resolution near the tight core is tied to R via r_A=(R/2)(ξ+η) — a GRID
  artifact (NOT BSSE: ghost-H basis makes zero difference; confirmed by refinement
  shrinking the spread), swamping the 0.088 Ha bond → inward collapse; (ii) a co-dominant
  **two-center-basis-completeness** half (Coulomb-Sturmian-like orbitals at fixed ζ span
  the true 1s² with R-dependent fidelity — no grid fixes it). Kind: SOFT (numerical,
  solvable) — a graded grid helped only ~15% (0.081→0.069 Ha over [2.6,3.75]); the real
  fix is a dedicated atom-centered/analytic core, a scoped future build.
- **Engine hardening (keeper):** fixed `_attraction_matrix` to the analytic
  (Z_A/r_A+Z_B/r_B)·J = (R/2)²[Z_A(ξ−η)+Z_B(ξ+η)] (the 1/r focus singularity cancels the
  Jacobian) — grid-robust; identical on the stock grid, so all prior results stand.
- **Falsifier for the numerical wall:** a dedicated atom-centered (or exact-integral)
  tight-core treatment in the prolate valence CI that reaches sub-% R_eq would breach it;
  a proof that the two-center basis cannot represent a tight core R-consistently at any
  cost would harden it. Guards `tests/test_prolate_allelectron_fci.py`; full chronicle
  `debug/track_logs/prolate_native_lih.md`; CHANGELOG v5.15.7.

### Delta 2026-09-21 (v5.15.8–9) — the "N=4 explicit-r₁₂ 4-body wall" is SOFTENED (scalar Coulomb term), with a Be energy behind it

**What was implicitly HARD.** This register (and CLAUDE.md §2, `memory/r12_generalization_boundary_n3_n4.md`) referred throughout to "the N=4 explicit-r₁₂ 4-body wall" as the boundary the exact-algebraic James-Coolidge program (r₁₂ in the basis, NO strong-orthogonality projector ⇒ NO resolution-of-identity) hits past N=3 — treated as the point where exact-no-RI *stops*, an RI-class obstruction like Gaussian-F12. That strong reading was never tested. It is now **REFUTED for the scalar Coulomb term** (PI-directed diagnostic + energy build).

- **Status → SOFTENED · Kind → SOFT (cost/encoding, not RI/decidability).** The only genuinely-4-body-connected term in `⟨Φ|F H F|Φ⟩` (F=Σf_pq multiplicative) is the scalar Coulomb **chain** `f₁₂·f₃₄·(1/r₁₃)` — kinetic gradients are pair-local (`∇_k f₁₂·∇_k f₃₄ = 0` for disjoint pairs) and V_ne is one-body, so only the two-body Coulomb can bridge two disjoint correlation edges. That term is **exact, RI-free, terminating** (the bridging-Coulomb multipole sum truncates at `L ≤ 2·l_bridge` for s-type leaves — measured s→0, p→2, d→4 — and in general at `L ≤ 2(l_bridge + l_leaf)`: for p-leaves the L=4 channel the s-leaf bound discards is nonzero at 24σ, audit 2026-09-22; still finite ⇒ still RI-free) and **reducible** (factorizes across the bridge via the Legendre addition theorem into a finite (L,M) contraction of 3-body vertex kernels — the chain analogue of the N=3 triangle contraction).
- **Confirmed with numbers:** reduced == brute-force 12-D MC at **0.02%** (2p bridge / 1s leaves, exercises L>0) and, with the antisymmetric Be determinant (exchange chains, **L=0 all-s leaves only**), rel 8.6e-5 (−0.2σ) at 4e7 samples on a 2026-09-22 re-run — the **0.25%** originally logged here was not reproduced and Paper 12 quotes 2e-4; `debug/r12ci_4e_be_integral.py`, `be_r12ci_4body_exchange.py`). And demonstrated **inside a Be R12-CI energy**: E_R12 = −14.5572 Ha, corr −18.2 mHa (19% of Be's 94 mHa), variational, with the 4-body content in `gV = E[(F−F̄)²V_ee]` (`debug/be_r12ci_full.py`; the ill-conditioned coupling h and σ² done analytically-exact).
- **What stays HARD (separate axes — do NOT read this softening as global):** (i) **quantum encoding** — an explicit 4-body correlation operator is a 4-body Pauli string (§3 TC 2nd-quant plateau 3.4%, TC angular-gradient 2.66× Pauli); real wall for the quantum-simulation product. (ii) **non-Hermitian TC** (vector legs, non-abelian Gaunt/6j; §3 2026-08-23) — a different, harder operator. (iii) **kinetic-vector 4-body** — SETTLED (corrected 2026-09-22; this line said UNTESTED after v5.15.13–14 had closed it): `∇_k f₁₂·∇_k f₃₄ = 0` is trivial for disjoint pairs, and the IBP identity `gT2+gT3 = −½⟨χ Σ∇²F⟩` (`geovac/lih_r12ci/gT.py:7-11`) is exact for any real Φ₀; validated g_T to 0.56% vs VMC (`tests/test_lih_r12ci.py`). Not a wall.
- **On this narrow axis GeoVac is ahead of Gaussian-F12** (which uses RI for exactly these 3-/4-electron integrals). The general resolution was recorded here as "the bridge addition theorem, likely extending to N≥5"; **re-tiered 2026-09-22 (audit): a VALIDATED REDUCTION, not a theorem** — no written derivation, no stated hypotheses; N≥5 is conjectural. At the linear-in-F ansatz nothing is ever >4-body-connected for any N; what N≥5 or a multi-determinant reference changes is the DENSITY (non-isotropic leaves, non-block-factorized reference). At audit time (2026-09-22) that path was checked only by a 5%-tolerance MC-vs-MC test (`debug/r12ci_4e_wall_q2_reducibility.py:117`); the marriage build's Phase 1 (G-leaf, 2026-09-23) then ran it directly and it **PASSED** on two real non-isotropic densities of the 74-det CI reference: signed p-like core×bond leaf (6-D route, |reduced−MC|=1.07e-8 vs tol 1.73e-8) and directional bond leaf (12-D, 1.8e-4). Both dressed leaves carry channels beyond the 2·l_bridge bound (l=4 at ~1e-3), confirming the general 2(l_bridge+l_leaf) bound is operative. So the generality path is validated-on-real-densities at the integral level (the assembled multi-det energy is Phase 2/3). Guard `tests/test_lih_r12ci_bridge.py`; memo `debug/sprint_lih_marriage_memo.md`; plan `debug/lih_marriage_build_plan.md`.
- **Falsifier / breach direction:** a proof that the kinetic-vector 4-body term (or the N≥5 case) forces RI would re-harden a sub-part; a full-analytic H₁₁ Be energy (making g exact) would firm the SOFT verdict quantitatively. Kind is SOFT because the residual is polynomial cost + a native 4-body angular kernel, not an RI/decidability obstruction. Memo `debug/sprint_n4_wall_diagnostic_memo.md` (§7b = energy); CHANGELOG v5.15.8–9.

## Architecture-swap audit: composed → prolate-native diatomic (2026-09-19, v5.14.6; PI-adjudicated)

**What this is.** Not a shared-mechanism cluster — an audit of the §3 ledger against a
specific architecture change (the corpus's move away from composed toward prolate-native
diatomics, Paper 12 v5.13+). It fills, retroactively, the `rests-on:` edge the
negative-side record has never had: which architectural component each wall's *proof*
used. Fuller spec: `debug/sprint_walls_composed_prolate_scope_memo.md` (transient; the
load-bearing rules are restated here so this record stands alone).

**Primary gate (PI caution).** The audit bites only where an architecture *choice*
exists. Atom/diatomic → prolate-native / NOCI / exact-integral are real alternatives →
in scope. **Triatomic+ (3+ centers): composed is the SOLE architecture (prolate has two
foci); its walls there are permanent constraints, not stale artifacts → OUT OF SCOPE,
STAND by default.** Verdicts are per (wall × system-class): the same PK wall is
MIS-SCOPED for diatomics and STANDING for triatomics.

**Integrity rule (non-negotiable).** MIS-SCOPED = "proven on machinery the prolate-native
target drops, so it no longer *informs* the new program." NOT "the approach works."
The entire MIS-SCOPED set gates behind ONE unbuilt artifact — a prolate-native
≥4-electron two-center CI (LiH 4e, NaH ~12e); only the 2-electron case (H₂) exists, and
nothing shows a ≥4e one binds or reaches accuracy. Reactivating any MIS-SCOPED approach
needs a fresh positive test, never the removal of the old wall
([[feedback_validate_before_reducing]]).

**`rests-on:` vocabulary:** COMPOSED / PK / NESTED / CONCAT / LÖWDIN / HYDROGENIC-PER-N /
ADIABATIC-HYPERSPHERICAL / GEOMETRY-2FOCI / TWO-BODY-CUSP / SYMMETRY-SPARSITY /
PROLATE-NATIVE.

### MIS-SCOPED — proven on dropped machinery; do NOT read as "now works" (25 rows)

| Row(s) | rests-on | guardrail | anti-laundering: fresh test to reactivate |
|:--|:--|:--|:--|
| W1e / second-row cluster: P1–P16, P21 (PK barrier, screened-Schrödinger, multi-zeta, 3-bucket, kernel-shape, rank-1 PK, mean-field J–K, Schmidt, [Ne] correlation, max_n=4, DMRG-FCIDUMP, LiH-qubit-binds, kwargs, explicit-core HF, off-diag cross-block h1, spectral-action, NaH Z_orb) | COMPOSED / PK / LÖWDIN | P15 adj. (Track DF) | Root (§3 lines 103/125): heteronuclear binding needs cross-center coupling → Löwdin across incompatible per-center coords → destroys Gaunt sparsity. Prolate has one shared coord + native orthogonality, so the dropped component is exactly what these rest on. **Does NOT mean NaH/LiH bind.** Test: build a prolate explicit-/all-electron LiH/NaH CI with Neumann V_ee and *measure* binding + R_eq. Unbuilt. |
| C2 — PK modifications (×6) | PK | — | Prolate 4e LiH is still a hard 4e two-center CI; PK negatives don't make it easy, only don't bear on the PK-free route. Test: prolate explicit-core LiH CI. |
| C4 — inter-group antisymmetry (×3) | COMPOSED | — | Prolate 2e has one coordinate system → antisymmetry native (shown for H₂). ≥3e needs a shared-coordinate antisymmetrized prolate CI; unbuilt. |
| C8 — l_max via 2D solver | COMPOSED / PK | — | Prolate angular convergence is structured (split-region Legendre, 3j-terminating; Paper 12 azimuthal 92.4→99.1%), not divergent. |
| C10 — TC Jastrow in adiabatic solver (disagreement #2: OPEN-LEANING → MIS-SCOPED) | ADIABATIC-HYPERSPHERICAL + TWO-BODY-CUSP | — | Row's own lesson: "needs a direct variational/FCI framework, not adiabatic." Prolate Neumann-CI *is* direct-CI. Test = P37/P38 (partially run; caps at basis size, needs the Neumann A_n/B_l/X_l moments extended for r₁₂). |
| C15 — coupled composition (cross-block ERIs) | COMPOSED / PK | — | Its Pauli/1-norm/29%-error came from a composed basis lacking two-center h1; prolate builds two-center integrals natively. Test: prolate ≥4e CI. |
| C16, C17 — single-center / charge-center nested LiH | NESTED | **8-9 + Track DF** | Single-center-basis limitations; prolate uses the two-foci coordinate natively. Papers 8-9's scoped theorem UNAFFECTED (still governs single-center Sturmian). Test: prolate 4e LiH CI. |
| C23 — balanced+frozen-core **NaH only** (MgH₂ excluded, triatomic; disagreement #4: row split) | COMPOSED / frozen-core | — | "Frozen [Ne] hides core screening" is a composed frozen-core artifact. Test: prolate explicit-core NaH CI (~12e — enormous, unbuilt). Strong caveat: nothing suggests tractable. |

### SUPERSEDED (1)
**C1 — LCAO / single-S³ / graph concatenation** (disagreement #3). MIS-SCOPED w.r.t.
prolate *and* architecture-locked corpus-wide (CLAUDE.md §2, v0.9.x). Guardrail
**FCI-M** — recorded superseded, **not retired**; FCI-M still governs any future
graph-concatenation proposal. (Prolate H₂/H₂⁺ already bind at machine precision /
99.8%, so the "R-independent kinetic energy → monotone PES" wall describes a dropped
architecture.)

### STANDING — bind the new program, or hold geometry-independently (23)

| Row(s) | rests-on | why it travels |
|:--|:--|:--|
| **Tier 1 — the real constraints:** P35 (exact-up-to-one-seed two-center ERI), P37 (r₁₂/FD on prolate 2e), P38 (r₁₂+Neumann), P39 (two-block radial exponent) | PROLATE-NATIVE | Native to the target; directly constrain it. P38: extending r₁₂ needs the Neumann moments extended, not the existing path. |
| **Two-foci ceiling / three-center genus jump** | GEOMETRY-2FOCI | The hard architectural cap — the prolate program stops at diatomics. |
| Cusp: C3, C11, C12, C14 + Group-C P26–P33 | TWO-BODY-CUSP | The e-e cusp is a two-body wavefunction feature, basis/geometry-independent (Paper 18 embedding tier). Re-scope only the hyperspherical (α,θ₁₂) parametrization; keep the substance. |
| Graph-native He: C20, C21; off-axis atom: C19, C22 | graph-native-S³ / atomic | Never touched composed; the swap has no purchase. |
| P22 (SO(4)-Wigner-D probe, falsified), P23 (overlap-slope tilt, curve-fit artifact) | diagnostic | Geometry-independent falsified diagnostics. |
| **P34 — k_n = Z/n completeness** (disagreement #1: OPEN-LEANING → STANDING) | HYDROGENIC-PER-N / completeness | A basis-completeness theorem ("non-orthogonality is the price of completeness"; bound per-function-scaled sets plateau ~60 mHa) that travels to *any* per-function-scaled orthonormal set incl. prolate; STANDING is the conservative call. Its per-BLOCK sibling (P39) already settled that per-block does not cost completeness. |

### STANDS-BUT-ORTHOGONAL — still true, but bind qubit-resource, not bond length (12)
C9 (Sturmian-CI 1-norm; guardrail **8-9**), C13, C18 (guardrail **Track DF**), P17–P20
(relativistic tapering — also selection-rule facts true in any basis), P24, P25 (Löwdin /
non-orthogonal encoding; guardrail-adjacent Track DF), P32, P36, P40.

### OPEN-LEANING — premise changed, untested (1)
**C5 — full N-electron radial solvers** (adiabatic / coupled-channel / 2D). Hyperspherical-
solver walls; prolate uses Neumann-CI. The "angular basis is the bottleneck" lesson may
or may not travel; needs a prolate re-test.

### HELD — STANDS-by-default, PI-directed not-classified (2)
**C6** (geometric elevation) and **C7** (diagnostic arcs — ledger says only "see
CHANGELOG"). Ledger text too thin to classify faithfully; not guessed. Revisit needs the
CHANGELOG detail.

### Guardrail adjudication (PI-confirmed 2026-09-19)
The 5 core + 3 adjacent guardrail rows (C1 FCI-M; C9, C16, C17 Papers 8-9 / Track DF;
C18 Track DF; adjacent P15, P24, P25) are recorded as **"does not constrain
prolate-native diatomics" — NOT retired.** Each guardrail still governs its own domain:
Papers 8-9's scoped Sturmian theorem still owns single-center encodings; FCI-M still owns
graph concatenation; Track DF still owns nested single-Hilbert-space. Re-scoping-for-
prolate ≠ weakening the guardrail.

### Tally
MIS-SCOPED 25 · SUPERSEDED 1 (C1) · STANDING 23 · STANDS-BUT-ORTHOGONAL 12 ·
OPEN-LEANING 1 (C5) · HELD 2 (C6, C7). Category rows bundle multiple sub-instances, so
these are *row* counts, not approach counts. **Excluded as triatomic+ (STAND by
default):** polyatomic coupling ×3, Poly-0/Poly-2 three-center, gerade lever, Woodbury
SW-penalty, MgH₂. **Excluded as non-chemistry:** the QED / gravity / periods / α /
nuclear / Lorentzian / MZV rows.

### Falsifier / what would move these
Any MIS-SCOPED row flips back to STANDING the moment a prolate-native ≥4e two-center CI
is built and *fails* for the same reason the composed row failed — i.e. if the
obstruction proves geometry-independent after all. The audit's integrity rests on that CI
being built before any MIS-SCOPED approach is treated as live. The `rests-on:` column
above is the first instance of the negative-side dependency edge; if it earns its keep,
promoting it to a standing field is the general fix for wall-staleness under future
component swaps.

### The load-bearing takeaway
The audit did not clear a path; it showed the composed/PK/nested ledger is about
machinery the prolate direction leaves behind, so it neither obstructs nor endorses that
direction. The walls that actually bind are the recent PROLATE-NATIVE ones (Tier 1) and
the two-foci ceiling. **Unobstructed by old evidence, entirely untested — gated on one
unbuilt prolate ≥4e CI.**

---

## Cluster: PERIODS-TRANSCENDENCE — what sets the transcendence class of a two-electron integral (worked 2026-10-16, v5.20.18, PI-invoked)

**Scope.** The walls on *what number class a GeoVac two-electron integral can land in*, and on the routes that tried to lower it. Current-state sources checked for every member: Paper 59 (`sec:genus`, `sec:obstruction`, `sec:f12`, `sec:scope`, `sec:beyond`), Paper 61 (`sec:modular`, the 64-digit guarded PSLQ), Paper 12 (`sec:r12`), CHANGELOG v5.1.10 → v5.20.17, `docs/predictions/2026-10-13_*.md`.

### Member walls

| Wall (ledger / owner) | Status | Kind | Evidence that decides it |
|:--|:--|:--|:--|
| Brute high-precision T2 by integrating the fibre over the family (ledger 2026-08-20) | **MIS-SCOPED** → route-specific | SOFT | The digits the row said were unreachable (≥30) were reached by a different route in the same release: the oscillation-free tail-analytic evaluator and then the (k,w) refactorization took T2 to 64 digits (CHANGELOG v5.1.10 "fibre precision wall removed"; Paper 61 `sec:modular`, "decided at 64 digits"). The wall is on *fibre-over-family quadrature*, not on the precision. Corrected scope: that route caps near 14 digits; use the Bessel-moment evaluators. |
| ℚ(i) seam (P56 Dirac ↔ P59 T2) as a shared mechanism (ledger 2026-08-21, 4 instances) | **STANDING** | HARD | Its stated upgrade condition was "β(2) measured in T2 at ≥32 digits". Paper 61's guarded PSLQ is decisively negative through weight three, including the Catalan-completed ring, at 64 digits and height ≤10: T2 carries no measured conductor-4 content at reachable heights. The seam stays CONVERGENT; the structural leg (a compact-but-pure-Tate scalar sector falsifies the Wick pun) is unchanged. |
| Seam involution: τ = i pinned by a physical ℤ₂ (ledger 2026-08-21) | **STANDING**, refined | HARD | Still no ℤ₂ pins τ = i. Refinement from Paper 61 `sec:modular`: τ = i *does occur*, as the in-domain CM fibre ρ = 1/2 of the Legendre family (period Γ(1/4)²/(4√π), Chowla–Selberg), i.e. as one fibre among a continuum, not as a fixed point of a physical involution. The row's "ρ ↦ 1/ρ, fixed locus the pure-Tate cusp" stands. |
| Elliptic/CM-adapted basis for chemical accuracy (ledger 2026-08-23, 3 instances) | **STANDING** | HARD | Every accuracy gain since has come from the electron-electron cusp, never from radial exponent placement: explicit-r₁₂ H₂ to 0.053 mHa with exact integrals (Paper 12, v5.14.10); LiH to −8.0585 via geminal shape and core e-n factors (v5.16.0, v5.20.6); the two-block radial exponent worth ~0 (v5.14.2). The period of the *integral* is not a *basis* principle. |
| Elliptic-dilogarithm route obstructed; the rank-4 Picard–Fuchs operator is irreducible (Paper 59 `sec:obstruction`, [SYMBOLIC + MEASURED]; not a §3 row — a paper-owned negative) | **STANDING** | HARD | Fourier–Laplace transform of an irreducible connection + exact-integer monodromy with no invariant coordinate subspace (Paper 59 abstract and `sec:obstruction`; `tests/test_routeC_momentum.py` `test_L4_monodromy_certificate_irreducible`). No factorization closes T2. |
| Partial closed form of the three-centre block is worse than none (Poly-2, ledger 2026-08-14) | **STANDING** | HARD | Paper 59 `sec:scope`: "a value that is part closed-form and part numerical is neither more accurate nor decidable"; the cancelling T1/T2 contributions are physics (8-fold symmetry bit-exact). Supports the all-or-nothing consequence below. |
| Generic three-electron chain fibre has no relation in K3 + elliptic + Tate bases (P3e-2, 2026-10-14; filed this run, ledger + §3) | **STANDING** | OPEN-LEANING | 320-digit PSLQ, one mass triple (1,2,4), basis omits dilogarithms of mass ratios and products of elliptic periods; the positive control recovers its exact relation in every basis (`docs/predictions/2026-10-13_three_electron_k3.md`; CHANGELOG v5.20.15; Paper 59 `sec:three_electron`, [OBSERVATION, evidence]). Named better attempt: widen the basis before widening the computation. **Update 2026-10-20 (pre-registered P3e-2W, `docs/predictions/2026-10-19_p3e2_widened.md`):** the named attempt was made. Five bases (adding dilogarithms, trilogarithms, Bailey-adapted products of elliptic periods, the vacuum sunrise, logarithms times periods) at four generic triples, a swap test, and two threshold triples: no relation, 17 of 20 generic cells valid; on the slice G(t;1,1) the fibre is a K3 function integrated once (Paper 59 `sec:density_ring`). Still STANDING; the kind (OPEN-LEANING vs SOFT) is left to the next `/walls` pass. |

**Adjacent, not members.** The Yukawa-PSLQ null (ledger 2026-06-03) and the W1e period-class null (2026-06-04) are "calibration data are not outer-factor periods", Paper 57's two-family finding, a different mechanism (inner-factor input data), left for an ALPHA/INNER-FACTOR pass. The Hain–Brown period-level negative (`memory/hain_brown_identification.md`) belongs with the P56 Tannakian arc. The S^(4)-stage numerics rows (Levin on log-modulated summands, nested nsum, sumem, float64 contamination, sympy rank; ledger 2026-06-11/13) are a METHOD family (how not to sum), not a transcendence family; listed here only so the next pass does not re-sort them.

**The l seam and m seam (owed by `docs/atlas.md`).** Neither needs a new row. The l-seam wall — cross-centre l-selection dies at the second nucleus — is the `two_kinds_of_sparsity` crystallization already inside CHEM-ACCURACY (Paper 58 `thm:abelian_residue`, seam J1). Bending's release of m is a *door* with prior art (Renner, Walsh; seam J5), not a wall. The water closed-form wall (seam J6) is this cluster's member 5/6 plus the GEOMETRY-2FOCI row in the architecture-swap audit.

### Crystallization (FORMING → CRYSTALLIZED)

**One sentence.** For an integral over Slater densities with rational kernels, the transcendence class is fixed by the density topology the kernels couple — zero or one distinct two-centre density gives genus zero, two give genus one (same-pair densities, the two-centre exchange class, build the same fibre family yet close in the two-centre seed set, so the fibre genus alone does not decide the class), and a second kernel meeting at a vertex lifts the curve to a surface (K3) — and not by the operator between the densities, the basis, the number of nuclei, or any factorization; so nothing short of changing what the kernels couple lowers the class, the three-centre radial core has no closed form through any factorization (whether the integrated integral has one is open), and a part-closed-form evaluation is worse than none. The rule is a pattern across the tested instances, not a derivation (seam J9, Obs; re-tiered 2026-10-16 by the /qa claim-search DELTA).

**Independent support (≥3).** (i) Paper 59 `sec:f12`: geminal kernels leave the curve unchanged (operator-independence), [SYMBOLIC] at one centre, test-backed at several. (ii) Four-centre P1 confirmed: no third root, same genus-1 family (`tests/test_paper59_four_centre_genus.py`); Paper 59 weighs it as a calibration, nearly a corollary of the momentum form, not an independent test. (iii) The CM-adapted basis null (the difficulty does not move with the basis). (iv) The irreducible PF operator (no factorization). (v) The K3 chain: the class rises exactly when a second kernel shares a vertex (`tests/test_paper59_three_electron_k3.py`). (vi) Poly-2: partial closed forms cancel against what they omit.

**Falsifier.** A genus-two fibre from a kernel coupling only two densities, or a third independent root at four centres (P1's own falsifiers); or a basis whose exponents placed by the CM fibre beats generic variational placement at matched size.

**Operational consequence (dispatch rule).** (A) *Class question* — count the distinct two-centre densities the kernel couples and the kernels sharing a vertex; read the rung off that, not off the molecule. (B) *Accuracy question* — evaluate (shared-centre channel sums reach ~10⁻⁹ by L ≈ 12); do not chase a closed form of the three-centre radial core by factorization, and never ship a partial one. (C) *Digits for an integer-relation search* — use the Bessel-moment / (k,w)-refactorization evaluators of Paper 61, not fibre-over-family quadrature. (D) *Next-rung question* — widen the PSLQ basis (dilogarithms of mass ratios, products of elliptic periods, more mass triples) before widening the computation.

**Promotion.** None needed: the statement is already a paper claim (Paper 59 `sec:f12`, extended by `sec:beyond`) and an Observation seam row (J9; re-tiered from F by the 2026-10-16 /qa claim-search DELTA); this pass records the cluster as its negative-side support and ties the walls to the guards.

### Delta 2026-10-27 (v5.21.19, `/walls PERIODS-TRANSCENDENCE`, PI-invoked) — the operator clause re-scoped

**Current state checked.** Paper 59 `sec:f12` (the geminal-master paragraph, v5.21.18), `sec:density_ring` (the γ = 1 mechanism; the Gaussian geminal; the off-slice kite and Y), `sec:three_electron` (the second slice); CHANGELOG v5.20.19 → v5.21.18; records `2026-10-22_geminal_shapes`, `2026-10-23_gamma1_mechanism`, `2026-10-23_gaussian_geminal_triangle`, `2026-10-27_h2_geminal_elliptic`, `2026-10-21_generic_slice`.

**The crystallized statement: MIS-SCOPED in one clause.** "Never by the operator" was too broad. Four results, each in the owning paper, show the kernel moving the class:

| Result | What the kernel does | Owner | Guard |
|:--|:--|:--|:--|
| GA2 (v5.21.7) | A Gaussian geminal's avatar has an irregular point at infinity: it leaves the Picard–Fuchs ladder for exponential periods; Slater avatars stay on it | Paper 59 `sec:density_ring` | `tests/test_paper59_geminal_shapes.py` |
| GY3 missed (v5.21.9) | With e^{−r²} on both edges the correlated triangle loses its K3 piece: "the K3 rung of the Slater triangle is not set by the densities and the Coulomb line alone" | Paper 59 `sec:density_ring` | `tests/test_paper59_gy_gaussian.py` |
| MG1–MG3 (v5.21.8) | The correlated triangle's order drops at γ = 1 because the geminal pair and the density pair then have equal masses; at γ = 1/2, 3/2, 2, 3 it does not | Paper 59 `sec:density_ring` | `tests/test_paper59_mg_mechanism.py` |
| HG1–HG4 (v5.21.18) | A Slater geminal (a rational kernel, 4π/(k²+w₁²) with the 1/r₁₂) lifts the two-centre master integral's operator from Bessel to genus one, though both densities sit on one centre pair. Operator level only; the physical integrals are its exponent derivatives | Paper 59 `sec:f12` | `tests/test_paper59_h2_geminal_master.py` |

What still holds, checked: at the level of the curve a rational kernel does not change the curve the densities build (Paper 59 `sec:f12`, [SYMBOLIC], `test_f12_multicenter_shares_eri_elliptic_genus`); the basis, nuclei-count and factorization clauses are untouched by anything since 2026-10-16.

**Re-crystallized sentence (replaces the 2026-10-16 one, which is kept below as history).** For integrals over Slater densities, the transcendence class is set by the integral read as a graph: the density topology the kernels couple, together with the masses on the kernel lines. A massless Coulomb kernel adds no threshold, so for Coulomb integrals the density topology alone sets the rung: zero or one distinct two-centre density gives genus zero, two give genus one, same-pair densities close in the two-centre seed set, and a second kernel at a vertex lifts curve to surface (K3). A massive kernel is a line with its own mass. A Slater geminal keeps the curve the densities build but can lift a same-pair master to genus one, and it degenerates when its mass equals a density pair's (γ = 1). A Gaussian geminal is not a propagator of this kind and leaves the Picard–Fuchs ladder. The class is never set by the basis, the number of nuclei, or a factorization. It is a pattern across the tested instances, not a derivation.

*Why this reading is not new physics.* Paper 59 `sec:density_ring` already reduces the kite as a two-loop D = 3 Feynman integral, with density lines, squared geminal lines and a massless Coulomb line (IB1–IB4, v5.21.10). The re-scope just reads the cluster through that dictionary: thresholds come from massive lines, and the Coulomb line has none.

**Falsifier (widened).** As before: a genus-two fibre from a kernel coupling two densities, a third root at four centres, or a CM-placed basis beating generic placement. Added: a Coulomb-only integral whose class depends on anything but the density topology, or a Slater-geminal integral whose class does not change at an equal-mass value of its exponent.

**Dispatch rule (A) amended.** *Class question:* count the distinct two-centre densities, the kernels sharing a vertex, **and the massive kernel lines with their masses**. Read the rung off the graph. Read it off the density topology alone only when every kernel is Coulomb.

**Member status updates.**
- Generic three-electron chain fibre (P3e-2 / P3e-2W): **STANDING, kind settled to SOFT (basis-limited, cause located).** The kind was left open on 2026-10-20. The second slice (t;1,2) runs through two of the null triples. Its operator has an order-four piece R, equivalent by an exact map to the operator R_V of ∫K₀(tx)K₀(x)K₀(2x)dx (GS2 held, CHANGELOG v5.21.4). R_V is irreducible and orthogonal (v5.21.14). None of the five P3e-2W bases contains R_V's periods, so the null is what a basis missing that piece predicts. The next attempt is named: put R_V's periods (the three-mass Bessel moment) in the basis.
- Every other member: unchanged, evidence as recorded on 2026-10-16. Nothing since has touched the irreducible PF operator, Poly-2, the seam involution, the ℚ(i) seam, or the CM-basis null.

**Loci re-scoped in the same pass (summary surfaces move with the claim).** Paper 59 `sec:f12`'s opening sentence ("governed not by the operator but by the density topology") contradicted the same section's master paragraph and `sec:density_ring`. It now reads "at the level of the curve", with pointers to both. The same cap was added to seam J9 (`docs/forced_free_seam.md`), Paper 57's J9 item and atlas row 19. The group 2 synthesis was already scoped in v5.21.18.

**Promotion.** None needed. The sharper statement is assembled from paper-owned results, and J9 stays an Observation. The PI decides whether Paper 59 should state the graph reading (massive vs massless lines) as its own [OBSERVATION] sentence. This pass did not add it.

**Promotion decided (2026-10-27, PI-delegated: "Can you make those decisions and enact them").** The graph reading is now its own [OBSERVATION + PRIOR ART] paragraph closing Paper 59 `sec:f12`, with the falsifier above and four guards. Prior art credited: internal masses set a Feynman integral's geometry (Bloch–Vanhove; Vanhove 2014). Tier stays Observation. The section's remaining "lives in the densities" sentence was restated at the curve level and points to it.

*Superseded 2026-10-16 one-sentence statement* (kept verbatim above under "Crystallization"; its "and not by the operator between the densities" clause is the one re-scoped here).

### This run's delta (2026-10-16)
- PERIODS-TRANSCENDENCE: *forming (unworked)* → CRYSTALLIZED; members 5 → 7.
- Brute-T2 row: STANDING → MIS-SCOPED (route-specific), with the 64-digit evidence.
- Seam-involution row: refined (τ = i occurs as the in-domain CM fibre, not as a pinned point).
- New member filed: P3e-2 generic chain fibre (ledger + CLAUDE.md §3 row), OPEN-LEANING.
- The atlas's owed "l-seam / m-seam rows" discharged by pointer (no new wall exists there).

---

## Note 2026-10-19 (v5.20.29, PI-delegated, not a `/walls` pass): HF-4 priced, kept as a multi-focal instance

The PI delegated one decision in session ("You choose"): whether the hyperfine door moves HF-4 out of the
multi-focal composition wall ([[multi_focal_wall_pattern]]; Paper 57 `sec:unity`).

- **Wall: STANDING**, evidence: Paper 23 `sec:finite_nucleus` ("A finite magnetization", 2026-10-19); record
  `docs/predictions/2026-10-19_zemach_door.md` (ZM1-ZM3 hold, ZM4 fails in part); seam J16.
- **HF-4 instance: MIS-SCOPED, re-stated.** Its text ("no magnetization-density distribution carried") names a
  limit of one code function. With the distribution carried, the electron's labels survive and the hyperfine
  factor is set by the Zemach moment of charge and magnetization together, one joint moment per order [superseded 2026-10-19: further joint constants as NRQED matching
  coefficients; see the addendum below], which no
  label fixes. That is the J14 shape (PI decision 2026-10-18: J14 is a door and the seventh instance), so HF-4
  stays an instance, priced. The wall's mechanism is supported, not breached: two spatial distributions compose
  at leading order only through a joint moment, the Zemach radius, that neither radius fixes alone.
- **Classification unchanged** (HARD at the level of the mechanism; the door is priced, not closed).
- **Addendum 2026-10-19 (EFT check, `debug/eft_multifocal_check_memo.md`, sources verified by the PM):** for
  the nuclear instances (J14, HF-4/J16) the price is the NRQED matching statement (Hill-Paz, PRL 107, 160402
  (2011); Karshenboim, Phys. Lett. A 225, 97 (1997)): one coefficient per independent contact operator at each
  order, logs fixed by lower orders. "One moment per order" was too strong (two at fourth order) and "universal
  ln + 2 gamma_E" holds only for the log's coefficient. These instances support the wall's existence, not that
  it is GeoVac-specific; the GeoVac-specific content is that the label skeleton has no native matching step
  (Paper 57 `sec:unity`, on the Paper 32 LS-8a template). X4-1 (pre-registered) confirmed the fourth-order
  n-dependence in our own code to 5e-9.

## Promotion candidates (PI-gated)

**CHEM-ACCURACY two-walls compounding.** One-sentence mechanism: *chemical accuracy (Wall A, breachable) and native sparsity (Wall B, hard) co-locate favorably only at one center, so breaching the cusp for molecules re-triggers the non-orthogonality wall.* Support: ≥3 independent negatives (R12-CI conditioning v5.0.7; flexible-Jastrow cap v5.0.8; elliptic-basis null v5.0.1) + 1 breach (atomic xTC v5.0.9) + 2 crystallized meta-findings it subsumes. Falsifier: above. **Proposed home:** a sharpening of [[two_kinds_of_sparsity]] and/or a Paper 14 remark tying the atomic xTC win to the molecular tension. **PI decision required** to write it into a paper.

---

## Backlog — clusters to work next
Run `/walls <cluster>` on any of: QED-SELECTION-RULES · ALPHA-COMBINATION · GRAVITY-SUBSTRATE · RELATIVISTIC-TAPERING · NUCLEAR-CONTINUUM-RESPONSE · MOLECULAR-BINDING-W1e. Each needs its members pulled from §3, re-verified against current state, tagged hard/soft, and a crystallization attempted.
