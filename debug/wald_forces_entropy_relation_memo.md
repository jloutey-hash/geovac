# Wald forces the entropy/action-$G$ relation; the G4-2 factor of 2 is not a cone coefficient

**Date:** 2026-05-30
**Path:** Conversational chase, gravity arc. Started from the "calibration = boundary condition" reframe: the sector-wise Mellin moment map (G4-5d) reads each gravitational constant as a moment of the cutoff/boundary profile $f$ — tip/entropy $\leftrightarrow \phi(0)$, Einstein--Hilbert/$G \leftrightarrow \phi(1)$, cosmological constant $\leftrightarrow \phi(2)$. The "overdetermination test" asked whether one boundary profile fits $G$, $\Lambda_{cc}$, and $S_{BH}$ consistently. It collapses to one question: the factor-of-2 between the two routes to Newton's $G$.
**Verdict:** **POSITIVE-STRUCTURAL (Wald) + ERRATUM (cone) + OPEN (convention audit).**

## 1. The two routes to $G$ and the factor of 2

- **Einstein--Hilbert route (G7), $\phi(1)$ channel:** $\Geff = 6\pi/\Lambda^2$ (Gaussian).
- **Entropy route (G4-2), $\phi(0)$ channel:** $S_{BH} = A\Lambda^2/(12\pi) \Rightarrow G_N = 3\pi/\Lambda^2$ via $S = A/(4G)$.
- Ratio $\Geff/G_N = 2$.

## 2. Overdetermination test: naive version is vacuous

Three observables pick out three different Mellin moments $\phi(0), \phi(1), \phi(2)$ of one function $f$. A function has infinitely many moments, so three numbers never overdetermine it — with $f$ free, "one $f$ fits all three" is automatic and empty (curve-fit trap, named explicitly). The test only has teeth where two observables share a constant. They do: $G$ is reachable both ways. So the entire content collapses to the factor of 2.

## 3. The factor of 2 is FORCED (Wald), not calibration

- G1--G2: the Dirac-sector spectral action is **two-term exact** — Einstein--Hilbert + cosmological constant, **no higher curvature**.
- Pure Einstein gravity $\Rightarrow$ **Wald's theorem** (Wald 1993): $S_{BH} = A/(4G)$ with the **same** $G$ as in $(1/16\pi G)R$. (The cosmological constant is a constant in the Lagrangian; $\partial/\partial R_{abcd}$ of a constant is zero, so it does not touch the horizon entropy.)
- Therefore the action-$G$ and entropy-$G$ are **identically equal**, by theorem. The factor of 2 between our two derivations is **forced to be a bookkeeping artifact** — not physical, not free calibration data.

**Consequence for the overdetermination test:** it PASSES, forced. $\phi(0)$ (entropy) and $\phi(1)$ ($G$) are not independent free data; Wald locks them. This sharpens G8: the *values* $\phi(k)$ are calibration (free $f$), but the *relations* among the gravitational constants are forced.

## 4. ERRATUM: the cone coefficient is NOT the factor of 2

The G4-2 memo and Paper 51 attributed the factor of 2 to "scalar vs Dirac conical-defect coefficient." **This is wrong.** The 2D Dirac and scalar cones share the conical coefficient $(1/12)(1/\alpha - \alpha)$, three independent ways:

1. **Central charge:** a free 2D Dirac fermion has $c = 1$; a free real scalar has $c = 1$. The 2D conical tip coefficient scales with $c$. Equal $c \Rightarrow$ equal cone coefficient.
2. **Fursaev--Miele 1996** ("spin $1/2$ resembles the scalar case") — already recorded as the gravity arc's own dead-end note (CLAUDE.md §3, task #26).
3. **Our own discrete substrate (G4-4c):** the spinor tip coefficient extracts to $-1/12$ bit-exact (and the anti-periodic half-integer spectrum is structurally essential; the scalar BC extraction was the one that came out garbage at 28--66%).

To lock $3\pi \to 6\pi$ the entropy route would need its cone coefficient *halved* (to $1/24$). It is not halved; it is $1/12$, same as scalar. So the cone is exonerated, and G4-2's attribution is internally inconsistent with the arc's own G4-4 data.

## 5. Where the factor of 2 actually lives (OPEN)

A normalization mismatch between the **4D-bulk** ($S^3 \times S^1$, G7) and **factorized 2D$\times$2D** ($D^2 \times S^2$, G4-2) spectral-action setups. Most plausible specific spot: the distribution of the 4D Dirac spinor trace $a_0 = 4 = 2 \otimes 2$ across the factorized product (G4-2 used a scalar-flavored cone $\times$ a 2-component Dirac sphere). Closing it is a **convention-alignment audit** (Str-vs-Tr, $D^2$ normalization, $a_0$ counting across both derivations), not a coefficient swap. Wald guarantees the audit must close $3\pi \to 6\pi$.

## 6. Honest scope

- **Banked:** Wald forces the entropy/action-$G$ relation (the overdetermination test passes, forced). Cone exonerated three ways. G4-2 attribution corrected.
- **NOT banked:** the explicit recomputation $3\pi \to 6\pi$. Doing it by tuning the cone coefficient would be curve-fitting the answer — and the cone is the one piece confirmed *not* to carry the 2. The honest close is the convention audit.
- **Does NOT compute the numbers:** passing the overdetermination test shows the boundary *description* is internally consistent; it does not determine $\phi(0), \phi(1), \phi(2)$, which still ride on the free $f$. The deep test (does the discrete edge *force* $f$) leans "reproduces, not forces" per L6 (replica-weight-harmless).

## 7. Files / edits

- This memo.
- `papers/group5_qed_gauge/paper_51_gravity_arc.tex` — §G4-2 parenthetical replaced with a Wald remark; summary one-liner corrected; `wald1993` bibitem added.
- `debug/g4_2_conical_replica_memo.md` — erratum pointer appended.

## 8. Cross-references

- **G4-2** (`debug/g4_2_conical_replica_memo.md`): the derivation whose factor-of-2 attribution this corrects.
- **G7** (`debug/g7_extremality_newton_memo.md`): the $\phi(1)$-route $G$.
- **G4-5d** (`debug/g4_5d_cutoff_dependence_memo.md`): sector-wise Mellin moment map.
- **G8** (`debug/g8_cutoff_dependence_memo.md`): cutoff-as-calibration; this memo sharpens it (relations forced, values free).
- **Wald 1993**, Phys. Rev. D 48, R3427.
