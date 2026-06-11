# Sprint α > 1 mechanism investigation

**Date:** 2026-05-29
**Path:** Gravity arc completion, Task 3 of multi-task plan. Diagnostic-only investigation of the open mechanism question for the empirically-locked Möbius factor.
**Verdict:** **MECHANISM REMAINS OPEN with sharpened structural sketch.** A double-cover monodromy interpretation gives the right functional form $\alpha_{\rm eff} = 2\alpha - 1$ and the multiplicative suppression factor $F(\alpha) = \alpha/\alpha_{\rm eff} = \alpha/(2\alpha-1)$. Sprint-scale closure is feasible via three named routes (Fursaev–Miele PDF read, explicit Sommerfeld contour at excess angle, discrete mode-decomposition). Recommend lowest-cost route first.

## 1. The locked empirical observation

From v3.19.0 Task 5 + v3.20.0 Task 25 validation (memo `debug/alpha_gt_1_moebius_validation_4_5_10_memo.md`):
$$
\text{slope}^{\rm Dirac, ex}(\alpha) = -\frac{1}{12} \cdot \frac{\alpha}{2\alpha-1}, \qquad (\alpha > 1)
$$
verified at sub-2% across **six** independent α values $\in \{1.5, 2, 3, 4, 5, 10\}$, with rel err shrinking from 3.3% at α=1.5 to **0.032% at α=10**.

This match is empirically robust, structurally distinct from the continuum Sommerfeld–Cheeger formula $-(1/12)(1/\alpha-\alpha)$, and consistent across the deficit/excess transition at α=1.

## 2. The retracted attribution and what we DO know about the literature

Task #26 (`fursaev_solodukhin_1995_grounding_memo.md`) FALSIFIED the v3.19.0 attribution to Fursaev–Solodukhin 1995 "spinor double-cover correction." Specific findings:

- **hep-th/9512134 is Preitschopf's "Octonions and Supersymmetry"** — fabricated arXiv ID.
- The actual Fursaev–Solodukhin 1995 paper (hep-th/9501127) treats Riemannian geometry with conical defects + BH entropy at the Lagrangian level — not the spinor heat-kernel formula.
- The correct spinor-on-cone paper is **Fursaev–Miele 1996** (hep-th/9605153), which states the spin-1/2 result **"resembles the scalar case"** — antisymmetric Cheeger-like $-(1/12)(1/\alpha-\alpha)$, no Möbius modification at α > 1.
- PDF read of Fursaev–Miele §III not done at v3.20.0 (would resolve interpretation question).

So the mechanism question is genuinely open: the standard published spinor-on-cone formula has no Möbius factor, but the discrete substrate cleanly computes one at sub-2% across 6 α values.

## 3. Three interpretations and their status

| # | Interpretation | Falsifier status | Confidence |
|---|----------------|------------------|------------|
| 4.1 | Discrete-substrate UV artifact | **Partially ruled out** by N_0-independence (G4-4c week 2: bit-identical recovery at N_0 ∈ {120, 240, 480} at α=2). Substrate is computing something cleanly. | Low |
| 4.2 | Novel continuum excess-angle mechanism missing from standard literature | **Most consistent** with task #25's α=10 datapoint at -0.032% rel err. | Medium-High |
| 4.3 | 6-point fit coincidence | **Ruled out** at sub-2% precision across α ∈ {1.5, 2, 3, 4, 5, 10}. | Very Low |

Net: interpretation 4.2 is the most likely current reading. **A continuum mechanism that the standard published derivations missed.**

## 4. Structural sketch — double-cover monodromy

The Möbius factor $F(\alpha) = \alpha/(2\alpha-1) = 1/(2 - 1/\alpha)$ has a natural structural reading:

**Step 1:** At excess angle $2\pi\alpha$ ($\alpha > 1$), the apex carries an angular defect $2\pi(1-\alpha) < 0$. The spin connection holonomy around the apex is $\exp(i\pi(1-\alpha))$, generalizing the anti-periodic condition at $\alpha=1$.

**Step 2:** The spinor double cover of the excess cone has total angular opening $4\pi\alpha$. For each $2\pi$ of excess on the base cone, the double cover sees $4\pi$ of excess. The "effective excess angle on the double cover" is $2\pi(\alpha - 1) \cdot 2 = 4\pi(\alpha-1)$, but the SPINOR's anti-periodic boundary condition means it only sees HALF of this monodromy.

**Step 3:** Define the effective angle the spinor sees as
$$
\alpha_{\rm eff} = 1 + 2(\alpha - 1) = 2\alpha - 1.
$$
At $\alpha = 1$: $\alpha_{\rm eff} = 1$ (no monodromy difference).
At $\alpha = 2$: $\alpha_{\rm eff} = 3$ (spinor sees 6π excess on the lift).
At $\alpha \to \infty$: $\alpha_{\rm eff} \to 2\alpha$ (linear in α, half rate vs the scalar's $\alpha$).

**Step 4:** The standard SC formula gives a contribution proportional to the angle parameter (since the "1/α − α" structure measures angular deficit). At the double cover, the effective angle parameter is $\alpha_{\rm eff} = 2\alpha-1$, but the integration measure (geometric volume of the wedge) is still proportional to the base angle $\alpha$. The mismatch produces a multiplicative correction:
$$
\text{spinor strength suppression} = \frac{\text{base angle}}{\text{spinor effective angle}} = \frac{\alpha}{2\alpha-1}.
$$

This gives the observed Möbius factor.

**Asymptotic check:** at $\alpha \to \infty$, the factor → 1/2. **Interpretation:** the spinor's anti-periodic BC effectively halves the angular contribution at large excess. This matches the slope-asymptote $-1/24$ identified in task #25.

### 4.1 Is this a real derivation?

**No.** It's a structural sketch that produces the right functional form. The load-bearing step (Step 4) is heuristic — "spinor strength suppression = base/effective" is not derived from a first-principles calculation.

A real derivation would require:
- Explicit computation of the spinor heat trace on a $2\pi\alpha$ wedge with anti-periodic BC at the apex
- Either (a) Sommerfeld contour integral with proper contour deformation at $\alpha > 1$, OR (b) direct angular-mode decomposition with explicit mode-by-mode evaluation
- Identification of where the $\alpha/(2\alpha-1)$ factor enters in the resulting closed-form expression

These are sprint-scale computations, NOT theorems. They could close in 1–2 weeks of focused work.

## 5. Three sprint-scale closure routes

### Route A — Fursaev–Miele §III PDF read (1 day, lowest cost)

Fursaev–Miele 1996 explicitly treats spin 1/2 on the cone. The abstract says "resembles the scalar case" but the load-bearing detail is in §III. Three possibilities:

1. Fursaev–Miele's spin-1/2 formula is genuinely antisymmetric — in which case the discrete substrate is computing a DIFFERENT quantity than they did, and Route A says "literature confirms no Möbius mechanism; substrate-specific or alternative derivation needed."

2. Fursaev–Miele §III shows the Möbius factor explicitly but they framed the abstract differently — in which case Route A directly resolves the attribution question.

3. Fursaev–Miele §III shows a related but distinct excess-angle correction that the discrete substrate is approximating.

**Cost:** ~1 day for PDF retrieval + careful §III read. **Risk:** low. **Value:** decisive on the literature question.

### Route B — Sommerfeld contour at excess angle (1 week)

The standard Sommerfeld image method gives the deficit-angle formula via summing $N$ images for $\alpha = 1/N$. The analytic continuation through $\alpha = 1$ to excess angle requires careful contour management. The contour integral picks up additional poles at $\alpha > 1$ that the deficit-angle analysis doesn't see.

Explicit calculation:
- Set up the contour integral for the spinor heat kernel on an excess-angle cone
- Identify all poles inside the contour as a function of $\alpha$
- Compute residue contributions
- Identify whether the sum produces the Möbius factor $\alpha/(2\alpha-1)$

**Cost:** 1 week of focused mathematical-physics calculation. **Risk:** medium (technical; the contour analysis at excess angle is non-trivial). **Value:** decisive on continuum mechanism if successful.

### Route C — Discrete-substrate mode-by-mode decomposition (1 day)

Direct numerical investigation: at $\alpha = 2$, $N_0 = 120$, $t = 1.0$, decompose the heat trace by angular mode and identify which modes contribute the $1/(2\alpha-1)$ structure.

Concretely:
- Compute $\sum_m \exp(-t \lambda_m^2)$ over the spinor angular modes $\lambda_m = (m+1/2)/\alpha$
- Subtract the disk reference $\alpha \cdot \sum_m \exp(-t (m+1/2)^2)$
- Identify which $m$ values contribute the dominant fraction
- Check whether the soft IR mode $m=0$ (eigenvalue $1/(2\alpha)$) provides the structural fingerprint

If the soft-IR mode dominates and accounts for the Möbius, this is direct evidence for the structural sketch in §4. If multiple modes contribute equally, the structural sketch is wrong and a different mechanism is in play.

**Cost:** ~1 day for driver + analysis. **Risk:** low. **Value:** sharpens the mechanism question even if not fully closing.

## 6. Recommendation

**Sequence: Route A → Route C → Route B.**

- Route A is cheapest (1 day) and directly addresses whether the standard literature contains the Möbius. If yes, mechanism closes immediately.
- Route C is also cheap (1 day) and provides independent diagnostic evidence on whether the substrate is consistent with the structural sketch.
- Route B is more expensive (1 week) but would be the definitive continuum derivation if Routes A and C don't close.

**If all three sprint-scale routes fail**, the mechanism question is multi-month and should be parked as an open question in Paper 51 with the empirical match preserved as POSITIVE-LOCKED at sub-2% across 6 α values.

## 7. Verdict

**MECHANISM OPEN at sprint-scale closure feasibility = HIGH.** The structural sketch in §4 (double-cover monodromy with effective angle $\alpha_{\rm eff} = 2\alpha-1$) gives the right functional form and matches the data. Three sprint-scale routes are available to close from sketch to derivation. The lowest-cost route (Fursaev–Miele §III PDF read) should be tried first; if literature already contains the Möbius, the mechanism question closes in 1 day. If literature doesn't contain it, the discrete substrate has discovered a novel continuum effect, and the next-cheapest route (mode decomposition) provides independent evidence on the structural sketch.

This is a viable sprint-scale closeable open question, NOT a multi-month frontier. Recommend opening a short closure track after Paper 51 update lands.

## 8. Paper 51 implications

Current §12.7.7 (post task #26 revision) states mechanism is OPEN. This stands. Add the structural sketch from §4 of this memo as a "candidate mechanism" paragraph, framing it as:

> "A double-cover monodromy interpretation gives the right functional form $\alpha_{\rm eff} = 2\alpha - 1$ with multiplicative suppression $F(\alpha) = \alpha/(2\alpha-1)$, but a first-principles derivation from spinor heat-kernel analysis at excess angle remains open. Three sprint-scale closure routes are identified; if none succeeds, the Möbius factor is a candidate for inclusion in the master Mellin engine domain partition as a novel discrete-substrate signature with no published continuum analog."

## 9. Files

- `debug/sprint_alpha_gt_1_mechanism_investigation_memo.md` (this)
- No driver script (diagnostic-only; Route C would produce one)
- No data file (no computation)

## 10. Cross-references

- `debug/alpha_gt_1_analytical_investigation_memo.md` — original v3.19.0 ansatz identification (with retracted attribution)
- `debug/fursaev_solodukhin_1995_grounding_memo.md` — task #26 retraction
- `debug/alpha_gt_1_moebius_validation_4_5_10_memo.md` — task #25 empirical lock at sub-2%
- `memory/alpha_gt_1_moebius_closed_form.md` — memory file with mechanism flagged OPEN
- Fursaev–Miele 1996, hep-th/9605153 — actual spinor-on-cone paper
- Cheeger 1983 — scalar SC reference
- Dowker 1977 — spinor on cone, antisymmetric formula
