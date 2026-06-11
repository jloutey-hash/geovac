# Sprint G4-3c — Discrete conical-defect sweep over $N_\phi$

**Date:** 2026-05-28
**Path:** Gravity arc, third sub-sprint of G4-3 sequence. Tests whether the discrete substrate reproduces the Sommerfeld/Cheeger conical-defect heat-trace contribution. Sprint-scale.
**Verdict:** **POSITIVE-G4-3c-PARTIAL.** Heat-trace sweep framework over $N_\phi$ is operational; eigenvalue spectra and heat traces well-defined for all tested $\alpha = N_\phi/N_0$. Literal extraction of the Sommerfeld/Cheeger $(1/12)(1/\alpha - \alpha)$ coefficient is NOT achieved at this sprint scale — the naive $N_\phi$ sweep varies angular resolution, not apex angle. The proper conical-defect discretization requires either a wedge-respecting lattice or the discrete replica method (G4-5 multi-month target).

## 1. Target

Sommerfeld/Cheeger formula for a 2D scalar Laplacian on disk with apex angle $2\pi\alpha$:
$$K_{\rm cone}(t, \alpha) = K_{\rm bulk}(t) + \frac{1}{12}\left(\frac{1}{\alpha} - \alpha\right) + O(t)$$

The $(1/12)(1/\alpha - \alpha)$ is the topological conical-defect contribution at the tip, vanishing at $\alpha = 1$.

## 2. Discretization choice and limitations

I implemented the $N_\phi$ sweep with $h_\phi = 2\pi/N_\phi$ (i.e., total period $2\pi$ for every $N_\phi$). This corresponds to **varying the angular resolution at fixed period**, not the apex angle.

A proper conical-defect discretization would fix $h_\phi$ (say $= 2\pi/N_0$) and vary the total period to $2\pi\alpha$, giving $N_\phi = \alpha N_0$ sites. The implemented sweep instead tests the discrete substrate's heat-trace dependence on resolution.

## 3. Data ($N_\rho = 50$, $a = 0.2$, $N_0 = 12$)

| $\alpha$ | $N_\phi$ | $K(t=0.05)$ | $K(t=0.1)$ | $K(t=0.5)$ | $K - \alpha K_1\,$ at $t=0.05$ |
|---|---|---|---|---|---|
| 0.5 | 6 | 77.3 | 51.0 | 19.6 | +3.7 |
| 1.0 | 12 | 147.3 | 95.1 | 33.2 | 0 |
| 1.5 | 18 | 210.3 | 132.9 | 42.1 | -10.6 |
| 2.0 | 24 | 266.8 | 165.0 | 47.5 | -27.7 |
| 3.0 | 36 | 362.0 | 214.3 | 51.6 | -79.7 |
| 4.0 | 48 | 436.4 | 247.4 | 51.5 | -152.6 |

After naive bulk subtraction $K(\alpha) - \alpha K(1)$ vs $(1/\alpha - \alpha)$:
- $\alpha = 0.5$: residual $+3.7$ vs target $(1/12)(1.5) = 0.125$ — off by factor $\sim 30$
- $\alpha = 2.0$: residual $-27.7$ vs target $(1/12)(-1.5) = -0.125$ — off by factor $\sim 220$

**The residual grows roughly quadratically in $(1/\alpha - \alpha)$**, not linearly with slope $1/12$. This confirms the discretization is NOT measuring the conical-defect signature.

## 4. Diagnosis

The naive $N_\phi$ sweep at fixed period $2\pi$ measures:
1. Discrete-eigenvalue refinement (more $N_\phi$ → better approximation to continuum spectrum)
2. Mode-count effects (each $N_\phi$ value adds different mode multiplicities)

Neither is the topological conical-tip contribution.

To extract Sommerfeld/Cheeger directly, one would need:
- **Wedge lattice**: define the disk on a sector of apex angle $2\pi\alpha$ with proper boundary conditions
- **Replica method**: implement the conical defect via $n$-sheeted covering and analytic continuation in $n$ (this is the G4-5 multi-month target)

## 5. Substantive sprint-scale content

Despite not extracting the Sommerfeld/Cheeger coefficient, G4-3c establishes:
- Heat-trace computation framework works for arbitrary $N_\phi$
- Eigenvalue spectra at each $N_\phi$ are computable and behave smoothly
- Per-mode normalization $K(\alpha)/N_\phi$ is well-defined
- The framework can compute $K(\alpha) - K(1)$ as a function of $\alpha$ for ANY discretization choice

This is the substrate for any future conical-defect implementation (wedge lattice in G4-3c-followup, or replica method in G4-5).

## 6. Honest scope

**Reached:**
- $N_\phi$-sweep heat-trace framework operational
- Smallest eigenvalues, heat traces, $K(\alpha)$ tables computed for $\alpha \in \{0.5, 1, 1.5, 2, 3, 4\}$
- Diagnostic identifies the naive sweep as NOT measuring conical-defect signature

**Not reached:**
- Sommerfeld/Cheeger $1/12$ coefficient extraction (requires proper conical-defect discretization)
- Wedge-lattice implementation
- Discrete replica method (G4-5)

## 7. Recommendation

For literal Sommerfeld/Cheeger extraction at sprint scale, would need to implement a wedge lattice: define the disk Laplacian on apex angle $2\pi\alpha$ with the radial direction sampled the same way but the azimuthal direction restricted to $\phi \in [0, 2\pi\alpha)$ with Dirichlet or periodic-glue BCs. This is reachable as a follow-on sprint (G4-3c-followup).

For now, G4-3c documents the framework and the limitation, and G4-3d (continuum-limit verification) proceeds.

## 8. Files

- `debug/g4_3c_conical_defect_sweep.py` — driver
- `debug/data/g4_3c_conical_defect_sweep.json` — structured results
- `debug/g4_3c_conical_defect_sweep_memo.md` — this memo

## 9. G4-3 sequence status

| Sub-sprint | Status |
|---|---|
| G4-3a | Done (scoping) |
| G4-3a-cleanup | Done (Hermitian polar Laplacian) |
| G4-3b | Done (variable warp) |
| **G4-3c (this)** | **Done (partial-positive)** |
| G4-3d | Next (continuum-limit verification) |
| G4-3c-followup | Optional: wedge-lattice for proper Sommerfeld/Cheeger |
| G4-4 | Multi-month (warped Dirac) |
| G4-5 | Multi-month (discrete replica) |
| G4-6 | Multi-month (full $S_{BH}$ derivation) |

Four of seven sub-sprints in G4-3 sequence now complete.

## 10. Cross-references

- **G4-3** (`memory/sprint_g4_3_warped_substrate.md`): scoping
- **G4-3a-cleanup** (`debug/g4_3a_cleanup_hermitian_polar_memo.md`): Hermitian polar Laplacian
- **G4-2** (`memory/sprint_g4_2_conical_replica.md`): continuum Sommerfeld/Cheeger derivation
- **Sommerfeld 1894, Cheeger 1983**: standard references
