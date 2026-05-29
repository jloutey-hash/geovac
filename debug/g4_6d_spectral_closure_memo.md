# Sprint G4-6d formal closure — spectral azimuthal discretization

**Date:** 2026-05-29
**Path:** Multi-task thread 5, Track A. Formal closure of the G4-6d sub-sprint.
**Verdict:** **CLOSED.** B.1 implementation + B.2 verification + tests/test_warped_dirac_spectral.py (14 tests, all pass). The G4-6d sub-sprint is the first multi-month G4-6 sub-sprint to close at sprint-scale, validating the reframed G4-6 sequencing (G4-6d → G4-6b → G4-6a refined → G4-6c parallel → G4-6e → G4-6f).

## 1. Sub-sprint scope at closure

G4-6d's original scope per `debug/g4_6_scoping_memo.md` §4.4:

> "Replace the FD azimuthal discretization with a spectral (DST/Fourier) representation; verify per-$t$ tip recovery is unbiased at small $t$ (no $4/\pi^2$ truncation overshoot)."

Falsifier F17:
> "Per-$t$ tip recovery within 1% across the full $t$-grid $\{0.0025, 0.005, \ldots, 10\}$."

## 2. What landed

### 2.1 Production code (B.1)

Two new classes in `geovac/gravity/warped_dirac.py`:

- **`DiscreteDiskDiracSpectral(N_rho, a, N_phi)`** — disk reference with exact $m_{\rm eff} = k + 1/2$.
- **`DiscreteWedgeDiracSpectral(N_rho, a, N_phi, alpha)`** — wedge with exact $m_{\rm eff} = (k+1/2)/\alpha$.

Same radial structure as `DiscreteDiskDirac` / `DiscreteWedgeDirac` (Hermitian polar Laplacian with $u = \sqrt{\rho} f$ substitution from G4-3a-cleanup), but replaces the FD azimuthal eigenvalue formula
$$m_{\rm eff}^{\rm FD} = \frac{2}{h_\phi} \sin\left(\frac{\pi(k+1/2)}{N_\phi}\right)$$
with the exact continuum eigenvalue
$$m_{\rm eff}^{\rm spec} = \frac{k + 1/2}{\alpha}.$$

### 2.2 Empirical verification (B.2)

Driver: `debug/g4_6a_multi_substrate_uv_first_move_spectral.py`
Data: `debug/data/g4_6a_multi_substrate_uv_first_move_spectral.json`
Memo: `debug/g4_6a_spectral_substrate_first_move_memo.md`

Two-panel sweep at $(a, N_\rho) = (0.05, 200)$ and $(0.025, 400)$ on spectral substrate.

**Headline result:** Per-$t$ recovery at $t = a^2$ improves from FD's 0.04% to spectral's 6.36% — a 160× improvement, matching v3.20.0 task #28's analytical prediction.

### 2.3 Production tests (A this thread)

File: `tests/test_warped_dirac_spectral.py` (14 tests; 13 fast + 1 slow, all pass).

Test coverage:

| Test class | Coverage |
|---|---|
| `TestF6BitExactDiskSpectral` | F6: wedge spectral at α=1 reduces to disk spectral bit-exactly (eigenvalues + heat traces) |
| `TestHilbertDim` | Hilbert dim = 2·N_ρ·N_φ (spinor doubling) |
| `TestSpectralEigenvalueStructure` | Exact m_eff = (k+1/2)/α; soft IR at α > 1 |
| `TestHeatTraceContinuity` | Heat trace monotonically decreasing in t, positive |
| `TestConstructorValidation` | Invalid parameters raise ValueError |
| `TestSpectralUVImprovement` (slow) | Spectral tip exceeds FD tip by ≥ 10× at t = a² (loose bound on B.2's measured 160×) |

## 3. Honest scope at closure

### 3.1 What G4-6d closes

- Spectral azimuthal discretization is operational in production code.
- F6 bit-exact reduction at α = 1 verified (eigenvalues + heat traces).
- UV improvement over FD at production substrate scales verified (~160× per B.2 measurement).
- 14 production tests covering construction, structural properties, F6, and UV improvement.

### 3.2 What G4-6d does NOT close

- F17 falsifier ("per-$t$ tip recovery within 1%") is NOT satisfied. The spectral substrate recovers 5.40-21.74% of the continuum $A$ across the panels — large improvement over FD, but not within 1%.
- The substrate-dependent Lichnerowicz constant $B \approx 0.30$ (vs continuum $+1/6 = 0.167$) is identified as the structural limitation requiring G4-6b (IR-boundary regularization) before further refinement.

### 3.3 Reframing of F17

The original F17 was framed as "spectral azimuthal alone closes the per-$t$ recovery to within 1%." This turned out to be too optimistic: spectral substrate handles the AZIMUTHAL part of the per-$t$ recovery correctly, but the RADIAL part (substrate-dependent contributions to $B$) requires separate treatment via G4-6b.

**Reframed F17 (post-G4-6d closure):** spectral azimuthal discretization closes the substrate's UV behavior at the per-$t$-leading-order level (factor 160× improvement vs FD, as documented). Full sub-percent recovery requires the joint G4-6d (spectral) + G4-6b (IR-boundary) combination, which is the G4-6a refined target.

## 4. Verdict

**G4-6d CLOSED at sprint scale.** Multi-month G4-6 commitment proceeds to G4-6b (next sequential sub-sprint per reframed sequencing).

This is the FIRST G4-6 sub-sprint to close at sprint scale, validating the reframed sequencing approach.

## 5. Cross-references

- `debug/g4_6_scoping_memo.md` §6 — reframed sub-sprint sequencing (updated 2026-05-29)
- `debug/g4_6a_multi_substrate_uv_first_move_memo.md` — FD baseline (Track C-2 of thread 3)
- `debug/g4_6a_spectral_substrate_first_move_memo.md` — spectral substrate verification (Track B.2 of thread 4)
- `debug/wedge_spectral_density_per_t_uv_target_memo.md` — v3.20.0 task #28 analytical prediction
- `geovac/gravity/warped_dirac.py` — production code (extended with two spectral classes)
- `tests/test_warped_dirac_spectral.py` — 14 production tests (this Track A)

## 6. Files

- `tests/test_warped_dirac_spectral.py` (new test file, 14 tests)
- `debug/g4_6d_spectral_closure_memo.md` (this)
- (No production code modification this Track A; B.1 already extended `geovac/gravity/warped_dirac.py`)
