# Sprint G4-4b scoping — Variable-warp Dirac on the discrete cigar substrate

**Date:** 2026-05-28
**Path:** Gravity arc, sub-sprint of the multi-month G4-4 commitment. **G4-4a is closed at the load-bearing falsifier level** (F1, F2, F3 all positive-verified in weeks 1-3 of the present session). G4-4b is the second sub-sprint of G4-4, extending the constant-warp Dirac to the **variable-warp** geometry that realizes the actual Euclidean Schwarzschild cigar.
**Verdict:** **POSITIVE-SCOPING-G4-4b.** Variable-warp architecture is well-defined as a structurally additive extension on top of G4-4a; four load-bearing falsifiers F4-tip, F5-asym, F6-Riemannian-limit, F7-factorization-loss named with operational test protocols; sub-sprint sequence G4-4b a/b/c/d sized at 4-8 weeks. Multi-month commitment scoping unchanged from G4-4 (~5-8 months for G4-4 end-to-end).

---

## §1. Context

The Euclidean Schwarzschild cigar's near-horizon geometry is
$$ds^2 = d\rho^2 + \rho^2\,d\phi^2 + r(\rho)^2\,d\Omega_2^2$$

G4-4a closed the **constant-warp** case $r(\rho) = r_h$ at sprint scale (this session): F1 (factorization $K_{\rm cigar} = K_{D^2} \cdot K_{S^2}$, bit-exact), F2 (chirality grading, operator-level), F3 (continuum Weyl-Selberg, 0.4% at sweet spot), lowest-mode $|\lambda_{\min}| \to \pi/R$ at 0.5%. Architecture in `geovac/gravity/warped_dirac.py`; 54 tests passing.

**The constant-warp case captures the near-horizon geometry at exact equality $r = r_h$.** It does NOT capture:
- The asymptotic-flat region where $r(\rho) \to \rho$ at large $\rho$ (Schwarzschild far field)
- The smooth-tip regularity at $\rho \to 0$ (no conical defect when $\alpha = 1$)
- The spin-connection coupling that arises from $r'(\rho) \ne 0$

These are the load-bearing physics of the cigar — *they are what makes it a Schwarzschild cigar rather than a $D^2 \times S^2$ product*. G4-4b builds the Dirac operator that respects them.

## §2. Target physics

### §2.1 Smooth-tip warp profile

The canonical smooth-tip warp profile is
$$r(\rho) = r_h\sqrt{1 + (\rho/r_h)^2}$$
with:
- **Tip behavior** ($\rho \ll r_h$): $r(\rho) \to r_h$, smooth (no conical defect at $\alpha = 1$). Derivative $r'(\rho) \to \rho/r_h$.
- **Asymptotic behavior** ($\rho \gg r_h$): $r(\rho) \to \rho$. Derivative $r'(\rho) \to 1$ (matches flat Schwarzschild far field).

This is the **G4-3b substrate** (already implemented in `debug/g4_3b_variable_warp.py`) extended to the spinor sector.

### §2.2 Warped product Dirac with spin connection

The 4D Dirac operator on the warped product $D^2 \times_w S^2$ with warp factor $r(\rho)$ (Camporesi 1996, standard form):
$$D_{\rm cigar}^{\rm var} = D_{D^2} \otimes I_{S^2} + \gamma^5_{D^2} \otimes \frac{D_{S^2}}{r(\rho)} + \frac{r'(\rho)}{r(\rho)}\,\gamma^\rho \otimes I_{S^2}$$

The **warp-derivative term** $\frac{r'(\rho)}{r(\rho)}\,\gamma^\rho$ is the spin connection on the warped product. It arises from the non-zero extrinsic curvature of the $S^2_{r(\rho)}$ slices in the cigar.

For smooth-tip warp:
$$\frac{r'(\rho)}{r(\rho)} = \frac{\rho/r_h^2}{1 + (\rho/r_h)^2} \cdot \frac{r_h^2}{r_h^2 + \rho^2} \cdot r_h = \frac{\rho}{\rho^2 + r_h^2}$$

- $\rho \to 0$: $\rho/r_h^2$ (linear, tip-regular)
- $\rho \to \infty$: $1/\rho$ (asymptotic-free, matches flat Schwarzschild)

### §2.3 Squared variable-warp Dirac

The squared operator picks up structural corrections beyond the constant-warp form:
$$(D_{\rm cigar}^{\rm var})^2 = (D_{D^2}^{\rm var})^2 \otimes I + I \otimes \frac{D_{S^2}^2}{r(\rho)^2} + \text{cross terms involving } r'(\rho)$$

where $D_{D^2}^{\rm var}$ inherits the warp-derivative term in the radial sector. **The factorization $K_{\rm cigar} = K_{D^2} \cdot K_{S^2}$ that closed F1 in G4-4a is BROKEN at variable warp**: cross terms in $D^2$ couple the disk and sphere sectors via $r'(\rho)$. This is a *structural feature* — it's what distinguishes the genuine cigar from the product geometry.

## §3. Inputs from G4-4a (what's already done)

All G4-4a infrastructure transfers cleanly:

| Piece | Source | Reuse |
|---|---|---|
| 2D Cl(2,0) gamma matrices | `geovac/gravity/warped_dirac.py` | Same $\gamma^1, \gamma^2, \gamma^5$ |
| `DiscreteDiskDirac` (squared, anti-periodic phi) | same | Constant-warp limit baseline |
| `DiscreteDirac2D` (explicit linear D) | same | Operator-level substrate |
| `S2DiracSpectrum` (Camporesi-Higuchi) | same | $|\lambda^{S^2}| = (n+1)/r$ — now $r = r(\rho)$ |
| Hermitian polar Laplacian | G4-3a-cleanup | Same $L_k$ blocks |
| Variable warp profile | G4-3b | Same $r(\rho)$ smooth-tip |
| F1, F2, F3 falsifier infrastructure | G4-4a | Adapted for variable case |

**No re-derivation of G4-4a results.** G4-4b is additive.

## §4. Architectural extension for variable warp

### §4.1 New computational object

`VariableWarpDirac` class extending the G4-4a `WarpedDiracConstant` architecture. Key difference: $r$ is no longer a scalar but a discrete function $r(\rho_k)$ evaluated at each radial site. The S² Dirac sub-block is now site-indexed:
$$\frac{D_{S^2}}{r(\rho)} \to \text{spectrum } \pm\frac{n+1}{r(\rho_k)} \text{ at radial site } k$$

This breaks the tensor-product structure: each radial site sees a *different* effective S² mass.

### §4.2 Three architectural levels

**Level 1 — Spectrum-level variable warp.** Compute the full cigar squared-Dirac spectrum by summing over radial sites:
$$K_{\rm cigar}^{\rm var}(t) = \sum_{\text{disk modes}} \sum_{S^2 \text{ modes}} \sum_k w_k \cdot e^{-t[\lambda_{\rm disk}^k + (n+1)^2/r(\rho_k)^2]}$$
where $w_k$ is the radial-site weight from the eigenvector amplitude. Cheap, but requires careful weighting.

**Level 2 — Operator-level variable warp.** Build the explicit $D_{\rm cigar}^{\rm var}$ as a sparse matrix on the full spinor Hilbert space, including the warp-derivative spin connection term. Expensive but operator-faithful.

**Level 3 — Effective-radial reduction.** Project the warp-dependent S² sector onto the lowest few $S^2$ modes and treat them as effective massive radial fields. Cheapest, but assumes a specific spectral truncation.

**Recommendation for G4-4b first move (week 1):** start with Level 1 spectrum-level, validate against G4-4a constant-warp at $r(\rho) = r_h$ (load-bearing reduction), then promote to Level 2 in week 2.

## §5. Load-bearing falsifiers

### F4-tip: Tip-regularity of warp-derivative term

$$r'(\rho)/r(\rho) = \rho/(\rho^2 + r_h^2) \to \rho/r_h^2 \text{ as } \rho \to 0$$

**Test:** the discrete warp-derivative array is finite at the smallest radial site (no $1/\rho$ singularity at the apex). Verify $r'(\rho_1)/r(\rho_1) \approx a/r_h^2$ at $\rho_1 = a$, small.

**Verdict expected:** PASS bit-exact (analytical identity at the discrete site).

### F5-asym: Asymptotic-free recovery at large $\rho$

For $\rho \gg r_h$, the warp $r(\rho) \to \rho$, and the variable-warp Dirac approaches the flat-space asymptotic form. The heat trace in the asymptotic region should match the Weyl law for a flat $S^2$ with radius $\rho$, integrated over the disk.

**Test:** at large $\rho$ truncation $\rho_{\rm IR} \gg r_h$, the asymptotic contribution to $K(t)$ approaches the analytical Weyl prediction for the spherical Schwarzschild radial geometry. Specifically, the $S^2$ sector at site $k$ contributes $\sim r(\rho_k)^2 / t$ to the leading Weyl term; summed over $k$ this should give the proper 4D Schwarzschild Weyl law $\sim \int \rho\, d\rho\, r(\rho)^2 / (4\pi t)^2$.

**Verdict expected:** PASS within a few % at the asymptotic-free panel ($R \gg r_h$).

### F6-Riemannian-limit: Reduction to constant warp at $r(\rho) = r_h$

If the warp profile is forced to be constant (e.g., set the implementation to $r(\rho) = r_h$ identically), the variable-warp Dirac must reduce bit-exact to the G4-4a constant-warp form.

**Test:** spectrum and heat trace of `VariableWarpDirac(disk, sphere, r_profile=constant r_h)` match `WarpedDiracConstant(disk, sphere)` at machine precision.

**Verdict expected:** PASS bit-exact (operator-level identity by construction).

### F7-factorization-loss: Quantification of factorization-loss

The constant-warp factorization $K_{\rm cigar} = K_{D^2} \cdot K_{S^2}$ is broken at variable warp. Quantify:
$$\Delta_{\rm fact}(t) := \left| K_{\rm cigar}^{\rm var}(t) - K_{D^2}^{\rm var}(t) \cdot K_{S^2_{r_h}}(t) \right|$$
and verify it scales as $r'(\rho)^2 \cdot$ (something) — the *structural prediction* is that the factorization-loss is controlled by the warp derivative squared, vanishing as $r'(\rho) \to 0$ (constant warp limit).

**Test:** monotonic in warp variation. At small warp variation (smooth-tip with $\rho \ll r_h$), $\Delta_{\rm fact}$ small; at large warp variation (large $\rho$, asymptotic free), $\Delta_{\rm fact}$ grows.

**Verdict expected:** PARTIAL — sign and scaling correct; quantitative form may require fitting.

## §6. G4-4b sub-sprint sequence

Four sub-sprints sized per the G4-4a cadence (week ≈ focused PM session + sub-agent dispatch):

| Sub-sprint | Scope | Effort |
|---|---|---|
| **G4-4b-a** (first move) | Spectrum-level variable warp; F6 Riemannian-limit recovery; F4 tip-regularity at the warp-derivative | 1 week |
| **G4-4b-b** | F7 factorization-loss quantification; structural pattern of $\Delta_{\rm fact}(t)$ vs warp variation | 1-2 weeks |
| **G4-4b-c** | F5 asymptotic-free recovery at large $\rho$; comparison to Schwarzschild Weyl law | 1-2 weeks |
| **G4-4b-d** | Explicit operator-level variable warp Dirac (Level 2); spin-connection term construction | 2-3 weeks |

**Total G4-4b commitment: 4-8 weeks**, consistent with the G4-4 scoping estimate.

## §7. First-move plan (G4-4b-a)

### §7.1 Module structure

Extension of `geovac/gravity/warped_dirac.py` (~150 new lines) with `VariableWarpDirac` class:

```
- VariableWarpDirac(disk, sphere, warp_profile, r_h)
  - warp_profile: callable rho -> r(rho), or precomputed array
  - radial_S2_blocks: list of S^2 squared eigenvalue arrays, one per radial site
  - squared_eigenvalues(): outer sum, but weighted by radial eigenvector amplitude at each site
  - heat_trace_variable(t): direct sum

- Smooth-tip warp function (matches G4-3b convention)
- Constant warp companion for F6 Riemannian-limit check

- Falsifier functions:
  - verify_F4_tip_regular(...)
  - verify_F5_asymptotic_free(...)
  - verify_F6_Riemannian_limit(...) — load-bearing reduction
  - verify_F7_factorization_loss(...)
```

### §7.2 Test architecture

New tests in `tests/test_warped_dirac.py` (~15 new tests + 2 slow):
- VariableWarpDirac construction and parameter validation
- F4-tip at small panel
- F6-Riemannian-limit bit-exact at constant $r(\rho)$
- F7 monotonic in warp variation
- Heat trace at variable warp positive and decreasing in $t$
- Variable-warp $S^2$ block construction (site-indexed eigenvalues)
- Warp profile evaluation correctness (smooth-tip + asymptotic)

### §7.3 Driver and memo

- `debug/g4_4b_a_first_move.py`: driver
- `debug/g4_4b_a_first_move_memo.md`: closure memo (target verdict POSITIVE-G4-4b-a-VERIFIED with F4, F6 passing)
- `debug/data/g4_4b_a_first_move.json`: numerical panel

### §7.4 Sequencing for G4-4b-a (one-week breakdown)

- **Day 1-2**: warp profile array construction; F4 tip-regularity test; tests for smooth-tip behavior
- **Day 3-4**: `VariableWarpDirac` class; spectrum-level variable warp; F6 Riemannian-limit at constant $r(\rho) = r_h$ bit-exact
- **Day 5-6**: heat-trace computation at variable warp; first-pass F7 factorization-loss diagnostic
- **Day 7**: closure memo, sub-sprint review, hand-off to G4-4b-b

## §8. Multi-month commitment scoping (unchanged)

G4-4b sits in the multi-month G4-4 commitment that opens the full discrete-substrate $S_{\BH}$ derivation. Sequence:

| Sub-sprint | Effort | Cumulative |
|---|---|---|
| G4-4 scoping (T3, done) | sprint | done |
| G4-4a (this session: weeks 1-3 done, F1/F2/F3 closed) | 4-8 weeks total | partial done |
| **G4-4b** (this scoping) | 4-8 weeks | in flight |
| G4-4c (conical defect spinor) | 4-8 weeks | queued |
| G4-4d (Seeley-DeWitt extraction) | 4-8 weeks | queued |
| G4-4e (BC sectors) | 2-4 weeks | queued |
| G4-4f (replica preparation) | 4-8 weeks | queued |
| **G4-4 total** | **5-8 months** | |
| G4-5 (discrete replica method) | 2-4 months | queued |
| G4-6 (full $S_{\BH}$ derivation) | 2-4 months | queued |

End-to-end G4-4 through G4-6: **9-16 months**. Unchanged from the G4-4 scoping memo (T3).

## §9. Risk-tier classification

**Risk tier 1 (low) — sprint-scale doable:**
- F4-tip: analytical identity, sprint-scale bit-exact
- F6-Riemannian-limit: load-bearing reduction, sprint-scale bit-exact (operator-level identity)
- VariableWarpDirac module construction: mechanical extension of G4-4a infrastructure

**Risk tier 2 (medium) — requires sustained attention:**
- F7-factorization-loss structural form: requires fitting or analytical derivation; structurally important
- F5-asymptotic-free at large $\rho$: needs careful Weyl-law comparison with the asymptotic Schwarzschild form
- Operator-level Level 2 construction with spin-connection terms

**Risk tier 3 (higher) — multi-month:**
- Variable-warp Seeley-DeWitt extraction (G4-4d, depends on G4-4b)
- Replica-method preparation at variable warp (G4-4f)
- Joint UV refinement with the radial site dependence (compute scaling)

## §10. Honest scope

This is a **scoping memo**. It does NOT contain new computational results. It:
- Confirms G4-4a outputs (F1, F2, F3 closed) are sufficient input
- Names the architectural extension G4-4b needs (`VariableWarpDirac` class + spin-connection terms)
- Names four load-bearing falsifiers F4, F5, F6, F7 with operational test protocols
- Names sub-sprint sequence G4-4b a/b/c/d with effort estimates
- Identifies risk-tier classification

Per the G4-4 scoping discipline (T3), the multi-month G4-4 commitment is **structurally additive on top of G4-4a's closure**. No re-derivation of G4-4a results.

## §11. Strategic positioning

The variable-warp Dirac is the **first substantive structural advance** beyond the constant-warp product geometry. Constant-warp $D^2 \times S^2$ at $r_h$ is mathematically a tensor product — interesting but structurally indistinguishable from any product Riemannian geometry. The cigar's actual content lives in **the variable-warp coupling** (F7 factorization-loss is its structural signature; F5 asymptotic-free recovery validates it captures the actual Schwarzschild far field; F4 tip-regularity validates it captures the actual smooth horizon).

G4-4b is the first sub-sprint where the framework produces *genuinely cigar-specific* content — content that distinguishes the Schwarzschild cigar from a $D^2 \times S^2$ product.

## §12. Cross-references

- **G4-4 scoping (T3)**: `debug/g4_4_warped_dirac_scoping_memo.md` (architectural blueprint)
- **G4-4a closure (this session, three weeks)**: F1, F2, F3 closed at sprint scale
  - Week 1: `debug/g4_4a_first_move_memo.md`
  - Week 2: `debug/g4_4a_week2_explicit_dirac_memo.md`
  - Week 3: `debug/g4_4a_week3_quantitative_f3_memo.md`
- **G4-3b variable warp scalar**: `debug/g4_3b_variable_warp_memo.md` (scalar substrate, transfers to spinor)
- **G4-3a-cleanup Hermitian polar Laplacian**: substrate for radial blocks
- **Camporesi 1996**: warped-product spin connection
- **Paper 23 §V**: Camporesi-Higuchi $S^2$ spinor spectrum

## §13. Files

- `debug/g4_4b_variable_warp_scoping_memo.md` — this memo
- (Future, when G4-4b-a launches): `geovac/gravity/warped_dirac.py` extended with `VariableWarpDirac`, `tests/test_warped_dirac.py` extended, `debug/g4_4b_a_first_move.py` + JSON + memo
