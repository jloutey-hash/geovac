# LiH r12 t^2 three-body geminal, Step 1 (diagnostic, 2026-10-05)

**Verdict: STOP on the energy lever.** The reductions themselves work: the separable route matches Monte Carlo on every piece tested. The energy is the problem. Of the 30.1 mHa core cusp, about 0.3 mHa is t-dependent content that orbitals cannot represent. On a marriage-type reference, an r12 t^2 geminal gains about 0.03 to 1.2 mHa per core pair. That cannot close a 21.5 mHa gap.

Scripts: `debug/r12t2_step1_{hylleraas,hylleraas_b,separable,triangle,laplacian_check}.py`.
Data and logs: `debug/data/r12t2_step1_*.{json,log}`.
Threads were pinned to 1. Nothing in geovac/, tests/, papers/, CHANGELOG or CLAUDE.md was modified.

## Task A: what do the t-terms carry?

### A1. Reproducing the 30.1 mHa
`geovac.prolate_allelectron.helike_cusp.energy` (called by `debug/lih_additive_f12.py`) gives these Li+ values:
- E({t2, s}) = -7.246855
- E({u, t2, s, u2, ut2}) = -7.276999
- The difference, the core cusp, is **30.144 mHa**. The run takes 0.9 s.

My re-implementation on a configurable grid reproduces this to 3.6e-15. On a finer grid (96x72, RMAX 10) the cusp is 30.03 mHa.

Grid caveat: absolute energies on these grids sit 1-3 mHa above the analytic single-zeta E0. For example, Li+ E0 is -7.21984 on the grid against an analytic -7.22266. So the "left vs exact" columns are overstated by about that amount. Every *increment* reported below agrees to within 0.05 mHa between the two grids.

### A2. Term-by-term decomposition (linear Hylleraas, Phi0 = 1s(zeta)^2 with frozen zeta, fine grid)
Notation: u = r12, s = r1 + r2, t = r1 - r2.

| basis | Li+ rec (mHa) | He rec | H- rec |
|:--|--:|--:|--:|
| ref (1s^2) | 0 | 0 | 0 |
| {s} | 0.00 | 0.00 | 0.00 |
| {u} | 28.67 | 27.51 | 23.89 |
| {u,u2} | 31.11 | 29.46 | 24.74 |
| {u,u2,s} (t-free) | 46.51 | 43.97 | 36.43 |
| {u,u2,s,t2} | 56.65 | 55.49 | 53.89 |
| {u,u2,s,ut2} | 53.29 | 52.01 | 51.65 |
| rich {u,u2,s,t2,ut2} | 57.00 | 55.70 | 54.03 |
| total corr (E0 - exact) | 59.56 | 56.95 | 55.34 |

Increments (mHa), Li+ / He / H- (the H- column uses a wide grid with RMAX 22):
- **t2 on top of the t-free set: 10.14 / 11.51 / 17.45.**
- ut2 on top of the t-free set: 6.78 / 8.03 / 15.21.
- t2 + ut2 together on top of the t-free set: 10.49 / 11.72 / 17.60.
- **ut2 on top of {u,u2,s,t2}: 0.35 / 0.21 / 0.15.** This is the only t-content beyond orbitals; see A3.
- t2 on top of {u,u2,s,ut2}: 3.71 / 3.69 / 2.39.

He control:
- Rich 5-term: -2.90246.
- {u,t2,s}: -2.90118. The classic optimized-zeta 3-term Hylleraas is -2.90244. Ours is lower-flexibility because zeta is frozen at 27/16.
- Exact: -2.903724. The rich set recovers 98% of the 57 mHa correlation.

### A3. Which terms are genuinely beyond orbitals
- t^2 = r1^2 + r2^2 - 2 r1 r2 times Phi0 is a sum of products of one-electron radial functions. It is **radial correlation, representable by orbitals**, and a core-enriched CI already contains it.
- u^2 = r1^2 + r2^2 - 2 r1.r2 is also orbital-representable, through p-wave products.
- Only **odd** powers of u (u, u t^2, ...) lie outside every finite orbital CI.

So the "10 mHa of t-content" in the plain decomposition is radial correlation the M=16 reference already has. The genuinely three-body increment, ut2 on top of everything else, is **0.15-0.35 mHa per pair**.

This also means the 30.1 mHa "cusp" (rich minus {t2, s}) is not three-body content. It is u and u^2 plus re-optimisation of the radial coefficients in the presence of u; ut2 contributes 0.35 of it. Paper 19 (around L976-980) says the rich basis captured e-e-n content "that a pure u(r12) cannot". **These numbers do not support that attribution.** This is flagged for the PI; I made no paper edit.

### A4. Marriage analog (the faithful molecular translation)
`debug/r12t2_step1_hylleraas_b.py` sets up the analog as follows:
- Reference: Psi_orb = P(s, t^2, u^2) Phi0, an even-u polynomial of degree ≤ 2, 3 or 4 (4, 7 or 13 primitives). This is the orbital-CI analog.
- Marriage: the linear span {Psi_orb, g Psi_orb} with g = e^{-gam u}.
- The r12 t^2 lever adds t^2 g Psi_orb and/or u t^2 g Psi_orb.

Note: g Psi_orb **already contains** t^2 g with a fixed ratio. The new terms buy only the freedom to change that ratio.

Li+ results (fine grid, mHa):

| ref | gam | marriage | +t2g | +ut2g | +both | relax ref coeffs | t-terms after relax | 2nd u-shape {g,ug} |
|:--|:--:|--:|--:|--:|--:|--:|--:|--:|
| deg2 | 1.5 | 6.23 | 0.41 | 0.03 | 2.88 | 9.06 | 0.67 | 5.09 |
| deg3 | 1.5 | 2.58 | 0.14 | 0.09 | 2.03 | 7.32 | 0.45 | 4.21 |
| deg4 | 1.5 | 1.25 | 0.04 | 0.10 | 1.15 | 5.71 | 0.24 | 2.84 |
| deg4 | 3.0 | 2.69 | 0.31 | 0.03 | 0.84 | 4.00 | 0.30 | 1.84 |
| deg4 | 5.4 | 3.85 | 0.28 | 0.29 | 0.30 | 1.87 | 0.38 | 0.13 |

He gives the same pattern: t-terms 0.26-1.16 mHa, relaxation 1.6-5.8, after relaxation 0.23-0.50.
H- (valence proxy, part a, with a radial reference): t-terms 0.25-1.3 mHa.

Reading: **the t-terms carry 0.03-1.2 mHa per core pair** on realistic references (up to 2.9 on the poorest degree-2 reference). Once the reference coefficients are relaxed, that falls to 0.2-0.7 mHa. "+both" is larger than the sum of its parts because the pair partly mimics relaxation of the reference.

### A5. Molecular range
The molecule has:
- the Li core pair, giving ≲ 0.1-1.2 mHa;
- the valence pair, ≲ 0.1-0.6 mHa;
- core-valence pairs, where t is measured from the wrong nucleus for one electron, so smaller still.

**Plausible gain from an r12 t^2 geminal: about 0.3-3 mHa (central about 1 mHa) against the 21.5 mHa gap**, i.e. 1-15%.

Caveats on that range:
- Transfer from atom to molecule is by analogy only. The isolated-pair test has no core-valence coupling.
- The upper end assumes the reference is as poor as the degree-2 analog.
- The molecular geminal already contains the u-only content and, through the product with Psi_CI, the fixed-ratio t^2-times-u content. The A4 numbers already account for that.
- None of this is variational information about the molecule.

### A6. What the remaining 21.5 mHa more plausibly is (hypotheses, ranked by the analog)
1. **Reference coefficients frozen in the 2x2.** The marriage (and the VMC-Padé twin) fix the CI coefficients. Relaxing them alongside the geminal is worth 1.9-9 mHa per pair in the analog, and the effect is largest at the compromise gam = 1.5. In QMC practice, CI-coefficient re-optimisation in the presence of the Jastrow is the standard companion. It is also what e-e-n Jastrow terms partly emulate, which reconciles this with the "e-e-n" framing. The molecular test is cheap: either an (N_det+1)-dimensional {Phi_I} ∪ {(F - F̄) Psi_CI} problem (one f-edge h-type elements per determinant), or a VMC linear-method run with CI coefficients unfrozen.
2. **Orbital-basis incompleteness of the reference.** The M>16 ceiling is -8.032, and a linear multiplier only partly repairs it. In the analog, the residual after marriage drops 7.5 → 1.5 mHa as the reference degree rises.
3. **Geminal range/shape (one gam for both core and valence).** In the analog a second u-shape gives 0-5 mHa, but the **molecular** VMC multi-geminal test already capped this at +1.4 mHa (ledger, 2026-09-26). So this is not the main lever.

## Task B: separability (does the RI-free machinery carry u t^2 g?)

### B1. Algebra
With t = r_1A - r_2A and k(u) = u g(u):

f(1,2) = k(u) [ r_1A^2 · 1 + 1 · r_2A^2 - 2 r_1A · r_2A ]

This is an exact rank-3 sum of (one-electron weight) × (pure two-body kernel) × (one-electron weight). Consequences:
- The f-dressing field is Psi^f_P = r_A^2 Psi^k_P + Psi^k_{r_A^2 P} - 2 r_A Psi^k_{r_A P}.
- Products stay separable:
  - f^2 = k^2 t^4 has rank 5, using t^{2n} = Σ_m C(2n,m)(-1)^m r_1A^{2n-m} r_2A^m.
  - f/u = g t^2 has rank 3.
  - A shared-vertex product f_ij f_ik has rank 9, with the weights multiplied at the shared vertex.
- With g = e^{-gam u}, the kernels needed are u e^{-gam u} (existing linexp), u^2 e^{-2 gam u}, e^{-gam u}, and the Yukawa e^{-gam u}/u.

Effect on each energy piece:
- **F̄, sigma^2, h_Vne, h_Vee, g_Vne:** reweighted densities only.
- **Triangle and bridge (g_Vee):** each (c,d) weight pair is the existing reduction, applied to reweighted μ (shared vertex) and edge densities. The r_A weights are axisymmetric (m = 0), so the m-structure of the π reference is unchanged. The η-polynomial degree rises by ≤ 4, so the termination bound becomes 2(l_bridge + l_leaf) with l_leaf raised by up to 4. That is finite, so the reduction stays RI-free; it needs no new structure, only a larger LMAX budget.
- **Kinetic h_T Part B and g_T (gT2+gT3):** these need ∇_i^2 f. It is separable:

  ∇_1^2 f = t^2 (k'' + 2k'/u) + 2k' t (r_1A^2 - r_2A^2 + u^2)/(r_1A u) + k (6 - 4 r_2A/r_1A)

  The r̂_1A·û_12 cross term closes through the law of cosines in the triangle (1, 2, A). The weights are {r^n, r_2A^n/r_1A} and the kernels {∇^2 k, k'/u, u k', k}, i.e. Yukawa, Slater, linexp and u^2-Slater. Verified symbolically to 2e-14 (`r12t2_step1_laplacian_check.py`). The new 1/r_A weight is integrable; the prolate Jacobian cancels it.
- **h_T Part A and gT1:** gradient densities G_pq times f-dressing fields, reweighted as above.

**No piece needs a reduction that does not exist.** What does *not* carry over:
- The 2-MO isotropic closed-form shortcuts (`f_dress_iso`, `yukawa_pot_iso`, `_hartree_1s` on reweighted densities): r_A^n times a B-centred or ab density is not isotropic. The general Neumann / grid paths and the Coulomb-minus-smooth Yukawa trick, already used by the marriage engine, cover these.
- **Cost:** about rank^2 ≈ 9-25× on g_Vee.

### B2. Numeric check on the 2-MO reference
The separable route feeds reweighted fields into the package's own `hVee.cov_FA_FB` and `triangle.triangle_raw`. It is compared against brute-force MC with exact i.i.d. rejection sampling of |Phi0|^2 (proposal ½[A²B² + B²A²]; no Metropolis bias). f = u t^2 e^{-gam u}; 1e8 configurations for the first table and 4e7 for the triangles.

| quantity | separable | MC ± err | rel dev | n σ |
|:--|--:|--:|--:|--:|
| F̄ (gam 1.5) | 0.803217 | 0.803189 ± 2.5e-5 | 3.5e-5 | 1.1 |
| sigma^2 (incl. 3-body share-vertex + 4-body disjoint) | 0.066405 | 0.066402 ± 8.9e-6 | 5.3e-5 | 0.4 |
| h_Vne | -0.936750 | -0.936783 ± 3.5e-4 | 3.5e-5 | 0.1 |
| h_Vee (incl. ⟨f12/r13⟩, ⟨f12/r34⟩), exact Neumann | 0.116757 | 0.116698 ± 4.4e-5 | 5.1e-4 | 1.3 |
| triangle T1 = ⟨f13 f23/r12⟩, exact Neumann | 0.0035608 | 0.0035627 ± 1.1e-6 | 5.6e-4 | 1.8 |
| triangle T2 = ⟨f12 f23/r13⟩, exact Neumann | 0.0230677 | 0.0230776 ± 9.0e-6 | 4.3e-4 | 1.1 |
| F̄ / sigma^2 (gam 0.5, long-range stress) | 20.33970 / 43.4613 | 20.34077 / 43.4667 | 5e-5 / 1.2e-4 | 1.6 / 1.1 |

Positive control: the plain e^{-0.5u} geminal through the same wrapper reproduces the banked F̄ = 1.8807 and sigma^2 = 0.1372 (vs MC: 0.9σ). The control's triangles pass at 0.4σ / 0.1σ.

All rows pass the target (≤ 1e-3 relative, within 3σ). φ-quadrature is converged: NPHI 28 → 64 leaves F̄ and sigma^2 unchanged to 6 digits.

**Side finding: legacy Neumann path.** On the legacy (default) Neumann path, T2 is off by 2.6e-3 (9.5σ) for the **unchanged** plain geminal and by 4.1e-3 for u t^2. `USE_EXACT_NEUMANN = True` removes the error. This is the known O(NXI^-2) kink error (Phase 0b), now measured on the inter-spin Coulomb-edge triangle of the 2-MO path at high MC precision. Not a separability defect. The original gate's "T2 0.4σ" only passed because its Metropolis error bars were loose.

The 4-body bridge was not run separately. Its inputs, reweighted leaf and bridge densities, are the same objects validated above, and the disjoint 4-body ⟨f13 f24⟩ inside sigma^2 passes.

## Dead-end check (docs/failed_approaches_ledger.md)
- No row covers r12 t^2 / e-e-n geminals specifically.
- Related rows:
  - "Per-pair γ / second correlation length" (2026-08-26): u-shape flexibility < 0.6 mHa.
  - "Multi-geminal / F^2" (2026-09-26): +1.4 mHa saturation in the molecule.
  - The memory DO-NOT: "a Gaussian e-e-n term only DISTORTS".
- This result belongs to the same family: correlation-factor flexibility beyond the Kato u-shape is not where the LiH residual lives. It is a candidate new ledger row (PI/PM to log).

## Verdict
**STOP** for the r12 t^2 build as the route to -8.062:
- The genuine (non-orbital) t-content is 0.15-0.35 mHa per pair in the linear decomposition, and 0.03-1.2 mHa per core pair in the marriage analog.
- The molecular range is about 0.3-3 mHa against 21.5 mHa.

Separability itself is **GO**. It is exact, finite rank, and needs no new reduction. Every tested piece matches MC to ≤ 5.6e-4, on the exact-Neumann path. So an e-e-n geminal could be added cheaply if ever wanted.

Recommended next diagnostic (cheap): unfreeze the CI coefficients alongside the geminal in the molecule. This is the leading candidate for the residual (analog: 2-9 mHa per pair). Test it first by VMC linear method on the cached M=16 wavefunction, then analytically as an (N_det+1)-dimensional problem with h-type elements only.
