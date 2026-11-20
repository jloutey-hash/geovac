# Sturmian Basis Attempts for LiH — Complete Arc (v0.9.20 to v0.9.34)

**Date:** March 11, 2026
**Status:** CLOSED — Structural theorem proven, arc complete
**Conclusion:** All Sturmian approaches failed for LiH. LCAO architecture vindicated.

---

## Summary

29 diagnostic versions (v0.9.8 through v0.9.36) were tested for LiH binding
in the GeoVac framework. The Sturmian basis arc (v0.9.20-v0.9.34) explored
every reasonable variant of the shared-p0 Sturmian construction. All failed.
The failure is not a bug but a structural theorem: the Sturmian Hamiltonian
is proportional to the overlap matrix, making eigenvalues geometry-independent.

---

## Version-by-Version Arc

### v0.9.20-v0.9.25: Atom-Centered Sturmian Variants
- **v0.9.20**: Basic Sturmian diagonal (eps_k = p0^2/2 - beta_k*Z*p0/n).
  H 1s diagonal POSITIVE (+1.84 Ha). LiH UNBOUND.
- **v0.9.21**: Added SW D-matrix off-diagonal. Still UNBOUND (H 1s deficit 2.35 Ha).
- **v0.9.22-25**: Various cross-nuclear corrections, angular weighting,
  shell-radius scaling. All UNBOUND or overbinding. The fundamental issue:
  shared p0 ~ 3.17 exceeds 2*Z_B = 2 threshold for H binding.

### v0.9.26-v0.9.27: SO(4) CG Cross-Center V_ee
- SO(4) Wigner D-matrix rotation of Sturmian Slater F0 integrals.
- J_cross(1s_A, 1s_B) = 1.98 Ha EXCEEDS same-center J = 1.88 Ha (artifact).
- Cannot rescue single-p0: V_ee adds repulsion, worsens deficit.
- **CONCLUSION**: Single-S3 shared-p0 fundamentally incompatible with
  heteronuclear molecules.

### v0.9.28: Molecular Beta from Bond Sphere Node Weights
- beta_k = (self_nuclear + cross_nuclear) / self_nuclear.
- beta(H 1s) = 1.312 (needed 1.582 for binding). 46% improvement but
  insufficient. LiH UNBOUND.

### v0.9.29: Prolate Spheroidal Molecular Sturmian Betas
- Complete rewrite with prolate spheroidal coordinate separation.
- H2+ verification: E = -1.096 Ha (0.6% error) — validates approach.
- LiH betas: 1sigma=1.033, 2sigma=1.950, 1pi=1.976.
- Problem: molecular betas are for MOs, not atom-centered orbitals.
  Mapping is ad hoc.

### v0.9.30: MO-to-Atom-Center Beta Projection
- Dominant-overlap projection of MO betas onto atomic orbitals.
- FCI: E = -42.5 Ha (UNPHYSICAL, 5.3x error).
- Root cause: beta scale mismatch (prolate spheroidal vs atomic Sturmian)
  + wrong MO assignment.

### v0.9.31: MO Sturmian FCI (First Success in Energy Scale)
- Full FCI in molecular Sturmian basis with Lowdin orthogonalization.
- E = -10.07 Ha (25% error). BOUND but overbinding.
- Without orthogonalization: E = -13.15 Ha (63% error).
- Population-weighted V_ee underestimates repulsion.

### v0.9.32: Canonical Orthogonalization + Exact J Integrals
- Poisson solve for direct Coulomb integrals on prolate spheroidal grid.
- E = -7.131 Ha (11.6% error). UNBOUND.
- Error halved vs v0.9.31. PES shape correct (minimum near R=3.5).
- Root cause: excess J repulsion from compact shared-p0 orbitals.

### v0.9.33: Exact Exchange K Integrals
- K via Poisson solve. K_ohno was OVERESTIMATED (0.269 vs 0.174 for 1sigma-2sigma).
- Removing accidental compensation worsens energy: E = -6.796 Ha (15.8% error).
- UNBOUND. p0 feedback does not self-correct compactness.

### v0.9.34: Dual-p0 MO Sturmian Basis
- Two MO sets: Li-scale (nmax=3, p0_Li) and H-scale (nmax=2, p0_H).
- Combined 14 spatial MOs, 28 spinorbs, 20,475 SDs.
- E ~ -10.1 Ha (25% error). BOUND but same overbinding as v0.9.31.
- p0_Li drifts from 3.845 to ~2.25, p0_H drifts from 1.0 to ~1.74.
  The two scales CONVERGE to a single geometric mean.
- **ROOT CAUSE**: Prolate spheroidal MOs are inherently delocalized.
  Cannot produce atom-specific scales. Dual-p0 eliminates itself.

---

## Structural Theorem (Proven, v0.9.37)

**Theorem**: In the molecular Sturmian basis with shared momentum scale p0,
the one-electron Hamiltonian matrix satisfies:
  H[i,j] = (p0^2/2)(1/beta_j - 1) * S[i,j]
The Hamiltonian is proportional to the overlap matrix. Consequently, the
generalized eigenvalues E_k = (p0^2/2)(1/beta_k - 1) are independent of
the bond angle gamma and hence independent of internuclear distance R.

**Corollary 1**: Binding requires beta_k(R) — the full R-dependent prolate
spheroidal eigenvalue spectrum. This requires solving the two-center
Schrodinger PDE at every R, reintroducing continuous geometry.

**Corollary 2**: For heteronuclear diatomics, no shared p0 exists that
simultaneously represents both atomic length scales. The dual-p0 extension
fails because molecular Sturmian functions are inherently delocalized.

**Numerical evidence**: v0.9.34 bond sphere Hamiltonian analysis confirms
E_total(R) is a pure 1/R curve with no bound state when beta is fixed.

---

## Additional Negative Results (v0.9.35-v0.9.36)

### v0.9.35: Cross-Atom V_ee All-L Extension
- Extended to all (n,l) pairs: 392 J entries (was 18 s-only).
- dE(all_l - s_only) = +0.449 mHa — negligible.
- Ground state is s-dominated (1s^2 2s^2). NOT the R_eq mechanism.

### v0.9.36: Orbital Exponent Relaxation
- Per-atom zeta scaling. zeta_B*(R=3.015) ~ 1.30.
- dE_relax nearly R-independent: -0.272 (R=2.5), -0.279 (R=3.015) Ha.
- R_eq unchanged at ~2.5. NOT the R_eq mechanism.

### Harmonic Phase Locking (v0.9.37 diagnostic)
- No critical points of D-matrix elements map to R_eq ~ 3.015 bohr.
- R_eq is an emergent variational balance, not a geometric fixed point.

---

## Final Assessment

The Sturmian bond sphere is geometrically elegant but computationally
intractable for heteronuclear diatomics without solving the continuous PDE
at every R. The atom-centered LCAO architecture (MolecularLatticeIndex)
preserves atom-specific length scales and achieves D_e_CP = 0.093 Ha
(1.0% error) at nmax=3. The LCAO path is the correct discrete framework
for GeoVac molecular calculations.
