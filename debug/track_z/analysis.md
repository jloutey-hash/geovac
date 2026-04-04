# Track Z: Geometric Elevation of Level 4 Piecewise Structure

## Problem Statement

The Level 4 angular Hamiltonian H_ang(rho) contains nuclear coupling V_nuc
with piecewise structure from `f_k(s, rho) = (min(s,rho)/max(s,rho))^k / max(s,rho)`.
The min/max creates a derivative discontinuity (kink) at s = rho, where
s = cos(alpha) or sin(alpha). Track S (v2.0.13) proved no global algebraic
curve P(rho, mu) = 0 exists -- linearity residuals 22.5% (l_max=1), 58.6%
(l_max=2).

Question: can geometric elevation to a higher-dimensional space resolve
this piecewise structure, analogous to how Fock's S3 projection resolves
the 1/r singularity of hydrogen?

## Avenue 1: Algebraic Blow-up

### Construction

The piecewise function min(s, rho) = (s + rho - |s - rho|) / 2 has a kink
at s = rho due to |s - rho|. In algebraic geometry, the absolute value
|x| = sqrt(x^2) is resolved by the real blow-up: introduce sigma^2 = (s-rho)^2,
so sigma = s - rho on one sheet and sigma = rho - s on the other.

On the blown-up space (alpha, sigma), define:
  - Sheet +: sigma = s - rho >= 0  (i.e., s >= rho)
    min(s,rho) = rho, max(s,rho) = s
    f_k = rho^k / s^{k+1}
  - Sheet -: sigma = rho - s >= 0  (i.e., s <= rho)
    min(s,rho) = s, max(s,rho) = rho
    f_k = s^k / rho^{k+1}

Within each sheet, f_k is a smooth rational function of (s, rho).

### Assessment: NEGATIVE

The blow-up resolves the kink but does NOT produce a global algebraic curve.
The obstruction is structural:

1. **Two electrons, four boundary surfaces.** Electron 1 has s1 = cos(alpha),
   creating a boundary at cos(alpha) = rho. Electron 2 has s2 = sin(alpha),
   creating a boundary at sin(alpha) = rho. For general l_max, ALL multipole
   terms k = 0, 1, ..., l1+l2 contribute, each with its own min/max
   switching at the same boundaries. The blow-up produces 2^2 = 4 sheets
   (one per combination of sign choices for the two electrons).

2. **Sheet boundaries are rho-dependent.** The boundary cos(alpha) = rho
   moves in alpha as rho changes. At each rho, the alpha-integration domain
   [0, pi/2] is split at alpha_1 = arccos(rho) and alpha_2 = arcsin(rho).
   After spectral projection (integrating against Jacobi polynomials),
   the matrix elements become integrals with rho-dependent limits. These
   produce transcendental functions of rho (inverse trig), not polynomials.

3. **Characteristic polynomial obstruction.** Even within a single sheet,
   P(rho, mu) involves integrals like:
     integral_0^{arccos(rho)} phi_k(alpha) * [rho^k / cos(alpha)^{k+1}] * phi_j(alpha) d(alpha)
   The upper limit arccos(rho) makes this a transcendental function of rho.
   No polynomial P(rho, mu) = 0 is possible.

4. **Cost comparison.** Computing on 4 sheets with rho-dependent integration
   limits is MORE expensive than the current 130 x 50x50 eigensolves, not
   less. The blow-up adds complexity without eliminating the per-rho
   diagonalization.

**Conclusion:** The algebraic blow-up correctly identifies the piecewise
structure but does not produce a computationally useful algebraic curve.
The obstruction is that spectral projection with rho-dependent integration
limits produces transcendental (not polynomial) matrix elements.


## Avenue 2: Lie Algebra Elevation

### Construction

The free angular problem (V_nuc = 0, V_ee = 0) has SO(6) symmetry, with
Casimir eigenvalues nu(nu+4)/2. The nuclear potential breaks this symmetry
differently in each region:

- Region s > rho: f_k ~ rho^k / s^{k+1}, which is a polynomial in rho
  times a power of 1/s. This preserves the homogeneity structure.
- Region s < rho: f_k ~ s^k / rho^{k+1}, which is a polynomial in s
  times a power of 1/rho. This has a different scaling.

Question: is there G containing SO(6) where both regions are restrictions
of a single representation?

### Assessment: NEGATIVE

The obstruction is fundamental:

1. **Different scaling dimensions.** In the s > rho region, V_nuc scales
   as rho^k * s^{-(k+1)}. In the s < rho region, it scales as
   s^k * rho^{-(k+1)}. These have DIFFERENT scaling weights under the
   SO(6) Casimir. No single irrep of a larger group can produce two
   different scaling behaviors in different regions of the SAME variable.

2. **The piecewise structure is NOT a symmetry breaking pattern.** Standard
   Lie algebra elevation works when a potential breaks G to H, and you embed
   both in G' where the full potential transforms as a single irrep. But
   min/max is not an algebraic function -- it is a C^0 but not C^1 function
   of (s, rho). No finite-dimensional Lie group has representations whose
   matrix elements are piecewise-smooth but not smooth.

3. **Comparison with the Fock precedent.** Fock's elevation works because
   1/r is the Green's function of the Laplacian on S3, which is a smooth
   object on S3. The transformation r -> chi(r) maps the singularity into
   a smooth function. The key: 1/r is singular in R3 but corresponds to a
   smooth geometric object (Green's function) on a compact manifold. In
   contrast, min(s, rho)/max(s, rho) is already defined on a compact domain
   [0, pi/2] and its non-smoothness is intrinsic -- it does not arise from
   a coordinate singularity that could be resolved by compactification.

4. **SO(6) to SO(7) or SO(8) analysis.** The natural candidates for
   G containing SO(6) are SO(7) and SO(8). SO(7) adds one dimension;
   the new coordinate would need to parametrize the sign of (s - rho).
   But this is exactly the blow-up from Avenue 1, and the same obstruction
   applies: the boundary moves with rho. SO(8) (relevant via triality for
   D4) does not help because the piecewise structure is not related to the
   D4 triality automorphism.

5. **Representation-theoretic impossibility.** For a Lie group G, matrix
   elements of irreducible representations are real-analytic functions.
   The function min(s, rho)/max(s, rho) is C^0 but not C^1 at s = rho.
   Therefore it CANNOT be a matrix element of any finite-dimensional
   representation of any Lie group. This is a clean obstruction.

**Conclusion:** The piecewise structure of V_nuc cannot be resolved by
Lie algebra elevation. The non-smoothness of min/max is intrinsic (not
a coordinate artifact) and is incompatible with matrix elements of
finite-dimensional Lie group representations.


## Avenue 3: S3 x S3 Hybrid

### Construction

Each electron has its own S3 via Fock projection:
- Electron 1 at r1: Fock-maps to S3_1, with 1/r1 smooth as Green's function
- Electron 2 at r2: Fock-maps to S3_2, with 1/r2 smooth as Green's function

The product S3 x S3 is 6-dimensional (matching the 6D configuration space
of two electrons in 3D). On S3 x S3:
- 1/r_i (nuclear attraction) is smooth on the corresponding S3 factor
- The basis is tensor products of S3 harmonics: Y_{n1,l1,m1} x Y_{n2,l2,m2}

### Assessment: NEGATIVE

Multiple obstructions:

1. **Nuclear attraction is NOT smooth on S3 x S3 for molecules.** The Fock
   map p0 = sqrt(-2E) depends on the energy eigenvalue. For a molecule with
   two nuclei at positions R_A, R_B, the nuclear potential is:
     V_nuc = -Z_A/|r - R_A| - Z_B/|r - R_B|
   Each 1/|r - R_X| maps to a smooth Green's function on S3 only if the
   S3 is centered at R_X. With two centers, there is no single S3 where
   BOTH 1/|r - R_A| and 1/|r - R_B| are simultaneously smooth. This is
   exactly the multi-center problem (Papers 8-9) that motivated the natural
   geometry hierarchy.

2. **V_ee representation.** On S3 x S3, the electron-electron repulsion
   1/|r1 - r2| couples the two factors. Track W proved: 1/r12 cannot be
   the Green's function on S5 (wrong singularity dimension). On S3 x S3,
   1/r12 is not the Green's function of any natural differential operator
   (the Laplacian on S3 x S3 is the sum of Laplacians, whose Green's
   function factorizes). So V_ee would need to be treated as an external
   charge function on S3 x S3, the same as in the current approach.

3. **The Fock map is energy-dependent.** The stereographic radius p0 depends
   on the energy being computed. For a two-electron molecule, there is no
   single p0 that simultaneously Fock-projects both electrons correctly
   (this is the structural theorem of Papers 8-9). The S3 x S3 construction
   inherits this problem: the Fock map for each electron depends on the
   total energy, creating a nonlinear self-consistency requirement.

4. **Equivalence to existing approach.** In molecule-frame hyperspherical
   coordinates, the alpha angle already parametrizes the relative distance
   of the two electrons (cos(alpha) ~ r1/R_e, sin(alpha) ~ r2/R_e). The
   channel structure (l1, l2) already represents the angular momenta on
   the two "S3 factors." The existing spectral solver IS effectively
   working on S3 x S3 restricted to a fixed hyperradius R_e. The min/max
   structure in V_nuc arises from the NUCLEAR positions being at finite
   distance from the origin, which no change of electron coordinates can
   remove.

5. **D_e ceiling.** Even if S3 x S3 worked, it would not exceed the 94.3%
   D_e ceiling because the limitation comes from partial-wave (channel)
   convergence of V_ee (the cusp), which is independent of the coordinate
   system for V_nuc.

**Conclusion:** S3 x S3 does not resolve the piecewise structure because:
(a) two-center nuclear attraction cannot be smooth on any single compact
manifold, (b) V_ee does not factorize on the product, and (c) the
approach is essentially equivalent to the existing one up to relabeling.


## Master Conclusion

All three avenues are NEGATIVE. The piecewise structure of the Level 4
angular Hamiltonian is IRREDUCIBLE -- it cannot be resolved by geometric
elevation to any higher-dimensional space. The obstruction is:

**The min/max structure arises from nuclear positions being at finite
distance from the coordinate origin.** Unlike the 1/r singularity (which
is a coordinate artifact removable by Fock projection to S3), the min/max
boundary s = rho is a PHYSICAL boundary: it is the surface where the
electron is equidistant from the origin and the nucleus. This surface
moves with rho = R/(2R_e), creating an intrinsically rho-dependent
partition of the angular domain. No coordinate transformation or group
embedding can make this partition rho-independent because the nuclear
positions themselves depend on rho.

**Classification:** The piecewise structure is an embedding exchange
constant in the Paper 18 taxonomy. It is the irreducible geometric content
of placing two nuclei at finite separation in the two-electron
hyperspherical space. Just as Track W showed 1/r12 is the irreducible
embedding exchange constant for electron-electron interaction, the min/max
boundary is the irreducible geometric consequence of two-center nuclear
attraction in hyperspherical coordinates.

**Interaction with Track AA (l_max convergence):** Geometric elevation
does NOT change l_max convergence, even if it had worked. The l_max
convergence rate is governed by the V_ee cusp (Track X, Track Y), which
is independent of V_nuc's smoothness. Resolving V_nuc's piecewise
structure would eliminate the need for per-rho diagonalization but would
not accelerate channel convergence.

**Computational implication confirmed:** The current approach (spectral
Jacobi angular solver, Track K, 130 x 50x50 eigensolves) is optimal
for the Level 4 angular sweep. The piecewise structure is irreducible,
and the spectral solver already achieves 269x speedup over FD. Further
improvements would come from caching/interpolation of eigenvalues across
rho-points, not from geometric restructuring.
