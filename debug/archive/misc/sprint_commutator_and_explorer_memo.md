# Sprint: commutator diagnostic (2) + Sturmian-addition-theorem Explorer (1)
**Date:** 2026-07-06/07
**Context:** follow-on to the /ahha A/B/C tests. PI approved "run (2)'s cheap
commutator diagnostic + dispatch (1)'s Explorer in parallel."
**Both lines independently converged on ONE object: the singular values of the
two-center overlap block S_AB.**

## Diagnostic (2): is the multi-focal wall a non-commuting-projections obstruction?
Model the two focal-length projections as the orthogonal projectors P_A, P_B onto
the single-center orbital subspaces of the LiH bond block (both Z=1). Same-Z
hydrogenic orbitals are orthonormal within a center, so the singular values of the
cross-overlap block S_AB are exactly cos(theta_k) (principal angles), and
||[P_A,P_B]|| = max_k cos·sin. Overlaps from the VALIDATED Topos-3
`overlap_two_center` (exact 1s-1s Slater to 12 digits). Driver
`debug/sprint_commutator_probe.py`, data `debug/data/commutator_probe.{json,log}`.

**Result (R_true=3.015, m=0 block {1s,2s,2p0}):**
```
  S_AB =  1s_B     2s_B    2p0_B
  1s_A  +0.3455  -0.2518  -0.4950
  2s_A  -0.2518  +0.7993  +0.1268
  2p0_A +0.4950  -0.1268  +0.4786
  singular values (cos theta) = 0.9913, 0.7107, 0.3856   (all < 1)
  principal angles            = 7.6, 44.7, 67.3 degrees
  ||[P_A,P_B]||               = 0.50   (NOT bit-zero)
```
**Verdict: claim (2) SURVIVES decisively.** The two center-projections do NOT
commute (three distinct principal angles, none 0 or 90 deg). So the multi-focal
wall's structural cause IS non-commutativity — there is no naive composed
projection. The structure is FINITE (3 angles per m-block) => the Jones basic
construction / intermediate subalgebra is finite and structured => "there is a
clean composed object we have not built" is viable, NOT a dense-continuum mess.
Also: S_AB is DENSE — large cross-n (<1s|2s>=-0.25) and cross-l (<1s|2p0>=-0.50)
entries => NOT block-diagonal in n; all sigma<1 => NOT an orthogonal congruence.

**Honest caveat:** this uses HYDROGENIC (physical, Z-scaled, position-space)
orbitals. Decisive for claim (2). SUGGESTIVE-not-decisive for Paper 8's specific
claim, which is about the shared-p0 STURMIAN (S^3) overlap — a different object
that needs its own computation (the convergent next step below).

## Explorer (1): orthogonalization-free two-center Sturmian integrals
**Verdict: GO** on the literal question. The Avery many-center Coulomb-Sturmian
apparatus gives all four two-center integrals (overlap S, kinetic T, V_ne, ERI) in
closed/recursive form on Fock's S^3, genuinely separate centers, and **without an
orthogonalization step to obtain them or to solve the classical generalized
eigenproblem** (metric = the SW overlap; Herbst-Avery-Dreuw run HF directly in the
non-orthogonal basis). Verified cites: Shibuya-Wulfman, Proc.R.Soc.A 286 (1965);
Avery, IJQC 100, 121 (2004); Avery, Adv.Quantum Chem. 67, 129 (2013);
Herbst-Avery-Dreuw, PRA 99, 012512 (2019, arXiv:1811.05777). **We already own a
fragment:** `geovac/shibuya_wulfman.py` (orthogonalization-free multipole cross-
center V_ne, `_lam` mismatched-exponent variant). Missing: the SW overlap metric +
hyperspherical-harmonic ERI.

**Honest caveats (Explorer):** (i) binding requires the scale k to FLOAT with R
(k(R) from the Sturmian secular equation, per-R algebraic — the guardrail's own
"R-dependent beta_k(R)" resolution, but closed-form not PDE-grid); the "graph"
is then R-parametrized, not one fixed lattice. (ii) The orthogonalization wall
REAPPEARS at qubit encoding — JW needs orthonormal spin-orbitals => Loewdin =>
Paper 8 BU-2's measured 2.8-4.5x Pauli-1-norm inflation (the encoding analog of
Track DF Sprint 5's 14x). So this helps CLASSICAL binding; the qubit product still
pays the orthogonalization tax.

## The flag (guardrail keystone — escalated, NOT acted on)
Paper 8's Sturmian Structural Theorem builds the cross-center block as the
**orthogonal SO(4) Wigner D-matrix** D^(n)(gamma), block-diagonal in n (line 324-
330: "the second nucleus acts as an SO(4) rotation, which cannot mix inequivalent
representations"), and asserts (line 335-337) "**No closed-form ... cross-nuclear
coupling exists in the literature, and the standard Shibuya-Wulfman formula is
exact only within a single n-manifold.**" Both lines here contradict that
dismissal: (a) the Explorer finds Avery's closed-form cross-n many-center SW
apparatus in the literature; (b) the probe shows the genuine two-center overlap is
dense/cross-n/sigma<1, NOT the block-diagonal orthogonal D-matrix.
**Precise re-pricing (narrow, fair):** the THEOREM is correct for its D-matrix
MODEL; what is in question is Paper 8's *dismissal of the cross-n escape hatch*
("no closed form exists / needs SO(4,2) / PDE grid"). If Avery's closed form is
real, the cross-n route to R-dependent binding is ALGEBRAICALLY accessible, which
Paper 8 rules out. **UNVERIFIED** against Avery primary text (paywalled; Explorer
flagged); touches a guardrail => PI call.

## Convergent decisive next step (both lines point here)
Compute the genuine **shared-p0 Coulomb-Sturmian** two-center overlap S_AB(kR) for
LiH and test: (i) does it mix n (refuting block-diagonal D-matrix)? (ii) sigma<1
(not an orthogonal congruence)? Plus verify Avery's primary text on whether the SW
overlap is single-n or cross-n. Cheap; it is the load-bearing hinge for BOTH the
Paper 8 re-pricing AND whether the Avery route is real. HOLD any code build
(extending shibuya_wulfman.py) until the hinge test resolves.

## Two prizes, different value (honest)
1. **Structural (high-confidence, essentially in hand):** multi-focal wall =
   non-commuting center-projections; Loewdin = forcing them to commute = the
   sparsity-destroying move. Clean characterization regardless of binding.
2. **Capability (uncertain, bigger):** algebraic Avery-Sturmian route to
   R-dependent binding that may re-price Paper 8's dismissal. Needs the hinge test
   + primary-text verification; and even if it binds classically, the qubit
   product still hits Loewdin — so its value for GeoVac's QC deliverable is open.

## HINGE TEST + VERIFY RESULTS (2026-07-07) — both CONFIRM cross-n
**PI decision was "hinge test + verify, hold build."** Both lines completed.

**(a) In-house: shared-p0 Coulomb-Sturmian two-center overlap MIXES n, strongly.**
Driver `debug/sprint_sturmian_hinge.py` (chi_{nl}(k) = R_nl_modern at Z=n*k; reuses
the validated Topos-3 two-center overlap). At k=1, R=3.015, m=0:
```
  within-n  <chi_1s|chi_1s>  = +0.3455
  CROSS-n   <chi_1s|chi_2s>  = -0.4711   (LARGER than the within-n diagonal)
  CROSS-n   <chi_1s|chi_2p0> = -0.5209   (LARGER still)
```
A block-diagonal SO(4) D-matrix gives EXACTLY 0 for the cross-n entries. They are
not just nonzero — they dominate. The genuine shared-p0 two-center Sturmian overlap
is emphatically NOT block-diagonal in n. (L2 metric; the reviewer's structural
argument shows the n-mixing comes from the TRANSLATION itself, metric-independent.)

**(b) Primary-text (citation-reviewer): Paper 8 lines 335-337 WRONG, flagged LARGE.**
- Q1 (n-mixing): CROSS-n. A second Coulomb center = real-space translation =
  e^{i p.R} in momentum space; on Fock's S^3 this is NOT an SO(4) isometry (would
  preserve lambda=n-1) but a plane-wave expansion over ALL hyperspherical degrees
  => couples all n. The SW integrals ARE the coefficients of that cross-n
  expansion; functions of s=kR (continuous), not a discrete D^(n)(gamma). "SW exact
  only within a single n-manifold" is BACKWARDS.
- Q2 (closed form): EXISTS. Avery, IJQC 100, 121 (2004); Red & Weatherford, IJQC
  (2004) "general formula ... no restrictions on the quantum numbers." So "no
  closed-form cross-nuclear coupling exists in the literature" is FALSE.
- **Careful scope (reviewer):** the guardrail THEOREM (shared-p0 single-center =>
  H,S share one SO(4) congruence => R-independent eigenvalues) is a DISTINCT result
  that does NOT depend on the false SW claim. The LARGE hit is on the literature-
  characterization PROSE (lines 335-337), NOT the structural theorem. The theorem
  STANDS; only Paper 8's DISMISSAL of the cross-n escape hatch is wrong.
- Sourcing caveat: paywalls blocked verbatim Avery-formula extraction; verdict
  rests on 4 mutually-consistent primary-grounded facts. PM re-verifies the LARGE
  hit against primary PDFs per no_sonnet_for_literature before any keystone edit.
- SMALL bib fixes found: DeFazio2003 initial P.->D.; Aquilanti1990 year 1990->1996
  + pages + title; Avery2006/avery_book duplicate; Aquilanti2001 key-year 2001->2003;
  Aquilanti1986 soft mis-cite.

**Three-line corroboration:** hydrogenic commutator probe (dense, cross-n, sigma<1)
+ shared-p0 Sturmian overlap (cross-n dominates) + primary literature (SW = cross-n
plane-wave expansion, closed forms exist) all agree.

**Honest deflation (must stay front-and-center):** this is NOT "we found how to
bind molecules." Herbst-Avery-Dreuw ALREADY bind (HF in non-orthogonal Coulomb-
Sturmians). What is new for GeoVac: (i) Paper 8's prose wrongly dismissed the cross-
n route -> correction needed; (ii) the cross-n SW coupling is the algebraic
composition object GeoVac was missing; (iii) the qubit-encoding Loewdin wall
(2.8-4.5x 1-norm) STILL bites -> whether the cross-n route helps GeoVac's actual
qubit product is the real open question. PI reserves the Paper 8 edit + build call.

## SESSION-2 CONTINUATION (2026-07-07): "verify+apply Paper 8 fix" + "probe the
## qubit-encoding wall". Both done. The Paper 8 issue is BIGGER than lines 335-337.

### Q2 verification (WebSearch, abstracts confirmed real, exact pages)
- Avery, IJQC 100, 121-130 (2004), DOI 10.1002/qua.10820: SW integrals via S^3
  hyperspherical harmonics "when the Sturmian basis functions are used in MOLECULAR
  calculations."
- Red & Weatherford, IJQC 100, 208-213 (2004): "an explicit and general formula for
  the Shibuya-Wulfman matrix, which may easily be programmed ... depends on the
  translation distance multiplied by the screening parameter" (continuous s=kR).
=> Paper 8 lines 335-337 both clauses refuted (Q1 by my computation, Q2 by these).

### Qubit-encoding wall diagnostic (`debug/sprint_encoding_wall_probe.py`)
- **m-selection (Hopf/azimuthal sparsity) is Loewdin-INVARIANT** (S m-block-diagonal
  => S^{-1/2} too; cross-m mass = 0 exactly) -> TRANSFERS to the qubit product.
- **within-m l-selection is DESTROYED** (Loewdin makes 2p0 46% l=0; ||S[l0,l1]||=0.72
  => no l-preserving orthogonalization exists; intrinsic, not a Loewdin artifact).
- **basic-construction (2x2 conjugate-pair) form does NOT rescue it** -- reaching it
  is itself a dense l-mixing rotation.
- **=> JW forces the l-mixing.** Prize #2 transfers ONLY the m-sparsity; the
  l-sparsity needs a NON-ORTHOGONAL fermionic encoding -- a research direction.

### The Paper 8 issue is THEOREM-SCOPE, not just stale prose
- **Backing test IMPOSES the orthogonal D-matrix cross-block** (`S_AB=f*D^(n)`,
  `V_AB propto D^(n)`, "hard-coded-and-asserted", test docstring l.43-62) => the
  R-independence ("corpse") is a **tautology of the imposed model**, not computed
  from the genuine two-center Sturmian overlap.
- **eq:sw_exact (l.613)** is the cross-nuclear POTENTIAL <S^A|-Z_B/r_B|S^A>, which
  legitimately collapses to a block-diagonal D-matrix via Sturmian 1/r-orthogonality.
  BUT the theorem ALSO needs the OVERLAP <S^A|S^B> to be that block-diagonal D-matrix,
  and the overlap has no 1/r to invoke that collapse -> it couples n. My shared-p0
  computation: <chi_1s(A)|chi_2s(B)> = -0.47 (cross-n). Theorem's OVERLAP premise
  FALSE for the genuine basis.
- **Paper already self-contradicts**: l.988-992 "different Sturmian n-shells GENUINELY
  DO couple ... making the D-matrix formula exact" -- opposite of l.335-337.
- **Net:** guardrail's "Sturmian molecular can't bind" force UNDERMINED (shows the
  D-matrix MODEL can't bind, not the genuine cross-n basis, which binds -- Herbst-
  Avery-Dreuw). cor:binding's "binding requires PDE / violates discrete premise" is
  wrong: cross-n SW gives R-dependence ALGEBRAICALLY (eq:sw_exact form, s=kR).
- => Fix is a theorem-framing reconciliation, bigger than the reviewer's l.335-337
  flag. PAUSED the edit; surfaced to PI (keystone scope grew).

### Non-orthogonal-encoding Explorer (PI: "scope the non-orthogonal encoding") = STOP
Does a fermion-to-qubit encoding keep S in {a,a^dag}=S_ij and preserve Gaunt l-sparsity
without a dense S^{-1/2}? **STOP for the literal object; BORDERLINE (NOCI/VB) only.**
- **Dual-basis theorem** (Artacho-del Bosch, PRA 43, 5770 (1991); Hu-Ratner-Seideman,
  arXiv:1511.06467): the metric never disappears -- push it into operators, states, or
  coefficients, it always lands as a dense S^{-1}/S^{-1/2} contraction. So "keep
  {a,a^dag}=S" and "keep the sparse Pauli Hamiltonian" are mutually exclusive for a
  single global second-quantized operator.
- **The one direct attempt** (Marruzzo et al., Adv. Quantum Chem. 92, 245 (2025),
  arXiv:2509.12680, non-orthogonal JW for VB): trades Loewdin S^{-1/2} for a
  biorthogonal S^{-1} -- which **densifies just as badly** (Pauli M^2/M^4 or M^3/M^6;
  loses Hermitian index symmetry). Confirms the wall from the qubit side; **S^{-1}
  densifies as badly as S^{-1/2}** -- the orthogonalizer choice is not the lever.
- **Only opening (BORDERLINE):** NOCI / VB subspace eigensolvers (Huggins et al., NJP
  22, 073009 (2020); Baek et al., PRX Quantum 4, 030307 (2023)) put non-orthogonality
  at the many-body STATE level, keep native integral sparsity in the classical
  matrix-element step, but pay with a device-measured many-body overlap S_IJ +
  generalized-eigenvalue conditioning. A different architecture, not a JW encoding.

**Arc closure:** prize #2 (cross-n Sturmian) helps CLASSICAL binding but does NOT
rescue the qubit sparsity -- confirmed from BOTH my encoding diagnostic (l-mixing
intrinsic) AND the literature (metric densifies regardless). The l-sparsity lives in
the BARE integral tensor; canonical anticommutation forces a dense-metric congruence;
the only way to keep it is to NOT transform the integrals (NOCI/VB), relocating the
cost to the many-body overlap. Prize #1 (multi-focal wall = non-commuting projections;
Loewdin = forced commutation = densification) stands as the solid, self-contained
deliverable.

## 6. Honest scope
- **Theorem grade:** no NEW theorem proved. One keystone **correction** landed:
  Paper 8 lines 335-337 retracted (SW is cross-n; closed forms exist), the Structural
  Theorem re-scoped to its imposed single-n D-matrix model (Remark `rem:overlap_imposed`),
  cor:binding softened. Backed by `tests/test_paper8_overlap_cross_n.py` (PASS) + verified
  primary abstracts.
- **Structural (well-supported, computed):** multi-focal/composition wall = non-commuting
  center-projections (‖[P_A,P_B]‖=0.50 at R_true; three principal angles 7.6/44.7/67.3°);
  the genuine two-center Coulomb-Sturmian overlap is cross-n (⟨χ_1s|χ_2s⟩=−0.47 shared-p0,
  bit-validated vs exact 1s-1s Slater); Löwdin is m-selection-invariant (cross-m mass=0
  exactly) but l-selection-destroying within m (2p0→46% l=0), intrinsically (no
  l-preserving orthogonalization; conjugate-pair form is itself a dense l-rotation).
- **Literature-grounded (verified to abstract level):** Avery cross-n SW closed forms exist
  (Avery IJQC 100/121; Red-Weatherford IJQC 100/208, both 2004); non-orthogonal qubit
  encoding = STOP (dual-basis theorem Artacho-del Bosch 1991; Marruzzo 2025 biorthogonal
  S⁻¹ densifies as badly as Löwdin S⁻¹ᐟ²); NOCI/VB is the only BORDERLINE opening.
- **Numerical observation (single-system, LiH; hydrogenic + shared-p0 k=1; mostly one R):**
  the overlap matrices, gradient deficit, contamination fractions above.
- **Named open follow-ons (none launched):** (a) NOCI/VB-on-GeoVac feasibility (the only
  literature opening — relocates metric cost to a device-measured many-body overlap S_IJ;
  scope before any build); (b) whether the cross-n *classical* binding demonstration is
  worth building given the qubit product does not benefit (PI call); (c) the A-vs-C
  structural-wall-vs-Pulay distinction from the sibling A/B/C sprint (needs balanced n_max≥4,
  out of reach).
- **Caveats:** probes are hydrogenic / shared-p0 at k=1 and mostly single-R; the Q2
  closed-form claim rests on verified *abstracts* (Avery/Red-Weatherford full text paywalled;
  primary formula not verbatim-extracted); the cross-n *overlap* refutation is decisive (my
  computation + the SW-is-a-translation structural argument), the "closed form exists" leg is
  literature-grounded not re-derived.

## 7. Follow-on: the sparsity-DESTROYING option (PI-directed 2026-07-07)
PI asked to explore actually paying the density to buy binding: restore the neglected
overlap S in balanced LiH, Loewdin-orthogonalize the integrals, measure geometry + Pauli.
Driver `debug/sprint_lowdin_tradeoff.py`, data `debug/data/lowdin_tradeoff.json`.

**Result: the cheap retrofit is STRUCTURALLY IMPOSSIBLE (clean negative).**
- Sparse baseline reproduces exactly (R_eq +6.9%, omega_e +45% / 2037 cm^-1, 878 Pauli) --
  machinery validated. Robust two-centre overlap (fixed Gauss-Legendre, 1s-1s exact to 1e-16)
  replaced the hanging topos3 adaptive quad.
- The restored overlap S is **non-PSD** (S_min=-0.037 at short R; 2 near-null directions):
  the balanced basis is near-linearly-dependent once overlap is honest (redundant s-functions:
  core-1s Z=3, bond-1s Z=1, partner-1s overlap 0.33-0.65).
- Loewdin (even floored) gives **garbage** (E=-54 Ha vs -15.2 correct; NO interior minimum;
  Pauli 878->15706 = 17.9x, consistent with Track DF's 14x).
- **Diagnosis:** the builder's h1/eri are NOT clean <i|h|j> matrix elements of a basis with a
  computable overlap (h1 diagonal = isolated-atom -Z^2/2n^2; cross terms = multipole /
  operator-approximation surrogates). Bolting an independent S + Loewdin onto them mixes
  incompatible objects. To do it correctly you must rebuild ALL integrals (S, h1, eri)
  CONSISTENTLY from the same basis functions = standard two-centre molecular Coulomb-Sturmian
  QC (Avery / Herbst-Avery-Dreuw) -> GeoVac sheds its only distinctive advantage (sparsity).
- **Honest scope:** did NOT prove a proper cross-n Sturmian build fails to bind (Herbst-Avery-
  Dreuw bind fine); showed (a) the retrofit is inconsistent/impossible, (b) the proper version
  IS standard QC at ~18x the Pauli cost, and (c) A/B/C already predicted overlap-restoration
  would not fix the geometry (defect = max_n orbital-basis completeness, a different axis).
- **Verdict: the sparsity-destroying option is a REPLACEMENT of GeoVac with standard QC, not a
  modification of it.** Not worth pursuing unless the goal is to reimplement molecular Sturmian
  QC. Confirms the session's through-line: sparsity and accuracy are the two ends of one trade;
  you cannot keep the sparse Pauli advantage and fix the geometry at the same time.
