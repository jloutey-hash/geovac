# Sprint: polyatomic composition = irreducible three-body operator structure (2026-08-25)

Canonical memo. PI-driven conceptual arc (a single conversation of structural reasoning)
that produced a concrete, firmed new structural result: **the polyatomic composition wall
is not reducible to the pairwise/diatomic picture — a third center opens a large,
irreducible three-body operator-algebra channel.**

## 1. The conceptual arc (PI reasoning → the object)
Working from the chemistry accuracy/γ-fragility thread, the PI reasoned bottom-up to the
non-commutative-geometry structure of the framework, in steps that each checked out:
must-project (the skeleton is not spatial) → *which* projection (metric vs density; γ(r)
is the density projection) → combine the projections (multiply them) → counter-rotate the
tilt (the commutator [P,Q] carries the ±i/2 and IS the rotation generator = the gauge
connection; the counter-rotation is Löwdin, dense) → sum over orderings and aggregate "like
Feynman" (= the Dyson series / spectral action / JLO cocycle). This is the WH1 / spectral-
triple program rederived from chemistry intuition. The open, undeveloped object it lands on:
the operator-algebra of the **multi-center** projections, and whether 3+ centers carry
irreducible three-body content.

## 2. Current-state check (before claiming novelty)
- Spectral action → chemistry: **thoroughly dead** (3 §3 dead-ends — wrong functional for
  chemistry observables at finite cutoff; doesn't bind; no Pauli reduction).
- JLO cocycle: **computed on the base spectral triple** and comes out **cocommutative**
  (S₃-invariant, order-independent) — because the base algebra is commutative. So "aggregate
  over orderings" is trivial on the discrete/graph side. This *confirms* the framework's own
  multi-focal split (discrete labels compose cleanly; spatial projections don't) — the
  non-commutativity is in the spatial center-projections, not the base algebra.
- The multi-center composition as a gauge/JLO object was **flagged but never built** (the
  composition-wall memo's "Jones basic construction / clean composed object we have not
  built"). So the 3-center question is genuinely undeveloped.

## 3. LiH (2 centers) — reproduced + rigorously cross-checked
The two center-projections (bond block, R=3.015) meet at three principal angles
**7.6° / 44.7° / 67.3°**, ‖[P_A,P_B]‖ = **0.500** (the 44.7° angle sits at the 45° ceiling —
maximal non-commutativity). Verified two ways: the SVD-of-cross-overlap shortcut AND building
the actual 6×6 projectors and computing [P_A,P_B] directly (both 0.500, bit-agree).
Key: for **two** projections the pairwise angles are the *complete* invariant (Halmos / CS
decomposition) — flat, tame, closed.

## 4. BeH₂ (3 centers) — the new result: large irreducible three-body content
Linear H–Be–H, σ pair {1s, 2p0} per center, three metric-orthogonal center-projections,
the nested commutator ‖[[P_A,P_B],P_C]‖ as the 3-way invariant (**exactly 0 iff the three
reduce to pairwise blocks**). Driver `debug/beh2_three_way_invariant.py`, data
`debug/data/beh2_three_way.json`. Overlaps from the validated Topos-3 two-center engine
(mismatched-Z supported, spot-checked); linear-geometry parity (2p0 odd) handled explicitly.

| | 3-way bonding (d=2.5) | far (d=15) |
|---|---|---|
| Structural (Z=1 all) | **0.421** (86% of pairwise 0.49) | 0.019 |
| Physical (Be Z=2, H Z=1) | **0.320** (65% of pairwise 0.49) | 0.0006 |

- **Large:** the three-body invariant is 0.32–0.42 at bonding — 65–86% of the *pairwise*
  magnitude, not a perturbative correction. For a diatomic this quantity cannot exist.
- **Real bonding structure:** decays cleanly to ~0 as the atoms separate (0.0006 at d=15 with
  compact physical Be).
- **Survives charge asymmetry:** Be Z=2 vs H Z=1 gives 0.32 (vs 0.42 at Z=1), still large.

## 5. Interpretation
**Diatomic = pairwise (Halmos, flat). Polyatomic = pairwise + a large irreducible
three-body term (curved).** The third center opens an operator-algebra channel the two-center
case cannot have — the concrete "curvature/holonomy" of the composition (the leading
non-trivial term of the PI's "aggregate over orderings": the nested commutator is exactly the
first obstruction to the orderings agreeing). This is a precise, structural, on-thesis reason
polyatomics are categorically harder than diatomics — a forced/free stake, not a chemistry
tool. Chemistry payoff remains dead (§2); this is a *characterization* of the wall.

## 5b. The FINISHED classification (2026-08-25 follow-on: full algebra)
Driver `debug/beh2_algebra_classification.py` (self-validating: the `commutant()` routine is
checked on three known cases — single rank-2 → dim 20, commuting pair → 12, generic 2-dim
pair → 1). Result for physical BeH₂ (Be Z=2, H Z=1, d=2.5), H dim = 6:

- **The triple generates the FULL matrix algebra M₆ and acts IRREDUCIBLY** — commutant
  dim(A′) = **1** (only scalars commute with all three). By von Neumann double-commutant,
  A′ = ℂI ⟺ A = M₆. The whole 6-dim space is a single irreducible three-body block
  (census: three-way 6, pairwise 0, trivial 0).
- **Every PAIR is REDUCIBLE** — dim(A′) = **6** for {Be,H1} and {H1,H2} (Halmos: decomposes
  into ≤2×2 blocks, pairwise angles complete).

**The categorical jump, finished:** diatomic (2 projections) = REDUCIBLE (always ⊕≤2×2, fully
classified by angles); polyatomic (3 projections) = IRREDUCIBLE (the full matrix algebra, no
decomposition exists). The three centers weld the entire space into one indecomposable object.
The sharpest form of "categorically harder": the three-center problem cannot be decomposed at
all. Bonding-distance phenomenon (far → commuting → reducible, per the §4 far control).

**Redirection of the JLO piece (what finishing the theorem revealed):** since the algebra is
just M₆, the JLO cocycle of the *bare algebra* is standard (M_n has trivial cyclic
cohomology). The depth is NOT the algebra type (settled: irreducible M₆) but the **three-
subspace configuration** — how the three projections sit inside M₆ (the angles) and how that
configuration **varies over molecular geometry** (the Bargmann/geometric-phase connection and
its curvature over the nuclear-geometry manifold). That geometric family is the genuine
remaining object, cleaner than "JLO of the algebra."

## 5c. The geometric family — a lattice of conical intersections + π Berry phase (2026-08-25, sub-agent-verified)
Since the algebra is settled (M₆), the live object is the three-subspace CONFIGURATION and its
variation over molecular geometry. For real orbitals the Bargmann invariant Δ=Tr(P_Be P_H1 P_H2)
is real ⇒ the geometric phase is Z₂. Three parallel sub-agents (Lines 1/2/3) took this from a
single interpolated CI to a verified, exact, bit-checked picture — including an honest deflation.
Drivers: `debug/beh2_bargmann_geometry.py`, `beh2_ci_refine.py` (linear survey);
`beh2_ci_exact_landscape.py` (Line 1, exact census + fast overlap engine);
`beh2_bending_ci*.py` + `fast_two_center_overlap.py` (Line 2, bent geometry); `beh2_onebody_probe.py`
(Line 3). Data in `debug/data/beh2_*.json`.

- **Exactly 3 CIs in the linear stretch region (Line 1, exact):** central on the C₂ᵥ line at
  **d\*≈2.445 bohr** (exact F-gap 2.6e-8) plus an H₁↔H₂ **mirror pair** off-axis at
  (2.698,2.266)/(2.266,2.698) (gap 3.3e-8). All three genuine diabolical points (linear dispersion
  slope ≈0.14; tight loops π; band-2/3 crossings at F≈0.71). A fast float64 prolate-spheroidal
  overlap evaluator (validated vs mpmath to 2.3e-14) replaced the slow/interpolated route — the
  earlier π→0→π radius-flip was a *grazing* artifact (loop skimming the off-axis pair at r≈0.31),
  but its inference "off-axis CIs exist" was correct.
- **π Berry phase (Z₂):** parallel-transported F-eigenvector sign-flips around any clean loop
  enclosing an odd number of CIs; radius-independent ⇒ a genuine topological invariant; far loops 0.
- **Bending enriches + a structural Renner–Teller mechanism (Line 2, GO, bit-validated):** the bent
  oriented-overlap engine (Slater–Koster) reproduces the linear F-spectrum at θ=180° to 1.7e-13.
  The central CI survives as a continuous curve **d\*(θ)** (2.445@180° → 2.02@110°), π at every
  angle/radius; the landscape grows to ≥3 CI branches. Renner–Teller appears structurally: at
  θ=180° a σ- and a π-level *touch* with the σ–π block exactly 0 (Berry 0, non-conical); bending
  switches on σ/π coupling (0→0.145 as θ:180°→150°) and **promotes it into a real CI** — a CI that
  cannot exist on the linear axis. (Correctness note: the linear driver's Cholesky whitening is not
  Σ-equivariant, Löwdin is — but eigenvalues/CIs/the algebra classification are whitening-invariant,
  so all numbers stand; Löwdin only sharpens even/odd labels.)

**Honest boundary (Line 3, DISTINCT):** this is the abstract CONFIGURATION operator F=ΣP_i, NOT the
electronic Hamiltonian. BeH₂'s genuine electronic CI is the classic Be+H₂ *insertion* intersection
(two-configuration ¹A₁/¹A₁) at a **bent, inserted** geometry (Be–H≈3.0, H–H≈2.55 bohr; Purvis–Shepard–
Brown–Bartlett 1983); our config-CI sits at the *linear equilibrium* where the real electronics are
benign/single-reference. Different geometry, deformation space, and operator; GeoVac's own BeH₂ builders
are linear-only and (composed) block-diagonal, so they cannot even reach the physical CI's geometry.
The shared π phase + antisymmetric-stretch coupling are **generic** to any symmetric-triatomic CI.

**Corrected framing:** the operator-algebra result (irreducible M₆) and this geometric phase are the
same three-center object two ways — the algebra says the centers weld into one irreducible whole; the
geometry says that whole carries a *rich lattice* of conical intersections with π holonomy, including a
structural Renner–Teller mechanism. This is a genuine geometric-phase structure of the composition —
the **structural/generic cousin** of molecular CI physics, NOT a claim of identity with BeH₂'s physical
electronic conical intersection (Line 3, DISTINCT).

## 6. Honest scope
- **Firmed:** the 3-body signal is large at bonding, decays cleanly to ~0 far, and survives
  physical charge asymmetry (two independent charge regimes agree).
- **Caveat (remaining):** the nested commutator is *one* defensible diagnostic (0 iff
  pairwise-reducible). The *full* operator-algebra classification (the three-projection
  indecomposables / the actual JLO cocycle of the multi-center algebra) is the natural next
  step — that would turn "large 3-body content" into a complete structural statement.
- **Basis:** σ pair {1s,2p0}, m=0; a richer basis (add 2s, higher l/m) would quantify but is
  not expected to remove the signal.
- **Paper candidate (PI-gated):** a remark in Paper 32 (spectral triple) or a short structural
  note — "the polyatomic composition wall carries irreducible three-body operator content" —
  NOT auto-written (fresh result + the one-diagnostic caveat).

Drivers: `debug/beh2_three_way_invariant.py` (BeH₂), `debug/sprint_commutator_probe.py` (LiH,
pre-existing). Data: `debug/data/beh2_three_way.json`.

## 7. DIAGNOSTIC follow-on (2026-08-25): is the reducible->irreducible jump generic?
Drivers `debug/diag_three_center_generality.py` (+ `diag_three_center_discriminator.py`),
data `debug/data/diag_three_center_generality.json`. Overlaps validated: fast engine vs mpmath
topos3 incl. the new 2s states, err <=4e-13. commutant() self-check (20/12/1) passes.

**Q1 basis robustness (linear, sigma menu enriched {1s,2p0}->{1s,2s,2p0}):** the triple stays
IRREDUCIBLE (M6 -> M9) at Z=1 and Z=2; pairs stay reducible (dimA' 6->12). Adding 2s does not
break it. GO.

**Q2 geometry generality (in-plane {s,px,pz}x3, 9-dim, Slater-Koster):**
- **New clean result:** at EXACTLY linear (180 deg) the in-plane space is REDUCIBLE
  (dimA'=2; sigma{s,pz} + pi{px} decouple by axial symmetry). ANY bend (160..90 deg) welds it
  IRREDUCIBLE M9. This is the operator-algebra dual of the memo's Renner-Teller bending
  mechanism -- bending mixes sigma/pi and welds the full in-plane algebra. Symmetry-protected,
  tolerance-independent.
- **Real H2O geometry (R_OH=1.809, 104.5 deg) is IRREDUCIBLE M9** at apex charge Z_O=1,2,4.
  Generality confirmed for a bent, heavier-apex triatomic.

**Important refinement (the paper remark slightly conflates two invariants):** commutant dim=1
("irreducible", the CATEGORY jump 2->3) is a KNIFE-EDGE / GENERIC property of >=3 coupled
projections -- the distance sweep shows it holds from d=2 out to d=15 and only flips to
reducible once the residual overlap drops below the numerical tol (d>=20). It is a statement
about center COUNT, not bonding. The BONDING-specific signal is the nested-commutator
MAGNITUDE ||[[P,P],P]|| (0.35 @ bonding -> 0.04 @ d=10 -> ~0 far), which decays smoothly.
Report BOTH: (A) categorical irreducibility (robust, basis/geometry-general) and (B) three-body
strength (physical, bonding-dependent).

**Theory-home CORRECTION (supersedes the session's earlier "D4 finite-type ladder"):** the
object is ORTHOGONAL PROJECTIONS (a *-representation question), not arbitrary linear subspaces
(GL/quiver). The correct dichotomy is Halmos: TWO projections = TAME / Type I (universal
C*-algebra = subalgebra of M2(C[0,1])); THREE-OR-MORE projections = *-WILD (unclassifiable).
So the 2->3 center jump is exactly the tame->wild threshold in operator algebras, and the
molecular triples robustly land on IRREDUCIBLE (Schurian) points of that wild family -- geometry
SELECTS a tame-looking point out of a wild problem. Verified vs literature (arXiv:1207.6890
"C*-algebras generated by three projections"; Halmos two-projections theorem). The linear-
subspace D4/tilde-D4 ladder (finite/tame/wild at 2/4/5 subspaces) is the WRONG lens for this
object and is retracted from the framing.

## 8. The Berry curvature over the nuclear-geometry manifold (2026-08-25, computed)
Driver `debug/beh2_berry_curvature.py`, data `debug/data/beh2_berry_curvature.json`. The memo's
named "genuine remaining object", now finished. F = sum_i P_i over the linear-stretch (d1,d2)
plane; sigma pair {1s,2p0} => F is 6x6 REAL-SYMMETRIC (validated exact overlaps). Fukui-
Hatsugai-Suzuki lattice Berry flux on a 141x141 grid.

**Result — the curvature is a flat Z2 connection with pi-flux deltas at the CIs:**
- Plaquette flux is EXACTLY in {0, pi} everywhere (max|flux-{0,pi}| = 0.0, machine-exact) =>
  the connection is FLAT; real-symmetric forces Z2, not a smooth U(1) field.
- Curvature = **3 pi-flux sources** at the known CIs: central (2.446,2.446) + antisymmetric
  mirror pair (2.268,2.696)/(2.696,2.268). Total 3*pi. Carried ONLY by the two crossing bands
  (2,3); the other four bands are flat/sourceless.
- Holonomy = Z2 rep of pi_1(M \ CI) -> {+-1}, loop phase = pi*(#enclosed mod 2). Cross-checked
  by the independent parallel-transport routine: pi around 1 CI AND pi around all 3 (odd). Agrees.
- Each source is a genuine diabolical point: gap ~ linear in radius (slope->0.141), anisotropic
  cone (min/max ~ 1/13).

**Each CI = the real (equatorial) section of a charge-1/2 Berry monopole (local 2x2 model):**
in a fixed real reference frame Q (crossing eigenvectors at the CI), F2 = a0 I + a_x sx + a_y sy
+ a_z sz with **a_y == 0 exactly** (real family sits on the sphere's equator). Branching vectors
**g = grad a_z = [0.0506,0.0506] (proportional to the SYMMETRIC stretch d1+d2 = tuning)** and
**h = grad a_x = [0.0038,-0.0038] (proportional to the ANTISYMMETRIC stretch d1-d2 = coupling)**;
Jacobian det(g,h) = -4e-4 != 0 => (a_x,a_z) winds once => pi holonomy. The missing sy (a_y) axis
is the reality/planarity-breaking direction that would lift Z2 -> a smooth U(1) monopole.

**Finished statement.** The configuration operator's geometric phase over nuclear geometry is a
FLAT Z2 line bundle punctured at its conical intersections, curvature = sum_CI pi*delta; it is
the Longuet-Higgins / Herzberg-LH molecular-Aharonov-Bohm structure realized for F (not the
electronic Hamiltonian). Z2 (not U(1)) is FORCED by TIME-REVERSAL (reality of the overlaps), NOT
planarity -- CORRECTED in section 9: a non-coplanar *real* center stays real-symmetric, so the
lift needs a magnetic/complex (T-breaking) phase, not geometry. g//symmetric, h//antisymmetric stretch ties the
sources to the Renner-Teller coupling coordinate. PI-gated Paper 32 capture candidate (not
auto-written): extend rem:multicenter_composition or a short new remark + a backing test.

## 9. The U(1) lift — each CI is a charge-+-1 Berry monopole (2026-08-25, computed)
Driver `debug/beh2_u1_lift.py`, data `debug/data/beh2_u1_lift.json`. Realizes Z2 -> U(1) and
corrects a wording error from section 8.

**Correction (section 8 said "planarity"):** the Z2 is protected by TIME-REVERSAL SYMMETRY = the
reality of the Coulomb overlaps, NOT by molecular planarity. Real orbitals give real overlaps in
ANY 3D geometry, so a non-coplanar *real* 4th center leaves F real-symmetric and does NOT lift
Z2. The "non-coplanar center" route named earlier is WRONG for the lift; only a genuine complex
(T-breaking) phase lifts it.

**The lift (magnetic/Peierls phase):** put e^{i phi} on the H1<->H2 overlap block (flux phi
through the H-Be-H loop) -> F(d1,d2,phi) complex-Hermitian. Verified:
- **F(phi=0) exactly real** (max|imag|=0); phi **opens the gap linearly** at the CI (slope ~0.55)
  = the empty sigma_y axis switched on.
- **3x3 branching Jacobian** (a_x,a_z,a_y over d1,d2,phi) det = 1.06e-4 != 0 => genuine 3D cone.
- **First Chern number over a sphere enclosing (d\*,d\*,0)** = **+1 (lower crossing band) / -1
  (upper)** = a genuine integer U(1) Berry monopole. The real phi=0 plane cuts it at the EQUATOR,
  recovering the pi (Z2) holonomy -> the CI is the T-symmetric section of a charge-+-1 monopole.
- Illustration: at phi=0 the (d1,d2) plaquette flux is exactly in {0,pi} (Z2); at phi!=0 it is no
  longer pinned there (smooth U(1) curvature).

**Captured:** Paper 32 `rem:config_berry_curvature` last paragraph corrected+extended (time-
reversal-forced Z2 + the magnetic-phase U(1) monopole, Chern +-1); backing test
`test_paper32_berry_curvature.py::test_u1_lift_is_integer_monopole` (slow, ~13s: F(0) real, gap
prop phi, Chern +-1). CHANGELOG v5.1.0 + claim_test_matrix updated. Paper compiles clean (87 pp).

**Further step (not done):** the 4-center (tame->wild / D~4 four-subspace) route carries a
continuous cross-ratio parameter classically; whether that is the SAME lambda as Paper 59's
Legendre/Gamma(2) elliptic frontier is the open cross-corpus lead (an /aha-flavored question).
