# The composed nuclear–electronic deuterium Hamiltonian as an explicit Connes-style real-space multi-particle spectral triple

**GeoVac Project Zenodo Memo, 2026-05-07**
*Track NI standalone observation memo, Observations status*
*Cross-references: Paper 23 §VI; Paper 32 §V, §VIII.D*

---

## Abstract

Track NI of the GeoVac project (Paper 23 §VI) constructs a composed
Hamiltonian on a 26-qubit register coupling a 16-qubit deuteron Hamiltonian
on the harmonic-oscillator basis with a 10-qubit hydrogenic electronic
Hamiltonian on the Fock-projected $S^3$ lattice. Cross-register coupling
reproduces the hydrogen 21 cm singlet–triplet hyperfine gap analytically at
$1.62 \times 10^{-7}$ Ha. We position this construction inside the
Connes–van Suijlekom spectral-truncation framework, audit it against the
standard Connes axioms following Paper 32, and observe that the canonical
Connes-style noncommutative-geometry literature (Standard Model
almost-commutative geometry; spectral-truncation convergence theorems;
type III $\sigma$-spectral triples) does not contain a directly comparable
real-space multi-particle construction. The empirical content is a
proof-of-concept; the structural observation is the gap.

---

## 1. Introduction

The GeoVac framework (Papers 0, 7, 14, 21, 32) treats the discrete graph
Laplacian on the Fock-projected three-sphere $S^3$ as a dimensionless
topological object that encodes hydrogenic spectral data via the conformal
projection of Fock 1935. Paper 32 reads this construction as a graph
spectral triple in the lineage of Marcolli–van Suijlekom 2014 + Connes–van
Suijlekom 2021. Working hypothesis WH1 (CLAUDE.md §1.7) was upgraded to
PROVEN in May 2026 (Paper 38) on the closure of a five-lemma
Gromov–Hausdorff convergence theorem for the truncated triple.

What this memo formalizes is a *separate* observation: that the Track NI
construction in Paper 23 §VI — built for entirely empirical reasons, to
demonstrate the composed-architecture machinery on a multi-species qubit
register — is structurally an *explicit Connes-style real-space
multi-particle spectral triple*, and that the published noncommutative
geometry / spectral-triple literature does not contain a directly comparable
construction. The Standard Model almost-commutative geometry of
Chamseddine–Connes couples a Riemannian manifold $M$ with a finite
*internal* factor $F$; the spectral-truncation convergence machinery of
Connes–van Suijlekom 2021 / Hekkelman–McDonald 2024 / Latrémolière 2026
covers single triples or one-infinite-Riemannian-plus-one-finite tensor
products. Track NI does neither: it constructs a tensor product of *two
infinite real-space spectral triples* for two physically distinct quantum
particles at categorically different focal lengths, with a cross-register
hyperfine coupling on the spin sectors.

This memo is an observation, not a theorem. It does not prove a convergence
theorem; it does not derive new physics; it does not claim the construction
is the *first* of its kind in the broader chemistry-physics literature
(see §6 for caveats). It earns its DOI for a structural literature
observation paired with a calibrated empirical validation.

---

## 2. The composed nuclear–electronic Hamiltonian

The Track NI Hamiltonian (Paper 23 §VI Eq. eq:composed-ne) acts on the
tensor-product Hilbert space

$$\mathcal{H}_{\mathrm{NE}} \;=\; \mathcal{H}_N \otimes \mathcal{H}_e,$$

where $\mathcal{H}_N$ is the 16-qubit nuclear deuteron register on the 3D
harmonic-oscillator basis with Minnesota NN potential (Paper 23 §IV) and
$\mathcal{H}_e$ is the 10-qubit hydrogenic electronic register on the
Fock-projected $S^3$ lattice at $n_{\max} = 2$ (Paper 14 §IV).

The full Hamiltonian decomposes into four blocks:

$$H_{\mathrm{NE}} \;=\; H_N \otimes I_e \;+\; I_N \otimes H_e
  \;+\; H_{\mathrm{fs}} \;+\; H_{\mathrm{hf}},$$

with $H_N$ the nuclear deuteron block (592 Pauli terms on 16 qubits),
$H_e$ the electronic hydrogen block (10 Pauli terms on 10 qubits),
$H_{\mathrm{fs}} = -2 Z \alpha m_e r_p$ a classical scalar shift
(absorbed into the identity component, evaluated at $r_p = 0.0316$ bohr),
and $H_{\mathrm{hf}}$ the cross-register hyperfine coupling
$A_{\mathrm{hf}} \mathbf{I}_p \cdot \mathbf{S}_e$ implemented via 12 Pauli
strings on the spin sectors of the proton 0s and electron 1s registers.

Resource counts (Paper 23 §VI.2 Table tab:ne-scales):

- 26 qubits total (16 nuclear + 10 electronic);
- 614 non-identity Pauli terms (592 nuclear + 10 electronic + 12
  cross-register hyperfine);
- 13-orders-of-magnitude coefficient hierarchy (~10⁵ MeV nuclear,
  ~10⁰ Ha electronic, ~10⁻⁷ Ha hyperfine);
- Reproducible from `geovac/nuclear/nuclear_electronic.py`.

Empirical validation: the singlet–triplet gap at the hyperfine block
evaluates analytically to

$$\Delta E_{\mathrm{hf}} \;=\; \frac{3 A_{\mathrm{hf}}}{4} \;=\;
1.62 \times 10^{-7} \mathrm{Ha},$$

matching the canonical hydrogen 21 cm transition wavelength to the
precision of the input $A_{\mathrm{hf}}$ value. This is not a prediction
—  it is a reproduction of $A_{\mathrm{hf}}$ as a Layer-2 calibrated input
from $|\psi_{1s}(0)|^2 = Z^3/\pi$ and the electron and proton magnetic
moments. The structural content of the empirical validation is that the
cross-register $\mathbf{I} \cdot \mathbf{S}$ operator, encoded as 12 Pauli
strings on a multi-species register, returns the physical answer when
evaluated.

---

## 3. Spectral-triple structure

Reading the Track NI Hamiltonian as a spectral triple in the sense of
Paper 32 §III–IV:

- **Algebra.** $\mathcal{A}_{\mathrm{NE}} = \mathcal{A}_N \otimes
  \mathcal{A}_e$, with $\mathcal{A}_N$ the multipliers on the 3D HO
  shell-basis register (functions diagonal in the $(n_x, n_y, n_z)$
  basis) and $\mathcal{A}_e$ the multipliers on the Fock-projected $S^3$
  lattice (diagonal in the $(n, l, m_l)$ basis, per Paper 32 §III). Both
  factors are infinite-dimensional in the continuum limit and are
  finite at each truncation $n_{\max}$.
- **Hilbert space.** $\mathcal{H}_{\mathrm{NE}} = \mathcal{H}_N \otimes
  \mathcal{H}_e$, with the nuclear factor carrying the deuteron's
  isospin and spin labels (proton + neutron, $\mathbf{I}_p$,
  $\mathbf{I}_n$ irreducibles) and the electronic factor carrying the
  KO-dim-3 Camporesi–Higuchi spinor bundle truncated to $n_{\max} = 2$.
- **Dirac operator.** $D_{\mathrm{NE}} = D_N \otimes I_e + I_N \otimes
  D_e + D_{\mathrm{cross}}$, with $D_N$ the HO-basis Dirac analog (HO
  shell label as eigenvalue), $D_e$ the truthful Camporesi–Higuchi
  spectrum on the electronic register, and $D_{\mathrm{cross}}$ the
  cross-register operator implementing the hyperfine
  $\mathbf{I}_p \cdot \mathbf{S}_e$ coupling on the spin sectors.
- **Real structure** (where applicable). The composed real structure
  $J_{\mathrm{NE}} = J_N \otimes J_e$ inherits Connes axiom compatibility
  from the factors: $J_e$ from the truthful Camporesi–Higuchi sector
  (Paper 32 §IV, KO-dim 3); $J_N$ from the HO-basis nuclear Dirac analog
  (KO-dim 0 or 8 in standard convention).

The structural feature that distinguishes this construction from the SM
almost-commutative geometry is $D_{\mathrm{cross}}$: in Connes' SM ACG,
the cross-block coupling lives on the *finite* internal factor (the
Yukawa Dirac block on $A_F$), while in Track NI it lives on the *spin
sectors of two real-space registers*. The KO-dimension structure also
differs: SM ACG combines KO-dim-4 manifold with KO-dim-6 internal
factor for total KO-dim 10 ≡ 2 (mod 8); Track NI combines two
real-space factors with their own KO-dim assignments.

The spectral-triple-level interpretation is *not* a re-derivation of the
empirical content of Paper 23 §VI; it is a re-reading of an existing
construction in the spectral-triple vocabulary. The empirical content
(qubit count, Pauli count, hyperfine validation) is unchanged.

---

## 4. Comparison to Connes' Standard Model construction

The Connes–Marcolli–Chamseddine SM construction (Chamseddine–Connes 2010;
Connes–Marcolli 2008; Lizzi 2018 review arXiv:1805.00411; Pati–Salam from
NCG arXiv:2511.07672) couples spacetime (Riemannian manifold $M$) with a
finite-dimensional internal factor $F$ to produce an almost-commutative
spectral triple $M \times F$. The internal factor $F$ encodes Yukawa
couplings, hypercharge, gauge group structure, and generation count.
Three structural differences distinguish Track NI:

- **Substrate.** SM ACG has $M \times F$ with $F$ finite-dimensional and
  internal (one Hilbert space per fermion species per generation);
  Track NI has $M_e \times M_N$ with both factors infinite-dimensional
  and real-space (one Hilbert space per real-space register, with each
  register itself containing internal labels for spin/isospin).
- **Coupling.** SM ACG's cross-block coupling is the Yukawa Dirac block
  on the finite factor (Paper 18 §IV.6 inner-factor Dirichlet ring);
  Track NI's cross-block coupling is the hyperfine
  $\mathbf{I} \cdot \mathbf{S}$ operator on the spin sectors of two
  real-space registers, which is a *spatial* spectral-triple
  construction in Paper 18's outer-factor sense.
- **Calibration.** SM ACG's parameters (Yukawas, hypercharge, mixing
  angles) are calibrated to particle-physics observables; Track NI's
  parameters ($A_{\mathrm{hf}}$, $r_p$) are calibrated to atomic-physics
  observables (hyperfine spectroscopy + electron–proton mass ratio).

These differences are structural, not just cosmetic. They place Track NI
in a different region of the spectral-triple landscape than Connes' SM
construction, even though both are *almost-commutative* in the broad
sense (one finite-dimensional or finite-truncated factor combined with
another). In Paper 18's six-tier exchange-constant taxonomy, Track NI's
parameters $A_{\mathrm{hf}}$, $r_p$ live in the *outer* calibration tier
(spatial scale of the proton, magnetic moments calibrated empirically),
while SM ACG Yukawas live in the *inner-factor input data* tier. The
two tiers tensor cleanly under the AC factorization theorem (Paper 18
§IV.6 Theorem thm:ac_factorization) but neither determines the other.

---

## 5. Comparison to spectral-truncation convergence theorems

The spectral-truncation convergence machinery — Connes–van Suijlekom 2021
(Commun. Math. Phys. 383, 2021); Hekkelman–McDonald 2024
(arXiv:2412.00628); Leimbach–van Suijlekom 2024 (Adv. Math. 439, 109496);
Farsi–Latrémolière 2024 (Adv. Math. 437, 109442); Latrémolière 2026
(arXiv:2603.19128) — proves convergence theorems for *single* spectral
triples (the Connes–vS truncated $C(S^1)$ on the circle, generalized) or
for almost-commutative products with *one* finite factor (e.g.
Latrémolière 2026 explicitly handles "products of the canonical spectral
triple of a compact connected spin manifold with a finite-dimensional
spectral triple"). Paper 38 closes the GeoVac single-factor case,
extending these theorems from the Toeplitz $S^1$ and torus cases to the
Camporesi–Higuchi spectral triple on $S^3$.

What is *not* yet handled in this published lineage is the convergence
of a tensor product of *two infinite metric* spectral triples. Track NI
is the empirical counterpart of that gap: it constructs a finite-truncation
composed object with two infinite real-space factors at $n_{\max} = 2$
that is internally consistent (qubit encoding works, Pauli decomposition
closes, hyperfine validation passes) and is the natural object whose
convergence is the open NCG question. We do *not* prove convergence of
the truncation $\mathcal{T}_{n_{\max}}^N \otimes \mathcal{T}_{n_{\max}}^e$
to the continuum tensor product as $n_{\max} \to \infty$; that
convergence theorem is the W2b cross-manifold blocker of Paper 32
§VIII.D, which is itself genuinely open in the published NCG literature.
Closing that convergence theorem is roadmapped as a possible Phase C
sprint (`debug/multifocal_phase_b_synthesis_memo.md` Section 3) but is
not within the scope of this memo.

The structural reading is therefore that Track NI is a finite-truncation
*example* of the kind of object whose convergence is open in NCG; it is
not itself a convergence proof.

---

## 6. Limitations and explicit non-claims

This memo is scope-honest about what it does and does not establish.

**Limitations.**

- *Known sign issue.* The off-diagonal Pauli encoding for $a^\dagger_i
  a_j$ on the deuteron register (Paper 23 §VI.4) carries a known sign
  issue using $XY - YX$ instead of the expected positive Hermitian
  combination. This affects only the off-diagonal one-body piece of
  the deuteron block, which in the HO basis is identically zero for
  the diagonal kinetic term, so the deuteron FCI spectrum and the
  hyperfine validation are unaffected. Documented but not yet patched.
- *Multi-scale precision problem.* The 13-orders-of-magnitude
  coefficient hierarchy (~10⁵ MeV nuclear, ~10⁰ Ha electronic, ~10⁻⁷ Ha
  hyperfine; Paper 23 §VI.3) makes single-pass quantum simulation
  impractical on near-term hardware. Block-partitioned solving is the
  operative regime.
- *No convergence proof.* Track NI is a finite-truncation construction
  at $n_{\max} = 2$; the tensor-product spectral-triple convergence
  theorem is W2b in Paper 32 §VIII.D and is open in the broader
  literature.
- *No spatial cross-register coupling.* The cross-register operator
  $D_{\mathrm{cross}}$ is hyperfine $\mathbf{I} \cdot \mathbf{S}$ only;
  no spatial $V_{eN}(\hat{\mathbf{r}}_e, \hat{\mathbf{R}}_N)$ recoil
  coupling is implemented. This is the W1a wall of Paper 32 §VIII.D and
  is not closed by the present construction.
- *No magnetization-density operator.* The proton register is point-like
  in its spatial coordinate; no operator on the proton register
  represents the Zemach radius $r_Z$. This is the W1b wall and is also
  not closed by the present construction.

**Explicit non-claims.**

- This memo does **not** claim Track NI is the *first* such construction
  in any literature. The published NCG ecosystem search (10 broad
  queries + 4 deep arXiv WebFetches in Phase B-position T1, May 2026)
  did not surface a directly comparable construction, but a paper in a
  chemistry-physics borderline venue (J. Math. Chem.; Few-Body Systems;
  Int. J. Quantum Chem.) cannot be ruled out. The title and abstract
  use "an explicit … construction" rather than "the first."
- This memo does **not** claim Track NI is *production-ready*. The
  multi-scale precision problem and the known sign issue both stand;
  Track NI is a proof-of-concept of an architecture, not a deployable
  quantum simulation.
- This memo does **not** prove a *convergence theorem*. W2b is open in
  the broader spectral-action / NCG literature.
- This memo does **not** close the *W1a recoil* or *W1b magnetization*
  walls. The Track NI cross-register operator is hyperfine-only.
- This memo does **not** *derive* any new physics. The empirical
  validation (21 cm gap = $3 A_{\mathrm{hf}}/4$) is a consistency check
  on the encoded Hamiltonian, not a prediction.

The originality observation is the literature gap and the empirical
validation is the consistency check. Both are honest.

---

## 6.1 Update (2026-05-08): operator-level closure of W1a and W1b

This memo was first written before the operator-level closures of the
two cross-register walls W1a (cross-register V_eN at distinct focal
lengths) and W1b (magnetization-density inner-fluctuation) had landed.
Phase~C of the multi-focal sprint sequence has since closed both at
the leading order. The memo's scope-honest framing — that Track NI as
constructed in Paper 23 §VI has hyperfine-only cross-register coupling
and does not by itself close the W1a or W1b walls — is unchanged. The
new operator-level constructions sit alongside Track NI's hyperfine
cross-register operator as architecturally compatible additions; they
do not retroactively change the present memo's empirical content
(26-qubit register, 614 non-identity Pauli terms, 21\,cm hyperfine
gap).

The two new operators are documented in Paper 23 §VII (added in the
same update that documents this paragraph) and are reproduced for
reference here:

- **W1a operator-level closure.** A cross-register Coulomb operator
  $V_{eN}(\hat{\mathbf{r}}_e, \hat{\mathbf{R}}_p)$ on the joint
  electron–proton Sturmian register, built on the Roothaan~1951 closed
  form $J_0(\lambda_e, \lambda_n) = \lambda_e \lambda_n (\lambda_e^2 +
  3 \lambda_e \lambda_n + \lambda_n^2)/(\lambda_e + \lambda_n)^3$ for
  the bilinear matrix element of $1/|\mathbf{r} - \mathbf{R}|$
  (`geovac/cross_register_vne.py`, ~1330 lines, 56 fast tests). At the
  calibrated quantum-motional focal length $\lambda_n = 2\sqrt{M_p/m_e}
  \approx 85.7$, the leading-order Bethe–Salpeter recoil shift is
  reproduced verbatim; the residual $2.86\%$ against the full
  Bethe–Salpeter expansion was subsequently identified by symbolic
  Taylor analysis as a basis-truncation artifact of the $1s \times 1s$
  Sturmian representation at $n_{\max} = 1$ (half-integer powers of
  $m_e/m_p$ that are structurally absent from the
  Pachucki–Patkóš–Yerokhin 2023 Foldy–Wouthuysen-reduced two-particle
  Hamiltonian). Multi-shell $n_{\max} \geq 2$ on both registers is the
  named structural follow-up.

- **W1b operator-level closure.** A magnetization-density
  inner-fluctuation $\omega_{\mathrm{magn}}$ on the joint register,
  parameterized by a unit-normalized radial profile $\rho_M(r;\,r_Z)$
  (`geovac/magnetization_density.py`, ~480 lines, 27 fast tests). The
  leading multipole $L = 0$ matrix element collapses to the first
  radial moment $M_1[\rho_M] = r_Z$ automatically, recovering the
  Eides Tab.~7.3 leading-order Zemach shift at $0.012\%$ residual
  ($-39.495$\,ppm vs.~$-39.500$\,ppm, residual consistent with rounding
  of the bohr–fermi conversion constant). The construction realizes
  the Phase~B-W1b-diag verdict ("W1b is downstream of W1a, sharing
  operator infrastructure with one additional Layer-2 calibration
  scalar"): the proton charge-density profile and the magnetization
  density both enter the same diagonal-density Pauli encoding; only
  the radial kernel changes ($-Z/|\mathbf{r}_e - \mathbf{R}_p|$ vs.
  $\rho_M$), with the calibration scalar shifting from $R_p$ to $r_Z$.

The structural realization is that $\omega_{\mathrm{magn}}$ is W1a's
inner-fluctuation sibling on the same composed triple, in direct
analogy with the Connes–Chamseddine decomposition $\omega =
\omega_{\mathrm{gauge}} + \omega_{\mathrm{Higgs}}$ on the SM
almost-commutative spectral triple. Each component is built from a
fixed parameterization of the proton's internal structure ($\rho_E$
for $\omega_{\mathrm{recoil}}$, $\rho_M$ for $\omega_{\mathrm{magn}}$);
the framework computes the operator action on the joint Hilbert space;
the choice of profile family and calibration scalars is Layer-2 input
from QCD/scattering data, in the sense of Paper 18 §IV's inner-factor
input data tier. This is also the structural-skeleton scope statement
of CLAUDE.md §3.5: the framework computes the operator structure, and
the calibration data ($r_Z$, $R_p$, $\lambda_n$ choice, $\rho_M$ shape)
enters as Layer-2 input.

The operator-level closure of W1a and W1b does **not** close the W2b
cross-manifold convergence theorem that surrounds Track NI. W2b is the
two-infinite-factor tensor-product propinquity convergence question
(Paper 32 §VIII.D); the present memo's structural observation — that
Track NI is a finite-truncation example of an object whose convergence
is open in the published NCG literature — is unchanged. A parallel
sprint (Paper 39, in preparation, 2026-05-08) is targeting that
convergence theorem directly.

---

## 7. Conclusion

Track NI is a proof-of-concept of an explicit Connes-style real-space
multi-particle spectral triple. The structural observation it grounds —
that the canonical noncommutative-geometry / spectral-triple lineage
does not contain a directly comparable real-space multi-particle
construction — is the durable contribution of this memo. The empirical
content (26 qubits, 614 non-identity Pauli terms, $1.62 \times 10^{-7}$
Ha hyperfine gap) is honest about its scope: it is a calibrated
consistency check, not a prediction.

The W2b cross-manifold convergence theorem that surrounds Track NI is
roadmapped as a possible Phase C sprint
(`debug/multifocal_phase_b_synthesis_memo.md` Section 3). If closed, it
would convert Track NI from "an explicit construction in a region of
the spectral-triple landscape that is otherwise empty" into "an
explicit instance of a tensor-product spectral triple whose convergence
is established." The present memo does the prior — it lands the
literature observation at Zenodo-memo scale, with all caveats spelled
out.

---

## Cross-references

- **Paper 23 §VI** (the empirical content): composed nuclear–electronic
  deuterium Hamiltonian, qubit construction, hyperfine validation,
  resource counts, known limitation.
- **Paper 32 §V** (sub-sector identification): Track NI is the explicit
  composed real-space construction discussed in the
  "Different $D$, same $\mathcal{A}$" row of Table tab:subsectors and
  the cross-reference paragraph following it.
- **Paper 32 §VIII.D** (frontier-of-field framing): the W2b
  cross-manifold blocker of which Track NI is the empirical example.
- **Paper 18 §IV.6** (inner-factor input data tier): the structural
  separation between outer-factor and inner-factor calibration content
  that distinguishes Track NI from the SM almost-commutative geometry.
- **Paper 38** (WH1 PROVEN): closes the single-factor case
  (Camporesi–Higuchi on $S^3$) of the convergence theorem whose
  tensor-product extension (W2b) is the open question Track NI lives
  inside.
- **Paper 14 §IV** (composed qubit architecture): the qubit-counting
  and Pauli-decomposition machinery that Track NI inherits.
- **Paper 17 §VI** (composed natural geometries): the architectural
  pattern (separate registers per particle group) that Track NI
  instantiates for nuclear + electronic species.
- **CLAUDE.md §1.7 WH1, WH4** (working hypotheses register): the
  spectral-triple framing context within which this observation lives.

## Sprint provenance

- **Phase A multi-focal-composition synthesis** (May 2026): identified
  W2b cross-manifold blocker as the structural setting for Track NI.
- **Phase B-position memo** (May 2026,
  `debug/multifocal_b_position_memo.md`): T1 literature search (10
  queries + 4 deep WebFetches) confirmed the published NCG / spectral-triple
  lineage does not contain a directly comparable real-space
  multi-particle construction; T2 publication-tier recommendation
  selected option (a) Zenodo memo over standalone Paper 39.
- **Phase B synthesis** (May 2026,
  `debug/multifocal_phase_b_synthesis_memo.md`): consolidated panel
  result + Phase C dispatch plan including this memo as a small
  parallel-track artifact.
- **Phase C-positioning** (May 2026): paste of Paper 32 §VIII.D
  frontier-of-field framing, Paper 18 §IV.6 second-packing-axiom
  open-question paragraph, Paper 23 §VI Track NI cross-reference,
  Paper 32 §V cross-reference, this Zenodo memo, CLAUDE.md WH4
  deflation update, and §2 sprint outcome update.

---

**End of Track NI Zenodo memo.**
