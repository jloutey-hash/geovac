# Multi-Focal Audit: Tool Census (Phase A)

**Date:** 2026-05-07
**Author:** Sub-agent (audit pass)
**Scope:** Per-tool catalogue of existing GeoVac machinery that could enable
multi-focal composition. Multi-focal here means: making operators couple across
different focal lengths $p_0 = \sqrt{-2E}$ — different orbital sectors,
different particle masses, different UV scales — within a single Hamiltonian or
spectral triple. The multi-focal seam is where the framework currently shows
five independent walls (Sprint H1, LS-8a, HF-3, HF-4, HF-5; see
`memory/multi_focal_wall_pattern.md`).

**Methodological note.** A tool that handles N focal lengths *in parallel*
(distinct, non-interacting blocks) is not the same as a tool that genuinely
*composes operators across focal lengths*. Several deployed pieces of
machinery do the former; almost none do the latter. The audit is honest about
this distinction.

---

## 1. Recoupling: 3j / 6j / 9j / 12j / Racah

**Lives in:** Paper 14 §III.D (recoupling-sparsity table); Paper 17 (composed
ERI builder); `geovac/composed_qubit.py::_wigner3j`; `geovac/dirac_matrix_elements.py`
(Wigner symbols for jj-coupled angular coefficient $X_k$); Track DD spin-tensor
identity in Paper 14 §V.

**What it currently does.** Wigner 3j enforces angular-momentum selection at
every vertex of every two-body operator; 6j performs basis transformations
between coupled and uncoupled angular bases (proven sparsity gain 18-26% for
H-set vs. uncoupled at 4 electrons, Paper 14 Table); 9j and Racah identities
underlie Drake's 6j-derived spin-tensor patterns (Sprint 4 DD, $f_{SS}/f_{SOO}$
J-dependence). The 12j shows up implicitly in the four-vertex two-loop QED
construction.

**Cross-focal redeployment possibility.** Recoupling is the natural language
for *cross-block angular composition* whenever the blocks share an angular
manifold. Concrete extension: define a "cross-focal 9j" that couples
$(L_A, L_B; L_{AB})$ at focal length $p_0^{(A)} = Z_A/n_A$ to
$(L_C, L_D; L_{CD})$ at focal length $p_0^{(C)} = Z_C/n_C$ via a recoupling
intertwiner whose radial part is the Shibuya-Wulfman or hypergeometric Slater
overlap (which already handles the radial focal-length mismatch). This is not
just "evaluating 6j at two values of l" — it is using the 9j *as the*
intertwiner from one focal sector's coupled basis to another's, with the radial
focal-length mismatch absorbed into the bracketing.

**Already deployed for cross-focal?** Yes-but-only-in-case-X. The composed
ERI builder uses 3j at each block individually; cross-block ERIs (Track CD
balanced coupled) use 3j inside multipole expansion at each center. Recoupling
is invoked for sparsity and selection, not as an inter-focal intertwiner.

**Additional infrastructure needed.** A small lemma generalizing the Drake
6j Racah identity to mixed-focal-length spin tensors: given two blocks at
$p_0^{(A)} \neq p_0^{(B)}$, write $V_{ee}$ in a doubly-coupled basis where
the angular recoupling factors out of the radial focal-length integral. This
is a direct algebraic exercise — not a research bet — and would unify Track CD
balanced coupled, Track NI hyperfine, and a hypothetical recoil V_ne under one
recoupling formalism.

**Speculative class:** Strong.

---

## 2. Hopf base / fiber decomposition

**Lives in:** Paper 25 (whole paper); `geovac/hopf_bundle.py` (`fock_chi`,
`state_to_hopf_angles`); Paper 2 K = π(B+F−Δ) construction.

**What it currently does.** Decomposes S³ into S² base × S¹ fiber via Hopf
fibration. Maps each lattice state $(n, l, m)$ to explicit (χ, ψ₁, ψ₂)
coordinates via Fock projection. Quotient by S¹ produces the 6-sector S²
graph used in Paper 2's α conjecture. U(1) gauge phase lives on the fiber.

**Cross-focal redeployment possibility.** The Hopf fibration is the cleanest
existing vocabulary for *separating fast and slow degrees of freedom* on the
same manifold. If two focal lengths can be identified with two natural
projections of the same Hopf bundle (say, base and fiber, or two different
sub-bundles), then operator composition across them is the standard
fibre-bundle pull-back/push-forward. Concretely: the electron-1s focal length
($p_0 = Z$) and the proton finite-size focal length ($p_0 \sim 1/r_p$) might
both project from a shared S³ on which one acts on the base and the other on
the fiber. This is speculative but structurally clean.

**Already deployed for cross-focal?** No. Currently used only as a
mathematical object inside one focal sector (Paper 2 α conjecture; Paper 25
gauge structure on a single graph).

**Additional infrastructure needed.** A construction identifying which Hopf
bundle naturally hosts a given two-focal observable; the matter sector of the
proton-electron tensor product must be recognized as living on a definite
Hopf base + fiber, not just two registers tensored.

**Speculative class:** Moderate.

---

## 3. Wigner-D matrix rotation (multi-center)

**Lives in:** `geovac/shibuya_wulfman.py::_build_rotation_matrix_real_sh`,
`_build_block_rotation_matrix`; Paper 19 §V; Track CM (l=2 D-matrix, v2.1.1);
multi-center Sprints CU-CX.

**What it currently does.** Builds (2l+1)×(2l+1) real-spherical-harmonic
rotation matrices via the standard Wigner small-d formula, converted from
complex SH basis through a fixed unitary $U$. Block-diagonal assembly across
mixed-l blocks. Used to rotate cross-center V_ne computed in the z-frame to
arbitrary 3D directions for non-collinear molecules (H₂O, CH₂O, etc.).

**Cross-focal redeployment possibility.** Wigner-D's *natural argument* is
an SO(3) rotation, but the same machinery generalizes to SO(4) via the SU(2)
× SU(2) double cover — and the Fock projection IS an SO(4) action on S³.
Concrete extension: replace the SO(3) rotation by an SO(4) rotation that maps
focal length $p_0^{(1)}$ to $p_0^{(2)}$ on S³ (a non-trivial isometry of the
Fock-projected sphere). This would let one express a "focal-length boost"
operator as an SO(4) rotation, with matrix elements computable by the same
Wigner-D framework but in the 4D rep.

**Already deployed for cross-focal?** Yes-but-only-in-case-X. The current
deployment handles cross-*center* (different spatial origins) but not
cross-*focal* (different binding energies). Spatial center and binding-energy
focal length are different geometric data.

**Additional infrastructure needed.** SO(4) D-matrices replacing SO(3), and
an explicit identification of the "focal-length boost" as a one-parameter
subgroup (probably the dilatation generator, which is part of the SO(2,1)
conformal subgroup of SO(4,2)). The SU(2) × SU(2) Avery-Wen-Avery 3-Y integral
in `so4_three_y_integral.py` already handles SO(4) 3-Y integrals — the rotation
infrastructure would parallel that.

**Speculative class:** Moderate.

---

## 4. Almost-commutative inner factor / D_F

**Lives in:** `geovac/almost_commutative.py`; Paper 32 §VIII.C; Sprint H1
memo `debug/h1_ac_extension_memo.md`; `inner_factor_mellin_engine.md` memory.

**What it currently does.** Tensors GeoVac with finite spectral triple
$(A_F, H_F, D_F, J_F)$ realizing electroweak data. $D = D_{GV} \otimes 1_F +
\gamma_{GV} \otimes D_F$. The off-diagonal $D_F$ block carries Yukawa data;
matter/antimatter doubling makes the order-zero condition hold automatically.
Three Connes axioms verified exactly (truthful CH).

**Cross-focal redeployment possibility.** This is the *closest existing
analog of a multi-focal composition theorem*. Tensoring $D_{GV}$ at one focal
length with an internal $D_F$ at a categorically different scale — the
electroweak / Yukawa scale — is exactly cross-focal composition, and the
spectral-action machinery handles it (see η-trivialization theorem +
factorization theorem in `inner_factor_mellin_engine.md`). The catch: the
inner factor's data $D_F$ is *prescribed*, not derived from $D_{GV}$, so the
framework does not autonomously generate Yukawas. But for any *known* second
focal length, the AC extension is the right machinery.

**Already deployed for cross-focal?** Yes — but only at the structural-
skeleton level (Pauli sparsity factorization, KO-dim addition, η-triviality
on the inner factor). Not deployed as a multi-focal operator on physical
observables (Higgs gap, etc.).

**Additional infrastructure needed.** A finite spectral triple representing
the *physical* second focal length one wants — for recoil, that would be a
small two-state $H_F$ with a $D_F$ encoding the proton kinetic operator (at
nuclear focal length $p_0 \sim m_p \cdot v$); for Zemach, $D_F$ encoding the
magnetization radial structure. The construction is mechanical once $D_F$ is
specified; the question is whether GeoVac structure picks $D_F$ uniquely or
demands it as input. Sprint H1 verdict says it is input, not derived.

**Speculative class:** Strong (machinery exists, deployment is mechanical;
but the *choice* of $D_F$ is the bottleneck identified by H1).

---

## 5. Coulomb-Sturmian basis at λ ≠ Z/n

**Lives in:** Papers 8-9 (Sturmian theorem); Paper 19 (Sturmian for
molecules); `geovac/sturmian_solver.py`; `geovac/molecular_sturmian.py`;
Track BU (multi-electron Sturmian CI, Paper 14 §IV.M).

**What it currently does.** Sturmian basis at a single fixed λ (the
"Sturmian momentum") gives orbitals at that λ a complete basis with H ∝ S
(structural theorem). The molecular Sturmian uses prolate spheroidal
separation. Track BU's negative result on multi-electron Sturmian in molecules
showed Lowdin orthogonalization inflates 1-norm 2.8-4.5×.

**Cross-focal redeployment possibility.** This is the *most natural existing
multi-focal tool*. The Sturmian basis at λ ≠ Z/n is *literally* an orbital at
the wrong focal length — Papers 8-9 build the entire framework around this
freedom. The Lamb shift sprint (LS-3) showed the Sturmian basis at λ = Z/n
*is* the Fock graph reparameterized for bound-state space, which means
choosing λ = different value gives the *same Fock graph* viewed from a
different focal length. Multi-focal composition would be: take two basis sets
at $\lambda_1$ and $\lambda_2$, both spanning the same Hilbert space, and
write a Hamiltonian whose blocks live at different λ.

**Already deployed for cross-focal?** No — single-λ at every deployment
point, even when the molecule has natural multi-scale structure (LiH does).
Track BU's negative result was for *constructing wavefunctions* with mixed-λ;
the multi-focal *composition* angle was not explored.

**Additional infrastructure needed.** A "transfer operator" between Sturmian
bases at different λ. This is an algebraic object (the matrix elements
$\langle n,l,m;\lambda_1 | n',l',m';\lambda_2\rangle$ have closed form in terms
of hypergeometric functions). Once that's available, composing operators
across λ is mechanical.

**Speculative class:** Strong. This may be the highest-yield tool in the
audit because the radial focal-length mismatch is its *raison d'être*.

---

## 6. Shibuya-Wulfman two-center integrals

**Lives in:** `geovac/shibuya_wulfman.py`; Paper 19 §III-IV; Track CD balanced
coupled.

**What it currently does.** Computes
$\langle \psi_{nlm}^A | -Z_B/|r-R_B| | \psi_{n'l'm'}^A\rangle$ via multipole
expansion of $1/|r-R_B|$ in $P_L(\cos\theta)$, terminating exactly at $L_{max}
= 2 l_{max}$ by Gaunt selection rules. Analytical evaluation via incomplete
gamma functions (machine precision).

**Cross-focal redeployment possibility.** The current deployment couples two
*spatial* centers at the *same* $Z_{orb}$ (the orbital wavefunctions are at
focal length $p_0 = Z_{orb}/n$, regardless of which nucleus the V_ne points
at). To make it multi-focal, generalize to bra and ket at *different* focal
lengths: $\langle n_A, l_A, m_A; p_0^{(A)} | -Z_B/|r-R_B| | n_C, l_C, m_C;
p_0^{(C)}\rangle$. The radial integral becomes a generalized hypergeometric
function (still elementary; the polynomial product just has two exponential
decay rates). The angular machinery is unchanged.

**Already deployed for cross-focal?** Yes-but-only-in-case-X. Current
deployment is single-focal on bra/ket but multi-center on potential. Two-focal
on bra/ket is a small extension.

**Additional infrastructure needed.** Replace `_hydrogenic_poly_coeffs` to
return $(c, \alpha)$ for two different effective charges, then track two
exponential decay rates through the polynomial product. The
`_split_integral_analytical` already handles arbitrary α via incomplete gamma.

**Speculative class:** Strong.

---

## 7. Multipole expansion of 1/|r-R_B|

**Lives in:** Paper 19 §III.C; balanced coupled `geovac/balanced_coupled.py`;
Sprint CD memos.

**What it currently does.** Expands a non-local two-center potential as a
sum of finitely many local multipole terms via the standard Legendre
expansion, and proves $L_{max} = 2 l_{max}$ termination from Gaunt selection.

**Cross-focal redeployment possibility.** Multipole expansion is one of the
*few* tools in the framework that is genuinely *manifold-agnostic*: as soon
as the kernel has a multipole representation, the angular machinery is fixed
by Gaunt/3j and only the radial part depends on focal length. Generalizes
to: any two-body Coulomb-like operator coupling two different focal sectors.

**Already deployed for cross-focal?** Yes-but-only-in-case-X. Used for
cross-center $V_{ne}$ with single-focal orbitals. The honest cross-focal
deployment would be e.g. $V(r_e, r_n)$ in the proton-electron problem, which
*has* a multipole expansion — but the framework currently treats $r_n$
classically (HF-3 negative result).

**Additional infrastructure needed.** The bottleneck is not multipole
machinery itself but the second-quantization side: representing $r_n$ as a
quantum operator on the proton register and expanding $|r_e - r_n|^{-1}$ in
multipoles where *both* coordinates are operator-valued. This is the
nuclear-electronic two-body Coulomb the HF-3 sprint flagged as missing.

**Speculative class:** Strong (machinery), but the missing piece is the
two-register operator structure, not the expansion.

---

## 8. Cross-register tensor product with bilinear coupling

**Lives in:** `geovac/nuclear/nuclear_electronic.py::build_deuterium_composed_hamiltonian`;
Track NI; Paper 23 Sec VI; HF-3 / HF-4 / Sprint HF.

**What it currently does.** Tensors a nuclear register
(O(MeV) HO basis) with an electronic register (O(Ha) hydrogen basis) into one
Hamiltonian; couples them via a *finite-size* one-body shift on the electron
1s number operator and a *hyperfine* I·S spin-spin term. Coefficient ratio
2×10¹³ across nuclear/electronic/hyperfine scales is real and works.

**Cross-focal redeployment possibility.** This IS the most explicit
cross-focal architecture in the framework — two registers at categorically
different energy scales (different focal lengths in the most concrete sense)
sharing a single qubit space. But: the ONLY cross-couplings present are
spin-spin (I·S, scalar parameter) and a classical scalar finite-size shift.
There is NO operator coupling spatial coordinates across registers. Multi-focal
extension: add genuine V_ne(r_e, R_n) as a two-body Coulomb between registers,
with R_n a real operator on the nuclear register (currently it's a classical
constant R_PROTON_BOHR).

**Already deployed for cross-focal?** Yes — for spin couplings. No — for
spatial couplings. This is exactly the wall HF-3 (recoil) and HF-4 (Zemach)
hit.

**Additional infrastructure needed.** A coordinate operator on the HO-basis
nuclear register and a multipole expansion of $|r_e - R_n|^{-1}$ where both
arguments are operator-valued. The HO basis has natural matrix elements
$\langle n_r, l, m_l | r^k | n_r', l', m_l'\rangle$ from Moshinsky-Talmi
(see `geovac/nuclear/moshinsky.py`); these would feed the multipole expansion
of Tool 7.

**Speculative class:** Strong. The architecture is *already* multi-focal
in the registers; the missing piece is wiring spatial coupling. This is
arguably the highest-leverage point for genuine multi-focal composition.

---

## 9. Bargmann transform / Bargmann-Segal lattice

**Lives in:** Paper 24; `geovac/nuclear/bargmann_graph.py`.

**What it currently does.** Maps $L^2(\mathbb{R}^3)$ to Fock-Bargmann
$F^2(\mathbb{C}^3)$ unitarily. HO eigenstates become solid-harmonic
polynomials of degree N. Restricted to $S^5 \subset \mathbb{C}^3$, gives the
Hardy-space sector with SU(3)-symmetric (N,0) irreps. The graph is bit-exactly
π-free.

**Cross-focal redeployment possibility.** Bargmann is a *transform*
(unitary equivalence) that takes one Hilbert space to another — and the natural
"second focal length" is exactly the holomorphic-vs-Riemannian dichotomy
(Paper 24 §IV-V). Cross-focal redeployment: use the Bargmann transform as a
basis-change between a Coulomb-S³ block (Riemannian, second-order) and a
HO-S⁵ block (holomorphic, first-order) inside the same Hamiltonian. This is
the G4b cross-manifold spectral-triple problem (Paper 32 §VIII open question)
in computational disguise.

**Already deployed for cross-focal?** No. Bargmann is used only in the
HO-only construction; never as a transform between Coulomb and HO sectors of
the same wavefunction.

**Additional infrastructure needed.** A formal Coulomb→HO transform (the
Cooper-Ginocchio-Wipf-style canonical map sending hydrogen radial states to
HO radial states via Levi-Civita regularization) and recognition that this
is structurally a Bargmann-like unitary intertwining two focal-length sectors.
Possibly out of scope at the spectral-triple level (Paper 32 §VIII.C, G4b).

**Speculative class:** Weak. Mathematically beautiful but the framework
already flags the cross-manifold gap as deferred.

---

## 10. Berezin reconstruction / Berezin-Toeplitz quantization

**Lives in:** `geovac/berezin_reconstruction.py`; Paper 38 (R2.5 L4); Sprint
WH1-R2.5.

**What it currently does.** Constructs $B_{n_{\max}}: C(S^3) \to O_{n_{\max}}$
as $B_{n_{\max}}(f) = P_{n_{\max}}(K_{n_{\max}} * f) P_{n_{\max}}$ — a positive
contractive approximate-identity map from continuous functions on S³ to the
truncated operator system. Used in WH1 PROVEN (GH convergence five-lemma
chain).

**Cross-focal redeployment possibility.** Berezin is a *quantization map*
that takes a classical observable to its quantum compression, and Berezin-
Toeplitz (the calling pair $(B, P)$) is the canonical bridge between two
spectral-triple resolutions. Cross-focal redeployment: at two different
focal lengths $p_0^{(1)}, p_0^{(2)}$ one obtains two different truncations of
the same continuum theory; the Berezin map transports operators between them
via the continuum as a common reference. This is exactly the
"projection-as-tunneling-pair" perspective from Paper 38.

**Already deployed for cross-focal?** No — single focal length with
$n_{max} \to \infty$. But the structural ingredient *is* a multi-resolution
bridge.

**Additional infrastructure needed.** A two-focal-length variant
$B_{n_{\max}^{(1)}, n_{\max}^{(2)}}$ identifying when the two truncations are
at *different* focal lengths rather than the same focal length at different
cutoffs. The existing Latrémolière propinquity machinery in
`geovac/gh_convergence.py` accommodates this, conceptually.

**Speculative class:** Moderate. Beautiful infrastructure; the multi-focal
deployment requires a non-obvious identification of "focal length" with a
metric-spectral-triple parameter.

---

## 11. Master Mellin engine k = 0, 1, 2

**Lives in:** Paper 32 §VIII (case-exhaustion + master mechanism + domain
partition remarks); Paper 18 §III.7; Sprint TS-E1, TS-E3, MR-A/B/C; memory
`mellin_engine_domain_partition.md`.

**What it currently does.** Classifies all π-sources in any finite chain of
Paper 34 projections as Mellin transforms $\mathcal{M}[\mathrm{Tr}(D^k
\cdot e^{-tD^2})]$ at $k \in \{0, 1, 2\}$ (k=0 Hopf-base, k=1 vertex-parity
Hurwitz, k=2 Seeley-DeWitt). Mechanism-and-domain partition: the index k
classifies BOTH the M-mechanism AND the natural class of observable.

**Cross-focal redeployment possibility.** The Mellin transform is *the*
standard operator-theoretic device for combining contributions across scales.
The k=0,1,2 partition naturally extends to a *focal-length-resolved* Mellin
where the trace runs over a multi-focal Hilbert space and the heat kernel
$e^{-tD^2}$ is replaced by a two-parameter heat kernel $e^{-t_1 D_1^2 -
t_2 D_2^2}$. The Mellin transform of this two-time heat kernel encodes
cross-focal mixing.

**Already deployed for cross-focal?** No — single-focal at every deployment.
But the inner-factor Mellin engine theorem
(`memory/inner_factor_mellin_engine.md`) is a step in this direction:
$D^2 = D_{GV}^2 \otimes 1 + 1 \otimes D_F^2$ produces a two-factor heat
kernel, and the Mellin engine on each factor is independent.

**Additional infrastructure needed.** A genuine two-time Mellin transform
on Hilbert spaces of differing focal-length. This is mathematically standard
(noncommutative integration with two regulators); the GeoVac-specific
ingredient would be identifying which physical observables produce which
two-time signatures.

**Speculative class:** Moderate-Strong (the inner-factor theorem already
half-deploys this).

---

## 12. Truncated operator system O_{n_max}

**Lives in:** `geovac/operator_system.py`; Paper 32 §III; WH1 Round 1+2
(prop = 2 Toeplitz match).

**What it currently does.** Realizes the Connes-vS spectral truncation
$O_{n_{\max}} = P_{n_{\max}} C^\infty(S^3) P_{n_{\max}}$ as a *-closed but
not multiplicatively-closed linear subspace; verifies prop = 2 (matching
Toeplitz S¹).

**Cross-focal redeployment possibility.** Operator systems are designed for
*compositions that are not algebra homomorphisms* — exactly what one expects
for cross-focal composition where the multiplicative structure does not
descend cleanly. A multi-focal operator system would be
$P_{n_{\max}^{(1)}} \otimes P_{n_{\max}^{(2)}} \cdot C^\infty(S^3 \times S^3)
\cdot P_{n_{\max}^{(1)}} \otimes P_{n_{\max}^{(2)}}$. Propagation number,
witness pairs, dim-fraction analyses all generalize.

**Already deployed for cross-focal?** No — single-focal on a single S³.

**Additional infrastructure needed.** A two-focal operator system class
extending the existing one. Computationally cheap; the question is whether
it adds physics or remains a categorical exercise.

**Speculative class:** Moderate.

---

## 13. Real structure J at finite n_max

**Lives in:** `geovac/real_structure.py`; Paper 32 §IV; WH1-Connes Step 2.

**What it currently does.** Constructs antilinear J = U·conj on the
Camporesi-Higuchi spinor bundle, verifies $J^2 = -I$, $JD = +DJ$ exactly at
finite $n_{\max}$.

**Cross-focal redeployment possibility.** $J$ is the *charge conjugation*,
which already commutes between matter and antimatter sectors in the AC
extension (Tool 4). Generalizing to "focal-length conjugation" would mean a
$J$ that interchanges two focal-length sectors with appropriate phase. The
$J_F = \sigma_x$ swap between matter and antimatter sectors is exactly such
a sector-swap.

**Already deployed for cross-focal?** Yes — internally to the AC factor
(matter ↔ antimatter sector swap is a categorical "two-focal" object). Not
deployed for two physical focal lengths.

**Additional infrastructure needed.** Identification of physical focal-length
pairs that admit a real-structure swap. This is a stronger requirement than
operator-system compatibility.

**Speculative class:** Weak. Real structure is a strong constraint and
restricts which focal-length pairs admit the symmetry.

---

## 14. Wilson lattice gauge

**Lives in:** Papers 25, 30; Sprint ST-SU3; `geovac/su2_wilson_gauge.py`,
`geovac/su3_wilson_s5.py`.

**What it currently does.** Three Wilson lattice gauge constructions on
GeoVac sub-manifolds: U(1) on S³ Hopf graph (Paper 25), SU(2) on S³ = SU(2)
Coulomb graph (Paper 30), SU(3) on S⁵ Bargmann graph (ST-SU3, gauge-only).

**Cross-focal redeployment possibility.** A Wilson link variable $U_e$ is
naturally a *transport* between the two endpoints of an edge. If two
endpoints sit at *different focal lengths*, the link variable is a focal-
length-mixing transport. This is a small extension of standard Wilson lattice
gauge (where endpoints are "the same kind of thing"). Concretely: dress the
GeoVac graph with edges connecting Coulomb-S³ shells to HO-S⁵ shells, and
put a Wilson link on each. The action is gauge-invariant by construction.

**Already deployed for cross-focal?** No — each of the three Wilson
constructions sits on a single sub-manifold.

**Additional infrastructure needed.** A graph that genuinely spans two focal
sectors (not just two replicas of the same sector). The Bargmann transform
(Tool 9) might supply the inter-focal edges.

**Speculative class:** Weak-Moderate. Conceptually clean but presupposes a
cross-focal graph that the framework does not yet construct.

---

## 15. Frozen-core machinery (screened Z_eff(r))

**Lives in:** `geovac/neon_core.py`; Paper 17 §III; Sprint 7b screened SO.

**What it currently does.** Bridges core electrons (compact, high-Z focal
length) and valence electrons (diffuse, low-Z_eff focal length) by analytical
$Z_{\text{eff}}(r)$ from Clementi-Raimondi exponents. Solves a single-electron
Schrödinger equation in the screened potential. Sprint 7b extension: $\langle
1/r^3\rangle$ for spin-orbit from the screened wavefunction.

**Cross-focal redeployment possibility.** This is the *operationally most
deployed* multi-focal tool in the framework — it bridges two focal lengths
(core $p_0 \sim Z$, valence $p_0 \sim 1$) by constructing a single radial
potential that interpolates between them. Generalizes by replacing the
Clementi-Raimondi profile with any radial function bridging two focal lengths
of the user's choice.

**Already deployed for cross-focal?** Yes — in the most direct sense of any
tool in the audit. But: it is a *mean-field* bridge (one-body screened
potential), not a *quantum operator* coupling two focal sectors. So it does
not address the multi-focal-composition wall.

**Additional infrastructure needed.** A quantum (operator-valued, not
classical) version of $Z_{\text{eff}}(r)$ that retains the core sector as
quantum DOF rather than collapsing it to a screening function. This is
exactly the explicit-correlation generalization that the literature calls
"explicit core" (vs. frozen core).

**Speculative class:** Strong (foundation), but the operator-valued upgrade
is a significant infrastructure project.

---

## 16. Composed-block architecture (different Z per block)

**Lives in:** `geovac/composed_qubit.py`; `geovac/molecular_spec.py`; Paper 17.

**What it currently does.** Each `OrbitalBlock` carries its own
$Z_{\text{center}}$ — the focal length of that block's orbitals — and the
total Hamiltonian is assembled by union. For LiH: Li-core block at Z=3,
LiH-bond block at Z_eff=1 (Li-side) + Z=1 (H-side). PK barrier separates
core/valence focal lengths.

**Cross-focal redeployment possibility.** The composed architecture is the
*most-deployed parallel-focal-length* tool: it handles N focal lengths *in
parallel* (different Z per block) but couples them only through (a) PK
projection, (b) cross-block ERIs, (c) cross-center V_ne. The PK is a one-body
operator; cross-block ERIs are the closest thing to cross-focal *two-body*
coupling, but they assume orbitals at different Z share spatial origin or
are linked by Shibuya-Wulfman (Tool 6).

**Already deployed for cross-focal?** Yes-but-only-in-case-X. Different
focal lengths in parallel; cross-focal coupling is restricted to (a) one-body
PK and (b) two-body Coulomb at the same spatial point or via two-center
multipole. There is no "kinetic mixing" between focal lengths at the
composed-block level — that is the wall HF-5 (multi-loop a_e renormalization)
hits.

**Additional infrastructure needed.** Genuine cross-focal *kinetic* terms,
i.e. an operator $T_{AB}$ in the Hamiltonian that takes a block-A orbital and
returns a block-B orbital. The framework has no such operator currently;
adding one would be a major architectural extension.

**Speculative class:** Strong (foundation; it IS the multi-focal architecture
for orbitals), but the cross-focal kinetic gap is real.

---

## 17. Other / supporting tools (briefer)

| Tool | Module | Cross-focal angle | Class |
|:-----|:-------|:------------------|:------|
| Hypergeometric Slater integrals | `hypergeometric_slater.py` | Already handles Slater integrals at arbitrary single Z; mixed-Z generalization (bra at Z₁, ket at Z₂) is a one-page lemma | Strong |
| Moshinsky-Talmi brackets | `nuclear/moshinsky.py` | Standard inter-focal HO machinery (CM ↔ relative coordinates at the same ω); ω-shift generalization handles two HO focal lengths | Moderate |
| Chirality grading γ_GV | `chirality_grading.py` | G3 NEGATIVE (γ_GV and γ_F commute independently). Confirms chirality is *not* a cross-focal handle | Weak |
| Drake-Swainson asymptotic subtraction | LS-4 sprint | Splits a sum at K and re-adds tail; the "intermediate K" is a regulator that generalizes to a "focal-length splitter" between low- and high-energy contributions | Moderate |
| Fock graph Hodge decomposition | `fock_graph_hodge.py` | Single-focal currently; co-exact / pendant-edge structures naturally extend to multi-focal | Moderate |
| Down-folding (DUCC) | `downfolding.py` | Already an inter-focal tool (active vs frozen at different scales) | Moderate |
| Coulomb-Sturmian molecular bond sphere (Papers 8-9) | `sturmian_solver.py` | The Sturmian theorem H ∝ S says eigenvalues at fixed λ are R-independent — a focal-length invariance property that *constrains* multi-focal composition | Weak (constraint, not enabler) |

---

## Summary

### Top 3 tools to test first

1. **Coulomb-Sturmian basis at $\lambda \neq Z/n$ (Tool 5).** This is the
   one tool in the framework whose *purpose* is to decouple basis focal
   length from physical focal length. Track BU's negative result was for
   wavefunctions; the *operator-composition* angle was never tested. The
   transfer operator between Sturmian bases at different λ has closed form
   and would unify several open questions (recoil, Drake regularization,
   Bethe log convergence) under a single multi-focal language.

2. **Cross-register tensor product with a genuine two-body coordinate
   coupling (Tool 8 + Tool 7 + Tool 14 above).** Track NI is *already*
   multi-focal at the register level (electron Ha + nucleus MeV); the wall
   is that V_ne(r_e, R_n) is currently a classical scalar, not a two-body
   operator. Adding an operator-valued $R_n$ on the HO register and using
   the multipole expansion of $|r_e - R_n|^{-1}$ across registers is a
   focused engineering task (machinery exists in pieces) and would close
   HF-3 / HF-4 directly.

3. **Shibuya-Wulfman with mixed-focal bra/ket (Tool 6).** Smallest extension
   per unit potential payoff. Replace the single $Z_{\text{orb}}$ with two
   exponential decay rates in `_hydrogenic_poly_coeffs`, propagate through
   the polynomial product, and the rest of the analytical incomplete-gamma
   pipeline is unchanged. Would let composed-block ERIs across blocks at
   different Z become genuinely cross-focal rather than approximate.

### Cross-tool combinations worth considering

- **Tool 5 + Tool 6**: Sturmian-at-different-λ as bra/ket on a Shibuya-Wulfman
  cross-center potential. Combines the orbital focal-length flexibility of
  Sturmians with the two-center reach of Shibuya-Wulfman. Plausibly the
  cleanest existing way to write down a cross-focal V_ne where bra, ket, and
  potential each carry independent focal data.

- **Tool 8 + Tool 11**: Cross-register tensor product fed into a two-time
  master Mellin engine. The two-factor heat kernel
  $e^{-t_1 D_e^2 - t_2 D_n^2}$ has a Mellin transform whose poles encode
  cross-register transcendental signatures, and the inner-factor Mellin
  theorem (memory: `inner_factor_mellin_engine.md`) is the prototype.

- **Tool 4 + Tool 13**: AC extension with a focal-length-conjugating $J$.
  If the second focal length comes with a natural sector-swap (matter ↔
  antimatter is the prototype), the spectral-triple machinery handles
  two-focal composition almost for free. Speculative but well-defined.

- **Tool 9 + Tool 14**: Wilson links along the Bargmann-induced edges
  between Coulomb-S³ and HO-S⁵ sectors. Closes G4b at the lattice level
  even if the spectral-triple cross-manifold construction (Paper 32 §VIII.C
  open) remains open.

### Honest scope and uncertainties

- **Did not cover** in depth: angular sparsity theorem (Paper 22, more a
  result than a tool); full two-loop QED machinery (Track LS-7 / LS-8a);
  details of the master Mellin engine's k=1 vertex-parity case beyond the
  Paper 28 vertex code; Track DF nested-hyperspherical sparsity result
  (relevant only if cross-focal composition lives on an S^(3N-1) manifold,
  which the failed-approaches table suggests it does not).

- **Surprising finding (one).** The Coulomb-Sturmian basis at $\lambda \neq
  Z/n$ is *exactly* a focal-length-shifted basis on the same Fock graph
  (Paper 36 LS-3 sprint identification: Sturmian basis at $\lambda = Z/n$
  IS the Fock graph reparameterized for bound-state space). I had filed
  Sturmians under "bond-sphere theorem" (Tool 5 felt like a niche
  molecular-method tool) and missed that LS-3 already names it as a
  focal-length-shifted reparameterization. This makes Tool 5 the highest-
  yield tool in the audit, and it is currently *single-focal* in every
  deployment despite its raison d'être being multi-focal.

- **Surprising finding (two).** Frozen-core machinery (Tool 15) is the most
  *deployed* multi-focal tool by a wide margin, but it is mean-field — it
  classicizes one focal sector to make composition tractable. The honest
  multi-focal-composition wall is that the framework *quantizes parallel
  focal lengths* (composed blocks, Tool 16) and *classicizes serial focal
  lengths* (frozen core, Tool 15). The genuine multi-focal seam is between
  these two — a regime where serial focal lengths must remain quantum and
  parallel coupling is non-trivial. None of the existing tools occupy that
  seam directly.

- **Uncertainty.** The classification "Strong / Moderate / Weak" is my
  best guess for engineering-effort and physical payoff combined; a tool
  classified Moderate could be Strong if it turns out the bottleneck for
  a specific multi-focal observable hits exactly its sweet spot. Tools 1,
  5, 6, 8 are the most uniformly Strong across plausible target observables;
  Tool 4 is Strong for spectral-triple targets but tied to the H1 verdict
  on Yukawa data being *input*; Tool 11 is Moderate-Strong but the
  computational cost is unclear.

- **Caveat on "deployed for cross-focal".** I was conservative: a tool
  flagged "Yes" must currently take operators at *physically different*
  focal lengths and produce a coupling. "Yes-but-only-in-case-X" means it
  handles a degenerate case (same Z, multiple centers; same n_max,
  different shells) but not the genuine cross-focal target. "No" means
  single-focal in every deployment.
