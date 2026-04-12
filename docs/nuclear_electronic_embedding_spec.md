# Nuclear-Electronic Embedding Architecture: Theory Specification

**Track ND — Working Document (April 2026)**

---

## 1. Physical Motivation

Standard computational quantum chemistry treats nuclei as classical point charges and solves only for electronic structure. This is the Born-Oppenheimer approximation, and it is excellent for most chemistry. However, several classes of observables require the nuclear wavefunction or nuclear-electronic coupling:

- **Isotope shifts**: The finite nuclear charge radius differs between isotopes ($^{1}$H vs $^{2}$H vs $^{3}$H, $^{3}$He vs $^{4}$He). The resulting shift in electronic energy levels scales as $(Z\alpha)^4 m_e c^2 (R_{\mathrm{nuc}}/a_0)^2$ and has been measured to sub-MHz precision.
- **Hyperfine structure**: The nuclear magnetic dipole moment couples to the electronic angular momentum via $H_{\mathrm{hf}} = A_{\mathrm{hf}} \, \mathbf{I} \cdot \mathbf{J}$, producing MHz-scale splittings (the 21-cm hydrogen line, atomic clock transitions).
- **Nuclear quadrupole coupling**: Non-spherical nuclei ($I \geq 1$) interact with the electric field gradient at the nucleus, producing kHz-scale splittings observable in NMR and microwave spectroscopy.
- **Recoil corrections**: The finite nuclear mass introduces a mass polarization term $\sim (m_e / M_{\mathrm{nuc}}) \, \mathbf{p}_i \cdot \mathbf{p}_j$ between electrons. For hydrogen, this is a $\sim 5 \times 10^{-4}$ relative correction.
- **Nuclear structure observables**: Binding energies, charge radii, and excited state spectra of the nucleus itself are active targets for quantum simulation (Roggero et al., 2020).

The energy scale hierarchy spans roughly nine orders of magnitude:

| Scale | Magnitude | Physics |
|:------|:----------|:--------|
| Nuclear binding | $\sim 1$--$10$ MeV/nucleon | Strong + Coulomb force between nucleons |
| Electronic binding | $\sim 1$--$10^3$ eV | Coulomb force, electron kinetic energy |
| Fine structure | $\sim 10^{-3}$--$10^{-1}$ eV | Spin-orbit coupling, relativistic corrections |
| Hyperfine structure | $\sim 10^{-7}$--$10^{-5}$ eV | Nuclear moment -- electron coupling |
| Finite nuclear size | $\sim 10^{-10}$--$10^{-8}$ eV | Proton charge distribution |

A composed Hamiltonian that encompasses both nuclear and electronic degrees of freedom within a single qubit register would, in principle, give direct access to all of these scales from a unified simulation. The practical question is whether the $\sim 10^6$-fold energy scale separation between nuclear and electronic physics can be handled within the composed block architecture without catastrophic precision loss.

**Connection to GeoVac's existing architecture.** The composed geometry already handles multi-scale physics: core electrons at $\sim 100$ eV are composed with valence electrons at $\sim 1$--$10$ eV, separated by PK pseudopotentials or balanced coupling. Extending this to include a nuclear block at $\sim 10^6$ eV is architecturally natural --- it adds one more tier to the existing hierarchy. The inter-block coupling (nuclear charge form factor) is a one-body operator on the electronic side, analogous to the existing $V_{ne}$ Coulomb terms.

---

## 2. Block Hierarchy

The composed nuclear-electronic Hamiltonian has three blocks, ordered by energy scale.

### Block A: Nuclear Shell

**Fermion content:** $Z$ protons $+ (A - Z)$ neutrons. Protons and neutrons are *distinct* fermion species --- antisymmetrization applies *within* each species but not between them. This is a critical architectural difference from the electronic blocks, where all particles are identical electrons.

**Single-particle basis:** The standard nuclear shell model uses the harmonic oscillator (HO) as the starting point, labeled by quantum numbers $(N, n_r, l, j, m_j)$ where $N = 2n_r + l$ is the HO major shell number. The physical single-particle potential is the Woods-Saxon form:

$$V_{\mathrm{WS}}(r) = \frac{-V_0}{1 + \exp\bigl((r - R_0)/a\bigr)}$$

with $R_0 \approx 1.2 \, A^{1/3}$ fm, $a \approx 0.65$ fm, $V_0 \approx 50$ MeV. A spin-orbit term $V_{\mathrm{so}}(r) \, \mathbf{l} \cdot \mathbf{s}$ splits each $(n_r, l)$ level into $j = l \pm 1/2$ doublets and is responsible for the nuclear magic numbers (2, 8, 20, 28, 50, 82, 126).

**Two-body interaction:** The nucleon-nucleon (NN) force is fundamentally different from the Coulomb interaction:
- Short-range: $\sim 1$--$2$ fm (mediated by pion exchange at long range, quark-gluon physics at short range)
- Spin-dependent: the $^1S_0$ (spin-singlet) and $^3S_1$ (spin-triplet) channels have different scattering lengths
- Tensor force: $V_T(r) \, S_{12}$ where $S_{12} = 3(\boldsymbol{\sigma}_1 \cdot \hat{r})(\boldsymbol{\sigma}_2 \cdot \hat{r}) - \boldsymbol{\sigma}_1 \cdot \boldsymbol{\sigma}_2$ mixes $L$ and $L \pm 2$
- Three-body forces (3NF) are significant for $A \geq 3$ (Fujita-Miyazawa, chiral EFT at N$^2$LO)

**Energy scale:** $\sim 8$ MeV/nucleon binding energy; single-particle level spacings $\sim 1$--$10$ MeV; $\hbar\omega_{\mathrm{HO}} \approx 41 A^{-1/3}$ MeV.

**Qubit encoding:** Each nucleon orbital $(n_r, l, j, m_j)$ maps to one qubit per spin-species via Jordan-Wigner (protons and neutrons encoded separately). For $N_{\mathrm{shells}}$ HO major shells, the number of single-particle states per species is $\sum_{N=0}^{N_{\mathrm{shells}}-1} (N+1)(N+2)/2$. At $N_{\mathrm{shells}} = 2$: 4 proton orbitals + 4 neutron orbitals = 8 nuclear qubits (excluding spin-orbit splitting, which increases this). With $j$-coupling: $N_{\mathrm{shells}} = 2$ gives $1s_{1/2}$ (2 states) + $0p_{3/2}$ (4) + $0p_{1/2}$ (2) = 8 states per species, so $Q_{\mathrm{nuc}} = 16$.

### Block B: Electronic Core

**Fermion content:** Inner-shell electrons (e.g., $1s^2$ for Li--Ne, $1s^2 2s^2 2p^6$ for Na--Ar).

**Single-particle basis:** Hydrogenic orbitals $(n, l, m, \sigma)$ in the field of nuclear charge $Z$, exactly as implemented in GeoVac's existing Level 3 / Level 5 solvers.

**Modification from nuclear embedding:** The point-charge Coulomb potential $V = -Z/r$ is replaced by the potential of the nuclear charge distribution $\rho_p(\mathbf{r})$:

$$V_{ne}(r) = -e^2 \int \frac{\rho_p(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|} \, d^3r'$$

For $r$ outside the nuclear charge radius ($r \gg R_{\mathrm{nuc}}$), this is identical to $-Z/r$. The correction is significant only for $s$-wave electrons that penetrate the nucleus, and is of order $(Z\alpha)^4 m_e c^2 (R_{\mathrm{nuc}}/a_0)^2$.

**Energy scale:** $\sim Z^2 \times 13.6$ eV for innermost shells. For He: 24.6 eV (1s); for Li: 64.4 eV (1s).

**GeoVac status:** Fully implemented. Core screening via `geovac/core_screening.py`, frozen-core treatment via `geovac/neon_core.py`, composed pipeline via `geovac/composed_qubit.py`.

### Block C: Electronic Valence

**Fermion content:** Outer electrons in screened orbitals.

**Single-particle basis:** Hydrogenic orbitals with effective nuclear charge $Z_{\mathrm{eff}}$, plus partner-side hydrogen orbitals for bond blocks. Basis defined by `OrbitalBlock` in `MolecularSpec`.

**Coupling to core:** Via PK pseudopotential or balanced coupling (cross-center $V_{ne}$ + cross-block ERIs). Already implemented.

**Energy scale:** $\sim 1$--$10$ eV. This is the standard domain of computational chemistry.

**GeoVac status:** Fully implemented for 30 molecules across three periodic rows.

---

## 3. Inter-Block Coupling Operators

### 3.1 Nuclear $\leftrightarrow$ Electronic Core Coupling

**Primary: nuclear charge form factor.** The proton charge distribution $\rho_p(\mathbf{r})$ within the nucleus creates a potential that deviates from $-Z/r$ at distances $r \lesssim R_{\mathrm{nuc}}$. For a uniform sphere of radius $R$:

$$V(r) = \begin{cases} -\dfrac{Ze^2}{r} & r > R \\[6pt] -\dfrac{Ze^2}{2R}\left(3 - \dfrac{r^2}{R^2}\right) & r \leq R \end{cases}$$

The finite-size correction is:

$$\Delta V(r) = V_{\mathrm{extended}}(r) - V_{\mathrm{point}}(r) = \begin{cases} 0 & r > R \\[4pt] \dfrac{Ze^2}{r} - \dfrac{Ze^2}{2R}\left(3 - \dfrac{r^2}{R^2}\right) & r \leq R \end{cases}$$

This operator is diagonal in the nuclear shell model basis (it depends on $\rho_p$, which is determined by the proton wavefunction) and acts as a one-body perturbation on electronic states. Its matrix elements require the overlap of the electronic wavefunction with the nuclear interior:

$$\langle n l m | \Delta V | n' l' m' \rangle = \delta_{ll'} \delta_{mm'} \int_0^{R_{\mathrm{nuc}}} R_{nl}(r) \, \Delta V(r) \, R_{n'l'}(r) \, r^2 \, dr$$

Only $s$-wave ($l = 0$) electronic orbitals have nonzero density at the origin, so this coupling is sparse: it affects only the $s$-electron block of the electronic Hamiltonian.

**Secondary: hyperfine coupling.** The nuclear magnetic dipole moment $\boldsymbol{\mu}_I = g_I \mu_N \mathbf{I}$ creates a magnetic field at the electron:

$$H_{\mathrm{hf}} = A_{\mathrm{hf}} \, \mathbf{I} \cdot \mathbf{J}$$

where $A_{\mathrm{hf}} \propto g_I \, |\psi_e(0)|^2$ for $s$-electrons (Fermi contact) and $A_{\mathrm{hf}} \propto g_I \, \langle r^{-3} \rangle$ for $l > 0$ (dipolar). The hyperfine constant for hydrogen 1s is $A_{\mathrm{hf}} = 1420.405$ MHz $\approx 5.88 \times 10^{-6}$ eV.

In a composed framework, this is a two-block operator: $\mathbf{I}$ acts on the nuclear block, $\mathbf{J}$ acts on the electronic block. This is structurally analogous to the existing cross-block ERI coupling, but with angular momentum operators replacing Coulomb integrals.

**Nuclear quadrupole coupling.** For nuclei with $I \geq 1$, the electric quadrupole moment $Q$ couples to the electronic field gradient:

$$H_Q = \frac{eQ}{4I(2I-1)} \sum_{q} (-1)^q \, V^{(2)}_{-q} \, T^{(2)}_q(I)$$

where $V^{(2)}_q$ is the rank-2 electric field gradient tensor at the nucleus and $T^{(2)}_q(I)$ is the nuclear quadrupole tensor. This is a rank-2 tensor coupling between blocks.

### 3.2 Nuclear $\leftrightarrow$ Electronic Valence Coupling

Effectively zero for valence electrons. The nuclear charge radius is $\sim 10^{-5} a_0$, and valence electron density at the origin is suppressed by centrifugal barriers ($l \geq 1$) or by core screening ($s$-valence wavefunctions are small at $r \sim R_{\mathrm{nuc}}$ because the core electrons expel them).

The coupling is *indirect*: nuclear structure $\to$ core electron wavefunction $\to$ core screening $\to$ valence energy levels. In the composed architecture, this flows automatically: Block A determines $\rho_p$, which modifies Block B's $V_{ne}$, which modifies Block B's screening of Block C.

### 3.3 Inter-Nucleon Coupling

The NN interaction is the dominant term in Block A. Unlike the electronic Coulomb interaction, it is:

1. **Short-range** ($\sim 1$--$2$ fm), so it is nonzero only between nucleons in overlapping spatial orbitals
2. **Spin-dependent**, with different strengths in spin-singlet and spin-triplet channels
3. **Contains a tensor force** $V_T \, S_{12}$ that mixes orbital angular momentum $L$ with $L \pm 2$
4. **Has significant three-body components** for $A \geq 3$ (Fujita-Miyazawa force from virtual $\Delta$ excitation)

Standard parameterizations include:
- **Skyrme interaction**: zero-range effective force with density dependence (nuclear DFT)
- **Gogny interaction**: finite-range Gaussian (shell model CI)
- **Chiral EFT**: systematic expansion organized by power counting (Weinberg, Epelbaum et al.)
- **JISP/Daejeon16**: NN potentials fitted to phase shifts, suitable for ab initio shell model

For the composed architecture, the NN interaction is treated as a two-body operator *within* Block A, analogous to the electron-electron Coulomb interaction within electronic blocks. The key structural question is whether the NN tensor force preserves the angular momentum selection rules that give GeoVac its sparsity advantage (see Section 8, Open Questions).

---

## 4. Proof of Concept: Deuterium ($^2$H + 1 electron)

### Nuclear Block

The deuteron consists of 1 proton + 1 neutron. Since protons and neutrons are distinguishable fermions, there is no antisymmetrization between them. The ground state is $^3S_1$ (total spin $S = 1$, orbital $L = 0$), with binding energy $B = 2.2246$ MeV. The deuteron also has a $\sim 5.7\%$ D-state ($L = 2$) admixture from the tensor force.

**Nuclear basis.** With $N_{\mathrm{shells}} = 2$ HO shells per species:
- Proton: $0s_{1/2}$ (2 states) + $0p_{3/2}$ (4) + $0p_{1/2}$ (2) = 8 states
- Neutron: same = 8 states
- $Q_{\mathrm{nuc}} = 16$ qubits

At $N_{\mathrm{shells}} = 3$, adding the $0d_{5/2}$, $1s_{1/2}$, $0d_{3/2}$ orbitals: 14 states/species, $Q_{\mathrm{nuc}} = 28$ qubits. The $N_{\mathrm{shells}} = 3$ basis is needed to capture the D-state admixture ($l = 2$ in relative coordinates).

**NN interaction.** The deuteron is the simplest nuclear bound state --- only one NN pair. The Hamiltonian in the nuclear block is:

$$H_{\mathrm{nuc}} = T_p + T_n + V_{NN}(|\mathbf{r}_p - \mathbf{r}_n|, \text{spin, isospin})$$

In the HO basis, this is a standard two-body matrix element calculation, well-studied in nuclear structure (see Suhonen, "From Nucleons to Nucleus," 2007).

### Electronic Block

One electron in the field of a (nearly) point charge $Z = 1$. This is hydrogen, GeoVac's Level 1 system. At $n_{\mathrm{max}} = 2$: $Q_{\mathrm{elec}} = 8$ qubits.

### Coupling

The finite nuclear size correction for deuterium is:

$$\Delta E_{1s} \approx \frac{2}{3} (Z\alpha)^4 m_e c^2 \left(\frac{R_d}{a_0}\right)^2 \approx 3.4 \times 10^{-10} \text{ eV}$$

using the deuteron charge radius $R_d = 2.1421$ fm. This is $\sim 10^{-10}$ of the electronic binding energy.

### Resource Estimate

| Component | Qubits | Dominant terms |
|:----------|:-------|:---------------|
| Nuclear ($N_{\mathrm{shells}} = 2$) | 16 | $V_{NN}$ two-body |
| Electronic ($n_{\mathrm{max}} = 2$) | 8 | $T_e + V_{ne}$ one-body |
| Coupling | 0 (modifies $V_{ne}$) | $\Delta V$ finite-size |
| **Total** | **24** | |

At $N_{\mathrm{shells}} = 3$: $Q_{\mathrm{nuc}} = 28$, total 36 qubits.

The coupling $\Delta V$ modifies only the electronic $s$-orbital one-body matrix elements. It adds no new Pauli terms --- it adjusts existing $V_{ne}$ coefficients by $\sim 10^{-10}$ relative amounts.

---

## 5. Proof of Concept: $^4$He + 2 Electrons

### Nuclear Block

Helium-4 has 2 protons + 2 neutrons. It is doubly magic ($Z = 2$, $N = 2$) with a $0^+$ ground state (spin-0, $L = 0$) and binding energy $B = 28.296$ MeV. The charge radius is $R = 1.6755$ fm.

**Nuclear basis.** Each species (proton, neutron) is a 2-fermion system, requiring antisymmetrization within species.
- $N_{\mathrm{shells}} = 2$: 8 states/species, $Q_{\mathrm{nuc}} = 32$ qubits
- $N_{\mathrm{shells}} = 3$: 14 states/species, $Q_{\mathrm{nuc}} = 56$ qubits

**NN interaction.** Six NN pairs: $\binom{2}{2}_p + 2 \times 2 + \binom{2}{2}_n = 1 + 4 + 1 = 6$ pairs (1 pp, 4 pn, 1 nn). The proton-proton pair also has a Coulomb repulsion contribution ($\sim 0.7$ MeV), which is a small perturbation on the $\sim 28$ MeV total binding.

### Electronic Block

Two electrons: GeoVac's Level 3 system (He atom), solved to 0.004% accuracy with the 2D variational solver + cusp correction. At $n_{\mathrm{max}} = 2$: $Q_{\mathrm{elec}} = 8$ qubits.

### Coupling

**Finite nuclear size.** The correction for $^4$He is:

$$\Delta E_{1s} \approx \frac{2}{3} Z^4 \alpha^4 m_e c^2 \left(\frac{R_{\mathrm{He}}}{a_0}\right)^2 \approx 2 \times 10^{-9} \text{ eV}$$

The isotope shift between $^3$He and $^4$He (different charge radii) produces a measurable $\sim 30$ GHz shift in the $2^3S_1 \to 2^3P$ transition, but this is dominated by the mass difference, not the charge radius difference.

**Nuclear quadrupole.** Zero: $^4$He has $I = 0$.

**Mass polarization.** The recoil correction $H_{\mathrm{MP}} = (1/M_{\mathrm{nuc}}) \, \mathbf{p}_1 \cdot \mathbf{p}_2$ between the two electrons contributes $\sim 1.3 \times 10^{-4}$ eV. This is larger than the finite nuclear size correction and can be expressed as a two-body operator in the electronic block (no nuclear qubits needed).

### Resource Estimate

| Component | Qubits ($N_{\mathrm{shells}} = 2$) | Qubits ($N_{\mathrm{shells}} = 3$) |
|:----------|:----------------------------------:|:----------------------------------:|
| Nuclear (p + n) | 32 | 56 |
| Electronic ($n_{\mathrm{max}} = 2$) | 8 | 8 |
| **Total** | **40** | **64** |

Both estimates are within reach of near-term quantum hardware. The nuclear block dominates the qubit count.

---

## 6. Formal Structure: Composed Nuclear-Electronic Hamiltonian

The complete Hamiltonian is:

$$H = H_{\mathrm{nuc}} + H_{\mathrm{elec}} + H_{\mathrm{coupling}}$$

### Nuclear Hamiltonian

$$H_{\mathrm{nuc}} = \sum_{i \in \text{protons}} T_i^{(\mathrm{nuc})} + \sum_{i \in \text{neutrons}} T_i^{(\mathrm{nuc})} + \sum_{i \in p} V_{\mathrm{WS}}^{(p)}(r_i) + \sum_{i \in n} V_{\mathrm{WS}}^{(n)}(r_i) + \sum_{i} H_{\mathrm{SO}}(r_i) + \sum_{i < j} V_{NN}(r_{ij})$$

where the sums run over proton and neutron indices separately, $V_{\mathrm{WS}}^{(p/n)}$ are the Woods-Saxon potentials (differing by Coulomb for protons), and $V_{NN}$ is the two-body nuclear force.

### Electronic Hamiltonian

$$H_{\mathrm{elec}} = \sum_{i \in \text{electrons}} T_i^{(\mathrm{elec})} + \sum_{i} V_{ne}(\mathbf{r}_i) + \sum_{i < j} \frac{e^2}{|\mathbf{r}_i - \mathbf{r}_j|}$$

This is the standard electronic Hamiltonian, already fully implemented in GeoVac.

### Coupling Hamiltonian

$$H_{\mathrm{coupling}} = V_{\mathrm{finite\text{-}size}} + H_{\mathrm{hyperfine}} + H_{\mathrm{recoil}}$$

**Finite nuclear size:**

$$V_{\mathrm{finite\text{-}size}} = \sum_{i \in \text{electrons}} \Bigl[ V_{ne}^{(\mathrm{extended})}(\mathbf{r}_i; \{\mathbf{r}_p\}) - V_{ne}^{(\mathrm{point})}(\mathbf{r}_i) \Bigr]$$

This depends on the proton positions $\{\mathbf{r}_p\}$ in the nuclear block, making it a genuine cross-block operator. In the composed architecture, it would be encoded as a set of Pauli strings that span both the nuclear and electronic qubit registers.

**Hyperfine coupling:**

$$H_{\mathrm{hyperfine}} = A_{\mathrm{hf}} \, \mathbf{I} \cdot \mathbf{J}$$

This is a product operator: $\mathbf{I}$ acts on nuclear qubits, $\mathbf{J}$ on electronic qubits. In Pauli representation: $\mathbf{I} \cdot \mathbf{J} = I_x J_x + I_y J_y + I_z J_z$, where each component is a tensor product of Pauli strings on the two registers. This adds $O(Q_{\mathrm{nuc}} \times Q_{\mathrm{elec}})$ cross-register Pauli terms.

**Recoil correction:**

$$H_{\mathrm{recoil}} = \frac{1}{2M_{\mathrm{nuc}}} \left(\sum_{i \in \text{electrons}} \mathbf{p}_i\right)^2$$

This is purely electronic (center-of-mass recoil) and requires no nuclear qubits. It is a two-body operator that can be absorbed into the electronic block as a correction to the kinetic energy.

### Implementation Status

| Term | GeoVac status | Notes |
|:-----|:-------------|:------|
| $T_i^{(\mathrm{elec})}$ | Implemented | Graph Laplacian (Level 1) or spectral (Levels 2--5) |
| $V_{ne}$ (point charge) | Implemented | One-body, diagonal in $(n,l,m)$ |
| $e^2/r_{ij}$ (e-e) | Implemented | Slater integrals, Gaunt selection rules |
| $T_i^{(\mathrm{nuc})}$ | New (uses HO infrastructure) | `harmonic_shell.py` from Track NA |
| $V_{\mathrm{WS}}$ | New | Angular integrals reuse Gaunt; radial integrals are new (Woods-Saxon radial form) |
| $H_{\mathrm{SO}}^{(\mathrm{nuc})}$ | New | $\mathbf{l} \cdot \mathbf{s}$ coupling; requires $j$-$j$ basis or $LS \to jj$ transformation |
| $V_{NN}$ | Genuinely new | Short-range, spin-dependent, tensor. No analog in electronic structure |
| $V_{\mathrm{finite\text{-}size}}$ | New (small) | Modifies existing $V_{ne}$ matrix elements for $s$-orbitals |
| $H_{\mathrm{hyperfine}}$ | New | Cross-register $\mathbf{I} \cdot \mathbf{J}$ product operator |
| $H_{\mathrm{recoil}}$ | New (electronic only) | Two-body correction to kinetic energy |

### Perturbative Hierarchy

The coupling terms satisfy $\|H_{\mathrm{coupling}}\| / \|H_{\mathrm{elec}}\| \sim 10^{-7}$--$10^{-4}$, meaning they can be treated as perturbations on top of separate nuclear and electronic solutions. This has two consequences:

1. **Validation strategy:** Solve $H_{\mathrm{nuc}}$ and $H_{\mathrm{elec}}$ independently first, verify against known results, then add $H_{\mathrm{coupling}}$ and check that perturbation theory agrees with the full diagonalization.
2. **Precision requirement:** Resolving $H_{\mathrm{coupling}}$ effects on a quantum computer requires energy resolution below $\sim 10^{-7}$ Ha ($\sim \mu$eV). This is below current VQE noise floors but accessible to QPE with sufficient circuit depth.

---

## 7. Market Positioning

### Existing Nuclear Structure Approaches

| Method | Scope | Scaling | Limitation |
|:-------|:------|:--------|:-----------|
| NCSM (Navratil, Barrett) | $A \leq 16$ | $O(N!)$ basis dimension | Exponential wall at $A \sim 12$--$16$ |
| Nuclear DFT (Skyrme, Gogny) | All nuclei | $O(N^3)$ | Mean-field; misses correlations |
| Coupled Cluster (Hagen et al.) | $A \leq 132$ ($^{132}$Sn) | Polynomial | Needs soft interaction; open-shell challenging |
| Lattice QCD | Fundamental | Extreme | Currently limited to $A \leq 4$ with physical pion mass |
| Quantum simulation (Roggero et al.) | Active research | $\sim Q^{2-4}$ gate count | Early stage; nuclear NN encoding is non-trivial |

### GeoVac's Potential Contribution

The GeoVac framework offers several structural advantages that *may* transfer to the nuclear domain:

1. **Angular momentum sparsity.** If the NN interaction preserves the angular momentum selection rules analogous to Gaunt integrals, the nuclear Hamiltonian would inherit GeoVac's $O(Q^{2-3})$ Pauli scaling rather than the generic $O(Q^4)$. This is the central open question (Section 8).

2. **Composed block architecture.** Separating protons and neutrons into distinct fermion blocks, with inter-species coupling handled as cross-block terms, is natural in the composed framework. The existing `MolecularSpec` / `OrbitalBlock` infrastructure already supports multi-block Hamiltonians.

3. **Unified nuclear-electronic description.** No existing code, to our knowledge, places nuclear and electronic degrees of freedom in a single qubit Hamiltonian for quantum simulation. The composed architecture's ability to handle $10^6$-fold energy scale separations (nuclear MeV vs electronic eV) via block decomposition is a potential differentiator.

4. **Ecosystem export.** GeoVac's existing OpenFermion/Qiskit/PennyLane pipeline (`ecosystem_export.py`) would immediately provide nuclear-electronic Hamiltonians in standard quantum computing formats.

**Honest caveats.** Nuclear structure theory is a mature field with deep domain expertise. GeoVac's graph Laplacian / S$^3$ formalism was derived for the Coulomb problem; its applicability to the short-range, spin-dependent NN interaction is unproven. The angular momentum basis is compatible (both use spherical harmonics and Clebsch-Gordan coefficients), but the radial structure is fundamentally different (HO vs hydrogenic, finite-range vs $1/r$). The specification below is a mathematical framework; whether it produces competitive results is an empirical question.

---

## 8. Open Questions for PI Review

The following questions require physics judgment and cannot be resolved by computation alone.

1. **Angular momentum selection rules for $V_{NN}$.** The NN central force preserves $l$ and thus respects Gaunt-like selection rules. The tensor force $S_{12}$, however, couples $L$ to $L \pm 2$. Does this break the angular sparsity, or does it merely expand the selection rule from $|\Delta l| \leq 1$ (Gaunt) to $|\Delta l| \leq 2$ (tensor)? If the latter, the sparsity structure is modified but not destroyed. Nuclear structure references (e.g., Suhonen Ch. 11) would clarify this.

2. **Two fermion species in composed blocks.** The existing `OrbitalBlock` assumes all particles are identical electrons. Protons and neutrons are distinct species: they share the same spatial orbitals but do not antisymmetrize with each other. The qubit encoding needs to handle $\mathcal{H}_p \otimes \mathcal{H}_n$ (tensor product, not wedge product) between species. Is this compatible with the existing Jordan-Wigner pipeline, or does it require a new encoding scheme?

3. **$j$-$j$ coupling vs $LS$ coupling.** Nuclear physics universally uses $j$-$j$ coupling ($j = l \pm 1/2$ labels single-particle states), while GeoVac's electronic structure uses $LS$ coupling (separate $l, m_l, \sigma$ labels). The $j$-$j$ basis diagonalizes the spin-orbit interaction, which is essential for reproducing magic numbers. Should the nuclear block use $j$-$j$ from the start, or should $LS$ be used for consistency with the electronic blocks, with spin-orbit treated as a perturbation?

4. **Energy scale separation and precision.** The nuclear binding energy ($\sim 10^6$ eV) exceeds the electronic binding energy ($\sim 10$ eV) by a factor of $\sim 10^5$. The coupling terms are another $\sim 10^{-7}$ smaller. On a quantum computer, the energy is encoded in the phase of the wavefunction. Can a single QPE circuit resolve energy differences spanning 12 orders of magnitude, or must the blocks be solved separately with their own phase estimation rounds?

5. **Nuclear force parameterization.** Which NN interaction is most appropriate for a composed-architecture proof of concept? Options range from simple (Minnesota potential, 3 Gaussians) to sophisticated (chiral N$^3$LO, hundreds of terms). The choice determines the two-body matrix element complexity and thus the Pauli term count in Block A.

6. **Three-nucleon forces.** For $A \geq 3$, three-body nuclear forces contribute $\sim 1$--$2$ MeV/nucleon. In second quantization, these become 6-index tensors $V_{ijk,lmn}$. GeoVac currently handles only one-body and two-body operators. Is a three-body extension feasible within the Pauli framework, and what is the scaling?

7. **Measurability of finite nuclear size on quantum hardware.** The finite-size correction for He is $\sim 10^{-9}$ eV $\sim 10^{-10}$ Ha. Current quantum hardware has energy resolution $\sim 10^{-3}$--$10^{-4}$ Ha (VQE) to $\sim 10^{-6}$ Ha (QPE with error correction). Is the nuclear-electronic coupling observable on any foreseeable quantum hardware, or is this purely a theoretical exercise?

8. **Proton charge radius puzzle.** The $\sim 4\%$ discrepancy between the proton charge radius measured from muonic hydrogen vs electronic hydrogen (partially resolved, 2019 CODATA) involves precisely the finite nuclear size correction. Could a composed nuclear-electronic Hamiltonian on a quantum computer provide any insight, or is the precision requirement ($\sim 10^{-14}$ eV) hopelessly beyond reach?

---

## References

- M. G. Mayer and J. H. D. Jensen, *Elementary Theory of Nuclear Shell Structure* (Wiley, 1955).
- I. Talmi, *Simple Models of Complex Nuclei* (Harwood, 1993).
- P. Ring and P. Schuck, *The Nuclear Many-Body Problem* (Springer, 1980).
- J. Suhonen, *From Nucleons to Nucleus* (Springer, 2007).
- A. Roggero, A. C. Y. Li, J. Carlson, R. Gupta, and G. N. Perdue, "Quantum computing for neutrino-nucleus scattering," Phys. Rev. D **101**, 074038 (2020).
- P. Navratil, J. P. Vary, and B. R. Barrett, "Properties of $^{12}$C in the ab initio nuclear shell model," Phys. Rev. Lett. **84**, 5728 (2000).
- G. Hagen, T. Papenbrock, M. Hjorth-Jensen, and D. J. Dean, "Coupled-cluster computations of atomic nuclei," Rep. Prog. Phys. **77**, 096302 (2014).
- E. Epelbaum, H.-W. Hammer, and U.-G. Meissner, "Modern theory of nuclear forces," Rev. Mod. Phys. **81**, 1773 (2009).
- GeoVac Papers 0, 7, 14, 17, 19 (composed architecture, Pauli scaling, balanced coupling).
