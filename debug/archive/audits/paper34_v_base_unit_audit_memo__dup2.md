# Paper 34 §V audit against base-unit decomposition (§IV.5, Table 2)

**Date:** 2026-05-09
**Source:** post-§IV.5 addition (this session); audit pass over §V matches catalogue and §V.B off-precision rows.

## Method

For each catalogue row I checked:
1. Are all base units required by the observable's compound dimension actually injected by some projection in the listed chain?
2. Is every projection name in the chain present in §III, or is it an informal sub-mechanism?
3. Are external calibration constants (g-factors, structural radii, multi-loop literature inputs) tagged explicitly as Layer-2 inputs, or treated as if they were native?

## Findings

### A. Clean rows — chain consistent with base-unit dictionary

Layer 1 (zero-projection) rows are all clean. Pauli count, magic numbers, Gaunt selection rules, π-free certificates, angular sparsity, discrete c_1 — all dimensionless integer/rational, no projections, consistent.

Single-projection rows that pass:
- Π = 1/(48π²) VP coefficient — spectral action, dimensionless (Λ-independent coefficient)
- β(α) = 2α²/(3π) — spectral action, dimensionless
- D² heat kernel SD coefficients — spectral action; the energy variable is Λ²
- F = ζ_R(2) = D_{n²}(d_max) — Fock + Dirichlet (graph-internal structural identification)
- Δ⁻¹ = g_3^Dirac = 40 — spinor lift, integer degeneracy
- B = 42 — Layer 1 graph trace
- E_Cas^{S³,Dirac} = 17/480 — Fock + spinor lift, dimensionless rational

Two-projection rows that pass:
- Uehling shift -27.13 MHz — spectral action ∘ Fock, energy
- 21 cm HF-2 (with Parker-Toms c₁) — chain explicit, dimensions consistent
- All Bethe log rows (LS-3 acceleration form, LS-4 Drake-Swainson) — dimensionless, chain explicit

Three-projection rows:
- Lamb shift LS-4 — Fock ∘ spectral action ∘ Sturmian ∘ Drake-Swainson, energy
- μH Lamb shift — same as above + rest-mass, energy

**Net clean count: ~22 of 35 catalogued rows.**

### B. Rows that omit Fock as foundational anchor

Some rows quote results in Hartree (energy in atomic units) without listing Fock. Fock is the foundational projection that anchors atomic units, so any energy-dimensioned result implicitly invokes it. Rows where Fock should arguably be explicit:

- 2p doublet α²/32 Ha (Z=1) — chain says "spinor lift" alone; should be "spinor lift ∘ Fock"
- Dirac fine-structure formula (n≤4) — same: "spinor lift" should be "spinor lift ∘ Fock"

Cosmetic, not a real inconsistency, but flagging for completeness.

### C. Rows where compositional rest-mass projections need explicit accounting

Observation 3 (Compositional projection): every distinct particle species requires its own rest-mass projection, and Breit retardation gates on two of them. Several rows have the right physics but undercount the rest-mass projections in the chain.

**21 cm hyperfine HF-1** — variables listed: α, g_p, m_e/m_p. The mass ratio implies rest-mass projection for the proton (m_p replacing the foundational m_e anchor). Chain currently lists "spinor ∘ Fermi contact ∘ Fock" — no rest-mass. The reduced-mass content is implicit in "Bohr-Fermi" but not chain-tagged.

**μH 1S Bohr-Fermi 182.4433 meV** — variables: α, g_p, m_μ/m_e, m_e/m_p. **Two** mass ratios (lepton-electron and electron-proton) imply two distinct particle species beyond electron. Chain lists "rest-mass projection ($m_e \to m_\mu$)" — only one. Should arguably list two: muon rest-mass (lepton swap) and proton rest-mass (nucleus mass). At equal-mass limit ($m_l/m_n = 1$) Breit retardation gates on, so the two-rest-mass structure here is the same machinery that surfaces in Ps 1S-2S.

**Mu 1S HFS** — same pattern: two rest-mass projections (m_μ for nucleus, m_e for lepton kept). Chain lists one.

**D 1S hyperfine** — α, g_d, m_e/m_d. One rest-mass for deuteron. But g_d is a Layer-2 calibration input.

**Mu 1S-2S** — m_p → m_μ swap. One rest-mass projection. Consistent with single particle change.

This isn't wrong physics — every row's numerical result is computed correctly. It's a chain-level taxonomy gap: Observation 3's structural counting (one rest-mass per species) isn't reflected in the chain notation.

### D. Informal projection names used in chains but not in §III

Audit yields a list of mechanism names that appear in chains but aren't in §III's sixteen-projection list. Each is either (a) a sub-feature of an existing §III projection, (b) a Layer-2 input, or (c) a candidate new §III projection. My classification:

| Informal name | Where it appears | Classification | Recommended action |
|---|---|---|---|
| "Fermi contact" | 21cm, μH, D, Mu HFS rows | Sub-feature: spinor + 3j + Fock at delta(r₁₂); the contact density is a Layer-1 graph quantity, the I·S coupling is spinor-derived | Replace with "spinor ∘ Fock" + note about delta-function vertex; OR add §III sub-projection "Fermi contact vertex" |
| "Conformal coupling" | E_Cas = 1/240 | Specific value m² = 1 in rest-mass projection (per Paper 35 KG-3) | Tag as "Fock ∘ rest-mass (m² = 1, conformal value)" |
| "Casimir/ζ-reg" | HO Casimir, Stefan-Boltzmann | Sub-feature of spectral action: heat kernel in Λ → ∞ limit, equivalently ζ-regularization of mode sum | Absorb into "spectral action" with parenthetical (ζ-reg) |
| "Discrete Green's" | Coulomb-like Green's function | Layer-1 operation: matrix inverse (L+I)⁻¹ on the bare graph | Reclassify as Layer-1 derived quantity, no projection |
| "Vertex parity" | D_even closed form | Sub-feature: vertex parity selection rule from spinor lift's CG decomposition + 3j coupling | Replace with "spinor lift" with vertex selection note |
| "Magnetization-density" | μH Zemach, HF-4 | **Candidate 17th §III projection** — operator-level construction in `geovac/magnetization_density.py`; treats Zemach radius r_Z as focal length | Promote to §III as new projection (see below) |
| "Wigner 6j" | He fine structure | Derived structure: 6j arises from coupling three angular momenta, computable from 3j | Note as "Wigner 3j (with 6j coupling)" |
| "Proper-time integration" | Heisenberg-Euler | Sub-feature of spectral action: Schwinger Mellin parametrization | Absorb into "spectral action" |
| "Internal multi-focal" | He 2³P (1s, 2p at distinct Z_eff) | **Candidate composition rule** — two Fock projections at different effective charges in one chain | See §F below |

### E. External calibration constants entering raw

g_p, g_d (proton/deuteron g-factors), r_Z (Zemach radius), Penin-Pivovarov α⁴ Breit input, Eides Tab. 7.3 multi-loop QED, Karshenboim corrections — these are Layer-2 inputs that close residuals but the chain doesn't always tag them.

The cleanest tagging convention in the catalogue is HF-4: "+ external Zemach r_Z" inline in the row title. This should be the standard for any Layer-2-input-augmented match.

### F. Two structural mechanisms that may deserve §III status

#### F.1. Magnetization-density projection (candidate §III.17)

Used in HF-4 Zemach and μH Zemach. The mechanism: replace point-magnetization with extended distribution at characteristic radius r_Z; modify the Fermi-contact term by a structural factor proportional to r_Z. Operator-level construction exists (`geovac/magnetization_density.py`).

Variables: r_Z (Zemach radius)
Dimension: length
Base-unit injection: new [L] scale (separate from internuclear or atomic; nuclear-internal)
Transcendental signature: ring-preserving over the framework's atomic-units ring at leading order; sub-leading corrections are profile-dependent (Gaussian vs exponential — verified independent at leading)
Role: scale-extension (new [L] instance)

This would make [L] a five-instance axis: foundational (a₀), stereographic coordinate (r), Wigner D position ($\vec{R}_{AB}$), mol-frame distance (R), and now magnetization radius (r_Z). All five are physically distinct length scales (atomic, coordinate-explicit, internuclear, intra-molecule, intra-nuclear).

#### F.2. Multi-focal Fock composition (candidate §III.18 OR composition-rule observation)

Used in He internal multi-focal: 1s and 2p see different effective Z (Z_eff(1s)=2, Z_eff(2p)=1) in the same atom. The same physical electron Hilbert space carries two Fock anchors. Currently §III treats Fock as one projection per system; multi-Z_eff at the orbital level is a composition pattern.

Two ways this could land:
- **As a 17th projection**: "multi-Z_eff Fock" with variables {Z_a, Z_b}, dimension energy via two Rydbergs.
- **As a composition rule**: noting that two Fock projections can compose at the orbital level, similar to how two rest-mass projections compose at the particle level via Breit retardation.

The composition-rule reading is structurally cleaner because it parallels Observation 3 (Breit composition of rest-mass). Specifically: just as a *ratio* of two rest-mass projections gates Breit retardation, a *ratio* of two Fock-Z values gates the internal multi-focal pattern. The He fine-structure rows are the empirical anchor for this composition.

If we go with the composition-rule reading, Observation 3 generalizes to:
> Variables that are functions of prior projection variables introduce new structural content. Currently identified composition rules: (a) ratio of two rest-mass [M] scales → Breit retardation at α⁴; (b) ratio of two Fock-Z [Q] labels → internal multi-focal at... (whatever order this surfaces).

The "internal multi-focal at α^p" question would then become a falsifiable structural prediction.

## Recommendations (decision matrix)

| Action | Cost | Benefit | Recommendation |
|---|---|---|---|
| Standardize Fock-explicit in atomic-units chains | low (search/replace ~5 rows) | cosmetic consistency | apply |
| Multi-particle rows list rest-mass per species | medium (need to audit ~10 rows) | structural taxonomy alignment | apply |
| Replace informal names per Table D | medium (~10 rows touched) | matches §III canonical list | apply |
| Promote magnetization-density to §III.17 | medium (write §III subsection) | extends [L] axis to five instances; closes Zemach Layer-2 | flag for PI decision |
| Promote multi-focal to §III.17/18 OR Observation 3 generalization | medium (write either §III subsection or extend Observation 3) | makes He fine structure structurally first-class; adds composition-rule axis | flag for PI decision; composition-rule reading is cleaner |
| Tag Layer-2 calibration constants explicitly in chain | low (notational convention + ~6 rows) | audit-ready provenance | apply |

## Hunting candidates (for the next phase, after audit applied)

The user's stated goal: "go hunting for new unit signatures and attempt some decompositions with them."

The audit surfaces three loose threads:

1. **Magnetization-density [L] scale (r_Z)** — if promoted, it asks: are there other intra-nuclear or intra-particle length scales that should enter as projections? Charge radius r_p ≠ Zemach r_Z; deuteron quadrupole moment Q_d adds a tensorial intra-nuclear length; nucleon-nucleon scattering length a_NN. Each is a candidate independent [L] scale.

2. **Multi-focal Fock composition** — He's two Z_eff values are the cleanest case. What about systems with three or more orbital-resolved Z_eff (LiH composed has core+valence with distinct Z; BeH₂, H₂O have more)? Is the composition rule a binary operation ("ratio") or is there a multi-focal generalization with N free parameters?

3. **Inverse temperature ↔ inverse mass duality** — the observation/temporal-window projection injects [T] (which equals 1/[E] at ℏ=1, equivalently inverse energy). The rest-mass projection injects [M] (which also equals [E] at c=1). Both project the framework into a 1D parameter space adjacent to [E] but distinct from it. Are there other "inverse-energy-adjacent" projections? Cosmological scale factor (Hubble length, equivalently 1/Hubble time)? Magnetic field strength B at e=ℏ=1 has dim energy/charge, which in atomic units is dim energy. An "external field" projection would inject a new [E] scale not currently in the framework.

If you want to go hunting, the most concrete short-loop probes are (a) charge-radius vs Zemach-radius distinction (does the framework predict their structural relationship?) and (b) external uniform magnetic field — does this fit as a 17th projection that injects [E] via the Zeeman gap?

