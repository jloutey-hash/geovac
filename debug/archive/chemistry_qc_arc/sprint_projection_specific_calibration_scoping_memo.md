# Sprint Projection-specific calibration scoping

**Date:** 2026-05-29
**Path:** Gravity arc completion, follow-on to W3 (Wolfenstein) sprint of 2026-05-08. Diagnostic-only.
**Verdict:** **NO-GO for a global "test every projection" sprint, with substantive structural finding.** Of 28 projections in Paper 34, only ~5 carry parametric Class-D calibration data, and of those only **2 (§III.14 rest mass and §III.16 Breit retardation) show known GeoVac-internal structure** in their calibration values (Koide cone and Paper 2 K-formula respectively). The remaining 26 projections either carry no parametric calibration (geometric/structural) or carry external system-specific values with no GeoVac-internal relations. The user's reframe of W3 is sharp and useful, but it sharpens *where* to look, not whether the global question has a positive answer.

## 1. Question

PI question 2026-05-29: the W3 sprint of 2026-05-08 tested the GLOBAL hypothesis (do Wolfenstein values appear in master Mellin engine M1/M2 rings?) and got a clean negative. Sharper reframe: does each of Paper 34's 28 projections carry its OWN calibration parameter set, with internal structure even if the global Mellin rings don't host it?

This is the projection-specific analog of the W3 question. Diagnostic-only.

## 2. Methodology

For each of the 28 projections (§III.1 – §III.28), identify:

(a) What calibration data — if any — is attached to this projection?
(b) Is the data parametric (Class-D Class-1 numerical input) or structural/conventional?
(c) Does the data set have known internal structure (relations, ratios, predictive content)?

Four classes of calibration content per projection:

- **Class G**: geometric/structural (e.g., 1/(4π) from Hopf measure) — derived, no freedom.
- **Class C**: convention-valued (e.g., gauge choice, KO-dim convention, Euler angles) — arbitrary representational.
- **Class E**: system-specific external input (e.g., bond length R, nuclear ⟨r²⟩_E) — measurable from experiment.
- **Class D**: universal calibration parameters (e.g., cutoff f, gauge couplings g, masses m, α) — Class-1 calibration in `external_input_three_class_partition`.

Only Class D is structurally interesting for the W3 question.

## 3. Projection-by-projection inventory

| # | Projection | Calibration data | Class | Known internal structure? |
|---|---|---|---|---|
| 1 | Fock conformal (§III.1) | nuclear Z (integer) | E | quantum number, no |
| 2 | Hopf bundle (§III.2) | 1/(4π) = Vol(S²)/(4π) | G | structural, not a parameter |
| 3 | Bargmann–Segal (§III.3) | ℏω per system | E | none beyond system |
| 4 | Stereographic (§III.4) | p₀ = √(−2E) | derived | structural, not free |
| 5 | Sturmian (§III.5) | exponent λ = Z/n | derived | structural, not free |
| 6 | **Spectral action (§III.6)** | **cutoff function f** | **D** | Mellin moment map (G4-5); no constraint on $\phi(s)$ values |
| 7 | Spinor lift (§III.7) | γ-matrix basis, KO-dim | C | convention |
| 8 | Wigner 3j (§III.8) | none | — | structural |
| 9 | Wigner D (§III.9) | Euler angles | C | convention |
| 10 | **Wilson loop (§III.10)** | **g₁, g₂, g₃** | **D** | RGE near-unification at GUT scale (standard SM) |
| 11 | Vector photon promotion (§III.11) | 1/(4π) per loop | G | derived geometric |
| 12 | Mol-frame (§III.12) | bond length R | E | experimental input |
| 13 | Drake–Swainson (§III.13) | regularization scale K, D denominator | derived | structural choice |
| 14 | **Rest mass (§III.14)** | **m_e, m_μ, m_τ, m_u, m_d, m_s, m_c, m_b, m_t, m_p, m_n** | **D** | **Koide cone at 1 arcsec (charged leptons)** |
| 15 | Observation / temporal window (§III.15) | τ (window) | C | convention |
| 16 | **Breit retardation (§III.16)** | **α = 1/137.0359990…** | **D** | **K = π(B + F − Δ) at 10⁻⁸ (Paper 2, conjectural)** |
| 17 | Nuclear charge density (§III.17) | ⟨r²⟩_E per nucleus | E | nuclear physics (NN scattering) |
| 18 | Nuclear magnetization density (§III.18) | Zemach r_Z per nucleus | E | nuclear physics; §V.D convention catalogue documents inter-compilation differences |
| 19 | Nuclear tensor multipole (§III.19) | Q_N per nucleus | E | nuclear physics |
| 20 | Phillips–Kleinman (§III.20) | core eigenvalues E_c | derived | structural |
| 21 | Multipole / Gaunt termination (§III.21) | none | — | structural |
| 22 | Bipolar harmonic / Drake combining (§III.22) | (3/50, −2/5, 3/2, −1) | G | pure rationals; structurally fixed |
| 23 | Symmetry / Young tableau (§III.23) | irrep | C | group theory |
| 24 | Adiabatic / Born–Oppenheimer (§III.24) | none | — | structural separation |
| 25 | Coupled-channel (§III.25) | none | — | derived |
| 26 | Gauge choice (§III.26) | gauge parameter ξ | C | convention |
| 27 | Wick rotation (§III.27) | β (KMS temperature, often 2π for BW) | C | convention or experimental |
| 28 | Apparatus identity (§III.28) | Born rule | D | external (Class 1 per memory `external_input_three_class_partition`) |

**Counts:**
- Class G (geometric/structural): §III.2, §III.11, §III.22 → 3 projections
- Class C (convention-valued): §III.7, §III.9, §III.15, §III.23, §III.26, §III.27 → 6 projections
- Class E (system-specific external): §III.1, §III.3, §III.12, §III.17, §III.18, §III.19 → 6 projections
- Class D (universal calibration parameters): §III.6, §III.10, §III.14, §III.16, §III.28 → 5 projections
- Derived/structural (no free input): §III.4, §III.5, §III.8, §III.13, §III.20, §III.21, §III.24, §III.25 → 8 projections

5 projections carry Class-D calibration data. Of these, **2 show known internal structure**.

## 4. The two Class-D projections with known internal structure

### 4.1 §III.14 rest mass — Koide cone

The Koide relation (Koide 1983):
$$
\frac{m_e + m_\mu + m_\tau}{(\sqrt{m_e} + \sqrt{m_\mu} + \sqrt{m_\tau})^2} = \frac{2}{3}
$$

verified at sub-arcsecond precision (~1 arcsec on the 3-vector cone angle, ~10⁻⁵ on the ratio).

This is an internal-structure observation among the charged-lepton mass values. Calibration data within a single projection (§III.14, rest mass) showing a structural ratio.

**Honest status:**
- Empirical observation, not a derived theorem in GeoVac.
- Koide has been catalogued for decades; many proposed mechanisms in the wider literature (Yukawa textures, family-symmetry models, no consensus).
- GeoVac does NOT currently derive Koide from any internal principle.
- Memory file: `memory/koide_cone_clean_negative.md` — third independent calibration-data NEGATIVE on Koide-specific PSLQ + structural-angle test. The Koide relation itself stands as empirical; the GeoVac-internal mechanism question is closed-NEGATIVE.

Sprint W3 of 2026-05-08 included a lepton-mass spectrum test as a sibling track. 2/9000 within-1σ hits at random-chance density. No internal structure beyond Koide itself was found.

### 4.2 §III.16 Breit retardation — Paper 2 K-formula

The conjectural relation:
$$
K \equiv \pi(B + F - \Delta) = 1/\alpha + O(10^{-8})
$$

with B = 42 (finite Casimir trace at m=3), F = π²/6 = $D_{n^2}(s=4)$ (Fock-degeneracy Dirichlet at packing exponent), Δ = 1/40 = $g_3^{\rm Dirac}$ (Dirac mode count at the natural cutoff).

This is calibration data (α) within a single projection (§III.16) showing a structural three-part decomposition.

**Honest status:**
- 12 mechanisms eliminated for K-formula derivation (Phases 4B-4I, Sprint A, Sprint K-CC; CLAUDE.md §2 WH5).
- Combination rule remains conjectural (Paper 2 in Observations; combination-rule "conjectural" label protected per CLAUDE.md §13.5).
- Three component identifications (B, F, Δ) are derived; the SUM matching α⁻¹ is the unexplained coincidence.

## 5. The three Class-D projections WITHOUT known internal structure

### 5.1 §III.6 spectral action — cutoff function f

CC scoping memo (Task 1 of this thread) covers this. The sector-wise Mellin moment map (G4-5) gives a structural articulation:
- $\phi(0)$ ↔ tip / $S_{\rm BH}$
- $\phi(1)$ ↔ Einstein–Hilbert
- $\phi(2)$ ↔ cosmological constant

These are independently tunable. No GeoVac-internal constraint between them. The CC problem requires $\phi(2)/\phi(1)^2 \approx 10^{-124}$ with $\phi(0), \phi(1)$ both $O(1)$, with no natural cutoff function satisfying it.

**No internal structure beyond moment-map decomposition.**

### 5.2 §III.10 Wilson loop — gauge couplings g₁, g₂, g₃

Standard SM observations:
- RGE running brings the couplings near unification at $M_{\rm GUT} \sim 10^{16}$ GeV (especially with MSSM matter content).
- Exact unification at GUT requires specific matter content (SU(5) needs 24, SO(10) needs 45/16 etc.).
- Weinberg angle $\sin^2\theta_W$ relates g₁, g₂ at electroweak scale.

These are all standard SM observations, not GeoVac-internal. Sprint W3 falsified the master Mellin engine hypothesis for Wolfenstein parameters. The GeoVac-internal question is whether any structure beyond standard SM relations exists.

**No GeoVac-internal structure identified beyond standard SM RGE.**

### 5.3 §III.28 apparatus identity — Born rule

Per memory `external_input_three_class_partition.md` updated 2026-05-26: Born rule is calibration data (Class 1) at the substrate level. Framework uses Gleason's theorem via Hilbert-space inheritance; no substrate-level derivation. NO internal structure to extract — the Born rule is the structure.

## 6. Five potential follow-on diagnostic angles

If the W3 question is reopened, the projection-specific framing identifies these as the natural targets:

| Angle | Cost | Risk | Value |
|---|---|---|---|
| Koide-like relations on quark mass triples (u,c,t; d,s,b) | 1 day | Low | Sprint W3 lepton-mass recheck partially covered this; likely NEGATIVE |
| Koide-like relations on Wilson gauge couplings (g₁, g₂, g₃) | 1 day | Low | Known SM territory; standard RGE constraints |
| Koide-like relations on nucleon mass triples (m_p, m_n, m_α) | 1 day | Low | Nuclear physics territory; not GeoVac-internal |
| Koide-like on charged-lepton CHARGES (e, μ, τ all carry charge -e) | 0 | trivial | All identical, no test |
| Paper 2 K-formula analog for OTHER calibration constants | 1 week | Medium | Sprint A eliminated 9+ mechanisms; specific projection-K targets exhausted |

**None of these is high-value.** The Koide observation is exhausted at the lepton level; the K-formula is exhausted at the α level; standard SM RGE covers gauge couplings.

## 7. Substantive structural finding

The projection-specific framing produces a sharp observation worth keeping in CLAUDE.md §1.7 (WH5 sibling):

> **Among GeoVac's 28 Paper-34 projections, exactly TWO carry parametric Class-D calibration data with known internal structure: §III.14 (Koide cone for charged leptons) and §III.16 (Paper 2 K-formula for α). Three additional projections (§III.6 spectral action, §III.10 Wilson loop, §III.28 apparatus identity) carry Class-D calibration data without GeoVac-internal structure. The remaining 23 projections carry no parametric calibration to test. This concentration — 2-of-28 projections showing internal calibration structure — supports the structural-skeleton scope reading: calibration data is genuinely external for the overwhelming majority of projections; the two known exceptions (Koide, Paper 2 K) are framework observations that have already been investigated and remain open as numerical observations without derivation.**

This is a stronger form of the structural-skeleton scope statement than what's currently in `geovac_structural_skeleton_scope_pattern.md`. It can sharpen the memory file and Paper 18 §IV.6 (inner-factor input data tier).

## 8. Verdict

**NO-GO for opening a projection-specific calibration parameter sprint.**

Substantive findings worth preserving:

1. **Projection-specific framing is structurally clean.** 28 projections → 5 with Class-D calibration data → 2 with known internal structure. This concentration is informative.

2. **The two structured calibration projections (§III.14, §III.16) are already framework-known observations** (Koide cone, Paper 2 K-formula). Neither has been derived; both remain numerical observations.

3. **The three remaining Class-D projections (§III.6, §III.10, §III.28) have no GeoVac-internal structure** beyond what standard physics provides.

4. **The structural-skeleton scope reading is REINFORCED at the projection level.** 26 of 28 projections have either no Class-D calibration or have calibration without GeoVac-internal structure. The framework's external-input partition holds at the projection-by-projection level.

## 9. Recommendations

1. **Do NOT open a projection-by-projection W3-style PSLQ sprint.** The W3 mechanical-basis sprint of 2026-05-08 already covered the global question, and the projection-specific reframe sharpens *where* to look but doesn't change the verdict.

2. **Add the 2-of-28 observation** from §7 to Paper 18 §IV.6 (inner-factor input data tier) or as a new structural reading in CLAUDE.md §1.7. This sharpens the structural-skeleton statement at the projection level.

3. **Park §III.14 (Koide) and §III.16 (Paper 2 K-formula) as the canonical W3 targets** — if the second packing axiom question is ever reopened, these are the two existing observations to address.

4. **No new sprints recommended.** The projection-specific scoping CLOSES the W3 reframe as a productive direction.

## 10. Cross-references

- `debug/multifocal_b_w3_diag_memo.md` — 14-catalog of "second packing axiom" speculations (Phase B-W3-diag, 2026-05-07)
- `debug/w3_mechanical_basis_memo.md` — Sprint W3 mechanical-basis falsification of Wolfenstein hypothesis
- `debug/w3_lepton_mass_recheck_memo.md` — lepton-mass W3 sibling track
- `memory/w3_spectral_zeta_candidate.md` — W3 falsification memory
- `memory/koide_cone_clean_negative.md` — Koide-specific PSLQ negative
- `memory/external_input_three_class_partition.md` — Class 1 / 2 / 3 calibration partition
- `memory/geovac_structural_skeleton_scope_pattern.md` — structural-skeleton scope pattern
- Paper 34 §III.1-§III.28 — projection inventory
- Paper 18 §IV.6 — inner-factor input data tier
- `papers/group5_qed_gauge/paper_2_alpha.tex` — Paper 2 K-formula conjecture

## 11. Files

- `debug/sprint_projection_specific_calibration_scoping_memo.md` (this)
- No driver script (diagnostic-only)
- No data file (no computation)
