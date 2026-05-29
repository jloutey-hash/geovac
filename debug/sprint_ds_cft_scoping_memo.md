# Sprint dS/CFT scoping — where does GeoVac sit relative to the dS/CFT correspondence?

**Date:** 2026-05-29
**Path:** Continuation thread, Task 9. Diagnostic-only.
**Verdict:** **GO — publishable framework-positioning contribution identified.** GeoVac is NOT dS/CFT (no bulk dS gravity in the framework), and NOT AdS/CFT (no anti-de Sitter bulk), but is structurally adjacent to BOTH at the boundary-CFT level. The S³ substrate sits exactly where Euclidean dS/CFT (Maldacena 2002) places its boundary CFT₃; Paper 50's bit-exact KPS matches reproduce standard dS-boundary CFT data via a discrete spectral-triple mechanism that is structurally distinct from any holographic-bulk construction. **The framework's actual position is: discrete spectral-triple realization of boundary CFT₃-on-S³, without bulk dual.** This is sui generis and worth articulating publicly.

## 1. The question

PI question from the conversation thread opening this multi-task work: *"could we be a DS/CFT correspondence? Like without the A."*

This scoping pass tests whether the question has a substantive answer beyond "no, we don't have a bulk." Specifically:
- What does Strominger 2001 / Maldacena 2002 dS/CFT formulation require?
- What does GeoVac actually have (Paper 50 + gravity arc + spectral triple machinery)?
- Where do the structural overlaps and gaps land?
- Is there a publishable framework-positioning contribution distinct from Paper 50's existing scope?

## 2. dS/CFT requirements (Strominger 2001, Maldacena 2002)

### 2.1 Strominger 2001 Lorentzian formulation

Strominger~\cite{strominger2001} proposed a duality between gravity on $d+1$-dimensional de Sitter space dS_{d+1} and a Euclidean CFT_d on its conformal boundary. The structural ingredients:

- **Bulk:** dS_{d+1} with metric $ds^2 = -d\tau^2 + \cosh^2(\tau)\,d\Omega_d^2$ in global coordinates. For $d=3$: dS_4 with spatial slices $S^3$ growing exponentially in $\tau$.
- **Boundaries:** future ($I^+$) and past ($I^-$) timelike infinities are both $S^d$.
- **Boundary CFT:** Euclidean CFT_d on $I^+ \cong S^d$.
- **Dictionary:** bulk fields $\phi$ correspond to boundary operators $\mathcal{O}_\Delta$ via $\phi(\tau \to \infty, \Omega) \sim e^{-\Delta \tau} \mathcal{O}_\Delta(\Omega)$.
- **Partition function identification:** $Z_{dS}[\phi_0] = \langle e^{-\int \phi_0 \mathcal{O}}\rangle_{CFT_d}$.

For our case ($d = 3$): the conjectured boundary CFT_3 lives on $S^3$.

### 2.2 Maldacena 2002 Euclidean sharpening

Maldacena~\cite{maldacena2002} sharpened via Wick rotation:

- **Euclidean dS_{d+1}:** $S^{d+1}$ (round sphere).
- **Boundary location:** equatorial $S^d$ at imaginary time $\tau = 0$.
- **Wavefunction of the universe:** $\Psi_{dS}[\phi_0] = Z_{CFT_d}[\phi_0]$ where the LHS is the Hartle-Hawking-type wavefunctional at $I^+$ and the RHS is the boundary CFT partition function.
- **For dS_4:** $\Psi_{dS_4}[\phi_0] = Z_{CFT_3-on-S^3}[\phi_0]$.

The structural identification is then: the CFT_3-on-S^3 partition function IS the Hartle-Hawking wavefunctional of dS_4 at I^+.

### 2.3 Requirements for a "dS/CFT correspondence" claim

| # | Requirement | Why it matters |
|---|---|---|
| R1 | Bulk dS_{d+1} or Euclidean S^{d+1} | The "dS" side of dS/CFT |
| R2 | Boundary CFT_d on the equator/I^+ S^d | The "CFT" side of dS/CFT |
| R3 | Holographic dictionary mapping bulk fields ↔ boundary operators | The "correspondence" |
| R4 | Partition function identification $Z_{dS} = Z_{CFT}$ | The operational content of the correspondence |
| R5 | Witten functional integration over bulk geometries | The dynamical generating function |

## 3. What GeoVac actually has

### 3.1 Boundary CFT_3-on-S^3 side (POSITIVE)

Paper 50 establishes:
- **R2 fully satisfied** at boundary side: CFT_3-on-S^3 partition functions reproduced bit-exactly for conformally coupled scalar and free Weyl Dirac (KPS~\cite{klebanov_pufu_safdi2011} match at 61+ digits with PSLQ integer relation $[8, 2, -3]$ for scalar; symbolic-exact match for Dirac).
- **F-theorem coefficient reproduced**: $F_s = (\log 2)/8 - 3\zeta(3)/(16\pi^2)$ and $F_D = (\log 2)/4 + 3\zeta(3)/(8\pi^2)$ — the standard KPS $F$-coefficients.
- **Master Mellin engine M2/M3 decomposition** (Paper 50 §3.4): scalar/Dirac combinations project orthogonally onto M2 (Seeley-DeWitt) and M3 (vertex-parity) sub-mechanisms.

The boundary CFT side is structurally present at theorem-grade rigor.

### 3.2 Bulk dS or Euclidean S^4 side (BLOCKED)

Paper 50 §6 Proposition `prop:bulk_blocked` explicitly states:
- No AdS_4 / H^4 machinery in `geovac/` codebase.
- Sprint RH-B 2026-04-17 closed Wick continuation $S^3 \to H^3$ as clean dead end (no natural discrete subgroup $\Gamma \subset SO(3,1)$ from framework invariants).
- Sprint L3e-P3 2026-05-23: Paper 38's $4/\pi$ rate does not transport to non-compact Coulomb setting.

For DE SITTER side specifically (extending Paper 50's AdS/H discussion):
- No Euclidean S^4 machinery in `geovac/`.
- The framework's $S^3$ is the Fock-projection of the Coulomb spectrum (Paper 7), not the spatial slice of any 4D bulk.
- The cigar geometry of the gravity arc (Euclidean Schwarzschild) is Euclidean black hole, not Euclidean de Sitter — a different geometry.
- Hartle-Hawking wavefunctional construction would require functional integration over 4D geometries; the framework integrates over spectral-triple labels (n, l, m, s), not over metrics.

R1, R3, R4, R5 are all BLOCKED at the framework's existing infrastructure.

### 3.3 Cigar + Hawking content (ADJACENT but DISTINCT)

The gravity arc has cigar geometry, Bekenstein-Hawking entropy, four-witness Wick-rotation theorem. Could this give the dS bulk side via continuation?

**No, structurally.** The cigar is Euclidean Schwarzschild (a Euclidean BLACK HOLE), which is distinct from Euclidean dS (which is $S^{d+1}$). The two Euclidean geometries share the property of having a thermodynamic temperature via cigar tip regularity (Hawking $T_H = 1/(8\pi M)$, dS $T_{dS} = 1/(2\pi \ell)$), but the bulk topologies are different ($\mathbb{R}^2 \times S^2$ for Schwarzschild, $S^4$ for dS). The four-witness theorem unifies them at the modular-flow / Tomita-Takesaki level (Paper 42, 43), but does not provide a single bulk geometry that could serve as the "dS side" of a dS/CFT correspondence.

## 4. The structural overlap and gap

### 4.1 Where GeoVac IS dS/CFT-adjacent

Three structural overlaps:

**O1: The S³ substrate sits where the dS/CFT boundary CFT would sit.** Euclidean dS_4 has equatorial S³ as the boundary CFT_3 host (per Maldacena 2002). GeoVac's S³ substrate is structurally the same manifold — same topology, same metric (with appropriate radius choice), same SO(4) symmetry. The Fock-projection identification (Paper 7) gives the same S³ that dS/CFT would name as its boundary.

**O2: Paper 50's bit-exact CFT_3 reproduction IS standard dS-boundary CFT data.** The F-coefficient is the standard KPS F-coefficient that monotonically decreases under RG (F-theorem). The scalar and Dirac partition function values match standard continuum CFT_3 results. If a holographic dS/CFT computation produced these CFT_3 values, the bit-exact match would be the verification. GeoVac produces them directly from spectral data.

**O3: The wedge KMS state (Paper 50 §4) plus four-witness theorem (Paper 42) gives the modular structure that BW + Hawking + Unruh + Sewell unify.** This modular structure IS the boundary side of any holographic correspondence involving thermal states.

### 4.2 Where GeoVac is NOT dS/CFT

Three categorical differences:

**D1: No bulk dS geometry.** The framework does NOT contain dS_4 (Lorentzian) or S^4 (Euclidean) as a manifold. There is no "bulk" to be holographic FROM.

**D2: No holographic dictionary.** GeoVac produces CFT_3 partition function data DIRECTLY from the discrete spectral-triple structure on S^3. No bulk fields are mapped to boundary operators because there are no bulk fields.

**D3: No Hartle-Hawking wavefunctional.** dS/CFT's operational content is $\Psi_{dS}[\phi_0] = Z_{CFT_3}[\phi_0]$ — a functional integral identification. GeoVac's $Z_{CFT_3}$ is a SPECTRAL trace on the discrete substrate, not a functional integral.

### 4.3 The structurally distinct position

GeoVac sits in a **third category** distinct from both AdS/CFT and dS/CFT:

- **AdS/CFT:** bulk anti-de Sitter gravity ↔ boundary CFT, holographic dictionary, Witten functional integral.
- **dS/CFT:** bulk de Sitter gravity ↔ I^+ boundary CFT (or equator CFT in Euclidean), holographic dictionary, Hartle-Hawking wavefunctional.
- **GeoVac:** discrete spectral triple on S^3 ↔ continuum CFT_3 on S^3 (Latrémolière propinquity convergence per Paper 38), spectral-triple discretization, propinquity-style convergence theorem.

The third category is structurally distinguished by:
- The "bulk" replaced by a discrete spectral-triple substrate (operator-system level, finite cutoff).
- The "correspondence" replaced by propinquity convergence in the continuum limit (Paper 38 5-lemma proof).
- The boundary CFT data emerges from spectral asymptotics (zeta-function derivatives at zero), not from holographic reconstruction.

## 5. Is there a publishable contribution?

**YES.** Specifically:

### 5.1 Framework-positioning paper candidate

Title (working): *"Discrete spectral realization of CFT_3-on-S^3 partition functions: a non-holographic alternative to dS/CFT for boundary CFT data."*

Headline content:
- (A) Paper 50's bit-exact CFT_3 matches REPRODUCE the boundary data that dS/CFT (Maldacena 2002) would predict via Hartle-Hawking, WITHOUT invoking any bulk dS geometry.
- (B) The mechanism is propinquity convergence (Paper 38) of the discrete spectral triple, NOT holographic reconstruction. This gives an operationally distinct realization.
- (C) The wedge KMS / Connes-Rovelli thermal-time structure (Papers 42, 49) corresponds to the boundary-side modular structure that dS/CFT, AdS/CFT, and rigorous-thermal-QFT all share.
- (D) The framework's structural-skeleton scope (this thread's Task 2 sharpening: 2-of-28 projections carry internal calibration) extends naturally to "no bulk geometry inferred from boundary CFT data" as a structural feature, not a defect.

This is a sui generis position. It explicitly distinguishes GeoVac from both AdS/CFT and dS/CFT while clarifying where the framework's boundary CFT_3 data fits.

### 5.2 What this paper would NOT claim

Honest scope boundaries:
- Not a new holographic correspondence (no bulk = no holography).
- Not a derivation of dS/CFT from GeoVac (we don't have the bulk to compare).
- Not a refutation of dS/CFT (the boundary data is consistent; we just produce it differently).
- Not a Theory of Everything candidate (per CLAUDE.md §1.5 rhetoric rule).

### 5.3 Estimated effort

- Drafting: 1-2 weeks at the Papers 47/48 cadence (math.OA style standalone, ~15-20 pages).
- Concurrent-work check: needed to verify no published paper has made the same "discrete CFT-on-S^3 without bulk" claim. Tentatively unique given Paper 50's CONCURRENT_WORK_CLEAR result (May 2026).
- Risk: medium. The contribution is positional/structural rather than producing new theorems; venues (math.OA, math-ph, hep-th) would receive it depending on framing.

### 5.4 Alternative: Paper 50 §8 extension

If a standalone paper is too much commitment, Paper 50 §8 already has a catalogue and open targets section. A new subsection §8.X "Structural distinction from dS/CFT" could capture the headline of §5.1 in 1-2 pages, contained within Paper 50.

Cost: ~1-2 days. Risk: low. This is the cheap closure.

## 6. Verdict

**GO at GO-LOW-COST.** The structural distinction from dS/CFT is real, articulable, and worth recording publicly. Recommend Paper 50 §8 extension (cheap closure, ~1-2 days) as the first move; full standalone paper as a follow-on if PI wants to invest more.

## 7. Cross-references

- Paper 50 — main paper for CFT_3-on-S^3 bit-exact matches
- Paper 50 §6 — bulk-side blocked statement (AdS analog)
- Paper 50 §8 — catalogue and open targets (natural extension venue)
- Paper 38 — Latrémolière propinquity convergence (the "correspondence" mechanism)
- Paper 42, 43 — four-witness Wick-rotation theorem (boundary modular structure)
- Paper 49 — Connes-Rovelli thermal-time stack
- Strominger 2001 (arXiv:hep-th/0106113) — dS/CFT proposal
- Maldacena 2002 (arXiv:astro-ph/0210603) — Euclidean dS/CFT sharpening
- Klebanov-Pufu-Safdi 2011 — F-theorem reference for the bit-exact CFT_3 data

## 8. Files

- `debug/sprint_ds_cft_scoping_memo.md` (this)
- No driver script (diagnostic-only)
- No data file (no computation)
