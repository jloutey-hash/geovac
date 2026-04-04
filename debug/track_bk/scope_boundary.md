# Track BK: GeoVac Atomic Scope Boundary

## Overview

Assessment of which atoms are reachable by the GeoVac composed architecture,
given the current natural geometry hierarchy and solver infrastructure.

## First Row (Z=1-10): Fully Supported

| Atom | Config | Core | Valence | GeoVac Status |
|------|--------|------|---------|---------------|
| H (1)  | 1s^1 | none | 1e | Level 1 (exact) |
| He (2) | 1s^2 | none | 2e | Level 3 (0.05%) |
| Li (3) | [He] 2s^1 | 1s^2 | 1e | Composed L3+L2 (R_eq 5.3%) |
| Be (4) | [He] 2s^2 | 1s^2 | 2e | Composed L3+L4 (R_eq 11.7% in BeH2) |
| B (5)  | [He] 2s^2 2p^1 | 1s^2 | 3e | Feasible but untested |
| C (6)  | [He] 2s^2 2p^2 | 1s^2 | 4e | Feasible but untested |
| N (7)  | [He] 2s^2 2p^3 | 1s^2 | 5e | Feasible but untested |
| O (8)  | [He] 2s^2 2p^4 | 1s^2 | 6e | Composed (H2O, R_eq 26%) |
| F (9)  | [He] 2s^2 2p^5 | 1s^2 | 7e | Feasible but untested |
| Ne (10)| [He] 2s^2 2p^6 | 1s^2 | 8e | No molecular interest |

**Core treatment:** All first-row atoms (Li-Ne) have a He-like 1s^2 core,
which the Level 3 solver handles exactly. PK parameters are derived ab initio
from the He-like core solution (Paper 17).

**Valence treatment:** The composed architecture groups valence electrons into
pairs (bond pairs, lone pairs). Each pair is treated as a Level 4 block with
Z_eff screening. The number of valence blocks scales linearly with valence
electron count.

**What works:**
- Core screening: Z_eff(r) from algebraic_zeff.py
- PK pseudopotential: Z^2-scaled from Li^2+ reference
- Bond pair coupling: Slater F^0 integrals (algebraic)
- Qubit encoding: O(Q^2.5) Pauli scaling confirmed for Li, Be, O

**What is missing (for B-F):**
- Orbital filling order: PK does not reproduce s-p splitting (see
  filling_order_investigation.md). The basis construction must manually
  assign valence configurations.
- Open-shell treatment: B, C, N have partially filled 2p subshells.
  The current composed architecture assumes closed-shell pairs. Open-shell
  atoms would need unpaired electron treatment (single-electron blocks).
- Valence correlation: at Z_eff >= 4 (C-Ne), the PK barrier becomes very
  large (see Track BK data: 67-284 Ha). The Z^2-scaled PK may need
  recalibration from direct Be^3+, B^4+, etc. core solutions rather than
  extrapolation from Li^2+.

**Assessment:** First-row atoms are fully within architectural scope. The
composed infrastructure handles any 2-electron core + N-valence-electron
system. The practical bottleneck is testing and validation, not fundamental
limitations. The PK Z^2 scaling should be validated against direct
CoreScreening solutions for Z >= 5.

## Second Row (Z=11-18): Partially Feasible

| Atom | Config | Core | Valence | Status |
|------|--------|------|---------|--------|
| Na (11) | [Ne] 3s^1 | 1s^2 2s^2 2p^6 | 1e | Needs 10e core |
| Mg (12) | [Ne] 3s^2 | 1s^2 2s^2 2p^6 | 2e | Needs 10e core |
| Al (13) | [Ne] 3s^2 3p^1 | 10e | 3e | Needs 10e core |
| Si (14) | [Ne] 3s^2 3p^2 | 10e | 4e | Needs 10e core |
| P (15)  | [Ne] 3s^2 3p^3 | 10e | 5e | Needs 10e core |
| S (16)  | [Ne] 3s^2 3p^4 | 10e | 6e | Needs 10e core |
| Cl (17) | [Ne] 3s^2 3p^5 | 10e | 7e | Needs 10e core |
| Ar (18) | [Ne] 3s^2 3p^6 | 10e | 8e | No molecular interest |

**The 10-electron core problem:**

The first-row composed architecture uses a 2-electron (He-like) core solved
by the Level 3 (hyperspherical) solver. Second-row atoms have a 10-electron
[Ne] core: 1s^2 + 2s^2 + 2p^6.

Options for handling this:

1. **Frozen-core tabulation.** Solve the 10-electron Ne-like core once (using
   the full N-electron solver or external data), tabulate its properties
   (energy, density, Z_eff(r), PK parameters), and use these as inputs to
   the composed valence solver. This is the most practical approach.
   - Pro: Reuses existing valence infrastructure unchanged.
   - Con: The 10-electron core solution itself is expensive. The current
     N-electron solver (Level 4N) handles 4 electrons at l_max=2 with
     difficulty (Track AK: 750 spectral dim). A 10-electron version would
     require SO(30) angular machinery.
   - Feasibility: **Yes, if core properties are tabulated externally** (e.g.,
     from NIST atomic data or Hartree-Fock calculations). The GeoVac valence
     solver only needs E_core, n_core(r), and the resulting Z_eff(r) and PK
     parameters — it does not need to solve the core quantum mechanically.

2. **Recursive composition.** Treat the 10-electron core as nested composed
   blocks: 1s^2 (Level 3) + 2s^2 (Level 4 with PK from 1s^2 core) + 2p^6
   (three Level 4 pairs with PK from 1s^2+2s^2).
   - Pro: Fully ab initio, no external data needed.
   - Con: The 2p^6 subshell has 3 pairs in 3 spatial orientations (m=-1,0,+1).
     These pairs interact strongly via exchange, which the current composed
     architecture handles poorly at high Z_eff (Track BK: lone pair coupling
     is unphysical at Z_eff >= 6 due to overly large Slater integrals).
   - Feasibility: **Unlikely to be rigorous.** The 2s-2p exchange coupling
     within the [Ne] core is a major correlation effect that recursive
     composition would struggle to capture.

3. **Effective core potential (ECP) from literature.** Use published ECPs
   (e.g., Stuttgart/Cologne) for the [Ne] core and focus GeoVac on the
   valence electrons only.
   - Pro: Standard approach in quantum chemistry; well-validated ECPs exist.
   - Con: Introduces external parameters, losing the ab initio character.
   - Feasibility: **Straightforward** if the goal is practical computation
     rather than foundational purity.

**Assessment:** Second-row atoms are feasible with approach (1) or (3) but
require external core data. A fully ab initio treatment via recursive
composition (approach 2) is unlikely to work without significant advances
in the many-electron solver. Near-term recommendation: tabulate [Ne] core
properties from NIST/HF data and use the existing composed valence
infrastructure.

## Transition Metals (Z=21-30): Out of Scope

| Atom | Config | Core | Issue |
|------|--------|------|-------|
| Sc (21) | [Ar] 3d^1 4s^2 | 18e | 3d/4s near-degeneracy |
| Ti-Zn   | [Ar] 3d^n 4s^m | 18e | Multi-reference, strong correlation |

**Why transition metals are out of scope:**

1. **18-electron core.** The [Ar] core requires solving (or tabulating) an
   18-electron system. The recursive composition challenges from the
   10-electron case are compounded.

2. **3d/4s near-degeneracy.** The Madelung rule predicts 4s fills before 3d,
   but 3d and 4s are nearly degenerate. As shown in this investigation, PK
   cannot reproduce the s-d energy ordering (wrong sign, wrong magnitude).
   The composed architecture would need explicit 3d-4s correlation, which
   is a multi-reference problem.

3. **Strong correlation.** Transition metal chemistry involves partially
   filled d-shells with strong electron correlation. This is one of the
   hardest problems in quantum chemistry and is far beyond the current
   composed architecture's capabilities.

4. **Spin-orbit coupling.** For heavier transition metals (Z > 30), spin-orbit
   effects become significant. The current framework is non-relativistic.

**Assessment:** Transition metals are clearly out of scope for the foreseeable
future. They would require: (a) a production-quality many-electron solver for
18+ electrons, (b) multi-reference treatment of 3d/4s near-degeneracy, and
(c) possibly relativistic corrections.

## Verdict: Near-Term Reachability

| Category | Atoms | Status | Bottleneck |
|----------|-------|--------|------------|
| Fully operational | H, He, Li, Be, O (in H2O) | Production | Validation |
| Architecturally ready | B, C, N, F | Need testing | PK validation at Z >= 5 |
| Feasible with external data | Na-Ar | Approach (1) or (3) | 10e core tabulation |
| Out of scope | Z > 20 | Fundamental | Multi-reference, many-electron core |

**Recommendation:** The near-term priority is completing first-row validation
(B-F atoms in molecules). This requires no architectural changes — only
testing the composed infrastructure at higher Z_eff and validating PK
parameters against direct CoreScreening solutions. Second-row support is a
medium-term goal requiring a design decision on core treatment (tabulation
vs. ECP).
