# -*- coding: utf-8 -*-
"""CLAUDE.md compaction round 2 (2026-06-10, v3.110.1, PI-authorized).

Moves verbatim content out of the always-loaded file into docs/ homes, splices in
compact replacements. Zero deletion: every removed byte lands in an archive doc.

Sections treated:
  S1.5  two chronicle paragraphs -> docs/development_frontier_archive.md (compact kept)
  S1.7  full WH register -> docs/wh_register_history.md (compact claim/falsifier/status kept)
  S2    legacy bullets + oversize prose -> docs/development_frontier_archive.md
  S6    full paper notes -> docs/paper_notes_archive.md (loading tiers + flags kept)
  S7    full entry-point catalogue -> docs/code_architecture.md (top-10 stub kept)
  S12   full algebraic registry -> docs/algebraic_registry.md (stub kept)
  S13.11 rule 9 added (replace-don't-append), PI-authorized
"""
import io, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CM = ROOT / "CLAUDE.md"
DOCS = ROOT / "docs"

lines = CM.read_text(encoding="utf-8").splitlines()


def find(prefix, start=0):
    for i in range(start, len(lines)):
        if lines[i].startswith(prefix):
            return i
    raise SystemExit("MARKER NOT FOUND: " + prefix)


i15 = find("## 1.5. Positioning")
i16 = find("## 1.6. Project Phase")
i17 = find("## 1.7. Working Hypotheses")
i18 = find("## 1.8. Multi-Observable")
i2 = find("## 2. Current Development Frontier")
i3 = find("## 3. Approaches That Failed")
i6 = find("## 6. Paper Series")
i7 = find("## 7. Code Architecture")
i8 = find("## 8. Coding Standards")
i12 = find("## 12. Algebraic Registry")
i13 = find("## 13. Multi-Agent Protocol")

sec15 = lines[i15:i16]
sec17 = lines[i17:i18]
sec2 = lines[i2:i3]
sec6 = lines[i6:i7]
sec7 = lines[i7:i8]
sec12 = lines[i12:i13]

HDR = ("> Extracted verbatim from CLAUDE.md on 2026-06-10 (v3.110.0 state) during "
       "compaction round 2. This file is the canonical archive; the CLAUDE.md "
       "section now holds only the compact working form.\n")

# ---- archives -------------------------------------------------------------
def write_doc(name, title, chunks):
    p = DOCS / name
    if p.exists():
        raise SystemExit("REFUSING TO OVERWRITE EXISTING: " + str(p))
    buf = ["# " + title, "", HDR]
    for label, body in chunks:
        buf.append("## " + label)
        buf.append("")
        buf.extend(body)
        buf.append("")
    p.write_text("\n".join(buf), encoding="utf-8")
    print("wrote", p.name, sum(len(x) for x in buf) // 1024, "KB")


# S1.5 chronicle paragraphs
qsp = find("**Quantum simulation positioning:**", i15)
pk = find("**PK limitation and research path:**", i15)
write_doc("development_frontier_archive.md",
          "Development frontier archive (CLAUDE.md S1.5 chronicle + S2, frozen at v3.110.0)",
          [("S1.5 Quantum simulation positioning (verbatim)", [lines[qsp]]),
           ("S1.5 PK limitation and research path (verbatim)", [lines[pk]]),
           ("S2 Current Development Frontier (verbatim)", sec2)])

write_doc("wh_register_history.md",
          "Working Hypotheses register — full chronicle (CLAUDE.md S1.7, frozen at v3.110.0)",
          [("S1.7 verbatim", sec17)])

write_doc("paper_notes_archive.md",
          "Paper Series detailed notes (CLAUDE.md S6, frozen at v3.110.0)",
          [("S6 verbatim", sec6)])

write_doc("code_architecture.md",
          "Code architecture — full entry-point catalogue (extracted from CLAUDE.md S7, live document)",
          [("S7 (update here, not in CLAUDE.md)", sec7)])

write_doc("algebraic_registry.md",
          "Algebraic registry (extracted from CLAUDE.md S12, live document)",
          [("S12 (update here, not in CLAUDE.md)", sec12)])

# ---- compact replacements -------------------------------------------------
QSP_NEW = ("**Quantum simulation positioning:** The classical solver investigation (v2.0.6-23) "
           "characterized what the natural-geometry hierarchy can and cannot achieve for ground-state "
           "PES; the primary computational value proposition is quantum simulation. The composed "
           "architecture produces qubit Hamiltonians with O(Q^2.5) Pauli scaling (vs O(Q^3.9-4.3) "
           "Gaussian), 51x-1,712x fewer Pauli terms across LiH/BeH2/H2O, structural Gaunt sparsity "
           "compatible with all downstream optimizations (tapering, grouping, tensor factorization), "
           "and a 40-molecule library across three periodic-table rows with N_Pauli = 11.11 x Q "
           "universal for composed builds. Position for VQE/NISQ (Pauli count, QWC groups) rather than "
           "fault-tolerant QPE; LiH electronic-only 1-norm matches STO-3G at 0.97x with 2.7x fewer "
           "Pauli terms. Full track-by-track chronicle: `docs/development_frontier_archive.md` + CHANGELOG.")

PK_NEW = ("**PK limitation and research path:** The composed framework's PK pseudopotential is the "
          "accuracy bottleneck (5.3-26% R_eq error); six modification attempts failed (Section 3). The "
          "balanced coupled Hamiltonian (Papers 19/20) is the PK-free alternative: cross-center V_ne "
          "via multipole expansion with exact Gaunt termination; LiH 4e FCI reaches 0.20% energy at "
          "n_max=3 with structural R_eq drift (~8.8%) — energy converges excellently, geometry drifts; "
          "optimal for single-point quantum simulation at fixed geometries. Resource anchor: composed "
          "LiH 334 Pauli vs STO-3G 907 vs cc-pVDZ 63,519. Chronicle: `docs/development_frontier_archive.md`.")

lines[qsp] = QSP_NEW
lines[pk] = PK_NEW

SEC17_NEW = """## 1.7. Working Hypotheses (Internal Register)

This section is a **bold-claim register**, distinct from the rhetoric of the papers. Papers remain cautious under §1.5 (dual-description framing, no ontological priority); this register is what sub-agents and the PM may *reason from* during synthesis work. Nothing here appears in papers unless promoted after its falsifier clears. Full status chronicles (the sprint-by-sprint evidence trail, April–June 2026) live in `docs/wh_register_history.md`; this section holds only claim, falsifier, and current status.

**Governance:**
- PM may update a WH's "Status" line based on sprint evidence — REPLACE the line and move superseded text to `docs/wh_register_history.md`; never append (§13.11 rule 9).
- Adding, retiring, or promoting a WH to paper-level claim requires explicit PI direction.
- Retired WHs move to the history doc with rationale; never silently deleted.
- No WH is a license to bypass the rhetoric rule in papers or the verification gates in §13.4.

---

**WH1 — GeoVac is an almost-commutative spectral triple.** A = functions on the Fock-projected S³ graph; H = scalar/spinor state space; D = Camporesi–Higuchi Dirac; non-abelian gauge structure enters as inner derivations of the almost-commutative extension A ⊗ M_n(ℂ) (Marcolli–van Suijlekom lineage; Papers 25/30/32).
*Falsifier:* a GeoVac observable demonstrably inconsistent with any spectral-action expansion, or a violated structural axiom (order-one, reality).
*Status:* **PROVEN — unconditional (2026-06-10).** Paper 38: the discrete truncations converge to the round-S³ spectral triple in van Suijlekom's state-space GH distance at rate (4/π + o(1))·log n_max/n_max, on the truthful CH substrate (translation-seminorm metrization; frozen falsifier `tests/test_p38_action_seminorm.py`). The Lorentzian extension (Papers 45–49) is DESCOPED (P45 annihilation theorem); repair path = Toeplitz temporal compressions (see WH7).

**WH2 — Paper 18 is the Seeley-DeWitt + ζ-invariant decomposition of this spectral triple.** The transcendental taxonomy is the structured output of spectral-action geometry, organized by operator order × bundle type.
*Falsifier:* a transcendental in a GeoVac observable that cannot be placed in the grid.
*Status:* three of four axis-quadrants filled; mixed-Tate sharpening POSITIVE (2026-06-03) — M2 on S³ sits in the pure-Tate sub-ring ⊕_k π^{2k}·ℚ (Fathizadeh–Marcolli inherited). See `debug/sprint_mixed_tate_test_memo.md`.

**WH3 — The lattice exists a priori; match to physics is evidence, not derivation.** The packing construction is independent of known physics; its persistent match (Fock S³, nuclear magic numbers, Dirac fine structure, Pauli sparsity) is evidence physics is hosted by a discrete spectral triple, not that the lattice was reverse-engineered.
*Falsifier:* the packing construction fails to force the (n, l, m, s) structure or the n²−1 spectrum; or a match turns out to depend on a hidden physics-informed parameter.
*Status:* origin-story framing permits strong ontological claim internally; papers stay under §1.5.

**WH4 (deflated 2026-05-07) — The four-way S³ unity is one Fock-projection statement plus three forced consequences.** Bertrand + SO(4) force S³; the Hopf base, the CH spinor bundle, and SU(2) Wilson all follow from S³ = SU(2) parallelizability / maximal-torus structure. Does NOT extend to inner-factor selection (Yukawas, generations remain unforced).
*Falsifier:* a construction forcing one of the four roles onto a different manifold; or a published argument decoupling a consequence from the Fock input.
*Status:* deflated to a single-input forcing statement; outer structural unity essentially closed (Sprint TS-D + Paper 38).

**WH5 — α is a projection constant, not a derivable number.** K = π(B + F − Δ) composes three structurally independent spectral objects (finite Casimir trace; Fock Dirichlet ζ(2); Dirac boundary count 1/40); the right open question is why the sum equals α⁻¹, not how to derive each piece.
*Falsifier:* a spectral-triple construction deriving K as a single coefficient of a well-defined functional.
*Status:* TWELVE mechanisms eliminated (Phases 4B–4I + Sprint A + Sprint K-CC, including the T9 algebraic obstruction: no single CC heat-kernel expansion contains B, F, Δ as terms). Standing reading: three-regime projection coincidence. Paper 2 stays in Observations; combination rule conjectural (§13.5 hard prohibition).

**WH6 — GeoVac's RH-adjacent object is the Dirac spectral zeta D(s), not classical ζ.** Internal GUE-like zero statistics (CV ≈ 0.35–0.40); the classical-RH bridge is closed by three independent walls (zeros not on one line; no spectral-triple-natural functional equation, 48 OoM; wrong Weyl class).
*Falsifier:* D(s) zeros on a single critical line at larger samples; or a natural functional equation closing the RH-O gap.
*Status:* paused (2026-04-18); if resumed, the target is D(s) and its spectral-action interpretation.

**WH7 — Time-discreteness is observer-compactification (registered 2026-06-10, PI direction).** The only temporal structure the framework can see metrically is compactified time, and compactified time is automatically discrete. Inputs: (i) Paper 35 — π enters exactly at temporal compactification (Matsubara 2πk/β); (ii) P45 annihilation theorem — ℝ-time is Lipschitz-invisible in the v1 architecture; (iii) Paper 47 — the three temporal carriers are spectrally indistinguishable. Reading: discreteness of time is supplied by the observer's compact integration window (KMS β = 2π, four-witness theorem, Connes–Rovelli thermal time) — time is the prototype free-side projection. Temporal restriction of the organizing observation below.
*Falsifier (primary):* the Toeplitz temporal-compression program — a metrically visible temporal algebra on a NON-compact carrier without compactification weakens WH7 to convention; a proven annihilation-type obstruction for the framework's whole non-compact temporal class forces it.
*Falsifier (secondary, from Paper 35):* a GeoVac observable containing π whose evaluation provably involves no temporal/spectral integration.
*Status:* REGISTERED (2026-06-10). Honest cap: input (iii) means discrete-vs-continuous time may be empirically undecidable at every currently computable level; papers stay under §1.5.

---

**Organizing observation — discreteness is compactness (established).** The discrete, bit-exact skeleton is the closed/compact regime: compact groups have discrete spectra (Peter-Weyl), and this is the mathematical content of the S³ graph. The continuum is what remains when compactness is released. Calibration data is continuum (un-packed) data, which is why the skeleton is rational/forced and calibration is transcendental/free. Well-captured in Paper 18 §III (compactness thesis) and Paper 35 (temporal compactification injects π). The "second packing axiom" framing is retired; the 2026-05-30 confinement reframing is archived as an organizing reading (2026-05-31).

---
""".splitlines()

# S2 compact: keep note + June bullets (through v3.97.0) + best-results/key-results
# + all post-bullet lines <= 1500 chars (drops the four RH-sprint monsters and
# the long QC/balanced status paragraphs, all archived verbatim above).
note_end = i2 + 1
while not lines[note_end].startswith("- **"):
    note_end += 1
bullets, k = [], note_end
while k < i3 and (lines[k].startswith("- **") or lines[k].strip() == ""):
    if lines[k].startswith("- **"):
        bullets.append(lines[k])
        if "v3.97.0" in lines[k]:
            k += 1
            break
    k += 1
idx_best = find("**Best results by system type:**", i2)
post = [ln for ln in lines[idx_best:i3] if len(ln) <= 1500]
dropped = sum(1 for ln in lines[idx_best:i3] if len(ln) > 1500)

SEC2_NEW = (lines[i2:note_end]
            + bullets
            + ["",
               "> Older sprint index (v2.x–v3.96.0), the long-form arc chronicles, and the RH "
               "sprint records moved verbatim to `docs/development_frontier_archive.md` "
               "(2026-06-10 compaction). CHANGELOG.md remains the canonical chronicle going forward.",
               ""]
            + post
            + [""])
print("S2: kept", len(bullets), "bullets; dropped", dropped, "oversize prose lines to archive")

SEC6_NEW = """## 6. Paper Series

> Newcomer-facing status map with one-liners: `papers/INDEX.md`. Detailed per-paper notes (the long key-result descriptions, frozen at v3.110.0): `docs/paper_notes_archive.md`. Topic → paper lookup: `docs/topic_to_paper_lookup.md`. This section keeps only what agents need at dispatch time: loading tiers, the folder map, and live status flags.

### Loading tiers

**Always load** (framework identity): Paper 0 (packing axiom, K = −1/16) · 1 (spectral graph methods) · 7 (S³ proof, 18 symbolic proofs) · 14 (qubit encoding headline) · 16 (S_N periodicity) · 22 (angular sparsity theorem) · 23 (nuclear hub, Fock rigidity) · 24 (Bargmann-Segal S⁵, Coulomb/HO asymmetry) · 27 (entropy as projection) · 31 (universal/Coulomb partition) · 32 (the spectral triple; §VIII theorems).

**Load on topic** (full list and statuses in `papers/INDEX.md`): chemistry solvers → 8–9, 11, 12, 13, 15, 17, 19, FCI-A/M; QC resources → 20; QED/gauge/gravity → 2, 25, 28, 30, 33, 36, 41, 51; math.OA arc → 29, 38, 39, 40, 42–50, 52, 53; foundations/periods → 18, 54, 55, 56, 57; precision → 26, 34, 35.

**GUARDRAIL papers** — MUST load before any investigation in their domain (trigger words and protocol in §3.5): Papers 8–9 (single-center / unified-basis / shared-exponent molecular — Sturmian structural theorem), FCI-M (graph-concatenation molecular), Track DF record (nested hyperspherical).

### Folder organization (audience groups, reorganized 2026-05-22)

| Folder | Audience | Papers (.tex) |
|:-------|:---------|:-------------:|
| `papers/group1_operator_algebras/` | math.OA / NCG | 16 |
| `papers/group2_quantum_chemistry/` | quantum chemists | 9 |
| `papers/group3_foundations/` | mathematical physicists | 11 |
| `papers/group4_quantum_computing/` | QC / NISQ / VQE | 4 |
| `papers/group5_qed_gauge/` | HEP / gauge theory | 8 |
| `papers/group6_precision_observations/` | precision AMO | 4 |
| `papers/synthesis/` | cross-group narratives + field guide | 3 |
| `papers/archive/` | historical (3, 4, 5, 6, 10, 18v1, 21) | 7 |

### Live status flags every agent must know

- **Paper 38 UNCONDITIONAL** (2026-06-10, translation-seminorm metrization; falsifier `tests/test_p38_action_seminorm.py`). The WH1 keystone.
- **Papers 45/46 DESCOPED in place** (2026-06-09: P45 main theorem withdrawn — K⁺ seminorm annihilates; falsifier `tests/test_p45_kplus_degeneracy.py`). Do NOT cite their pre-descope claims. **Papers 47/48/49 PARTIAL** (in-paper Status notes; the norm-resolvent arrow and the TICI/cocycle algebra survive).
- **Paper 2 is an Observation**; the combination rule K = π(B + F − Δ) stays labeled conjectural (§13.5 hard prohibition).
- **Paper 34** is the living projection catalogue (28 projections); **Paper 18 §III.7** is the master Mellin engine; tag every transcendental against both (memory rule).
- Papers are corrected **in place** (de-versioning directive 2026-06-10); git/Zenodo are the version record. No splinter files.
"""
SEC6_NEW = SEC6_NEW.splitlines() + ["", "---", ""]

SEC7_NEW = """## 7. Code Architecture

> Full entry-point catalogue (~120 rows) and solver-method table: `docs/code_architecture.md` (live document — update there, not here). Most-used entry points:

| Task | Module · Entry point |
|:-----|:---------------------|
| Atomic lattice / Hamiltonian | `geovac/lattice.py` `GeometricLattice(Z, max_n)` · `geovac/hamiltonian.py` `GraphHamiltonian(lattice)` |
| Multi-electron FCI | `geovac/lattice_index.py` `LatticeIndex(Z, n_electrons, max_n)`; direct CI at N_SD ≥ 5000 |
| Molecular spec + composed builder | `geovac/molecular_spec.py` `MolecularSpec` · `geovac/composed_qubit.py` `build_composed_hamiltonian(spec)` + `*_spec()` factories |
| Balanced coupled builder | `geovac/balanced_coupled.py` `build_balanced_hamiltonian(spec, nuclei)` |
| Ecosystem export | `geovac/ecosystem_export.py` `hamiltonian(name, tapered=None/'global'/'per_block'/'extended'/'full')` → `.to_qiskit()/.to_openfermion()/.to_pennylane()` |
| Z₂ tapering | `geovac/z2_tapering.py` `apply_hopf_tapering()` · `geovac/extended_tapering.py` |
| Frozen cores / cross-center V_ne | `geovac/neon_core.py` `FrozenCore(Z)` · `geovac/shibuya_wulfman.py` `compute_cross_center_vne()` |
| Slater integrals (exact) | `geovac/hypergeometric_slater.py` `compute_rk_float()` (threshold-dispatched n ≥ 5 → exact Fraction) |
| Operator system / Connes distance / GH | `geovac/operator_system.py` `TruncatedOperatorSystem(n_max)` · `geovac/connes_distance.py` · `geovac/gh_convergence.py` |
| Physical constants | no central module; `-1/16` may be used directly (§8); `ALPHA`/`C_LIGHT` live next to their modules |
"""
SEC7_NEW = SEC7_NEW.splitlines() + ["", "---", ""]

SEC12_NEW = """## 12. Algebraic Registry

Tracks which matrix elements at each level are computed algebraically vs numerically. **Full registry** (Levels 2 / 3 / 4 / 4N / 5 + spin-ful Tier-2 tables): `docs/algebraic_registry.md` (live document — update statuses there, not here). Status vocabulary: **algebraic** (closed-form from quantum numbers) / **algebraic (implicit)** (defined by P = 0 with known coefficient ring; pointwise diagonalization is convenience, not necessity) / **algebraic-pending** (route identified, production still uses quadrature) / **numerical-required** (no known algebraic replacement). The §4 prime-directive test governs changes: anything touching quantum-number labels or selection rules is prohibited; improving radial-amplitude evaluation within a channel is legitimate.

---
"""
SEC12_NEW = SEC12_NEW.splitlines()

# Splice from the bottom up so earlier indices stay valid.
lines[i12:i13] = SEC12_NEW
lines[i7:i8] = SEC7_NEW
lines[i6:i7] = SEC6_NEW
lines[i2:i3] = SEC2_NEW
lines[i17:i18] = SEC17_NEW

# S13.11 rule 9 (re-find: indices shifted)
r8 = find("8. **One canonical record per fact.**")
lines[r8 + 1:r8 + 1] = ["",
    "9. **Status updates replace, never append (added 2026-06-10, PI-authorized).** When a WH status, "
    "a §6 status flag, or a paper-state description changes, REPLACE the existing text and move the "
    "superseded version to its history home (`docs/wh_register_history.md`, `docs/paper_notes_archive.md`, "
    "CHANGELOG.md). Chronicling-by-appending inside CLAUDE.md is the failure mode that regrew the file "
    "from 1,263 lines (2026-05-31 compaction) to 1,400 lines / 320 KB by 2026-06-10."]

out = "\n".join(lines) + "\n"
CM.write_text(out, encoding="utf-8")
print("CLAUDE.md:", len(lines), "lines,", len(out.encode('utf-8')) // 1024, "KB")
