# Sprint outreach polish (2026-06-04)

## 1. Setup

User context: after the v3.49.0 Q5'-scoping arc (period-ring deflated to MT(ℤ[i,1/2],4); cosmic-Galois U^* identified as the viable multi-year target) and the v3.48.0 identity-paper drafting (Field Guide), the project is at the point where outreach to external collaborators (Marcolli/Fathizadeh group, Latrémolière, Mondino, Hekkelman) is the natural next move. Before that wave goes out, the corpus needs polish: visible blemishes that would erode reviewer confidence (Paper 14's `\fbox` TODO placeholders), undefined citations/references that compile as `[?]`, and §13.4a equation-verification coverage gaps the followon register has been carrying.

This sprint closed five polish items from `debug/followon_register.md` that the user explicitly scoped to "not reach-out":

| Item | Source | Description |
|------|--------|-------------|
| C1 | v3.45.3 corpus cleanup | Paper 14 `\fbox` placeholders for two figures (`pauli_scaling.png`, `eri_density.png`) |
| A9 | Verification cleanup `wf_39fc232e-906` (2026-06-04) | 22 of 28 Paper 34 §III projections without numerical spot-check tests |
| C2 | Phase 3 physics audit (v3.45.2) | Paper 14 §sec:ft_gaussian BLISS-THC fault-tolerant comparison not built |
| C3 (NEW, surfaced by C1) | C1 compile audit | 8 missing bibitems + 3 missing section labels in Paper 14 |
| C4 (NEW, surfaced by A9 batch 1) | Batch 1 test ModuleNotFoundError | CLAUDE.md §7 claims `geovac/constants.py` exists; it doesn't |

The sprint is mechanical polish, not new physics. The discipline is reverse-direction: surface what the existing claims-on-paper say and verify each, fix what disagrees, document what's verified.

## 2. What was done

### 2.1 C1 — Paper 14 figures

The pre-existing `benchmarks/qubit_encoding/pauli_comparison.py` was originally written for an H₂-like Z=1 study and only emits `pauli_scaling.png`; Paper 14's Tables I and II use mixed Z=2 He data plus published Gaussian baselines, and the paper needs `eri_density.png` as a companion plot. Rather than retrofit a benchmark script with different physical scope, a dedicated figure-generator was created at `benchmarks/qubit_encoding/generate_paper14_figures.py` that hard-codes Paper 14's canonical Table I/II data points and produces both figures so the figures match the paper exactly. The extracted GeoVac scaling exponent from the table data is 3.147, matching Paper 14's claimed α=3.15 — a useful sanity check.

Both figures written to `papers/group4_quantum_computing/paper_14_figures/{pauli_scaling,eri_density}.png`. The two `\fbox` placeholders replaced with `\includegraphics`. Paper 14 three-pass clean compile: 21 pages.

### 2.2 A9 — Paper 34 §13.4a equation verification

Three batches of test files following the format of the original `tests/test_paper34_projection_spot_checks.py` (6 of 28 covered).

**Batch 1** (`test_paper34_projection_spot_checks_batch1.py`, 8 load-bearing projections, 31 tests):
- §III.1 Fock conformal — κ=-1/16 rational prefactor, Vol(S³)=2π², integer spectrum -(n²-1)
- §III.5 Sturmian — relabeling preserves Layer 1 rationality
- §III.11 Vector-photon — 1/(4π) = Vol(S²)/(4·4π²) closed form
- §III.13 Drake–Swainson — D_drake(n,l) = 2(2l+1)Z⁴/n³ combinatorial-rational
- §III.16 Two-body Dirac/Breit — closed-form Z=1 matrix elements via production `compute_radial`
- §III.17 Foldy/Friar — (2π/3) prefactor, 1s contact density π-cancellation
- §III.18 Zemach — linear-in-r_Z LO, profile independence at LO
- §III.19 Tensor multipole — rank-2 triangle inequality, m-conservation

**Batch 2** (`..._batch2.py`, 7 gauge/symmetry rows, 22 tests):
- §III.9 Wigner D — unitarity of d^l(β) at β ∈ {π/4, π/3, π/2}; Q[√2] ring at β=π/2
- §III.10 Wilson plaquette — maximal-torus diagonal SU(2) → U(1) action bit-exact reduction
- §III.20 Phillips-Kleinman — projector idempotency P_c² = P_c, no transcendental injected
- §III.21 Multipole/Gaunt termination — L > l₁+l₂ vanishes; LiH n_max=2 count check
- §III.22 Bipolar harmonic — triangle |k₁-k₂| ≤ K ≤ k₁+k₂ via Wigner 3j
- §III.23 Young tableau — integer characters of S_N, hook-length formula on S₄/S₅
- §III.26 Gauge choice — Coulomb-gauge identity 1/(4π) = Vol(S²)/(4·4π²), no variable introduced

**Batch 3** (`..._batch3.py`, 7 remaining rows, 11 tests):
- §III.3 Bargmann-Segal — π-free at finite N_max; Vol(S⁵) = π³ in continuum
- §III.4 Stereographic — chordal distance conformal factor identity
- §III.12 Mol-frame hyperspherical — Gaunt rationality at angular level via squared identity
- §III.15 Observation/temporal-window — Matsubara modes 2πk/β bosonic, (2k+1)π/β fermionic; ζ(4) = π⁴/90 for Stefan-Boltzmann
- §III.24 Adiabatic/Born-Oppenheimer — BO small parameter m_e/M_n, orthonormality of fast eigenvectors
- §III.25 Coupled-channel — linear matrix pencil H₀ + R·V^coupling gives polynomial P(R,μ) of degree l_max+1
- §III.28 Apparatus identity — von Neumann entropy log N at maximally mixed, Gibbs closed form for 2-level

**Coverage now 28 of 28 named projections.** Combined test count 108 (pytest-expanded with parametrize); raw test function count 64.

**Substantive correction surfaced inside batch 1.** Paper 34 §III.16 listed `R⁰_BP(1s,2s;1s,2s) = -4 log 2 - 19/9 + 9 log(3)/2` (numerically ~0.060). The production module `geovac.breit_integrals.compute_radial(1,0,2,0,1,0,2,0, k=0, kernel_type='breit', Z=1)` returns the pure rational `4/81 ≈ 0.0494`, and the pre-existing regression test `tests/test_breit_integrals.py:131-132` explicitly tests for this value with the comment "BP-retarded integral is a pure rational." So Paper 34's stated closed form was wrong, and the origin of the error traces to `debug/ps_1s2s_autopsy_track4_memo.md` line 81 (Sprint Calc-Ps-1S2S Track 4 autopsy, 2026-05-09 vintage), which itself listed the same wrong formula and inconsistent numerical value (0.060055).

Fixed per CLAUDE.md §13.8 PM authority to correct claims contradicted by computational evidence:
- Paper 34 `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` §III.16 — closed form rewritten as `R^0_BP(1s,2s;1s,2s) = 4/81` with annotation that the value is a pure rational, no log content, production-verified by the two named test files
- `debug/ps_1s2s_autopsy_track4_memo.md` line 81 corrected with a "CORRECTION 2026-06-04 (A9 batch 1)" provenance note preserving the original wrong entry's mention so the audit trail is intact
- Paper 34 three-pass clean compile after correction: 125 pages (unchanged)

### 2.3 C2 — Paper 14 BLISS-THC contextualization

The followon register's C2 description estimated 1-2 days because it required reading the external Caesura et al. 2025 paper (arXiv:2501.06165) carefully. The paper PDF was fetched via WebFetch and extracted via `pdftotext`. Caesura's Table II (P450, 63e/58o, Q=116) gives:

| Method | Memory qubits | 1-norm λ (Ha) | Toffolis |
|--------|---|---|---|
| DF [Goings] | 4,922 | — | — |
| SCDF | 1,706 | — | — |
| THC [Goings] | 1,434 | 388.9 | 7.79×10⁹ |
| THC re-estimated | 1,357 | — | — |
| **BLISS-THC** | **999** | **~133** | **1.71×10⁹** |
| Theoretical lower bound | — | 69.3 | — |

Speedup breakdown vs THC baseline: ~25× from AV photonic compilation + ~8× from BLISS within THC + ~1.1× from circuit improvements = ~233× total runtime speedup.

A literal head-to-head table is structurally impossible because the molecule sets and circuit-family targets don't overlap: Caesura is P450 at Q=116 only; GeoVac's largest composed library entry is C₂H₆ at Q=160, but its operating range is LiH/H₂O at Q=30/70. Pauli count and 1-norm/Toffoli are different metrics targeting different algorithmic regimes (VQE measurement overhead vs FT qubitization-LCU). Per CLAUDE.md §1.5 positioning: "Published DF/THC lambda values for Q<100 do not exist — GeoVac's competitive landscape is defined by raw JW baselines where it wins decisively."

The honest fix is a contextualizing paragraph, not a comparison table. Added "BLISS-THC at the production-scale ceiling" paragraph to Paper 14 §sec:ft_gaussian after "Computed 1-norm comparison." The paragraph quotes Caesura's actual numbers, names the 233× speedup breakdown, and explicitly states why a literal head-to-head can't be constructed at current library scope. Paper 14: 21 → 22 pages, three-pass clean, zero undefined refs.

### 2.4 C3 — Paper 14 missing bibitems and section labels

C1's compile audit surfaced 8 cite keys used but missing from Paper 14's inline `\thebibliography` block, plus 3 cross-paper section refs that don't resolve in per-paper compile. Strategy: add bibitems where verifiable from corpus / well-known, mark uncertain ones as `[VERIFY-REF]` placeholders for PI rather than fabricate.

**Section labels (3, all converted to explicit prose):**
- `sec:pk_partitioning` → "implemented in `geovac/pk_partitioning.py` per CLAUDE.md §7"
- `sec:frozen_core` → "implemented in `geovac/neon_core.py` (extended to [Ar] in v2.2.0; equivalently, the Z=20 entry in the atomic classifier)"
- `sec:paper22_spinor` (×2 sites) → "Paper 22 §5 'Restatement in the spinor basis'" — verified that the section actually exists in Paper 22 (line 583 of `paper_22_angular_sparsity.tex`)

**Bibitems (8, with confidence-marked status):**
- Verified additions (4): `Dyall` (Oxford 2007 Relativistic QC textbook, cross-checked against Paper 31's entry); `reiher2017` (PNAS 114:7555, well-known nitrogenase QPE paper); `Szmytkowski2007` (J. Math. Chem. 42:397, also mentioned in CLAUDE.md §12); `BreitPauli` (Bethe-Salpeter 1957 canonical textbook reference)
- Placeholder for PI (4): `BJL`, `ChildsBerry`, `MartinezYRomero2004`, `Sunaga2025` — added as bibitems with `[VERIFY-REF: PI please supply ...]` content and `% [VERIFY]` comments giving line number, context, and best-guess candidate. The PDF compile now has **zero** undefined-citation warnings (the placeholders show as text in the bibliography slot, not as `[?]` markers).

This converts a silent compile warning into a visible source-level TODO that the PI can either confirm or replace before outreach goes out.

### 2.5 C4 — CLAUDE.md §7 constants drift

A9 batch 1's first test fired `ModuleNotFoundError: No module named 'geovac.constants'` because the CLAUDE.md §7 Code Architecture table listed:

> | Physical constants | `geovac/constants.py` | `HBAR`, `C`, `ALPHA`, etc. |

The file does not exist. Survey of the actual code state:
- `C_LIGHT = 137.036` lives at `geovac/dirac_hamiltonian.py:30`
- `ALPHA = mpmath.mpf(1) / mpmath.mpf("137.035999084")` lives at `geovac/two_loop_self_energy.py:107`
- `KAPPA_SCALAR = Rational(-1, 16)` lives at `geovac/graph_qed_propagator.py:103`
- No central constants module exists; the various worktree-debug scripts that referenced `geovac.constants` were echoing the CLAUDE.md claim

Two options: (a) create `geovac/constants.py` and migrate definitions + update all callers; (b) update CLAUDE.md §7 to reflect reality. Chose (b) per CLAUDE.md §13.5 access control (PM may edit §7 for module-path updates). The §7 row now reads:

> | Physical constants | scattered next to their modules | `C_LIGHT` in `geovac/dirac_hamiltonian.py`; `ALPHA` in `geovac/two_loop_self_energy.py`; `KAPPA_SCALAR = Rational(-1,16)` in `geovac/graph_qed_propagator.py`. No central `geovac/constants.py` module. `-1/16` per CLAUDE.md §8 may be used directly. |

Option (a) is the structurally cleaner outcome but touches multiple modules' imports and is a separate refactor. The CLAUDE.md fix is the lower-disruption, immediately-correct move; a future sprint can create the central module if/when there's a forcing reason.

## 3. Net counts

| Artifact | Before | After |
|---|---|---|
| Paper 14 pages, three-pass clean | 21 (with 2 `\fbox` placeholders) | 22 |
| Paper 14 undefined citations | 10 (8 unique keys) | 0 |
| Paper 14 undefined references | 5 (3 unique labels) | 0 |
| Paper 34 §III projection spot-check coverage | 6 of 28 | **28 of 28** |
| Paper 34 spot-check test files | 1 | 4 |
| Paper 34 spot-check tests (pytest-expanded) | 12 | 108 |
| Paper 34 pages, three-pass clean | 125 | 125 (unchanged by §III.16 correction) |
| Followon-register active items (sections A + C) | 6 (A9, A10–A14, C1, C2) | 5 (A10–A14 — all are sprint-scale research follow-ons, not polish) |

## 4. Files touched

**New:**
- `benchmarks/qubit_encoding/generate_paper14_figures.py` — dedicated Paper 14 figure generator
- `papers/group4_quantum_computing/paper_14_figures/pauli_scaling.png`
- `papers/group4_quantum_computing/paper_14_figures/eri_density.png`
- `tests/test_paper34_projection_spot_checks_batch1.py` — 8 projections, 31 tests
- `tests/test_paper34_projection_spot_checks_batch2.py` — 7 projections, 22 tests
- `tests/test_paper34_projection_spot_checks_batch3.py` — 7 projections, 11 tests

**Modified:**
- `papers/group4_quantum_computing/paper_14_qubit_encoding.tex` — figures, 8 bibitems, 3 section-label conversions, BLISS-THC paragraph (lines ~282, ~326, ~1249, ~1594, ~1694, ~1789, ~2316–2350, ~2669)
- `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` — §III.16 closed-form correction (lines ~641–660)
- `tests/test_paper34_projection_spot_checks.py` — header updated to reference batches 1/2/3
- `debug/ps_1s2s_autopsy_track4_memo.md` — line 81 + correction provenance note added
- `debug/followon_register.md` — A9, C1, C2 marked closed; C3, C4 added then C4 closed
- `CLAUDE.md` — §7 Physical constants row updated to reflect actual code state

## 5. Tests

All affected suites pass. Spot check on the verification-affected modules:

```
tests/test_paper34_projection_spot_checks*.py  108 passed, 1 skipped
tests/test_breit_integrals.py                  38 passed
Total                                          146 passed, 1 skipped
```

No regressions. The Paper 34 §III.16 correction in particular is verified by both the new `test_paper34_III16_breit_R0_1s2s_1s2s` (matches `4/81` exactly) and the pre-existing `test_breit_pauli_rational_values` parametrize at the same orbital quartet — two independent test sites now agree.

## 6. Honest scope

**Closed at theorem grade:** Nothing. This sprint added no new theorems.

**Closed at structural level:** Nothing structural was added; the §III.16 correction is a numerical-fact correction, not a structural change.

**Closed at numerical-observation / verification level:**
- 108 spot-check tests covering all 28 named Paper 34 §III projections; each test verifies one stated transcendental signature, load-bearing identity, or algebraic-ring claim from the paper, against either symbolic sympy reduction, production-module output, or analytical-limit check per §13.4a protocol
- `R^0_BP(1s,2s;1s,2s) = 4/81` (pure rational) — verified bit-exact at two independent test sites and against the production module's exact Fraction/sympy arithmetic
- Paper 14 figure data points reproduce the paper's claimed α=3.15 scaling exponent to three significant figures (3.147 extracted)

**What remains as named open follow-ons:**
- Four `[VERIFY-REF]` placeholder bibitems in Paper 14 (`BJL`, `ChildsBerry`, `MartinezYRomero2004`, `Sunaga2025`) — PI must supply verified citations before outreach to the quantum-computing audience. Marked inline with line-number context and best-guess candidates.
- Reach-out queue (B1–B4 in followon register) is unchanged; this sprint did not touch it per user scope ("not reach-out").
- Substantive research follow-ons A10–A14 (Paper 40 log-power ansatz, Paper 42 Pythagorean formal proof, Paper 50 Maxwell extension, Paper 53 Q1–Q4 cluster, gravity-arc Q1–Q3) are unchanged — these are research items, not polish gaps, and the followon register correctly carries them across this sprint.

**What was NOT verified in this sprint:** The 4 placeholder bibitems' actual identities; the substance of Caesura et al. 2025 beyond the Table II/III/IV numbers extracted via pdftotext (no independent verification of their numerical claims; relayed as quoted with proper attribution); the constants-module decision (option (b) chosen; option (a) deferred).

**Verification-protocol observation worth recording.** The Paper 34 §III.16 correction is a clean example of the verification protocol working as intended (CLAUDE.md §13.4a, §13.4 consistency gate). The original closed form had been on paper since the Sprint Calc-Ps-1S2S Track 4 autopsy in May 2026 — a month — without anyone noticing the discrepancy with the production module. The discrepancy surfaced the moment the spot-check test pattern was applied. This is exactly what the discipline `feedback_audit_numerical_claims` is for: when a paper says "X = closed form", run X through the production code path and check. The cost of running the test (5 minutes) is small; the cost of an outreach email containing a wrong closed form would have been very high.

## 7. Ready for /release

- CHANGELOG.md entry: ready to add as `## [v3.51.0]`
- CLAUDE.md §1 version bump: v3.50.0 → v3.51.0
- CLAUDE.md §2 one-liner: drafted (≤30 words)
- All affected papers three-pass clean compile verified
- 146 affected tests pass, no regressions
- No hard-prohibitions violated (no fitted parameters, no geometry-hierarchy changes, no negative-result deletions, no combination-rule label change in Paper 2)
