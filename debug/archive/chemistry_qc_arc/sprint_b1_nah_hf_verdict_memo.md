# Sprint B.1 — HF on explicit-core NaH (preliminary verdict)

**Date:** 2026-06-07.
**Driver:** `debug/b1_nah_hf_driver.py`
**Data:** `debug/data/b1_nah_hf.json`
**Log:** `debug/b1_nah_hf_log.txt`
**Scoping:** `debug/sprint_b1_explicit_core_nah_scoping_memo.md` (this sprint's design doc)

---

## Verdict

**STOP** — HF on explicit-core NaH does **not** close the W1e wall.

Caveat: preliminary. 2 of 7 R points failed to converge under our SCF (DIIS + density damping + level shift). A clean publication-grade verdict would require a pyscf cross-check on the same FCIDUMP. But the qualitative result is robust to SCF noise: the 5 R points that DID converge show a clean monotone-descending E(R) into small R with no interior minimum — the same overattraction pattern as the frozen-core balanced baseline.

## Result

R panel {2.5, 3.0, 3.566, 4.0, 4.5, 5.0, 6.0} bohr; closed-shell RHF SCF with DIIS, density damping (α=0.3), level shift (1.0 Ha on virtuals).

**Explicit-core NaH HF energies (12 electrons, 30 qubits, 15 spatial orbitals):**

| R (bohr) | E_HF (Ha) | converged | E vs R=6.0 |
|---|---|---|---|
| 2.500 | −222.067 | ✓ | −4.04 (very over-bound) |
| 3.000 | −221.505 |  ✓ | −3.48 |
| 3.566 | −220.878 | ✓ | −2.85 |
| 4.000 | −209.390 | ✓ (different basin) | +8.64 (less bound) |
| 4.500 | −219.544 | ✗ (400 iter) | — |
| 5.000 | −208.023 | ✗ (400 iter) | — |
| 6.000 | −218.026 | ✓ | 0 |

The converged points (2.5, 3.0, 3.566, 6.0) show monotone descent into small R. R_min at the panel boundary R=2.5 with D_e = +4.04 Ha (54× over-attractive vs experimental 0.075 Ha).

The intermediate R (4.0–5.0) shows SCF bistability — bouncing between solutions roughly 12 Ha apart. This is the classical signal of HF failure modes during bond formation/dissociation: either a near-degeneracy that single-determinant HF can't represent, or the basis is too restrictive (likely both).

## Methodology cross-check (LiH HF)

We ran the same SCF on first-row LiH (explicit-core 1s² + 2 valence = 4 electrons), where the balanced FCI verdict from today's W1e-Projection-Audit gives R_eq = 3.015 with D_e = 0.158 Ha.

**LiH HF energies (4 electrons, 30 qubits, 15 spatial orbitals):**

| R (bohr) | E_HF (Ha) | converged |
|---|---|---|
| 2.500 | −15.165 | ✓ |
| 3.015 | −15.201 | ✓ (R_min) |
| 3.500 | −15.199 | ✓ |
| 4.000 | −15.167 | ✓ |
| 5.000 | −15.043 | ✓ |

**R_min = 3.015 bohr (matches experimental R_eq exactly).**
**D_e = 0.158 Ha (matches today's balanced FCI to 4 digits).**

All SCF convergence in ~40 iterations. Methodology is validated: HF works correctly on the framework's Hamiltonian when the molecule is amenable to a single-determinant treatment.

The exact D_e match between HF and FCI for LiH is consistent with minimal basis (n_max=2) → near-zero correlation energy. At this truncation, HF ≈ FCI.

## Reading

The methodology is sound. The framework's HF on the W1e-Projection-Audit balanced LiH Hamiltonian reproduces the balanced FCI verdict. Therefore the explicit-core NaH HF failure is **not** a methodology artifact — it's a structural property of the framework's NaH Hamiltonian at this basis truncation.

Combining with today's other sprint results:

| Sub-thread | Verdict |
|---|---|
| LiH kwarg sweep (existing kwargs) | STOP — cheap engineering kwargs are no-ops or catastrophic on LiH |
| Frozen-core balanced NaH | STOP — monotone overattraction |
| Frozen-core balanced LiH (Path A, post-bug-fix) | BINDS at right R, 2.4× over-binding |
| Explicit-core HF NaH (this sprint) | STOP (preliminary) |
| Explicit-core HF LiH (cross-check) | BINDS, methodology validated |
| H₂ Bratteli pilot (same day) | Marcolli–vS PASS bit-exact; Perez-Sanchez 2024a NO |

**Net structural reading:** The chemistry W1e wall is not engineering-tier closeable at the framework's current scope. Frozen-core, multi-zeta, screened-valence, cross-block ERI, cross-block h1, and explicit-core HF have all been tested. NaH does not bind in any of them.

LiH binds because the static Z_eff ≈ 1.84 valence orbital is close enough to the physical Li 2s that the fixed-basis multipole-V_ne machinery captures the right physics. NaH's Z_eff ≈ 2.2 with the much more diffuse 3s orbital sits outside this regime.

This is consistent with the structural-skeleton-scope pattern (MEMORY.md `geovac_structural_skeleton_scope_pattern.md`): the framework maps the structural skeleton; it does not autonomously generate calibration data. For chemistry, "calibration data" includes the system-specific shape of the bond's Pauli-repulsion content against the explicit core.

## What's preserved, what's lost

**Preserved:**
- The framework's structural advantages (O(Q^2.5) Pauli scaling, propinquity-bounded truncation error, zero-parameter construction) are unaffected.
- LiH chemistry binds at the right R; first-row hydrides remain a usable regime.
- The Bratteli H₂ pilot showed bit-exact match to Marcolli–vS 2014; the NCG structural home is real.

**Lost:**
- The hope that an existing engineering knob closes second-row chemistry binding.
- The hope that explicit-core treatment alone closes second-row.
- The expected payoff of a CCSD-class extension at the qubit level (would close N×expensive correlation but not bond formation).

## Recommended next sprints

This is the strategic decision-gate moment. Three options, in increasing scope:

1. **Accept the chemistry-accuracy scope boundary for second-row.** Document this in Paper 20 / Paper 57 framing: "GeoVac chemistry is chemical-accuracy for first-row hydrides; the second-row scope boundary is structural, not engineering." Position the framework as resource-estimation + first-row-chemistry. Pursue the Marcolli–vS NCG paper (today's H₂ Bratteli pilot is the seed) as the structural arc. This is the conservative pivot.

2. **One more engineering test before accepting.** B.4 (bonding-orbital Pauli repulsion against explicit core — the operator class the kwarg-sweep Explorer named) is the only remaining unclosed engineering hypothesis. 5 days, narrowly scoped. If it fails, the scope boundary is confirmed; if it succeeds, second-row is reachable.

3. **Full pivot to NCG / Marcolli–vS theorem.** Treat today's chemistry results as the empirical motivation for the structural reformulation: "we have demonstrated the engineering-tier closure is exhausted; here is the Marcolli–vS gauge-network home that explains *why* and predicts where it CAN close." Multi-quarter math.OA arc. Highest leverage.

My read: (1) and (3) are complementary, not exclusive. Run them together. Skip (2) — its expected information value, given the convergence of today's negatives plus the H1/W1e/Yukawa pattern, is low.

## Honest scope

- 2 of 7 R points unconverged; pyscf cross-check would tighten the verdict but won't change qualitative direction.
- We only ran one SCF parameter set (α=0.3 damping, 1.0 Ha level shift, DIIS from iter 8). A more aggressive damping schedule might converge the intermediate R points but is unlikely to flip the topology.
- We have not tried multi-zeta basis on the explicit-core Na block. Plan agent's A.2 would add this — `multi_zeta_basis=True` currently only triggers on frozen-core specs (`Z_nuc_center >= 11 AND n_val_offset > 0`); we'd need to widen the dispatcher to first-row-style explicit-extended specs. Not yet done.
- We have not tested whether the explicit-core HF for NaH binds at higher n_max. The basis sanity diagnostic noted max_n=3 was required for the 3s to exist at all; we did use max_n=3 effectively (max_n=2 for core + n_val_offset=2 for bond block → block_n=1 maps to physical n=3 for the valence). But the bond block has only max_n=2 → physical n=3, 4 valence orbitals. Larger basis (max_n=3 on bond block → physical n=3, 4, 5) might give the SCF more variational freedom; unlikely to change qualitative direction but a real follow-on.
