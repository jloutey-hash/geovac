# Prolate (Level 2) Pauli count: first measurement (2026-09-26)

Analysis #5 of `debug/sprint_product_basis_tradeoff_scoping_memo.md` (the G7 gap). Driver:
`debug/prolate_pauli_count.py`. Data: `debug/data/prolate_pauli_count.json`, logs
`debug/data/prolate_pauli_count_*.log`, H2 integrals `debug/data/prolate_pauli_h2_M*_ints.npz`.

## Method

- **Prolate integrals.** H2 at R = 1.4011: Paper 12 product functions used as one-electron
  orbitals (`prolate_allelectron_c4.h2_orbitals_c4`, the validated G-H2PI engine): σ4, σ6,
  and σ6 plus 1–4 π shells, so M = 4, 6, 8, 10, 12, 14. LiH at R = 3.015: the cached Route C
  10-orbital set (`prolate_separation_constant_stageC_ints.npz`, chemist (pq|rs), definite m).
  **It reproduces E_FCI = −8.007942 (log −8.00794).** Orthogonalization and every orbital
  rotation are done block-wise in signed m, so each MO keeps a definite m. Three orbital
  choices: bare-h₁ MOs, RHF MOs, FCI natural orbitals. A "real-π" variant rotates each
  π± pair to cos φ / sin φ.
- **Qubit operator.** `build_fermion_op_from_integrals` + openfermion `jordan_wigner`.
  N_Pauli and λ are **non-identity**, the corpus convention. The registry aliases 838 to
  "including identity" for composed LiH; this run measures 837 non-identity and 838 with
  identity. Identity-included values are in the JSON. QWC uses
  `measurement_grouping.count_qwc_groups`. ERI density is the fraction of the M⁴ MO-basis
  (pq|rs) with magnitude > 1e-12.
- **Energy check.** Every JW Pauli string is applied to the (N, S_z = 0) computational-basis
  sector and the result diagonalized. The full 2^Q space is never used. **Every row
  reproduces its FCI energy to ≤ 1e-8 Ha.**
- **Gaussian baselines.** PySCF has no wheel for this Windows / Python 3.14 build, so the
  driver carries its own McMurchie–Davidson engine, with basis data from
  `basis_set_exchange` (installed to the scratchpad). On openfermion's shipped PySCF data
  (H2 STO-3G, H2 6-31G, LiH STO-3G) it matches RHF and FCI to ≤ 1e-8, and N_Pauli, λ and
  QWC exactly. It also reproduces the literature H2/cc-pVDZ FCI (−1.16342) and LiH/cc-pVDZ
  RHF (−7.98373). Degenerate π/δ pairs are cleaned to real px/py orbitals. A complex-m
  variant, (px ± i py)/√2, is also reported.

## Results

H2 (R = 1.4011 bohr; exact E = −1.17447 Ha). Prolate rows use bare-h₁ MOs. RHF and NO rows
have identical N_Pauli and ERI density, with λ within ±10%; see JSON.

| basis | Q | N_Pauli | λ (Ha) | QWC | ERI dens. | E (Ha) |
|:--|--:|--:|--:|--:|--:|--:|
| STO-3G | 4 | 14 | 1.89 | 5 | 0.500 | −1.13727 |
| **prolate σ4** | 8 | **192** | 14.47 | 51 | 0.531 | −1.13087 |
| 6-31G | 8 | 184 | 11.46 | 69 | 0.500 | −1.15168 |
| cc-pVDZ, lowest 4 MOs | 8 | 184 | 8.68 | 69 | 0.500 | −1.14954 |
| **prolate σ6** | 12 | **922** | 43.01 | 312 | 0.506 | −1.16036 |
| cc-pVDZ, lowest 6 (real / complex-m π) | 12 | 326 / 522 | 25.99 / 26.58 | 77 / 124 | 0.179 | −1.15907 |
| **prolate σ6+1π (complex / real π)** | 16 | **1740 / 1256** | 66.1 / 65.1 | 611 / 382 | 0.220 | −1.16941 |
| cc-pVDZ, lowest 8 (real / complex) | 16 | 1004 / 1572 | 56.6 / 58.6 | 353 / 626 | 0.174 / 0.211 | −1.16144 |
| **prolate σ6+2π (complex / real π)** | 20 | **4710 / 2442** | 104.7 / 101.0 | 1746 / 788 | 0.171 | −1.17091 |
| cc-pVDZ full (real / complex) | 20 | 2238 / 4362 | 101.3 / 105.2 | 840 / 1946 | 0.158 | −1.16342 |
| 6-31G** full (real) | 20 | 2238 | 124.8 | 851 | 0.158 | −1.16515 |
| prolate σ6+3π (complex / real) | 24 | 9240 / 4200 | 150.6 / 143.3 | 3978 / 1569 | 0.138 | −1.17151 |
| prolate σ6+4π (complex / real) | 28 | 17082 / 7562 | 209.6 / 197.5 | 6656 / 2683 | 0.131 | −1.17180 |
| GeoVac library H2, max_n = 2 | 10 | 279 | 6.40 | 64 | 0.171 | +0.19906* |
| GeoVac library H2, max_n = 3 | 28 | 14046 | 37.79 | 5552 | 0.099 | +0.19812* |

\*The library "H2" is a single-centre Z_eff = 1 bond-pair block plus 1/R. It does not bind,
so it serves as a resource-count reference only.

LiH (R = 3.015 bohr; exact E ≈ −8.0705 Ha; all-electron, 4e).

| basis | Q | N_Pauli | λ (Ha) | QWC | ERI dens. | E (Ha) |
|:--|--:|--:|--:|--:|--:|--:|
| STO-3G | 12 | 630 | 12.34 | 179 | 0.352 | −7.88240 |
| **prolate Route C, bare-h₁ MOs (complex π)** | 20 | **8370** | 61.90 | 2810 | 0.487 | **−8.00794** |
| prolate Route C, RHF MOs | 20 | 8370 | 53.94 | 2740 | 0.487 | −8.00794 |
| prolate Route C, FCI NOs | 20 | 8370 | 55.11 | 2752 | 0.487 | −8.00794 |
| prolate Route C, real π | 20 | 6830 | 61.13 | 2101 | 0.487 | −8.00794 |
| 6-31G, lowest 10 RHF MOs (real / complex π) | 20 | 4426 / 8550 | 30.22 / 32.32 | 1428 / 3115 | 0.315 / 0.312 | −7.98873 |
| cc-pVDZ, lowest 10 (real) | 20 | 4426 | 28.39 | 1420 | 0.315 | −7.98822 |
| 6-31G full | 22 | 6861 | 41.98 | 2168 | 0.333 | −7.99828 |
| GeoVac composed LiH, max_n = 2 | 30 | 837 | 27.71 | 64 | 0.006 | (−14.476)† |
| cc-pVDZ, lowest 15 | 30 | 17565 | 80.76 | 6067 | 0.243 | −7.99430 |
| cc-pVDZ full | 38 | 47041 | 175.56 | n/c | 0.251 | n/c |

†The exported composed operator excludes PK (`to_openfermion`), so its sector ground state
is not a LiH energy. Q38 QWC and energy were not computed (cost).

## What the numbers say

1. **The scoping memo's prediction holds. Prolate orbitals carry no ERI sparsity beyond
   symmetry.** The σ-only sets have no m selection. Their ERI density (0.53, 0.51) matches the
   g/u-parity floor of a dense Gaussian set (0.50), and the N_Pauli ratio is 192/184 = 1.04 at
   Q8. The only extra zeros come from M_L (π shells) and parity, exactly as in any
   axially-symmetric Gaussian basis. There is no Neumann or product-structure sparsity.
2. **At matched Q and matched convention, prolate is slightly worse than, or level with, a
   Gaussian basis.**
   - H2, Q20: 4710 vs 4362 (complex, +8%) and 2442 vs 2238 (real, +9%).
   - H2, Q16: +11% (complex), +25% (real).
   - LiH, Q20, complex-m: 8370 vs 8550 (−2%, parity).
   - LiH, Q20, real: 6830 vs 4426 (+54%). The Route C set has one π pair among 10 orbitals,
     whereas the lowest-10 Gaussian MOs contain two. The difference is set composition, not
     prolate structure.
   - H2 at Q12 (σ6, 922) vs cc-pVDZ lowest 6 (326, which includes a π pair): 2.8× worse,
     for the same reason.
3. **The orbital convention matters as much as the basis.** Rotating π± to real cos/sin
   halves the count (H2 Q20: 4710 → 2442; Gaussian 4362 → 2238). The composed library uses
   complex m, so any cross-basis comparison must fix the convention. The real form is the
   cheaper one for both.
4. **λ: prolate is up to 2× worse for LiH and level for H2 at Q20.**
   - LiH, Q20: 54–62 Ha vs 28–32 Ha (Gaussian).
   - H2, Q20: 101–105 Ha, level with cc-pVDZ.
   - H2 at small Q: 14.5 vs 8.7–11.5 Ha (Q8); 43 vs 26 Ha (Q12).

   A likely source is near-linear dependence in the prolate set (cond(S) 6.5e4 for LiH),
   which yields high-energy virtuals. RHF MOs cut LiH λ by 13%; NOs do not help further.
5. **Energy at matched Q is where prolate wins.**
   - H2, Q20: −1.17091 vs −1.16342 (cc-pVDZ), 7.5 mHa closer to exact, i.e. within 3.6 mHa
     of exact vs 11.1 mHa.
   - H2, Q16: 8.0 mHa better.
   - LiH, Q20: −8.00794 vs −7.98873 / −7.98822, 19–20 mHa better. That also beats 6-31G
     full at Q22 (−7.99828) and cc-pVDZ lowest-15 at Q30 (−7.99430).
   - At Q8 prolate σ4 is 21 mHa **worse** than 6-31G. The advantage appears once π
     functions enter.

   So per qubit, prolate buys accuracy, not sparsity.
6. **Against GeoVac composed, prolate loses the resource contest outright.** Composed LiH
   (Q30) has 837 terms, λ 27.7 and 64 QWC groups, against 8370 / 61.9 / 2810 for prolate
   at Q20: 10× the terms, 2.2× the λ and 44× the QWC groups. Composed's ERI density is 0.006,
   against prolate's 0.49. Composed pays for that sparsity with the PK and accuracy
   bottleneck (§1.5).

## Verdict

**Prolate loses on resources, as predicted.** A Level-2 prolate qubit Hamiltonian is dense
within M_L. At matched Q and orbital convention its N_Pauli is within −2% to +25% of a
same-size Gaussian set, and its λ is equal (H2) or about 2× worse (LiH). Composed GeoVac is
an order of magnitude sparser. There is **no G7 advantage.** The one positive is accuracy per
qubit: at Q20 it beats the Gaussian bases by 7.5 mHa (H2) and 19 mHa (LiH). For Paper 14
this is a benchmarked negative: Level-2 geometry buys energy per qubit, not sparsity.

## Caveats

- The Gaussian baselines use my McMurchie–Davidson engine, not PySCF. It is validated
  against PySCF-generated openfermion data for s and p functions (bit-level agreement in
  N_Pauli, λ and QWC) and against literature cc-pVDZ energies. d functions (Li cc-pVDZ)
  are checked only through the RHF energy.
- Truncating Gaussian sets to the lowest RHF MOs is a choice. Frozen-core truncation was not
  used, because the prolate LiH set is all-electron.
- Trenev et al.'s 276 / 5851 / 63519 include a 2-qubit reduction, so they are not
  comparable to the unreduced 630 / 6861 / 47041 here. Their cc-pVDZ count is nonetheless
  larger, which suggests a different geometry or convention. Not resolved.
- The prolate orbital sets are the engine's existing hand-chosen sets, not optimized for Q.
