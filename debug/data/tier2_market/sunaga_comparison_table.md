# GeoVac (Tier 2) vs Sunaga et al. 2025 -- Head-to-head resource comparison

Source: Sunaga et al., PRA 111, 022817 (2025); arXiv:2406.04992.

GeoVac: relativistic composed pipeline, n_max=2, Dirac (kappa, m_j) basis.

Sunaga: VQE-UCCSD, cv2z Dyall basis, DC Hamiltonian -> JW, 3 occ + 15 unocc spinorbitals = 18q active space.


## Part A: GeoVac raw resource metrics (reproduced from T3 memo)

| Molecule | Q | N_pauli (rel) | N_pauli (scalar) | rel/scalar | lambda_ni (rel, Ha) | lambda_ni (scalar, Ha) | QWC (rel) | QWC (scalar) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| LiH | 30 | 805 | 333 | 2.42x | 35.90 | 37.23 | 55 | 21 |
| BeH | 30 | 805 | 333 | 2.42x | 141.32 | 139.12 | 52 | 21 |
| CaH | 20 | 534 | 222 | 2.41x | 13.87 | 16.60 | 52 | 21 |

## Part B: Head-to-head with Sunaga 2025

Sunaga publishes explicit per-molecule Pauli counts only for **RaH-18q**.  We therefore compare each GeoVac molecule against the RaH-18q baseline as the single calibrated Sunaga number.  Q mismatch is flagged.

| Molecule | GeoVac Q | GeoVac N_pauli | Sunaga target | Sunaga Q | Sunaga N_pauli | GeoVac/Sunaga | Flag |
|:---|:---:|:---:|:---|:---:|:---:|:---:|:---|
| LiH | 30 | 805 | RaH-18q | 18 | 47099 | 0.017x | Q mismatch |
| BeH | 30 | 805 | RaH-18q | 18 | 47099 | 0.017x | Q mismatch |
| CaH | 20 | 534 | RaH-18q | 18 | 47099 | 0.011x | Q mismatch |

## Part C: 1-norm vs Sunaga (DEFERRED)

Sunaga does not report per-molecule 1-norms for BeH/MgH/CaH/SrH/BaH/RaH in the published figures/tables (Section III/IV).  Only ground-state energies (-3178.63 Ha for SrH STO-6G) and PDM values are reported.  Direct 1-norm comparison is a **fetchable-via-direct-paper-read** deferred item.


## Part D: QWC groups vs Sunaga (DEFERRED)

Sunaga reports **1 dominant QWC clique** for the 6-qubit PDM operator (19 Pauli terms) and for the 12-qubit PDM (53 terms).  QWC groups for the full Hamiltonian are not reported.  Comparison deferred.


## Part E: 2-qubit gate count (UCCSD/VQE ansatz)

GeoVac Tier 2 deliverable is the Hamiltonian, not a VQE ansatz, so this cell is not directly comparable.  For reference, Sunaga 12q SrH ES-VQE starts at 9,148 2qg and compresses to 13 after RL-ZX+Cflow pipeline.  The GeoVac 30-qubit relativistic Hamiltonian would need an analogous UCCSD ansatz + pipeline-optimisation pass for head-to-head gate comparison, which is out of scope for T4.


## Part F: Matched-Z family comparison

The only Sunaga molecule in scope for GeoVac (Z <= 36) is CaH (Z=20).  For that molecule only:

- GeoVac CaH Q=20, N_pauli=534, lambda_ni=13.87 Ha, QWC=52.

- Sunaga CaH: 18q active space, but no per-molecule Pauli count published.  The Sunaga 18q baselines are all calibrated on RaH (47,099 Pauli).

- GeoVac CaH at Q=20 is 1/18 the qubit count advantage over Sunaga.  GeoVac composed architecture has structural asymmetry with frozen-core reduction (Ca [Ar] core is inert) that Sunaga's cv2z Dyall basis does not have.


## Deferred items (fetchable-via-direct-paper-read)

- Per-molecule Pauli counts for BeH/MgH/CaH/SrH/BaH at 18q
- Per-molecule 1-norms (lambda) at 18q or 12q
- Per-molecule QWC group counts

Retrieving the Sunaga Supplemental Material (Tables S1-S3) would populate the per-molecule Pauli/1-norm/QWC cells for all six molecules.

