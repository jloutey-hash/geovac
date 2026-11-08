# Sprint memo — xTC PoC on Li (Coulomb-Sturmian)

**Date:** 2026-08-23  **Verdict:** _TBD (fill after study)_

**Question.** Does xTC — the 3-body TC operator L3 contracted to an effective
2-body operator via the reference 1-RDM (Christlmaier–Kats–Alavi, JCP 159 014113
(2023)) — recover accuracy at a *smaller / sparser* two-body qubit Hamiltonian than
plain FCI, on GeoVac's Coulomb-Sturmian basis? Li (1s²2s, doublet) is the minimal
system where the genuine 3-body term is live (He/H₂ have only the j=k two-body piece).

Code: `debug/xtc_poc_li.py` (engine), `debug/xtc_poc_li_validate.py` (gates),
`debug/xtc_poc_li_study.py` (study) → `debug/data/xtc_poc_li_study.json`.

## Recipe implemented
- **Basis.** Coulomb-Sturmian s-radials R_{n0}(r; nk) (shared decay k), built on one
  clustered radial grid, Löwdin-orthonormalized. Everything (S, h1, 2-body, 3-body)
  on the same grid → convention-free.
- **Geminal.** Slater u(r) = −(1/2γ)e^{−γr}, u'(0)=1/2 (singlet Kato cusp),
  u'(r)=½e^{−γr}. Single spin-independent Jastrow (caveat below).
- **TC Hamiltonian.** H̃ = H + D + K + L3.
  - D (Hermitian 2-body) ≡ the finite kernel w(r)=(1−e^{−γr})/r + (γ/2)e^{−γr} − ¼e^{−2γr}
    (matches Track 2 `ctf12_validate`).
  - K (non-Hermitian 2-body convective) = −u'(r₁₂) r̂₁₂·(∇₁−∇₂).
  - L3 (3-body) = −½ Σ_i Σ_{j≠i,k≠i,j≠k} u'(r_ij)u'(r_ik)(r̂_ij·r̂_ik).
    For s-orbitals the angular average **factorizes** into a shared-vertex product
    A(a,b,c)=g(a,b)g(a,c), g(a,b)=½∫₋₁¹ u'(r_ab)(a−bx)/r_ab dx  (MC-validated to ~1e-3;
    g ≡ the convective vertex kernel kA). This makes V3 a cheap triple radial integral.
- **xTC contraction.** L3 second-quantized as −½ Σ V3 a†a†a† aaa; fully antisymmetrized
  W; effective bare 0/1/2-body from single/double/triple γ-contraction over the reference
  occupation (aufbau 1s↑1s↓2s↑), converted normal-ordered→bare via Wick with diagonal γ.
- **Diagonalization.** Convention-free 2nd-quantized operator application, particle-number
  projected FCI. Non-Hermitian (K, xTC) → `scipy.linalg.eig`, real ground state.

## Validation gates (all PASS)
- Grid S, h1 vs analytic Sturmian: 4e-16 / 3e-11. Coulomb ERI matches to grid
  precision (the 5e-2 vs `_slater_rk` is the latter's coarse 500-pt linspace; He
  s-limit −2.87859 vs −2.879029 confirms our ERIs).
- **xTC contraction exact on ref+singles+doubles:** effective-operator ref-row matches
  the exact 3-body H3 to **1.7e-18** on all ≤2-excitation determinants; ref diagonal
  collapses to v0 exactly. Dropped residual = pure-3-body normal-ordered part (genuine).
- **geminal→0 (γ large):** TC2, exactTC, xTC all → plain FCI monotonically
  (|xTC−plain| = 4.8e-3 @ γ=4 → 4.4e-4 @ γ=10).

## Results

### (1) γ-scan (ns=3, k=1.5) — xTC vs exact 3-body
| γ | plain | TC2 | exactTC | xTC | 3-body exact / xTC (mHa) | xTC fidelity |
|---|-------|-----|---------|-----|--------------------------|--------------|
| 0.60 | −7.24237 | −7.29552 | −7.31425 | −7.31382 | −18.72 / −18.29 | 97.7% |
| 0.80 | −7.24237 | −7.28795 | −7.29627 | −7.29597 | −8.32 / −8.02 | 96.4% |
| 1.00 | −7.24237 | −7.28124 | −7.28503 | −7.28484 | −3.79 / −3.60 | 95.0% |
| 1.20 | −7.24237 | −7.27542 | −7.27718 | −7.27706 | −1.76 / −1.64 | 93.2% |
| 1.50 | −7.24237 | −7.26830 | −7.26886 | −7.26880 | −0.55 / −0.50 | 91% |

**xTC reproduces the exact 3-body TC energy to 0.1–0.4 mHa (93–98% of the genuine
3-body shift) across the whole geminal range.** The 3-body term is real and grows for
diffuse geminals (up to ~19 mHa at γ=0.6). TC (2-body + xTC 3-body) recovers 26–72 mHa
of correlation beyond plain FCI *at the same basis*. (TC is non-variational — energy
decreases with smaller γ; γ should be fixed by stationarity in production, not by min-E.)

### (2) accuracy vs basis (s-only, k optimized per ns for plain)
_TBD from JSON — E_plain(ns) vs E_xTC(ns); does xTC(ns) reach plain(ns')?_

### (3) qubit / sparsity of the effective 2-body operator
_TBD from JSON — 1-norm (LCU λ proxy) and nnz, plain vs xTC at matched basis._

## Verdict & caveats
_TBD._

Caveats already known:
- **Single common-k Coulomb-Sturmian caps absolute accuracy** for Li (can't span the
  1s core vs 2s valence scales) — affects plain and xTC equally; not an xTC failure.
- **s-only cannot exercise Track 1's angular Gaunt sparsity** (all l=0). The 1-norm
  comparison is measured; the "does the contracted 2-body inherit the ~5% angular
  density" question needs p-orbitals (the p-inclusive L3 four_Y machinery of Track 1).
- Single spin-independent geminal (u'(0)=½ for the singlet 1s pair; over-corrects the
  parallel-spin cusp). A spin-resolved Jastrow is a standard refinement.
- Classical PoC: the effective operator is non-Hermitian → a quantum algorithm needs
  QEVE/QITE (noted, not implemented).
