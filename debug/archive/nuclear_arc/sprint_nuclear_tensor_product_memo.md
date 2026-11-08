# Sprint: Nuclear Model + Tensor-Product Spectral Action
**Date:** 2026-06-01
**Status:** CLOSED

## 1. Goal

Two-part sprint: (A) push the nuclear model to production grade (binding, observables, tensor force), (B) attack the inter-particle interaction problem head-on via the tensor-product spectral action.

## 2. Deliverables

### A. Nuclear model improvements

**Minnesota V_S/V_T swap (BUG FIX).** The singlet (V_S = -178 MeV) and triplet (V_T = -91.4 MeV) Minnesota potential parameters in `geovac/nuclear/minnesota.py` were swapped. The ground state was ¹S₀ (J=0) instead of the correct ³S₁ (J=1). Fixed: V_T = -178 MeV (triplet, stronger), V_S = -91.4 MeV (singlet, weaker). Eigenvalue spectrum unchanged (values swap channels); all Pauli counts, 1-norms, and resource estimates invariant. Paper 23 parameter ordering corrected.

**Sturmian (exponential) basis for deuteron.** Replaced the HO basis (Gaussian tails, B_d ≈ 0 at hw=8) with a Coulomb-Sturmian-like basis with exponential tails exp(-αr/2). At n_basis=16, α=1.1: B_d = 2.197 MeV (-1.2% vs exp 2.225), r_d = 2.118 fm (-0.5%). Conditioning < 10² throughout. The HO basis at N_shells=2 gives B_d ≈ 0 for the relative motion — the Sturmian basis is categorically superior for binding.

**Tensor force (coupled ³S₁ + ³D₁).** Added a Gaussian tensor V_T0·exp(-κ_T·r²)·S₁₂ to Minnesota. S₁₂ matrix elements corrected to standard values: ⟨S|S₁₂|S⟩=0, ⟨D|S₁₂|S⟩=√8, ⟨D|S₁₂|D⟩=−2. Best fit at κ_T=0.7: V_T0=−16.1 MeV, Q_d=0.286 fm² (exact), B_d=2.395 MeV (+7.7%), P_D=0.37%, μ_d=0.878 n.m. (+2.4%). The B_d/Q_d/P_D tension is physical — Gaussian tensor cannot simultaneously fit all three (requires OPEP 1/r³ core for P_D ~ 5%).

**Nuclear structure observables (M_J-projected).** Complete catalogue at hw=8 with M_J=+1 projection: r_d=2.119 fm (−0.4%), μ_d=0.880 n.m. (+2.6%, S-wave limit), Q_d=0 (no tensor), α_E=0.694 fm³ (+9.7%). He-4 r_pp=1.527 fm (+4.5%) at hw=20. Zemach radius (point-nucleon) r_Z=2.569 fm (+13.7%), Friar moment 26.6 fm³ (+51%) — Gaussian tail artifact, not basis size.

**LIT method (diagnostic negative).** Lorentz Integral Transform validated (bit-exact vs eigendecomposition) but unhelpful for finite-basis effective interaction. The +9.7% α_E error is from Minnesota (EWSR=152% > TRK, too much low-energy strength), not from the method. N_shells=2→3 gives BIT-IDENTICAL α_E (Minnesota is basis-converged at N_shells=2). Tikhonov inversion adds positive bias.

**NaH Z_orb scan (diagnostic negative for chemistry bridge).** Scanned Z_orb ∈ {0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0} for the NaH bond block at R_eq=3.566. PES monotonically descending at every Z_orb. Descent magnitude similar (0.32–0.36 Ha). The NaH overattraction is structural to the balanced-coupled architecture (W1e multi-determinant), not basis extent. The nuclear Sturmian improvement does NOT transfer to the chemistry solver.

### B. Tensor-product spectral action (HEADLINE)

**Step 1 — Free tensor-product Dirac.** Built D_total = D₁⊗I + γ₁⊗D₂ on the chirality-doubled scalar S³ basis at n_max=2,3. Key finding: {D, γ} = 0 exactly for our construction, so D²_total = D²₁⊗I + I⊗D²₂ with ZERO cross term. Heat trace factorizes to machine precision. Physically correct: free particles don't interact.

**Step 2 — Inner fluctuations.** Computed [D_total, a₁⊗a₂] for 3-Y multiplier algebra elements. Both terms nonzero (‖term1‖ = ‖term2‖). The algebra generates cross-particle coupling.

**Step 3 — Gauged spectral action.** Built A = (a₁⊗a₂)·[D_total, b₁⊗b₂] from the strongest commutator pair. Computed {D, A} and decomposed into factorizable + connected parts via partial-trace subtraction. **Connected fraction = 75.7% (n_max=2), 74.3% (n_max=3).** Stable across n_max — structural, not finite-size. Heat trace departs from factorized by 0.28% (n_max=2) to 0.66% (n_max=3) at ε=1.

**Step 4 — Angular decomposition.** Extracted the (l₁,l₂)→(l₁',l₂') channel structure of the connected part in the (+,+) chirality sector.

| Property | n_max=2 | n_max=3 |
|:---------|--------:|--------:|
| m-conservation | 100% | 100% |
| Gaunt-compatible | 100% | 65.3% |
| k=0 monopole | 100% | 96.5% |
| k=2 quadrupole | 0% | 3.5% |
| Hermiticity | exact | exact |
| Exchange symmetry | exact | 10⁻¹⁷ |

The monopole + quadrupole hierarchy matches the Coulomb multipole expansion. The k=1 dipole is absent (correct for identical particles by parity). The 34.7% Gaunt-incompatible at n_max=3 comes from using a single-multiplier gauge field; summing over all algebra generators (the full rotationally invariant gauge field) should recover full Gaunt compatibility.

## 3. Files Created

- `debug/deuteron_lit_polarizability.py` — LIT method implementation
- `debug/deuteron_lit_diagnostic.py` — LIT diagnostic (EWSR, N_shells convergence)
- `debug/nuclear_structure_observables.py` — v1 observables (M_J bug)
- `debug/nuclear_structure_observables_v2.py` — v2 M_J-projected observables
- `debug/deuteron_sturmian_nuclear.py` — Sturmian basis deuteron solver
- `debug/deuteron_tensor_force.py` — Coupled ³S₁+³D₁ with Gaussian tensor
- `debug/tensor_product_dirac.py` — Free tensor-product Dirac diagnostic
- `debug/tensor_product_gauged.py` — Gauged spectral action + factorization test
- `debug/tensor_product_decompose.py` — Angular channel decomposition
- `debug/data/deuteron_lit_polarizability.json`
- `debug/data/deuteron_lit_diagnostic.json`
- `debug/data/nuclear_structure_observables.json`
- `debug/data/nuclear_structure_observables_v2.json`
- `debug/data/deuteron_sturmian_nuclear.json`
- `debug/data/deuteron_tensor_force.json`

## 4. Files Modified

- `geovac/nuclear/minnesota.py` — V_S/V_T swap fix (lines 53-59)
- `papers/group4_quantum_computing/paper_23_nuclear_shell.tex` — parameter ordering (line 322)

## 5. Decisions

- LIT does not improve over SOS for finite-basis effective interaction — CLOSED as diagnostic negative
- NaH overattraction is W1e (multi-determinant), not basis extent — chemistry bridge hypothesis RULED OUT
- The gauged tensor-product spectral action IS the head-on path for deriving inter-particle interactions
- Paper 23 numerical results invariant under V_S/V_T fix — only quantum number labels corrected

## 6. Honest Scope

**Closed at production grade:**
- Minnesota V_S/V_T fix (bit-exact eigenvalue invariance verified)
- Paper 23 parameter ordering correction (compiles clean, 12 pages)
- M_J-projected deuteron observables (r_d −0.4%, μ_d +2.6%)
- NaH Z_orb scan (9 values, all monotonically descending)

**Numerical observation (not theorem-grade):**
- Sturmian B_d = 2.197 MeV at n=16, α=1.1 (−1.2%)
- Tensor force Q_d = 0.286 fm² at V_T0=−16.1 MeV
- Connected fraction 75% of {D, A} at n_max=2,3
- k=0 monopole 96.5–100%, k=2 quadrupole 3.5% at n_max=3
- 100% m-conservation, 65–100% Gaunt compatibility

**Structural sketch (needs theorem-grade proof):**
- {D, γ} = 0 for the chirality-doubled scalar construction (verified numerically, obvious analytically)
- D²_total factorizes (direct consequence of {D, γ} = 0)
- The connected fraction measures genuine two-body content (partial-trace subtraction is standard)

**Named open follow-ons:**
1. Sum over ALL multiplier pairs in the gauge field — should recover 100% Gaunt compatibility
2. Extract radial weights (k=0 vs k=2 ratio) and compare with Paper 12 Slater integrals
3. OPEP tensor (1/r³ core) for proper P_D ~ 5% and μ_d correction
4. Write up the tensor-product spectral action as a paper (potentially Paper 54)
5. Extend to SU(2) Wilson gauge (Paper 30) for the nuclear force channel
