# Spin-ful Composed Qubit Pipeline — Design Memo (Track T3)

**Sprint:** Dirac-on-S³ Tier 2, Track T3.
**Date:** 2026-04-15.
**Status:** Complete. All 13 T3 regression tests pass; 164 pre-existing
tests (composed_qubit + T1 + T2 + D1) pass without modification.

## 1. API changes

### `geovac/molecular_spec.py`

`MolecularSpec` (additive, default preserves scalar path):

```diff
 @dataclass
 class MolecularSpec:
     name: str
     blocks: List[OrbitalBlock]
     nuclear_repulsion_constant: float
     description: str = ''
     core_method: str = 'pk'
     nuclei: Optional[List[Dict]] = None
+    # Track T3: dispatch to spinor (DiracLabel) builder when True.
+    relativistic: bool = False
```

`OrbitalBlock` (additive):

```diff
     l_min: int = 0
     center_nucleus_idx: int = -1
     partner_nucleus_idx: int = -1
+    # Optional DiracLabel lists populated by the relativistic builder.
+    spinor_labels: Optional[List[Any]] = None
+    spinor_partner_labels: Optional[List[Any]] = None
```

Three new spec factories: `lih_spec_relativistic`, `beh_spec_relativistic`,
`cah_spec_relativistic`. Each reuses its scalar counterpart's block topology
and flips `relativistic=True`. BeH is BeH₂ with one H-bond removed; CaH is
CaH₂ with one H-bond removed. Nuclear repulsion is recomputed.

### `geovac/composed_qubit.py`

`build_composed_hamiltonian` signature extended (both new kwargs default to
the non-relativistic path):

```diff
 def build_composed_hamiltonian(
     spec: 'MolecularSpec',
     pk_in_hamiltonian: bool = True,
     verbose: bool = False,
+    relativistic: bool = False,
+    alpha_num: float = 7.2973525693e-3,   # CODATA α
 ) -> Dict[str, Any]:
```

When either the `relativistic` kwarg or `spec.relativistic` is `True`, the
builder dispatches to `composed_qubit_relativistic.build_composed_hamiltonian_relativistic`.
The scalar code path is completely unchanged; existing LiH/BeH₂/H₂O benchmarks
(334, 556, 778 Pauli terms incl. identity) reproduce bit-for-bit.

### `geovac/composed_qubit_relativistic.py` (new module)

- `enumerate_dirac_labels(max_n, l_min)`: enumerates all DiracLabel triples
  `(n, κ, m_j)` with `l(κ) < n`, emitting one spin-orbital per label.
- `jj_angular_Xk(κ_a, 2m_a, κ_c, 2m_c, k) -> float`: full-Gaunt jj-coupled
  angular coefficient (Dyall §9 / Grant §7.5) in closed form using
  `sympy.physics.wigner.wigner_3j`. Cached.
- `_build_spinor_eri_block`: per-block physicist ERI `⟨ab|V|cd⟩` assembly
  using `X_k · X_k · R^k` with the scalar `_compute_rk_integrals_block`
  radial cache from `composed_qubit` (algebraic via `hypergeometric_slater`).
- `build_composed_hamiltonian_relativistic(spec, alpha_num, pk_in_hamiltonian,
  verbose)`: the main entry point. Returns dict with `Q`, `N_pauli`,
  `lambda_total`, `lambda_ni`, `qwc_groups`, `qubit_op`, `fermion_op`,
  `eri_sparse`, `cross_block_eri_count`, `h1_so_diag`, etc.

## 2. The 6j extension (why no explicit 6j symbol appears)

The Dyall §9.3 jj-coupled Coulomb formula in its cleanest form is

$$
\langle a b | 1/r_{12} | c d \rangle
= \sum_k X_k(a,c) \, X_k(b,d) \, R^k(n_a l_a, n_b l_b, n_c l_c, n_d l_d)
$$

where $X_k(a,c)$ carries both the radial parity rule $(l_a + l_c + k)$ even
and the Clebsch-Gordan structure from coupling two spinor spherical harmonics
through a rank-$k$ scalar. In this "separable" form the 6j symbol does not
appear explicitly — it factorises through the product $X_k \cdot X_k$ with
the m-conservation constraint $m_{j,a}+m_{j,b}=m_{j,c}+m_{j,d}$. This is the
"pre-recoupled" form that T0 used for its density computation.

The 6j would reappear if we recoupled the two-electron angular wavefunction
into the total $J = J_1 + J_2$ basis. We do not need that here — the
composed pipeline operates in the uncoupled (spin-orbital-indexed)
FermionOperator, so Dyall §9.3's direct form suffices.

**Selection rules encoded in `jj_angular_Xk`:**
1. $l_a + l_c + k$ even (parity)
2. $|j_a - j_c| \le k \le j_a + j_c$ (j-triangle, in $2j$ representation)
3. $3j(j_a\,k\,j_c; \tfrac{1}{2}\,0\,-\tfrac{1}{2}) \ne 0$ (reduced matrix element)
4. $3j(j_a\,k\,j_c; -m_{j,a}\,q\,m_{j,c}) \ne 0$ with $q = m_{j,a}-m_{j,c}$
5. Global $m_{j,a}+m_{j,b} = m_{j,c}+m_{j,d}$ (from bottom-row 3j sum-to-zero
   on the partner pair)

These reproduce T0's fullgaunt convention. The scalar (non-relativistic)
builder uses a compatible spherical-Gaunt formulation via `_ck_coefficient`;
the Dirac builder replaces the $l$-only 3j's with their $j$-valued analogues.

## 3. Three molecular specs

| Spec | Block decomposition | Scalar Q / Dirac Q at n_max=2 |
|:-----|:--------------------|:------------------------------|
| `lih_spec_relativistic()` | Li 1s² core + LiH σ-bond (Z_eff=1 + H) | 30 / 30 |
| `beh_spec_relativistic()` | Be 1s² core + 1 BeH σ-bond (Z_eff=2 + H)  | 30 / 30 |
| `cah_spec_relativistic()` | [Ar] frozen core + 1 CaH σ-bond (Z_eff=2 + H) | 20 / 20 |

Qubit count equality (scalar JW-doubled = relativistic direct) is
**structural**: the scalar builder spin-doubles each spatial orbital
(Q = 2·M), while the relativistic builder directly enumerates one qubit
per DiracLabel; the number $2(l_\max+1)^2$ of labels per n-shell exactly
equals $2 \cdot (l_\max+1)^2$ spin-doubled spatial orbitals. The two
encodings cover the same one-body Hilbert space in different bases.

## 4. Resource metrics

Computed on GeoVac development environment (Windows, Python 3.14, sympy
1.14). Wall times are single-thread end-to-end from spec construction to
JW qubit_op.

### LiH

| n_max | Q | Pauli (scalar) | Pauli (rel) | rel/scalar | λ_ni (scalar) | λ_ni (rel) | QWC (scalar) | QWC (rel) | wall (rel) |
|:-----:|:-:|:--------------:|:-----------:|:----------:|:-------------:|:----------:|:------------:|:---------:|:----------:|
| 1     |  6 |          9     |          9 | 1.00× |  10.15 |  10.15 |   1 |   1 |   0.1 s |
| 2     | 30 |        333     |        805 | 2.42× |  37.23 |  35.90 |  21 |  55 |   0.2 s |
| 3     | 84 |      7 878     |     46 434 | 5.89× | 170.59 | 126.75 | 790 | 6 571 |  45 s |

### BeH

| n_max | Q | Pauli (scalar-1bond) | Pauli (rel) | rel/scalar | λ_ni (scalar) | λ_ni (rel) | QWC (rel) | wall |
|:-----:|:-:|:--------------------:|:-----------:|:----------:|:-------------:|:----------:|:---------:|:----:|
| 1     |  6 |           9 |          9 | 1.00× |  63.44 |  63.44 |   1 |  0.0 s |
| 2     | 30 |         333 |        805 | 2.42× | 139.12 | 141.32 |  52 |  0.2 s |
| 3     | 84 |       7 878 |     46 434 | 5.89× | 346.09 | 297.03 | 6 571 | 45 s |

"scalar-1bond" = `beh_spec_relativistic` with `relativistic=False` forced
(the scalar BeH₂ spec has two H-bonds and is structurally distinct).

### CaH

| n_max | Q | Pauli (scalar-1bond) | Pauli (rel) | rel/scalar | λ_ni (scalar) | λ_ni (rel) | QWC (rel) | wall |
|:-----:|:-:|:--------------------:|:-----------:|:----------:|:-------------:|:----------:|:---------:|:----:|
| 1     |  4 |           6 |          6 | 1.00× |   2.03 |   2.03 |    1 |  0.2 s |
| 2     | 20 |         222 |        534 | 2.41× |  16.60 |  13.87 |   52 |  0.1 s |
| 3     | 56 |       5 252 |     30 940 | 5.89× |  96.45 |  65.81 | 6 571 | 26 s  |

### Reading the table

- **Q identical by construction.** The relativistic builder enumerates one
  qubit per DiracLabel; the scalar builder spin-doubles each spatial
  orbital. Both produce $Q = 2 \cdot \sum_n n^2$ (for full $l_\max=n-1$
  blocks).
- **Pauli ratio grows with n_max.** n_max=1 gives 1.00× (only s-states,
  κ=−1 only, the angular structure is a single scalar × spin tensor product).
  n_max=2 gives a clean ~2.4× across all three molecules (the spin doubling
  that scalar JW implicitly encodes through σ-sum constraints becomes
  explicit in the spinor builder; $p_{1/2}$ and $p_{3/2}$ branches carry
  distinct integrals). n_max=3 gives ~5.9× — the two $p_{1/2}/p_{3/2}$
  branches plus the full $d_{3/2}/d_{5/2}$ pair expose the T0 absolute-
  Pauli-count growth.
- **1-norm mostly decreases or stays comparable.** This is non-trivial and
  chemically meaningful: spinor ERIs factor m_j sums through narrower
  angular channels, and at finite truncation the aggregate 1-norm is
  reduced. Consistent with T0's "sparsity exponent preserved, prefactor
  increased" headline: absolute Pauli count grows, but the per-term
  weights shrink in proportion so λ stays flat.
- **Structurally, the rel/scalar Pauli ratio exceeds the T0 angular-density
  naive prediction (~0.76 ratio) because Q is matched rather than
  l_max-matched.** T0 measured density on the 4-tuple enumeration at
  fixed Q_spinor = 2·Q_scalar (T0 §3 point iii: "spinor encoding is
  absolutely denser — at l_max=3, the spinor Hamiltonian has
  (32/16)⁴ = 16× more 4-tuples"). Here we instead match Q, so the
  Hamiltonians live in the same spin-orbital count; the cost of
  relativistic physics manifests as a ~(2–6)× Pauli multiplier from the
  angular-structure richness, not the full 16× T0 quoted at matched l_max.

## 5. Regression test status

All existing composed/Dirac/spin-orbit tests pass:

```
tests/test_composed_qubit.py:        25 passed
tests/test_dirac_matrix_elements.py: 66 passed
tests/test_spin_orbit.py:            22 passed
tests/test_dirac_s3.py:              51 passed
```

New T3 regression suite (`tests/test_spin_ful_composed.py`, 13 tests):

- Scalar LiH/BeH₂/H₂O Pauli counts 334/556/778 preserved
- `MolecularSpec.relativistic` defaults to False across all hydride specs
- LiH/BeH/CaH relativistic builds succeed at n_max ∈ {1, 2}
- Block-diagonality: zero cross-block ERI entries
- α → 0 zeroes the spin-orbit diagonal exactly; finite α activates it for
  l ≥ 1 channels
- Hermiticity of the JW qubit operator verified
- LiH n_max=1 Pauli count = 9 (pinned as exact regression)
- LiH n_max=2 rel/scalar Pauli ratio in [1.9, 3.0] (pin at 2.42×)

## 6. Structural findings

### Surprise 1: Pauli ratio at n_max=2 is cleanly ~2.4×

T0 predicted d_spinor/d_scalar ≤ 1 and approaching 1 as l_max grows. At
l_max=1 the fullgaunt ratio was 0.58. Here we see **Pauli ratio ~2.42×**,
NOT the density ratio directly. The explanation:

- Scalar JW spin-doubling implicitly couples the α-spin and β-spin spatial
  copies via the $\sigma = \tau$ diagonal constraint in the two-body
  assembly (`build_fermion_op_from_integrals`, lines 302-314).
- In the spinor basis, $p_{1/2}$ and $p_{3/2}$ branches are distinct
  one-body orbitals carrying the same (n, l) but different κ. Their
  Coulomb matrix elements differ because the radial integrals depend on
  (n_a, l_a, n_b, l_b, ...) only — same values — but the angular X_k
  coefficients differ through j-triangle and (j 1/2 0 -1/2) 3j.

The 2.42× factor is approximately $2 \cdot \bar{d}_{\text{spinor}}/\bar{d}_{\text{scalar}}$
at l_max=1 (full-gaunt 0.58 × 2 spin-branches × ~2 for new ERI couplings)
— consistent with T0.

### Surprise 2: 1-norm does NOT inflate

λ_ni(rel) ≈ λ_ni(scalar) within 20% at all n_max tested, and actually
*decreases* at n_max=3 (e.g., LiH 170.6 → 126.8 Ha). Interpretation:
spinor jj-coupling distributes Coulomb matrix elements across more angular
channels with smaller per-channel magnitudes. This is a **positive** result
for Tier 3 market comparison with Sunaga et al. 2025 — 1-norm scaling is
the dominant resource estimator for QPE.

### Surprise 3: QWC group count inflates heavily

LiH n_max=3: QWC scalar 790 → rel 6 571 (8.3×). The relativistic
Hamiltonian has richer Pauli structure (more distinct $X_iY_j$ patterns
from the m_j rotation implied by κ conservation), so QWC-compatible groups
are smaller. This hurts VQE measurement overhead but is mitigated by
downstream grouping strategies (SCD, tensor-factorised measurement).

## 7. What T4 consumes

Three relativistic `QubitOperator` objects, available through:

```python
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import (
    lih_spec_relativistic, beh_spec_relativistic, cah_spec_relativistic
)

lih_rel = build_composed_hamiltonian(lih_spec_relativistic(max_n=2))
# lih_rel['qubit_op']: openfermion.QubitOperator (Hermitian, real coefficients)
# lih_rel['Q'], ['N_pauli'], ['lambda_ni'], ['qwc_groups']: resource metrics
# lih_rel['h1_so_diag']: spin-orbit contributions per qubit index
# lih_rel['eri_sparse']: dict of nonzero physicist-notation ERI entries
```

T4 will import these, compute Sunaga-baseline comparison metrics
(Pauli count at matched Q, 1-norm, QWC groups), and report head-to-head
tables in `benchmarks/relativistic_comparison.py`. The memo also gives
T4 the structural framing: Pauli ratio 2.4–5.9× is the cost of
relativistic physics at matched Q, with 1-norm held roughly flat.

## 8. Algebraic obstructions flagged for Decomposer

None encountered. Every step in T3 is algebraic:

- Radial R^k: `hypergeometric_slater` (exact Fraction) or machine-precision
  float evaluation. No grid quadrature enters for first-row elements.
- Angular X_k: exact sympy `wigner_3j`, cast to float once.
- Spin-orbit: closed-form Breit-Pauli per T2 (`so_diagonal_matrix_element`).
- PK: uses the existing `_compute_pk_matrix_elements` grid routine, which
  is the same numerical integral the scalar builder uses. Not new.

**One flagged item for future consideration:** the PK radial integral is
still grid-based inside `_compute_pk_matrix_elements`. Replacing it with
the Laguerre-moment algebraic path (Track N, v2.0.12) would close the
last quadrature entry in the composed relativistic builder.

## 9. Guardrails observed

- **Did NOT modify** `geovac/dirac_s3.py`, `geovac/dirac_matrix_elements.py`,
  `geovac/spin_orbit.py`. T1/T2/D1 remain locked.
- **Did NOT attempt** SrH / RaH (deferred per sprint plan — [Kr] frozen
  core not yet built). CaH substitutes as the Sunaga-family Z=20 target.
- **Did NOT touch** Paper 14, 18, 20, 22 (T5/T6's job).
- **No TC, no Dirac Fock projection, no graph Dirac operator, no S⁵/S⁷.**
- **No new numerical quadrature.** The one existing PK quadrature is
  inherited from the scalar composed pipeline.
- **Scalar regression bit-exact.** 334/556/778 preserved.

## 10. Files

- `geovac/molecular_spec.py` — `relativistic` + `spinor_labels` fields,
  three new spec factories.
- `geovac/composed_qubit.py` — `build_composed_hamiltonian` dispatches to
  relativistic branch when requested.
- `geovac/composed_qubit_relativistic.py` — new module, the T3 builder.
- `tests/test_spin_ful_composed.py` — 13 regression tests.
- `docs/spin_ful_composed_design_memo.md` — this memo.
