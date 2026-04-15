# Dirac-on-S³ Infrastructure — Design Memo (Track D1)

**Status:** Complete, v1 — 2026-04-15.
**Module:** `geovac/dirac_s3.py`
**Tests:** `tests/test_dirac_s3.py` (51/51 passing under exact sympy arithmetic).
**Scope:** Tier 1 sprint, Track D1 only. Matrix-element engines, Dirac
Fock-projection theorems, and graph-edge discretizations are explicitly
out of scope (deferred to later tracks or Tier 2).

---

## 1. Purpose

Provide the minimal infrastructure that Tracks D2, D3, D4 need to
investigate whether the Dirac-on-S³ sector promotes Paper 2's
K = π(B + F − Δ) combination rule from a three-tier coincidence into a
common-generator theorem. Concretely:

- A π-free rational spectrum on unit S³.
- A labeling of spinor spherical harmonics compatible with existing
  scalar (n, l, m) Fock graph nodes, augmented by spin / chirality σ.
- A certification that the spectrum and degeneracies live in ℚ — direct
  analog of Paper 24 §III's Bargmann–Segal certificate.

Nothing in D1 computes a single floating-point number.

---

## 2. API surface

### Spectrum
- `dirac_eigenvalue_abs(n, convention=...)` → `sympy.Rational` = `n + 3/2`.
- `dirac_degeneracy(n, sector=..., convention=...)` → `int`
  (2(n+1)(n+2) Dirac, (n+1)(n+2) Weyl).
- `delta_inverse_identity()` → `(40, Rational(1, 40))` — the Phase 4H
  SM-D identity Δ⁻¹ = g₃^Dirac.

### Spinor harmonics (labels only; no analytic expressions)
- Dataclass `SpinorHarmonicLabel(n_fock, l, m, sigma, chirality)`.
  Stored integer `m` holds **2·m_j** (twice the total angular momentum
  projection), so that all dataclass fields except `sigma` are plain
  `int`, and `sigma`/`m_j` retain their natural half-integer algebra
  through `sympy.Rational`. Documented in the docstring.
- `spinor_labels_at_n(n, sector=..., convention=...)`.
- `iter_spinor_labels(n_max, ...)` / `count_spinor_labels(n_max, ...)`.

### Certification
- `verify_pi_free(n_max, sector=...)` → `bool`. Verifies every
  eigenvalue is a `sympy.Rational` and every degeneracy is a positive
  `int` up to n_max. Also checks the label generator produces exactly
  `g_n` labels per level — this is a stronger invariant than the bare
  π-free claim, and catches label-generator bugs immediately.

### Convention conversion
- `fock_to_ch`, `ch_to_fock`. Invertible, tested.

---

## 3. Label mapping: (n, l, m, σ) ↔ Camporesi–Higuchi

The existing GeoVac scalar graph (`geovac/lattice.py`) uses

```
n_Fock ∈ {1, 2, ...},   l ∈ {0, ..., n_Fock - 1},   m ∈ {-l, ..., l}.
```

Camporesi–Higuchi's natural index is `n_CH ∈ {0, 1, 2, ...}`, with
Spin(4) = SU(2)_L × SU(2)_R irreps

- ψ₊ sector: (j_L, j_R) = ((n_CH+1)/2, n_CH/2), dimension
  (n_CH+1)(n_CH+2).
- ψ₋ sector: (n_CH/2, (n_CH+1)/2), same dimension.

**Convention adopted:** n_Fock = n_CH + 1. This preserves compatibility
with the scalar graph while letting D2/D3/D4 use either index.

Under the diagonal SU(2) (angular momentum on an embedded S²), each
chirality irrep decomposes as j = 1/2, 3/2, …, n_CH + 1/2. Writing
j = l + 1/2 with l ∈ {0, …, n_CH} gives the orbital label l. Each j
contributes 2j + 1 values of m_j; summing yields (n_CH+1)(n_CH+2),
confirming the count per chirality.

**Forgetful map to scalar graph:** dropping σ and chirality takes a
spinor label (n_Fock, l, m_j, σ, c) to a "nearby" scalar node, but the
map is not a strict 2-to-1 lift onto scalar (n, l, m), because m_j is
half-integer. D2/D3 do not need a forgetful map; D4 (Hopf-equivariant)
will use the fiber-charge q = m_L − m_R (via the SU(2)_L × SU(2)_R
decomposition), which is well-defined at the label level here. Any
track that wants a strict orbital-to-spinor lift should add a `j`
quantum number and use Clebsch–Gordan coefficients — that is a D2+
concern and intentionally not implemented in D1.

---

## 4. Weyl vs. full Dirac

Both sectors are exposed; the caller chooses via `sector="weyl"` or
`sector="dirac"`. Use cases:

- **D2 (Dirac analog of B = 42):** use `sector="dirac"`. Paper 2's
  B is a *full spectral* invariant for the scalar case; the natural
  Dirac analog sums over both chiralities.
- **D3 (Dirac analog of F):** can try both. The Dirichlet series
  D_{g_n^{Dirac}}(s) at s = 4 doubles D_{g_n^{Weyl}}(s) by
  construction, so they differ by a factor of 2 only.
- **D4 (Orduz / Hopf-equivariant):** `sector="dirac"` — the 40 states
  at n = 3 that must be decomposed over fiber charges q.

The Weyl sector is the "minimal" version used when only one chirality
makes sense (e.g., Witten-index computations, or if D3 finds that a
single-chirality sum gives a cleaner ζ identity than the Dirac one).

---

## 5. Convention choices and the alternatives rejected

| Choice | Alternative considered | Reason for the choice |
|:---|:---|:---|
| n_Fock = n_CH + 1 | Keep CH's n ≥ 0 and retire Fock's n ≥ 1 | Existing scalar graph nodes use Fock n ≥ 1; a global reindex would churn every Paper 7/13/18 reference. The `convention=` kwarg lets both live side-by-side. |
| σ ∈ ±1/2 as `sympy.Rational` | `+1`/`-1` integer | Matches physical half-integer spin and the natural σ = chirality/2 identity within a single-chirality sector. |
| Store 2·m_j as int in dataclass field `m` | Store m_j as sympy.Rational | Keeps hashability and frozen-dataclass semantics clean; a half-integer m_j is recoverable as `Rational(m, 2)`. Documented in the class docstring. |
| Label orbital quantum number by `l` with j = l + 1/2 implicit | Add explicit `j` field | D2/D3/D4 do not need j at the label layer; adding it forces a decision on how to handle the two σ values within a fixed j, which Bär's and the CH conventions handle differently. Deferred to D2+. |
| No matrix elements in D1 | Provide a placeholder spinor-harmonic callable | Explicit out-of-scope per the sprint plan. D2/D3 build matrix elements if and only if they need them, using the labeling committed here. |
| π-free certificate on *spectrum* only | Certify normalized harmonics too | The harmonics carry conventional √π normalization; that is a calibration in Paper 18's taxonomy, not a physical transcendental. Paper 24 §III does the same: Bargmann–Segal edge weights are certified rational, even though the Bargmann kernel itself carries π through dμ. |

---

## 6. Transcendental content (Paper 18 taxonomy)

**Every quantity returned by `dirac_s3.py` is in ℚ.**

- Intrinsic: none.
- Calibration: the √π normalization of spinor harmonics is a
  calibration (choice of inner-product measure). It does **not** appear
  in any routine here.
- Embedding: none. The Dirac operator on unit S³ has no embedding
  exchange constant — it is the natural spinor Laplace–Beltrami object.
- Flow: none.

This is the same status as Paper 24's Bargmann–Segal lattice (bit-exact
π-free in rational arithmetic), and is why D2/D3/D4 can proceed in pure
symbolic sympy without any floating-point fallback.

If D2 or D3 finds that a particular Dirichlet combination equals a
specific value of ζ or π^k, that transcendental enters through the
*analysis of the infinite sum*, not through the spectrum itself — the
same pattern as Phase 4F α-J (F = D_{n²}(d_max) = ζ(2)).

---

## 7. Guardrails respected

- **No local discrete Dirac on the Fock graph.** The Dirac operator in
  this module lives on the continuum spinor bundle over S³, not on a
  simplicial complex. Ginsparg–Wilson's theorem (no local
  chirality-preserving discrete Dirac without fermion doublers on a
  finite graph) is not violated because nothing here claims to be a
  graph-edge discretization. D1 provides spectrum + labels only.
- **No Paper 2 edits.** Paper 2 §IV rewrite is D5's responsibility.
- **No matrix element engines.** Deferred.
- **No Dirac Fock-projection theorem.** Deferred (Tier 2).

---

## 8. What D2/D3/D4 will need to know about this API

1. `dirac_degeneracy` and `dirac_eigenvalue_abs` are the canonical
   entry points for any symbolic computation. They never return float.
2. `count_spinor_labels(n_max, sector="dirac")` is the exact truncated
   Hilbert-space dimension for the full Dirac sector; use it when
   comparing to Paper 2's B = 42 (where the scalar truncation is at
   m = 3, meaning CH levels 0..2 — see D2 memo).
3. The `SpinorHarmonicLabel` dataclass is frozen and hashable, so it
   can be used as a dict key or set element for bookkeeping
   (D4 Hopf-charge decomposition is the likely first consumer).
4. If a track wants to iterate the chirality-+ sector only, filter
   Dirac labels by `chirality == +1` — this is identical to iterating
   the Weyl sector.
5. If a track wants the stored integer `m` as a half-integer, use
   `Rational(label.m, 2)`.
6. **The module does NOT yet commit to a particular j labeling inside
   the (l, σ) basis.** A D2+ track that builds Wigner–Eckart matrix
   elements will need to add an explicit j field and decide whether to
   use the ( (n+1)/2, n/2 ) SU(2)×SU(2) basis directly or the (j, m_j)
   diagonal-SU(2) basis. This design decision is deferred, not
   pre-empted, by D1's labeling.

---

## 9. References

- R. Camporesi and A. Higuchi, "On the eigenfunctions of the Dirac
  operator on spheres and real hyperbolic spaces", *J. Geom. Phys.*
  20 (1996) 1–18, arXiv:gr-qc/9505009. (Spectrum on S^d.)
- C. Bär, "The Dirac operator on space forms of positive curvature",
  *J. Math. Soc. Japan* 48 (1996) 69–83. (Explicit spinor harmonics
  on S³, (l, j, m_j) labeling.)
- GeoVac Paper 24 §III (π-free rational certification pattern for
  Bargmann–Segal lattice on S⁵).
- GeoVac Paper 18 (exchange constants taxonomy).
- GeoVac CLAUDE.md §2 Phase 4H SM-D: Δ⁻¹ = g_3^Dirac(S³) = 40.
- GeoVac `agents/DECOMPOSER.md` and `docs/dirac_s3_tier1_sprint_plan.md`
  (governing algebraic-first philosophy for the sprint).
