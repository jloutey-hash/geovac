# L1 Infrastructure Reusability Audit — Modular Hamiltonian on T_{n_max}

**Date:** 2026-05-16
**Sprint:** L1-B (scoping; parallel to L1-A architecture, L1-C witness specs)
**Status:** Audit complete; reusability matrix + new-module plan + cost estimate below.
**No production code modified.**

---

## §1. Executive summary

The Bisognano–Wichmann modular Hamiltonian on the truncated Camporesi–Higuchi
spectral triple $\mathcal{T}_{n_\max}$ is **a fundamentally new object** that
the framework does not currently host, but the supporting infrastructure for
its construction, verification, and continuum-lift is **substantially in place**
(roughly 70% reusable directly, 25% reusable with thin extensions, 5% genuinely
new). The single new load-bearing module is `geovac/modular_hamiltonian.py`
(estimated 600–900 lines), wrapping a wedge restriction, a KMS state
constructor, and a modular flow $\sigma_t = e^{itK}$.

**Critical finding (the one structural distinction that drives the audit):**
the existing real structure `RealStructure` (KO-dim 3 charge conjugation,
$J^2=-I$, $JD=+DJ$) is **NOT** the Tomita–Takesaki antilinear conjugation
$J_{T\!T}$ that enters $S = J_{T\!T} \Delta^{1/2}$. They are categorically
different antilinear objects: $J_{GV}$ is a kinematic spinor charge-conjugation
fixed once at construction time, intrinsic to KO-dim; $J_{T\!T}$ is
state-dependent (each KMS state defines its own $J_{T\!T}$ via Tomita's theorem).
Both can coexist on the same Hilbert space, and the L1 module must not conflate
them. The existing `J.apply()` API can be reused as an antilinear-operator
template, but the `RealStructure` data type itself is the wrong physical object.

---

## §2. Module-by-module reusability matrix

| Module | Lines | Tests | Reusable for L1? | Load-bearing functions |
|:-------|------:|------:|:-----------------|:-----------------------|
| `operator_system.py` | 777 | 342 | **Direct YES** (host for K) | `TruncatedOperatorSystem`, `multiplier_matrices`, `contains`, `_extract_matrix_basis` |
| `full_dirac_operator_system.py` | 767 | 511 | **Direct YES** (Dirac D, Hilbert space) | `full_dirac_basis`, `camporesi_higuchi_full_dirac_matrix`, `FullDiracTruncatedOperatorSystem` |
| `spinor_operator_system.py` | 527 | 316 | **Direct YES** (Weyl reduction) | `spinor_basis`, `build_spinor_multiplier_matrix`, `camporesi_higuchi_dirac_matrix` |
| `dirac_s3.py` | 455 | 207 | **Direct YES** (spectrum API) | `dirac_eigenvalue_abs`, `dirac_degeneracy`, `spinor_labels_at_n` |
| `thermal_tensor_triple.py` | 726 | 388 | **Partial YES** (KMS via Matsubara) | `matsubara_spectrum`, `scalar_thermal_free_energy_S3_x_S1`, `tensor_modes` |
| `real_structure.py` | 671 | 566 | **Partial YES** (antilinear operator template ONLY) | `RealStructure.apply`, `RealStructure.apply_to_operator` |
| `berezin_reconstruction.py` | 591 | 627 | **Partial YES** (states, not just observables) | `BerezinReconstruction`, `PlancherelSymbol`, `apply` |
| `gh_convergence.py` | 883 | 545 | **Partial YES** (extend to modular convergence) | `TunnelingPair`, `compute_propinquity_bound`, `gamma_rate` |
| `central_fejer_su2.py` | 1006 | 814 | **Direct YES** (rate estimates for modular lift) | `gamma_rate`, `central_multiplier_cb_norm`, `normalization_constant` |
| `connes_distance.py` | 623 | 207 | **NO** (different optimization) | (none directly; SDP machinery is for sup distances, not modular flow) |
| `circulant_s3.py` | 541 | 219 | **NO** (falsification comparator) | (none for L1) |
| `almost_commutative.py` | 897 | 603 | **NO** (H1 Higgs scope, unrelated) | (none for L1) |
| `gh_convergence_tensor.py` | 1881 | (large) | **Partial YES** (tensor-extension template for L1 follow-on) | (none for L1-B core; relevant only if Hartle–Hawking witness needs a tensor structure) |

Net: ~5,300 lines directly reusable; ~3,000 lines reusable with thin extensions; ~2,000 lines irrelevant. About 28 functions are load-bearing.

---

## §3. Per-module assessment

### 3.1 `operator_system.py` — DIRECT YES (host for K)

**Provides:** `TruncatedOperatorSystem(n_max)` with `.multiplier_matrices`,
`.contains(M)`, `.identity_in_O()`, `.is_star_closed()`, propagation-number
machinery. The Hilbert-space basis is the scalar Fock basis (`HyperLabel(n, l, m)`).

**Reusability for L1.** The modular Hamiltonian $K$ on a wedge sub-system is
naturally a Hermitian element of the operator system on the spinor truncation
(see §3.2), so the `contains` test machinery (Frobenius lstsq projection) is
load-bearing for *verifying* that any candidate $K$ lies in the spanned
operator-system span, and for projecting external trial operators onto $O_{n_\max}$.

**Reuse plan.** No modification. Extend by adding a *wedge-restriction* projector
$P_W$ as a new sparse Hermitian operator (probably built as a new function
`wedge_projection(op_sys, axis, ...)` inside `modular_hamiltonian.py`); the
resulting sub-system $P_W O_{n_\max} P_W$ is itself an operator system whose
self-adjoint elements include $K$. No new API on `operator_system.py` itself.

**Tests.** 342 lines existing; no extension needed.

### 3.2 `full_dirac_operator_system.py` — DIRECT YES (D and the Hilbert space)

**Provides:** `FullDiracTruncatedOperatorSystem(n_max)` with `dim_H =
2*spinor_dim(n_max)` (40 at $n_\max=3$, matching Paper 2 $\Delta^{-1}$); the
truthful CH full Dirac matrix $D = \chi \cdot (n_\text{fock} + 1/2)$ block-diagonal
in chirality; and an offdiag-CH option for SDP boundedness only.

**Reusability for L1.** This IS the Hilbert space the modular Hamiltonian
lives on. The truthful CH Dirac has the property $JD = +DJ$ (Track 2 of WH1-J
sprint), making it the *only* admissible Dirac for any spectral-action-like
construction. The `camporesi_higuchi_full_dirac_matrix` returns
`np.diag([chirality * (n_fock + 0.5) for ...])`, which gives direct
diagonalization access — important because if we can write
$K = 2\pi \cdot \text{wedge}(D)$ for some wedge-resriction prescription, then
$e^{itK}$ is a single `np.exp(1j*t*diag)` on diagonal blocks.

**Reuse plan.** No modification. The wedge restriction adds a $P_W$ projector
on top, and the modular flow $\sigma_t(O) = e^{itK} O e^{-itK}$ is computed
in the full Hilbert space using the existing `D`.

**Tests.** 511 lines (including R3.5 chirality-grading verification). No extension
needed; the L1 module gets its own test file.

### 3.3 `thermal_tensor_triple.py` — PARTIAL YES (KMS state via Matsubara)

**Provides:** `matsubara_spectrum(β, k_max, fermionic)` returning the temporal
$S^1_\beta$ spectrum $\omega_k = 2\pi(k + 1/2)/\beta$ (fermionic) or $2\pi k/\beta$
(bosonic); `scalar_thermal_free_energy_S3_x_S1`; `stefan_boltzmann_factorization`;
and the headline reproduction of Hawking $T_H = 1/(8\pi M)$ on the Euclidean
Schwarzschild cigar (Sprint TD Track 4) at the spectrum level.

**Reusability for L1.** The KMS state at $\beta = 2\pi/\kappa_g$ is structurally
identical to the Matsubara construction for a finite-T QFT. The
`matsubara_spectrum` is the exact temporal-mode quantization implied by
$\sigma_{i\beta}(O) = O$. The L1 KMS state $\omega_\beta(O) = \text{Tr}(e^{-\beta K} O)/Z$
is the *operator-level* analog of what this module already supplies at the
spectrum level — the bridge from "thermal partition function from spatial
spectrum" to "thermal state on the truncated operator system" is short.

**Reuse plan.** Two thin additions needed (in `modular_hamiltonian.py`, not here):
(a) a `KMSState(beta, K)` class that exposes `omega(O) = trace(exp(-beta*K) @ O) / Z`
with an explicitly-computed partition function; (b) a `verify_kms_condition`
helper that checks $\omega(\sigma_t(A) B) = \omega(B \sigma_{t+i\beta}(A))$ at
test pairs $(A, B)$. The Matsubara apparatus stays inside `thermal_tensor_triple`;
the L1 module imports and calls it for sanity checks (e.g., confirming
that the KMS state's $\omega(D^k e^{-\beta D^2})$ trace agrees with the
master-Mellin-engine evaluation at $k=0,1,2$).

**Tests.** 388 lines existing (37 tests); extend with cross-validation tests in the new module.

### 3.4 `real_structure.py` — PARTIAL YES (antilinear API template ONLY)

**Critical structural caveat.** The existing `RealStructure` is the Connes
KO-dim 3 charge conjugation $J_{GV}$ with $J_{GV}^2 = -I$ and $J_{GV} D = +D J_{GV}$
(R2.5 Step 2 closure, 2026-05-05). The Tomita–Takesaki antilinear
conjugation $J_{T\!T}$ is **a different object**:

- $J_{T\!T}^2 = +I$ always (not $-I$).
- $J_{T\!T}$ depends on the choice of cyclic vector / state $|\Omega\rangle$;
  there is one $J_{T\!T}$ per KMS state.
- $J_{T\!T}$ is defined via the polar decomposition $S = J_{T\!T} \Delta^{1/2}$
  where $S O |\Omega\rangle = O^\dagger |\Omega\rangle$ on the cyclic-separating
  vector.

So `RealStructure` is *the wrong kind of antilinear conjugation* for Tomita's
theorem. Mixing them up in the modular-Hamiltonian construction would produce a
silent algebraic error. The L1 module must construct $J_{T\!T}$ separately.

**Reusability for L1.** The mechanical realization of an antilinear operator as
$J(\psi) = U \cdot \text{conj}(\psi)$ with a stored unitary matrix $U$ (and
`apply_to_operator`: $J O J^{-1} = U \cdot \text{conj}(O) \cdot U^T$) is exactly
the right *implementation template*. The L1 module should define a
`TomitaConjugation` class that mirrors this API (`apply`, `apply_to_operator`,
`J_squared_matrix`) but with $J_{T\!T}^2 = +I$ and the cyclic-vector dependence
made explicit. The `_spinor_J_phase` and `_full_dirac_J_target` helpers do
*not* transfer — they encode KO-dim 3 charge conjugation, not modular conjugation.

**Reuse plan.** Borrow the antilinear-operator dataclass *pattern*; do not import
or call `RealStructure` directly in L1. Document in L1's module docstring
that $J_{GV}$ and $J_{T\!T}$ both exist on $\mathcal{H}_{n_\max}$ and are
independent commuting (or non-commuting; this is a sub-question for L1
verification) antilinear involutions.

**Tests.** 566 lines existing on $J_{GV}$; the new $J_{T\!T}$ gets its own test file.

### 3.5 `berezin_reconstruction.py` — PARTIAL YES (extends naturally to states)

**Provides:** `BerezinReconstruction(n_max)` realizing
$B_{n_\max}: C(S^3) \to O_{n_\max}$ via convolution with the L2 central Fejér
kernel followed by truncation; positivity, contractivity, approximate identity,
L3 Lipschitz compatibility. This is the load-bearing map for *observables* in
the continuum-lift direction.

**Reusability for L1.** For lifting the modular Hamiltonian to the continuum
we need the *dual* map for states: a Berezin map $B_*^{n_\max}: \text{States}(O_{n_\max}) \to \text{States}(C(S^3))$ which on KMS states should send
$\omega_\beta^{(n_\max)} \to \omega_\beta^{\text{continuum}}$ in some appropriate
norm. The existing `BerezinReconstruction` directly implements this on
observables; by duality (since $B_{n_\max}$ is positive and unital up to
$\gamma_{n_\max}$), $B_*^{n_\max}(\omega)(f) := \omega(B_{n_\max}(f))$ is the
right dual construction.

**Reuse plan.** Add a `BerezinStateDual(n_max)` class to `modular_hamiltonian.py`
that wraps an existing `BerezinReconstruction` instance and exposes
`apply_to_state(omega)`. The properties (positivity preservation,
approximate-identity rate $\gamma_{n_\max}$) inherit automatically from the
observable-side properties.

**Tests.** 627 lines existing. Extension: 15–20 new tests on the dual map
inside the new L1 test file.

### 3.6 `gh_convergence.py` and `central_fejer_su2.py` — DIRECT YES for rate estimates

**Provides (gh_convergence):** Paper 38 five-lemma chain assembly. `TunnelingPair`
packages L1'–L4 ingredients. `compute_propinquity_bound(n_max)` returns
$\Lambda(\mathcal{T}_{n_\max}, \mathcal{T}_{S^3}) \le C_3 \cdot \gamma_{n_\max}$.

**Provides (central_fejer_su2):** $\gamma_n = O(\log n / n)$ rate with closed-form
small-n values (sympy-exact), $4/\pi$ asymptote (M1 Hopf-base signature).

**Reusability for L1.** The propinquity convergence theorem of Paper 38 is for
the *metric-spectral-triple* propinquity, which controls observables and their
Lipschitz structure. The modular Hamiltonian convergence requires a *modular
propinquity* which is, as of 2026-05-16, an open object in NCG literature
(no Lorentzian extension of Latrémolière 2026 is published). However, the
*rate constant* $\gamma_{n_\max}$ that controls the metric-spectral-triple
convergence almost certainly also controls the modular-Hamiltonian convergence
up to a different (but bounded) prefactor — the M1 mechanism that produces
$4/\pi$ on the metric side is the same M1 mechanism that produces $2\pi$ on
the modular side (per Track D / Unruh four-witness theorem).

**Reuse plan.** Add `modular_propinquity_estimate(n_max, kappa_g)` to the new
L1 module. It calls `gamma_rate(n_max)` for the rate, multiplies by a
prefactor accounting for the wedge-restriction and the choice of surface
gravity, and reports a *qualitative-rate* bound (matching Paper 38's
qualitative-rate-only status). This is honest about scope: L1 produces a
plausibility estimate, not a theorem. The theorem-level lift is the named
follow-up (L2-style sprint on modular propinquity).

**Tests.** 545 + 814 lines existing. No modification.

### 3.7 `connes_distance.py` — NO (different optimization problem)

The SDP machinery in `connes_distance.py` solves a sup-distance optimization
$\sup \{|\phi(x) - \psi(x)| : \|[D, x]\|_\text{op} \le 1\}$ to extract the Connes
distance between two pure states. The modular Hamiltonian construction is not
an optimization at all — $K$ is *defined* by the wedge restriction +
Tomita–Takesaki theorem, not extracted by maximization. The Tomita–Takesaki
theorem itself uses *polar decomposition* of $S$ on a cyclic-separating vector,
which is a numerical linear algebra task (SVD), not a semidefinite program.

**Caveat.** A different SDP application *may* arise as a sanity check (e.g.,
"is the candidate $K$ self-adjoint of unit operator norm on the wedge?"), but
none of `connes_distance.py`'s current entry points implement that test.

**Tests.** 207 lines. Not relevant.

### 3.8 `circulant_s3.py` and `almost_commutative.py` — NO

`circulant_s3.py` is the WH1 R3 falsification comparator; it deliberately
breaks the operator-system structure and is a comparator object only.

`almost_commutative.py` is the Sprint H1 Higgs scoping infrastructure; its
$\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$ + matter/antimatter doubling
is unrelated to wedge restriction / modular flow.

Neither is relevant for L1.

### 3.9 `gh_convergence_tensor.py` — Partial YES (deferred to L1 follow-on)

The tensor-product propinquity machinery (Paper 39) is *not* needed for the
single-factor BW landing on a wedge of $\mathcal{T}_{n_\max}$. It becomes
relevant if a witness specifies the bipartite structure "wedge + complement"
on $\mathcal{T}_{n_\max}$ via a tensor decomposition of the Hilbert space
$\mathcal{H}_{n_\max} = \mathcal{H}_W \otimes \mathcal{H}_{W'}$ (which is the
natural decomposition for entanglement-entropy formulations of modular flow,
e.g. CFT modular Hamiltonians). Defer to L1-C witness specs.

---

## §4. Specifically-asked-for assessments

### 4.1 Wedge restriction as a sub-system of $O_{n_\max}$

Yes, naturally. $P_W O_{n_\max} P_W$ for any orthogonal projection $P_W$ on
$\mathcal{H}_{n_\max}$ is itself a $*$-closed but multiplicatively-non-closed
linear subspace of $M_{N(n_\max)}(\mathbb{C})$, i.e., an operator system. The
existing `TruncatedOperatorSystem.contains` immediately gives a membership
test for the *full* $O_{n_\max}$; restricting it to the wedge sub-system
requires a thin sub-operator-system class (probably ~30 lines).

The non-trivial design question is what defines $P_W$ on $\mathcal{H}_{n_\max}$:
the wedge is a Lorentzian-causal region of S³, but $\mathcal{H}_{n_\max}$ is
parametrized by SO(4) labels $(n, l, m)$ (Riemannian). The standard prescription
(BW theorem on a sphere) is a hemispherical wedge "above an equator" relative
to a chosen polar axis; on $\mathcal{H}_{n_\max}$ this corresponds to a
projection onto states with a parity property under the chosen reflection
(e.g., $P_W = \frac{1}{2}(I + R)$ where $R$ is the reflection-through-equator
operator on $S^3$). The L1-A architecture track should fix the wedge definition;
L1-B notes only that the existing infrastructure supports it.

### 4.2 Real structure J — KO-dim 3 vs. Tomita J

Categorically different. See §3.4. Both can coexist; the L1 module must give
each its own class. Falsifiable sub-question for the L1 sprint: do
$J_{GV}$ and $J_{T\!T}$ commute or anti-commute on the Hartle–Hawking KMS state?
If they commute, the L1 module benefits from a unified antilinear-operator
abstract base class; if they anti-commute or have a non-trivial relation, the
relation itself is publishable structural content.

### 4.3 Truncated Dirac for $e^{itK}$ computation

If $K = 2\pi \cdot D|_W$ (a wedge-restricted Dirac) — which is the
"most-naive" BW prescription — then $e^{itK}$ is a single matrix exponential
applied to a Hermitian operator with known sparse structure (D is diagonal
on the truthful CH basis). The full diagonalization is free: $e^{itK}$ at any
$t$ is `np.diag(np.exp(1j * t * 2*np.pi * K_diag))`. Modular flow on operators
is then $\sigma_t(O) = e^{itK} O e^{-itK}$, a single conjugation per $t$ — O($N^3$)
in $N = \dim \mathcal{H}_{n_\max}$, completely tractable at $n_\max \le 5$.

If the BW prescription involves a non-trivial wedge projection that breaks
diagonality (e.g., $K = \log(\rho_W)$ for the reduced density matrix
$\rho_W = \text{tr}_{W'} |\Omega\rangle\langle\Omega|$), then $K$ is no
longer diagonal in the CH basis and requires a separate Hermitian-matrix
diagonalization. Still O($N^3$); still tractable.

### 4.4 Thermal apparatus as a thin wrapper

Mostly yes. `ModularHamiltonian(beta=2π)` can be built as a class that takes the
Dirac, the wedge projection, and the surface gravity, and exposes
`flow(O, t)`, `kms_state()`, `verify_kms_period()`. The mechanical thermal
spectrum computation in `matsubara_spectrum` is reused for sanity checks
(e.g., the thermal trace $\text{Tr}(e^{-\beta K})$ should equal the product
over Matsubara modes weighted by spatial degeneracies). The wrapper is *not*
purely thin — it has its own logical content (the wedge prescription, the
Tomita conjugation, the KMS verification) — but the existing thermal
apparatus prevents duplicating the Matsubara-sum infrastructure.

### 4.5 GH-convergence extension

Yes, as a rate-estimate adapter. The full theorem-level lift of
$\Lambda_\text{modular}(\mathcal{T}_{n_\max}, \mathcal{T}_{S^3}) \to 0$ requires
a *modular propinquity* (open NCG question). The L1 module can produce a
plausibility estimate by:
1. Computing $\gamma_{n_\max}$ via `gamma_rate(n_max)`.
2. Computing the wedge-restriction operator norm bound on $\|K - K_{S^3}|_W\|$.
3. Reporting a heuristic $\Lambda_\text{modular} \le C_\text{wedge} \cdot \gamma_{n_\max}$
   (qualitative-rate, mirroring Paper 38's scope).

The full theorem-level lift is a named multi-month follow-on.

### 4.6 Connes-distance SDP

Not directly relevant. The polar decomposition of $S$ is plain SVD on a
$N \times N$ matrix; no SDP. The wedge "size" parameter (one of the natural
choices: the spectral width of $P_W$, the trace of $P_W D P_W$ over
$\text{tr}(P_W)$) might be optimized via a small LP, but cvxpy is overkill
for that.

### 4.7 Berezin reconstruction for states

Yes, by duality. See §3.5. The natural definition
$B_*^{n_\max}(\omega)(f) := \omega(B_{n_\max}(f))$ inherits positivity,
$\|B_*^{n_\max}(\omega)\| = \omega(B_{n_\max}(1)) = \omega(I)$ (so it preserves
unital states), and approximate-identity rate $\gamma_{n_\max}$. The L1 module
adds ~50 lines of dual-map adapter code on top of existing
`BerezinReconstruction`.

---

## §5. New modules to build

### 5.1 `geovac/modular_hamiltonian.py` — load-bearing new module (estimated 600–900 lines)

Single module. Hosts all the new content. Logical structure:

```python
# Tomita conjugation (J^2 = +I, state-dependent)
class TomitaConjugation:
    """Antilinear J_TT from polar decomposition of S on a cyclic vector."""
    def __init__(self, dirac, cyclic_vector): ...
    def apply(self, psi): ...                 # antilinear
    def apply_to_operator(self, O): ...       # O -> J_TT O J_TT^{-1}
    def J_squared_matrix(self): ...           # should equal +I

# Wedge restriction
class WedgeProjection:
    """Orthogonal projection P_W onto a Lorentzian-causal wedge sub-Hilbert."""
    def __init__(self, op_sys, axis, kind='hemispherical'): ...
    @property
    def P(self) -> np.ndarray: ...
    def restrict(self, O: np.ndarray) -> np.ndarray:  # P_W O P_W
        ...

# Modular Hamiltonian
class ModularHamiltonian:
    """K such that sigma_t = e^{itK} is the wedge-modular flow."""
    def __init__(self, op_sys, full_dirac, wedge, surface_gravity_kappa_g): ...
    @property
    def beta(self) -> float: return 2 * np.pi / self.kappa_g
    @property
    def K(self) -> np.ndarray: ...                       # Hermitian
    def flow(self, O: np.ndarray, t: complex) -> np.ndarray:
        """sigma_t(O) = e^{itK} O e^{-itK}, real or imaginary t."""
    def verify_kms_period(self, O: np.ndarray, tol=1e-9) -> tuple[bool, float]:
        """sigma_{i*2*pi}(O) == O up to tol."""

# KMS state
class KMSState:
    """omega_beta(O) = Tr(exp(-beta K) O) / Z."""
    def __init__(self, modular_hamiltonian): ...
    def omega(self, O: np.ndarray) -> complex: ...
    def verify_kms_condition(self, A, B, t, tol=1e-9) -> tuple[bool, float]: ...

# Berezin dual map for states
class BerezinStateDual:
    """B_*^{n_max}(omega)(f) := omega(B_{n_max}(f))."""
    def __init__(self, berezin: BerezinReconstruction): ...
    def apply_to_state(self, omega: KMSState) -> Callable[[np.ndarray], complex]: ...

# Modular convergence rate estimate
def modular_propinquity_estimate(n_max, kappa_g) -> dict:
    """Qualitative-rate L_mod(T_{n_max}, T_{S^3}) bound."""
    ...
```

### 5.2 Other new modules — none required at L1-B scope

A separate `geovac/wedge_restriction.py` could host the wedge projection if
it grows; if the wedge is a simple algebraic object (projection onto specific
$(n, l, m)$ patterns), it fits inside `modular_hamiltonian.py`. Recommend keeping
the L1 build as **one module** to avoid premature decomposition.

`geovac/kms_state.py` similarly fits inside the parent module — KMS state is
50–80 lines of explicit Hermitian-matrix exponential code.

---

## §6. Test extension plan

| Existing test | Extend? | Action |
|:--------------|:--------|:-------|
| `test_operator_system.py` | No | No new tests; wedge sub-system is tested in new file |
| `test_real_structure.py` | No | Existing tests on $J_{GV}$ stay; do not conflate |
| `test_full_dirac_operator_system.py` | No | Existing D verification stays |
| `test_thermal_tensor_triple.py` | Lightly (2-3 tests) | Cross-check: KMS state from modular Hamiltonian reproduces `matsubara_spectrum` partition function |
| `test_berezin_reconstruction.py` | Lightly (3-5 tests) | Cross-check that `BerezinStateDual` inherits positivity / contractivity |
| `test_gh_convergence.py` | Lightly (2 tests) | Cross-check that modular rate estimate is bounded by metric rate |
| `test_modular_hamiltonian.py` | **NEW** (~400-600 lines) | All of L1's verification: $J_{T\!T}^2=+I$, wedge $P_W^2=P_W$, KMS period $\sigma_{i\beta}=\text{id}$, KMS condition $\omega(\sigma_t(A) B) = \omega(B \sigma_{t+i\beta}(A))$, four-witness checks (Hawking-T, Sewell, BW, Unruh) |

Total new test code: roughly 500 lines; total extensions to existing tests:
roughly 100 lines.

---

## §7. Bottlenecks (existing modules needing API changes vs additions)

**No API changes needed in any existing module.** This is the clean finding
of the audit. The L1 module is purely additive: it imports from
`operator_system`, `full_dirac_operator_system`, `berezin_reconstruction`,
`thermal_tensor_triple`, and `central_fejer_su2`, and exposes its own classes.

Three soft bottlenecks (additions, not changes):

1. `BerezinReconstruction` does not currently expose its underlying convolution
   kernel as a callable object. The dual-state construction is cleaner if
   `BerezinReconstruction` gains a `kernel(f) -> Callable[[g], float]` method
   returning the convolved function. Not strictly needed — the dual can also
   call `BerezinReconstruction.apply` repeatedly — but cleaner.

2. `FullDiracTruncatedOperatorSystem` does not currently expose an explicit
   `chirality_grading()` operator (the matrix $\gamma$ with $\gamma^2 = I$
   and $\gamma D = -D\gamma$). The Tomita conjugation construction may benefit
   from it. Adding it is ~10 lines.

3. `thermal_tensor_triple.matsubara_spectrum` returns a sympy `Expr` symbol
   list; the L1 module needs numerical evaluation, which is one `mp.mpf`
   conversion call. Not a bottleneck.

---

## §8. Implementation cost estimate

Assuming one full-session agent (single dispatch, all-tracks-in-parallel
discipline, ~6-8 hours focused work):

| Phase | Estimated effort |
|:------|-----------------:|
| Architecture scaffolding (class definitions, docstrings) | 1.5 h |
| `TomitaConjugation` (polar decomposition + verification) | 1 h |
| `WedgeProjection` (one canonical wedge, e.g. hemispherical-in-$n$) | 1 h |
| `ModularHamiltonian.K` computation + `.flow(O, t)` | 1.5 h |
| `KMSState` (trace, partition function, KMS verification) | 1 h |
| `BerezinStateDual` (50 lines + 5 tests) | 0.5 h |
| Test file (verification of all four witnesses at $n_\max = 2, 3$) | 1.5 h |
| Documentation + memo | 1 h |

**Total: 8–9 agent-hours**, fitting one full session with light overrun risk.
The risk concentrates in the wedge-projection-on-$S^3$ definition (§4.1)
and in the four-witness verification (Hawking-T, Sewell, BW, Unruh require
slightly different surface-gravity conventions; getting all four to land at
$\sigma_{i\cdot 2\pi} = \text{id}$ simultaneously is the load-bearing test).

**Critical-path note for L1-A architecture track.** L1-A should pin the
wedge definition (axis choice, parity prescription, signature of $P_W$ on the
truthful CH basis) before L1 implementation starts. If L1-A leaves this open,
the implementation agent must pick a default (recommended: hemispherical
in the *first* basis label coordinate, $P_W = \frac{1}{2}(I + R_{e_1})$ where
$R_{e_1}$ is reflection through the equator orthogonal to the first SO(4)
direction); the choice is not unique but is structurally informative either way.

---

## §9. Honest scope statement

What the L1 build will produce:
- A working `ModularHamiltonian` class on $\mathcal{T}_{n_\max}$ at $n_\max \le 5$.
- KMS state, modular flow, KMS-period verification at $\sigma_{i \cdot 2\pi}$.
- Four-witness reproduction at the *finite-cutoff structural* level (numerical
  agreement to machine precision on quantities the construction is built to
  reproduce; failure on quantities it is not).
- A qualitative-rate estimate of modular-propinquity convergence.

What L1 will *not* produce:
- A theorem-level modular-propinquity convergence (open NCG question; multi-month).
- The full Tomita–Takesaki spectrum at the continuum limit (L1-C may sketch).
- Any claim about cross-manifold modular structure (Paper 24 §V W2b
  blocker still applies).
- The literal identification of the framework's $J_{GV}$ with $J_{T\!T}$
  (they are different objects per §3.4).

The audit confirms the L1 build is **scoped at a single agent-session** and
that the necessary infrastructure is **70%+ in place**. The remaining 30%
is one new module of well-defined logical content.

---

## §10. Files referenced

Production modules audited:
- `geovac/operator_system.py`
- `geovac/full_dirac_operator_system.py`
- `geovac/spinor_operator_system.py`
- `geovac/dirac_s3.py`
- `geovac/thermal_tensor_triple.py`
- `geovac/real_structure.py`
- `geovac/berezin_reconstruction.py`
- `geovac/gh_convergence.py`
- `geovac/central_fejer_su2.py`
- `geovac/connes_distance.py` (negative — included for completeness)
- `geovac/circulant_s3.py` (negative)
- `geovac/almost_commutative.py` (negative)
- `geovac/gh_convergence_tensor.py` (partial — deferred to L1 follow-on)

Memos / context:
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` §III, §IV, §VIII
- `papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex`
- `debug/bisognano_wichmann_track_d_memo.md`
- `debug/unruh_pendant_memo.md`
- `CLAUDE.md` §1.7 WH1 PROVEN, §2 R3.5 / R2.5 / Berezin / J entries

No production code or papers were modified by this audit.
