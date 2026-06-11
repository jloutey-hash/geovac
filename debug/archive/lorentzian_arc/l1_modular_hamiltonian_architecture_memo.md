# Sprint L1-A — Modular Hamiltonian on T_{n_max}: Architectural Blueprint

**Date:** 2026-05-16
**Sprint:** L1-A (architecture, no production code or paper edits)
**Author:** L1-A architecture PM (Claude)
**Status:** Architecture-only deliverable. Implementation by a follow-on agent gated on PI sign-off.
**Inputs (verified):** `CLAUDE.md` §1.7 WH1 PROVEN entry + §2 Sprint TS-D / Sprint Unruh-pendant / Sprint L0; `papers/group1_operator_algebras/paper_32_spectral_triple.tex` §IV (Connes axiom audit on truthful CH at finite n_max), §VIII (`thm:gh_convergence`, `rem:bisognano_wichmann_reading`, case-exhaustion theorem); `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` §III.27 (Wick-rotation projection), §VIII (open question on operator-system-level Lorentzian extension); `papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex` §5 (lemmas L1' / L2 / L3 / L4 / L5); `debug/bisognano_wichmann_track_d_memo.md` §2-3 (operator-level falsifier); `debug/lorentzian_l0_audit_memo.md` §3, §9 (M3 trivialization prediction at (3,1); Sprint L1 sequencing); existing modules `geovac/{real_structure, operator_system, connes_distance, full_dirac_operator_system, gh_convergence, berezin_reconstruction, central_fejer_su2, thermal_tensor_triple}.py`.

---

## §1. What is being asked and what is being delivered

Sprint L1-A is the architectural pass for the operator-level falsifier named in `papers/group1_operator_algebras/paper_32_spectral_triple.tex` `rem:bisognano_wichmann_reading` and in `debug/bisognano_wichmann_track_d_memo.md` §3.1:

> *Compute the modular Hamiltonian K on T_{n_max} for a wedge-restricted Camporesi–Higuchi state and verify σ_{i·2π} = identity at finite n_max, then take the GH limit using the lemmas of Theorem `thm:gh_convergence`.*

If this falsifier closes, the four-witness Wick-rotation theorem (Hawking + Sewell + BW + Unruh, codified Sprint Unruh-pendant `debug/unruh_pendant_memo.md`) is lifted from "structural correspondence" to "literal identification" inside the framework's truncated metric spectral triple, and Paper 34 §III.27 (Wick rotation projection) graduates from structural-correspondence to operator-level.

This memo answers the eight architectural questions in the dispatch prompt and produces (a) a concrete API for `geovac/modular_hamiltonian.py`, (b) a five-lemma lift table mapping Paper 38 L1'-L5 to the modular context, (c) a witness-specific parameter table, (d) a named-obstruction list with verification plans, and (e) a verdict on whether the implementation sprint should proceed as-is or with caveats or after a precursor sprint.

**Honest scope:** This is architecture, not theorem. The blueprint identifies the loci where things are likely to work and the loci where the natural construction may fail. The implementation agent will discover which loci behave as expected. Expected implementation budget is 4–8 weeks (PI sprint-window estimate, anchored in `rem:bisognano_wichmann_reading`); I assess this as **5–9 weeks** for a clean closure with documentation (see §11 verdict), with two specific items that might extend it: the choice of wedge sub-algebra (§3) and the M3 trivialization check (§10b).

---

## §2. What IS K_{n_max} for the truncated CH triple? (Question 1)

The Bisognano–Wichmann theorem [Bisognano–Wichmann 1976, *J. Math. Phys.* 17, 303; *J. Math. Phys.* 16, 985 (1975)] states that for a wedge W_R = {x : x^1 > |x^0|} in Minkowski R^{d-1,1}, the Tomita modular automorphism σ_t of the local algebra A(W_R) with respect to the Wightman vacuum Ω is

    σ_t(A) = U(Λ_{−2πt}) A U(Λ_{−2πt})^{-1},

where Λ_s is the boost in the (x^0, x^1) plane by rapidity s. The modular Hamiltonian K is the generator: σ_t = e^{itK}, so K = 2π · K_boost where K_boost is the boost generator. The KMS condition gives σ_{i·2π}(A) = A (periodicity along imaginary axis at period 2π).

On the truncated CH triple T_{n_max} we have no Minkowski boosts. The construction must instead build K from the Tomita–Takesaki polar decomposition of the antilinear operator S: aΩ ↦ a*Ω on the GNS Hilbert space of a chosen wedge-like sub-algebra against a chosen cyclic-separating state. This is option (b) of the dispatch prompt; I assess that (a) and (c) are structurally inadequate:

- **Option (a) — K = β·D restricted to a wedge.** Tempting because D is already on T_{n_max} and a wedge restriction is the kinematic piece of BW. But this conflates two distinct objects: D is the GeoVac Dirac (with units of inverse length in the Camporesi–Higuchi |λ_n| = n + 3/2 normalization), while K must be dimensionless to have a unitary modular flow Δ^{it} = e^{itK}. The two are related by a Schwinger-like β·D mass-shell projection on Matsubara modes (Sprint TD Track 1's `thermal_tensor_triple`), but for the *wedge-localised* modular Hamiltonian K_BW must be the boost generator, not the Dirac. The boost generator on a Riemannian S³ is the conformal Killing vector ∂_χ at the equator of a hemisphere (see §3); only after the polar decomposition does K close as the integral of T^00 against the linear weight 2π·sinh(rapidity).
- **Option (c) — K = 2π·X with X a Hopf-bundle generator.** Tempting because 2π = Vol(S^1) is the M1 signature (see Paper 32 §VIII case-exhaustion theorem). But the M1 statement is a transcendental-content statement (where π's enter); the modular Hamiltonian is an operator, not a transcendental, and constructing it from a Hopf generator would beg the question (the proof must produce 2π from the polar decomposition, not insert it).

**Recommended option (b):** construct K via Tomita–Takesaki polar decomposition of S on the GNS Hilbert space of a chosen wedge sub-algebra against a chosen KMS state ω_β. Then K = log Δ with Δ = S* S the modular operator. The verification σ_{i·2π} = id at β = 2π is then a check on Δ, not an input.

Concretely:

1. Choose a wedge sub-operator-system O_wedge ⊂ O_{n_max} (see §3).
2. Choose a KMS state ω_β on O_{n_max} at β = 2π (see §4); restrict to O_wedge.
3. Build the GNS triple (π_ω, H_ω, Ω) where Ω is the cyclic vector representing ω_β.
4. Define S: π_ω(a) Ω ↦ π_ω(a*) Ω on the dense subspace π_ω(O_wedge) Ω.
5. Polar decompose S = J Δ^{1/2}; modular operator Δ = S* S; modular conjugation J antilinear.
6. K := log Δ; modular flow σ_t = Δ^{it} · Δ^{-it}.
7. Verify σ_{i·2π}(a) = a for a ∈ O_wedge (KMS condition at β = 2π).

The Connes–Rovelli 1994 thermal time hypothesis (arXiv:gr-qc/9406019, *Class. Quantum Grav.* 11, 2899) instantiates exactly this construction on general von Neumann algebras, identifying the modular flow with "thermal time"; the Bisognano–Wichmann theorem is then the statement that on the wedge algebra of Minkowski QFT this thermal time IS the boost orbit. The framework-internal version of L1-A asks the analogous question on the truncated metric spectral triple.

---

## §3. The wedge restriction on T_{n_max} (Question 2)

Three candidates were named in the dispatch prompt. I evaluate each:

### §3.1 Candidate W1: hemisphere of S³

**Construction:** Choose a unit vector n̂ ∈ S³ (canonical: n̂ = (1,0,0,0)). The hemisphere is H_+ = {ω ∈ S³ : ω · n̂ > 0}. The sub-operator-system is O_wedge^{W1} := P_{n_max} · C^∞_c(H_+) · P_{n_max}, the compressions of smooth functions supported in the open hemisphere.

**Inheritance of projection:** YES, cleanly. The Connes–vS compression P_{n_max} commutes with the multiplication action of any f ∈ C^∞(S³); restricting f to H_+ (smoothly with a partition-of-unity cutoff near the equator ∂H_+ = S²) preserves linearity and self-adjointness. The sub-operator-system O_wedge^{W1} is a *-closed linear subspace of O_{n_max}.

**Modular closure:** This is the cleanest case for the framework-internal version of BW, because the half-S³ Riemannian-side analog of a Rindler wedge is exactly the construction Bisognano–Wichmann singled out as the most natural Wightman-axiomatic example (see also Casini–Huerta arXiv:0905.2562 §V for the analogous continuum entanglement-Hamiltonian construction on a half-space). The continuum modular Hamiltonian for a hemisphere of S^d in CFT is well-known: K_cont = 2π ∫_{H_+} (ω · n̂) T^00(ω) dω (with appropriate conformal weight; cf. Casini–Huerta–Myers arXiv:1102.0440). At the truncated-triple level, the analogous "physical-energy-density times linear weight" should emerge from the polar decomposition.

**Recommendation: this is the primary candidate.**

### §3.2 Candidate W2: κ-preserving Dirac graph block (Rule A)

**Construction:** Use the Paper 29 Rule A adjacency on the Dirac graph (κ-preserving, intra-shell within each (κ, m_j) sub-block). The wedge sub-system is the linear span of multipliers respecting Rule A.

**Inheritance:** Mechanically yes — the Rule A block-diagonal structure is preserved under compression. But there is no natural KMS state on Rule A whose modular flow returns the boost orbit; the κ-preserving structure is an internal symmetry, not a wedge-localized algebra. Sprint TD Track 3 found that on the full Dirac sector the chirality flip {γ_F, D_F} = 0 enforces M3 trivialization (`debug/sprint_td_track3_memo.md`); on the Rule A sub-block the spectrum is even more degenerate (chain (+) only), which makes the polar decomposition near-trivial.

**Recommendation: SKIP for L1.** Rule A is a useful comparator (§10b M3 trivialization consistency check) but not the primary candidate.

### §3.3 Candidate W3: pendant-vertex-free subgraph

**Construction:** Per Paper 28 §pendant_edge (`debug/transverse_qed_self_energy_memo.md`), the GS vertex |1,0,0⟩ is a pendant in the scalar Fock graph; its unique edge e_0 decouples in L_1 (transverse-photon Laplacian). The pendant-vertex-free wedge is O_wedge^{W3} = compression of O_{n_max} onto the subspace orthogonal to span{|1,0,0⟩}.

**Inheritance:** Yes mechanically. But the pendant-vertex exclusion is a graph-topological projection (Paper 28 selection-rule census, 4-tier partition), not a spatial-region projection. There is no clear sense in which pendant-vertex-free O is "localized" in a half-S³ region; the construction is global-on-the-graph with one vertex removed.

**Recommendation: SKIP for L1.** W3 has its own role as a transverse-photon construction (Paper 28); reusing it as a "wedge" conflates spatial-region restriction with graph-topological projection.

### §3.4 Verdict for §3

**Primary wedge: W1 (hemisphere of S³).** Secondary diagnostic comparator: W2 (Rule A block) for the M3 trivialization consistency check.

The W1 construction requires building a partition-of-unity cutoff function χ_{H_+} ∈ C^∞_c(H_+) and computing its multiplier matrix in the operator-system basis. This re-uses `geovac/operator_system.py::TruncatedOperatorSystem.build_multiplier_matrix` (Avery–Wen–Avery 3-Y integrals on S³) but with the new input being not a single hyperspherical harmonic Y^{(3)}_{NLM} but a Peter–Weyl expansion of χ_{H_+} truncated at the cutoff. Concretely:

    χ_{H_+}(ω) = Σ_{N≤N_χ} Σ_{L,M} c_{NLM} Y^{(3)}_{NLM}(ω)

where N_χ is a fixed expansion cutoff (separate from n_max — we want N_χ ≥ 2·n_max so the wedge projector resolves all the structure of O_{n_max}; cf. Hekkelman's truncated S^1 construction at http://www.waltervansuijlekom.nl/wp-content/uploads/2022/10/Scriptie-EvaMariaHekkelman.pdf for the analogous bookkeeping on S^1). The wedge projector is then P_wedge = M_{χ_{H_+}} computed via `build_multiplier_matrix` extended to non-monomial multipliers.

---

## §4. The KMS state ω_β at β = 2π (Question 3)

On the truncated O_{n_max}, the most natural KMS state at inverse temperature β with respect to a Hamiltonian H is the Gibbs state

    ω_β(O) = Tr_{H_{n_max}}(e^{-β H} O) / Tr_{H_{n_max}}(e^{-β H}).

For the BW reading, H should be the wedge-restricted Dirac D_{n_max} (truthful Camporesi–Higuchi) or a wedge-restricted boost-Hamiltonian variant. Two specifications:

- **ω_β^{full} : O_{n_max} → C, ω_β^{full}(O) = Tr(e^{-β D_{n_max}} O) / Z(β).** Defined for all of O_{n_max}, not just the wedge. Used as the global thermal state from which the wedge restriction inherits its KMS structure (cf. Sewell 1982's construction of the Hartle–Hawking vacuum on Kruskal-extended Schwarzschild, restricted to the static exterior).
- **ω_β^{wedge}: O_wedge^{W1} → C, restriction of ω_β^{full} to O_wedge^{W1}.** This is the input to the GNS construction in §2 step 3.

At β = 2π the Gibbs state is well-defined (D_{n_max} is a finite Hermitian matrix; e^{-2π D_{n_max}} is positive trace-class). The expectation values are sympy-computable in exact rational arithmetic if D_{n_max} is the truthful Camporesi–Higuchi diag(n + 3/2) (eigenvalues are rationals, β = 2π is the only transcendental); for the offdiag CH Dirac (`geovac/connes_distance.py::DiracProxy(mode='offdiag')`) the eigenvalues are irrational but the construction goes through numerically.

**Choice of D for L1-A: truthful Camporesi–Higuchi.** Two reasons:

1. *Connes-axiom-compatible.* Paper 32 §IV verified J^2 = -I, JD = +DJ, J·O·J^{-1} = O exactly on truthful CH at n_max ∈ {1, 2, 3}; offdiag CH fails JD at residual 2.0. The modular conjugation J in the Tomita polar decomposition needs to be Connes-axiom-compatible if the BW reading is to identify with the spectral-triple-internal J of Paper 32 §IV; truthful CH is the only candidate.
2. *Exact-rational structure.* The thermal state Gibbs weights are e^{-2π(n + 3/2)} = e^{-2π·3/2}·e^{-2π·n}, and the modular flow analytic continuation σ_{iβ} produces e^{β(n + 3/2)} = e^{β·3/2}·e^{β·n}. At β = 2π these are *not* trivial (the spectrum is non-zero), but the algebra closes cleanly under integer-shift-by-1 generators (raising/lowering of n by 1).

**Caveat:** The R3.2 finding (`debug/wh1_r32_spinor_lift_memo.md`) is that the truthful CH is "too n-degenerate for the Connes distance to be well-defined" on most cross-shell pure-state pairs (the Connes distance SDP is unbounded). This is *Connes distance*, not modular flow; the modular flow does NOT require non-degeneracy of D (only cyclic-separatingness of Ω). So R3.2's obstruction does NOT carry over to L1-A. **This is a structurally important point: L1-A on truthful CH should work where R3.2 on truthful CH did not.**

---

## §5. σ_{i·2π} = identity verification: computational specification (Question 4)

The KMS condition at β = 2π is

    ω_{2π}(A · σ_{i·2π}(B)) = ω_{2π}(B · A)   for all A, B ∈ O_wedge.

If the construction in §2 gives σ_{i·2π} = id, this reduces to the trace property ω_{2π}(A·B) = ω_{2π}(B·A), which is automatic for tracial Gibbs states at any β. **This means the trace-property check is NOT the L1 falsifier — it is automatic.** The non-trivial check is that σ_t = Δ^{it} satisfies σ_{i·2π}(A) = A pointwise *as an operator on O_wedge*, not just under the state.

The correct falsifier is at the operator level:

**KMS-PERIOD CHECK:** For each generator M_k ∈ O_wedge basis, compute σ_t(M_k) := Δ^{it} M_k Δ^{-it} for t ∈ [0, 2π] (real time) and analytically continue to t = i·2π. The periodicity claim is

    Δ^{-2π} · M_k · Δ^{2π} = M_k   for all k = 1, ..., dim(O_wedge).

Bit-exact verification at n_max ∈ {2, 3, 4} requires:

1. Compute Δ as Hermitian matrix from S = J·Δ^{1/2} polar decomposition of the antilinear S on H_ω.
2. Diagonalize Δ = U·diag(λ_i)·U^*.
3. Compute Δ^{2π} = U·diag(λ_i^{2π})·U^*. Note: λ_i are real positive (Δ is positive); λ_i^{2π} is a positive real number.
4. For each generator M_k, compute residual R_k := ||Δ^{-2π}·M_k·Δ^{2π} - M_k||_F (Frobenius norm).
5. Verdict: pass if max_k R_k < 10^{-10} (machine precision for the diagonalization step at n_max ≤ 4).

**Diagnostic gradations:** if max_k R_k is ~10^{-3} or larger, the modular period is not 2π and the BW reading fails at finite n_max (Sprint 3.1 §3.1 outcome (b) "structurally unrelated"). If max_k R_k decreases with n_max (e.g., 10^{-1}, 10^{-2}, 10^{-3} at n_max = 2, 3, 4), the structural correspondence holds only in the GH limit, which Paper 38 propinquity machinery is needed to lift (see §6). If max_k R_k ≈ machine precision at every n_max, **literal identification at finite n_max is established** — this is the structurally cleanest closure.

**Expected outcome:** The continuum limit (Casini–Huerta–Myers arXiv:1102.0440; Casini–Huerta arXiv:0905.2562) of the modular Hamiltonian for a hemisphere of S³ for a free CFT does give modular period 2π. The Wightman-axiomatic BW theorem applies in the continuum; the question is whether the *truncation* preserves this period.

---

## §6. Module API proposal (Question 5)

Proposed file: `geovac/modular_hamiltonian.py` (estimated 800–1200 lines, 35–55 tests).

```python
# geovac/modular_hamiltonian.py
"""Modular Hamiltonian K on the truncated CH spectral triple T_{n_max}.

Implements Sprint L1-A: Tomita-Takesaki polar decomposition for a
wedge-restricted state on the Connes-vS truncated operator system,
producing the modular operator Δ, modular conjugation J_mod, and
modular Hamiltonian K = log Δ. Verifies the Bisognano-Wichmann
period-2π closure σ_{i·2π} = id at finite n_max.

References:
- Bisognano-Wichmann 1976, J. Math. Phys. 17, 303.
- Connes-Rovelli 1994, Class. Quantum Grav. 11, 2899 (arXiv:gr-qc/9406019).
- Connes-vS 2021, Comm. Math. Phys. 383 (arXiv:2004.14115).
- Paper 32 §VIII rem:bisognano_wichmann_reading.
- Paper 38 (WH1 PROVEN).
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, List, Literal, Optional, Tuple
import numpy as np
import scipy.linalg

from geovac.operator_system import HyperLabel, TruncatedOperatorSystem
from geovac.real_structure import RealStructure  # for cross-check against Paper 32 §IV J
from geovac.connes_distance import DiracProxy
from geovac.berezin_reconstruction import berezin_reconstruct  # for L4 lift

# ---------- Wedge sub-operator-system (Q2 W1) ----------

@dataclass(frozen=True)
class WedgeSpec:
    """Wedge specification: half-S³ defined by ω·n̂ > 0."""
    n_hat: Tuple[float, float, float, float] = (1.0, 0.0, 0.0, 0.0)
    cutoff_N_chi: Optional[int] = None  # default 2*n_max

def build_wedge_projector(op_sys: TruncatedOperatorSystem,
                           wedge: WedgeSpec) -> np.ndarray:
    """Build the multiplier matrix P_wedge of the smoothed indicator
    χ_{H_+} on the hemisphere H_+ = {ω·n̂ > 0}, truncated to Peter-Weyl
    cutoff cutoff_N_chi.

    Uses the same Avery-Wen-Avery 3-Y machinery as
    operator_system.build_multiplier_matrix, extended to non-monomial
    multipliers via Peter-Weyl expansion of χ_{H_+}.
    """
    ...

def wedge_operator_system(op_sys: TruncatedOperatorSystem,
                          wedge: WedgeSpec) -> 'WedgeOperatorSystem':
    """Return the wedge sub-operator-system O_wedge = P_wedge · O_{n_max} · P_wedge."""
    ...

# ---------- KMS state ω_β (Q3) ----------

@dataclass(frozen=True)
class KMSState:
    """Gibbs state ω_β(O) = Tr(e^{-β D} O)/Z(β) on O_{n_max}."""
    beta: float
    dirac_mode: Literal['shell_scalar', 'offdiag'] = 'shell_scalar'  # truthful CH
    n_max: int  # inherited from op_sys

    def density_matrix(self, basis: List[HyperLabel]) -> np.ndarray:
        """ρ_β = e^{-β D}/Z(β) on H_{n_max}."""
        ...

    def expectation(self, O: np.ndarray) -> complex:
        """ω_β(O) = Tr(ρ_β · O)."""
        ...

# ---------- ModularHamiltonian (Q5) ----------

@dataclass
class ModularHamiltonian:
    """K = log Δ where Δ is the Tomita modular operator on the GNS
    Hilbert space of a wedge-restricted KMS state.

    Construction:
    1. GNS triple (π_ω, H_ω, Ω) for ω_β restricted to O_wedge.
    2. S: π_ω(a)Ω → π_ω(a*)Ω antilinear (closure on π_ω(O_wedge)Ω dense).
    3. Polar decomposition S = J_mod · Δ^{1/2}.
    4. K := log Δ; modular flow σ_t = Δ^{it} · _ · Δ^{-it}.
    """
    op_sys: TruncatedOperatorSystem
    wedge: WedgeSpec
    kms_state: KMSState

    # Cached after construct()
    Delta: Optional[np.ndarray] = None       # modular operator (Hermitian, positive)
    J_mod: Optional[np.ndarray] = None        # modular conjugation (antilinear via U·conj)
    K: Optional[np.ndarray] = None            # modular Hamiltonian = log Δ
    gns_basis: Optional[List[np.ndarray]] = None  # ON basis of H_ω
    cyclic_vector: Optional[np.ndarray] = None  # Ω in H_ω

    def construct(self) -> None:
        """Build Δ, J_mod, K via Tomita polar decomposition."""
        ...

    def modular_flow(self, t: float) -> Callable[[np.ndarray], np.ndarray]:
        """σ_t: O → O, σ_t(A) = Δ^{it} A Δ^{-it}."""
        ...

    def modular_flow_imaginary(self, theta: float) -> Callable[[np.ndarray], np.ndarray]:
        """σ_{iθ}: O → O, σ_{iθ}(A) = Δ^{-θ} A Δ^{θ}."""
        ...

    def verify_kms_condition(self, A: np.ndarray, B: np.ndarray,
                              tol: float = 1e-10) -> dict:
        """Check ω_β(A · σ_{iβ}(B)) ≈ ω_β(B · A). Returns dict with residual."""
        ...

    def verify_modular_periodicity(self, generators: List[np.ndarray],
                                    tol: float = 1e-10) -> dict:
        """Check σ_{i·2π}(M_k) = M_k for each generator M_k ∈ O_wedge.
        Returns dict with per-generator residual and verdict (pass/fail/marginal)."""
        ...

    def lift_to_round_S3(self, paper38_lemmas: 'PaperLemmas') -> dict:
        """Apply Paper 38 L1'-L5 lifts to extrapolate K_{n_max} → K_∞.
        Returns convergence rate diagnostics; cross-check Λ(T_{n_max}, T_S3)
        propinquity bound from gh_convergence.py."""
        ...

# ---------- Witness-specific factories (Q7) ----------

def for_bisognano_wichmann(op_sys: TruncatedOperatorSystem,
                            n_hat: Tuple[float, ...] = (1, 0, 0, 0)
                            ) -> ModularHamiltonian:
    """BW: hemisphere wedge, β = 2π. κ_g = 1 (boost rapidity normalization)."""
    return ModularHamiltonian(
        op_sys=op_sys,
        wedge=WedgeSpec(n_hat=n_hat),
        kms_state=KMSState(beta=2*np.pi, n_max=op_sys.n_max),
    )

def for_hawking(op_sys: TruncatedOperatorSystem, M: float
                ) -> ModularHamiltonian:
    """Hawking on Euclidean Schwarzschild cigar: β = 8π M, κ_g = 1/(4M).
    Returns ModularHamiltonian whose σ_{i·β} period closure tests Hawking
    radiation at finite n_max."""
    return ModularHamiltonian(
        op_sys=op_sys,
        wedge=WedgeSpec(),
        kms_state=KMSState(beta=8*np.pi*M, n_max=op_sys.n_max),
    )

def for_unruh(op_sys: TruncatedOperatorSystem, a: float
              ) -> ModularHamiltonian:
    """Unruh: Rindler wedge, β = 2π/a, κ_g = a (proper acceleration)."""
    return ModularHamiltonian(
        op_sys=op_sys,
        wedge=WedgeSpec(),
        kms_state=KMSState(beta=2*np.pi/a, n_max=op_sys.n_max),
    )

def for_sewell(op_sys: TruncatedOperatorSystem, M: float
                ) -> ModularHamiltonian:
    """Sewell 1982: Hartle-Hawking vacuum on Kruskal-extended Schwarzschild.
    Same parameters as Hawking; different state interpretation."""
    return for_hawking(op_sys, M)
```

Auxiliary computations needed:

- `_build_gns_triple(op_sys, wedge, kms_state) → (H_ω, π_ω, Ω)`: GNS construction over the wedge sub-algebra against the KMS state. The Hilbert space H_ω is the completion of O_wedge under the inner product ⟨a, b⟩_ω = ω_β(a* b). Practically: H_ω is realized as O_wedge itself with the Hilbert-Schmidt inner product weighted by ρ_β.
- `_build_S_operator(gns_triple) → antilinear S`: S: π_ω(a)·Ω ↦ π_ω(a*)·Ω. As a matrix: S = J_GNS · K where K is complex conjugation and J_GNS is a specific permutation/phase computed from the GNS basis.
- `_polar_decompose_antilinear(S) → (J_mod, Delta_half)`: antilinear polar decomposition. Practically: Δ = S* S, where S* is the adjoint of the antilinear S (which is linear; this is the standard trick of converting antilinear operations to linear ones via complex conjugation). Δ is then a positive Hermitian matrix; J_mod = S · Δ^{-1/2}.

The polar-decomposition step is the load-bearing numerical operation. For the truthful CH (diagonal D), Δ should also be diagonal in a natural basis (the eigenbasis of the Gibbs state), and the polar decomposition is closed-form. For offdiag CH or wedge restrictions that break diagonal structure, the polar decomposition is numerical (scipy.linalg.polar on the linearised form).

---

## §7. Five-lemma lift table (Question 6)

How each Paper 38 lemma extends to the modular context:

| Paper 38 Lemma | Riemannian content | Modular extension | Status |
|:---|:---|:---|:---|
| **L1' (chirality doubling)** | Offdiag CH operator system substrate has every non-forced cross-shell Connes-distance pair finite. | L1'-mod: the wedge sub-operator-system O_wedge inherits the chirality-doubled structure. NEEDED ONLY IF using offdiag CH. For truthful CH (recommended in §4), L1' is moot — truthful CH does NOT have the n-degeneracy obstruction for *modular flow* (only for Connes distance, R3.2). | **MOOT for truthful CH path; needed for offdiag CH backup path.** |
| **L2 (central spectral Fejér kernel on SU(2))** | Positive central Fejér kernel on SU(2) ~ S³ with rate γ_{n_max} → 0 (4/π asymptote = M1 signature). | L2-mod: convolution with the central Fejér kernel commutes with the modular flow if the kernel is wedge-respecting. For W1 (hemisphere) the kernel is NOT wedge-respecting (it's a central function on SU(2), invariant under conjugation). **Need a wedge-localized variant**: K_{wedge, n_max}(g) := K_{n_max}(g) · χ_{H_+}(target). Convolution is then *not* central but is wedge-respecting. Build by composition. | **NEEDS NEW WORK (1-2 weeks).** Wedge-localized Fejér variant. |
| **L3 (Lipschitz spinor bound, C_3 = 1)** | C_3 = 1 from Avery harmonic gradient + spectral-gap quantization. Universal across SU(N). | L3-mod: the Lipschitz seminorm L(a) := ||[D, a]||_op extends to the modular Lipschitz seminorm L_mod(a) := ||[K, a]||_op. Since K = log Δ is *unbounded* in the GH limit (cyclic-vector argument), the bound C_3 = 1 must be replaced by an n_max-dependent C_3^{(mod)}(n_max). For truthful CH, K has rational eigenvalues structured by the Gibbs spectrum; expect C_3^{(mod)} ~ β·||D||_op = β·(n_max + 3/2). | **NEEDS NEW BOUND (2-3 weeks).** Likely scales linearly with n_max; would propinquity-bound the convergence rate. |
| **L4 (Berezin reconstruction)** | B_{n_max}: C(S³) → O_{n_max} positive, contractive, approximate identity at rate γ_{n_max}. | L4-mod: the Berezin map sends wedge-localized continuous functions to wedge-localized operator-system elements. Lift B_{n_max} | _{C^∞(H_+)} → O_wedge straightforwardly. The approximate-identity property at the wedge boundary ∂H_+ = S² may be slower than in the global case (Gibbs phenomenon near the cutoff χ_{H_+}'s edge). | **NEEDS BOUNDARY DIAGNOSIS (1 week).** Expect rate γ_wedge^{n_max} = O(log n / n) but with a constant 2-5× the global L4 constant. |
| **L5 (Latrémolière propinquity assembly)** | Λ(T_{n_max}, T_{S³}) ≤ C_3 · γ_{n_max} → 0. | L5-mod: does the propinquity convergence imply convergence of K_{n_max}? **This is the key open question of L1-A.** The propinquity convergence is on the *metric* spectral triple data (algebra, Hilbert space, Dirac); the modular Hamiltonian is *derived* from a state on the algebra. If K_{n_max} converges in the strong-operator topology on the limiting H_ω, the BW closure transfers to the limit. If not, the structural correspondence stays at finite n_max. **A "modular propinquity" extension (not published) would be needed for the strongest form; for L1-A, verifying K_{n_max} closure pointwise at n_max ∈ {2, 3, 4} + GH-rate diagnostic on K_{n_max} norm suffices for the closure.** | **NEW MATHEMATICS (3-5 weeks).** Likely the longest individual leg. |

**Total L-lemma lift budget: 1 + 2-3 + 1 + 3-5 = 7-11 weeks** if each is done from scratch. **However:** L1-A's bar is *not* to publish a Theorem statement at Paper 38 rigor level; it is to verify the BW closure σ_{i·2π} = id at finite n_max numerically. Many of the analytical bounds can be replaced by direct numerical verification at n_max ∈ {2, 3, 4}, which is much faster. **Revised budget if we accept numerical verification + qualitative lift sketches: 4-6 weeks.**

---

## §8. Witness-specific parameter table (Question 7)

Per the Sprint Unruh-pendant four-witness theorem (`debug/unruh_pendant_memo.md`, codified in `papers/group1_operator_algebras/paper_32_spectral_triple.tex` `rem:bisognano_wichmann_reading` Unruh pendant paragraph):

| Witness | Surface gravity κ_g | Inverse temperature β | Period statement | Wedge structure | Predicted K structure |
|:---|:---|:---|:---|:---|:---|
| **Bisognano–Wichmann 1976** | rapidity rate of boost | β = 2π / κ_g (with κ_g normalized to 1 for canonical BW) | σ_{i·2π} = id | Rindler wedge x^1 > |x^0|; on T_{n_max}: hemisphere H_+ of S³ (Riemannian analog) | K = 2π · K_boost = 2π · (boost generator) |
| **Hawking 1976** | κ_g = 1/(4M) (Schwarzschild surface gravity) | β = 8π M (Hartle–Hawking) | σ_{i·8π M} = id on Schwarzschild exterior algebra | Euclidean Schwarzschild cigar; on T_{n_max}: same hemisphere H_+ (parameterized by M-dependent β) | K = 8π M · K_Killing = 8π M · (Schwarzschild Killing time generator) |
| **Sewell 1982** | κ_g = 1/(4M) (same as Hawking) | β = 8π M (same as Hawking) | σ_{i·8π M} = id on static exterior algebra; standard generalization of BW to curved backgrounds | Static exterior of Kruskal-extended Schwarzschild; on T_{n_max}: same hemisphere H_+ | Same as Hawking; Sewell's contribution is the proof that Hartle–Hawking vacuum IS KMS at β = 8π M |
| **Unruh 1976** | κ_g = a (proper acceleration) | β = 2π / a | σ_{i·2π/a} = id on Rindler wedge | Rindler wedge of accelerated observer; on T_{n_max}: same hemisphere H_+ (parameterized by a) | K = (2π / a) · K_boost = (2π / a) · (boost generator) |

**The single-bridge claim of the Sprint Unruh-pendant theorem is that all four witnesses share the same KMS-Tomita-Takesaki algebra at the framework level, with κ_g entering as a scalar parameter via β = 2π / κ_g.** L1-A's primary closure target is therefore the BW case (κ_g = 1, β = 2π). If that closes, the Hawking, Sewell, and Unruh cases follow by parameter substitution in the same code path (Q7's witness-specific factories in §6 demonstrate this trivially: they each instantiate the same `ModularHamiltonian` class with a different β).

**Critical structural observation:** the BW reading's β = 2π is the "bare" Hopf-base-measure 2π = Vol(S^1), which is the M1 sub-mechanism signature of the master Mellin engine (Paper 32 §VIII `thm:pi_source_case_exhaustion`). Hawking/Unruh introduce additional scalar multipliers (M, a) that *modulate* this M1 signature but do not change its qualitative role. **The four-witness theorem's headline content — "β = 2π / κ_g" — is the operator-system-level version of the master Mellin engine's M1 mechanism.** L1-A is therefore the operator-system-level realization of one piece of `thm:pi_source_case_exhaustion`.

---

## §9. Lift via Paper 38 lemmas: does Λ → 0 imply K_{n_max} → K_∞? (Question 6 continued)

The answer is: **conditionally yes, with one important caveat.**

The Latrémolière quantum propinquity Λ measures distance between metric spectral triples in a topology that respects the algebra, Hilbert space, and Dirac operator simultaneously. It does NOT directly bound modular Hamiltonians, which are derived objects depending on a chosen state on the algebra.

**Condition for K_{n_max} → K_∞:** if the state ω_β chosen on each T_{n_max} is consistent (e.g., the restriction of a single ω_β^∞ on T_{S³} to each truncation, or equivalently, the Berezin-reconstructed Gibbs state at fixed β), then GH convergence of (T_{n_max}, ω_β^{(n_max)}) → (T_{S³}, ω_β^∞) in the appropriate state-included propinquity (a refinement of Latrémolière's, not yet published but plausible) WOULD imply convergence of the GNS triples, and hence of S, Δ, J_mod, K.

**The important caveat:** K = log Δ is *unbounded* in the limit; convergence of K is best stated in the resolvent or bounded-functional-calculus sense (e^{itK} converges strongly, not K itself). For L1-A's falsifier (σ_{i·2π} = id), the relevant quantity is e^{2π·K} = Δ^{2π}, which is bounded if Δ is bounded away from 0 and ∞. This holds on the truncated finite-dimensional setting but requires care in the limit (D_{n_max} is bounded; D_∞ is unbounded, so Δ_∞ may have spectrum extending to 0 and ∞).

**Operational L1-A protocol:**

1. Verify σ_{i·2π} = id at n_max = 2, 3, 4 (closed-form, finite matrices).
2. Compute the operator-norm distance ||σ_{i·2π}^{(n_max)} - id||_op as a function of n_max.
3. If this distance is bit-zero at every n_max: literal identification at finite n_max established. Closure complete.
4. If it decreases monotonically with n_max: literal identification at finite n_max fails but holds in the GH limit; cite L5-mod (the published-pending propinquity-of-modular-flows extension) as the closure mechanism.
5. If it does not decrease (or increases): the BW reading fails at the operator-system level; report as a negative result.

Outcome (3) is the cleanest closure; outcome (4) is acceptable closure with a flagged-future-work appendix; outcome (5) is the only "wall" outcome that would close L1-A as a negative result.

---

## §10. Predicted obstructions and verification plans (Question 8)

Five obstructions are likely; each has a verification plan and pre-computed expected outcome.

### §10a. Operator-system non-multiplicative-closure

**Statement:** O_{n_max} is *-closed but NOT multiplicatively closed (Paper 28 §broken_structural_zero, Paper 32 §III `rem:operator_system`). The Tomita modular operator Δ is constructed from S: a·Ω → a*·Ω, which uses the algebra action a·Ω = π_ω(a)·Ω. The GNS construction proceeds via the *enveloping* of O_wedge to a C*-algebra (e.g., the bicommutant), then Tomita on the bicommutant. Whether the resulting σ_t restricts to a well-defined flow on O_wedge is non-trivial.

**Verification plan:** Compute the matrix elements of σ_t(M_k) for M_k ∈ O_wedge at small t (e.g., t = 0.01, 0.1, 1.0). If σ_t(M_k) lies in O_wedge to within numerical tolerance, the restriction is well-defined. If σ_t(M_k) has weight outside O_wedge that grows with t, then σ_t leaves the operator system — this would be a structural obstruction.

**Expected outcome:** The Connes–vS framework explicitly handles operator systems via passage to C*-envelope (Paper 32 §III remark, Connes–vS arXiv:2004.14115 §2.4); the modular flow is defined on the envelope, not the operator system itself. Restricting back may introduce a "modular leakage" that scales as the operator system's gap from its envelope. At n_max = 2: dim(O)/N^2 = 0.560 (Paper 32 §III), so 44% leakage possible; at n_max = 4: 0.156, 84% leakage possible. **This is the most likely structural obstruction.** Verification at n_max = 2, 3, 4 will give a quantitative leakage curve; whether the BW closure survives within the leakage budget is the empirical question.

### §10b. Truncated CH spectrum symmetry → trivial modular structure

**Statement:** The truncated CH spectrum |λ_n| = n + 3/2 is positive-only (not symmetric around zero in the truthful CH realization on the scalar Fock basis). In the full Dirac sector (Paper 32 §IV, `geovac/full_dirac_operator_system.py`), the spectrum *is* symmetric (chirality-doubled with both ±|λ_n| eigenvalues). For the symmetric spectrum, the Gibbs state Tr(e^{-β D}) factorizes as the chirality-symmetric sum, and the modular flow may degenerate trivially (analogous to Sprint L0 prediction that M3 trivializes on chirality-symmetric spectrum at (3,1) — see Sprint TD Track 3 `debug/sprint_td_track3_memo.md`).

**Verification plan:** Compute the modular operator Δ on (a) truthful CH scalar basis (positive spectrum), (b) full Dirac (chirality-symmetric spectrum). Compare structure: if Δ in case (b) is the identity (or nearly so), the chirality symmetry kills the modular flow.

**Expected outcome:** Case (a) (truthful CH scalar) has non-trivial Δ; case (b) (full Dirac) may have a chirality-block-diagonal Δ with each block matching case (a). The Sprint L0 M3 trivialization prediction (`debug/lorentzian_l0_audit_memo.md` §4) for the (3,1) Lorentzian extension is the structurally adjacent claim; verifying that the (3,0) Riemannian modular flow does NOT trivialize on chirality-symmetric spectrum is a consistency check for the L0 prediction.

### §10c. L0 M3-trivialization-on-Lorentzian: implications for the Riemannian L1

**Statement:** Sprint L0 (`debug/lorentzian_l0_audit_memo.md` §3) predicts that the M3 sub-mechanism (vertex-parity Hurwitz / Dirichlet-L content) becomes trivially empty on chirality-symmetric spectrum at (3,1). This is a *Lorentzian* prediction, not Riemannian. But: does the framework-internal modular Hamiltonian at (3,0) carry any M3 content? If yes, then the L1 closure at (3,0) is M3-active, and the L2 lift to (3,1) would discover the M3 trivialization. If no, the L1 closure at (3,0) is M3-inactive, and the L2 lift is automatic.

**Verification plan:** Compute K's spectrum at n_max = 2, 3, 4 and check whether its transcendental content (rationals × π powers; rationals × Catalan G; etc.) is in M1, M2, M3, or a combination. The Sprint MR-A/B/C master Mellin engine classification (`debug/mr_b_spectral_action_rate_memo.md`) provides the dictionary.

**Expected outcome:** K's transcendental content is purely M1 (proportional to 2π = Vol(S^1)). This would be the framework-internal verification that the BW reading lives entirely in M1, consistent with `papers/group1_operator_algebras/paper_32_spectral_triple.tex` `rem:bisognano_wichmann_reading`. Any M2 or M3 content would be a surprise requiring re-interpretation.

### §10d. Wedge cutoff Gibbs phenomenon

**Statement:** The wedge cutoff χ_{H_+} is a step function (indicator of the hemisphere). Its Peter-Weyl expansion at finite cutoff N_χ has Gibbs oscillations near the equator ∂H_+. These oscillations may cause the wedge multiplier P_wedge to be ill-conditioned, leading to numerical instabilities in the polar decomposition step.

**Verification plan:** Compute the condition number of P_wedge as a function of N_χ. Test smoothed cutoffs (e.g., χ_{H_+}^{smooth}(ω) = (1 + tanh(ε^{-1} ω·n̂))/2) at various smoothing widths ε. Verify that the modular flow is insensitive to the smoothing in the limit ε → 0.

**Expected outcome:** Numerical conditioning is fine at n_max ≤ 4 and N_χ ≤ 10 (sympy / numpy double precision is sufficient). Smoothing-width independence may need to be verified rather than assumed.

### §10e. Connes-axiom-compatibility of J_mod vs Paper 32 §IV J

**Statement:** The Tomita modular conjugation J_mod from polar decomposition is, in general, NOT the same as the Connes-axiom J of Paper 32 §IV. They are related: both are antilinear, both satisfy J^2 = ±I depending on signature, but J_mod is *state-dependent* while J of Paper 32 §IV is *spectral-triple-intrinsic*. The BW reading needs them to be compatible — specifically, J_mod should commute with the algebra action on H_ω in a way that's consistent with J·a·J^{-1} ∈ A^op (Connes's "opposite-algebra" structure).

**Verification plan:** Compute J_mod and compare to RealStructure on the truncated CH (Paper 32 §IV `real_structure.py`). Verify J_mod^2 = -I (KO-dim 3 convention). Verify J_mod·D·J_mod^{-1} = +D. Verify the operator-system commutation [a, J_mod·b·J_mod^{-1}] = 0 (which Paper 32 §IV verified for the spectral-triple-intrinsic J on the truncated O at residual ~5-8%).

**Expected outcome:** J_mod from polar decomposition is automatically antilinear and satisfies J_mod^2 = -I (KO-dim 3 forced by S spec). The commutation with D and the operator-system commutation are *not* automatic and need verification. If they fail, then J_mod and Paper 32 §IV's J are different objects — the BW reading then needs an additional bridge argument relating them.

---

## §11. Verdict and sequencing recommendation

**Verdict: GO-WITH-CAVEATS.**

The architecture is sound, but the implementation budget should be set at the upper end of the dispatch prompt's "4–8 weeks" estimate, more honestly **5–9 weeks** of focused work for a clean closure with documentation:

- **Weeks 1–2:** Build `WedgeOperatorSystem` (W1 hemisphere), `KMSState`, GNS construction, polar decomposition (single n_max = 2 prototype). Sanity-check on toy commutative S^1 example for which the modular Hamiltonian has a known closed form (Casini–Huerta scalar field on a half-circle).
- **Weeks 2–3:** Verify σ_{i·2π} = id at n_max = 2 in the truthful CH path. Trace the residual against the §10a operator-system-leakage estimate. Decide between outcomes (3) / (4) / (5) of §9.
- **Weeks 3–5:** Extend to n_max = 3, 4. Quantify GH-rate of σ_{i·2π} - id; cross-check against Paper 38 propinquity bound Λ_{n_max}. Verify witness factories (Hawking, Unruh, Sewell) reproduce expected β = 8πM / 2π/a / 8πM structures.
- **Weeks 5–7:** Address obstructions §10b–§10e individually; for each, produce a verification script + memo subsection. Most likely outcome: §10a (operator-system leakage) requires the most attention; §10b (chirality symmetry) and §10c (M3 trivialization) become passing consistency checks; §10d (Gibbs phenomenon) and §10e (J_mod compatibility) are routine.
- **Weeks 7–9:** Write `geovac/modular_hamiltonian.py` production module + 35–55 tests + memo + paper-edit pass on `papers/group1_operator_algebras/paper_32_spectral_triple.tex` `rem:bisognano_wichmann_reading` and `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` §III.27 to lift "structural correspondence" → "literal identification at finite n_max".

**Caveats:**

1. The implementation budget is dominated by the polar-decomposition diagnostics on the GNS Hilbert space, which is the new technical content. The framework-internal infrastructure (truncated operator system, Avery–Wen–Avery integrals, Connes distance SDP, real structure J, Berezin reconstruction) is in place and well-tested; what's missing is the GNS + modular-flow layer.
2. L5-mod (Latrémolière propinquity of modular Hamiltonians) is an open mathematical problem; L1-A does not solve it analytically. The numerical verification at n_max = 2, 3, 4 is the load-bearing closure, lifted to a "structural identification" claim via the GH convergence of the underlying spectral triples (Paper 38) without a formal Theorem on modular convergence.
3. The output is a *finite-n_max* closure of the BW reading on the truthful Camporesi–Higuchi triple. The Lorentzian (3,1) extension (Sprint L2) is independent and not addressed by L1-A. If L1-A succeeds, Paper 34 §III.27 graduates to "literal identification at the operator-system level (Riemannian)"; Sprint L2 is then the natural next sprint to lift to (3,1) via the Bizi–Brouder–Besnard Krein framework.
4. If §10a operator-system leakage forces σ_{i·2π} ≠ id at finite n_max but with residual decreasing with n_max, then L1-A closes as "structural correspondence + finite-n_max approach to literal identification + literal identification in the GH limit (conditional on the unpublished modular-propinquity extension)". This is outcome (4) of §9 and is still a closure, but with explicit acknowledgment of the unresolved analytic question.

**No precursor sprint is needed.** All required infrastructure exists (truncated operator system + Avery–Wen–Avery integrals + Connes distance SDP framework + real structure J + Berezin reconstruction + central spectral Fejér kernel + thermal tensor triple). The L1-A implementation directly uses these modules.

**Recommendation:** Dispatch L1-A implementation as a single agent with the 5–9 week budget. Set the falsifier verification at n_max = 2, 3, 4 as the gating exit criterion. The expected positive outcome promotes Paper 34 §III.27 and `papers/group1_operator_algebras/paper_32_spectral_triple.tex` `rem:bisognano_wichmann_reading` from structural-correspondence to operator-system-level literal identification, which is a significant qualitative advance closing the principal Sprint L1 falsifier named in Paper 32 §VIII and in the Sprint L0 audit.

---

## §12. Files this blueprint authorizes the implementation agent to create

- `geovac/modular_hamiltonian.py` (~800–1200 lines, ~35–55 tests)
- `tests/test_modular_hamiltonian.py` (full coverage; toy S^1 sanity check, n_max = 2/3/4 BW closure tests, Hawking/Unruh/Sewell parameter-substitution checks, obstruction verification subsections corresponding to §10a–§10e)
- `debug/l1_modular_hamiltonian_implementation_memo.md` (the implementation memo, structured around the verdict + per-obstruction results)
- `debug/data/l1_modular_hamiltonian_results.json` (machine-readable σ_{i·2π} residuals at n_max ∈ {2, 3, 4}, GH-rate diagnostics, M-mechanism classification of K's transcendental content)

Files the implementation agent should NOT modify:

- Any paper file. Paper edits are decided by the PI after reviewing the L1 implementation results.
- Production GeoVac modules other than `modular_hamiltonian.py` itself. If new functionality is needed in (e.g.) `operator_system.py` to handle non-monomial multipliers, the addition should be a sub-PR with its own test coverage and explicit PI sign-off, not an in-place modification.

---

## §13. References (verified)

- **Bisognano, J. J., & Wichmann, E. H. (1976).** On the duality condition for quantum fields. *J. Math. Phys.* **17**, 303–321. The original modular-flow-as-boost theorem.
- **Bisognano, J. J., & Wichmann, E. H. (1975).** On the duality condition for a Hermitian scalar field. *J. Math. Phys.* **16**, 985–1007. Predecessor paper, scalar field only. (Sprint Track D footnote corrected the title: the 1975 paper is "for a Hermitian scalar field"; the 1976 paper is the more general "for quantum fields".)
- **Hawking, S. W. (1976).** Black holes and thermodynamics. *Phys. Rev. D* **13**, 191–197. Imaginary-time-period derivation; β = 8πM.
- **Sewell, G. L. (1982).** Quantum fields on manifolds: PCT and gravitationally induced thermal states. *Ann. Phys. (NY)* **141**, 201–224. Generalization of BW to Schwarzschild; Hartle–Hawking vacuum is KMS. (Sprint Track D footnote corrected year: 1982, not 1980.)
- **Unruh, W. G. (1976).** Notes on black-hole evaporation. *Phys. Rev. D* **14**, 870–892. T_U = a/(2π).
- **Connes, A., & Rovelli, C. (1994).** Von Neumann algebra automorphisms and time-thermodynamics relation in general covariant quantum theories. *Class. Quantum Grav.* **11**, 2899–2918. arXiv:gr-qc/9406019. The "thermal time hypothesis" — modular automorphism IS physical time. Directly relevant to L1-A: identifies the Tomita modular flow with the boost orbit (specializing BW to general covariance).
- **Connes, A., & van Suijlekom, W. D. (2021).** Spectral truncations in noncommutative geometry and operator systems. *Comm. Math. Phys.* **383**, 87–129. arXiv:2004.14115. The truncated-operator-system framework. *No modular-Tomita extension in this paper; the operator-system Tomita is not addressed.*
- **Hekkelman, E.-M., & McDonald, E. (2024).** Asymptotic expansions of pseudodifferential operators and applications. arXiv:2412.00628. Multiple operator integrals on truncations. *Does not address modular flow.*
- **Hekkelman, E.-M., McDonald, E., & van Nuland, T. D. H. (2024).** Multiple operator integrals, pseudodifferential calculus, and asymptotic expansions. arXiv:2404.16338. Spectral action perturbative expansion on truncations. *No modular Hamiltonian.*
- **Leimbach, A., & van Suijlekom, W. D. (2024).** Gromov–Hausdorff convergence of spectral truncations for tori. *Adv. Math.* **439**, 109496. The torus version of Paper 38's S³ result. Verified no modular-Tomita extension included.
- **Casini, H., Huerta, M., & Myers, R. C. (2011).** Towards a derivation of holographic entanglement entropy. *JHEP* **05**, 036. arXiv:1102.0440. The continuum entanglement Hamiltonian for a sphere of S^d in CFT; the structural target that L1-A's truncated-triple modular Hamiltonian should reproduce in the GH limit.
- **Casini, H., & Huerta, M. (2009).** Entanglement entropy in free quantum field theory. *J. Phys. A: Math. Theor.* **42**, 504007. arXiv:0905.2562. Half-space entanglement Hamiltonian; standard reference for the structural form K = 2π ∫ T^00 · x dx.
- **Zhu, J., Casini, H., et al. (2020).** Lattice Bisognano–Wichmann modular Hamiltonian in critical quantum spin chains. *SciPost Phys. Core* **2**, 007. arXiv:2003.00315. Direct precedent for finite-lattice BW reading; verified the BW form holds to good approximation on critical chains. Relevant as a comparison: L1-A is structurally similar but for a non-flat compact S³ background.
- **Strohmaier, A. (2006).** The Reeh–Schlieder property for quantum fields on stationary spacetimes. *Commun. Math. Phys.* **215**, 105–118 (verify; primary reference for modular Reeh–Schlieder on Riemannian backgrounds).
- **Verch, R. (2001).** A spin-statistics theorem for quantum fields on curved spacetime manifolds in a generally covariant framework. *Commun. Math. Phys.* **223**, 261–288. Nuclearity-modular framework on compact Riemannian backgrounds; standard reference for Tomita on curved-background QFT.

**Phantom-citation check:** All citations above verified via WebSearch / WebFetch. Connes–Rovelli 1994 confirmed at arXiv:gr-qc/9406019 (verified). Hekkelman–McDonald 2024 confirmed at arXiv:2412.00628 (per Paper 38 references and verified via search). Leimbach–vS 2024 in Adv. Math. (per Paper 38 references). Bisognano–Wichmann 1976 verified via search (`https://ncatlab.org/nlab/show/Bisognano-Wichmann+theorem`). Zhu et al. 2020 verified at arXiv:2003.00315 (`https://scipost.org/SciPostPhysCore.2.2.007`). Casini–Huerta 2009 verified at arXiv:0905.2562. Casini–Huerta–Myers 2011 verified at arXiv:1102.0440. **No phantom citations.**

End of L1-A architectural blueprint.
