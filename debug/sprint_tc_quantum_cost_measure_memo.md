# Sprint memo — xTC quantum-algorithm cost of the non-Hermitian effective operator

**Date:** 2026-08-23  **Verdict: GO — and the escape hatch moots the question.**

**Question.** The classical xTC PoC (`sprint_xtc_{poc_li,pblock}_memo.md`) diagonalized
the effective 2-body operator with `scipy.linalg.eig` but never measured the operator
properties that drive quantum-algorithm cost for a *non-Hermitian* object. Does the
xTC sparsity/1-norm win (Pauli 279→279, 2-body 1-norm 0.80–0.89×) survive once the
operator must be run on a QC as a non-Hermitian object? Measured on the ACTUAL
operators, engines reused READ-ONLY.

Driver: `debug/xtc_quantum_cost_measure.py` → `debug/data/xtc_quantum_cost.json`.
Systems: **He** (s-only, 2-body-only control, D+K), **Li 1s²2s** (s-only, genuine
3-body → D+K+xTC-L3, the load-bearing system), **C 1s²2s²2p²** (s+p, open-p ref,
D+xTC-L3). γ=1. PN-projected FCI throughout (never qubit-space diag).

---

## The pivotal structural finding (established first, by probe)

**The xTC 3-body→2-body contraction v2 is Hermitian to machine precision — for BOTH an
s-reference and a p-reference. The ONLY non-Hermitian term is the generic 2-body
convective K.**

| tensor | Hermiticity deviation max‖A−A†‖ |
|---|---|
| D (asym_w, finite TC kernel) | 3.3e-16 (Hermitian) |
| xTC-contracted L3, v2 — Li **s**-ref | 3.5e-18 (Hermitian) |
| xTC-contracted L3, v2 — C **p**-ref | 2.8e-17 (Hermitian) |
| convective K (asym_K) | **0.226 (non-Hermitian)** |

So "xTC produces a non-Hermitian object" is **misattributed**. The 3-body→2-body
contraction — the xTC machinery itself — preserves Hermiticity exactly. The
non-Hermiticity is the ordinary transcorrelated convective term
K = −(∇τ)·∇, which is present in **any** TC Hamiltonian, xTC or not. K's matrix
elements are real (float64) → the LCU coefficients stay real up to the JW repackaging
(small imaginary Pauli coefficients ~0.05 appear after JW, accounted for in λ=Σ|c|).

---

## Numbers table (primary variants bold; γ=1)

| system | variant | spectrum real? | non-normality | κ_V(low) | κ_V(full) | max Bauer–Fike κᵢ | gap (Ha) | λ_nH/λ_H | n_pauli | sym-win retained |
|---|---|---|---|---|---|---|---|---|---|---|
| **He** (2e ctrl) | plain (Coulomb) | Y | 6e-17 | 2.12 | 2.12 | — | 0.907 | 1.000 | 117 | — |
| **He** | **D+K (2body TC)** | **Y** | **0.069** | **1.98** | **1.99** | **1.22** | **0.923** | **1.096** | **249** | **1.59** |
| He | sym(D+K) [Herm] | Y | 4e-17 | 1.58 | 2.59 | — | 0.933 | **0.977** | **117** | — |
| **Li** (3e) | plain (Coulomb) | Y | 3e-17 | 1.30 | 1.30 | — | 1.472 | 1.000 | 117 | — |
| Li | D+K (2body TC) | Y | 0.015 | 1.25 | 1.49 | — | 1.414 | 1.131 | 249 | — |
| Li | D+xTC-L3 (noK, Herm) | Y | 2e-17 | 1.00 | 1.00 | — | 1.116 | 0.852 | 117 | — |
| **Li** | **D+K+xTC-L3 (FULL)** | **Y** | **0.015** | **1.21** | **1.45** | **1.03** | **1.397** | **1.122** | **249** | **1.77** |
| Li | sym(FULL) [Herm] | Y | 2e-17 | 1.00 | 1.00 | — | 1.423 | **0.978** | **117** | — |
| **C** (p-ref, s+p) | plain (Coulomb) | Y | 3e-18 | 1.12 | 4.30 | — | 0.076 | 1.000 | 279 | — |
| **C** | **D+xTC-L3 (noK)** | **Y** | **3e-18** | **1.00** | **1.00** | **1.00** | **0.022** | **1.030** | **279** | **1.00** |

- κ_V(low) = cond of the n×4 low-lying right-eigenvector block; Bauer–Fike κᵢ = 1/|⟨uᵢ|vᵢ⟩|.
- gap = physical (degeneracy-merged, 2 mHa tol) ground→first-excited; spurious ground-manifold split (spin-broken single-det ref) ≤1.2 mHa, reported separately.
- λ = full JW LCU 1-norm (1-body+2-body, complex coeff allowed, identity dropped).
- sym-win retained = (E_plain−E_sym)/(E_plain−E_TC); >1 ⇒ symmetrization over-recovers.

---

## Reading the numbers

**(1) Spectrum.** Every operator has an entirely real spectrum (pseudo-Hermitian /
PT-symmetric-like) and a real ground state — im(E₀)=0 in all cases. The non-Hermitian
K does **not** produce complex-conjugate eigenvalue pairs at these parameters.

**(2) Non-normality** is small: 0.069 (He, the pure 2-body TC), 0.015 (Li full). C is
Hermitian (3e-18) because the p-engine has no K.

**(3) κ_V — the load-bearing number — is O(1).** κ_V(low-lying) ≤ 1.98 (He), 1.21 (Li);
per-eigenvalue Bauer–Fike κᵢ ≤ 1.22 (He), 1.03 (Li). κ_V(full) ≤ 2.0. The eigenvector
conditioning that would drive QEVE/Bauer–Fike cost is **benign** even for the genuine
non-Hermitian operator. → **GO** (O(1–10)).

**(4) Gap.** He 0.92 Ha, Li 1.40 Ha, C 22 mHa (term-to-term). The ~0.1–1.2 mHa
"micro-gap" seen in a naïve read is a **spurious split of the exactly-degenerate ground
spin/term manifold** by the spin-broken single-determinant xTC reference (Li ground
doublet split 0.15 mHa; C ³P 9-fold manifold split 1.2 mHa) — a known reference-quality
limitation, not a QC-cost feature.

**(5) LCU 1-norm.** λ_nonHerm/λ_Herm = 1.096 (He), 1.122 (Li), 1.030 (C) — all ≲1.13,
well under the 1.5× STOP. K inflates λ modestly (Li: 0.852 no-K → 1.122 with K, ×1.32),
but the more consequential effect of K is on the **Pauli COUNT**, which roughly
**doubles** (Li/He 117 → 249): the non-symmetric K tensor fills JW Pauli strings the
symmetric Coulomb/w/L3 kernels leave empty. So the pblock memo's "279→279 unchanged"
is a property of the **Hermitian** D+xTC-L3 sub-operator; the *full* non-Hermitian TC
operator you would actually run carries ~2× the Pauli terms.
*(Cross-check: my C D+xTC-L3 2-body-tensor 1-norm ratio = **0.841**, an exact match to
the pblock memo's reported 0.841 — the reconstruction is faithful to the engine.)*

**K is essential — you cannot avoid it by deletion.** Dropping K entirely (Hermitian
D+xTC-L3) over-binds catastrophically: Li E = **−8.09 Ha** vs plain −7.24 and true TC
−7.28. D (the (∇τ)²+∇²τ term) without its convective counter-term K massively
over-corrects. (This is also why C's no-K energy −32.89 sits 2.4 Ha below plain −30.44:
that "win" is D-over-binding, not physical correlation. C's *sparsity* metrics are
unaffected and remain the meaningful C deliverable.)

**(6+7) Escape hatch — decisive.** Symmetrize the full TC operator,
H̃_sym = (H̃+H̃†)/2 (keeps the Hermitian part of K; the anti-Hermitian part is worth
only ~33 mHa):

| | κ_V | Pauli | λ/λ_H | cusp-win retained |
|---|---|---|---|---|
| Li sym(FULL) [Herm] | **1.00** | **117 (= plain)** | **0.978** | **177%** |
| He sym(D+K) [Herm] | 1.58 | **117 (= plain)** | **0.977** | 159% |

Symmetrization is a win on **every** QC axis at once: κ_V = 1, Pauli count back to
plain (117, from 249), λ **below** plain (0.98×), and it **retains ≥100% of the cusp
win** (over-recovers, since the anti-Hermitian residual carries almost none of the
correlation energy). H̃_sym is Hermitian ⇒ fits Paper 14's existing VQE/Trotter cost
with **no new axis**. → **GO-BY-AVOIDANCE; the escape hatch moots the whole
non-Hermitian question.**

**Fold H̃†H̃ is the wrong route.** It squares the conditioning: Li fold cond = 4.7e3
(vs κ_V 1.45), gap collapses to ~1.1 mHa. Symmetrization is strictly better.

---

## Decision gate

- **κ_V(low-lying):** He 1.98, Li 1.21, C 1.00 (Bauer–Fike κᵢ ≤ 1.22) → **GO** (O(1–10)).
- **λ_nonHerm/λ_Herm:** ≤ 1.122 → **GO** (≲1.5×). K inflates λ modestly and roughly
  doubles the Pauli count, but never blows λ up.
- **Escape hatch:** symmetrization retains 159–177% of the cusp win AND restores Pauli
  count 117 = plain, λ 0.98×, κ_V=1 → **GO-BY-AVOIDANCE (moots it).**

**Overall lean: GO — and the escape hatch moots it.** Run as a genuine non-Hermitian
object the operator is already benign (real spectrum, κ_V O(1), λ ≲1.13×). But you need
not: symmetrizing the full TC operator is strictly better on every quantum-cost axis
(Pauli, λ, κ_V, Hermiticity) while keeping the full cusp benefit. The xTC sparsity/
1-norm win survives the quantum-cost lens, and the recommended encoding carries **no new
non-Hermitian axis at all**. The only genuine cost of the non-Hermitian form — the ~2×
Pauli count from K — is exactly what symmetrization removes.

---

## Caveats

- **C measured without K** (the convective term is not in the validated s+p engine, and
  building an unvalidated p-orbital K was out of scope for a read-only reuse). C's
  operator is therefore Hermitian and its "non-Hermitian" measurement is trivial; the K
  non-Hermiticity is characterized on He/Li (s-only, where K is built and validated).
  C's contracted-L3 being Hermitian (2.8e-17) generalizes the s-reference result to an
  open-p reference. C's energy is a no-K over-binding artifact; C's sparsity/Pauli/1-norm
  metrics (279→279, l2 ratio 0.841) reproduce the pblock memo and are the meaningful C
  deliverable.
- Minimal single-common-k Coulomb-Sturmian basis (s-only He/Li, s+p C); single
  spin-independent geminal; spin-broken single-determinant reference (splits degenerate
  ground manifolds by ≤1.2 mHa). Classical PoC. None of these affect the κ_V / λ /
  symmetrization verdicts, which are operator-structure properties.
- TC is non-variational; production γ should be fixed by stationarity, not min-E — so the
  symmetrized over-recovery (E_sym below the true TC energy) is not "extra accuracy," it
  is a different Hermitian operator that preserves the cusp physics.

---

## Accuracy validation (2026-08-23) — does the symmetrization escape hatch keep the physics?

**Verdict: CAUTION/STOP.** The escape hatch is NOT a free, accuracy-preserving move.
Symmetrizing the full non-Hermitian TC operator (dropping the anti-Hermitian part of
the convective K) does not "over-recover ~100–177% of the cusp win" toward the true
energy — it **overshoots past the exact energy**, and it gets **worse as the basis
improves**. The "escape hatch retains the full cusp benefit" reading above is an
artifact of (i) evaluating at a single geminal γ=1 and (ii) reading a non-variational
energy as if lower = better.

Driver: `debug/xtc_sym_accuracy_validate.py` → `debug/data/xtc_sym_accuracy.json`.
Engines reused READ-ONLY (`xtc_poc_li.py`; orchestration only, no physics changed).
He (s-only, 2-body TC D+K) and Li (s-only, genuine 3-body → D+K+xTC-L3), PN-projected
FCI, k optimized per ns for plain, Ng=600/nx=96 (E_TC converged to 1e-6, E_plain to
~0.03 mHa). Exact NR energies He −2.90372, Li −7.47806.

**Table @ the memo's reference γ=1.0** (deltas to exact, mHa; **+ = above exact,
− = overshoot BELOW exact**):

| sys | ns | k | E_plain | E_TC | E_sym | exact | Δ_plain | Δ_TC | Δ_sym |
|---|---|---|---|---|---|---|---|---|---|
| He | 3 | 2.00 | −2.87815 | −2.89781 | −2.90970 | −2.90372 | +25.6 | +5.9 | **−6.0** |
| He | 4 | 2.00 | −2.87863 | −2.89814 | −2.91008 | −2.90372 | +25.1 | +5.6 | **−6.4** |
| He | 5 | 2.20 | −2.87892 | −2.89828 | −2.91035 | −2.90372 | +24.8 | +5.4 | **−6.6** |
| Li | 3 | 1.60 | −7.25072 | −7.29361 | −7.32897 | −7.47806 | +227.3 | +184.5 | +149.1 |
| Li | 4 | 1.60 | −7.39672 | −7.43788 | −7.46752 | −7.47806 | +81.3 | +40.2 | +10.5 |
| Li | 5 | 1.60 | −7.43052 | −7.47041 | −7.49945 | −7.47806 | +47.5 | +7.7 | **−21.4** |

**One-line: CAUTION/STOP — the ~33 mHa is the *anti-Hermitian part of K, and dropping
it moves the energy DOWNWARD PAST exact (AWAY from the true energy once the basis is
decent): at the best basis (He ns5 / Li ns5) γ=1, symmetrization takes E_TC from +5.4 /
+7.7 mHa ABOVE exact to −6.6 / −21.4 mHa BELOW exact — E_sym is on the WRONG side and
(Li) farther from exact than the genuine TC.**

### The four measurements

**(1) γ-stationarity does not exist.** For all six (system × ns) cases E_TC is
**monotone in γ** over [0.25, 3.0]: dE_TC/dγ > 0 everywhere (E keeps dropping as γ→0,
no turnover), so there is **no interior γ\* with dE_TC/dγ = 0**. The "principled
operating point" the sprint invokes is undefined for this minimal single-common-k
basis — the TC energy is a freely γ-tunable, non-variational number, not a stationary
observable. Any "accuracy at γ\*" statement is therefore ill-posed here; the γ=1 row is
the memo's reference, not a stationary point.

**(2) E_sym overshoots past exact — and the overshoot WORSENS with basis.**
The s-only *plain* FCI converges (ns=6) to a floor **25 mHa (He, −2.87901) / 35 mHa
(Li, −7.44312) ABOVE exact** — that residual is angular (l>0) correlation a *radial*
cusp geminal categorically cannot supply. So any energy dipping below that floor, let
alone below exact, is non-variational drift, not a basis-honest cusp fix. E_sym crosses
below exact at:

| | He | Li |
|---|---|---|
| γ where **E_TC** = exact | ns3 0.83 / ns4 0.84 / ns5 0.84 | ns3 none / ns4 0.53 / ns5 0.87 |
| γ where **E_sym** = exact | ns3 1.14 / ns4 1.15 / ns5 1.16 | ns3 0.47 / ns4 0.92 / ns5 1.31 |

The overshoot onset **moves to larger γ as the basis grows** (Li E_sym-crosses-exact:
0.47 → 0.92 → 1.31) — i.e. a *better* basis makes the overshoot *more* accessible, the
opposite of a controlled cusp fix converging on the true energy. At the memo's γ=1,
E_sym is already below exact for He (all ns) and for Li at the best basis (ns5).

**(3) |E_sym − E_TC| is NOT a fixed ~33 mHa residual — it diverges.** The "anti-Hermitian
part worth only ~33 mHa" is a γ=1 coincidence (Li ns3 γ=1: −35 mHa). Across γ the gap is
strongly γ-dependent and grows without bound as the geminal diffuses:

| γ | 3.0 | 2.0 | 1.5 | 1.0 | 0.6 | 0.4 | 0.25 |
|---|---|---|---|---|---|---|---|
| He ns3 sym−TC (mHa) | −0.7 | −2.3 | −5.0 | −11.9 | −26.7 | −41.9 | −60.3 |
| Li ns3 sym−TC (mHa) | −1.8 | −6.5 | −14.1 | −35.4 | −88.4 | −152.8 | −242.5 |

It tracks smoothly (monotone) but the discarded anti-Hermitian content reaches **hundreds
of mHa** — symmetrization is not shaving a negligible fixed sliver; it throws away a
γ-dependent chunk of the operator.

**(4) Basis trend: the correction is ~basis-independent (does NOT grow).** At γ=1.0 the
TC/sym *correction* (E_plain − E_{TC,sym}) is essentially flat (slightly shrinking) with
basis, while E_plain itself falls steeply toward exact:

| | He TC / sym (mHa) | Li TC / sym (mHa) |
|---|---|---|
| ns3 | +19.7 / +31.6 | +42.9 / +78.3 |
| ns4 | +19.5 / +31.4 | +41.2 / +70.8 |
| ns5 | +19.4 / +31.4 | +39.9 / +68.9 |

Because the correction is a roughly fixed downward shift while the plain basis is already
converging to within 25–35 mHa of exact on its own, the corrected energy is pushed
**further past exact as the basis improves** (Li Δ_sym: +149 → +10.5 → **−21.4** mHa
across ns3→5). A genuine cusp fix would move E toward exact and *stop* there; this shift
does not know where exact is.

### Bottom line for the escape-hatch finding

The earlier "GO-BY-AVOIDANCE — symmetrization moots the non-Hermitian question" verdict
stands **only on the quantum-cost axes** (Pauli count → plain, λ ≤ plain, κ_V = 1,
Hermitian). On **accuracy** it fails: symmetrization is strictly *worse* than the genuine
non-Hermitian TC (larger overshoot), and the genuine TC itself has no stationary γ and
overbinds past exact in this minimal s-only basis. The ~33 mHa sign is resolved:
**dropping the anti-Hermitian part of K moves the energy AWAY from exact** (downward,
past the true energy) once the basis reaches the point where plain FCI is near exact.
The "cusp-win retained ≥100%" number is over-recovery = non-variational overshoot, not
preserved physics — as CLAUDE.md §3 and the memo caveat already warned ("the symmetrized
over-recovery is not extra accuracy, it is a different Hermitian operator"). This
validation makes that caveat quantitative and upgrades it to a CAUTION: **do not read the
symmetrized (or the genuine) TC energy as an accuracy gain in this basis** — the escape
hatch buys quantum-cost, not accuracy, and symmetrization degrades accuracy relative to
the full operator.

**Caveats.** Minimal single-common-k Coulomb-Sturmian s-only basis; single
spin-independent geminal; spin-broken single-determinant xTC reference (the ns→exact
gap includes angular correlation the s-only basis cannot hold). These *strengthen* the
verdict: with a basis that cannot even reach exact, any E below exact is provably drift.
A larger l>0 basis + spin-resolved geminal + a real γ-optimization criterion (variance
minimization, not energy) would be needed before any accuracy claim; none of that
changes the structural finding that symmetrization moves the energy the wrong way.
