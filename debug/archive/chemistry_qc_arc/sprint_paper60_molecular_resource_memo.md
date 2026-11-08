# Sprint memo — Paper 60 molecular resource estimate (2026-08-18)

**Task.** Turn the measured molecular Shibuya–Wulfman conditioning (`cond(S)~N^1.85`) and the
gerade/ungerade lever (Paper 60 §molecular) into an *actual* qubit / T-gate (block-encoding)
resource count for the one-electron H₂⁺ isoenergetic secular equation `[W − kS]C=0`, head-to-head
against the same problem in a Gaussian LCAO metric.

**Cost model (grounded in the paper's already-cited prior art).** Whiten the generalized
eigenproblem to the standard one `M̃ = S^{-1/2} W S^{-1/2}` (Liang 2112.02554); the metric enters
only through `S^{-1/2}`, applied by QSVT with a degree-`d_inv ~ κ·ln(κ/ε)` polynomial for `x^{-1/2}`
on `[1/κ,1]` (Gilyén–Su–Low–Wiebe 1806.01838; Childs–Kothari–Somma 1511.02306), `κ=cond(S)`. This
`d_inv` *is* the metric penalty. QPE queries `Q ~ (π/2)λ_eff/ε_k` (Low–Chuang 1610.06546),
`λ_eff = k_max`. Each query = 1 BE(W) + `d_inv` BE(S).

## Results (H₂⁺, R=2 bohr, chemical accuracy ε_E=1.6 mHa)

**Conditioning (exact momentum-space, `debug/sturmian_molecular_resource.py` Part 1):**

| metric | growth | cond @ N=8 | cond @ N=16 |
|---|---|---|---|
| SW full | N^1.81 (R²=1.000) | 25.8 | 91.8 |
| SW gerade | **flat ~2.0–2.3** | 2.03 | 2.28 |
| SW ungerade | N^1.86 | 14.1 | 52.5 |
| Gaussian (even-temp ratio 2) | N^3.52 (R²=0.989) | 2.0e3 | 2.0e4 |

**Resource table (Part 2), N=16 (n_max=8/center), n_phase=10, λ_eff=1.72:**

| metric | κ | d_inv | qubits |
|---|---|---|---|
| atomic (no metric) | 1 | 0 | 23 |
| **SW gerade (σ_g ground)** | **2.3** | **18** | **26** |
| SW full | 91.8 | 1050 | 31 |
| Gaussian LCAO (ratio 2) | 1.99e4 | 3.3e5 | 31 |

**Head-to-head (Part 3):** Gauss/SW-full d_inv ratio 20–335×; Gauss/SW-gerade 100–27000×.
**Large-R lever (Part 4):** SW cond → 2.2 at R=10 (d_inv 862→17 over R=1.4→10); Gaussian stays ~2.3e3.

## The headline — the gerade lever is the molecular payoff

The H₂⁺ ground state σ_g is **gerade**. The gerade sector of the SW metric block-diagonalizes off
and is **flat, κ≈2 across N=4→20**, so the ground-state `S^{-1/2}` costs `d_inv≈16–18` *independent
of basis size* — a ~60× reduction versus the full metric at N=16, on a register halved to
⌈log₂(N/2)⌉ qubits. This is *more* than the ~2× symmetry factor already noted in the paper: the
metric penalty for the chemically relevant state **does not grow**. Growth is confined to the
excited/ungerade states.

## Honest caveats (kept in the paper)

1. **Resource-MODEL tier, not measured.** d_inv carries an O(1) constant → *absolute* counts are
   order-of-magnitude; the O(1) cancels in cross-metric *ratios*, which are the robust deliverable.
2. **Gaussian comparison is ratio-dependent** (`sturmian_h2plus_validate.py` Part 2): ratio 1.6
   (near-complete) → N^6.1, κ~1e5; ratio 3.0 (sparse) → N^1.4 but ~10× the SW prefactor + coverage
   gaps. This is the classic Gaussian coverage-vs-linear-dependence trade-off; the Coulomb-Sturmians
   evade it (complete at one scale, polynomially conditioned, small prefactor). Framed as trade-off,
   NOT "Gaussian always worse."
3. **Production Gaussian QC absorbs the overlap classically** in a prior mean-field step. The metric
   penalty priced here is the cost of the *on-device angular generation* that motivates the Sturmian
   encoding in the first place — the honest scope of the comparison.
4. **One-electron only.** Matches the paper's SW §; the many-electron molecular case (T′ on top of
   the SW metric) is not measured and not claimed.

## Validation (`sturmian_h2plus_validate.py` Part 1; test `test_paper60_h2plus_isoenergetic_binds`)

Solving `[W − kS]C=0` for H₂⁺ by isoenergetic self-consistency (build at scale k, find k′=k) gives
σ_g at **k=1.475 → E_elec=−1.088 Ha** vs reference −1.1026 Ha (1.3% s-only underbinding, no
p-polarization, stays above exact). Pins λ_eff = k_max ≈ 1.72. The algorithm genuinely binds the
molecule — the resource estimate rests on a working method, not an abstraction.

## Artifacts

- `debug/sturmian_molecular_resource.py` — resource driver (Parts 1–4).
- `debug/sturmian_h2plus_validate.py` — H₂⁺ isoenergetic solve + Gaussian ratio robustness.
- `papers/group2_quantum_chemistry/paper_60_*.tex` — new `\subsection{Resource estimate}` (sec:resource)
  + `tab:resource`; abstract + conclusion updated; +2 cites (gslw2019, cks2017). 5-pass clean, 5 pp.
- `tests/test_paper60_sturmian.py` — +3 tests (gerade flat; SW beats Gaussian metric; H₂⁺ binds).
  Full file 11/11 pass incl. `--slow`.

## Status

Certified→**Phase-4 re-review OWED compounds** (consistent with 58/59). Paper edits autonomous
(§13.8). Git commit / version tag / release deferred to PI (close-out phase: merge-to-main + Release
are PI-only; a Release mints a Zenodo DOI).

---

## Follow-on (same day) — many-electron molecular case: the metric does not compound

**Question (PI):** can we do the many-electron molecular case?

**Result [MEASURED].** For a many-electron molecule the generalized-Sturmian secular equation is
still `[T⁰ + T' − pκ 1]B = 0` — a *standard* (metric-free) eigenproblem (Avery BK6 eq 6.19 general
form); only the atomic conveniences (T⁰ diagonal, T′ pure-number) are lost. The resource question is
whether the SW metric penalty **compounds** with electron number. **It does not.**

A k-electron configuration overlap is the **k-th compound matrix** of the one-electron overlap, with
conditioning `(λ₁···λ_k)/(λ_{N−k+1}···λ_N)` in the one-electron overlap spectrum {λ_i}:

| basis for the orbitals | k=1 | k=2 | k=3 (half-fill) |
|---|---|---|---|
| RAW L² many-center CS (N=6, R=1.4) | 39.1 (1×) | 320 (8.2×) | 365 (9.3×) |
| SW molecular Sturmians (S_SW-orthonormal) | **1.000** | **1.000** | **1.000** |

(cond values cross-checked bit-exactly against the compound-matrix eigenvalue product formula.)
Raw grows with electron number (worst at half-filling, bounded by cond(S_L²)^min(k,N−k), driven by
the near-linear-dependence eigenvalue). The SW route — diagonalize the one-electron SW problem ONCE
(cost = cond(S_SW), priced above, g/u-improvable) — gives molecular Sturmians that are S_SW-orthonormal,
so determinants of them have **identity** config overlap at every electron number.

**⇒ The molecular metric penalty is one-electron, not N-electron.** Paid once at the one-electron
molecular-orbital construction; the many-electron secular equation inherits an identity metric. This
is the on-device analogue of the once-executed classical mean-field step. The many-electron extension
of the one-electron result.

**Honest boundary [OPEN].** What remains at the N-electron level is the block-encoding 1-norm of
`T⁰ + T'`. T⁰ is one-electron/bounded; T′ is the two-center ERI tensor, which closes in this framework
weight-1 π-free {E₁,ln,γ} (Papers 58/59; all four classes closed-form on THIS branch, `two_center_eri.py`
`aabb_value` + hybrid/exchange). Whether the config-space 1-norm stays sublinear (as it does for atoms,
eq:sublinear) is the concrete next measurement — needs the shared-scale Sturmian two-center ERI tensor
assembled (the production engine is *hydrogenic* a=Z/n, so this is a wiring task, not a new integral).
Unlike atoms, molecular T′ is not pure-number (retains geometry dependence) but stays pκ-independent
(isoenergetic) → assembled once per geometry, not per energy.

**Not attempted this pass (honest scope):** the full interacting isoenergetic H₂ solve (all 4 ERI
classes at shared scale + 2-electron self-consistency). The metric-compounding result is the rigorous,
self-contained deliverable and it answers the resource question; the interacting binding + T′ 1-norm is
the next sprint (engine exists, wiring + basis-scale reconciliation needed).

**Artifacts:** `debug/sturmian_manyelectron_metric.py`; Paper 60 new `\subsection{Many-electron
molecules}` (sec:manyelectron); abstract + conclusion updated; +1 test
`test_paper60_manyelectron_metric_does_not_compound` (12/12 pass). 3-pass clean, 6 pp.

---

## Follow-on 2 (same day) — enabling ERI validated + interacting H₂ binds; 1-norm still OPEN

**Push:** wire the two-center ERIs into an interacting many-electron molecule and measure the
config-space 1-norm (the sec:manyelectron [OPEN]).

**Enabling ERI — VALIDATED [MEASURED].** Built a self-contained numerical two-center shared-scale
Coulomb-Sturmian s-ERI (multipole expansion about nucleus A; every center class uniform). Validated
against the framework's EXACT closed form:
- one-center (1s1s|1s1s) = 5k/8: **err 3.5e-6**
- two-center (AA|BB) vs `two_center_eri.aabb_value`: **err ~1e-4** at R=1.4/2.0/3.0.
So the shared-scale Sturmian two-center ERIs (the molecular T′) are computable and validated. (The
production `aabb_value` is hydrogenic a=Z/n; shared-scale reached by Z=k·n per orbital — but the
numerical multipole route sidesteps per-class engine mapping and is what's validated here.)

**Interacting H₂ BINDS [MEASURED].** Minimal basis (1s on each center, common scale k=1), metric-free
molecular-Sturmian basis, eps_g computed on the SAME grid/scale as J (consistent):

| R (bohr) | eps_g | J(σ_g²) | E_elec | E_tot |
|---|---|---|---|---|
| 1.4 | −1.183 | 0.566 | −1.800 | −1.086 |
| **1.6** | −1.136 | 0.553 | −1.719 | **−1.094** |
| 2.0 | −1.052 | 0.527 | −1.576 | −1.076 |

**Interior minimum at R≈1.6, E_tot≈−1.094** = the textbook single-ζ minimal-basis H₂ value (Szabo–Ostlund;
underbinds exact −1.174 as s-only single-ζ must, stays above exact). The interacting many-electron
molecular method binds a real molecule. (My first pass used a made-up too-shallow eps table giving
E_tot~−0.80 — WRONG; fixed by computing eps rigorously on the grid. Lesson: no hard-coded reference
numbers, compute consistently.)

**Config-space 1-norm — still OPEN (honest).** Measured a PROXY config 1-norm (K^2.72 over K=3,10,21)
but it is NOT a faithful test of the atomic sublinearity (eq:sublinear): it omits the isoenergetic
−1/pκ weighting and the Goscinskian pure-number structure that MADE the atomic T′ decay with config
index, h_diag was zeroed, and ~K² mostly reflects the K² entry count. A faithful test needs the full
isoenergetic secular matrix at basis > s-only (which reaches too few configs for a clean exponent).
The enabling integrals are now validated; the scaling verdict is the next step. **NOT put in the paper
as an exponent.**

**Paper:** sec:manyelectron [OPEN]→ split into [MEASURED] (ERI validated + H₂ binds) + a sharpened
[OPEN] (the full isoenergetic secular-matrix 1-norm). +2 tests (ERI-vs-exact, H₂-binds). 14/14 pass.

**Artifacts:** `debug/sturmian_h2_eri.py` (validated multipole two-center s-ERI + binding + proxy
1-norm). 3-pass clean, 6 pp.

---

## Follow-on 3 (same day) — the isoenergetic-1-norm sublinearity is a SINGLE-CENTER property (obstruction)

**Push:** assemble the actual isoenergetic molecular secular matrix (−1/pκ weighting) at p-inclusive
basis, measure whether its 1-norm is sublinear.

**Result: OBSTRUCTION identified — no trustworthy number.** The first assembly (fixed-scale
molecular-Sturmian basis, T⁰ = diag(Σ_p λ_p)) produced garbage (self-consistency pκ~3.5, E~−6 Ha).
Root cause, verified exactly:

- The isoenergetic collective scale is the **root-sum-of-squares** √(Σ k_p²), NOT a plain sum. For
  He 1s²: √(k²+k²) = 2√2 = **Z R_ν** → bare E = −4 (paper, correct); a plain sum gives −8 (my bug).
- So a fixed-scale shortcut (T⁰=Σλ) FAILS the non-interacting limit: pκ = 2k_σg where the correct
  value is √2 k_σg (E off by 2×). Confirmed numerically both ways.

**Why the atomic sublinearity does not transfer.** eq:sublinear rides on the clean **diagonal
T⁰ = Z R_ν**, which exists only because (i) the Goscinskian orbitals scale WITH pκ (Q_ν=pκ/R_ν) so pκ
cancels the −1/pκ prefactor → M is pκ-independent linear, and (ii) on ONE center ⟨1/r⟩=Q/n² collapses
the collective scale to √(Σk²)=ZR_ν. On TWO centers ⟨1/r_A⟩+⟨1/r_B⟩ gains off-center, Q·R-dependent
terms → T⁰ is non-diagonal AND pκ-dependent (the same kR nonlinearity as the one-electron SW case).
So the single-center structure that PRODUCED the sublinearity is gone.

**Faithful measurement requires** the pκ-scaled molecular Goscinskian: 2-electron configs of atomic
Sturmians on the two centers at config charge Q_ν=pκ/R_ν, two-center nuclear T⁰ + two-center ERIs at
config-specific scales, solved self-consistently. **Validation gate: R→∞ = 2 H atoms (−1.0 Ha),
R→0 = He (−2.847).** That's the defined next step — NOT the fixed-scale shortcut.

**Discipline note.** Caught the wrong T⁰ via its own non-interacting-limit sanity check before it
reached the paper. The K^0.89 the shortcut produced is discarded (built on the doubled diagonal). The
honest finding — sublinearity is single-center-specific, obstruction identified — is more informative
than a (wrong) number would have been.

**Paper:** sec:manyelectron [OPEN] rewritten as [OPEN, obstruction identified]; abstract [OPEN]
sharpened. +1 test `test_paper60_collective_scale_is_root_sum_of_squares` (15/15). Driver
`debug/sturmian_iso_secular_1norm.py` is now the obstruction diagnostic. 3-pass clean, 6 pp.

---

## Follow-on 4 (same day) — validated two-center CI + the 1-norm verdict: NOT sublinear for molecules

**Push:** build the p_kappa-scaled molecular Goscinskian and measure the isoenergetic 1-norm.

**What I built + validated.**
- **Mixed-charge two-center integral layer** (`debug/sturmian_goscinskian_integrals.py`): s-orbital
  overlap, nuclear attraction (1/r_A spherical + 1/r_B Legendre-expanded about A), ERIs, kinetic
  (via the shared-scale Sturmian ODE −½∇²χ_n = nk/r_c χ_n − ½k²χ_n, no gradients). Validated:
  ⟨1s|1/r_A|1s⟩=a (1e-4), ⟨1s_A|1/r_B|1s_A⟩ vs closed form (7e-5), (1s1s|1s1s)=5a/8 (1e-5),
  (AA|BB) vs `aabb_value` (1.5e-4), mixed overlap (5e-10), H-atom E=−0.5 exact.
- **Two-center Sturmian 2-electron H₂ CI** (`debug/sturmian_h2_ci_1norm.py`): **dissociates correctly**
  — E_tot→−1.0 Ha (2 H atoms) at R=6 (bit at R=20 needs higher Lmax = single-center-expansion
  anisotropy limit, not physics); binds E_tot≈−1.15 at R=1.4 (above exact −1.174). Validated solver.

**Two bugs found + fixed via the validation gates** (the gates earned their keep):
1. ERI multipole grouping (i,k)(j,l) → (i,j)(k,l) [chemist density pairing]. Caught by (AA|BB) at a≠1.
2. Löwdin ERI transform einsum `'pi,...'` → `'ip,...'` (transformed with X^T not X). Caught by FCI
   landing above RHF / R→∞ not giving −1.0.

**VERDICT [MEASURED]: the molecular 1-norm is NOT sublinear.** The standard block-encoding
λ = Σ|h| + Σ|(pq|rs)| in the Löwdin-orthonormalized Sturmian basis grows **polynomially,
λ ~ n_orb^2.2** (R=1.4), no sublinear behavior. Combined with the obstruction (follow-on 3): the
atomic sublinearity (eq:sublinear, K^0.84) is a property of the **isoenergetic secular matrix**, not
of the basis, and it does not survive the loss of the single-center diagonal T⁰=ZR_ν. So molecules
pay the ordinary polynomial block-encoding cost; the real savings are the metric levers (gerade
sector, large-R, one-electron-only metric), not a sublinear matrix.

**Not delivered (honest):** the exact isoenergetic molecular *secular matrix* M B=pκB with the −pκ 1
form for general V₀ was NOT assembled (its identity-RHS derivation for non-orthogonal molecular
configs needs the primary source). But the resource QUESTION is answered via the standard λ, which is
what actually sets qubitization cost.

**Paper:** sec:manyelectron + abstract [OPEN]→[MEASURED] (answered negative). +1 test
`test_paper60_h2_ci_dissociates_and_binds` (16/16). Integral layer + CI drivers reusable. 3-pass clean.
