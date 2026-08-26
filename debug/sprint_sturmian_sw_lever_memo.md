# Sprint memo — Molecular Shibuya–Wulfman conditioning: momentum-space cross-check + levers

**Date:** 2026-08-18 · **Branch:** work/sparsity-boundary · **Diagnostic only** (Paper 60 §molecular backing).
**Drivers:** `debug/sturmian_sw_momentum.py` (this sprint) vs `debug/sturmian_sw_conditioning.py` (position-space baseline).

## Question
Can the molecular SW metric conditioning — which grows with basis size in the position-space
calc — be tamed or better characterized, via the momentum-space (Fock) evaluation and/or physical
levers? Two-center `s`-orbital Coulomb–Sturmians, common scale `k=1`.

## Method — independent momentum-space evaluation
The SW matrix in momentum space [PhD eq 10.5.13] is `S_{μ'μ}=∫d³p e^{ip·R}(2k/(k²+p²))³ Y*_{μ'}Y_μ`.
Under the Fock map `p→u∈S³` the measure `d³p=((k²+p²)/2k)³ dΩ₃` and the SW weight `(2k/(k²+p²))³`
**cancel the Jacobian exactly**, leaving an integral over the unit 3-sphere. For `s`-orbitals the
θ,φ plane-wave integral gives a sinc and the χ (polar-`S³`) integral collapses via
`sin²χ·U_{n'−1}U_{n−1}=sin(n'χ)sin(nχ)` to a clean **1-D integral**:

```
SW :  S_{n',n}(R) = (2/π) ∫₀^π sin(n'χ) sin(nχ) · sinc(kR·cot(χ/2)) dχ
L2 :  m_{n',n}(R) = (2/π) ∫₀^π sin(n'χ) sin(nχ) · (1−cosχ) · sinc(kR·cot(χ/2)) dχ
```
(`sinc(x)=sin x/x`). Full two-center matrix = `[[intra, inter(R)],[inter(R), intra]]`. This is a
*completely different quadrature* from the position-space 2-D grid + numerical gradients — different
variables, integrand, and singularity structure — so it is a genuine independent check. A **third**
route (direct 3-D momentum quadrature over `|p|`) confirms the 1-D χ-reduction to ~1e-10.

## Part (1) — VALIDATION GATE: **PASS**
| R | nmax | cond(L2) mom | cond(L2) pos | cond(SW) mom | cond(SW) pos |
|--:|--:|--:|--:|--:|--:|
|1.4|3|39.11|39.12|30.07|30.09|
|2.0|3|20.92|20.93|15.43|15.44|
|4.0|3|11.13|11.14| 4.43| 4.43|

All within 0.1%. Intra-center blocks:
- **SW intra = EXACT identity** (momentum: diag `[1,1,1]`, off-diag 7e-17, cond 1.0). The position
  grid's `{1,0.999,0.998}` / off-diag 2e-3 were **grid artifacts** — the identity is exact, confirming
  Paper 60 claim (i) at machine precision.
- **L2 intra** = tridiagonal `[[1,−½,0],[−½,1,−½],[0,−½,1]]`, cond 5.828 (position gave 5.83).

## Part (2) — GROWTH is POLYNOMIAL, not exponential
cond vs basis size `N=2·nmax`, nmax=1..12, log-log and log-linear fits:

| R | SW fit | L2 fit |
|--:|:--|:--|
|1.4| `cond ~ N^1.85` (R²=0.9999) | `N^1.70` (R²=0.9994) |
|2.0| `cond ~ N^1.82` (R²=0.9998) | `N^1.66` (R²=0.9994) |
|4.0| `cond ~ N^1.73` (R²=0.9990) | `N^1.83` (R²=0.9999) |

Power-law beats exponential decisively at every R (exp fits R²≈0.94). **Both SW and L2 grow at the
same polynomial rate ~N^1.7–1.85 (≈N²).** SW is a strict improvement over L2 (smaller prefactor at
bonding R, much better large-R limit), **not** a different asymptotic class. Spectral driver (L4):
`λ_max(SW)` saturates to **< 2** (bounded above); the growth is *entirely* `λ_min→0` as ≈N^−1.8
(basis approaching linear dependence). No plateau through nmax=12 (cond 396 at N=24, R=1.4).

## Part (3) — LEVERS
**L1 — large separation R:** cond(SW) → 1 as R→∞, and fast (nmax=8: 355→24→6.7→1.6 over R=1→4→8→20).
**Crucially L2 does NOT** — it plateaus at the intra-center L2 conditioning, which itself grows with
basis (5.8, 14.1, 38.6 for nmax=3,5,8 at R=20). SW's floor is 1 (intra=identity); L2's floor grows.
So SW's large-R advantage over L2 is *qualitative and unbounded*. But chemistry lives at bonding
R~1.4–2, where cond is largest.

**L2 — momentum vs position give the SAME matrix ⇒ SAME cond.** Momentum evaluation removes grid
error and is cheap/exact (1-D integral or 3-term recurrence), but the condition number is *intrinsic
to S* and invariant to how it is computed. **Momentum-space does NOT give a better-conditioned
matrix** — it is an *evaluation* win, not a *conditioning* win. (Honest correction to the paper's
"momentum-space evaluation" lever framing.)

**L3/L5 — gerade/ungerade symmetry adaptation (homonuclear):** the orthogonal transform
`(A_n±B_n)/√2` **block-diagonalizes S exactly** (off-block 1e-16), splitting it into a
**perfectly-conditioned gerade sector** (cond ≈ 1.5–2.3, essentially flat in both N and R) and an
ungerade sector carrying all the ill-conditioning (cond N^2.00). Because the molecular secular
problem `[W−kS]C=0` respects the same g/u symmetry (identical centers), the two sectors can be solved
**independently**, so the relevant cost is `max(cond_g,cond_u) ≈ cond_u`, not
`cond(full)=λ_max(g)/λ_min(u)`. This is a **constant-factor ~1.75–2× win plus a free gerade sector**
— genuine in the symmetric regime — but it does *not* change the ungerade sector's polynomial growth.

## HONEST TWO-WAY VERDICT
**Tamed, partially — and much better characterized than "grows with basis size":**
1. The growth is **POLYNOMIAL (~N^1.8), not exponential** (R²>0.999). This is the headline: a poly(N)
   conditioning multiplier is a *manageable* block-encoding cost, not a showstopper. `λ_max<2` bounded;
   only `λ_min→0` polynomially.
2. **Momentum-space is conditioning-neutral** (same matrix) but an evaluation win — exact, grid-free,
   3-term recurrence. It buys cheap construction, not better conditioning.
3. **Two real mitigation levers:** (a) larger effective separation R drives cond(SW)→1 (SW-specific;
   L2 cannot follow); (b) gerade/ungerade splitting isolates a perfectly-conditioned sector and roughly
   halves the effective cost in the symmetric case.

**Genuine residual cost the block-encoding must pay:** at fixed bonding geometry the ungerade sector's
`λ_min→0` as ~N^−2, so a **poly(N) conditioning multiplier is genuinely paid** and cannot be
transformed away (orthogonal/unitary maps leave cond invariant; momentum evaluation reproduces the
same matrix). The levers reduce the prefactor and isolate a good sector; they do not remove the
polynomial growth.

**Net:** the molecular metric is *mitigated and now quantified* — polynomial not exponential, with SW
strictly better than L2 and two concrete levers — but *not dissolved*. The paper's "mitigated, not
dissolved" framing is correct; this sprint sharpens it to "polynomially-growing, momentum-neutral,
lever-reducible."

## Suggested Paper 60 §molecular refinements (PI integrates)
- The SW intra-block identity is **exact** (machine precision in momentum form), not `{1,0.999,0.998}`
  (those are position-grid artifacts).
- Add the growth rate: **cond(SW) ~ N^1.8 (polynomial, R²>0.999)**, same rate as L2, `λ_max<2`.
- Correct the momentum-space lever: it is an **evaluation** win (exact/cheap), **conditioning-neutral**
  (same matrix, same cond) — not a route to a better-conditioned matrix.
- Add the two working levers: large-R (cond→1, L2 cannot) and gerade/ungerade splitting
  (perfectly-conditioned gerade sector; ~2× effective-cost reduction in the symmetric case).
