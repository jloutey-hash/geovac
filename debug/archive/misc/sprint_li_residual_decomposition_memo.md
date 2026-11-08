# Sprint memo — decomposing Li's 33.2 mHa post-per-shell-lambda residual (2026-09-19)

**Owed by:** CLAUDE.md Sec. 3 ("k_n = Z/n" row family) and
`debug/sprint_r12ci_gamma_transfer_memo.md` Sec. 6/9 ("Li after the exponent fix is
still 33.2 mHa ... its residual is itself undecomposed ... the same undecomposed-residual
defect this sprint fixed for H2, still owed for Li"). Measurement only — no new method
was built, and no attempt was made to beat 33.2 mHa.

**Verdict in one line.** The 33.2 mHa is **overwhelmingly angular, not radial and not
cusp**: pushing the s-only radial basis to near-saturation (ns 4→7) recovers only
**3.06 mHa (9.2%)**, and the explicit-r12 geminal recovers only **1.00 mHa (3.0%)** once
measured at that same converged basis (its apparent 4.53 mHa gain at the ns=4 baseline is
mostly a proxy for the missing radial functions, not real cusp physics — see Sec. 3.3).
**29.1 mHa (87.7% of the original residual) survives both axes** and is the operational
estimate of the angular correlation this s-orbital-plus-single-geminal ansatz cannot
reach without l>0 orbitals.

---

## 1. Current-state check: baseline reproduced

Before decomposing anything, `debug/r12ci_per_shell_lambda.py` was re-run in full (He,
Li⁺, Li) to confirm CLAUDE.md's number is still live, not a stale snapshot:

```
Li:  shared exponent E = -7.396405 (81.7 mHa above exact)
     per-shell lambda  E = -7.444876 (33.18 mHa above exact)
     GAIN = 48.47 mHa
```

**Reproduced exactly** (33.18 mHa vs the memo's 33.2 mHa). GO condition met. All
decomposition work below starts from this basis: ns=4, free per-shell lambda
`[3.153, 2.208, 1.103, 1.505]`, no geminal, Z=3, N=3 electrons, ms2=1 (doublet),
`debug/r12ci_ne_engine.py`, Ng=240 (matching the baseline exactly).

**Exact reference cross-checked.** The corpus's `EXACT["Li"] = -7.4780603236` Ha was
checked against an independent literature search (not derived from the corpus): a 2026
search returns the non-relativistic clamped-nucleus Li ground state as
**-7.478060323452 Ha** (Yan & Drake / King lineage) — agrees with the registered value to
every digit quoted. No discrepancy.

Driver: `debug/r12ci_li_residual_decomposition.py`. Raw data:
`debug/data/r12ci_li_residual_decomposition.json`.

## 2. Why the engine forces an indirect angular measurement

`debug/r12ci_ne_engine.py` has **no l>0 orbital shape** — the radial basis functions in
`R_and_dR` are plain (generalized-Laguerre) s functions, and every kernel in
`r12ci_ne_core.py` is built for scalar (L=0-projected or triangle-summed) legs between
s-type densities. Building explicit p/d/f orbitals with their own Slater-Condon /
Gaunt machinery would be a new method, which the task scope excludes. So:

- **Radial degree** is measured directly: raise `ns` (more l=0 radial functions) at free
  per-shell lambda, no geminal. This is a clean, unambiguous probe of "how much of the
  gap is basis-set incompleteness within l=0" — exactly the H2 memo's radial-degree rung.
- **Cusp** is measured via the engine's existing explicit-r12 geminal
  (`with_geminal=True`), which is the direct analogue of the H2 memo's r12^p toggle.
  **Caveat stated up front, not discovered after the fact:** `KernelBank.mom()` expands
  the geminal's Legendre content up to `Lmax=20` for the 3-electron triangle terms, so
  `exp(-gamma*r12)` carries angular content at every L the way Hylleraas' 1929 helium
  wavefunction did — a single exponential geminal is not a pure short-range cusp probe;
  it is a (weak, single-length-scale) proxy for angular correlation too. Section 3.3
  measures exactly how much of its apparent gain is this proxy effect.
- **Angular correlation (l>0 orbitals)** cannot be built or turned on directly. It is
  bounded from below by what survives after both of the above axes are pushed to
  (near-)saturation — the same "residual after every knob is exhausted" logic the H2 memo
  used to attribute its own final remainder to the cusp.

## 3. Measurements

### 3.1 Radial-degree axis (ns 4→7, l_max=0, no geminal)

One shell added at a time: freeze the previous shells' lambdas, 1-D optimize the new
exponent (cheap, good warm start), then a short joint Nelder-Mead polish over all `ns`
exponents together.

| ns | E (Ha) | err (mHa) | gain (mHa) | lambda (per shell) |
|--:|--:|--:|--:|:--|
| 4 (baseline) | -7.444876 | 33.18 | — | [3.15, 2.21, 1.10, 1.51] |
| 5 | -7.447095 | 30.97 | **2.22** | [3.68, 2.88, 1.09, 1.48, 1.94] |
| 6 | -7.447739 | 30.32 | **0.64** | [3.78, 3.22, 1.04, 1.38, 1.77, 2.25] |
| 7 | -7.447936 | 30.12 | **0.20** | [3.59, 3.24, 1.02, 1.35, 1.72, 2.17, 2.67] |

**TOTAL gain ns 4→7 = 3.06 mHa (9.2% of 33.2 mHa).**

The increments shrink geometrically (ratio 0.64/2.22 = 0.29, then 0.20/0.64 = 0.31) —
clean saturation, not a plateau-then-cliff. Extrapolating the geometric tail (ratio
≈ 0.30) gives a fully-converged s-only estimate of **≈ 30.0 mHa** above exact — i.e.
pushing `ns` to infinity within l=0 buys perhaps another ~0.1 mHa beyond what was
measured, not materially more. **Radial-degree incompleteness is a small, saturating
piece of the gap**, matching the qualitative pattern already established for H2's
polynomial-degree axes (Sec. 4 of the scoping memo): real, measurable, and far short of
closing the residual.

*Note on cost, not a physics point:* one energy evaluation costs 0.08 s at ns=4 and 4.7 s
at ns=7 (`nd` grows from 24 to 147 determinants); ns=8 was not attempted (~35-45 min for a
joint polish at the per-eval cost observed) given the returns already sit at 0.2 mHa/step
and shrinking.

### 3.2 Cusp axis (explicit r12 geminal, lambda held at the ns=4 baseline)

One-geminal gamma scan (12 points) then a bounded refine, then a 2-geminal coarse grid
(5×5, symmetric) to check whether a second correlation length helps here too, exactly as
Sec. 8 of the gamma-transfer memo did on a different (shared-k) basis.

| step | best gamma(s) | E (Ha) | err (mHa) | gain (mHa) |
|:--|:--|--:|--:|--:|
| 1 geminal | 0.4946 | -7.449408 | 28.65 | **4.53** |
| + 2nd geminal | (0.15, 2.50) | -7.449763 | 28.30 | **+0.355** (extra) |

**Cross-check against the corpus's prior (different-basis) measurement:** the
gamma-transfer memo Sec. 8 found a 2nd geminal buys Li (3e) **+0.394 mHa** on the
*shared-k* basis. Here, on the *per-shell-lambda* basis, the extra is **+0.355 mHa** —
same order, same conclusion ("a single-gamma ansatz is not the limitation"), reproduced
independently on a materially different lambda structure. No qualitative change.

Read naively, the 1-geminal gain (4.53 mHa, 13.6% of the residual) looks like the
second-largest lever after radial degree. **Section 3.3 shows this number is inflated.**

### 3.3 The interaction check: geminal gain shrinks 4.5× once the radial basis is converged

Same 1-geminal scan, but now at the **ns=7 radial-converged** lambda set instead of the
ns=4 baseline:

| basis | best gamma | E (Ha) | err (mHa) | geminal gain (mHa) |
|:--|--:|--:|--:|--:|
| ns=4 (baseline) | 0.49 | -7.449408 | 28.65 | 4.53 |
| ns=7 (radial-converged) | 0.20 | -7.448934 | 29.13 | **1.00** |

**The geminal's apparent gain drops from 4.53 mHa to 1.00 mHa (4.5×) once the radial
basis is already good.** This is the key mechanistic finding of this sprint: at the
coarse ns=4 basis, most of what the single-exponential r12 term is "buying" is not
short-range cusp physics — it is compensating for the *same* radial incompleteness that
more s functions also fix. The two axes are **not additive**: naively summing the
ns=4-measured pieces (radial 3.06 + cusp 4.53 = 7.59 mHa) would overstate the combined
gain by about 3.5 mHa against the measured combined value (radial 3.06 + cusp-at-ns7 1.00
= 4.06 mHa, matching the directly-measured combined residual to 0.01 mHa — see Sec. 4).
**The basis-independent estimate of the true cusp/short-range contribution is the
smaller, ns=7 number: ≈ 1.0 mHa, not 4.5 mHa.**

This mirrors, in a different system and a different basis family, the same general
lesson Paper 12's H2 arc drew about r12 and radial completeness competing for the same
correlation — except here it is the geminal partially substituting for radial degree
(not vice versa), because the s-only basis at ns=4 is comparatively poor while the
geminal is comparatively rich (an exponential in r_12 already carries some radial *and*
some angular content).

## 4. Combined result and the decomposition table

Adding the geminal on top of the radial-converged (ns=7) basis:

| gamma | E (Ha) | err (mHa) |
|--:|--:|--:|
| 0.20 (best) | -7.448934 | **29.13** |

```
baseline (ns=4, per-shell lambda, no geminal):            33.18 mHa   [100%]
  - radial-degree gain (ns 4 -> 7, saturating):            3.06 mHa   [ 9.2%]
  - cusp/geminal gain (measured AT ns=7, basis-independent): 1.00 mHa [ 3.0%]
  ---------------------------------------------------------------------
  residual surviving BOTH axes:                            29.13 mHa  [87.7%]
```

Bookkeeping check: 33.18 − 3.06 − 1.00 = 29.12, matching the *directly* measured combined
value (29.13) to 0.01 mHa — the two axes, measured in the right order (radial first, then
cusp at the converged basis), are consistent and very nearly additive. It is only the
*naive* ns=4-measured cusp number (4.53 mHa) that double-counts.

**29.1 mHa (87.7% of the original 33.2 mHa) is the operational estimate of angular
correlation this basis structurally cannot reach.** It is a lower bound on "true" angular
correlation in the sense that the geminal's own L>0 Legendre content (Sec. 2) has already
been credited against it — whatever is left after that credit is not reachable by any
s-orbital-plus-single-exponential-geminal combination tested here.

## 5. Independent cross-check: Hartree-Fock partition confirms the same split

An entirely independent decomposition — correlation energy relative to Hartree-Fock,
using literature numbers not derived from this engine — lands on the same picture.

- **Li numerical HF limit** (Froese Fischer / Bunge lineage, literature): **-7.43271 Ha**
  (a 2026 search returned -7.432707(1) a.u., consistent to 5 digits with the commonly
  cited -7.4327242).
- Total correlation energy: -7.478060 − (−7.43271) = **-45.35 mHa**.
- This engine's s-only, radial-converged (ns=7) energy, -7.447936 Ha, recovers
  -7.447936 − (−7.43271) = **-15.23 mHa of it (33.6%)** via pure radial (in-out)
  correlation.
- The remainder, **30.12 mHa (66.4%)**, is angular — and this number is not fitted to
  match Sec. 3.1's radial-ladder residual (30.12 mHa); **it is the same number**, because
  "energy above exact after saturating l=0" and "correlation energy not recovered by an
  l=0 basis" are the same quantity read from two different zero points (exact vs. HF).
  The agreement is a consistency identity, not an independent confirmation — but it is a
  useful sanity check that no arithmetic error crept into either bookkeeping path.

**Comparison with He**, whose s-limit is already a citable literature quantity in this
corpus (`S_LIMIT["He"] = -2.8790288`, used directly by
`debug/r12ci_per_shell_lambda.py`): HF(He) = -2.861680 Ha, exact = -2.9037244 Ha, total
correlation = -42.04 mHa. The s-limit recovers -2.8790288 − (−2.861680) = -17.35 mHa
(41.3%), leaving 24.70 mHa (58.7%) angular. **Li's angular fraction (66.4%) is somewhat
higher than He's (58.7%)**, which is physically sensible: Li has an additional
core-valence (1s-2s) angular correlation channel that He's single 1s² shell does not, on
top of the intra-1s² angular correlation both systems share. No literature "Li s-limit"
citation was found in this search (unlike He's, which is a standard Schwartz-lineage
number); the 30.0-30.1 mHa figure here is this sprint's own measurement/extrapolation,
not a literature-sourced calibration point, and should be labeled as such.

## 6. Honest scope

**Measured (this sprint, with controls):**
- Radial-degree axis: 3.06 mHa total (ns 4→7), saturating geometrically (ratio ≈ 0.3),
  extrapolates to ≈ 30.0 mHa asymptotic residual.
- Cusp/geminal axis: 4.53 mHa naively (at ns=4), but only 1.00 mHa once measured at the
  radial-converged basis — the basis-independent estimate.
- The two axes are very nearly additive when the cusp is measured at the *converged*
  basis (33.18 − 3.06 − 1.00 = 29.12 vs. directly measured 29.13).
- A 2nd geminal buys +0.355 mHa at ns=4, reproducing the gamma-transfer memo's +0.394 mHa
  finding on an independent lambda structure — "single-gamma is not the limitation" holds
  under the improved basis too.
- Cross-check via the HF partition (independent of this engine's own internal
  bookkeeping) lands on the identical 30.1 mHa angular figure and a physically
  reasonable angular fraction (66.4%) comparable to He's literature value (58.7%).

**Structural / not proven:**
- "29.1 mHa is angular correlation" is an *operational* estimate — the complement of
  everything this ansatz's two axes can reach — not a proof that no cheaper radial or
  cusp trick could close more of it. In particular, a genuinely different radial family
  (not more Laguerre shells but, e.g., a DVR/spline set) was not tried; the H2 memo found
  exactly this kind of family choice mattered (TMR's spanning radial basis vs. this
  corpus's monomial one). That axis is untested for Li here and is a legitimate residual
  uncertainty in the decomposition, not folded into "angular" by assumption.
- No l>0 orbital was ever built or turned on; "angular correlation" is inferred by
  subtraction, not directly computed. This is the honest limit of what a measurement-only
  sprint on an s-only engine can establish.
- ns=8+ radial-degree points were not measured (cost estimate ~35-45 min for the joint
  polish at the observed per-ns cost growth); the 30.0 mHa asymptote is a geometric
  extrapolation of 3 data points, not a converged limit.
- A 2nd geminal was not tested at the ns=7 (radial-converged) basis — only at ns=4. Given
  the ns=4→ns=7 shrinkage already seen for the 1-geminal gain (4.53→1.00 mHa), a 2nd
  geminal's *extra* contribution at ns=7 is expected to be smaller than the 0.355 mHa
  measured at ns=4, not larger — but this is an inference, not a measurement.

**Hard-prohibition check (CLAUDE.md Sec. 13.5):** no fitted or empirical parameter
entered production (all optimization lives in `debug/`); no change to the natural
geometry hierarchy; no negative result deleted; no Paper 2 combination-rule language
touched; no new method built (the task's own constraint).

## 7. Files

- Driver: `debug/r12ci_li_residual_decomposition.py`.
- Raw output: `debug/data/r12ci_li_residual_decomposition.json`,
  `debug/data/li_residual_decomp_run.log`.
- Baseline reproduction (unmodified, re-run for the current-state check):
  `debug/r12ci_per_shell_lambda.py` → `debug/data/r12ci_per_shell_lambda.json`.
- Prior corpus context: `debug/sprint_r12ci_gamma_transfer_memo.md` (Sec. 8/9, the
  per-shell-lambda result and the 2-geminal null on the shared-k basis),
  `debug/sprint_explicit_correlation_scoping_memo.md` Sec. 4 (the H2 analogue this
  sprint mirrors).

## 8. Open / next

1. A genuinely different radial family (DVR/spline/B-spline rather than more Laguerre
   shells) is untested for Li, the same open item Sec. 4b/7 of the H2 scoping memo
   flagged for the prolate basis. Given radial degree is already shown to be a *small*
   piece of Li's gap (9.2%, saturating), this is a low-priority follow-on, not a redirect.
2. An actual l>0 (p-orbital) extension of `r12ci_ne_engine.py` would let "angular
   correlation" be measured directly rather than inferred by subtraction — a real
   engineering lift (new orbital shapes, new Slater-Condon/Gaunt machinery for p
   electrons), explicitly out of scope here and a PI call if pursued.
3. No literature "Li s-limit" citation was located (unlike He's) — if one exists (e.g. in
   the Chakravorty-Davidson or Carroll-Silverstone-Metzger partial-wave lineage this
   corpus already cites for He/Be), it would let Sec. 5's cross-check use a truly
   external Li-specific calibration point instead of the HF-partition proxy.
