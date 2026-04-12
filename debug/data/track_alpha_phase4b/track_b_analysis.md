# Track alpha-B: Packing pi as Hopf exchange constant

**Date:** 2026-04-10
**Sprint:** Phase 4B alpha, Track alpha-B
**Status:** Mixed. Subtask 1 clarifies but does NOT close the question; Subtask 2 produces one new clean identity (Delta = 1/(B-2)).

## Target

From Paper 2:

```
K = pi (B + F - Delta) = 137.036064...
B     = 42                        # degeneracy-weighted SO(3) Casimir trace, shells 1..3
F     = zeta(2) = pi^2/6          # S^1 fiber spectral zeta at s=1
Delta = 1/40 = 1/(|lambda_3|*N(2))# S^3 boundary correction
K/pi  = 43.619934066848226...
```

The question: is the outer `pi` geometrically forced by the axioms, or is
it an empirical choice?

## Subtask 1 verdict: AMBIGUOUS / PAPER 18 WORDING IS IMPRECISE

Paper 0 derives sigma_0 = pi d_0^2 / 2 as "area per state" from:

* `N_init = 2` (Axiom 1)
* isotropic circular shell 1 of radius `d_0` (Axiom 2)
* area of shell-1 disk = pi d_0^2
* sigma_0 = (pi d_0^2) / 2 = (omega_2 / N_init) * d_0^2

where `omega_2 = pi` is the **area of the unit disk in R^2** (the
2-dimensional ball volume, not the 2-sphere surface area).

Paper 18 (Sec V, lines 740-810) already claims "the overall factor pi is the
S^2 Weyl exchange constant for K." **This wording is imprecise.** Paper 18's
own catalog (lines 422-430) lists the S^2 Weyl constant as

    omega_2 / (2 pi)^2 = pi / (4 pi^2) = 1 / (4 pi)

which is *inverse* pi, not pi. The outer factor in K is `pi^{+1}`, not
`pi^{-1}`. Therefore the naive identification

    outer pi (in K) =?= S^2 Weyl constant

is numerically wrong by a factor of 4 pi^2.

### What is actually true

The defensible identification is:

    outer pi in K = omega_2 = area of the unit disk in R^2
                 = the NUMERATOR of the S^2 Weyl constant
                 = the packing-plane prefactor in sigma_0 = (pi/2) d_0^2

i.e. the same pi that enters Paper 0's fundamental area enters Paper 2's K
as an overall factor. This is a **structural** match: both are
"pi-type-2-area" constants — the ball volume omega_d at d=2. It is NOT the
Weyl *count-density* constant 1/(4 pi).

### Is the outer pi "forced"?

Partially. Paper 0's axioms force a pi to appear in sigma_0 via the
choice of circular shells (Axiom 2) and the ratio of annular area to
state count. That pi is `omega_2`. If one accepts the packing
construction as canonical and then asks "what geometric prefactor
converts the discrete Casimir combination B + F - Delta into a physical
coupling strength," the only pi^1 prefactor naturally available from
the construction is the same `omega_2 = pi` that already appears in
sigma_0.

What is NOT shown (and remains open):

1. Why exactly `pi^1` rather than `pi^2`, `1/pi`, or `pi/2`. The
   combinatorial search in Paper 2 (Sec "Statistical Validation")
   uses 14 pi-prefactors and pi^1 is just one of them.
2. Why the prefactor applies to `(B + F - Delta)` rather than to
   individual summands. The R^2 packing picture does not force
   the exchange to be on the combined sum.
3. Whether the `/2` factor (i.e., why sigma_0 = pi d_0^2 / 2 rather
   than pi d_0^2) should appear or be absorbed in K's definition.

**Net:** `outer pi` being a "pi^1 geometric exchange" is defensible.
Identifying it with the S^2 Weyl constant, as Paper 18 currently does,
is **numerically incorrect**. Identifying it with `omega_2` (unit disk
area) is consistent with both Paper 0's sigma_0 and the R^2 packing
picture. This is a refinement of the current Paper 18 wording, not a
derivation.

**Action flagged for plan mode:** Paper 18 Sec V (lines 764-767, 804)
should be corrected: "the S^2 Weyl exchange constant" -> "the 2D ball
volume omega_2 = pi, which is the numerator of the S^2 Weyl constant
and the same pi that enters Paper 0's fundamental area sigma_0."

## Subtask 2 verdict: PARTIAL POSITIVE

### B = 42 is algebraically identified (already in Paper 2)

```
B = sum_{n=1}^{3} sum_{l=0}^{n-1} (2l+1) l(l+1)
  = sum_{n=1}^{3} n^2 (n^2 - 1) / 2
  = sum_{n=1}^{3} dim(V_n) * c_2^{SO(4)}(V_n)
  = 0 + 6 + 36 = 42.
```

This is a **degeneracy-weighted SO(4) Casimir trace** over the first
three Peter-Weyl shells of S^3. B is not an empirical input: it is
forced once the selection principle `B/N = dim(S^3) = 3` is imposed,
which by Paper 2 Sec 3 has the unique positive-integer solution
`n_max = 3`. Two structural layers determine B:

1. The selection rule `n_max = 3` (algebraic, derived from
   `(m+4)(m-3) = 0`).
2. The Casimir trace identity (representation theory on S^3).

So B is **not empirical**; it is the Peter-Weyl Casimir trace through
the unique cutoff selected by the second algebraic principle. This
was already stated in Paper 2 and Paper 18; we add nothing here.

### Delta = 1/40 is a clean Casimir-trace inverse (NEW observation)

The computation reveals:

```
Delta = 1/40 = 1/(|lambda_3| * N(2)) = 1/(8 * 5)       # Paper 2 form
             = 1/(T_c2 - 2) = 1/(B - 2) = 1/40         # NEW
             = 1/(B - dim(S^3) + 1)                     # = 1/(42 - 2)
```

The two forms are equivalent because

    |lambda_3| * N(2) = 8 * 5 = 40
    B - 2             = 42 - 2 = 40

i.e. the boundary correction denominator happens to equal B - 2. The
`-2` coincides with `N_init = 2` from Paper 0 Axiom 1 (the initial two
states on shell 1, whose Casimir contribution is zero). This suggests
the reading:

    Delta = 1/(B - N_init)

but with only one data point (n_max = 3) we cannot distinguish this from
the Paper 2 form 1/(|lambda_3| * N(2)). The match is a single-point
coincidence, not a proof. It is classified as a **near-miss
identification** pending additional structural data.

### F = pi^2 / 6 remains an unexplained injection

The computation:

    F = zeta(2) = pi^2/6

is not produced by any Casimir trace on S^3. It comes from the S^1 fiber
spectral zeta function. In the Peter-Weyl picture of `S^1 -> S^3 -> S^2`
the fiber S^1 contributes its own spectral data that is independent of the
base Casimir trace. No truncation of Casimir powers through shells 1..3
reproduces pi^2/6: the closest integer traces give 42, 84, 153, 589.5
(Table 1 below) — none of these have a pi^2/6 signature.

F is therefore the **single transcendental injection** in K that is not
determined by packing/Casimir data. It enters as a separate Weyl-Selberg
exchange constant for the S^1 fiber.

### Near-miss table (K/pi = 43.6199340668...)

| Candidate                      | Value                 | rel diff   |
|--------------------------------|-----------------------|------------|
| 42 + pi^2/6 - 1/40             | 43.6199340668...      | 0          |
| 42 + pi^2/6 - 1/50             | 43.6249340668...      | 1.1e-4     |
| 42 + pi^2/6 - 1/30             | 43.6116007335...      | 1.9e-4     |
| 42 + pi^2/6 - 1/60             | 43.6282674001...      | 1.9e-4     |
| 42 + pi^2/6 - 1/20             | 43.5949340668...      | 5.7e-4     |
| 42 + zeta(2) (no Delta)        | 43.6449340668...      | 5.7e-4     |
| 42 + pi^2/6 - 1/10             | 43.5449340668...      | 1.7e-3     |

No non-Paper-2 candidate comes within 1e-5. The Paper 2 form is
numerically unique among the near-miss candidates.

### Casimir traces over shells 1..3

| Trace                                 | Value |
|---------------------------------------|-------|
| T_c2    = sum n^2 c_2(n)              | 42    |
| T_deg   = sum n^2                     | 14    |
| T_|lam| = sum n^2 |lambda_n| = 2 T_c2 | 84    |
| T_c22   = sum n^2 c_2(n)^2            | 153   |
| T_c23   = sum n^2 c_2(n)^3            | 589.5 |

The only trace that matches a component of K/pi is T_c2 = 42 = B. No
trace of higher Casimir powers reproduces pi^2/6 or 43.62.

## Subtask 3: Positive results to flag for plan-mode review

1. **Paper 18 wording correction (high-confidence).** The claim that
   the outer pi is "the S^2 Weyl exchange constant" is numerically
   wrong (1/(4 pi) != pi). The correct statement is: outer pi =
   omega_2, the 2D ball volume, which enters as the numerator of the
   S^2 Weyl constant and is the same pi that fixes sigma_0 = (pi/2) d_0^2
   in Paper 0. This is a documentation fix, not new physics. Plan mode
   should review and apply.

2. **New near-identity Delta = 1/(B - 2) (low-confidence).** The
   boundary correction 1/40 admits the alternative form 1/(B - N_init)
   with N_init = 2. This collapses Delta into the same Casimir-trace
   data as B. It is a single-point coincidence pending a second data
   point (which does not exist at n_max = 3 by the selection rule).
   Plan mode should decide whether this is worth mentioning in Paper 2
   as an observation.

3. **F = pi^2/6 remains a SECOND independent transcendental
   injection.** The derivation program for K must account for at least
   two independent pi-type transcendentals: `omega_2 = pi` (from
   packing / S^2 base) and `zeta(2) = pi^2/6` (from S^1 fiber). These
   are not reducible to each other through any Casimir-trace identity
   tested here. The naive "everything is one S^2 Weyl constant" story
   is therefore inadequate.

## Open questions

1. Is there a single spectral-geometric object on the Hopf bundle
   whose trace yields the **entire** combination B + F - Delta (rather
   than three separate pieces glued by an outer pi)? Heat kernel,
   analytic torsion, and spectral determinant attempts have all
   failed (Paper 2 Sec VII); the current finding does not change that.
2. Is the `-2 = -N_init` in Delta = 1/(B - 2) principled, or does the
   `|lambda_3| * N(2)` form have independent meaning?
3. Why does the *combination* (B + F - Delta) carry a single outer pi
   rather than three independent pi-factors (one per summand)?

## Recommendation

Pi is *partially* forced: Paper 0's packing plane forces a pi^1 of
type `omega_2` to exist somewhere in any coupling-strength exchange
built on this construction. That pi matches the outer pi in K
*structurally* but not via a rigorous chain of equalities. The sprint
target "prove outer pi is the packing exchange constant" is NOT
achieved; what is achieved is a refinement of the Paper 18 wording
and one new minor identity for Delta.

The three historic derivation attempts for K (spectral determinants,
analytic torsion, heat kernels, APS — Paper 2 Sec VII) remain negative.
No new derivation pathway is suggested by this investigation.

**Verdict: document the Paper 18 wording fix and the Delta =
1/(B-2) observation; do not modify Paper 2; do not claim the outer pi
is derived.**
