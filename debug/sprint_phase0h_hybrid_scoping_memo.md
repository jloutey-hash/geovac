# Phase 0-h — scoping the hybrid class (AA|AB), (AB|BB)

**Date:** 2026-08-11 · branch `work/sparsity-boundary`
**Driver:** `debug/phase0h_hybrid_scoping.py`
**Called for by:** build plan §8.3 ("need their own scoping pass before any code")

## Verdict: **GO**, with the seed named in advance and one gate re-priced

The hybrid class is a **generalization of machinery already in hand**, not the
from-scratch Ruedenberg Part II build §8.3 assumed. It is elementary when the
one-center pair is s-type and carries the **already-classified Stieltjes seed**
otherwise. No new transcendental class appears.

This does **not** extend to the exchange class (AB|AB) — see "What this does not
cover".

---

## HQ1 — the reduction holds, and §8.3's stated obstruction does not apply

Hybrid means one distribution is one-center and the other is two-center:

    rho_1 = conj(chi_a^A) chi_b^A          one-center at A
    rho_2 = conj(chi_c^A) chi_d^B          two-center overlap density

Because `rho_1` is one-center, its potential is closed-form (increment 1), so

    (ab|cd) = int rho_2 V_1 d3r
            = sum g_1 * gaunt * int d3r F(r_A) Y_{L'M'}(Om_A) G(r_B) Y_{ld md}(Om_B)

with `F = conj(R_c) * V_L` and `G = R_d`. That is **exactly the master integral
increment 1c already evaluates**, only with a more general A-side radial
function (1c's was a shell potential). The angular sum is finite: two nested
terminating Gaunt couplings.

Verified two ways on `(1s_A 1s_A | 1s_A 1s_B)`, Z_A=3, Z_B=1, R=3:

| reference | result |
|---|---|
| pointwise evaluation bypassing the Gaunt re-coupling | **2.8e-17** |
| `eri_md`, 12-Gaussian fit | 2.3e-07 |

**Correction to §8.3's premise.** §8.3 blocked this class on the two-center
density's angular content "not terminating about either nucleus." That is true
of the density but **irrelevant**, because the reduction never expands it about
either nucleus: `chi_c` stays on A, `chi_d` stays on B, and the (r_A, r_B)
domain absorbs the rest with no angular expansion at all. The non-termination is
avoided by not performing the expansion, not by defeating it. This is why the
class is cheaper than the plan assumed.

---

## HQ2 — power counting, and the exact rule

The r_A integral runs over `[|r_B - R|, r_B + R]`; `int r^p e^{-a r} dr` is
elementary for `p >= 0` and carries E_1 for `p <= -1`. Measured minimum power:

| l_a | l_b | l_c | l_d | min power of r_A | r_A integral |
|:---:|:---:|:---:|:---:|:---:|---|
| 0 | 0 | 0 | 0 | +0 | elementary |
| 0 | 0 | 1 | 1 | +0 | elementary |
| 0 | 0 | 2 | 2 | +0 | elementary |
| 1 | 0 | 1 | 0 | −2 | E_1 |
| 1 | 1 | 0 | 0 | −4 | E_1 |
| 1 | 1 | 1 | 1 | −4 | E_1 |
| 2 | 1 | 1 | 0 | −6 | E_1 |
| 2 | 2 | 0 | 0 | −8 | E_1 |

**Rule, confirmed by every row:** `min power = l_c − L − L'` with `L <= l_a+l_b`
and `L' <= l_c + L`, so the floor is

    min power of r_A  =  -2 (l_a + l_b)

**It depends only on the one-center pair.** `l_c` and `l_d` cancel out — the
two-center side is irrelevant to whether the seed appears. The s-type case is
exactly marginal (`p = 0`, not comfortably positive), which is why it is
elementary but only just.

Note what does *not* happen here. In (AA|BB) the seed was blocked because a
one-center orbital product starts at `k = l1 + l2` while Gaunt caps `L` at the
same value. The hybrid class has no such protection: the r_A-dependence is
`R_c * V_L`, and `R_c` starts at `l_c`, which is the wrong quantity to fight
`-(L+1)`.

---

## HQ3 — the seed is the known one, and the integral closes on it

Probe `(2p0_A 2p0_A | 1s_A 1s_B)`, Z_A=3, Z_B=1. E_1 survives the r_A integral,
at decay rates **{3, 6} = {a, a+b}** — the orbital exponent of `chi_c`, and that
plus the decay of `rho_1`. Sums of orbital exponents, nothing else.

The endpoints are `|r_B − R|` and `r_B + R`, so the outer integral must be
checked too. By parts:

    int_0^inf e^{-ct} E_1(a(t+R)) dt  =  E_1(aR)/c  -  e^{cR} E_1((c+a)R)/c

verified numerically at three (c, a, R) points, **worst deviation 2.3e-18**.

So the class closes on elementary terms plus `E_1(lambda R)` constants with
`lambda` a sum of orbital exponents — i.e. `e^{+-a} E_1(a * shift)`, the
**Stieltjes seed of Phase 0 Q2 / Paper 18 "Level 2"**.

> **CORRECTED 2026-08-11 by the increment-2 pre-build diagnostic
> (`debug/inc2_prebuild_diagnostic.py`). The seed set above is INCOMPLETE.**
>
> The r_A range is `[|r_B - R|, r_B + R]`, so there are TWO endpoints. This leg
> checked only `r_B + R`, which is indeed E_1-closed. The other endpoint,
> `|r_B - R|`, passes through **zero** at the coincidence `r_B = R`, where E_1 is
> logarithmically singular — and that endpoint is not E_1-closed:
>
>     int_0^R  e^{-cu} E_1(au) du = (1/c)[ln((a+c)/a) + E_1((a+c)R)
>                                          - e^{-cR} E_1(aR)]
>     int_0^inf e^{-ct} E_1(at) dt = ln((a+c)/a)/c
>
> both verified to ~5e-16. Euler gamma cancels, but **a logarithm survives**.
>
> Corrected seed set: **{E_1(lambda R)} ∪ {ln(rate ratio)}**. The log's argument
> is a ratio of decay rates and is **R-independent**, which still distinguishes
> it from the exchange class's `ln a` (argument scales with R) and from that
> class's explicit gamma. The monotone-growth reading across classes survives in
> refined form; the claim "hybrid carries E_1 alone" does not.
>
> Also settled by the same diagnostic: the E_1 coefficients **survive the sum
> over (L, L')** — checked on three quartets — so they are not a per-term
> artifact and the builder must carry them rather than simplify them away.

A pleasing consistency: Phase 0 Q2's seed prediction was right all along. It was
increment 1 that put it in the wrong class — (AA|BB) never had it. It belongs
here. The E_1 branch of `upper_integral`, built during the pre-1c cleanup and
pinned by a test that had no live consumer, turns out to have been built for
this class.

---

## Byproduct: Gate (d) is mis-priced for this class, by ~70x

The plan's Gate (d) says `eri_md` agreement "to ~1e-6 is the expectation."
That holds for (AA|BB), where both densities are one-center. It is **wrong for
hybrids**, because a two-center overlap density samples the exponential tail
between the nuclei — exactly where a Gaussian fit is worst:

| n_gauss | ⟨fit\|STO⟩ | deviation from the exact reduction |
|:---:|---|---|
| 6 (default) | 0.999999381 | **6.9e-05** |
| 8 | 0.999999973 | 1.4e-06 |
| 10 | 0.999999998 | 4.7e-07 |
| 12 | 1.000000000 | 2.3e-07 |

At the default 6 Gaussians the gate is ~70x looser than advertised. Anyone
refereeing increment 2 against it would either reject a correct implementation
or accept a 1e-5-level error. **Use `n_gauss >= 10` for classes 2 and 3, and
treat `eri_md` as a coarse gate there, not a precision one.**

This was nearly a false negative in this very sprint: HQ1 first read as a FAILED
reduction at 6.9e-05 before the fit sweep and the pointwise route located the
error in the reference rather than in the derivation.

---

## What this does not cover

**The exchange class (AB|AB) is untouched by this result.** The whole reduction
rests on `rho_1` being one-center so that `V_1` is closed-form. In the exchange
class both distributions are two-center, neither has a closed-form potential,
and the reduction simply does not start. That is the genuine Ruedenberg Part II
problem and it needs its own scoping pass. Do not read this memo as scoping it.

Weights: hybrid and exchange together are ~80% of the census tensor; §8.1 does
not split that figure between them.

---

## What increment 2 would actually be

Reused, already validated: `angular_factor` (1c, checked to 4e-14 against direct
angular quadrature), the (r_A, r_B) domain and region logic (1c), `V_L_radial`
and `multipole_decomposition` (increment 1), and `upper_integral` **including
its E_1 branch** (pre-1c cleanup).

Genuinely new: an r_A integrator that admits negative powers and carries E_1
endpoint terms, and an r_B integrator that consumes E_1 integrands via the
closure identity verified above.

Not reusable: 1c's shell-kernel trick. It works because the A-side is a shell
potential, which makes the r_A integral the *angular* one (all-even powers, so
no logarithm). In the hybrid class the r_A integral is a genuine radial
integral, so that particular protection is gone. The 1c lesson still applies in
its general form — prefer a formulation that never generates spurious
transcendentals over one that generates and cancels them — but the specific
device does not carry over.

Open before coding, in order:
1. Do the E_1 coefficients survive the *sum* over (L, L'), or cancel as they did
   in 1c? HQ3 measured survival per-term on one quartet; per-sum cancellation
   across the whole term set was not tested.
2. Is there a formulation that keeps the r_A powers non-negative? Worth one
   pass before accepting the seed, given how 1c went.

---

## Files

- `debug/phase0h_hybrid_scoping.py` — the three probes
- build plan §8.3 → replaced by §8.4 with this result
