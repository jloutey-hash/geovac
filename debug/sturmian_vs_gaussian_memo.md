# Do Coulomb Sturmians converge faster per basis function than Gaussians?

**2026-08-22, exploratory (PI-directed). Driver:
`debug/sturmian_vs_gaussian_convergence.py`. H2 at R = 1.4 a0, all-electron
FCI, exact = −1.174476 Ha. No paper claim.**

## Why this was run

It is the premise under "replace the basis wholesale". Sturmians carry the
correct nuclear cusp and the correct exponential tail; Gaussians have neither.
If that buys accuracy *per basis function*, then the corpus's closed-form
two-center integrals remove the one obstacle (integral cost) keeping
Sturmian-basis correlated methods out of production. If it does not, the
program has no foundation.

QC-1 (v4.77.0) had already measured something adjacent and got a negative, but
at H2, s-only, M=2 — small enough that it might not generalize. This widens it.

## Validation gates — all PASS

| gate | result |
|---|---|
| **V1** Gaussian pipeline vs stored literature | STO-3G H2 FCI **−1.137276** vs stored **−1.1373** (Szabo & Ostlund), Δ **0.02 mHa** |
| **V1 soft** hardcoded basis data | 6-31G −1.151679 (lit −1.1516, Δ 0.1 mHa); cc-pVDZ −1.163399 (lit −1.1636, Δ 0.2 mHa) |
| **V2** Sturmian definition | ⟨χ_nl\|1/r\|χ_n′l⟩ = 0 for n≠n′ to **~5e-4** relative, all pairs, l = 0 and 1 |
| **V2 discrimination** | a *hydrogenic* set (exponent Z/n) gives **2.1e-1** — 400× larger, so the gate genuinely separates Sturmian from hydrogenic |
| **V3** conditioning / variational | s_min > 3e-3 everywhere; E monotone in M within each family |

## Results (H2, R = 1.4)

**Gaussian family (real published contractions)**

| basis | M | E (Ha) | err |
|---|---|---|---|
| STO-3G | 2 | −1.137276 | 0.037200 |
| 6-31G | 4 | −1.151679 | 0.022797 |
| 6-31G** | 10 | −1.165153 | 0.009323 |
| cc-pVDZ | 10 | −1.163399 | 0.011077 |
| cc-pVTZ (d dropped) | 18 | −1.170853 | 0.003623 |

**Sturmian family (shared exponent k, optimized at each size)**

| shells | M | k* | E (Ha) | err |
|---|---|---|---|---|
| 1s | 2 | 1.20 | −1.147479 | 0.026997 |
| 1s+2s | 4 | 1.20 | −1.152300 | 0.022176 |
| 1s+2s+2p | 10 | 1.35 | −1.165288 | 0.009188 |
| 1s+2s+2p+3s | 12 | 1.50 | −1.168781 | 0.005695 |

## The M=2 row is void — and it is the whole lesson in miniature

At M=2 **both families are a 1s Slater function**. The Sturmian 1s *is* a
Slater 1s. There is no shape difference to measure.

Control (`STO6G_1S` scaled, ζ scanned):

| | E (Ha) |
|---|---|
| STO-3G (3 primitives, ζ = 1.24 fixed) | −1.137276 |
| STO-6G (6 primitives, ζ = 1.20 optimized) | −1.147481 |
| **Sturmian 1s (6 primitives, k = 1.20 optimized)** | **−1.147479** |

**Sturmian − STO-6G(opt) = +0.002 mHa.** Identical functions. The apparent
10.2 mHa "Sturmian advantage" at M=2 is entirely *primitive count* (6 vs 3)
plus *exponent optimization* — **zero shape-class content**.

Primitives are cheap and do not enter basis-function count. This is exactly
QC-1's finding reproduced from the opposite direction: **Slater/Sturmian shapes
buy accuracy per PRIMITIVE, which is free, and the question is whether they buy
anything per BASIS FUNCTION, which is what costs.** Letting the M=2 row stand as
a win would have repeated the error QC-1 already documented.

Informative matched points are therefore M ≥ 10, where the families genuinely
differ in composition (shared-k 1s+2s+2p vs independently tuned 2s+1p).

## Matched-M comparison

| M | Sturmian | Gaussian | diff | note |
|---|---|---|---|---|
| 2 | −1.147479 | −1.137276 | 10.2 mHa | **void** (same function; fit + ζ only) |
| 4 | −1.152300 | −1.151679 | 0.6 mHa | below chemical accuracy |
| 10 | −1.165288 | −1.165153 (6-31G**) | **0.14 mHa** | tie |
| 10 | −1.165288 | −1.163399 (cc-pVDZ) | 1.9 mHa | marginal |
| 18 | −1.170705 | −1.170853 (cc-pVTZ, d dropped) | **+0.15 mHa** | **Gaussian ahead** |

Chemical accuracy is 1.6 mHa. Excluding the void row, the two families sit
**within ~2 mHa of each other at every size**, and the difference **decays
monotonically and changes sign**:

    M=2   −10.20 mHa   (void: identical functions)
    M=4    −0.62 mHa
    M=10   −1.89 mHa (vs cc-pVDZ) / −0.14 mHa (vs 6-31G**)
    M=18   +0.15 mHa   <- Gaussian ahead

## Verdict: BORDERLINE, trending negative. The premise is not supported.

The driver printed BORDERLINE, but its arithmetic counted the M=2 row, which
the control proved void. Applying the gate correctly:

* GO required **> 1 mHa at ≥ 2 matched M with a non-shrinking advantage**.
  Among informative rows exactly **one** clears 1 mHa (M=10 vs cc-pVDZ), and
  even that depends on which M=10 Gaussian set is chosen — against 6-31G** at
  the same M the gap is 0.14 mHa. The advantage also shrinks at every step and
  inverts by M=18. **Not GO on either clause.**
* Strict STOP required the Gaussian to be equal-or-lower at *every* matched M,
  which is not literally true. So the honest label is **BORDERLINE**.

Practical reading: **there is no per-basis-function advantage to build a
program on.** QC-1 generalizes from M=2 to M=18.

## Two ways this test was CONSERVATIVE toward Sturmians

Both cut the same way, which strengthens the negative:

1. **Sturmians were given far more primitives.** Each Sturmian here concatenates
   6 Gaussians *per power term* — up to 18 primitives for a 3s — against 1–3
   for the standard Gaussian functions. Under the premise's own logic
   (primitives are free, basis functions cost) that is legitimate, and the
   Sturmians still lost at M=18.
2. **The M=18 Gaussian was a stripped design.** cc-pVTZ is built with a d
   shell; removing it to match composition leaves a set that is not an optimal
   18-function Gaussian basis. A properly designed one would likely be better.

## What Sturmians *did* do well, stated fairly

They reach parity using **one** variational parameter — the shared exponent k —
against basis sets whose every exponent and contraction coefficient has been
independently optimized over decades. That is a real structural economy. It is
just not an *accuracy-per-function* advantage, which is what the program needed.

## Where this leaves "replace the basis wholesale"

Unsupported on its central premise for H₂. Options, in order of cost:

* Re-examine whether the diatomic payoff was ever accuracy-per-function, or
  something else (the metric/conditioning levers of Paper 60, or F12).

---

# LiH follow-on: what the shared exponent costs (completed)

Driver `debug/sturmian_vs_free_lih.py`, R = 3.015 a₀, 4-electron all-electron
FCI. Compares (A) Sturmian, one shared exponent **per center**, against (B) the
same shells and shapes with every ζ optimized independently. Same renderer, same
integral engine — the only difference is whether exponents are tied.

| M | (A) Sturmian | (B) free-ζ | B − A | k_Li | k_H |
|---|---|---|---|---|---|
| 3 | −7.954707 | −7.968275 | 13.6 mHa | 2.70 | 0.80 |
| 6 | −7.973665 | −7.991028 | 17.4 mHa | 2.70 | 0.80 |
| 7 | −7.986686 | −8.008550 | 21.9 mHa | 2.70 | 0.80 |
| 8 | −8.006870 | −8.022093 | 15.2 mHa | 3.30 | 0.80 |
| 10 | −8.001743 | −8.009360 | 7.6 mHa | 3.00 | 0.80 |

**VERDICT: COSTLY** — free-ζ wins at 5 of 5 sizes, every gap above chemical
accuracy (1.6 mHa), smallest 7.6 mHa ≈ 4.8×.

**Correction to the interim reading.** From the first three rungs I reported the
cost as *growing with M*. Over the full ladder it **peaks at M=7 and declines**
(13.6 → 17.4 → 21.9 → 15.2 → 7.6). Plausible reading: more functions supply
flexibility that partly compensates for tied exponents, so the constraint may
matter less asymptotically. Costly at every size tested; not monotone.

**Caveat that limits the trend, not the verdict.** The ladder is **not nested
past M=7**: M=8 is Li{1s,2s,2p,3s}+H{1s,2s} while M=10 is
Li{1s,2s,2p,3s,3p}+H{1s} — different compositions, not supersets. So the energy
rising from M=8 to M=10 is a composition artifact, not a variational failure,
and the 15.2/7.6 values are not part of a convergence sequence. Only M = 3, 6, 7
are nested, and across those the gap does grow.

**Edge check:** k_Li ∈ {2.70, 3.00, 3.30} and k_H = 0.80 are interior to the
scanned grids (1.8–3.6 and 0.35–1.4), so the Sturmian side is not grid-limited.

**Frozen core was abandoned here.** Freezing by *orbital energy* — the intended
fix for the NaH index-freezing bug — was measured **112–130 mHa** off
all-electron, because the lowest eigenvector of the bare one-electron
Hamiltonian sees an unscreened Z=3 nucleus and is far too compact. A correct
frozen core needs a self-consistent Fock orbital. The ladder is therefore capped
at M ≤ 10 where dense 4-electron FCI is exact.

## Combined standing (H₂ + LiH)

* **H₂:** no per-function advantage over real contracted Gaussian sets; the gap
  decays and inverts by M=18.
* **LiH:** the shared exponent — the property that *makes* a set Sturmian —
  costs 7.6–21.9 mHa at every size tested.

Consistent with Papers 8–9: no shared exponent serves Z=3 and Z=1. Here even
sharing *within* a center, across n, is expensive. **"Replace the basis
wholesale" has no accuracy-per-function foundation on either system.**
