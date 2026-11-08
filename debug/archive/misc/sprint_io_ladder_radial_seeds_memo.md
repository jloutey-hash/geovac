# Sprint memo — I/O ladder accounting, Rung 2 (radial-seed compactness)

**Date:** 2026-08-17
**Branch:** work/sparsity-boundary
**Status:** DIAGNOSTIC ONLY. Driver `debug/io_ladder_radial_seeds.py`. No paper /
CHANGELOG / version edits.
**Axis under test:** Avery's device-I/O thesis — for a hyperspherical /
Coulomb-Sturmian encoding, the ANGULAR Hamiltonian structure is π-free and
generated on-device from `(n,l,m)` labels via Gaunt/3j/6j (nothing loaded); only
the 1D RADIAL factors need classical tabulation and shipping. This rung measures
the size of that radial table, sharply: the count of DISTINCT closed-form
radial-seed INSTANCES vs the naive per-entry count.

**Not the same axis as** N3b (tensor-density sparsity, dead past one centre) or
QC-1 (contracted Gaussian matches Slater at matched M, no qubit/Pauli win). This
counts distinct scalar arguments fed to `{E1, ln, exp, gamma}`, not tensor
nonzeros and not Pauli terms.

**Relation to Rung 1** (`debug/io_ladder_accounting.py`, same date): Rung 1
established the compression at SHELL-QUARTET granularity (`(n,l)` per centre,
angular `m`-degeneracy generated not loaded) and named this rung's job in its
§8: "replace `unique8(S)` with the genuine count of distinct closed-form seed
instances... and add the three-centre elliptic-Bessel-moment seed class." Both
done here.

---

## 1. The question

For the actual closed-form engine (`geovac/two_center_eri.py`, Papers 58/59, not
the composed-builder abstraction Rung 1 used), how many DISTINCT scalar
arguments does the radial machinery actually evaluate — as opposed to the naive
per-tensor-entry count?

## 2. The mechanism (grounded in the live code, not asserted)

`radial_product`'s combined decay rate for a same-centre orbital pair
`(n1,l1,m1)-(n2,l2,m2)` is `b = Z/n1 + Z/n2` — a function of `(Z,n1,n2)` ONLY.
`l`, `m` never enter. Every transcendental argument the engine builds — every
`E1(mu R)`, `ln(rate ratio)`, `exp(-lambda R)` — is built from sums of these
`Z/n` rates. So many `(n,l,m)`-quartets sharing the same `n`-tuple share the
SAME transcendental argument, differing only in the (algebraically generated,
Gaunt/Wigner) polynomial coefficients multiplying it.

**Part 0 confirms this against the live code**, not just the docstrings: two
`(AA|BB)` quartets with identical `n`'s but different `l` (2s-2s vs 2p-2p self
density) give BIT-IDENTICAL `exp`-rate sets; two HYBRID quartets with identical
`n`'s but different `m` (2p₀ vs 2p₊₁ trio) give BIT-IDENTICAL `E1`-rate sets
(`{2,4,5,7}` both times, generated from the same `(mu,ad)=(6,1)` pair). PASS.

## 3. Classes and seed keys (two-centre engine)

Same-atom-only quartets (all 4 orbitals on one centre) are out of this engine's
scope — a different, already-understood machinery (same-centre Slater `R_k`
integrals) — and are excluded from both the entry and seed counts.

| class | split | seed key |
|:--|:--|:--|
| `AA\|BB` | bra 2×A, ket 2×B (or reverse) | `(bA, bB)`, `bA=Z_A(1/n1+1/n2)` — E1-free, structurally (file header) |
| hybrid | 3 on one centre + 1 on other | `(trio_centre, mu, ad)`, `mu=Z_trio·Σ(1/n)` over 3 trio orbitals, `ad=Z_lone/n_lone` |
| exchange | 2+2, one A + one B in EACH of bra/ket | `(alpha_bra,beta_bra,alpha_ket,beta_ket)`, raw single-orbital rates, canonicalized under bra↔ket swap |

Entry count = the RAW (unreduced) `M^4`-style count, matching the convention
already established elsewhere in this repo (p58 census "dense" `M^4`,
Poly-0/Poly-2's `2401`) rather than an 8-fold-symmetry-reduced count — the more
naive of the two, and it lets the compression ratio and the water T1/T2 counts
be cross-checked directly against existing repo numbers.

## 4. PART A results — two-centre engine (Li/H, Z_A=3, Z_B=1, matches p58's own census system)

| n_max | M | M^4 (naive) | class | entries | distinct seeds | ratio |
|--:|--:|--:|:--|--:|--:|--:|
| 2 | 10 | 10,000 | aabb | 1,250 | 9 | 138.9× |
| 2 | 10 | | hybrid | 5,000 | 16 | 312.5× |
| 2 | 10 | | exchange | 2,500 | 10 | 250.0× |
| 2 | 10 | | **TOTAL** | **8,750** | **35** | **250.0×** |
| 3 | 28 | 614,656 | aabb | 76,832 | 36 | 2,134× |
| 3 | 28 | | hybrid | 307,328 | 60 | 5,122× |
| 3 | 28 | | exchange | 153,664 | 45 | 3,415× |
| 3 | 28 | | **TOTAL** | **537,824** | **141** | **3,814×** |

(Same-atom quartets excluded: 1,250 at n_max=2, 76,832 at n_max=3 — matches
`M_A^4+M_B^4` exactly, an independent hand-checkable identity.)

**Growth exponents** (M: 10→28): entries `~ M^4.00` (as expected, raw `M^4`
convention), **distinct seeds `~ M^1.35`** — comfortably sub-`M^2`, so the
compression ratio itself grows `~ M^2.65` and the win *widens* with basis size,
mirroring Rung 1's finding on the angular sector (§4b there: the win grows
without bound along the basis-richness axis).

Every count in this table is independently hand-verifiable from the
combinatorics of the {A,B}⁴ centre-assignment patterns (2 same-atom, 8 hybrid, 2
aabb, 4 exchange patterns out of 16, each carrying `M_A^{#A}·M_B^{#B}` raw
entries) — done by hand during this sprint and matched the driver's numbers
exactly, an independent check on the classification logic.

## 5. PART B results — three-centre engine (H₂O, matches Poly-0/Poly-2's geometry+basis exactly)

Cross-check: **140 T1 + 280 T2 = 420 three-centre entries out of 2401 total —
matches Poly-0/Poly-2's established figures exactly** (not re-derived on a
different convention; same STO-shape minimal basis: O 1s/2s/2p×3 + 2×H 1s).

| class | status | entries | distinct instances | ratio |
|:--|:--|--:|--:|--:|
| T1 `(XX\|YZ)` | **CLOSED**, weight-1, γ-free (v4.81.0) | 140 | 12 | 11.7× |
| T2 `(XY\|XZ)` | **OPEN**, elliptic Bessel moment (Paper 59) | 280 | 15 | 18.7× |

T1 seed key: `(mu_X, {(Y,zeta_Y,D(X,Y)), (Z,zeta_Z,D(X,Z))})`, `mu_X` the
X-pair's combined rate — structurally the hybrid+exchange rate pattern one
level up, same `(Z,n)`-only dependence, canonicalized under the `(ab|cd)=(cd|ab)`
bra↔ket symmetry.

T2 seed key: same shape, `(za,zb,zc,zd,D(X,Y),D(X,Z))` per shared-centre choice
X, canonicalized the same way.

**The ratios here are modest compared to Part A** — because this is a MINIMAL
basis (M=7, only 3 distinct ζ values on O, 1 on each H) borrowed unchanged from
Poly-0/Poly-2 for direct comparability, not a richer `n_max`-graded hydrogenic
basis. Part A's own n_max=2→3 comparison shows the ratio should grow sharply
with basis richness; this was not re-measured for the 3-centre case (cap: one 3c
system, no basis sweep), so the 3-centre growth-with-richness claim is
STRUCTURAL (same `Z/n` mechanism), not separately measured here.

## 6. The T2 caveat — stated plainly, not folded into the ratio

The 18.7× count-compression for T2 is real but **qualitatively weaker than
every other row in this memo**, and must not be read the same way:

- `c1=s(1-s)`, `c2=t(1-t)` (the elliptic curve's parameters) are UNIVERSAL
  functions of the Feynman parameters `(s,t)∈[0,1]²` — the SAME continuous
  family for every quartet in every molecule. This part is genuinely shared
  library content, like the `E1`/`ln` special functions themselves.
- What varies per quartet — the geometric distances `D1,D2` and the orbital
  "mass" content `(za,zb,zc,zd)` — DOES compress the same `(Z,n)`-only way as
  the 2-centre engine (confirmed structurally, not yet by an independent
  symbolic spot-check the way Part 0 did for the 2-centre classes).
- **BUT each T2 "instance" is not a lookup-able scalar.** Unlike `E1(mu R)`
  (one number), a T2 instance is an entire 2D-Feynman elliptic-period integral
  that must be numerically quadratured over `(s,t)` — no closed-form /
  elliptic-polylog evaluator exists yet (this is Paper 59's named open
  frontier, "THE Avery-call topic," `debug/sprint_routeC_momentum_memo.md`).
  So while the DISCRETE parameter count is small, the cost of *shipping* one
  instance is qualitatively higher than shipping a 2-centre seed. **Negative,
  stated plainly: the T2 compression claim is currently a count-of-parameters
  result, not a count-of-shippable-numbers result** — the gap between those two
  is exactly the open elliptic-polylog closed form.

## 7. Verdict against the decision gate

**GO for the two-centre engine (Part A), unambiguously.** Distinct-seed count
is 2–3 orders of magnitude below the naive entry count at both tested `n_max`,
growing sub-`M^2` (measured `M^1.35`) while entries grow `M^4` — the compression
ratio widens with basis size, the regime chemistry accuracy actually needs
(bigger `n_max` per centre).

**GO, with an explicit qualitative caveat, for the three-centre engine (Part
B).** T1 (closed, weight-1, γ-free) is a clean GO by the same mechanism as
Part A, just measured on a small basis (11.7×, expected to widen with `n_max`
the same way Part A does). T2 (elliptic, open) compresses in PARAMETER COUNT
(18.7×) but the per-instance object is not yet a shippable scalar — the
elliptic-polylog closed form is the missing piece, not merely absent from this
driver but genuinely unbuilt in the corpus. Report this as the honest limit of
Rung 2's reach on the three-centre class.

## 8. Files

**Created:** `debug/io_ladder_radial_seeds.py` (driver), this memo.
**Modified:** none (`geovac/`, papers, CHANGELOG, tests all untouched — pure
diagnostic per the task cap).

## 9. Honesty caveats / scope

1. **Cap respected:** one 2-centre system (Li/H, n_max ∈ {2,3}), one 3-centre
   system (H₂O, single basis/geometry, reused verbatim from Poly-0/Poly-2, not
   swept). BeH₂ was NOT re-run for Part B; the H₂O cross-check against
   established 140/280/420 figures is the validation instead of a second
   system.
2. **Same-atom quartets excluded from Part A** (a different, already-understood
   engine) and **1-centre/2-centre-only quartets excluded from Part B's
   3-centre classification** — both explicitly, both counted and reported, not
   silently dropped.
3. **T1/T2/exchange seed canonicalization**: all three are canonicalized under
   the `(ab|cd)=(cd|ab)` bra↔ket symmetry (sorting the two "sides" of the seed
   tuple together, not independently) — an actual bug caught and fixed during
   this sprint. Switching the entry loop from an 8-fold-symmetry-reduced form
   (one representative per physical integral) to the raw `M^4` form (§9.6)
   exposed it directly: raw entries visit BOTH bra↔ket index orderings of the
   same physical integral, and the *exchange* class's seed key was not sorted
   across that swap, so its distinct-seed count at n_max=2 measured 16 before
   the fix and 10 after (10 is the physically correct count — the two bra/ket
   orderings of each seed were being recorded as two separate tuples). The
   same latent bug was caught and fixed in the T1/T2 seed keys before they
   were ever reported; the values in this memo (12 and 15) are the post-fix,
   canonicalized counts.
4. **T2's "distinct instances" are a NECESSARY-PARAMETER characterization, not
   a proof that 15 numbers suffice to reconstruct the tensor** — see §6. This is
   the sprint's one explicit negative-flavored finding, stated per instructions.
5. **PART 0's spot-check is two examples per class (aabb, hybrid), not
   exhaustive** — chosen to demonstrate `l,m`-blindness at fixed `n`, which is
   the load-bearing mechanism; the counting in Parts A/B relies on this holding
   generally, which is what the file's own `radial_product`/`radial_poly` code
   structurally guarantees (no `l`,`m` argument ever reaches the decay-rate
   computation) — not re-derived from scratch per quartet in the driver, for
   tractability (614,656 raw quartets at n_max=3 would be infeasible to run
   through full sympy symbolic derivation).
6. **RAW (`M^4`) entry convention chosen over 8-fold-symmetry-reduced** for
   direct comparability with existing repo numbers (p58 census, Poly-0/Poly-2).
   Using the reduced convention instead would divide entry counts by roughly
   8× uniformly without changing the seed counts, i.e. the reported ratios
   would shrink by roughly that same factor — the qualitative verdict (GO,
   growing sub-`M^2`) is convention-independent; the specific ratio numbers
   are not.
