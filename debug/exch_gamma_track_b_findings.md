# Track B — the EXCHANGE two-centre ERI class `{E₁, ln, γ}`: resurgent structure

Date: 2026-08-21 | Branch: `work/sparsity-boundary` | **NO paper edits made**

Drivers: `debug/exch_gamma_census.py`, `debug/exch_gamma_borel.py`,
`debug/exch_gamma_general.py`, `debug/exch_gamma_assembly.py`,
`debug/exch_gamma_unify.py`
Logs: `debug/data/exch_gamma_{census,borel,general,assembly,unify}.txt`

Predecessor: `debug/aha_track3_findings.md` (4 objects, pattern 4/4). This track is
its named next falsifier — §4.4 there: *"the exchange class `(AB|AB)` with its `ln R`
… the first corpus object where the Borel singularity is not a pole and not a pure
algebraic branch. If its Stokes constant carries a period … the general law breaks
exactly at the point where the transcendental becomes `R`-dependent."*

---

## 0. GATE VERDICT: **GO** — pattern stands at **5/5**

The Stokes data of the exchange class is `2πi × ℚ(rates)` — **integer charges**, no
radical at all. No period appears anywhere in the connection data. The boundary
datum is **forced** by the Borel skeleton, and the forcing law is *stronger* than
the hybrid class's (it is the charge vector itself, not merely a group membership).

**The one-line statement.** *The `R`-dependent log does not break the skeleton — it
moves one Borel singularity onto the sector's own origin (`λ → 0`, the degenerate
member of the same four-point family), turning a pole into a `u ln u` branch point;
the monodromy is still `2πi ×` integer, and Euler's `γ` turns out never to exist in
the Borel plane at all — it is the bookkeeping cost of writing that origin
singularity in the `R` variable.*

---

## 1. The object, and why it is the right representative

`geovac.two_center_eri.ordered_xi_closed(p₁, p₂)` — the `τ = 0`, `σ = 0`,
`j₁ = j₂ = 0` ordered-ξ kernel (Increment 3c, build plan §8.5.3), with
`p_i = rate_i · R`. Generalisation `ordered_xi_general(τ, σ, H₁, H₂, j₁, j₂, p₁, p₂)`
(Increment 3e, §8.5.5) is analysed alongside throughout.

Both are already validated against nested quadrature by the tracked suite
(`tests/test_two_center_eri_aabb.py::test_ordered_xi_*`, 10 passed, ≤ 1.0e-13), so
the resurgence below is derived from a *validated* closed form.

**Representative check** (`exch_gamma_assembly.py`). The assembled exchange ERI is
`Σ_τ [η-half] × [ordered-ξ-half]`.

| # | statement | result |
|:--|:--|:--|
| A1 | the η half `∫₋₁¹ ηᵏ P_τ(η) e^{−qRη} dη` is **elementary** — `coeff · R^{−n} · e^{±qR}` only, no `E₁`, no `log`, no `γ` | PASS 9/9 `(k,τ)` |
| A2 | one assembled τ term keeps a γ-free Borel transform and Stokes polynomials over `ℚ(rates)`; positions are the ξ-half positions **rigidly translated** by `±q` | PASS 6/6 |

So all transcendental content of the exchange class lives in the ξ half, and the ξ
half is the right object. **Scope kept honest:** for heteronuclear centres the τ sum
is *infinite* (build plan EQ1 — it terminates iff the two centres share an
exponent). Every term is classified here; the resurgent structure of the infinite
**sum** is a separate question this track does not settle.

---

## 2. STEP 1 — exact term census, and the reduced normal form

Parser: `exch_gamma_census.parse` — the T3 parser with the `R`-dependent-log *raise*
replaced by an explicit `lnR` channel and `EulerGamma` pulled out. It still raises on
anything outside `coeff·R^k·e^{−cR}·X`, `X ∈ {1, E₁(λR), ln R, ln q, γ}`. Nothing raised.

Sweep: 6 rate pairs `(a, b) = (p₁/R, p₂/R)` including the degenerate `a = b`:
`(3/2,1), (2,5/2), (1/2,7/2), (1,1), (5/2,3/2), (7/3,4/5)`.

| # | statement | result |
|:--|:--|:--|
| L1 | every term parses inside the census | PASS 6/6 |
| L2 | **the whole object sits in ONE trans-series sector**, action `A = a + b` | PASS 6/6 |
| L3 | `coeff(γ) = coeff(ln R)` identically (the Ein bundle) | PASS 6/6 |
| L4 | the reduced normal form below is an **exact** identity (residual `0`) | PASS 6/6 |
| L5 | `κ = (2a)(2b)/(2A)` | PASS 6/6 |
| L6 | the general `(τ,σ,H,j)` kernel parses inside the same census | PASS 6/6 |
| L7 | the γ/`ln R` bundle survives general `(τ,σ,H,j)` | PASS 6/6 |

**Reduced normal form** (derived, then verified as an exact symbolic identity against
`ordered_xi_closed`; `A = a + b`, `κ = 2ab/(a+b)`):

```
F(R) = 1/(2ab R²) [ e^{−AR}(γ + ln(κR))
                    + e^{(a−b)R} E₁(2aR)
                    + e^{(b−a)R} E₁(2bR)
                    − e^{(a+b)R} E₁(2(a+b)R) ]
```

All four terms sit at the **same** exponential action `A = a+b` — the exchange class
is a **single-sector** object, unlike the hybrid class's three sectors
`{a_d, a_c, μ}`. The `E₁` rates are exactly `{2a, 2b, 2A}`.

`ln R` appears at **exactly one order** (`R^{−2}`, the leading one) and never
squared: a single weight-1 log trans-monomial (checked, all 6 pairs).

---

## 3. STEP 2 — the Borel transform, in closed form

Write `F(R) = e^{−AR} φ(R)`, `φ(R) = ∫₀^∞ e^{−ξR} 𝔅(ξ) dξ`. Then **exactly**

```
𝔅(ξ) = (1/2ab) [ −ξ ln(ξ/κ)
                 + (ξ+2a) ln((ξ+2a)/(2a))
                 + (ξ+2b) ln((ξ+2b)/(2b))
                 − (ξ+2A) ln((ξ+2A)/(2A)) ]

       = (1/2ab) Σ_j q_j [ (ξ−ξ_j) ln(ξ−ξ_j) − (−ξ_j) ln(−ξ_j) ]
```

verified by direct Laplace transform at `mp.dps = 90`:

| check | rows | worst agreement |
|:--|:--:|:--|
| L8 `φ(R)` (closed form) vs `∫₀^∞ e^{−ξR}𝔅(ξ)dξ` | 18 (6 rate pairs × `R ∈ {1, 2.5, 7}`) | **90.4 digits** |
| G5 same, generic construction from the census, general `(τ,σ,H,j)` | 18 | **57.9 digits** (dps 60) |

---

## 4. STEP 3 — Borel singularity table and Stokes data

### 4.1 The table

| local `ξ_j` | global `ζ_j = A + ξ_j` | charge `q_j` | type | `S_j/(2πi)` |
|:--|:--|:--:|:--|:--|
| `0` | `+(a+b)` | `−1` | `u ln u` log branch — **at the sector's own origin** | `−1/(2ab)` |
| `−2a` | `−(a−b)` | `+1` | `u ln u` log branch | `+1/(2ab)` |
| `−2b` | `+(a−b)` | `+1` | `u ln u` log branch | `+1/(2ab)` |
| `−2A` | `−(a+b)` | `−1` | `u ln u` log branch | `−1/(2ab)` |

`Σ_j q_j = 0` (charge neutrality). The four global positions are `±(a+b)` and
`±(a−b)` — the **sum and difference rates**, symmetric under `ζ → −ζ`.

Monodromy `u → e^{2πi}u` on `q_j u ln u` gives `Disc 𝔅 = 2πi q_j (ξ−ξ_j)/(2ab)`, i.e.
the alien derivative attaches the **single term** `2πi q_j /(2ab R²)` — the
trans-series terminates at one instanton, exactly as for the `e^a E₁(a)` seed.

### 4.2 Three independent routes, all agreeing

**Route (i) — closed form.** Read `q_j` off `𝔅` above.

**Route (ii) — BLIND large order.** Uses only the exact rational coefficient
sequence `a_m` of `φ` taken from the *parsed* closed form. Two stages, both blind:

* *type/scale*: three consecutive coefficients at `m₀ = 34` give the shift `s` and
  `λ_min` with no prior assumption (`a_{m+1}/a_m = −(m+1−s)/λ`);
* *positions and charges*: exact **Prony/Hankel over ℚ** on
  `v_m = a_m(−1)^{m−1}/(m−s)! = Σ_j q_j x_j^{m−s+1}`, minimal order auto-detected.

**Synthetic validation first** (required — log-corrected large-order fits are the
classic trap):

| synthetic model | true `s` | blind `s` | recovered `(λ_j, q_j)` | verdict |
|:--|:--:|:--:|:--|:--|
| `q(ξ+λ)ln(ξ+λ)`, `λ={2,5,9}`, `q={1, −3/7, √2}` | 3 | 3.0 | `(2,1), (5,−3/7), (9,√2)` | **OK** (irrational charge recovered exactly) |
| `q ln(ξ+λ)`, `λ={3,8}`, `q={2,−5}` | 2 | 2.0 | `(3,2), (8,−5)` | **OK** |
| `q/(ξ+λ)`, `λ={4,11}`, `q={1,7/5}` | 1 | 1.0 | `(4,1), (11,7/5)` | **OK** |

The extractor therefore distinguishes **pole (s=1) / bare log (s=2) / `u ln u`
(s=3)** blind. On the real object:

| `(a,b)` | blind `s` | blind `λ_min` | blind `(λ_j, q_j)` | Prony residual |
|:--|:--|:--|:--|:--:|
| `(3/2,1)` | 3.000272 | 1.9999845 | `(2, 1/3), (3, 1/3), (5, −1/3)` | **0** |
| `(2,5/2)` | 3.0333887 | 3.9964592 | `(4, 1/10), (5, 1/10), (9, −1/10)` | **0** |
| `(1/2,7/2)` | 3.0 | 1.0 | `(1, 2/7), (7, 2/7), (8, −2/7)` | **0** |
| `(1,1)` | 3.0 | 2.0 | `(2, 1), (4, −1/2)` (confluent) | **0** |
| `(5/2,3/2)` | 3.0000134 | 2.9999988 | `(3, 2/15), (5, 2/15), (8, −2/15)` | **0** |
| `(7/3,4/5)` | 3.0 | 1.6 | `(8/5, 15/56), (14/3, 15/56), (94/15, −15/56)` | **0** |

`s = 3` blind ⇒ the singularity is `(ξ−ξ_j)ln(ξ−ξ_j)`, **not** a pole and **not** a
bare log. (Consequence: the coefficients grow like `(m−3)!`, i.e. the exchange class
is *two factorial orders less divergent* than a pole-type object at the same action.)

**Route (iii) — the `R`-plane branch cut.** Prediction
`Disc F = −(2πi)/(2ab R²) Σ_j q_j e^{−ζ_j R}` from the table alone, tested at
`R = −x ± i·10⁻⁶⁰`, `x ∈ {1.3, 2.7}`, dps 90:

| check | rows | worst relative residual |
|:--|:--:|:--|
| L11 | 12 | **2.80e-60** (contour-offset limited) |

**Three-way agreement: ≥ 25 digits on every rate pair, not just one** — 90 digits
(route i), exact over ℚ (route ii), 60 digits (route iii).

---

## 5. STEP 4 — classification

### (a) Field of the Stokes data — **ALGEBRAIC, up to `2πi`**

`S_j = 2πi · q_j/(2ab)` with `q_j ∈ {−1,+1,+1,−1} ⊂ ℤ`. So `S_j ∈ 2πi·ℚ(rates)`;
**no radical is required at all** (the hybrid class needed `ℚ(√d)`). π-power `π^{+1}`.

At general `(τ,σ,H,j)` the charge becomes a **Stokes polynomial** `P_j(ξ)` with
`Disc 𝔅 = 2πi P_j(ξ)`; e.g. at `(a,b)=(3/2,1)`, `τ=1`:

```
ξ_j = −5 : −(ξ−1)(ξ+5)(2ξ+7)/54      ξ_j = −3 : (ξ+3)(2ξ²+9ξ−9)/54
ξ_j = −2 : (ξ+2)(2ξ²+11ξ−4)/54       ξ_j =  0 : −ξ(ξ+6)(2ξ+3)/54
```

| # | statement | result |
|:--|:--|:--|
| G4 | every `P_j` has **rational** coefficients, all 21 cases (3 rate pairs × 7 `(τ,σ,H,j)`) | PASS |
| G1c | **sum rule** `Σ_j P_j(ξ) ≡ 0` — charge neutrality generalised | PASS 21/21 |
| G3/L13 | singularity positions stay at `{0, −2a, −2b, −2A}` at every `(τ,σ,H,j)` | PASS 21/21 |

### (b) Is the boundary datum forced? — **YES, and more tightly than in the hybrid class**

Two forcings, which are one statement:

1. **`κ` is the charge-weighted product of the other singularity positions**
   `κ = Π_{j≠0} λ_j^{q_j} = (2a)(2b)/(2A)` — L5/L12, PASS 6/6, exact.
2. **`γ + ln κ` is a single forced number.** `γ` and `ln κ` enter `𝔅` only through
   the combination `γ + ln κ`, and charge neutrality `Σ_j q_j = 0` is exactly what
   removes the regular `ξ`-linear remainder, fixing it. Checked two ways per case:
   `d/dξ [𝔅 − (pure charge sum)] = 0`, and the remainder takes the *same* value at
   two rational `ξ` (exact). L12 PASS 6/6.

**Sharper still, from the Borel side (G2):** because `coeff(γ) = coeff(ln R)`,
**`γ` cancels identically from `𝔅`** — checked symbolically at 21/21 general cases.
`γ` never appears in the Borel plane. It is not boundary *data*; it is the price of
writing the origin log-branch in the `R` variable, since
`L⁻¹[R^{−m}(γ + ln R)] = ξ^{m−1}(H_{m−1} − ln ξ)/(m−1)!` is the γ-free object.

**General-`(τ,σ,H,j)` forcing (G1/G1b).** Normalise
`𝔅 = P₀(ξ)ln ξ + Σ_{j≠0} P_j(ξ) ln(1+ξ/λ_j) + Q(ξ)` (canonical: every log factor
vanishes at `ξ=0`). `Q` is then the boundary datum. Result: at **every** order of
`Q` and **every** case, `Q`'s log content is a rational multiple of the *same*
vector `(1,1,−1)` in the basis `{ln 2a, ln 2b, ln 2A}` — i.e. it lies on the
**one-dimensional ℚ-line spanned by `ln κ`**, with `κ` the `τ=0` charge vector.
G1 PASS 21/21, G1b PASS 21/21.

### (c) Cut-freeness

`𝔅(ξ)` is analytic on the whole Laplace ray `(0,∞)` — the three non-trivial
singularities sit at `ξ = −2a, −2b, −2A < 0` and the fourth is the *endpoint* `ξ=0`.
So the exchange kernel is **Borel summable along `R > 0` with no lateral
ambiguity** — verified numerically by rotating the contour to `±π/6, ±π/4`
(L12b: worst relative gap between lateral resummations **4.6e-89**).

This is a *different* mechanism from the hybrid class, where cut-freeness came from
sector constants cancelling pairwise. And it comes at a price the hybrid did not pay:

> `Disc_{arg R = π} F ≠ 0` (L11) — **the exchange kernel is genuinely multivalued in
> complex `R`**, whereas the hybrid ERI was single-valued. That is exactly what the
> `R`-dependent log buys. But the discontinuity it creates is `2πi × rational ×
> elementary` — no new transcendental crosses the cut.

---

## 6. STEP 5 — decoy controls

Per the numerical-coincidence audit rule, every structural law was re-tested with
deliberately wrong data.

**On `κ = Π λ_j^{q_j}`** (`exch_gamma_census.py`, all 6 rate pairs):

| decoy | matches `κ` |
|:--|:--:|
| `(2a)(2b)(2A)` (charge `+1` on `2A`) | 0/6 |
| `(2a)(2A)/(2b)` (charges permuted) | 0/6 |
| `(2a+2b)/2` (additive not multiplicative) | 0/6 |
| `(2a)(2b)/(2A)²` (wrong charge weight) | 0/6 |
| `√((2a)(2b))` (geometric mean) | 0/6 |

**On the charge vector** (`exch_gamma_borel.py`, all 6 pairs):

| decoy charge vector on `{2a, 2b, 2A}` | matches `κ` |
|:--|:--:|
| `all +1` | 0/6 |
| `swap sign on 2b` | 0/6 |
| `charge 2 on 2A` | 0/6 |
| `swap sign on 2a` | **1/6** |
| `charges (2, 1, −1)` | **1/6** |
| `charges (0, 1, −1)` | **1/6** |

> **Honest note on the three 1/6 hits.** All three land on the *same* rate pair,
> `(a,b) = (1/2, 7/2)`, and for the same reason: there `2a = 1`, so `1^e = 1` for
> every exponent `e` and the `2a` slot carries no information. That is a degenerate
> rate, not a law — each fails on the other five pairs. This is precisely why the
> sweep runs over six rate pairs rather than one, and why the `κ` law itself
> (§5 (b), L5/L12) is checked as an exact symbolic identity rather than a match.

**On the general-case forcing** (`exch_gamma_general.py`): `ln(11/13)`, `ln(101/2)`,
`ln 13` are all *not* in the span of `{ln λ_j}` for every rate pair (9/9 fail as
required). Note also that G1b's direction test is *computed*, not assumed — any
order whose log direction differed from `(1,1,−1)` would have been flagged, and none
was.

---

## 7. Cross-class check: is it the *same* boundary law as the hybrid? (`exch_gamma_unify.py`)

| # | statement | result |
|:--|:--|:--|
| U1 | the hybrid class's boundary log lies in the multiplicative group generated by **its** Borel singularity positions | **7/7** |
| U1b | the `Λ` extracted here (by log-content vectors — a route completely different from T3's) agrees with T3's published cross-ratio up to inversion: `10/7, 35/11, 10/7, 4/3, 95/143, 9/5, 99/19` | **7/7** |
| U2 | but the *exponent pattern* is **not** the same rule | — |

**Verdict: partial unification.** Both classes obey the weak law *"the boundary log
is a `{0,±1}` monomial in the Borel singularity positions"* — so the boundary period
is skeleton-forced in both. They do **not** obey a single exponent rule:

* hybrid: `(+1,−1,−1,+1)` on `{a_c∓a_d, μ∓a_d}` — pairs singularities **across two
  divergent sectors** by the sign of `a_d`; the log itself sits in a *third*,
  `E₁`-free sector that has no Borel singularity of its own;
* exchange: `(−1,+1,+1,−1)` on `{0, 2a, 2b, 2A}` — all four in **one** sector, and
  the exponents *are* the Stokes charges.

The exchange class's stronger statement is therefore specific to its single-sector
geometry, where the log-carrying term and the singularities live together. (U1b is
also a useful by-product: an independent reproduction of the T3 cross-ratio law.)

---

## 8. Where this leaves the pattern

| object | singularity type | Stokes data | π-power | boundary period | forced by skeleton? |
|:--|:--|:--|:--|:--|:--|
| `e^a E₁(a)` (P18 seed) | simple pole | `2πi · 1` | `π^{+1}` | none | — |
| hybrid ERI `{E₁, ln}` (P58) | pole, order `k+1` | `2πi · ℚ(√d)` | `π^{+1}` | `ln Λ`, `Λ ∈ ℚ` | yes (cross-sector monomial) |
| **exchange ERI `{E₁, ln, γ}` (P58, this track)** | **`u ln u` log branch, one of them AT the origin** | **`2πi · ℚ(rates)`, charges in ℤ** | `π^{+1}` | **`γ + ln κ`, `κ = Π λ^q ∈ ℚ`** | **yes — exponents = the Stokes charges** |
| `N(D)` one-mass fibre (P59) | √-branch | algebraic over `ℚ(ρ)` | `π^{0}` | `K(1−ρ)` (elliptic) | no |
| `T2` second cusp (P59) | `(z*−w)^{3/2}` | algebraic `× 1/√π` | `π^{−1/2}` | Γ(2) MMV — [OPEN] | no |

**5/5 algebraic up to a π-power.** The refinement this track adds to T3 §4.1's
"π-power is fixed by the singularity type": it is fixed by whether the local exponent
is **integral** (pole *or* log branch → `2πi`) or **half-integral** (`√` → `π⁰`;
`3/2` → `1/√π`). The log branch is a *new* type and it lands on the integral side,
so the rule survives with a sharper statement of what it depends on.

And on T3 §4.2's boundary-period ladder — *skeleton-forced at weight ≤ 1,
skeleton-independent from genus ≥ 1* — the exchange class is the **strongest**
weight-1 rung yet: not merely forced, but forced with the Stokes charges as literal
exponents, and with `γ` shown to be no boundary datum at all.

---

## 9. Honest limitations

1. **Derived from a validated closed form, not measured on the integral.** As in T3
   §2.8: the closed form is independently validated against nested quadrature by the
   tracked suite (≤ 1.0e-13); the Borel/Stokes work is an analysis *of that form*,
   and the large-order route is a blind cross-check of the bookkeeping, not a second
   independent determination of the ERI. Measuring exponentially small sectors off
   the quadrature is not feasible.
2. **The infinite τ sum is not classified.** Each τ term is (§1, A2); the resurgent
   structure of the heteronuclear infinite sum is a separate question.
3. **Rational rates only.** The sweep uses rational `(a,b)`, as the physical
   hydrogenic exponents are. The coefficient field is `ℚ` here; the invariant claim
   is that it is algebraic over the rate field.
4. **`ordered_xi_general` sampled, not proved.** 7 parameter tuples × 3 rate pairs
   (`τ ≤ 3`, `σ ≤ 2`, `H ≤ 2`, `j ≤ 1`). Uniform across all of them, but a general
   proof in `(τ,σ,H,j)` is not attempted.
5. **Genuinely new territory checked, but only one exponent order deep.** At general
   `(τ,σ,H,j)` the Borel branches reach order `(ξ+λ)^{13}ln(ξ+λ)`. All are
   `2πi × ℚ(rates)`; no simple poles (`k = −1`) were observed in the tested set,
   which is a fact about the sample, not a theorem.

---

## 10. Suggested next falsifier

Every corpus object tested so far has its Borel singularities at **rate-linear**
positions (sums and differences of orbital exponents), which is what keeps the
connection data rational. The two objects whose skeletons are *not* algebraic
(`N(D)`, `T2`) are exactly the ones where the singularity positions move on a curve.
The sharp next question is therefore not another chemistry class but:

> is there a corpus object whose Borel singularity **positions** are periods
> (rather than algebraic in the parameters) while its **charges** stay rational?

`T2`'s Γ(2) resurgent Lambert series (Paper 59 `sec:modular`) is the standing
candidate, and it is already the specialist hand-off — so the chemistry side of this
axis is, on present evidence, closed at algebraic.
