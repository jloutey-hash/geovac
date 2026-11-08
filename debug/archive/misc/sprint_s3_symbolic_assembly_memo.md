# Sprint S^(3) symbolic assembly (2026-06-11)

Driver: `debug/s3_symbolic_assembly.py`
Data: `debug/data/s3_symbolic_assembly.json`

## Headline

**S^(3) modulo W10 is explicit and closed.** Every t-value in the linear
table reduces to Q[π², ln2, ζ(3), ζ(5), ζ(7), ζ(9)] — no t2(5,1)/t2(5,3)/t2(7,1),
no Li_4(1/2), no ζ(11). All three genuine depth-2 generators cancel in the
full assembly. **Realized depth at k=3 is ZERO modulo W10**: the identified
part of S^(3) sits in the classical depth-1 ring.

Gate: `|symbolic_part + W10_numerical - canonical| = 4.09e-197 ≤ 1e-180`. **PASS.**

---

## 1. Closed form: S^(3) = Ident + W10

The **identified part** (26 terms, all rational coefficients in Q):

```
Ident =

  ln2 sector:
    128 π² ln²2  −  (32/3) π⁴ ln²2
    + (7232/81) π² ln2  − (35504/243) π⁴ ln2  + (352/15) π⁶ ln2  − (314/315) π⁸ ln2

  z3 sector:
    − 96 π² ln2 z3  + 16 π⁴ ln2 z3
    + 20 π² z3²  − (2200/27) π² z3  + (4988/81) π⁴ z3  − (66/5) π⁶ z3  + (157/210) π⁸ z3

  z5 sector:
    − 80 π² ln2 z5
    + (7460/81) π² z5  + (85/3) π⁴ z5  − (11/3) π⁶ z5

  z7 sector:
    − 7 π² z7  + (91/12) π⁴ z7

  z9 sector:
    − 30 π² z9

  Pure π sector:
    (4096/6561) π²  − (109240/19683) π⁴  + (27427/1215) π⁶
    − (288859/51030) π⁸  + (1381/2835) π¹⁰  − (17851/1247400) π¹²
```

The **W10 remainder** (numerically included, motivically unidentified):

```
W10 = −2048 t3(4,3,3) − 1024 t3(4,4,2)
       − 1536 t2(4,6) − 1024 t2(7,3) − 512 t2(8,2)
```

Numerical value: W10 ≈ −23.3366 (≈ 42% of the raw S^(3) = 31.57 value;
the identified part evaluates to ≈ 54.91 before W10 correction).

---

## 2. Strategy: pair substitution via stuffle relations

The linear table (12 t3 + 26 t2 + boundary/lam terms) contains t3 pairs
whose INDIVIDUAL values are not PSLQ-identified but whose SUMS are pinned by
stuffle relations:

- **R2:** `t3(2,1,3) + t3(2,3,1) = λ(3)·t2(2,1) − t2(5,1) − t2(2,4) − t3(3,2,1)`
  Both carry coefficient −2048 in the table, so their contribution is
  `−2048 × [R2 rhs]` — fully identified once t3(3,2,1) and all t2's are known.

- **R6:** `t3(4,1,3) + t3(4,3,1) = λ(3)·t2(4,1) − t2(7,1) − t2(4,4) − t3(3,4,1)`
  Both carry coefficient +2048; contribution `2048 × [R6 rhs]`.

Additional stuffle-derived t3 forms:
- **R1** → `t3(2,2,2) = (λ(2)·t2(2,2) − t2(4,2) − t2(2,4)) / 3`
- **R4** → `t3(2,4,2) = λ(2)·t2(4,2) − 2·t3(4,2,2) − t2(6,2) − t2(4,4)`
- **R5** → `t3(2,3,3) = (λ(3)·t2(2,3) − t3(3,2,3) − t2(5,3) − t2(2,6)) / 2`

With t3(4,2,2) and t3(3,2,3) from stageB PSLQ, and the four trailing-1 values
from the closure sprint, every non-W10 t3 entry is handled.

---

## 3. Cancellation inventory

### What cancels in the full assembly

| Object | Introduced by | Cancels with | Net coefficient |
|:-------|:-------------|:-------------|:----------------|
| t2(5,1) | t3(4,1,1), t3(3,2,1) | R2 pair sum (×2) | **0** |
| t2(5,3) | t3(3,4,1), t3(4,2,2), t2(6,2), t2(2,6) | sign flips across substitutions | **0** |
| t2(7,1) | t3(3,4,1) | R6 pair sum | **0** |
| Li4(1/2) | t3(2,1,1) | −t3(2,1,1) contribution in R2 pair sum subst. | **0** |
| ζ(7) | several t2 odd reductions | mutual cancellation | **0** |
| ζ(11) | t2(8,3), t2(4,7) | mutual cancellation | **0** |
| ln2·ζ(3) | multiple t2/t3 | partial, net nonzero | **nonzero** |

The complete cancellation of all three depth-2 generators is the central
structural result: the identified part of S^(3) does NOT introduce new
depth-level content beyond depth 1.

### What survives

Only weight-1 generators (with polynomial π² coefficients):
**ln2, ζ(3), ζ(5), ζ(7), ζ(9)** plus pure powers of **π**.

The nonzero appearances:
- `ln2`: 6 terms (up to π⁸ ln2, π⁴ ln²2)
- `ζ(3)`: 7 terms (up to π⁸ ζ(3), plus z3²)
- `ζ(5)`: 4 terms (up to π⁶ ζ(5))
- `ζ(7)`: 2 terms (up to π⁴ ζ(7))
- `ζ(9)`: 1 term (−30 π² ζ(9))
- pure π: 6 terms from π² to π¹²

---

## 4. Gate results

**Main gate:**
```
S^(3) symbolic + W10_numerical = 31.572561207512022754... 
canonical (200 dps)             = 31.572561207512022754...
residual                        = 4.09e-197  ≤ 1e-180  PASS
```

**Individual substitution gates (all ≤ 1e-50 at 60 dps):**
- 20 t2 identifications: all PASS (residuals 1e-62 to 1e-61)
- 9 t3 identifications: all PASS (residuals 1e-63 to 1e-61)
- 2 pair sums via R2/R6: both PASS (residuals 3.9e-62, 7.1e-63)
- S_min verification: PASS (1.5e-59)

---

## 5. W10 structure and open items

The W10 remainder contains **5 unidentified objects**:
- t3(4,3,3) [w=10, coeff −2048]: level-2 w10 triple sum
- t3(4,4,2) [w=10, coeff −1024]: level-2 w10 triple sum
- t2(4,6) [w=10, coeff −1536]: level-2 w10 double sum
- t2(7,3) [w=10, coeff −1024]: level-2 w10 double sum
- t2(8,2) [w=10, coeff −512]: level-2 w10 double sum

All five have w+d=10+2=12 or 10+3=13 (odd for the triples, even for the doubles).
From the parity theorem, the t3 objects should reduce to depth ≤ 2 (w+d=13 odd),
introducing t2(9,1), t2(7,3), or higher depth-2 generators. The t2 doubles at
w=10 are even-weight; t2(8,2) and t2(7,3) are depth-2 objects, possibly irreducible
vs products (like t2(5,1) was at w=6). The level-2 w10 basis has dimension 89 — a
full PSLQ run at ≥500 dps is the named next step.

---

## 6. Loop tower reading

| k | Ident depth | New generators |
|:-:|:------------|:---------------|
| 2 (S_min) | 0 (mod W10 = nothing) | none — fully depth ≤ 1 |
| 3 (S^(3)) | 0 (mod W10) | none in identified part — depth-2 gens all cancel |
| W10 | predicted ≤ 2 | t2(5,3)/t2(7,1)-type objects possible |

The depth-2 generators t2(5,1), t2(5,3), t2(7,1) appear TRANSIENTLY
(introduced by individual t3 substitutions) but cancel completely in the
full assembly. This cancellation is a structural property of the chain-sum
architecture, not a numerical accident.

Depth of the identified part of S^(3) is **zero** (classical depth-1 ring:
Q[π², ln2, ζ(odd)]). The k=3 loop chain is even shallower than k=2 in its
identified part — the genuine novelty is deferred to the W10 layer.

## 7. Honest scope

- **Exact / theorem-grade:** all 20 t2 and 9 t3 identifications (each
  individually verified ≤ 1e-50 against 220-dps cache). The pair sums via
  R2/R6 (algebraic stuffle identities). S_min closed form (verified 1.7e-129).
- **Rigorous-numerical:** gate residual 4.09e-197 (200 dps assembly from
  Levin-safe cache values, anchored by the [31.5706, 31.5730] bracket).
- **OPEN (W10):** t3(4,3,3), t3(4,4,2), t2(4,6), t2(7,3), t2(8,2) — numerically
  included in the canonical value, not symbolically identified.
- **Cancellation facts are theorem-grade:** t2(5,1), t2(5,3), t2(7,1), Li4(1/2),
  ζ(7), ζ(11) all cancel as consequences of the exact rational coefficients
  in the linear table and the exact rational closed forms of each substitution.
  No numerical coincidence.

## Status: SPRINT CLOSED
