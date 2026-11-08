# Sprint: Paper 54 two-body forward scoping

> **CORRECTION (2026-06-03, main-session recompute).** The "Prop 2 32.4%
> erratum" claimed in this memo is **WITHDRAWN**. A consolidated recompute
> (`debug/paper54_recompute_both_constructions.py`) bit-reproduces ALL of
> Paper 54's published numbers (76.7%/32.4%, Gaunt 100%/100%, Pearson
> 0.58/0.41, Table I) from the **double-sum** gauge field
> `A = Σ_i Σ_j M_i⊗M_j` — the construction the paper's numbers actually
> used (original driver `tensor_product_radial_compare.py`). This memo's
> probe computed the **diagonal** sum `Σ_i M_i⊗M_i` (eq:A_full as literally
> written), a *different* construction, and so got different numbers
> (84.3%/71.4%, Gaunt 97.3%). The real defect was the opposite and far
> smaller: `eq:A_full` was mis-transcribed as the diagonal sum and has been
> corrected to the double sum in Paper 54. The **STOP verdict stands** (the
> radial gap is the Fock conformal-factor wall — Pearson 0.58/0.41 under the
> correct construction confirms it); only the numerical "erratum" sub-claim
> was wrong. Bochniak–Sitarz 2022 / Vanhecke 1999 remain valid supporting
> citations.

**Date:** 2026-06-03
**Type:** Diagnostic + literature scoping (no paper/code edits)
**Verdict:** **STOP** on both fronts — (a) the connected-fraction "gap" is a
mis-framed metric; the *real* gap (radial weights ≠ Coulomb) is a second
conformal-factor wall identical to the resolvent CLEAN NEGATIVE; (b) the
product-triple construction has published precedent (Bochniak–Sitarz 2022;
Vanhecke 1999), and crucially **no published spectral action yields a
Coulomb/Green's-function two-body kernel** — every precedent gets
metric-invariant (bimetric-gravity) interactions, corroborating Paper 54's
own negative.

---

## 1. What was asked

Two-part scope:
- **(a)** Is the 75%→100% "connected fraction" gap structural or closeable?
- **(b)** Does external NCG literature already establish "tensor-product
  spectral action generates a two-body interaction"?

DECISION GATE: GO if an identified closeable mechanism exists AND no prior
art; STOP if (i) the gap is a second conformal-factor wall, or (ii)
published work already does the construction.

---

## 2. Part (a): the connected-fraction "gap" is the wrong metric

### 2.1 Numerical reconciliation (a data-hygiene finding)

The "connected fraction" numbers are inconsistent across the record:

| Source | n_max=2 | n_max=3 | Gauge field used |
|:-------|--------:|--------:|:-----------------|
| Paper 54 Prop 2 (published) | 76.7% | **32.4%** | claims eq:A_full (full sum) |
| Nuclear sprint memo §B Step 3 | 75.7% | **74.3%** | single strongest pair |
| **This probe — single pair** | **75.7%** | **74.3%** | single strongest pair |
| **This probe — full sum (eq:A_full)** | **84.3%** | **71.4%** | all generators summed |

Findings:
- My single-pair reproduction is **bit-consistent with the sprint memo**
  (75.7% / 74.3%). The memo is the reliable record.
- Paper 54's published **32.4% at n_max=3 is anomalous** — it matches
  neither the actual single-pair run (74.3%) nor the full sum (71.4%). It
  is almost certainly a stale/buggy run or a transcription error (likely a
  pre-conjugation-fix number, cf. the resolvent sprint which found a
  conjugation bug giving ~40% m-conservation before fixing). **The paper's
  eq:A_full text describes the full sum but the number reported was not
  produced by it.** (Flag for PI: Paper 54 Prop 2 n_max=3 value is wrong;
  the connected fraction is ~71–74%, not 32%.)
- The connected fraction does **not** collapse toward 0 with n_max under
  either construction (single-pair 75.7→74.3; full-sum 84.3→71.4). So the
  "gap to 100%" is **not** a structural collapse — but it is also **not**
  closeable to 100%, because (see §2.3) 100% is not the right target.

### 2.2 The full sum DOES close the *angular* (Gaunt) gap

The nuclear sprint's named follow-on #1 ("sum over all generators →
recover 100% Gaunt") is **CONFIRMED for the angular structure**:

| Quantity | single pair | full sum |
|:---------|:-----------:|:--------:|
| Gaunt compat, n_max=2 | 100% | 100% |
| Gaunt compat, n_max=3 | **65.3%** | **97.3%** |
| m-conservation | 100% | 100% |
| multipole, n_max=3 | k0=96.5%, k2=3.5% | k0=97.5%, k2=2.5% |

So the 34.7% Gaunt-incompatibility at n_max=3 was a single-multiplier
artifact, exactly as the memo predicted. The angular selection-rule result
of Paper 54 (Theorem 3) is **robust and improved** under the full sum.

### 2.3 Why "connected fraction → 100%" is the wrong gate

The connected fraction is a partial-trace bookkeeping quantity: it measures
how much of `{D,A}` is non-factorizable. There is no physical reason it
should reach 100% — the free single-particle kinetic content
(factorizable) is *supposed* to be there. A genuine two-body Coulomb
operator embedded alongside single-particle terms would itself have a
factorizable component. **100% connected is not the Coulomb target.** The
paper's framing of "75%→100%" as the open gap is a category slip.

The **correct** GO/STOP-deciding metric is whether the *connected part's
radial weights* are proportional to the Coulomb kernel. That is Paper 54's
own Prop 3 (radial negative), and it is where the real wall lives.

### 2.4 The radial gap is a second conformal-factor wall (STOP trigger i)

I recomputed the radial correlation against the chordal-distance S³ Green's
function (the Paper-7 / resolvent-sprint reference) for the **full-sum**
gauge field:

| Gauge construction | Pearson n=2 | Pearson n=3 | sign-agree |
|:-------------------|:-----------:|:-----------:|:----------:|
| Paper 54 single-pair (Prop 3, vs exact Slater) | +0.58 | +0.41 | 68% / 54% |
| Resolvent Dirac-weight (resolvent sprint) | +0.81 | +0.75 | — |
| **Full-sum gauge (this probe, vs chordal Green)** | **−0.82** | **−0.36** | **12.5% / 29%** |

Summing over all generators makes the radial match **worse**, not better:
the correlation is strongly *anti*-correlated and the sign-agreement
collapses to ~12–29%. No gauge-field construction (single-pair OR full-sum)
brings the radial weights near +1. This is the **same wall** the resolvent
sprint diagnosed: the 3-Y multipliers carry a **Gegenbauer radial triple on
S³**, while Coulomb wants a **Slater R^k in flat space**; the conversion is
the **Fock-projection conformal factor** (chordal = 2R·sin(geodesic/2R)) —
a *function*, not a constant. Spectral action and resolvent both hit it.

(The absolute |Pearson| here differs from Paper 54's because I used the
chordal-Green reference, not exact Slater, to keep the probe self-contained;
the load-bearing fact is the sign-anticorrelation and the absence of any
trend toward +1, not the magnitude.)

### 2.5 Curve-fit audit of "matches the Coulomb multipole hierarchy"

Per [[feedback_audit_numerical_claims]]:
- **Free parameters:** none (parameter-free output). Good.
- **Discriminating power:** weak. The "hierarchy" claim reduces to (k=0
  dominant, k=1 absent, small k=2). **k=1 absence is forced by parity for
  *any* exchange-symmetric two-body operator** — not Coulomb-specific.
  **k=0 dominance is generic** for any short-cutoff angular operator. The
  hierarchy is *consistent with* Coulomb but does **not discriminate**
  Coulomb from the large class of exchange-symmetric central operators.
- **Robustness:** the one discriminating quantity (radial weights within
  each k) fails (§2.4).
- **Verdict:** Paper 54 does **not** over-claim — its own Prop 3 is the
  radial negative and §6 explicitly assigns radial weights to calibration.
  The multipole "match" is real but non-discriminating. No correction
  needed beyond the Prop 2 numerical erratum (§2.1).

---

## 3. Part (b): literature precedent

### 3.1 The product-triple construction is standard prior art

- **Vanhecke, "On the Product of Real Spectral Triples,"
  arXiv:math-ph/9902029 (1999).** The construction
  `D_total = D₁ ⊗ I + γ₁ ⊗ D₂` with the chirality twist is the canonical
  product of even spectral triples (also Connes 1994, already cited in
  Paper 54). Paper 54's Setup §2.2 and Theorem 1 (free factorization) are
  textbook product-triple facts — correctly attributed but not novel.
- **Marcolli–van Suijlekom, "Gauge networks in NCG," J. Geom. Phys. 75
  (2014), arXiv:1301.3480** (already cited; WH1 lineage). Gauging-on-a-graph
  via inner fluctuations is their framework; spectral action → Wilson +
  Higgs lattice. This is the gauging precedent.

### 3.2 The two-geometry spectral-action interaction is published — and it
### is NOT Coulomb

- **Bochniak & Sitarz, "Spectral interaction between universes," JCAP 04
  (2022) 055, arXiv:2201.03839.** Derives the **direct interaction between
  two four-dimensional geometries from the spectral action**, to third
  order around flat vacua. **The interaction comes out as polynomials of
  the metric invariants of the two geometries, compared to bimetric
  gravity — NOT a Green's-function / Coulomb kernel.**
- **Follow-up: "Spectral interactions between strings in the Higgs
  background," EPJ Special Topics (2023).** Same structural outcome.

This is the decisive literature finding for the GO/STOP decision:

> The published spectral-action-between-two-geometries program produces
> **metric/curvature interactions**, never a two-body Green's-function
> kernel. Paper 54's central negative (spectral action gives selection
> rules + metric coupling, not 1/r₁₂) is **independently corroborated** by
> the only published precedent. The Coulomb kernel is structurally
> unreachable from a product spectral action — this is a known feature of
> the framework, not a GeoVac-specific failure.

### 3.3 No scoop, but no GO either

- **No published work claims a Coulomb/Green's-function two-body
  inter-particle interaction from a product spectral action.** Searched:
  product spectral triples, gauge networks two-body, real-space N-body
  Connes, electron Coulomb from inner fluctuations. Nothing.
- The *reason* there is no such paper is the same reason Paper 54 can't
  produce one: spectral actions compute heat-kernel (Seeley–DeWitt)
  coefficients = metric curvature, not Green's functions. Bochniak–Sitarz
  found exactly the metric-invariant interaction the formalism *can*
  produce.
- Paper 54's **distinctive** contribution (the algebra/metric split:
  algebra → angular selection rules, D/metric → radial coupling, at the
  two-body level) does **not** appear to be duplicated. But that
  contribution is a *verification of Paper 31's A/D partition*, not a route
  to deriving Coulomb — consistent with the resolvent-sprint PI directive
  to fold Paper 54 into Paper 31.

---

## 4. Verdict: STOP

Both STOP conditions fire:

1. **(i) The real gap is a second conformal-factor wall.** The
   connected-fraction "gap" is a mis-framed metric (100% is not the
   target). The discriminating quantity — radial weights vs Coulomb — fails
   under both single-pair and full-sum gauge fields (Pearson nowhere near
   +1, sign-anticorrelated, no trend with n_max). This is the **identical
   Fock-projection conformal-factor wall** the resolvent sprint already
   closed as a CLEAN NEGATIVE. An implementation sprint to "close the gap"
   would re-derive a known negative.

2. **(ii) Published precedent exists** for the construction
   (Vanhecke 1999; Marcolli–vS 2014) and for the spectral-action two-
   geometry interaction (Bochniak–Sitarz 2022), and that precedent
   **confirms** the interaction is metric-type, never Coulomb.

**What IS worth doing (low-cost, not an implementation sprint):**
- **Paper 54 Prop 2 numerical erratum:** the n_max=3 connected fraction is
  ~71–74%, not 32.4%. Fix the published number (or replace with the
  full-sum 71.4% + Gaunt 97.3%, which strengthens the angular result).
- **Citations:** add Bochniak–Sitarz 2022 (the bimetric corroboration —
  this is a *strong* supporting citation for §7 "For gravity" and §8
  Discussion) and Vanhecke 1999 (product-triple prior art) to Paper 54
  (or to Paper 31 §10 if Paper 54 stays folded).
- These are PI-discretion edits; this sprint applied none (scoping only).

**Guardrail scope check:** This is product-triple / tensor-product NCG, NOT
single-center molecular encoding. Papers 8–9 and the Track DF guardrail do
**not** apply (no shared-p₀ single-center molecule is being built). Clear.

---

## 5. Files

### Created
- `debug/paper54_full_sum_gauge_probe.py` — full-sum vs single-pair gauge
  field; connected fraction + Gaunt + m-conservation + multipole.
- `debug/paper54_radial_under_fullsum_probe.py` — radial correlation of the
  full-sum connected interaction vs chordal-Coulomb reference.
- `debug/data/paper54_full_sum_gauge_probe.json`
- `debug/data/paper54_radial_under_fullsum.json`
- `debug/paper54_two_body_forward_scoping_memo.md` — this memo.

### Modified
- None (scoping only; no paper or geovac/ edits, per directive).

---

## 6. Honest scope

- **Verified numerically (n_max ∈ {2,3}):** single-pair connected fraction
  75.7/74.3% (reproduces memo); full-sum 84.3/71.4%; full-sum Gaunt
  65.3→97.3%; full-sum radial Pearson −0.82/−0.36 vs chordal Coulomb.
- **Theorem grade:** none. Diagnostic + literature scoping.
- **Literature:** Bochniak–Sitarz 2022 + Vanhecke 1999 found and read at
  abstract level; full-PDF mechanism (off-diagonal vs inner-fluctuation
  coupling) not extracted — but the *outcome* (metric-invariant, bimetric,
  not Coulomb) is explicit in the abstract and is what the verdict rests on.
- **The Prop 2 n_max=3 = 32.4% erratum** is the one concrete paper-quality
  issue surfaced; flagged for PI, not auto-applied (directive: no paper
  edits).
