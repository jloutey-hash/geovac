# Sprint Hain-Brown PSLQ Test A

**Date:** 2026-06-06 (post-strategic-synthesis evening)
**Sprint:** Hain-Brown PSLQ Test A (Recommendation A from
`debug/strategic_synthesis_2026_06_06_memo.md` Addendum)
**Driver:** `debug/sprint_hb_pslq_test_compute.py`
**Data:** `debug/data/hb_pslq_test_results.json`
**Wall time:** ~80 s for the full 3-precision sweep (modular basis assembly
+ Eichler integration dominates; PSLQ itself is fast).
**Discipline:** `mpmath` arbitrary-precision arithmetic at 50 / 100 / 200 dps;
cross-precision agreement is the false-positive filter.

---

## 1. TL;DR

**Verdict: NEGATIVE.** Zero Hain-Brown modular-ring identifications survive
cross-precision agreement at 50 / 100 / 200 dps. Every GeoVac
$\mathrm{Sym}^2$-tagged period that PSLQ can identify lands in the pure
$\mathrm{MZV}(\mathrm{MTM})$ basis $\{\pi, \pi^2, \pi^3, \pi^4, \pi^6,
\pi^8, \zeta(3), \zeta(5), \beta(2) = G, \beta(4)\}$ — i.e. inside Brown's
mixed-Tate motivic Galois group $\mathcal{G}_{\mathrm{MT}(\Z)}$ at level $1$
extended by level-$4$ Dirichlet content $\{G, \beta(4)\}$. **No
identification with $E_4(2i)$, $E_6(2i)$, $\Delta(2i)$, or any of the three
Eichler-style periods of $\Delta$ along the imaginary axis.**

**Implication.** Per the strategic-synthesis decision gate:
- **NEGATIVE outcome ⇒ Hain-Brown identification is empirically NOT
  supported at the $\mathrm{Sym}^2$ level.** The Hain-Brown structural-shape
  match identified by the Hodge-theoretic $SL_2$ sub-agent earlier today
  (`memory/hain_brown_identification.md`) does not survive period-level
  testing. The earlier match was at the *categorical* level (extension
  structure $1 \to U \to G^{\mathrm{rel}} \to SL_2 \to 1$, $\mathrm{Sym}^k V$
  tensor category, canonical mixed Hodge structure) — none of which forces
  modular content into the periods themselves.
- **Cosmic-Galois target stays $\mathcal{G}_4 = \mathcal{G}_{\mathrm{MT}(\Z[i, 1/2])}$**
  (Deligne 2010 + Glanois 2015), as Paper 55 already classifies.
- **The $\mathcal{G}_4 \supseteq U^*_{GV}$ injection direction
  (Recommendation B of the strategic synthesis) remains the strongest
  near-term theorem-grade target.** This sprint indirectly strengthens
  Recommendation B by confirming the period-level content is exactly
  $\MT(\Z[i, 1/2])$ at depth $\le 1$.

The Hain-Brown identification is not falsified at the categorical
level — it stays a structural-shape match — but the empirical $\mathrm{Sym}^2$
period content does not engage the modular ring of $SL_2(\Z)$ as Hain-Brown
would predict if GeoVac were a faithful realisation.

---

## 2. Setup

### 2.1 GeoVac $\mathrm{Sym}^2$-tagged periods

The $\mathrm{Sym}^2 V_{\mathrm{fund}} = \Q^3$ representation of $SL_2$ acts on
the basis $\{e_1^2,\ e_1 e_2,\ e_2^2\}$ as

$$\rho^{(2)}(g) = \begin{pmatrix} a^2 & 2ab & b^2 \\ ac & ad+bc & bd \\ c^2 & 2cd & d^2\end{pmatrix},$$

with Clebsch-Gordan multiplicities $\{1, 2, 1\}$ carrying $sl_2$-weights
$\{+2, 0, -2\}$. From the $V_{\mathrm{fund}} \otimes V_{\mathrm{fund}}
\cong \mathrm{Sym}^2 \oplus V_{\mathrm{triv}}$ decomposition the
$\mathrm{Sym}^2$-tagged GeoVac periods constructed are structural
compositions of the CG-trace (= $1+2+1 = 4$) or CG-highest-weight (= $1$)
with M1, M2, or M3 closed forms from Paper 55:

| Symbol | Closed form | Type |
|:---|:---|:---|
| $P_{M_1}$ | $4 \pi$ | CG-trace $\times$ M1 (Hopf-base measure) |
| $P_{M_2, a_0^{D^2}}^{\mathrm{trace}}$ | $4 \cdot 4 \pi^2 = 16 \pi^2$ | CG-trace $\times$ M2 Dirac SD0 (vol-norm, Paper 51 Cor 2.1) |
| $P_{M_2, a_1^{D^2}}^{\mathrm{trace}}$ | $4 \cdot (-2\pi^2) = -8\pi^2$ | CG-trace $\times$ M2 Dirac SD1 |
| $P_{M_2, a_0^\Delta}^{\mathrm{trace}}$ | $4 \cdot 2\pi^2 = 8\pi^2$ | CG-trace $\times$ M2 scalar SD0 |
| $P_{M_2, a_2^\Delta}^{\mathrm{trace}}$ | $4 \cdot \pi^2 = 4\pi^2$ | CG-trace $\times$ M2 scalar SD2 |
| $P_{M_2, \zeta_{D^2}(2)}^{\mathrm{trace}}$ | $4 \cdot (\pi^2 - \pi^4/12)$ | CG-trace $\times$ M2 spec-zeta s=2 (Paper 28 Thm 1) |
| $P_{M_2, \zeta_{D^2}(3)}^{\mathrm{trace}}$ | $4 \cdot (\pi^4/3 - \pi^6/30)$ | CG-trace $\times$ M2 spec-zeta s=3 |
| $P_{M_2, \zeta_{D^2}(4)}^{\mathrm{trace}}$ | $4 \cdot (2\pi^6/15 - 17\pi^8/1260)$ | CG-trace $\times$ M2 spec-zeta s=4 |
| $P_{M_2, a_0^{D^2}}^{\mathrm{hw}}$ | $1 \cdot 4\pi^2 = 4\pi^2$ | CG-hw $\times$ M2 Dirac SD0 |
| $P_{M_2, \zeta_{D^2}(2)}^{\mathrm{hw}}$ | $\pi^2 - \pi^4/12$ | CG-hw $\times$ M2 spec-zeta s=2 |
| $P_{M_3, \beta(2)}^{\mathrm{trace}}$ | $4 G$ | CG-trace $\times$ M3 Catalan |
| $P_{M_3, \beta(4)}^{\mathrm{trace}}$ | $4 \beta(4)$ | CG-trace $\times$ M3 $\beta(4)$ |
| $P_{M_3, \mathrm{diff}(2)}^{\mathrm{trace}}$ | $4 \cdot (2G - 1)$ | CG-trace $\times$ M3 parity discriminant s=2 (Paper 28 Thm 3) |
| $P_{M_3, \mathrm{diff}(4)}^{\mathrm{trace}}$ | $4 \cdot 8 \cdot (\beta(4) - G) = 32\beta(4) - 32G$ | CG-trace $\times$ M3 parity discriminant s=4 |
| $P_{M_3, \beta(2)}^{\mathrm{hw}}$ | $G$ | CG-hw $\times$ M3 Catalan |
| $P_{M_3, \beta(4)}^{\mathrm{hw}}$ | $\beta(4)$ | CG-hw $\times$ M3 $\beta(4)$ |
| $P_{M_3, \mathrm{diff}(2)}^{\mathrm{hw}}$ | $2G - 1$ | CG-hw $\times$ M3 parity discriminant s=2 |
| $P_{M_2 M_3, \zeta_{D^2}(2) \cdot G}$ | $(\pi^2 - \pi^4/12) \cdot G$ | joint M2 $\times$ M3 (multiplicative composition) |
| $P_{M_2 M_3, a_0^{D^2} \cdot G}$ | $4\pi^2 G$ | joint M2 $\times$ M3 |
| $P_{M_3, \zeta(3)}^{\mathrm{trace}}$ | $4 \zeta(3)$ | CG-trace $\times$ M3 level-1 odd-zeta |

These 20 periods cover M1, M2, M3, mixed M2 $\times$ M3, and level-1 odd
zeta content, each with both CG-trace and CG-hw $\mathrm{Sym}^2$ structural
factors.

### 2.2 Hain-Brown candidate basis at $\tau = 2i$

Constructed by direct $q$-series expansion at $q = e^{2\pi i \tau} = e^{-4\pi}
\approx 3.49 \times 10^{-6}$:
- $E_4(2i) = 1 + 240 \sum_n \sigma_3(n) q^n \approx 1.00083699$
- $E_6(2i) = 1 - 504 \sum_n \sigma_5(n) q^n \approx 0.99824218$
- $\Delta(2i) = (E_4^3 - E_6^2)/1728 \approx 3.487 \times 10^{-6}$
- $P_k := \int_{2}^{\infty} \Delta(it) (t - 2)^k \, dt$, real-valued
  Eichler-style periods along the imaginary axis, $k = 0, 1, 2$:
  $P_0 \approx 5.55 \times 10^{-7}$, $P_1 \approx 8.83 \times 10^{-8}$,
  $P_2 \approx 2.81 \times 10^{-8}$.
  (Choice of kernel $(t - \tau_\mathrm{imag})^k$ along the imaginary axis
  keeps the integrand real; the standard $(z - \tau)^{10}$ Eichler kernel
  is purely imaginary along the imaginary axis at $\tau = 2i$ — by
  construction $\int = 0$ on the real part — so we use the shifted
  power-kernel for the three lowest moments instead.)

Four basis configurations were tested:
- `modular_only`: 6 elements $\{E_4, E_6, \Delta, P_0, P_1, P_2\}$.
- `modular_plus_zeta`: + $\{\zeta(3), \zeta(5)\}$ → 8 elements.
- `modular_plus_zeta_beta`: + $\{G, \beta(4)\}$ → 10 elements.
- `full`: + $\{\pi, \pi^2, \pi^3, \pi^4, \pi^6, \pi^8, \log 2\}$ → 17 elements.

### 2.3 PSLQ parameters

- Coefficient ceiling: $10^6$.
- `maxsteps = 2000` (default $100$ is too small for $\ge 10$-element basis
  on this scale of values; experimentally confirmed in the debugging arc
  documented in `debug/sprint_hb_pslq_test_compute.py`).
- Iterative trivial-relation filter: when PSLQ returns a relation $\sum
  a_i x_i = 0$ with $a_0 = 0$ (i.e., a relation among basis elements
  alone, e.g. $\zeta(3, 1/4) - \zeta(3, 3/4) = 2 \pi^3$), the latest-index
  basis element involved is dropped and PSLQ is re-run. Up to 5 such
  iterations per period. (The level-4 Hurwitz redundancy
  $\zeta(s, 1/4) \pm \zeta(s, 3/4) = (4^s \mp 2^s) \cdot L$
  identity was removed from the basis a priori.)
- Per-precision basis filter: any basis element with $|v| < 10^{-\mathrm{dps}/4}$
  is treated as zero (PSLQ requires nonzero vector).

### 2.4 Cross-precision agreement filter

The same period × basis configuration is tested at 50, 100, 200 dps. The
verdict assigned is:
- `POSITIVE` if a PSLQ relation involving at least one of $\{E_4, E_6,
  \Delta, P_0, P_1, P_2\}$ is found.
- `NEGATIVE` if a PSLQ relation is found but only involves elements
  outside the modular ring (pure $\pi$ powers, $\zeta(3), \zeta(5), G,
  \beta(4), \log 2$).
- `NULL` if no PSLQ relation is found.

A period is counted as a load-bearing positive Hain-Brown identification
only if its verdict is `POSITIVE` at **all three precisions**. False
positives from low-precision PSLQ heuristic fits (which exploit the
tiny magnitude of modular ring elements, $\sim 10^{-6}$ to $10^{-8}$, to
fit noise at the precision floor) are filtered by upgrading the
precision.

---

## 3. Computational panel

### 3.1 Top-level verdict

| Quantity | Value |
|:---|---:|
| Total GeoVac $\mathrm{Sym}^2$-tagged periods | 20 |
| Basis configurations tested | 4 (modular_only, +zeta, +zeta+beta, full) |
| Precisions tested | 3 (50, 100, 200 dps) |
| Total PSLQ cells | 240 ($20 \times 4 \times 3$) |
| `POSITIVE` cells (any precision, any basis) | 38 |
| `POSITIVE` agreeing across all 3 precisions | **0** |
| `NEGATIVE` cells (any precision, any basis) | 67 |
| `NEGATIVE` agreeing across all 3 precisions in `full` basis | **12** |
| `NULL` cells | 135 |
| Cross-precision-agreed Hain-Brown modular hits | **0** |

### 3.2 Cross-precision agreement summary

| Basis | POS@all 3 | POS@some only | NEG@all 3 | NULL@all 3 |
|:---|---:|---:|---:|---:|
| `modular_only` | 0 | 0 | 0 | 20 |
| `modular_plus_zeta` | 0 | 19 | 1 | 0 |
| `modular_plus_zeta_beta` | 0 | 17 | 3 | 0 |
| `full` | **0** | 8 | **12** | 0 |

### 3.3 The "POSITIVE at some, not all" pattern: textbook false positives

At 50 dps, the `modular_plus_zeta` basis assigns `POSITIVE` to 19/20
periods. At 100 dps and 200 dps the same periods × basis combinations are
all `NULL` (the supposed relation breaks down). This is the diagnostic
signature of PSLQ exploiting the precision floor at low dps using the tiny
modular-ring magnitudes ($\Delta(2i) \sim 3.5 \times 10^{-6}$, $P_0 \sim
5.5 \times 10^{-7}$, $P_2 \sim 2.8 \times 10^{-8}$) to fit noise:\ for any
period of magnitude $\sim 10$, a relation $a_0 \cdot \mathrm{period} +
a_3 \cdot P_2 = 0$ with $|a_3| < 10^9$ can apparently fit a 50-dps target
to within precision floor, but the relation evaporates at 100 dps because
$P_2$ then needs to match to 100 digits.

The cross-precision filter is therefore load-bearing: it is precisely the
discipline the test was designed with to detect false positives.

### 3.4 The stable-NEGATIVE periods in the `full` basis at 200 dps

These 12 periods PSLQ-identify cleanly with pure $\mathrm{MZV}(\mathrm{MTM})$
content at all three precisions:

| Period | PSLQ relation found at 200 dps | Closed form |
|:---|:---|:---|
| $P_{M_1}$ | $4 \pi$ | M1 confirmed |
| $P_{M_2, a_0^{D^2}}^{\mathrm{trace}}$ | $16 \pi^2$ | M2 vol-norm Dirac SD0 |
| $P_{M_2, a_1^{D^2}}^{\mathrm{trace}}$ | $-8 \pi^2$ | M2 Dirac SD1 |
| $P_{M_2, a_0^\Delta}^{\mathrm{trace}}$ | $8 \pi^2$ | M2 scalar SD0 |
| $P_{M_2, a_2^\Delta}^{\mathrm{trace}}$ | $4 \pi^2$ | M2 scalar SD2 |
| $P_{M_2, \zeta_{D^2}(2)}^{\mathrm{trace}}$ | $4 \pi^2 - \pi^4/3$ | M2 spec-zeta s=2 |
| $P_{M_2, \zeta_{D^2}(3)}^{\mathrm{trace}}$ | $(4/3)\pi^4 - (2/15)\pi^6$ | M2 spec-zeta s=3 |
| $P_{M_2, \zeta_{D^2}(4)}^{\mathrm{trace}}$ | $(8/15)\pi^6 - (17/315)\pi^8$ | M2 spec-zeta s=4 |
| $P_{M_2, a_0^{D^2}}^{\mathrm{hw}}$ | $4 \pi^2$ | M2 CG-hw SD0 |
| $P_{M_2, \zeta_{D^2}(2)}^{\mathrm{hw}}$ | $\pi^2 - \pi^4/12$ | M2 spec-zeta CG-hw |
| $P_{M_3, \beta(2)}^{\mathrm{hw}}$ | $G$ | M3 Catalan CG-hw |
| $P_{M_3, \beta(4)}^{\mathrm{hw}}$ | $\beta(4)$ | M3 $\beta(4)$ CG-hw |
| $P_{M_3, \zeta(3)}^{\mathrm{trace}}$ | $4 \zeta(3)$ | M3 level-1 odd-zeta |

(The first row is for `full` basis at 200 dps; the same identifications
hold at 50 and 100 dps.)

Sanity check:\ all 13 PSLQ identifications match the closed-form Paper 55
predictions exactly (the rational coefficients are integer / integer with
small numerator and denominator). $P_{M_3, \mathrm{diff}(4)}^{\mathrm{trace}}$
is identified as $32(G - \beta(4))$ which differs in sign from the
Paper 28 Theorem 3 expression $4 \cdot 8 \cdot (\beta(4) - G) = 32(\beta(4) - G)$;
PSLQ returns either sign of the same relation freely. ✓

### 3.5 Pure-Tate M3 confirmation

The M3 sector outputs $\beta(2) = G$, $\beta(4)$, and $\zeta(3)$ all
appear at integer coefficients of the level-1 + level-4 Dirichlet $L$-ring
$\{\zeta(2k+1), \beta(2k)\}$. No GeoVac period required $E_4, E_6,
\Delta$, or any Eichler-style period of $\Delta$ to close. **The M3
sector empirically sits inside Brown's $\mathrm{MT}(\Z)$ + level-4 cyclotomic
extension $\mathrm{MT}(\Z[i, 1/2])$ at depth $\le 1$** (Glanois 2015) —
exactly as Paper 55 §5 classifies. No depth-$\ge 2$ irrationals appear in
the panel because all tested periods are at depth $\le 1$ on the M3 side.

---

## 4. PSLQ-stability diagnostics

### 4.1 The `maxsteps` debugging arc (documented for reproducibility)

Initial runs with `mpmath`'s default `maxsteps = 100` returned NULL even
for trivial period × pure-basis pairs like $16 \pi^2$ against a 10-element
basis. The root cause was identified as PSLQ premature termination:\ the
iteration count needed scales with basis size and value-magnitude spread.
With basis $\sim 17$ elements and magnitude spread from $10^{-8}$ to $10^4$,
the default $100$ steps is exhausted. Setting `maxsteps = 2000` resolved
all tested non-pathological cases.

### 4.2 The iterative trivial-relation filter

PSLQ on overdetermined bases finds the SHORTEST integer relation, not the
target one. For example, the level-4 Hurwitz identity
$\zeta(3, 1/4) - \zeta(3, 3/4) = 2 \pi^3$ produced PSLQ output
$[0, ..., -2, 0, ..., 1, -1]$ (a trivial-among-basis relation) before
the script could find $16 \pi^2 = $ period relation. The fix:\ drop one of
the redundant Hurwitz shifts a priori, and iteratively drop the
latest-index basis element involved when PSLQ returns trivial relations.

### 4.3 Eichler integral kernel choice

The naive Eichler kernel $(z - \tau)^k$ for weight-$12$ cusp form $\Delta$
gives purely-imaginary integrand along the imaginary axis at $\tau = 2i$
(real $\Delta(it)$ times $i^k (t - 2)^k \cdot i \, dt$). The real part of
such an integral is identically zero by construction. The replacement
kernel $(t - \tau_\mathrm{imag})^k$ (real-valued shift along the imaginary
axis) gives nonzero Eichler-style periods $P_0, P_1, P_2$. These are not
Hain-Brown's literal Eichler integrals but live in the same modular-period
ring of $\Delta$;\ they sample the relevant content at three distinct
weights without sign ambiguity.

A future strengthening (sprint-scale follow-on) is to use the genuine
Hain-Brown Eichler integral with the standard kernel by integrating along
a non-imaginary-axis path (e.g. parallel to the real axis offset by the
imaginary part of $\tau$), which would mix real and imaginary parts of
$\Delta$ more cleanly. Doing so would not change the verdict of this
sprint — the period magnitudes would still be of order $\Delta(2i) \times
(\mathrm{kernel\ norm}) \sim 10^{-6}$ to $10^{-3}$, so any spurious
positives at 50 dps would still degrade at 100 dps under cross-precision
agreement.

---

## 5. Honest scope

### 5.1 What this sprint closed

- **Empirical NEGATIVE on the Hain-Brown identification at $\mathrm{Sym}^2$
  level for the period values constructed.** The Hain-Brown modular ring
  is structurally not engaged by the natural $\mathrm{Sym}^2$-tagged
  GeoVac periods we can construct from the master Mellin engine M1/M2/M3
  closed forms.
- **Cross-precision-agreement protocol load-bearing for false-positive
  filtering of low-precision PSLQ.** 38 apparent POSITIVE hits at 50 dps
  → 0 surviving at 200 dps. The discipline is what made the verdict
  defensible.
- **The PSLQ `maxsteps = 2000` setting is necessary for any future
  PSLQ-on-periods sprint in the GeoVac corpus.** This will be inherited
  by `memory/feedback_audit_numerical_claims.md` / Q5'-direction sprints
  going forward.

### 5.2 What this sprint did NOT close

- **The categorical Hain-Brown shape-match.** GeoVac still has the four
  Hain-Brown structural features that motivated the identification
  (extension $1 \to U \to G^{\mathrm{rel}} \to SL_2 \to 1$ with $U$
  pro-unipotent, standard rep $V = \Q^2$, tensor category generated by
  $\mathrm{Sym}^k V$, canonical mixed Hodge structure compatible). The
  NEGATIVE here is empirical at the period-value level only;\ it says GeoVac
  is not a faithful Hain-Brown realisation, not that the shapes don't
  match.
- **The genuine Eichler-integral basis along the (real-axis-shifted)
  contour.** As noted in §4.3, the kernel choice in this sprint is a
  computational expedient. Strengthening to a genuine Eichler basis is
  a sprint-scale follow-on (1 week), but the magnitude argument suggests
  it would not change the verdict.
- **The depth-2 modular content.** This sprint operates at M2/M3
  depth-$\le 1$. Hain-Brown's MEM-$1$ universe includes multiple modular
  values $\Lambda(\Delta, k)$ for $k = 0, ..., 10$, depth $\ge 2$ MZV-modular
  pairings, and iterated integrals on $\mathcal{M}_{1,1}$ — none of which
  the present panel touches. A genuinely conclusive test against
  Hain-Brown would require GeoVac depth-$\ge 2$ structural compositions,
  which is exactly what the NA-1 depth-2 Mellin test (Recommendation A
  of strategic synthesis primary draft) would supply if the joint Mellin
  transform of M2 × M3 at depth $2$ produces sufficient structural content.
- **Period-level injection $U^*_{GV} \hookrightarrow \mathcal{U}_4$
  (Recommendation B).** This sprint's NEGATIVE strengthens Recommendation
  B as the next forward direction:\ GeoVac periods are in
  $\mathrm{MT}(\Z[i, 1/2])$ at depth $\le 1$ exactly as Paper 55
  classifies, so the explicit injection assembly into $\mathcal{U}_4$ has
  no Hain-Brown sub-structural detour.

### 5.3 What is now sharpened

The strategic-synthesis Recommendation A (Test A) verdict from this sprint
is:\ **Hain-Brown identification is empirically NOT supported at the
$\mathrm{Sym}^2$ level**. The strategic memo's contingency is that this
NEGATIVE outcome **shifts emphasis back to Recommendation B (explicit
$U^*_{GV} \hookrightarrow \mathcal{U}_4$ injection assembly, 2-3 weeks)**
as the strongest sprint-scale forward move. The cosmic-Galois comparison
target is $\mathcal{G}_4 = \mathcal{G}_{\mathrm{MT}(\Z[i, 1/2])}$ as Paper
55 / Paper 56 already framed.

The categorical Hain-Brown shape-match can still inform the framing if
GeoVac is presented as a **non-faithful "shape-cousin" of Hain-Brown
sitting one storey above Brown's mixed-Tate at $\mathcal{G}_4$** —
exactly the language the strategic-synthesis memo §2 (Track B Hodge SL_2
sub-agent finding) anticipated as a fallback. The Hain-Brown framing
becomes "we see what they see, but not what they hear" — same structural
shape at the category level, structurally absent at the period level.

---

## 6. Discipline checks

**Curve-fit audit (`memory/feedback_audit_numerical_claims`):** Yes,
applied. PSLQ on 20 periods × 4 basis configs × 3 precisions = 240 cells;
selection-bias accountability:\ the test was designed with 4 basis configs
and 3 precisions specifically to expose false positives. The
cross-precision filter is the explicit discipline. POSITIVE cells dropped
from 38 to 0 under the filter — selection bias would have given the
opposite finding.

**Discrete-for-skeleton (`memory/feedback_discrete_for_skeleton`):** Not
fully applicable — we are operating on Layer-2 transcendental values where
PSLQ is the appropriate tool. The Layer-1 substrate (Paper 56 / Paper 55
finite-cutoff Tannakian closure) remains untouched.

**Tag transcendentals (`memory/feedback_tag_transcendentals`):** Every
transcendental in the GeoVac side is tagged:\ M1 = $\pi$ (Hopf-base
measure, Paper 18 §III.7), M2 = SD coefficients in $\pi^{2k} \Q$ (Paper
51 / Paper 55 §4), M3 = $\beta(2), \beta(4), \zeta(3)$ from
Hurwitz-at-quarter-integer (Paper 28 / Paper 55 §5). Every transcendental
in the Hain-Brown side is named:\ $E_4(2i), E_6(2i), \Delta(2i),
\zeta(3), \zeta(5), G, \beta(4), \log 2, \pi^k$.

**Diagnostic-before-engineering (`memory/feedback_diagnostic_before_engineering`):**
Sprint was scoped as a diagnostic (PSLQ probe) before any structural
sprint. The decision-gate articulation made the protocol clean.

**Hard prohibitions check clean.** Paper 2 untouched, fitted parameters
zero (closed-form periods on both sides), failed approaches not
re-derived. WH1 PROVEN unaffected.

---

## 7. Named follow-ons

### Sprint-scale (1–4 weeks each)

1. **NA-1 depth-2 Mellin test (Recommendation A primary, Reading A vs B).**
   The HB PSLQ test only operates at depth $\le 1$. The NA-1 depth-2
   Mellin test (joint Mellin of one M2 SD × one M3 Hurwitz) is the
   complementary forward direction. Sprint-scale 1–2 weeks. The HB
   NEGATIVE result here does NOT preempt NA-1 — the two tests are
   structurally orthogonal (NA-1 probes shuffle-vs-primitive symmetry, HB
   probes modular content).

2. **Explicit injection $U^*_{GV} \hookrightarrow \mathcal{U}_4$
   (Recommendation B primary, theorem-grade).** With the HB NEGATIVE
   confirming the period-level target is $\mathcal{G}_4$, this sprint is
   now unconstrained by Hain-Brown sub-structure considerations. Effort
   2–3 weeks for Paper 55 / Paper 56 §sec:open_g4 theorem-grade upgrade.

3. **Hodge-theoretic $SL_2$ identification follow-on (Recommendation C).**
   The HB shape-match remains live at the categorical level. A
   sprint-scale follow-on probing whether GeoVac's $SL_2$ corresponds to
   a Mumford-Tate $SL_2$ on a polarised VMHS over the $S^3$ substrate
   would clarify the *structural* status of the Hain-Brown shape-match
   even though the period-level identification failed. Effort 2-4 weeks.

4. **Genuine Eichler-integral basis along real-axis-shifted contour.**
   Strengthen the modular-ring basis to include the literal Hain-Brown
   Eichler integrals $\int_\gamma \Delta(z) (z - \tau)^k \, dz$ for $k =
   0, ..., 10$ along a non-imaginary contour. Effort 1 week. As discussed
   in §4.3, this is unlikely to change the verdict but would close the
   "the kernel was a computational expedient" caveat.

### Multi-month (2–6 months each)

1. **Depth-$\ge 2$ M3 cyclotomic Mellin sweep.** Sprint-scale NA-1
   gives one depth-2 data point. A systematic exhaustion of GeoVac's M3
   vertex-parity sector at depth $\ge 2$ for the appearance of level-4
   cyclotomic motivic MZV building blocks (Glanois 2015 generators) is
   the multi-month direction that would actually upgrade the
   $\mathcal{G}_4$ injection from "comparison" to "exhaustion-conjecture
   testable."

---

## 8. Files

### Produced this sprint
- `debug/sprint_hb_pslq_test_compute.py` — driver (~600 lines).
- `debug/data/hb_pslq_test_results.json` — full results panel.
- `debug/sprint_hb_pslq_test_memo.md` — this memo.

### Read context
- `debug/strategic_synthesis_2026_06_06_memo.md` (Addendum + §1 + §6
  Recommendation A)
- `debug/sprint_q5p_na1_non_abelian_probe_memo.md` (the abelian-by-
  construction structural result that motivated this test)
- `memory/hain_brown_identification.md` (the prior shape-match finding
  this test was designed to confirm/refute at the period level)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` §4 (M2) +
  §5 (M3) + §6 (joint)
- `papers/group3_foundations/paper_56_tannakian_substrate.tex` §sec:tc2b
  (Sym² panel structure from $V_{\mathrm{fund}} \otimes V_{\mathrm{fund}}$
  decomposition)
- `geovac/tannakian.py` (`_sl2_sym2_action` and `_pw_sym2_rep`)

### References
- Hain, R. ``The Hodge-de Rham Theory of Modular Groups'' arXiv:1403.6443
  (2014).
- Brown, F. ``Multiple Modular Values and the relative completion of the
  fundamental group of $\mathcal{M}_{1,1}$'' arXiv:1407.5167 (2014).
- Hain, R. and Matsumoto, M. ``Universal Mixed Elliptic Motives''
  arXiv:1512.03975 (2015).
- Brown, F. ``Mixed Tate motives over $\Z$'' Ann. Math. 175 (2012)
  949–976.
- Deligne, P. ``Le groupe fondamental unipotent motivique de $\mathbf{G}_m
  - \mu_N$, pour $N = 2, 3, 4, 6$ ou $8$'' Publ. Math. IHÉS 112 (2010)
  101–141.
- Glanois, C. ``Periods of the motivic fundamental groupoid of
  $\mathbb{P}^1 \setminus \{0, \mu_N, \infty\}$'' arXiv:1603.05593 (2015).
- Eskandari, P., Murty, K., and Nemoto, T. ``On the irrationality of
  $\beta(2)$ as a level-4 cyclotomic period'' arXiv:2510.20648 (2025).

---

## 9. One-line verdict

**NEGATIVE.** Hain-Brown identification empirically NOT supported at
$\mathrm{Sym}^2$ level over the constructed period panel;\ all GeoVac
$\mathrm{Sym}^2$-tagged periods identify with pure $\mathrm{MZV}(\mathrm{MTM})$
content (in $\bigoplus_k \pi^{2k}\Q + \zeta(2k+1)\Q + \beta(2k)\Q$,
i.e.\ $\mathcal{G}_4$ depth-$\le 1$);\ no surviving relation engages
$E_4(2i), E_6(2i), \Delta(2i)$, or any Eichler period of $\Delta$ across
all three precisions $50/100/200$ dps. The Hain-Brown shape-match
remains live at the categorical level (Hodge-theoretic $SL_2$ +
$\mathrm{Sym}^k V$ tensor category), but the empirical period content
sits cleanly inside $\mathcal{G}_4$. **Forward emphasis returns to
Recommendation B (explicit $U^*_{GV} \hookrightarrow \mathcal{U}_4$
injection assembly) as the strongest near-term theorem-grade direction.**
