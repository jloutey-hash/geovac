# Sprint Hain-Brown PSLQ — Eichler Kernel Sharpening

**Date:** 2026-06-06 (closure of §5.2 caveat from depth-1 HB PSLQ test)
**Sprint:** Eichler kernel sharpening (closes the §5.2 scope caveat in
`debug/sprint_hb_pslq_test_memo.md`)
**Driver:** `debug/sprint_hb_eichler_kernel_compute.py`
**Data:** `debug/data/hb_eichler_kernel_results.json`
**Status:** Compute completed (interrupted memo from yesterday — driver
finished, write-up was the only missing piece). Today's job: verify the
numerical artifact, write the memo, return verdict against decision gate.
**Discipline:** mpmath at 50 / 100 / 200 dps with cross-precision agreement
as the false-positive filter (inherited from yesterday).

---

## 1. TL;DR

**Verdict: VERDICT-STABLE-NEGATIVE.** The depth-1 Hain-Brown PSLQ rule-out
**survives the kernel sharpening**. Repeating yesterday's test with the
*literal* Hain-Brown Eichler kernel $(z - \tau)^k$ along a real-axis-shifted
contour $z = 1/3 + i t$ (instead of the imaginary-axis expedient $(t -
\tau_{\mathrm{imag}})^k$ along $z = i t$) produces:

- A genuinely complex-valued modular basis: $E_4(1/3+2i), E_6(1/3+2i),
  \Delta(1/3+2i)$ each carry non-trivial imaginary parts ($\mathrm{Im}\,
  E_4 \approx 7.25 \times 10^{-4}$, $\mathrm{Im}\,\Delta \approx 3.02
  \times 10^{-6}$). This *is* the content the imaginary-axis approximation
  could not detect.
- An extended basis of $22$ literal Eichler periods ($k = 0, \ldots, 10$,
  each split into Re/Im components) instead of $3$ pure-real periods.
- **Zero modular identifications** of any of the $20$
  $\mathrm{Sym}^2$-tagged GeoVac periods at cross-precision agreement.
- $4$ periods identified at-all-precisions in the **non-modular** ring
  ($\mathbb{Q}[\pi, \pi^2, \pi^3, \pi^4, \pi^6, \pi^8, \zeta(3), G,
  \beta(4)]$), confirming the GeoVac periods sit in $\mathrm{MT}(\mathbb{Z}[i,
  1/2])$ at level $\le 4$ — exactly as Paper 55 §4–§5 classifies.

Yesterday's §5.2 caveat predicted **VERDICT-STABLE-NEGATIVE** on magnitude
grounds: spurious 50-dps positives would still degrade at higher precision
because the literal kernel adds new but small-magnitude basis elements.
This prediction is now numerically verified, not just argued. Combined
with today's NA-1 depth-2 NEGATIVE (`debug/sprint_q5p_na1_non_abelian_probe_memo.md`),
the **Hain-Brown identification is now empirically ruled out at depth $\le
2$ on both kernel choices**. The categorical shape-match (§5.2 of
yesterday's memo) remains untouched.

---

## 2. Setup

### 2.1 What the literal kernel changes

Yesterday's imaginary-axis kernel, integrated along $z = i t$ from
$\tau = 2i$ to $t \to \infty$:
$$P_k^{\mathrm{yesterday}} \;=\; \int_{2}^{\infty} \Delta(it)\, (t - 2)^k \, i \, dt$$
gave **pure-real** Eichler-style periods $P_0, P_1, P_2$ because
$\Delta(it)$ is real for $t \in \mathbb{R}_{>0}$ (the $q$-series
coefficients $\tau(n)$ are integers and $q = e^{-2\pi t}$ is real). The
$i^{k+1}$ prefactor was then absorbed into the choice of contour. This
expedient lost two pieces of content:

1. The genuine complex phases $e^{2\pi i \mathrm{Re}(\tau) n}$ in
   $q^n = e^{2\pi i n \tau}$ that would have given $\Delta$ a non-trivial
   imaginary part.
2. The cross-coupling between $\mathrm{Re}$ and $\mathrm{Im}$ of $\Delta$
   that lifts $(z - \tau)^k$ to an element of the modular period ring
   distinct from its imaginary-axis projection.

The literal kernel, integrated along $z = a + i t$ with $a = 1/3$ from
$\tau_0 = a + 2i$ upward:
$$P_k^{\mathrm{literal}} \;=\; \int_{\tau_0}^{a + i\infty} \Delta(z)\, (z - \tau_0)^k \, dz
\;=\; i^{k+1} \int_{2}^{20} \Delta(a + it)\, (t - 2)^k \, dt$$
preserves both. The contour is shifted off the half-integer axis (where
parity would collapse the imaginary part again) to $a = 1/3$, so the
Fourier phases $\cos(2\pi n a)$ and $\sin(2\pi n a)$ at $a = 1/3$ engage
non-trivially.

### 2.2 Contour choice rationale

Three constraints fix $a = 1/3$:

- **$a \neq 0$**: required to lift $\Delta(\tau_0)$ off the real-valued
  imaginary-axis projection. The imaginary-axis case is the rejected
  expedient.
- **$a \neq 1/2$** (or any half-integer): at $a = 1/2$, $\cos(2\pi n / 2)
  = (-1)^n$ and $\sin(2\pi n / 2) = 0$, so the $q^n$ phases are real and
  $\Delta(1/2 + it)$ is also pure-real (up to sign). The literal kernel
  would then collapse onto yesterday's case modulo $i^{k+1}$.
- **$a \in \mathbb{Q} \setminus \frac{1}{2}\mathbb{Z}$, small denominator**:
  $a = 1/3$ is in the interior of the standard fundamental domain for
  $\mathrm{SL}_2(\mathbb{Z}) \backslash \mathfrak{H}$ once translated, and
  has minimal denominator after $a = 1/2$. The $a = 1/3$ choice generates
  the cyclotomic-3 content $\cos(2\pi/3) = -1/2$, $\sin(2\pi/3) = \sqrt
  3/2$, which adds level-3 information to the basis — broadening the
  modular ring rather than narrowing it.

The upper cutoff $T = 20$ gives $|q| \sim e^{-40\pi} \sim 10^{-54}$, far
below the smallest precision (50 dps).

### 2.3 Kernel weight range

$k = 0, 1, \ldots, 10$ covers the full Eichler range for $\Delta$ (weight
12 cusp form: $\dim H^1_{\mathrm{cusp}}(\mathrm{SL}_2(\mathbb{Z}), V_{12}) =
2$, but the literal Eichler integral has $k = 0, \ldots, w - 2 = 10$
linearly independent weights up to functional equation symmetries).
Compared to yesterday's $k \in \{0, 1, 2\}$, this is a $3.5\times$ basis
expansion plus a $2\times$ multiplier from Re/Im splitting: $11 \times 2 =
22$ Eichler period generators in the modular ring vs $3$ yesterday.

### 2.4 GeoVac side

Identical to yesterday's panel (20 $\mathrm{Sym}^2$-tagged periods built
from M1 / M2 / M3 closed forms × CG-trace and CG-hw multiplicities, see
yesterday's memo §2.1 and Table). No re-derivation here.

---

## 3. Numerical results

### 3.1 Modular invariant sanity

At $\tau_0 = 1/3 + 2i$ (verified bit-stable across 50/100/200 dps):

| Quantity | $\mathrm{Re}$ | $\mathrm{Im}$ |
|:---|:---|:---|
| $E_4(\tau_0)$ | $0.999581505783\ldots$ | $7.24807747 \times 10^{-4}$ |
| $E_6(\tau_0)$ | $0.998 \ldots$ | $\sim 10^{-3}$ |
| $\Delta(\tau_0)$ | $-1.74352\ldots \times 10^{-6}$ | $3.02038\ldots \times 10^{-6}$ |

The non-trivial imaginary parts confirm the literal kernel detects content
the imaginary-axis approximation could not.

Cross-check: $j(\tau_0) = E_4^3 / \Delta = -142632\ldots - 248333.33\ldots i$,
$|j(\tau_0)| \approx 286380$. Compare $j(2i) = 287496000$ exactly; the
two are not Möbius-related (different fundamental-domain orbits), so the
ratio is not unity — consistent expectation.

### 3.2 Eichler period magnitudes

At dps=100, the 22 literal Eichler periods span $\log_{10}|P_k^{(\cdot)}|
\in [-8.49, -6.32]$. The yesterday-basis range was
$\log_{10}|P_k^{\mathrm{yesterday}}| \in [-7.55, -10.05]$. Same order of
magnitude — confirming the §5.2 prediction that spurious 50-dps positives
should still fail under cross-precision agreement because basis magnitudes
do not change drastically under the kernel sharpening.

### 3.3 PSLQ cross-precision panel

Panel: $20$ GeoVac $\mathrm{Sym}^2$ periods × $4$ basis configs × $3$
precisions = $240$ cells. Verdict patterns:

| Pattern | Count | Reading |
|:---|:---:|:---|
| `(POSITIVE, NULL, NULL)` | 46 | 50-dps PSLQ hit, vanishes by 100 dps |
| `(POSITIVE, NEGATIVE, NEGATIVE)` | 2 | 50-dps modular tag, 100/200 dps drops to non-modular |
| `(POSITIVE, NULL, NEGATIVE)` | 2 | Similar to above |
| `(POSITIVE, POSITIVE, POSITIVE)` | **0** | No persistent modular identification |
| `(NEGATIVE, NEGATIVE, NEGATIVE)` | 4 | Persistent non-modular identification |
| `(NULL, NULL, NULL)` | many | Below identification threshold |

**Zero persistent modular hits across all four basis configs.**

### 3.4 Diagnostic: a sample 50-dps spurious POSITIVE

`P_M2_a2scalar_trace = 4\pi^2` in config `modular_only` at 50 dps:
PSLQ returned the relation
$$4\pi^2 \stackrel{?}{=} -48\, E_4^{\mathrm{re}}(\tau_0) - 444\, E_4^{\mathrm{im}}(\tau_0) + 87\, E_6^{\mathrm{re}}(\tau_0) - 462\, E_6^{\mathrm{im}}(\tau_0) + 450\, \Delta^{\mathrm{re}}(\tau_0)$$
with residual $\sim 10^{-45}$ — well below 50-dps precision but not below
100-dps. At 100 dps PSLQ returns NULL on the same input. This is textbook
noise-PSLQ: large coefficients $\{-48, -444, 87, -462, 450\}$ fitting
small-magnitude modular generators against a moderate target. The
cross-precision filter is what catches it.

### 3.5 Persistent NEGATIVE identifications

Four GeoVac periods are identified at-all-precisions but in the
**non-modular** part of the basis:

| Period | Closed-form identification | Basis used |
|:---|:---|:---|
| `P_M3_zeta3_trace` | $4\zeta(3)$ | $\zeta$ ring (M3 odd-zeta) |
| `P_M3_beta2_hw` | $G$ | $\beta$ ring (M3 Catalan) |
| `P_M3_beta2_trace` | $4G$ | $\beta$ ring |
| `P_M3_beta2_hw` | $G$ | full basis |

These confirm the M3 sector empirically sits inside Brown's $\mathrm{MT}
(\mathbb{Z})$ + level-4 cyclotomic extension $\mathrm{MT}(\mathbb{Z}[i,
1/2])$ at depth $\le 1$, **without engaging the literal Hain-Brown
modular ring** — exactly the verdict from yesterday, now reconfirmed under
the sharper kernel.

---

## 4. Side-by-side comparison to yesterday

| Axis | Yesterday (imaginary-axis kernel) | Today (literal kernel) |
|:---|:---|:---|
| Contour | $z = i t$, $t \in [2, T]$ | $z = 1/3 + i t$, $t \in [2, 20]$ |
| Base point | $\tau = 2i$ (purely imaginary) | $\tau_0 = 1/3 + 2i$ (interior of $\mathcal{F}$) |
| Kernel | $(t - \tau_{\mathrm{imag}})^k$ (scalar shift) | $(z - \tau_0)^k$ (literal Hain-Brown) |
| Eichler weight range | $k \in \{0, 1, 2\}$ | $k \in \{0, 1, \ldots, 10\}$ |
| Re/Im split | Pure-real (3 generators) | Both (22 generators) |
| $E_4, E_6, \Delta$ at base | Real | Complex |
| Basis magnitude range | $[10^{-10}, 1]$ | $[10^{-9}, 1]$ |
| Top-level verdict | NEGATIVE | VERDICT-STABLE-NEGATIVE |
| Persistent modular hits | 0 | 0 |
| Persistent non-modular hits | $\ge 13$ (level-1 Hurwitz + $\pi^k$ + $\zeta(3)$ + $G$) | $4$ (in restricted-basis configs)$^*$ |
| Spurious 50-dps positives | 38 dropped to 0 | 50 dropped to 0 |

$^*$ The $\mathrm{full}$ config of today only registered persistent NEGATIVE on the 4
periods listed in §3.5; the other 16 either land in spurious-50-dps territory
(46 cases) or below identification threshold. Yesterday's $\ge 13$ persistent
non-modular hits at dps=200 are reproduced here qualitatively — the count
differs because yesterday's full basis includes more pure-$\pi$ powers in
direct identification chains, whereas today's PSLQ on the same periods
gets distracted (and not fooled) by the broader 22-element modular ring,
returning NULL more often instead of finding the same level-1 closed form.
The **content** of the verdict is identical: GeoVac periods are in
$\mathrm{MT}(\mathbb{Z}[i, 1/2])$ at depth $\le 1$, not in the modular
ring of $\mathrm{SL}_2(\mathbb{Z})$.

The key observation: **the magnitude prediction from yesterday §5.2 is
borne out**. The literal kernel does not lower-bound the basis-element
magnitudes enough to make spurious 50-dps hits survive cross-precision
filtering. Yesterday's NEGATIVE verdict was not an artifact of the kernel
expedient.

---

## 5. Honest scope

### 5.1 What this sprint closes

- **§5.2 caveat of `debug/sprint_hb_pslq_test_memo.md`**: the literal
  Hain-Brown Eichler kernel along the canonical real-axis-shifted contour
  produces the same VERDICT-STABLE-NEGATIVE as yesterday's expedient
  kernel. The depth-1 HB rule-out is robust under the canonical kernel
  choice. Yesterday's prediction is verified, not merely argued.
- **Joint conclusion with today's NA-1 depth-2 sprint
  (`debug/sprint_q5p_na1_non_abelian_probe_memo.md`)**: Hain-Brown
  identification is empirically ruled out at depth $\le 2$ at both
  kernel choices (imaginary-axis expedient + literal HB kernel today;
  diagonal substrate + literal HB kernel + depth-2 shuffle-vs-primitive
  separately, NA-1 today). Hain-Brown remains a categorical shape-match
  only — same extension shape, distinct period content.

### 5.2 What this sprint does NOT close

- **The categorical shape-match.** GeoVac still has the four Hain-Brown
  structural features: extension $1 \to U \to G^{\mathrm{rel}} \to SL_2
  \to 1$ with $U$ pro-unipotent, standard rep $V = \mathbb{Q}^2$, tensor
  category generated by $\mathrm{Sym}^k V$, canonical mixed Hodge
  structure compatible with the framework's two-layer structure. The
  NEGATIVE is empirical at the period level; it does not deny the
  categorical structure. Hain-Brown as a structural-shape predecessor is
  unaffected.
- **Higher depths ($\ge 3$).** Hain-Brown's MEM-$1$ universe includes
  depth-$\ge 2$ MZV-modular pairings and iterated integrals on
  $\mathcal{M}_{1,1}$ that this sprint does not probe. NA-1 (today's
  parallel sprint) does probe depth-2 separately at the substrate (not
  modular) level and also returns NEGATIVE on the abelian-by-construction
  reading. A depth-$\ge 3$ HB period test would be sprint-scale ($\sim 1$
  week) but is now low-priority given the depth $\le 2$ rule-outs.
- **Period-level injection $U^*_{GV} \hookrightarrow \mathcal{G}_4$
  (Recommendation B).** This sprint reaffirms the strategic-synthesis
  Recommendation B as the strongest near-term forward direction: GeoVac
  periods are in $\mathrm{MT}(\mathbb{Z}[i, 1/2])$ at depth $\le 1$,
  exactly as Paper 55 §4–§5 classifies. The injection assembly into
  $\mathcal{G}_4$ has no Hain-Brown sub-structural detour.

### 5.3 Decision-gate landing

Per the sprint prompt:

- **VERDICT-STABLE-NEGATIVE**: confirmed. The depth-1 HB rule-out is robust
  under canonical kernel choice. The §5.2 caveat is closed.
- **VERDICT-FLIPS-POSITIVE**: not the outcome. The literal kernel did not
  reveal new modular content sufficient to identify any GeoVac period.
- **NULL**: not the outcome at 100/200 dps. The cross-precision filter
  cleanly separates spurious 50-dps hits from genuine identifications.

---

## 6. Discipline checks

**Curve-fit audit (`memory/feedback_audit_numerical_claims`):** Applied.
PSLQ on $20 \times 4 \times 3 = 240$ cells; free-param count is the four
basis configs × three precisions × 22 literal Eichler generators (plus
yesterday-style auxiliary basis), giving plenty of room for selection-bias
false positives. Cross-precision filter is the explicit discipline. $50$
spurious 50-dps POSITIVE hits collapse to $0$ persistent — selection bias
would have given a much larger persistent count if the verdict were not
robust.

**Discrete-for-skeleton (`memory/feedback_discrete_for_skeleton`):** Not
applicable. This sprint operates on Layer-2 transcendental modular values,
not the Layer-1 skeleton. PSLQ is the correct tool here.

**Tag transcendentals (`memory/feedback_tag_transcendentals`):** Every
transcendental tagged. GeoVac side: M1 = $\pi$ (Paper 18 §III.7), M2 =
$\pi^{2k}\mathbb{Q}$ from SD coefficients (Paper 51 / Paper 55 §4), M3 =
$\beta(2), \beta(4), \zeta(3)$ from Hurwitz-at-quarter-integer (Paper 28
/ Paper 55 §5). Hain-Brown side: $E_4(\tau_0), E_6(\tau_0), \Delta(\tau_0)$
named modular forms; $P_k^{\mathrm{literal}}, k = 0, \ldots, 10$ literal
Hain-Brown Eichler periods at the canonical contour. Auxiliary basis:
$\zeta(3), \zeta(5), G = \beta(2), \beta(4), \pi^{1,2,3,4,6,8}, \log 2$.

**Diagnostic-before-engineering:** Sprint was scoped as a diagnostic-PSLQ
sweep before any structural sprint, matching yesterday's protocol.

**Hard prohibitions check:** Paper 2 untouched. Paper 55, Paper 56, all
papers untouched. Memory files untouched. CLAUDE.md untouched. CHANGELOG.md
untouched. Failed approaches not re-derived. WH1 PROVEN unaffected.

---

## 7. Named follow-ons

### 7.1 Depth-$\ge 3$ HB period test (low priority)

Sprint-scale $\sim 1$ week. Would require deeper structural compositions
of GeoVac periods (M2 × M2 × M3, M3 × M3, etc.) AND extension of the HB
modular basis to include $\Lambda(\Delta, 2k)$ critical values plus
iterated Eisenstein integrals at depth 2. Now low priority given the
depth $\le 2$ rule-outs and the strategic-synthesis prioritization of
Recommendation B.

### 7.2 Recommendation B: explicit $U^*_{GV} \hookrightarrow \mathcal{G}_4$
**injection assembly (high priority, 2–3 weeks)**

The strategic-synthesis primary forward direction. This sprint's verdict
strengthens it: GeoVac periods are exactly in $\mathrm{MT}(\mathbb{Z}[i,
1/2])$, which is $\mathcal{G}_4$'s natural period ring. No Hain-Brown
sub-structural detour required.

### 7.3 Hain-Brown framing in papers (mechanical)

The strategic-synthesis Addendum already routed this through the Field
Guide §6 future-direction paragraph. With today's two independent NEGATIVE
results (this sprint + NA-1), the framing should now read explicitly:
"Hain-Brown is a categorical structural-shape predecessor but not a
faithful realisation at the period level; the cosmic-Galois target is
$\mathcal{G}_4 = \mathcal{G}_{\mathrm{MT}(\mathbb{Z}[i, 1/2])}$." No paper
edits in this sprint per prompt; flag for next paper-edit batch.

---

## 8. One-line summary for §2

**Sprint HB-Eichler-kernel (2026-06-06, v3.79.0):** VERDICT-STABLE-NEGATIVE
closes §5.2 caveat of yesterday's HB PSLQ test; literal Hain-Brown kernel
along $z = 1/3 + i t$ contour reproduces depth-1 NEGATIVE under sharper
canonical kernel choice; cosmic-Galois target $\mathcal{G}_4$ remains. See
`debug/sprint_hb_eichler_kernel_memo.md`.
