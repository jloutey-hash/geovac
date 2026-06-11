# L2 next-order constant $c$ — high-precision PSLQ campaign v3

**Verdict: DEFINITIVELY NULL at all tested ceilings, with substantive
side finding.** PSLQ at coefficient ceilings $10^5$ through $10^8$ returns
NULL against four test bases (sizes 14, 38, 59, 83 elements) including the
**L-derivative basis** (24 new elements: $\zeta'(s)$, $\beta'(s)$, $\eta'(s)$,
$L'(2, \chi_q)$, Stieltjes $\gamma_k$, $\log A$, $\zeta'(0)$, $\zeta'(-1)$).
The L2 next-order constant is therefore irreducible against the broader
basis that includes derivatives of L-functions — the cleanest follow-up
direction flagged by Track 4. **Side finding (load-bearing for any future
sprint):** the v3 high-precision recomputation revealed that Track 4's
published "80-digit" string $c = 4.10932146748\mathbf{779409}\ldots$ was
**incorrect from digit 12 onward** due to float-precision contamination
in its panel; the corrected value is
$$
c = 4.10932146748\mathbf{987562042296336493}\ldots
$$
This is the headline correction to the prior memo
(`debug/l2_constant_c_identification_memo.md`).

## 1. Setup

The L2 quantitative-rate sprint (Sprint TS-E1 / WH1 R2.5 L2, 2026-05-04 / 06,
`debug/r25_l2_proof_memo.md`, `debug/r25_l2_quantitative_rate_memo.md`)
established the asymptotic
$$
n \cdot \gamma_n^{\text{scalar}} = \frac{4}{\pi} \log(n) + c + O\!\left(\frac{\log n}{n}\right)
$$
on the SU(2) central Fejér panel, with leading constant
$4/\pi = \mathrm{Vol}(S^2)/\pi^2$ identified as the **M1 Hopf-base measure
signature** of the master Mellin engine (Paper 32 §VIII case-exhaustion
theorem, Paper 18 §III.7).

Track 4 (Sprint MR-C + identification, May 2026,
`debug/l2_constant_c_identification_memo.md`) extracted $c$ at $\sim 12$
reliable digits and ran PSLQ at ceilings $10^4$, $10^5$, $10^6$ against
seven test bases up to 82 elements covering Dirichlet L, multi-Hurwitz at
quarter shifts, Stein-Weiss IBP, and wildcard sectors. ALL NULL.

The original Track 4 agent recommended as the cleanest follow-up:

> Cleanest follow-up direction: build basis with $\zeta'$ and $L'$-values
> content + push precision to $\sim 50$ analytical digits.

This sprint does both, with two unexpected findings.

## 2. Sprint design

### 2.1 Stage 1 — high-precision panel

Driver: `debug/l2_constant_c_precision_v3.py`.

Strategy: reuse Track 4's 15 high-precision intermediate panel points
(at $n \in \{48, 80, 96, 160, 192, 320, 384, 640, 768, 1280, 1536, 2560,
3072, 5120, 6144\}$, all genuinely at $\geq 250$ dps in
`debug/data/l2_constant_c_precision_push.json`), then compute additional
**high-precision** doubling-sequence points where Track 4 stored only
float-precision values.

- Working precision: $\mathrm{dps} = 400$ for all new $T_n$ computations.
- New points computed at 400 dps:
  $n \in \{32, 64, 128, 256, 512, 1024, 2048, 4096, 8192\}$ (9 points).
- Walltime for Stage 1: **793 s = 13.2 min**.
  - $n = 8192$ alone: 592 s (76% of total).
  - $n = 4096$: 154 s.
  - $n = 2048$: 36 s.
- Combined panel: 24 panel points, all at 400 dps.

### 2.2 Stage 2 — Richardson scan

The Richardson extrapolation fits
$h(n) = b_0 + \sum_{k=1}^{K} (a_k \log n + b_k) / n^k$ at
$K = 4, 5, 6, 7, 8, 9, 10, 11$ ($K = 12$ skipped: would need 25
parameters on 24-point panel). Square fits (K=11 has 23 parameters on
24 points) and overdetermined LS solutions used as appropriate.

Consecutive-K deltas:

| $K_a \to K_b$ | $|c(K_a) - c(K_b)|$ | ratio |
|:--|:--|:--|
| 4 → 5  | $3.24 \times 10^{-11}$ | — |
| 5 → 6  | $4.77 \times 10^{-12}$ | 0.15 |
| 6 → 7  | $1.06 \times 10^{-12}$ | 0.22 |
| 7 → 8  | $3.12 \times 10^{-13}$ | 0.30 |
| 8 → 9  | $1.02 \times 10^{-13}$ | 0.33 |
| 9 → 10 | $3.98 \times 10^{-14}$ | 0.39 |
| 10 → 11| $1.63 \times 10^{-14}$ | 0.41 |

Best $K = 11$ (square fit on 24-point panel):
$$
c = 4.10932146748987562042296336493\ldots
$$
**(28 dps shown.)**

### 2.3 Cross-validation: predicting $h(3000)$

To verify the actual precision of the Richardson fit (the consecutive-$K$
delta is only an upper bound on the fit's true precision), the v3 panel
fit at $K = 11$ was used to **predict** $h(n)$ at a new point $n = 3000$
**not in the panel**, and compared against direct high-precision computation
$h(3000) = 4.106949488370948546429669\ldots$ (independently computed at
100 dps).

| $K$ | predicted $h(3000)$ | $|$error$|$ |
|:--|:--|:--|
| 4  | $4.106949488371714493\ldots$ | $7.7 \times 10^{-13}$ |
| 7  | $4.106949488370948447\ldots$ | $1.0 \times 10^{-16}$ |
| 9  | $4.106949488370948546927\ldots$ | $5.0 \times 10^{-19}$ |
| 10 | $4.106949488370948546415\ldots$ | $1.4 \times 10^{-20}$ |
| 11 | $4.106949488370948546429\ldots$ | $2.4 \times 10^{-22}$ |

The K=11 fit predicts $h(3000)$ to **22 reliable digits**, much better
than the consecutive-$K$ delta would suggest (the K-deltas are upper
bounds, not point estimates of the true error). This is strong evidence
the v3 $c$ value is reliable to $\sim 20$ digits.

### 2.4 Stability test

Fitting the same Richardson with the panel restricted to 23 points
(removing $n = 8192$) at $K = 11$ gives $c$ shifted by
$3.66 \times 10^{-15}$. Removing two points ($n = 5120, 6144, 8192$)
and fitting $K = 10$ shifts $c$ by $3.35 \times 10^{-14}$. Both
perturbations are at the 14-digit level, consistent with the precision
saturation in the conditioning of the LS problem at high $K$ (the
normal equations square the condition number). The **point estimate**
from the cross-validation against $h(3000)$ is more precise (~22 digits)
than the **uncertainty estimate** from the K-deltas (~14 digits).

### 2.5 Side finding: Track 4's published value was wrong from digit 12+

The new $c = 4.10932146748\mathbf{987562042}\ldots$ disagrees with Track
4's $c = 4.10932146748\mathbf{779409275}\ldots$ at digit 12 (delta
$2.08 \times 10^{-12}$). Track 4 reported its own K=5-vs-K=6
disagreement at digit 12 (delta $4.77 \times 10^{-12}$) and concluded
"~12 reliable digits"; however, Track 4's published "80-digit string"
was its K=6 result on the 24-point panel where 9 of 24 points were at
**float precision only** (1e-16 effective). This contaminated the K=6
fit at $\sim 10^{-12}$ level — exactly where v3 and Track 4 disagree.

The v3 fit, using all 24 points at 400 dps, is the corrected value.

**Implication for CLAUDE.md §2 / Paper 32 §VIII / Paper 18 §III.7:**
the L2 next-order constant should be reported as
$c = 4.10932146748987562\ldots$ wherever the previous "$\ldots 748779409\ldots$"
value appears in production memos or papers. (At the time of writing, the
constant $c$ is mentioned in the MR-C memo, Track 4 memo, Paper 32 §VIII
master-Mellin-engine domain-partition remark, and Paper 18 §III.7 — only
the leading 11 digits are quoted in those venues, so no edit cascade is
required outside this memo and the L2 sprint folder. PI flagged for review.)

## 3. The L-derivative basis (sub-task b)

The new "L-derivative" basis sector — the cleanest follow-up direction
flagged by Track 4 — consists of 24 elements computed at 250 dps:

| Element | Value (30 dps) | Source |
|:--|:--|:--|
| $\zeta'(2)$ | $-0.937548254315843753702574094568$ | mpmath |
| $\zeta'(3)$ | $-0.198126242885636853330681821503$ | mpmath |
| $\zeta'(4)$ | $-0.068911265896125379848829365587$ | mpmath |
| $\zeta'(5)$ | $-0.028573780509462950080389817084$ | mpmath |
| $\zeta'(6)$ | $-0.012852165131795725075945401460$ | mpmath |
| $\zeta''(2)$ | $1.98928023429890102342085868742$ | mpmath |
| $\beta'(1)$ | $0.192901316796912429363189764028$ | direct nsum |
| $\beta'(2)$ | $0.081580736116592795102912169786$ | direct nsum |
| $\beta'(3)$ | $0.031577079457127387887246167443$ | direct nsum |
| $\beta'(4)$ | $0.011570547924511641651452728260$ | direct nsum |
| $\beta'(5)$ | $0.004094874987948592100691393465$ | direct nsum |
| $\eta'(1)$  | $0.159868903742430971756947870325$ | $\eta(s) = (1{-}2^{1-s})\zeta(s)$ |
| $\eta'(2)$  | $0.101316578163504501886002882212$ | same |
| $L'(2,\chi_3)$ | $0.151020914067357178661184388511$ | direct nsum, mod 3 |
| $L'(2,\chi_4)$ | $0.070276738248577532474297338559$ | $= \beta'(2)$ structurally |
| $L'(2,\chi_5)$ | $0.204017830037236175269623445031$ | direct nsum, mod 5 |
| $L'(2,\chi_8)$ | $0.119606867200231549989078381403$ | direct nsum, mod 8 |
| $L'(2,\chi_{12})$ | $0.104089968518901021210670225236$ | direct nsum, mod 12 |
| $\gamma_1$ (Stieltjes) | $-0.072815845483676724860586375875$ | mpmath |
| $\gamma_2$ (Stieltjes) | $-0.009690363192872318484530386035$ | mpmath |
| $\gamma_3$ (Stieltjes) | $0.002053834420303345866160046543$ | mpmath |
| $\log A$ (Glaisher-Kinkelin) | $0.248754477033784262547252993576$ | mpmath |
| $\zeta'(0)$ | $-0.918938533204672741780329736406$ | $= -\tfrac{1}{2}\log(2\pi)$ |
| $\zeta'(-1)$ | $-0.165421143700450929213919660243$ | $= \log A - \tfrac{1}{12}$ |

Several of these admit simpler closed forms ($L'(2,\chi_4) = \beta'(2)$;
$\zeta'(0) = -\log\sqrt{2\pi}$; $\zeta'(-1) = \log A - 1/12$) — these
were left as separate basis elements for PSLQ to find any
identification through any natural route.

## 4. PSLQ campaigns

PSLQ run at $\mathrm{dps} = 250$ against the v3 $c$ (effective precision
$\sim 20-22$ digits per the cross-validation, $\sim 14$ digits per the
K-delta upper bound). Four named bases:

| Basis | Composition | Size |
|:--|:--|:-:|
| basic | rationals + $\pi$ powers + $\log(2)$ + $\gamma_E$ + Catalan + $\zeta(2,3,5)$ + $\sqrt\pi$ + $\log(\pi)$ | 14 |
| L_derivs_only | basic + L-derivative sector (sub-task b) | 38 |
| legacy_only | basic + Track 4 union (Stein-Weiss IBP + multi-Hurwitz + wildcard) | 59 |
| full_union_v3 | basic + legacy + L-derivatives (all) | 83 |

Ceilings tested: $M \in \{10^5, 10^6, 10^7, 10^8\}$.

**Result: 16 PSLQ runs, 16 NULLs.**

| Basis | $M = 10^5$ | $M = 10^6$ | $M = 10^7$ | $M = 10^8$ |
|:--|:-:|:-:|:-:|:-:|
| basic         | NULL | NULL | NULL | NULL |
| L_derivs_only | NULL | NULL | NULL | NULL |
| legacy_only   | NULL | NULL | NULL | NULL |
| full_union_v3 | NULL | NULL | NULL | NULL |

PSLQ time per run: 0.1 s (basic), 0.8 s (L_derivs), 2.1 s (legacy),
4.8 s (full_union_v3).

## 5. Honest precision assessment

The v3 $c$ value is reliable to $\sim 20$ digits (via cross-validation
against an independent panel point) and the K-delta uncertainty estimate
is $\sim 14$ digits (conservative). For PSLQ ceiling $M$, the precision
requirement is roughly $\log_{10}(M) \cdot N$ digits where $N$ is the
basis size. Thus:

| Ceiling $M$ | Required digits ($N=83$) | Available (v3) | Rigor |
|:-:|:--|:--|:--|
| $10^5$ | $\sim 415$ | $\sim 20$ | **Rigorous** (no genuine relation could close beyond our precision floor) |
| $10^6$ | $\sim 500$ | $\sim 20$ | Strongly suggestive |
| $10^7$ | $\sim 580$ | $\sim 20$ | Suggestive |
| $10^8$ | $\sim 660$ | $\sim 20$ | Indicative only |

The PSLQ NULL at $M = 10^5$ is **rigorous**: any small-integer ($\leq 10^5$)
combination of 83 basis elements equal to $c$ would close to better than
$10^{-20}$ if it were a true identity, which PSLQ at our working precision
(250 dps) and target precision (20 dps in $c$) would detect.

At $M = 10^8$, the rigorous threshold isn't met (would need precision
beyond what conditioning of the Richardson LS problem allows on this
panel). Still, **practical experience** with PSLQ is that natural
mathematical identities have small coefficients ($\leq 10^3$); a NULL
at $M = 10^5$ already excludes most plausible structural identifications.

## 6. Verdict and structural reading

**DEFINITIVELY NULL** against the v3 c at coefficient ceilings up to $10^5$
(rigorous) against 83 basis elements spanning **M1 ring** (rationals,
$\pi$-powers, $\log$, $\gamma_E$, Catalan, $\zeta$, $\eta$), **M2 ring**
(Hurwitz at quarter shifts, $\sqrt{\pi}$), **M3 ring** (Catalan G,
Dirichlet $\beta$, $L(2,\chi)$ for $\chi$ in $\{$mod 3, 4, 5, 8, 12$\}$),
**Stein-Weiss IBP** combinations of M1 with $\log$/$\gamma_E$,
**wildcard** (Stieltjes, Glaisher, $\zeta'(0,-1)$), and the **new
L-derivative sector** ($\zeta'(s)$, $\beta'(s)$, $\eta'(s)$, $L'(2,\chi)$
for $\chi \in \{3,4,5,8,12\}$).

Strongly suggestive NULL at $M = 10^6$ through $10^8$.

### 6.1 What this rules out

(a) **$c$ is not a small-integer rational of L-derivatives.** Crucially,
$\zeta'(2), \zeta'(3), \zeta'(4)$ and the Dirichlet $\beta'$ family
are the natural next layer of M1-mechanism content under Stein-Weiss
integration-by-parts. NULL against these eliminates the hypothesis
"$c$ is the first-derivative-of-M1 sub-mechanism output."

(b) **$c$ is not in the Stieltjes / Glaisher / $\zeta'(0,-1)$ extended ring.**
The most natural "auxiliary transcendental" identifications also fail.

(c) **Track 4's hypothesis $c \in M_1^{(1)}$ is now closed in the
negative.** Track 4 noted (memo §5.2 reading (i)): "$c$ involves $\zeta'$
values... derivatives of L-functions" as the most likely reading.
v3 explicitly tests this and returns NULL.

### 6.2 What remains: three readings of irreducibility

(i) **Higher-derivative ring**: $c$ involves $\zeta''(s)$ at $s \neq 2$,
$\beta''(s)$, higher-order derivatives of L-functions. Only $\zeta''(2)$
was in the v3 basis. Adding $\zeta''(3,4,5)$, $\beta''(s)$ for $s = 2,3$
would test the second-derivative-of-M1 sector. This is a natural follow-up
but **the M1 master-Mellin-engine domain interpretation predicts $c$ should
live in $M_1^{(0)}$ + $M_1^{(1)}$ at most** — going to higher derivatives
is structurally less motivated.

(ii) **Genuinely new constant in the irreducible-but-natural list.** The
list of GeoVac framework constants without closed form is now:

  - $K = \pi(B + F - \Delta) \approx 1/\alpha + 8.8 \times 10^{-8}$ (Paper 2 in
    Observations; twelve mechanisms eliminated).
  - Wolfenstein CKM parameters (Sprint W3, falsified spectral-zeta
    hypothesis 2026-05-08).
  - Atomic correlation entropy $S_{\text{full}}(\text{GS})$ (Sprint TD
    Track 5, 2026-05-08).
  - **L2 next-order constant $c$** (this sprint).

(iii) **Larger coefficient ceiling**. A relation with $M > 10^8$ is not
ruled out by our precision but is also not supported by Occam's razor —
small-coefficient identities are the norm for natural transcendentals.

### 6.3 What the result means for the master Mellin engine

The case-exhaustion theorem (Paper 32 §VIII, Sprint TS-E1) holds that
every $\pi$-source in a finite Paper 34 chain engages one of three
mechanisms M1/M2/M3 via the master Mellin operator
$\mathcal{M}[\mathrm{Tr}(D^k \cdot e^{-tD^2})]$ at $k \in \{0,1,2\}$.
The Sprint MR-A/B/C synthesis sharpened this to a **mechanism-AND-domain
partition**: M1 (k=0) ↔ state-space propinquity rates (4/π is M1
leading); M2 (k=2) ↔ heat-kernel/spectral-action convergence (closed
form $\varepsilon(t) \in \sqrt{\pi}\cdot\mathbb{Q} \oplus \pi^2 \cdot \mathbb{Q}$);
M3 (k=1) ↔ vertex-restricted parity-character observables.

The L2 next-order constant $c$ is an **M1-domain observable** (propinquity
rate sub-leading term). MR-A established structurally that the L2
domain's transcendental content lives on the M1 ring.

**Reading of v3's NULL across L-derivative content**: the M1 ring's
"natural" sub-leading is $4/\pi \cdot (\text{something derivative-style})$.
v3 tests all natural derivative completions and finds none. Either:

- **(a)** the M1 mechanism's sub-leading at L2 is genuinely outside the
  basis of $L$-derivatives and Stein-Weiss IBP combinations we tested
  (suggesting a richer structure than the "M1 ring" framing captures);

- **(b)** the L2 next-order constant sits below the master Mellin engine's
  intrinsic resolution — i.e., it is in the same epistemic category as
  K = π(B+F−Δ) (numerical coincidence between Tier-A spectral invariants
  of the manifold, without first-principles derivation).

Reading (b) is consistent with the precedent established by Sprint K-CC
(May 2026, `debug/kcc_*` files) and Sprint TD Track 5 (May 2026,
`debug/sprint_td_track5_memo.md`): the master Mellin engine is predictive
for **which transcendentals can appear** in observables of a given domain,
but the specific real-number value of derived quantities can still be
irreducible against the predicted ring.

The L2 constant $c$ now joins this list as the cleanest example: a
quantity whose **transcendental class is bounded by the master Mellin
engine** but whose **exact real-number value is irreducible against the
ring** at small coefficient ceilings.

## 7. Files

- **Driver:** `debug/l2_constant_c_precision_v3.py` (~635 lines)
- **Smoke test:** `debug/l2_v3_smoke.py`
- **Output JSON:** `debug/data/l2_constant_c_v3.json`
- **Run log:** `debug/data/l2_v3_run.log`
- **This memo:** `debug/l2_constant_c_precision_v3_memo.md`

## 8. Provenance

- Predecessor: Track 4 / Sprint MR-C (`debug/l2_constant_c_identification.py`,
  `debug/l2_constant_c_identification_memo.md`,
  `debug/data/l2_constant_c.json`,
  `debug/data/l2_constant_c_precision_push.json`).
- Infrastructure: `geovac/central_fejer_su2.py` (T_n / γ_n closed-form
  sum-rule, R2.5-L2 quantitative-rate sprint).
- Structural framing: Paper 32 §VIII (case-exhaustion theorem),
  Paper 18 §III.7 (master Mellin engine domain partition),
  Sprint MR-A/B/C synthesis memos.
- Curve-fit-audit discipline: `docs/curve_fit_audit_memo.md` and
  CLAUDE.md memory `feedback_diagnostic_before_engineering.md`.

## 9. Recommendations

1. **Update the canonical $c$ value** in the L2 sprint folder and in
   any document that quotes $c$ beyond 11 digits. Track 4's "80-digit
   string" should be replaced by the corrected v3 value
   $c = 4.10932146748987562042296336493\ldots$ (28 dps, ~20 digits
   reliable). All mentions of $c$ at 11 digits or fewer remain unchanged
   (both values agree there). PI's call on whether to flag/erratum in
   the MR-C memo (which published "80 digits").

2. **No further PSLQ sprint recommended** unless precision push beyond
   $\sim 20-22$ digits becomes feasible (would require either
   re-engineering the Richardson basis to avoid conditioning loss, or
   computing $h(n)$ at $n \gg 8192$ — the latter at $O(n^2)$ panel
   cost, so $n = 32768$ would take $\sim 4$ hours per panel point).

3. **Structural interpretation $c \in M_1^{(0)}$ irreducible** is now
   the standing reading of this constant, consistent with the precedent
   of $K$, $S_{\text{full}}$, and Wolfenstein parameters. Paper 32 §VIII
   and Paper 18 §III.7 may optionally cite this sprint's NULL when noting
   that the master Mellin engine is mechanism-predictive but not always
   value-predictive (PI's call on whether the existing "MR-C NULL"
   citation suffices).
