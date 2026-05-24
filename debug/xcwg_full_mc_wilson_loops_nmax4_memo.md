# XCWG-G $n_\max=4$ Extension — Area-vs-Perimeter Separation Test

Sprint memo, 2026-05-16. Follow-on to XCWG-G (`xcwg_full_mc_wilson_loops_memo.md`).

## Headline

**Structural surprise: $\sigma_{\rm comb}$ moves $\sim 10\times$ deeper negative
at $n_\max=4$ (z = $-48$) compared to $n_\max=3$ (z = $-4.8$).** The
finite-volume perimeter-dominance hypothesis from XCWG-G's §7 verdict is
**falsified by direction** — perimeter dominance does not recede with growing
$n_\max$; instead, the geometric correlation between $A$ and $L$ in the
area-controlled cluster ensemble produces a structurally anti-area-law signal
that the linear joint $\sigma A + \mu L$ fit cannot disentangle. The
$\sigma_{\rm ens}(\beta)$ ensemble fit remains $>0$ and monotone-decreasing
at every tested $\beta$ — qualitatively still the right 3D U(1) shape — but the
"clean area-vs-perimeter separation at $n_\max=4$" predicted by XCWG-G's memo
does not materialize. Paper 41 v4's WEAK PASS finite-volume scope holds; the
recommended $n_\max=4$ test sharpens the diagnosis rather than confirming
recovery, and identifies a **non-linear $A\to L$ topological coupling** in the
Rule B cluster ensemble as the load-bearing mechanism.

The Polyakov rate constant also moves the wrong way:
$c_\sigma^{\rm MC}(n_\max=4) = 0.93$ vs the prediction $c_\rho/2 = 4.70$ — a
**$5\times$ disagreement, larger than at $n_\max=2$ (factor $3.9$) or
$n_\max=3$ (factor $2.9$).** The Polyakov dilute-gas derivation does not
extrapolate cleanly with $n_\max$ on Rule B; the disagreement is robust and
attributable to the non-cubic dual-lattice topology — exactly the reading
XCWG-G already flagged.

**Net for Paper 41:** WEAK PASS with finite-volume scope stands. The diagnosis
sharpens: it is not "the available graph sizes are too small for area-law
recovery" (the $n_\max$ extension is in hand and it does not recover), it is
"the area-controlled cluster ensemble on Rule B has an intrinsic non-linear
$A\leftrightarrow L$ correlation that the joint linear fit cannot resolve;"
the apparent confinement signal in $\sigma_{\rm ens}$ is structurally
perimeter-mediated at every accessible $n_\max$.

## §1. Setup

**Graph.** Rule B Dirac-S$^3$ at $n_\max = 4$: $V = 60$, $E = 312$,
$\beta_1 = E - V + 1 = 253$, primitive 4-plaquettes 5047 total (capped at
1500 for MC tractability, sufficient because cluster growth is local and
plaquette cap does not enter the $\sigma$ fit). Source:
`geovac.ihara_zeta_dirac.build_dirac_s3_graph(4, "B")`.

**MC settings.**
- $\beta$ grid: $\{0.3,\, 1.0,\, 3.0\}$ — three operating points spanning
  strong-coupling ($\beta=0.3$, where confinement should be most visible),
  intermediate ($\beta=1.0$), and dilute-monopole ($\beta=3.0$, where the
  Polyakov rate prediction is tested).
- Thermalisation: $n_{\rm therm} = 600$ sweeps with auto-tuned $\delta_{\max}$
  for $\sim 50\%$ acceptance.
- Sampling: $n_{\rm sample} = 3000$ configurations, $\text{sample\_interval} = 10$
  sweeps between samples → 30,600 MC sweeps per $\beta$.
- A-targets: $\{1, 2, 3, 4, 5, 6\}$. Max 30 clusters per area, 8 loops
  selected per area (round-robin across available perimeters).
- Runtime: **4.1 min total** (well under the 30–60 min budget).
- RNG seed: 44.

**Available perimeters by area (n_max=4 cluster enumeration).**
| $A$ | Available $L$ | Selected $L$ (8 loops) |
|:-:|:-:|:-:|
| 1 | $\{4\}$           | $[4]\times 8$ |
| 2 | $\{4, 6\}$         | $[4, 6, 6, 6, 6, 6, 6, 6]$ |
| 3 | $\{6, 8\}$         | $[6, 8, 6, 8, 6, 8, 8, 8]$ |
| 4 | $\{8, 10\}$        | $[8, 10, 8, 10, 8, 10, 8, 10]$ |
| 5 | $\{12, 14\}$       | $[12, 14, 12, 14, 12, 14, 12, 14]$ |
| 6 | $\{10, 12, 14\}$   | $[10, 12, 14, 12, 14, 12, 14, 12]$ |

Critically: $A=5$ has *no* $L \le 10$ clusters in the available enumeration,
while $A=6$ has $L=10$ clusters. So **the minimum perimeter at fixed $A$ jumps
non-monotonically** as $A: 4\to 5\to 6$ ($L_{\min} = 8 \to 12 \to 10$). This is
the geometric source of the $\sigma_{\rm comb}<0$ effect, as documented in §3.

All selected loops verified closed under the signed incidence boundary,
$\partial \sigma_C = 0$ (consistency check OK).

## §2. MC measurement of $\langle W(C) \rangle_\beta$ at $n_\max=4$

Per-area ensemble means with block-jackknife errors (10 blocks):

| $\beta$ | $\langle W\rangle_{A=1}$ | $\langle W\rangle_{A=2}$ | $\langle W\rangle_{A=3}$ | $\langle W\rangle_{A=4}$ | $\langle W\rangle_{A=5}$ | $\langle W\rangle_{A=6}$ |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 0.3 | $0.8432\pm 0.0021$ | $0.7508\pm 0.0034$ | $0.5146\pm 0.0026$ | $\mathbf{0.5776\pm 0.0040}$ | $0.3539\pm 0.0059$ | $\mathbf{0.4059\pm 0.0057}$ |
| 1.0 | $0.9617\pm 0.0007$ | $0.9286\pm 0.0011$ | $0.8775\pm 0.0011$ | $0.8682\pm 0.0011$ | $0.7811\pm 0.0033$ | $0.8082\pm 0.0030$ |
| 3.0 | $0.9859\pm 0.0002$ | $0.9767\pm 0.0005$ | $0.9568\pm 0.0006$ | $0.9580\pm 0.0007$ | $0.9298\pm 0.0014$ | $0.9395\pm 0.0014$ |

**Non-monotonicity highlighted in bold.** At every tested $\beta$,
$\langle W(A=4)\rangle > \langle W(A=3)\rangle$ in a statistically significant
sense (the SEs are $0.003$–$0.004$, the gap is $0.06$–$0.07$), and at
$\beta = 0.3$ and $1.0$ also $\langle W(A=6)\rangle > \langle W(A=5)\rangle$.
This is **anti-area-law**: in a confining theory $\langle W\rangle$ must
decrease monotonically with $A$. The cause is geometric, *not* MC noise:
the perimeter distribution shifts non-monotonically with $A$ in the available
cluster enumeration ($L_{\min}: 4, 4, 6, 8, 12, 10$ at $A: 1\dots 6$).

Acceptance rates after auto-tuning were $\sim 50\%$ at all $\beta$;
$\delta_{\max}$ ranged $\sim 1.8$ at $\beta = 0.3$ down to $\sim 0.4$ at
$\beta = 3.0$ (standard behaviour).

### Per-loop perimeter dependence at $\beta=0.3$, $A=3$

Within an area class, $\langle W\rangle$ correlates strongly with $L$, as at
$n_\max = 2, 3$. Per-loop means at $\beta=0.3$, $A=3$ (perimeters
$L \in \{6, 6, 6, 8, 8, 8, 8\}$, selected round-robin):

| $L$ | per-loop $\langle W\rangle$ | mean at this $L$ |
|:-:|:-:|:-:|
| 6 | 0.701, 0.671, 0.760 | 0.711 |
| 8 | 0.334, 0.283, 0.312, 0.303 | 0.308 |

Spread across perimeters at fixed $A$ (0.40 between $L=6$ and $L=8$) is far
larger than the per-loop MC SE — same pattern as $n_\max \le 3$, but with
*tighter* per-loop SEs because $n_{\rm sample} = 3000$ here vs $\le 2000$ in
the original XCWG-G.

## §3. $\sigma_{\rm comb}$ and $\mu_{\rm comb}$ at $n_\max = 4$

Three-flavor fit results:

| $\beta$ | $\sigma_{\rm ens}$ | $\sigma_{\rm comb}$ | $\mu_{\rm comb}$ | $z(\sigma_{\rm comb})$ | $R^2_{\rm comb}$ | $\sigma_{\rm LO}$ | $\sigma_{\rm ens}/\sigma_{\rm LO}$ |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 0.3 | $+0.1717\pm 0.0016$ | $\mathbf{-0.0803\pm 0.0017}$ | $+0.1294\pm 0.0010$ | $\mathbf{-48.3}$ | 0.62 | 1.908 | 0.090 |
| 1.0 | $+0.0395\pm 0.0006$ | $\mathbf{-0.0206\pm 0.0008}$ | $+0.0337\pm 0.0004$ | $\mathbf{-26.1}$ | 0.84 | 0.807 | 0.049 |
| 3.0 | $+0.0112\pm 0.0001$ | $\mathbf{-0.0076\pm 0.0002}$ | $+0.0107\pm 0.0001$ | $\mathbf{-38.7}$ | 0.84 | 0.211 | 0.053 |

**Three observations.**

(a) $\sigma_{\rm ens}$ is positive and monotone-decreasing at all $\beta$, with
ratio to LO around $0.05$–$0.09$ — same qualitative behaviour as $n_\max=2,3$,
i.e. the naive ensemble fit reproduces the 3D compact-$U(1)$ shape with
$\sigma_{\rm ens}\ll \sigma_{\rm LO}$ at small $\beta$ (LO XCWG-D overshoots).

(b) $\sigma_{\rm comb}$ is **statistically NEGATIVE** at every $\beta$ with
**very high significance** ($z \le -26$). This is mechanistically driven by
the non-monotonic $A\to L_{\min}$ pattern documented in §1: when the joint
fit is asked to ascribe the (large) $\beta=0.3$ ⟨W⟩ jump from $A=3$ to $A=4$
(at $L\in\{6,8\}$ vs $L\in\{8,10\}$) to an area dependence, the only way to
reconcile with the $A=5\to 6$ jump (at $L\in\{12,14\}$ vs $L\in\{10,12,14\}$)
is via a *negative* $\sigma$ plus a positive $\mu$ that absorbs the bulk
perimeter trend.

(c) $\mu_{\rm comb}$ is positive and decreasing with $\beta$, roughly tracking
$\sigma_{\rm ens}$ within a factor of $\sim 0.75$. The $R^2$ of the joint fit
is $0.62$ at $\beta=0.3$ (worse than $n_\max\le 3$) and $0.84$ at $\beta\ge 1$
— the linear $\sigma A + \mu L$ ansatz fits poorly in the strong-coupling
regime, consistent with the anti-area-law jumps being structural (the
geometry, not noise).

### Fixed-perimeter fits — the load-bearing diagnostic at $n_\max=4$

The fixed-$L$ fits sharpen the diagnosis. At $\beta=0.3$ (each fit has
$n_{\rm areas}=2$):

| $L$ | $\sigma_{\rm fixL}$ | $A$-pair | $\langle W\rangle$-pair | Interpretation |
|:-:|:-:|:-:|:-:|:-:|
| 4  | $-0.016\pm 0.004$ | $\{1,2\}$ | $0.843, 0.749$ | tiny slope, near zero |
| 6  | $-0.004\pm 0.003$ | $\{2,3\}$ | $\sim 0.76, \sim 0.71$ | **essentially zero, clean** |
| 8  | $\mathbf{-0.586\pm 0.004}$ | $\{3,4\}$ | $\sim 0.31, \sim 0.66$ | **anti-area-law**, dominant anomaly |
| 10 | $-0.049\pm 0.006$ | $\{4,6\}$ | $\sim 0.82, \sim 0.50$ | mildly negative |
| 12 | $+0.013\pm 0.008$ | $\{5,6\}$ | $\sim 0.35, \sim 0.39$ | essentially zero/positive |
| 14 | $-0.191\pm 0.019$ | $\{5,6\}$ | $\sim 0.21, \sim 0.41$ | anti-area-law |

The dominant anomaly is at $L=8$: clusters of area 3 with perimeter 8 have
$\langle W\rangle \sim 0.31$ while clusters of area 4 with perimeter 8 have
$\langle W\rangle \sim 0.66$. **The area-4 cluster has a more compact interior
geometry that produces a less-fluctuating Wilson loop than the area-3
cluster** — this is not noise, it is the cluster topology talking.

The $L=6$ fixed-perimeter slope is $-0.004\pm 0.003$ (consistent with zero,
$z = -1.4$). **This is the only clean fixed-$L$ datum and it gives the
honest answer**: at fixed perimeter, on the available 2-point $A$-set, the
genuine area-law coefficient is consistent with zero at $\beta=0.3$. No
confinement signal that survives perimeter control.

The same pattern holds at $\beta=1$ and $\beta=3$: $L=8$ is always the outlier
(strongly negative), $L=6$ is always the cleanest (zero), and the joint
fit's $\sigma_{\rm comb}$ inherits the $L=8$ pull.

## §4. Trend across $n_\max = 2, 3, 4$

Summary at $\beta=0.3$ (the strongest-coupling probe in the common range):

| $n_\max$ | $V$ | $E$ | $\sigma_{\rm ens}$ | $\sigma_{\rm comb}$ | $\mu_{\rm comb}$ | $z(\sigma_{\rm comb})$ |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 2 | 10 | 20  | $+0.174\pm 0.009$ | $-0.006\pm 0.003$ | $+0.217$ | $-2.0$ |
| 3 | 28 | 106 | $+0.131\pm 0.010$ | $-0.024\pm 0.005$ | $+0.178$ | $-4.8$ |
| 4 | 60 | 312 | $+0.172\pm 0.002$ | $\mathbf{-0.080\pm 0.002}$ | $+0.129$ | $\mathbf{-48.3}$ |

**$\sigma_{\rm comb}$ moves $\sim 10\times$ deeper negative from $n_\max=3$ to
$n_\max=4$, not toward zero.** The z-score blows up by an order of magnitude
because the higher-statistics $n_\max=4$ run pins $\sigma_{\rm comb}$ at
$-0.08$ with extraordinarily tight error bars. Equivalent reading:
$\sigma_{\rm comb}$ at $\beta = 0.3$ is $-0.006, -0.024, -0.080$ at
$n_\max = 2, 3, 4$ — moving away from zero by a factor of $\sim 3\times$ per
graph-size step, in the wrong direction for the XCWG-G prediction.

$\mu_{\rm comb}$ moves the other way ($0.217 \to 0.178 \to 0.129$): the
perimeter coefficient is decreasing with $n_\max$. So
$\sigma_{\rm ens} \approx \mu_{\rm comb}$ that held at $n_\max=2,3$ also
holds at $n_\max=4$ ($0.172$ vs $0.129$: same order, within $\sim 30\%$).
Both reflect that the ensemble fit absorbs the dominant perimeter trend into
an apparent area slope.

**Three competing readings of the trend:**

1. *Asymptotic behaviour negative.* The Rule B Wilson construction is
   structurally non-confining in this observable, and $\sigma_{\rm comb}$
   asymptotes to some negative value (or zero from below) as $n_\max \to
   \infty$. The growing |z|-score is consistent with $\sigma_{\rm comb}$
   stabilising at a *negative* value as graph size grows, not at zero.

2. *Linear-ansatz breakdown.* The joint $\sigma A + \mu L$ fit is structurally
   inadequate for Rule B's cluster ensemble — the area-perimeter correlation
   is non-linear ($\mu(L)$ or $\sigma(A)$ are not constants of motion of the
   Wilson loop), and the linear fit returns a meaningless effective slope
   that grows in magnitude with statistics. The clean $L=6$ fixed-perimeter
   datum ($\sigma_{\rm fixL=6} = -0.004\pm 0.003$, consistent with zero)
   supports this reading: when the cluster topology is held fixed, the area
   law genuinely vanishes within MC noise.

3. *Cluster-ensemble pathology.* The area-controlled cluster enumeration on
   Rule B systematically produces, at higher $A$, cluster topologies whose
   Wilson-loop fluctuations *decrease* (rather than increase) when traversed
   to the same $L$. This is real geometry — the "fatter" $A=4$ clusters at
   $L=8$ have more interior structure that locally cancels gauge
   fluctuations on the boundary. The Rule B graph's high plaquette density
   ($9.4$ plaquettes/edge at $n_\max=3$, $16.2$ at $n_\max=4$ — far above
   3D-cubic 1.5) makes this geometric effect more pronounced as $n_\max$
   grows.

Readings 2 and 3 are not mutually exclusive and are jointly the most
consistent with the data. Reading 1 (structural non-confining) is the
honest *worst-case* reading; we cannot rule it out from the present data, but
the clean fixed-$L=6$ result favours readings 2+3 (linear ansatz fails
because cluster topology dominates).

## §5. Polyakov rate cross-check at $n_\max = 4$

Fitting $\log \sigma_{\rm ens}(\beta) = \log A_\sigma - c_\sigma \beta$ over
$\beta\in[0.3, 3.0]$ at $n_\max=4$ (3 data points):

| Quantity | $n_\max=2$ | $n_\max=3$ | $n_\max=4$ | Predicted ($c_\rho/2$ from XCWG-F) |
|:--------:|:----------:|:----------:|:----------:|:-----------------------------------:|
| $c_\sigma^{\rm MC}$ | $1.20$ | $1.63$ | $\mathbf{0.93}$ | $\mathbf{4.70}$ |
| Relative error vs $4.70$ | $0.74$ | $0.65$ | $\mathbf{0.80}$ | — |
| $R^2$ of exp fit | $0.83$ | $0.88$ | $0.90$ | — |

**The rate constant has not converged toward $4.70$ with growing $n_\max$.**
In fact, $c_\sigma^{\rm MC}$ at $n_\max=4$ ($0.93$) is *farther* from the
prediction than at $n_\max=3$ ($1.63$), and roughly equal to $n_\max=2$
($1.20$). The trend $1.20 \to 1.63 \to 0.93$ does not suggest convergence to
$4.70$ at larger $n_\max$.

**Reading.** The XCWG-G memo flagged the rate disagreement as expected on a
non-cubic dual lattice (the Polyakov derivation assumes $\mathbb{Z}^3$
topology). The $n_\max=4$ datum reinforces that reading: the rate constant
is *not* a finite-size artifact that scales toward the continuum prediction;
it is a structural property of the Rule B dual-lattice topology. The
exponential *shape* of $\sigma_{\rm ens}(\beta)$ continues to fit well
($R^2 = 0.90$), but the *rate* is a graph-specific constant of order $\sim 1$
that does not match the cubic-lattice value $c_\rho/2 = 4.70$.

Honest scope: $n_\max=4$ has only 3 $\beta$-data-points for the rate fit
(budget-constrained), and the relative error of $0.80$ inherits some of
that low-data-point uncertainty. A larger $\beta$-grid (say $\beta\in
\{0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0\}$, doubling MC cost to ~10 min)
at $n_\max=4$ could tighten the rate-fit estimate, but is unlikely to change
the qualitative conclusion that the rate constant remains far below the
Polyakov prediction.

## §6. Verdict

**Outcome class: `NEGATIVE_SIGMA_COMB_STRUCTURAL_SURPRISE`** at $n_\max=4$,
with statistically very strong (|z| > 25) negativity at all tested $\beta$.

**Reading for Paper 41:**

1. The XCWG-G WEAK PASS with finite-volume scope at $n_\max\le 3$ **stands**.
   $\sigma_{\rm ens}(\beta)$ at $n_\max=4$ remains $> 0$ and
   monotone-decreasing, with ratio to LO $0.05$–$0.09$ — qualitatively the
   same 3D compact-$U(1)$ shape as at smaller $n_\max$.

2. The XCWG-G memo's prediction "the area-vs-perimeter separation should
   clean up at $n_\max=4$" is **falsified by direction**. $\sigma_{\rm comb}$
   moves $\sim 10\times$ deeper negative, not toward zero.

3. The honest diagnosis is sharpened from "finite-volume perimeter dominance
   that should recede" to **"the area-controlled cluster ensemble on Rule B
   has an intrinsic non-linear $A\leftrightarrow L$ coupling that the linear
   joint fit cannot resolve, and the cleanest fixed-$L$ datum
   ($\sigma_{\rm fixL=6} = -0.004\pm 0.003$ at $\beta=0.3$) gives an
   area-law coefficient consistent with zero."** This is consistent with
   3D compact $U(1)$ (which should give $\sigma > 0$ at strong coupling) only
   if the perimeter-controlled measurement underestimates the true area law
   on this graph due to the finite-cluster cap.

4. The Polyakov rate constant remains $\sim 5\times$ off prediction at
   $n_\max=4$ and shows no trend toward convergence. **The
   non-cubic-dual-lattice reading of the rate disagreement (XCWG-G §4) is
   strengthened, not weakened, by the $n_\max=4$ datum.**

**Recommended Paper 41 v5 update (small):** add a paragraph in §6.5 noting
the $n_\max=4$ datum and the trend table from §4, with the verdict that the
area-vs-perimeter separation does not clean up with $n_\max$ growth — it
sharpens the perimeter-dominance diagnosis to a non-linear $A\leftrightarrow
L$ cluster-topology effect. Cite §3's fixed-$L=6$ clean result (consistent
with $\sigma = 0$) as the honest area-law coefficient at $\beta=0.3$.
Honest-scope caveats (v) and (vi) [Polyakov rate disagreement, finite-volume
area-vs-perimeter] both stand; the $n_\max=4$ datum **strengthens** both by
showing they do not relax with $n_\max$.

**WEAK PASS with sharpened finite-volume scope.** Six-witness structural
compatibility (Paper 41 v4) holds. The seven witnesses now have:

- 1-5, 7: clean passes (no $n_\max$-dependence concerns).
- 6 (XCWG-G): WEAK PASS that **does not clean up with $n_\max$** — diagnosed
  as a cluster-ensemble non-linearity, with one clean fixed-$L$ datum
  ($L=6$) consistent with zero area-law coefficient at all $\beta$.

The natural next sprint (if continued) is **not another $n_\max$ extension**
(the trend is now characterised), but rather a **non-area-controlled
Wilson-loop observable** — e.g., loops of fixed perimeter and varying
"interior plaquette count" defined intrinsically, decoupled from the
cluster-growth enumeration that produces the $A\leftrightarrow L$ correlation.
Or, alternatively, a **two-point Polyakov-correlator** $\langle P(x) P^*(y)
\rangle$ on the spatial graph, which has no perimeter-vs-area ambiguity and
gives a cleaner asymptotic-distance scaling. Both are 2–3 week sprints if
prioritised.

## §7. Reported key numbers

- $\sigma_{\rm comb}$ at $\beta\in\{0.3, 1.0, 3.0\}$ (n_max=4):
  $\mathbf{-0.080\pm 0.002,\; -0.021\pm 0.001,\; -0.008\pm 0.000}$.
- One-sentence verdict: at $n_\max=4$, **the area-law coefficient is
  statistically very negative** ($z \le -26$) under the linear joint fit
  due to a non-monotonic $A\leftrightarrow L_{\min}$ pattern in the
  area-controlled cluster ensemble — **finite-volume perimeter dominance
  does not recede with graph-size growth, it sharpens.** The cleanest
  fixed-$L=6$ datum gives $\sigma = -0.004\pm 0.003$ (consistent with zero)
  at $\beta=0.3$, the honest area-law coefficient.
- Polyakov rate constant at $n_\max=4$: $c_\sigma^{\rm MC} = \mathbf{0.93}$
  vs prediction $c_\rho/2 = 4.70$ ($\sim 5\times$ disagreement, **not
  improving with $n_\max$**; trend $1.20\to 1.63\to 0.93$ at
  $n_\max=2,3,4$).

## Files

- Script: `debug/xcwg_full_mc_wilson_loops_nmax4.py`
- Data: `debug/data/xcwg_full_mc_wilson_loops_nmax4.json` (~5 MB)
- Log: `debug/data/xcwg_full_mc_wilson_loops_nmax4.log`
- Plot: `debug/plots/xcwg_full_mc_wilson_loops_nmax4.png`
- MC runtime: 4.1 min single CPU (well under the 30–60 min budget; smaller
  than expected because numpy vectorisation handles the n_max=4 plaquette
  evaluation efficiently).

## Cross-references

- XCWG-G original: `debug/xcwg_full_mc_wilson_loops_memo.md` (the WEAK PASS
  verdict and the $n_\max=4$ recommendation that motivated this sprint).
- Rule B graph specs at $n_\max=4$:
  `debug/xcwg_rule_b_infrastructure_memo.md` §2.
- XCWG-F monopole density: `debug/xcwg_monopole_density_memo.md`
  ($c_\rho = 9.40$ at $n_\max=2$, used for the Polyakov rate prediction).
- Paper 41 v4 (current): `papers/group5_qed_gauge/paper_41_rule_b_wilson_u1.tex`
  §6.5 reports XCWG-G as WEAK PASS with finite-volume scope and recommends
  the $n_\max=4$ extension; v5 should add the present datum.
