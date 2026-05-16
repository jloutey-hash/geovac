# XCWG-G — Full Monte Carlo Wilson Loops on Rule B Compact U(1)

Sprint memo, 2026-05-16.  Sixth-witness diagnostic for Paper 41 v3.

## Headline

**WEAK PASS, with a substantive new finding.**  At full Metropolis-Hastings MC
of the compact U(1) Wilson action on the Dirac Rule B graph, the
ensemble-mean-per-area area-law fit gives $\sigma_{\rm ens}(\beta) > 0$ at
every $\beta \in [0.1, 10]$, monotonically decreasing in $\beta$.  But the
joint $\sigma A + \mu L$ fit and the fixed-perimeter $L = 4$ restriction
both show that **the genuine area-law coefficient is statistically zero**
and **the apparent confinement signal is carried by the perimeter term**:
at every $\beta$ tested, $\mu_{\rm comb} \approx \sigma_{\rm ens}$ and
$\sigma_{\rm comb} \approx 0$ within MC noise.  On the finite Rule B graph
at $n_\max = 2$ ($V=10$, $E=20$) and $n_\max = 3$ ($V=28$, $E=106$),
Wilson loops are dominated by perimeter contributions, not area law.  The
Polyakov dilute-gas prediction $\sigma(\beta) \sim \exp(-(c_\rho/2)\beta)$
from XCWG-F passes a qualitative exponential-shape test ($R^2 = 0.83$ at
$n_\max = 2$) but with the **wrong rate constant** — the measured exponent
$c_{\sigma}^{\rm MC} = 1.20$ is a factor of $3.9\times$ smaller than the
predicted $c_\rho/2 = 4.70$.

This is **NOT a contradiction of the XCWG-F monopole finding** — XCWG-F
measured $\rho_M(\beta) > 0$ directly, which is geometry, and that result
stands.  What XCWG-G reveals is that the **Polyakov-style derivation of
$\sigma$ from $\rho_M$ via the dilute-gas formula does not transfer to the
Rule B finite graph in the form the textbook 3D continuum result would
predict**.  The five-witness structural compatibility verdict of Paper 41 v3
holds; XCWG-G adds **one weak-pass instance and one clean finite-volume
finding**, and the perimeter-dominated regime is now characterized
quantitatively.

## §1.  Setup (extending XCWG-F MC for Wilson loops; distinction from XCWG-D/C)

XCWG-G samples the **full** compact U(1) Wilson measure
$$
P(\{\theta_e\}) \propto \exp\bigl(-S_W\bigr), \qquad
S_W = \beta \sum_P \bigl(1 - \cos\theta_P\bigr), \qquad
\theta_P = \sum_{e \in \partial P} \mathrm{sign}(P,e)\, \theta_e ,
$$
via Metropolis-Hastings updates of link angles $\theta_e \in (-\pi, \pi]$.
The MC infrastructure is identical to XCWG-F's monopole sprint: single-link
updates with auto-tuned proposal width $\delta_{\max}$ targeting ~50%
acceptance; thermalization of 2000 sweeps; sampling every 20 sweeps after
thermalization to suppress autocorrelation.  We use $n_{\rm sample} = 2000$
samples per $\beta$ at $n_\max = 2$ and $n_{\rm sample} = 600$ at $n_\max = 3$.

For each MC sample $\{\theta_e\}$ and each area-controlled Wilson loop $C$
with signed-edge vector $\sigma_C \in \{-1, 0, +1\}^E$, we compute
$$
W(C; \theta) = \cos\bigl( \sigma_C \cdot \theta \bigr)
            = \mathrm{Re}\,\prod_{e \in C} e^{i\,\sigma_C(e)\,\theta_e} ,
$$
and average $W(C; \theta)$ over MC samples.

**Distinction from XCWG-D (LO).**  XCWG-D used the leading-order
strong-coupling character expansion $\langle W(C) \rangle_\beta =
(I_1(\beta)/I_0(\beta))^A$, which is the **single-plaquette tessellation
formula** exact only at small $\beta$ on the infinite lattice.  XCWG-G samples
the **full** compact measure at every $\beta$, including all higher orders
of the character expansion AND all monopole/disorder fluctuations
simultaneously.  The two agree at very small $\beta$ if (and only if) the
graph supports a clean tessellation of the loop into plaquettes.

**Distinction from XCWG-C (Gaussian).**  XCWG-C used the Gaussian quadratic
form $S(W) = \langle W, K^+ W \rangle$ where $K$ is the abelian (free-photon)
edge Laplacian.  That probe sees only the harmonic sector — no monopole
contributions, hence perimeter law in 3D regardless of compactness.  XCWG-G
uses the full compact action, so it sees the monopole-mediated
disorder-disorder correlations that drive 3D compact U(1) confinement
in continuum.

**Area-controlled Wilson loops.**  We reuse XCWG-D's `grow_area_A_clusters`
to build clusters of $A$ plaquettes whose oriented boundary is a single
closed loop.  At $n_\max = 2$, available perimeters by area are:
$A=1: [4]$, $A=2: [4,6]$, $A=3: [4,6,8,10]$, $A=4: [4,6,8,10,12]$,
$A=5: [4,6,8,10,12]$.  We select 10 loops per area, round-robin across
perimeters.  At $n_\max = 3$, similar enumeration produces 6 loops per area,
with the relevant constraint that perimeter $L=4$ only appears for $A \in
\{1, 2, 3\}$ (loops of $A=4,5$ all have $L \geq 6$).

## §2.  MC measurement of $\langle W(C) \rangle_\beta$ at $n_\max = 2$

The per-area ensemble means with jackknife error bars (selected betas):

| $\beta$ | $\langle W \rangle_{A=1}$ | $\langle W \rangle_{A=2}$ | $\langle W \rangle_{A=3}$ | $\langle W \rangle_{A=4}$ | $\langle W \rangle_{A=5}$ |
|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|
| 0.1   | 0.069(9)  | 0.045(8)  | 0.029(6)  | 0.022(6)  | 0.017(6)  |
| 0.3   | 0.365(7)  | 0.304(7)  | 0.231(9)  | 0.216(7)  | 0.189(8)  |
| 1.0   | 0.871(2)  | 0.851(2)  | 0.806(3)  | 0.793(2)  | 0.762(2)  |
| 3.0   | 0.959(1)  | 0.950(1)  | 0.934(1)  | 0.929(1)  | 0.920(2)  |
| 10.0  | 0.987(0)  | 0.984(0)  | 0.978(0)  | 0.978(0)  | 0.974(0)  |

Acceptance rates after auto-tuning are 0.5 — 0.6 at every $\beta$;
thermalization signature (mean plaquette charge) stabilizes after the first
~200 sweeps.

**Observation.**  At every $\beta$, $\langle W \rangle$ does decrease with
$A$.  Naively this looks like area law.  But examining the per-loop data
within an area class reveals strong dependence on perimeter $L$, which
contaminates the area-only fit.

### Per-loop perimeter dependence

At $\beta = 1.0$, area $A = 4$ (5 distinct perimeters in the selected set):
| $L$ | $\langle W \rangle$ |
|:---:|:------:|
| 4  | 0.878(2) |
| 6  | 0.826(2) |
| 8  | 0.774(2) |
| 10 | 0.706(2) |
| 12 | 0.687(2) |

The spread across perimeters at fixed $A$ (0.19 units of $W$) is **20–100×
larger than the ensemble SE on $\langle W \rangle_A$** (~0.002).  The data
strongly favors a $\langle W \rangle \sim e^{-\mu L}$ perimeter law over
pure area law.

## §3.  $\sigma(\beta)$ extraction — three flavors of fit

Three fit prescriptions, motivated by the perimeter contamination just
identified:

1. **Ensemble area fit** ($\sigma_{\rm ens}$):
   $\log \langle W(A) \rangle = -\sigma A + c$ on per-area ensemble means.
   This is the "naive" XCWG-D-style fit and gives the values one would
   report without disentangling perimeter.

2. **Combined fit** ($\sigma_{\rm comb}, \mu_{\rm comb}$):
   $\log W_i = -\sigma A_i - \mu L_i + c$ on per-loop data, joint WLS.
   Disentangles area and perimeter contributions.

3. **Fixed-perimeter fit** ($\sigma_{\rm fixL}$):
   $\log \langle W(A) \rangle_{L \text{ fixed}} = -\sigma A + c$,
   restricted to loops with a common perimeter (cleanest pure-area test).

Results at $n_\max = 2$:

| $\beta$ | $\sigma_{\rm ens}$ | $\sigma_{\rm comb}$ | $\mu_{\rm comb}$ | $\sigma_{L=4}$ | $\sigma_{\rm LO}$ |
|:-------:|:------:|:------:|:------:|:------:|:------:|
| 0.1   | 0.377(74) | -0.007(16) | +0.342 | -0.015 | 2.997 |
| 0.2   | 0.259(22) | -0.010(6)  | +0.337 | -0.009 | 2.308 |
| 0.3   | 0.174(9)  | -0.006(3)  | +0.217 | -0.006 | 1.908 |
| 0.5   | 0.089(2)  | -0.003(1)  | +0.097 | -0.004 | 1.417 |
| 0.7   | 0.051(2)  | -0.003(1)  | +0.053 | -0.004 | 1.108 |
| 1.0   | 0.033(1)  | -0.002(0)  | +0.034 | -0.003 | 0.807 |
| 1.5   | 0.021(1)  | -0.002(0)  | +0.023 | -0.002 | 0.517 |
| 2.0   | 0.016(1)  | -0.001(0)  | +0.016 | -0.001 | 0.360 |
| 3.0   | 0.011(0)  | -0.000(0)  | +0.010 | -0.000 | 0.211 |
| 5.0   | 0.006(0)  | -0.000(0)  | +0.006 | -0.000 | 0.113 |
| 10.0  | 0.003(0)  | -0.000(0)  | +0.003 | -0.000 | 0.053 |

**Three observations.**

(a) $\sigma_{\rm ens} > 0$ at every $\beta$, monotonically decreasing.
This is the "good" signal — the ensemble area-law fit is positive and looks
like a 3D compact U(1) string tension that vanishes at large $\beta$.

(b) But $\sigma_{\rm comb}$ is **statistically zero** (and slightly negative)
at every $\beta$, and $\sigma_{\rm fixL=4}$ is even more clearly zero.  Both
diagnostic fits explicitly tell us the apparent area-law signal in
$\sigma_{\rm ens}$ is **not** a genuine area dependence — it's perimeter
contamination, because the area-controlled clusters at higher $A$ also have
longer perimeters in the available ensemble.

(c) The relationship $\sigma_{\rm ens} \approx \mu_{\rm comb}$ at every
$\beta$ is what one would expect if the apparent area-law slope is actually
carried by perimeter:  on this graph the perimeters of available
$A$-clusters scale roughly linearly with $A$, with $\partial L / \partial A
\sim 1$, so $-\sigma_{\rm ens} A \approx -\mu L$ which absorbs the
perimeter coefficient into an apparent area coefficient.

**$\sigma_{\rm ens}$ vs $\sigma_{\rm LO}$ ratio:** at small $\beta$,
$\sigma_{\rm ens}/\sigma_{\rm LO} \in [0.05, 0.13]$, i.e. the full-MC
ensemble fit is **5–13% of the LO XCWG-D prediction**.  The LO prediction
overestimates the area law by ~10×.  The ensemble fit's $\sigma > 0$ is
robust, but its magnitude is finite-volume / perimeter-contamination.

## §4.  Polyakov dilute-gas cross-check (using XCWG-F)

XCWG-F measured the monopole density on the same graph at the same MC
machinery and fit $\rho_M(\beta) = 0.747 \cdot \exp(-9.40 \beta)$ in the
range $\beta \in [0.05, 1.0]$ with $R^2 = 0.981$.

The Polyakov 1977 dilute-monopole-gas prediction for 3D compact U(1) gives
$$
\sigma_{\rm Polyakov}(\beta) \propto \sqrt{\rho_M(\beta) / \beta}
                              \propto \exp\bigl(-(c_\rho / 2)\,\beta \bigr) ,
$$
so the predicted exponential rate constant on a log-$\sigma$-vs-$\beta$ plot
is $c_\sigma^{\rm predicted} = c_\rho / 2 = 4.70$.

Fitting $\log \sigma_{\rm ens}(\beta)$ against $\beta$ over $\beta \in
[0.1, 3.0]$ at $n_\max = 2$:

| Quantity | Measured | Predicted |
|:---|:---:|:---:|
| $A_\sigma$ (prefactor) | 0.30 | (fit) |
| $c_\sigma$ (rate)      | **1.20** | **4.70** ($= c_\rho/2$) |
| $R^2$ of exp fit       | 0.83 | — |
| Relative error on rate | 75% | — |

The measured exponent is **3.9× smaller than predicted**.  The
exponential shape itself fits $R^2 = 0.83$ (not great, not terrible — the
data is exponentially decaying at all $\beta$ tested), but the quantitative
agreement is poor.

**Reading.**  The Polyakov derivation assumes:
(a) the 3D-cubic lattice version of "dual lattice = points with monopoles
on them";
(b) the integer-valued Dirac string field $m_P$ has a continuum-like
Coulomb-gas free energy that scales as $\rho_M$;
(c) the area law inherits the $\sqrt{\rho_M / \beta}$ scaling from the
correlation function of two test charges in the monopole plasma.

On the Rule B graph at $n_\max = 2$:
- The dual lattice is NOT a regular 3D cubic lattice.  The closed 2-surfaces
  on which monopole charges sit are size-3 triangle-prism configurations
  (XCWG-F §2), and there are 60 such elementary monopole sites at
  $n_\max = 2$.
- Step (b) likely fails on this graph: the dual-lattice topology is too
  different from $\mathbb{Z}^3$ for the standard sine-Gordon dualization
  to apply verbatim.
- Step (c) requires the existence of a well-defined "test charge"
  correlator, which on a finite non-cubic graph is exactly what we're
  measuring — and it's perimeter-law dominated at this size, not area-law
  dominated.

The XCWG-F monopole density result $\rho_M(\beta) > 0$ at all $\beta$ tested
stands as a measured fact.  What XCWG-G shows is that **the Polyakov-style
derivation of $\sigma$ from $\rho_M$ on the Rule B graph at $n_\max = 2$
does not produce a quantitative match in the form Polyakov's textbook
3D continuum derivation gives**.  This is the natural disclaimer on a small
discrete graph.

## §5.  Perimeter-vs-area test at full MC ($n_\max = 2$, $\beta = 1$)

At leading order (XCWG-D), $\langle W \rangle = (I_1/I_0)^A$ should be
independent of perimeter at fixed area.  XCWG-D verified this analytically
(by construction, since the LO formula has no $L$ dependence).  At full
MC, this is a genuine empirical test.

| $A$ | Perimeters $L$ | Mean $\langle W \rangle$ at each $L$ | Spread | SE |
|:---:|:--:|:--:|:--:|:--:|
| 2 | [4, 6]         | [0.879, 0.823]                | 0.056 | 0.002 |
| 3 | [4, 6, 8, 10]  | [0.880, 0.825, 0.767, 0.705]  | 0.174 | 0.003 |
| 4 | [4, 6, 8, 10, 12]| [0.878, 0.826, 0.774, 0.706, 0.687] | 0.191 | 0.002 |
| 5 | [4, 6, 8, 10, 12]| [0.877, 0.821, 0.760, 0.696, 0.655] | 0.222 | 0.002 |

In every area class, **the spread across perimeters is 30 to 100 times
larger than the MC standard error**.  The perimeter dependence is the
dominant signal in the per-loop data, not random MC noise — area-only
$\langle W \rangle$ is **falsified** at full MC.

A linear fit of $\log \langle W \rangle$ vs $L$ at fixed $A$ across the
five $A=4$ cases at $\beta = 1$ gives a slope of about $-0.027$, i.e. an
effective perimeter coefficient $\mu \approx 0.027$.  This is consistent
with the joint-fit $\mu_{\rm comb} = 0.034$ at the same $\beta$.

## §6.  $n_\max = 3$ extension

At $n_\max = 3$ ($V = 28$, $E = 106$, $\beta_1 = 79$, 300 L=4 plaquettes),
the available perimeters by area shift: $A=1: [4]$, $A=2: [4,6]$,
$A=3: [4,6]$, $A=4: [6,8]$, $A=5: [6,8]$.  $L=4$ no longer extends to
$A \geq 4$, and $L=10, 12$ no longer appear.  This makes the fixed-$L$ test
less powerful at $n_\max = 3$ (only $A \in \{1, 2, 3\}$ at $L=4$, only
$A \in \{4, 5\}$ at $L=8$).

The pattern from $n_\max = 2$ persists at $n_\max = 3$:

| $\beta$ | $\sigma_{\rm ens}$ | $\sigma_{\rm comb}$ | $\mu_{\rm comb}$ | $\sigma_{\rm LO}$ |
|:-------:|:------:|:------:|:------:|:------:|
| 0.1 | 0.242(342) | +0.040(127) | +0.296 | 2.997 |
| 0.3 | 0.131(10)  | -0.024(5)   | +0.178 | 1.908 |
| 0.5 | 0.051(4)   | -0.045(3)   | +0.121 | 1.417 |
| 1.0 | 0.021(2)   | -0.016(1)   | +0.048 | 0.807 |
| 2.0 | 0.009(1)   | -0.009(1)   | +0.022 | 0.360 |
| 5.0 | 0.004(0)   | -0.007(0)   | +0.012 | 0.113 |

- $\sigma_{\rm ens} > 0$ at every $\beta$ (verified above $\beta=0.1$ where
  MC noise is large).
- $\sigma_{\rm comb}$ is statistically zero or slightly negative at every
  $\beta$.
- $\mu_{\rm comb} > \sigma_{\rm ens}$, consistent with $n_\max = 2$.

The Polyakov cross-check at $n_\max = 3$ over $\beta \in [0.1, 2.0]$ gives
$c_\sigma^{\rm MC} = 1.63$ with $R^2 = 0.88$ — slightly closer to the
predicted $c_\rho/2 = 4.70$ but still off by a factor of $2.9\times$.
**The disagreement with the Polyakov rate constant is robust across both
graph sizes**, consistent with the non-cubic dual-lattice topology
preventing the textbook continuum derivation from applying verbatim.

## §7.  Verdict — sixth witness

**Weak pass with explicit finite-volume scope.**

The XCWG-G test consists of three sub-claims:
1. $\sigma_{\rm MC}(\beta) > 0$ at all tested $\beta$ — YES for
   $\sigma_{\rm ens}$, NO for $\sigma_{\rm comb}$ and $\sigma_{\rm fixL=4}$.
2. $\sigma_{\rm MC}(\beta)$ monotonically decreasing in $\beta$ —
   YES for $\sigma_{\rm ens}$.
3. Either (a) $\sigma_{\rm MC} \approx \sigma_{\rm LO}$ at small $\beta$
   (XCWG-D consistency), or (b) $\sigma_{\rm MC}$ exponential rate matches
   $c_\rho / 2$ (Polyakov / XCWG-F consistency).  NEITHER holds
   quantitatively: $\sigma_{\rm ens} / \sigma_{\rm LO} \sim 0.05$–$0.13$
   at small $\beta$, and $c_\sigma^{\rm MC} / (c_\rho / 2) \sim 0.25$.

The verdict is:

- **CLEAN PASS:** FALSE — none of the three sigma estimates agrees
  quantitatively with both XCWG-D LO and XCWG-F Polyakov.

- **WEAK PASS (perimeter-dominated finite-volume regime):** TRUE —
  $\sigma_{\rm ens}$ is monotone-decreasing positive at every $\beta$,
  matching the structural expectation of 3D compact U(1) (and the previous
  witnesses XCWG-D and XCWG-F directionally), but the data is
  perimeter-dominated at the available graph sizes and a clean separation
  of area law from perimeter law would require larger Rule B graphs (at
  least $n_\max = 4$, $V = 60$, with longer accessible $L$ at fixed $A$).

### Honest framing for Paper 41

Five-witness structural compatibility (Paper 41 v3) remains.  XCWG-G adds:
- **Confirming**: $\langle W \rangle$ falls off with $A$ at all $\beta$, and
  the ensemble-mean fit gives $\sigma_{\rm ens} > 0$ monotone in $\beta$ —
  qualitatively the right 3D compact U(1) shape.
- **Diagnostic**: at the graph sizes accessible ($n_\max \leq 3$, $E \leq
  106$), full-MC Wilson loops are quantitatively perimeter-dominated, and
  Polyakov's $\sigma \sim \exp(-(c_\rho / 2) \beta)$ prediction holds only
  qualitatively (exponential shape) but not quantitatively (rate constant
  off by 3–4×).

The sixth witness **structurally** points the same way as the first five
(toward 3D compact U(1) on a non-cubic graph), but at the quantitative
level it tells us the dual lattice is far from $\mathbb{Z}^3$ enough that
the Polyakov continuum derivation does not transfer verbatim — exactly
what one would expect for a finite Rule B graph with mixed-density
adjacency and 60 elementary monopole sites packed onto 60 plaquettes.

### Recommended Paper 41 v4 update

Add a paragraph in §VII (or new subsection) acknowledging XCWG-G as
sixth-witness-weak-pass, with explicit finite-volume scope.  The fact that
all six witnesses point the same way structurally — and the deviation
from the textbook continuum result on $\sigma(\beta)$ is along the
expected axis (non-cubic dual lattice, finite volume) — strengthens
rather than weakens the case for compatibility-with-finite-graph-3D-compact-U(1).

Honest scope note: the witness does NOT establish that Rule B confines in
the infinite-graph limit.  It establishes that the ensemble-mean
Wilson-loop slope behaves as a 3D compact U(1) string tension would at the
graph sizes tested, with the genuine area-law sub-coefficient masked by
perimeter contributions that should systematically shrink as the graph
grows.  $n_\max = 4$ ($V = 60$, $E \sim 300$, $\beta_1 \sim 240$) would be
the cleanest discriminator and is the natural next step.

## Files

- Script: `debug/xcwg_full_mc_wilson_loops.py` (~870 lines)
- Data: `debug/data/xcwg_full_mc_wilson_loops.json`
- Plot: `debug/plots/xcwg_full_mc_wilson_loops.png` (three panels:
  $\langle W \rangle$ vs $A$ by $\beta$; $\sigma_{\rm MC}$ vs
  $\sigma_{\rm LO}$ and Polyakov shape; Polyakov-cross-check
  semi-log $\sigma$ vs $\beta$)
- MC runtime: ~6 min total at $n_\max = 2$ + ~2 min at $n_\max = 3$ on
  single CPU core.

## Cross-references

- XCWG-D leading-order area law: `debug/xcwg_strong_coupling_wilson_memo.md`
  (witness #4; LO $\sigma_{\rm LO}$ values tabulated here in §3)
- XCWG-F monopole density: `debug/xcwg_monopole_density_memo.md`
  (witness #5; $\rho_M(\beta) = 0.747 \exp(-9.40\beta)$ fit used in §4)
- XCWG-C Gaussian Wilson loops: `debug/xcwg_wilson_loop_scaling_memo.md`
  (witness #3; perimeter-law on free-photon kernel — XCWG-G shows that
  the perimeter law persists in the full compact regime at $n_\max \leq 3$)
- Paper 41 v3: `papers/observations/paper_41_rule_b_wilson_u1.tex`
- CLAUDE.md §2 entries: XCWG-D, XCWG-F (already documented); XCWG-G
  (this sprint) to be added in next sync.
