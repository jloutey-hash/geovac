# XCWG Polyakov-rate refinement on Rule B compact U(1) — resolution by structural reading

**Sprint:** XCWG-H (post-XCWG-G follow-up on Paper 41 v4 caveat (vi))
**Date:** 2026-05-16
**Files:**
- `debug/xcwg_polyakov_rate_refinement.py` (driver, ~580 lines)
- `debug/data/xcwg_polyakov_rate_refinement.json` (tables)

## Headline

**RESOLVED — but by a structural reading, not by a rate-constant correction to the
dilute-gas formula.** The 3–4× quantitative disagreement between XCWG-F's
$c_\rho = 9.40$ (monopole-density rate) and XCWG-G's $c_\sigma^{\rm MC} = 1.20$
(measured Wilson-loop "string tension" rate) at $n_{\max}=2$ is a **category
mismatch between observables, not a missing structural correction to the
dilute-gas formula**. The "measured $c_\sigma$" extracted by XCWG-G is
statistically identical to the **perimeter-coefficient rate** $c_{\mu_{\rm comb}}$
on the same MC data ($c_{\sigma_{\rm ens}} = 1.196$ vs $c_{\mu_{\rm comb}} = 1.242$;
ratio $\sigma_{\rm ens} / \mu_{\rm comb} \to 1$ within MC noise across the fit
window). The TRUE area-law coefficient $\sigma_{\rm comb}$ — extracted from the
joint area-plus-perimeter fit on the same data — is **statistically zero at every
$\beta$** ($|\sigma_{\rm comb}| < 0.011$ Ha at every $\beta$, mean $-0.0036$),
which is **exactly what Polyakov's $\sigma \sim \sqrt{\rho_M / \beta_V}$ predicts**
given that the predicted $\sigma_{\rm Polyakov}(\beta = 1.5) \sim 5 \times 10^{-9}$
sits four orders below the XCWG-G MC detection floor.

The honest reading therefore is: **the Polyakov dilute-gas formula on Rule B is
not in conflict with the data; it just predicts an unmeasurably small area
coefficient at the graph sizes tested, leaving the dominant Wilson-loop falloff
signal in the perimeter sector**, which is governed by an unrelated mechanism
(strong-coupling link self-energy, $c_\mu \sim 1$ at small $\beta$).

That said, the non-cubic dual-lattice topology on Rule B **does** produce real
structural corrections relative to $\mathbb{Z}^3$: the dual graph $G_B^*$ has
average degree $9.60$ vs $\mathbb{Z}^3$'s $6.0$ (60% denser), and the dual
Laplacian spectral gap is $\lambda_2 = 2.138$ vs $\mathbb{Z}^3$'s $3.0$ (29%
smaller). In the refined dilute-gas formula these factors enter the
**prefactor** $\sigma_0$, not the rate constant $c_\sigma$ — which remains
$c_\rho/2$ in any sine-Gordon dualization. The finite-volume floor effect
predicts $\sigma \to {\rm const} \sim 0.66$ at large $\beta$, but this is
irrelevant for the rate-constant comparison since $\sigma_{\rm comb}$ stays
zero anyway.

Paper 41 v4 caveat (vi) can be reframed accordingly: the disagreement is not a
gap in the framework; it is a difference between the **ensemble area-law slope**
(which is perimeter-contaminated and bears no relation to Polyakov) and the
**Polyakov-predicted area coefficient** (which is consistent with measurement
modulo the MC floor).

---

## §1. Standard Polyakov derivation: recap

Polyakov 1977 (NPB 120, 429) derives confinement in 3D compact U(1) lattice
gauge theory via the **sine-Gordon dualization** of monopole-instanton
configurations. The compact gauge field $\theta_e$ on links of a 3D cubic lattice
$\mathbb{Z}^3$ admits a dual scalar field $\chi_x$ on dual lattice sites
$x \in \mathbb{Z}^{3*}$ (which on $\mathbb{Z}^3$ are the 3-cells = unit cubes).
Standard derivation:

$$
S_{\rm dual}[\chi] \;=\; \frac{1}{2\beta_V}\sum_{\langle xy\rangle}(\chi_x - \chi_y)^2 \;-\; 2z\sum_x \cos(2\pi \chi_x) ,
$$

where $z = \tfrac{1}{2}\rho_M^{(\rm LO)}(\beta)$ is the **monopole fugacity**.
Expanding the cosine to quadratic order yields an effective scalar mass

$$
m_D^2 \;=\; (2\pi)^2 \cdot 2z / \beta_V \;=\; \frac{4\pi^2 \rho_M}{\beta_V},
$$

the **Debye screening mass**, which absorbs to the Polyakov 1977 form
$m_D^2 = \rho_M / \beta_V$ under a $4\pi^2$ rescaling of $\rho_M$. Either
convention gives the same physical content. The string tension between two test
charges at large separation comes from the Coulomb-gas propagator of $\chi$ at
finite Debye screening:

$$
\boxed{\;\sigma(\beta) \;=\; \frac{4}{\pi}\sqrt{\rho_M / \beta_V}\;}
$$

(Polyakov 1977 eq. 4.12 in conventions where $\beta_V$ is the Villain coupling).
The $4/\pi$ prefactor is the **M1 Hopf-base measure** $\text{Vol}(S^2)/\pi^2$
signature of the Sprint TS-E1 master Mellin engine (Paper 32 §VIII). Inserting
the dilute-gas form $\rho_M = A\,e^{-c_\rho \beta}$ gives

$$
\log \sigma(\beta) \;=\; \log\bigl((4/\pi)\sqrt{A/\beta_V}\bigr) \;-\; (c_\rho/2)\,\beta - \tfrac{1}{2}\log\beta_V .
$$

The **rate constant** for $\log\sigma$ vs $\beta$ is therefore

$$
c_\sigma^{\rm Polyakov} \;=\; c_\rho / 2 .
$$

Standard derivation. No deviation possible without modifying the underlying
sine-Gordon mechanism. For XCWG-F's $c_\rho = 9.40$ at $n_{\max}=2$, this gives
$c_\sigma^{\rm predicted} = 4.70$. XCWG-G measured $c_\sigma^{\rm MC} = 1.20$ —
**a 3.9× disagreement on the rate constant**, the subject of Paper 41 v4 caveat
(vi).

The question we test: does Rule B's non-cubic dual-lattice topology give a
structural correction to this formula that brings the predicted rate constant
down to $\approx 1.20$?

## §2. Dual graph $G_B^*$ of Rule B

### §2.1 Construction

On a cubic lattice, dual vertices are the elementary closed 2-surfaces (the
3-cubes, each bounded by 6 plaquettes), and dual edges connect 3-cubes that
share a 2-face (i.e., a plaquette). This is the **plaquette-sharing dual**.

On Rule B the analog is unambiguous (XCWG-F established this):

- **Dual vertices of $G_B^*$**: the elementary closed 2-cycles, which are the
  size-$k_\delta = 3$ triangle-prism configurations of three $L=4$ plaquettes
  sharing common edge axes. There are 60 such "monopole sites" at $n_{\max} = 2$
  and 200 at $n_{\max} = 3$ (the latter capped by the
  `enumerate_closed_2cycles_up_to_size` `max_count = 200` budget).

- **Dual edges**: two monopole sites are dual-adjacent if and only if they share
  at least one primal plaquette. We tabulate both the **unweighted** (binary
  shared-or-not) and **weighted** (count of shared plaquettes) adjacencies; on
  Rule B every two-site pair sharing any plaquette shares exactly one, so the
  two coincide.

### §2.2 Tabulation at $n_{\max} \in \{2, 3\}$ and comparison to $\mathbb{Z}^3$

| Quantity | Rule B $n_{\max}=2$ | Rule B $n_{\max}=3$ | $\mathbb{Z}^3_3$ (3×3×3 with PBC) |
|:---|:--:|:--:|:--:|
| primal $V, E, P$ | 10 / 20 / 44 | 28 / 106 / 500 (capped) | 27 / 81 / 81 |
| monopole sites $\lvert S\rvert$ | 60 | 200 (capped) | 27 |
| dual graph $V^*, E^*$ | 60 / 288 | 200 / 1213 | 27 / 81 |
| dual avg degree | **9.60** | **12.13** | **6.00** |
| dual min/max degree | 8 / 12 | 3 / 27 | 6 / 6 (regular) |
| dual connected components | 1 | 1 | 1 |

The $\mathbb{Z}^3$ entry serves as the standard Polyakov baseline. On $\mathbb{Z}^3$
each dual vertex (3-cube) is adjacent to exactly six others (sharing one of its
six faces), so the dual is the 6-regular graph $\mathbb{Z}^3$ with $\lvert S\rvert$
sites and $\lvert E^*\rvert = 3\lvert S\rvert$ edges (each shared face counted
once). At $\lvert S\rvert = 27$ the $3\cdot 27 = 81 = E^*$ entry confirms the
construction is correct, and the cube boundary closure $\sum_{P\in S}\text{sign}(P,S)\,
d_1[P,:] = 0$ for every cube is verified to machine precision.

**Structural reading.** Rule B's dual graph has 60% higher average connectivity
than $\mathbb{Z}^3$'s, and at $n_{\max}=3$ the dual graph is heterogeneous
(degree spread 3 to 27, no regular structure). This is the "non-cubic" topology
the caveat (vi) names.

## §3. Dual Laplacian spectrum

We compute the eigenvalues of the unweighted dual Laplacian $L_{\rm dual} =
D_{\rm dual} - A_{\rm dual}$ on $G_B^*$ and on $\mathbb{Z}^3_3$. Key statistics:

| Spectrum statistic | Rule B $n_{\max}=2$ | Rule B $n_{\max}=3$ | $\mathbb{Z}^3_3$ |
|:---|:--:|:--:|:--:|
| spectral gap $\lambda_2$ | **2.1379** | 1.4563 | **3.0000** |
| max eigenvalue $\lambda_{\max}$ | 14.2749 | 28.8709 | 9.0000 |
| mean eigenvalue $\bar\lambda$ | 9.6000 | 12.1300 | 6.0000 |
| zero eigvals (= components) | 1 | 1 | 1 |

The spectral gap quantifies the **smallest non-trivial eigenmode** of the dual
graph, i.e., the slowest decay mode of correlations along the dual structure.

- **Rule B at $n_{\max}=2$**: $\lambda_2 = 2.138$, **29% smaller** than $\mathbb{Z}^3_3$'s
  $\lambda_2 = 3.0$. Smaller gap means the dual graph has a "longest-wavelength
  fluctuation mode" that decays more slowly — the analog of a longer "effective
  Debye length" set by the finite graph geometry.

- **Rule B at $n_{\max}=3$**: $\lambda_2 = 1.456$, even smaller (the graph is
  larger so longest wavelengths fit). This is the expected $\lambda_2 \to 0$
  trend in the thermodynamic limit.

The full eigenvalue sequence at $n_{\max}=2$ (first 10 sorted):

`[0.0, 2.138, 5.226, 5.226, 5.561, 7.027, 7.027, 7.439, 7.439, 8.156]`

— a clean spectrum showing one zero (component count) and a positive gap, no
spurious near-zero modes.

For comparison, $\mathbb{Z}^3_3$'s spectrum is the discrete cosine transform of
the cubic lattice: eigenvalues $\lambda_{(k_1,k_2,k_3)} = \sum_{i=1}^3 2(1-\cos(2\pi
k_i/3))$ for $k_i \in \{0, 1, 2\}$, with the minimum nonzero value $2(1 - \cos(2\pi/3))
= 3$. We verified this to machine precision in the script.

## §4. Refined dilute-gas formula on a non-cubic dual graph

### §4.1 Where the dual Laplacian enters

The sine-Gordon kinetic term $\tfrac{1}{2\beta_V}\sum_{\langle xy\rangle}(\chi_x -
\chi_y)^2$ on the dual lattice is, in matrix notation, $\tfrac{1}{2\beta_V}
\chi^\top L_{\rm dual}\chi$. After Fourier transform on the discrete graph
($\chi = \sum_k a_k \psi_k$ with $L_{\rm dual}\psi_k = \lambda_k\psi_k$), the
free propagator becomes

$$
\langle a_k a_{k'}\rangle = \beta_V\frac{\delta_{kk'}}{\lambda_k + m_D^2},
$$

with $m_D^2 = 4\pi^2 \rho_M / \beta_V$ as before. The Wilson-loop area law then
arises (Polyakov 1977 §4) from the dual-field propagator evaluated between two
test-charge insertions separated by transverse distance $L$:

$$
\sigma \sim \frac{1}{\pi} \sqrt{\frac{2 m_D^2}{\beta_V}}.
$$

The same formula on a general graph holds for the *continuum* (zero-momentum)
component of the propagator, i.e., the contribution from modes with $\lambda_k$
much smaller than $m_D^2$. The **rate constant** $c_\sigma = c_\rho/2$ is thus
**topology-independent** at the level of the standard sine-Gordon derivation:
the dual graph's spectrum only enters the prefactor (which mode contributes
most), not the rate.

### §4.2 Finite-volume floor

On a finite graph $\lambda_2 > 0$, which means the dual Laplacian has a finite
spectral gap. The propagator at $k$ corresponding to $\lambda_2$ is
$1/(\lambda_2 + m_D^2)$, and the analog of $\sigma$ saturates at

$$
\sigma_{\rm floor}^{(\lambda_2)} \;\sim\; \frac{1}{\pi}\sqrt{\frac{2\lambda_2}{\beta_V}}.
$$

For Rule B at $n_{\max}=2$ with $\lambda_2 = 2.138$ this gives $\sigma_{\rm floor}
\approx 0.658$. The crossover between the Polyakov regime ($m_D^2 \gg \lambda_2$)
and the saturated-floor regime ($m_D^2 \ll \lambda_2$) happens at $4\pi^2 \rho_M
\approx \lambda_2$, which using the XCWG-F fit $\rho_M = 0.747\,e^{-9.40\beta}$
gives **$\beta_{\rm cross} \approx 0.279$**.

If the **measured slope** of $\log\sigma$ vs $\beta$ averaged the two regimes
linearly, we would predict a slope intermediate between $0$ (saturated floor)
and $c_\rho/2 = 4.70$ (pure Polyakov). Fitting the model $\sigma_{\rm refined}^2
= \sigma_{\rm floor}^2 + \sigma_{\rm Polyakov}^2$ (Pythagorean soft floor) over
$\beta \in [0.1, 3.0]$ at $n_{\max}=2$ gives slope $-0.18$, **way smaller** than
the measured $-1.20$.

Fitting the alternative model $\sigma_{\rm refined}^2 = (1/\pi^2)(2/\beta_V)
(\lambda_2 + m_D^2)$ (quadratic-in-mass softer floor) over $\beta \in [0.1, 3.0]$
gives slope $-0.21$, **also too small**.

**Conclusion.** Neither floor-based refinement of the dilute-gas formula matches
the measured slope. The refined formula slopes are *too flat*; the measured
slope $1.20$ sits between the standard Polyakov $4.70$ and the floor-saturated
$\sim 0.18$, suggesting that **neither model describes what is actually being
measured**.

## §5. What XCWG-G actually measured: $c_{\sigma_{\rm ens}} = c_{\mu_{\rm comb}}$

This is the load-bearing structural finding. XCWG-G's full-MC sigma analysis
fits three different observables to the per-loop Wilson data:

- **$\sigma_{\rm ens}$**: ensemble area-law slope $\log\langle W(A)\rangle$ vs
  area $A$, using per-area ensemble means. This is the "naive XCWG-D-style"
  fit and is **perimeter-contaminated** because at fixed area class the
  available loops have varying perimeters.

- **$\sigma_{\rm comb}$, $\mu_{\rm comb}$**: joint area-plus-perimeter fit
  $\log W_i = -\sigma A_i - \mu L_i + c$ on per-loop data. $\sigma_{\rm comb}$
  is the **clean area coefficient**; $\mu_{\rm comb}$ is the perimeter
  coefficient.

- **$\sigma_{\rm fixL}$**: pure area fit restricted to loops at a single fixed
  perimeter — cleanest pure-area test.

XCWG-G's headline numbers at $n_{\max}=2$ from `xcwg_full_mc_wilson_loops.json`:

| $\beta$ | $\sigma_{\rm ens}$ | $\sigma_{\rm comb}$ | $\mu_{\rm comb}$ | $\sigma_{\rm LO}$ |
|:--:|:--:|:--:|:--:|:--:|
| 0.10 | 0.377 | **-0.007** | 0.342 | 2.997 |
| 0.30 | 0.174 | **-0.006** | 0.217 | 1.908 |
| 0.50 | 0.089 | **-0.003** | 0.097 | 1.417 |
| 1.00 | 0.033 | **-0.002** | 0.034 | 0.807 |
| 2.00 | 0.016 | **-0.001** | 0.016 | 0.360 |
| 3.00 | 0.011 | **-0.000** | 0.011 | 0.211 |

The XCWG-G v1 memo explicitly states: *"$\sigma_{\rm ens} \approx \mu_{\rm comb}$
at every $\beta$"*, identifying this as a direct empirical signature that the
apparent area-law signal is carried by perimeter contributions on this graph.

### §5.1 Numerical verification

We compute the rate constants (slope of $\log(\cdot)$ vs $\beta$) for all four
quantities on the XCWG-G data at $n_{\max}=2$, fit over $\beta \in [0.1, 3.0]$:

| Quantity | Rate $c$ | $R^2$ | XCWG-G v1 reported |
|:---|:--:|:--:|:--:|
| $c_{\sigma_{\rm ens}}$ | **1.196** | 0.831 | 1.20 |
| $c_{\mu_{\rm comb}}$ | **1.242** | 0.838 | (not reported as rate) |
| $c_{\sigma_{\rm LO}}$ | 0.900 | 0.954 | (not a rate; algebraic decay) |
| $c_{\sigma_{\rm comb}}$ | (undefined, $\sigma_{\rm comb}$ is statistically zero) | — | "statistically zero" |

**$c_{\sigma_{\rm ens}} = 1.196$ and $c_{\mu_{\rm comb}} = 1.242$ agree to within
4%**. The point-by-point ratio is $\sigma_{\rm ens}/\mu_{\rm comb} = 1.034 \pm
0.115$ across the fit window (CV $\approx 11\%$).

**This is the structural fact**: the "measured $c_\sigma = 1.20$" of Paper 41 v4
caveat (vi) is empirically the **perimeter coefficient decay rate**, and only
nominally the "string-tension rate" — the genuine area coefficient is
statistically zero in this MC.

### §5.2 What does $c_{\mu_{\rm comb}} \approx 1.2$ actually correspond to?

$\mu_{\rm comb}$ is the **self-energy per unit perimeter** of a Wilson loop in
the full compact theory. The natural strong-coupling expansion gives (XCWG-D)

$$
\langle W(C)\rangle \;=\; \bigl(I_1(\beta)/I_0(\beta)\bigr)^{A_{\min}} \cdot \bigl(I_0(\beta)\bigr)^{-L} \cdot (\text{NLO corrections})
$$

at the level of the tessellation factor, where $I_0(\beta) \to e^\beta/\sqrt{2\pi\beta}$
at large $\beta$. The strong-coupling perimeter coefficient is then $\mu_{\rm comb}
\sim -\log(I_0)^{-1} \cdot (\text{some O}(1)\text{ factor})$. We can also
write $\mu_{\rm comb} \sim c \cdot \log(I_1/I_0) \cdot (\text{some small
coefficient})$ from the leading-tadpole resummation.

We measure $c_{\sigma_{\rm LO}} = 0.900$ where $\sigma_{\rm LO} = -\log(I_1(\beta)/
I_0(\beta))$ — this is the **strong-coupling area-law character-expansion rate**
from XCWG-D (which is algebraic in $\beta$ at large $\beta$, $\sigma_{\rm LO} \to
1/(2\beta)$, not exponential, but the fit over $[0.1, 3.0]$ extracts an apparent
exponential rate from the curve segment). The measured $c_{\mu_{\rm comb}} = 1.242$
is **within a factor of 1.4 of $c_{\sigma_{\rm LO}}$** — both numbers are of order 1,
controlled by the same character-expansion structure $I_n(\beta)/I_0(\beta)$.

In contrast, $c_\rho/2 = 4.7$ comes from the **monopole instanton density** which
decays much faster ($\rho_M \propto e^{-2\pi^2 v_0 \beta}$ with $v_0 \approx 0.47$).
The monopole sector and the perimeter self-energy sector are **completely
separate** dynamical mechanisms.

## §6. Sanity check on $\mathbb{Z}^3$

On a 3×3×3 cubic lattice with periodic boundary conditions:

- 27 vertices, 81 edges, 81 plaquettes, **27 monopole sites** (one per unit cube),
  cube boundary closure verified to machine precision.
- Dual graph: 27 vertices, 81 edges, **avg degree exactly 6** (regular).
- Dual Laplacian spectrum (analytic from DFT): smallest nonzero $\lambda_2 = 3.0$,
  max $\lambda_{\max} = 9.0$, mean $\bar\lambda = 6.0$ — matches the computation
  to machine precision.

Thus the refined-formula machinery reduces correctly to the standard
Polyakov-on-cubic structure when applied to $\mathbb{Z}^3$. The non-trivial
structural corrections (avg degree, spectral gap) are visible specifically
because Rule B is non-cubic, but they affect the **prefactor** of $\sigma$, not
the rate.

## §7. Verdict

### §7.1 Three readings, three verdicts

The Polyakov rate-constant disagreement of Paper 41 v4 caveat (vi) is RESOLVED
under one of three readings, depending on what observable is meant by "$\sigma$":

**(A) Resolution by structural reading (favored).** The XCWG-G "measured
$c_\sigma = 1.20$" is the rate of $\sigma_{\rm ens}$, which equals the rate of
$\mu_{\rm comb}$ to within 4% — these are the same observable up to MC noise.
$\mu_{\rm comb}$ is the **perimeter self-energy coefficient** in the strong-coupling
character expansion, related to $I_n(\beta)/I_0(\beta)$ at the link level; the
$O(1)$ rate $c \approx 1$ is the natural scale of this expansion. The Polyakov
prediction $\sigma \sim \sqrt{\rho_M / \beta_V}$ refers to the **area-law string
tension** $\sigma_{\rm comb}$, which on the XCWG-G data is **statistically zero
at every $\beta$**. So Polyakov's formula is in fact consistent with the
measurement — both predict $\sigma_{\rm area} \approx 0$ at the graph sizes
tested. The 3.9× rate disagreement is a category mismatch between two different
observables ($\sigma_{\rm ens} \neq \sigma_{\rm comb}$), not a structural failure.

**(B) "Refined dilute-gas formula" verdict (negative).** If we insist on
comparing the Polyakov formula to the XCWG-G $\sigma_{\rm ens}$ observable
(which is what Paper 41 v4 caveat (vi) literally states), the standard formula
predicts $c_\sigma^{\rm Polyakov} = c_\rho/2 = 4.7$ — disagreeing by 3.9×.
Including a finite-volume floor from the dual Laplacian spectral gap $\lambda_2
= 2.14$ predicts a slope $\approx -0.18$ to $-0.21$ depending on the soft-floor
convention — disagreeing in the OTHER direction by 5–7×. Neither structural
correction recovers the measured rate. **By this reading, the disagreement is
NOT resolved by non-cubic dual-lattice topology.**

**(C) Hybrid honest reading.** The non-cubic dual-lattice topology DOES produce
real structural corrections to the Polyakov dilute-gas analysis:

- Average dual degree $9.6$ vs $\mathbb{Z}^3$'s $6.0$ — 60% denser. Enters
  $\sigma$ prefactor as $\sqrt{(\bar\lambda/6)} \approx 1.27$.
- Dual spectral gap $\lambda_2 = 2.14$ vs $\mathbb{Z}^3$'s $3.0$ — 29% smaller.
  Sets a finite-volume floor for $\sigma$ at $\sigma_{\rm floor} \approx 0.66$.
- Floor crossover at $\beta_{\rm cross} \approx 0.28$, exactly where rho_M starts
  to drop below the saturation scale.

But these corrections **do not change the rate constant** $c_\sigma = c_\rho/2$
in the standard sine-Gordon derivation; they only modify the **prefactor**
$\sigma_0$ and the **finite-volume floor** $\sigma_{\rm floor}$. The genuine
3.9× rate disagreement of caveat (vi) — if it were a real disagreement between
the same observable measured one way and predicted another — would not be
resolved by these topological corrections.

The honest synthesis is reading (A) + reading (C): the rate-constant
"disagreement" is an artifact of comparing different observables; the
non-cubic dual-lattice corrections ARE real but they would only matter if the
Polyakov prediction were being applied to the right observable in the right
regime (which on the available graph sizes it isn't, because $\sigma_{\rm
Polyakov}(\beta = 1.5) \approx 5\times 10^{-9}$ is far below the MC detection
floor).

### §7.2 Predicted $c_\sigma$ and verdict statement

| Quantity | $n_{\max} = 2$ | $n_{\max} = 3$ |
|:---|:--:|:--:|
| $c_{\sigma}^{\rm MC}$ (XCWG-G measured = $c_{\sigma_{\rm ens}}$) | 1.196 | 1.633 |
| $c_{\mu_{\rm comb}}$ (perimeter rate from same MC) | **1.242** | **1.346** |
| $c_{\sigma_{\rm comb}}$ (genuine area-law rate) | $\sigma_{\rm comb} \approx 0$ | $\sigma_{\rm comb} \approx 0$ |
| $c_{\sigma}^{\rm Polyakov}$ (standard, $c_\rho/2$) | 4.701 | 4.476 |
| $c_{\sigma}^{\rm refined-floor}$ (Pythagorean) | 0.177 | 0.361 |
| $|c_\rho/2 - c_{\sigma}^{\rm MC}|$ | 3.50 (3.9×) | 2.85 (2.7×) |
| $|c_\mu - c_{\sigma}^{\rm MC}|$ | **0.046 (4%)** | **0.287 (18%)** |

**At $n_{\max}=2$: the disagreement between $c_\sigma^{\rm MC}$ and $c_{\mu_{\rm
comb}}$ is 4% — the same observable, not a Polyakov-prediction failure.**

At $n_{\max}=3$ the agreement is only 18% (the MC error is larger because
$n_{\rm sample} = 600$ vs $2000$ at $n_{\max}=2$), but the same structural
ordering $c_{\sigma_{\rm ens}} > c_{\mu_{\rm comb}}$ persists.

The dual-lattice topology effects on the prefactor — $\sigma_{\rm Polyakov}^0
\sim (\bar\lambda/6)^{1/2} \times (4/\pi)\sqrt{A/\beta_V}$ — would modify the
PREDICTED $\sigma_{\rm comb}$ values by a factor $\sim 1.27$ on Rule B vs Z^3,
which is well within the MC noise of the genuine $\sigma_{\rm comb}$
measurement.

### §7.3 One-sentence verdict

**The 3-4× quantitative disagreement of Paper 41 v4 caveat (vi) is RESOLVED
under the structural reading that XCWG-G's "measured Polyakov rate" $c_\sigma =
1.20$ is in fact the perimeter coefficient $c_{\mu_{\rm comb}}$ — matching the
true area-law $c_{\sigma_{\rm comb}}$ ($\approx 0$ within MC noise) and the
Polyakov formula's prediction simultaneously — while the non-cubic dual-lattice
topology DOES produce real structural corrections (60% higher avg dual degree,
29% smaller spectral gap $\lambda_2$) that enter the PREFACTOR of $\sigma$, not
the rate constant.**

## §8. Recommendation for Paper 41 v4 caveat (vi)

Replace the current text of caveat (vi):

> *"Polyakov rate constant quantitative disagreement (3-4×) between XCWG-F
> prediction and XCWG-G direct measurement — reflects non-cubic dual-lattice
> topology."*

with:

> *"The XCWG-G ensemble area-law slope $c_{\sigma_{\rm ens}} = 1.20$ does not
> compare to the Polyakov-XCWG-F prediction $c_\rho/2 = 4.70$: it is
> structurally identical to the perimeter coefficient rate
> $c_{\mu_{\rm comb}}$ on the same MC data (agreement to 4% at $n_{\max}=2$),
> reflecting strong-coupling perimeter self-energy ($O(1)$ rate from
> $I_n(\beta)/I_0(\beta)$ character expansion) rather than the area-law string
> tension. The true area-law coefficient $\sigma_{\rm comb}$ is statistically
> zero at every $\beta$ tested, which is consistent with the Polyakov formula
> given that $\sigma_{\rm Polyakov}(\beta \geq 1.5) \sim 5 \times 10^{-9}$ is
> four orders below MC detection floor. The non-cubic dual-lattice topology of
> Rule B does produce real structural corrections (avg dual degree 9.6 vs $\mathbb{Z}^3$'s
> 6.0; dual spectral gap $\lambda_2 = 2.14$ vs $\mathbb{Z}^3$'s 3.0) but these
> affect the prefactor of $\sigma$, not the rate constant in the standard
> sine-Gordon derivation."*

This converts a "honest-scope caveat" into a "structural finding plus honest
scope on prefactor".

---

## §9. Reproducibility

```bash
cd Project_Geometric
python debug/xcwg_polyakov_rate_refinement.py
```

Outputs:
- `debug/data/xcwg_polyakov_rate_refinement.json` — full numerical tables
- This memo

Runtime: ~20 s on single CPU (dominated by $n_{\max}=3$ plaquette enumeration
at the size-3 closed-cycle level).

**Bibliography note (in-memo).** Polyakov 1977 NPB 120 429; Banks-Myerson-Kogut
1977 NPB 129 493; Guth 1980 PRD 21 2291; DeGrand-Toussaint 1980 PRD 22 2478;
Paper 41 v4 §VII (`papers/observations/paper_41_rule_b_wilson_u1.tex`); XCWG-F
memo `debug/xcwg_monopole_density_memo.md`; XCWG-G memo
`debug/xcwg_full_mc_wilson_loops_memo.md`; Paper 32 §VIII (master Mellin engine,
$4/\pi$ as M1 Hopf-base measure $\text{Vol}(S^2)/\pi^2$).
