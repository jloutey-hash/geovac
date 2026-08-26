# Sprint memo -- the T2 fibre as a 2D Euclidean two-propagator correlator (kinematics)

Date: 2026-08-21 | Branch: work/sparsity-boundary (Paper 59 UNCERTIFIED, in flight)
Owning paper: `papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex`
  (`eq:K0`, `eq:besselmoment`, `eq:laplace`, `sec:genus`, `sec:modular`)
Drivers (new): `debug/t2_euclidean_dictionary.py` (sections A,B,D,E,F),
  `debug/_t2_dict_secC.py` (C0-C3, fixed-GL), `debug/_t2_dict_secC3.py` (C3 light grid)
Data: `debug/data/_t2_dict_sec{C,C3,C3b,D}.out`, `_t2_dict_ABEF.out`
Predecessors (built on, NOT re-derived): `sprint_T2_zhou_bd_adaptation_memo.md`
  (closed-form Stokes triple), `sprint_routeC_irregular_resurgent_memo.md` (N(D) resurgence),
  `sprint_T2_eichler_lambert_memo.md` (co-area reduction).

**NO PAPER WAS MODIFIED.** Section 5 holds a *proposed* paragraph for PI decision.

---

## 0. Task (PI standing question: "what can we say about kinematics?")

Make the 2D-Euclidean-QFT reading of the Paper 59 fibre *exact* (prefactors, measures,
what every symbol is physically), validate it numerically, then re-read the closed-form
Stokes data and the corner obstruction kinematically, and test honestly whether this
sharpens Paper 35.

---

## 1. Step (i) -- the dictionary, made exact

### 1.1 `eq:K0` *is* the two-dimensional Euclidean propagator [MEASURED, 1.3e-30]

The Green's function of $-\Delta + a^2$ on $\mathbb{R}^2$ is $G_a(\tau,x) =
\frac{1}{2\pi}K_0\!\big(a\sqrt{\tau^2+x^2}\big)$. Paper 59's `eq:K0` reads

```
int_0^inf cos(k b) e^{-D sqrt(c k^2 + m^2)} / sqrt(c k^2 + m^2) dk
      = (1/sqrt c) K_0( a sqrt(p^2 + b^2) ),      a = m/sqrt(c),  p = D sqrt(c)
      = (2 pi / sqrt c) G_a(p, b)                 [ = 2 pi a G_a(p,b) at m = 1 ]
```

so the chemistry object is the 2D propagator with **mass** $a$, **Euclidean time** $p$ and
**space** $b$. Validated at five $(c,m,D,b)$ points, worst relative deviation
**1.28e-30** (section A). *Note the prefactor is $2\pi/\sqrt c$, which equals $2\pi a$ only
at $m=1$; the general-$m$ form is needed and was the one arithmetic slip caught by the gate.*

### 1.2 `eq:besselmoment` is a two-propagator correlator [MEASURED, 2.0e-31]

Parseval for cosine transforms turns the radial core into

```
M = int_0^inf dk  prod_i e^{-D_i Delta_i}/Delta_i
  = (2 / (pi sqrt(c1 c2))) int_0^inf db K_0(a1 sqrt(p1^2+b^2)) K_0(a2 sqrt(p2^2+b^2))
  = 4 pi a1 a2  int_R db  G_{a1}(p1, b) G_{a2}(p2, b) .
```

Validated at five $(c_1,c_2,D_1,D_2)$ points, worst **2.04e-31** (section B). The physical
reading: **two Euclidean packets of unequal mass, each a point source evolved for Euclidean
time $p_i$, overlapped on a common spatial slice.** The $\int db$ is the slice; the two
propagators do *not* share a time.

The $D\to0$ slice reproduces `eq:period` to **0.0** at four $(c_1,c_2)$ (section B) once the
prefactor is written symmetrically, $\int_0^\infty dk/\sqrt{(c_1k^2+1)(c_2k^2+1)} =
K(m)/\sqrt{c_{\max}}$, $m = 1-c_{\min}/c_{\max}$. **NIT for the PI:** as printed, `eq:period`'s
$1/(a_1\sqrt{c_1c_2})$ equals $1/\sqrt{c_{\min}}$, which is the correct value only under the
implicit convention $c_1 = c_{\min}$; the symmetric form $K(m)/\sqrt{c_{\max}}$ is
convention-free. Not an error, a reading hazard.

**Second NIT for the PI (internal inconsistency).** `sec:reduction` says "as $\rho\to0$ (coincident scales) $N\to K_0(D)$", but $\rho=c_2/c_1\to0$ is the limit in which the *second*
scale collapses (mass $a_2\to\infty$), i.e. the one-mass limit; the coincident-scale locus is
$\rho=1$, as `sec:modular` correctly states ("the cusp $\rho=1$ is the coincident-scale
diagonal"). Both $\rho=0$ and $\rho=1$ are genus-zero cusps of $\Gamma(2)$, but for different
reasons (branch points running off to infinity versus colliding), so the parenthetical in
`sec:reduction` should read "the one-mass limit". Two sections currently disagree.

### 1.3 What the $\partial_\zeta$'s are: a Euclidean source-time superposition [exact algebra]

The physical fibre factor is $P(x,k) = c\,e^{-\Delta}(\Delta^{-3}+3\Delta^{-4}+3\Delta^{-5})$,
$c = x(1-x)$, $\Delta=\sqrt{ck^2+1}$. Using $e^{-D\Delta}\Delta^{-n} =
\frac{1}{(n-1)!}\int_D^\infty (u-D)^{n-1}e^{-u\Delta}du$,

```
P(x,k)/c = int_D^inf w_D(u) e^{-u Delta} du ,   w_D(u) = (u-D)^2/2 + (u-D)^3/2 + (u-D)^4/8
```

(at the physical $D=1$: lower limit $u=1$). Since $\int_0^\infty\cos(kb)e^{-u\Delta}dk =
-\partial_p\big[K_0(a\sqrt{p^2+b^2})\big]_{p=u\sqrt c}$, the position-space packet is

```
f_i(b) = sqrt(c_i) int_{u >= D_i} w(u) u K_1(R)/R du ,   R = sqrt(u^2 + b^2/c_i)
       = c_i int  w(u) [ -d/dp  2 pi G_{a_i}(p,b) ]_{p = u sqrt(c_i)} du .
```

So a **Slater** density (as against a Yukawa) is a weighted superposition of point sources
emitted at Euclidean times $p \ge p_i = D_i\sqrt{c_i}$, and the *endpoint* $u=D_i$ of that
superposition is the physical time offset. Validated against the direct cosine transform at
five $(x,b)$ points, worst **5.35e-18** (C1).

### 1.4 The $j_0$ factor is a uniform box smearing of the shared space coordinate

$j_0(k|\bm W|) = \frac{1}{2|\bm W|}\int_{-|W|}^{|W|}\cos(kv)\,dv$ (exact, |d| = 0.0), hence

```
J(s,t) = int_0^inf dk j0(k|W|) P(s,k) P(t,k)
       = (2/pi) int_0^inf db  f_1(b) <f_2>_{|W|}(b),
   <f_2>_{|W|}(b) = (1/2|W|) int_{-|W|}^{|W|} f_2(|b+v|) dv .
```

Two validations. (C2) the **unsmeared** case ($|W|\to0$), both factors in propagator form:
$\int_0^\infty P(s,k)P(t,k)dk = (2/\pi)\int_0^\infty f_1f_2\,db$ at four $(s,t)$, worst
**2.1e-14**. (C3) the **full** fibre including the smearing, both factors in propagator form:
at $(s,t)=(0.3,0.2)$, $J_{\rm momentum}=0.208486308690628474$ vs
$J_{\rm position}=0.208486308727350167$, rel. dev. **1.8e-10**; at
$(s,t)=(0.4,0.4)$, $0.3052821202670761$ vs $0.305282120296606$, rel. dev. **9.7e-11**.
Both are nested-quadrature-limited, not identity-limited -- C1/C2 carry the precision.

### 1.5 The Feynman parameters are Kallen-Lehmann spectral masses [exact algebra]

$c = s(1-s)$ gives $a = \zeta/\sqrt{c} \ge 2\zeta$, with the minimum at $s=\tfrac12$. This is
exactly the two-propagator Feynman combination: **each two-center density is dispersed over
its own two-particle cut, starting at $2\zeta$**, and the "Fock scale" is the point on that
cut. Consequences:

* $\rho = c_2/c_1 = (a_1/a_2)^2$ is the **squared ratio of the two spectral masses**. The
  elliptic modulus is a mass-splitting datum, nothing else.
* The three $\Gamma(2)$ cusps are: $\rho=1$ equal masses (the coincident-scale diagonal),
  $\rho\to0,\infty$ one mass at the far (heavy) end of its cut.
* The diagonal cusp datum $A$ of `sec:modular` is the **double-threshold** value: verified
  $J(s{=}\tfrac12,t{=}\tfrac12) = 0.31621612846726750686 = 4A$ exactly (ratio 4.0 to 12
  digits) -- both spectral masses sitting at $2\zeta$, with $|\bm W| = s+t = 1$.

### 1.6 The exponential twist is mass x time, and it is family-invariant [exact algebra]

$a_i p_i = (\zeta/\sqrt{c_i})(D_i\sqrt{c_i}) = \zeta D_i$, independent of the Feynman
parameter. So as the parameter integral sweeps the whole family of masses, **the
dimensionless Euclidean-time displacement stays pinned at $\zeta D_i$**; the physical
$D=1$ (bohr, $\zeta=1$) twist is *one Compton wavelength of Euclidean-time displacement*.
$D\to0$ is coincident times = the pure period = the regular $\Gamma(2)$ MMV shadow.

### 1.7 The collinear three-center geometry, mapped

| chemistry | kinematics |
|---|---|
| bond lengths $D_1=\lvert X-Y\rvert$, $D_2=\lvert X-Z\rvert$ | the two **Euclidean-time offsets** $p_i = D_i\sqrt{c_i}$ |
| Slater exponent $\zeta$ | sets the cut threshold $2\zeta$ and the twist $\zeta D$ |
| Feynman parameters $s,t$ | positions on the two spectral cuts (the masses) |
| composite vertex $\bm W(s,t) = (t-s)X + sY - tZ$ | the **spatial** smearing range $\lvert\bm W\rvert$ |
| reference geometry $X{=}0,\,Y{=}(0,0,1),\,Z{=}(0,0,-1)$ | $p_1=p_2=\sqrt c$, $\lvert\bm W\rvert = s+t \in[0,2]$ |

The sharp point: **the third center enters only through the space coordinate.** The time
offsets are per-density and geometry-fixed; the third center changes only $\bm W$.

**Flagged for the PI (a care point, not an error).** `sec:genus` grades the transcendence by
center count ("genus zero at two centers, genus one at three"), but the curve
$y^2=(c_1k^2+1)(c_2k^2+1)$ knows only about the two Feynman scales, and the two-center
*exchange* class $(AB\vert AB)$ also carries two independent Feynman parameters -- yet it
closes at weight one. The kinematic discriminator is $\bm W$, not the curve: for
$(AB\vert AB)$, $\bm W = (t-s)(X-Y)$ is *slaved to the same bond vector* and vanishes on the
degeneration locus $s=t$; for $(XY\vert XZ)$, $\bm W$ is generically nonzero at $s=t$ and
carries the independent angle at $X$. So the genus-one statement is a statement about the
**fibre**, and the two-center family integral degenerates. One clarifying clause in
`sec:genus` would remove a real reader trap. (Not tested here; flagged only.)

---

## 2. Step (ii) -- the Stokes data is a complexified light cone

### 2.1 Derivation (exact, from 1.3-1.4)

$f_i(b)$ is analytic in $b^2$; its branch points are where $R = a_i\sqrt{p^2+b^2}$ vanishes at
the **endpoint** $u=D_i$ of the source-time superposition, i.e. where

```
b^2 + p_i^2 = 0          (the Euclidean interval of propagator i closes)
```

-- in Euclidean signature reachable only at complex $b = \pm i p_i$: the **complexified light
cone**. Under the box smearing, $\langle f_2\rangle(\beta) = \frac{1}{2|W|}\int_{-|W|}^{|W|}
f(\beta+v)dv$ has branch points in $v$ at $v = -\beta \pm i p_1$; an *interior* crossing is
contour-deformable, so a singularity of the smeared function occurs only at the **endpoint
pinch** $v=\pm|W|$:

```
beta* = +- |W| +- i p_1  ,        z* = (beta*)^2 ,        |z*| = p_1^2 + |W|^2 .
```

This is exactly the Stokes location of `sec:modular`, $z^*=-(\sqrt{c_1}\mp ib)^2$, once one
notes $p_1 = D_1\sqrt{c_1} = \sqrt{c_1}$ **at the physical $D=1$**. So:

> $|z^*|$ is the squared Euclidean interval between the source and the far end of the
> smearing range, and $\sqrt{z^*} = |W| + i p_1$ is the light-cone coordinate $x + i\tau$.
> The trans-series action $a_2\sqrt{z^*}$ is $m_2 \times$ (complexified interval) -- the
> exchange of the heavy packet across the light cone of the light one.

### 2.2 A refinement the derivation forces: the general-$D$ form

The memo/paper closed form is the $D=1$ slice. The derivation gives the general twist:

```
z*(D) = -( D sqrt(c_1)  -+ i |W| )^2  =  -( p_1 -+ i |W| )^2 ,     |z*| = p_1^2 + |W|^2 .
```

i.e. **literally the squared Euclidean interval**, for any twist. This is a genuinely new,
sharply testable statement (it differs from $c_1+b^2$ by a factor $D^2$ on the time leg).

### 2.3 Numerical validation (section D) [MEASURED]

Witness: the Taylor coefficients of the smeared position-space fibre at $\beta=0$ are
$(-1)^n m_n/(2n)!$ with $m_n = \int k^{2n} j_0(k|W|) P_D(c,k)dk$. Removing the exact algebraic
prefactor ($S_n = m_n/(2n-4)!$ for $|W|>0$, $m_n/(2n-3)!$ for $|W|=0$) gives
$|S_n|^{1/n}\to 1/|z^*|$; the **sign-flip spacing** independently fixes
$\arg(p_1-i|W|)$. Twelve coefficients each ($n \le 16$, dps 30):

| $c$ | $\lvert W\rvert$ | $D$ | $p_1$ | predicted $\lvert z^*\rvert$ | fitted | rel.d | phase pred | phase meas | rel.d |
|---|---|---|---|---|---|---|---|---|---|
| 0.2 | 0 | 1 | 0.4472 | 0.20 | 0.20126 | 0.63% | -- | -- | -- |
| 0.05 | 0 | 1 | 0.2236 | 0.05 | 0.050292 | 0.58% | -- | -- | -- |
| 0.2 | 0.5 | 1 | 0.4472 | 0.45 | 0.45827 | 1.8% | 0.8411 | 0.8727 | 3.8% |
| 0.1 | 0.5 | 1 | 0.3162 | 0.35 | 0.36795 | 5.1% | 1.0069 | 1.0472 | 4.0% |
| 0.2 | 1.0 | 1 | 0.4472 | 1.20 | 1.21691 | 1.4% | 1.1503 | 1.0996 | 4.4% |
| 0.2 | 0.5 | **0.5** | 0.2236 | **0.30** | 0.290585 | 3.1% | 1.1503 | 1.0996 | 4.4% |
| 0.2 | 0.5 | **2** | 0.8944 | **1.05** | 1.087523 | 3.6% | 0.5097 | 0.5236 | 2.7% |
| 0.1 | 0.8 | **1.5** | 0.4743 | **0.865** | 0.857695 | 0.85% | 1.0356 | 1.0472 | 1.1% |

The table also decides *endpoint* against *interior* pinch, which is the one step of 2.1 that
is a judgement rather than algebra: an interior-crossing reading would give
$|z^*| = p_1^2$ alone (0.20 at row 3, 0.05 at row 6), against measured 0.458 and 0.291 --
off by factors 2.3 and 5.8, while the endpoint form lands within a few per cent.

The three general-$D$ rows are decisive for 2.2: at $D=0.5$ the naive $c+b^2$ would predict
0.45 (measured 0.291, a 35% miss) and the interval form predicts 0.30 (3.1%); at $D=2$ the
naive form predicts 0.45 against a measured 1.088 and the interval form 1.05 (3.6%). The few
per cent residual is the finite-order fit over 12 coefficients (the $b=0$ rows, free of the
oscillation, sit at 0.6%). This is also an **independent confirmation of the memo's
Borel-Pade determination of $z^*$** by an unrelated method (Taylor radius of the
position-space fibre rather than Borel-Pade of the $c_t$-series).

### 2.4 Why there is no median ambiguity

In Euclidean signature a nonzero interval never vanishes for real $(p_1,|W|)$, so $z^*$ is
strictly off the positive real Borel axis and the Laplace contour is unobstructed. **The
absence of a Stokes ambiguity in the two-scale resummation is Euclidean positivity of the
interval** -- and it must fail exactly where the interval closes, which is the corner below.

### 2.5 The $(0,0)$ corner is the coincidence (ultraviolet) point

At $(s,t)\to(0,0)$ both $p_i\to0$ and $|W|\to0$, so $|z^*|\to0$: the Borel singularities
accumulate at the origin. Kinematically the interval closes -- **the UV / short-distance
point** -- and simultaneously the masses $a_i = 1/\sqrt{c_i}$ diverge. Both readings agree:
short distance = heavy mass = UV.

Duffy at the corner, $s=\sigma^2 a$, $t=\sigma^2(1-a)$: $p_i \sim \sigma$ but
$|W| = s+t = \sigma^2$, so the spatial separation closes *faster* than the time offsets and
$k|W| \sim \sigma \to 0$: **the oscillation switches off**. Power counting on
$J = \int dk\,j_0 P_1P_2$ with $k\sim 1/\sigma$: each density contributes $c_i = m_i^{-2}
\sim\sigma^2$, the loop range contributes $\sigma^{-1}$, hence

```
J = sigma^3 A(a) + O(sigma^5),   A(a) = a(1-a) int_0^inf G(D_a) G(D_{1-a}) dkap
```

with $G(\Delta)=e^{-\Delta}(\Delta^{-3}+3\Delta^{-4}+3\Delta^{-5})$, $D_a=\sqrt{a\kappa^2+1}$.
Validated: $A(0.5)=0.96981365880091$, $A(0.3)=0.82209185304782$, residual
$|J/\sigma^3 - A|$ falling exactly as $\sigma^2$ ($7.7\!\times\!10^{-3}\to1.9\!\times\!10^{-3}
\to4.8\!\times\!10^{-4}\to9.8\!\times\!10^{-5}$ over $\sigma=0.1\ldots0.0125$). **So the
$\rho^{3/2}$ non-analyticity of the outer integral (with $\rho=s+t=\sigma^2$) is just the odd
ultraviolet mass dimension $m^{-2}m^{-2}m$ of the two-propagator overlap.** $A(a)$ carries
**no $\pi$**, and its curve $y^2=(a\kappa^2+1)((1-a)\kappa^2+1)$ is genus one for $a\ne\tfrac12$:
the mass *ratio* survives the UV limit even though both masses diverge, so **the Duffy angle
$a$ at the corner IS the elliptic modulus** ($\rho_{\rm mod} = (1-a)/a$) and the radial $\sigma$
is the overall scale. That is precisely the co-area split (modulus outer, scale inner) seen
locally at the UV point -- an independent confirmation that the co-area reduction is the
kinematically natural one.

### 2.6 New: the three oscillatory corners, in closed form [MEASURED, 1e-4 and converging]

At the corners with $|W|\ne0$ the same power counting applies but the oscillation survives and
$\int_0^\infty \sin(\Lambda\kappa)h(\kappa)\,d\kappa/\kappa \to \tfrac{\pi}{2}h(0)$ (Dirichlet)
supplies one extra power of $\sigma$ and an explicit $\pi$:

```
J = sigma^4 * (pi/2) G(1)^2 alpha(1-alpha) / |W| + O(sigma^5),     G(1) = 7/e .
```

Validated at four (corner, $\alpha$) combinations -- $(1,0)$ with $\alpha=0.5,0.3$ and
$(1,1)$ with $\alpha=0.5,0.35$, the latter testing the $1/|W|$ dependence at $|W|\to2$ --
each over four $\sigma$, with the relative deviation falling as $\sigma^2$:
$6.4\!\times\!10^{-3} \to 1.6\!\times\!10^{-3} \to 4\!\times\!10^{-4} \to 1\!\times\!10^{-4}$
uniformly. Example: $(1,1)$, $\alpha=0.5$, $\sigma=0.01$: $J/\sigma^4 = 1.30201291159$ vs
predicted $1.30214312264$.

This closes, in closed form, the paper's previously-measured statement that "the three
oscillatory corners are integer-leading (leading order $\rho^2$)": the coefficient is now
explicit, and the contrast with the $(0,0)$ corner is the whole point --
**oscillation on $\Rightarrow$ integer power and an explicit $\pi$; oscillation off (coincidence)
$\Rightarrow$ half-integer power and no $\pi$.**

---

## 3. Step (iii) -- the Paper 35 tie, honestly

### 3.1 The 2x2 table (validated, section F)

With $L(D,\rho) = \int_1^\infty e^{-Dx}dx/\sqrt{(x^2-1)(\rho x^2+1-\rho)}$ (`eq:laplace` up to
$\sqrt{c_1}$):

| | $\rho \in\{0,1\}$ (equal or infinitely-split masses) | $\rho$ generic (split masses) |
|---|---|---|
| **$D=0$** (coincident Euclidean times) | $L=\pi/2$ at $\rho=1$ (exact) -- genus 0, **regular**, elementary | $L = K(1-\rho)$; at $\rho=\tfrac12$, $1.854074677301372 = K(\tfrac12)$ (exact) -- genus 1, **regular** period |
| **$D\ne0$** (Euclidean time displacement) | $L = K_0(D)$ at $\rho=0$ (rel.d $\le 3\!\times\!10^{-18}$ at $D=0.7,1,2$) -- genus 0, **irregular** (classical Bessel resurgence) | the new object; $L(1,\tfrac12)=0.3568617578026771$ -- genus 1 **and** irregular |

Supporting checks: the $D=0$ shadow is annihilated by the Legendre operator
$\rho(1-\rho)\partial_\rho^2+(1-2\rho)\partial_\rho-\tfrac14$ to $\le10^{-31}$ (Fuchsian =
regular); the large-$D$ Watson series is Gevrey-1 at **every** $\rho$ including $\rho=0$, with
Borel radius $R = \min(2,\,1/\sqrt\rho) = \min(2,\, m_2/m_1)$ -- measured (root-test fit over
80 exact series coefficients) $R = 2.02, 1.36, 1.07$ at $\rho = 0.1, 0.5, 0.9$ against
predicted $2, 1.4142, 1.0541$ (1--4%, finite-order).

### 3.2 Is the partition clean? -- YES at the on/off level, TRIANGULAR at the data level

* **Genus depends on $\rho$ only.** The curve $y^2=(c_1k^2+1)(c_2k^2+1)$ contains no $D$ at
  all. Exact, not measured.
* **Irregularity depends on $D$ only.** At $D=0$ the object is a period of a Fuchsian
  (Legendre) equation; at $D\ne0$ it is Gevrey-1 for every $\rho$, including the genus-zero
  slice $\rho=0$ where it is literally $K_0(D)$.
* **But the Stokes data is not factorized.** The Borel radius $\min(2, m_2/m_1)$ and the
  singularity locations depend on the mass ratio; the periods do not depend on $D$. The
  dependence is therefore *triangular*, not a product: (genus $\leftarrow \rho$),
  (Stokes data $\leftarrow \rho$ and $D$).

So the candidate sharpening should be stated as: **mass split turns on the elliptic level;
Euclidean-time displacement turns on the irregular level; the elliptic modulus then also
grades the irregular data.**

### 3.3 Where $\pi$ enters -- a corroboration of WH7/Paper 35, not a refinement of it

Every $\pi$ in the object is traceable, and none of them comes from the Euclidean-time axis:

* $2\pi$ / $2/\pi$ in `eq:K0`, `eq:besselmoment` and Parseval: the momentum-space **measure**
  $d^2k/(2\pi)^2$ of the 2D propagator (equivalently, the $d^3k/k^2$ Coulomb kernel collapsing
  the 3D measure to a 1D one). Spectral-integration $\pi$.
* $\pi/2$ at the genus-zero cusp $\rho=1$, $D=0$: $\int_0^\infty dk/(ck^2+1)$ -- integration
  over the whole momentum half-line.
* $\pi/2$ in the oscillatory-corner coefficient (2.6): the Dirichlet limit of the *same*
  half-line integration, in the ultraviolet.
* $\sqrt\pi$ in $K(\tfrac12)=\Gamma(\tfrac14)^2/4\sqrt\pi$: Chowla-Selberg, already tagged in
  Paper 59.
* $1/\sqrt\pi$ in the Stokes amplitude $(c_1+b^2)^{3/2}/(4\sqrt\pi\sqrt{c_1}b)$: a
  $\Gamma(\tfrac52)$, i.e. the master Mellin engine at half-integer argument.

The Euclidean-time displacement $D$ contributes **no $\pi$** -- it contributes $e^{-D}$, the
Gevrey-1 divergence and the Stokes triple. That is exactly what WH7 predicts for a *non-compact*
temporal direction: $D$ here is a displacement, not a compactified circle, so no $2\pi$
appears; the compactness/discreteness -> $\pi$ implication is untouched. Honest scope: this
**corroborates** the WH7 reading and adds a concrete instance in which the non-compact
temporal axis is shown to inject irregular (exponential/resurgent) content instead of $\pi$;
it is **not** a new mechanism and it is **not** the Matsubara $2\pi$ of Paper 35 (which is a
statement about compactified time). Calling it a "refinement of Paper 35" would overclaim; a
one-line corroborating cross-reference is what the evidence supports.

**And one caveat that cuts the other way.** Euclidean signature does not distinguish
time from space: the assignment "$D$ = Euclidean time, $b$ = space" is forced by the
*structure* of the correlator ($D$ is the fixed shift carried by each propagator; $b$ is
the coordinate the two share and that is integrated out) but not by any metric. So the
statement above is really "the *shift* axis injects exponentials, the *integrated* axis
injects $\pi$" -- a statement about which coordinate is summed over, dressed in
temporal language. That is itself consistent with the project's own signature-blindness
findings (WH7, the K$^+$/Toeplitz probe): nothing here can see a signature, and the
temporal reading should be offered as a reading, not as a metric fact.

### 3.4 Transcendental tags (house rule)

| constant | where | Paper 18 tier | Paper 34 chain |
|---|---|---|---|
| $\pi$ (measure) | 2D propagator normalisation, Parseval | embedding, genus-graded; M1-type measure $\pi$ | momentum-space measure / spectral integration |
| $\pi/2$ (Dirichlet) | oscillatory-corner coefficient (**new**) | embedding, genus 0 | half-line spectral integration in the UV limit |
| $1/\sqrt\pi$ | Stokes amplitude | master Mellin engine at half-integer argument ($\Gamma(\tfrac52)$) | same family as the half-integer Hurwitz entries |
| $e^{-\zeta D}$, Gevrey-1, Stokes triple | the twist | irregular / exponential-period tier (Paper 59 `sec:obstruction`) | **new candidate row**: non-compact Euclidean-time displacement -> exponential period |
| $K(\tfrac12)=\Gamma(\tfrac14)^2/4\sqrt\pi$, $G$ | period / quasiperiod ring | elliptic (genus 1), CM $\Gamma$-value | already carried by Paper 59 `sec:modular` |

The only genuinely *new* tag is the fourth row: an exponential/irregular entry keyed to a
non-compact Euclidean-time displacement, distinct from the compactification $2\pi$ entry.

---

## 4. Honest scope

**Exact / derived (high confidence).** The dictionary of Section 1 is exact algebra, and the
three identities that could hide an error (D1, D2, the position-space packet) are validated at
$10^{-30}$, $10^{-31}$ and $10^{-18}$. The endpoint-pinch derivation of $\beta^*$ (2.1) is exact
given 1.3-1.4. The corner power counts (2.5, 2.6) are exact asymptotics.

**Measured.** $|z^*| = p_1^2+|W|^2$ including the general-$D$ refinement, at eight
$(c,|W|,D)$ points, 0.6-5.1% on the modulus and 1.1-4.4% on the phase from 12 coefficients;
the oscillatory-corner coefficient to $10^{-4}$ with a clean $\sigma^2$ approach; the 2x2 table.

**Not established.** (a) No theorem: $\beta^*$ is derived from the leading endpoint behaviour
plus contour deformability, then validated -- not proved. (b) Nothing here produces the T2
closed form; the family-integration obstruction of the predecessor memo stands, and 2.5 only
*explains* it (the interval closes at $(0,0)$, so the Borel singularities collide at the origin
and the "no median ambiguity" argument of 2.4 fails exactly there). (c) The $\pi$ analysis is a
corroboration of WH7, not a new mechanism. (d) The `sec:genus` care point of 1.7 is flagged,
not tested.

**Negative / caution.** The naive $D$-independent reading $|z^*| = c_1+b^2$ is *wrong away from
$D=1$* (35% at $D=0.5$); anyone reusing the closed form off the physical twist must use the
interval form.

---

## 5. PROPOSED Paper 59 paragraph -- NOT APPLIED, for PI decision

Suggested placement: end of `sec:reduction` (the dictionary half) with the second and third
paragraphs in `sec:modular` after the existing Stokes paragraph. Tier tags as marked.
**Application note:** the final sentence cites `loutey_paper35`, which is NOT currently in
Paper 59's bibliography (only `loutey_paper18` and `loutey_paper34` are) -- either add the
bibitem or retarget the cite to Paper 18/34.

> \textbf{[MEASURED]} \emph{The kinematics.} Equation~\eqref{eq:K0} is the frequency
> representation of the two-dimensional Euclidean propagator: with
> $G_{a}(\tau,x)=(2\pi)^{-1}K_{0}(a\sqrt{\tau^{2}+x^{2}})$ the Green's function of
> $-\Delta+a^{2}$ on $\mathbb{R}^{2}$, its left-hand side is exactly
> $(2\pi/\sqrt{c})\,G_{a}(p,b)$ with mass $a=\zeta/\sqrt{c}$, Euclidean time $p=D\sqrt{c}$ and
> space $b$, and the fibre~\eqref{eq:besselmoment} is exactly the two-propagator correlator
> $M=4\pi a_{1}a_{2}\int_{\mathbb{R}}db\,G_{a_{1}}(p_{1},b)\,G_{a_{2}}(p_{2},b)$: two Euclidean
> packets of unequal mass, evolved from a point source for Euclidean times $p_{1},p_{2}$ and
> overlapped on one spatial slice \textbf{[MEASURED, $10^{-30}$ and $10^{-31}$ at five parameter
> points each]}. Three things become visible. First, the Feynman parameters are
> Källén--Lehmann spectral variables: $c=s(1-s)$ forces $a=\zeta/\sqrt{c}\ge2\zeta$, so each
> two-centre density is dispersed over its own two-particle cut and the elliptic modulus
> $\rho=c_{2}/c_{1}=(a_{1}/a_{2})^{2}$ is the squared ratio of the two spectral masses---the
> genus jump is a mass-splitting phenomenon, and the three cusps are equal masses ($\rho=1$)
> and one mass at the heavy end of its cut ($\rho\to0,\infty$). Second, the exponential twist
> is $a_{i}p_{i}=\zeta D_{i}$, \emph{independent of the Feynman parameter}: the family is swept
> at fixed mass$\times$time, and the physical point is one Compton wavelength of Euclidean-time
> displacement, the $u=D$ endpoint of the source-time superposition
> $e^{-D\Delta}\Delta^{-n}=\frac{1}{(n-1)!}\int_{D}^{\infty}(u-D)^{n-1}e^{-u\Delta}\,du$ that
> converts the Yukawa into the Slater. Third, $j_{0}(k|\bm W|)$ is a uniform smearing of the
> shared spatial coordinate over the composite separation, so the third centre enters
> \emph{only} through space while the bond lengths set the two Euclidean-time offsets
> \textbf{[MEASURED, $1.8\times10^{-10}$ for the full position-space fibre]}.
>
> \textbf{[MEASURED]} \emph{The Stokes location is a complexified light cone.} In position
> space each packet is analytic in $b^{2}$ with branch points at the endpoint of its
> source-time superposition, where the Euclidean interval closes, $b^{2}+p_{i}^{2}=0$; the
> $j_{0}$ smearing moves the singularity to the endpoint pinch
> $\beta^{*}=\pm|\bm W|\pm i p_{1}$, so the Borel location above is the squared complexified
> interval, $z^{*}=(\beta^{*})^{2}$, $|z^{*}|=p_{1}^{2}+|\bm W|^{2}$, of which
> $-(\sqrt{c_{1}}\mp ib)^{2}$ is the $D=1$ specialisation $p_{1}=\sqrt{c_{1}}$; the trans-series
> action $a_{2}\sqrt{z^{*}}$ is the heavy mass times that interval. Moment growth of the
> position-space fibre confirms modulus and phase at eight $(c,|\bm W|,D)$ points, including
> $D\neq1$, where the interval form and the $D$-independent form differ by $35\%$
> \textbf{[MEASURED, $0.6$--$5\%$ from twelve coefficients]}. Two consequences. The median
> resummation is unambiguous because a Euclidean interval never vanishes at real separation:
> the absence of a Stokes ambiguity is Euclidean positivity. And the $(s,t)\to(0,0)$ corner is
> the coincidence (ultraviolet) point, where the interval closes and the Borel singularities
> collide at the origin---which is why that corner, and only that corner, obstructs the family
> integration. There the masses diverge at fixed ratio, so the modulus survives and the Duffy
> angle \emph{is} the modulus; the leading behaviour is the ultraviolet power count
> $J=\sigma^{3}A(a)+O(\sigma^{5})$, $\sigma^{3}=m_{1}^{-2}m_{2}^{-2}m$, so the $\rho^{3/2}$
> non-analyticity is an odd mass dimension and carries no $\pi$. At the three corners with
> $|\bm W|\neq0$ the oscillation survives, a Dirichlet limit supplies one further power and an
> explicit $\pi$, and the integer-leading behaviour is now closed:
> $J=\sigma^{4}\,\frac{\pi}{2}G(1)^{2}\alpha(1-\alpha)/|\bm W|+O(\sigma^{5})$ with $G(1)=7/e$
> \textbf{[MEASURED, $10^{-4}$ with $O(\sigma^{2})$ approach, four corner/angle combinations]}.
>
> \textbf{[OBSERVATION]} Read against the two-layer taxonomy, the object separates on two
> kinematic axes: the mass split turns on the elliptic (genus-one) level and the Euclidean-time
> displacement turns on the irregular level---the curve carries no $D$, and the large-$D$ series
> is Gevrey-one at every $\rho$, including the genus-zero slice $\rho=0$ where the fibre is
> $K_{0}(D)$---while the modulus then also grades the irregular data through the Borel radius
> $\min(2,\,a_{2}/a_{1})$, the smallest difference of the four exponential rates. Every $\pi$
> in the object comes from the momentum-space measure or from integration over the momentum
> half-line; the non-compact Euclidean-time direction contributes no $\pi$ but the exponential
> and its Stokes data, consistent with the compactification reading of
> Paper~35~\cite{loutey_paper35}.

---

## 6. Files

* `debug/t2_euclidean_dictionary.py` -- sections A (propagator identity), B (two-propagator
  correlator + period), D (light-cone/moment growth incl. general $D$), E (corner scalings),
  F (2x2 table, Borel radius, Legendre check). `python debug/t2_euclidean_dictionary.py ABDEF 30`
* `debug/_t2_dict_secC.py` -- C0 kernel identity, C1 position-space packet, C2 unsmeared
  overlap, C3 full smeared fibre (fixed Gauss-Legendre so the nesting is affordable).
* `debug/_t2_dict_secC3.py` -- C3 on a lighter grid.
* `debug/data/_t2_dict_sec{C,C3,C3b,D}.out`, `_t2_dict_ABEF.out` -- run logs.
* Reproduce: `python debug/t2_euclidean_dictionary.py ABDEF 30` (about 30 s) and
  `python debug/_t2_dict_secC.py 22` (C0-C2 about 4 min; C3 about 8 min per point).
* Nothing committed; no paper modified; PI controls commit/release.
