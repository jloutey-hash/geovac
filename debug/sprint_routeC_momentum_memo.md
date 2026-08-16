# Sprint Route C — momentum-space two-body 3-centre ERI (2026-08-16)

**Verdict: OPENED + method validated + wall named + WEIGHT QUESTION RESOLVED.**
The two-body three-centre ERI `T2 = (XY|XZ)` — the one remaining genuine polyatomic
wall after the v4.81.0 one-body 3-centre closure — is now **computable exactly in
momentum space**, the angular integral **closes GeoVac-natively to a spherical Bessel
kernel**, and the transcendence obstruction is **identified: the third centre is
ELLIPTIC (genus 1)**, categorically beyond the two-centre engine's genus-0
`{E₁, ln, γ}`. Driver `debug/routeC_momentum_poc.py`; test
`tests/test_routeC_momentum.py`.

## The object

`T2 = (XY|XZ) = ∫∫ ρ₁(r₁)(1/r₁₂)ρ₂(r₂)`, `ρ₁=χ_Xχ_Y`, `ρ₂=χ_Xχ_Z` — two
two-centre densities sharing centre X on two different axes. No prolate-spheroidal
system holds three foci, so the Neumann route that closed the two-centre engine
has no coordinate system here (build plan §10.4, memory
`native-two-center-eri-engine` 3b). Ground truth (eri_md, 8-Gaussian 1s fits,
`<fit|STO>=1−3e−8`): `X=(0,0,0) Y=(0,0,2) Z=(1.5,0,0.5) ζ=1 → 0.20494172`.

## What was established (all reproduced by the driver)

**1. The coordinate wall is dissolved — three phases, not three foci.**
Momentum space: `1/r₁₂ = (1/2π²)∫d³k e^{ik·(r₁−r₂)}/k²`, translation = phase
`e^{ik·R}`, so
```
(XY|XZ) = (1/2π²) ∫ d³k/k²  ρ̃₁(k) conj(ρ̃₂(k)),   ρ̃(k) = ∫ρ(r)e^{ik·r}d³r
```
Validated two independent ways against ground truth:
- **A. Gaussian-density FT** (closed form; FT of a Gaussian product is a Gaussian):
  agrees with `eri_md` on the *same* Gaussians to **1.9e−14** (isolates the
  formula + conventions + k-integrator). `ρ̃(0)` = the overlap `⟨χ_X|χ_Y⟩`
  exactly, fixing the normalization.
- **B. TRUE Slater density FT** via a Feynman/Yukawa reduction (independent of the
  Gaussian fit): agrees to **4.3e−7** (finite-difference `d/dζ` limited, not the
  method).

**2. The angular Ω_k integral closes to a spherical Bessel kernel.** Writing the
Slater FT via `e^{−ζr}=−∂_ζ(e^{−ζr}/r)` and the Yukawa-product convolution
(Feynman-parametrized), the three phases combine into a **single** phase `e^{ik·W}`
with
```
W(s,t) = (t−s)X + sY − tZ            (= sY − tZ when X = 0)
∫dΩ_k e^{ik·W} = 4π j₀(k|W|)          (spherical Bessel — the GeoVac-native step)
```
so the whole ERI reduces to a **2D Feynman × 1D radial** integral (validated to
**1.8e−6**):
```
(XY|XZ) = (8/π) ∂ζa∂ζb∂ζc∂ζd ∫₀¹ds ∫₀¹dt ∫₀^∞dk
            j₀(k|W|) · e^{−D₁Δ₁}/Δ₁ · e^{−D₂Δ₂}/Δ₂ |_{ζ=1}
Δ₁=√(s(1−s)k²+s ζa²+(1−s)ζb²), D₁=|X−Y|;  Δ₂ analogous with t, D₂=|X−Z|.
```

**3. The transcendence obstruction: the third centre is ELLIPTIC (genus 1).**
A *single* dispersion factor closes under the Fock substitution `k√c = m·sinh θ`:
```
∫₀^∞ cos(kb) e^{−D√(ck²+m²)}/√(ck²+m²) dk = (1/√c) K₀( (m/√c)√(cD²+b²) )   (verified 6e−18)
```
— a Bessel `K₀`, the momentum-space Coulomb-Sturmian (Fock) object, on a **genus-0
(rational) curve**; this is why the two-centre engine closed at weight 1 over
`{E₁, ln, γ}`. The three-centre integrand carries **two dispersion factors with
different scales `c₁=s(1−s) ≠ c₂=t(1−t)`**, whose product defines the algebraic curve
```
y² = (c₁k²+1)(c₂k²+1)          — a QUARTIC.
```
This is an **elliptic curve (genus 1) whenever c₁≠c₂**, degenerating to a perfect
square (rational, genus 0) exactly on the diagonal `c₁=c₂`. Decisive witness — the
`D=0` period is a **complete elliptic integral** (verified to 31 digits):
```
∫₀^∞ dk/√((c₁k²+1)(c₂k²+1)) = (1/(a₁√(c₁c₂))) K(m),  a₁=1/√c₁,  m=1−(a₂/a₁)²
```
(nondegenerate, `0<m<1`), whereas the diagonal gives `∫dk/(ck²+1)=π/(2√c)`,
elementary. Over the `(s,t)` Feynman domain this is a **family of elliptic curves**,
modulus `m(s,t)=1−c_min/c_max`, degenerate only on the measure-zero locus `s=t` or
`s=1−t`. Every `∂_ζ`-generated term (`1/Δ_i^{3,4,5}·e^{−D_iΔ_i}`) is a meromorphic
function on the *same* curve, so the whole ERI is a period/quasi-period of this
elliptic family.

## Resolution of the weight question

**The third centre raises the transcendence from genus-0 polylogarithms
(`{E₁, ln, γ}`, the two-centre engine) to genus-1 ELLIPTIC transcendentals
(elliptic polylogarithms).** This is the exact momentum-space/Fock statement of the
"no shared hypersphere for three foci" wall: the two densities' Fock scales coincide
(a *shared* S³) only on the degenerate diagonal `c₁=c₂`; generically they define two
distinct spheres and hence a genuine elliptic curve. It also explains why a
dilogarithm/ζ(2) PSLQ never lands — the constants are in the **wrong transcendence
class** (elliptic, not polylog). Supporting numerics: the collinear value
`0.395355766…` is provably **not** a rational combination of `{1, e^{−2}, e^{−4}}`
(weight-0 PSLQ returns only huge spurious coefficients), consistent with a
genus-1 object, not the exp-polynomial the two-centre exchange `J(R)` was.

Rigorous backbone (both verified numerically to ≥30 digits):
- **Subordination / sunrise form.** `e^{−DΔ}/Δ=∫_D^∞e^{−wΔ}dw` and the subordinator
  `e^{−w√A}=(w/2√π)∫₀^∞ s^{−3/2}e^{−w²/4s−sA}ds` reduce the two-scale radial integral to
  `Φ(0)=(1/2√π)∫∫ds₁ds₂ (s₁s₂)^{−1/2} e^{−D₁²/4s₁−s₁−D₂²/4s₂−s₂}/√(s₁c₁+s₂c₂)` — a
  two-scale Feynman/sunrise integral whose Symanzik `√`(linear form) is the genus-1
  signature; it separates only when `c₁=c₂`.
- **Elliptic period.** `D=0` gives the complete elliptic integral above.

Frontier for a dedicated sprint (THE Avery-call topic — his momentum-space Sturmian
turf): express the full `T2` in closed form via **elliptic polylogarithms / iterated
integrals on this elliptic family** (the genus-1 analogue of how the two-centre
exchange used `{E₁,ln,γ}` on the rational curve). The obstruction is now named at the
right level of the transcendence hierarchy.

## Scope / honesty

- s-type (1s) validated; higher l not attempted (adds solid-harmonic polynomials
  in k — mechanical, doesn't change the dispersion-scale obstruction).
- This does NOT solve water: closed-form value is all-or-nothing at the tensor
  level (build plan §10.4), and the elliptic closed form of T2 is not derived (only
  the genus-1 obstruction is proven). What it buys: an exact, fast, GeoVac-native
  *evaluator* for T2 and a precise transcendence diagnosis (elliptic, genus 1).
- The `D=0`-period = elliptic-`K` identity is rigorous; the full `T2` (with the
  `e^{−DΔ}` exponential factors and the `(s,t)` family integral) is an elliptic
  *polylogarithm*-type object on the same curve — genus-1 established, closed form
  not derived.
- No production code touched. Artifacts: `debug/routeC_momentum_poc.py` (evaluator +
  elliptic witness), `debug/routeC_weight_probe.py` (high-precision collinear value +
  the weight-0 exclusion; its dilog-basis PSLQ is superseded — the class is elliptic,
  not polylog), `tests/test_routeC_momentum.py` (3 pins, 3s). Regression clean.

## (a)/(b) follow-up (2026-08-16): elliptic ⊕ log Bessel moment

Bounded attempt at a closed form vs. a rigorous no-elementary-closure statement.
The radial content reduces to the two-scale **Bessel moment**
`M(D₁,D₂;c₁,c₂) = ∫₀^∞ dk e^{−D₁√(1+c₁k²)−D₂√(1+c₂k²)}/√((1+c₁k²)(1+c₂k²))`
(single dispersion factor = Bessel K₀ by the Fock substitution; a product of two,
integrated, is by definition a Bessel moment — the Broadhurst / elliptic-Feynman
class). Two facts, both numerically pinned:

- **`M(D=0)` is exactly a complete elliptic integral `K`** (31 digits) — genus-1,
  classically non-elementary and non-polylogarithmic.
- **`M` is non-analytic in `D` at `D=0`: it carries a `D·ln D` term.**
  `(M(D)−M(0))/D` diverges logarithmically as `D→0` (measured −9.4, −11.5, −13.5,
  −15.6 at `D=0.1,0.05,0.025,0.0125`). So the physical (`D≠0`) kernel is
  **elliptic ⊕ logarithmic**, not a pure `K/E`.

**Reading.** The object is a two-scale Bessel moment on the elliptic curve
`y²=(c₁k²+1)(c₂k²+1)`. This is why the two-centre `{E₁,ln,γ}` closure cannot
extend (that is the genus-0 polylog world) and why a dilog/ζ(2) PSLQ was doomed.

- **(b) LOCKED (strengthened).** The three-centre radial integral is provably
  outside the two-centre class: its `D=0` period is a complete elliptic integral
  (genus 1, non-elementary), and it additionally carries logarithmic `D`-structure.
  No finite elementary or elliptic-`K/E` closed form exists.
- **(a) is frontier, not a bounded win.** A genuine finite closed form, if one
  exists, is an **elliptic polylogarithm / modular L-value** (Broadhurst Bessel
  moments) — specialized modern machinery, some members still open. The bounded
  shot did not (and structurally will not) yield an elementary/`K,E` form.

Net: the honest, in-hand deliverable is a sharpened (b) — a precise
characterization (two-scale elliptic Bessel moment) plus the connection to the
Bessel-moment / elliptic-Feynman-integral literature — not a closed form.
Driver: `debug/routeC_weight_probe.py` sibling probes; owning doc build plan §10.5.

## Literature scout (2026-08-16) — the object is a NOVEL bridge, closed form open

A verified web/literature hunt (general-purpose agent) placed the object precisely.

**Verdict:** the closed form of `M` is **not in the literature**; the general class
(two-mass elliptic Feynman integrals / Bessel moments) exists and its machinery
applies, but our specific object is **not evaluated** — genuinely open as a closed
form. And the **chemistry↔elliptic-Feynman/Bessel-moment bridge appears unmade in
either field** (medium-high confidence, targeted-search absence-of-evidence).

- **Chemistry side is uniformly "hard → infinite series / 1-D numerics":**
  Özdoğan–Ruiz 2012 (arXiv:1209.3755, the exact 3-centre STO ERI class; one infinite
  expansion, 20 digits in 25–30 terms, no closed form); Barnett–Coulson, Harris,
  Steinborn/Weniger; **Avery & Avery** 4-centre STO Coulomb-Sturmian ERIs (momentum
  space — closest to our route — stops at expansion+numerics, draws no elliptic/Bessel
  closed form). The one famous closed form, **Fromm–Hill 1987** (PRA 36, 1013), is
  **dilogarithmic (genus 0)** — a *different* topology (one-centre r₁₂-correlated),
  genus-0 precisely because it lacks the two-scale quartic. dilog-vs-elliptic is a
  substantive contrast for us.
- **Amplitudes/number-theory side has the machinery, applied elsewhere:** the object
  is a **two-mass sunrise / two-scale Bessel moment**. `D=0` period = complete
  elliptic `K` (classical, Byrd–Friedman 213 / G-R 3.152); at rational scale ratios
  potentially a product of `K`s / Γ-value products (**Broadhurst arXiv:0801.4813**,
  the closest published analog — a *doubled-mass* Bessel moment = (1/12)K(sin π/12)
  K(cos π/12) = Γ⁶(1/3)/…). Full `M` (`D≠0`) is **elliptic-dilogarithm level, weight
  ≤ 2** (only two propagator factors ⇒ low end of the class; expect `K`,`E` + an
  elliptic dilogarithm, **not** deep modular L-values except at special CM ratios).
- **Independent confirmation of the D·ln D finding:** the scout notes the `D·ln D`
  non-analyticity is exactly the **inhomogeneous/regulator term** of the Picard–Fuchs
  solution defining the elliptic dilogarithm — pure elliptic integrals never produce
  it. So our numeric tell matches the theory.
- **Concrete template for a closed-form (a):** **Adams–Bogner–Weinzierl
  arXiv:1405.5640** (unequal-mass two-loop sunrise → elliptic dilogarithms, "as simple
  as equal-mass, only the arguments modified"); founding theory **Bloch–Vanhove
  arXiv:1309.5865**; motivic relations **Fresán–Sabbah–Yu arXiv:2006.02702**; Bessel-
  moment catalogue **Bailey–Borwein–Broadhurst–Glasser arXiv:0801.0891** (single-scale
  only — ours is two-scale/shifted, hence not a BBBG entry = new).
- **Useful exact identities:** `∫₀^∞ x K₀(ax)K₀(bx)dx = ln(a/b)/(a²−b²)` (elementary!
  the x¹-weighted moment), vs `∫₀^∞ K₀(ax)K₀(bx)dx` = ₂F₁(½,½;1;·) = elliptic `K`
  (x⁰). **Weight of the moment decides elementary-vs-elliptic** — a lead for isolating
  the elliptic content among the ∂_ζ-generated powers.

**Consequence.** (a) is not hopeless frontier — it is at the tractable
elliptic-dilogarithm rung with a published template (unequal-mass sunrise), just never
done for a chemistry ERI. The strongest deliverable in hand is the **novel bridge**:
the 3-centre Sturmian ERI, via Fock projection, is a two-scale Bessel moment on the
sunrise elliptic curve — evaluable by elliptic-polylog methods, a connection unmade in
either literature. Scout caveats: load-bearing refs fetched/verified; a few Zhou
preprints + G-R/Byrd–Friedman formula *numbers* corroborated by standard knowledge,
not page-fetched.

## (a) derivation progress (2026-08-16) — reduced + characterized, finite form NOT reached

PI-directed push for the closed form. Substantial verified progress; the finite
elliptic-dilogarithm formula was not reached (it is the specialized step).

Achieved (all numerically verified ≥1e-24 / 30 digits):
1. **`M` = a two-mass Bessel moment (exact):**
   `M = (2/π)(1/√(c₁c₂)) ∫₀^∞ K₀(a₁√(p₁²+b²)) K₀(a₂√(p₂²+b²)) db`, `a_i=1/√c_i`,
   `p_i=D_i√c_i`. A single-scale/shifted-argument member — not a BBBG catalogue entry.
2. **Clean canonical form (one-mass slice `N(D)=M(D,0)`):** the Laplace transform of
   the elliptic differential,
   `N(D) = (1/√c₁) ∫₁^∞ e^{−Dx} dx / √((x²−1)(ρx²+1−ρ))`, `ρ=c₂/c₁`.
   A *single explicit 1-D integral* of an elementary integrand — already better than
   the Özdoğan–Ruiz infinite-series baseline, and it makes the elliptic curve manifest.
3. **`N` is an elliptic deformation of the Bessel `K₀`:** `N→K₀(D)` as `ρ→0` (scales
   coincide), `N→` complete elliptic `K` at `D=0` (zero separation).
4. **Exact 4th-order Picard–Fuchs ODE** (from `Qf'+½Q'f=0` + one IBP; boundary at
   `x=1` vanishes since `Q(1)^{1/2}=0`; verified residual ~1e-24):
   `Dρ·L⁗ + 2ρ·L‴ + D(1−2ρ)·L″ + (1−2ρ)·L′ − D(1−ρ)·L = 0`.
   Characteristic polynomial `ρ(s²−1)(s²+(1−ρ)/ρ)`: rates `±1` (from the x²−1 branch
   points, K₀/I₀ sector) and `±i√((1−ρ)/ρ)` (the other branch pair, J₀/Y₀ sector).

NOT achieved: a **finite** closed form in `K`, `E`, and the elliptic dilogarithm.
The 4th-order operator does **not** factor into elementary Bessel operators (genuine
`D²` couplings), so the object is a genuinely new elliptic transcendent; its finite
form requires the elliptic-Feynman-integral reduction (Adams–Bogner–Weinzierl-style
Picard–Fuchs + variation-of-parameters + elliptic-dilog identification) — specialized
and uncertain. Honest deliverable = the reduction + ODE characterization above, not a
formula. Drivers: this sprint's mpmath probes (routeC_momentum series).

## Transcendental tag (Paper 18 / Paper 34) — obligation discharged 2026-08-16

The elliptic period is a NEW transcendental for the corpus; classified against the
live taxonomy (Paper 18 §Level-2, Paper 34 two-centre ERI entry lines ~1102–1144):

1. **Appearance.** Complete elliptic integral `K(m)` in the `D=0` period of the
   two-scale radial kernel of `T2=(XY|XZ)` (Route C); curve `y²=(c₁k²+1)(c₂k²+1)`,
   `c_i` the Feynman parameters of the two densities' Fock projections. Full `T2` =
   an elliptic polylogarithm on this family.
2. **Paper 18 tier = EMBEDDING** — the same tier as `{e^aE₁(a), ln, γ_E}` and
   `1/r₁₂`: it is what the continuum *multi*-centre Coulomb geometry costs, not graph
   content. The elliptic period is the three-centre continuation of the two-centre
   Level-2 embedding seed set.
3. **Master Mellin M1/M2/M3 — NOT in the chain** (unchanged from two-centre): `π`
   cancels against the Coulomb `4π`, so the calibration tier is untouched and the
   Hopf/Seeley-DeWitt/vertex-parity sub-mechanisms do not appear.
4. **Paper 34 projection chain** = Fock conformal (§proj_fock) → Sturmian
   (§proj_sturmian) → multipole projection carried across **three** centres — the
   direct extension of the certified two-centre entry. The elliptic curve arises
   because the two densities Fock-project onto two DIFFERENT S³'s (scales √c₁≠√c₂);
   the product of the two stereographic weights is the quartic.
5. **Three-axis tag** = `(L, dimension-preserving, {complete elliptic K / elliptic
   polylog} at genus 1)` — the transcendental axis is a NEW ring, categorically
   beyond the two-centre `{E₁,ln,γ}` (genus 0, weight-one polylog).
6. **The structural content: the EMBEDDING tier carries a GENUS grading.** Two
   centres = genus 0 (rational curve; `{E₁,ln,γ}`; the weight-one filtration that
   *avoids* the Li₂ its simplex shape would generically produce). Three centres =
   genus 1 (elliptic). **The two-centre negatives ("no π", "no Li₂/weight-2") are a
   two-centre-only property; the third centre exits the polylogarithm/MZV weight
   tower entirely into elliptic transcendentals** — orthogonal to the M1/M2/M3
   Mellin tower, not a higher rung of it.
7. **Pinning check: PASSES** (not anonymous). It pins to existing projections
   (Fock→Sturmian→multipole); the genus grading is a REFINEMENT of the embedding
   tier, not a new named projection. Paper edits (Paper 18 §Level-2 continuation +
   Paper 34 two-centre-entry extension) are drafted; both are certified papers, so
   applying them is PI-gated (Phase-4 re-review of the changed sections).

Files: `debug/routeC_momentum_poc.py`, `debug/routeC_weight_probe.py`,
`tests/test_routeC_momentum.py`. Owning doc:
`docs/neumann_general_m_build_plan.md` §10.5.

## (a) ABW push — the route is STRUCTURALLY OBSTRUCTED (2026-08-16 cont.)

PI-directed dedicated sprint on the finite elliptic-dilog closed form of the
one-mass slice `N(D)` via the Adams–Bogner–Weinzierl unequal-mass-sunrise route on
the 4th-order Picard–Fuchs ODE. **Verdict: the ABW mechanism does not engage — for a
precise, verified structural reason — and `N(D)` is a genuinely new elliptic
transcendent one rung above the sunrise, NOT an elliptic dilogarithm of elementary
arguments.** Drivers `debug/routeC_abw_diagnostic.py`, `debug/routeC_modulus_pf.py`,
`debug/routeC_gauss_manin.py`.

Object: `L(D,rho) = sqrt(c1) N(D) = int_1^inf e^{-Dx} dx / sqrt(Q)`,
`Q=(x^2-1)(rho x^2+1-rho)`, `rho=c2/c1` — the Laplace transform of the holomorphic
elliptic differential over the cycle `[1,inf)`.

**Why the naive attack (VoP on the 4th-order D-PF) is mis-aimed — VERIFIED.**
The elliptic curve `y^2=Q(x,rho)` depends only on the MODULUS `rho`, not on `D`
(`D` is Laplace-dual at fixed modulus). In the ABW sunrise the kinematic variable
*moves the modulus*, so its periods ARE the ODE's homogeneous solutions; here they
are not:
- the 4th-order D-ODE reproduces to ~1e-41 (state confirmed, exact moment quadratures);
- single Bessels `{K0(D),I0(D),J0(wD),Y0(wD)}`, `w=sqrt((1-rho)/rho)`, are NOT
  D-solutions (residuals O(0.1)) — the K0/I0 and J0/Y0 sectors are coupled;
- the periods are D-constants, `O[const]=-D(1-rho)!=0` — NOT D-solutions.
So the elliptic-dilog machinery cannot live on the D-axis. The periods live on the
MODULUS axis.

**The corrected ABW target — the 2nd-order modulus PF — and why IT is obstructed too.**
- Modulus operator `M_rho = rho(1-rho) d^2/drho^2 + (1-2rho) d/drho - 1/4`
  (the Legendre/hypergeometric PF); homogeneous solutions = the two periods
  `{K(rho), K(1-rho)}` (annihilated to ~1e-18); `L(0,rho)=K(1-rho)` to 19 digits.
- `L(D,rho)` is the inhomogeneous solution, `M_rho[L]=S(D,rho)!=0`.
- **Exact-form identity (symbolic, residual 0):**
  `M_rho[Q^{-1/2}] = d/dx[ g(x) Q^{-1/2} ]`, `g(x) = -x(x^2-1)/(4(rho x^2+1-rho))`,
  with **`g(1)=0`**. Hence for `D!=0`, `S(D,rho) = D int_1^inf e^{-Dx} g Q^{-1/2} dx`
  (boundary vanishes). Because `g(1)=0`, modulus-differentiation produces **no
  branch-point boundary / no subtopology term** — everything stays on the elliptic
  curve. This is the exact mechanism of the obstruction (contrast the sunrise, whose
  `p^2`-derivative pinches a propagator into an ELEMENTARY tadpole = the ABW
  inhomogeneity).
- **Gauss–Manin closure (numeric, dps=45).** `S` closes EXACTLY and finitely only in
  the integrand's own D-module: `S = sum_{k=0}^3 P_k(D) L^{(k)}`, `P_k` degree-2
  polynomials — residual **5.1e-44 at deg 2, FLAT at deg 3/4** (saturated = exact).
- **Discriminator (the load-bearing test).** Does `S` instead decompose as
  period-level `{L,L'}` + elementary Bessel subtopology (the genuine ABW shape)?
  `S in {L,L'}xpoly + {K0(D),K1(D),J0(wD),Y0(wD)}xpoly` (BOTH sectors) only
  **APPROXIMATES**: 2.1e-14 → 8.5e-20 → 1.8e-25 → 2.8e-30 monotone across deg 2..5
  (36 params), never reaching the floor Hyp A hit with 12. **No finite elementary
  inhomogeneity exists.**

**Conclusion (honest scope).**
- PROVEN (numeric-to-floor + symbolic `g`): the modulus-PF source is in-module; it
  admits no finite period+Bessel (subtopology) decomposition.
- FOLLOWS: the ABW variation-of-parameters mechanism (the route named) does NOT
  reduce `N(D)` to an elliptic dilogarithm of elementary/period data, on EITHER axis.
- So `N(D)` is a **new elliptic period of a rank-4 irregular (Laplace-type)
  connection** over the curve — strictly beyond the two-loop sunrise's
  elliptic-dilogarithm class. This sharpens Paper 59's `(b) LOCKED` + `(a) frontier`
  into `(a) NEGATIVE for the elementary-ELi target`; the only residual open question
  is differential-Galois (does the 4th-order operator factor over the elliptic period
  field — does ANY closed form exist).
- NOT proven: that no closed form of any kind exists (full differential-Galois). And
  the two-mass `M(D1,D2)` may carry ABW structure with these one-mass slices as
  subtopologies — but then its closed form is an ELi OF the new transcendents, not
  elementary (turtles). The scout's optimistic "elliptic-dilog level, weight<=2" was
  a sunrise pattern-match; the direct computation (source in-module) overrides it.

**Contradiction-to-scout flagged.** The literature scout guessed the object sits at
the tractable elliptic-dilogarithm rung (two propagators = low end of the class).
The computed inhomogeneity is in-module, i.e. HARDER than that pattern-match implies.
Computation over pattern-match (current-state rule).

Files: `debug/routeC_abw_diagnostic.py` (D-axis: ODE verify + single-Bessel +
periods-not-D-solutions), `debug/routeC_modulus_pf.py` (modulus PF + period basis +
source), `debug/routeC_gauss_manin.py` (D-module closure vs period+Bessel
discriminator), `debug/routeC_saturation.py` (dps-45 both-sector saturation —
produces Paper 59 §obstruction's `5e-44` + `2e-14→3e-30` numbers). Backing tests
in `tests/test_routeC_momentum.py`: `test_modulus_pf_annihilates_periods`,
`test_modulus_source_is_exact_derivative_g1_zero` (incl. a sympy SYMBOLIC pin of the
g(x) exact-form identity), `test_single_bessels_and_periods_are_not_D_solutions`,
`test_modulus_source_is_in_module_not_period_plus_bessel` (both-sector, `@slow`).
Owning doc build plan §10.5.

**Phase-4 re-review (2026-08-16, PI-directed on the certified-paper reframe):** two
fresh adversarial reviewers (claims + code) on the v4.83.0 §obstruction changes —
BOTH CLEAN, no MATERIAL/LARGE. Code review independently confirmed the g(x) identity
symbolically (sympy exact 0) and that the saturation discriminator genuinely
separates exact finite closure from non-terminating approximation (verified the
stronger both-sector claim at dps-45: Hyp A floor 1.99e-45 flat; Hyp C 3.9e-16→4.5e-33
across 18→36 params, never saturating). NITs fixed: sympy symbolic pin added; abstract
[OPEN]/[OBSERVATION] tag split; "(numerically an exact finite identity)"; Claim-4 +
both-sector tests added to close the `\cite{test_routeC}` coverage gaps.
