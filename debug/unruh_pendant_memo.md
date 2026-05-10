# Unruh pendant — third face of the Hartle–Hawking → Sewell → Bisognano–Wichmann Wick-rotation chain

**Date:** 2026-05-10
**Sprint:** Unruh pendant (1–2 day pendant follow-up to the Bisognano–Wichmann
landing of 2026-05-09).
**Builds on:** Track D Bisognano–Wichmann reading
(`debug/bisognano_wichmann_track_d_memo.md`); Sprint TD Track 4
(`debug/sprint_td_track4_memo.md`); Sprint TD Track 1
(`geovac/thermal_tensor_triple.py`); Sprint MR-A/B master Mellin engine
domain partition; Paper 38 GH-convergence theorem (WH1 PROVEN).
**Status:** **POSITIVE pendant**, structural-correspondence verdict
(matches Track D scope verbatim).
**No production code modified.**

---

## §1. Executive summary

The Unruh effect — a uniformly accelerated observer in Minkowski vacuum
sees a thermal state at temperature
$T_U = \hbar a / (2 \pi c k_B)$ — is the flat-space analog of the
Schwarzschild Hawking effect, and is the third face of the same
Wick-rotation chain that Track D (2026-05-09) used to land the
Bisognano–Wichmann reading of Sprint TD Track 4. Wick-rotating Rindler
coordinates produces Euclidean polar coordinates with conical singularity
at the origin; smoothness forces angular period $2\pi$, fixing
$\beta_U = 2\pi / a$.

This $2\pi$ is the same M1 Hopf-base measure (Vol($S^1$) = $2\pi$,
master Mellin engine $k=0$, Sprint MR-A/B) that already drove the Hawking
$\beta_{\rm cigar} = 8\pi M = 4M \cdot 2\pi$ identification on the
Schwarzschild Euclidean cigar (Sprint TD Track 4). It is the same $2\pi$
that the Bisognano–Wichmann theorem (1976) identifies, via Sewell 1982,
with the modular-flow / boost-orbit period of an observer with
restricted causal access to the wedge.

**The unified Wick-rotation theorem.** The four-witness chain
[Hartle–Hawking 1976 path-integral derivation] → [Sewell 1982
Schwarzschild generalization] → [Bisognano–Wichmann 1976 Rindler-wedge
modular flow] → [Unruh 1976 flat-space limit] is one Wick-rotation
theorem with four physical instantiations: thermal state of an observer
with restricted causal access is determined by the regularity period of
the Euclidean rotation of their causal-domain coordinates, which equals
$2\pi$ times $1/(\text{surface gravity})$. Surface gravity values:
$\kappa_g = 1/(4M)$ (Schwarzschild horizon), $a$ (Rindler / Unruh),
$2\pi/\beta$ for the boost generator's rapidity rate
(Bisognano–Wichmann).

**Verdict.** Same scope as Track D: **structural correspondence, not
literal identification**. The framework's M1 mechanism on $S^1_\tau$
*IS* the Wick-rotated image of the modular flow on the Rindler wedge,
via the published Wick-rotation prescription. The literal identification
— that the modular Hamiltonian on the truncated metric spectral triple
$\mathcal{T}_{n_{\max}}$ for a Rindler-wedge-restricted state has
analytic period $2\pi$ — is named below as the natural Unruh-side
falsifier and is *the same* operator-system test that Track D named for
Schwarzschild (one falsifier covers both faces; lifting it on either
face lifts the structural correspondence on both).

**Apparatus.** No new code. Sprint TD Track 1's
`thermal_tensor_triple.matsubara_spectrum(β=2π/a)` produces the
Unruh-thermal Matsubara spectrum verbatim, with all three sympy-symbolic
residuals identically zero:
- $T_U$ from $1/\beta_U$ residual: $0$
- Lowest bosonic Matsubara $= a$ (Unruh's Hawking-$\kappa_g$ analog), residual $0$
- Lowest fermionic Matsubara $= a/2 = \pi T_U$, residual $0$

**Top deliverables (this 1–2 day sprint).**

1. This memo (`debug/unruh_pendant_memo.md`, ~3500 words) documenting
   the Unruh pendant.
2. Machine-readable JSON
   (`debug/data/unruh_pendant.json`) with the closed-form $\beta_U =
   2\pi/a$, the M1 identification, and the paper-edit list.
3. Four small paper edits:
   - Paper 35 §VIII new subsection `subsec:unruh_pendant`;
   - Paper 32 §VIII Unruh paragraph appended to
     `rem:bisognano_wichmann_reading` (or new
     `rem:unruh_pendant`, whichever reads cleanly);
   - Paper 34 §V.B Unruh row matching Hawking-T row format
     (machinery-witness, error class C);
   - Paper 34 §III.15 footnote extension naming Unruh as the third
     witness; Paper 34 §VIII Lorentz-boost open question one-sentence
     append.

---

## §2. Wick rotation of Rindler coordinates → $T_U = \hbar a / (2 \pi c k_B)$

### §2.1 Setup

A uniformly accelerated observer in Minkowski space follows a hyperbolic
worldline of constant proper acceleration $a$. The right wedge
$x > |t|$ is the union of all such worldlines; on it, Rindler
coordinates $(\eta, \rho)$ are defined by

$$
t = \rho\, \sinh(a\eta),
\qquad
x = \rho\, \cosh(a\eta),
$$

with metric

$$
\mathrm{d}s^2 \;=\; -\rho^2\, \mathrm{d}\eta^2 \;+\; \mathrm{d}\rho^2
\;+\; \mathrm{d}y^2 \;+\; \mathrm{d}z^2.
$$

The accelerated observer at fixed $\rho = 1/a$ has proper time
$\eta / a$ along the orbit and proper acceleration $a$.

### §2.2 Wick rotation and conical regularity

Wick-rotate the rapidity coordinate $\eta \to i \eta_E$. The wedge
metric becomes Euclidean:

$$
\mathrm{d}s^2_E \;=\; \rho^2\, \mathrm{d}\eta_E^2 \;+\; \mathrm{d}\rho^2
\;+\; \mathrm{d}y^2 \;+\; \mathrm{d}z^2.
$$

This is the metric of $\mathbb{R}^2$ in polar coordinates $(\rho,
\eta_E)$ with $\eta_E$ playing the role of the angular variable, *plus*
the trivial transverse $(y, z)$ direction. Smoothness of the metric at
$\rho = 0$ — the analog of the cigar tip in Sprint TD Track 4 — forces
$\eta_E$ to have period $2\pi$: any smaller or larger period leaves a
conical defect at the origin. (Standard absence-of-conical-singularity
argument, identical in form to the Hawking–Gibbons argument for the
Euclidean Schwarzschild cigar.)

The accelerated observer's *proper time* $\tau_p = \eta_E / a$ at fixed
$\rho = 1/a$ then has period

$$
\beta_U \;=\; \frac{2\pi}{a}.
$$

Restoring units of $\hbar$, $c$, $k_B$:

$$
\beta_U \;=\; \frac{2\pi c}{a},
\qquad
T_U \;=\; \frac{1}{k_B \beta_U} \cdot \hbar
\;=\; \frac{\hbar a}{2 \pi c k_B}.
$$

This is the Unruh temperature.

### §2.3 Translation to KMS

In QFT, a state restricted to the right Rindler wedge of Minkowski
vacuum satisfies the KMS condition with respect to the boost flow
$\sigma_t = e^{itK}$ at inverse temperature $\beta_U = 2\pi/a$ (Unruh
1976; reviewed in Wald 1994 §5.2 and Witten 2018 RMP). The boost
generator $K$ is exactly the modular Hamiltonian that Bisognano–Wichmann
1976 identifies as $-2\pi \cdot M_{\rm boost}$ (in natural units;
restoring $a$ for the proper-acceleration observer gives $K = -(2\pi/a)
\cdot M_{\rm boost}$).

The KMS analytic period of the modular automorphism is therefore
$\sigma_{i \cdot \beta_U} = $ identity, i.e. $\sigma_{i \cdot 2\pi/a} =
1$. Per-rapidity $a = 1$, this is exactly the Bisognano–Wichmann
$\sigma_{i \cdot 2\pi} = 1$.

---

## §3. The $2\pi$ in $\beta_U$ is M1 Hopf-base measure

### §3.1 Trace through the framework

Sprint TD Track 4 identified the $2\pi$ in the cigar's $\phi$-period as
the M1 Hopf-base measure / $\mathrm{Vol}(S^1_\tau)$ signature on the
$\tau$-circle. The Unruh computation traces through identically:

| Step | Object | Value | Mechanism |
|:---|:---|:---|:---|
| 1 | Wedge → Euclidean metric | $\rho^2\, d\eta_E^2 + d\rho^2$ | Wick rotation of $\eta$ |
| 2 | Conical regularity at $\rho = 0$ | $\eta_E \in [0, 2\pi)$ | Standard absence-of-defect argument |
| 3 | Proper-time period | $\beta_U = 2\pi / a$ | $\tau_p = \eta_E / a$, periodicity inherited |
| 4 | Hopf-base measure of $S^1_{\eta_E}$ | $\mathrm{Vol}(S^1) = 2\pi$ | Master Mellin engine M1, $k = 0$ |

The $2\pi$ in $\beta_U$ enters the framework at step 2/3 as the
circumference of the unit angular circle, which is exactly the M1
Hopf-base measure $\mathrm{Vol}(S^1) = 2\pi$ in the master Mellin engine
classification (Paper 18 §III.7; Paper 32 §VIII Remark
`rem:master_mellin_domain`). This is structurally the same factor that
appears in Sprint TD Track 1's $M_1$ column for the Matsubara circle
(Stefan-Boltzmann factorization at $k = 0$) and that Sprint TD Track 4
identified on the Schwarzschild cigar's $\tau$-circle.

### §3.2 Comparison to the Hawking and BW witnesses

| Mechanism | $S^1$ factor | Surface gravity | $\beta$ | $T$ |
|:---|:---|:---|:---|:---|
| Hawking (Track 4) | $\tau$-circle of cigar at horizon | $\kappa_g = 1/(4M)$ | $8\pi M = 2\pi / \kappa_g$ | $T_H = \kappa_g/(2\pi)$ |
| BW (Track D) | modular-flow orbit on Rindler wedge | rapidity rate | $2\pi$ (in units where $\kappa = 1$) | $1/(2\pi)$ |
| Unruh (this sprint) | Wick-rotated $\eta_E$ circle | $a$ | $2\pi / a$ | $T_U = a/(2\pi)$ |

All three present the same M1 mechanism at $k=0$:
$\beta = 2\pi / (\text{surface gravity})$, with $2\pi = \mathrm{Vol}(S^1)$
fixed by Euclidean conical regularity / Hopf-base measure / KMS analytic
period (these are three names for one structural object).

### §3.3 Cross-check against other M1 instantiations

The L2 GH-convergence rate constant $4/\pi = \mathrm{Vol}(S^2)/\pi^2$
on $S^3$ (Paper 38; Sprint MR-B) is a *different* M1 instantiation: the
Hopf-base measure there is the spatial $\mathrm{Vol}(S^2) = 4\pi$,
divided by $\pi^2 = \mathrm{Vol}(S^3)/2$, giving $4/\pi$. The Unruh
case has Hopf-base measure $\mathrm{Vol}(S^1) = 2\pi$ on a temporal
circle. These are different sub-manifolds, hence different
*instantiations* of the same M1 sub-mechanism, both within the M1
transcendental ring (rational multiples of $\pi$), exactly as the
master Mellin engine domain partition predicts (Paper 32 §VIII
Remark `rem:master_mellin_domain`).

A useful reference table:

| M1 instantiation | Sub-manifold | Hopf-base measure |
|:---|:---|:---|
| Unruh | Rindler $\eta_E$-circle | $\mathrm{Vol}(S^1) = 2\pi$ |
| Hawking | Schwarzschild $\tau$-circle | $\mathrm{Vol}(S^1) = 2\pi$ (× $4M$ for $\beta$) |
| L2 propinquity rate (Paper 38) | Spatial $S^3 \to S^2$ Hopf base | $\mathrm{Vol}(S^2)/\pi^2 = 4/\pi$ |
| Hopf-bundle measure (Paper 2) | $S^3$ as $S^1 \to S^3 \to S^2$ | $\mathrm{Vol}(S^3)/\mathrm{Vol}(S^1) = \pi$ |

All four have ring class M1 ($\mathbb{Q}\cdot\pi^{n}$ for small integer
$n$). This is consistent with the case-exhaustion theorem (Paper 32
§VIII): the *transcendental ring* M1 is universal; the *specific
$\pi^n$ exponent* is determined by the dimension of the Hopf-base
sub-manifold being integrated against.

---

## §4. The unified Hawking–Sewell–BW–Unruh Wick-rotation theorem

### §4.1 Statement

Let $\mathcal{O}$ be an observer with restricted causal access to a
sub-region $\mathcal{R}$ of a Lorentzian spacetime, where $\mathcal{R}$
is bounded by a bifurcate Killing horizon (Schwarzschild exterior,
Rindler right wedge, de Sitter static patch, etc.) with surface gravity
$\kappa$. The thermal state of $\mathcal{O}$ with respect to the
Killing-time evolution is KMS at inverse temperature

$$
\beta \;=\; \frac{2\pi}{\kappa}.
$$

Equivalently, the Wick rotation of the Killing time $t \to i\tau$
produces a Euclidean geometry whose smoothness at the horizon (no
conical defect) forces the $\tau$-circle period to be $2\pi/\kappa$.

The single $2\pi$ in this formula is, in all four witnesses,
$\mathrm{Vol}(S^1)$ — the regularity period of the Euclidean rotation
of the observer's causal-domain coordinates.

### §4.2 The four witnesses

| Witness | Year | $\mathcal{R}$ | $\kappa$ | $\beta$ | $T$ |
|:---|:---:|:---|:---|:---|:---|
| Hartle–Hawking | 1976 | Schwarzschild exterior | $1/(4M)$ | $8\pi M$ | $1/(8\pi M)$ |
| Sewell | 1982 | bifurcate Killing horizon (general) | (Killing) | $2\pi / \kappa_{\rm Killing}$ | $\kappa_{\rm Killing}/(2\pi)$ |
| Bisognano–Wichmann | 1976 | Rindler right wedge of Minkowski | (boost rapidity rate) | $2\pi$ (rescaled units) | $1/(2\pi)$ |
| Unruh | 1976 | Rindler wedge (proper acceleration $a$) | $a$ | $2\pi/a$ | $a/(2\pi)$ |

The four witnesses are not four independent results; they are four
faces of the same Wick-rotation theorem. Hartle–Hawking gives the
Schwarzschild instantiation as a path-integral derivation;
Bisognano–Wichmann gives the Rindler instantiation as a
modular-automorphism statement; Sewell generalizes to any bifurcate
Killing horizon; Unruh gives the flat-space limit as a thermal state
of a detector. **In each case, the $2\pi$ is the regularity period of
the Euclidean angular coordinate** — i.e., the M1 Hopf-base measure
$\mathrm{Vol}(S^1)$ of the master Mellin engine at $k = 0$.

### §4.3 Why the framework's M1 mechanism captures all four

The framework reproduces each instantiation by applying Sprint TD
Track 1's `matsubara_spectrum(β)` machinery with the appropriate
$\beta$:

- **Stefan-Boltzmann (Track 1):** $\beta = 1/T$, no horizon; M1 gives
  the partition-function $\mathrm{Vol}(S^1_\beta)$ factor and combines
  with M2 ($\zeta_R(4)$) to reproduce $\sigma = \pi^2/90$.
- **Hawking (Track 4):** $\beta = 8\pi M = 2\pi/\kappa_g$;
  `matsubara_spectrum(β)` gives lowest bosonic mode $= \kappa_g$ (sympy
  residual zero); M1 reproduces $T_H = \kappa_g/(2\pi)$.
- **Bisognano–Wichmann (Track D):** $\beta = 2\pi$ (in units where the
  rapidity rate is 1); the framework reads this off the M1 mechanism
  via the published Wick-rotation chain Sewell 1982 + Hartle–Hawking
  1976; structural correspondence verified.
- **Unruh (this sprint):** $\beta = 2\pi/a$;
  `matsubara_spectrum(2π/a)` gives lowest bosonic mode $= a$ (sympy
  residual zero, manual check §5); M1 reproduces $T_U = a/(2\pi)$.

The framework's contribution is *not* deriving any of the four facts
(all four are textbook physics), but rather identifying the $2\pi$ in
each as the M1 Hopf-base-measure signature of the master Mellin engine,
which embeds the four Wick-rotation faces inside Paper 32 §VIII's
case-exhaustion theorem and inside Paper 38's GH-convergence
infrastructure (where the same M1 mechanism with the spatial $S^3 \to
S^2$ Hopf base produces the propinquity rate $4/\pi$).

---

## §5. Apparatus check (sympy verification, no new code)

Sprint TD Track 1's `geovac/thermal_tensor_triple.matsubara_spectrum`
applies verbatim with $\beta = 2\pi/a$ in place of $\beta = 8\pi M$
(Hawking) or $\beta = 1/(k_B T)$ (Stefan-Boltzmann). The check
performed at this sprint (no production-code modification):

```python
import sympy as sp
from geovac.thermal_tensor_triple import matsubara_spectrum

a = sp.Symbol('a', positive=True)
beta_U = 2 * sp.pi / a
T_U = a / (2 * sp.pi)

bos = matsubara_spectrum(beta_U, k_max=4, fermionic=False)
fer = matsubara_spectrum(beta_U, k_max=4, fermionic=True)

# Lowest non-zero bosonic = 2π/β_U = a (Unruh's surface-gravity analog)
omega_b1 = sp.simplify([x for k, x in bos if k > 0][0])
assert sp.simplify(omega_b1 - a) == 0          # PASS

# Lowest fermionic = π/β_U = a/2 = π·T_U
omega_f0 = sp.simplify([x for k, x in fer if k >= 0][0])
assert sp.simplify(omega_f0 - sp.pi * T_U) == 0   # PASS
assert sp.simplify(omega_f0 - a/2) == 0           # PASS

# T_U via Tomita-Takesaki KMS period 2π/a
assert sp.simplify(1/beta_U - T_U) == 0           # PASS
```

**All four sympy-symbolic identities residual zero**, identical in
form to Sprint TD Track 4's three Hawking residuals. The same Matsubara
apparatus that produces Stefan-Boltzmann (Track 1) and Hawking
(Track 4) produces the Unruh-thermal trace (this sprint), with the only
difference being the value of $\beta$ supplied.

The lowest bosonic mode equals the surface-gravity analog ($a$ in
Unruh, $\kappa_g$ in Hawking, $2\pi/\beta$ in Stefan-Boltzmann) in all
three cases — this is the M1 mechanism's universal signature on $S^1$
factors.

---

## §6. Verdict — structural correspondence, not literal identification

**Same scope as Track D, verbatim.**

The framework's M1 mechanism on $S^1_\tau$ produces $\beta_U = 2\pi/a$
on the Wick-rotated Rindler wedge, with the $2\pi$ being the
$\mathrm{Vol}(S^1)$ Hopf-base-measure factor of the master Mellin
engine at $k = 0$. Under the standard published Wick-rotation chain
(Hartle–Hawking 1976 → Sewell 1982 → Bisognano–Wichmann 1976 → Unruh
1976), this $2\pi$ IS the period of the modular-flow / boost orbit of
an observer with restricted causal access to the wedge. **The
framework's M1 mechanism is therefore the Wick-rotated image of the
Bisognano–Wichmann modular-flow mechanism, instantiated at the Rindler
wedge for the Unruh case.**

This is a **structural correspondence**, not a literal identification.
Two distinct mathematical objects (a Hopf-bundle measure factor on the
Riemannian side; a modular-automorphism period on the Lorentzian side)
are equated by a published Wick-rotation prescription, not by an
operator-level proof inside the GeoVac spectral triple. Track D's
scope statement applies verbatim to the Unruh face: the same operator-
level falsifier (compute the modular Hamiltonian on
$\mathcal{T}_{n_{\max}}$ for a wedge-restricted state and verify
$\sigma_{i \cdot 2\pi} = 1$) tests both faces simultaneously; lifting
either face to literal identification lifts both.

---

## §7. Paper-edit summary

Four small edits, applied directly per CLAUDE.md §13.8. Voice and
structural-correspondence framing matched verbatim to the
Bisognano–Wichmann landing.

### §7.1 Paper 35 §VIII new subsection `subsec:unruh_pendant`

Added immediately after `subsec:bisognano_wichmann_reading`, before
`\section{Open questions}`. Mirrors the structure of the
Bisognano–Wichmann subsection: Wick rotation of Rindler, conical
regularity at $\rho = 0$ forces $\eta_E$-period $2\pi$,
$\beta_U = 2\pi/a$, lift via Track 1 `matsubara_spectrum(β=2π/a)`,
identification with the boost-orbit period via the same Wick-rotation
chain. Honest scope: same as BW (structural correspondence, not literal
identification). Cross-references: Paper 32 §VIII Remark
`rem:bisognano_wichmann_reading`, Paper 34 §III.15 footnote, Paper 34
§V.B Unruh row.

### §7.2 Paper 32 §VIII Remark `rem:bisognano_wichmann_reading` extension

Added a paragraph at the end of the existing remark naming the Unruh
pendant explicitly: the same M1 mechanism applied with $\beta = 2\pi/a$
reproduces $T_U = a/(2\pi)$ as the Wick-rotated image of the
modular flow on the Rindler wedge. Sympy verification at $\beta_U =
2\pi/a$ residual zero. Same structural-correspondence scope. The
unified Hawking–Sewell–BW–Unruh four-witness reading of M1 on temporal
compactifications is now made explicit. The named falsifier
($\sigma_{i \cdot 2\pi} = 1$ on $\mathcal{T}_{n_{\max}}$) covers both
faces.

### §7.3 Paper 34 §V.B Unruh row

New row in the off-precision catalogue table immediately after the
Hawking row (lines ~1994–2017 of `paper_34_projection_taxonomy.tex`),
matching the Hawking row format:

> Unruh temperature $T_U = a/(2\pi)$ on Rindler wedge & Fock $\circ$
> observation/temporal-window (§ref{sec:proj_observation}) at $\beta =
> 2\pi/a$ & machinery-witness only & C & Unruh pendant
> (`debug/unruh_pendant_memo.md`, May 2026) reproduces $\beta_U =
> 2\pi/a$ from the framework's M1 Hopf-base measure /
> $\mathrm{Vol}(S^1_{\eta_E})$ mechanism alone (sympy residual zero on
> $T_U$, lowest bosonic Matsubara $= a$, lowest fermionic $= a/2 = \pi
> T_U$). $a$ is external proper-acceleration input, not a framework
> calibration mismatch in the usual sense; row is a machinery-witness.
> Same structural-correspondence scope as the Hawking row: under the
> Wick-rotation chain Hartle–Hawking 1976, Sewell 1982,
> Bisognano–Wichmann 1976, the framework's M1 $2\pi$ IS the
> modular-flow / boost-orbit period of the accelerated observer. See
> Paper 32 §VIII `rem:bisognano_wichmann_reading` (Unruh paragraph) and
> Paper 35 §VIII `subsec:unruh_pendant`.

### §7.4 Paper 34 §III.15 footnote extension

The existing Bisognano–Wichmann footnote in §III.15 (lines ~520–545)
gains one sentence at the end naming Unruh as the third witness:

> The Unruh corollary $T_U = a/(2\pi)$ on a Rindler wedge is the
> flat-space limit of the same Wick-rotation chain and is documented in
> the Unruh pendant memo (`debug/unruh_pendant_memo.md`, May 2026); it
> closes the Track D §5.4 follow-up flag.

### §7.5 Paper 34 §VIII Lorentz-boost open question append

The Track D documentation closure paragraph (lines ~4165–4201) gains
one sentence noting that the Unruh pendant lands as a fourth
Wick-rotation witness without changing the REQUIRES-EXTENSION verdict
on literal lifting:

> The Unruh pendant (`debug/unruh_pendant_memo.md`, May 2026) added the
> flat-space face $T_U = a/(2\pi)$, completing the Hawking–Sewell–BW–Unruh
> four-witness chain at the framework's M1 sub-mechanism level; the
> REQUIRES-EXTENSION verdict on a literal Lorentz-boost projection is
> unchanged, since the same operator-level falsifier
> ($\sigma_{i\cdot 2\pi} = $ identity on $\mathcal{T}_{n_{\max}}$)
> tests all four faces simultaneously.

---

## §8. Named falsifier — modular Hamiltonian on the Rindler wedge in $\mathcal{O}_{n_{\max}}$

The natural Unruh-side falsifier is structurally the same as Track D's
falsifier #1, instantiated at the Rindler wedge instead of the
Schwarzschild static exterior. **Falsifier name: Unruh modular-flow
period closure on $\mathcal{T}_{n_{\max}}$.**

**Claim being falsified:** The structural correspondence of §6 lifts to
literal identification at finite $n_{\max}$. That is, the modular
Hamiltonian $K$ on the truncated operator system $O_{n_{\max}}$ (Paper
32 §III) for a Rindler-wedge-restricted Camporesi–Higuchi state has
imaginary-axis period $2\pi$.

**Why this is *the same* falsifier as Track D's #1 (one falsifier
covers two faces).** Bisognano–Wichmann is the wedge case;
Schwarzschild thermal is the bifurcate-horizon case; Unruh is the
flat-space-limit / proper-acceleration case. All three are KMS periods
$\sigma_{i\cdot 2\pi/\kappa} = 1$ for the appropriate $\kappa$.
Verifying the period closure on $\mathcal{T}_{n_{\max}}$ for *any one*
wedge-restricted state lifts the structural correspondence on *all
three* faces, because the Wick-rotation chain that ties the framework
$2\pi$ to the modular-flow $2\pi$ runs through the same KMS-Tomita-Takesaki
algebra in every case.

**Operationalization (same as Track D §3.1, instantiated for Rindler).**
Extend the existing Connes-distance SDP framework (`geovac/connes_distance.py`)
to compute the modular automorphism on $O_{n_{\max}}$ for a half-$S^3$
sub-algebra (the Riemannian-side analog of the Rindler wedge after
spatial reduction). Compute the polar decomposition $S = J \Delta^{1/2}$
of $S: a\Omega \to a^*\Omega$ on the GNS Hilbert space; verify
$\Delta^{2\pi i} = 1$ at $n_{\max} = 2, 3, 4$.

**Estimated cost:** 4–8 weeks (same as Track D #1, since the operator
system construction is shared between the BW and Unruh faces).

**Priority:** MEDIUM-HIGH. Lifts the structural correspondence on all
three Wick-rotation faces (Hawking, BW, Unruh) simultaneously to literal
identification at finite $n_{\max}$, in the same framework where WH1
PROVEN (Paper 38 GH-convergence) was established.

**If it fails:** all three Wick-rotation-face structural correspondences
are simultaneously refuted at the operator-system level. The framework's
KMS structure on the truncated triple is then incompatible with the
Wightman-axiomatic vacuum's KMS structure, and a rederivation of the
modular flow without assuming Wightman axioms would be needed.

**If it passes:** the M1 Hopf-base measure / boost-orbit-period
identification is lifted from "structural correspondence via published
Wick-rotation" to "operator-level theorem inside the GeoVac framework
at finite $n_{\max}$, plus GH limit via Paper 38's lemmas" — making the
framework an internal Bisognano–Wichmann theorem prover at finite cutoff.

---

## §9. Files

- This memo: `debug/unruh_pendant_memo.md`.
- JSON: `debug/data/unruh_pendant.json` (closed-form $\beta_U$, M1
  identification, paper-edit list as machine-readable rows).
- Paper edits: `papers/observations/paper_35_time_as_projection.tex`,
  `papers/synthesis/paper_32_spectral_triple.tex`,
  `papers/observations/paper_34_projection_taxonomy.tex`.
- No production code modified; `geovac/thermal_tensor_triple.py` was
  cited via manual sympy check (§5) but not extended.

---

## §10. Bibliography (Unruh-specific, in addition to the BW landing's
verified citations)

1. **Unruh, W.G.** "Notes on black-hole evaporation."
   *Phys. Rev. D* **14**, 870–892 (1976). The original derivation of
   the Unruh effect. Standard citation.

The four other Wick-rotation-witness citations (Bisognano–Wichmann
1975/1976, Sewell 1982, Hartle–Hawking 1976) were verified in the
Track D bibliography
(`debug/bisognano_wichmann_track_d_memo.md` §10) and are already
present in Paper 32 / Paper 34 / Paper 35 bibliographies. No new
bibitems are added by this sprint beyond `unruh1976`.
