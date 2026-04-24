# Sprint A Memo — Sign Pattern Check on K = π(B + F − Δ)

**Author:** Sub-agent (PM dispatch, α-program reframe Sprint A)
**Scope:** WH1 / WH5 (CLAUDE.md §1.7); structural check on the EXISTING
decomposition. NO new common-generator proposal.
**Date:** 2026-04-18

---

## 1. Summary verdict

**NO MATCH (with one suggestive partial along the eta-invariant axis).**

The standard Connes–Chamseddine spectral-action expansion
$\mathrm{Tr}\,f(D/\Lambda) \;\sim\; \sum_k f_k\,\Lambda^{d-k}\,a_k$
on a closed Riemannian manifold puts every Seeley–DeWitt coefficient
$a_k$ in front of a positive moment $f_k = \int_0^\infty f(u)\,u^{k-1}\,du$
of a *positive* cutoff function $f$ (equivalently, $f$ a positive
Schwartz function). The bulk expansion is therefore a **monotone all-positive
sum**: there is no structural rule in the Connes–Chamseddine framework that
forces alternation of sign across $a_0, a_2, a_4$.

The **sign pattern $(+, +, -)$** of Paper 2's $K = \pi(B + F - \Delta)$
does *not* match the bulk Seeley–DeWitt rule. The single rule in the
broader heat-kernel / index-theory literature that *does* generate a
sign-flipped boundary correction sitting next to bulk coefficients is the
**Atiyah–Patodi–Singer (APS) eta-invariant correction**, and we identify
that as the *only* structural channel under which a $(+, +, -)$ pattern
could be a matched feature rather than ad hoc.

The match-strength is **weak/suggestive**. The map
$B \leftrightarrow a_2$-type bulk Casimir trace, $F \leftrightarrow$
zeta-regularized fiber spectral invariant, $\Delta^{-1} = g_3^\text{Dirac}
\leftrightarrow$ APS-type single-level boundary mode count, is internally
consistent in shape but is *not* derivable from any one functional.

---

## 2. Heat-kernel sign rules surveyed

Sources consulted: Vassilevich's heat-kernel manual (arXiv hep-th/0306138);
Connes & Chamseddine "The Spectral Action Principle" (CMP 186 (1997) 731,
arXiv hep-th/9606001); Chamseddine & Connes "The Uncanny Precision of the
Spectral Action" (CMP 293 (2010) 867); van Suijlekom *Noncommutative Geometry
and Particle Physics* (2nd ed., 2024); Estrada–Gracia-Bondía–Várilly
distributional approach (arXiv 1101.4804 Renormalization of the Spectral
Action). All citations to specific equation numbers below are paraphrases
sufficient for the structural check; the project memo is not a literature
review.

### Rule (R1) — Connes–Chamseddine bulk expansion (cosmological/EH/curvature^2)

$\mathrm{Tr}\,f(D/\Lambda) \simeq f_4\,\Lambda^4\,a_0 + f_2\,\Lambda^2\,a_2
+ f_0\,a_4 + O(\Lambda^{-2})$

with $f_k := \int_0^\infty f(u)\,u^{k-1}\,du$ for $k>0$ and $f_0 := f(0)$.
For the standard cutoff function (positive even Schwartz function, e.g.
the Gaussian $\exp(-x^2)$), every $f_k > 0$, every bulk $a_k$ is a
real spectral invariant whose sign is determined by curvature
(positive on $S^d$). On unit $S^3$:

| Term | $a_k$ | Sign |
|------|-------|------|
| Cosmological | $a_0 = (4\pi)^{-3/2} \cdot \dim_S \cdot \mathrm{Vol}(S^3) = (4\pi)^{-3/2}\cdot 4 \cdot 2\pi^2$ | $+$ |
| Einstein–Hilbert | $a_2 \propto \int R/6 = +$ | $+$ |
| Curvature$^2$ | $a_4$ on a constant-curvature space form: combination of $5R^2 - 2|\mathrm{Ric}|^2 + 2|\mathrm{Riem}|^2 - 30 R E + 60 E^2$, **net positive on $S^3$** with Lichnerowicz $E = R/4$ | $+$ |

Conclusion (R1): all three bulk coefficients enter with the SAME sign on $S^3$.

### Rule (R2) — Bochner–Weitzenböck / Lichnerowicz curvature correction

For the squared Dirac operator $D^2 = \nabla^*\nabla + R/4$, the Lichnerowicz
endomorphism $E = R/4$ enters $a_4$ with definite sign. On unit $S^3$:
$E = 3/2 > 0$. There is no Bochner sign-flip relative to the bulk.

### Rule (R3) — APS eta-invariant boundary correction

When the manifold has boundary, or when one works on a finite-cutoff
truncation (Estrada–Gracia-Bondía distributional setting), the **APS
correction** appears. The Atiyah–Patodi–Singer index theorem on a manifold
$M$ with boundary $\partial M$ has the structure

$\mathrm{ind}(D) = \int_M \widehat{A} - \frac{\eta(D_{\partial M}) + h}{2}$

where $\widehat{A}$ is the bulk Pontryagin/$\widehat{A}$-roof contribution
and the boundary term is **subtracted**. In a finite-spectral-truncation
setting, the analog is a "boundary mode count" living on the cutoff edge.
This is the only standard heat-kernel / index-theory rule we surveyed in
which a finite-mode degeneracy enters with the **opposite sign** to the
bulk contributions.

### Rule (R4) — Estrada–Gracia-Bondía moment asymptotic expansion

The distributional MAE for $\mathrm{Tr}\,f(D)$ has structure
$\sum_n c_n M_n[f]/n!$ with $c_n$ Cesàro coefficients of the spectral
density. Signs of $c_n$ are not constrained to be positive in general; the
distributional approach permits negative $c_n$ when the spectral density
has subtractions. However, this generality does not by itself produce a
specific $(+, +, -)$ pattern — it just relaxes R1.

### Rule (R5) — Spectral truncation à la van Suijlekom

The "spectral truncations" framework (van Suijlekom 2nd ed. Ch. 13;
Connes–van Suijlekom CMP 2020) replaces $\mathrm{Tr}\,f(D)$ by the trace on
the finite operator system cut off at eigenvalues $|\lambda| \le N$. This
introduces a **Tolerance term** that is genuinely a finite-cutoff
correction, but its sign is not pinned positive or negative — it depends on
the specific cutoff and is treated as a small remainder, not as a separate
additive object competing with bulk $a_k$.

---

## 3. Mapping B, F, Δ to operator sectors under each rule

### Under (R1) bulk Connes–Chamseddine

| Paper-2 piece | Matches a$_k$? | Comment |
|---|---|---|
| $B = 42$ | shape of a finite-truncation Casimir trace; closer to a *truncated $a_2$ with a Casimir weight* than to any specific $a_k$ | NOT a heat-kernel SD coefficient; it is a finite Peter–Weyl Casimir trace at $m=3$. Lives on the scalar Laplace–Beltrami sector. |
| $F = \pi^2/6 = \zeta_R(2)$ | shape of a $\zeta$-regularized invariant | Spectral-zeta-at-2 of fiber $S^1$ (or, per Phase 4F, $D_{n^2}(d_\text{max}) = \zeta_R(2)$ on the scalar $S^3$ Fock degeneracy). Lives in the spectral-zeta family; structurally akin to a regularized-determinant or $\zeta'(0)$ contribution, NOT a Seeley–DeWitt $a_k$. |
| $\Delta = 1/40$ | none under R1 | $\Delta^{-1} = g_3^\text{Dirac} = 40$ is a single-level Dirac mode degeneracy, not a curvature integral. |

Mapping under R1: B is *truncated bulk*, F is *fiber zeta-regularization*,
Δ is *boundary mode count*. **All three are different object types**; this
already echoes Paper 2 §IV.A "Theorem (Three homes, structural mystery)".
R1 does not pin the signs; it only says bulk $a_k$ contributions are
all-positive — which $B$ and $F$ both respect (both enter with $+$).

### Under (R3) APS / boundary correction

The APS rule says: "subtract the boundary mode contribution". With
$\Delta^{-1} = g_3^\text{Dirac} = 40$ identified as the Dirac single-level
mode count at the cutoff, the sign $-\Delta$ in $K = \pi(B + F - \Delta)$
matches the APS sign convention for the eta-invariant boundary correction
in spirit:
- Bulk pieces ($B$, $F$) enter with $+$.
- Boundary/cutoff piece ($\Delta$) enters with $-$.

This is the only rule under which $(+, +, -)$ is structurally motivated
rather than ad hoc. The match is **shape only, not derivational**: APS
applies to manifold-with-boundary or to spectral-flow on a one-parameter
family, neither of which is literally what Paper 2 has — Paper 2 has a
finite Fock cutoff $n \le 3$ on a closed manifold, which the truncation
framework (R5) handles, but R5 does not pin the sign.

### Under (R4) Estrada distributional MAE

Sign-permissive but not predictive. Cannot adjudicate.

### Under (R5) van Suijlekom spectral truncation

Sign-permissive but not predictive. Cannot adjudicate.

---

## 4. Fit / mismatch discussion

**What MATCHES (suggestive partial):**
1. Three structurally distinct object types in the formula: this is exactly
   the Estrada / van Suijlekom expectation when one mixes a curvature
   integral (bulk), a zeta-regularized invariant (fiber), and a finite
   cutoff term (boundary). Paper 2's "three homes" theorem is the same
   observation phrased without the spectral-action language.
2. The minus sign on $\Delta$ is consistent with APS-style boundary
   subtraction. No other heat-kernel / spectral-action rule we surveyed
   would put a $(-)$ in front of a finite mode degeneracy.
3. The leading $\pi$ overall prefactor is consistent with R1's
   $(4\pi)^{-d/2}$ measure factor on a $d=3$ manifold (modulo coefficient).

**What MISMATCHES (decisive structural):**
1. Paper 2 has **no cutoff function** $f$. There are no moments $f_k$. The
   formula $K = \pi(B + F - \Delta)$ is a single number with single weights
   $(+1, +1, -1)$, not the leading three terms of an expansion in $\Lambda$.
   The shape-match to the spectral-action expansion is at the level of "three
   pieces of distinct type" only.
2. The $(+, +, -)$ sign pattern is NOT the SD pattern $(a_0, a_2, a_4)$
   on $S^3$, which is $(+, +, +)$.
3. The three Paper 2 pieces live in **different operator sectors** (scalar
   Laplace–Beltrami for $B$ and $F$; spinor for $\Delta$) and at **different
   degrees of regularization** (finite trace, zeta-regularized infinite sum,
   single-level degeneracy). Connes–Chamseddine bulk expansion lives in a
   single sector at a single degree. A genuine spectral-action expansion
   would not mix these.
4. Phases 4B–4I (CLAUDE.md §3) eliminated *nine* common-generator
   mechanisms including the closest spectral-action analogs (Phase 4E α-I:
   S$^5$ spectral geometry produces $\zeta_R(\text{odd}) + \log 2$ in
   determinants and pure rationals in SD coefficients, but **never** the
   Paper 2 $(+, +, -)$ pattern). The α-I track explicitly checked SD
   coefficients $a_0, a_1, a_2$ on $S^3$ and $S^5$; all are positive,
   none reproduce $\Delta = 1/40$ or its negative.
5. The APS match (R3) is **shape only**: the actual APS eta-invariant of
   the Camporesi–Higuchi Dirac operator on $S^3$ is well-defined and is not
   $1/40$; the eta-invariant of the APS-style boundary operator at a finite
   Fock cutoff $n \le 3$ would be a sum of half-integer eigenvalues, not a
   reciprocal mode count. So $\Delta$ is not literally an APS eta-invariant.

**Verdict refined:** the sign pattern $(+, +, -)$ is **shape-compatible**
with an APS-style heat-kernel decomposition into (bulk curvature integral)
$+$ (regularized zeta) $-$ (boundary/cutoff correction) but is **not
derivable** from any specific spectral-action functional on $S^3$, and the
three pieces are not the three leading SD coefficients of any single
expansion. Phases 4B–4I have already eliminated the closest candidate
mechanisms.

---

## 5. Recommendation for Paper 2 reframe §III

Two compatible recommendations, neither requiring modification of the
existing Theorem (Three homes, structural mystery) in §IV.A:

### (A) Add an honest one-paragraph framing remark to §III

Replace nothing in §III's spectral-invariant derivations; simply add a new
final paragraph noting:

> The combination rule $K = \pi(B + F - \Delta)$ has the **shape** of an
> APS-style spectral decomposition into a finite curvature/Casimir bulk
> trace (sign $+$), a zeta-regularized fiber invariant (sign $+$), and a
> single-level boundary mode correction (sign $-$). However, no
> spectral-action functional on $S^3$ in the Connes–Chamseddine framework
> generates this combination as its leading three terms: in the standard
> expansion $\mathrm{Tr}\,f(D/\Lambda) \sim \sum_k f_k\,\Lambda^{d-k}\,a_k$
> the bulk Seeley–DeWitt coefficients $a_0, a_2, a_4$ on the round $S^3$
> are all positive, and the cutoff moments $f_k$ are positive for any
> positive Schwartz cutoff. The $(+, +, -)$ sign pattern in $K$ is
> **shape-compatible** with the APS eta-invariant boundary subtraction but
> is not derivable from APS as written, since $\Delta$ is a Dirac
> single-level degeneracy rather than the eta-invariant of any boundary
> operator. We record this as a structural observation; it sharpens the
> Phase 4I "three homes" theorem by attaching the APS-shape interpretation
> to the minus sign on $\Delta$.

### (B) Optional: add an "open question" closing remark to §IV

> Open question: is there a finite-resolution spectral truncation in the
> sense of Connes–van Suijlekom (CMP 2020) on the $S^3$ Hopf bundle whose
> "tolerance" remainder term reproduces $-\Delta = -1/40$ when truncated at
> $n = 3$? The Phase 4 program did not test this specific truncation
> functional. This is a candidate Sprint B target if the α-program reframe
> is reopened.

Both edits are **conservative** and consistent with the WH1/WH5 framing in
CLAUDE.md §1.7: WH5 says α is best read as a projection constant between
three regimes; WH1 says GeoVac has the *shape* of an almost-commutative
spectral triple. Section A above is the operational consequence: name the
shape-match honestly, attribute the $-$ sign to APS-flavor boundary
correction *qualitatively*, do not claim derivation.

---

## 6. References (consulted)

1. A.\ Connes, A.H.\ Chamseddine, "The Spectral Action Principle,"
   CMP 186 (1997) 731, [hep-th/9606001](https://arxiv.org/abs/hep-th/9606001).
2. A.H.\ Chamseddine, A.\ Connes, "The Uncanny Precision of the Spectral
   Action," CMP 293 (2010) 867.
3. D.V.\ Vassilevich, "Heat kernel expansion: user's manual,"
   Phys. Rep. 388 (2003) 279,
   [hep-th/0306138](https://arxiv.org/abs/hep-th/0306138).
4. W.D.\ van Suijlekom, *Noncommutative Geometry and Particle Physics*,
   2nd ed., Mathematical Physics Studies, Springer (2024).
5. A.\ Connes, W.D.\ van Suijlekom, "Spectral Truncations in Noncommutative
   Geometry and Operator Systems," CMP (2020),
   [arXiv:2004.14115](https://arxiv.org/abs/2004.14115).
6. A.H.\ Chamseddine, A.\ Connes, "Renormalization of the Spectral
   Action," [arXiv:1101.4804](https://arxiv.org/abs/1101.4804).
7. M.F.\ Atiyah, V.K.\ Patodi, I.M.\ Singer, "Spectral asymmetry and
   Riemannian geometry I," Math. Proc. Camb. Phil. Soc. 77 (1975) 43.
8. GeoVac Phase 4E Track α-I:
   `debug/data/track_alpha_phase4e/track_i_analysis.md`.
9. GeoVac qed_vacuum_polarization.py: explicit SD coefficients
   $a_0, a_1, a_2$ on $S^3$, all positive.
10. GeoVac Paper 28 §4: $D_\mathrm{even}(s) - D_\mathrm{odd}(s) = 2^{s-1}\,
    (\beta(s) - \beta(s-2))$ — the only place in the framework where a
    *signed* Dirichlet difference is structural.

---

## 7. Honest scope statement

This memo is a **structural shape-check**, not a derivation. It does *not*
propose a new common-generator mechanism (would violate the sprint
constraint). It does *not* modify any existing files. The headline finding
— APS-shape-compatibility of the minus sign on $\Delta$ — is suggestive
but **weak evidence**, on par with the Phase 4 partial identifications
($\kappa \leftrightarrow B$, $F = D_{n^2}(d_\text{max})$, $\Delta^{-1} =
g_3^\text{Dirac}$): each is a clean structural observation, none is a
derivation.

The match strength on the Sprint A question itself is therefore:
**WEAK PARTIAL, APS-shape only** — sufficient to motivate Recommendation A
(one paragraph in Paper 2 §III), insufficient to motivate any change to the
"conjectural" status of Eq. (K).

