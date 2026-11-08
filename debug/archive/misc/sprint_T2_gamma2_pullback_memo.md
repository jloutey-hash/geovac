# Sprint memo — T2 as a Γ(2) length-2 wt-3 iterated-Eisenstein integral? (PI-directed structural test)
Date: 2026-08-20 | Branch: work/sparsity-boundary (uncommitted)
Owning paper: papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex (sec:modular, sec:bessel_algebra)
Drivers: debug/routeC_gamma2_eisenstein_basis.py, debug/routeC_gamma2_pullback.py
Data: debug/data/gamma2_eisenstein_basis.json

## Task
Determine whether the integrated collinear observable T2 = 0.3953557659017139641…(~19 dig,
frozen) — established weight-3, shown to REQUIRE G=Catalan — is a length-2 weight-3 Γ(2)
iterated-Eisenstein integral / Broadhurst–Dorigoni resurgent Lambert series. Build the basis,
test the pullback. DIGIT-INDEPENDENT structural work.

## The object (exact)
Collinear geometry X=0, Y=(0,0,1), Z=(0,0,−1), 1s ζ=1 ⇒ D₁=D₂=1, |W|=s+t.
  T2 = (8/π) ∫₀¹ds ∫₀¹dt J(s,t),  J symmetric,
  J(s,t) = ∫₀^∞ dk j₀(k(s+t)) P(s,k) P(t,k),
  P(x,k) = c e^{−Δ}(Δ⁻³+3Δ⁻⁴+3Δ⁻⁵),  c=x(1−x),  Δ=√(ck²+1).
The fibre curve is y²=(x²−1)(ρx²+1−ρ), modulus ρ=c_t/c_s, modular via λ(τ)=1−ρ (X(2)/Legendre).

## VERDICT — BORDERLINE (rigorous determination; NO finite candidate)
**T2 is NOT a finite length-2 weight-3 Γ(2) iterated-Eisenstein integral. It is a Γ(2)
RESURGENT LAMBERT SERIES (Broadhurst–Dorigoni class). The finite length-2 weight-3 Γ(2) MMV
in the ring {π, ϖ, 1/ϖ, G} is the regular D→0 SHADOW of T2, not T2 itself.**
This CONFIRMS the paper's/memo's hand-off (the twist is obstructive), now with a proof of the
discriminator rather than an assertion. The residual open question is whether the physical
(D=1-twisted) value nonetheless COLLAPSES to a finite ring element — undecided at ~19 dig, needs
~32 dig OR the explicit BD Eichler summation (specialist step).

## Deliverable A — the basis (built, verified)
Driver routeC_gamma2_eisenstein_basis.py. All classical identities checked at TWO precisions
(dps 30/45) with independent evaluators (theta-series vs mpmath ellipk/ellipe/quad).

- **Weight-2 Eisenstein space:** M₂(Γ(2)) = span{θ₂⁴, θ₄⁴}, **dim 2**, Jacobi θ₃⁴=θ₂⁴+θ₄⁴
  (verified |·|≤4e-16). S₂(Γ(2))=0 (genus 0) ⇒ ALL weight-2 forms are Eisenstein. These are the
  weight-2 Eisenstein generators the two Feynman integrations pull back against.
- **Length-1 Eichler dictionary → {π, ϖ, 1/ϖ, G}** (all err ≤ 2e-31 at dps30, ≤2e-46 at dps45):
    · ∫₀¹ K(k) dk = 2G          (Catalan G = L(2,χ₋₄), the weight-2 Eisenstein L-value — NATIVE)
    · ∫₀¹ E(k) dk = G + 1/2
    · E(1/2) = π/(4ϖ) + ϖ/2      (the second-kind period ⇒ quasiperiod 1/ϖ NATIVE)
    · ϖ = K(1/2) = Γ(1/4)²/(4√π)  (CM period at τ=i, disc −4; K_theta=(π/2)θ₃² matches ellipk 0.0)
  So the length-1 layer produces exactly the four ring generators π, ϖ, 1/ϖ, G.
- **Length-2 layer lands in the ring** (dps40, guarded PSLQ): the length-2 cusp-to-cusp
  double-Eichler ∫₀¹K(∫₀^k K)dk = **2G² exactly** (h=2 closure, residual 2.3e-41) — a length-2
  iterated integral of the weight-2 Γ(2) form producing a graded-ring value (here G², weight 4).
- **Ring dimension:** {π:1, ϖ:1, 1/ϖ:−1, G:2} weight-graded to ≤3 (cap: transc. degree b+c≤2,
  G-degree ≤1, no simultaneous ϖ & 1/ϖ) → **dim 24** (weight 0/1/2/3 = 4/6/7/7). Convention-
  dependent; consistent with the memo's dim-20 (tighter cap). The weight-3 generators (no 1/ϖ):
  {ϖG, πG, πϖ², π²ϖ, π³}, plus the 1/ϖ-carrying ones {π²G/ϖ, π³G/ϖ²}. **G enters only at weight ≥2**,
  so a closed form REQUIRING G is necessarily weight-3 — matching the sharpened target.

## Deliverable B — the pullback / the DECISIVE discriminator
Driver routeC_gamma2_pullback.py.

The heuristic that makes T2 "look like" a length-2 wt-3 Γ(2) MMV: the (s,t)→τ change of
variables supplies two weight-2 Eisenstein Jacobians λ'(τ) (one per Feynman integration), and
the fibre supplies the period ⇒ value in {π,ϖ,1/ϖ,G}. **This heuristic has a gap: it treats the
fibre as a period. It is not.**

**Test (1) — the fibre is a resurgent/irregular period, NOT a modular period (PROVEN):**
The one-mass fibre N(D) = (1/√c₁)∫₁^∞ e^{−Dx}/√Q dx (T2's fibre as the second mass →0) has:
  · N(0) = K(1−ρ)  = a Γ(2) modular period (Fuchsian);
  · N(D>0) = a period of the rank-4 IRREGULAR connection eq:pf (corpus: Poincaré rank 1 at ∞,
    four exponential sectors λ∈{±1,±i√((1−ρ)/ρ)}; test_L4_irregular_at_infinity).
Independent reconfirmation this session (rho=1/5, D=8,12; dps60): the large-D Watson series
N ~ e^{−D}Σ b_k D^{−(k+1/2)} is **asymptotic (divergent)** — the truncation error DECREASES to a
minimum at k* then RISES (k*=14 at D=8, k*=25 at D=12; k*∝D), the unmistakable Gevrey-1 signature
(a convergent q-series would keep improving). The Borel radius is **S=2 exactly**:
|b_{k+1}/b_k|/(k+½)→0.5027→1/S, matching the corpus's dominant Borel singularity at ζ=−2 (branch
points {−2,−1±iω}). **A finite iterated-Eisenstein integral (MMV) is a FUCHSIAN period — a
convergent q-series with NO essential singularity at the cusp.** The fibre is not that. The D=1
physical fibre carries the SAME exponential (Bessel) twist.

**The structural conclusion:** T2 = ∫∫ (irregular-period fibre) over the modulus X(2).
Integrating irregular-period fibres over the modular curve does not restore Fuchsian-ness — the
exponential twist survives the ρ-integration (it lives in the SCALE direction, orthogonal to the
modular direction). Hence T2 is an **integrated irregular period = a Γ(2) resurgent Lambert
series** (Broadhurst–Dorigoni arXiv:2607.14020, at level Γ(2) — the universal Legendre level,
BELOW BD's Γ₁(6) sunrise). The pure Γ(2) MMV in {π,ϖ,1/ϖ,G} is its regular D→0 shadow.

**Test (2) — guarded PSLQ reconfirms the numerics track is underpowered:** frozen T2 (~19 dig)
vs the corrected weight-3 ring (dim 24) at tol 1e-16/1e-18, decoy-controlled: the only relations
found are high-height (h=118) and matched by the decoy (h=175), or V-coeff-0 basis-internal
identities. No resolvable low-height weight-3 closure at ~19 dig; decisive test needs ~32 dig.
(Duplicates/agrees with the parallel numerics track; included as an independent check.)

## Why the length/weight/level is "same numbers, different class"
- length = 2 (two Feynman integrations) — unchanged.
- level = Γ(2) (Legendre universal family, λ=1−ρ) — unchanged.
- weight = 3 (needs G) — unchanged FOR THE SHADOW.
- CLASS: finite Fuchsian MMV (what "length-2 wt-3 Γ(2) iterated-Eisenstein integral" names) →
  **irregular/exponential resurgent Lambert series** (Bessel moment of the Legendre family,
  Fresán–Sabbah–Yu irregular period). T2 lives in the larger class; the finite MMV is its D→0 tip.

## The sharp hand-off (what closes it)
Two routes, both specialist:
1. **BD Eichler summation for the twisted integrand.** Pull the D=1 fibre's Bessel twist through
   the λ=1−ρ modular map; Eichler-integrate the weight-2 Γ(2) Eisenstein series against the twist
   to get the explicit resurgent Lambert series Σ a(n) qⁿ/(1−qⁿ) (+ Stokes/trans-series data), as
   BD do for sunrise/banana at Γ₁(6). The Γ(2) analog is length-2. This IS the closed form; it is
   resurgent, not a finite ring element (unless it collapses — see 2). Machinery: Broedel–Duhr
   1803.10256 / 1912.00077 (iterated-Eisenstein numeric engine), BD 2607.14020 (the twist).
2. **32-digit PSLQ vs {π,ϖ,1/ϖ,G} weight-3.** Only decides the residual "does the resurgent
   series collapse to a finite ring element at the collinear/CM point" question. Blocked by the
   OUTER-integral corner convergence (~1e-19 wall; see sprint_T2_modular_derivation_memo.md);
   needs a cleaner outer quadrature or the modular/q-series rep (= route 1).
The natural closing step is route 1 (Brown/Kleinschmidt for the MMV/Eisenstein machinery; Avery
for the momentum-space Sturmian side). Route 2 only fences the collapse question.

## Transcendental tags (CLAUDE.md rule)
- π: Paper 18 Layer-2 M1 pure-Tate (π^{2k}·ℚ on S³); Paper 34 Hopf-measure / temporal-compactif.
- ϖ=K(1/2)=Γ(1/4)²/(4√π): elliptic genus-1 layer (Paper 18 §Level-2 genus grading, v4.82.0);
  disc-4 CM period; cosmic-Galois elliptic Rung 1 (memory cosmic_galois_elliptic_rung1).
- G=Catalan=β(2)=L(2,χ₋₄): weight-2 Eisenstein L-value, native to the Γ(2)/Legendre family
  (∫₀¹K=2G). Paper 18 "Catalan G via vertex parity" home; here the Eisenstein-L-value instance.

## Files
- debug/routeC_gamma2_eisenstein_basis.py — (A) weight-2 space, length-1 dictionary, length-2
  landing, ring dimension. debug/data/gamma2_eisenstein_basis.json.
- debug/routeC_gamma2_pullback.py — (B) fibre-resurgence discriminator + guarded T2 PSLQ +
  structural conclusion.
- NOTHING committed; NO paper edits (PI integrates both tracks). Frozen T2 headline untouched.
