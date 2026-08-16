# Cosmic-Galois probe — RUNG 1 = decisive GO (2026-08-16)

**Question (PI, this session).** Does the Paper 59 elliptic Bessel moment carry a
*cosmic Galois* structure — i.e., does GeoVac's aspirational **elliptic** cosmic
Galois (Paper 56's Hain–Brown mixed-elliptic-motive target, one storey above the
established mixed-Tate injection into `G_4`) get *populated* by a concrete chemistry
period? Rung 1 (the modular-level test) is the decisive GO/STOP. **Verdict: GO, and
cleaner than expected.** Driver `debug/routeC_cosmic_galois_rung1.py`.

## The object
The Route C elliptic family (normalized c1=1, c2=ρ): `E_ρ : y² = (x²−1)(ρx²+1−ρ)`,
base `ρ = c₂/c₁ = t(1−t)/[s(1−s)]`. Its two periods are the verified `K(1−ρ)`, `K(ρ)`
(solutions of the modulus PF `M_ρ`), so the modular parameter is
`τ(ρ) = i K(ρ)/K(1−ρ)`.

## Results (all exact numerics, `debug/routeC_cosmic_galois_rung1.py`)

**A. Modular identification — CONFIRMED to 1e-41.** The theta-series modular lambda of
τ(ρ) equals the modulus exactly: `λ(τ(ρ)) = 1 − ρ` (ρ=0.37→0.63, 0.6→0.4, 0.15→0.85,
0.85→0.15, all |diff| ≤ 1e-41). So the family is **the Legendre / Γ(2) universal
family** (the universal elliptic curve over the modular curve X(2)), pulled back along
the **rational-linear** modulus map `λ = 1 − ρ`. The periods `K(1−ρ), K(ρ)` ARE the
Γ(2) modular periods.

**B. Non-isotrivial.** `j(ρ)` varies (69194, 968452, 6221 at ρ=0.37, 0.6, 0.15) — a
genuine family, not a constant curve.

**C. Three cusps.** Singular fibers over ρ∈{0,1,∞} (↔ λ∈{1,0,∞}) — matching Γ(2)'s
three cusps. ρ=1 is the c₁=c₂ diagonal (genus-0 nodal degeneration; τ→i∞).

**D. Physical domain.** `ρ ∈ [0.36, 2.78]` over (s,t)∈[0,1]², **includes ρ=1** (the
cusp) and **ρ=1/2** (a CM fiber, see E).

**E. Rung-2 CM data point — EXACT.** ρ=1/2 is in the physical range and gives **τ=i**
(CM by ℤ[i], lemniscatic). The period there is a Γ-value (Chowla–Selberg):
`K(1/2) = Γ(1/4)²/(4√π)` to 30 digits (|diff|=0). **A molecular integral's period, at a
physically-visited fiber, is a Γ-value** — the Route C period genuinely populates the
modular/Γ-value (cosmic-Galois) ring.

## What this means (honest scope)
- **GO for the elliptic cosmic Galois at the FAMILY/fiber-period level.** GeoVac's
  Paper 56 aspirational elliptic target (Hain–Brown MEM over `M_{1,1}`) is populated by
  a concrete chemistry period: the family is modular (Γ(2)), its fiber periods are
  modular periods, and at the CM fiber ρ=1/2 the value is a Γ-value. The earlier
  abstract-panel "Test A" negative (Sym² panel periods, depth 1–2) is NOT contradicted —
  that tested a different object; this is a concrete Feynman-period object and it lands.
- **Level = Γ(2)** — the classical *universal* Legendre level, LOWER than the sunrise's
  Γ₁(6). Honest caveat: landing on the universal Legendre family is expected for any
  4-branch-point family (any elliptic curve appears); the *informative* content is that
  the modulus map is **rational and linear** (λ=1−ρ) and the physical domain contains a
  **CM fiber** (ρ=1/2) and a **cusp** (ρ=1). It is a real, concrete instantiation, not
  an exotic-level discovery.
- **Rung 3 remains open (= the frontier).** Whether the *integrated* T2 (the (s,t)
  integral over the family) is a **Γ(2) multiple modular value / elliptic
  polylogarithm** — the cosmic-Galois statement about the *observable*, not just the
  fiber periods — is the same open object as the ABW closed form (memo
  `sprint_routeC_momentum_memo.md` §"(a) ABW push"). Rung 1 identifies its modular home
  (Γ(2)); realizing it as an MMV is the collaboration/frontier piece.

## Rung 2 — periods are Γ-values across discriminants (CONFIRMED)
`debug/routeC_cosmic_galois_rung2.py`. Since ρ sweeps all of (0,∞) = all of X(2), the
domain hits every CM fiber. Confirmed the period is a Γ-value at TWO distinct
fundamental discriminants (so it tracks the discriminant, not a lemniscatic fluke),
each by a **direct closed-form residual** (~1e-51, independent of PSLQ):
- **disc −4** (τ=i, ρ=1/2): `K = Γ(1/4)²/(4√π)` — PSLQ log-relation `[2,1,4,−4]`,
  two-precision-stable; residual 2.7e-51.
- **disc −8** (τ=i√2, ρ=2√2−2≈0.828): `K = (1+√2)^{1/2}Γ(1/8)Γ(3/8)/(2^{13/4}√π)` —
  residual 2.7e-51. **PSLQ discipline note:** the disc−8 PSLQ leg (`[4,2,13,−2,−4,−4]`,
  the correct relation) was NOT two-precision-stable — a genuine PSLQ trap flagged by the
  cross-precision check; the confirmation rests on the direct closed-form residual, not
  the PSLQ. (Followed the corpus PSLQ discipline: cross-precision caught it.)

Verdict: the periods populate the modular/Γ-value (cosmic-Galois) ring **systematically**.

## Rung 3 — the integrated T2 as a Γ(2) multiple modular value (FRONTIER)
`debug/routeC_cosmic_galois_rung3.py`. The observable-level statement: the integrated
T2 pulls back (via λ=1−ρ) to an iterated integral over X(2) of period-weighted forms =
a candidate **Γ(2) MMV / elliptic polylogarithm** — the SAME open object as the ABW
closed form. Four independent lines converge on "genuine Γ(2) elliptic MMV, not
reducible": (a) modular home Γ(2) [Rung 1]; (b) fiber periods are Γ-values [Rung 2];
(c) NOT in the weight-0/1/2 polylog ring [`routeC_weight_probe.py`]; (d) the single
fiber does not close it [ABW in-module obstruction]. Bounded probe: the integrated
collinear value is not a low-height single-fiber period product (precision-limited,
illustrative). **Two frontier obstacles, both flagged honestly:** (i) explicit
realization needs the Γ(2) iterated-Eisenstein/elliptic-polylog basis; (ii) even a
PSLQ-grade *value* of the integrated observable needs a dedicated evaluator — naive
nested tanh-sinh times out (>260 s at dps≥18). = the collaboration/frontier piece
(Avery track, HELD).

## Files
`debug/routeC_cosmic_galois_rung1.py` (A–E), `debug/routeC_cosmic_galois_rung2.py`
(CM-fiber Γ-value confirmation), `debug/routeC_cosmic_galois_rung3.py` (frontier setup
+ bounded probe). Cross-refs: Paper 56
`sec:period_map` (the mixed-Tate period map this extends to genus 1), Paper 59
(the elliptic object), memory `hain_brown_identification` (the target + source list
incl. Tapušković 2303.17534, the sunrise cosmic Galois).
