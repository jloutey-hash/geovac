# Sprint B3 Phase 1 — boost/modular-flow seminorm probe (2026-06-10, v3.113.0)

**Goal:** first probe of the named Lorentzian candidate from the sharpened Paper 45 Q1:
can the wedge boost — the modular flow of the four-witness theorem (Paper 42,
K = J_polar in the doubled `two_m_j` convention) — serve as the TEMPORAL generator of
a quantum-metric structure, so that signature enters through the generator?

**Verdict: POSITIVE-STRUCTURED.** Four exact results on the j ≤ 1 Peter–Weyl window
with the full reachable multiplier system (rank 55, the n_max = 3 system dimension),
boost K = 2J_z (`debug/wh7_b3_boost_seminorm_probe.py`, JSON in `debug/data/`,
falsifier `tests/test_wh7_b3_boost.py`, 6/6; Paper 45 Q1 paragraph applied, 24 pp
GATE: PASS):

| Check | Result |
|:--|:--|
| T0 system rank | 55 = Σ_{b≤2}(2b+1)² exactly |
| T1 flow compatibility | σ_{2π}(F) = F at 1.3×10⁻¹⁵ for ALL 55 elements (four-witness closure on this substrate); **half-period grading σ_π(F) = (−1)^{2b}F bit-exact (6.6×10⁻¹⁶)** — a finite-cutoff spin–statistics shadow: the modular flow grades the operator system by band parity (Bose/Fermi) at β/2 |
| T2 boost-alone kernel | dim ker(ad_K) = **9 = Σ_{b∈{0,1,2}}(2b+1)** exactly — the boost-invariant (axial) multipliers. The structured middle between the two known failure modes: NOT the P45 annihilation (dim 55 = everything), NOT the metric kernel condition (dim 1) |
| T3 frame completion | joint kernel of {J_x, J_y, J_z} = ℂ1 (dim 1) — transversal legs restore the metric condition |
| T4 causal classifier | Q(F) = ‖[J_z,F]‖² − ‖[J_x,F]‖² − ‖[J_y,F]‖²: sign agrees with the symbol classifier sign(2m′² − b(b+1)) in **9/9 classes**. m′ = 0 purely spacelike (ratio −1 exact); top weights timelike from b = 3/2; **b = 1 top weight ON the cone (q_min = 0 to 10⁻⁶, class straddles)** — a genuine null ray in the operator algebra at the exact integer the symbol algebra predicts |

## Reading

The geometry demanded T2 and delivered it: a single boost flow cannot metrize the
transversal directions (boost orbits don't reach them), and the kernel is exactly the
boost-invariant subalgebra with the right dimension count — by contrast with June 9's
annihilation, where the kernel was *everything* because the algebra was wrong. The
boost is a sound temporal leg. T4 is the first genuinely Lorentzian-flavored structure
in the corpus: a signed leg comparison that classifies operators as
timelike/null/spacelike, with the cone passing through an exact integer locus
(2m′² = b(b+1) at b = 1, |m′| = 1). T1's half-period grading ties the candidate to the
four-witness/KMS machinery: the same flow that closes at 2π grades by spin at π.

**Caveats (honest):** the off-extremal Q-ratios deviate from the idealized symbol
values (compression effects); several measured ratios look like simple rationals
(e.g. +3/5 at b = 2 top, −1/3 at b = 1/2) — flagged as unverified observations per the
numerical-coincidence rule, closed-form identification is follow-up work. The
classifier convention (sum-of-squares legs) is one of several; sign pattern stability
under convention changes should be checked in Phase 2. No convergence statement is
made or implied yet.

## Follow-ons (B3 Phase 2)
1. Closed forms for the Q-ratios (symbol calculus on compressions; PSLQ only as
   post-identification per the discrete-for-skeleton rule).
2. State-level reverse-triangle structure along boost orbits (the twin-paradox
   super-additivity of Paper 49's TICI algebra, now with a non-degenerate metric leg).
3. Wedge restriction (HemisphericWedge of `geovac/modular_hamiltonian.py`) — repeat
   T2/T4 on the wedge-compressed system; edge effects at the equator.
4. The convergence question: a Mondino–Sämann-shaped statement (Lorentzian pre-length
   convergence) with the boost leg as temporal generator — the actual B3 prize.
