# Sprint memo — Neumann-r12 build (explicit correlation in the accurate engine)

**Goal (PI-approved 2026-09-19):** "explicit correlation done right, at physical
accuracy" = port explicit r12 into the ACCURATE algebraic prolate engine
(neumann_vee / neumann_vee_general_m / prolate_recondition), not the crude 5D
quadrature engine. Target: push H2 below Paper 12's 0.32 mHa toward the exact
value (physical accuracy, [[feedback-physical-accuracy-target]]).

Status: **DESIGN DE-RISKED, build not started.** Reduction identity proven at
machine precision (below). No new special functions needed — smaller than the
scoping memo's "5D IBP or Sack recurrences" estimate.

---

## 0. Corrections that reset the starting point (2026-09-18/19)

- The analytical kinetic-energy fix for r12^p is ALREADY built and is the only
  path solve_hylleraas can reach (geovac/hylleraas.py; tests/test_hylleraas.py
  asserts it). It was NOT unbuilt — the scoping memo's "named-but-unbuilt" line
  predated the v5.13.10 retraction.
- March 2026's "94.7% at 9 functions" was a FINITE-DIFFERENCE ARTIFACT (v5.13.10
  retraction). Done correctly (analytical kinetic), the crude r12 engine reaches
  ~78% at N=9, climbing to ~90% at N=18.

## 1. Head-to-head diagnostic (the go/no-go) — GO

Same crude 5D-quadrature engine + grid (N_xi=20,N_eta=14,N_phi=16, xi_max=12,
R=1.4011), analytical kinetic, D_e% vs function count. Driver:
`debug/data` via scratchpad hyl_headtohead.py (numbers below reproduced there).

| N | plain CI (p=0) | explicit r12 |
|--:|--:|--:|
| 12 | 48.3% | 75.0% (p=1) |
| 18 | — | 85.9% (p=2) |
| 24 | — | 80.9% (p=1) |
| 27 | 67.7% | — |
| 36 | — | 87.3% (p=2) |
| 72 | 69.4% | — |
| 110 | 69.7% (plateaued) | — |

Plain CI plateaus ~70% in this engine (N=27->110 barely moves) — a grid wall,
and the p=0 baseline here is handicapped, so the raw pp gap OVERSTATES the
advantage. But the direction is unambiguous and matches the physics: r12 reaches
far higher accuracy at every matched N and is still climbing where p=0 stalled.
This is the real signal the FD artifact was faking. The crude engine cannot hit
physical accuracy either way, so the absolute proof requires the Neumann engine —
i.e. the build. GO.

## 2. Design: odd r12 powers reduce to the EXISTING general-m Neumann engine

Sigma / Dphi-averaged, prolate spheroidal:
    r12^2 = (R/2)^2 (A - B cosDphi),  A = (z1-z2)^2+rho1^2+rho2^2, B = 2 rho1 rho2
    => r12^(2k+1) = (R/2)^(2k+1) (A - B cosDphi)^(k+1) * (1/r12) * (R/2)... [units]

Binomial-expand (A - B cosDphi)^(k+1); cos^j Dphi -> sum of cos(m Dphi). So every
ODD power's Dphi-average is a FINITE combination of <cos(m Dphi)/r12>_phi, i.e.
the m-th Neumann term. EVEN powers are polynomial in (A, B^2) -> existing cheap
xi/eta moments.

**Parity closes cleanly.** B = 2 rho1 rho2 carries sqrt((xi^2-1)(1-eta^2)); B^j
multiplies the m=j Neumann term whose P_l^j carry the MATCHING
sqrt(xi^2-1) sqrt(1-eta^2) factors, so B^j x (m=j term) is polynomial. No
half-integer powers survive. => the whole extension is binomial algebra on top of
neumann_vee_general_m.py (v5.12.7). NO new special functions; the "Sack /
generalized Neumann" concern is dissolved.

## 3. PoC validation (2026-09-19) — machine precision

Driver: `debug/neumann_r12_reduction_poc.py` (5 geometries, R=2, dense Dphi
quadrature ground truth). RHS averages by direct quadrature, so R1/R2 test the
algebraic reduction itself (constant-free).

| check | tests | worst rel err |
|:--|:--|--:|
| R1 | r12^1 = A<1/r12> - B<cosDphi/r12> | 1.8e-16 |
| R2 | r12^3 -> m=0,1,2 Neumann terms | 1.8e-16 |
| R3 | corpus m=0 Neumann vs direct <1/r12> | 5.3e-7 (l-trunc, worst when xi_< ~ xi_>) |
| E2/E4 | even powers polynomial (<r12^2>=A, <r12^4>=A^2+B^2/2) | 2.6e-14 |

## 4. Remaining work (the actual build)

**Integration surface located (2026-09-19) — neumann_vee_general_m.py::vee_matrix.**
That engine ALREADY carries a general azimuthal-channel structure: `mus` = basis
azimuthal quantum numbers, `m_set` = {a+b, |a-b|} over mu-pairs = the Neumann
m-channels, and `ms_pairs = (m, s)` with s the (xi^2-1)^s power and the eta
selection cap `l <= q_max + 2s - m`. This is EXACTLY the substrate the r12
reduction needs: the reduction's cos(m Dphi) terms route to channel m, and B^j's
sqrt((xi^2-1)(1-eta^2)) pairs with the m=j block's (xi^2-1)^{m/2} to land on an
integer s. So the r12 wiring is a modification of `vee_matrix`'s term loop, not a
new engine.

1. **Wire the reduction into `vee_matrix`.** For a basis carrying r12^p, the V_ee
   element is <phi_i | r12^{P-1} | phi_j>, P = p_i + p_j:
   - **P-1 even** -> r12^{P-1} = (A - B cosDphi)^{(P-1)/2} is polynomial: an
     OVERLAP-type integral (no 1/r12). Needs a poly-moment path (even powers of B
     are polynomial; <cos^j Dphi> are Beta constants). NEW small path.
     **VALIDATED 2026-09-19** at the r12^2 case (p=1 overlap): the algebraic
     poly-moment overlap (products of xi-moments A_n(2a) and eta-moments 2/(q+1))
     matches the crude 5D-quadrature compute_overlap_matrix on a pure-p=1 basis
     to max rel 5.4e-6 (grid-limited ground truth; algebraic side exact). Driver
     `debug/neumann_r12_evenpower_poc.py`.
   - **P-1 odd** -> r12^{P-1} = (A - B cosDphi)^{P/2} * (1/r12): expand the
     binomial into A^a B^b cos^b Dphi, map cos^b Dphi -> sum cos(m Dphi) (m<=b),
     and for each (a,b,m) add the general-m Neumann term with (p1,q1,p2,q2)
     shifted by the A^a B^b polynomial powers and routed to channel m. Reuses
     `build_Xtab` / `_A_moment` / the Ytab eta moments unchanged.
     **m=0 and m=1 kernels PINNED 2026-09-19** (the two channels r12^1 needs):
       <1/r12>_phi   = +(2/R) sum_l (2l+1) P_l(xi_<)Q_l(xi_>)P_l(eta1)P_l(eta2)   [validated R3]
       <cosDphi/r12>_phi = -(2/R) sum_l (2l+1)[(l-1)!/(l+1)!]^2
                            P_l^1(xi_<)Q_l^1(xi_>)P_l^1(eta1)P_l^1(eta2)
     with EXPLICIT derivative-based associated functions -- P_l^1(x)=sqrt|x^2-1|*P_l'(x),
     Q_l^1(xi)=sqrt(xi^2-1)*Q_l'(xi) -- because scipy lpmv/lqmn return NaN/garbage
     for xi>1 at higher l. Constant/sign fitted empirically: ratio = -1.00000 at
     all 4 test points (driver `debug/.../m1_constant.py`; to be moved to debug/).
     PARITY: B=2 rho1 rho2 carries sqrt((xi^2-1)(1-eta^2)) on each electron; the
     m=1 element's P_l^1 carry the MATCHING sqrt, so B*(m=1 term) is polynomial ->
     integer (xi^2-1)^s, i.e. the neumann_vee_general_m (m,s) structure exactly.
     **MATRIX VALIDATED 2026-09-19**: the full odd-power V_ee assembly (A K0 - B K1,
     ordered-xi integrals + m=0 C_l / m=1 D_l eta-moments, the -1 constant, parity
     cancellation) reproduces the 5D-quadrature V_ee for a sigma p=1 basis to max
     rel 9.2e-6. Driver `debug/neumann_r12_oddpower_vee_poc.py`. So BOTH new integral
     families -- even-power overlap moment and odd-power m=0/m=1 V_ee -- are validated.
     Ordered-xi here is 2D Gauss; production swaps in compute_Xl (m=0) /
     neumann_vee_general_m.build_Xtab (m=1,s=1), both corpus-validated.

## 4b. What remains for a first ENERGY (updated 2026-09-19)

- **Kinetic -- VALIDATED, and it is NOT hard.** For phi = g r12 the Green-form
  <T> = (1/2)∫(grad g_i.grad g_j) r12^2 + (vector cross) + g_i g_j |grad r12|^2.
  |grad_k r12|^2 = 1 exactly, and the vector cross term
  (1/2)∫(r1-r2).(grad1-grad2)(g_i g_j) collapses under IBP (div1(r1-r2)=3,
  div2=-3) to -3<g_i|g_j>; the |grad r12|^2 term adds +<g_i|g_j>, giving
     <phi_i|T|phi_j> = (1/2)∫(grad g_i.grad g_j) r12^2 dV - 2 <g_i|g_j>_{p=0}
  BOTH terms r12-EVEN. VALIDATED vs the crude analytical kinetic to max rel
  9.4e-5 (`debug/neumann_r12_kinetic_poc.py`). **The prior "last hard piece
  (James-Coolidge)" framing was WRONG -- the vector terms cancel.**
- **V_ne -- the "two-Coulomb hybrid" DISSOLVES for homonuclear H2 (VALIDATED).**
  I feared the mixed p0xp1 V_ne = <g_i g_j r12/r1A> was a Roothaan-Ruedenberg
  hybrid integral. It is NOT, for H2: the homonuclear symmetry gives
  1/r1A + 1/r1B = (4/R) xi1/(xi1^2-eta1^2), so (1/r1A+1/r1B)*J1 = (R^2/2) xi1 --
  the (xi1^2-eta1^2) CANCELS the Jacobian and V_ne*J1J2 is POLYNOMIAL:
  -(R^2/2)(R/2)^3 [xi1(xi2^2-eta2^2)+xi2(xi1^2-eta1^2)]. So p0xp1 V_ne is just the
  ODD r12^1 machinery with that V_ne polynomial (P_Vne) in place of Jp1 Jp2.
  VALIDATED vs 5D quadrature to max rel 8.3e-6 (`debug/neumann_r12_vne_p0p1_poc.py`).
  p1xp1 V_ne = r12^2 V_ne = EVEN (A*V_ne poly). (The hybrid integral only bites
  HETERONUCLEAR, where Z_A != Z_B leaves a genuine 1/(xi+-eta) -- not our case.)
- **NO HARD INTEGRAL REMAINS for the mixed p={0,1} homonuclear H2 energy.** Every
  block reduces to validated even/odd machinery:
    overlap: p0p0 std, p0p1 r12^1 odd (validated integral), p1p1 r12^2 even (validated)
    V_ee:    p0p0 r12^-1 (existing Neumann), p0p1 r12^0 = <g_i g_j> std, p1p1 odd (validated)
    V_ne:    p0p0 homonuclear-poly std, p0p1 odd*P_Vne (validated), p1p1 even
    kinetic: p0p0 std g-kinetic, p1p1 even-collapse (validated), p0p1 VALIDATED (below)
- **p0p1 kinetic: VALIDATED 2026-09-19 (max rel 2.4e-5).** Prolate form: the
  azimuthal term VANISHES (p0 bra phi-independent, d_phi phi_i=0), leaving
  <d_xi1 phi_j>_phi = d_xi1 g_j <r12>_phi + g_j (R/2)^2 [d_xi1 A K0 - d_xi1 B K1]/2
  (and same for eta1, xi2, eta2), with <r12>_phi=(R/2)^2[A K0 - B K1] and K0,K1 the
  pointwise Neumann sums. d_xi1 B carries the sqrt but routes through K1 (m=1) and
  cancels, as predicted. Validated vs the crude analytical kinetic p0xp1 subblock
  (`debug/neumann_r12_kinetic_p0p1_poc.py`). **ALL integral blocks now validated.**
- **First-energy plan:** the fully-algebraic mixed p={0,1} energy is the target and
  is now pure ASSEMBLY -- wire all blocks (each even/odd, all validated except the
  p0p1 kinetic which needs implementing) into one (H, S) at consistent prefactors
  and diagonalize vs exact -1.174476 Ha. No hard integral gates it any longer.

## 6. UNIFIED ENGINE built + VALIDATED (2026-09-19)

`scratchpad/unified_engine.py`: one (H,S) builder for mixed p={0,1}, Delta-phi
average via the validated Neumann kernels (K0,K1) + geometry A,B, xi/eta by Gauss
quadrature, kinetic block-dispatched (p0p0 std / p1p1 even-collapse / p0p1 odd).
S unifies via <r12^(pi+pj)>, V_ne via <r12^(pi+pj)>*Vne, V_ee via <r12^(pi+pj-1)>.

**Validated -- and it exposed that the CRUDE engine is the grid-limited one.**
At grid (20,14), N=6: my S matches crude S to 7e-6; my H differs by 18 mHa, my E
= -1.151301 vs crude -1.132645. FIRST mis-step: I mis-attributed this to a bug in
my engine after an N_phi test showed crude stable -- but N_phi was the WRONG axis.
Arbitration on xi,eta: the algebraic Neumann V_ee (analytic, no grid; two
independent Neumann impls agree at 1.22671) and the crude ELLIPTIC V_ee disagree
by 8%; refining the elliptic's xi,eta grid CONVERGES it to the Neumann value
(9.2% -> 3.1% -> 1.4% -> 0.6% at nx=90). So the algebraic Neumann V_ee is EXACT
and the crude was xi,eta-grid-limited by the 1/r12 cusp. Confirmed at the energy
level: the crude MIXED-basis E converges straight to my unified value as xi,eta
refines -- +18.7 -> +6.0 -> +2.0 mHa at (20,14)/(36,26)/(56,40). **My unified
engine is the accurate one.**
LESSON (recorded): the crude engine is NOT a trustworthy V_ee reference at coarse
xi,eta; validate against a grid-CONVERGED crude or the analytic Neumann, not the
default grid. Also: my pointwise-K0 V_ee actually SMOOTHS the cusp (truncated
Legendre sum) so it integrates well on a coarse grid -- closer to exact than the
singular elliptic at the same grid.
CONSEQUENCE: the earlier "first r12 energy" 82.68% (crude, N=24, Sec.5) was itself
xi,eta-under-converged and UNDERSTATED r12's benefit; the accurate engine gives
more (N=6 already 86.7% D_e).

**Remaining for physical accuracy (unchanged in kind, now the only work):**
scale to a large mixed basis; make xi fully algebraic (X_l / build_Xtab) rather
than Gauss for speed+exactness; handle conditioning (Paper-12-style re-basing +
precision). No new physics -- the engine is validated.

## 7. ALGEBRAIC CORE + CONDITIONING WALL measured (2026-09-19)

`debug/neumann_r12_algebraic_core.py`: grid-free S+Vee+Vne via mono-moments +
ordered-xi X_l + analytic eta moments (the scalable form; kinetic not yet in it).
- **S validated: max rel 1.6e-7 vs grid-converged crude (nx=90).** Core correct.
- **cond(S) vs basis (float64): 1.1e5 (N=6) -> 7.6e9 (N=24) -> 9.0e11 (N=92).**
  Extrapolates past float64's ~1e16 by N~150-200. So a large-basis physical-accuracy
  energy CANNOT be float64 -- **mpf + re-basing is MANDATORY** (confirms Paper 12's
  route, now measured for the r12 case).
- **Build speed:** my quick core is O(N^2)-pairs x expensive-inner (Cl/Dl recompute
  Legendre per call): 2.4s (N=6) -> 46s (N=24) -> 941s (N=92). NOT fundamental --
  the production engine must precompute Cl/Dl tables (cf. compute_Cl_table) and use
  the algebraic X_l recurrence. A perf task, separate from the mpf/conditioning one.

## 8. MPF ENGINE (physical-accuracy build; PI-committed 2026-09-19)

`debug/prolate_r12_mpf.py` — extends Paper 12's mpf engine (`prolate_recondition.py`
+ `neumann_vee_general_m.py`) with r12. Mapping: EVEN r12 powers factorize
(reuse `_ov`/`_vne`/`_kin` with A/P_Vne power shifts); ODD powers reuse the
`_build_Xtab_mp` Neumann X-table (X_l^{0,0}=build_X(0), X_l^{1,1}=build_X(1)).
Convention = prolate UNSYM ProductFn (symmetric ground state via the eigensolve).
Each block validated vs a FIXED unsymmetrized quadrature (the first quad reference
had a missing factor-2 in the exponent — the overlap of two e^{-a(xi1+xi2)} is
e^{-2a...}; fixed 2026-09-20).

ALL FIVE integral blocks validated (mpf, vs quad, all ~1e-4 grid-limited):
- overlap even (r12^2): 1.2e-4
- V_ee odd (r12^1): 6.2e-5
- V_ne mixed (p0p0/p1p1 even + p0p1 odd, homonuclear P_Vne): 1.4e-5
- kinetic p1p1 (even-collapse (R/2)^2 KG_A - 2<g|g>): 1.4e-4
- kinetic p0p0: reuse `one_body_mp` `_kin` (already Paper-12-validated)
- **kinetic p0p1 (odd): 9.05e-05** -- the last block (2026-09-20)
- overlap odd (p0p1) + V_ee even (p0p1=<gg>) + V_ee p0p0 (r12^-1, K0=`vee_mp`):
  covered by the same machinery (odd = the V_ee r12^1 core; rest are p=0 paths)

**p0p1 kinetic collapsed CLEANLY (the anticipated hard block was easy).** The
IBP identity  T = 1/2<grad phi_u . grad phi_v>  with phi_u=g_u (p0), phi_v=g_v r12
(p1) splits into Piece A = 1/2<(grad g_u.grad g_v) r12> and Piece B =
1/2<g_v (grad g_u . grad r12)>. IBP on Piece B moves grad off r12:
Piece B = -1/2<r12 (grad g_u.grad g_v)> - 1/2<r12 g_v grad^2 g_u>. **Piece A
cancels the first term of Piece B exactly**, leaving
   **T_{p0xp1} = -1/2 int r12 g_v (grad^2 g_u) dV**
-- a PURE odd-r12^1 element with NO grad_r12 term and NO sqrt-cancellation fold.
The Laplacian NUMERATOR L[g_u] = d_xi((xi^2-1)d_xi g_u)+d_eta((1-eta^2)d_eta g_u)
(the metric denominator cancels the electron's Jacobian) is a polynomial, so it
routes through the *same* A K0 - B K1 machinery as the validated odd V_ee, with
base = g_v * L[g_u] * (Jacobian of the other electron). Prefactor -1/2 h6 (2pi)^2
(2/R) (one (R/2)^2 less than V_ee, from grad^2's 1/(R/2)^2). Impl:
`_L_polys`, `_kin_odd_base`, `_odd_K0`/`_odd_K1`, `kinetic_mixed_mpf`.

ENGINE COMPLETE. Assembly (`assemble_mixed`) + first energy driver
(`debug/r12ci_first_energy.py`) built and run.

**FIRST ENERGY from the exact mpf engine (2026-09-20, sigma-only, R=1.4011):**
  p=0 only (no r12), n=9 :  E_tot=-1.118024  D_e=67.6%  err=-56.5 mHa
  p={0,1} (with r12), n=18: E_tot=-1.153819  D_e=88.2%  err=-20.7 mHa
r12 lifts D_e 67.6%->88.2% at matched radial basis, cutting error 2.7x. The
residual at 88% is ANGULAR correlation (l=m=0 here), not conditioning
(cond(S)=4.8e8, all kept). Next lever = angular channels (l_max).

**VERIFIED RESULT (2026-09-20) -- BEATS Paper 12's headline at a modest basis.**
Convergence + mpf cross-check ladder (`debug/r12ci_convergence_ladder.py`,
`assemble_mixed` + `solve_canonical_mpf`; gerade sigma product basis, mu=0):

  (j=2,l=0) n= 18  cond=4.8e8   E_tot=-1.1538193  D_e=88.16%  err=-20.66 mHa
  (j=2,l=2) n= 90  cond=1.2e11  E_tot=-1.1743084  D_e=99.905% err=-0.1666 mHa

  * mpf solve == float64 solve to **0.00 uHa** at BOTH points -> the -0.167 mHa
    is NOT a conditioning artifact (float64 downcast is clean at cond<=1.2e11).
  * both variational (E_tot ABOVE exact -1.174475); monotone decrease.
  * **-0.167 mHa vs Paper 12's re-based-CI 0.32 mHa** -- ~2x better, and a
    DIFFERENT method (explicit r12 vs CI expansion). p=0 control at n=90 = 92.25%
    (-13.5 mHa); r12 cuts the error 80x.
  * full-grid engine validation (N_xi=24,eta=18,phi=28): all 5 blocks pass,
    p0xp1 kinetic 1.955e-5, full mixed matrix 9.35e-5.

Held OUT of Paper 12 / CLAUDE.md S5 / CHANGELOG pending (a) the (3,2)/(2,4)
convergence points confirming continued descent toward exact, and (b) a PI call
(new result beating the headline = PI decides capture, not a PM re-measure).

Convergence CONTINUES (3rd point, old code):
  (j=3,l=2) n=160 cond=1.2e15  float err=-0.0789 mHa (kept 158/160)
                               mpf   err=-0.0787 mHa (kept 160/160, dfloat=0.22uHa)
  -> monotone -20.66 -> -0.167 -> -0.079 mHa, all ABOVE exact; at n=160
     conditioning bites (float64 drops 2 vecs) and the mpf solve pulls ahead.
     BOTH verification checks PASS; result is solid.

**ODD ROUTING OPTIMIZED (2026-09-20) -- 65x.** The odd-r12^1 value of a unit
monomial (1,P1,Q1,P2,Q2) is a fixed scalar G_odd[P1,Q1,P2,Q2]=(A K0)+(B K1),
memoized once (`_make_godd`/`_odd_value`/`_odd_context`); the X-table is built
ONCE per assembly and shared by `_kern_odd1`, `vne_mpf`, `kinetic_mixed_mpf`
(was 3 separate builds). n=90 assembly **653s -> 10s**; energies reproduce the
quad-validated old code to <0.02 uHa at (2,0) and (2,2) -> refactor correct.

**SCALING LADDER (optimized engine + mpf-orthogonalized solve, 2026-09-20):**
  (3,2) n=160 cond=1.2e15  mpf err=-0.0787 mHa
  (2,4) n=234 cond=7.4e13  mpf err=-0.1403 mHa
  (3,3) n=256 cond=2.6e16  mpf err=-0.0620 mHa
  (4,3) n=400 cond=2.7e19  mpf err=-0.0528 mHa (float -0.0688; mpf +16 uHa better)
  (3,4) n=416 cond=5.0e18  mpf err=-0.0535 mHa
=> PLATEAU at ~-0.05 mHa (99.97% D_e), still 6x past Paper 12's 0.32 mHa.
   At n=400 the mpf solve is decisively needed (float64 drops 97/400 vectors).

**alpha is NOT the lever (scan at (3,3), 2026-09-20):**
  a=1.0 -0.0625 | a=1.2 -0.0603(best) | a=1.4 -0.2066 | a=1.6 -0.895 | a=1.8 -2.86
  Optimum ~1.0-1.2 (2 uHa gain); Paper 12's a=1.40 HURTS -- explicit r12 already
  supplies the correlation a contracted alpha was compensating for.

**FLOOR DIAGNOSED = the r12-power truncation (p<=1).** cond(S) reaching 1e19 with
mpf keeping ~all vectors => the p in {0,1}, single-alpha, sigma space is nearly
SPANNED at ~0.05 mHa. The microhartree/physical-accuracy regime needs genuinely
new basis directions = **higher r12 powers p>=2** (Kolos-Wolniewicz uses r12 up
to ~5). Build path:
  - even r12^2/r12^4 blocks: ALREADY handled by `_kern_even` (any even power via A^s).
  - odd r12^3: extend G_odd to m=0,1,2 (r12^3 = (R/2)^4 (A-BcosDphi)^2 / r12;
    (A-BcosDphi)^2 = A^2 - 2A B cosDphi + B^2 cos^2Dphi -> K0(A^2,B^2), K1(AB), K2(B^2)).
  - kinetic cross-blocks p0xp2, p1xp2, p2xp2 (same IBP-collapse + G_odd routing).
This is the next substantial sub-project (PI call).

Solve note: mpf `eigsy` O(n^3) is the large-n bottleneck now (assembly is fast).
For n in the many-hundreds, Paper 12's re-basing (`_factored_cob`/
`_normalized_solve`) keeps a float64 solve viable; raw mpf eigsy caps ~n=400.

## 9. B-PROBE: does exact-algebraic r12 generalize? (2026-09-20)

PI-directed (option B): test whether the exact-algebraic explicit-r12 survives a
second center or a third electron. First probe = **HeH+** (2e HETERONUCLEAR),
the cleanest: the whole r12 engine (overlap, V_ee, kinetic) is CHARGE-INDEPENDENT
and transfers unchanged; only V_ne changes.

**Machinery generalizes -- VERDICT YES.** Per electron
Z_A/r_iA + Z_B/r_iB = (2/R)[(Z_A+Z_B)xi + (Z_B-Z_A)eta]/(xi^2-eta^2), still
POLYNOMIAL after the Jacobian. New `_pvne_hetero_terms` / `vne_hetero_mpf` /
`assemble_hetero` (focus A=(xi+eta) carries Z_A). The (Z_B-Z_A) eta term is new:
it breaks gerade<->ungerade, so a heteronuclear basis needs BOTH angular parities
and V_ne couples them (`build_basis_full` = all l,m). Validated vs quad: 6.1e-5.
Reduces to `vne_mpf` at Z_A=Z_B=1.

**HeH+ energies (R=1.4632, E_ref=-2.97869 Ha, nuclear=Z_A Z_B/R=1.3669):**
  (2,1) n= 72  best a=1.6  E_tot=-2.84544  err=-133.2 mHa
  (2,2) n=162  best a=1.6  E_tot=-2.96704  err= -11.7 mHa   (11x better, +angular)
Variational throughout; CONVERGING toward the reference.

**Accuracy bottleneck moved and is named = single-exponent scale mismatch.** 11.7
mHa (n=162) vs H2's 0.05 mHa at comparable size: one alpha cannot serve both the
contracted He (Z=2) and diffuse H (Z=1). A BASIS limitation (Paper 12's
single-alpha strain, worse here), not a machinery failure -- separately fixable
by center-specific / two-block exponents. Driver `debug/heh_probe.py`.

### 9a. (a) HeH+ two-block -- DIAGNOSIS REFUTED (2026-09-20)

Two-exponent (two-block) support built: `_kin2_grad` (the only new piece -- the
kinetic gradient hits each function's OWN exponent; S/V_ne/V_ee cross-blocks are
exponent-symmetric, c=a_a+a_b, and reuse existing machinery) + `assemble_hetero_2block`
(p=0; V_ee stitched from 3 vee_mp calls). Sanity: two-block(1.5,1.5)==single(1.5)
to 0.00 uHa.

HeH+ p=0, single-alpha vs two-block:
  (2,1) n=36/72 : single -246.6 mHa, two-block -246.3 -> gain 0.365 mHa
  (2,2) n=81/162: single  -39.5 mHa, two-block  -39.1 -> gain 0.326 mHa
Basis SIZE (2,1)->(2,2): -246.6 -> -39.5 mHa (84% of the error).

Convergence (p={0,1}, r12, a=1.6, `debug/heh_converge.py`):
  (2,2) n=162 -11.65 mHa | (3,2) n=288 -11.41 (radial+1: -0.25) |
  (2,3) n=288 -0.88 mHa (ANGULAR+1: -10.8) -> HeH+ reaches sub-mHa, variational.
Lever = ANGULAR basis size, definitively (l_max +1 = -10.8 mHa; radial +1 = -0.25;
two-block exponent = +0.3). Machinery generalizes QUANTITATIVELY.

**My "single-alpha scale mismatch" diagnosis for HeH+ is REFUTED.** Two-block
barely helps (~0.3 mHa) and radial +1 barely helps (-0.25); the lever is ANGULAR
basis size (-10.8 mHa for l_max 2->3), not the exponent.
Consistent with Paper 12's H2 two-block plateau and the physics: HeH+'s 2
electrons share ONE sigma bond = ONE length scale, so a second exponent has no
work to do. Multi-exponent is a MULTI-SHELL lever (Li core+valence, v5.14.9:
81.7->33.2 mHa), NOT a heteronuclear-charge one. Driver `debug/heh_2block.py`.

### 9b. (b) 3-electron -- ALREADY ANSWERED in the corpus; boundary is N=4

Per the current-state rule, checked existing work before building: the corpus has
a 3e (and Ne) R12-CI infrastructure (`debug/r12ci_3e_triangle_kernel.py`,
`r12ci_3e_vertex_rules.py`, `r12ci_per_shell_lambda.py`, promoted
`geovac/transcorrelated_sturmian.py`).

**N=3: exact-algebraic explicit-r12 closes with NO resolution-of-identity.** For
N=3 the operator product <Phi|F H F|Phi> stays at most 3-body (all indices in
{1,2,3}). The term inventory (verified 2026-09-20, `r12ci_3e_vertex_rules.py`):
all 6 shapes close with closed-form angular rules -- RULE A (<P_a P_b>=d_a0 d_b0,
shared vertex -> L=0), RULE B (vec-vec kA factorization), and the TRIANGLE rule
<P_a(12)P_b(13)P_c(23)> = delta_abc/(2a+1)^2 (only all-equal multipoles survive;
a matrix contraction, not a 6D quadrature). No vector-leg triangle exists.
"Remaining work is ASSEMBLY, not new angular derivation." A working Li engine
exists (s-only, v5.14.9). Multi-lambda keeps integrals closed-form
(`shibuya_wulfman._hydrogenic_poly_coeffs_lam`).

**N=4 (Be): the WALL -- first case that genuinely needs 4-body operators**
(disjoint-pair term f_ij (1/r_kl) with 4 distinct indices; no <=3-body
reduction). Documented in the triangle-kernel derivation + the TC/Be arc (CLAUDE
S3: TC 3-body collapse, TC second-quantization plateau ~3.4%).

**COMBINED B VERDICT:** exact-algebraic explicit-r12 generalizes across a second
center (HeH+, this session) AND a third electron (Li, existing) -- exact, no RI,
up to N=3. The boundary/wall is N=4 (4-body operators). The heteronuclear accuracy
gap is basis-size; the multi-exponent lever is multi-shell, not heteronuclear.

## 5. FIRST r12 ENERGY (grid-limited crude engine, 2026-09-19)

Matched (j,l)=(2,1) truncation, converged-ish grid (22x16x18), 3-point alpha scan,
crude 5D-analytical engine (`debug/data/.../first_r12_energy`):
  p=0 (plain CI)  N=12  E=-1.092369  D_e=52.93%  err=82.1 mHa
  p={0,1} (+r12)  N=24  E=-1.144270  D_e=82.68%  err=30.2 mHa
So explicit r12 nearly DOUBLES the captured binding (+29.75 pp) at matched
truncation -- the arc's premise (r12 helps H2) confirmed in an actual energy.
**But this is NOT physical accuracy and NOT competitive with Paper 12 (99.81%,
0.32 mHa):** it is grid-limited (crude 5D quadrature) AND small-basis (N=24),
and the p=1 build took 490 s -- which is exactly WHY the algebraic (grid-free,
fast) engine matters. The algebraic blocks are validated to reproduce this engine
to ~1e-5; the remaining work is assembling them at scale (grid-free, large basis),
where r12's per-function efficiency (head-to-head) is the path to the 0.32 mHa
residual. The 30.2 mHa here is a PoC that r12 works, not the load-bearing result.
   A and B in engine variables: A = (z1-z2)^2 + rho1^2 + rho2^2 (z=xi*eta,
   rho^2=(xi^2-1)(1-eta^2)) -> polynomial in (xi1,eta1,xi2,eta2); B = 2 rho1 rho2.
   The (xi1 eta1 - xi2 eta2)^2 cross term mixes electrons' xi,eta -> expands into
   a small sum of separable (p1,q1)x(p2,q2) monomials. Book-keeping, not new math.
2. **Kinetic for r12^p.** d(r12^p)/dx = (p/2) r12^(p-2) d(r12^2)/dx, and
   d(r12^2)/dx is polynomial, so IBP kinetic terms reduce to r12^(even/odd)
   moments by the SAME identity. To validate (not yet done).
3. **Selection rules / truncation.** Confirm the eta selection rule and Neumann
   l-truncation carry over per m-channel (the general-m code already imposes
   `l > Q + 2s - m` + parity; v5.12.9 gotcha).
4. **Validate energy** vs 5D-quadrature ground truth on small bases, then vs the
   exact H2 (Kolos-Wolniewicz -1.174476 Ha) — target the ~0.19 mHa cusp residual
   that p=0 basis growth cannot reach.

## 5. Scope verdict

Smaller than the scoping memo (`debug/sprint_explicit_correlation_scoping_memo.md`)
estimated: no new special functions, the odd-power engine already exists
(general-m, v5.12.7), and the reduction is exact algebra validated to 1e-16. The
real work is wiring + the kinetic reduction + validation.
