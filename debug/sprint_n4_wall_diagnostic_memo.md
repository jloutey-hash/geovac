# Sprint memo — the N=4 explicit-r12 "wall" diagnostic (is it RI-forced, or soft?)

**Date:** 2026-09-20/21. **Type:** diagnostic (no new energy). **Verdict:** the documented
"N=4 = THE WALL, needs 4-body operators, no <=3-body reduction" is **SOFT for the scalar
Coulomb bridging term** — that term is exact, RI-free, terminating, and reducible;
confirmed with a number on the Be-relevant integral (reduced closed form == 12-D brute
Monte-Carlo, rel 2.5e-4). The genuine hard walls are separate axes (quantum encoding;
non-Hermitian TC; kinetic-vector 4-body untested).

## 0. Why this sprint

Standing state before: exact-algebraic explicit-r12 (James-Coolidge style — r12 in the
BASIS, **no** strong-orthogonality projector, hence **no** resolution-of-identity) closes
with no RI through N=3 via three angular rules (RULE A shared-vertex->L=0; RULE B vec-vec
factorization; TRIANGLE `<P_a(12)P_b(13)P_c(23)>=d_abc/(2a+1)^2`). The corpus recorded N=4
(Be) as "THE WALL — first case that genuinely needs 4-body operators, no <=3-body
reduction" (`memory/r12_generalization_boundary_n3_n4.md`,
`debug/sprint_neumann_r12_build_memo.md` sec.9b). PI's open hope: the N<=3 mapping might
"point toward a general resolution."

The strong reading of that wall — that N=4 is where exact-no-RI *stops*, an RI-class
obstruction like Gaussian F12 — was never tested. This sprint tests it.

**Distinct from the documented dead end** (ledger 2026-08-23, "TC three-body operator
collapse via Gaunt/6j"): that tried to COLLAPSE a 3-body operator to 2-body via an abelian
plane-wave momentum trick and FAILED (GeoVac's Y.Y=sum_Lambda is non-abelian → "on the
Gaussian 3-body side"). This sprint does the OPPOSITE — it accepts no collapse and asks
whether the genuine 4-body integral is still finite/closed-form. Hermitian variational
scalar object, not the non-Hermitian TC commutator with vector legs.

## 1. The object

In `<Phi| F H F |Phi>` (F = sum_{p<q} f_pq multiplicative), the ONLY genuinely-4-body-
connected term at N=4 is the scalar Coulomb CHAIN

    T4 = f(r12) * f(r34) * (1/r13)        graph 2-1-3-4, four distinct electrons

**Term-enumeration argument (why this is the only one):** kinetic gradients are pair-local
(`grad_k f12 . grad_k f34 = 0` for disjoint pairs since {1,2}∩{3,4}=∅) and V_ne is one-body,
so NEITHER can bridge two disjoint correlation edges. Only the two-body Coulomb `1/r_kl` can
connect electrons one-from-each-pair. The linear-in-F terms `<Phi|f12(1/r34)|Phi>` with
disjoint indices FACTORIZE (2-body x 2-body); the genuine 4-body needs both f's AND a
bridging Coulomb, which first appears at second order (F H F). The N=3 case never produced
a chain because with 3 electrons every pair shares a vertex (→ triangle, closes via d_abc).

## 2. Q1 — TERMINATION (`debug/r12ci_4e_wall_diagnostic.py`)

Expand `1/r13 = sum_L g_L(r1,r3) P_L(u1.u3)`. Angular factor
`J(a,c,L; l1..l4) = < w1 w2 w3 w4 P_a(u1.u2) P_c(u3.u4) P_L(u1.u3) >`, w_i=|Y_li|^2. MC scan
of the largest nonzero L:

| bridge orbital | largest nonzero L | 2*l_bridge |
|:--|:--:|:--:|
| s | 0 | 0 |
| p | 2 | 2 |
| d | 4 | 4 |

**The Coulomb multipole sum TRUNCATES at L = 2*l_bridge** — bounded by the finite orbital
angular content, NOT the (unbounded) Coulomb. Gaussian F12 needs RI precisely to avoid
these 3-/4-electron integrals; GeoVac's terminating angular sum computes them exactly.
**→ the "no RI" property carries through N=4.**

(Odd correlation multipoles vanish against even diagonal densities |Y_l|^2 — a toy
artifact, not a selection rule; the termination bound is physical and general: everything
meeting a bridge vertex has finite multipole content.)

## 3. Q2 — REDUCIBILITY (`debug/r12ci_4e_wall_q2_reducibility.py`)

By the Legendre addition theorem `P_L(u1.u3)=4pi/(2L+1) sum_M Y_LM(u1) conj(Y_LM(u3))`, the
chain FACTORIZES across the bridge:
`J_chain(a,c,L) = 4pi/(2L+1) sum_M A_LM conj(C_LM)`, each vertex a 3-body-ish object (center
+ leaf + one free bridge index). Confirmed numerically on a fully-active case (leaves l=1,
bridge l=1, a=c=2) at every channel Q1 says carries signal: L=0 rel 6.1e-3, L=2 rel 2.8e-3.
(L=4 is identically 0 by Q1; its 6.8% is noise/noise.)

**→ the connected 4-body angular object is a finite (L,M) contraction of THREE-body vertex
kernels** — the chain analogue of the N=3 `triangle_contract`, NOT an irreducible 4-index
blow-up. Radial cost is a separable O(n^2) chain matmul.

## 4. THE Be NUMBER — reduced closed form vs brute force (`debug/r12ci_4e_be_integral.py`)

The genuinely-4-body Be-relevant integral (bridge e1,e3 = 2p, Z=1.0, Be's angular
correlating space; leaves e2,e4 = 1s, Z=3.7, the core; f=exp(-0.5 r)):

    I = INT rho1 rho2 rho3 rho4  f(r12) f(r34) (1/r13)  d3r1..d3r4

computed two independent ways:
- **BRUTE:** full 12-D importance-sampled Monte-Carlo (no reduction used).
- **REDUCED:** 1s LEAVES integrate out into a 1D radial DRESSING of their bridge partner
  (`Phi_leaf(r) = INT |phi_1s|^2 f(r12) d3r2` = spherical avg of f), leaving a standard
  TWO-electron Slater-Condon Coulomb integral of the dressed 2p densities,
  `I = sum_{L=0,2} Theta^L R^L`. Deterministic; L-sum terminates at L=2; no RI.

| | value |
|:--|:--|
| REDUCED (deterministic) | **5.08756e-2** (L=0: 4.7208e-2 ; L=2: 3.6675e-3) |
| BRUTE N=2e7 | 5.0874e-2 ± 1.9e-5 → **0.1 sigma** (rel 2.7e-5) |
| BRUTE N=8e7 | 5.0863e-2 ± 9.5e-6 → **1.4 sigma** (rel 2.5e-4) |

Agreement at the MC-limited ~2e-4 level. **The 4-electron integral reduces exactly to a 1D
leaf dressing + a standard 2-electron Coulomb integral, no RI, no truncation.** The bridge
L=2 channel carries 7.2% of the total → the angular bridge is genuinely exercised (not a
trivial monopole product). This is the physically dominant Be 4-body class (2p correlating
space ↔ 1s core via a bridging Coulomb), not a toy.

**Process note (independent-route cross-check, per `feedback_independent_route_crosscheck`):**
the first run disagreed 22% (1540 sigma). Brute was tight/converged, so the REDUCED side had
the bug — a dropped `(2L+1)` in a hand-derived `Theta^L`. Replaced the hand-derived constant
with a deterministic angular quadrature → agreement. The cross-check caught a false
"confirmation" before it landed; the number is trustworthy because two unrelated methods meet.

## 5. Verdict, scoped honestly

**The N=4 scalar Coulomb bridging term is NOT an RI/decidability wall.** It is exact,
RI-free, its Coulomb sum terminates at L=2*l_bridge, and it reduces to a finite bridge-(L,M)
contraction of vertex kernels (for s-leaves: literally a 2-electron Slater-Condon integral of
dressed densities). The corpus's "N=4 needs 4-body operators / no <=3-body reduction" is TRUE
as a term statement but its *implication* (exact-no-RI stops here) is **refuted**. The general
resolution the PI hoped for is concrete: **the bridge addition theorem** — every disjoint-pair
coupling factorizes across its bridging two-body operator into vertex kernels; likely extends
to N>=5 (more bridges, same per-bridge factorization), testable.

**What this is NOT (do not overclaim):**
1. A diagnostic on angular structure + ONE validated integral, NOT a Be R12-CI energy. It
   removes the "RI wall" belief and shows the path is a polynomial build; it delivers no
   spectroscopic number.
2. Scoped to the **scalar Coulomb** chain. The variational `<Phi|F H F|Phi>` also has
   KINETIC 4-body pieces — argued above not to bridge (pair-local gradients), but NOT
   numerically tested. The non-Hermitian TC route (vector legs, the 2026-08-23 non-abelian
   failure) is a different, harder operator.
3. The **QUANTUM-ENCODING wall is separate and real**: an explicit 4-body correlation
   operator is a 4-body Pauli string (ledger: TC 2nd-quant plateau 3.4%, angular-gradient
   2.66x Pauli) — the wall for the quantum-simulation product, independent of the classical
   closed form.

**Ahead vs catch-up (the PI's framing question):** on this narrow axis GeoVac is genuinely
AHEAD of Gaussian-F12 (which uses RI for exactly these 4-electron integrals) — the angular
closed-forms keep them exact and RI-free. It is a classical-integral advantage; the
product-relevant wall (quantum encoding) is unmoved.

## 6. Next step (owed, PI call)

Convert the structural verdict to a full result: build the Be (or a 4-electron) R12-CI matrix
element end-to-end (including exchange chains — same structural class, permuted labels) via
the reduced form, vs a brute reference, and quote an energy. Then the "exact, RI-free, N=4"
claim carries a spectroscopic number, not just an integral.

## 7. Files

- `debug/r12ci_4e_wall_diagnostic.py` — term inventory + Q1 termination.
- `debug/r12ci_4e_wall_q2_reducibility.py` — Q2 bridge factorization.
- `debug/r12ci_4e_be_integral.py` — the Be number (reduced == brute).

## 7b. Be R12-CI ENERGY (the "4-body inside an energy" build, 2026-09-21)

PI-directed follow-on: turn the validated 4-body reduction into a Be R12-CI correlation
energy. 2x2 variational ansatz {Phi_0, F Phi_0}, Phi_0 = Be 1s^2 2s^2 (minimal Slater,
E0 = -14.539, 34 mHa from HF limit), F = sum f_ij explicit correlation.

**Three findings, each load-bearing:**
1. **The linear-geminal energy is ill-conditioned.** With f=1-e^{-r} (->1 at large r),
   Fbar~5 and {Phi_0, FPhi_0} are 99.8% parallel; the coupling h = H_01 - Fbar*E0 is a
   ~0.03 residual of two ~71-magnitude numbers, so h needs H_01 to ~1e-4 relative --
   beyond any MC. Every pure-MC route (VMC nodes; matrix-element MC heavy Coulomb tails +
   generalized-eigensolver noise bias) failed for this reason.
2. **A short-range geminal fixes conditioning.** f = r e^{-2r} (same cusp f'(0)=1, ->0
   large r): parallelism 0.994->0.86, Fbar 4.9->0.3. Now h is a moderate difference.
3. **The ill-conditioned pieces must be analytic; the rest can be MC.** Built the noise-free
   analytic block-reduction engine (`debug/be_r12ci_full.py`): all-s orbitals -> every term
   is a radial integral over up-block {1,2} x down-block {3,4} with monopole kernels.
   Validated: S_00/S_01/S_11 exact (Slater-Condon to 3e-13); H_01/S_00 = -4.4642 matches MC
   -4.458 (0.15%), grid-converged 0.013%. So **h = hT+hV = 0.135-0.167 = -0.03232 is exact**
   (both halves analytic), and **sigma2 = 0.015167 is exact**. g = <G|H|G> is
   well-conditioned (no cancellation) -> its 0.5% MC value suffices.

**RESULT: E_R12 = -14.5572 Ha, correlation captured -18.2 mHa (19% of Be's 94 mHa true
correlation), variational (above exact -14.6674).** A single short-range explicit-r12
geminal, with the 4-body RI-free content inside gV = E[(F-Fbar)^2 V_ee], computed with the
ill-conditioned coupling h and sigma2 done ANALYTICALLY (exact) and the well-conditioned g
by MC. This is the "4-body RI-free evaluation works inside a Be energy" demonstration with a
trustworthy number.

**NOT overclaimed:** PoC-level (single geminal, minimal basis; not spectroscopic -- capturing
19% of correlation). The full analytic H_11 (F^2 V_ee three-pair) would make g exact too, but
g is well-conditioned so it changes the answer < MC noise; the engine's purpose (kill the
ill-conditioning by making h, sigma2 exact) is achieved. Files: `debug/be_r12ci_{reference,
matelem,ortho,analytic,engine,full,4body_exchange}.py`.

## 8. Records touched

CHANGELOG v5.15.8; CLAUDE.md §2 one-liner; `memory/r12_generalization_boundary_n3_n4.md`
updated (boundary reframed N=4-hard-wall → N=4-scalar-soft-wall). **Flagged to PI (not
self-applied):** the walls register (`docs/walls/register.md`) N=4 entry and any Paper 12
explicit-correlation-section note — reframing a registered wall is a /walls / PI matter.
