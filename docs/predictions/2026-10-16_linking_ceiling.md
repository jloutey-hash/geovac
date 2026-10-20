# Pre-registered prediction: the chemistry period ladder climbs with linked densities, not with nuclei

**Registered:** 2026-10-16, before any computation in this file (PI-directed: "let's test that").
**Extends:** `2026-10-13_four_centre_genus.md` (P1) and `2026-10-13_three_electron_k3.md` (P3a-P3e);
Paper 59 `sec:four_centre`, `sec:three_electron`.
**Not committed at registration** (no commit authorized this session); computations below start after this file is written.

## The question

Does the transcendence class of a molecular integral climb without limit, or is it capped by
how the electrons are linked?

## Correction recorded first

An electron repulsion integral (ERI) involves at most four orbitals, so at most four centres, and it
couples exactly two charge densities. So for Hamiltonian matrix elements the rung is capped at the
elliptic one (Paper 59 + P1) whatever the molecule. A climb can only come from integrals that link
n >= 3 densities through n - 1 or more kernels, as explicitly correlated (r12, geminal) methods use.

## Setup

A tree integral: n two-centre 1s transition densities rho_v, joined by n - 1 rational kernels
(Coulomb or geminal). At the fingerprint fibre (D = 0, phases trivial) each density carries one
dispersion root Delta_v = sqrt(c_v L_v^2 + 1), L_v the signed sum of the momenta on its kernels.
The holomorphic period of a tree is P_n = int_{R^{n-1}} prod_v 1/Delta_v(L_v(k)) d^{n-1}k
(the convention of P3b/P3d: one-dimensional momenta).

## Predictions

**PL1 (tree identity, symbolic).** For every tree with n densities,
P_n = 2^n / (pi sqrt(prod c_v)) * int_0^inf prod_v K0(a_v x) dx,  a_v = c_v^{-1/2},
independent of the tree's shape (star, chain, any). Reason: Fourier-transforming each 1/Delta_v
to K0(|y|/sqrt c_v)/(pi sqrt c_v) turns each of the n - 1 momentum integrals into a delta function
that identifies the two endpoint positions, and a tree is connected. n = 3 reproduces P3d's
8/(pi sqrt(ca cb cc)) V3. So the rung depends on how many densities are linked, not on the pattern.

**PL2 (n = 4 opens a new rung).** The one-mass family F_4(t) = int_0^inf K0(x)^3 K0(t x) dx satisfies
a linear ODE in t whose part carrying F_4 is NOT built from elliptic (order-two) pieces: its solutions
satisfy no constant-coefficient quadratic relation (a symmetric cube of an order-two operator would
give 3; a tensor product of two would give 1).
Controls, which must come out as stated or the test is void:
- F_2(t) = int K0(x) K0(t x) dx: order two, elliptic, 0 relations.
- F_3(t) = int K0(x)^2 K0(t x) dx: elliptic squared (Bailey/BBBG eq. 33), so at least 1 relation.

**PL3 (consequence, not separately tested).** If PL1 and PL2 hold: Hamiltonian integrals sit at the
elliptic rung for any molecule; three-density correlated integrals sit at elliptic squared (no new
class at the holomorphic level, by Bailey); four-density integrals are the first to climb. The
ceiling is set by the method (how many densities one integral links), not by the system size.

**PL4 (hunch, low confidence).** The climb continues: F_5 is also not built from lower rungs
(consistent with Broadhurst's five-Bessel modular L-values). Tested only if cheap.

## Falsifiers

1. PL1: a numerical evaluation of a four-density tree (star or chain) disagrees with the
   Bessel-moment formula, or star and chain differ.
2. PL2: the irreducible part of F_4's operator carries one or more quadratic relations among its
   solutions, i.e. it factors into elliptic pieces as F_3 does.
3. Void: a control fails (F_3 shows no relation, or F_2 is not order two).

## Honest caps, stated before the result

- Everything here is a fibre statement at D = 0, like every rung in Paper 59. Whether the Feynman
  integrations lower the class of the integrated value is the separate open PI item 1
  (the exchange class shows they can).
- The holomorphic period is a lower bound on the class: P3e showed the true fibre carries more.
- One-dimensional momenta are the P3b/P3d convention. The real three-dimensional integrals reduce
  to Bessel moments of the same family with other x-weights (P3e's control: x^5 K1K1K1 reduces to
  the x^0 periods at n = 3); that this carries to n = 4 is assumed, not shown.

## Result (2026-10-16, same session; scripts `debug/linking_ceiling_{A_tree,B_ode,C_quadrics,D_yukawa,D2_local,D4_irred_n}.py`, data in `debug/data/linking_ceiling_*`)

**PL1 HOLDS. PL2 HOLDS. Controls pass. PL4: partly (see below).**

| Check | Outcome |
|:--|:--|
| PL1, four-density tree, star and chain (fast-decay variant, 3-D sinh-trapezoid) | both equal the Bessel-moment formula to 4e-9 to 2e-7 relative, at the grid's own convergence level, at three scale sets including equal scales. Shape does not matter. The identity itself is a one-line Fourier argument (delta functions on a connected tree), so this checks the bookkeeping, not a conjecture. |
| Step B: exact ODE in t for F_n, n = 2..5 (Mellin route) | orders 2, 3, 4, 5. Each checked on F_n by direct quadrature: relative residual 1e-28 to 1e-32 at t = 0.3, 0.7, 1.6 (30 digits). |
| Control F_2 | theta^2 F = t^2 (theta+1)^2 F: the Legendre (elliptic K) equation. 0 quadratic relations. |
| Control F_3 | 4 theta^3 F = t^2 (theta+1)^3 F, i.e. 3F2(1/2,1/2,1/2;1,1;t^2/4) (argument corrected 2026-10-17 from t^2: the recurrence of the operator gives t^2/4, singular at t = 2), which is Clausen's square of 2F1(1/4,1/4;1;.): symmetric square of an elliptic equation. 1 quadratic relation. Matches Bailey. |
| **F_4** | order 4, **0 quadratic relations** (exact, at t0 = 1/2 and 7/5, 60 Taylor terms). In x = 1/(4t^2) it is the MUM operator th^4 - x(2th+1)^2(10th^2+10th+3) + 9x^2(2th+1)^2(2th+3)^2 (symbolic identity); in z = x/4 its holomorphic solution has integer coefficients B_n = C(2n,n)^2 sum_k C(n,k)^2 C(2k,k) (checked n < 40; integral to n = 300): the Hadamard product of C(2n,n)^2 with the equal-mass sunrise sequence. Integral mirror map. Yukawa coupling K(q) = 1 + 16q + 1152q^2 + 99520q^3 + ... (integer coefficients), **not constant**, so not a symmetric cube. Conifold points x = 1/36, 1/4 with exponents (0,1,1,2). Irreducible: the holomorphic solution satisfies no ODE of order 1, 2 or 3 with polynomial coefficients of degree <= 40 (rank mod a 61-bit prime; the same search finds the operator itself at order 4, degree 2). |
| F_5 (PL4) | order 5, MUM at t = inf, irreducible (same search, degree <= 40), **1 quadratic relation**: an orthogonal (rank-5, K3-type) structure. A symmetric fourth power of an elliptic equation would give 6; so F_5 is not built from elliptic pieces either. Whether it is new relative to F_4 is not settled: SO(5) is the exterior-square image of Sp(4), so F_5 could be tied to the F_4 rung. |

**Reading.** An irreducible order-4 local system with no invariant quadric and a non-constant
Yukawa coupling has differential Galois group Sp(4) or SL(4). Neither is a quotient of a product
of SL(2)'s, so F_4 is not a sum of products of elliptic periods with algebraic moduli. The same
then holds for the generic four-mass family, since restricting a Bailey-type elliptic product
formula to masses (1,1,1,t) would give one. So the ladder, at the holomorphic-period level, is:

| densities linked by one integral | where it lands |
|:--|:--|
| 2 (every Hamiltonian ERI) | elliptic |
| 3 (three-electron correlated integrals) | K3, collapses to elliptic squared (Bailey) |
| 4 (four-electron, e.g. the bridge class) | Calabi-Yau-type, Sp(4): the first genuinely new rung |
| 5 | orthogonal rank 5, not elliptic-built; relation to the n = 4 rung open |

**The ceiling is the method, not the molecule.** Ordinary Hamiltonian integrals stay elliptic for
any molecule. The class climbs only with the number of densities one integral links.

**Caps.** (1) Fibre and holomorphic-period statements only; the integrated-value question (PI item
1) is untouched and applies to every rung. (2) The four-Bessel object is the one-dimensional
moment family; that the real three-dimensional four-electron integral reduces into the same
family is assumed (P3e checked the analogue at n = 3 only). (3) The instanton numbers from K(q) are 16, 142, 11056/3, 121470, ...; with Broadhurst's normalisation N_k = N0 n_k, N0 = 3, they are integral. Prior-art scan (same day, after the result): the operator is AESZ #58 (CYDB 2.9, a Hadamard product), with the same holomorphic coefficients and Yukawa coupling; see Paper 59 `sec:four_density` and CHANGELOG v5.20.20 for the full prior art. (4) Irreducibility is shown up to the stated
degree bound, not for unbounded degree.

**Follow-up (2026-10-17; memos `debug/sprint_{banana_dictionary,minimality_audit}_memo.md`).**
The operators of F_n for n = 2 to 5 are the Hadamard products of the Legendre operator with the
equal-mass banana operators with n - 1 propagators, in y = 1/(16 t^2), checked exactly; n = 4 is
AESZ #58 and n = 5 is the fifth-order operator the tables print with entry #39 (prior art:
Almkvist-Zudilin, math/0402386; AESZ, math/0507430). So five densities give no new type, and the
open question of this record about n = 5 is answered by the literature. Irreducibility for n = 2
to 5 now follows from local exponents and Fuchs' relation, with no degree bound. Paper 59
`sec:four_density` carries this.

## Novelty classification (appended 2026-10-19)

Classes per `docs/predictions/README.md`, set after the fact by a literature audit (`debug/prediction_novelty_audit_memo.md`; this record predates the novelty field). FINDING = held, not found in the literature; CHECK = held, already known; MISS = failed.

- PL1: FINDING, routine (the standard position-space form of banana integrals)
- PL2: CHECK (AESZ #58)
- PL3: UNSCORED
- PL4: CHECK
