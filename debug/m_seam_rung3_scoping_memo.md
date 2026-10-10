# m-seam rung 3 — scoping (2026-10-10, PI-directed)

**Question.** Can a multi-centre Coulomb-Sturmian basis (Avery's world) capture the
bending relay that a one-centre expansion could not converge (rung 2S: a2 at R = 2 moves
0.7..4.5 with n_max and l_max), before the l > 0 conditioning wall (Paper 60, 2026-10-09)
forces canonical orthogonalisation to stall?

Context: rungs 1/2/2S in `tests/test_renner_teller_{m_splitting,rung2,rung2_sturmian}.py`;
chronicle in CHANGELOG [Unreleased]. Model throughout: one electron, central Z_c = 1, two
protons at distance R in directions (s, 0, ±c), bend s = sin(half-angle).

## 1. Inventory

| Piece | Exists? | Where |
|---|---|---|
| One-centre Sturmian, any l, bent protons, fixed-S and fixed-R | yes | rung 2S tests |
| Two-centre SW metric, sigma basis (m = 0) along the axis, prolate quadrature | yes | `geovac/sturmian_sigma_basis.py` |
| Three-centre (water) SW metric, s-only (distance-only blocks), momentum chi-reduction | yes (debug) | `debug/sturmian_sw_water_conditioning.py`, `debug/sturmian_sw_momentum.py` |
| Two-centre Sturmian s-orbital integral engine (overlap, V, T, ERI) | yes, s only | `geovac/sturmian_integrals.py` |
| **Two-centre SW block for general (l, m) in arbitrary orientation** | **no** | to build |
| Canonical orthogonalisation (Avery's remedy) | yes | `CANON_TOL` in `sturmian_sigma_basis.py` |

## 2. The structural fact that makes rung 3 feasible

In the **fixed-S (Shibuya-Wulfman kernel) posing** the one-electron many-centre problem is
the standard Hermitian eigenproblem K C = k C, K = D^1/2 S D^1/2, D = diag(Z_a/n), with
S_{t't} = (n_t/k)<chi_t'|1/r_{a(t)}|chi_t> = <chi_t'|(p^2+k^2)/(2k^2)|chi_t>. Every entry
involves only the two centres the functions sit on. **Three nuclei need no three-centre
integral.** (Shibuya-Wulfman 1965; Avery & Avery; module header of
`sturmian_sigma_basis.py`.) The fixed-R Galerkin posing would need <chi_a|1/r_c|chi_b>,
c != a, b: a three-centre one-electron integral. Avoid it.

Caveats, all recorded on 2026-10-09:
- The fixed-S posing converges more slowly than Galerkin (H2+ at R = 2: 8.5e-2, 2.7e-2,
  8.2e-3 Ha at 3/6/15 functions per centre vs Galerkin 8.5e-4, 9.1e-5, 5.4e-6).
- It is non-variational.
- Each nucleus's potential is carried by its own functions, so every nucleus needs a
  basis. Rung 3 is multi-centre by construction; it is not a one-centre basis plus a few
  proton functions.
- Fixed R by self-consistency k = lambda(kR), exactly as in rung 2S (shape fixed, S
  scales all distances).

## 3. The missing piece and how to build it

Two-centre SW block for general (n l m) at separation vector d:

    S_{nlm, n'l'm'}(d) = sum_L 4 pi i^L Gaunt(l m, l' m', L M) Y_LM(d^) R^L_{nl,n'l'}(k d),
    R^L = int f_nl(p) f_n'l'(p) (p^2+k^2)/(2k^2) j_L(p d) p^2 dp,

where f_nl are the Sturmian momentum radials (rational in p). R^L is a closed form,
polynomial(kd) e^{-kd} (residues), so it is exact in sympy; mpmath quadrature is the
referee. Then convert to real harmonics.

Validation gates (all to 1e-10):
1. s-s blocks reproduce the water chi-reduction;
2. l > 0, m = 0 along z reproduce `sturmian_sigma_basis.py` (independent prolate route);
3. d = 0 gives the intra-centre identity;
4. rotation covariance: the block at d^ equals the Wigner-rotated z-axis block;
5. S_ab = S_ba^T.

## 4. Plan with decision gates

- **Step 1: general SW block** (the bulk of the work). Gate = §3 validations.
- **Step 2: rung 3 pi energy.** Same model, R = 2; central n <= 3 (l <= 2) and protons
  n <= 2 (1s, 2s, 2p) as the smallest meaningful basis, then grow. Continuity check: the
  true model energy is the common limit of rung 2S (l -> infinity) and rung 3, so their
  trends must meet. The rung 2S references are only that: the variational referee gives
  -0.727050 (l <= 2), and the secular value is -0.73681 (l <= 3; l = 4 leaves it
  unchanged). Gate: the pi energy changes by less than 1e-3 between consecutive bases.
- **Step 3: a2 vs basis and canonical threshold** (1e-8..1e-5). Reuse rung 2S's
  self-consistency + Richardson machinery.
  - **GO:** consecutive-basis change < 10% and threshold-insensitive.
  - **BORDERLINE:** converges only at a fixed threshold.
  - **STOP:** the threshold dependence is at least as large as the signal.
  - Prediction: a2 is a small second-order quantity, the one most exposed to the
    directions that canonical orthogonalisation drops, so the wall is likely to show here
    first.

## 5. Guardrails and dead ends checked

- Papers 8-9: not applicable (multi-centre basis).
- §3 rows nearby, NOT re-attempted: "point-group gerade lever for polyatomic conditioning
  (water C2v)" and "low-rank/Woodbury flattening of the SW metric penalty". Rung 3
  measures conditioning; it does not try to cure it.
- Composition wall (Löwdin densification) concerns qubit encodings; not relevant to this
  classical one-electron calculation.

## 6. What either outcome means

The question is whether the Shibuya-Wulfman fixed-S posing, in a basis containing pi
functions, resolves the bending relay (the |m| = 1 Renner-Teller coefficient) before
canonical orthogonalisation stalls. A positive answer is the first m != 0 multi-centre
Sturmian observable. A negative answer locates the l > 0 wall on a physical quantity
rather than on cond(S). Reference material on the SW posing, if needed: Avery's work in
`C:\Users\jlout\Desktop\Avery` (outside the repo).

## 7. Outcome (same day)

- **Step 1 GO:** `geovac/sturmian_general_sw.py` + `tests/test_sturmian_general_sw.py`.
  The radial factor is an exact closed form (one residue at p = i; no pi).
- **Step 2:** the pi energy converges with the PROTON basis: -0.490, -0.657, -0.7182,
  -0.7192. Growing the centre basis barely moves it.
- **Step 3: STOP for a converged a2 in this posing.** a2 swings (+1.12, -0.40, -0.43,
  -1.39, -1.55, +1.11) because one sigma state sits 0.022-0.027 from the pi state in k,
  with bend coupling ~0.18, and it crosses the pi level between the two largest bases.
- **The predicted conditioning failure did not occur.** The fixed-S kernel posing never
  inverts S. The obstacle is one slowly converging sigma level.
- **Next lever:** pin the sigma-pi gap at linear geometry with a more accurate method,
  then treat sigma-pi as a 2x2 (pseudo-Jahn-Teller / Walsh) problem in s, which stays
  valid at physical bends, where a2 does not. *(Premise corrected by §8: the physical
  gap is not small.)*

## 8. Rung 4 (same day): the converged referee

`tests/test_renner_teller_rung4.py`. A complete radial grid (FE-DVR) plus partial waves
to l = 55, with exact second-order perturbation theory at linear geometry:
**a2 = 0.44941, E_pi = -0.7385753 at R = 2.**
- The near-degeneracy in §7 was a property of the fixed-S posing: there, levels are
  compared at one shared S, not one shared R. At fixed R the dominant sigma partner, the
  first excited z-even sigma, sits 0.093 Ha BELOW pi. It carries 0.742 of the 0.848
  sigma relay, with coupling 0.263.
- The |m| = 2 relay cancels exactly from the splitting.
- A two-state model gives 0.344: the right sign and 77% of the size.
- Rung 3 (fixed-S Sturmian) is therefore not the route to bending coefficients. A posing
  that compares levels at one shared R is, and the one-centre grid is already one.
