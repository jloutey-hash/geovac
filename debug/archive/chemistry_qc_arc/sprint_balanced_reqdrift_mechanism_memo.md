# Sprint: balanced-coupled LiH R_eq-drift MECHANISM PROOF (diagnostic)

**Date:** 2026-09-20
**Mode:** diagnostic-only (documented wall area). No production fix; this pins the
mechanism of an established wall, per the mission's WALL deliverable ("a proven
obstruction with the mechanism pinned", CLAUDE.md §1.7).
**Builds on:** v5.14.8, which localized the drift to ONE term (cross-center V_ne
R-slope, 8.8% too weak) via an exact Hellmann--Feynman force decomposition, and
*named* the mechanism ("the fixed, unpolarized hydrogenic bond orbital Z_orb=1
carries no R-adaptive relaxation") without demonstrating it. Owning paper: 19,
Sec. "Fixed-geometry energy versus well shape". Memo predecessor:
`debug/sprint_balanced_coupled_reqdrift_memo.md`.

**Verdict:** the named mechanism is **DEMONSTRATED** (analytically, artifact-free)
and the wall is **STRENGTHENED** — direct variational relaxation is structurally
inexpressible in the zero-parameter balanced construction.

---

## 1. The mechanism, demonstrated analytically (artifact-free)

The balanced LiH well is a strict two-term balance: V_NN (outward,
dV_NN/dR = -Z_Li Z_H/R^2 = -0.3300 Ha/bohr at R_true=3.015) vs the cross-center
V_ne (inward). At the true minimum the cross-V_ne R-slope delivers only +0.3018,
leaving the residual outward tilt -0.0290 (= 8.8% of the required force). The
question v5.14.8 left open: *is that deficit really the bond orbital's finite
extent (lack of contraction)?*

**Closed-form model.** For a 1s Slater orbital of exponent z centered at A,
attracted to a nucleus of charge Z at B (separation R):

    <phi_z | -Z/r_B | phi_z> = -(Z/R) [ 1 - (1 + zR) e^{-2zR} ]

Its R-slope (the Hellmann--Feynman cross force) as a function of contraction z,
at R = 3.015, Z = Z_Li = 3 (the H-side bond electron attracted to Li):

| z    | dV/dR (Ha/bohr) | % of point-charge limit |
|-----:|----------------:|------------------------:|
| 1.00 |        +0.31001 |                   93.9% |
| 1.30 |        +0.32488 |                   98.4% |
| 1.50 |        +0.32804 |                   99.4% |
| 2.00 |        +0.32986 |                  100.0% |
| 3.00 |        +0.33002 |                  100.0% |
| inf  |        +0.33002 |                  100.0% |

**The point-charge limit is +Z/R^2 = +0.3300 = the required force Z_Li Z_H/R^2
exactly** (Z_H = 1). This is not a coincidence: a fully-localized (point) density
attracts the partner nucleus with exactly the classical Coulomb force that
cancels V_NN — the Hellmann--Feynman force-balance identity. The finite-extent
factor [1 - (1+zR)e^{-2zR}] < 1 both weakens the attraction and *softens its
R-slope*; contraction (larger z) removes the softening.

**So the mechanism is pinned:** the R_eq drift is the finite-extent screening of
the bond orbital's cross-attraction; full contraction recovers the exact
point-charge force = exact force balance = correct bond length. The fixed z=1
orbital delivers 93.9% of it — a **6.1% deficit** from a single 1s orbital, same
sign and direction as the measured **8.8%**. The residual 2.7 pp is the more
diffuse 2s/2p admixture (larger relative screening) and the second, weaker cross
term (Li-side electron -> H, Z=1); an exact occupation-weighted two-sided
reproduction is available but not required — mechanism and limit are exact.

Driver: `debug/balanced_reqdrift_relaxation_probe.py` (analytic block inline in
the sprint transcript / CHANGELOG). Fast, closed-form.

## 2. Why direct relaxation is inexpressible here (honest negative + structural finding)

The obvious test — let the bond exponent relax R-adaptively and watch the drift
collapse — was run first and **failed by construction**, which is itself the
finding.

Setting the bond block's shared exponent z (both Z_center and Z_partner) and
re-solving the balanced FCI on an R grid gave a **variational catastrophe**:

    E_relax(R_true) = -9.565 Ha   (exact LiH = -8.071; below the true floor)
    z*(R) pegged at the grid ceiling (2.30) at every R; gain "1636 mHa"

Energies below exact are impossible for a valid variational method, so this is
not physics. **Verified cause:** the balanced construction's one-body term is the
hydrogenic baseline -Z_orb^2/(2 n^2) (documented at `balanced_coupled.py:503` as
the W1c-residual issue), i.e. **Z_orb is locked to the assumed nuclear charge**.
Measured: the bond-block h1_no_pk trace deepens -11.00 -> -19.58 Ha as z: 1 -> 2.3
(the Li-core diagonal -4.50 is untouched; the deepening is entirely the bond
orbitals following -z^2/2n^2), while cross-V_ne only moves -6.38 -> -8.08. So
cranking z makes the code believe the H nucleus (physical Z=1) has charge 2.3,
deepening the one-body energy unphysically.

**Structural consequence (strengthens the wall):** orbital contraction is *not a
free variational parameter* in the zero-parameter balanced construction — the
recipe ties orbital exponent to nuclear charge by design. A genuine contracting
bond orbital requires computing h1 as physical kinetic + nuclear integrals
(decoupling exponent from Z_nuc), which is a different, non-zero-parameter
Hamiltonian class. **This is the concrete mechanism behind "healing breaks the
zero-parameter construction"** (v5.14.8): it is not merely that a knob would be
added — the knob cannot even be expressed without changing the construction's
one-body term.

## 3. Cross-check against the corpus (H2, He): the drift is a recipe artifact

The finding is consistent with, and explained by, the accurate systems:

- **He** (Level 3, atom): no bond -> no bond-length observable. Its residual is
  the e-e cusp (a two-body correlation feature), a different axis; done in a
  variational hyperspherical basis, so no fixed-orbital drift.
- **H2** (Level 2, prolate spheroidal natural geometry): 99.97% D_e, no
  balanced-style R_eq drift, *because* its basis exponent is variationally
  optimized (alpha=1.40 optimum; explicit-r_12 for correlation) — the orbital IS
  free to contract. This is the control: where the orbital can contract (H2, He),
  energy and geometry both come out right; where it is frozen and Z-locked
  (balanced/composed LiH), the geometry drifts.

The drift therefore is a price of the cheap composed/balanced recipe, not a law
of nature — the "exact != accurate" / two-kinds-of-sparsity picture
([[two_kinds_of_sparsity]], CHEM-ACCURACY cluster). It worsens with the recipe's
reach: LiH 5.3-8.8% -> BeH2 11.7% -> H2O 19.4% (Paper 17/19).

## 4. Consequence for the path forward (PI-directed 2026-09-20)

The accuracy route for LiH is the **prolate-native (H2-style) treatment** with a
variational, contraction-capable basis — the "unbuilt prolate >=4e two-center CI"
the walls audit repeatedly names. Trick to stay inside the proven regime: **freeze
the Li 1s^2 core** (as balanced already does) -> a 2-valence-electron two-center
problem = the N=2 sweet spot where explicit-r_12 works (H2, and TODAY's HeH+,
v5.14.11). The one genuinely new piece is representing the frozen-core screening
of the valence *exactly* (frozen-density Hartree--exchange in prolate
coordinates), NOT via the crude PK shortcut that was the composed accuracy
bottleneck. Scoping this is the next diagnostic step; if it lands, build from the
`assemble_hetero` / `vne_hetero_mpf` HeH+ engine and measure R_eq + energy against
the drift proven here.

## 5. Honest scope / caveats

- **Grade:** analytic closed-form demonstration (mechanism + limit exact) +
  numerical structural finding (h1 baseline artifact). Not a theorem.
- **1s model vs full:** the 1s-to-Li term gives 6.1 of the ~8.8 pp; the remainder
  is multi-orbital + the second cross term (not separately measured this sprint).
- **n_max=2** throughout (fast, live). The v5.14.8 n_max=3 numbers (0.20% E /
  8.8% R_eq) remain the paper values (documented ~2.3 h/pt coverage gap).
- The point-charge-limit identity is geometry/charge-general (dV/dR -> Z/R^2 as
  z -> inf for any Z, R); the *required* force being Z_A Z_B/R^2 makes the two
  coincide for any diatomic — a general statement, stated for LiH here.

Files: this memo; `debug/balanced_reqdrift_relaxation_probe.py`;
`debug/data/balanced_reqdrift_relaxation.json`.
