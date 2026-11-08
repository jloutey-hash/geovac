# Sprint: Chemistry-error projection decomposition (diagnostic)
**Date:** 2026-07-05
**Mode:** diagnostic-before-engineering (PK has 6+ honest negatives; we MODEL the
residual, we do not try to reduce it — a genuinely new move, prior-art-cleared).
**Goal (idea a, primacy):** characterize the chemistry solver's error as projection
content; (idea b, fallback) a corrector. PI worry driving it: "runway on (a)?"

## Prior-art clearance (two read-only sweeps)
NEW as an integrated program. Ingredients scattered (accuracy_target_analysis.md
3-row category table; cusp productionized; overlap-correction prototypes; Track AA
CBS on H2 D_e only; composed-vs-balanced LiH juxtaposed once). No signed
cross-system projection-attributed atlas; no structural model of the PK residual;
§1.8 / Paper 34 §V.C are atomic-precision-only (no chemistry-solver analog). §3
dead-ends all try to REDUCE error, none MODEL it.

## Deliverables
- `debug/data/chem_error_atlas.md` — verified two-column (energy vs R_eq) signed
  atlas, 15 rows, from owning papers. Caught a live doc bug: **H2O composed is
  19.4% (Paper 19), not 26% (stale in CLAUDE.md §5 + accuracy_target_analysis.md).**
- `debug/sprint_pk_amplification_lih.py` + `debug/data/pk_amplification_lih.json`
  — tilt/curvature autopsy, balanced (n=2,3) + composed/PK (l=2,3,4), from cached
  current-convention curves.
- `debug/sprint_pk_amplification_finegrid.py` + `..._finegrid.json` — live n_max=2
  fine grid; validates cache == current solver and gives robust derivatives.

## Method
The "error" is TWO objects. For computed PES E_c(R), residual eps=E_c-E_true:
- energy error  = eps(R_true)             [R_true=3.015 exp]
- geometry error dR = R_eq,c - R_true = -eps'(R_true)/E_c''(R_true)
Since E_true'(R_true)=0, **eps'(R_true) = E_c'(R_true) exactly (reference-free)**,
and E_c'' is offset-immune, so the whole decomposition is convention-proof (does
NOT need a reference PES or the correct absolute-energy convention — which matters,
because the live balanced E_coupled carries an R-independent -7.28 Ha core
double-count). True stiffness from spectroscopy: k = mu*w_e^2 = 0.0659 Ha/bohr^2
(LiH, w_e=1405.65 cm^-1, mu(7Li1H)).

## Results (LiH)
| curve | energy err eps | tilt eps'(3.015) Ha/bohr | E_c'' / k_true | R_eq err |
|-------|---------------:|-------------------------:|---------------:|---------:|
| balanced n=2 | 1.75% | -0.030 (live fine-grid -0.0298) | 2.24x | 6.9% |
| balanced n=3 | 0.20% | -0.036 (coarse 5-pt fit)         | 2.22x | 8.8% |
| composed l=2 | 1.63% (over) | -0.016 | 1.51x | 5.4% |
| composed l=3 | 1.91% (over) | -0.037 | 1.45x | 14.5% |
| composed l=4 | 2.38% (over) | -0.061 | 1.57x | 23.3% |

Consistency -eps'/E_c'' reproduces dR to ~5% (parabolic identity holds).

## Four findings
1. **R_eq error = residual TILT / computed curvature** — governed by the *slope*
   of the error at R_true, not its magnitude. Verified both architectures.
2. **Balanced: the tilt is FROZEN (~-0.03) while energy converges 8.8x (n=2->3).**
   This IS the "energy converges but geometry drifts" paradox: energy probes
   eps(R_true) (heals), geometry probes eps'(R_true) (structural, does not heal).
   [n=3 tilt from a coarse 5-pt fit; the qualitative "stays O(0.03), does not
   trend to 0" is robust, exact value uncertain. n=3 live infeasible: ~2.3 h/pt.]
3. **Composed/PK: the tilt GROWS ~linearly with l_max** (-0.016, -0.037, -0.061 at
   l=2,3,4; ~-0.022/l_max) AND energy diverges (over-binds more). This linearly
   growing tilt is the mechanistic SOURCE of the famous non-saturating +0.15-0.30
   bohr/l_max R_eq drift (guardrail fact): R_eq drift = tilt/curv, curv~const,
   tilt~linear in l_max.
4. **Both wells are TOO STIFF (1.5-2.2x k_true), not soft.** So the original
   "shallow-well amplifies a small error" framing is WRONG-SIGNED at the mechanism
   level; the stiff well SUPPRESSES R_eq error (if true stiffness, LiH balanced
   R_eq err would be ~15-18%). Cross-system atlas trend (floppier=worse) still
   reconciles: R_eq err ~ tilt/(c*k_true), computed curv tracks k_true by a
   ~constant factor, so floppier real systems (small k_true) => bigger R_eq err
   for comparable tilt.

## Direction sign (physical)
Tilt eps'(R_true) < 0 everywhere => residual favors LONGER bonds => R_eq pushed
out. Balanced: single-center-per-block basis resolves separated atoms (large R)
better than the overlapping bonding region (short R). Composed: PK over-repulsion
that worsens as added channels dilute the (0,0) barrier weight.

## Verdict on (a) runway
YES, sharper than the entry hypothesis. Clean structural reframe: chemistry error
= two independent objects; the headline 5-26% is a geometry-TILT defect, largely
orthogonal to the (convergent, QC-relevant) single-geometry electronic accuracy.
Next (make-or-break for a genuine (a) result or a (b) corrector): does the tilt
eps'(R_true) have a zero-parameter structural form? Two concrete leads:
- composed tilt ~linear in l_max/n_channels (already visible: -0.022/l_max) —
  predict slope from PK (0,0)-channel dilution.
- balanced tilt ~ R-derivative of block-block density overlap at R_true — test
  eps'(R_true) vs d/dR<block_A|block_B>(R_true) across first-row hydrides.

## TILT-MODEL step (chase the overlap law) — 2026-07-05
Deliverables: `debug/sprint_tilt_overlap_law.py` + `..._law.json`;
`debug/sprint_tilt_stiffness_honest.py`.

**Overlap-slope law FALSIFIED.** LiH bond block = Li-valence(Z_eff=1.0)+H(Z=1),
both centers hydrogen-like Z=1, so the leading bond overlap is the parameter-free
1s-1s Slater S(R)=e^{-R}(1+R+R^2/3). Tested eps'(R) ∝ S'(R) (eps~S) and ∝ 2SS'
(eps~S^2). Regressions give R^2=0.996-0.9998 — but that is an **intercept-driven
curve-fit artifact** (audit rule): the pointwise ratio eps'(R)/S'(R) runs +0.207 →
0 → -0.065 across the window (114% spread, sign change). Root cause: eps'(R) is
~LINEAR in R while any two-center overlap slope is exponential — a linear residual
cannot be proportional to an exponential. The residual near equilibrium is a
low-order (curvature-mismatch) polynomial, NOT overlap-shaped. Clean negative.

**Stiffness finding #4 CORRECTED.** The "1.5-2.2x too stiff" was evaluated at
R_true=3.015, which sits on the inner repulsive wall (inflates curvature).
Re-evaluated at each curve's OWN minimum (omega_e-relevant point):
- balanced n=2: omega_e=2040 cm^-1, **+45%** (curv 2.1x) — GENUINELY too stiff.
- composed l=2/3/4: omega_e=1639/1479/1409 (+16.6/+5.2/+0.3%) — mildly high,
  DECREASING with l_max; reconciles with Paper 17's reported 1471 (+4.6%). The
  composed "too stiff" was the inner-wall artifact; at its own min it's ~right.

**Architecture-specific geometry defects (corrected picture):**
- balanced: stiff well (omega_e +45%) + moderate outward tilt (R_eq +7-9%).
  Energy 0.20% but BOTH shape derivatives (position, curvature) wrong. New: the
  best-energy solver gives a 45%-too-high vibrational frequency.
- composed: ~correct stiffness, R_eq error is pure outward displacement from the
  PK barrier, growing ~linearly with l_max (the tilt ∝ l_max mechanism).

## Verdict after the tilt step
Diagnosis is rich and the reframe is solid (error = well-SHAPE defect: position +
curvature, architecture-specific; NOT depth — energy converges to 0.20%). But the
one clean structural LAW we chased (overlap-slope) is FALSIFIED. No zero-parameter
tilt model in hand. Remaining leads, both un-cheap or un-clean:
- balanced stiffness universality (is omega_e ~+45% across hydrides?) — needs
  BeH2 balanced (~600k-det FCI, hours/point) or balanced omega_e from papers
  (not currently reported).
- composed displacement ∝ l_max — a drift-rate model, not zero-parameter.
(b)-corrector: fixing balanced geometry needs a 2-handle (curvature rescale +
shift) correction = fitted unless 2.1x and the shift prove universal/derivable.

## Honest caveats
- balanced n=3 tilt is a coarse fit (finding 2 qualitative not quantitative there).
- composed absolute energies from an approximate (non-variational) solver.
- cross-system stiffness/tilt only computed for LiH; universality untested (BeH2
  expensive). pyscf unavailable on Windows (no external reference PES this sprint).

## 6. Honest scope
- **Theorem grade:** nothing new; no proofs.
- **Structural sketch:** the error decomposition itself -- fixed-geometry energy
  (residual magnitude), R_eq (residual slope at the true min), and omega_e
  (curvature) as three independent read-outs of ONE residual. Offset-immune and
  well-defined, but a framing/decomposition, not a theorem.
- **Numerical observation (MEASURED, single-system LiH):** balanced n_max=2
  omega_e ~= 2040 cm^-1 (+45%), curvature 2.11x k_true, R_eq +6.9%; composed
  omega_e at own minimum +4.6-17% (mildly high); balanced energy 0.20% @ n_max=3
  (pre-existing). Backed by tests/test_paper19_well_shape.py (QA-passed).
- **Falsified (clean negative):** the overlap-slope law eps'(R) proportional to
  S'(R). R^2=0.996 was an intercept-driven curve-fit artifact; the pointwise
  ratio test (114% spread, sign change) is the discriminator. Residual is a
  low-order polynomial, not overlap-exponential.
- **Corrected in-sprint:** the "1.5-2.2x too stiff" at R_true was an inner-wall
  evaluation artifact; re-anchored at each curve's own minimum.
- **Named open follow-ons (none launched -- diagnostic banked):** (a) cross-
  hydride universality of the +45% balanced stiffness -- needs BeH2 balanced
  (~600k-det FCI); (b) composed displacement proportional to l_max as a drift-
  rate model (not zero-parameter); (c) transfer of the cusp / two-body conformal
  factor into the qubit representation (known hard wall, TC dead-ends).
