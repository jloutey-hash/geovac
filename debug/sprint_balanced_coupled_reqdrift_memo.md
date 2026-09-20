# Sprint: balanced-coupled LiH R_eq-drift term localization (DIAGNOSTIC)

**Date:** 2026-09-19
**Mode:** diagnostic-only (documented wall area: §3 "Balanced + frozen-core PES",
"Overlap-slope law for the balanced-solver R_eq tilt"). No fix attempted.
**Goal:** localize which Hamiltonian term/approximation drives the balanced-coupled
LiH bond-length (R_eq) drift — the "0.20 % energy / ~9 % R_eq" split, energy
converges while geometry drifts.
**Verdict:** **GO — drift localized to ONE term: the cross-center V_ne
(Shibuya–Wulfman). Its R-slope is the sole electronic restoring force and it is
8.8 % too weak. Every other term is bit-exactly R-independent (force = 0).**

Driver: `debug/sprint_balanced_reqdrift_termdecomp.py`
Data:   `debug/data/balanced_reqdrift_termdecomp.json`, `debug/data/reqdrift_run.log`

---

## 1. Current-state confirmation (not trusting a memo)

**Papers 19 & 20 (authoritative source, re-read this sprint) both state the current
split:** LiH balanced-coupled, n_max=3, analytical integrals:
E(3.015) = −8.055 Ha = **0.20 % energy error**, R_eq = 3.28 bohr = **8.8 % error**,
structural outward drift +0.053 bohr per n_max, crossing *past* exact at n_max=4
(over-binds). Paper 19 already localizes the residual to "the one-particle orbital
basis" via the measured statement: *at the true minimum the computed electronic
energy gradient is 8.8 % weaker than the +Z_A Z_B/R² required to balance nuclear
repulsion*, and its n_max=4 direct-CI (`test_balanced_direct_ci.py`) shows an
irreducible orbital-basis-limit wall (tilt worsens to −0.0399, R_eq → +10.1 %,
curvature 2.2× too stiff).

**I reproduced the split LIVE at n_max=2** (n_max=3 is ~2.3 h/point and has no CI
regression test — a documented deferred coverage gap, CHANGELOG M-B; not re-run):

| quantity | measured (n_max=2) | paper (n_max=2) |
|---|---|---|
| E(R_true) [TC] | −7.928 Ha, **1.77 %** | −7.93, 1.7–1.8 % |
| R_eq (parabolic min) | **3.218 bohr, 6.7 %** | 3.226, 7.0 % |
| tilt ε'(R_true)=dE/dR | **−0.0290 Ha/bohr** | −0.030 |

Match confirms the current code reproduces the "energy converges / geometry drifts"
structure. STOP condition not triggered.

---

## 2. Method — rigorous per-term force decomposition (term freezing)

The 2026-07-05 memo established R_eq error = residual TILT ε'(R_true) / curvature,
and that the tilt is frozen (~−0.03) while energy converges. The well-shape test
split the tilt into dV_NN/dR (exact) + electronic gradient (8.8 % too weak). **Neither
decomposed the electronic gradient BY TERM.** That is this sprint's contribution.

For any term T, define H_frozenT(R) by holding T at T(R_ref=R_true) while all else
varies with R. At R_ref, H_frozenT(R_ref)=H(R_ref) so the ground state ψ is
identical, and Hellmann–Feynman gives the term's exact force contribution:

    F_T = tilt_full − tilt_frozenT = ⟨ψ| dT/dR |ψ⟩,   with   Σ_T F_T = tilt_full.

Implementation: `build_balanced_hamiltonian` returns the component matrices
separately (`h1_no_pk`, `h1_cross_vne`, `eri`, `nuclear_repulsion`); `eri_within`
comes from `build_composed_hamiltonian`; cross-block ERI = eri_balanced − eri_within;
V_NN = nuclear_repulsion − E_core, E_core = −7.2799 (R-independent). Each variant is
re-diagonalized (`coupled_fci_energy`) on an 8-point R grid (2.70–3.80 bohr).

---

## 3. Results

### 3a. R-dependence probe (max|Δ(matrix)| vs R_true, over R = 2.70–3.80)

| term | R-dependence |
|---|---|
| within-block h1 (`h1_no_pk`) | **0.00e+00 — bit-exactly R-independent** |
| within-block ERI (`eri_within`) | **0.00e+00 — bit-exactly R-independent** |
| cross-block ERI (`eri_balanced−eri_within`) | **0.00e+00 — bit-exactly R-independent** |
| cross-center V_ne (`h1_cross_vne`) | 0.04–0.20 Ha (R-dependent) |
| V_NN (3/R) | 0.04–0.21 Ha (R-dependent) |

**Only two operators in the entire balanced Hamiltonian depend on R: the nuclear
repulsion V_NN and the cross-center V_ne.** Everything else is frozen.

The cross-block ERI being *bit-exactly* R-independent is a structural fact, not a
grid artifact: `cross_block_mp2._compute_rk_integrals_cross` evaluates a
**single-center Slater integral** (both radial coordinates from a common origin,
standard 1/r_> multipole; the only "two-center" feature is a different Z_eff per
orbital set). No internuclear separation R enters, so the inter-block
electron–electron repulsion carries **zero** geometry.

### 3b. Force decomposition at R_true (the localization)

    tilt_full = dE/dR|R_true                     = −0.0290 Ha/bohr   (<0 ⇒ outward)
    required electronic force (= +Z_Li Z_H/R²)   = +0.3300 Ha/bohr

| term T | F_T = ⟨ψ|dT/dR|ψ⟩ (Ha/bohr) | role |
|---|---:|---|
| V_NN | **−0.3308** | = exact d(3/R)/dR = −3/R²; outward push |
| **cross-center V_ne** | **+0.3018** | **the SOLE electronic restoring force (inward)** |
| cross-block ERI | −0.0000 | zero (R-independent) |
| within-block ERI | −0.0000 | zero (R-independent) |
| within-block h1 | −0.0000 | zero (R-independent) |
| **Σ_T F_T** | **−0.0290** | **= tilt_full exactly** (consistency) |

**The balanced LiH well is a two-term R-balance: V_NN (out, −0.3308) vs cross-center
V_ne (in, +0.3018).** For a zero net force at R_true the cross-V_ne slope would have
to be +0.3308; it delivers +0.3018 — **a deficit of +0.0290 = 8.8 %**, which IS the
residual outward tilt and thus the entire R_eq drift. The 8.8 % reproduces Paper 19's
independently-measured "electronic gradient 8.8 % too weak" — now attributed to a
single operator.

### 3c. R_eq of each frozen variant (corroboration)

| variant | R_eq (bohr) | note |
|---|---|---|
| FULL | 3.218 (6.7 %) | baseline |
| freeze cross-block ERI | 3.218 (6.7 %) | **identical** — no effect |
| freeze within ERI | 3.218 (6.7 %) | **identical** — no effect |
| freeze within h1 | 3.218 (6.7 %) | **identical** — no effect |
| freeze V_NN | no interior min | removes the repulsive wall's slope |
| freeze cross-center V_ne | no interior min | **removes the entire restoring force** |
| freeze ALL electronic | 3.218 (6.7 %) | identical to FULL ⇒ electronic R-force lives entirely in cross-V_ne |

Freezing the three R-independent terms leaves R_eq bit-unchanged; freezing *either*
of the two R-dependent terms destroys the minimum. This confirms the well is built by,
and only by, the V_NN ↔ cross-V_ne balance.

---

## 4. Verdict and mechanism

**The R_eq drift is localized, with a term-exact force budget, to the cross-center
V_ne (Shibuya–Wulfman) integral.** It supplies 100 % of the electronic restoring
force, and its R-derivative is 8.8 % too weak to balance the exact nuclear-repulsion
slope. No other term contributes to the force at all (all bit-exactly R-independent).

**Why the slope is too weak (energy-accurate but geometry-wrong).** Cross-V_ne is
⟨ψ_A|−Z_B/|r−R_B||ψ_A⟩ evaluated in **fixed hydrogenic orbitals** (Z_orb = 1 on the
bond block, Z = 3 on the core). As R changes, only the operator's R_AB moves; the
orbitals never polarize or contract toward the bond. The *magnitude* of cross-V_ne at
a single geometry is accurate to the basis (⇒ energy converges to 0.20 %), but its
*R-derivative* — which physically encodes how the bonding density redistributes as the
nuclei move (orbital relaxation / bond polarization) — is under-delivered because that
density response is absent from a fixed, un-polarized, too-diffuse Z=1 basis. Hence the
internuclear attractive (Hellmann–Feynman) force is too weak and R_eq drifts outward.
This is precisely Paper 19's "localized to the one-particle orbital basis" and the
2026-07-05 memo's "single-center-per-block basis resolves the separated atoms better
than the overlapping bond region," now pinned to the specific operator carrying the
deficit. It is consistent with the n_max=4 direct-CI wall (adding radial shells at
fixed Z_orb does not heal it — the tilt worsens — because the missing ingredient is
R-adaptive polarization, not more shells).

**Secondary structural finding.** The cross-block electron–electron repulsion is
geometry-blind (single-center Slater form), so the balanced model contains *no*
R-dependent electron–electron physics whatsoever. This does not offer a repair lead:
a true two-center V_ee_cross has d⟨V⟩/dR < 0 (repulsion grows as nuclei approach), so
restoring its R-dependence would push R_eq *further outward*, worsening the deficit —
and adding genuine two-center R-dependence is exactly the Löwdin/sparsity-destroying
route already walled (§3 "sparsity-destroying option", "non-orthogonal fermionic
encoding").

---

## 5. Recommended next step

**No new fix is warranted** — this closes to a documented wall. The drift is not an
integral-evaluation, angular-truncation, or basis-*size* defect; it is the missing
R-relaxation of the fixed-Z_orb hydrogenic basis inside the one term (cross-V_ne) that
provides binding. Closing it requires an R-dependent (polarizing/relaxing) bond-block
orbital, which breaks the zero-parameter fixed-basis construction and re-enters the
PK / Löwdin / non-orthogonal-encoding walls of §3. The diagnostic value is the sharper
statement now available for Paper 19 (should the PI want it): the "8.8 %-too-weak
electronic gradient" is entirely the cross-center V_ne R-slope, and the well is a
strict two-term V_NN ↔ cross-V_ne balance — every other term is exactly R-independent.

## 6. Honest scope / caveats

- **Theorem grade:** none; a numerical force decomposition + structural code reading.
- **Measured, live, n_max=2 (LiH):** tilt −0.0290, R_eq 3.218 (6.7 %), E 1.77 %;
  F_cross_vne = +0.3018 vs required +0.3300 (8.8 % deficit); Σ F_T = tilt_full to
  1e-4. Force decomposition is exact (Hellmann–Feynman) and consistent.
- **n_max=3 (0.20 %/8.8 %) taken from Papers 19/20 + the cached curve** (live n_max=3
  ~2.3 h/pt, no CI test — deferred coverage gap). The localization is basis-structural
  and the tilt is frozen/worsening with n_max, so it transfers; but the term-force
  numbers above are the n_max=2 values.
- **BeH₂ / H₂O (context, not recomputed):** balanced polyatomic PES is infeasible
  (~600 k-det FCI, hours/point). The composed/PK R_eq errors are 11.7 % (BeH₂) and
  19.4 % (H₂O) (Paper 17), a *different* architecture whose drift is the PK
  (0,0)-channel dilution growing ∝ l_max (2026-07-05 memo finding 3), not the
  cross-V_ne slope diagnosed here for balanced LiH.
