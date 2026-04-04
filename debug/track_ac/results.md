# Track AC: Channel Diagnostic — Even-Odd Staircase in Level 4 l_max Convergence

## VERDICT: PHYSICS (Selection Rule Effect), NOT a Bug

The even-odd staircase in D_e convergence is a **selection rule effect** intrinsic to
the sigma channel structure, amplified by frozen pi channels but not caused by them.

---

## Part 1: Channel Diagnostic

### Channel Counts by l_max (homonuclear, gerade l1+l2=even)

| l_max | N_sigma | N_pi | N_delta | N_total (m_max=2) | N_frozen (AA config) |
|------:|--------:|-----:|--------:|-------------------:|---------------------:|
|     1 |       2 |    2 |       0 |                  4 |                    4 |
|     2 |       5 |    4 |       2 |                 11 |                    9 |
|     3 |       8 |   10 |       4 |                 22 |                   12 |
|     4 |      13 |   16 |      10 |                 39 |                   17 |
|     5 |      18 |   26 |      16 |                 60 |                   22 |
|     6 |      25 |   36 |      26 |                 87 |                   29 |
|     7 |      32 |   50 |      36 |                118 |                   36 |

### New Sigma Channels at Each l_max

| l_max | New channels | Count | Key feature |
|------:|:------------|------:|:------------|
|     1 | (0,0), (1,1) | 2 | Ground + first mixed |
|     2 | (0,2), (2,0), (2,2) | 3 | First (0,L) and (L,0) — strong nuclear coupling to (0,0) |
|     3 | (1,3), (3,1), (3,3) | 3 | Cross-terms only, weaker coupling |
|     4 | (0,4), (4,0), (2,4), (4,2), (4,4) | 5 | New (0,L) and (L,0) + cross-terms |
|     5 | (1,5), (5,1), (3,5), (5,3), (5,5) | 5 | Cross-terms only |
|     6 | (0,6), (6,0), (2,6), (6,2), (4,6), (6,4), (6,6) | 7 | New (0,L) and (L,0) + cross-terms |

**Root cause of staircase:** The nuclear coupling is a ONE-ELECTRON operator: it
is diagonal in one electron's angular momentum and changes the other via Gaunt
integrals. Therefore:
- (0,L) couples directly to (0,0) via electron 2 (l2: L -> 0, l1: 0 = 0 diagonal)
- (L,0) couples directly to (0,0) via electron 1 (l1: L -> 0, l2: 0 = 0 diagonal)
- (l1,l2) with BOTH l1 > 0 and l2 > 0 CANNOT couple directly to (0,0)

Channels (0,L) and (L,0) appear ONLY at even l_max (since l1+l2 must be even for
gerade symmetry, and 0+L=L must be even). At odd l_max, all new channels have
both l1 >= 1 and l2 >= 1 (e.g., (1,3), (3,1), (3,3)) and couple to (0,0) only
INDIRECTLY through intermediate channels. This indirect coupling is much weaker,
explaining the near-zero gain at odd l_max steps.

### Sigma-Pi Coupling: DECOUPLED in Angular Equation

Nuclear coupling selection rule: diagonal in (m1, m2).
- sigma (m1=0, m2=0) couples ONLY to sigma
- pi (m1=1, m2=-1) couples ONLY to pi

V_ee coupling: also diagonal in (m1, m2).

Therefore sigma and pi are **completely decoupled** in the angular eigenvalue problem.
They interact only through the hyperradial 2D tensor product solver.

### Gerade Constraint

All channels at all l_max values satisfy l1+l2 = even. No violations found.

---

## Part 1e: Sigma-Only vs Frozen-Pi Comparison (2D solver)

| l_max | N_ch_sigma | D_e% sigma-only | D_e% frozen-pi (AA) | Diff (pi contribution) |
|------:|-----------:|----------------:|--------------------:|-----------------------:|
|     2 |          5 |           79.78 |               86.44 |                  +6.66 |
|     3 |          8 |           80.43 |               87.09 |                  +6.66 |
|     4 |         13 |           87.31 |               93.87 |                  +6.56 |
|     5 |         18 |           87.39 |               93.95 |                  +6.56 |
|     6 |         25 |        (running) |               95.64 |              (running) |

**Key findings:**

1. Pi channels contribute a **constant ~6.6 pp** to D_e across all l_max values.
   This is because pi is frozen at l_max_per_m[1]=2 (always 4 pi channels).

2. The even-odd staircase is **identical** in sigma-only:
   - l_max 2->3: +0.65 pp (sigma-only) vs +0.66 pp (frozen-pi)
   - l_max 3->4: +6.88 pp (sigma-only) vs +6.78 pp (frozen-pi)
   - l_max 4->5: +0.08 pp (sigma-only) vs +0.08 pp (frozen-pi)

3. This proves the staircase is **intrinsic to sigma channel structure**, not from
   frozen pi. The even-l_max jumps come from adding (0,L) and (L,0) channels.

---

## Part 2: Dissociation Limit

### E(R) at Large R (l_max=4, adiabatic, R_e_max=15.0)

| R (bohr) | E(sigma+pi) | E-E_atoms | E(sigma) | E-E_atoms |
|----------:|------------:|----------:|---------:|----------:|
|       5.0 |   -0.941442 |  +0.05856 | -0.940255 | +0.05975 |
|      10.0 |   -0.728866 |  +0.27113 | -0.728711 | +0.27129 |
|      15.0 |   -0.576441 |  +0.42356 | -0.576409 | +0.42359 |
|      20.0 |   -0.381060 |  +0.61894 | -0.380994 | +0.61901 |
|      30.0 |   -0.092335 |  +0.90767 | -0.090147 | +0.90985 |
|      50.0 |   +0.024098 |  +1.02410 | +0.031325 | +1.03133 |

**The solver does NOT converge to E_atoms = -1.0 Ha at large R.**

### Root Cause: R_e_max Too Small

The adiabatic potential U(R_e) minimum shifts to larger R_e as R increases:
- R=1.4: U_min at R_e=1.256 (well within grid)
- R=10: U_min at R_e=6.947 (near grid edge)
- R=30: U_min at R_e=15.000 (AT grid boundary)
- R=50: U_min at R_e=15.000 (AT grid boundary)

The R_e_max=15.0 default truncates the wavefunction for R >= ~10 bohr.
Even with R_e_max scaled to 1.5*R, convergence is very slow — this is a
fundamental limitation of hyperspherical coordinates at dissociation
(two distant atoms require R_e ~ R/sqrt(2)).

### Impact on D_e

Track AA correctly used **E_atoms = -1.0 Ha** (exact analytical value) as the
dissociation reference, NOT E(large R). Therefore the D_e% values in Track AA
are correct. The dissociation limit convergence issue does NOT affect the
reported results.

The adiabatic solver's D_e > 100% at l_max >= 4 is real: the adiabatic
approximation adds ~10 pp of artificial binding compared to the variational 2D
solver (confirmed constant across l_max in Track AA).

---

## Part 3: Delta Channel Test

Single-point test at R=1.4, l_max=4, 2D solver:

| Configuration | N_ch | D_e% | Wall time |
|:-------------|-----:|-----:|----------:|
| sigma-only (m_max=0) | 13 | 87.31 | ~44s/pt |
| sigma + frozen-pi (AA) | 17 | 93.87 | ~44s/pt |
| sigma + pi + delta (m_max=2) | 39 | 94.52 | 384s/pt |
| sigma + full-pi (AA, l_max=6) | 61 | 96.00 | 754s/pt |

Delta channels at l_max=4 contribute **+0.65 pp** beyond frozen-pi (93.87% -> 94.52%).
This is a modest improvement at 8.7x the cost (384s vs 44s). The full-pi result at
l_max=6 (96.00%) is more cost-effective.

l_max=6 with m_max=2 would give 87 channels — estimated 1500+ seconds per point.
Not cost-effective compared to increasing l_max with sigma+pi.

---

## Summary

1. **Staircase: PHYSICS.** The even-odd pattern is an intrinsic selection rule effect.
   Even l_max adds (0,L) and (L,0) sigma channels with strong nuclear coupling to (0,0);
   odd l_max adds only weaker cross-terms. Pi channels are completely decoupled from
   sigma in the angular equation and contribute a constant ~6.6 pp.

2. **Dissociation: KNOWN LIMITATION.** The solver cannot converge E(R) to -1.0 Ha at
   large R because R_e_max=15.0 is too small. This does not affect D_e when using the
   correct E_atoms=-1.0 reference.

3. **Delta channels: MARGINAL.** +0.65 pp at l_max=4 for 8.7x cost increase.
   Higher sigma l_max with full pi is more cost-effective.

4. **No bugs found.** Channel generation is correct, gerade constraint is enforced,
   selection rules are properly implemented.
