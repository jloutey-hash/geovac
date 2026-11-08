# Lithium-7 2²S₁/₂ HFS Roothaan Autopsy — Track 5

**Date:** 2026-05-18
**Sprint:** Multi-track Roothaan-autopsy continuation, Track 5
**Status:** CLOSED with substantive structural finding

## Reference

$\nu_\text{HFS}(^7\text{Li},\, 2^2S_{1/2},\, F=2 \leftrightarrow F=1) = 803{,}504{,}086.0(34)$ Hz
(Beckmann, Boklen, Elke 1974; confirmed Walls 2003).
Nuclear: $I=3/2$, $g_{^7\text{Li}} = 2.170951(4)$ (CODATA $\mu/I/\mu_N$),
$m_{^7\text{Li}} \approx 12{,}786.4\, m_e$.

## Architecture under test

First **multi-electron HFS** in the catalogue at $I=3/2$, exercising
THREE projection-dictionary entries SIMULTANEOUSLY:

- **§III.18** nuclear magnetization-density / Zemach ($r_Z(^7\text{Li}) \approx 3.71$ fm)
- **§III.20** Phillips-Kleinman / core-valence orthogonality ($1s^2$ core screens $2s$ valence)
- **§III.22** bipolar harmonic / Drake combining (multi-electron Fermi contact via 1s-2s-2p CI mixing)

Previous catalogue entries are all single-electron HFS observables
(H 21cm, μH 1S HFS, D 1S HFS, Mu 1S HFS, Ps 1S HFS); Li 7 is the
first 3-electron HFS.

## Five components

1. **Bohr-Fermi** at $I=3/2$, $J=1/2$ with hydrogenic $|\psi_{2s}(0)|^2 = Z_\text{eff}^3 / (8\pi)$
   at $Z_\text{eff} = 1.279$ (Clementi-Raimondi 1963).
   Multiplicity factor $\Delta\nu(F=2 \leftrightarrow F=1) = 2 A_\text{hf}$
   (worked out from $E_F = A/2 [F(F+1) - I(I+1) - J(J+1)]$:
   $E_2 = 3A/4$, $E_1 = -5A/4$, $\Delta E = 2A$).
2. **Schwinger $a_e \approx \alpha/(2\pi)$**.
3. **Reduced-mass / cross-register recoil $(1 + m_e/m_{^7\text{Li}})^{-3}$**
   via §III.14 rest-mass projection.
4. **§III.18 Zemach** with $r_Z(^7\text{Li}) = 3.71$ fm at operator level.
5. **§III.20 PK + §III.22 bipolar harmonic** — DOCUMENTED, not numerically applied.

## Results

| Component | $\nu$ [MHz] | Residual |
| --- | ---: | ---: |
| 1. Bohr-Fermi (Z_eff=1.279, hydrogenic 2s, I=3/2 mult 2) | 82.977 | -89.673% |
| 2. + Schwinger $a_e$ | 83.073 | -89.661% |
| 3. + reduced-mass recoil | 83.054 | -89.664% |
| 4. + §III.18 Zemach ($r_Z=3.71$ fm, Z_eff-scaled) | 83.039 | -89.665% |
| **Experimental** | **803.504** | — |

The framework-native chain (BF + Schwinger + recoil + Zemach) sits at
**$-89.7\%$**. The dominant gap is not 1% — it is a factor of ~9.7×.

## Headline structural finding

**The Li-7 HFS observable exposes a multi-electron contact-density cliff
that is an order of magnitude deeper than the 1.10% framework atomic
accuracy floor (Paper FCI-A).**

The hydrogenic-$Z_\text{eff}$ approximation underestimates $|\psi_{2s}(0)|^2$
by a factor of $\sim 10\times$ for Li 2s valence — far outside any
correction this catalogue has surfaced before. The cliff is a structural
property of the alkali valence-orbital contact density, not a literature-
convention exposure.

**Mechanism (well-documented in atomic physics, here mapped to GeoVac
projections):**

- **§III.20 Phillips-Kleinman / core-valence orthogonality:**
  the 2s wavefunction must satisfy $\langle 2s | 1s_\text{core} \rangle = 0$
  exactly; this constraint pulls the 2s radial node inward and dramatically
  enhances the 2s amplitude at the origin. Clementi-Roetti 1974 SCF
  $|\psi_{2s}(0)|^2$ is enhanced $\sim 8\times$ over the hydrogenic
  $Z_\text{eff} = 1.279$ approximation.
- **§III.22 bipolar harmonic / Drake combining:** multi-electron Fermi
  contact via 1s-2s-2p CI mixing. The unpaired 2s electron's spin density
  at the nucleus is modified by the 1s² core via spin-spin and core-
  polarization couplings (Lindgren 1985; Yan & Drake 2003); contributes
  ~1-2% on top of the PK enhancement.

The Z_eff-sensitivity grid in the driver confirms this:
- $Z_\text{eff} = 1.279$ (CR67): $\nu = 83.0$ MHz ($-89.7\%$)
- $Z_\text{eff} = 3.000$ (no screening): $\nu = 1070.8$ MHz ($+33.3\%$)
- Experimental $\nu = 803.5$ MHz corresponds to effective $Z_\text{eff} \approx 2.73$.

The effective $Z_\text{eff}$ for Li 2s contact density sits closer to 2.7
than to 1.3 — the PK orthogonality enhancement pulls 2s amplitude back
toward the unscreened-nuclear value at the origin, while the orbital's
spatial extent remains screened. This is the well-known *Fermi contact
enhancement* in alkali atoms, which standard atomic-physics frameworks
handle via SCF + CI but the GeoVac framework currently does not.

## Operator-level verifications (machine precision)

- **$\hat{I}\cdot\hat{S}$ at $I=3/2$, $J=1/2$:** eigenvalues $\{+3/4, -5/4\}$,
  multiplicity (max-min) = 2.0 to machine precision (target: 2.0 by CG).
- **Pauli encoding:** $Q_\text{total} = 3$ qubits (2 nuclear binary for $I=3/2$
  in 4-dim space + 1 electronic), 8 non-identity Pauli terms — same minimal
  structure as D HFS (10 terms at $I=1$) and H 21cm (3 terms at $I=1/2$).
- **§III.18 operator I-independence:** Pauli count = 4 (same as H/D); confirms
  §III.18 depends only on spatial register, not nuclear-spin qubits.

## Convention exposure

Unlike the D HFS Layer-2 itemization split (§V.D.1 Eides vs PY 2010,
~25 mfm in $r_Z(p)$), **no major Li-7 HFS literature-convention mismatch
exists**. The experimental value Beckmann-Boklen-Elke 1974 = 803.504086(34) MHz
is unchallenged; Yan & Drake 2003 ab initio MCDHF reproduces it to ~0.001%
using $r_Z(^7\text{Li}) = 3.71$ fm.

## New §V.D-class observation (class iii)

This autopsy surfaces a NEW §V.D-class exposure: **multi-electron
contact-density cliff**. Distinct from the four existing §V.D entries
(literature convention mismatches Eides/PY/Karshenboim/Antognini):

- §V.D.1 D HFS (Eides Tab 7.3 vs PY 2010): $\sim 25$ mfm in $r_Z(p)$
- §V.D.2 μH VP (Antognini vs Krauth): $\sim 100$ ppm
- §V.D.3 μH 21cm (Karshenboim recoil aggregation)
- §V.D.4 Friar profile (RMS vs first-moment, $4/\pi$ factor)
- **§V.D.5 (NEW, this autopsy): multi-electron contact-density cliff**
  — for any alkali ns valence orbital, hydrogenic-$Z_\text{eff}$ misses
  $|\psi_{ns}(0)|^2$ by $\sim 5\text{--}10\times$ unless PK orthogonality is
  explicitly enforced. Magnitude: $\sim 10\times$ on $A_\text{hf}$ for Li-7;
  scales similarly for Na, K, Rb, Cs valence HFS.

This is a class (iii) exposure: not literature convention, not framework
PRECISION ceiling (the 1.10% atomic accuracy floor), but framework
ARCHITECTURE ceiling — the choice of hydrogenic-$Z_\text{eff}$ basis at
the first-pass level cannot represent core-valence orthogonality at the
amplitude-at-origin level.

## Coverage matrix extension

New axis added: **electron count** (single vs multi-electron HFS).

| System | $I$ | $N_e$ | Residual |
| --- | --- | --- | ---: |
| H 21cm | 1/2 | 1 | +18 ppm |
| D 1S HFS | 1 | 1 | +286 ppm |
| μH 1S HFS | 1/2 | 1 | +2 ppm |
| Mu 1S HFS | 1/2 | 1 | +199 ppm |
| Ps 1S HFS | 1/2 | 0 | +4900 ppm |
| Mu 1S-2S | 1/2 | 1 | -0.11 ppm |
| Mu Lamb | 1/2 | 1 | +130 ppm |
| μH Lamb | 1/2 | 1 | +1000 ppm |
| **Li-7 2²S HFS (this autopsy)** | **3/2** | **3** | **-89.7%** |

The catalogue now spans five orthogonal axes (mass hierarchy / nuclear
spin / observable type / QCD content / electron count). Li-7 adds the
first multi-electron HFS entry and exposes the multi-electron
contact-density cliff.

## Named follow-on sprint

**Hylleraas-Eckart double-zeta for Li 2²S HFS:**
- Module: extend `geovac/hylleraas_r12.py` from He singlet ground state
  (currently +0.0006% at $\omega=4$) to Li 2-electron doublet + 1s-pair
  reference, then to 3-electron Hylleraas with $r_{12}, r_{13}, r_{23}$.
- Expected improvement: framework Li atomic accuracy 1.10% → sub-0.1%
  (Drake-class), AND closure of the $|\psi_{2s}(0)|^2$ enhancement at the
  same time (Hylleraas explicitly satisfies the PK orthogonality at the
  amplitude level, no separate Z_eff fitting required).
- Downstream impact: unlocks separable testing of §III.18 (Zemach via
  $r_Z$), §III.20 (PK orthogonality), §III.22 (bipolar core polarization),
  and multi-loop QED in multi-electron HFS. The alkali HFS frontier
  (Li, Na, K, Rb, Cs) becomes a testbed for the framework's projection-
  chain dictionary at $Z \in [3, 55]$.

## Connection to the Cs HFS Z>20 cliff

The §V.C.6 Cs HFS Track-2 diagnostic identified that the heavy-atom
screening cliff is structurally due to **missing radial nodes in the
all-positive-coefficient single-zeta hydrogenic basis** (`debug/z_cliff_diagnostic_track2_memo.md`).
The CR67 single-zeta for Cs 6s peaks at $r \approx 0.63$ bohr instead of
physical $r \approx 2$ bohr; the same mechanism with the same direction
operates for Li 2s — but the magnitude is much larger (factor $\sim 10\times$
on contact density vs factor $\sim 1.6\times$ on Cs HFS).

**Two cliffs, one mechanism:** the Li-7 HFS multi-electron contact-density
cliff (this autopsy) and the Cs HFS Z>20 cliff (§V.C.6) are both
manifestations of the same structural limitation — single-zeta hydrogenic
basis cannot represent the PK-orthogonalized SCF radial wavefunction.
The Z>20 cliff diagnostic (BBB93 multi-zeta closes ~20% of Cs residual)
and the Hylleraas-Eckart closure named here are two routes to the same
structural fix.

## Verdict

- §III.18 magnetization-density operator at $I=3/2$ reproduces leading-order
  Eides Zemach at hydrogenic $Z=1$ to machine precision; profile
  independence preserved; Z_eff scaling applied for Li 2s valence.
- $\hat{I}\cdot\hat{S}$ Hamiltonian operator at $I=3/2$ has the correct CG
  eigenvalues $\{+3/4, -5/4\}$ with multiplicity 2 at machine precision.
- **Multi-electron contact-density cliff** identified as the dominant
  Layer-2 contribution: $\sim 10\times$ on $A_\text{hf}$, beyond what
  hydrogenic-$Z_\text{eff}$ can capture. The §III.20 PK enhancement and
  §III.22 bipolar core polarization together account for the cliff.
- §III.18 cannot be separately tested in multi-electron HFS until §III.20
  is engineered correctly in the framework (Hylleraas-Eckart route).
- This is a structural-skeleton-scope finding for multi-electron HFS:
  the framework provides the skeleton (selection rules, transcendental
  class, Pauli encoding); the multi-electron *amplitude-at-origin*
  requires explicit core-valence orthogonality that the single-zeta
  hydrogenic basis cannot supply.

## Files

- `debug/Li7_HFS_autopsy_track5.py` (~700 lines)
- `debug/Li7_HFS_autopsy_track5_memo.md` (this file)
- `debug/data/Li7_HFS_autopsy_track5.json` (~14 KB structured outputs)
- Paper 34 §V.C subsection added (`sec:autopsy_li7_hfs`)
- Paper 34 §V.B row added (Li-7 HFS off-precision, error code = framework
  precision cliff / multi-electron architecture wall)
