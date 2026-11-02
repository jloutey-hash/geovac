# Pre-registered prediction: the hydrogen torus, one exact winding and one leaking winding (TR)

**Registered 2026-11-02** (PI: "Ok. Do that doughnut."), before any scored case is computed.
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** record HR (`docs/predictions/2026-11-01_hydrogen_ring.md`: site vs ring-momentum labels on H6; its
B1-B3 Z2 and locality findings); the torus prior-art memo `debug/lit_scan/prior_art_tori_quantum_numbers_memo.md`.
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope. One atom-centred 1s function per hydrogen; no
single-centre or shared-exponent molecular basis; no graph concatenation.

## Why

The PI asked whether a torus molecule carries more than one set of quantum numbers, used differently. A torus has
two loops: around the hole (toroidal) and around the tube (poloidal). A ring of N small rings of M hydrogens has an
exact N-fold rotation about the axis, so the toroidal winding k is an exact label. The poloidal winding m is not:
bending a straight tube into a torus squeezes the inner side and stretches the outer side, and that difference
mixes m. The predictions say how the two labels behave and how the m leak shrinks as the torus gets thinner (more
rings at fixed spacing, so the tube looks locally like a straight cylinder, where m would be exact).

## Method (fixed before scoring)

- **Driver:** `debug/torus_study.py`, sha256 `7ad717097a28874d5a8ed18fc9ee1f1e97ab14f47ba9144a74272f2d9ad5ba16`.
  It loads `debug/hring_study.py` (HR's driver, sha256 unchanged,
  `afd9102800943aa8c45e4c6718209aad5f4823d9ff45f592e9c7f4cc97534f28`) for shared helpers and replaces only its Z2
  counter, which builds bit masks from numpy int64 indices and overflows past 64 qubits. The replacement uses
  Python integers. Checked before registration: identical Z2 counts to the original on all nine (geometry,
  labelling) pairs of N x M = 6x3, 8x3, 6x4; completes on 12x3 (72 qubits), where the original raises
  OverflowError. HR never exceeded 60 qubits, so its scored counts are unaffected.
- **Geometry:** atom (j, a) at toroidal angle 2 pi j/N, poloidal angle 2 pi a/M (a = 0 on the outer equator);
  R = s / (2 sin(pi/N)), r = d / (2 sin(pi/M)). **Scored family: s = 1.4 A, d = 1.0 A** (the controls used s = 1.2).
- **Basis:** one Slater 1s per H, zeta 1.24, fitted by 3 Gaussians; Lowdin site orbitals (commute with every
  symmetry of the molecule).
- **Labellings:** site; tor (real DFT over j per a); full (real DFT over j and over a).
- **Measures:** leakage shares in the complex (k, m) basis (share of the Frobenius norm of h and of (pq|rs) carried
  by k-breaking entries, and by m-breaking entries per |delta m| mod M); per real labelling, the Z2 count, JW Pauli
  counts at truncations 1e-10 ... 1e-2 and the 1-norm; the Pauli 1-norm of the m-breaking part (a rigorous bound on
  how far dropping it can move any energy), and its share of the full 1-norm.
- **Scored cases** (`python debug/torus_study.py case N M 1.4 1.0 [--light]`; rows to `debug/data/torus_study.jsonl`):
  - M = 3: N = 10, 11, 14 (full labels); N = 18, 22, 28 (`--light`, leakage only).
  - M = 4: N = 10, 11 (full labels); N = 14, 18, 22 (`--light`).
  - r/R spans 0.255 to 0.092 (M = 3) and 0.312 to 0.144 (M = 4); closest inner-side H-H gap 0.96 A.
- **Slopes:** local log-log slope of a leak share against r/R between consecutive scored cases at the same M,
  ordered by decreasing r/R.

## What I know before comparing (disclosed)

Sixteen non-scored control cases at s = 1.2, d = 1.0, in `debug/data/torus_controls.jsonl` (sha256
`ddaa587ad3ecfd463b220896551b0ca699eef8780506141939fc70658d032ad0`): M = 3, N = 8-12, 16, 20, 24, 32; M = 4,
N = 8, 9, 10, 12, 16, 20, 24. Nine have full labels (M = 3: N = 8-12; M = 4: N = 8, 9, 10, 12).

| Measure | Control result |
|:--|:--|
| k-breaking share, max over h, g, all 16 | 8.0e-15 |
| Z2 site/tor/full | 2/4/5 at every even N, 2/3/4 at every odd N (9 of 9) |
| Delta m = 1 local slopes, M = 3 | h 0.675 -> 0.886, g 0.727 -> 0.940, rising every step |
| Delta m = 1 local slopes, M = 4 | h 0.790 -> 0.917, g 0.929 -> 0.997, rising every step |
| Delta m = 2 local slopes, M = 4 | h 1.724 -> 1.865, g 1.779 -> 1.977, rising every step |
| full/tor Pauli count at 1e-10 | M = 3: 0.50649-0.50683; M = 4: 0.53072-0.53090 |
| 1e-2 truncation | site count below both tor and full (9 of 9) |
| m-breaking share of full 1-norm | falls with N at each M (M = 3: 0.137 -> 0.098; M = 4: 0.203 -> 0.131) |

The Z2 counts were derived before the controls (site: the two spin-number parities; tor: plus the k-parity
pieces, one more at even N where k = N/2 is real; full: one more from m); the controls confirmed them, including
both odd-N cases. The slope predictions (1 for |delta m| = 1, 2 for |delta m| = 2) come from expanding the
geometry in r/R: the inner-outer difference first enters at order (r/R)^1 with a cos(theta) shape, which shifts m
by one, and a shift by two needs it twice. The full/tor ratio was noticed on the controls and has no derivation.

## Predictions

- **TR1 (exact label).** The k-breaking share is at most 1e-13 for h and for g in every scored case.
- **TR2 (Z2 counts).** Site/tor/full = 2/4/5 at N = 10 and 14; 2/3/4 at N = 11; for both M where computed.
- **TR3 (leak exponent).** For each M and each of h, g:
  (a) the Delta m = 1 local slopes rise at every step toward thinner tori, and every one lies in [0.60, 1.05];
  (b) for M = 4, the Delta m = 2 local slopes rise at every step, and every one lies in [1.55, 2.10];
  (c) at the thinnest step, the g slope is within 0.15 of 1 (Delta m = 1, M = 3), within 0.10 of 1
  (Delta m = 1, M = 4) and within 0.15 of 2 (Delta m = 2, M = 4).
- **TR4 (counting ratio).** full/tor Pauli count at 1e-10 lies in 0.5066 +- 0.002 (M = 3) and 0.5308 +- 0.002
  (M = 4) in every full-label case. If it is a selection-rule count, a different spacing does not move it.
- **TR5 (locality reverses the order).** At the 1e-2 truncation, the site count is below both tor and full in every
  full-label case (HR-B3 on the torus).
- **TR6 (the leak's energy weight falls).** The m-breaking share of the full 1-norm falls strictly with N at each M
  over the full-label cases.

## Novelty status (set before computing)

- Curvature coupling of toroidal and poloidal modes is known: Encinosa and Williamson, "Electron gases in toroidal
  shells: mode coupling and state functions", arXiv:2207.01745 (Physica E, doi 10.1016/j.physe.2022.115558); the
  cylinder-folded-into-torus reading appears in the plasmonics literature ("broken poloidal symmetry and plasmonic
  eigenmodes on a torus", Monash). TR3 is a check of that mechanism in a molecule, not new physics.
- Momentum labels and translation symmetries for term reduction and qubit tapering are known: Manrique et al.,
  arXiv:2008.08694 (momentum-space UCC, point-group tapering for periodic systems); Picozzi, arXiv:2606.05777
  (periodic symmetry-adapted encoding with translation generators); Bravyi et al., arXiv:1701.08213. TR2 and TR5
  apply them to a molecule with one exact and one approximate winding.
- Not found in these sources: the fixed full/tor count ratio (TR4), and the m-breaking 1-norm as a no-SCF bound on
  the cost of keeping m as a label (TR6). Scope of the scan: three web searches and these abstracts; not
  exhaustive.

## Controls

The sixteen control cases above, at a different ring spacing, set every band. No scored case was run before this
record was committed.

## Decision

- **GO:** TR1, TR2 and TR4 held, and TR3 held for g at both M.
- **BORDERLINE:** TR1 and TR2 held, but TR3 or TR4 missed in part.
- **STOP:** TR1 or TR2 missed (the exact label leaks, or the Z2 derivation is wrong).
TR5 and TR6 are reported but do not enter the decision.

## Reading if GO, stated in advance

A torus molecule carries two winding labels that do different jobs. The winding around the hole is exact and
gives the qubit symmetries. The winding around the tube is approximate; its leak shrinks like (r/R)^|Delta m|,
so it becomes useful as the torus thins, and its cost as a label is bounded without any SCF. This confirms a known
geometric mechanism in a molecule. It does not by itself say anything about gauge groups.

## Disposition (2026-11-02, appended after scoring; predictions above unedited)

**Decision: GO.** All six predictions held. Scored rows `debug/data/torus_study.jsonl` (11 cases, s = 1.4,
d = 1.0; sha256 `c2de2d26ac0006e65055b11169ca16057132ed8d4f377e381bdc8092bb88f3d1`), driver hash as registered.

| Prediction | Result | Verdict |
|:--|:--|:--|
| TR1 k-breaking <= 1e-13 | max 8.4e-15 (h and g, all 11) | HELD |
| TR2 Z2 site/tor/full | M = 3: N = 10, 14 2/4/5, N = 11 2/3/4; M = 4: N = 10 2/4/5, N = 11 2/3/4 | HELD |
| TR3a Delta m = 1 slopes rising, in [0.60, 1.05] | M = 3: h 0.765, 0.799, 0.838, 0.866, 0.887; g 0.834, 0.867, 0.902, 0.926, 0.943. M = 4: h 0.862, 0.883, 0.906, 0.920; g 0.978, 0.989, 0.996, 0.998 | HELD |
| TR3b Delta m = 2 slopes rising, in [1.55, 2.10] | M = 4: h 1.766, 1.803, 1.844, 1.871; g 1.912, 1.942, 1.967, 1.980 | HELD |
| TR3c thinnest g slope | 0.943 (to 1, bar 0.15); 0.998 (to 1, bar 0.10); 1.980 (to 2, bar 0.15) | HELD |
| TR4 full/tor at 1e-10 | M = 3: 0.50657, 0.50649, 0.50643; M = 4: 0.53081, 0.53081 | HELD |
| TR5 site lowest at 1e-2 | M = 3: 2131 / 8235 / 7689, 2542 / 8400 / 8862, 3991 / 14063 / 11957; M = 4: 3681 / 12901 / 11553, 4357 / 12629 / 10701 | HELD (not in decision) |
| TR6 m-breaking 1-norm share falls | M = 3: 0.0964, 0.0878, 0.0714; M = 4: 0.1318, 0.1168 | HELD (not in decision) |

Imaginary residue after dropping the m-breaking part: at most 2.1e-14 (the back-transform is consistent).

**Reading.** As stated in advance: the winding around the hole is exact and carries the extra Z2; the winding
around the tube leaks as (r/R)^|Delta m| in the thin limit, approached from below, and its cost as a label (the
m-breaking 1-norm share) falls as the torus thins. This confirms a known geometric mechanism in a molecule.

**Caps.** One minimal Slater 1s per H, no correlation energy computed; two tube sizes (M = 3, 4); slopes are
local finite differences, so "approaches 1 and 2" is read from the trend, not from a fitted limit. TR4's ratio held
at a second spacing to 3e-4, which supports a counting origin, but no derivation exists; it stays an observation
until one does.
