# NaH basis extension — does the in-basis FCI ceiling rise?

**2026-08-22, exploratory (PI: "wire that up as a debug for now and see where we
land"). No paper claim made. Drivers: `debug/noci_nah_basis_extension.py`,
`debug/noci_nah_b2_diagnostic.py`.**

## Question

Paper 58 attributes NaH's D_e shortfall to minimal-basis incompleteness, naming
three missing ingredients: polarization, diffuse functions, BSSE correction.
None had been tried. NOCI-3 already recovers 91.1% of in-basis FCI, so NOCI is
not the constraint — the basis is. The number to move is the **in-basis FCI
ceiling**: published 1.175 eV against experiment 1.961 eV.

## Answer: yes, but modestly — and two prior numbers were flattered

All values all-electron FCI (no frozen core), experiment = 1.961 eV @ 3.566 a0.

| basis | M | raw D_e | CP D_e | BSSE | CP % exp | CP R_eq |
|---|---|---|---|---|---|---|
| base | 7 | 1.173 @ 3.58 | **1.092** @ 3.73 | 0.081 | 55.7% | +4.6% |
| + diffuse s (ζ=0.30) | 8 | 1.286 @ 3.65 | **1.181** @ 3.80 | 0.104 | 60.2% | +6.6% |
| + diffuse s (ζ=0.75) | 8 | 1.471 @ 3.61 | **1.290** @ 3.77 | 0.181 | 65.8% | +5.7% |
| + H 2pz | 8 | 1.645 @ 3.39 | **1.207** @ 3.70 | 0.438 | 61.6% | +3.8% |

Counterpoise-corrected, ONE added function buys **+0.11 to +0.20 eV** against a
0.868 eV deficit — 13–23%. Real, in the right direction, not transformative. At
~0.15 eV per function this is a basis-saturation program, not a quick win.

NOCI compactness reproduced: **91.3%** of in-basis FCI at M=7 (published 91.1%).

## Two Paper 58 numbers are BSSE-flattered

Both are raw (atom references in atom-only bases). Counterpoise-corrected,
all-electron, verified on a fine R grid (0.05 a0 spacing, 3.30–4.10):

| | published (raw) | counterpoise-corrected |
|---|---|---|
| in-basis FCI R_eq | +0.8% | **+4.6%** (3.583 → 3.729 a0) |
| in-basis FCI D_e | 1.175 eV (59.9%) | **1.092 eV (55.7%)** |

The paper discloses "no basis-set-superposition correction, the atom references
being computed in atom-only bases" — but attaches it to the **D_e** discussion.
The R_eq claim ("the in-basis FCI landing at +0.8%") carries no such caveat, and
it depends on the same thing: BSSE shortens the bond artificially, so the close
geometric agreement is partly borrowing. **PI decision owed** on whether to
scope that sentence.

## Trap 1 — frozen core by orbital INDEX is unsafe (the ladder it produced is void)

Freezing the 5 lowest *basis-function-indexed* Löwdin orbitals validated at
M=7 (0.21% of D_e) and was then applied at M=8–14, where it was never
validated. Discriminating test on one basis (base + H 2pz, M=8):

| | D_e | R_eq |
|---|---|---|
| all-electron | 1.645 eV | 3.39 |
| frozen core | 0.447 eV | 5.16 |

**1.2 eV apart.** Mechanism: Löwdin orbital 4 derives from Na 2pz; H 2pz sits on
the same axis with the same symmetry and mixes strongly at short R. Freezing it
removes the bonding flexibility and pushes the minimum outward. Freezing by
index is only safe when no added function shares symmetry with a core function.

Everything in the `B1`–`B5` rungs of `debug/data/noci_nah_basis_extension.json`
used it. **Treat that file's ladder as void.** A correct frozen core selects by
orbital *energy*, not index.

## Trap 2 — BSSE grows with basis size, so uncorrected ladders cannot answer this

BSSE went 0.081 eV (M=7) → 0.438 eV (M=8 + H 2pz), a 5× rise from one function.
Raw, adding H 2pz looks like +0.47 eV (1.173 → 1.645). Counterpoise-corrected it
is +0.11 eV. **Roughly three quarters of the apparent gain was borrowing.**
Since basis size is the axis under test and BSSE grows along it, an uncorrected
D_e ladder measures mostly its own artifact.

Counterpoise (Boys–Bernardi ghost-basis references at the same R) is implemented
in `atom_energy_cp` / `scan_cp`.

## A rejected explanation, recorded

The first collapse was attributed to linear dependence — the ζ=0.75 "diffuse"
function overlaps the H 1s by 0.9695, which looks damning. It is not: the
smallest overlap eigenvalue stays at 1.5e-2 (and 4.0e-1 for the H 2pz basis),
three to five orders above the ill-conditioning threshold. The Löwdin transform
was well conditioned throughout and the energies were variational. **A 97%
pair overlap is not by itself ill-conditioning** — check the spectrum, not the
pair.

## Where it stops

M=9+ all-electron needs a sparse/Davidson FCI (the dense solver caps near 2000
determinants) or an energy-ordered frozen core. Both are real builds. Not
started.
