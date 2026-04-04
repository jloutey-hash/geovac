# Track BK: Orbital Filling Order (Madelung Rule) Investigation

## Summary

**Result: NEGATIVE.** The l-dependent PK pseudopotential does not reproduce
the correct s-p orbital energy splitting. PK pushes s-orbitals UP (less
bound), while in real atoms s-orbitals are MORE bound than p due to core
penetration. The effect has the wrong sign and is orders of magnitude too
large.

## Background

In hydrogenic atoms, orbitals with the same principal quantum number n are
l-degenerate: E(ns) = E(np) = E(nd) = -Z_eff^2/(2n^2). In multi-electron
atoms, electron-electron screening breaks this degeneracy. The empirical
Madelung rule (n+l ordering) predicts filling order.

The GeoVac composed architecture uses a Phillips-Kleinman (PK) pseudopotential
V_PK(r) = A * exp(-B*r^2) / r^2 to enforce core-valence orthogonality. With
the delta_{l,0} prescription (only l=0 orbitals get PK), this creates an
l-dependent energy shift. The question: does this produce physically correct
s-p splitting?

## Computed PK s-p Splittings

PK parameters use Z^2-scaled He-like defaults from composed_qubit.py:
- Li (Z=3): A=6.93, B=7.00 (Paper 17 Table 1)
- C (Z=6): A=27.72, B=28.00
- N (Z=7): A=37.73, B=38.11
- O (Z=8): A=49.28, B=49.78

### n=2 orbital energies

| Atom | Z_eff | E_hydro(n=2) | <2s\|V_PK\|2s> | GV E(2s) | GV E(2p) | GV gap | HF gap | Sign |
|------|-------|-------------|-----------------|----------|----------|--------|--------|------|
| Li   | 1     | -0.125 Ha   | 0.716 Ha        | +0.591   | -0.125   | -0.716 | +0.067 | WRONG |
| C    | 4     | -2.000 Ha   | 67.77 Ha        | +65.77   | -2.000   | -67.77 | +0.272 | WRONG |
| N    | 5     | -3.125 Ha   | 148.3 Ha        | +145.1   | -3.125   | -148.3 | +0.372 | WRONG |
| O    | 6     | -4.500 Ha   | 284.2 Ha        | +279.7   | -4.500   | -284.2 | +0.612 | WRONG |

Convention: GV gap = E(2p) - E(2s). Positive means 2p is higher (less bound).
HF gap is always positive (2p above 2s). GV gap is always negative (2s pushed
above 2p by PK).

### Key observations

1. **Wrong sign.** PK pushes s-orbitals UP because it is a repulsive barrier.
   In real atoms, s-orbitals are MORE bound (lower energy) than p because
   they penetrate the core and see more nuclear charge. The PK effect is
   opposite to the physics.

2. **Wrong magnitude.** The PK matrix elements are enormous (67-284 Ha for
   C-O), while the physical s-p gap is 0.07-0.61 Ha. Even if the sign were
   correct, the magnitude is off by 2-3 orders of magnitude.

3. **Rapid Z_eff scaling.** <2s|V_PK|2s> scales roughly as Z_eff^5 (from
   Z_eff^2 in PK parameters times Z_eff^3 in wavefunction compactness).
   This is far steeper than the physical gap scaling (~Z_eff).

### n=3 splittings (Madelung-relevant)

| Atom | <3s\|V_PK\|3s> | <3p\|V_PK\|3p> | <3d\|V_PK\|3d> |
|------|-----------------|-----------------|-----------------|
| N    | 43.59 Ha        | 1.365 Ha        | 0.008 Ha        |
| O    | 83.50 Ha        | 2.854 Ha        | 0.019 Ha        |

The l-dependence of PK matrix elements (s >> p >> d) is qualitatively
correct (higher l means less core penetration, so less PK overlap).
But this ordering is irrelevant because the overall effect pushes
s-orbitals UP, not down.

### Cross-n PK coupling

PK has significant off-diagonal elements: for N, <2s|V_PK|3s> = 80.4 Ha,
comparable to the diagonal elements. This means PK mixes different n values,
which is a separate (and potentially useful) effect for the composed
architecture. But it does not help with s-p splitting.

## Why PK Cannot Reproduce s-p Splitting

The physics is clear:

1. **Core penetration effect (real physics):** s-orbitals have nonzero
   density at the nucleus. Near r=0, the effective nuclear charge is Z
   (unscreened), not Z_eff = Z - N_core. This extra attraction pulls s
   below p. The effect requires r-dependent Z_eff(r): Z_eff -> Z as r -> 0,
   Z_eff -> Z - N_core as r -> infinity.

2. **PK effect (GeoVac):** The PK barrier V_PK ~ A*exp(-Br^2)/r^2 is
   concentrated near r=0, exactly where s-orbitals have the most density.
   But PK is REPULSIVE (positive), pushing s UP rather than pulling it DOWN.

3. **Structural impossibility:** PK enforces core-valence orthogonality by
   RAISING the energy of valence orbitals that overlap with the core. This
   is the opposite of the physical mechanism that makes s-orbitals more
   bound. No adjustment of PK parameters can fix this — the sign is wrong
   by construction.

## What WOULD Be Needed

To reproduce the correct s-p splitting within the composed architecture,
one would need:

1. **r-dependent Z_eff(r).** The `algebraic_zeff.py` module already computes
   Z_eff(r) from the core density. If the valence Hamiltonian used Z_eff(r)
   instead of the constant Z_eff = Z - N_core, the s-orbitals would see
   more attraction near the nucleus, naturally producing s < p.

2. **Explicit core-valence exchange.** The exchange integral K between core
   and valence orbitals contributes a negative (stabilizing) correction that
   is l-dependent. This is exactly what the LCAO approach captured (pre-v1.0)
   but the composed architecture replaces with PK.

3. **Self-consistent field.** A mean-field treatment where each electron sees
   the effective potential of all others (like HF) naturally produces the
   correct orbital ordering. The composed architecture's ab initio approach
   trades this for simplicity.

## mu_free and Filling Order

Paper 16 defines mu_free = 2(N-2)^2 as the Pauli centrifugal cost for N
electrons. This is a total-system invariant from the SO(3N) Casimir eigenvalue
on S^{3N-1}. It depends only on total electron count N and the spatial
S_N irrep (determined by spin state), NOT on how those electrons fill
individual (n,l) orbitals. Therefore mu_free contains no information about
orbital filling order.

This is structurally inevitable: mu_free characterizes the free (unperturbed)
angular kinetic energy of the N-electron system, before any inter-electron
or electron-nuclear interaction. The filling order is a property of the
INTERACTING system.

## Conclusion

The l-dependent PK pseudopotential in the GeoVac composed architecture
does not capture orbital filling rules (Madelung rule). This is a structural
limitation: PK is a repulsive barrier that pushes s-orbitals up, while the
physical effect (core penetration) pulls them down. The framework correctly
captures orbital degeneracy (from the hydrogenic basis) and its breaking
(from PK + Z_eff screening), but the splitting has the wrong sign.

This does not affect the framework's value for quantum simulation — the
composed architecture targets structural sparsity of the qubit Hamiltonian,
not orbital energy ordering. For the systems currently studied (LiH, BeH2,
H2O), the valence basis uses Z_eff = Z - N_core and PK to approximate the
correct physics. The PK barrier is too large for single-particle orbital
energies but its role in the many-body Hamiltonian (enforcing orthogonality)
is different from reproducing orbital energies.
