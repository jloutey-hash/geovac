# T7 Memo: gamma = sqrt(1-(Za)^2) Radial Corrections

**Track:** T7 (Tier 3 Dirac sprint)
**Date:** 2026-04-15
**Status:** Complete with known limitations.

## 1. Formulas implemented

### s = -1 (all states, any n_r)

Derived via the Hellmann-Feynman theorem applied to the exact Dirac-Coulomb energy

    E = (1/alpha^2) * [(n_r + gamma) / N_D - 1]

differentiating w.r.t. Z (total derivative, accounting for gamma's Z-dependence
through dg/dZ = -Z*alpha^2/gamma):

    <r^{-1}>_{n,kappa} = Z * (gamma*n_r + kappa^2) / (gamma * N_D^3)

where gamma = sqrt(kappa^2 - Z^2*alpha^2), N_D^2 = n_r^2 + 2*n_r*gamma + kappa^2.

**NR limit:** gamma -> |kappa|, N_D -> n. Numerator -> |kappa|*n_r + kappa^2 = |kappa|*n.
Denominator -> |kappa|*n^3. Result: Z/n^2. Matches T1.

### s = -2, -3 (n_r = 0 states only)

For n_r = 0, both large and small Dirac radial components have the same functional
form r^gamma * exp(-lambda*r) (single confluent hypergeometric term). The
expectation values reduce to Pochhammer ratios:

    <r^{-s}>_{n_r=0} = (2*Z/|kappa|)^s / [(2gamma)(2gamma-1)...(2gamma-s+1)]

Explicitly:

    <r^{-2}> = 4*Z^2 / (kappa^2 * 2*gamma * (2gamma-1))
    <r^{-3}> = 8*Z^3 / (|kappa|^3 * 2*gamma * (2gamma-1) * (2gamma-2))

**NR limits verified:**
- <r^{-2}> -> Z^2/(n^3*(l+1/2)) for l = kappa_to_l(kappa). Confirmed for 1s, 2p3/2, 3d5/2.
- <r^{-3}> -> Z^3/(n^3*l*(l+1/2)*(l+1)). Confirmed for 2p3/2, 3d5/2.

### s = -3, l = 0

Raises ValueError (divergence), matching the non-relativistic behavior.

### s = 1, 2

Non-relativistic hydrogenic forms (relativistic corrections are O(alpha^2) and
not implemented at this stage).

## 2. Limitations

**n_r >= 1 with s = -2, -3:** raises NotImplementedError. The Pochhammer formula
does not apply because the Dirac radial wavefunction involves a confluent
hypergeometric function (Laguerre polynomial of degree n_r). The exact closed forms
require the Kramers-Pasternak three-term recursion with verified coefficients.

Multiple attempts to derive the recursion coefficients from first principles during
this session failed due to unit-conversion ambiguities between atomic and natural
relativistic units. The known recursion (Blanchard 1974, Shabaev 2002) uses natural
units (hbar = m = c = 1) where the conversion to atomic units is non-trivial.

**Impact on T8:** The spin-orbit coupling H_SO = xi(r)*L.S requires <r^{-3}> for
l >= 1 states. For the n_r = 0 states (2p3/2, 3d5/2, ...), the formula is available.
For n_r >= 1 states (2p1/2, 3p1/2, 3p3/2 at n_r=1, ...), T8 can use the
non-relativistic <1/r^3> from T1 as the leading term, with the understanding
that relativistic corrections to this are O(alpha^2).

## 3. Transcendental content

All outputs are algebraic over Q(Z, alpha, gamma_kappa):
- gamma_kappa = sqrt(kappa^2 - Z^2*alpha^2) satisfies gamma^2 + Z^2*alpha^2 = kappa^2.
- No pi, no zeta, no log, no Gamma functions.
- This confirms membership in the R_sp spinor ring from T5.
- The relativistic correction factor gamma is a new "spinor-intrinsic" algebraic
  quantity per Paper 18's taxonomy — it arises from the first-order Dirac operator
  (not the second-order Laplacian).

## 4. What T8 can consume

- `radial_expectation_relativistic(n, kappa, -1, Z, alpha)`: exact <1/r> for all states.
- `radial_expectation_relativistic(n, kappa, -3, Z, alpha)`: exact <1/r^3> for n_r=0
  states with l >= 1 (the dominant spin-orbit correction term).
- `dirac_principal_quantum_number(n, kappa, Z, alpha)`: the effective N_D.
- For n_r >= 1 spin-orbit: use T1's `inverse_r_cubed_hydrogenic(n, l, Z)` as the
  leading-order value. The full relativistic <1/r^3> for these states requires
  the Kramers-Pasternak recursion (deferred to future work).

## 5. Test summary

97/97 tests pass (66 T1 regression + 31 new T7 tests).

## 6. References

- Hellmann-Feynman for <1/r>: direct dE/dZ of the Dirac eigenvalue.
- Pochhammer ratios: from the r^gamma * exp(-lambda*r) single-term structure of
  the n_r = 0 Dirac wavefunction.
- BLP "Quantum Electrodynamics" §36 (Dirac-Coulomb expectation values).
- Rose "Relativistic Electron Theory" (1961) Ch. 7.
- Grant "Relativistic Quantum Theory of Atoms and Molecules" (2007) Ch. 4-5.
- Blanchard (1974) "Coulomb Green's function" (Kramers-Pasternak recursion).
- Martinez-y-Romero et al. (2004) arXiv:physics/0402061 (three-term recursions).
