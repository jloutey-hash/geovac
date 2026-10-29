# Pre-registered check: the one-electron bond-angle operator derived bottom-up (BU)

**Registered:** 2026-10-28, before any derivation. PI-directed ("Maybe we should put a sub agent on this? Bottom-up
derivation ... build the order-8 equation from the label rules. That would confirm the guessed equation
independently.").
**Commit status:** committed on its own, before the derivation is run.
**Rests on:** NA (`2026-10-28_nuclear_attraction_angle_operator.md`): the generic sequence of case N1
(D1 = 1.8, D2 = 2.6, zX = 1.3, zH = 1.0) is annihilated by a unique operator of order 8 and degree 31 in the ansatz
sum a_ij x^j T^i, found by exact guessing mod two primes.

## The check

Derive an operator annihilating V(x) for N1 **without guessing**: from the three-term rules in L of i_L, k_L and
P_L, the exact radial integrations, and closure/elimination algorithms (creative telescoping, Zeilberger-type
summation, or elimination in the radial Poisson structure). No null vector may be fitted to coefficient data; data
may be used only to check. The guessed operator enters only at the final comparison.

**Outcomes.**
- **GO:** the derived operator, reduced mod both primes, equals the guessed (8, 31) operator up to a scalar, or is
  a left multiple of it (then the guessed one right-divides it exactly, checked mod both primes).
- **BORDERLINE:** a derived operator annihilates the generic sequence (checked on all L <= 1000, both primes) but
  the relation to the guessed one cannot be settled (for example, it is too large to divide).
- **STOP:** no operator derived within the agent's budget, or a derived operator that does not annihilate the
  sequence (that would point at an error on one side; both are then re-examined).

Probability of GO: 0.6.

## Novelty

The methods are standard (Zeilberger, Chyzak, Koutschan); a bottom-up derivation of a bond-angle operator for a
three-centre integral was not found (AO novelty check).
