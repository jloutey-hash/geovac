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

## Disposition (appended 2026-10-28, after the derivation; nothing above edited)

**Outcome: GO** (the left-multiple clause). Derivation by a sub-agent (memo `debug/bu_bottom_up_memo.md`; scripts
`debug/bu_{ibp,system,derive,data,check,compare,compare_ctrl}.py`); one restart after an authentication error, same
context. Route: integration-by-parts identities for each radial piece's Bessel family (sympy, L symbolic; one master per
piece), giving a 14-dimensional first-order system in L with coefficients in Q(L), U_{L+1} = A(L) U_L,
F_L/(2 pi N) = c(L) . U_L (checked against quadrature to 6e-48); a recurrence in L by the cyclic-vector method (order 9,
degree 40); and the x-operator by requiring sum a_ij x^j T^i to annihilate c . U for every solution of the system, at
random points of F_p and at the low-L boundary rows. Result: order 8 none up to degree 256; order 9 at degree 30
(nullity 1), then (10, 22), (11, 18), (12, 17): NA's trade table shifted up by one order, on both primes.

**Comparison** (`debug/data/bu_compare.json`): on both primes NA's guessed (8, 31) operator G right-divides the derived
operator D exactly, D = M G with M of order 1 (remainder 0). Negative controls fail the division (G with one
coefficient changed; a random operator of G's shape; D's order-8 truncation).

**PM verification.** Independence: `debug/bu_derive.py` reads only the derived system (`debug/data/bu_system_v2.json`);
its null-space search acts on operator images through the system at random L, not on coefficient data. Annihilation,
with NA's own code (`debug/na_modguess.py`) on the generic sequence rebuilt from `debug/na_exact.py`, prime 2^61 - 1:
D leaves 0 nonzero rows out of 962; D with one coefficient changed leaves 962.

**Reading.** The bottom-up route reaches the same operator up to one extra first-order left factor: the system's
solution space holds the generic sequence's plus one more solution, which only initial values (data) would remove.
The guessed order-8 operator of NA is confirmed by an independent route.
