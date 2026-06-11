# Track DP: PSLQ Identification of Drake 1971 Combining Coefficients

**Sprint:** Sprint 6 Track DP.
**Date:** 2026-04-15.
**Status:** **HONEST NEGATIVE on PSLQ; POSITIVE on independent confirmation.**

## 1. Sprint question

Can the Drake combining coefficients (3/50, -2/5, 3/2, -1) be identified
via PSLQ integer relation algorithm on high-precision Breit-Pauli matrix
elements?

## 2. Result summary

| Test | Status | Finding |
|:-----|:------:|:--------|
| Z^3 scaling of M^k integrals | **VERIFIED** | Machine precision across Z=1,2,3,4,10 |
| NIST extraction A_SS | **CONFIRMED** | 0.09% match to Drake (3/50, -2/5) |
| NIST extraction A_SOO | **CONFIRMED** | 0.07% match to Drake (3/2, -1) |
| PSLQ (minimal M^2 basis) | **FAILED** | NIST ~6 sig figs insufficient for PSLQ |
| PSLQ (wide 6-integral basis) | **INCONCLUSIVE** | Found relation with leading zero; artifact |
| C^(2)_0(hat_12)/r_12^3 exchange angular | **ZERO** | Exchange integral vanishes for (1s)(2p) |
| Angular kernel identification | **STRUCTURAL** | Direct integral is NOT a single M^k kernel |

## 3. Why PSLQ failed

PSLQ requires the input values to have more precision than the complexity
of the target relation. For rationals with denominator <= 50, PSLQ needs
at least ~4 significant digits (log10(50^2) ~ 3.4). But the NIST splittings
have only ~6 significant figures, and the linear system solving to extract
A_SS and A_SOO amplifies errors. The PSLQ target has only ~4-5 reliable
digits, which is borderline for identifying denominator-50 rationals.

The PSLQ approach would work if we had:
(a) NIST data to 10+ significant figures, OR
(b) A theoretical computation of the Breit-Pauli matrix element at high
    precision (not going through NIST)

Option (b) requires the correct bipolar expansion coefficients for
C^(K)(hat_12)/r_12^{K+1} (the Sprint 5 DV obstruction). The PSLQ
approach was intended to BYPASS this, but it cannot because NIST is the
only independent source of the matrix element value.

## 4. Key structural finding

The angular integral of C^(2)_0(hat_12)/r_12^3 on the (1s)(2p) exchange
density is **exactly zero**:

```
<Y_00(1) Y_10(2)| C^(2)_0(hat_12)/r_12^3 |Y_10(1) Y_00(2)> = 0
```

verified at 6 test points (r1,r2) to machine precision (~1e-18).

This means the rank-2 tensor part of the SS operator contributes to the
DIRECT Slater integral only. The "M^2_exch" term in Drake's formula
comes from a DIFFERENT origin: it involves the Slater exchange integral
(swapping orbital labels in the ket), not the quantum-mechanical exchange
from wavefunction antisymmetrization.

In other words, Drake's "direct" and "exchange" M^k integrals refer to
the orbital-index structure:
- M^k_dir = R^k(1s,2p; 1s,2p) -- same orbitals on each electron
- M^k_exch = R^k(1s,2p; 2p,1s) -- swapped orbitals on electron 2

NOT to the direct/exchange parts of the antisymmetrized matrix element.

## 5. Independent confirmation of (3/50, -2/5, 3/2, -1)

The NIST extraction gives:
- A_SS(NIST) / A_SS(Drake) = 0.999089 (0.09% match)
- A_SOO(NIST) / A_SOO(Drake) = 1.000732 (0.07% match)

The ~0.1% residual is consistent with higher-order corrections
(mass polarization, QED) not included in the Breit-Pauli approximation.

Z^3 scaling verified to machine precision at Z=1,2,3,4,10, confirming
that the angular combining coefficients are Z-independent.

## 6. Files

| File | Purpose |
|:-----|:--------|
| `debug/dp_drake_pslq.py` | Main computation script |
| `debug/data/dp_drake_pslq.json` | Numerical results |
| `debug/dp_drake_pslq_memo.md` | This memo |

## 7. Conclusion

The Drake coefficients (3/50, -2/5, 3/2, -1) are **independently confirmed**
to 0.09% via NIST extraction + exact radial integrals. PSLQ identification
failed due to insufficient NIST precision. The structural derivation remains
open (Sprint 5 DV obstruction: bipolar expansion coefficients for
C^(K)(hat_12)/r_12^{K+1}).

**New structural insight:** the rank-2 tensor operator C^(2)(hat_12)/r_12^3
has zero exchange angular integral for (1s)(2p), meaning Drake's M^k_exch
comes from the Slater integral orbital ordering, not from wavefunction
antisymmetrization.

## 8. Paper updates

None needed. The existing docstring in `breit_integrals.py` correctly states
that the coefficients were "identified by rational search (Track BF-D)" and
subsequent Sprints (DD, DV) characterized the structure. This sprint adds the
NIST extraction confirmation at 0.09% but does not change any claims.

## 9. Tests added

One new test in `tests/test_breit_integrals.py`:
- `test_drake_nist_extraction_confirms_coefficients`: verifies that NIST-extracted
  A_SS and A_SOO match Drake predictions to <0.2% via independent linear-system
  extraction from the three NIST splittings.
