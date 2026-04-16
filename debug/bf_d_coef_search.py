"""BF-D supporting: coefficient search for A_SS, A_SOO combining.

Given the 6 candidate combinations in bf_d_he_2P_drake.py all miss NIST by
>100%, this script does a 2-dimensional search over rational coefficients
(c, d) such that A_SS = α² · c · (M² combination) and A_SOO = α² · d ·
(M¹ combination), to find which simple rational forms produce <20% error.

This is a DIAGNOSTIC search, not a fit. We enumerate small rationals
c ∈ {±1/1, ±1/2, ±1/3, ±1/5, ±1/10, ±1/20, ±3/10, ±3/50, ...} and same for d,
and report which achieve <20% on the span.

If some (c, d) gives exact NIST, that's a "numerological" result and tells
us what the correct Drake 1971 coefficients are.
"""

from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fractions import Fraction
from geovac.breit_integrals import breit_ss_radial

ALPHA_CODATA = 7.2973525693e-3
HA_TO_MHZ = 6.5796839204e9
NIST_MHZ = {"P0-P1": 29616.951, "P1-P2": 2291.178, "P0-P2": 31908.129}
F_SS = {0: -2.0, 1: 1.0, 2: -0.2}
F_SOO = {0: 2.0, 1: 1.0, 2: -1.0}
X_J = {J: J * (J + 1) - 4 for J in (0, 1, 2)}


def eval_span_err(A_SS, A_SOO, zeta):
    E_SO = {J: (zeta/2) * X_J[J] for J in (0, 1, 2)}
    E_SS = {J: A_SS * F_SS[J] for J in (0, 1, 2)}
    E_SOO = {J: A_SOO * F_SOO[J] for J in (0, 1, 2)}
    E = {J: E_SO[J] + E_SS[J] + E_SOO[J] for J in (0, 1, 2)}
    split = {"P0-P1": (E[0]-E[1])*HA_TO_MHZ,
             "P1-P2": (E[1]-E[2])*HA_TO_MHZ,
             "P0-P2": (E[0]-E[2])*HA_TO_MHZ}
    return {k: (split[k]-NIST_MHZ[k])/NIST_MHZ[k]*100 for k in split}


def main():
    Z = 2
    M0d = float(breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 0, Z=Z))
    M1d = float(breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 1, Z=Z))
    M2d = float(breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 2, Z=Z))
    M0e = float(breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 0, Z=Z))
    M1e = float(breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 1, Z=Z))
    M2e = float(breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 2, Z=Z))

    alpha = ALPHA_CODATA
    zeta = alpha**2 * Z * 1**3 / 24.0  # Z_eff=1

    # BR-C fit target:
    A_SS_target = -1.202e-6
    A_SOO_target = +5.333e-6

    print(f"At Z=2, alpha=CODATA, Z_eff=1:")
    print(f"  zeta = {zeta:.4e} Ha")
    print(f"  NIST-fit target: A_SS = {A_SS_target:+.4e}, A_SOO = {A_SOO_target:+.4e}\n")

    # Enumerate rational combinations
    rationals = [Fraction(s*p, q) for p in (1, 2, 3) for q in (1, 2, 3, 5, 6, 10, 15, 30, 50)
                                     for s in (-1, +1)]
    rationals = sorted(set(rationals), key=lambda r: (abs(r.numerator) + abs(r.denominator), r))

    # A_SS from two coefficients: A_SS = α² (c_d M²_d + c_e M²_e)
    # Search over small rationals c_d, c_e such that this matches target within 10%
    print("=== Search for A_SS ≈ α² (c_d M²_dir + c_e M²_exch) to ±10% of target ===")
    ss_candidates = []
    for c_d in rationals:
        for c_e in rationals:
            A_SS = alpha**2 * (float(c_d) * M2d + float(c_e) * M2e)
            if abs(A_SS - A_SS_target) < abs(A_SS_target) * 0.10:
                ss_candidates.append((c_d, c_e, A_SS, (A_SS - A_SS_target)/A_SS_target * 100))
    ss_candidates.sort(key=lambda x: abs(x[3]))
    print(f"  Found {len(ss_candidates)} with <10% match. Top 10:")
    for c_d, c_e, A_SS, err in ss_candidates[:10]:
        print(f"    c_d={c_d}, c_e={c_e}:  A_SS = {A_SS:+.4e} ({err:+.2f}% from target)")

    # Same for A_SOO with k=1 integrals
    print("\n=== Search for A_SOO ≈ α² (d_d M¹_dir + d_e M¹_exch) to ±10% of target ===")
    soo_candidates = []
    for d_d in rationals:
        for d_e in rationals:
            A_SOO = alpha**2 * (float(d_d) * M1d + float(d_e) * M1e)
            if abs(A_SOO - A_SOO_target) < abs(A_SOO_target) * 0.10:
                soo_candidates.append((d_d, d_e, A_SOO, (A_SOO - A_SOO_target)/A_SOO_target * 100))
    soo_candidates.sort(key=lambda x: abs(x[3]))
    print(f"  Found {len(soo_candidates)} with <10% match. Top 10:")
    for d_d, d_e, A_SOO, err in soo_candidates[:10]:
        print(f"    d_d={d_d}, d_e={d_e}:  A_SOO = {A_SOO:+.4e} ({err:+.2f}% from target)")

    # Also try with M^0, M^1, M^2 mixed for A_SOO (since SOO has k=0 and k=1 pieces)
    print("\n=== Search for A_SOO ≈ α² (d_0 M⁰_comb + d_1 M¹_comb) mixing orders ===")
    mixed_candidates = []
    for d_0d in rationals:
        for d_1d in rationals:
            for d_0e in (Fraction(0), ):  # Simplify: k=0 exchange only at d_0e=0
                for d_1e in rationals:
                    A_SOO = alpha**2 * (
                        float(d_0d) * M0d + float(d_0e) * M0e +
                        float(d_1d) * M1d + float(d_1e) * M1e
                    )
                    if abs(A_SOO - A_SOO_target) < abs(A_SOO_target) * 0.05:
                        mixed_candidates.append(
                            (d_0d, d_0e, d_1d, d_1e, A_SOO,
                             (A_SOO - A_SOO_target)/A_SOO_target*100))
    mixed_candidates.sort(key=lambda x: abs(x[5]))
    print(f"  Found {len(mixed_candidates)} with <5% match. Top 10:")
    for d_0d, d_0e, d_1d, d_1e, A_SOO, err in mixed_candidates[:10]:
        print(f"    [d_0d, d_0e, d_1d, d_1e] = [{d_0d}, {d_0e}, {d_1d}, {d_1e}]: "
              f"A_SOO = {A_SOO:+.4e} ({err:+.2f}% from target)")

    # Now test combining the best A_SS and A_SOO
    if ss_candidates and soo_candidates:
        best_ss = ss_candidates[0]
        best_soo = soo_candidates[0]
        c_d, c_e, A_SS, _ = best_ss
        d_d, d_e, A_SOO, _ = best_soo
        errs = eval_span_err(A_SS, A_SOO, zeta)
        max_err = max(abs(e) for e in errs.values())
        print(f"\n=== Combined best A_SS × best A_SOO ===")
        print(f"  A_SS = α² ({c_d} M²_d + {c_e} M²_e) = {A_SS:+.4e}")
        print(f"  A_SOO = α² ({d_d} M¹_d + {d_e} M¹_e) = {A_SOO:+.4e}")
        print(f"  Splitting errors: {errs}")
        print(f"  max |rel err| = {max_err:.2f}%")

    # Brute force the combined 4-coefficient search for <20% span error
    print("\n=== Brute-force 4-coef search: coefficients where |span err| < 20% ===")
    found_20pct = []
    rationals_small = [Fraction(s*p, q) for p in (1,) for q in (1, 2, 3, 5, 10)
                                              for s in (-1, +1)] + [Fraction(0)]
    for c_d in rationals_small:
        for c_e in rationals_small:
            for d_d in rationals_small:
                for d_e in rationals_small:
                    A_SS = alpha**2 * (float(c_d) * M2d + float(c_e) * M2e)
                    A_SOO = alpha**2 * (float(d_d) * M1d + float(d_e) * M1e)
                    errs = eval_span_err(A_SS, A_SOO, zeta)
                    max_err = max(abs(e) for e in errs.values())
                    if max_err < 20:
                        found_20pct.append((c_d, c_e, d_d, d_e, errs, max_err))
    found_20pct.sort(key=lambda x: x[5])
    print(f"  Found {len(found_20pct)} 4-tuples with max|err| < 20%. Top 10:")
    for c_d, c_e, d_d, d_e, errs, max_err in found_20pct[:10]:
        print(f"    (c_d={c_d}, c_e={c_e}, d_d={d_d}, d_e={d_e}):  max_err = {max_err:.2f}%")


if __name__ == "__main__":
    main()
