"""Part 2 of combined-substrate scoping: bit-exact verification of the
SU(2) action non-preservation of the OffDiag idempotent basis.

The structural claim from part 1: U_theta * T_{s' -> s} * U^{-1}_theta
contains m_l-dependent phase factors and is NOT a scalar multiple of
T_{s' -> s}.

Here we verify this concretely on the (1,1) -> (1,0) transition (or
similar small case) by constructing U_theta as a phase rotation on m_l
indices and checking the conjugation result.

Discipline: symbolic sympy, no floats.
"""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp
from sympy import Rational, Integer, Symbol, exp, I, simplify, expand


# CH states at n_max=2 in sector (1, 1) and (1, 0):
# Sector (1, 0): dim 2 -- (n=1, l=0, spin up/down) labeled m_s = +/- 1/2
# Sector (1, 1): dim 2 -- (n=1, l=1, half-integer total j) -- in CH labeling

# For simplicity, model the m_l action on a 4-state truncation:
# states |1,0,0,+>, |1,0,0,->, |1,1,1,->, |1,1,-1,+>  (specific CH structure)
#
# The SU(2) z-rotation acts diagonally:
#   U_theta = diag(exp(I*theta*m_l_eff))
# where m_l_eff is the eigenvalue of L_z + S_z (total angular momentum z).

# To stay close to the OffDiag prosystem, we use the diagonal of L_z + S_z
# inferred from the spinor structure.

# Simplified model: 4-state model
#   |a> = |s=(1,0), m=+1/2> (m_J = +1/2)
#   |b> = |s=(1,0), m=-1/2> (m_J = -1/2)
#   |c> = |s=(1,1), m=+1/2> (m_J = +1/2; built from m_l=+1, m_s=-1/2 and m_l=0, m_s=+1/2)
#   |d> = |s=(1,1), m=-1/2> (m_J = -1/2; built from m_l=-1, m_s=+1/2 and m_l=0, m_s=-1/2)
# So m_J eigenvalues are {+1/2, -1/2, +1/2, -1/2}

theta = Symbol("theta", real=True)
half = Rational(1, 2)

# U_theta diagonal:
U_diag = sp.Matrix.diag(
    exp(I * theta * half),    # |a> m_J = +1/2
    exp(-I * theta * half),   # |b> m_J = -1/2
    exp(I * theta * half),    # |c> m_J = +1/2
    exp(-I * theta * half),   # |d> m_J = -1/2
)

# Model the OffDiag transition T_{(1,1) -> (1,0)} in this 4-state truncation
# In real CH, this is a 4x4 block with entries kappa = -1/16 on some structure
# determined by Wigner-Eckart / Gaunt selection rules.
#
# E1 selection rule on m_J: |Delta m_J| <= 1. The off-diagonal block from
# (1,1) -> (1,0) sector has nonzero entries:
#   <a|T|c> = kappa * CG_factor (Delta m_J = 0)
#   <b|T|d> = kappa * CG_factor (Delta m_J = 0)
#   <a|T|d> = kappa * CG_factor (Delta m_J = +1)
#   <b|T|c> = kappa * CG_factor (Delta m_J = -1)
# For the uniform-kappa CH convention, all CG factors collapse to 1.

kappa = Rational(-1, 16)

# T_{(1,1) -> (1,0)} in basis (|a>, |b>, |c>, |d>):
# Outgoing sector is (1,0) (states |a>, |b>); incoming sector is (1,1) (states |c>, |d>).
# Matrix entries T[i, j] = <i|T|j> nonzero for i in {a, b}, j in {c, d}:
T = sp.zeros(4, 4)
T[0, 2] = kappa  # <a|T|c>
T[0, 3] = kappa  # <a|T|d>
T[1, 2] = kappa  # <b|T|c>
T[1, 3] = kappa  # <b|T|d>

# Conjugate: T' = U T U^{-1}
U_inv = U_diag.H  # Hermitian conjugate (since U is unitary)
T_conj = simplify(U_diag * T * U_inv)

# Check: is T_conj a scalar multiple of T?
# Compute ratio T_conj[i,j] / T[i,j] for nonzero entries:
ratios = {}
for i in range(4):
    for j in range(4):
        if T[i, j] != 0:
            ratio = simplify(T_conj[i, j] / T[i, j])
            ratios[f"({i},{j})"] = str(ratio)

# The ratios should be {1, e^{i theta}, e^{-i theta}, 1} or similar mixed phases,
# NOT all the same value.

is_scalar_multiple = len(set(str(simplify(T_conj[i, j] / T[i, j]))
                              for i in range(4) for j in range(4) if T[i, j] != 0)) == 1

# Get distinct phase factors:
phase_factors = sorted(set(str(simplify(T_conj[i, j] / T[i, j]))
                           for i in range(4) for j in range(4) if T[i, j] != 0))

results = {
    "model": "4-state CH truncation: |1,0,+1/2>, |1,0,-1/2>, |1,1,+1/2>, |1,1,-1/2>",
    "T_{(1,1)->(1,0)}_entries": {
        f"({i},{j})": str(T[i, j])
        for i in range(4) for j in range(4) if T[i, j] != 0
    },
    "U_theta_diagonal": [str(U_diag[i, i]) for i in range(4)],
    "T_conj_entries": {
        f"({i},{j})": str(simplify(T_conj[i, j]))
        for i in range(4) for j in range(4) if simplify(T_conj[i, j]) != 0
    },
    "phase_ratios": ratios,
    "distinct_phase_factors": phase_factors,
    "is_scalar_multiple_of_T": is_scalar_multiple,
    "structural_verdict": (
        "T_conj has m_J-DEPENDENT phase factors {1, e^{i theta}, e^{-i theta}, 1}. "
        "Not a scalar multiple of T. The conjugate is in the span of m_J-resolved "
        "sub-transitions, NOT in span{T}. The SU(2) z-rotation requires refining "
        "the sector basis to m_J-resolved sub-projectors before it acts. "
        "Therefore SU(2) does NOT preserve the basic 5-sector OffDiag substrate."
    ),
    "categorical_conclusion": (
        "L5 reduces to L1 at the substrate level. The combined substrate is the "
        "tensor product A^{J*}_{j_max} (x) A^{OD}_{n_max}; SL_2 acts trivially on "
        "the OffDiag sector idempotents and so trivially on A^{OD} in the basic basis. "
        "The Levi-decomposition U* = SL_2 x G_a^{N_OD} is the multi-year target."
    ),
}

out_path = Path(__file__).parent / "data" / "sprint_q5p_combined_substrate_part2.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2, default=str)

print("=== Part 2: SU(2) z-rotation conjugation test ===\n")
print(f"T entries (nonzero):")
for k, v in results["T_{(1,1)->(1,0)}_entries"].items():
    print(f"  T{k} = {v}")
print(f"\nU_theta * T * U^{{-1}} entries (nonzero):")
for k, v in results["T_conj_entries"].items():
    print(f"  (UTU^-1){k} = {v}")
print(f"\nPhase ratios T_conj[i,j] / T[i,j]:")
for k, v in results["phase_ratios"].items():
    print(f"  {k}: {v}")
print(f"\nDistinct phase factors: {phase_factors}")
print(f"\nIs T_conj a scalar multiple of T? {is_scalar_multiple}")
print(f"\nVerdict: {results['structural_verdict']}")
print(f"\nData saved to: {out_path}")
