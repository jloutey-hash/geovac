"""
Track alpha-C: Avery-Aquilanti Sturmian Connection to kappa = -1/16 and K.

Three subtasks:
  (1) Compute Coulomb-Sturmian / hydrogenic normalization ratios for (n,l)
      up to n=3 in exact rational arithmetic and check whether kappa = -1/16
      appears.
  (2) Compute the Fock-weight matrix element <1|w|1>_R3 / <1|1>_R3 at
      p0 = Z/n = 1, with w(p) = 1/(p^2 + p0^2)^2 (position-space kernel
      equivalent), and check whether kappa = -1/16 appears.
  (3) Hopf restriction: fiber-averaged vs base-only Sturmian normalizations,
      ratios that might equal 2*pi / pi / B = 42.

All computations use sympy with exact rational (and, where needed, symbolic pi)
arithmetic.
"""

from __future__ import annotations

import json
import os
from fractions import Fraction

import sympy as sp

# --------------------------------------------------------------------------
# Symbolic objects
# --------------------------------------------------------------------------
r, p, p0, k, Z, chi, alpha_var = sp.symbols(
    "r p p0 k Z chi alpha", positive=True, real=True
)
n_sym, l_sym = sp.symbols("n l", integer=True, positive=True)

pi = sp.pi


def hydrogenic_radial(n: int, l: int, Zval: sp.Expr):
    """
    Standard hydrogenic radial wavefunction R_{nl}(r) with
    int_0^inf |R|^2 r^2 dr = 1.
    Returns R_{nl}(r) as a sympy expression in r.
    """
    # rho = 2 Z r / n
    rho = 2 * Zval * r / n
    # Associated Laguerre L_{n-l-1}^{2l+1}(rho)
    Lag = sp.assoc_laguerre(n - l - 1, 2 * l + 1, rho)
    # Normalization
    N = sp.sqrt(
        (2 * Zval / n) ** 3
        * sp.factorial(n - l - 1)
        / (2 * n * sp.factorial(n + l))
    )
    R = N * sp.exp(-rho / 2) * rho ** l * Lag
    return sp.simplify(R)


def radial_norm_sq(expr):
    """int_0^inf |expr|^2 r^2 dr, symbolic."""
    return sp.integrate(expr ** 2 * r ** 2, (r, 0, sp.oo))


# --------------------------------------------------------------------------
# SUBTASK 1: Coulomb Sturmian <-> hydrogenic normalization ratios
# --------------------------------------------------------------------------
#
# The Coulomb Sturmian basis S_{nlm}(r;k) is defined by a common exponent k
# (independent of n). The standard Avery/Aquilanti convention is
#
#    S_{nlm}(r;k) = N_nl^S (2 k r)^l e^{-kr} L_{n-l-1}^{2l+1}(2 k r) Y_{lm}
#
# normalized against the Sturmian weight k/r:
#
#    int S_{n'l'm'}*(r;k) (k/r) S_{nlm}(r;k) d^3 r = delta_{nn'} delta_{ll'} delta_{mm'}
#
# This is the potential-weighted inner product.  Equivalently, without the
# weight,
#
#    int |S_{nlm}|^2 d^3 r = n / k  (the "Sturmian overlap" relation).
#
# The corresponding hydrogenic orbital (energy Z^2/(2 n^2)) uses
# exponent k_n = Z/n (different for each n).  Writing the hydrogenic orbital
# as
#
#    psi_{nlm}(r) = N_nl^H (2 k_n r)^l e^{-k_n r} L_{n-l-1}^{2l+1}(2 k_n r) Y_{lm}
#
# with int |psi|^2 d^3 r = 1, the ratio
#
#    R_{nl} := (N_nl^H / N_nl^S)  (at k = k_n = Z/n)
#
# measures the relative calibration of the two normalizations at the same
# energy shell.  We compute these ratios exactly using sympy.

def sturmian_radial(n: int, l: int, kval):
    """
    Coulomb Sturmian normalized via ⟨S|S⟩_{R3} = n/k (Avery's convention).

    The radial functional form is identical to hydrogenic with Z -> k*n,
    i.e. S_{nl}(r;k) = hydrogenic_radial(n, l, Z=k*n) up to an overall
    rescaling that enforces ⟨S|S⟩_{R3} = n/k.
    """
    # Start from the unnormalized form
    rho = 2 * kval * r  # Sturmians use 2kr, not 2Zr/n
    Lag = sp.assoc_laguerre(n - l - 1, 2 * l + 1, rho)
    unnorm = sp.exp(-kval * r) * rho ** l * Lag
    # Compute <unnorm|unnorm>
    nsq = sp.integrate(unnorm ** 2 * r ** 2, (r, 0, sp.oo))
    nsq = sp.simplify(nsq)
    # Enforce <S|S> = n/k
    N = sp.sqrt(sp.Rational(n) / kval / nsq)
    return sp.simplify(N * unnorm), nsq, N


def subtask_1():
    results = {}
    cases = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]
    Z_val = sp.Integer(1)
    kappa_geovac = sp.Rational(-1, 16)

    for (n, l) in cases:
        k_val = Z_val / n  # energy shell at the hydrogenic eigenvalue
        # Hydrogenic radial (normalized to 1)
        psi = hydrogenic_radial(n, l, Z_val)
        # Sturmian radial at k = Z/n (normalized to n/k)
        chi_s, chi_nsq_raw, chi_N = sturmian_radial(n, l, k_val)

        # Overlap in r^3
        overlap = sp.integrate(psi * chi_s * r ** 2, (r, 0, sp.oo))
        overlap = sp.simplify(overlap)

        # Sturmian norm (should be n/k = n^2 at k=1/n)
        chi_norm_sq = sp.integrate(chi_s ** 2 * r ** 2, (r, 0, sp.oo))
        chi_norm_sq = sp.simplify(chi_norm_sq)

        # Ratio R_nl = psi/chi_s at the origin-regular leading behavior.
        # Equivalently, since both have the same functional form up to
        # an overall constant, the ratio is just
        #     R_nl = (hydrogenic prefactor) / (Sturmian prefactor).
        # Compute at a fiducial point, say r = 1/Z (Bohr-like scale).
        r_test = sp.Rational(1)
        R_ratio = sp.simplify(psi.subs(r, r_test) / chi_s.subs(r, r_test))

        results[f"n={n},l={l}"] = {
            "k": str(k_val),
            "hydrogenic_norm_sq_check": str(
                sp.simplify(sp.integrate(psi ** 2 * r ** 2, (r, 0, sp.oo)))
            ),
            "sturmian_norm_sq_nk": str(chi_norm_sq),
            "expected_nk": str(sp.Rational(n) / k_val),
            "overlap_psi_sturmian": str(overlap),
            "ratio_R_nl (psi/sturmian at r=1)": str(R_ratio),
        }

    # --- Check whether kappa = -1/16 appears ---
    # Well-known identity:
    # int S_{nlm}^* (-1/2 nabla^2 - k^2/2) S_{nlm} d^3 r = 0
    # Hydrogenic kinetic + potential = -k^2/2 (total energy).
    # In GeoVac units the graph Laplacian eigenvalues l_n = -(n^2 - 1)
    # map to E_n = (1/16) * l_n /... NO: the mapping is
    #     E_n = -1/(2 n^2),   lambda_n = -(n^2 - 1)
    # so E_n = -1/(2 n^2) and lambda_n/n^2 = -(1 - 1/n^2). The kappa
    # factor arises from the distinction between the graph Laplacian
    # (D - A) and the continuous Laplace-Beltrami (Paper 7 notes
    # kappa = -1/16 comes from a 16 = 4 * (d_max^2 = 16) / 2 factor).
    #
    # Check: do any Sturmian normalization ratios equal +- 1/16, +- 1/4,
    # +- 1/2, or are related to them by simple powers of p0?
    targets = {
        "kappa = -1/16": sp.Rational(-1, 16),
        "+1/16": sp.Rational(1, 16),
        "1/8": sp.Rational(1, 8),
        "1/4": sp.Rational(1, 4),
        "1/2": sp.Rational(1, 2),
    }
    matches = []
    for case_label, info in results.items():
        for name, tval in targets.items():
            r_expr = sp.sympify(info["ratio_R_nl (psi/sturmian at r=1)"])
            if sp.simplify(r_expr - tval) == 0:
                matches.append((case_label, "R_nl", name))
            nsq_expr = sp.sympify(info["sturmian_norm_sq_nk"])
            if sp.simplify(nsq_expr - tval) == 0:
                matches.append((case_label, "sturmian_norm_sq", name))

    results["_matches_kappa_targets"] = matches
    results["_kappa_geovac"] = str(kappa_geovac)
    return results


# --------------------------------------------------------------------------
# SUBTASK 2: Fock weight matrix element <1|w|1>_R3
# --------------------------------------------------------------------------
#
# In momentum space, the Fock projection introduces a weight
#   w(p) = [(p0^2 + p^2) / (2 p0)]^2
# so that the momentum-space hydrogenic wavefunction phi(p) is related to a
# hyperspherical harmonic Y on S^3 by
#   phi(p) = (w(p))^{-1} Y(Omega(p))
# and the inner product becomes
#   int |Y|^2 dOmega_{S^3} = int |phi|^2 w(p)^{-2} dOmega_{S^3}.
#
# In position space the analogue is the "Coulomb weight" ⟨psi|(p0/r)|psi⟩:
# this is the Sturmian potential-weight inner product at k = p0.  For the
# n=1 ground state,
#   ⟨1s|(1/r)|1s⟩_{R3} = Z.
# At p0 = k = Z/n = Z, this gives
#   ratio = ⟨1s|(p0/r)|1s⟩ / ⟨1s|1s⟩ = Z * Z / 1 = Z^2.
#
# More relevant: the S^3 normalization introduces a factor of
#   (2/p0) (Jacobian of the Fock map at the south pole times p0^3 volume).
# At p0 = 1 this is just a factor of 2.

def subtask_2():
    results = {}
    Z_val = sp.Integer(1)
    # n=1 hydrogenic
    psi1 = hydrogenic_radial(1, 0, Z_val)
    nsq = sp.simplify(sp.integrate(psi1 ** 2 * r ** 2, (r, 0, sp.oo)))
    # <1|1/r|1>
    inv_r = sp.simplify(sp.integrate(psi1 ** 2 * r, (r, 0, sp.oo)))
    # <1|r|1>
    exp_r = sp.simplify(sp.integrate(psi1 ** 2 * r ** 3, (r, 0, sp.oo)))

    # Momentum-space Fock weight.  The 1s momentum wavefunction is
    #   phi_1s(p) = (1/pi) * (2 p0)^{5/2} / (p^2 + p0^2)^2
    # (see, e.g., Bethe-Salpeter Eq. 8.8).  The normalization check is
    #   int |phi_1s(p)|^2 d^3 p = 1.
    # At p0 = Z = 1:
    #
    #   <phi|(p^2 + p0^2)^(-2)|phi> / <phi|phi>
    #
    # is the Fock-weighted average.

    p_var, p0_var = sp.symbols("p p0", positive=True)
    phi_1s_p = (1 / sp.pi) * (2 * p0_var) ** sp.Rational(5, 2) / (
        p_var ** 2 + p0_var ** 2
    ) ** 2
    # radial integral in 3D: 4*pi * int phi^2 p^2 dp
    norm_phi = sp.simplify(
        4 * sp.pi * sp.integrate(phi_1s_p ** 2 * p_var ** 2, (p_var, 0, sp.oo))
    )

    # Weight w(p) = (p^2 + p0^2)^(-2)  (what remains of the Fock factor)
    w = 1 / (p_var ** 2 + p0_var ** 2) ** 2
    weighted = sp.simplify(
        4 * sp.pi * sp.integrate(phi_1s_p ** 2 * w * p_var ** 2, (p_var, 0, sp.oo))
    )
    ratio = sp.simplify(weighted / norm_phi)

    # Specialize at p0 = 1
    ratio_p0_1 = sp.simplify(ratio.subs(p0_var, 1))

    # And also the "conformal factor squared at the south pole",
    # Omega(p=0)^2 = (2/p0)^2 = 4/p0^2.
    conformal_sq = sp.Rational(4)  # at p0 = 1

    results = {
        "norm_phi_1s": str(norm_phi),
        "inv_r_ground": str(inv_r),
        "exp_r_ground": str(exp_r),
        "Fock_weighted_ratio (p0 free)": str(ratio),
        "Fock_weighted_ratio at p0=1": str(ratio_p0_1),
        "conformal_factor_squared at p0=1": str(conformal_sq),
        "kappa_target": "-1/16",
        "any_match_to_kappa": (
            sp.simplify(ratio_p0_1 - sp.Rational(-1, 16)) == 0
            or sp.simplify(ratio_p0_1 - sp.Rational(1, 16)) == 0
        ),
    }
    return results


# --------------------------------------------------------------------------
# SUBTASK 3: Hopf restriction - fiber vs base normalizations
# --------------------------------------------------------------------------
#
# The Hopf bundle S^1 -> S^3 -> S^2.  A hyperspherical harmonic Y_{nlm}
# on S^3 decomposes into S^2 spherical harmonics Y_{lm} tensored with
# an S^1 phase e^{i m phi_fiber}.  Averaging over the S^1 fiber gives
# (2*pi) times the S^2 base component with m = 0 only.
#
# The paper 2 alpha formula uses
#   B = sum_{n=1}^{3} sum_{l=0}^{n-1} (2 l + 1) l (l + 1)
# which is 42 for n=1..3.  Let's verify this and then construct the
# analogous Sturmian-weighted sums.

def subtask_3():
    results = {}

    # 1. Verify B = 42
    B = sum((2 * l + 1) * l * (l + 1) for n in range(1, 4) for l in range(n))
    results["B (sum of (2l+1) l(l+1), n=1..3)"] = B

    # Check some other "natural" sums:
    sums = {
        "sum (2l+1)": sum(
            (2 * l + 1) for n in range(1, 4) for l in range(n)
        ),
        "sum (2l+1) l": sum(
            (2 * l + 1) * l for n in range(1, 4) for l in range(n)
        ),
        "sum (2l+1) l^2": sum(
            (2 * l + 1) * l ** 2 for n in range(1, 4) for l in range(n)
        ),
        "sum (2l+1) l(l+1)": sum(
            (2 * l + 1) * l * (l + 1) for n in range(1, 4) for l in range(n)
        ),
        "sum (2l+1) n^2": sum(
            (2 * l + 1) * n ** 2 for n in range(1, 4) for l in range(n)
        ),
        "sum (2l+1) (n^2 - 1)": sum(
            (2 * l + 1) * (n ** 2 - 1) for n in range(1, 4) for l in range(n)
        ),
        "sum n^2 (from n=1..3)": sum(n ** 2 for n in range(1, 4)),
    }
    results["angular_sums_n1to3"] = sums

    # 2. Compute Sturmian normalization factors and try weighted sums
    #    Using <S_nl|S_nl> = n/k = n^2 at k = 1/n.
    sturmian_norm_sq = {}
    for n in range(1, 4):
        for l in range(n):
            # <S|S> = n/k = n^2 (at k = 1/n, Z=1)
            sturmian_norm_sq[(n, l)] = n ** 2
    results["sturmian_norm_sq (n^2 at k=1/n)"] = {
        f"n={n},l={l}": val for (n, l), val in sturmian_norm_sq.items()
    }

    # 3. Fiber vs base
    #
    # The S^3 volume is 2 pi^2; the S^2 volume is 4 pi; the S^1 volume is
    # 2 pi.  Ratio S^3 / S^2 = pi/2.  Ratio S^1 = 2 pi.
    results["S3_volume"] = "2*pi^2"
    results["S2_volume"] = "4*pi"
    results["S1_volume"] = "2*pi"
    results["ratio_S3_to_S2"] = "pi/2"
    results["ratio_S1"] = "2*pi"

    # Compute "fiber-averaged" sums: average over m, i.e. divide by (2l+1)
    # fiber_avg = sum_{n,l} sum_m (l(l+1)) / (2 l+1) = sum (2l+1)(l(l+1))/(2l+1)
    #           = sum l(l+1) = 0 + 2 + 0 + 2 + 6 = 10 for n=1..3
    fiber_avg = sum(l * (l + 1) for n in range(1, 4) for l in range(n))
    results["fiber_avg_l(l+1)"] = fiber_avg
    # Base-only: m=0 only, weight 1
    base_only = sum(l * (l + 1) for n in range(1, 4) for l in range(n))
    results["base_only_l(l+1)"] = base_only

    # Ratio B / fiber_avg
    results["B_over_fiber_avg"] = B / fiber_avg if fiber_avg else None
    results["B_equals_pi_times_something"] = False  # 42 is exact integer

    # 4. Hopf-weighted Sturmian sum
    # The Sturmian norm n^2 at k=1/n, times (2l+1), summed:
    sturm_weighted = sum(
        (2 * l + 1) * n ** 2 for n in range(1, 4) for l in range(n)
    )
    results["sturm_weighted_n^2_(2l+1)"] = sturm_weighted
    # weighted by l(l+1):
    sturm_weighted_B = sum(
        (2 * l + 1) * l * (l + 1) * n ** 2
        for n in range(1, 4) for l in range(n)
    )
    results["sturm_weighted_n^2_(2l+1)l(l+1)"] = sturm_weighted_B
    # Hopf B = 42 (does NOT get multiplied by n^2 in Paper 2). So the
    # Sturmian normalization does NOT yield 42 directly; the Paper 2
    # formula uses the *unweighted* Casimir trace of the base.

    # 5. Check if any ratio gives K / pi = 42 + pi^2/6 - 1/40 = 43.619934
    K_over_pi_numeric = 42 + float(sp.pi) ** 2 / 6 - 1 / 40
    results["K/pi_numeric"] = float(K_over_pi_numeric)

    return results


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    out = {}
    print("=" * 60)
    print("SUBTASK 1: Sturmian <-> hydrogenic normalization")
    print("=" * 60)
    s1 = subtask_1()
    out["subtask_1"] = s1
    for k, v in s1.items():
        print(f"{k}: {v}")

    print()
    print("=" * 60)
    print("SUBTASK 2: Fock weight <1|w|1>")
    print("=" * 60)
    s2 = subtask_2()
    out["subtask_2"] = s2
    for k, v in s2.items():
        print(f"{k}: {v}")

    print()
    print("=" * 60)
    print("SUBTASK 3: Hopf fiber/base restriction")
    print("=" * 60)
    s3 = subtask_3()
    out["subtask_3"] = s3
    for k, v in s3.items():
        print(f"{k}: {v}")

    data_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data",
        "track_alpha_phase4b",
    )
    os.makedirs(data_dir, exist_ok=True)
    json_path = os.path.join(data_dir, "track_c_sturmian.json")
    with open(json_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {json_path}")

    return out


if __name__ == "__main__":
    main()
