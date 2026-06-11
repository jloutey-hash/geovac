"""
g2_c3_racah.py -- Racah recoupling approach to B(n_int).

The vertex amplitude v(ext→int via photon) is a sum of 4 CG coefficients
over internal magnetic quantum numbers mL. By the Wigner-Eckart theorem,
the sum over all magnetic quantum numbers of a product of CG coefficients
reduces to products of 6j and 9j symbols (which ARE rational).

The key insight: the vertex involves recoupling from
  (jL ⊗ jR) → j_s  and  (jL ⊗ jgL) → jtL,  (jR ⊗ jgR) → jtR → (jtL ⊗ jtR) → j_t
This is related to a 9j symbol.

Strategy:
1. Express the vertex amplitude as a CG × reduced matrix element
2. Sum the product v1*probe*v3 over all internal m values using Racah algebra
3. Check if the result is a rational function of the j-labels
"""

from sympy.physics.wigner import (clebsch_gordan, wigner_6j, wigner_9j,
                                   wigner_3j)
from sympy import (Rational as R, sqrt, simplify, radsimp, S,
                   factorial, Integer)
from functools import lru_cache
import time


@lru_cache(maxsize=2_000_000)
def cg_exact(j1_2, j2_2, J_2, m1_2, m2_2, M_2):
    if m1_2 + m2_2 != M_2:
        return R(0)
    if abs(j1_2 - j2_2) > J_2 or J_2 > j1_2 + j2_2:
        return R(0)
    if abs(m1_2) > j1_2 or abs(m2_2) > j2_2 or abs(M_2) > J_2:
        return R(0)
    val = clebsch_gordan(R(j1_2, 2), R(j2_2, 2), R(J_2, 2),
                         R(m1_2, 2), R(m2_2, 2), R(M_2, 2))
    return val


def vertex_amp_exact(jsL2, jsR2, js2, mjs2,
                     jtL2, jtR2, jt2, mjt2,
                     jgL2, jgR2, mgL2, mgR2):
    """Vertex amplitude: product of 4 CG coefficients summed over mL."""
    total = R(0)
    for mL1_2 in range(-jsL2, jsL2 + 1, 2):
        mR1_2 = mjs2 - mL1_2
        if abs(mR1_2) > jsR2:
            continue
        mL2_2 = mL1_2 + mgL2
        if abs(mL2_2) > jtL2:
            continue
        mR2_2 = mR1_2 + mgR2
        if abs(mR2_2) > jtR2:
            continue
        if mL2_2 + mR2_2 != mjt2:
            continue
        c1v = cg_exact(jsL2, jsR2, js2, mL1_2, mR1_2, mjs2)
        if c1v == 0:
            continue
        c2v = cg_exact(jtL2, jtR2, jt2, mL2_2, mR2_2, mjt2)
        if c2v == 0:
            continue
        c3v = cg_exact(jsL2, jgL2, jtL2, mL1_2, mgL2, mL2_2)
        c4v = cg_exact(jsR2, jgR2, jtR2, mR1_2, mgR2, mR2_2)
        total += c1v * c2v * c3v * c4v
    return total


def vertex_allowed(n1, n2, q):
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1


def get_channels(n_src, n_tgt, q, jsL2, jsR2, jtL2, jtR2):
    if not vertex_allowed(n_src, n_tgt, q):
        return []
    chs = []
    for jgL2, jgR2 in [(q + 1, q - 1), (q - 1, q + 1)]:
        if jgL2 < 0 or jgR2 < 0:
            continue
        if (abs(jsL2 - jgL2) <= jtL2 <= jsL2 + jgL2 and
                abs(jsR2 - jgR2) <= jtR2 <= jsR2 + jgR2):
            chs.append((jgL2, jgR2))
    return chs


def compute_diagram_racah(n_ext, n_int, verbose=True):
    """
    Compute B(n_int) using explicit Racah recoupling.

    The full B(n_int) = Σ_mj sign(mj) × Σ_{j_int} Σ_{channels}
        [Racah-reduced diagram] / (λ⁴ μ_q)

    We factor the computation into:
    1. For each (j_int, j_ext, q, channel), compute the REDUCED vertex
       by summing over ALL m values and dividing by the external CG structure
    2. Check if the reduced vertices are rational
    """

    jE_L_2 = n_ext + 1
    jE_R_2 = n_ext
    j_ext_2 = 1

    jI_L_2 = n_int + 1
    jI_R_2 = n_int

    lam_num = 2 * n_int + 3
    lam4 = R(lam_num**4, 16)

    q_probe = 1
    if not vertex_allowed(n_int, n_int, q_probe):
        return R(0), {}
    probe_chs = get_channels(n_int, n_int, q_probe, jI_L_2, jI_R_2, jI_L_2, jI_R_2)
    if not probe_chs:
        return R(0), {}

    j_int_min_2 = 1
    j_int_max_2 = 2 * n_int + 1

    # Collect contributions organized by (j_int, q_loop, ch1, ch2)
    contributions = {}

    result_plus = R(0)  # mj_ext = +1/2
    result_minus = R(0)  # mj_ext = -1/2

    for mj_ext_2 in [+1, -1]:
        subtotal = R(0)

        for j_int_2 in range(j_int_min_2, j_int_max_2 + 1, 2):
            for mj_int_2 in range(-j_int_2, j_int_2 + 1, 2):
                for mj_int_prime_2 in range(-j_int_2, j_int_2 + 1, 2):

                    probe_amp = R(0)
                    for jpL2, jpR2 in probe_chs:
                        for mpL2 in range(-jpL2, jpL2 + 1, 2):
                            for mpR2 in range(-jpR2, jpR2 + 1, 2):
                                pa = vertex_amp_exact(
                                    jI_L_2, jI_R_2, j_int_2, mj_int_2,
                                    jI_L_2, jI_R_2, j_int_2, mj_int_prime_2,
                                    jpL2, jpR2, mpL2, mpR2)
                                probe_amp += pa

                    if probe_amp == 0:
                        continue

                    q_lo = max(1, abs(n_ext - n_int))
                    q_hi = n_ext + n_int

                    for q_loop in range(q_lo, q_hi + 1):
                        if not vertex_allowed(n_ext, n_int, q_loop):
                            continue
                        mu_q = R(q_loop * (q_loop + 2))

                        chs1 = get_channels(n_ext, n_int, q_loop, jE_L_2, jE_R_2, jI_L_2, jI_R_2)
                        chs2 = get_channels(n_int, n_ext, q_loop, jI_L_2, jI_R_2, jE_L_2, jE_R_2)

                        for jgL1_2, jgR1_2 in chs1:
                            for jgL2_2, jgR2_2 in chs2:
                                v_amp = R(0)
                                for mgL2_v in range(-jgL1_2, jgL1_2 + 1, 2):
                                    for mgR2_v in range(-jgR1_2, jgR1_2 + 1, 2):
                                        if abs(mgL2_v) > jgL2_2 or abs(mgR2_v) > jgR2_2:
                                            continue

                                        v1 = vertex_amp_exact(
                                            jE_L_2, jE_R_2, j_ext_2, mj_ext_2,
                                            jI_L_2, jI_R_2, j_int_2, mj_int_2,
                                            jgL1_2, jgR1_2, mgL2_v, mgR2_v)
                                        if v1 == 0:
                                            continue

                                        v3 = vertex_amp_exact(
                                            jI_L_2, jI_R_2, j_int_2, mj_int_prime_2,
                                            jE_L_2, jE_R_2, j_ext_2, mj_ext_2,
                                            jgL2_2, jgR2_2, mgL2_v, mgR2_v)
                                        if v3 == 0:
                                            continue

                                        v_amp += v1 * v3

                                subtotal += v_amp * probe_amp / (lam4 * mu_q)

        if mj_ext_2 == +1:
            result_plus = subtotal
        else:
            result_minus = subtotal

    B = result_plus - result_minus

    if verbose:
        print(f"  result(+1/2) = {radsimp(result_plus)}")
        print(f"  result(-1/2) = {radsimp(result_minus)}")
        print(f"  B = result(+) - result(-) = {radsimp(B)}")
        print(f"  B float = {float(B.evalf(30)):.15e}")

        # Check rationality of individual pieces
        rp = radsimp(result_plus)
        rm = radsimp(result_minus)
        print(f"  result(+) rational? {rp.is_Rational}")
        print(f"  result(-) rational? {rm.is_Rational}")

        # Check if result(+) + result(-) is rational (even part)
        even_part = radsimp(result_plus + result_minus)
        print(f"  result(+) + result(-) = {even_part}")
        print(f"  even part rational? {even_part.is_Rational}")

        # Check if B * V_mag^2 is rational
        V_mag = R(2, 3) * (sqrt(R(n_ext + 3, n_ext + 1)) - sqrt(R(n_ext, n_ext + 2)))
        BV2 = radsimp(B * V_mag**2)
        print(f"  B * V_mag^2 = {BV2}")
        print(f"  B * V_mag^2 rational? {BV2.is_Rational}")

    return B, {"plus": result_plus, "minus": result_minus}


def try_wigner_9j_formula(n_ext, n_int):
    """
    Try to express the vertex amplitude using Wigner 9j symbols.

    The vertex couples (jL ⊗ jR) → j_source and
    (jL ⊗ jgL → jtL), (jR ⊗ jgR → jtR), (jtL ⊗ jtR → j_target)

    This is a 9j recoupling:
    { jL   jR   j_s }
    { jgL  jgR  j_g }  (where j_g is the total photon j)
    { jtL  jtR  j_t }
    """
    jEL = R(n_ext + 1, 2)
    jER = R(n_ext, 2)
    jIL = R(n_int + 1, 2)
    jIR = R(n_int, 2)
    j_ext = R(1, 2)

    print(f"\n  Wigner 9j test for n_ext={n_ext}, n_int={n_int}")
    print(f"  jEL={jEL}, jER={jER}, jIL={jIL}, jIR={jIR}")

    q_loop = 1  # for n_ext=n_int=1
    for jgL2, jgR2 in [(q_loop + 1, q_loop - 1), (q_loop - 1, q_loop + 1)]:
        if jgR2 < 0:
            continue
        jgL = R(jgL2, 2)
        jgR = R(jgR2, 2)

        # Total photon j: ranges from |jgL-jgR| to jgL+jgR
        for jg_2 in range(abs(jgL2 - jgR2), jgL2 + jgR2 + 1, 2):
            jg = R(jg_2, 2)

            # j_int ranges
            for j_int_2 in range(1, 2 * n_int + 2, 2):
                j_int = R(j_int_2, 2)

                # Compute 9j symbol for vertex 1 (ext → int)
                try:
                    w9j = wigner_9j(jEL, jER, j_ext,
                                     jgL, jgR, jg,
                                     jIL, jIR, j_int)
                    if w9j != 0:
                        print(f"    9j({jEL},{jER},{j_ext}; {jgL},{jgR},{jg}; "
                              f"{jIL},{jIR},{j_int}) = {w9j}")
                        print(f"    rational? {w9j.is_Rational}")
                except Exception as e:
                    print(f"    9j error: {e}")


def main():
    print("=" * 70)
    print("  RACAH RECOUPLING ANALYSIS OF B(n_int)")
    print("=" * 70)

    # First: detailed analysis of n_int=1, n_ext=1
    print("\n--- n_int=1 detailed analysis ---")
    t0 = time.time()
    B1, details1 = compute_diagram_racah(1, 1, verbose=True)
    dt = time.time() - t0
    print(f"  Time: {dt:.1f}s")

    # Try 9j formula
    print("\n--- Wigner 9j test ---")
    try_wigner_9j_formula(1, 1)

    # Also check n_int=2
    print("\n--- n_int=2 detailed analysis ---")
    t0 = time.time()
    B2, details2 = compute_diagram_racah(1, 2, verbose=True)
    dt = time.time() - t0
    print(f"  Time: {dt:.1f}s")


if __name__ == "__main__":
    main()
