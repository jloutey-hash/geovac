"""
Extend the g-2 vertex correction from n_int=15 to n_int=25.
Goal: improve precision on c2 from ~2 digits to ~4 digits.

Loads all existing data (n_int=0..12 from ratio_investigation,
n_int=13..15 from g2_extended_nint) and continues.

Expected runtime: several hours (n_int=16 ~400s, scaling ~n^2.5).
"""

from sympy.physics.wigner import clebsch_gordan
from sympy import Rational, S, sqrt as ssqrt
import json, time, sys
import numpy as np

sys.path.insert(0, '.')


def half_ints(j):
    result = []
    m = -j
    while m <= j + Rational(1, 100):
        result.append(m)
        m += 1
    return result


def vertex_amp_pol(j_sL, j_sR, j_s, mj_s,
                   j_tL, j_tR, j_t, mj_t,
                   jgL, jgR, mgL, mgR):
    total = S.Zero
    for mL1 in half_ints(j_sL):
        mR1 = mj_s - mL1
        if abs(mR1) > j_sR:
            continue
        mL2 = mL1 + mgL
        if abs(mL2) > j_tL:
            continue
        mR2 = mR1 + mgR
        if abs(mR2) > j_tR:
            continue
        if mL2 + mR2 != mj_t:
            continue
        c1 = clebsch_gordan(j_sL, j_sR, j_s, mL1, mR1, mj_s)
        if c1 == 0:
            continue
        c2 = clebsch_gordan(j_tL, j_tR, j_t, mL2, mR2, mj_t)
        if c2 == 0:
            continue
        c3 = clebsch_gordan(j_sL, jgL, j_tL, mL1, mgL, mL2)
        c4 = clebsch_gordan(j_sR, jgR, j_tR, mR1, mgR, mR2)
        total += c1 * c2 * c3 * c4
    return total


def vertex_allowed(n1, n2, q):
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1


def get_channels(n_src, n_tgt, q, j_sL, j_sR, j_tL, j_tR):
    if not vertex_allowed(n_src, n_tgt, q):
        return []
    chs = []
    for jgL, jgR in [(Rational(q+1,2), Rational(q-1,2)),
                      (Rational(q-1,2), Rational(q+1,2))]:
        if jgL < 0 or jgR < 0:
            continue
        if (abs(j_sL - jgL) <= j_tL <= j_sL + jgL and
            abs(j_sR - jgR) <= j_tR <= j_sR + jgR):
            chs.append((jgL, jgR))
    return chs


def compute_vertex_3pt_single_nint(n_ext, j_ext, mj_ext, n_int, q_probe=1):
    jE_L = Rational(n_ext + 1, 2)
    jE_R = Rational(n_ext, 2)
    jI_L = Rational(n_int + 1, 2)
    jI_R = Rational(n_int, 2)

    lam = Rational(2*n_int + 3, 2)
    lam4 = lam**4

    j_int_min = abs(jI_L - jI_R)
    j_int_max = jI_L + jI_R

    if not vertex_allowed(n_int, n_int, q_probe):
        return S.Zero

    probe_chs = get_channels(n_int, n_int, q_probe,
                              jI_L, jI_R, jI_L, jI_R)
    if not probe_chs:
        return S.Zero

    total = S.Zero

    for j_int in half_ints(j_int_max):
        if j_int < j_int_min:
            continue
        for mj_int in half_ints(j_int):
            for mj_int_prime in half_ints(j_int):
                probe_amp = S.Zero
                for jpL, jpR in probe_chs:
                    for mpL in half_ints(jpL):
                        for mpR in half_ints(jpR):
                            pa = vertex_amp_pol(
                                jI_L, jI_R, j_int, mj_int,
                                jI_L, jI_R, j_int, mj_int_prime,
                                jpL, jpR, mpL, mpR)
                            probe_amp += pa
                if probe_amp == 0:
                    continue

                q_lo = max(1, abs(n_ext - n_int))
                q_hi = n_ext + n_int

                for q_loop in range(q_lo, q_hi + 1):
                    if not vertex_allowed(n_ext, n_int, q_loop):
                        continue
                    mu_q = Rational(q_loop * (q_loop + 2))

                    chs1 = get_channels(n_ext, n_int, q_loop,
                                        jE_L, jE_R, jI_L, jI_R)
                    chs2 = get_channels(n_int, n_ext, q_loop,
                                        jI_L, jI_R, jE_L, jE_R)

                    for jgL1, jgR1 in chs1:
                        for jgL2, jgR2 in chs2:
                            for mgL in half_ints(jgL1):
                                for mgR in half_ints(jgR1):
                                    if abs(mgL) > jgL2 or abs(mgR) > jgR2:
                                        continue

                                    v1 = vertex_amp_pol(
                                        jE_L, jE_R, j_ext, mj_ext,
                                        jI_L, jI_R, j_int, mj_int,
                                        jgL1, jgR1, mgL, mgR)
                                    if v1 == 0:
                                        continue

                                    v3 = vertex_amp_pol(
                                        jI_L, jI_R, j_int, mj_int_prime,
                                        jE_L, jE_R, j_ext, mj_ext,
                                        jgL2, jgR2, mgL, mgR)
                                    if v3 == 0:
                                        continue

                                    total += v1 * probe_amp * v3 / (lam4 * mu_q)

    return total


def main():
    n_ext = 1
    j_ext = Rational(1, 2)
    alpha = 7.2973525693e-3
    schwinger = alpha / (2 * 3.14159265358979)

    # Load all existing data
    with open('debug/data/alpha_g_minus_2_ratio_investigation.json') as f:
        base_data = json.load(f)

    with open('debug/data/g2_extended_nint.json') as f:
        ext_data = json.load(f)

    V_mag = base_data['V_magnetic_float']

    # Collect all existing B levels
    all_B = {}
    for lev in base_data['per_level']:
        all_B[lev['n_int']] = lev['B_level']
    for lev in ext_data['extended_levels']:
        all_B[lev['n_int']] = lev['B_level']

    n_int_done = max(all_B.keys())
    B_cumul = sum(all_B.values())

    print(f"Loaded existing data: n_int=0..{n_int_done}")
    print(f"Cumulative B = {B_cumul:.12e}")
    print(f"F2/Schwinger = {B_cumul / V_mag / schwinger:.12f}")
    print(f"Schwinger = {schwinger:.12e}")
    print(f"V_magnetic = {V_mag:.12e}")
    print()

    # Extend to n_int=25
    n_int_target = 25
    new_levels = []

    for n_int in range(n_int_done + 1, n_int_target + 1):
        print(f"Computing n_int={n_int}...", flush=True)
        t0 = time.time()

        contrib_up = compute_vertex_3pt_single_nint(
            n_ext, j_ext, Rational(1, 2), n_int, q_probe=1)
        contrib_dn = compute_vertex_3pt_single_nint(
            n_ext, j_ext, Rational(-1, 2), n_int, q_probe=1)

        B_level = float(contrib_up - contrib_dn)
        B_cumul += B_level
        F2_S = B_cumul / V_mag / schwinger

        elapsed = time.time() - t0

        all_B[n_int] = B_level

        print(f"  n_int={n_int}: B_level={B_level:.6e}  "
              f"cum F2/S={F2_S:.12f}  ({elapsed:.1f}s)", flush=True)

        new_levels.append({
            'n_int': n_int,
            'B_level': B_level,
            'F2_over_schwinger': F2_S,
            'time': elapsed,
        })

        # Save checkpoint after each level
        _save_checkpoint(all_B, new_levels, V_mag, schwinger)

    # Final analysis
    _final_analysis(all_B, V_mag, schwinger)


def _save_checkpoint(all_B, new_levels, V_mag, schwinger):
    """Save intermediate results after each n_int level."""
    B_cumul = sum(all_B.values())

    # Power-law fits on all data
    ns = sorted(n for n in all_B.keys() if n >= 3)
    even_n = [n for n in ns if n % 2 == 0 and n >= 4]
    even_B = [all_B[n] for n in even_n]
    odd_n = [n for n in ns if n % 2 == 1 and n >= 3]
    odd_B = [all_B[n] for n in odd_n]

    if len(even_n) >= 2 and len(odd_n) >= 2:
        log_ne = np.log(np.array(even_n, dtype=float))
        log_Be = np.log(np.abs(np.array(even_B)))
        Ae = np.column_stack([np.ones_like(log_ne), -log_ne])
        fit_e, _, _, _ = np.linalg.lstsq(Ae, log_Be, rcond=None)
        Ce, pe = np.exp(fit_e[0]), fit_e[1]

        log_no = np.log(np.array(odd_n, dtype=float))
        log_Bo = np.log(np.abs(np.array(odd_B)))
        Ao = np.column_stack([np.ones_like(log_no), -log_no])
        fit_o, _, _, _ = np.linalg.lstsq(Ao, log_Bo, rcond=None)
        Co, po = np.exp(fit_o[0]), fit_o[1]

        # Tail estimate
        n_start = max(all_B.keys()) + 1
        tail = 0.0
        for n in range(n_start, 100001):
            if n % 2 == 0:
                tail += Ce * n**(-pe)
            else:
                tail += Co * n**(-po)

        tail_f2s = tail / V_mag / schwinger
        corrected = B_cumul / V_mag / schwinger + tail_f2s
    else:
        Ce, pe, Co, po = 0, 0, 0, 0
        tail_f2s = 0
        corrected = B_cumul / V_mag / schwinger

    output = {
        'new_levels': new_levels,
        'all_B': {str(k): v for k, v in sorted(all_B.items())},
        'cumulative_B': B_cumul,
        'cumulative_F2_over_S': B_cumul / V_mag / schwinger,
        'tail_estimate': tail_f2s,
        'corrected_F2_over_S': corrected,
        'corrected_delta': corrected - 1.0,
        'power_law_even': {'C': float(Ce), 'p': float(pe)},
        'power_law_odd': {'C': float(Co), 'p': float(po)},
    }

    with open('debug/data/g2_extended_nint_v2.json', 'w') as f:
        json.dump(output, f, indent=2)


def _final_analysis(all_B, V_mag, schwinger):
    """Run tail estimate and PSLQ after all levels computed."""
    B_cumul = sum(all_B.values())
    F2_S = B_cumul / V_mag / schwinger

    ns = sorted(n for n in all_B.keys() if n >= 3)
    even_n = [n for n in ns if n % 2 == 0 and n >= 4]
    even_B = [all_B[n] for n in even_n]
    odd_n = [n for n in ns if n % 2 == 1 and n >= 3]
    odd_B = [all_B[n] for n in odd_n]

    log_ne = np.log(np.array(even_n, dtype=float))
    log_Be = np.log(np.abs(np.array(even_B)))
    Ae = np.column_stack([np.ones_like(log_ne), -log_ne])
    fit_e, _, _, _ = np.linalg.lstsq(Ae, log_Be, rcond=None)
    Ce, pe = np.exp(fit_e[0]), fit_e[1]

    log_no = np.log(np.array(odd_n, dtype=float))
    log_Bo = np.log(np.abs(np.array(odd_B)))
    Ao = np.column_stack([np.ones_like(log_no), -log_no])
    fit_o, _, _, _ = np.linalg.lstsq(Ao, log_Bo, rcond=None)
    Co, po = np.exp(fit_o[0]), fit_o[1]

    print(f"\n{'='*70}")
    print(f"  FINAL ANALYSIS (n_int=0..{max(all_B.keys())})")
    print(f"{'='*70}")
    print(f"\n  Power-law fits:")
    print(f"    Even: B ~ {Ce:.4e} * n^(-{pe:.2f})")
    print(f"    Odd:  B ~ {Co:.4e} * n^(-{po:.2f})")

    # Tail
    n_start = max(all_B.keys()) + 1
    tail = 0.0
    for n in range(n_start, 100001):
        if n % 2 == 0:
            tail += Ce * n**(-pe)
        else:
            tail += Co * n**(-po)

    tail_f2s = tail / V_mag / schwinger
    corrected = F2_S + tail_f2s
    delta = corrected - 1.0

    print(f"\n  Tail estimate (n_int={n_start}..100000):")
    print(f"    Tail B = {tail:.4e}")
    print(f"    Tail F2/S = {tail_f2s:.4e}")
    print(f"    Corrected F2/S = {corrected:.12f}")
    print(f"    delta = {delta:.12f}")
    print(f"    Precision: ~{-np.log10(abs(tail_f2s) / abs(delta)):.1f} digits")

    # Curvature coefficients
    lam2 = 6.25
    lam4 = lam2**2
    c1 = 0.5
    residual = delta - c1/lam2
    c2_apparent = residual * lam4

    print(f"\n  Curvature expansion:")
    print(f"    delta = c1/lam^2 + c2/lam^4 + ...")
    print(f"    c1 = 1/2 (Parker-Toms)")
    print(f"    Residual after c1: {residual:.12f}")
    print(f"    c2 (apparent) = {c2_apparent:.10f}")

    # PSLQ on c2
    try:
        import mpmath
        mpmath.mp.dps = 50
        c2_mp = mpmath.mpf(str(c2_apparent))

        print(f"\n  PSLQ identification of c2 = {float(c2_mp):.10f}:")

        bases = [
            (['c2', '1'], [c2_mp, mpmath.mpf(1)]),
            (['c2', '1', 'pi^2'], [c2_mp, mpmath.mpf(1), mpmath.pi**2]),
            (['c2', '1', 'pi'], [c2_mp, mpmath.mpf(1), mpmath.pi]),
            (['c2', '1', 'log2'], [c2_mp, mpmath.mpf(1), mpmath.log(2)]),
            (['c2', '1', 'G'], [c2_mp, mpmath.mpf(1), mpmath.catalan]),
            (['c2', '1', 'zeta3'], [c2_mp, mpmath.mpf(1), mpmath.zeta(3)]),
            (['c2', '1', 'pi^2', 'log2'], [c2_mp, mpmath.mpf(1), mpmath.pi**2, mpmath.log(2)]),
            (['c2', '1', 'pi^2', 'zeta3'], [c2_mp, mpmath.mpf(1), mpmath.pi**2, mpmath.zeta(3)]),
            (['c2', '1', 'pi', 'log2'], [c2_mp, mpmath.mpf(1), mpmath.pi, mpmath.log(2)]),
        ]

        for labels, vals in bases:
            try:
                rel = mpmath.pslq(vals, tol=1e-5, maxcoeff=500)
                if rel is not None and rel[0] != 0:
                    terms = [f"{r}*{l}" for r, l in zip(rel, labels) if r != 0]
                    ratio = sum(-r*float(v) for r, v in zip(rel[1:], vals[1:])) / rel[0]
                    err = abs(ratio - float(c2_mp))
                    print(f"    Basis {labels[1:]}: {' + '.join(terms)} = 0")
                    print(f"      => c2 = {ratio:.12f} (err {err:.2e})")
                else:
                    print(f"    Basis {labels[1:]}: no relation")
            except Exception as e:
                print(f"    Basis {labels[1:]}: error ({e})")

        # Also try PSLQ on delta directly
        delta_mp = mpmath.mpf(str(delta))
        print(f"\n  PSLQ on delta = {float(delta_mp):.12f}:")
        delta_bases = [
            (['delta', '1', 'pi^2'], [delta_mp, mpmath.mpf(1), mpmath.pi**2]),
            (['delta', '1', 'pi', 'pi^2'], [delta_mp, mpmath.mpf(1), mpmath.pi, mpmath.pi**2]),
        ]
        for labels, vals in delta_bases:
            try:
                rel = mpmath.pslq(vals, tol=1e-6, maxcoeff=500)
                if rel is not None and rel[0] != 0:
                    terms = [f"{r}*{l}" for r, l in zip(rel, labels) if r != 0]
                    ratio = sum(-r*float(v) for r, v in zip(rel[1:], vals[1:])) / rel[0]
                    err = abs(ratio - float(delta_mp))
                    print(f"    Basis {labels[1:]}: {' + '.join(terms)} = 0")
                    print(f"      => delta = {ratio:.12f} (err {err:.2e})")
                else:
                    print(f"    Basis {labels[1:]}: no relation")
            except Exception as e:
                print(f"    Basis {labels[1:]}: error ({e})")

    except ImportError:
        print("  (mpmath not available for PSLQ)")

    # Save final
    _save_checkpoint(all_B, [], V_mag, schwinger)
    print(f"\nFinal results saved to debug/data/g2_extended_nint_v2.json")


if __name__ == '__main__':
    main()
