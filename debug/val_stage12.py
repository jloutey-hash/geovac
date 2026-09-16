r"""Stage 1-2 validation: recurrence X table vs the driver's grid-based Xtab, mu<=1.

The driver (prolate_ci_general_m.vee_matrix) builds Xtab by grid quadrature with
q_deriv (differentiation), which is STABLE at mu<=1.  We rebuild the exact same
grid Xtab for a chosen (m_set, s_set, l_neumann, p_max) and compare it element by
element against general_m_engine.build_Xtab.
"""
import sys; sys.path.insert(0, 'debug')
import numpy as np
from numpy.polynomial import polynomial as P
from numpy.polynomial import legendre as L
import prolate_ci_general_m as drv
import general_m_engine as eng


def driver_Xtab(m_set, s_set, l_neumann, p_max, alpha):
    """Reproduce prolate_ci_general_m.vee_matrix's Xtab block build (grid route)."""
    grid = drv.XiGrid(alpha)
    xi = grid.xi
    decay = np.exp(-2.0 * alpha * xi)
    fcache = {}
    for s in s_set:
        xp = drv.poly_xi2m1(s)
        for Pp in range(p_max + 1):
            fcache[(s, Pp)] = P.polyval(xi, drv.shift(xp, Pp)) * decay
    Xtab = {}
    for l in range(l_neumann + 1):
        for m in m_set:
            if m > l:
                continue
            pv = P.polyval(xi, drv.legendre_deriv_poly(l, m))
            qv = drv.q_deriv(l, m, xi)
            for s in s_set:
                lo = np.empty((p_max + 1, len(xi)))
                hi = np.empty((p_max + 1, len(xi)))
                for Pp in range(p_max + 1):
                    f = fcache[(s, Pp)]
                    lo[Pp] = grid.cumulative(f * pv)
                    cq = grid.cumulative(f * qv)
                    hi[Pp] = cq[-1] + (grid.integral(f * qv) - cq[-1]) - cq
                mat = np.empty((p_max + 1, p_max + 1))
                for P1 in range(p_max + 1):
                    g1 = fcache[(s, P1)]
                    for P2 in range(p_max + 1):
                        mat[P1, P2] = grid.integral(g1 * (qv * lo[P2] + pv * hi[P2]))
                Xtab[(l, m, s)] = mat
    return Xtab


if __name__ == "__main__":
    alpha = 1.0
    # mu<=1 sector: m in {0,1,2}, s in {0,1,2}
    ms_pairs = sorted(set(
        (m, (a + b + m) // 2)
        for a in (0, 1) for b in (0, 1)
        for m in (a + b, abs(a - b))
        if (a + b + m) % 2 == 0
    ))
    m_set = sorted(set(m for m, _ in ms_pairs))
    s_set = sorted(set(s for _, s in ms_pairs))
    l_neumann, p_max = 8, 8
    print("ms_pairs (m,s) =", ms_pairs)
    print("m_set =", m_set, " s_set =", s_set)

    Xd = driver_Xtab(m_set, s_set, l_neumann, p_max, alpha)
    Xr = eng.build_Xtab(ms_pairs, l_neumann, p_max, alpha)

    print("\n block | max|abs diff| | max|rel diff| (over |X|>1e-8) | ||X_recur||")
    worst_rel = 0.0
    for (l, m, s) in sorted(Xr.keys()):
        A = Xr[(l, m, s)]
        B = Xd[(l, m, s)]
        num = np.abs(A - B)
        den = np.maximum(np.abs(B), 1e-8)
        mask = np.abs(B) > 1e-8
        rel = (num[mask] / den[mask]).max() if mask.any() else 0.0
        worst_rel = max(worst_rel, rel)
        print(f"  ({l},{m},{s}) | {num.max():.3e} | {rel:.3e} | {np.linalg.norm(A):.4e}")
    print(f"\n WORST rel diff over all mu<=1 X blocks = {worst_rel:.3e}")
