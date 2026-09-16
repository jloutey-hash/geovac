r"""Fire tests for tests/test_paper12_general_m_neumann.py.

Each guard in that file must REJECT the specific wrong answer it names.  Here we
construct that wrong answer and confirm the guard's assertion fails (fires).  A
guard whose rejected answer cannot be produced is not a guard (Sec. 9).

Run: python debug/firetest_p12_general_m.py
"""
import sys
import numpy as np

sys.path.insert(0, 'tests')
from geovac import neumann_vee_general_m as gm
from geovac import prolate_general_m as pg

E_EXACT = -1.174475
PASS, FAIL = "FIRES (guard works)", "DID NOT FIRE (guard is inert)"
results = []


def check(name, fired):
    results.append((name, fired))
    print(f"  [{'OK' if fired else 'XX'}] {name}: {PASS if fired else FAIL}",
          flush=True)


# 1. mu=0 reduction guard: a corrupted V_ee (one element scaled) must trip rel<1e-7
def fire_mu0():
    from geovac.hylleraas import HylleraasBasisFunction
    from geovac.neumann_vee import compute_vee_matrix_neumann
    a = 1.0
    mine = pg.generate_basis(2, 2, 0, a)
    theirs = [HylleraasBasisFunction(b.j, b.k, b.l, b.m, 0, a) for b in mine]
    v_mine = gm.vee_matrix(mine, pg.R_DEFAULT, l_neumann=14).copy()
    v_them = compute_vee_matrix_neumann(theirs, pg.R_DEFAULT, l_max=20)
    v_mine[3, 3] *= 1.01                       # inject a 1% error
    rel = np.abs(v_mine - v_them) / np.maximum(np.abs(v_them), 1e-12)
    check("mu0_reduction", not (rel.max() < 1e-7))


# 2. X-reference guard: a perturbed engine value must trip rel<1e-9
def fire_Xref():
    ref = 87.70947318323068                     # X[4,4,4](0,0)
    Xtab = gm.build_Xtab([(4, 4)], l_neumann=4, p_max=0, basis_alpha=1.0)
    got = Xtab[(4, 4, 4)][0, 0] * (1 + 1e-6)     # perturb by 1e-6
    check("X_reference", not (abs(got - ref) / abs(ref) < 1e-9))


# 3. mu=2 stability guard: the differentiation blow-up (grid engine) must trip
#    the variational assertion E > E_EXACT
def fire_mu2():
    basis = pg.generate_basis(2, 2, 2, 1.0)
    mom = pg.Moments(2.0, 6 * 2 + 6 * 4 + 20)
    s, h1 = pg.one_body(basis, pg.R_DEFAULT, 1.0, mom)
    v = pg.vee_matrix(basis, pg.R_DEFAULT, pg.XiGrid(1.0), 14)   # differentiation
    h = h1 + v + (1.0 / pg.R_DEFAULT) * s
    e, _k, _t = pg.solve_generalized(h, s)
    check(f"mu2_stability (grid E={e:.3f})", not (e > E_EXACT))


# 4. intact-weight guard: a monomial-expanded (bare) delta moment is divergent;
#    emulate the wrong path and confirm it is non-finite / non-positive
def fire_intact_weight():
    import mpmath as mp
    with mp.workdps(30):
        # the WRONG path: expand (xi^2-1)^s -> the bare xi^p d^m Q_l moment, which
        # diverges for m>=1 (the pieces the intact weight must keep together)
        m, l = 4, 4
        def bare(u):
            x = 1 + u
            return x ** 0 * mp.e ** (-mp.mpf(2.0) * x) * gm._RQ_mp(l, m, x)
        try:
            val = mp.quad(bare, [0, mp.mpf('0.001'), mp.mpf('0.1'), 1, mp.inf])
            fired = not (mp.isfinite(val) and abs(val) < mp.mpf(10) ** 8)
        except Exception:
            fired = True
        check("intact_weight (bare moment divergent)", fired)


if __name__ == "__main__":
    print("=== fire tests: each guard must reject its named wrong answer ===",
          flush=True)
    fire_mu0()
    fire_Xref()
    fire_mu2()
    fire_intact_weight()
    n_fired = sum(f for _, f in results)
    print(f"\n {n_fired}/{len(results)} guards fired against their wrong answers.",
          flush=True)
    sys.exit(0 if n_fired == len(results) else 1)
