r"""Decisive gate (b): mu<=2 (delta channel) at (2,2), recurrence vs grid,
plus an X-block accuracy re-check after the seed-breakpoint change."""
import time
import numpy as np
from geovac import neumann_vee_general_m as gm
from geovac import prolate_general_m as pg

E_EXACT = -1.174475
DE = pg.DE_EXACT
dep = lambda e: 100.0 * (-1.0 - e) / DE
P = lambda *a: print(*a, flush=True)

# --- X-block accuracy (semi-independent high-precision references) ---
REF = {(2, 2, 2, 0, 0): 0.0100142205252448,
       (4, 4, 4, 0, 0): 87.70947318323068,
       (6, 4, 4, 0, 0): 2144.07536568145,
       (8, 4, 4, 0, 2): 114762.1627845131}
Xt = gm.build_Xtab([(2, 2), (4, 4)], l_neumann=8, p_max=2, basis_alpha=1.0)
P("=== X-block accuracy after seed change ===")
worst = 0.0
for (l, m, s, P1, P2), r in REF.items():
    got = Xt[(l, m, s)][P1, P2]
    rel = abs(got - r) / abs(r)
    worst = max(worst, rel)
    P(f"  X[{l},{m},{s}]({P1},{P2}) rel err = {rel:.2e}")
P(f"  worst = {worst:.2e}")


def energy(engine, jm, lm, mm, alpha=1.0, ln=14):
    basis = pg.generate_basis(jm, lm, mm, alpha)
    mom = pg.Moments(2.0 * alpha, 6 * max(jm, lm) + 6 * (mm + 2) + 20)
    s, h1 = pg.one_body(basis, pg.R_DEFAULT, 1.0, mom)
    t = time.time()
    v = gm.vee_matrix(basis, pg.R_DEFAULT, ln) if engine == "recur" \
        else pg.vee_matrix(basis, pg.R_DEFAULT, pg.XiGrid(alpha), ln)
    dt = time.time() - t
    h = h1 + v + (1.0 / pg.R_DEFAULT) * s
    e, nk, nt = pg.solve_generalized(h, s)
    return e, len(basis), nk, nt, dt


P("\n=== GATE (b): mu<=2 (2,2), recurrence vs grid ===")
e1, *_ = energy("recur", 2, 2, 1)
P(f"  mu<=1 recur: E={e1:.6f}  {dep(e1):.2f}%")
e2, n2, nk2, nt2, dt2 = energy("recur", 2, 2, 2)
P(f"  mu<=2 recur: E={e2:.6f}  {dep(e2):.2f}%  ({nk2}/{nt2})  build {dt2:.0f}s "
  f" variational={e2>E_EXACT}")
eg, ng, nkg, ntg, _ = energy("grid", 2, 2, 2)
P(f"  mu<=2 grid:  E={eg:.4f}  {dep(eg):.1f}%  variational={eg>E_EXACT}")
P(f"  delta gain (recur, mu<=1 -> mu<=2): {1e3*(e1-e2):.3f} mHa")
