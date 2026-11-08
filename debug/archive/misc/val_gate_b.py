r"""Gate (b) + gate (a) recur-vs-grid, using the fast cached geovac engine."""
import time
import numpy as np
from geovac import neumann_vee_general_m as gm
from geovac import prolate_general_m as pg

E_EXACT = -1.174475
DE = pg.DE_EXACT


def dep(e):
    return 100.0 * (-1.0 - e) / DE


def energy(engine, jm, lm, mm, alpha=1.0, ln=14):
    basis = pg.generate_basis(jm, lm, mm, alpha)
    mom = pg.Moments(2.0 * alpha, 6 * max(jm, lm) + 6 * (mm + 2) + 20)
    s, h1 = pg.one_body(basis, pg.R_DEFAULT, 1.0, mom)
    t = time.time()
    if engine == "recur":
        v = gm.vee_matrix(basis, pg.R_DEFAULT, ln)
    else:
        v = pg.vee_matrix(basis, pg.R_DEFAULT, pg.XiGrid(alpha), ln)
    dt = time.time() - t
    h = h1 + v + (1.0 / pg.R_DEFAULT) * s
    e, nk, nt = pg.solve_generalized(h, s)
    return e, len(basis), nk, nt, dt


P = lambda *a: print(*a, flush=True)

P("=== GATE (a) recur vs grid at matched alpha=1.0 (grid is sound at mu<=1) ===")
for (jm, lm, mm) in [(2, 2, 1), (3, 3, 1)]:
    er, n, nkr, nt, dtr = energy("recur", jm, lm, mm)
    eg, _, nkg, _, dtg = energy("grid", jm, lm, mm)
    P(f"  ({jm},{lm}) mu<={mm}: recur E={er:.6f} ({dep(er):.2f}%, {nkr}/{nt}, {dtr:.0f}s)  "
      f"grid E={eg:.6f} ({dep(eg):.2f}%)  |dE|={1e6*abs(er-eg):.1f} uHa")

P("\n=== GATE (b) mu<=2 (delta) via recurrence -- grid blows up to -21 Ha ===")
for (jm, lm) in [(2, 2), (3, 3)]:
    er, n, nkr, nt, dtr = energy("recur", jm, lm, 2)
    eg, _, nkg, _, dtg = energy("grid", jm, lm, 2)
    P(f"  ({jm},{lm}) mu<=2: recur E={er:.6f} ({dep(er):.2f}%, {nkr}/{nt}, {dtr:.0f}s)  "
      f"grid E={eg:.4f} ({dep(eg):.1f}%)  variational(recur)={er>E_EXACT}")

P("\n=== delta gain at fixed basis (recur): mu<=1 -> mu<=2 ===")
for (jm, lm) in [(2, 2), (3, 3)]:
    e1, *_ = energy("recur", jm, lm, 1)
    e2, *_ = energy("recur", jm, lm, 2)
    P(f"  ({jm},{lm}): mu<=1 {dep(e1):.2f}%  mu<=2 {dep(e2):.2f}%  delta gain {1e3*(e1-e2):.2f} mHa")
