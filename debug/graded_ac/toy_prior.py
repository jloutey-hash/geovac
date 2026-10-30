"""Toy prior for the pre-registration (flat 2-torus only; no S^3 value).

Measures, at the toy analogue of the scored points (window R^2 = (n + 1/2)^2,
n = 2..5, t = tau / R^2, tau in {1, 4}): c, |c/4 - 1| and the pure-gauge
companion ratio rho_pg = |k_A^pg / k_A| (a = 1, same b).  These numbers are the
measured controls that justify the bounds 1e-3 and 0.1 in the draft.
Log: debug/graded_ac/toy_prior.log
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from graded_product import (eps2_heat_coefficient, kinetic_duhamel, torus_continuum_gauge,  # noqa: E402
                            torus_continuum_higgs_kin, torus_window)

LOG = []


def log(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def ops(R2, p, q):
    W = torus_window(R2)
    D = W.D_exact().to_numpy()
    a = W.shift_exact(*q).to_numpy()
    b = W.shift_exact(*p).to_numpy()
    om = a @ (D @ b - b @ D)
    A = om + om.conj().T
    phi = a @ b
    C = D @ phi - phi @ D
    pg = D @ b - b @ D
    Apg = pg + pg.conj().T
    return W, D, A, C, Apg


def main() -> None:
    p, q = (1, 0), (0, 1)
    log("toy_prior: flat T^2, pair a = e^{i q.x}, b = e^{i p.x}, p=(1,0), q=(0,1); pure gauge a = 1, same b")
    log("  tau  n  R^2     t          dim   c              |c/4-1|     rho_pg")
    for tau in (1.0, 4.0):
        for n in (2, 3, 4, 5):
            R2 = (n + 0.5) ** 2
            t = tau / R2
            W, D, A, C, Apg = ops(R2, p, q)
            kA = eps2_heat_coefficient(D, A, t)
            kP = kinetic_duhamel(D, C, t)
            kpg = eps2_heat_coefficient(D, Apg, t)
            c = 4 * (kA / torus_continuum_gauge(p, q, t)) / (kP / torus_continuum_higgs_kin(p, q, t))
            log(f"  {tau:<4} {n}  {R2:<6}  {t:.6f}   {W.dim:<4}  {c:+.8f}   {abs(c / 4 - 1):.3e}   {abs(kpg / kA):.3e}")
    (HERE / "toy_prior.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
