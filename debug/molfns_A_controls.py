"""MF step A: the mixed finite-size x Breit-Pauli energy Q(rn), and the hydrogen controls C1, C2.

Record: docs/predictions/2026-10-20_molecular_fns.md (committed 01f79fa5 before computing).

E_mixed(rn) = 2 <H_BP psi | eta> + (1/(8c^2)) <lap dV>,  eta = G' (dV psi),  dV = sphere - point at d.
  2 <H_BP psi|eta> = 2 (pi/(2c^2)) sum_a Z_a psi(R_a) eta(R_a)
                     - (1/(4c^2)) [ 4 <(E-V)psi|(E-V)eta> - 4 <(E-V)psi| Q dV psi> ]
  (1/(8c^2)) <lap dV> = (1/(8c^2)) 4 pi Z_d [ int rho_sphere psi^2 - psi(0)^2 ]
Q(rn) = c^2 E_mixed / <dV>  (c cancels).  S difference between two states = -(Q_a - Q_b).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from mpmath import mp, digamma, euler, log, mpf

sys.path.insert(0, str(Path(__file__).parent))
import molfns_solver as S  # noqa: E402

OUT = Path(__file__).parent / "data" / "molfns_A.json"


def mixed_Q(sysm: S.System, E: float, psi: np.ndarray, rn: float) -> float:
    g = sysm.grid
    r = g.r
    Zd = sysm.Zd
    dV = np.where(r < rn, -Zd * (3 * rn**2 - r**2) / (2 * rn**3) + Zd / r, 0.0)
    d = sysm.split(psi)
    s = sysm.join({l: dV * w for l, w in d.items()})
    eta = S.reduced_solve(sysm, E, psi, s)
    psi0 = S.value_at_origin(sysm, psi, rn / 100, rn / 5)
    eta0 = S.value_at_origin(sysm, eta, rn / 100, rn / 5)
    T1 = np.pi * Zd * psi0 * eta0
    for (Z, R) in sysm.others:
        T1 += np.pi * Z * S.value_on_axis(sysm, psi, R) * S.value_on_axis(sysm, eta, R)
    EVpsi = E * psi - sysm.apply_V(psi)
    EVeta = E * eta - sysm.apply_V(eta)
    Qs = s - psi * sysm.inner(psi, s)
    T2 = -(sysm.inner(EVpsi, EVeta) - sysm.inner(EVpsi, Qs))
    # int_0^rn psi^2 d^3r = sum_l int_0^rn u_l^2 dr, trapezoid in t ending on the grid point rn
    k = int(np.argmin(abs(r - rn)))
    e2t = np.exp(2 * g.t[: k + 1])
    wts = np.full(k + 1, g.h); wts[-1] *= 0.5
    inside = sum(float(np.sum(wts * e2t * w[: k + 1] ** 2)) for w in d.values())
    T3 = 0.5 * np.pi * Zd * (3 / (4 * np.pi * rn**3) * inside - psi0**2)
    dVexp = sysm.inner(psi, s)
    return float((T1 + T2 + T3) / dVexp)


def S_exact(n: int) -> float:
    mp.dps = 30
    return float(-log(n) + digamma(n) + euler - mpf(5 * n + 9) * (n - 1) / (4 * n * n))


def hydrogen(Z: float, h: float, rns):
    g = S.Grid(t0=-40.0, t1=np.log(150.0), h=h)
    sysm = S.System(grid=g, lmax=0, Zd=Z)
    vals, vecs = S.ground_state(sysm, sigma=-Z**2 / 2 - 0.05, k=3)
    out = {}
    for rn0 in rns:
        rn = g.snap(rn0)
        Qs = [mixed_Q(sysm, vals[j], vecs[:, j], rn) for j in range(3)]
        out[f"{rn:.3e}"] = {"Q": Qs, "S21": -(Qs[1] - Qs[0]) / Z**2, "S31": -(Qs[2] - Qs[0]) / Z**2}
    return vals.tolist(), out


if __name__ == "__main__":
    rec = {"S21_exact": S_exact(2) - S_exact(1), "S31_exact": S_exact(3) - S_exact(1)}
    print("exact", rec)
    for Z in (1.0, 2.0):
        for h in (0.02, 0.01):
            E, res = hydrogen(Z, h, [1e-3, 1e-4, 1e-5])
            print(f"Z={Z} h={h} E={E}")
            for k, v in res.items():
                print(f"   rn={k}  S21={v['S21']:.7f}  S31={v['S31']:.7f}  Q1s={v['Q'][0]:.5f}", flush=True)
            rec[f"Z{Z}_h{h}"] = {"E": E, "res": res}
    OUT.write_text(json.dumps(rec, indent=1))
