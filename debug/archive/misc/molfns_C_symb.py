"""MF step C: the far-nucleus Darwin term through the symmetric reduced Green's function.

Record: docs/predictions/2026-10-20_molecular_fns.md (committed 01f79fa5 before computing).
In step B the term pi Z_b psi(R_b) eta(R_b) uses point values at nucleus b from a single-centre sum,
which converges like 1/l at a cusp.  Since eta = G'(dV psi) and dV lives within rn of the origin,
eta(R_b) -> G'(R_b, 0) int dV psi, and G' is symmetric, so this term contributes pi Z_b G'(0, R_b) to
Q, independent of rn.  G'(0, R_b) = chi(0), chi = G' delta_b, read at the expansion centre.
Everything else is step B's mixed_Q with the b-term removed.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import molfns_solver as S  # noqa: E402
import molfns_A_controls as MA  # noqa: E402

OUT = Path(__file__).parent / "data"


def delta_b_source(sysm: S.System, R: float) -> np.ndarray:
    """Unit delta at R z-hat as w-values: channel l weight sqrt((2l+1)/(4 pi)) delta(r-R)/R^2."""
    g = sysm.grid
    k = int(np.argmin(abs(g.r - R)))
    q = {}
    for l in sysm.channels():
        v = np.zeros(g.n)
        v[k] = np.sqrt((2 * l + 1) / (4 * np.pi)) * np.exp(0.5 * g.t[k]) / (R * np.exp(2 * g.t[k]) * g.h)
        q[l] = v
    return sysm.join(q)


def Q_sym(sysm: S.System, E: float, vs, rn: float, Gb0: float) -> float:
    """Q(rn) with G'' (vs = [g, u] projected out) and the far-nucleus Darwin term as pi Z_b G''(0, R_b)."""
    psi = vs[0]
    g = sysm.grid
    r = g.r
    Zd = sysm.Zd
    dV = np.where(r < rn, -Zd * (3 * rn**2 - r**2) / (2 * rn**3) + Zd / r, 0.0)
    d = sysm.split(psi)
    s = sysm.join({l: dV * w for l, w in d.items()})
    eta = S.reduced_solve_multi(sysm, E, vs, s)
    psi0 = S.value_at_origin(sysm, psi, rn / 100, rn / 5)
    eta0 = S.value_at_origin(sysm, eta, rn / 100, rn / 5)
    T1 = np.pi * Zd * psi0 * eta0
    EVpsi = E * psi - sysm.apply_V(psi)
    EVeta = E * eta - sysm.apply_V(eta)
    Qs = s.copy()
    for v in vs:
        Qs = Qs - v * sysm.inner(v, s)
    T2 = -(sysm.inner(EVpsi, EVeta) - sysm.inner(EVpsi, Qs))
    k = int(np.argmin(abs(r - rn)))
    e2t = np.exp(2 * g.t[: k + 1])
    wts = np.full(k + 1, g.h); wts[-1] *= 0.5
    inside = sum(float(np.sum(wts * e2t * w[: k + 1] ** 2)) for w in d.values())
    T3 = 0.5 * np.pi * Zd * (3 / (4 * np.pi * rn**3) * inside - psi0**2)
    dVexp = sysm.inner(psi, s)
    Zb = sum(Z for Z, _ in sysm.others)
    return float((T1 + T2 + T3) / dVexp + np.pi * Zb * Gb0)


def run(R0: float, lmax: int, h: float, rns, t0: float = -30.0, rmax: float = 45.0, project_u: bool = True):
    g = S.Grid(t0=t0, t1=np.log(rmax), h=h, anchor=R0)
    R = g.snap(R0)
    mol = S.System(grid=g, lmax=lmax, Zd=1.0, others=[(1.0, R)])
    t = time.time()
    vals, vecs = S.ground_state(mol, sigma=-1.0 / R - 0.75, k=2)
    E, psi, u = float(vals[0]), vecs[:, 0], vecs[:, 1]
    # u should be odd under the swap of the nuclei: u(R_b) = -u(0); psi gerade: psi(R_b) = psi(0)
    par = (S.value_on_axis(mol, u, R) / S.value_at_origin(mol, u, 1e-6, 1e-4),
           S.value_on_axis(mol, psi, R) / S.value_at_origin(mol, psi, 1e-6, 1e-4))
    vs = [psi, u] if project_u else [psi]
    chi = S.reduced_solve_multi(mol, E, vs, delta_b_source(mol, R))
    Gb0 = S.value_at_origin(mol, chi, 1e-6, 1e-4)
    atom = S.System(grid=g, lmax=0, Zd=1.0)
    av, avec = S.ground_state(atom, sigma=-0.55, k=1)
    out = {"R": R, "lmax": lmax, "h": h, "E": E, "Eu": float(vals[1]), "parity_ratios": par,
           "Gb0": Gb0, "project_u": project_u, "rows": {}}
    for rn0 in rns:
        rn = g.snap(rn0)
        qm = Q_sym(mol, E, vs, rn, Gb0)
        qa = MA.mixed_Q(atom, float(av[0]), avec[:, 0], rn)
        out["rows"][f"{rn:.3e}"] = {"Qmol": qm, "QH": qa, "S": -(qm - qa)}
    out["sec"] = time.time() - t
    return out


if __name__ == "__main__":
    tag = sys.argv[1]
    plan = {"mfb": [(2.0, [32, 48, 64]), (4.0, [32, 48, 64]), (8.0, [32, 48, 64])]}[tag]
    rec = {}
    for R0, lms in plan:
        for lmax in lms:
            for pu in ([True, False] if (R0 == 2.0 and lmax == 32) else [True]):
                o = run(R0, lmax, 0.005, [1e-3, 3e-4, 1e-4], project_u=pu)
                Ss = [v["S"] for v in o["rows"].values()]
                print(f"R={o['R']:.4f} lmax={lmax} proj_u={pu} E={o['E']:.10f} Eu={o['Eu']:.6f} "
                      f"parity={o['parity_ratios'][0]:+.3f},{o['parity_ratios'][1]:+.3f} Gb0={o['Gb0']:.6f} "
                      f"S(rn)={['%.5f' % x for x in Ss]} {o['sec']:.0f}s", flush=True)
                rec[f"R{R0}_l{lmax}_pu{int(pu)}"] = o
                (OUT / f"molfns_C_{tag}.json").write_text(json.dumps(rec, indent=1))
