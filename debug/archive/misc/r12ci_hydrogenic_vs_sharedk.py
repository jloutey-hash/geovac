"""Part 2, decided cheaply: the base kappa wall is a property of the SHARED-k
Sturmian basis, not of the R12-CI.  Coulomb-Sturmians with k_n = Z/n ARE the
hydrogenic radials (same functional form), and hydrogenic radials with a common Z
are mutually orthonormal.  So swap the orbital scale rule and re-measure.

Reports, at matched function count:  kappa(S_1e), kappa(pair block), and the
He energy with and without a geminal.  (debug/ only, READ-ONLY use of the engine.)
"""
import importlib.util
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)

_spec = importlib.util.spec_from_file_location(
    "ctf12_r12ci_he", os.path.join(HERE, "ctf12_r12ci_he.py"))
ctf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ctf)

EXACT_HE, Z = ctf.EXACT_HE, 2.0
NG, NX = 700, 96            # wider grid: hydrogenic n=5 has k=Z/5=0.40
R_MAX = 42.0 / 0.40
R, WR = ctf.make_grid(0.40, NG, r_max=R_MAX)
W = R * R * WR


def tabs_shared_k(ns, k):
    Rt, dRt = {}, {}
    for n in range(1, ns + 1):
        Rt[n], dRt[n] = ctf.R_and_dR(n, R, k)
    return Rt, dRt


def tabs_hydrogenic(ns, Zc=Z):
    """Coulomb-Sturmian at k_n = Z/n  ==  the hydrogenic radial R_{n0}."""
    Rt, dRt = {}, {}
    for n in range(1, ns + 1):
        Rt[n], dRt[n] = ctf.R_and_dR(n, R, Zc / n)
    return Rt, dRt


def kappa_norm(M):
    d = np.sqrt(np.diag(M)); S = M / np.outer(d, d)
    ev = np.linalg.eigvalsh(S)
    return float(ev.max() / ev.min())


def s1e(Rt, ns):
    Rm = np.array([Rt[n] for n in range(1, ns + 1)])
    S = (Rm * W) @ Rm.T
    d = np.sqrt(np.diag(S))
    return S / np.outer(d, d)


_KCACHE = {}


def kernels(gamma):
    if gamma not in _KCACHE:
        _KCACHE[gamma] = ctf.build_kernels(R, gamma, nx=NX)
    return _KCACHE[gamma]


def run(Rt, dRt, ns, n_gem, gamma):
    K = kernels(gamma)
    bfs = []
    for i in range(1, ns + 1):
        for j in range(i, ns + 1):
            bfs.append(ctf.orbital_pair(i, j, Rt, dRt))
    for ref in range(1, n_gem + 1):
        bfs.append(ctf.geminal(ref, Rt, dRt))
    H, M = ctf.assemble(bfs, R, WR, K, Z=Z)
    ev, U = np.linalg.eigh(M)
    keep = ev > 1e-10
    X = U[:, keep] / np.sqrt(ev[keep])
    E = float(np.linalg.eigvalsh(X.T @ H @ X)[0])
    return E, kappa_norm(M), len(bfs)


out = {"exact": EXACT_HE, "grid": dict(Ng=NG, r_max=R_MAX), "rows": []}
print("=" * 92)
print("shared-k Coulomb-Sturmian (k=1.7)   vs   hydrogenic scale rule (k_n = Z/n)")
print("=" * 92)
print(f"{'ns':>3} | {'kap(S1e)':>9}{'kap(pair)':>10}{'E(no gem)':>12}{'err_mHa':>9}"
      f" | {'kap(S1e)':>9}{'kap(pair)':>10}{'E(no gem)':>12}{'err_mHa':>9}")
for ns in range(2, 6):
    Rs, dRs = tabs_shared_k(ns, 1.7)
    Rh, dRh = tabs_hydrogenic(ns)
    ks = kappa_norm(s1e(Rs, ns)); kh = kappa_norm(s1e(Rh, ns))
    Es, kps, _ = run(Rs, dRs, ns, 0, 0.7)
    Eh, kph, _ = run(Rh, dRh, ns, 0, 0.7)
    print(f"{ns:>3} | {ks:>9.2f}{kps:>10.2f}{Es:>12.6f}{(Es-EXACT_HE)*1000:>9.2f}"
          f" | {kh:>9.4f}{kph:>10.4f}{Eh:>12.6f}{(Eh-EXACT_HE)*1000:>9.2f}")
    out["rows"].append(dict(ns=ns, shared_k=dict(kappa_S1e=ks, kappa_pair=kps, E=Es,
                                                 err_mHa=(Es - EXACT_HE) * 1000),
                            hydrogenic=dict(kappa_S1e=kh, kappa_pair=kph, E=Eh,
                                            err_mHa=(Eh - EXACT_HE) * 1000)))

print()
print("with one geminal (gamma scan), ns=3:")
print(f"{'gamma':>6} | {'shared-k  E':>12}{'err':>8}{'kappa':>9}"
      f" | {'hydrogenic E':>13}{'err':>8}{'kappa':>9}")
out["geminal_ns3"] = []
for gamma in (0.4, 0.5, 0.7, 1.1, 2.0):
    Rs, dRs = tabs_shared_k(3, 1.7)
    Rh, dRh = tabs_hydrogenic(3)
    Es, ks, _ = run(Rs, dRs, 3, 1, gamma)
    Eh, kh, _ = run(Rh, dRh, 3, 1, gamma)
    print(f"{gamma:>6.2f} | {Es:>12.6f}{(Es-EXACT_HE)*1000:>8.2f}{ks:>9.1f}"
          f" | {Eh:>13.6f}{(Eh-EXACT_HE)*1000:>8.2f}{kh:>9.1f}")
    out["geminal_ns3"].append(dict(gamma=gamma,
                                   shared_k=dict(E=Es, err_mHa=(Es - EXACT_HE) * 1000, kappa=ks),
                                   hydrogenic=dict(E=Eh, err_mHa=(Eh - EXACT_HE) * 1000, kappa=kh)))

with open("debug/data/r12ci_hydrogenic_vs_sharedk.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/r12ci_hydrogenic_vs_sharedk.json")
