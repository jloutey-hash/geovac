"""Is the R12-CI geminal width DETERMINED?  (Avery "solve-and-tabulate", variational side.)

The v5.0.12 arc tested determined/tabulated gamma on the CHEAP non-Hermitian dressing
(negative).  Its own diagnosis: F12's fixed gamma works only VARIATIONALLY, and that route
was set aside as gated by kappa(S)~200 -- which the 2026-08-25 probe re-scoped to an
ENCODING cost (classically Lowdin is free: same span, same energy).  So test it here.

Transferability probe on the He ISOELECTRONIC SERIES.  Under Z-scaling,
    H = Z^2 [ h(rho) + (1/Z) 1/rho_12 ],   rho = Z r,
so a Z-independent *scaled* geminal width means  gamma_opt ∝ Z  EXACTLY.
  gamma_opt/Z constant  =>  determined law (physics), not a fitted table.
  gamma_opt/Z scatters  =>  a fit; stop.

Self-validation gate: Z=2 must reproduce the recorded He R12-CI operating point
(0.31-0.80 mHa at ns=3, n_gem=1).
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

# nonrelativistic infinite-nuclear-mass ground states, 1s^2 isoelectronic series
# (Pekeris / Drake high-accuracy values)
EXACT = {1: -0.5277510165, 2: -2.9037243770, 3: -7.2799134126, 4: -13.6555662384}

NG, NX = 800, 96
KMIN = 0.80                      # most diffuse case is Z=1 (k ~ 0.85)
R_MAX = 42.0 / KMIN
R, WR = ctf.make_grid(KMIN, NG, r_max=R_MAX)

_KC = {}


def kernels(gamma):
    if gamma not in _KC:
        _KC[gamma] = ctf.build_kernels(R, gamma, nx=NX)
    return _KC[gamma]


def energy(ns, k, gamma, Z, n_gem=1):
    Rt, dRt = {}, {}
    for n in range(1, ns + 1):
        Rt[n], dRt[n] = ctf.R_and_dR(n, R, k)
    bfs = []
    for i in range(1, ns + 1):
        for j in range(i, ns + 1):
            bfs.append(ctf.orbital_pair(i, j, Rt, dRt))
    for ref in range(1, n_gem + 1):
        bfs.append(ctf.geminal(ref, Rt, dRt))
    H, M = ctf.assemble(bfs, R, WR, kernels(gamma), Z=float(Z))
    ev, U = np.linalg.eigh(M)
    keep = ev > 1e-10
    X = U[:, keep] / np.sqrt(ev[keep])
    return float(np.linalg.eigvalsh(X.T @ H @ X)[0])


NS = 3
GAMMAS = [0.10, 0.15, 0.20, 0.30, 0.40, 0.55, 0.70, 0.90, 1.20, 1.60, 2.10, 2.80]
KFACS = [0.75, 0.85, 0.95]        # k = kfac * Z  (He optimum was k=1.7 at Z=2 -> 0.85)

out = {"ns": NS, "exact": EXACT, "grid": dict(Ng=NG, nx=NX, r_max=R_MAX),
       "gammas": GAMMAS, "kfacs": KFACS, "rows": [], "surface": {}}

print("=" * 84)
print(f"R12-CI variational gamma optimum across the He isoelectronic series (ns={NS}, 1 geminal)")
print("=" * 84)
print(f"{'Z':>2}{'k_opt':>8}{'k/Z':>7}{'gam_opt':>9}{'gam/Z':>8}"
      f"{'E':>14}{'err_mHa':>10}{'d2E/dg2 curv':>14}")
for Z in (1, 2, 3, 4):
    best = None
    surf = {}
    for kf in KFACS:
        k = kf * Z
        es = []
        for g in GAMMAS:
            try:
                E = energy(NS, k, g, Z)
            except Exception:
                E = np.nan
            es.append(E)
            if np.isfinite(E) and (best is None or E < best[0]):
                best = (E, k, g, kf)
        surf[f"kfac={kf}"] = es
    E, k, g, kf = best
    # local curvature in gamma at the optimum (flatness = transferability tolerance)
    gi = GAMMAS.index(g)
    if 0 < gi < len(GAMMAS) - 1:
        gm, gp = GAMMAS[gi - 1], GAMMAS[gi + 1]
        Em, Ep = energy(NS, k, gm, Z), energy(NS, k, gp, Z)
        curv = (Ep - 2 * E + Em) / ((gp - g) * (g - gm))
    else:
        curv = float("nan")
    err = (E - EXACT[Z]) * 1000
    print(f"{Z:>2}{k:>8.2f}{kf:>7.2f}{g:>9.2f}{g/Z:>8.3f}{E:>14.6f}{err:>10.2f}{curv:>14.3f}")
    out["rows"].append(dict(Z=Z, k_opt=k, kfac=kf, gamma_opt=g, gamma_over_Z=g / Z,
                            E=E, err_mHa=err, curvature=float(curv)))
    out["surface"][str(Z)] = surf

gz = [r["gamma_over_Z"] for r in out["rows"]]
gz_noH = [r["gamma_over_Z"] for r in out["rows"] if r["Z"] > 1]
print()
print(f"gamma_opt/Z  all Z : {['%.3f' % v for v in gz]}"
      f"   spread {max(gz)/min(gz):.2f}x")
print(f"gamma_opt/Z  Z>=2  : {['%.3f' % v for v in gz_noH]}"
      f"   spread {max(gz_noH)/min(gz_noH):.2f}x")
out["gamma_over_Z_spread_allZ"] = max(gz) / min(gz)
out["gamma_over_Z_spread_Zge2"] = max(gz_noH) / min(gz_noH)

# --- transferability: fix ONE scaled width, apply to every Z, price the loss ---
print()
print("TRANSFER TEST -- one fixed scaled width c = gamma/Z applied to every Z:")
print(f"{'c':>6} | " + "".join(f"{'Z=%d err' % Z:>12}" for Z in (1, 2, 3, 4)) + f"{'max err':>10}")
out["transfer"] = []
for c in (0.10, 0.15, 0.20, 0.25, 0.30, 0.40):
    errs = []
    for Z in (1, 2, 3, 4):
        kf = [r["kfac"] for r in out["rows"] if r["Z"] == Z][0]
        g = c * Z
        gg = min(GAMMAS, key=lambda x: abs(x - g))     # snap to cached kernel grid
        E = energy(NS, kf * Z, gg, Z)
        errs.append((E - EXACT[Z]) * 1000)
    print(f"{c:>6.2f} | " + "".join(f"{e:>12.2f}" for e in errs)
          + f"{max(errs):>10.2f}")
    out["transfer"].append(dict(c=c, errs_mHa=errs, max_err=max(errs)))

os.makedirs("debug/data", exist_ok=True)
with open("debug/data/r12ci_gamma_isoelectronic.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/r12ci_gamma_isoelectronic.json")
