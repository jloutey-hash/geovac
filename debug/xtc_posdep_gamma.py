"""DIAGNOSTIC #1: does a DETERMINED, position-dependent geminal width gamma(r) beat the
global-gamma fragility?  (Josh 2026-08-24; sequel to xtc_gamma_determined_probe.)

Literature target (explorer): the physical correlation length is position-dependent,
gamma(r) proportional to n(r)^{1/3} (Wigner-Seitz / electron-avoidance radius; local
range-separation mu(r)).  We bake that SHAPE in (no per-system fitting of the shape) and
scan only ONE overall scale gamma0.

Determined shape (physics-only, from the reference density -- NO exact-energy peek):
  gamma(rbar) = gamma0 * ( n(rbar) / n(r_typ) )^{1/3},   rbar = (R1+R2)/2,
where n(r) is the s-only reference (occupied-Lowdin-MO) density and r_typ is its
density-weighted mean radius, so gamma(r_typ)=gamma0 (comparable to the global gamma).
The width is TIGHT (large gamma) in the dense core, WIDE (small gamma) in the diffuse
tail -- exactly the shell-varying "diameter" picture.

DECISIVE METRIC (not "is E stationary" -- the non-Herm E_TC is monotone by nature):
  Is the accurate scale gamma0 TRANSFERABLE across He/Li/Be (values close together)
  and is the crossing LESS SHARP (|dE/dgamma0| smaller, esp. Be) than the global-gamma
  oracle?  GO if position-dependence collapses the scatter (0.83/0.97/0.60) toward one
  universal scale AND tames Be's slope (253); STOP if gamma0 stays scattered/sharp.

Reuses the tracked engine's kernels/FCI; only the kernel gamma is made an (R1,R2) array
(build_kernels is elementwise, so it accepts the array directly).
"""
from __future__ import annotations
import json
import numpy as np
from geovac import transcorrelated_sturmian as TC

EXACT = {"He": -2.903724, "Li": -7.478060, "Be": -14.667356}
ATOMS = [("He", 2, 2, 1.70, False), ("Li", 3, 3, 1.30, True), ("Be", 4, 4, 1.50, True)]
NS, NG, NX = 6, 500, 96
GAMMA0 = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
# global-gamma oracle from xtc_gamma_determined.json (for side-by-side):
GLOBAL = {"He": (0.83, 40), "Li": (0.97, 43), "Be": (0.60, 253)}


def occ_density(ns, k, Z, n_elec, r, wr):
    """Reference (occupied-Lowdin-MO) radial density n(r) on the grid, and <r>_n."""
    S, h1s, Rtab, W = TC.build_one_body(ns, r, wr, k, Z)
    X = TC.lowdin(S)                                   # S^{-1/2}, Lowdin MOs
    diag_e = np.diag(TC.transform_1(h1s, X))
    order = np.argsort(diag_e)
    n_occ_spatial = 1 if n_elec == 2 else 2            # 1s^2 (He) or 1s^2 2s.. (Li/Be)
    # radial function of Sturmian mu (mu=0..ns-1) is Rtab[mu+1][0]; MO = sum_mu X[mu,i] R_mu
    Rmat = np.array([Rtab[mu + 1][0] for mu in range(ns)])  # (ns, Ng)
    dens = np.zeros_like(r)
    for i in order[:n_occ_spatial]:
        mo = X[:, i] @ Rmat                             # (Ng,)
        occ = 2.0 if (n_elec != 3 or i != order[1]) else 1.0  # Li 2s singly occ
        dens += occ * mo * mo
    dens = np.maximum(dens, 1e-30)
    W = r * r * wr                                      # 3D radial measure
    r_typ = float(np.sum(W * dens * r) / np.sum(W * dens))
    return dens, r_typ


def build_posdep(ns, k, gamma0, Z, n_elec, with_L3, r, wr, dens, r_typ):
    """build_atomic_system, but with gamma(R1,R2) = gamma0 * (n(rbar)/n(r_typ))^{1/3}."""
    S, h1s, Rtab, W = TC.build_one_body(ns, r, wr, k, Z)
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    rbar = 0.5 * (R1 + R2)
    n_typ = float(np.interp(r_typ, r, dens))
    # determined n^{1/3} shape on the 1D grid, BOUNDED to [0.35,2.5]x (finite tail hole),
    # then NORMALIZED to density-weighted mean 1 so gamma0 IS the effective scale
    s1d = np.clip((dens / n_typ) ** (1.0 / 3.0), 0.35, 2.5)
    Wr = r * r * wr
    s1d = s1d / (np.sum(Wr * dens * s1d) / np.sum(Wr * dens))   # <s1d>_n = 1
    s_bar = np.interp(rbar.ravel(), r, s1d).reshape(rbar.shape)
    G = gamma0 * s_bar                                   # (Ng,Ng) position-dependent gamma
    Km = TC.build_kernels(r, G, nx=NX)                  # elementwise -> accepts the array
    eri_coul, eri_w, eri_K = TC.two_body(ns, Rtab, W, Km)
    V3 = TC.three_body(ns, Rtab, W, Km) if with_L3 else None
    X = TC.lowdin(S)
    h1o = TC.transform_1(h1s, X)
    eri_w_o = TC.transform_2(eri_w, X); eri_K_o = TC.transform_2(eri_K, X)
    V3o = TC.transform_3(V3, X) if with_L3 else None
    nso = 2 * ns
    dets, didx = TC.make_dets(nso, n_elec)
    hso = TC.h_spin(h1o, nso)
    asym_w = TC.asym_from_phys(eri_w_o, nso); asym_K = TC.asym_from_phys(eri_K_o, nso)
    order = np.argsort(np.diag(h1o))
    if n_elec == 2:
        ref = (2 * order[0], 2 * order[0] + 1)
    else:
        ref = (2 * order[0], 2 * order[0] + 1, 2 * order[1])
    # assemble the non-Herm TC operator (D + K [+ xTC-L3]) exactly as sys.xtc_full/tc2
    hso_eff = hso.copy(); asym_eff = asym_w + asym_K; v0 = 0.0
    if with_L3:
        v2, v1, v0 = TC.xtc_contract(V3o, nso, ref)
        asym_eff = asym_eff + v2
        hso_eff = hso_eff + v1
    H = TC.build_H(dets, didx, hso_eff, asym_eff, nso, v0=v0)
    E, im = TC.ground(H, hermitian=False)
    return float(E)


def main():
    print("PREDICTION: determined gamma(r)~n^{1/3} should collapse the accurate scale "
          "toward ONE transferable gamma0 and tame Be's slope (253) if position-dependence "
          "is the fix.\n", flush=True)
    print(f"{'atom':>4} | {'posdep: g0*':>10} {'best_mHa':>8} {'|dE/dg0|':>8} | "
          f"{'global g*':>9} {'gl.slope':>8}", flush=True)
    print("-" * 62, flush=True)
    out = []
    for label, Z, ne, k, L3 in ATOMS:
        ex = EXACT[label]
        r, wr = TC.make_grid(k, Ng=NG)
        dens, r_typ = occ_density(NS, k, Z, ne, r, wr)
        ds = []
        for g0 in GAMMA0:
            E = build_posdep(NS, k, g0, Z, ne, L3, r, wr, dens, r_typ)
            ds.append((E - ex) * 1e3)
        ds = np.array(ds)
        ib = int(np.argmin(np.abs(ds)))
        slope = float(abs(np.polyfit(GAMMA0, ds, 1)[0]))
        g_star, g_slope = GLOBAL[label]
        print(f"{label:>4} | {GAMMA0[ib]:>10.2f} {ds[ib]:>+8.2f} {slope:>8.1f} | "
              f"{g_star:>9.2f} {g_slope:>8.0f}", flush=True)
        out.append(dict(atom=label, r_typ=r_typ, best_gamma0=GAMMA0[ib],
                        best_mHa=float(ds[ib]), posdep_slope=slope,
                        global_gamma=g_star, global_slope=g_slope,
                        rows=[dict(g0=g, dmHa=float(d)) for g, d in zip(GAMMA0, ds)]))
    with open("debug/data/xtc_posdep_gamma.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nread: if the posdep best_gamma0 cluster (transferable) AND posdep_slope << "
          "global_slope (esp. Be), position-dependence dissolves the fragility.", flush=True)
    print("wrote debug/data/xtc_posdep_gamma.json", flush=True)


if __name__ == "__main__":
    main()
