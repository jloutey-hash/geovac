"""Paper 60 secular-matrix 1-norm: two-axis (R_MAX x N_GRID) convergence study.

Monkeypatches the module grid (r, dr, r2 must move together) and reports every
norm leg.  Writes debug/data/p60_twoaxis.json.
"""
import json, sys, time
import numpy as np
import geovac.sturmian_secular as S


def set_grid(rmax: float, npts: int) -> None:
    S.R_MAX = float(rmax)
    S.N_GRID = int(npts)
    S.r = np.linspace(1e-7, S.R_MAX, S.N_GRID)
    S.dr = S.r[1] - S.r[0]
    S.r2 = S.r * S.r
    S._GAUNT_CACHE.clear()
    S.reset_caches()


def metrics(cts, Z=2.0):
    cfgs = S.build_configs(cts)
    Tp = S.build_Tprime(cfgs)
    Rnu = np.array([c.Rnu for c in cfgs])
    M = Tp.copy()
    M[np.diag_indices_from(M)] += Z * Rnu
    dTp = np.abs(np.diag(Tp)).sum()
    tot = float(np.abs(Tp).sum())
    return dict(
        K=len(cfgs),
        M_total=float(np.abs(M).sum()),
        M_diag=float(np.abs(np.diag(M)).sum()),
        M_off=float(tot - dTp),
        Tp_full=tot,
        Tp_diag=float(dTp),
        T0=float(Z * Rnu.sum()),
        E=float(-np.sort(np.linalg.eigvalsh(M))[-1] ** 2 / 2),
    )


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    cts = S.gen_configs(3, {l: nmax for l in range(4)})
    out = []
    for rmax in (60, 120, 240, 480, 960, 1920):
        for npts in (9000, 18000, 36000, 72000):
            set_grid(rmax, npts)
            t0 = time.time()
            m = metrics(cts)
            m.update(rmax=rmax, npts=npts, dr=S.dr, wall=time.time() - t0, nmax=nmax)
            out.append(m)
            print(f"nmax={nmax} K={m['K']:4d} rmax={rmax:5d} npts={npts:6d} dr={S.dr:.5f} "
                  f"tot={m['M_total']:.5f} Mdiag={m['M_diag']:.5f} Moff={m['M_off']:.5f} "
                  f"Tpfull={m['Tp_full']:.5f} Tpdiag={m['Tp_diag']:.5f} T0={m['T0']:.5f} "
                  f"E={m['E']:.6f} [{m['wall']:.0f}s]")
            sys.stdout.flush()
    with open(f"debug/data/p60_twoaxis_n{nmax}.json", "w") as fh:
        json.dump(out, fh, indent=1)
