"""Quick timing + monkeypatch sanity for the Paper-60 secular grid study."""
import time, sys, math
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
    K = len(cfgs)
    off = np.abs(Tp).sum() - np.abs(np.diag(Tp)).sum()
    return dict(
        K=K,
        M_total=float(np.abs(M).sum()),
        M_diag=float(np.abs(np.diag(M)).sum()),
        M_off=float(off),
        Tp_full=float(np.abs(Tp).sum()),
        Tp_diag=float(np.abs(np.diag(Tp)).sum()),
        Tp_off=float(off),
        T0=float(Z * Rnu.sum()),
        E=float(-np.sort(np.linalg.eigvalsh(M))[-1] ** 2 / 2),
    )


if __name__ == "__main__":
    for (rmax, npts) in [(60, 18000), (240, 18000)]:
        for nmax in (7, 8):
            set_grid(rmax, npts)
            cts = S.gen_configs(3, {l: nmax for l in range(4)})
            t0 = time.time()
            m = metrics(cts)
            print(f"rmax={rmax} npts={npts} nmax={nmax} K={m['K']} "
                  f"t={time.time()-t0:.1f}s total={m['M_total']:.4f} "
                  f"diag={m['M_diag']:.4f} off={m['M_off']:.4f} "
                  f"Tpfull={m['Tp_full']:.4f} T0={m['T0']:.4f} E={m['E']:.5f}")
            sys.stdout.flush()
