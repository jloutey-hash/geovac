"""Paper 60 sec:molecular -- the Shibuya-Wulfman (SW) metric on a FULL sigma basis (l < n, m = 0),
versus the s-only sine-basis family every earlier Paper 60 conditioning measurement used.

Built 2026-10-09 after reading Avery & Avery, "Hyperspherical Harmonics and Their Physical
Applications" (World Scientific 2017): the Shibuya-Wulfman-orbitals section (fixed-S posing) and
the two-centre density expansion (ill-conditioned by overcompleteness at small kR).

The builder lives in geovac/sturmian_sigma_basis.py (position-space prolate quadrature,
independent of the momentum-space sine-basis code).  Anchors recorded here:
  * s-only Bessel deficit reproduces tests/test_paper60_one_direction.py's 0.380/0.696/0.907;
  * the fixed-S kernel posing agreed to every printed digit with an independent momentum-space
    4D-harmonic implementation (2026-10-09 scratchpad cross-check).

Output: debug/data/p60_sigma_basis_sw.json.  Diagnostic driver; backing tests in
tests/test_paper60_sigma_basis.py.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geovac.sturmian_sigma_basis import (  # noqa: E402
    NE, NX, basis, bessel_deficit, build, galerkin_energy, gerade_ungerade, kernel_energy)

H2P_EXACT = {1.0: -1.451786313, 2.0: -1.102634214, 4.0: -0.796084884, 8.0: -0.627570389}


def metric_stats(N: int, s: float, sonly: bool) -> dict:
    S, _, labels = build(N, 1.0, s, sonly)
    ev = np.linalg.eigvalsh(S)
    Sg, Su = gerade_ungerade(S, labels)
    return {
        "fns_per_centre": len(labels),
        "lam_min": float(ev[0]), "lam_max": float(ev[-1]), "cond": float(ev[-1] / ev[0]),
        "cond_gerade": float(np.linalg.cond(Sg)), "cond_ungerade": float(np.linalg.cond(Su)),
        "n_eig_below_1e-2": int((ev < 1e-2).sum()), "n_eig_below_1e-4": int((ev < 1e-4).sum()),
        "bessel_deficit_1sB": bessel_deficit(N, s, sonly),
    }


def main() -> None:
    t0 = time.time()
    out = {"meta": {"source": "debug/sturmian_sigma_basis_sw.py + geovac/sturmian_sigma_basis.py",
                    "date": "2026-10-09",
                    "quadrature": f"prolate Gauss-Laguerre {NX} x Gauss-Legendre {NE}",
                    "precision_floor": ("float64 entries are accurate to ~2e-12 absolute; any cond above "
                                        "~1e13 here is at the float floor and unreliable (e.g. kR=1, N>=5 "
                                        "gives negative cond; true 2.8e19 / 5.7e25 at N=5/6). Exact values: "
                                        "debug/p60_sigma_precision_crosscheck.py (60-digit moment route).")}}
    S, _, _ = build(3, 1.3, 2.0)
    out["checks"] = {"max_intra_dev_from_identity": float(np.max(np.abs(S[:6, :6] - np.eye(6))))}
    out["h2plus"] = []
    for N in (2, 3, 4, 5):
        for R, Eex in H2P_EXACT.items():
            out["h2plus"].append({"N": N, "fns_per_centre": len(basis(N)), "R": R, "exact": Eex,
                                  "galerkin_err": galerkin_energy(N, R) - Eex,
                                  "kernel_err": kernel_energy(N, R) - Eex})
    out["metric"] = []
    for s in (1.0, 1.5, 2.0, 3.0, 4.0, 6.0):
        for sonly, Ns in ((True, (4, 8, 12)), (False, (2, 3, 4, 5, 6))):
            for N in Ns:
                rec = {"kR": s, "family": "s-only" if sonly else "sigma", "N": N}
                rec.update(metric_stats(N, s, sonly))
                out["metric"].append(rec)
    out["meta"]["wall_s"] = time.time() - t0
    path = Path(__file__).parent / "data" / "p60_sigma_basis_sw.json"
    path.write_text(json.dumps(out, indent=1))
    for r in out["h2plus"]:
        print(f"N={r['N']} ({r['fns_per_centre']:2d}/centre) R={r['R']}: "
              f"Galerkin {r['galerkin_err']:+.2e}  kernel {r['kernel_err']:+.2e}")
    for r in out["metric"]:
        print(f"kR={r['kR']} {r['family']:6s} N={r['N']:2d}: cond={r['cond']:.3e} "
              f"g={r['cond_gerade']:.3e} u={r['cond_ungerade']:.3e} soft(<1e-2)={r['n_eig_below_1e-2']} "
              f"deficit={r['bessel_deficit_1sB']:.4f}")
    print(f"checks {out['checks']}  wall {out['meta']['wall_s']:.0f}s -> {path}")


if __name__ == "__main__":
    main()
