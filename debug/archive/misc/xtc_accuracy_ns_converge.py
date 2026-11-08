"""Separate basis-incompleteness from the cusp in the per-atom xTC accuracy test.

The ns=3 s-only run showed He clean (+25 mHa s-limit) but Li/Be huge (+227/+555 mHa) --
the small s-only basis cannot represent two shells (1s AND 2s), so that error is RADIAL
INCOMPLETENESS, not the electron-electron cusp.  This converges ns per atom: E_plain(ns)
should flatten at the s-limit (residual = angular/cusp), and only THEN is the xTC
improvement a clean cusp measurement.

Fixed physical-ish common k per atom (no per-point k-opt); unbuffered.
"""
from __future__ import annotations
import json, sys
import numpy as np
from geovac import transcorrelated_sturmian as TC

EXACT = {"He": -2.903724, "Li": -7.478060, "Be": -14.667356}
ATOMS = [("He", 2, 2, 1.70, False), ("Li", 3, 3, 1.30, True), ("Be", 4, 4, 1.50, True)]
NS_LIST = [3, 4, 5, 6]
NG, NX = 500, 96
GAMMAS = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]


def build(ns, k, g, Z, ne, L3):
    return TC.build_atomic_system(ns, k, g, Z=Z, n_elec=ne, Ng=NG, nx=NX, with_L3=L3)


def main():
    out = []
    for label, Z, ne, k, L3 in ATOMS:
        ex = EXACT[label]
        print(f"\n=== {label} (Z={Z}, {ne}e, k={k}) exact={ex:.5f} ===", flush=True)
        print(f"{'ns':>3} {'E_plain':>11} {'dplain(mHa)':>12}", flush=True)
        plain_by_ns = {}
        for ns in NS_LIST:
            s = build(ns, k, 1.0, Z, ne, L3)
            Ep = TC.plain_energy(s)
            plain_by_ns[ns] = Ep
            print(f"{ns:>3} {Ep:>11.5f} {(Ep-ex)*1e3:>+12.2f}", flush=True)
        # xTC gamma sweep at the largest ns (converged radial basis)
        ns = NS_LIST[-1]
        Ep = plain_by_ns[ns]
        print(f"  -- xTC gamma sweep at ns={ns} (E_plain {(Ep-ex)*1e3:+.2f} mHa) --", flush=True)
        print(f"  {'gamma':>5} {'E_TC':>11} {'dTC(mHa)':>9} {'imag':>9}", flush=True)
        rows = []
        for g in GAMMAS:
            s = build(ns, k, g, Z, ne, L3)
            op = s.xtc_full() if s.v2 is not None else s.tc2()
            E, im = TC.ground(TC.build_fci_matrix(s, op), hermitian=False)
            rows.append((g, E, (E-ex)*1e3, im))
            print(f"  {g:>5.2f} {E:>11.5f} {(E-ex)*1e3:>+9.2f} {im:>9.1e}", flush=True)
        dTC = np.array([r[2] for r in rows])
        ibest = int(np.argmin(np.abs(dTC)))
        crosses = bool(dTC.min() < 0 < dTC.max())
        print(f"  best xTC {dTC[ibest]:+.2f} mHa @ g={rows[ibest][0]:.2f}; "
              f"crosses exact (gamma-fragile)? {crosses}", flush=True)
        out.append(dict(atom=label, Z=Z, n_elec=ne, k=k, exact=ex,
                        plain_by_ns={str(n): plain_by_ns[n] for n in NS_LIST},
                        dplain_conv_mHa=(Ep-ex)*1e3,
                        best_xtc_mHa=float(dTC[ibest]), best_gamma=rows[ibest][0],
                        crosses_exact=crosses,
                        rows=[dict(gamma=g, E=E, dmHa=d, imag=im) for g,E,d,im in rows]))
    with open("debug/data/xtc_accuracy_ns_converge.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nwrote debug/data/xtc_accuracy_ns_converge.json", flush=True)


if __name__ == "__main__":
    main()
