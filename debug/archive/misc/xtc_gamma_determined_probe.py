"""DIAGNOSTIC: can the geminal width be DETERMINED (physics-only, no exact energy) so the
CHEAP native-TC operator is accurate -- i.e. does Avery's "solve-and-tabulate" instinct
beat the gamma-fragility?  (Josh's question, 2026-08-24.)

Framing.  "Determine-and-tabulate gamma, don't tune it" is not new -- it is literally F12
quantum chemistry: fix the geminal exponent to a TRANSFERABLE constant (gamma~1.0-1.3) and
never tune it.  But F12 uses the geminal VARIATIONALLY (R12-CI, CI coeffs adapt), where a
fixed gamma is fine -- at the cost of the overlap conditioning wall (kappa(S)~200, v5.0.7).
This probe asks the complementary, cheap-side question: does a fixed/physics-determined
gamma work for the CHEAP non-Hermitian DRESSING (the sparsity-preserving route)?

Physics-only (no-exact-energy) gamma candidates tested per atom:
  gF12 = 1.0   (F12 standard transferable exponent)
  gk   = k     (correlation length = orbital scale)
  gk2  = k/2
For each we report E_TC error.  Then a fine local sweep pins the *oracle* crossing gamma
and its SLOPE dE/dgamma (the fragility): how sharply E passes through exact.

PRE-REGISTERED PREDICTION (written before running):
  - The fixed/physics gammas give INCONSISTENT accuracy across He/Li/Be (esp. Be many x10
    mHa), so a tabulated cheap-dressing gamma FAILS.
  - The oracle crossing is non-transferable (varies per atom, not a clean scale multiple)
    and SHARPENS with electron count (|dE/dgamma| grows He<Li<Be).
  => cheap-side STOP; determine-and-tabulate works only on the variational (conditioning-
     limited) F12 side, confirming the two-walls tension.
"""
from __future__ import annotations
import json
import numpy as np
from geovac import transcorrelated_sturmian as TC

EXACT = {"He": -2.903724, "Li": -7.478060, "Be": -14.667356}
ATOMS = [("He", 2, 2, 1.70, False, 0.83),
         ("Li", 3, 3, 1.30, True, 1.02),
         ("Be", 4, 4, 1.50, True, 0.60)]   # last = coarse crossing from the ns-converge run
NS, NG, NX = 6, 500, 96


def E_TC(Z, ne, k, g, L3):
    s = TC.build_atomic_system(NS, k, g, Z=Z, n_elec=ne, Ng=NG, nx=NX, with_L3=L3)
    op = s.xtc_full() if s.v2 is not None else s.tc2()
    E, im = TC.ground(TC.build_fci_matrix(s, op), hermitian=False)
    return float(E)


def main():
    print("PREDICTION: fixed/physics gamma -> inconsistent (Be worst); oracle crossing "
          "non-transferable + sharpens with N_elec.\n", flush=True)
    out = []
    print(f"{'atom':>4} {'k':>4} | {'gF12=1.0':>9} {'g=k':>9} {'g=k/2':>9} | "
          f"{'oracle_g':>8} {'best_mHa':>8} {'|dE/dg|':>8}", flush=True)
    print("-" * 74, flush=True)
    for label, Z, ne, k, L3, g0 in ATOMS:
        ex = EXACT[label]
        dF12 = (E_TC(Z, ne, k, 1.0, L3) - ex) * 1e3
        dgk = (E_TC(Z, ne, k, k, L3) - ex) * 1e3
        dgk2 = (E_TC(Z, ne, k, k / 2, L3) - ex) * 1e3
        # fine local sweep to pin the crossing + slope
        gs = np.round(np.arange(g0 - 0.15, g0 + 0.16, 0.05), 3)
        ds = np.array([(E_TC(Z, ne, k, float(g), L3) - ex) * 1e3 for g in gs])
        ib = int(np.argmin(np.abs(ds)))
        g_oracle, d_best = float(gs[ib]), float(ds[ib])
        # local slope dE/dgamma (mHa per unit gamma) by finite difference near crossing
        slope = float(abs(np.polyfit(gs, ds, 1)[0]))
        print(f"{label:>4} {k:>4.2f} | {dF12:>+9.2f} {dgk:>+9.2f} {dgk2:>+9.2f} | "
              f"{g_oracle:>8.2f} {d_best:>+8.2f} {slope:>8.1f}", flush=True)
        out.append(dict(atom=label, Z=Z, n_elec=ne, k=k, exact=ex,
                        dF12_mHa=dF12, dgk_mHa=dgk, dgk2_mHa=dgk2,
                        oracle_gamma=g_oracle, best_mHa=d_best, slope_mHa_per_gamma=slope,
                        fine_sweep=[[float(g), float(d)] for g, d in zip(gs, ds)]))
    with open("debug/data/xtc_gamma_determined.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nread: dF12/dgk/dgk2 = error at physics-DETERMINED gamma (no exact-energy peek);",
          flush=True)
    print("      |dE/dg| in mHa per unit gamma = how sharply E crosses exact (the fragility).",
          flush=True)
    print("wrote debug/data/xtc_gamma_determined.json", flush=True)


if __name__ == "__main__":
    main()
