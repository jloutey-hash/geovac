"""Why does adding H 2pz collapse the NaH well?

B2 (base + diffuse-s + H 2pz, frozen core) gave D_e = 0.599 eV at R_eq = 5.0,
against B1's 1.339 eV at 3.62. Adding a basis function cannot legitimately do
that. Two explanations survive:

  (a) FROZEN CORE. The n_core=5 approximation was validated only at M=7
      (0.21% of D_e). At M=8 it already drifts 0.13 eV. It may break outright
      once p functions are present on H.
  (b) BSSE. Borrowing grows with basis size and is longest-ranged at large R,
      which would both inflate and outwardly shift the well.

Linear dependence is already excluded -- s_min = 1.3e-2 at B2, three orders
above the ill-conditioning threshold.

DISCRIMINATING TEST. Build B2' = base + H 2pz with NO diffuse function: M = 8,
so 12-electron all-electron FCI is still tractable (dim 1820) and can be
compared directly against the frozen-core result on the SAME basis.

  * if all-electron B2' is sane and frozen-core B2' collapses  -> (a)
  * if both collapse                                           -> (b), or real
  * if both are sane -> the collapse needs the diffuse+p COMBINATION

Counterpoise is reported alongside so (b) is measured, not assumed.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.noci_engine import HARTREE_TO_EV

from noci_nah_basis_extension import (      # noqa: E402
    D_E_EXP_EV, R_GRID, Z_H, Z_NA,
    atom_energy, atom_energy_cp, base_spec, build_shapes,
    molecule_energies, well,
)


def curve(spec, shapes, n_core, label):
    e_na = atom_energy(spec, shapes, "Na", 11.0, 11, n_core)
    e_h = atom_energy(spec, shapes, "H", 1.0, 1, 0)
    e_diss = e_na + e_h

    raw, cp, smin = [], [], []
    for r in R_GRID:
        e_mol, _, sm = molecule_energies(spec, shapes, r, n_core, None)
        raw.append(e_mol)
        smin.append(sm)
        cp.append(e_mol
                  - atom_energy_cp(spec, shapes, "Na", r, n_core)
                  - atom_energy_cp(spec, shapes, "H", r, n_core)
                  + e_diss)

    wr, wc = well(R_GRID, raw, e_diss), well(R_GRID, cp, e_diss)
    print(f"  {label:<34} raw {wr['D_e_Ha']*HARTREE_TO_EV:6.3f} eV @ "
          f"{wr['R_eq']:5.2f} | CP {wc['D_e_Ha']*HARTREE_TO_EV:6.3f} eV @ "
          f"{wc['R_eq']:5.2f} | s_min {min(smin):.2e}")
    return wr, wc, raw, cp


def main():
    print("=== B2 diagnostic: why does H 2pz collapse the well? ===")
    print(f"    experiment D_e = {D_E_EXP_EV} eV @ 3.566 a0\n")
    shapes = build_shapes()

    b0 = base_spec()
    b2p = b0 + [("H2pz", "H", "2p", (0, 0, 1), Z_H)]          # M=8, no diffuse
    b1 = b0 + [("Hdiff", "H", "1s", (0, 0, 0), 0.30)]          # M=8, diffuse
    print()

    print("[M=7 control]")
    curve(b0, shapes, 0, "base, ALL-ELECTRON")
    curve(b0, shapes, 5, "base, frozen core")

    print("\n[M=8, +H 2pz only -- the discriminating rung]")
    ae = curve(b2p, shapes, 0, "base+H2pz, ALL-ELECTRON")
    fc = curve(b2p, shapes, 5, "base+H2pz, frozen core")

    print("\n[M=8, +diffuse only -- for comparison]")
    curve(b1, shapes, 0, "base+diffuse, ALL-ELECTRON")
    curve(b1, shapes, 5, "base+diffuse, frozen core")

    print("\n=== verdict ===")
    d = abs(ae[0]["D_e_Ha"] - fc[0]["D_e_Ha"]) * HARTREE_TO_EV
    print(f"  frozen-core vs all-electron on base+H2pz: {d:.3f} eV apart")
    if d > 0.2:
        print("  -> FROZEN CORE is unreliable once p functions sit on H;")
        print("     the ladder's frozen-core magnitudes are not trustworthy.")
    else:
        print("  -> frozen core is NOT the culprit; look to BSSE / real physics.")
    print("  (compare each raw vs CP column above for the BSSE magnitude)")


if __name__ == "__main__":
    main()
