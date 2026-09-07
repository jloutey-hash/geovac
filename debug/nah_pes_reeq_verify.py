"""Verify the NaH well-minimum R_eq = 3.736 a0 (Paper 58 Table II) by an
actual PES scan of the 3-determinant fragment-native NoCI curve, reusing the
exact integral/ladder machinery from tests/test_paper58_nah_ladder.py.

Purpose: the R_eq/D_e figures were driver-backed and flagged (a prior driver
serialized a corrupted energy once). This re-runs the scan from the pinned
setup and reports the minimum + a parabolic refinement.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tests"))
import test_paper58_nah_ladder as T
_shapes, _nah_integrals, _ladders = T._shapes, T._nah_integrals, T._ladders
E_NA_REF, E_H_REF = T.E_NA_REF, T.E_H_REF
from geovac import noci_engine as E

HARTREE_EV = 27.211386245988

sh = _shapes()
e_diss = E_NA_REF + E_H_REF
dets3 = _ladders()["cov+ionH (3 dets)"]

Rgrid = np.round(np.arange(2.8, 5.21, 0.15), 3)
curve = []
for R in Rgrid:
    s, h, g, vnn = _nah_integrals(float(R), sh)
    e3, _ = E.noci_ground_gensc(dets3, s, h, g)
    tot = e3 + vnn
    curve.append((float(R), tot))
    print(f"R={R:6.3f}  E3+Vnn={tot:.6f}  bind_eV={(e_diss-tot)*HARTREE_EV:+.4f}")

Rs = np.array([c[0] for c in curve])
Es = np.array([c[1] for c in curve])
imin = int(Es.argmin())
print(f"\ngrid minimum: R={Rs[imin]:.3f}  E={Es[imin]:.6f}")

# parabolic refinement using the min and its two neighbours
if 0 < imin < len(Rs) - 1:
    x0, x1, x2 = Rs[imin-1], Rs[imin], Rs[imin+1]
    y0, y1, y2 = Es[imin-1], Es[imin], Es[imin+1]
    # vertex of parabola through 3 points
    denom = (x0-x1)*(x0-x2)*(x1-x2)
    A = (x2*(y1-y0) + x1*(y0-y2) + x0*(y2-y1)) / denom
    B = (x2*x2*(y0-y1) + x1*x1*(y2-y0) + x0*x0*(y1-y2)) / denom
    R_eq = -B/(2*A)
    C = y1 - A*x1*x1 - B*x1
    E_min = A*R_eq*R_eq + B*R_eq + C
    D_e = (e_diss - E_min) * HARTREE_EV
    print(f"parabolic R_eq = {R_eq:.4f} a0   E_min = {E_min:.6f} Ha")
    print(f"D_e = {D_e:.4f} eV   (vs Table II: R_eq=3.736, D_e=1.071 eV)")
    print(f"R_eq error vs 3.5667 exp = {100*(R_eq-3.5667)/3.5667:+.1f}%")
