r"""ADDITIVE F12 estimate for LiH (2026-09-22, PI-approved): the pure-orbital ceiling
(-8.03170, v5.15.19 M>16 push) + the r12 CUSP corrections the orbital basis structurally
cannot reach.  Standard F12 methodology (an explicitly-correlated correction added to a
conventional energy).

The cusp correction per electron pair = E(radial-converged) - E(radial + r12), computed with
the validated He-like 2e Hylleraas machinery (debug/lih_r12_ceiling_probe.py `energy`):
  radial-converged pair  = {t2=(r1-r2)^2, s=r1+r2} basis (the in-out correlation orbitals give)
  radial + r12           = {u=r12, t2, s, u2, ut2}   (adds the cusp)
The difference is the pure cusp -- the piece the geminal adds ON TOP of a radial-converged
orbital pair, so it is (approximately) additive to the -8.032 orbital energy with no double count.

Pairs: Li core (Z=3, tight, dominant), H- valence proxy (Z=1, UPPER bound -- the bonded LiH
valence is less correlated than free H-).  Core-valence (2-center, electrons on different
centers) cusp is small and neglected (budget-checked).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lih_r12_ceiling_probe import energy    # noqa: E402

RADIAL = ['t2', 's']
RICH = ['u', 't2', 's', 'u2', 'ut2']

# (name, Z, single-zeta ref, exact non-rel, note)
SYS = [
    ("Li core   (Z=3)", 3.0, 2.6875, -7.27991, "tight, dominant, transferable"),
    ("H- valence(Z=1)", 1.0, 0.6875, -0.52775, "UPPER bound (bonded valence less correlated)"),
    ("He control(Z=2)", 2.0, 27 / 16, -2.90372, "cross-check vs known Hylleraas"),
]

E_ORB = -8.03170          # pure-orbital LiH ceiling (v5.15.19)
E_EXACT = -8.070

print("=" * 82)
print("Per-pair r12 CUSP correction (He-like 2e machinery)")
print("=" * 82)
cusp = {}
for name, Z, zeta, exact, note in SYS:
    E0, _, _ = energy(Z, zeta, ['u'])
    E_rad, _ = energy(Z, zeta, RADIAL)[1], None
    E_rich = energy(Z, zeta, RICH)[1]
    c = E_rad - E_rich
    cusp[name] = c
    print(f"  {name}: ref E0={E0:.5f} (exact {exact})  radial={E_rad:.5f}  rich={E_rich:.5f}")
    print(f"      cusp = radial - rich = {c*1e3:5.1f} mHa   [{note}]")

print("\n" + "=" * 82)
print("ADDITIVE F12 estimate for LiH")
print("=" * 82)
core = cusp["Li core   (Z=3)"]
val = cusp["H- valence(Z=1)"]
gap = E_ORB - E_EXACT
print(f"  pure-orbital ceiling            E = {E_ORB:.5f}   ({(E_EXACT-E_ORB)*1e3:+.1f} mHa from exact)")
print(f"  exact-orbital gap (all cusps)     = {gap*1e3:+.1f} mHa   <- budget the cusps must fill")
print(f"  + core-core cusp (solid)          = {-core*1e3:+.1f} mHa   -> E = {E_ORB-core:.5f}")
print(f"  + valence cusp (H- upper bound)   = {-val*1e3:+.1f} mHa   -> E = {E_ORB-core-val:.5f}")
print(f"\n  >>> LiH additive-F12: core cusp alone -> {E_ORB-core:.4f} "
      f"({(E_EXACT-(E_ORB-core))*1e3:+.1f} mHa, near-chemical)")
print(f"      with valence (upper bnd) -> {E_ORB-core-val:.4f}; true value between these, "
      f"budget caps total cusp at {-gap*1e3:.0f} mHa.")
print(f"  Honest: additive estimate (not variational); core cusp transferable, valence approximate.")
