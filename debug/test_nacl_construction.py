"""Test NaCl LockedShellMolecule construction at nmax=3."""
import time
import sys
sys.path.insert(0, '.')

from geovac.locked_shell import LockedShellMolecule

# BH first (should be fast now)
print("=" * 60)
print("BH (Z=5, Z=1) at nmax=3")
print("=" * 60)
t0 = time.perf_counter()
bh = LockedShellMolecule(
    Z_A=5, Z_B=1, nmax_A=3, nmax_B=3,
    R=2.329, n_electrons=6,
    active_nmax=2,
)
t_bh = time.perf_counter() - t0
E_bh, _ = bh.solve()
print(f"BH total time: {t_bh:.2f}s, E = {E_bh[0]:.6f} Ha")
print(f"  Target: < 5s -> {'PASS' if t_bh < 5 else 'FAIL'}")

# NaH (Z=11, Z=1) — intermediate
print("\n" + "=" * 60)
print("NaH (Z=11, Z=1) at nmax=3")
print("=" * 60)
t0 = time.perf_counter()
nah = LockedShellMolecule(
    Z_A=11, Z_B=1, nmax_A=3, nmax_B=3,
    R=3.566, n_electrons=12,
    active_nmax=3,
)
t_nah = time.perf_counter() - t0
E_nah, _ = nah.solve()
print(f"NaH total time: {t_nah:.2f}s, E = {E_nah[0]:.6f} Ha")

# NaCl (Z=11, Z=17)
print("\n" + "=" * 60)
print("NaCl (Z=11, Z=17) at nmax=3")
print("=" * 60)
t0 = time.perf_counter()
nacl = LockedShellMolecule(
    Z_A=11, Z_B=17, nmax_A=3, nmax_B=3,
    R=4.461, n_electrons=28,
    active_nmax=3,
)
t_nacl = time.perf_counter() - t0
E_nacl, _ = nacl.solve()
print(f"NaCl total time: {t_nacl:.2f}s, E = {E_nacl[0]:.6f} Ha")
print(f"  Target: completes -> PASS")

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"BH:   {t_bh:.2f}s, E = {E_bh[0]:.6f} Ha")
print(f"NaH:  {t_nah:.2f}s, E = {E_nah[0]:.6f} Ha")
print(f"NaCl: {t_nacl:.2f}s, E = {E_nacl[0]:.6f} Ha")
