"""Test Track K integration: Level 4 with spectral angular solver.

Compares H2 D_e results between FD and spectral angular methods.
"""
import sys
sys.path.insert(0, '.')

import time
from geovac.level4_multichannel import solve_level4_h2_multichannel

R = 1.4  # H2 equilibrium

print("=== FD angular (baseline) ===")
t0 = time.perf_counter()
res_fd = solve_level4_h2_multichannel(
    R, l_max=2, n_alpha=200, angular_method='fd', verbose=True,
)
t_fd = time.perf_counter() - t0

print(f"\n=== Spectral angular (Track K) ===")
t0 = time.perf_counter()
res_sp = solve_level4_h2_multichannel(
    R, l_max=2, angular_method='spectral', n_basis_angular=10,
    verbose=True,
)
t_sp = time.perf_counter() - t0

print(f"\n=== Comparison ===")
print(f"FD:       E_total={res_fd['E_total']:.6f}  D_e%={res_fd['D_e_pct']:.1f}%  time={t_fd:.2f}s")
print(f"Spectral: E_total={res_sp['E_total']:.6f}  D_e%={res_sp['D_e_pct']:.1f}%  time={t_sp:.2f}s")
print(f"Speedup:  {t_fd/t_sp:.1f}x")
print(f"E_total delta: {abs(res_fd['E_total'] - res_sp['E_total']):.2e} Ha")
print(f"D_e delta:     {abs(res_fd['D_e'] - res_sp['D_e']):.2e} Ha")
