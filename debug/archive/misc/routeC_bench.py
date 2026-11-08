"""Micro-benchmark: mpf multiply-add throughput, and matrix bilinear-form cost."""
import time
import mpmath as mp

mp.mp.dps = 30

N = 4000000
a = mp.mpf('1.234567890123')
b = mp.mpf('9.876543210987')
t0 = time.time()
s = mp.mpf(0)
for _ in range(N):
    s += a * b
dt = time.time() - t0
print(f"plain loop: {N} mults in {dt:.3f}s -> {N/dt:.0f}/s")

t0 = time.time()
s2 = mp.fsum(a * b for _ in range(N))
dt = time.time() - t0
print(f"fsum: {N} mults in {dt:.3f}s -> {N/dt:.0f}/s")
