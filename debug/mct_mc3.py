"""MC3 of docs/predictions/2026-10-27_master_count_tabulation.md: scaling identity of the slice integral.

k-form   I = int_0^inf dk exp(-D Delta1) / (Delta1 Delta2),  Delta_i = sqrt(c_i k^2 + m_i^2)
scaled   I = m2^{-1} c1^{-1/2} int_1^inf exp(-Dhat x) (x^2-1)^{-1/2} (sigma x^2 + 1 - sigma)^{-1/2} dx,
         Dhat = D m1, sigma = c2 m1^2 / (c1 m2^2).
The two sides use different integration variables (k on [0,inf) vs x on [1,inf) with the endpoint singularity left
to tanh-sinh), so agreement is a check of the substitution, not an identity of the code.
"""
from __future__ import annotations

import random

import mpmath as mp


def k_form(c1, c2, m1, m2, D):
    f = lambda k: mp.exp(-D * mp.sqrt(c1 * k * k + m1 * m1)) / (
        mp.sqrt(c1 * k * k + m1 * m1) * mp.sqrt(c2 * k * k + m2 * m2))
    return mp.quad(f, [0, 1, 4, 16, mp.inf])


def x_form(c1, c2, m1, m2, D):
    dh = D * m1
    sg = c2 * m1**2 / (c1 * m2**2)
    f = lambda x: mp.exp(-dh * x) / mp.sqrt((x * x - 1) * (sg * x * x + 1 - sg))
    return mp.quad(f, [1, 2, 8, mp.inf]) / (m2 * mp.sqrt(c1))


if __name__ == "__main__":
    mp.mp.dps = 45
    rng = random.Random(20261029)
    worst = mp.mpf(0)
    for i in range(3):
        c1, c2 = [mp.mpf(rng.randint(1, 25)) / 100 for _ in range(2)]
        m1, m2 = [mp.mpf(rng.randint(5, 30)) / 10 for _ in range(2)]
        D = mp.mpf(rng.randint(1, 40)) / 10
        a, b = k_form(c1, c2, m1, m2, D), x_form(c1, c2, m1, m2, D)
        rel = abs(a - b) / abs(a)
        worst = max(worst, rel)
        print(f"pt{i} c1={c1} c2={c2} m1={m1} m2={m2} D={D}\n  k-form {mp.nstr(a, 32)}\n  x-form {mp.nstr(b, 32)}"
              f"\n  rel diff {mp.nstr(rel, 5)}", flush=True)
    print("worst rel diff", mp.nstr(worst, 5), "PASS(<1e-25)" if worst < mp.mpf(10)**-25 else "FAIL", flush=True)
