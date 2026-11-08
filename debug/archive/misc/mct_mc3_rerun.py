"""MC3 rerun (disclosed deviation): same three points and forms as mct_mc3.py, at mp.dps = 70 with quadrature
error estimates, to see which side limited the first run at ~2.7e-25."""
import random
import mpmath as mp
from mct_mc3 import k_form, x_form

if __name__ == "__main__":
    mp.mp.dps = 70
    rng = random.Random(20261029)
    worst = mp.mpf(0)
    for i in range(3):
        c1, c2 = [mp.mpf(rng.randint(1, 25)) / 100 for _ in range(2)]
        m1, m2 = [mp.mpf(rng.randint(5, 30)) / 10 for _ in range(2)]
        D = mp.mpf(rng.randint(1, 40)) / 10
        a, b = k_form(c1, c2, m1, m2, D), x_form(c1, c2, m1, m2, D)
        rel = abs(a - b) / abs(a); worst = max(worst, rel)
        print(f"pt{i} k {mp.nstr(a, 40)}\n    x {mp.nstr(b, 40)}\n    rel {mp.nstr(rel, 5)}", flush=True)
    print("worst rel diff", mp.nstr(worst, 5), "PASS(<1e-25)" if worst < mp.mpf(10)**-25 else "FAIL", flush=True)
