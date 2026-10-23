"""Fire test for tests/test_paper59_c2b_feynman_g2.py: each planted wrong answer must make its test fail."""
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent / "tests"))
import test_paper59_c2b_feynman_g2 as T  # noqa: E402

PLANTS = [
    ("gamma = 3/2 operator", dict(GAMMA=sp.Rational(3, 2)), T.test_order_two_piece_is_the_image_of_the_legendre_operator),
    ("real-momentum argument", dict(Z_SIGN=1), T.test_order_two_piece_is_the_image_of_the_legendre_operator),
    ("real-momentum argument, j", dict(Z_SIGN=1), T.test_j_invariant_at_gamma_2),
    ("gamma = 5/2 in j", dict(GAMMA=sp.Rational(5, 2)), T.test_j_invariant_at_gamma_2),
]
default = dict(GAMMA=T.GAMMA, Z_SIGN=T.Z_SIGN)
fired = 0
for name, patch, test in PLANTS:
    for k, v in default.items():
        setattr(T, k, v)
    for k, v in patch.items():
        setattr(T, k, v)
    try:
        test()
        print("SILENT (bad):", name, flush=True)
    except AssertionError:
        fired += 1
        print("fired:", name, flush=True)
for k, v in default.items():
    setattr(T, k, v)
print(f"{fired}/{len(PLANTS)} fired")
