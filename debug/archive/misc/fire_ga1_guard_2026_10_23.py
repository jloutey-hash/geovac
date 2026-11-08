"""Fire test for tests/test_paper59_ga1_geminal.py: each planted wrong answer must make its test fail."""
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent / "tests"))
import test_paper59_ga1_geminal as T  # noqa: E402

DEFAULT = dict(GAMMA=T.GAMMA, Z_SIGN=T.Z_SIGN, FACTOR_KEY=T.FACTOR_KEY)
PLANTS = [
    ("gamma = 2 Legendre operator", dict(GAMMA=sp.Integer(2)), T.test_factor_is_the_image_of_the_legendre_operator),
    ("real-momentum argument", dict(Z_SIGN=1), T.test_factor_is_the_image_of_the_legendre_operator),
    ("the exponential case's C2b as the factor", dict(FACTOR_KEY="c2b"), T.test_factor_is_the_image_of_the_legendre_operator),
]
fired = 0
for name, patch, test in PLANTS:
    for k, v in DEFAULT.items():
        setattr(T, k, v)
    for k, v in patch.items():
        setattr(T, k, v)
    try:
        test()
        print("SILENT (bad):", name, flush=True)
    except AssertionError:
        fired += 1
        print("fired:", name, flush=True)
for k, v in DEFAULT.items():
    setattr(T, k, v)
print(f"{fired}/{len(PLANTS)} fired")
