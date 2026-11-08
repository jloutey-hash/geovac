"""Fire test for tests/test_paper59_c2b_feynman.py: each planted wrong answer must make its test fail."""
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent / "tests"))
import test_paper59_c2b_feynman as T  # noqa: E402

DEFAULT = dict(NU_NU1=T.NU_NU1, Z_SIGN=T.Z_SIGN, TWIST=T.TWIST, KERNEL=T.KERNEL, PREF=T.PREF)
PLANTS = [
    ("real-momentum argument", dict(Z_SIGN=1), T.test_twisted_legendre_is_the_ec_period_operator),
    ("real-momentum argument, closed form", dict(Z_SIGN=1), T.test_closed_form_is_the_legendre_function),
    ("degree off by 1/200", dict(NU_NU1=-sp.Rational(1, 4) * sp.Rational(201, 200)), T.test_twisted_legendre_is_the_ec_period_operator),
    ("no twist", dict(TWIST=False), T.test_twisted_legendre_is_the_ec_period_operator),
    ("exp for cosh", dict(KERNEL="exp"), T.test_cosh_integral_solves_it_and_the_exp_candidate_does_not),
    ("prefactor pi^2/2", dict(PREF=2), T.test_closed_form_is_the_legendre_function),
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
