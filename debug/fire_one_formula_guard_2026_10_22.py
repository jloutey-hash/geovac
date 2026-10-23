"""Fire test for tests/test_paper59_one_formula.py: each planted wrong answer must make its test fail."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tests"))
import test_paper59_one_formula as T  # noqa: E402

DEFAULT = dict(Z_SIGN=T.Z_SIGN, PREF=T.PREF, HERON=T.HERON, RIGHT=T.RIGHT)
PLANTS = [
    ("prefactor pi^2/2", dict(PREF=2), T.test_zero_momentum_is_the_three_centre_fibre),
    ("real-momentum sign flipped", dict(Z_SIGN=-1), T.test_real_momentum_cosine_form),
    ("imaginary momentum taken as real", dict(Z_SIGN=-1), T.test_imaginary_momentum_is_a_triangle),
    ("Heron constant 8", dict(HERON=8), T.test_imaginary_momentum_is_a_triangle),
    ("right angle gamma^2 = t^2 - 1", dict(RIGHT=-1), T.test_right_angle_is_the_cm_fibre),
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
