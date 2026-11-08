"""Fire test for tests/test_paper59_geminal_shapes.py: each planted wrong answer must make its test fail."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tests"))
import test_paper59_geminal_shapes as T  # noqa: E402

DEFAULT = dict(LG_INNER=T.LG_INNER, COMPANION=T.COMPANION, ALT_MOMENTUM=T.ALT_MOMENTUM)
PLANTS = [
    ("L_G with 4 theta + 3", dict(LG_INNER=3), T.test_gaussian_avatar_is_annihilated),
    ("H_G offered as the growing companion", dict(COMPANION="K0"), T.test_gaussian_companion_has_gaussian_growth),
    ("the piece's own family as an alternative slot", dict(ALT_MOMENTUM=("gamma",)), T.test_momentum_slot_is_the_geminals),
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
