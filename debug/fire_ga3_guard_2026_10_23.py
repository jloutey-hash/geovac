"""Fire test for tests/test_paper59_ga3_geminal.py (record 2026-10-23_geminal_shapes_g2.md): each planted wrong
answer must make its test fail.  Runs debug/qa/fire_test.py four times; exit 0 only if all four fire."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
T = "tests/test_paper59_ga3_geminal.py"
D = "tests/data/paper59_ga3_2026_10_23.json"
PLANTS = [
    ("a wrong map denominator (data file)",
     ["--plant-in", D, "--plant", "(x - 9)**5*(x - 1)**5*(x + 3)=>(x - 9)**5*(x - 1)**5*(x + 4)", "-k", "test_order_two_piece_is_the_image"]),
    ("the real-momentum sign in the Legendre source",
     ["--plant", "Z_SIGN = -1=>Z_SIGN = 1", "-k", "test_order_two_piece_is_the_image"]),
    ("a wrong exponent at t = 4",
     ["--plant", "sp.Rational(-11, 2)=>sp.Rational(-9, 2)", "-k", "test_order_one_pieces"]),
    ("a wrong order-one solution (data file)",
     ["--plant-in", D, "--plant", "12715)/x**2=>12716)/x**2", "-k", "test_order_one_pieces"]),
]
fired = 0
for name, args in PLANTS:
    r = subprocess.run([sys.executable, "debug/qa/fire_test.py", T] + args, cwd=ROOT, capture_output=True, text=True)
    ok = r.returncode == 0
    fired += ok
    print(("fired:" if ok else "SILENT (bad):"), name, flush=True)
    if not ok:
        print(r.stdout[-600:])
print(f"{fired}/{len(PLANTS)} fired")
sys.exit(0 if fired == len(PLANTS) else 1)
