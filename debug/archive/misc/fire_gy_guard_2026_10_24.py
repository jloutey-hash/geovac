"""Fire test for tests/test_paper59_gy_gaussian.py (record 2026-10-23_gaussian_geminal_triangle.md): each planted wrong
answer must make its test fail.  Runs debug/qa/fire_test.py once per plant; exit 0 only if all fire."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
T = "tests/test_paper59_gy_gaussian.py"
D = "tests/data/paper59_gy_2026_10_24.json"
PLANTS = [
    ("a wrong operator coefficient (data file)",
     ["--plant-in", D, "--plant", "-421312132546560*u**9=>-421312132546561*u**9", "-k", "test_operator_annihilates"]),
    ("a wrong right factor F1 (data file)",
     ["--plant-in", D, "--plant", "(-636)*x^3=>(-637)*x^3", "-k", "test_two_pieces"]),
    ("a wrong left factor F2 (data file)",
     ["--plant-in", D, "--plant", "(-829830000)*x^18=>(-829830001)*x^18", "-k", "test_two_pieces"]),
    ("a wrong growth rate, exp(t^2/6) for exp(t^2/8)",
     ["--plant", "(C - sp.Rational(1, 4)) ** 2)=>(C - sp.Rational(1, 3)) ** 2)", "-k", "test_irregular"]),
    ("x = 4 claimed ordinary",
     ["--plant", "fuchs1, roots1 = _indicial_roots(L, sp.Integer(1))=>fuchs1, roots1 = _indicial_roots(L, sp.Integer(4))",
      "-k", "test_x_equals_one"]),
]
fired = 0
for name, args in PLANTS:
    r = subprocess.run([sys.executable, "debug/qa/fire_test.py", T] + args, cwd=ROOT, capture_output=True, text=True)
    ok = r.returncode == 0
    fired += ok
    print(("fired:" if ok else "SILENT (bad):"), name, flush=True)
    if not ok:
        print(r.stdout[-800:], r.stderr[-400:])
print(f"{fired}/{len(PLANTS)} fired")
sys.exit(0 if fired == len(PLANTS) else 1)
