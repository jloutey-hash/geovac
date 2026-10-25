"""Fire test for tests/test_paper59_ib_kite.py (record 2026-10-24_ibp_kite.md): each planted wrong answer must make
its test fail.  Runs debug/qa/fire_test.py once per plant; exit 0 only if all fire."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
T = "tests/test_paper59_ib_kite.py"
D = "tests/data/paper59_ib_2026_10_24.json"
PLANTS = [
    ("a wrong coefficient of the PF operator (data file)",
     ["--plant-in", D, "--plant", "142507253410806448147230747593=>142507253410806448147230747594", "-k", "test_pf_operator"]),
    ("a wrong Taylor coefficient of the PF series (data file)",
     ["--plant-in", D, "--plant", "-1.1787426104210813587879484339=>-1.1787426104210813587879484338", "-k", "test_pf_operator"]),
    ("a direct-pairing point claimed among the crossed cuts",
     ["--plant", '"-49/4", "-16"=>"-25/4", "-16"', "-k", "test_pf_singular"]),
    ("a crossed point claimed ordinary",
     ["--plant", 'DIRECT = ["-25/4", "-25"]=>DIRECT = ["-16", "-25"]', "-k", "test_pf_singular"]),
    ("the genuine point z = -1 claimed apparent",
     ["--plant", "red(F(-1)), p) == 8=>red(F(-1)), p) == 11", "-k", "test_genuine"]),
    ("a wrong coefficient of the gamma = 3/2 IBP operator (data file)",
     ["--plant-in", D, "--plant", "1807606699553104784=>1807606699553104785", "-k", "test_ibp_operator"]),
    ("the bubble reduced in D = 4",
     ["--plant", "assert _bubble_check(3) <=>assert _bubble_check(4) <", "-k", "test_reduction_method"]),
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
