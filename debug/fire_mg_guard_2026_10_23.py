"""Fire test for tests/test_paper59_mg_mechanism.py (record 2026-10-23_gamma1_mechanism.md): each planted wrong
answer must make its test fail.  Runs debug/qa/fire_test.py three times; exit 0 only if all three fire."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
T = "tests/test_paper59_mg_mechanism.py"
PLANTS = [
    ("a wrong map denominator at gamma = 3 (data file)",
     ["--plant-in", "tests/data/paper59_mg_2026_10_23.json",
      "--plant", "(x - 16)**4*(x - 4)**4*(x + 8)=>(x - 16)**4*(x - 4)**4*(x + 9)", "-k", "test_elliptic_piece"]),
    ("the real-momentum sign in the Legendre source",
     ["--plant", "def _source(g, z_sign=-1):=>def _source(g, z_sign=1):", "-k", "test_elliptic_piece"]),
    ("a wrong pair-mass set",
     ["--plant", "-4 * gam**2, -((1 + gam) ** 2)=>-9 * gam**2, -((1 + gam) ** 2)", "-k", "test_kite_singular"]),
]
fired = 0
for name, args in PLANTS:
    r = subprocess.run([sys.executable, "debug/qa/fire_test.py", T] + args, cwd=ROOT, capture_output=True, text=True)
    ok = r.returncode == 0
    fired += ok
    print(("fired:" if ok else "SILENT (bad):"), name, flush=True)
print(f"{fired}/{len(PLANTS)} fired")
sys.exit(0 if fired == len(PLANTS) else 1)
