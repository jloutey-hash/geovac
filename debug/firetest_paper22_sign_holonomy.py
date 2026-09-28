"""Fire tests for tests/test_paper22_sign_holonomy.py (Paper 22 sign-problem subsection).

Each plant breaks the subject in the way a guard claims to exclude; the guard must go RED.
Runner: debug/qa/fire_test.py (exit 0 = every plant fired).

  A  spanning-tree gauge fix dropped (BFS assigns s[j] = s[i])      -> frozen frustration must fire
  A2 bipartite gauge dropped ((-1)^(n+m) -> 1)                        -> adjacency gauge test must fire
  B  one V_ee off-diagonal sign flipped (1s^2 <-> 1s2s element)       -> frozen frustration/Delta must fire
  C  the retracted (1/16)A hop Hamiltonian swapped in for physical h1 -> frozen values + builder test must fire
  D  identity "rotation" (Gaunt zeros kept)                          -> native-vs-rotated inequality must fire
  E  random-sign control made vacuous (random signs = native signs)   -> not-better-than-random test must fire

Run from root:  python debug/firetest_paper22_sign_holonomy.py > debug/data/firetest_paper22_sign_holonomy.log 2>&1
"""
import os
import subprocess
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

TEST = "tests/test_paper22_sign_holonomy.py"
RUNNER = [sys.executable, "debug/qa/fire_test.py", TEST]
H_DIAG = "    h = np.diag([-Z * Z / (2.0 * n * n) for (n, l, m) in orbs])"
PLANTS = [
    ("A gauge fix dropped",
     ["                    s[j] = -np.sign(x) * s[i]=>                    s[j] = s[i]"],
     "test_frozen_gs_sector_frustration"),
    ("A2 bipartite gauge dropped",
     ["G = np.diag([(-1.0) ** (n + m) for=>G = np.diag([(1.0) ** (n + m) for"],
     "test_fock_adjacency_gauge_free"),
    ("B one V_ee sign flipped",
     ["    H = 0.5 * (H + H.T)\n    return H=>    H = 0.5 * (H + H.T)\n"
      "    H[0, 1] = -H[0, 1]\n    H[1, 0] = -H[1, 0]\n    return H"],
     "test_frozen_gs_sector_frustration"),
    ("C hop Hamiltonian swapped in",
     [H_DIAG + "=>" + H_DIAG + " + lat.adjacency.toarray() / 16.0"],
     "test_frozen_gs_sector_frustration or test_physical_hamiltonian_is_builder"),
    ("D identity rotation",
     ["    U = _random_orth(len(orbs), np.random.default_rng(1000))=>    U = np.eye(len(orbs))"],
     "test_gaunt_zeros_contain_sign_problem_n3"),
    ("E random-sign control vacuous",
     ["np.abs(w) * rng.choice([-1.0, 1.0], size=len(w))=>w + 0.0 * rng.choice([-1.0, 1.0], size=len(w))"],
     "test_native_signs_not_better_than_random_n3"),
]


def main() -> int:
    bad = 0
    for name, plants, sel in PLANTS:
        cmd = RUNNER + sum((["--plant", p] for p in plants), []) + ["-k", sel]
        print(f"\n=== {name} ===", flush=True)
        rc = subprocess.call(cmd)
        print(f"--- {name}: {'FIRED' if rc == 0 else 'DID NOT FIRE (rc %d)' % rc}", flush=True)
        bad += (rc != 0)
    print(f"\n{len(PLANTS) - bad}/{len(PLANTS)} plants fired")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
