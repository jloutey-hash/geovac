"""Fire tests for tests/test_graph_s3_nonconvergence.py (graph-Laplacian non-convergence + SO(4) Casimir).

Each plant breaks the subject in the specific way a guard claims to exclude; the guard must go RED.
Runner: debug/qa/fire_test.py (exit 0 = every plant fired).

  A  cross-l edge planted in the production lattice (geovac/lattice.py): the T+ radial edge
     (n,l,m) -> (n+1,l,m) is redirected to (n+1,l+1,m), joining l-blocks.
       -> (a) the block/kernel/component test and (b) the heat-kernel test must fire,
          for the binary AND topological_weights lattices.
  B  wrong Runge-Lenz reduced element, denominator: (2l+1)(2l+3) -> (2l+1)(2l+2).
       -> (c) the Casimir test must fire.
  C  wrong Runge-Lenz reduced element, radial factor: (n^2 - (l+1)^2) -> (n^2 - l^2).
       -> (c) the Casimir test must fire.

Run from root:  python debug/firetest_graph_s3.py > debug/data/firetest_graph_s3.log 2>&1
"""
import subprocess
import sys

TEST = "tests/test_graph_s3_nonconvergence.py"
RUNNER = [sys.executable, "debug/qa/fire_test.py", TEST]
PLANTS = [
    ("A cross-l edge (binary+topological) -> (a)", "geovac/lattice.py",
     ["                target = (n + 1, l, m)=>                target = (n + 1, l + 1, m)"],
     "test_graph_splits_into_nmax_l_blocks"),
    ("A cross-l edge (binary+topological) -> (b)", "geovac/lattice.py",
     ["                target = (n + 1, l, m)=>                target = (n + 1, l + 1, m)"],
     "test_heat_kernel_vanishes_between_l_blocks"),
    ("B wrong Runge-Lenz denominator -> (c)", TEST,
     ["((2 * l + 1) * (2 * l + 3))=>((2 * l + 1) * (2 * l + 2))"],
     "test_so4_casimir_is_n2_minus_1_on_every_shell"),
    ("C wrong Runge-Lenz radial factor -> (c)", TEST,
     ["np.sqrt((n * n - (l + 1) ** 2)=>np.sqrt((n * n - l ** 2)"],
     "test_so4_casimir_is_n2_minus_1_on_every_shell"),
]


def main() -> int:
    bad = 0
    for name, target, plants, sel in PLANTS:
        cmd = RUNNER + ["--plant-in", target] + sum((["--plant", p] for p in plants), []) + ["-k", sel]
        print(f"\n=== {name} ===", flush=True)
        rc = subprocess.call(cmd)
        print(f"--- {name}: {'FIRED' if rc == 0 else 'DID NOT FIRE (rc %d)' % rc}", flush=True)
        bad += (rc != 0)
    print(f"\n{len(PLANTS) - bad}/{len(PLANTS)} plants fired")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
