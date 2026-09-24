"""Fire tests for tests/test_lih_r12ci_neumann_exact.py (the exact prolate-Neumann operator guard).

Each plant breaks the subject in the specific way a guard claims to exclude; the guard must go RED.
Recorded here so the checks are reproducible, per the Sec. 9 rule that a guard whose rejected
answer cannot be named is not a guard.  Runner: debug/qa/fire_test.py (exit 0 = every plant fired).

  A  dead switch in hVee.neumann_potential (exact=True falls through to the legacy sum)
       -> 5 zeta/8 guards, (aa|bb) guard, default-switch guard must fire
  B  dead switch in triangle.coul_mode_potential
       -> the m=1 solid-harmonic guard and the default-switch guard must fire
  C  P and Q interchanged in the operator assembly (ordered integral built backwards)
       -> 5 zeta/8, (aa|bb), m=1 guards must fire
  D  sub-rule too coarse (n_gauss 60 -> 3: the P-side polynomial is no longer integrated exactly)
       -> 5 zeta/8, (aa|bb), m=1 guards must fire
  E  default switch flipped to True
       -> the default-switch guard must fire; the legacy-is-wrong guards pass an explicit exact=False
          and correctly stay green (observed: 1 failed, 2 passed)

Result 2026-09-23 (debug/data/firetest_lih_r12ci_neumann_exact.log): 5/5 plants FIRED.

Run from root:  python debug/firetest_lih_r12ci_neumann_exact.py > debug/data/firetest_lih_r12ci_neumann_exact.log 2>&1
"""
import subprocess
import sys

TEST = "tests/test_lih_r12ci_neumann_exact.py"
RUNNER = [sys.executable, "debug/qa/fire_test.py", TEST]
PLANTS = [
    ("A dead switch (hVee)", "geovac/lih_r12ci/hVee.py",
     "    if exact:\n        op = _KG.exact_neumann(LMAX, 0)=>    if False:\n        op = _KG.exact_neumann(LMAX, 0)",
     "test_exact_operator_reproduces_5zeta_over_8 or test_two_centre_aabb or test_default_switch"),
    ("B dead switch (triangle)", "geovac/lih_r12ci/triangle.py",
     "    if exact:\n        op = _KG.exact_neumann(LMAX, m)=>    if False:\n        op = _KG.exact_neumann(LMAX, m)",
     "test_general_m_mode_potential or test_default_switch"),
    ("C P/Q interchanged", "geovac/lih_r12ci/neumann_exact.py",
     "Qnode[:, :, i][:, :, None] * IP + Pnode[:, :, i][:, :, None] * IQ=>Qnode[:, :, i][:, :, None] * IQ + Pnode[:, :, i][:, :, None] * IP",
     "test_exact_operator_reproduces_5zeta_over_8 or test_two_centre_aabb or test_general_m_mode_potential"),
    ("D coarse sub-rule", "geovac/lih_r12ci/neumann_exact.py",
     "n_gauss: int = 60=>n_gauss: int = 3",
     "test_exact_operator_reproduces_5zeta_over_8 or test_two_centre_aabb or test_general_m_mode_potential"),
    ("E default flipped", "geovac/lih_r12ci/kernels.py",
     "USE_EXACT_NEUMANN = False=>USE_EXACT_NEUMANN = True",
     "test_default_switch or test_legacy_operator_is_wrong"),
]


def main() -> int:
    bad = 0
    for name, target, plant, sel in PLANTS:
        # the dead-switch plants span two lines: hVee/triangle keep `op` defined so the exact branch below
        # (`if exact: radial = op.radial(...)`) still runs -> fire_test needs the SECOND `if exact:` dead too
        plants = [plant]
        if name.startswith("A"):
            plants.append("        if exact:\n            radial = op.radial(g_l, l, 0)=>        if False:\n            radial = op.radial(g_l, l, 0)")
        if name.startswith("B"):
            plants.append("        if exact:\n            radial = op.radial(gB, l, m)=>        if False:\n            radial = op.radial(gB, l, m)")
        cmd = RUNNER + ["--plant-in", target] + sum((["--plant", p] for p in plants), []) + ["-k", sel]
        print(f"\n=== {name} ===", flush=True)
        rc = subprocess.call(cmd)
        print(f"--- {name}: {'FIRED' if rc == 0 else 'DID NOT FIRE (rc %d)' % rc}", flush=True)
        bad += (rc != 0)
    print(f"\n{len(PLANTS) - bad}/{len(PLANTS)} plants fired")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
