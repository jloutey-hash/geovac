"""Fire tests for tests/test_lih_r12ci_bridge.py (the Phase 1 / gate G-leaf bridge guard).

Each plant breaks the subject in the specific way a guard claims to exclude; the guard must go RED.
Recorded here so the checks are reproducible, per the Sec. 9 rule that a guard whose rejected answer
cannot be named is not a guard.  Runner: debug/qa/fire_test.py (exit 0 = every plant fired).

Requires the banked artifact debug/data/lih_marriage_phase1.npz (run debug/lih_marriage_phase1.py first);
otherwise the guard's fixture SKIPs and fire_test reports NOT RUN (not a dead guard).

  A  dead switch in hVee.neumann_potential (exact=True falls through to the legacy cumulative-GL sum)
       -> the reduction recomputed by the guard now uses the wrong (kinked-kernel) bridge, so it no
          longer matches the frozen exact-Neumann values: the C0/C1 production and C1 converged
          regressions must fire.
  B  P and Q interchanged in the exact operator (neumann_exact assembly built backwards)
       -> same reduction regressions must fire.
  C  the isotropic-shortcut fire test made vacuous: replace the (wrong) 1-D isotropic dressing of the
       signed p-like leaf with the CORRECT directional dressing (art C1_Psi_hi).  The shortcut then
       "agrees" with the correct reduction, so test_c1_isotropic_shortcut_is_wrong must fire (it asserts
       DISAGREEMENT > gate) -- proving that assertion is not vacuously true.
  D  the L<=2 fire test made vacuous: sum the Neumann series to L=34 instead of L=2, so the "truncated"
       value equals the full one.  test_c1_L2_truncation_is_wrong must fire.

Run from root:  python debug/firetest_lih_r12ci_bridge.py > debug/data/firetest_lih_r12ci_bridge.log 2>&1
"""
import subprocess
import sys

TEST = "tests/test_lih_r12ci_bridge.py"
RUNNER = [sys.executable, "debug/qa/fire_test.py", TEST]
REDUCTION_K = "test_c0_production_reduction or test_c1_production_reduction or test_c1_converged"
PLANTS = [
    ("A dead switch (hVee)", "geovac/lih_r12ci/hVee.py",
     ["    if exact:\n        op = _KG.exact_neumann(LMAX, 0)=>    if False:\n        op = _KG.exact_neumann(LMAX, 0)",
      "        if exact:\n            radial = op.radial(g_l, l, 0)=>        if False:\n            radial = op.radial(g_l, l, 0)"],
     REDUCTION_K),
    ("B P/Q interchanged (neumann_exact)", "geovac/lih_r12ci/neumann_exact.py",
     ["Qnode[:, :, i][:, :, None] * IP + Pnode[:, :, i][:, :, None] * IQ"
      "=>Qnode[:, :, i][:, :, None] * IQ + Pnode[:, :, i][:, :, None] * IP"],
     REDUCTION_K),
    ("C isotropic fire test vacuous", TEST,
     ['    Psi_iso = _isotropic_dressing(KG.rA.ravel(), rs, wr, rho_sph, gam)=>    Psi_iso = art["C1_Psi_hi"]'],
     "test_c1_isotropic_shortcut_is_wrong"),
    ("D L<=2 fire test vacuous", TEST,
     ["    I_L2 = float(cl[:3].sum())            # L <= 2 = the 2*l_bridge (s-leaf) bound"
      "=>    I_L2 = float(cl[:35].sum())            # L <= 2 = the 2*l_bridge (s-leaf) bound"],
     "test_c1_L2_truncation_is_wrong"),
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
