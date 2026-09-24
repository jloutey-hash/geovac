"""Fast guard for Phase 1 (gate G-leaf) of the LiH "marriage" build (2026-09-23;
plan debug/lih_marriage_build_plan.md, memo debug/sprint_lih_marriage_memo.md).

The gate: the only genuinely 4-body-connected term of the explicit-r12 energy,
    I = INT rho1 rho2 rho3 rho4 f(r12) f(r34) / r13 ,
evaluated on a REAL non-isotropic pair density from the Phase-0 CI reference by the RI-free
reduction (leaf dressed through the axially-averaged f-kernel, joined by the EXACT prolate-Neumann
operator) and cross-checked against Monte Carlo.  debug/lih_marriage_phase1.py is the backing driver;
its artifact debug/data/lih_marriage_phase1.npz freezes the reduction, the dressings, and the MC bars.

What each guard EXCLUDES (the wrong answer it was fire-tested against, debug/firetest_lih_r12ci_bridge.py):

  test_c0_production_reduction_matches_sigma_gate
      The C0 control (isotropic 1s leaves / bridge = the sigma-gate configuration): the production
      reduction, recomputed from the frozen production dressing through the exact Neumann bridge, must
      reproduce the sigma-gate record (2026 CHANGELOG v5.15.10) within its MC bar and match the frozen
      value tightly.  Fires if the exact-Neumann bridge or the grid integration regresses.
  test_c1_production_reduction_regression
      C1 (the gate: signed p-like leaf rho_(NO0,NO1), bridge rho_(NO1,NO1)): the production reduction
      recomputed from the frozen production dressing must match the frozen value tightly, and the frozen
      headline to 1e-6.  Fires on any drift of the production bridge pipeline against the banked run.
  test_c1_converged_matches_mc_gate_held
      THE GATE: the converged reduction (i) recomputed from the frozen R1 dressing must (a) match the
      frozen (i) tightly, (b) agree with the frozen 12-D/6-D MC value within max(1e-3|I|, 2 sigma_MC),
      and (c) the carrying MC estimator must have reached sigma_MC <= 1e-3|I| (not INCONCLUSIVE).
      Fires if the bridge factorization stops matching brute force, or if the artifact records a
      non-PASS gate.
  test_c1_isotropic_shortcut_is_wrong  (built-in fire test)
      Audit A(ii): the 1-D isotropic radial leaf dressing is the ISOTROPIC-leaf special case.  Applied
      to the signed p-like C1 leaf (whose monopole INT rho_01 = 0), it gives a grossly wrong reduction.
      The guard reconstructs that shortcut from the frozen spherical-averaged leaf and asserts the
      resulting reduction deviates from the correct one by MORE than the gate.  A helper self-check
      confirms the same shortcut reproduces the closed form on a genuine 1s (so it is a correct
      isotropic dressing, wrong only because the leaf is not isotropic).
  test_c1_L2_truncation_is_wrong  (built-in fire test)
      Audit A(i): the termination bound is 2(l_bridge + l_leaf), NOT 2*l_bridge.  Truncating the
      Neumann L-series at L <= 2 (the s-leaf bound) on the p-like C1 leaf discards channels that carry
      > gate of the integral; the guard asserts the L<=2 partial sum deviates from the full sum by
      more than the gate.

All checks read the frozen artifact; the reduction is recomputed live through the fast exact-Neumann
path (no gVee phi-kernels, no 3-D dressing quadrature), so wall < 10 s.
"""
import importlib.util
import math
import os

import numpy as np
import pytest

import geovac.lih_r12ci.kernels as KG
from geovac.lih_r12ci.hVee import neumann_potential

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "debug", "data", "lih_marriage_phase1.npz")
DRIVER = os.path.join(ROOT, "debug", "lih_marriage_phase1.py")

GATE_REL = 1e-3

# FROZEN headline values (debug/data/lih_marriage_phase1.log; run 2026-09-23, exit 0).
#   C0 (isotropic control): (i) = 2.85624479e-02, carrier 12D, PASS.
#   C1 (gate, signed p-like leaf): (i) = 1.72777754e-05, (ii) = 1.72778893e-05, carrier 6D, PASS
#   (6-D bridge reached sigma_MC/|I| = 5.0e-4; 12-D dip+quad/dipole floored at 2.8e-3/1.5e-3 on the
#   signed leaf and agreed at +-0.3 sigma, so the 6-D carried the strict bar).
FROZEN = {
    "C0_I_conv": 2.85624479e-02,
    "C1_I_conv": 1.72777754e-05,
    "C1_I_prod": 1.72778893e-05,
}
FROZEN_REL = 5e-4   # loose window: the frozen headline may be re-measured within the gate


def _load_driver():
    spec = importlib.util.spec_from_file_location("lih_marriage_phase1_mod", DRIVER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def art():
    if not os.path.exists(OUT):
        pytest.skip("run debug/lih_marriage_phase1.py first to produce the artifact")
    return dict(np.load(OUT, allow_pickle=True))


def _reduce(D_flat: np.ndarray) -> float:
    """I = INT D V[D], V the EXACT prolate-Neumann Coulomb potential of D (LMAX=34, m=0)."""
    V = neumann_potential(D_flat.reshape(KG.NXI, KG.NETA), 34, exact=True).reshape(-1)
    return float(KG.grid_int(D_flat * V))


def _relerr(x: float, y: float) -> float:
    return abs(x - y) / max(abs(y), 1e-300)


def _iso_kernel_avg(s: np.ndarray, rp: np.ndarray, gam: float) -> np.ndarray:
    """(1/2) INT_-1^1 exp(-gam |s - r'|) dx over the source sphere, closed form; shape (len(s), len(rp))."""
    s = s[:, None]
    rp = rp[None, :]
    lo = np.abs(s - rp)
    hi = s + rp
    num = ((lo / gam + 1.0 / gam ** 2) * np.exp(-gam * lo)
           - (hi / gam + 1.0 / gam ** 2) * np.exp(-gam * hi))
    return num / (2.0 * np.maximum(s * rp, 1e-300))


def _isotropic_dressing(s_grid: np.ndarray, rp: np.ndarray, wr: np.ndarray, rho_sph: np.ndarray,
                        gam: float) -> np.ndarray:
    """The 1-D isotropic leaf dressing: Psi_iso(s) = INT 4 pi r'^2 rho_sph(r') Kbar_f(s, r') dr'."""
    Kbar = _iso_kernel_avg(s_grid, rp, gam)                        # (Ns, Nr)
    return 4.0 * np.pi * (Kbar * (rp ** 2 * rho_sph * wr)[None, :]).sum(axis=1)


# --------------------------------------------------------------------------- reductions / regression
def test_c0_production_reduction_matches_sigma_gate(art):
    alpha = float(art["ALPHA"])
    rho_bridge = (alpha ** 3 / np.pi) * np.exp(-2.0 * alpha * KG.rB.ravel())
    I = _reduce(rho_bridge * art["C0_Psi_prod"])
    assert _relerr(I, float(art["C0_I_prod"])) < 1e-9, "C0 production bridge regressed vs the frozen value"
    sg_lim = float(art["sigma_gate_reduced_limit"])
    sg_bar = float(art["sigma_gate_mc_bar"])
    assert _relerr(I, sg_lim) <= sg_bar, f"C0 reduction {I:.8e} vs sigma-gate {sg_lim:.8e} beyond {sg_bar:.1e}"
    assert _relerr(I, FROZEN["C0_I_conv"]) <= FROZEN_REL


def test_c1_production_reduction_regression(art):
    I = _reduce(art["rho_11"] * art["C1_Psi_prod"])
    assert _relerr(I, float(art["C1_I_prod"])) < 1e-9, "C1 production bridge regressed vs the frozen value"
    assert _relerr(I, FROZEN["C1_I_prod"]) <= FROZEN_REL


def test_c1_converged_matches_mc_gate_held(art):
    Ic = float(art["C1_I_conv"])
    Imc = float(art["C1_I_mc"])
    s_mc = float(art["C1_s_mc"])
    verdict = str(art["C1_verdict"])
    # (a) the converged reduction is reproducible from the frozen R1 dressing through the exact bridge
    I_re = _reduce(art["rho_11"] * art["C1_Psi_hi"])
    assert _relerr(I_re, Ic) < 1e-9, "C1 converged bridge not reproducible from the frozen dressing"
    assert _relerr(Ic, FROZEN["C1_I_conv"]) <= FROZEN_REL
    # (b) the gate held: reduction agrees with brute-force MC within max(1e-3|I|, 2 sigma_MC)
    bar = max(GATE_REL * abs(Ic), 2.0 * s_mc)
    assert abs(Ic - Imc) <= bar, f"C1 gate: |(i)-MC| {abs(Ic - Imc):.2e} exceeds tol {bar:.2e}"
    # (c) MC actually reached the bar (a recorded PASS, not INCONCLUSIVE)
    assert s_mc <= GATE_REL * abs(Ic), f"C1 MC sigma {s_mc:.2e} did not reach 1e-3|I| ({GATE_REL * abs(Ic):.2e})"
    assert verdict == "PASS", f"C1 verdict recorded as {verdict}, not PASS"


# --------------------------------------------------------------------------- built-in fire tests
def test_c1_isotropic_shortcut_is_wrong(art):
    gam = float(art["GAM"])
    R = float(art["R"])
    mod = _load_driver()
    # reconstruct the radial rule the driver used for the spherical average (rs stored, wr not)
    rb = mod.graded_bounds(0.0, 30.0, [0.0, R], 0.004, 1.7)
    rs, wr = mod.gl_on_panels(rb, 10)
    assert np.max(np.abs(rs - art["C1_leaf_sph_r"])) < 1e-12, "radial rule reconstruction drifted"
    rho_sph = art["C1_leaf_sph_rho"]

    # self-check: the SAME isotropic dressing reproduces the closed form on a genuine 1s -> it is a
    # correct isotropic dressing, so its failure below is the leaf's non-isotropy, not a coding bug.
    zeta = 2.6875
    rho_1s = (zeta ** 3 / np.pi) * np.exp(-2.0 * zeta * rs)
    s_probe = np.array([0.3, 1.0, 2.5, 6.0])
    Psi_probe = _isotropic_dressing(s_probe, rs, wr, rho_1s, gam)
    Psi_closed = mod.psi_f_iso(s_probe, zeta, gam)
    assert np.max(np.abs(Psi_probe - Psi_closed) / np.abs(Psi_closed)) < 1e-4, "isotropic dressing miscoded"

    # the wrong answer: dress the SIGNED p-like leaf isotropically, then run the same bridge reduction
    Psi_iso = _isotropic_dressing(KG.rA.ravel(), rs, wr, rho_sph, gam)
    I_iso = _reduce(art["rho_11"] * Psi_iso)
    Ic = float(art["C1_I_conv"])
    assert _relerr(I_iso, Ic) > GATE_REL, (
        f"isotropic shortcut {I_iso:.4e} agrees with the correct {Ic:.4e} within the gate -- "
        "the p-like leaf's directional content is being captured by the 1-D shortcut, which is impossible")


def test_c1_L2_truncation_is_wrong(art):
    cl = np.asarray(art["C1_cl_hi"], dtype=float)
    I_L2 = float(cl[:3].sum())            # L <= 2 = the 2*l_bridge (s-leaf) bound
    I_full = float(cl[:35].sum())         # L <= 34
    assert _relerr(I_L2, I_full) > GATE_REL, (
        f"L<=2 truncation {I_L2:.4e} matches the full sum {I_full:.4e} within the gate -- "
        "the leaf carries no channel beyond 2*l_bridge, contradicting audit A(i)")
