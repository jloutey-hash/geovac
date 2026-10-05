"""CI coefficients relaxed alongside the geminal -- analytic ansatz (A) (Paper 19, LiH marriage).

Backs geovac/lih_r12ci/ci_relaxed.py + transition_tensor.py: the (N+1)-dimensional generalized
eigenproblem in {chi_I} U {(F - Fbar) psi} on the truncated M=16 pi reference, e^{-1.5 r12} geminal.

Fast: the non-symmetric transition tensor reproduces the direct 2x2-minor expansion of
<chi_I| v1(1) v2(2) v3(3) v4(4) |psi> for four DIFFERENT non-symmetric one-body matrices (so the
bra/ket order and the electron slots are both tested), and its overlap trace is d_I.

Slow (subprocess, fresh process, threads pinned to 1; ~6-12 min):
  * G1: the frozen 2-dimensional restriction of the new matrices reproduces the banked 2x2
    (golden.json pi_marriage_exp1.5_1e-2: -8.048140609752963) to < 1e-9 Ha, and the new-code
    contractions collapse onto the frozen pieces (aa = E0_grid, d.t = h, d.s = 0).
  * the nine-det relaxed energy is pinned (-8.0504035977, agrees with the step-2 VMC linear-method
    value -8.05047 +- 0.00024), is variational (> exact -8.0705) and lies below the frozen value.
Rejected wrong answers: a bra/ket-swapped or sign-flipped transition term (fast test); a dropped or
mis-weighted antisymmetric kinetic current term (pinned energy: G1 cannot see it -- it vanishes on the
reference); a mis-weighted Fbar subtraction in the kinetic column (G1 identity).
"""
from __future__ import annotations

import os
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
E_GOLDEN_FROZEN_9 = -8.048140609752963      # golden.json pi_marriage_exp1.5_1e-2 E_R12
E_RELAXED_9 = -8.050403597683575            # this module, nine-det (pinned 2026-10-05)
E_EXACT_LIH = -8.0705                       # nonrelativistic exact LiH (variational floor)


# --------------------------------------------------------------------------- fast
def _expected(ia, ib, d, apairs, v, I):
    """0.25 sum_J d_J [alpha 2x2-minor product expansion] [beta ...] with v1..v4 on electrons 1..4."""
    def pairterm(a, a2, va, vb):
        p, q = apairs[a]; p2, q2 = apairs[a2]
        return (va[p, p2] * vb[q, q2] - va[p, q2] * vb[q, p2]
                - va[q, p2] * vb[p, q2] + va[q, q2] * vb[p, p2])
    return 0.25 * sum(d[J] * pairterm(ia[I], ia[J], v[0], v[1]) * pairterm(ib[I], ib[J], v[2], v[3])
                      for J in range(len(d)))


def test_transition_tensor_matches_minor_expansion() -> None:
    from geovac.lih_r12ci.transition_tensor import (determinants, transition_stack, bra_det_ket_psi,
                                                    split_owners)
    Mk = 5
    rng = np.random.default_rng(3)
    n_ap = len(list(combinations(range(Mk), 2)))
    ct = rng.standard_normal(n_ap * n_ap)
    ct[rng.random(ct.size) < 0.8] = 0.0                    # a sparse "truncated" vector
    ct /= np.linalg.norm(ct)
    ia, ib, d, apairs = determinants(ct, Mk)
    assert len(d) >= 10
    st = transition_stack(*bra_det_ket_psi(ia, ib, d), apairs, Mk)
    parts = split_owners(st, len(d))
    v = [rng.standard_normal((Mk, Mk)) for _ in range(4)]   # non-symmetric, one per electron
    eye = np.eye(Mk)
    for I, (iA, iB, iC, iD, val) in enumerate(parts):
        vf = [x.reshape(-1) for x in v]                     # ordered pair x*Mk+y -> v[x, y]
        got = float(np.sum(val * vf[0][iA] * vf[1][iB] * vf[2][iC] * vf[3][iD]))
        assert abs(got - _expected(ia, ib, d, apairs, v, I)) < 1e-12
        e = eye.reshape(-1)
        assert abs(float(np.sum(val * e[iA] * e[iB] * e[iC] * e[iD])) - d[I]) < 1e-13


# --------------------------------------------------------------------------- slow
@pytest.fixture(scope="module")
def nine_det(tmp_path_factory) -> dict:
    out = tmp_path_factory.mktemp("cirelax") / "r.npz"
    env = dict(os.environ)
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[k] = "1"
    run = subprocess.run([sys.executable, str(ROOT / "debug" / "ci_relaxed_run.py"), "1e-2", "--save",
                          str(out)], cwd=ROOT, env=env, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    z = np.load(out, allow_pickle=True)
    import json
    return dict(E_A=float(z['E_A']), E_frozen=float(z['E_frozen']), ident=json.loads(str(z['identity'])))


@pytest.mark.slow
def test_frozen_restriction_reproduces_banked_2x2(nine_det: dict) -> None:
    assert abs(nine_det['E_frozen'] - E_GOLDEN_FROZEN_9) < 1e-9
    idt = nine_det['ident']
    for k in ('aa_vs_E0grid', 'h_vs_pieces', 'hT_vs', 'hVne_vs', 'hVee_vs', 'Fbar_vs', 's_dot_d'):
        assert abs(idt[k]) < 1e-10, (k, idt[k])


@pytest.mark.slow
def test_relaxed_nine_det_energy_pinned(nine_det: dict) -> None:
    assert abs(nine_det['E_A'] - E_RELAXED_9) < 2e-8


@pytest.mark.slow
def test_relaxed_is_variational_and_below_frozen(nine_det: dict) -> None:
    assert E_EXACT_LIH < nine_det['E_A'] < nine_det['E_frozen'] - 1e-3
