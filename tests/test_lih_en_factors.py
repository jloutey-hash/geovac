"""Electron-nucleus (core-region) extension of the CI-relaxed LiH eigenproblem (Paper 19, LiH marriage).

Backs geovac/lih_r12ci/en_factors.py + en_algebra.py + core_switch.py: the (A) problem
{chi_I} U {(F - Fbar) psi} extended by r12-free core-switched columns W1 psi = sum_i w_i psi,
W2 psi = sum_{i<j} w_i w_j psi (phase 1) and core-weighted geminals Fcc psi, Fc1 psi (phase 2), with the
parameter-free switch w = |phi_core|^2 / |phi_core(R_Li)|^2 (core natural orbital of the reference).

Fast (pure algebra, no grid): the symbolic operators evaluated pointwise in real space for a synthetic
non-spherical axial switch -- W1, W2 and their Laplacians, and the closed-form Laplacian of the
weighted geminals (g_lap: the grad w . grad f cross term decomposed into one-electron tokens incl.
the circulation +-1 tokens) against a finite-difference Laplacian.
Slow (subprocess, threads pinned to 1, ~70 min: phase 2 batched triangles dominate): nine-det G1 identity (dropping the new columns
reproduces ci_relaxed's pinned -8.0504035977 to < 1e-9) + internal consistency (transition kinetic
current vanishes on the reference; transition vs symmetric kinetic forms agree), and the pinned
phase-1 / phase-2 nine-det energies (variational, below (A)).
Rejected wrong answers: a sign / factor slip or a swapped circulation token in the cross term (fast);
a mis-weighted lap w or dropped current term in the transition kinetic (pinned energies); the
new-column embedding disturbing the (A) block (G1).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from geovac.lih_r12ci import en_algebra as A

ROOT = Path(__file__).resolve().parent.parent
E_RELAXED_9 = -8.050403597683575            # ci_relaxed nine-det (A), tests/test_lih_ci_relaxed.py
E_EXACT_LIH = -8.0705
E_PH1_9 = -8.051861170421194              # A + W1 + W2, nine-det (pinned 2026-10-05)
E_PH2_9 = -8.052394144646039              # A + W1 + W2 + Fcc + Fc1, nine-det (pinned 2026-10-05)
G = 1.5
_a, _b, _z0 = 0.7, 1.3, -0.4                # synthetic axial, non-spherical switch


# --------------------------------------------------------------------------- pointwise evaluator
def _tokfun(tok: str, r: np.ndarray) -> np.ndarray:
    x, y, z = r[..., 0], r[..., 1], r[..., 2]
    rho = np.sqrt(x * x + y * y); eip = (x + 1j * y) / rho
    w = np.exp(-_a * rho ** 2 - _b * (z - _z0) ** 2)
    drw = -2 * _a * rho * w; dzw = -2 * _b * (z - _z0) * w
    lapw = w * (4 * _a ** 2 * rho ** 2 - 4 * _a + 4 * _b ** 2 * (z - _z0) ** 2 - 2 * _b)
    return dict(w=w, lapw=lapw, u=rho * drw + z * dzw, dzw=dzw, z=z, drwP=drw * eip,
                drwM=drw * np.conj(eip), rhoP=rho * eip, rhoM=rho * np.conj(eip))[tok]


def _edgefun(key: str, d: np.ndarray) -> np.ndarray:
    return {'f': np.exp(-G * d), 'Y': np.exp(-G * d) / d, 'coul': 1.0 / d}[key]


def _eval(op, r: np.ndarray) -> np.ndarray:
    """op: Op (toks4 -> c) or GOp ((toks4, edges) -> c), at points r (Nw, 4, 3)."""
    out = np.zeros(r.shape[0], complex)
    for key, c in op.items():
        toks, edges = (key, ()) if isinstance(key[0], tuple) and (len(key) == 4) else key
        v = np.full(r.shape[0], c, complex)
        for e in range(4):
            for t in toks[e]:
                v = v * _tokfun(t, r[:, e])
        for (i, j, k) in edges:
            v = v * _edgefun(k, np.linalg.norm(r[:, i - 1] - r[:, j - 1], axis=-1))
        out += v
    return out


def _fd_lap(op, r: np.ndarray, h: float = 1e-4) -> np.ndarray:
    lap = np.zeros(r.shape[0], complex); v0 = _eval(op, r)
    for e in range(4):
        for dd in range(3):
            rp = r.copy(); rp[:, e, dd] += h; rm = r.copy(); rm[:, e, dd] -= h
            lap += (_eval(op, rp) - 2 * v0 + _eval(op, rm)) / h ** 2
    return lap


@pytest.fixture(scope="module")
def pts() -> np.ndarray:
    rng = np.random.default_rng(11)
    return rng.normal(scale=0.8, size=(40, 4, 3))


# --------------------------------------------------------------------------- fast
def test_W_operators_and_laplacians(pts) -> None:
    w = np.stack([_tokfun('w', pts[:, e]) for e in range(4)], axis=1).real
    W1 = w.sum(1); W2 = sum(w[:, i] * w[:, j] for i in range(4) for j in range(i + 1, 4))
    assert np.allclose(_eval(A.W1_op(), pts), W1, rtol=1e-13, atol=0)
    assert np.allclose(_eval(A.W2_op(), pts), W2, rtol=1e-13, atol=0)
    assert np.allclose(_eval(A.op_mul(A.W1_op(), A.W2_op()), pts), W1 * W2, rtol=1e-13, atol=0)
    for op in (A.W1_op(), A.W2_op()):
        lap = _fd_lap(op, pts)
        assert np.max(np.abs(_eval(A.lap_op(op), pts) - lap)) < 2e-5 * np.max(np.abs(lap))


def test_weighted_geminal_laplacian_closed_form(pts) -> None:
    """g_lap (incl. the u / dzw z / circulation-token cross term) == finite-difference Laplacian,
    and the result is real (the e^{+-i phi} tokens pair into cos(phi_k - phi_l))."""
    for kind in ('F', 'cc', 'c1'):
        op = A.g_geminal(kind)
        lap = _fd_lap(op, pts)
        got = _eval(A.g_lap(op, G), pts)
        assert np.max(np.abs(got.imag)) < 1e-12
        assert np.max(np.abs(got - lap)) < 2e-5 * np.max(np.abs(lap)), kind
    lapW = _fd_lap(A.g_from(A.W2_op()), pts)
    assert np.max(np.abs(_eval(A.g_lap(A.g_from(A.W2_op()), G), pts) - lapW)) < 2e-5 * np.max(np.abs(lapW))


def test_class_geminal_span(pts) -> None:
    """{F, Fcc, Fc1} spans the step-3 classes: F_cv = Fc1 - 2 Fcc, F_vv = F - Fc1 + Fcc."""
    w = np.stack([_tokfun('w', pts[:, e]) for e in range(4)], axis=1).real
    F = _eval(A.g_geminal('F'), pts).real; cc = _eval(A.g_geminal('cc'), pts).real
    c1 = _eval(A.g_geminal('c1'), pts).real
    vv = sum((1 - w[:, i - 1]) * (1 - w[:, j - 1]) * np.exp(-G * np.linalg.norm(pts[:, i - 1] - pts[:, j - 1], axis=-1))
             for (i, j) in A.PAIRS6)
    assert np.allclose(F - c1 + cc, vv, rtol=1e-12, atol=1e-14)


# --------------------------------------------------------------------------- slow
@pytest.fixture(scope="module")
def nine_det(tmp_path_factory) -> dict:
    out = tmp_path_factory.mktemp("en") / "r.npz"
    env = dict(os.environ)
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[k] = "1"
    run = subprocess.run([sys.executable, str(ROOT / "debug" / "en_factors_run.py"), "1e-2", "--phase2",
                          "--save", str(out)], cwd=ROOT, env=env, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    z = np.load(out, allow_pickle=True)
    z2 = np.load(str(out)[:-4] + "_phase2.npz", allow_pickle=True)
    return dict(E=json.loads(str(z['E'])), ident=json.loads(str(z['identity'])),
                E2=json.loads(str(z2['E'])), sw=json.loads(str(z['switch'])))


@pytest.mark.slow
def test_G1_identity_and_internal_consistency(nine_det: dict) -> None:
    assert abs(nine_det['E']['A'] - E_RELAXED_9) < 1e-9
    assert abs(nine_det['E2']['A'] - E_RELAXED_9) < 1e-9
    idt = nine_det['ident']
    for k in ('G1_A_vs_ci_relaxed', 'S_P_vs_dS_I', 'LA_ref_W1', 'LA_ref_W2', 'HPW_vs_kinsym_W1',
              'HPW_vs_kinsym_W2'):
        assert abs(idt[k]) < 1e-10, (k, idt[k])
    assert nine_det['sw']['core_no'] == 0 and nine_det['sw']['kind'] == 'dens'


@pytest.mark.slow
def test_nine_det_energies_pinned(nine_det: dict) -> None:
    assert abs(nine_det['E']['A_W12'] - E_PH1_9) < 2e-8
    assert abs(nine_det['E2']['A_W12_G'] - E_PH2_9) < 2e-8


@pytest.mark.slow
def test_variational_and_below_A(nine_det: dict) -> None:
    E1 = nine_det['E']['A_W12']; E2 = nine_det['E2']['A_W12_G']
    assert E_EXACT_LIH < E2 <= E1 < E_RELAXED_9 - 1e-3
