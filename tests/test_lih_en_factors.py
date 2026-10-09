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
Option b (PI 2026-10-09): + amplitude-switch columns A1 psi = sum_i a_i psi, A2 psi = sum_{i<j} a_i a_j psi,
a = phi_core / phi_core(R_Li) (signed; core_switch.AmpSwitch): fast real-space checks of the switch
(normalisation, a^2 = w, closed-form derivatives vs FD, node handling) and of the 'a' token algebra;
slow nine-det G1 (amp columns dropped -> density-only phase 1 to < 1e-9) + pinned combined energy.
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
E_PHB_9 = -8.057432973804504              # A + W1 + W2 + A1 + A2 (option b, dens + amp), nine-det (2026-10-09)
G = 1.5
_a, _b, _z0 = 0.7, 1.3, -0.4                # synthetic axial, non-spherical switch
_c, _e, _z1 = 0.35, 0.55, 0.25              # a second synthetic switch (amplitude token 'a')


# --------------------------------------------------------------------------- pointwise evaluator
def _tokfun(tok: str, r: np.ndarray) -> np.ndarray:
    x, y, z = r[..., 0], r[..., 1], r[..., 2]
    rho = np.sqrt(x * x + y * y); eip = (x + 1j * y) / rho
    w = np.exp(-_a * rho ** 2 - _b * (z - _z0) ** 2)
    drw = -2 * _a * rho * w; dzw = -2 * _b * (z - _z0) * w
    lapw = w * (4 * _a ** 2 * rho ** 2 - 4 * _a + 4 * _b ** 2 * (z - _z0) ** 2 - 2 * _b)
    # second, DIFFERENT synthetic switch for the option-b amplitude token 'a'
    am = np.exp(-_c * rho ** 2 - _e * (z - _z1) ** 2)
    lapa = am * (4 * _c ** 2 * rho ** 2 - 4 * _c + 4 * _e ** 2 * (z - _z1) ** 2 - 2 * _e)
    return dict(w=w, lapw=lapw, a=am, lapa=lapa, u=rho * drw + z * dzw, dzw=dzw, z=z, drwP=drw * eip,
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


def test_amp_token_operators_and_laplacians(pts) -> None:
    """Option b: W1('a'), W2('a') (amplitude token) evaluate to sum a_i, sum_{i<j} a_i a_j; lap_op maps
    'a' -> 'lapa' (NOT 'lapw'); en_ops('both') = [W1, W2, A1, A2]."""
    am = np.stack([_tokfun('a', pts[:, e]) for e in range(4)], axis=1).real
    A1 = am.sum(1); A2 = sum(am[:, i] * am[:, j] for i in range(4) for j in range(i + 1, 4))
    assert np.allclose(_eval(A.W1_op('a'), pts), A1, rtol=1e-13, atol=0)
    assert np.allclose(_eval(A.W2_op('a'), pts), A2, rtol=1e-13, atol=0)
    for op in (A.W1_op('a'), A.W2_op('a')):
        lap = _fd_lap(op, pts)
        assert np.max(np.abs(_eval(A.lap_op(op), pts) - lap)) < 2e-5 * np.max(np.abs(lap))
    # mixed same-electron product (W1_den * W1_amp carries w_i a_i) -- the 'a*w' variant
    mix = A.op_mul(A.W1_op('w'), A.W1_op('a'))
    w = np.stack([_tokfun('w', pts[:, e]) for e in range(4)], axis=1).real
    assert np.allclose(_eval(mix, pts), w.sum(1) * A1, rtol=1e-13, atol=0)
    assert any(t[e] == ('a', 'w') for t in mix for e in range(4))


def _amp_switches():
    sys.path.insert(0, str(ROOT / "debug"))
    import en_factors_stage0_lm as E0
    from geovac.lih_r12ci import pi_m16 as PM
    from geovac.lih_r12ci.core_switch import CoreSwitch, AmpSwitch
    T_no, occ, m_no, R = E0.switch_reference()
    prims = PM.build_primitives(PM.CFG, PM.PRIM_CACHE_PI, verbose=False)['prims']
    return (AmpSwitch.from_reference(prims, T_no, occ, m_no, R),
            CoreSwitch.from_reference(prims, T_no, occ, m_no, R), R)


def _fd_vgl(f, r, h=1e-4):
    v = f(r); g = np.zeros(r.shape); lap = np.zeros(len(r))
    for k in range(3):
        rp = r.copy(); rp[:, k] += h; rm = r.copy(); rm[:, k] -= h
        vp, vm = f(rp), f(rm)
        g[:, k] = (vp - vm) / (2 * h); lap += (vp - 2 * v + vm) / h ** 2
    return v, g, lap


def test_amp_switch_realspace_and_node_handling() -> None:
    """The real amplitude switch a = phi_c/phi_c(Li) of the M=16 pi reference: a(Li) = 1, a^2 == the
    density switch w, closed-form grad / lap == finite differences in the core AND on a line crossing
    the outer nodal surface, where the SIGNED a changes sign (an |a| implementation -- kinked, with a
    surface delta in its Laplacian -- fails the straddle assert), and the node is >= 2.5 bohr from Li
    with |a| < 1e-3 wherever a < 0 (so signed == |phi_c|/|phi_c(Li)| on the core region)."""
    amp, den, R = _amp_switches()
    li = np.array([[0.0, 0.0, -R / 2]])
    assert abs(amp.value(li)[0] - 1.0) < 1e-12
    rng = np.random.default_rng(5)
    core = rng.normal(scale=0.8, size=(30, 3)) + li
    assert np.allclose(amp.value(core) ** 2, den.value(core), rtol=1e-10, atol=1e-14)
    # points straddling the node surface along the bond axis (sign change near z ~ 1.34 bohr)
    zz = np.linspace(1.1, 1.6, 11)
    ring = np.stack([0.05 + 0 * zz, 0 * zz, zz], axis=-1)
    a_ring = amp.value(ring)
    assert a_ring.min() < 0 < a_ring.max(), "test points must straddle the node"
    for r in (core, ring):
        v, g, l = amp.vgl(r)
        fv, fg, fl = _fd_vgl(amp.value, r)
        assert np.max(np.abs(fg - g)) < 1e-6 * max(np.max(np.abs(g)), 1e-3)
        assert np.max(np.abs(fl - l)) < 1e-5 * max(np.max(np.abs(l)), 1e-3)
    dnode = amp.node_distance()
    assert 2.5 < dnode < np.inf
    rho = np.linspace(0, 6, 121); zg = np.linspace(-R / 2 - 6, R / 2 + 6, 241)
    Rg, Zg = np.meshgrid(rho, zg, indexing='ij')
    pts = np.stack([Rg, 0 * Rg, Zg], -1).reshape(-1, 3)
    av = amp.value(pts)
    assert np.max(np.abs(av[av < 0])) < 1e-3


# --------------------------------------------------------------------------- slow
def _run(tmp_path_factory, phase2: bool, switch: str = "dens") -> dict:
    out = tmp_path_factory.mktemp("en") / "r.npz"
    env = dict(os.environ)
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[k] = "1"
    cmd = [sys.executable, str(ROOT / "debug" / "en_factors_run.py"), "1e-2", "--save", str(out),
           "--switch", switch]
    if phase2:
        cmd.append("--phase2")
    run = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    z = np.load(out, allow_pickle=True)
    r = dict(E=json.loads(str(z['E'])), ident=json.loads(str(z['identity'])), sw=json.loads(str(z['switch'])))
    if phase2:
        r['E2'] = json.loads(str(np.load(str(out)[:-4] + "_phase2.npz", allow_pickle=True)['E']))
    return r


@pytest.fixture(scope="module")
def nine_det(tmp_path_factory) -> dict:
    """phase 1 only (~10 min)."""
    return _run(tmp_path_factory, False)


@pytest.fixture(scope="module")
def nine_det_both(tmp_path_factory) -> dict:
    """option b: phase 1 with the density AND amplitude switch columns (~20 min)."""
    return _run(tmp_path_factory, False, "both")


@pytest.fixture(scope="module")
def nine_det_p2(tmp_path_factory) -> dict:
    """phase 1 + phase 2 (~70 min; the batched triangle class dominates)."""
    return _run(tmp_path_factory, True)


@pytest.mark.slow
def test_G1_identity_and_internal_consistency(nine_det: dict) -> None:
    assert abs(nine_det['E']['A'] - E_RELAXED_9) < 1e-9
    idt = nine_det['ident']
    for k in ('G1_A_vs_ci_relaxed', 'S_P_vs_dS_I', 'LA_ref_W1', 'LA_ref_W2', 'HPW_vs_kinsym_W1',
              'HPW_vs_kinsym_W2'):
        assert abs(idt[k]) < 1e-10, (k, idt[k])
    assert nine_det['sw']['core_no'] == 0 and nine_det['sw']['kind'] == 'dens'


@pytest.mark.slow
def test_phase1_nine_det_energy_pinned(nine_det: dict) -> None:
    assert abs(nine_det['E']['A_W12'] - E_PH1_9) < 2e-8
    assert E_EXACT_LIH < nine_det['E']['A_W12'] < E_RELAXED_9 - 1e-3


@pytest.mark.slow
def test_phase2_nine_det_energy_pinned(nine_det_p2: dict) -> None:
    assert abs(nine_det_p2['E2']['A'] - E_RELAXED_9) < 1e-9
    assert abs(nine_det_p2['E2']['A_W12_G'] - E_PH2_9) < 2e-8
    assert E_EXACT_LIH < nine_det_p2['E2']['A_W12_G'] <= nine_det_p2['E']['A_W12']


@pytest.mark.slow
def test_option_b_nine_det_G1_and_pin(nine_det_both: dict) -> None:
    """Option b (PI 2026-10-09): [W1, W2] (density switch) + [A1, A2] (amplitude switch a = phi_c/phi_c(Li)).
    G1: dropping the amp columns reproduces the density-only phase 1 (and dropping all W, ci_relaxed)
    to < 1e-9; the amp columns pass the same internal-consistency identities; the combined energy is
    pinned, variational, and below both single-switch sub-ansatzes."""
    E = nine_det_both['E']
    assert abs(E['A'] - E_RELAXED_9) < 1e-9
    assert abs(E['A_W12'] - E_PH1_9) < 1e-9
    for k, v in nine_det_both['ident'].items():
        assert abs(v) < 1e-10, (k, v)
    assert nine_det_both['sw']['kind'] == 'both' and nine_det_both['sw']['node_dist_Li'] > 2.5
    assert abs(E['A_W12_A12'] - E_PHB_9) < 2e-8
    assert E_EXACT_LIH < E['A_W12_A12'] < min(E['A_W12'], E['A_A12']) - 1e-3
