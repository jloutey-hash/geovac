r"""VMC-over-FCI for the rigorous variational LiH energy (2026-09-22, PI-directed).

Build plan: debug/lih_rigorous_vmc_build_plan.md.  GOAL: a fully rigorous, variational
LiH energy E ~ -8.06 +/- few mHa that beats the pure-orbital ceiling -8.032, stays above
exact -8.0705 (variational; below = a BUG), double-counting-free by construction.

METHOD: E_VMC[Psi_CI * Jastrow] >= exact.  Psi_CI = the analytic-orbital FCI ground state
(the -8.032 ceiling wavefunction); Jastrow J = exp(sum_{i<j} u(r_ij)) supplies the e-e cusp
correlation the finite orbital basis structurally cannot reach (the ~30 mHa Li 1s^2 core
cusp + valence).  Sampled by Metropolis on |Psi_CI J|^2; gradient-form (bounded) kinetic
local energy (first derivatives only, no Laplacian cusp spikes).

The analytic orbital (prolate_allelectron_c4.OrbitalM):
  phi = N * xi^j * (xi^2-1)^{mu/2} * g(eta) * (1-eta^2)^{mu/2} * e^{-alpha xi} * e^{i s phi_azi}
KEY IDENTITY (this module):  (xi^2-1)^{mu/2}(1-eta^2)^{mu/2} e^{i s phi_azi}
                           = (2/R)^mu (x + i*sign(s)*y)^mu,   mu=|s|,
so phi = N * (2/R)^mu * xi^j e^{-alpha xi} * g(eta) * (x + i*sign(s)*y)^mu -- smooth in
(xi,eta) times a Cartesian polynomial, giving EXACT gradients with no on-axis singularity.

Coordinates (nucleus A=Li Z=3 at z=-R/2, B=H Z=1 at z=+R/2):
  rA=|r-(0,0,-R/2)|, rB=|r-(0,0,+R/2)|, xi=(rA+rB)/R, eta=(rA-rB)/R.

Run:  python debug/lih_vmc.py gate_orb     # step-2 orbital value+grad FD gate
      python debug/lih_vmc.py <more as built>
"""
from __future__ import annotations

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                    # debug/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # root


# ===========================================================================
# Step 2: arbitrary-point analytic orbital evaluator + real-space gradient
# ===========================================================================
def extract_prim(orb) -> dict:
    """Pull the float real-space parameters out of an OrbitalM (mpf fields -> float)."""
    return dict(
        N=float(orb.norm),
        j=int(orb.xi_power),
        alpha=float(orb.alpha),
        eta=np.array([float(c) for c in orb.eta_poly], dtype=float),  # g(eta)=sum c_k eta^k
        mu=int(orb.mu),
        s=int(orb.msign),
    )


def prolate_coords(r: np.ndarray, R: float):
    """r: (...,3) Cartesian -> (xi, eta, rA, rB) and the Cartesian gradients of xi, eta.

    Returns xi, eta (shape r.shape[:-1]) and gxi, geta (shape r.shape) = d xi/dr, d eta/dr.
    A (Li) at z=-R/2, B (H) at z=+R/2."""
    x = r[..., 0]; y = r[..., 1]; z = r[..., 2]
    dzA = z + R / 2.0
    dzB = z - R / 2.0
    rA = np.sqrt(x * x + y * y + dzA * dzA)
    rB = np.sqrt(x * x + y * y + dzB * dzB)
    xi = (rA + rB) / R
    eta = (rA - rB) / R
    # gradients of rA, rB
    grA = np.stack([x, y, dzA], axis=-1) / rA[..., None]
    grB = np.stack([x, y, dzB], axis=-1) / rB[..., None]
    gxi = (grA + grB) / R
    geta = (grA - grB) / R
    return xi, eta, rA, rB, gxi, geta


def orbital_value_grad(prim: dict, r: np.ndarray, R: float):
    """phi(r) and grad phi(r) for one analytic orbital, over a batch r: (...,3).

    Returns (val, grad): val shape r.shape[:-1] (complex), grad shape r.shape (complex).
    """
    N, j, alpha, eta_c, mu, s = (prim['N'], prim['j'], prim['alpha'],
                                 prim['eta'], prim['mu'], prim['s'])
    x = r[..., 0]; y = r[..., 1]
    xi, eta, rA, rB, gxi, geta = prolate_coords(r, R)

    # A(xi) = xi^j e^{-alpha xi};  A'(xi) = (j xi^{j-1} - alpha xi^j) e^{-alpha xi}
    exp = np.exp(-alpha * xi)
    xij = xi ** j
    A = xij * exp
    dA = (j * xi ** (j - 1) - alpha * xij) * exp if j > 0 else (-alpha * xij * exp)

    # B(eta) = g(eta) = sum c_k eta^k ;  B'(eta) = sum k c_k eta^{k-1}
    g = np.polynomial.polynomial.polyval(eta, eta_c)
    if len(eta_c) > 1:
        dg = np.polynomial.polynomial.polyval(eta, np.polynomial.polynomial.polyder(eta_c))
    else:
        dg = np.zeros_like(eta)

    pref = N * (2.0 / R) ** mu
    if mu == 0:
        W = np.ones_like(x, dtype=complex)
        gW = np.zeros(r.shape, dtype=complex)
    else:
        sign = 1.0 if s > 0 else -1.0
        w = x + 1j * sign * y                    # (x + i sign y)
        W = w ** mu
        dWdw = mu * w ** (mu - 1)
        gW = np.zeros(r.shape, dtype=complex)
        gW[..., 0] = dWdw * 1.0                   # d/dx (w) = 1
        gW[..., 1] = dWdw * (1j * sign)           # d/dy (w) = i sign
        # d/dz = 0

    # phi = pref * A(xi) * B(eta) * W
    val = pref * A * g * W
    # grad phi = pref [ A'(xi) B W gxi + A B'(eta) W geta + A B gW ]
    grad = (pref * (dA * g * W)[..., None] * gxi
            + pref * (A * dg * W)[..., None] * geta
            + pref * (A * g)[..., None] * gW)
    return val, grad


def orbital_vgl(prim: dict, r: np.ndarray, R: float):
    """phi, grad phi, AND laplacian(phi) analytically.  Uses: W=(x+i s y)^mu is
    holomorphic => lap_xy W = 0; grad(xi).grad(eta)=0 (orthogonal prolate gradients).
    Returns (val, grad, lap): val/lap shape r.shape[:-1] (complex), grad r.shape."""
    N, j, alpha, eta_c, mu, s = (prim['N'], prim['j'], prim['alpha'],
                                 prim['eta'], prim['mu'], prim['s'])
    x = r[..., 0]; y = r[..., 1]; z = r[..., 2]
    cAz = -R / 2.0; cBz = R / 2.0
    dzA = z - cAz; dzB = z - cBz
    rA = np.sqrt(x * x + y * y + dzA * dzA)
    rB = np.sqrt(x * x + y * y + dzB * dzB)
    xi = (rA + rB) / R
    eta = (rA - rB) / R
    uA = np.stack([x, y, dzA], axis=-1) / rA[..., None]
    uB = np.stack([x, y, dzB], axis=-1) / rB[..., None]
    gxi = (uA + uB) / R
    geta = (uA - uB) / R
    dotAB = np.sum(uA * uB, axis=-1)
    gxi2 = (2.0 + 2.0 * dotAB) / R ** 2          # |grad xi|^2
    geta2 = (2.0 - 2.0 * dotAB) / R ** 2         # |grad eta|^2
    lap_xi = (2.0 / rA + 2.0 / rB) / R
    lap_eta = (2.0 / rA - 2.0 / rB) / R

    # F(xi)=xi^j e^{-a xi}, F', F''
    exp = np.exp(-alpha * xi)
    xij = xi ** j
    F = xij * exp
    Fp = (j * xi ** (j - 1) - alpha * xij) * exp if j > 0 else (-alpha * xij * exp)
    if j >= 2:
        Fpp = (j * (j - 1) * xi ** (j - 2) - 2 * alpha * j * xi ** (j - 1)
               + alpha ** 2 * xij) * exp
    elif j == 1:
        Fpp = (-2 * alpha + alpha ** 2 * xi) * exp
    else:
        Fpp = alpha ** 2 * exp

    # G(eta)=poly, G', G''
    G = np.polynomial.polynomial.polyval(eta, eta_c)
    d1 = np.polynomial.polynomial.polyder(eta_c) if len(eta_c) > 1 else [0.0]
    d2 = np.polynomial.polynomial.polyder(eta_c, 2) if len(eta_c) > 2 else [0.0]
    Gp = np.polynomial.polynomial.polyval(eta, d1)
    Gpp = np.polynomial.polynomial.polyval(eta, d2)

    pref = N * (2.0 / R) ** mu
    if mu == 0:
        W = np.ones_like(x, dtype=complex)
        gW = np.zeros(r.shape, dtype=complex)
    else:
        sign = 1.0 if s > 0 else -1.0
        w = x + 1j * sign * y
        W = w ** mu
        dWdw = mu * w ** (mu - 1)
        gW = np.zeros(r.shape, dtype=complex)
        gW[..., 0] = dWdw
        gW[..., 1] = dWdw * (1j * sign)

    val = pref * F * G * W
    grad = (pref * (Fp * G * W)[..., None] * gxi
            + pref * (F * Gp * W)[..., None] * geta
            + pref * (F * G)[..., None] * gW)
    # laplacian: grad(xi).grad(eta)=0, lap_xy W = 0
    gxi_dot_gW = np.sum(gxi * gW, axis=-1)       # only x,y comps of gW nonzero
    geta_dot_gW = np.sum(geta * gW, axis=-1)
    lap = pref * (
        G * W * (Fpp * gxi2 + Fp * lap_xi)
        + F * W * (Gpp * geta2 + Gp * lap_eta)
        + 2.0 * Fp * G * gxi_dot_gW
        + 2.0 * F * Gp * geta_dot_gW
    )
    return val, grad, lap


def orbital_value(prim: dict, r: np.ndarray, R: float):
    """phi(r) value ONLY (no gradient) -- for the Metropolis accept step."""
    N, j, alpha, eta_c, mu, s = (prim['N'], prim['j'], prim['alpha'],
                                 prim['eta'], prim['mu'], prim['s'])
    x = r[..., 0]; y = r[..., 1]; z = r[..., 2]
    dzA = z + R / 2.0; dzB = z - R / 2.0
    rA = np.sqrt(x * x + y * y + dzA * dzA)
    rB = np.sqrt(x * x + y * y + dzB * dzB)
    xi = (rA + rB) / R; eta = (rA - rB) / R
    A = xi ** j * np.exp(-alpha * xi)
    g = np.polynomial.polynomial.polyval(eta, eta_c)
    pref = N * (2.0 / R) ** mu
    if mu == 0:
        return pref * A * g
    sign = 1.0 if s > 0 else -1.0
    W = (x + 1j * sign * y) ** mu
    return pref * A * g * W


def eval_values_batch(prims, r: np.ndarray, R: float):
    """Values of ALL primitive orbitals over a batch r:(...,3), computing the shared
    prolate coordinates ONCE.  Returns Phi (Nprim, *batch) complex."""
    x = r[..., 0]; y = r[..., 1]; z = r[..., 2]
    dzA = z + R / 2.0; dzB = z - R / 2.0
    rho2 = x * x + y * y
    rA = np.sqrt(rho2 + dzA * dzA)
    rB = np.sqrt(rho2 + dzB * dzB)
    xi = (rA + rB) / R; eta = (rA - rB) / R
    logxi = np.log(xi)
    out = np.empty((len(prims),) + r.shape[:-1], dtype=complex)
    for a, p in enumerate(prims):
        j = p['j']; alpha = p['alpha']; mu = p['mu']
        A = np.exp(j * logxi - alpha * xi) if j else np.exp(-alpha * xi)
        g = np.polynomial.polynomial.polyval(eta, p['eta'])
        val = (p['N'] * (2.0 / R) ** mu) * A * g
        if mu:
            sign = 1.0 if p['s'] > 0 else -1.0
            val = val * (x + 1j * sign * y) ** mu
        out[a] = val
    return out


def eval_orbitals(prims, r: np.ndarray, R: float):
    """Evaluate a list of primitive orbitals at electron positions r: (Ne,3).

    Returns Phi (Nprim, Ne) complex, gPhi (Nprim, Ne, 3) complex."""
    Nprim = len(prims)
    Ne = r.shape[0]
    Phi = np.empty((Nprim, Ne), dtype=complex)
    gPhi = np.empty((Nprim, Ne, 3), dtype=complex)
    for a, p in enumerate(prims):
        v, g = orbital_value_grad(p, r, R)
        Phi[a] = v
        gPhi[a] = g
    return Phi, gPhi


# ===========================================================================
# Gate: analytic gradient vs finite difference (step 2)
# ===========================================================================
def gate_orb():
    """FD gate on orbital_value_grad for sigma-core, sigma-valence, pi orbitals."""
    import mpmath as mp
    from prolate_mixed_eri import sto_orbital, valence_prolate_orbital, ZC_LI
    from prolate_allelectron_analytic_fci import sto_orbital_B
    from prolate_allelectron_c4 import from_sigma, valence_pi_orbital, _orb_on_grid

    R = 3.015
    orbs = {
        "Li 1s core (z=2.6875)": from_sigma(sto_orbital(ZC_LI, R, is_core=True)),
        "Li core2 (z=4.5)":      from_sigma(sto_orbital(mp.mpf('4.5'), R, is_core=True)),
        "H 1s (z=1.0)":          from_sigma(sto_orbital_B(mp.mpf('1.0'), R)),
        "bond sigma (0,0)":      from_sigma(valence_prolate_orbital(0, 0, mp.mpf('1.0'))),
        "bond sigma (2,1)":      from_sigma(valence_prolate_orbital(2, 1, mp.mpf('1.0'))),
        "pi +1 (0,0)":           valence_pi_orbital(0, 0, mp.mpf('1.0'), +1),
        "pi -1 (1,0)":           valence_pi_orbital(1, 0, mp.mpf('1.0'), -1),
    }
    rng = np.random.default_rng(1)
    # random points away from nuclei/axis
    pts = rng.uniform(-2.0, 2.0, size=(6, 3))
    pts[:, 2] += rng.uniform(-1.0, 1.0, size=6)
    h = 1e-6
    print("=" * 74)
    print("GATE step-2: orbital value + analytic gradient vs finite difference")
    print("=" * 74)
    ok = True
    for name, orb in orbs.items():
        prim = extract_prim(orb)
        val, grad = orbital_value_grad(prim, pts, R)
        # FD gradient
        gfd = np.empty_like(grad)
        for d in range(3):
            dp = np.zeros(3); dp[d] = h
            vp, _ = orbital_value_grad(prim, pts + dp, R)
            vm, _ = orbital_value_grad(prim, pts - dp, R)
            gfd[:, d] = (vp - vm) / (2 * h)
        err = np.max(np.abs(grad - gfd)) / max(np.max(np.abs(grad)), 1e-300)
        # cross-check the VALUE against the existing grid evaluator _orb_on_grid
        # (build a 1-point "grid" is awkward; instead check value magnitude is finite)
        finite = np.all(np.isfinite(val)) and np.all(np.isfinite(grad))
        tag = "PASS" if (err < 1e-6 and finite) else "FAIL"
        ok &= (err < 1e-6 and finite)
        print(f"  {name:24s}: max|grad-FD|/|grad| = {err:.2e}   {tag}")
    print(f"\n  GATE step-2: {'PASS' if ok else 'FAIL'}")
    return ok


def gate_value_vs_grid():
    """Cross-check orbital_value_grad VALUE against the engine's own _orb_on_grid,
    proving the real-space evaluator IS the same function the ERIs are built from."""
    import mpmath as mp
    from prolate_mixed_eri import sto_orbital, valence_prolate_orbital, ZC_LI
    from prolate_allelectron_c4 import from_sigma, valence_pi_orbital, _orb_on_grid, _grid_template

    R = 3.015
    tmpl = _grid_template(R, N_grid=24, xi_max=10.0)
    xi = tmpl['xi']; eta = tmpl['eta']
    # pick a few interior grid nodes, convert (xi,eta,phi=0) -> Cartesian, compare
    orbs = {
        "Li 1s core": from_sigma(sto_orbital(ZC_LI, R, is_core=True)),
        "bond (2,1)":  from_sigma(valence_prolate_orbital(2, 1, mp.mpf('1.0'))),
        "pi +1 (0,0)": valence_pi_orbital(0, 0, mp.mpf('1.0'), +1),
    }
    print("=" * 74)
    print("GATE step-2b: real-space value vs engine _orb_on_grid (phi_azi=0 slice)")
    print("=" * 74)
    ok = True
    ii, jj = len(xi) // 2, len(eta) // 2
    xi0, eta0 = xi[ii], eta[jj]
    # phi_azi = 0 -> x = rho, y = 0, with rho=(R/2)sqrt((xi^2-1)(1-eta^2)), z=(R/2)xi eta
    rho = (R / 2) * np.sqrt(max((xi0**2 - 1) * (1 - eta0**2), 0.0))
    zc = (R / 2) * xi0 * eta0
    rc = np.array([[rho, 0.0, zc]])
    for name, orb in orbs.items():
        prim = extract_prim(orb)
        val, _ = orbital_value_grad(prim, rc, R)
        grid = _orb_on_grid(orb, tmpl)['psi'][ii, jj]   # psi(xi,eta) part (phi stripped)
        # at phi_azi=0, e^{i s phi}=1, so val should equal grid (real) for this slice
        err = abs(val[0] - grid) / max(abs(grid), 1e-300)
        tag = "PASS" if err < 1e-9 else "FAIL"
        ok &= err < 1e-9
        print(f"  {name:14s}: real-space={val[0].real:+.8e}  grid={grid:+.8e}  rel={err:.1e} {tag}")
    print(f"\n  GATE step-2b: {'PASS' if ok else 'FAIL'}  (xi0={xi0:.3f}, eta0={eta0:.3f})")
    return ok


# ===========================================================================
# Step 1: FCI ground-state eigenVECTOR (fci_fast returns eigenvalue only)
# ===========================================================================
def fci_ground_vector(h1, eri, M, nelec):
    """Sparse connected-pair FCI (mirrors fci_fast.fci_energy_fast) returning the
    ground eigenvector + the determinant list, in the SAME _dets ordering.
    Returns (E0, c[nd], dets)."""
    import scipy.sparse as sp
    from scipy.sparse.linalg import eigsh
    from prolate_allelectron_fci import _dets, _matel
    from fci_fast import _popcount64, _FOUR

    na = nb = nelec // 2
    dets = _dets(M, na, nb)
    nd = len(dets)
    if nd == 1:
        e = _matel(h1, eri, dets[0], dets[0])
        return e, np.array([1.0]), dets
    masks = np.array([int(sum(1 << so for so in D)) for D in dets], dtype=np.uint64)
    rl, cl, vl = [], [], []
    diag = np.empty(nd)
    for a in range(nd):
        diag[a] = _matel(h1, eri, dets[a], dets[a])
        if a + 1 < nd:
            pc = _popcount64(masks[a + 1:] ^ masks[a])
            conn = np.nonzero(pc <= _FOUR)[0] + (a + 1)
            Da = dets[a]
            for b in conn.tolist():
                v = _matel(h1, eri, Da, dets[b])
                if v != 0.0:
                    rl.append(a); cl.append(b); vl.append(v)
    r = np.array(rl, dtype=np.int64); c_ = np.array(cl, dtype=np.int64); v = np.array(vl)
    ri = np.concatenate([np.arange(nd), r, c_])
    ci = np.concatenate([np.arange(nd), c_, r])
    vi = np.concatenate([diag, v, v])
    H = sp.csr_matrix((vi, (ri, ci)), shape=(nd, nd))
    ev, evec = eigsh(H, k=1, which='SA', return_eigenvectors=True)
    return float(ev[0]), evec[:, 0], dets


# ===========================================================================
# Step 1+3 support: reconstruct the analytic-orbital LiH FCI wavefunction
# ===========================================================================
class LiHWavefunction:
    """Everything needed to evaluate Psi_CI(r1..r4) in real space.
    prims  : list of Mk orthonormal-MO primitive-expansion data? NO -- we keep the
             M raw primitives + the T=(C.T@X) map, MO_p = sum_a T[a,p] prim_a.
    """
    def __init__(self, prims, T, c, apairs, bpairs, CS, R, Z_A, Z_B, Vnn,
                 E_fci, Mk, na, nb):
        self.prims = prims        # list of M primitive dicts (extract_prim)
        self.T = T                # (M, Mk) primitive -> orthonormal MO
        self.c = c                # (nd,) CI coefficients
        self.apairs = apairs      # list of alpha spatial-orbital tuples
        self.bpairs = bpairs      # list of beta  spatial-orbital tuples
        self.CS = CS              # (n_ap, n_bp) CI matrix folded with det sign sigma
        self.R = R; self.Z_A = Z_A; self.Z_B = Z_B; self.Vnn = Vnn
        self.E_fci = E_fci; self.Mk = Mk; self.na = na; self.nb = nb


def _perm_sign_inv(seq):
    """(-1)^{number of inversions} of an integer sequence."""
    n = len(seq); inv = 0
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                inv += 1
    return -1.0 if (inv & 1) else 1.0


def build_lih_wavefunction(Jb=2, Lb=1, npi=1, Jpi=1, Lpi=0, alpha=1.0,
                           core2=(4.5, 1.6), tol=1e-11, verbose=True):
    """Reconstruct the assemble_rebased pipeline, extract the CI eigenvector, and
    package the real-space wavefunction. Uses the float64 ERI engine (fast)."""
    import mpmath as mp
    import prolate_energy_ladder as L
    import lih_core2exp_probe as P
    from prolate_allelectron_c4 import one_body_general
    from prolate_float_eri import build_eri_tensor_m_f
    from itertools import combinations

    R = P.R; Z_A, Z_B = 3.0, 1.0; nelec = 4; Vnn = Z_A * Z_B / R
    c2 = list(core2) if core2 else None
    orbs, tags = P.lih_orbs(Jb, Lb, npi, Jpi, Lpi, alpha, c2)
    M = len(orbs)
    # rebasing (C) + one-body + ERI, exactly as assemble_rebased
    C = L.build_C(tags, alpha)
    S, h1 = one_body_general(orbs, R, Z_A, Z_B)
    eri = build_eri_tensor_m_f(orbs, R, verbose=verbose)
    So = C @ S @ C.T
    h1o = C @ h1 @ C.T
    erio = np.einsum('pa,qb,rc,sd,abcd->pqrs', C, C, C, C, eri, optimize=True)
    w, U = np.linalg.eigh(So)
    keep = w > tol * w[-1]
    X = U[:, keep] / np.sqrt(w[keep])
    Mk = int(keep.sum())
    h1_f = X.T @ h1o @ X
    eri_f = np.einsum('ap,bq,cr,ds,abcd->pqrs', X, X, X, X, erio, optimize=True)
    T = C.T @ X                                    # (M, Mk) primitive -> MO

    E0, c, dets = fci_ground_vector(h1_f, eri_f, Mk, nelec)
    na = nb = nelec // 2
    # Rayleigh cross-check (step-1 gate): c must have E0 as Rayleigh quotient
    E_tot = E0 + Vnn

    # determinant (alpha-pair, beta-pair) enumeration in the SAME order as _dets
    apairs = list(combinations(range(Mk), na))
    bpairs = list(combinations(range(Mk), nb))
    n_ap, n_bp = len(apairs), len(bpairs)
    assert n_ap * n_bp == len(dets), "determinant count mismatch"
    Cmat = c.reshape(n_ap, n_bp)                   # _dets loops a-outer, b-inner
    # sign sigma_I for each (a,b): grouped spin-orbital order [2*a..., 2*b+1...]
    sigma = np.empty((n_ap, n_bp))
    for ia, a in enumerate(apairs):
        for ib, b in enumerate(bpairs):
            grouped = [2 * p for p in a] + [2 * q + 1 for q in b]
            sigma[ia, ib] = _perm_sign_inv(grouped)
    CS = Cmat * sigma

    prims = [extract_prim(o) for o in orbs]
    if verbose:
        print(f"  build_lih_wavefunction: M={M} kept Mk={Mk} nd={len(dets)} "
              f"E_elec={E0:.5f} E_tot={E_tot:.5f}")
        print(f"    top |c|: {np.sort(np.abs(c))[::-1][:5].round(4)}  sum c^2={np.sum(c**2):.6f}")
    return LiHWavefunction(prims, T, c, apairs, bpairs, CS, R, Z_A, Z_B, Vnn,
                           E0, Mk, na, nb)


# ===========================================================================
# Step 3: multi-determinant Psi_CI evaluator (value + per-electron gradient),
# walker-batched.  Positions r: (Nw, Ne, 3), electrons ordered [alpha..., beta...].
# ===========================================================================
def _minors_and_grads(MOval, gMO, pairs, elec_idx):
    """For na=2: det of the 2x2 [MO_p(e0) MO_p(e1); MO_q(e0) MO_q(e1)] per pair, and its
    gradient wrt each of the 2 electrons.  MOval:(Mk,Nw,Ne) gMO:(Mk,Nw,Ne,3).
    elec_idx: the 2 electron indices for this spin.  Returns:
      dets (npairs,Nw), gdets (2, npairs, Nw, 3)  [grad wrt local electron 0,1]."""
    e0, e1 = elec_idx
    npairs = len(pairs)
    Nw = MOval.shape[1]
    dets = np.empty((npairs, Nw), dtype=complex)
    gdets = np.empty((2, npairs, Nw, 3), dtype=complex)
    for k, (p, q) in enumerate(pairs):
        Pp0 = MOval[p, :, e0]; Pp1 = MOval[p, :, e1]
        Pq0 = MOval[q, :, e0]; Pq1 = MOval[q, :, e1]
        dets[k] = Pp0 * Pq1 - Pq0 * Pp1
        gPp0 = gMO[p, :, e0, :]; gPp1 = gMO[p, :, e1, :]
        gPq0 = gMO[q, :, e0, :]; gPq1 = gMO[q, :, e1, :]
        # d/d(elec e0): p1,q1 fixed
        gdets[0, k] = gPp0 * Pq1[:, None] - gPq0 * Pp1[:, None]
        # d/d(elec e1): p0,q0 fixed
        gdets[1, k] = Pp0[:, None] * gPq1 - Pq0[:, None] * gPp1
    return dets, gdets


def _minors_full(MOval, gMO, lMO, pairs, elec_idx):
    """na=2 minors + per-electron gradient AND laplacian.
    Returns dets(npairs,Nw), gdets(2,npairs,Nw,3), ldets(2,npairs,Nw)."""
    e0, e1 = elec_idx
    npairs = len(pairs); Nw = MOval.shape[1]
    dets = np.empty((npairs, Nw), dtype=complex)
    gdets = np.empty((2, npairs, Nw, 3), dtype=complex)
    ldets = np.empty((2, npairs, Nw), dtype=complex)
    for k, (p, q) in enumerate(pairs):
        Pp0 = MOval[p, :, e0]; Pp1 = MOval[p, :, e1]
        Pq0 = MOval[q, :, e0]; Pq1 = MOval[q, :, e1]
        dets[k] = Pp0 * Pq1 - Pq0 * Pp1
        gdets[0, k] = gMO[p, :, e0, :] * Pq1[:, None] - gMO[q, :, e0, :] * Pp1[:, None]
        gdets[1, k] = Pp0[:, None] * gMO[q, :, e1, :] - Pq0[:, None] * gMO[p, :, e1, :]
        ldets[0, k] = lMO[p, :, e0] * Pq1 - lMO[q, :, e0] * Pp1
        ldets[1, k] = Pp0 * lMO[q, :, e1] - Pq0 * lMO[p, :, e1]
    return dets, gdets, ldets


def psi_ci_full(wf: LiHWavefunction, r: np.ndarray):
    """Psi_CI value, per-electron gradient, per-electron Laplacian (all analytic).
    Returns psi(Nw,), gpsi(Nw,Ne,3), lpsi(Nw,Ne)  [lpsi_i = lap_i Psi_CI]."""
    R = wf.R; Nw, Ne, _ = r.shape
    M = len(wf.prims)
    Phi = np.empty((M, Nw, Ne), dtype=complex)
    gPhi = np.empty((M, Nw, Ne, 3), dtype=complex)
    lPhi = np.empty((M, Nw, Ne), dtype=complex)
    for a, p in enumerate(wf.prims):
        v, g, l = orbital_vgl(p, r, R)
        Phi[a] = v; gPhi[a] = g; lPhi[a] = l
    MOval = np.einsum('ap,aWe->pWe', wf.T, Phi)
    gMO = np.einsum('ap,aWed->pWed', wf.T, gPhi)
    lMO = np.einsum('ap,aWe->pWe', wf.T, lPhi)
    na, nb = wf.na, wf.nb
    a_e = list(range(na)); b_e = list(range(na, na + nb))
    dA, gdA, ldA = _minors_full(MOval, gMO, lMO, wf.apairs, a_e)
    dB, gdB, ldB = _minors_full(MOval, gMO, lMO, wf.bpairs, b_e)
    CS = wf.CS
    tmpB = np.einsum('ab,bW->aW', CS, dB)          # sum_b CS[a,b] dB[b]
    tmpA = np.einsum('ab,aW->bW', CS, dA)
    psi = np.einsum('aW,aW->W', dA, tmpB)
    gpsi = np.zeros((Nw, Ne, 3), dtype=complex)
    lpsi = np.zeros((Nw, Ne), dtype=complex)
    for le, ge in enumerate(a_e):
        gpsi[:, ge, :] = np.einsum('aWd,aW->Wd', gdA[le], tmpB)
        lpsi[:, ge] = np.einsum('aW,aW->W', ldA[le], tmpB)
    for le, ge in enumerate(b_e):
        gpsi[:, ge, :] = np.einsum('bWd,bW->Wd', gdB[le], tmpA)
        lpsi[:, ge] = np.einsum('bW,bW->W', ldB[le], tmpA)
    return psi, gpsi, lpsi


def local_energy_analytic(wf, r, jas=None):
    """STANDARD local energy with ANALYTIC Laplacian.  Returns (E_L real(Nw,), psi)."""
    psi, gpsi, lpsi = psi_ci_full(wf, r)
    lap_over_psi = np.sum(lpsi, axis=1) / psi          # sum_i lap_i Psi_CI / Psi_CI
    V = potential(wf, r)
    if jas is None:
        return (V - 0.5 * lap_over_psi).real, psi
    glnpsi = gpsi / psi[:, None, None]
    glnJ = jas.grad_lnJ(r)
    lapJ = jas.lap_lnJ(r)
    cross = 2.0 * np.sum(glnpsi * glnJ, axis=(1, 2))
    gJ2 = np.sum(glnJ ** 2, axis=(1, 2))
    lap_tot = lap_over_psi + cross + lapJ + gJ2
    return (V - 0.5 * lap_tot).real, psi


def psi_ci(wf: LiHWavefunction, r: np.ndarray):
    """Psi_CI value and per-electron gradient. r:(Nw,Ne,3) -> (psi(Nw,), gpsi(Nw,Ne,3))."""
    R = wf.R
    Nw, Ne, _ = r.shape
    # orbital values at all electrons: Phi (M, Nw, Ne), gPhi (M, Nw, Ne, 3)
    M = len(wf.prims)
    Phi = np.empty((M, Nw, Ne), dtype=complex)
    gPhi = np.empty((M, Nw, Ne, 3), dtype=complex)
    for a, p in enumerate(wf.prims):
        v, g = orbital_value_grad(p, r, R)          # v:(Nw,Ne) g:(Nw,Ne,3)
        Phi[a] = v; gPhi[a] = g
    # MO values: MOval (Mk,Nw,Ne)
    MOval = np.einsum('ap,aWe->pWe', wf.T, Phi)
    gMO = np.einsum('ap,aWed->pWed', wf.T, gPhi)
    na, nb = wf.na, wf.nb
    a_e = list(range(na)); b_e = list(range(na, na + nb))
    dA, gdA = _minors_and_grads(MOval, gMO, wf.apairs, a_e)   # (n_ap,Nw),(2,n_ap,Nw,3)
    dB, gdB = _minors_and_grads(MOval, gMO, wf.bpairs, b_e)
    CS = wf.CS
    # psi = sum_ab CS[a,b] dA[a] dB[b] = dA^T CS dB   (per walker)
    tmp = np.einsum('ab,bW->aW', CS, dB)             # (n_ap,Nw)
    psi = np.einsum('aW,aW->W', dA, tmp)             # (Nw,)
    # gradients
    gpsi = np.zeros((Nw, Ne, 3), dtype=complex)
    # alpha electrons: grad wrt local a-electron le (0..na-1)
    for le, ge in enumerate(a_e):
        # d psi / d r_ge = (gdA[le])^T CS dB
        gpsi[:, ge, :] = np.einsum('aWd,aW->Wd', gdA[le], tmp)
    tmpB = np.einsum('ab,aW->bW', CS, dA)            # (n_bp,Nw)
    for le, ge in enumerate(b_e):
        gpsi[:, ge, :] = np.einsum('bWd,bW->Wd', gdB[le], tmpB)
    return psi, gpsi


# ===========================================================================
# Potential + gradient-form local energy
# ===========================================================================
def potential(wf: LiHWavefunction, r: np.ndarray):
    """V(r): (Nw,) real.  -Z_A/rA - Z_B/rB per electron + sum 1/r_ij + Vnn."""
    R = wf.R; Nw, Ne, _ = r.shape
    zc = np.array([0.0, 0.0, -R / 2.0]); zh = np.array([0.0, 0.0, R / 2.0])
    rA = np.linalg.norm(r - zc, axis=-1)             # (Nw,Ne)
    rB = np.linalg.norm(r - zh, axis=-1)
    Vne = -np.sum(wf.Z_A / rA + wf.Z_B / rB, axis=1)
    Vee = np.zeros(Nw)
    for i in range(Ne):
        for j in range(i + 1, Ne):
            Vee += 1.0 / np.linalg.norm(r[:, i] - r[:, j], axis=-1)
    return Vne + Vee + wf.Vnn


def _wf_to_dict(wf):
    return dict(prims=wf.prims, T=wf.T, c=wf.c, apairs=wf.apairs, bpairs=wf.bpairs,
                CS=wf.CS, R=wf.R, Z_A=wf.Z_A, Z_B=wf.Z_B, Vnn=wf.Vnn,
                E_fci=wf.E_fci, Mk=wf.Mk, na=wf.na, nb=wf.nb)


def _wf_from_dict(d):
    return LiHWavefunction(d['prims'], d['T'], d['c'], d['apairs'], d['bpairs'],
                           d['CS'], d['R'], d['Z_A'], d['Z_B'], d['Vnn'],
                           d['E_fci'], d['Mk'], d['na'], d['nb'])


def cached_wavefunction(tag, **kw):
    """Build (or load) a LiHWavefunction, cached as a MODULE-INDEPENDENT plain dict."""
    import pickle
    DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, f"lih_vmc_wf_{tag}.pkl")
    if os.path.exists(path):
        with open(path, "rb") as f:
            obj = pickle.load(f)
        return _wf_from_dict(obj) if isinstance(obj, dict) else obj  # dict or legacy
    wf = build_lih_wavefunction(**kw)
    with open(path, "wb") as f:
        pickle.dump(_wf_to_dict(wf), f)
    return wf


# ===========================================================================
# Step 4: Jastrow  J = exp(sum_{i<j} u(r_ij)),  u = cusp-correct linexp
# u_ss(r) = b_ss * r * exp(-gamma r);  b = 1/4 (parallel), 1/2 (antiparallel spin).
# For LiH electrons ordered [a-up, a-up, b-dn, b-dn]: pairs (0,1)=UU parallel,
# (2,3)=DD parallel, the four cross pairs = UD antiparallel.
# ===========================================================================
class Jastrow:
    def __init__(self, gamma, na, nb, cusp=True, amp=1.0):
        """u(r) = amp * b_spin * r * e^{-gamma r};  b_spin = 1/2 (antiparallel), 1/4
        (parallel).  amp=1 gives the exact e-e cusp; amp>1 deepens the hole (variational,
        cusp = amp*b_spin -- allowed for a trial wavefunction)."""
        self.gamma = gamma
        self.na = na; self.nb = nb
        self.Ne = na + nb
        self.b = np.zeros((self.Ne, self.Ne))
        for i in range(self.Ne):
            for j in range(self.Ne):
                same = (i < na) == (j < na)
                self.b[i, j] = 0.25 if same else 0.5
        if not cusp:
            self.b[:] = 0.5
        self.b = self.b * amp
        self.amp = amp

    def _u(self, rij, b):
        g = self.gamma
        return b * rij * np.exp(-g * rij)

    def _du(self, rij, b):
        g = self.gamma
        return b * np.exp(-g * rij) * (1.0 - g * rij)

    def lnJ(self, r):
        """sum_{i<j} u(r_ij).  r:(Nw,Ne,3) -> (Nw,)."""
        Nw, Ne, _ = r.shape
        tot = np.zeros(Nw)
        for i in range(Ne):
            for j in range(i + 1, Ne):
                d = np.linalg.norm(r[:, i] - r[:, j], axis=-1)
                tot += self._u(d, self.b[i, j])
        return tot

    def grad_lnJ(self, r):
        """grad_i lnJ = sum_{j!=i} u'(r_ij) (r_i-r_j)/r_ij.  -> (Nw,Ne,3) real."""
        Nw, Ne, _ = r.shape
        g = np.zeros((Nw, Ne, 3))
        for i in range(Ne):
            for j in range(Ne):
                if i == j:
                    continue
                dvec = r[:, i] - r[:, j]
                d = np.linalg.norm(dvec, axis=-1)
                d = np.where(d < 1e-12, 1e-12, d)
                g[:, i] += (self._du(d, self.b[i, j]) / d)[:, None] * dvec
        return g

    def lap_lnJ(self, r):
        """sum_i lap_i lnJ.  Per pair, lap_i u + lap_j u = 2(u'' + 2u'/r); the 2b/r
        cusp (r->0) cancels V's +1/r_ij when b=1/2 (antiparallel).  -> (Nw,) real."""
        Nw, Ne, _ = r.shape
        g = self.gamma
        tot = np.zeros(Nw)
        for i in range(Ne):
            for j in range(i + 1, Ne):
                d = np.linalg.norm(r[:, i] - r[:, j], axis=-1)
                d = np.where(d < 1e-12, 1e-12, d)
                e = self.b[i, j] * np.exp(-g * d)
                lap_u = e * (2.0 / d - 4.0 * g + g * g * d)   # u'' + 2u'/r
                tot += 2.0 * lap_u
        return tot


class PadeJastrow:
    """Standard 2-body Pade Jastrow  u(r) = A*r/(1+gamma r),  A = b_spin (1/2 anti-,
    1/4 parallel).  EXACT cusp (u'(0)=A) => finite-variance E_L; deeper, tunable hole
    (plateau A/gamma) than the cusp-fixed linexp.  This is the correct way to add hole
    depth without breaking the cusp."""
    def __init__(self, gamma, na, nb):
        self.gamma = gamma; self.na = na; self.nb = nb; self.Ne = na + nb
        self.A = np.zeros((self.Ne, self.Ne))
        for i in range(self.Ne):
            for j in range(self.Ne):
                same = (i < na) == (j < na)
                self.A[i, j] = 0.25 if same else 0.5

    def lnJ(self, r):
        Nw, Ne, _ = r.shape; tot = np.zeros(Nw); g = self.gamma
        for i in range(Ne):
            for j in range(i + 1, Ne):
                d = np.linalg.norm(r[:, i] - r[:, j], axis=-1)
                tot += self.A[i, j] * d / (1.0 + g * d)
        return tot

    def grad_lnJ(self, r):
        Nw, Ne, _ = r.shape; g = self.gamma; G = np.zeros((Nw, Ne, 3))
        for i in range(Ne):
            for j in range(Ne):
                if i == j:
                    continue
                dvec = r[:, i] - r[:, j]; d = np.linalg.norm(dvec, axis=-1)
                d = np.where(d < 1e-12, 1e-12, d)
                up = self.A[i, j] / (1.0 + g * d) ** 2      # u'(r)
                G[:, i] += (up / d)[:, None] * dvec
        return G

    def lap_lnJ(self, r):
        Nw, Ne, _ = r.shape; g = self.gamma; tot = np.zeros(Nw)
        for i in range(Ne):
            for j in range(i + 1, Ne):
                d = np.linalg.norm(r[:, i] - r[:, j], axis=-1)
                d = np.where(d < 1e-12, 1e-12, d)
                A = self.A[i, j]
                upp = -2.0 * A * g / (1.0 + g * d) ** 3      # u''
                uop = A / (1.0 + g * d) ** 2 / d             # u'/r
                tot += 2.0 * (upp + 2.0 * uop)               # both electrons
        return tot


class TwoPadeJastrow:
    """Two-scale Pade: u(r) = A1 r/(1+g1 r) + A2 r/(1+g2 r), A1=f*b_spin (short, g1
    large = CORE hole), A2=(1-f)*b_spin (long, g2 small = VALENCE hole).  A1+A2=b_spin
    => EXACT cusp preserved (finite variance).  Serves the tight core pair AND the
    diffuse valence pair with one translationally-invariant u (they sit at different r)."""
    def __init__(self, g1, g2, f, na, nb):
        self.g1 = g1; self.g2 = g2; self.f = f
        self.na = na; self.nb = nb; self.Ne = na + nb
        self.bs = np.zeros((self.Ne, self.Ne))
        for i in range(self.Ne):
            for j in range(self.Ne):
                same = (i < na) == (j < na)
                self.bs[i, j] = 0.25 if same else 0.5

    def lnJ(self, r):
        Nw, Ne, _ = r.shape; tot = np.zeros(Nw)
        for i in range(Ne):
            for j in range(i + 1, Ne):
                d = np.linalg.norm(r[:, i] - r[:, j], axis=-1)
                A1 = self.f * self.bs[i, j]; A2 = (1 - self.f) * self.bs[i, j]
                tot += A1 * d / (1 + self.g1 * d) + A2 * d / (1 + self.g2 * d)
        return tot

    def grad_lnJ(self, r):
        Nw, Ne, _ = r.shape; G = np.zeros((Nw, Ne, 3))
        for i in range(Ne):
            for j in range(Ne):
                if i == j:
                    continue
                dvec = r[:, i] - r[:, j]; d = np.linalg.norm(dvec, axis=-1)
                d = np.where(d < 1e-12, 1e-12, d)
                A1 = self.f * self.bs[i, j]; A2 = (1 - self.f) * self.bs[i, j]
                up = A1 / (1 + self.g1 * d) ** 2 + A2 / (1 + self.g2 * d) ** 2
                G[:, i] += (up / d)[:, None] * dvec
        return G

    def lap_lnJ(self, r):
        Nw, Ne, _ = r.shape; tot = np.zeros(Nw)
        for i in range(Ne):
            for j in range(i + 1, Ne):
                d = np.linalg.norm(r[:, i] - r[:, j], axis=-1)
                d = np.where(d < 1e-12, 1e-12, d)
                A1 = self.f * self.bs[i, j]; A2 = (1 - self.f) * self.bs[i, j]
                lap = 0.0
                for A, g in ((A1, self.g1), (A2, self.g2)):
                    upp = -2.0 * A * g / (1 + g * d) ** 3
                    uop = A / (1 + g * d) ** 2 / d
                    lap = lap + upp + 2.0 * uop
                tot += 2.0 * lap
        return tot


class PadePlus3:
    """2-body Pade (exact e-e cusp) + a 3-body electron-electron-nucleus term:
      U3 = -c * sum_{i<j} sum_I exp(-beta(r_iI^2 + r_jI^2)) exp(-delta r_ij^2)
    All-Gaussian => smooth at BOTH the e-e and e-n coalescences (no cusp modification,
    finite variance).  It deepens the correlation hole specifically when both electrons
    are near the SAME nucleus (the core pair) -- the e-e-n correlation a pure u(r12)
    Jastrow structurally cannot represent (= what the additive-F12 rich basis captured).
    3-body grad/lap by finite difference (the 3-body part is smooth; the 2-body cusp
    stays analytic)."""
    def __init__(self, gamma, c, beta, delta, na, nb, nuclei, h=1e-4):
        self.pade = PadeJastrow(gamma, na, nb)
        self.c = c; self.beta = beta; self.delta = delta
        self.na = na; self.nb = nb; self.Ne = na + nb
        self.nuclei = [np.asarray(n, float) for n in nuclei]
        self.h = h

    def _U3(self, r):
        Nw, Ne, _ = r.shape
        tot = np.zeros(Nw)
        G = []
        for I in self.nuclei:
            rI2 = np.sum((r - I) ** 2, axis=-1)          # (Nw,Ne)
            G.append(np.exp(-self.beta * rI2))
        for i in range(Ne):
            for j in range(i + 1, Ne):
                dij2 = np.sum((r[:, i] - r[:, j]) ** 2, axis=-1)
                H = np.exp(-self.delta * dij2)
                gij = np.zeros(Nw)
                for GI in G:
                    gij = gij + GI[:, i] * GI[:, j]
                tot += -self.c * gij * H
        return tot

    def lnJ(self, r):
        return self.pade.lnJ(r) + self._U3(r)

    def _fd_grad_lap(self, r):
        """central-diff grad and 2nd-diff lap of U3 in one pass."""
        Nw, Ne, _ = r.shape; h = self.h
        g = np.zeros((Nw, Ne, 3)); lap = np.zeros(Nw)
        U0 = self._U3(r)
        for i in range(Ne):
            for d in range(3):
                rp = r.copy(); rp[:, i, d] += h
                rm = r.copy(); rm[:, i, d] -= h
                Up = self._U3(rp); Um = self._U3(rm)
                g[:, i, d] = (Up - Um) / (2 * h)
                lap += (Up + Um - 2 * U0) / h ** 2
        return g, lap

    def grad_lnJ(self, r):
        g3, _ = self._fd_grad_lap(r)
        return self.pade.grad_lnJ(r) + g3

    def lap_lnJ(self, r):
        _, l3 = self._fd_grad_lap(r)
        return self.pade.lap_lnJ(r) + l3


class JastrowMulti:
    """Two-scale cusp Jastrow: u(r) = b*r*e^{-g1 r} + a*r*e^{-g2 r}.
    The cusp is carried entirely by the first term (u'(0)=b*b_spin, spin cusp b_spin);
    the second term (a, g2) adds a diffuse valence hole WITHOUT changing the cusp only
    if we keep its contribution to u'(0) accounted -- here u'(0) = (b + a) * scale.  To
    keep the exact cusp we fix the TOTAL short-range slope: b + a = 1 (times spin factor),
    so a is the single free shape parameter (fraction moved to the long scale)."""
    def __init__(self, gamma1, gamma2, a_long, na, nb):
        self.g1 = gamma1; self.g2 = gamma2
        self.na = na; self.nb = nb; self.Ne = na + nb
        # per-pair spin cusp factor: 0.5 antiparallel, 0.25 parallel
        self.bs = np.zeros((self.Ne, self.Ne))
        for i in range(self.Ne):
            for j in range(self.Ne):
                same = (i < na) == (j < na)
                self.bs[i, j] = 0.25 if same else 0.5
        # split the cusp slope between the two scales: b1 + b2 = bs, b2 = a_long*bs
        self.a = a_long

    def _terms(self, d, bs):
        b1 = (1.0 - self.a) * bs; b2 = self.a * bs
        e1 = b1 * np.exp(-self.g1 * d); e2 = b2 * np.exp(-self.g2 * d)
        return e1, e2

    def lnJ(self, r):
        Nw, Ne, _ = r.shape; tot = np.zeros(Nw)
        for i in range(Ne):
            for j in range(i + 1, Ne):
                d = np.linalg.norm(r[:, i] - r[:, j], axis=-1)
                e1, e2 = self._terms(d, self.bs[i, j])
                tot += d * (e1 + e2)          # u = r*(b1 e^{-g1 r} + b2 e^{-g2 r})
        return tot

    def grad_lnJ(self, r):
        Nw, Ne, _ = r.shape; g = np.zeros((Nw, Ne, 3))
        for i in range(Ne):
            for j in range(Ne):
                if i == j:
                    continue
                dvec = r[:, i] - r[:, j]; d = np.linalg.norm(dvec, axis=-1)
                d = np.where(d < 1e-12, 1e-12, d)
                e1, e2 = self._terms(d, self.bs[i, j])
                du = e1 * (1 - self.g1 * d) + e2 * (1 - self.g2 * d)   # u'(r)
                g[:, i] += (du / d)[:, None] * dvec
        return g

    def lap_lnJ(self, r):
        Nw, Ne, _ = r.shape; tot = np.zeros(Nw)
        for i in range(Ne):
            for j in range(i + 1, Ne):
                d = np.linalg.norm(r[:, i] - r[:, j], axis=-1)
                d = np.where(d < 1e-12, 1e-12, d)
                e1, e2 = self._terms(d, self.bs[i, j])
                # (u'' + 2u'/r) for each term b*r*e^{-g r}: e*(2/r - 4g + g^2 r)
                lap_u = (e1 * (2.0 / d - 4 * self.g1 + self.g1 ** 2 * d)
                         + e2 * (2.0 / d - 4 * self.g2 + self.g2 ** 2 * d))
                tot += 2.0 * lap_u
        return tot


# ===========================================================================
# Step 5+6: Metropolis VMC with the gradient-form (bounded) kinetic local energy
# ===========================================================================
def psi_ci_val(wf, r):
    """Psi_CI value only (no gradients), for FD Laplacians.  r:(Nw,Ne,3)->(Nw,) complex."""
    R = wf.R; Nw, Ne, _ = r.shape
    Phi = eval_values_batch(wf.prims, r, R)           # (M, Nw, Ne), coords computed once
    MOval = np.einsum('ap,aWe->pWe', wf.T, Phi)
    na, nb = wf.na, wf.nb
    a_e = list(range(na)); b_e = list(range(na, na + nb))
    dA = np.empty((len(wf.apairs), Nw), dtype=complex)
    for k, (p, q) in enumerate(wf.apairs):
        dA[k] = MOval[p, :, a_e[0]] * MOval[q, :, a_e[1]] - MOval[q, :, a_e[0]] * MOval[p, :, a_e[1]]
    dB = np.empty((len(wf.bpairs), Nw), dtype=complex)
    for k, (p, q) in enumerate(wf.bpairs):
        dB[k] = MOval[p, :, b_e[0]] * MOval[q, :, b_e[1]] - MOval[q, :, b_e[0]] * MOval[p, :, b_e[1]]
    return np.einsum('aW,ab,bW->W', dA, wf.CS, dB)


def local_energy_fd(wf, r, jas=None, h=1e-4):
    """STANDARD local energy E_L = V - 1/2 lap(Psi_T)/Psi_T, finite-variance at nodes.
    Psi_T = Psi_CI * J.  lap(Psi_CI)/Psi_CI by FD (smooth); Jastrow analytic (cusp exact,
    the 1/r in -1/2 lap(lnJ) cancels V's 1/r).  Returns (E_L real(Nw,), psi(Nw,))."""
    psi, gpsi = psi_ci(wf, r)
    Nw, Ne, _ = r.shape
    glnpsi = gpsi / psi[:, None, None]                    # (Nw,Ne,3) complex

    # FD Laplacian of Psi_CI / Psi_CI (total, complex)
    lap = np.zeros(Nw, dtype=complex)
    psi0 = psi
    for i in range(Ne):
        for d in range(3):
            rp = r.copy(); rp[:, i, d] += h
            rm = r.copy(); rm[:, i, d] -= h
            lap += (psi_ci_val(wf, rp) + psi_ci_val(wf, rm) - 2.0 * psi0) / (h * h)
    lap_over_psi = lap / psi0                             # complex

    V = potential(wf, r)
    if jas is None:
        EL = V - 0.5 * lap_over_psi
        return EL.real, psi
    # Jastrow analytic 2nd-derivative terms
    glnJ = jas.grad_lnJ(r)                                # (Nw,Ne,3) real
    lapJ = jas.lap_lnJ(r)                                 # (Nw,) real = sum_i lap_i lnJ
    cross = 2.0 * np.sum(glnpsi * glnJ, axis=(1, 2))      # complex
    gJ2 = np.sum(glnJ ** 2, axis=(1, 2))                  # real
    lap_tot = lap_over_psi + cross + lapJ + gJ2
    EL = V - 0.5 * lap_tot
    return EL.real, psi


# analytic local energy is the production estimator
local_energy = local_energy_analytic


def _init_walkers(wf, nwalk, rng):
    """Place na electrons near Li (z=-R/2), nb near H (z=+R/2)."""
    R = wf.R
    r = rng.normal(0, 0.8, size=(nwalk, wf.na + wf.nb, 3))
    for e in range(wf.na):
        r[:, e, 2] -= R / 2.0
    for e in range(wf.nb):
        r[:, wf.na + e, 2] += R / 2.0
    return r


def vmc(wf, nwalk=2000, nsweep=4000, nburn=800, step=0.45, jas=None, seed=0,
        thin=5, adapt=True, verbose=True):
    """Metropolis on |Psi_CI * J|^2 with single-electron moves; accumulate the
    gradient-form local energy.  Returns (E_mean, E_err, accept, E_tot_mean)."""
    rng = np.random.default_rng(seed)
    Ne = wf.na + wf.nb
    r = _init_walkers(wf, nwalk, rng)

    def logp(rr):
        psi = psi_ci_val(wf, rr)
        lp = 2.0 * np.log(np.abs(psi) + 1e-300)
        if jas is not None:
            lp = lp + 2.0 * jas.lnJ(rr)
        return lp, psi

    lp, _ = logp(r)
    acc = 0; ntry = 0
    est = []             # local-energy samples (electronic, no Vnn double? we include Vnn)
    sd = step
    for sweep in range(nsweep):
        for i in range(Ne):
            rprop = r.copy()
            rprop[:, i] += rng.normal(0, sd, size=(nwalk, 3))
            lpp, _ = logp(rprop)
            a = np.log(rng.uniform(size=nwalk)) < (lpp - lp)
            r[a] = rprop[a]; lp[a] = lpp[a]
            acc += int(a.sum()); ntry += nwalk
        if adapt and sweep < nburn and sweep % 50 == 49:
            ar = acc / max(ntry, 1)
            if ar > 0.6: sd *= 1.1
            elif ar < 0.4: sd /= 1.1
            acc = ntry = 0
        if sweep >= nburn and (sweep - nburn) % thin == 0:
            EL, _ = local_energy(wf, r, jas)
            est.append(EL.copy())
    est = np.array(est)                       # (nsamp, nwalk)
    # block over sweeps to reduce autocorrelation in the error bar
    per_sweep = est.mean(axis=1)              # (nsamp,)
    E_mean = per_sweep.mean()
    E_err = per_sweep.std(ddof=1) / np.sqrt(len(per_sweep))
    accept = acc / max(ntry, 1) if ntry else 0.0
    if verbose:
        print(f"    VMC: E={E_mean:.5f} +/- {E_err:.5f}  (nsamp={len(per_sweep)} "
              f"nwalk={nwalk} step={sd:.2f})")
    return E_mean, E_err, accept


def gate6(tag="small", Jb=1, Lb=0, npi=0, core2=None, nwalk=2000, nsweep=3000,
          nburn=600, seed=0):
    """MANDATORY GATE: VMC of Psi_CI WITHOUT Jastrow reproduces E_fci within stats."""
    wf = cached_wavefunction(tag, Jb=Jb, Lb=Lb, npi=npi, core2=core2, verbose=True)
    print("=" * 74)
    print(f"GATE-6: VMC(Psi_CI, J=1) vs FCI  [config {tag}]  E_fci_tot={wf.E_fci+wf.Vnn:.5f}")
    print("=" * 74)
    E, err, acc = vmc(wf, nwalk=nwalk, nsweep=nsweep, nburn=nburn, jas=None, seed=seed)
    dev = (E - (wf.E_fci + wf.Vnn)) / max(err, 1e-9)
    print(f"  VMC E_tot = {E:.5f} +/- {err:.5f}   FCI = {wf.E_fci+wf.Vnn:.5f}   "
          f"dev = {dev:+.1f} sigma   accept={acc:.2f}")
    ok = abs(dev) < 4.0
    print(f"  GATE-6: {'PASS' if ok else 'FAIL'} (VMC reproduces FCI within 4 sigma)")
    return ok


def gate_lap():
    """Analytic orbital Laplacian vs finite difference."""
    import mpmath as mp
    from prolate_mixed_eri import sto_orbital, valence_prolate_orbital, ZC_LI
    from prolate_allelectron_analytic_fci import sto_orbital_B
    from prolate_allelectron_c4 import from_sigma, valence_pi_orbital
    R = 3.015
    orbs = {
        "Li 1s core": from_sigma(sto_orbital(ZC_LI, R, is_core=True)),
        "Li core2 4.5": from_sigma(sto_orbital(mp.mpf('4.5'), R, is_core=True)),
        "H 1s": from_sigma(sto_orbital_B(mp.mpf('1.0'), R)),
        "bond (0,0)": from_sigma(valence_prolate_orbital(0, 0, mp.mpf('1.0'))),
        "bond (2,1)": from_sigma(valence_prolate_orbital(2, 1, mp.mpf('1.0'))),
        "pi +1 (0,0)": valence_pi_orbital(0, 0, mp.mpf('1.0'), +1),
        "pi -1 (1,0)": valence_pi_orbital(1, 0, mp.mpf('1.0'), -1),
    }
    rng = np.random.default_rng(3)
    pts = rng.uniform(-2.0, 2.0, size=(6, 3)); pts[:, 2] += rng.uniform(-1, 1, 6)
    h = 1e-5
    print("=" * 74)
    print("GATE: analytic orbital Laplacian vs finite difference")
    print("=" * 74)
    ok = True
    for name, orb in orbs.items():
        prim = extract_prim(orb)
        val, grad, lap = orbital_vgl(prim, pts, R)
        lfd = np.zeros_like(lap)
        v0, _, _ = orbital_vgl(prim, pts, R)
        for d in range(3):
            dp = np.zeros(3); dp[d] = h
            vp, _, _ = orbital_vgl(prim, pts + dp, R)
            vm, _, _ = orbital_vgl(prim, pts - dp, R)
            lfd += (vp + vm - 2 * v0) / h ** 2
        err = np.max(np.abs(lap - lfd)) / max(np.max(np.abs(lap)), 1e-300)
        ok &= err < 1e-4
        print(f"  {name:14s}: max|lap-FD|/|lap| = {err:.2e}  {'PASS' if err<1e-4 else 'FAIL'}")
    print(f"\n  GATE lap: {'PASS' if ok else 'FAIL'}")
    return ok


def gate_le(tag="small", **kw):
    """Analytic local energy vs FD local energy (must agree; validates the many-body
    analytic Laplacian end-to-end), with and without a Jastrow."""
    wf = cached_wavefunction(tag, **kw)
    rng = np.random.default_rng(5)
    r = _init_walkers(wf, 40, rng)
    # equilibrate a little so configs are physical
    print("=" * 74)
    print(f"GATE: analytic LE vs FD LE  [config {tag}]")
    print("=" * 74)
    for label, jas in [("J=1", None),
                       ("Jastrow g=0.8", Jastrow(0.8, wf.na, wf.nb))]:
        Ea, _ = local_energy_analytic(wf, r, jas)
        Ef, _ = local_energy_fd(wf, r, jas)
        d = np.max(np.abs(Ea - Ef))
        print(f"  {label:16s}: max|E_analytic - E_FD| = {d:.2e}  "
              f"mean|E|={np.mean(np.abs(Ea)):.2f}  {'PASS' if d < 1e-3 else 'FAIL'}")
    return True


def gate_cusp(wf, seed=0):
    """Unit test: the STANDARD local energy with the cusp Jastrow stays FINITE as two
    electrons approach (the -1/2 lap(lnJ) ~ -b/r cancels V's +1/r_ij for b=1/2)."""
    rng = np.random.default_rng(seed)
    jas = Jastrow(gamma=0.8, na=wf.na, nb=wf.nb)
    print("=" * 74)
    print("GATE cusp: E_L finite as an antiparallel pair (e0 up, e2 down) coalesces")
    print("=" * 74)
    base = _init_walkers(wf, 1, rng)
    ok = True
    for dr in (0.5, 0.1, 0.02, 0.005):
        r = base.copy()
        r[:, 2] = r[:, 0] + np.array([dr, 0, 0])      # put e2 (down) near e0 (up)
        EL, _ = local_energy(wf, r, jas)
        finite = np.all(np.isfinite(EL))
        ok &= finite
        print(f"  r_02={dr:.3f}: E_L={EL[0]:+.4f}  finite={finite}")
    print(f"  GATE cusp: {'PASS (bounded)' if ok else 'FAIL (blows up -> cusp mismatch)'}")
    return ok


def run_full(tag="m16", Jb=2, Lb=1, npi=1, Jpi=1, Lpi=0, core2=(4.5, 1.6),
             gammas=(0.3, 0.5, 0.8, 1.2), nwalk=1500, nsweep=3000, nburn=600,
             seed=0):
    """Full pipeline: build/cache the wavefunction, GATE-6 (J=1 vs FCI), cusp gate,
    then a Jastrow gamma-scan for the correlated energy."""
    wf = cached_wavefunction(tag, Jb=Jb, Lb=Lb, npi=npi, Jpi=Jpi, Lpi=Lpi,
                             core2=core2, verbose=True)
    Efci = wf.E_fci + wf.Vnn
    print("=" * 74)
    print(f"CONFIG {tag}: M={len(wf.prims)} Mk={wf.Mk} nd={len(wf.c)}  FCI E_tot={Efci:.5f}")
    print(f"  anchors: beat {Efci:.5f} (this basis), land near -8.062, stay > -8.0705")
    print("=" * 74)
    print("\n[GATE-6] VMC(Psi_CI, J=1) vs FCI:")
    E0, e0, a0 = vmc(wf, nwalk=nwalk, nsweep=nsweep, nburn=nburn, jas=None, seed=seed)
    dev = (E0 - Efci) / max(e0, 1e-9)
    print(f"  E={E0:.5f}+/-{e0:.5f}  FCI={Efci:.5f}  dev={dev:+.1f}sig  "
          f"{'PASS' if abs(dev) < 4 else 'FAIL'}")
    print("\n[cusp gate]")
    gate_cusp(wf)
    print("\n[Jastrow gamma-scan]  E_VMC[Psi_CI * exp(sum u)], u=b r e^{-g r}:")
    results = []
    for g in gammas:
        jas = Jastrow(gamma=g, na=wf.na, nb=wf.nb)
        E, err, acc = vmc(wf, nwalk=nwalk, nsweep=nsweep, nburn=nburn, jas=jas,
                          seed=seed, verbose=False)
        dcusp = (E - Efci) * 1e3
        results.append((g, E, err))
        flag = "  <-- BELOW EXACT (BUG)" if E < -8.0705 else ""
        print(f"  gamma={g:.2f}: E={E:.5f}+/-{err:.5f}  dE(vs orbital)={dcusp:+.1f}mHa"
              f"  vs exact -8.0705={(E+8.0705)*1e3:+.1f}mHa{flag}", flush=True)
    best = min(results, key=lambda t: t[1])
    # parabolic fit for the optimal gamma (if the minimum is interior)
    gs = np.array([r[0] for r in results]); Es = np.array([r[1] for r in results])
    g_opt = best[0]
    if 0 < np.argmin(Es) < len(Es) - 1:
        cpar = np.polyfit(gs, Es, 2)
        if cpar[0] > 0:
            g_opt = -cpar[1] / (2 * cpar[0])
    print(f"\n  best grid: gamma={best[0]:.2f}  E={best[1]:.5f}+/-{best[2]:.5f}  "
          f"({(best[1]+8.0705)*1e3:+.1f} mHa from exact)  | parabola g_opt~{g_opt:.2f}")
    return results, g_opt, wf


def gscan(tag="m16", gammas=(0.5, 0.8, 1.2, 1.8, 2.5), nwalk=800, nsweep=1500,
          nburn=350, seed=7, jclass=Jastrow, **kw):
    """Lean Jastrow gamma-scan (gate-6 already validated separately)."""
    wf = cached_wavefunction(tag, **kw)
    Efci = wf.E_fci + wf.Vnn
    print("=" * 74)
    print(f"GAMMA-SCAN  {tag}  FCI={Efci:.5f}  exact=-8.0705  ({jclass.__name__})")
    print("=" * 74, flush=True)
    res = []
    for g in gammas:
        jas = jclass(gamma=g, na=wf.na, nb=wf.nb)
        E, err, acc = vmc(wf, nwalk=nwalk, nsweep=nsweep, nburn=nburn, jas=jas,
                          seed=seed, verbose=False)
        flag = "  <-- BELOW EXACT (BUG)" if E < -8.0705 else ""
        print(f"  g={g:.2f}: E={E:.5f}+/-{err:.5f}  dE(cusp)={(E-Efci)*1e3:+.1f}mHa  "
              f"exact-gap={(E+8.0705)*1e3:+.1f}mHa  acc={acc:.2f}{flag}", flush=True)
        res.append((g, E, err))
    best = min(res, key=lambda t: t[1])
    print(f"\n  best: g={best[0]:.2f}  E={best[1]:.5f}+/-{best[2]:.5f}  "
          f"({(best[1]+8.0705)*1e3:+.1f} mHa from exact)")
    return res


def scan_list(tag, jlist, nwalk=800, nsweep=1600, nburn=350, seed=7, **kw):
    """Scan a list of (label, jastrow) pairs; jastrow=None means J=1 (gate-6)."""
    wf = cached_wavefunction(tag, **kw)
    Efci = wf.E_fci + wf.Vnn
    print("=" * 78)
    print(f"JASTROW SCAN  {tag}  FCI={Efci:.5f}  exact=-8.0705")
    print("=" * 78, flush=True)
    res = []
    for label, jas in jlist:
        E, err, acc = vmc(wf, nwalk=nwalk, nsweep=nsweep, nburn=nburn, jas=jas,
                          seed=seed, verbose=False)
        flag = "  <-- BELOW EXACT (BUG)" if E < -8.0705 else ""
        print(f"  {label:28s}: E={E:.5f}+/-{err:.5f}  dE={(E-Efci)*1e3:+.1f}mHa  "
              f"exact-gap={(E+8.0705)*1e3:+.1f}mHa{flag}", flush=True)
        res.append((label, E, err))
    best = min(res, key=lambda t: t[1])
    print(f"\n  best: {best[0]}  E={best[1]:.5f}+/-{best[2]:.5f}  "
          f"({(best[1]+8.0705)*1e3:+.1f} mHa from exact)")
    return res


def pade3_scan(tag="m16", gamma=1.4, combos=None, nwalk=800, nsweep=1600,
               nburn=350, seed=7, **kw):
    """Scan the 2-body Pade (fixed gamma) + 3-body e-e-n Gaussian (c,beta,delta)."""
    wf = cached_wavefunction(tag, **kw)
    R = wf.R
    nuclei = [np.array([0.0, 0.0, -R / 2.0]), np.array([0.0, 0.0, R / 2.0])]
    na, nb = wf.na, wf.nb
    if combos is None:
        # (c, beta, delta)
        combos = [
            (0.0, 0.0, 0.0),                 # ref = pure 2-body Pade
            (0.5, 3.0, 4.0),
            (1.0, 3.0, 4.0),
            (1.5, 3.0, 4.0),
            (1.0, 2.0, 3.0),
            (1.0, 4.0, 6.0),
        ]
    jlist = []
    for (c, beta, delta) in combos:
        if c == 0.0:
            jlist.append((f"Pade g={gamma} (ref)", PadeJastrow(gamma, na, nb)))
        else:
            jlist.append((f"+3body c={c} b={beta} d={delta}",
                          PadePlus3(gamma, c, beta, delta, na, nb, nuclei)))
    return scan_list(tag, jlist, nwalk=nwalk, nsweep=nsweep, nburn=nburn, seed=seed, **kw)


def two_pade_scan(tag="m16", na=2, nb=2, **kw):
    """Scan the two-scale Pade (g1 core, g2 valence, f split)."""
    combos = [
        ("2Pade g1=2.5 g2=0.8 f=0.6", (2.5, 0.8, 0.6)),
        ("2Pade g1=3.0 g2=0.8 f=0.6", (3.0, 0.8, 0.6)),
        ("2Pade g1=3.0 g2=0.8 f=0.7", (3.0, 0.8, 0.7)),
        ("2Pade g1=2.5 g2=1.0 f=0.5", (2.5, 1.0, 0.5)),
        ("2Pade g1=3.5 g2=0.6 f=0.7", (3.5, 0.6, 0.7)),
        ("1Pade g=1.5 (ref)",         None),
    ]
    jlist = []
    for label, p in combos:
        if p is None:
            jlist.append(("1Pade g=1.5 (ref)", PadeJastrow(1.5, na, nb)))
        else:
            g1, g2, f = p
            jlist.append((label, TwoPadeJastrow(g1, g2, f, na, nb)))
    return scan_list(tag, jlist, **kw)


def final_run(tag="m16", gamma=1.5, nwalk=4000, nsweep=8000, nburn=1500,
              seeds=(10, 11, 12), jclass=PadeJastrow, **kw):
    """High-statistics final energy at a fixed gamma, averaged over independent seeds."""
    wf = cached_wavefunction(tag, **kw)
    Efci = wf.E_fci + wf.Vnn
    jas = jclass(gamma=gamma, na=wf.na, nb=wf.nb)
    print("=" * 74)
    print(f"FINAL RUN  config {tag}  gamma={gamma}  (FCI {Efci:.5f}, exact -8.0705)")
    print("=" * 74)
    Es, errs = [], []
    for s in seeds:
        E, err, acc = vmc(wf, nwalk=nwalk, nsweep=nsweep, nburn=nburn, jas=jas,
                          seed=s, verbose=False)
        Es.append(E); errs.append(err)
        print(f"  seed {s}: E={E:.5f} +/- {err:.5f}  acc={acc:.2f}", flush=True)
    Es = np.array(Es); errs = np.array(errs)
    Emean = Es.mean(); Eerr = (Es.std(ddof=1) / np.sqrt(len(Es))
                               if len(Es) > 1 else errs[0])
    print(f"\n  FINAL E_tot = {Emean:.5f} +/- {Eerr:.5f} Ha")
    print(f"    vs orbital ceiling {Efci:.5f}: {(Emean-Efci)*1e3:+.1f} mHa")
    print(f"    vs exact -8.0705:            {(Emean+8.0705)*1e3:+.1f} mHa")
    variational = Emean > -8.0705
    print(f"    variational (E > -8.0705): {variational}"
          f"{'' if variational else '  <-- BUG'}")
    return Emean, Eerr


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "gate_orb"
    if arg == "gate_orb":
        gate_orb()
    elif arg == "gate_val":
        gate_value_vs_grid()
    elif arg == "gate2":
        a = gate_orb(); print(); b = gate_value_vs_grid()
        print(f"\nSTEP-2 OVERALL: {'PASS' if (a and b) else 'FAIL'}")
    elif arg == "buildwf":
        wf = build_lih_wavefunction(Jb=1, Lb=0, npi=0, core2=None, verbose=True)
        print(f"  Rayleigh gate: E_elec={wf.E_fci:.6f}  sum c^2={np.sum(wf.c**2):.6f}")
    elif arg == "gate6":
        gate6()
    elif arg == "gate_lap":
        gate_lap()
    elif arg == "gate_le":
        gate_lap(); print(); gate_le(tag="small", Jb=1, Lb=0, npi=0, core2=None)
    elif arg == "resave":
        # convert legacy class-pickle caches (built as __main__) to module-independent dicts
        import pickle, glob
        DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        for path in glob.glob(os.path.join(DATA, "lih_vmc_wf_*.pkl")):
            with open(path, "rb") as f:
                obj = pickle.load(f)
            if not isinstance(obj, dict):
                with open(path, "wb") as f:
                    pickle.dump(_wf_to_dict(obj), f)
                print(f"  resaved {os.path.basename(path)} (class -> dict)")
            else:
                print(f"  {os.path.basename(path)} already dict")
    elif arg == "buildm16":
        # build+cache the full M=16 wavefunction (the -8.029 basis)
        wf = cached_wavefunction("m16", Jb=2, Lb=1, npi=1, Jpi=1, Lpi=0,
                                 core2=(4.5, 1.6), verbose=True)
        print(f"  cached m16: Mk={wf.Mk} nd={len(wf.c)} E_tot={wf.E_fci+wf.Vnn:.5f}")
    elif arg == "full":
        run_full()
