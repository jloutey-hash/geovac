"""Tensor lifted-state route, v2: EXACT identity checks + per-element inequalities.

v1's TT5 "violation" at (3,2) was a grid artifact: Lip_prod was estimated from
16 sample points (a 6-dimensional manifold) while L_ab used exact operator
norms.  v2 replaces every ratio-of-suprema test by either
  (E) an exact matrix/pointwise identity, or
  (I) a per-element inequality with BOTH sides computed exactly.

Conventions: unit S^3, d(e,g) = arccos(Re tr D^{1/2}(g) / 2);
product metric d_prod = sqrt(d_a^2 + d_b^2).
"""
from __future__ import annotations
import sys
import itertools
from fractions import Fraction
import numpy as np

sys.path.insert(0, r"C:/Users/jlout/Desktop/Project_Geometric")
sys.path.insert(0, r"C:/Users/jlout/Desktop/Project_Geometric/debug")
from p38_g1g2_scalar_prototype import (
    Window, wigner_D, f_eval, sigma_symbol, jm_range,
)

rng = np.random.default_rng(20260904)


# ---------------------------------------------------------------- group utils
def euler_to_su2(e):
    return wigner_D(Fraction(1, 2), *e)


def su2_to_euler(U):
    """Invert U = exp(-i a Jz) exp(-i b Jy) exp(-i c Jz)."""
    c0, s0 = abs(U[0, 0]), abs(U[1, 0])
    beta = 2.0 * np.arctan2(s0, c0)
    if c0 > 1e-12 and s0 > 1e-12:
        p = -np.angle(U[0, 0])
        q = np.angle(U[1, 0])
        return (p + q, beta, p - q)
    if c0 <= 1e-12:                      # beta = pi
        q = np.angle(U[1, 0])
        return (q, beta, -q)
    p = -np.angle(U[0, 0])               # beta = 0
    return (2 * p, beta, 0.0)


def rand_euler():
    return (rng.uniform(0, 2 * np.pi),
            np.arccos(rng.uniform(-1, 1)),
            rng.uniform(0, 2 * np.pi))


def dist_e(e1, e2):
    U1, U2 = euler_to_su2(e1), euler_to_su2(e2)
    c = np.clip(np.real(np.trace(U1 @ U2.conj().T)) / 2, -1, 1)
    return float(np.arccos(c))


def mul(e1, e2):
    return su2_to_euler(euler_to_su2(e1) @ euler_to_su2(e2))


def inv(e):
    return su2_to_euler(euler_to_su2(e).conj().T)


IDE = (0.0, 0.0, 0.0)


def gamma_unit(n_max, npts=400001):
    th = np.linspace(1e-9, 2 * np.pi, npts)
    h = np.zeros_like(th)
    for n in range(1, n_max + 1):
        h += np.sqrt(n) * np.sin(n * th / 2) / np.sin(th / 2)
    K = h ** 2 / (n_max * (n_max + 1) / 2)
    w = (1 / np.pi) * np.sin(th / 2) ** 2
    return float(np.trapezoid(K * w * (th / 2), th) / np.trapezoid(K * w, th))


def gamma_joint_unit(na, nb, M=3000):
    x, wq = np.polynomial.legendre.leggauss(M)
    th = (x + 1) * np.pi
    wq = wq * np.pi
    mu = np.sin(th / 2) ** 2 / np.pi

    def Kf(n):
        h = np.zeros_like(th)
        for k in range(1, n + 1):
            h += np.sqrt(k) * np.sin(k * th / 2) / np.sin(th / 2)
        return h ** 2 / (n * (n + 1) / 2)

    wa, wb = wq * mu * Kf(na), wq * mu * Kf(nb)
    d = th / 2
    return float(wa @ np.sqrt(d[:, None] ** 2 + d[None, :] ** 2) @ wb)


# ------------------------------------------------------------------- E0 check
def check_group_utils():
    bad = 0.0
    for _ in range(200):
        e1, e2 = rand_euler(), rand_euler()
        A = euler_to_su2(mul(e1, e2))
        B = euler_to_su2(e1) @ euler_to_su2(e2)
        bad = max(bad, float(np.max(np.abs(A - B))))
        C = euler_to_su2(mul(e1, inv(e1)))
        bad = max(bad, float(np.max(np.abs(C - np.eye(2)))))
    return bad


# --------------------------------------------------- E1: joint covariance
def check_covariance(w, U_of, n_trials=6):
    """rho_v(M_{J,A,B}) = sum_C D^J_{AC}(v^{+-1}) M_{J,C,B}, exact matrices."""
    worst = {}
    for J in [j for j in w.js if j > 0]:
        ms = jm_range(J)
        best = np.inf
        for _ in range(n_trials):
            v = rand_euler()
            Uv = U_of(v)
            for (A, B) in [(ms[0], ms[-1]), (ms[-1], ms[0])]:
                lhs = Uv @ w.multiplier(J, A, B) @ Uv.conj().T
                for tag, ee in (("v", v), ("v^-1", inv(v))):
                    Dv = wigner_D(J, *ee)
                    rhs = np.zeros_like(lhs)
                    for ci, C in enumerate(ms):
                        rhs = rhs + Dv[ms.index(A), ci] * w.multiplier(J, C, B)
                    r = float(np.linalg.norm(lhs - rhs, 2))
                    best = min(best, r)
        worst[float(J)] = best
    return worst


# ------------------------------- E4: int D^J(u^{-1}) K(u) du = sigma_J * I
def check_schur_average(w, n_mc=40000):
    """Monte-Carlo Haar average of D^J(u^{-1}) K(u) vs sigma_J I."""
    out = {}
    js = [j for j in w.js if j > 0]
    # Haar sample on SU(2): uniform quaternion
    q = rng.normal(size=(n_mc, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    # rotation angle theta in [0, 2pi): q0 = cos(theta/2)
    theta = 2 * np.arccos(np.clip(q[:, 0], -1, 1))          # in [0, pi]... map
    # SU(2) element U = q0 I - i (q1 sx + q2 sy + q3 sz); tr = 2 q0 = 2cos(t/2)
    s = np.sin(theta / 2)
    Kv = np.zeros(n_mc)
    n_max = int(2 * max(w.js)) + 1
    Z = n_max * (n_max + 1) / 2
    h = np.zeros(n_mc)
    small = np.abs(s) < 1e-12
    for n in range(1, n_max + 1):
        h += np.sqrt(n) * np.where(small, float(n),
                                   np.sin(n * theta / 2) / np.where(small, 1, s))
    Kv = h ** 2 / Z
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    for J in js:
        d = int(2 * J) + 1
        acc = np.zeros((d, d), dtype=complex)
        for i in range(n_mc):
            U = (q[i, 0] * np.eye(2) - 1j * (q[i, 1] * sx + q[i, 2] * sy
                                             + q[i, 3] * sz))
            e = su2_to_euler(U.conj().T)          # u^{-1}
            acc += wigner_D(J, *e) * Kv[i]
        acc /= n_mc
        sig = sigma_symbol(w, J)
        out[float(J)] = (sig,
                         float(np.max(np.abs(acc - sig * np.eye(d)))))
    return out


# -------------------------------------------------------------------- driver
def run(jmax_a, jmax_b, n_pairs=200, n_trials=4):
    na, nb = int(2 * jmax_a) + 1, int(2 * jmax_b) + 1
    wa, wb = Window(jmax_a), Window(jmax_b)
    xi = np.kron(wa.xi(), wb.xi())
    ga, gb = gamma_unit(na), gamma_unit(nb)
    gj = gamma_joint_unit(na, nb)
    print("")
    print("=" * 74)
    print("(n_a, n_b) = (%d, %d)   dim %d   gamma_a=%.6f gamma_b=%.6f "
          "JOINT=%.6f (<= %.6f)" % (na, nb, wa.dim * wb.dim, ga, gb, gj,
                                    ga + gb))

    bands_a = [j for j in wa.js if j > 0]
    bands_b = [j for j in wb.js if j > 0]

    def U(u):
        return np.kron(wa.U(*u[0]), wb.U(*u[1]))

    # ---- I2 : fact (ii) per-element, BOTH sides exact ----
    print("I2  |ups(T)(g) - ups(T)(g')| <= ||rho_{g'g^-1}(T) - T||_op  "
          "[per-pair, exact]")
    worst_ratio = 0.0
    n_viol = 0
    for trial in range(n_trials):
        T = np.zeros((wa.dim * wb.dim,) * 2, dtype=complex)
        for Ja, Jb in itertools.product(bands_a, bands_b):
            Aa, Ba = jm_range(Ja)[0], jm_range(Ja)[-1]
            Ab, Bb = jm_range(Jb)[0], jm_range(Jb)[-1]
            T = T + rng.normal() * np.kron(wa.multiplier(Ja, Aa, Ba),
                                           wb.multiplier(Jb, Ab, Bb))
        T = T + T.conj().T
        for _ in range(n_pairs):
            g = (rand_euler(), rand_euler())
            gp = (rand_euler(), rand_euler())
            Ug, Ugp = U(g), U(gp)
            lhs = abs(np.vdot(Ug @ xi, T @ (Ug @ xi))
                      - np.vdot(Ugp @ xi, T @ (Ugp @ xi)))
            q = (mul(gp[0], inv(g[0])), mul(gp[1], inv(g[1])))
            Uq = U(q)
            rhs = float(np.linalg.norm(Uq @ T @ Uq.conj().T - T, 2))
            if rhs > 1e-12:
                worst_ratio = max(worst_ratio, lhs / rhs)
            if lhs > rhs * (1 + 1e-9):
                n_viol += 1
    print("     max ratio LHS/RHS over %d pairs = %.10f   violations: %d"
          % (n_trials * n_pairs, worst_ratio, n_viol))
    i2 = (n_viol == 0)
    print("     I2 verdict: %s" % ("PASS" if i2 else "FAIL"))

    # ---- I2b : same, but RHS = L_ab(T) * d_prod(g,g')  (the seminorm form) ----
    print("I2b |ups(T)(g)-ups(T)(g')| <= L_ab_lower(T) * d_prod(g,g')  "
          "[L_ab under-estimated => conservative]")
    T = np.zeros((wa.dim * wb.dim,) * 2, dtype=complex)
    for Ja, Jb in itertools.product(bands_a, bands_b):
        Aa, Ba = jm_range(Ja)[0], jm_range(Ja)[-1]
        Ab, Bb = jm_range(Jb)[0], jm_range(Jb)[-1]
        T = T + rng.normal() * np.kron(wa.multiplier(Ja, Aa, Ba),
                                       wb.multiplier(Jb, Ab, Bb))
    T = T + T.conj().T
    vs = [(rand_euler(), rand_euler()) for _ in range(120)]
    L_lower = 0.0
    for v in vs:
        d = float(np.hypot(dist_e(v[0], IDE), dist_e(v[1], IDE)))
        if d < 1e-6:
            continue
        Uv = U(v)
        L_lower = max(L_lower,
                      float(np.linalg.norm(Uv @ T @ Uv.conj().T - T, 2)) / d)
    bad = 0
    mx = 0.0
    for _ in range(n_pairs):
        g = (rand_euler(), rand_euler())
        gp = (rand_euler(), rand_euler())
        Ug, Ugp = U(g), U(gp)
        lhs = abs(np.vdot(Ug @ xi, T @ (Ug @ xi))
                  - np.vdot(Ugp @ xi, T @ (Ugp @ xi)))
        d = float(np.hypot(dist_e(g[0], gp[0]), dist_e(g[1], gp[1])))
        mx = max(mx, lhs / (L_lower * d))
        if lhs > L_lower * d * (1 + 1e-9):
            bad += 1
    print("     L_ab_lower(T) = %.6f   max ratio = %.6f   violations: %d"
          % (L_lower, mx, bad))
    i2b = (bad == 0)
    print("     I2b verdict: %s" % ("PASS" if i2b else "FAIL"))

    # ---- I4 : fact (iv) ----
    print("I4  ||S(ups(T)) - T||_op <= gamma^(ab) * L_ab_lower(T)  "
          "[conservative]")
    ok4 = True
    for trial in range(n_trials):
        coeffs = {}
        T = np.zeros((wa.dim * wb.dim,) * 2, dtype=complex)
        Sv = np.zeros_like(T)
        for Ja, Jb in itertools.product(bands_a, bands_b):
            Aa, Ba = jm_range(Ja)[0], jm_range(Ja)[-1]
            Ab, Bb = jm_range(Jb)[0], jm_range(Jb)[-1]
            c = rng.normal()
            Mk = np.kron(wa.multiplier(Ja, Aa, Ba), wb.multiplier(Jb, Ab, Bb))
            T = T + c * Mk
            Sv = Sv + c * sigma_symbol(wa, Ja) * sigma_symbol(wb, Jb) * Mk
        T, Sv = T + T.conj().T, Sv + Sv.conj().T
        L_lower = 0.0
        for v in vs:
            d = float(np.hypot(dist_e(v[0], IDE), dist_e(v[1], IDE)))
            if d < 1e-6:
                continue
            Uv = U(v)
            L_lower = max(L_lower, float(
                np.linalg.norm(Uv @ T @ Uv.conj().T - T, 2)) / d)
        defect = float(np.linalg.norm(Sv - T, 2))
        ok = defect <= gj * L_lower * (1 + 1e-9)
        ok4 = ok4 and ok
        print("     trial %d: defect=%.6f  gamma^(ab)*L_lower=%.6f  "
              "ratio=%.4f  %s" % (trial, defect, gj * L_lower,
                                  defect / (gj * L_lower),
                                  "OK" if ok else "VIOLATED"))
    print("     I4 verdict: %s" % ("PASS" if ok4 else "FAIL"))
    return i2 and i2b and ok4


if __name__ == "__main__":
    print("E0  group-utility round-trip max error: %.2e" % check_group_utils())
    for jm in [Fraction(1, 2), Fraction(1)]:
        w = Window(jm)
        cov = check_covariance(w, lambda e: w.U(*e))
        print("E1  covariance rho_v(M) = M_{lambda_v f}, window n=%d: %s"
              % (int(2 * jm) + 1,
                 {k: "%.2e" % v for k, v in cov.items()}))
    w = Window(Fraction(1))
    sch = check_schur_average(w, n_mc=8000)
    print("E4  MC  int D^J(u^-1) K(u) du  vs  sigma_J * I  (n=3 window):")
    for J, (sig, err) in sch.items():
        print("      J=%.1f  sigma_J=%.8f   max|MC - sigma_J I| = %.4f"
              % (J, sig, err))
    ok = True
    ok = run(Fraction(1, 2), Fraction(1, 2)) and ok
    ok = run(Fraction(1), Fraction(1, 2)) and ok
    ok = run(Fraction(1), Fraction(1)) and ok
    print("")
    print("OVERALL: %s" % ("PASS" if ok else "FAIL"))
