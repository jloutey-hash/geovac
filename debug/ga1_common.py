"""GA1 (docs/predictions/2026-10-22_geminal_shapes.md): shared helpers for the cusp-correct geminal
f = r exp(-gamma r) on both edges of the correlated triangle (apex placement, fingerprint fibre).

Derivation and constants: debug/ga1_0_kernel.py (sympy; the exp case reproduces ring_ext_QG_common's
kernel exactly).  Nothing in ring_ext_QG_common.py is modified; its node tables are reused.

Shapes: 'exp' f~ = 8 pi g/(k^2+g^2)^2 = -d/dg [4 pi/(k^2+g^2)];  'lin' f~ = 8 pi (3 g^2 - k^2)/(k^2+g^2)^3
= d^2/dg^2 [4 pi/(k^2+g^2)].  Both are gamma-derivatives of the Yukawa line, taken analytically.

Loop (route A, QG1 style): I_111(m1, m2, m3; s12, s13, s23) = pi^2 arctan(sqrt(det M)/D)/sqrt(det M) in
closed form; the lin loop int d^3 l g(l-x1) g(l-x2)/|l-x3|^2 with g(q) = (3 g^2 - q^2)/(q^2+g^2)^3
= (1/2) d^2/dg^2 (q^2+g^2)^-1 is loop_lin = (1/4) d^2/dg1^2 d^2/dg2^2 I_111(g1, g2, 0)
= (d_A + 2 A d_A^2)(d_B + 2 B d_B^2) I_111, A = g1^2, B = g2^2, by a truncated bivariate Taylor algebra
(degree <= 2 in each nilpotent), never by finite differences.

Kite (route B): Phi~(w) = int_0^1 int_0^1 du dv sqrt(u(1-u)) M^-3 sum_n q_n J_n(s, k2),
M = sqrt(1 + u(1-u) z), N = sqrt(g^2 + v(1-v) z), s = M + N, k2 = (u-v)^2 z, z = w^2, V = v(1-v),
q_n = Q_n + M Q_{n-1}, Q(r) = sum_n Q_n r^n:
 exp: Q0 = 6 g^2 V/N^5, Q1 = 6 g^2 V/N^4, Q2 = 2 g^2 V/N^3
 lin: Q0 = 6V/N^5 - 30 g^2 V/N^7 + 210 g^4 V^2/N^9, Q1 = 6V/N^4 - 30 g^2 V/N^6 + 210 g^4 V^2/N^8,
      Q2 = 2V/N^3 - 12 g^2 V/N^5 + 90 g^4 V^2/N^7, Q3 = -2 g^2 V/N^4 + 20 g^4 V^2/N^6, Q4 = 2 g^4 V^2/N^5
 J_0 = 1/D, J_1 = 2s/D^2, J_2 = 2(3s^2-k2)/D^3, J_3 = 24 s(s^2-k2)/D^4, J_4 = 24(k2^2 - 10 k2 s^2 + 5 s^4)/D^5,
 J_5 = 240 s (s^2 - 3 k2)(3 s^2 - k2)/D^6,  D = s^2 + k2.
Then Y(t) = (1/(2 pi^2)) int_0^oo w^2 Phi~(w) t^3 (w^2+t^2)^(-3/2) dw as in QG."""
from __future__ import annotations
import sys
from math import factorial
from pathlib import Path
import mpmath as mp
from flint import arb, acb
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm


# ---------------------------------------------------------------- truncated bivariate Taylor algebra
class TD2:
    """v + sum_{(i,j) != (0,0), i,j <= 2} c_ij e1^i e2^j, with e1^3 = e2^3 = 0.  c[i][j]."""
    __slots__ = ("c",)

    def __init__(self, v=None, c=None):
        if c is not None:
            self.c = c
        else:
            self.c = [[v, 0, 0], [0, 0, 0], [0, 0, 0]]

    @staticmethod
    def var(v, which):
        """v + e1 (which = 0) or v + e2 (which = 1)"""
        t = TD2(v)
        if which == 0:
            t.c[1][0] = 1
        else:
            t.c[0][1] = 1
        return t

    @staticmethod
    def lift(x):
        return x if isinstance(x, TD2) else TD2(x)

    @property
    def v(self):
        return self.c[0][0]

    def __add__(self, o):
        o = TD2.lift(o)
        return TD2(c=[[self.c[i][j] + o.c[i][j] for j in range(3)] for i in range(3)])
    __radd__ = __add__

    def __neg__(self):
        return TD2(c=[[-self.c[i][j] for j in range(3)] for i in range(3)])

    def __sub__(self, o):
        return self + (-TD2.lift(o))

    def __rsub__(self, o):
        return TD2.lift(o) + (-self)

    def __mul__(self, o):
        if not isinstance(o, TD2):
            return TD2(c=[[self.c[i][j] * o for j in range(3)] for i in range(3)])
        a, b = self.c, o.c
        out = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        for i in range(3):
            for j in range(3):
                acc = 0
                for k in range(i + 1):
                    for l in range(j + 1):
                        if a[k][l] != 0 and b[i - k][j - l] != 0:
                            acc += a[k][l] * b[i - k][j - l]
                out[i][j] = acc
        return TD2(c=out)
    __rmul__ = __mul__

    def fn(self, derivs):
        """f(self) from f, f', ..., f'''' at self.v (list of 5)"""
        d = TD2(c=[[self.c[i][j] if (i, j) != (0, 0) else 0 for j in range(3)] for i in range(3)])
        out = TD2(derivs[0])
        pw = TD2(1)
        for n in range(1, 5):
            pw = pw * d
            out = out + pw * (derivs[n] / factorial(n))
        return out

    def inv(self):
        x = self.v
        return self.fn([1 / x, -1 / x**2, 2 / x**3, -6 / x**4, 24 / x**5])

    def __truediv__(self, o):
        return self * TD2.lift(o).inv()

    def __rtruediv__(self, o):
        return TD2.lift(o) * self.inv()

    def sqrt(self):
        x = self.v
        s = mp.sqrt(x)
        return self.fn([s, s / (2 * x), -s / (4 * x**2), 3 * s / (8 * x**3), -15 * s / (16 * x**4)])

    def d(self, i, j):
        """d^i/dA^i d^j/dB^j (unmixed Taylor coefficient times i! j!)"""
        return self.c[i][j] * factorial(i) * factorial(j)


def A_fun4(y):
    """A(y) = arctan(sqrt y)/sqrt y and derivatives up to order 4 (analytic in y > -1).
    2 y A' + A = 1/(1+y)  =>  A^(n+1) = [(-1)^n n!/(1+y)^(n+1) - (2n+1) A^(n)] / (2y)."""
    if abs(y) < mp.mpf("0.1"):
        out = [mp.mpf(0)] * 5
        k = 0
        while True:
            base = (-1) ** k / mp.mpf(2 * k + 1)
            term = base * y**k
            for n in range(5):
                if k >= n:
                    out[n] += base * (factorial(k) // factorial(k - n)) * y ** (k - n)
            if k > 6 and abs(term) < mp.eps * mp.mpf(10) ** -5:
                break
            k += 1
        return out
    if y > 0:
        s = mp.sqrt(y)
        a0 = mp.atan(s) / s
    else:
        s = mp.sqrt(-y)
        a0 = mp.atanh(s) / s
    out = [a0]
    for n in range(4):
        out.append(((-1) ** n * factorial(n) / (1 + y) ** (n + 1) - (2 * n + 1) * out[n]) / (2 * y))
    return out


def tri_sq_TD2(a, b, c, s12, s13, s23):
    """I_111 with squared masses a, b (TD2 or numbers), c a number (the massless/third line)."""
    A_, B_ = TD2.lift(a), TD2.lift(b)
    m1, m2 = A_.sqrt(), B_.sqrt()
    m3 = mp.sqrt(c)
    M12 = (A_ + B_ + s12) * mp.mpf(0.5)
    M13 = (A_ + c + s13) * mp.mpf(0.5)
    M23 = (B_ + c + s23) * mp.mpf(0.5)
    det = A_ * (B_ * c - M23 * M23) - M12 * (M12 * c - M23 * M13) + M13 * (M12 * M23 - B_ * M13)
    D = m1 * m2 * m3 + M12 * m3 + M13 * m2 + M23 * m1
    y = det / (D * D)
    Ay = y.fn(A_fun4(y.v))
    return (Ay / D) * mp.pi ** 2


def loop_lin(m1, m2, m3, s12, s13, s23):
    """int d^3 l g_1(l-x1) g_2(l-x2) / (|l-x3|^2 + m3^2), g_i(q) = (3 m_i^2 - q^2)/(q^2+m_i^2)^3
    = (1/4) d^2/dm1^2 d^2/dm2^2 I_111 = (d_A + 2A d_A^2)(d_B + 2B d_B^2) I_111."""
    A = mp.mpf(m1) ** 2
    B = mp.mpf(m2) ** 2
    I = tri_sq_TD2(TD2.var(A, 0), TD2.var(B, 1), mp.mpf(m3) ** 2, s12, s13, s23)
    return I.d(1, 1) + 2 * B * I.d(1, 2) + 2 * A * I.d(2, 1) + 4 * A * B * I.d(2, 2)


def loop_exp_TD2(m1, m2, m3, s12, s13, s23):
    """I_221 = d_A d_B I_111 through the TD2 algebra (consistency check against Cm.I221)."""
    A = mp.mpf(m1) ** 2
    B = mp.mpf(m2) ** 2
    I = tri_sq_TD2(TD2.var(A, 0), TD2.var(B, 1), mp.mpf(m3) ** 2, s12, s13, s23)
    return I.d(1, 1)


# ---------------------------------------------------------------- route B: the kite
def kernel(M, N, V, g2, g4, k2, shape):
    """sum_n q_n J_n at one (u, v) node; M, N, k2 acb; V, g2 = g^2, g4 = g^4 arb."""
    N2 = N * N
    iN = 1 / N
    iN2 = iN * iN
    iN3 = iN2 * iN
    iN4 = iN2 * iN2
    iN5 = iN4 * iN
    if shape == "exp":
        Q = [6 * g2 * V * iN5, 6 * g2 * V * iN4, 2 * g2 * V * iN3]
    else:
        iN6 = iN3 * iN3
        iN7 = iN6 * iN
        iN8 = iN4 * iN4
        iN9 = iN8 * iN
        V2 = V * V
        Q = [6 * V * iN5 - 30 * g2 * V * iN7 + 210 * g4 * V2 * iN9,
             6 * V * iN4 - 30 * g2 * V * iN6 + 210 * g4 * V2 * iN8,
             2 * V * iN3 - 12 * g2 * V * iN5 + 90 * g4 * V2 * iN7,
             -2 * g2 * V * iN4 + 20 * g4 * V2 * iN6,
             2 * g4 * V2 * iN5]
    s = M + N
    s2 = s * s
    D = s2 + k2
    iD = 1 / D
    iD2 = iD * iD
    iD3 = iD2 * iD
    J = [iD, 2 * s * iD2, 2 * (3 * s2 - k2) * iD3, 24 * s * (s2 - k2) * iD2 * iD2]
    if shape != "exp":
        J.append(24 * (k2 * k2 - 10 * k2 * s2 + 5 * s2 * s2) * iD3 * iD2)
        J.append(240 * s * (s2 - 3 * k2) * (3 * s2 - k2) * iD3 * iD3)
    n = len(Q)
    tot = Q[0] * J[0]
    for i in range(1, n + 1):
        qi = (Q[i] if i < n else 0) + M * Q[i - 1]
        tot += qi * J[i]
    return tot


def Qcoeffs(N, V, g2, g4, shape):
    """the Q_n (functions of the v-node only)"""
    iN = 1 / N
    iN2 = iN * iN
    iN3 = iN2 * iN
    iN4 = iN2 * iN2
    iN5 = iN4 * iN
    if shape == "exp":
        return [6 * g2 * V * iN5, 6 * g2 * V * iN4, 2 * g2 * V * iN3]
    iN6 = iN3 * iN3
    iN7 = iN6 * iN
    iN8 = iN4 * iN4
    iN9 = iN8 * iN
    V2 = V * V
    return [6 * V * iN5 - 30 * g2 * V * iN7 + 210 * g4 * V2 * iN9,
            6 * V * iN4 - 30 * g2 * V * iN6 + 210 * g4 * V2 * iN8,
            2 * V * iN3 - 12 * g2 * V * iN5 + 90 * g4 * V2 * iN7,
            -2 * g2 * V * iN4 + 20 * g4 * V2 * iN6,
            2 * g4 * V2 * iN5]


def kernel_pair(M, N, Q, k2):
    """sum_n (Q_n + M Q_{n-1}) J_n(s, k2), s = M + N; Q from Qcoeffs."""
    s = M + N
    s2 = s * s
    D = s2 + k2
    iD = 1 / D
    iD2 = iD * iD
    iD3 = iD2 * iD
    n = len(Q)
    J = [iD, 2 * s * iD2, 2 * (3 * s2 - k2) * iD3, 24 * s * (s2 - k2) * iD2 * iD2]
    if n > 3:
        J.append(24 * (k2 * k2 - 10 * k2 * s2 + 5 * s2 * s2) * iD3 * iD2)
        J.append(240 * s * (s2 - 3 * k2) * (3 * s2 - k2) * iD3 * iD3)
    tot = Q[0] * J[0]
    for i in range(1, n + 1):
        qi = (Q[i] if i < n else 0) + M * Q[i - 1]
        tot += qi * J[i]
    return tot


def phi_z(z, g, nu, nv, shape="lin", symmetric=True):
    """phi(z) = Phi~(sqrt z) for the geminal kite (apex), z complex (acb), gamma = g (arb).
    symmetric: use the (u, v) -> (1-u, 1-v) symmetry of the integrand (nv must be even)."""
    z = acb(z)
    g = arb(g)
    g2 = g * g
    g4 = g2 * g2
    Mu = []
    for u, Wu in Cm.cheb2_nodes(nu):
        M = (1 + u * (1 - u) * z).sqrt()
        Mu.append((u, Wu / (M * M * M), M))
    vnodes = Cm.gl_nodes01(nv)
    if symmetric:
        assert nv % 2 == 0
        vnodes = [(v, 2 * Wv) for v, Wv in vnodes if v < arb(1) / 2]
        assert len(vnodes) == nv // 2
    Nv = []
    for v, Wv in vnodes:
        V = v * (1 - v)
        N = (g2 + V * z).sqrt()
        Nv.append((v, Wv, N, Qcoeffs(N, V, g2, g4, shape)))
    tot = acb(0)
    for v, cv, N, Q in Nv:
        acc = acb(0)
        for u, cu, M in Mu:
            k2 = (u - v) * (u - v) * z
            acc += cu * kernel_pair(M, N, Q, k2)
        tot += cv * acc
    return tot


def phi_ts(w, g, h="1/16", eps=40, shape="lin"):
    """phi at real w with tanh-sinh in u and v (for the Y grid, large real w)."""
    z = acb(w) * acb(w)
    g = arb(g)
    g2 = g * g
    g4 = g2 * g2
    nodes = Cm.ts_nodes01(arb(h), eps)
    Mu = []
    for u, um, Wu in nodes:
        U = u * um
        M = (1 + U * z).sqrt()
        Mu.append((u, Wu * U.sqrt() / (M * M * M), M))
    Nv = []
    for v, vm, Wv in nodes:
        V = v * vm
        N = (g2 + V * z).sqrt()
        Nv.append((v, Wv, N, Qcoeffs(N, V, g2, g4, shape)))
    tot = acb(0)
    for v, cv, N, Q in Nv:
        acc = acb(0)
        for u, cu, M in Mu:
            k2 = (u - v) * (u - v) * z
            acc += cu * kernel_pair(M, N, Q, k2)
        tot += cv * acc
    return tot


# ---------------------------------------------------------------- route A: the momentum triangle
def phi_routeA(w, g, shape="lin"):
    """Phi~(w) = (2pi)^-3 (2pi/w) int int_{|p-q|<=w<=p+q} p q rho~(p) rho~(q) T(p,q,w),
    T = (2pi)^-3 (4 pi) x [(8 pi g)^2 I_221 (exp) | (8 pi)^2 loop_lin (lin)]."""
    w = mp.mpf(w)
    g = mp.mpf(g)
    rt = lambda k: (k * k + 1) ** mp.mpf(-1.5)
    if shape == "exp":
        pref = (2 * mp.pi) ** -3 * (8 * mp.pi * g) ** 2 * 4 * mp.pi
        loop = lambda p, q: Cm.I221(g, g, 0, w * w, p * p, q * q)
    else:
        pref = (2 * mp.pi) ** -3 * (8 * mp.pi) ** 2 * 4 * mp.pi
        loop = lambda p, q: loop_lin(g, g, 0, w * w, p * p, q * q)

    def inner(p):
        return mp.quad(lambda q: q * rt(q) * loop(p, q), [abs(p - w), p + w])
    val = mp.quad(lambda p: p * rt(p) * inner(p), [0, w, 2 * w, 4 * w + 2, mp.inf])
    return (2 * mp.pi) ** -3 * 2 * mp.pi / w * pref * val
