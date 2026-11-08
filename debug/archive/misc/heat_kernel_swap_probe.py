"""Heat-kernel 'swap' probe: GeoVac binary graph vs continuum S^3.

Diagnostic only (no production code). Sections:
  A  block structure of the production GeometricLattice (binary weights)
  B  exact Bessel-image heat kernel per l-block vs expm
  C  heat trace, spectral dimension, Weyl expansion vs S^3 theta function
  D  kernel multiplicity + same-block pair fraction (heat-kernel support)
  E  limit spectral measure vs the Z^2 square-lattice DOS
  F  top edge: Dirichlet-rectangle infimum = Paper 0's C; accumulation
  G  Rydberg counting under kappa = -1/16
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import sys
import numpy as np
from scipy.linalg import expm
from scipy.special import ive
from scipy.sparse import diags

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from geovac.lattice import GeometricLattice  # noqa: E402

PI = np.pi


def path_lap(N: int) -> np.ndarray:
    L = np.zeros((N, N))
    for i in range(N - 1):
        L[i, i] += 1; L[i + 1, i + 1] += 1
        L[i, i + 1] = L[i + 1, i] = -1
    return L


def p_Z(n: np.ndarray, t: float) -> np.ndarray:
    """Continuous-time SRW on Z, generator = discrete Laplacian (rate 1 each way)."""
    return ive(np.abs(n), 2.0 * t)


def path_kernel_bessel(N: int, t: float, K: int = 60) -> np.ndarray:
    """Neumann (free-end) path P_N heat kernel by the method of images.
    Reflection planes at -1/2 and N-1/2 -> period 2N, mirror y -> -1-y."""
    x = np.arange(N)[:, None]; y = np.arange(N)[None, :]
    out = np.zeros((N, N))
    for k in range(-K, K + 1):
        out += p_Z(x - y + 2 * k * N, t) + p_Z(x + y + 1 + 2 * k * N, t)
    return out


def path_trace_bessel(N: int, t: float, K: int = 200) -> float:
    ks = np.arange(-K, K + 1)
    return N * p_Z(2 * ks * N, t).sum() + 0.5 * (1.0 - np.exp(-4.0 * t))


def path_trace_spec(N: int, t: float) -> float:
    j = np.arange(N)
    return np.exp(-t * (2 - 2 * np.cos(j * PI / N))).sum()


def lattice_L(nmax: int):
    lat = GeometricLattice(max_n=nmax)
    A = lat.adjacency.toarray()
    L = np.diag(A.sum(1)) - A
    return lat, L


print("=" * 72)
print("A. Block structure (binary production lattice)")
for nmax in range(2, 9):
    lat, L = lattice_L(nmax)
    cross = 0
    for i, (n1, l1, m1) in enumerate(lat.states):
        for j, (n2, l2, m2) in enumerate(lat.states):
            if L[i, j] != 0 and l1 != l2:
                cross += 1
    maxdev = 0.0
    for l in range(nmax):
        idx = [lat._state_index[(n, l, m)] for n in range(l + 1, nmax + 1)
               for m in range(-l, l + 1)]
        N, M = nmax - l, 2 * l + 1
        G = np.kron(path_lap(N), np.eye(M)) + np.kron(np.eye(N), path_lap(M))
        maxdev = max(maxdev, np.abs(L[np.ix_(idx, idx)] - G).max())
    ncomp = int(np.sum(np.abs(np.linalg.eigvalsh(L)) < 1e-10))
    print(f"  nmax={nmax}: cross-l entries={cross}, max|L_l - grid(P_{{nmax-l}} x P_{{2l+1}})|={maxdev:.1e},"
          f" dim ker L={ncomp}")

print("=" * 72)
print("B. Bessel-image heat kernel vs expm (full lattice)")
for nmax in (3, 5, 7):
    lat, L = lattice_L(nmax)
    for t in (0.01, 0.3, 1.0, 3.0, 10.0, 40.0):
        E = expm(-t * L)
        Kb = np.zeros_like(E)
        for l in range(nmax):
            idx = [lat._state_index[(n, l, m)] for n in range(l + 1, nmax + 1)
                   for m in range(-l, l + 1)]
            N, M = nmax - l, 2 * l + 1
            Kb[np.ix_(idx, idx)] = np.kron(path_kernel_bessel(N, t), path_kernel_bessel(M, t))
        tr_img = sum(path_trace_bessel(nmax - l, t) * path_trace_bessel(2 * l + 1, t) for l in range(nmax))
        print(f"  nmax={nmax} t={t:6.2f}: max|K_bessel - expm|={np.abs(Kb - E).max():.2e}"
              f"  |Tr_img - Tr|={abs(tr_img - np.trace(E)):.2e}")


def graph_trace(nmax: int, t: float) -> float:
    return sum(path_trace_spec(nmax - l, t) * path_trace_spec(2 * l + 1, t) for l in range(nmax))


def graph_weyl(nmax: int, t: float) -> float:
    s = 1.0 / np.sqrt(4 * PI * t)
    return sum(((nmax - l) * s + 0.5) * ((2 * l + 1) * s + 0.5) for l in range(nmax))


def s3_theta(tau: float, nn: int = 4000) -> float:
    n = np.arange(1, nn + 1, dtype=float)
    return float((n ** 2 * np.exp(-tau * (n ** 2 - 1))).sum())


def dlog(f, t, h=1e-4):
    return -2 * (np.log(f(t * (1 + h))) - np.log(f(t * (1 - h)))) / (np.log(1 + h) - np.log(1 - h))


print("=" * 72)
print("C. Heat trace / spectral dimension d_s = -2 dlnZ/dlnt")
for nmax in (100, 400):
    V = nmax * (nmax + 1) * (2 * nmax + 1) // 6
    print(f"  graph nmax={nmax} (V={V}):")
    for t in (0.5, 2, 10, 50, 200, 1000, 1e4, 1e5, 1e6):
        Z = graph_trace(nmax, t)
        ds = dlog(lambda tt: graph_trace(nmax, tt), t)
        W = graph_weyl(nmax, t)
        print(f"    t={t:9.1f}  Z={Z:14.4f}  d_s={ds:6.3f}  Weyl-rel-err={(W - Z) / Z:+.2e}")
print("  S^3 theta Z(tau)=sum n^2 exp(-tau(n^2-1)):")
for tau in (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 5.0):
    Z = s3_theta(tau)
    ds = dlog(s3_theta, tau)
    print(f"    tau={tau:8.4f}  Z={Z:14.4f}  d_s={ds:6.3f}  (sqrt(pi)/4)tau^-1.5 e^tau={np.sqrt(PI)/4*tau**-1.5*np.exp(tau):.4f}")
print("  Weyl coefficient: graph leading = V/(4 pi t) [2D]; S^3 leading = Vol(S^3)/(4pi tau)^{3/2} = 2pi^2/(4pi)^{1.5}"
      f" = {2*PI**2/(4*PI)**1.5:.6f} = sqrt(pi)/4 = {np.sqrt(PI)/4:.6f}")

print("=" * 72)
print("D. Kernel multiplicity and heat-kernel support")
for nmax in (10, 30, 100, 1000, 10000):
    Vl = np.array([(nmax - l) * (2 * l + 1) for l in range(nmax)], dtype=float)
    V = Vl.sum()
    frac = (Vl ** 2).sum() / V ** 2
    print(f"  nmax={nmax:6d}: dim ker = {nmax} (components); same-block pair fraction = {frac:.6f},"
          f" nmax*frac = {nmax*frac:.4f}  (asymptote 6/5 = 1.2)")

print("=" * 72)
print("E. Normalized spectral CDF vs Z^2 square-lattice DOS (KS distance)")
th = (np.arange(20000) + 0.5) * PI / 20000


def z2_cdf(x: np.ndarray) -> np.ndarray:
    a = 2 - 2 * np.cos(th)  # arcsine law on [0,4]
    y = x[:, None] - a[None, :]
    y = np.clip(y, 0, 4)
    F1 = np.arccos(1 - y / 2) / PI
    return F1.mean(axis=1)


xs = np.linspace(0, 8, 801)
Fz = z2_cdf(xs)
for nmax in (25, 50, 100, 200, 400):
    ev = []
    for l in range(nmax):
        N, M = nmax - l, 2 * l + 1
        a = 2 - 2 * np.cos(np.arange(N) * PI / N)
        b = 2 - 2 * np.cos(np.arange(M) * PI / M)
        ev.append((a[:, None] + b[None, :]).ravel())
    ev = np.sort(np.concatenate(ev))
    Fg = np.searchsorted(ev, xs, side="right") / ev.size
    print(f"  nmax={nmax:4d}: sup|F_graph - F_Z2| = {np.abs(Fg - Fz).max():.4f}")

print("=" * 72)
print("F. Top edge: n^2(8 - lambda) vs Dirichlet-rectangle infimum")
nu = np.linspace(1e-4, 1 - 1e-4, 2000001)
Cmin = (PI ** 2 * (1 / (1 - nu) ** 2 + 1 / (4 * nu ** 2))).min()
Cpaper = PI ** 2 / 4 * (2 + 2 ** (1 / 3)) ** 2 * (1 + 2 ** (-2 / 3))
print(f"  min_nu pi^2[1/(1-nu)^2 + 1/(2nu)^2] = {Cmin:.6f};  Paper-0 C = {Cpaper:.6f}")
for nmax in (50, 100, 200, 400, 800):
    tops = []
    for l in range(nmax):
        N, M = nmax - l, 2 * l + 1
        a = 2 - 2 * np.cos(np.arange(N) * PI / N)
        b = 2 - 2 * np.cos(np.arange(M) * PI / M)
        v = (a[:, None] + b[None, :]).ravel()
        tops.append(v[v > 7.0])
    top = np.sort(np.concatenate(tops))[::-1]
    s = nmax ** 2 * (8 - top)
    c15 = int((s <= 1.5 * Cpaper).sum()); c2 = int((s <= 2 * Cpaper).sum())
    # multiplicities among top 12 (group at 1e-12)
    t12 = top[:12]
    print(f"  nmax={nmax:4d}: n^2(8-lmax)={s[0]:.4f}  #{{<=1.5C}}={c15:4d} (/nmax={c15/nmax:.3f})"
          f"  #{{<=2C}}={c2:5d} (/nmax={c2/nmax:.3f})  top-5 scaled: {np.round(s[:5],3)}")

print("=" * 72)
print("G. Rydberg counting, kappa=-1/16: hydrogen needs exactly 1 eigenvalue with lambda>2 (E<-1/8),")
print("   5 with lambda>8/9 (E<-1/18), 14 with lambda>1/2 (E<-1/32)")
for nmax in (5, 10, 30, 100):
    ev = []
    for l in range(nmax):
        N, M = nmax - l, 2 * l + 1
        a = 2 - 2 * np.cos(np.arange(N) * PI / N)
        b = 2 - 2 * np.cos(np.arange(M) * PI / M)
        ev.append((a[:, None] + b[None, :]).ravel())
    ev = np.concatenate(ev)
    print(f"  nmax={nmax:4d}: #lambda>2 = {(ev > 2).sum():7d}, #>8/9 = {(ev > 8/9).sum():7d},"
          f" #>1/2 = {(ev > 0.5).sum():7d}  (total {ev.size})")
print("=" * 72)
print("H. Per-site trace limit Z(t)/V -> [e^{-2t} I_0(2t)]^2 (Z^2 return probability);")
print("   d_s^{Z2}(t) = 8t(1 - I_1(2t)/I_0(2t)) -> 2; DOS = K(1-(lam-4)^2/16)/(2 pi^2)")
from scipy.special import ellipk  # noqa: E402
for t in (0.5, 1.5, 5.0, 20.0):
    ref = ive(0, 2 * t) ** 2
    dsz = 8 * t * (1 - ive(1, 2 * t) / ive(0, 2 * t))
    row = []
    for nmax in (100, 400, 1600):
        V = nmax * (nmax + 1) * (2 * nmax + 1) / 6
        row.append(f"n{nmax}:{(graph_trace(nmax, t) / V - ref) / ref:+.2e}")
    print(f"  t={t:5.1f}: [e^-2t I0]^2={ref:.6e}, d_s^Z2={dsz:.4f}; rel dev  " + "  ".join(row))
lat15_ds = dlog(lambda tt: graph_trace(15, tt), 1.5)
print(f"  archived Paper 4 check: graph d_s(t=1.5) at nmax=15 (full spectrum) = {lat15_ds:.4f}; Z^2 limit ="
      f" {8*1.5*(1-ive(1,3.0)/ive(0,3.0)):.4f}")
nmax = 400
ev = []
for l in range(nmax):
    N, M = nmax - l, 2 * l + 1
    a = 2 - 2 * np.cos(np.arange(N) * PI / N)
    b = 2 - 2 * np.cos(np.arange(M) * PI / M)
    ev.append((a[:, None] + b[None, :]).ravel())
ev = np.concatenate(ev)
for lam in (1.0, 2.5, 3.5, 5.5, 7.0):
    w = 0.05
    emp = ((ev > lam - w) & (ev < lam + w)).sum() / ev.size / (2 * w)
    th_ = ellipk(1 - (lam - 4) ** 2 / 16) / (2 * PI ** 2)
    print(f"  DOS at lam={lam}: graph(nmax=400, bin 0.1)={emp:.4f}  Z^2 formula={th_:.4f}")
print("DONE")
