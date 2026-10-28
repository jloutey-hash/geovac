"""Guard for Paper 59 sec:density_ring, a production geminal fit priced on Y (record
docs/predictions/2026-10-26_production_geminal_fit.md, PG1-PG3).

Y at density scales (1, 2) and apex t, with exp(-gamma r) geminals on the two edges (gamma = 3/2, 3), replaced by a
six-Gaussian fit sum_i c_i exp(-alpha_i gamma^2 r^2): Y_fit = sum_ij c_i c_j Y_GG(alpha_i g1^2, alpha_j g2^2).
Claims guarded:
  * Psi4's shipped fit of exp(-r) (libmints/integralparameters.cc) makes Y too small by 12.9%, 10.6%, 7.9%, 5.5% at
    t = 1/2, 1, 2, 4, the error falling with t;
  * control: the record GF W1 six-term fit reproduces the exact Y to a few 1e-5 on the same engine, so the engine
    and the exact Y agree.
Y_EXACT is the 60-digit quadrature grid's Y (debug/data/pg_results.json, from pf_grid_d60_h0.03125.json), an
independent route.  The Gaussian engine below is a copy of debug/gf_engine.py (no debug/ import).
"""
import numpy as np
import pytest
from scipy.special import dawsn

PI = np.pi
TS = [0.5, 1.0, 2.0, 4.0]
G1, G2 = 1.5, 3.0
Y_EXACT = np.array([0.0002606666452394625, 0.0010263578351971288, 0.0029871508170413487, 0.006058044020756876])
PSI4 = {"alpha": [0.2209, 1.004, 3.622, 12.16, 45.87, 254.4],
        "c": [0.3144, 0.3037, 0.1681, 0.09811, 0.06024, 0.03726]}
W1_6 = {"alpha": [0.0651095395441957129586558886434, 0.158088415115097225441155054492,
                  0.407098898179791786557287609138, 1.1850565186745461807094223365,
                  4.23591553413811181432977214208, 23.1030314881090752031677966778],
        "c": [0.0212215075497409931599584377511, 0.131906430346820846630666087473,
              0.238573127640927500025176460228, 0.241817473442845131479913566066,
              0.184112601688493755571594437381, 0.121983877649037298665530760295]}


def rho(c, k):
    return c ** 3 * (k * k + c * c) ** -1.5


def G_loop(R, mu):
    x = R / (2 * np.sqrt(mu))
    small = x < 1e-6
    out = np.empty_like(R)
    out[~small] = 8 * mu * PI ** 1.5 * dawsn(x[~small]) / R[~small]
    out[small] = 8 * mu * PI ** 1.5 / (2 * np.sqrt(mu))
    return out


def gl(n, lo, hi):
    x, w = np.polynomial.legendre.leggauss(n)
    return (hi - lo) / 2 * x + (hi + lo) / 2, (hi - lo) / 2 * w


def half_line(n, lo, scale):
    u, wu = gl(n, 0.0, 1.0)
    return lo + scale * (u / (1 - u)) ** 2, wu * scale * 2 * u / (1 - u) ** 3


def geo(lo, hi, n, h0=0.5):
    br, h = [lo], h0
    while br[-1] + h < hi:
        br.append(lo + h)
        h *= 2
    br.append(hi)
    xs, ws = zip(*[gl(n, br[j], br[j + 1]) for j in range(len(br) - 1)])
    return np.concatenate(xs), np.concatenate(ws)


def two_sided(lo, hi, n, h0=0.5):
    m = (lo + hi) / 2
    x1, w1 = geo(lo, m, n, h0)
    x2, w2 = geo(-hi, -m, n, h0)
    return np.concatenate([x1, -x2]), np.concatenate([w1, w2])


def Y_GG(A, B, ts, a=1.0, b=2.0, n=24):
    cut = 4.0 + 12 * np.sqrt(A + B)
    x1, w1 = geo(0.0, cut, n)
    x2, w2 = half_line(2 * n, cut, cut)
    w, ww = np.concatenate([x1, x2]), np.concatenate([w1, w2])
    mu, c = A * B / (A + B), A / (A + B)
    pref_T = (2 * PI) ** -3 * 4 * PI * (PI * PI / (A * B)) ** 1.5
    ph = np.empty_like(w)
    for i, wi in enumerate(w):
        cu = wi + 8.0 + 12 * np.sqrt(mu)
        s1, ws1 = geo(wi, cu, n)
        s2, ws2 = half_line(2 * n, cu, cu)
        sg, wsg = np.concatenate([s1, s2]), np.concatenate([ws1, ws2])
        dl, wdl = two_sided(-wi, wi, n, min(0.5, wi / 2))
        S, D = np.meshgrid(sg, dl, indexing="ij")
        p, q = (S + D) / 2, (S - D) / 2
        R = np.sqrt(np.maximum((1 - c) * p * p + c * q * q - c * (1 - c) * wi * wi, 0.0))
        T = pref_T * np.exp(-wi * wi / (4 * (A + B))) * G_loop(R, mu)
        ph[i] = (2 * PI) ** -3 * 2 * PI / wi * np.sum(np.outer(wsg, wdl) * (S * S - D * D) / 8 * rho(a, p) * rho(b, q) * T)
    return np.array([1 / (2 * PI * PI) * np.sum(ww * w * w * rho(t, w) * ph) for t in ts])


def Y_fit(fit):
    al, c = fit["alpha"], fit["c"]
    return sum(c[i] * c[j] * Y_GG(al[i] * G1 ** 2, al[j] * G2 ** 2, TS) for i in range(6) for j in range(6))


@pytest.mark.slow
def test_psi4_fit_misses_y_and_the_engine_control_holds():
    E_ctrl = Y_fit(W1_6) / Y_EXACT - 1
    assert np.all(np.abs(E_ctrl) < 1e-4), E_ctrl
    E = Y_fit(PSI4) / Y_EXACT - 1
    assert np.allclose(E, [-1.285e-1, -1.063e-1, -7.90e-2, -5.49e-2], atol=1e-3), E
    assert np.all(np.diff(E) > 0)
