"""Paper 63, Proposition prop:tree_reduction: every tree of radial densities reduces to one dimension.

At the coincident-centre fibre, a tree of radial densities rho_v joined by radial kernels Y_e(|r_u - r_w|) equals
4 pi int W_root r^2 dr, with W_u = rho_u prod_children Phi_e W_w and
    (Phi_e f)(r) = int_0^oo k_e(r, s) f(s) s^2 ds,   k_e(r, s) = (2 pi/(r s)) int_{|r-s|}^{r+s} Y_e(q) q dq
(k_e = 4 pi / max(r, s) for Coulomb).  Fingerprint weights rho_v = a^3 K0(a r) / (2 pi^2) (Fourier Delta^-3).

Independent routes (the recursion is never checked only against itself):
- the paper's own n = 3, 4 formulas (chain, star; derived and checked against momentum space in their records),
  built from chi = int_0^x s K1 computed separately;
- the angular integral k_e against direct quadrature in mu = cos(angle), Coulomb and Yukawa;
- Monte Carlo of the full three-dimensional integral (no reduction) at n = 5.

Wrong answers rejected (named at each test):
- kernel 4 pi / min(r, s) (inner/outer swapped) or a missing s^2 measure;
- leaf potential chi(r)/r without the scale (chi(a r) is required);
- "the reduction depends on where the tree is rooted";
- a Yukawa line treated as Coulomb (the massive line changes the value).
"""
import math

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid, quad
from scipy.special import k0, k1

U = np.linspace(np.log(1e-10), np.log(90.0), 200001)
R = np.exp(U)
COUL = ("C",)


def _cum0(f):
    return cumulative_trapezoid(f * R, U, initial=0.0)


def _tail(f):
    c = _cum0(f)
    return c[-1] - c


def _rho(a):
    return a**3 * k0(a * R) / (2 * np.pi**2)


def _phi(W, kern):
    if kern[0] == "C":
        return 4 * np.pi * (_cum0(W * R**2) / R + _tail(W * R))
    g = kern[1]
    return 4 * np.pi / (g * R) * (np.exp(-g * R) * _cum0(W * R * np.sinh(g * R))
                                  + np.sinh(g * R) * _tail(W * R * np.exp(-g * R)))


def _recursion(scales, edges, root=0):
    adj = {v: [] for v in range(len(scales))}
    for x, y, k in edges:
        adj[x].append((y, k))
        adj[y].append((x, k))

    def W(v, parent):
        out = _rho(scales[v])
        for w, k in adj[v]:
            if w != parent:
                out = out * _phi(W(w, v), k)
        return out
    return 4 * np.pi * np.trapezoid(W(root, None) * R**3, U)


def _chi():
    return _cum0(R * k1(R))


T = 1.6


def test_reproduces_the_papers_chain_and_star():
    """n = 3, 4: the recursion equals C = (8/pi^3) t^3 int chi^2 K0(t x) dx and
    S = (16/pi^4) t^3 int chi^3 K0(t x) dx/x (sec:density_ring, sec:four_density). Rejects 1/min, a missing
    s^2 and an unscaled chi (each moves the value by far more than 1e-7)."""
    chi = _chi()
    chain = 8 / np.pi**3 * T**3 * np.trapezoid(chi**2 * k0(T * R) * R, U)
    star = 16 / np.pi**4 * T**3 * np.trapezoid(chi**3 * k0(T * R), U)
    assert _recursion([T, 1, 1], [(0, 1, COUL), (0, 2, COUL)]) == pytest.approx(chain, rel=1e-7)
    assert _recursion([T, 1, 1, 1], [(0, 1, COUL), (0, 2, COUL), (0, 3, COUL)]) == pytest.approx(star, rel=1e-7)


@pytest.mark.parametrize("a", [0.6, 1.3, 2.5])
def test_coulomb_leaf_potential_is_two_over_pi_chi_ar_over_r(a):
    """A Coulomb leaf of scale a contributes (2/pi) chi(a r)/r. Rejects chi(r)/r (scale dropped)."""
    pot = _phi(_rho(a), COUL)
    chi_a = np.interp(a * R, R, _chi())
    sel = (R > 0.05) & (R < 20)
    assert np.allclose(pot[sel], (2 / np.pi) * chi_a[sel] / R[sel], rtol=1e-6)


@pytest.mark.parametrize("kind", ["coulomb", "yukawa"])
def test_angular_kernel_against_direct_quadrature(kind):
    """k_e(r, s) = 2 pi int_{-1}^{1} Y(sqrt(r^2 + s^2 - 2 r s mu)) dmu, checked by quadrature in mu."""
    g = 1.5
    Y = (lambda q: 1 / q) if kind == "coulomb" else (lambda q: math.exp(-g * q) / q)
    for r, s in [(0.3, 1.7), (2.0, 0.4), (1.1, 1.3)]:
        direct = 2 * math.pi * quad(lambda m: Y(math.sqrt(r * r + s * s - 2 * r * s * m)), -1, 1, limit=200)[0]
        if kind == "coulomb":
            closed = 4 * math.pi / max(r, s)
        else:
            closed = 2 * math.pi / (g * r * s) * (math.exp(-g * abs(r - s)) - math.exp(-g * (r + s)))
        assert closed == pytest.approx(direct, rel=1e-9)


SC = [1.0, 0.7, 1.3, 2.0, 0.9]
TREES = {
    "star5": [(0, 1, COUL), (0, 2, COUL), (0, 3, COUL), (0, 4, COUL)],
    "path5": [(0, 1, COUL), (1, 2, COUL), (2, 3, COUL), (3, 4, COUL)],
    "fork5": [(0, 1, COUL), (1, 2, COUL), (1, 3, COUL), (3, 4, COUL)],
    "fork5_yukawa": [(0, 1, COUL), (1, 2, ("Y", 1.5)), (1, 3, COUL), (3, 4, COUL)],
}


@pytest.mark.parametrize("name", list(TREES))
def test_root_independence(name):
    """The fibre is the same from every root (consistency of the peeling). Rejects a root-dependent reduction."""
    vals = [_recursion(SC, TREES[name], root=r0) for r0 in range(5)]
    assert max(vals) - min(vals) < 1e-9 * abs(vals[0])


def _mc(scales, edges, nchunk=20, n=1_000_000, seed=7):
    rng = np.random.default_rng(seed)
    cdfs = []
    for a in scales:
        c = _cum0((2 / np.pi) * a**3 * R**2 * k0(a * R))
        cdfs.append(c / c[-1])
    acc = []
    for _ in range(nchunk):
        X = []
        for c in cdfs:
            r = np.interp(rng.random(n), c, R)
            v = rng.normal(size=(n, 3))
            X.append(v / np.linalg.norm(v, axis=1)[:, None] * r[:, None])
        p = np.ones(n)
        for x, y, k in edges:
            d = np.linalg.norm(X[x] - X[y], axis=1)
            p *= 1 / d if k[0] == "C" else np.exp(-k[1] * d) / d
        acc.append(p.mean())
    acc = np.array(acc)
    return acc.mean(), acc.std(ddof=1) / math.sqrt(len(acc))


@pytest.mark.parametrize("name", ["star5", "fork5_yukawa"])
def test_full_three_dimensional_monte_carlo(name):
    """n = 5, no reduction: 2e7 samples r_v ~ rho_v in 3D, mean of prod_e Y_e. Within 4 standard errors of the
    recursion. Rejects treating the Yukawa line as Coulomb (the fork's value moves by about 7x)."""
    rec = _recursion(SC, TREES[name])
    m, se = _mc(SC, TREES[name])
    assert abs(m - rec) < 4 * se
    assert se < 5e-3 * abs(rec)
