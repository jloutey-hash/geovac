"""Paper 32, scope paragraph after thm:h1_yukawa_non_selection (2026-10-29): regression guard.

The paragraph says: on the true Camporesi-Higuchi window the compressed gauge fluctuation
a[D,b] + h.c. is not gauge invariant; its pure-gauge companion (a = 1) is an edge-crossing
boundary sum; and that sum dominates the gauge functional whenever the heat cutoff sits at
the window edge (fixed tau = t (n_max + 1/2)^2), so the CCM coefficient 4 cannot be read at
finite cutoff by this instrument.  Record: docs/predictions/2026-10-29_ccm_yukawa_coefficient_window.md
(YC, GO); mechanism derived in the gauge-truncation literature memo, CHANGELOG v5.21.27.

What this file checks, all recomputed here (no debug/ import; the construction is ported):
  1. in the continuum the pure-gauge companion vanishes: [chi, [chi, D]] = 0 on every
     row whose neighbours all exist (order zero plus order one for a function);
  2. on the window the companion's eps^2 heat coefficient EQUALS the edge-crossing sum
     -sum_{i in W, j not in W} |chi_ij|^2 f'(lambda_i) (lambda_i - lambda_j), f = exp(-t x^2);
  3. at fixed tau the companion ratio rho_pg = |k_pg / k_A| reproduces the record's scored
     table and is of order one or larger at every point; at fixed t = 0.1 it is small and
     falls with n_max.

Construction (as in the record).  L^2(SU(2)) (x) C^2 with basis sqrt(2j+1) D^j_{m m'} (x) e_s;
the Dirac eigenbasis is the exact j (x) 1/2 Clebsch-Gordan coupling, J = j + 1/2 ->
lambda = +(2j + 3/2), J = j - 1/2 -> lambda = -(2j + 1/2).  a = D^{1/2}_{-1/2,+1/2},
b = D^{1/2}_{+1/2,+1/2}, acting by Clebsch-Gordan on (m, m'), identity on spin.
Written as a separate step after the paragraph (CLAUDE.md Sec. 9, guard-writing rule).
"""
import numpy as np
import scipy.sparse as sp


# ----------------------------------------------------------------------------- construction
def cg_half(j1x2, m1x2, sx2, Jx2):
    """<j1 m1; 1/2 s | J, m1+s> (Condon-Shortley), doubled-integer inputs."""
    j1x2, m1x2, sx2, Jx2 = (np.asarray(x) for x in (j1x2, m1x2, sx2, Jx2))
    j1 = j1x2 / 2.0
    M = (m1x2 + sx2) / 2.0
    den = 2 * j1 + 1
    plus = sx2 > 0
    v_up = np.where(plus, np.sqrt(np.maximum(j1 + M + 0.5, 0) / den),
                    np.sqrt(np.maximum(j1 - M + 0.5, 0) / den))
    v_dn = np.where(plus, -np.sqrt(np.maximum(j1 - M + 0.5, 0) / den),
                    np.sqrt(np.maximum(j1 + M + 0.5, 0) / den))
    out = np.where(Jx2 == j1x2 + 1, v_up, np.where(Jx2 == j1x2 - 1, v_dn, 0.0))
    return np.where(np.abs(M) <= Jx2 / 2.0 + 1e-12, out, 0.0)


class Sphere:
    """Ambient basis up to 2j = J2max and its exact coupled Dirac eigenbasis."""

    def __init__(self, J2max):
        self.J2max = J2max
        self.offs, cnt = {}, 0
        for j2 in range(J2max + 1):
            self.offs[j2] = cnt
            cnt += 2 * (j2 + 1) ** 2
        self.n = cnt
        lam, lev, rows, cols, vals, col = [], [], [], [], [], 0
        for j2 in range(J2max + 1):
            for Jc2 in (j2 + 1, j2 - 1):
                if Jc2 < 0:
                    continue
                for m2 in range(-j2, j2 + 1, 2):
                    for MJ2 in range(-Jc2, Jc2 + 1, 2):
                        for sx2 in (1, -1):
                            mp2 = MJ2 - sx2
                            if abs(mp2) > j2:
                                continue
                            c = float(cg_half(j2, mp2, sx2, Jc2))
                            if c != 0.0:
                                rows.append(self.idx(j2, m2, mp2, sx2))
                                cols.append(col)
                                vals.append(c)
                        lam.append(j2 + 1.5 if Jc2 == j2 + 1 else -(j2 + 0.5))
                        lev.append(j2 if Jc2 == j2 + 1 else j2 - 1)
                        col += 1
        assert col == self.n
        self.V = sp.csr_matrix((vals, (rows, cols)), shape=(self.n, self.n))
        self.lam = np.array(lam)
        self.level = np.array(lev)

    def idx(self, j2, m2, mp2, sx2):
        return self.offs[j2] + (((j2 - m2) // 2) * (j2 + 1) + (j2 - mp2) // 2) * 2 + (0 if sx2 > 0 else 1)

    def multiplier(self, M2, Mp2):
        """Multiplication by D^{1/2}_{M, M'} in the Dirac eigenbasis."""
        rows, cols, vals = [], [], []
        for j1 in range(self.J2max + 1):
            ms = np.arange(-j1, j1 + 1, 2)
            m1, mp1 = (x.ravel() for x in np.meshgrid(ms, ms, indexing='ij'))
            for j2 in (j1 + 1, j1 - 1):
                if j2 < 0 or j2 > self.J2max:
                    continue
                m2, mp2 = m1 + M2, mp1 + Mp2
                c = (np.sqrt((j1 + 1) / (j2 + 1)) * cg_half(np.full_like(m1, j1), m1, np.full_like(m1, M2), np.full_like(m1, j2))
                     * cg_half(np.full_like(m1, j1), mp1, np.full_like(m1, Mp2), np.full_like(m1, j2)))
                ok = (np.abs(m2) <= j2) & (np.abs(mp2) <= j2) & (c != 0)
                for s in (1, -1):
                    for k in np.nonzero(ok)[0]:
                        rows.append(self.idx(j2, int(m2[k]), int(mp2[k]), s))
                        cols.append(self.idx(j1, int(m1[k]), int(mp1[k]), s))
                        vals.append(c[k])
        X = sp.csr_matrix((vals, (rows, cols)), shape=(self.n, self.n))
        return (self.V.T @ X @ self.V).tocsr()


def comm_D(lam, X):
    L = sp.diags(lam)
    return (L @ X - X @ L).tocsr()


def _phi(x, y):
    """int_0^1 exp(-s x - (1-s) y) ds, stable on the diagonal."""
    d = x - y
    small = np.abs(d) < 1e-12
    safe = np.where(small, 1.0, d)
    return np.where(small, np.exp(-0.5 * (x + y)), np.exp(-y) * (-np.expm1(-safe)) / safe)


def eps2(lam, A, t):
    """eps^2 coefficient of Tr exp(-t (D + eps A)^2), D = diag(lam), exact (Duhamel)."""
    Ac = A.tocoo()
    w = np.abs(Ac.data) ** 2
    mu = lam ** 2
    first = -t * float(np.sum(w * np.exp(-t * mu[Ac.row])))
    i, j = Ac.row, Ac.col
    second = 0.5 * t * t * float(np.sum((lam[i] + lam[j]) ** 2 * w * _phi(t * mu[i], t * mu[j])))
    return first + second


_SPHERES = {}


def sphere(J2max):
    if J2max not in _SPHERES:
        S = Sphere(J2max)
        _SPHERES[J2max] = (S, S.multiplier(-1, +1), S.multiplier(+1, +1))
    return _SPHERES[J2max]


def window_k(n_max, t):
    """(k_A, k_pg) on the true window |lambda| <= n_max + 1/2, multipliers compressed first."""
    S, a, b = sphere(2 * n_max + 2)
    w = np.nonzero(np.abs(S.lam) <= n_max + 0.5)[0]
    aw, bw, lw = a[w][:, w].tocsr(), b[w][:, w].tocsr(), S.lam[w]
    Db = comm_D(lw, bw)
    om = (aw @ Db).tocsr()
    A = (om + om.getH()).tocsr()
    Apg = (Db + Db.getH()).tocsr()
    return eps2(lw, A, t), eps2(lw, Apg, t)


# ----------------------------------------------------------------------------- tests
def test_pure_gauge_companion_vanishes_in_the_continuum():
    """With chi = i(b - b^dag) self-adjoint, i[chi, D] = [D, b] + h.c. is the pure-gauge
    perturbation, and [chi, [chi, D]] = 0 on every interior row: a function commutes with
    Clifford multiplication by its own differential.  This is why the continuum companion
    is zero and the window one is not.

    Wrong answers rejected: D replaced by |D| (the sign(D) confusion; the double commutator
    is then of order one), and a chi whose i[chi, D] is not the companion.
    """
    S, a, b = sphere(12)
    lam = S.lam
    chi = (1j * (b - b.getH())).tocsr()
    B1 = (1j * (chi @ sp.diags(lam) - sp.diags(lam) @ chi)).tocsr()
    Db = comm_D(lam, b)
    assert abs(B1 - (Db + Db.getH())).max() < 1e-13
    dc = (chi @ comm_D(lam, chi) - comm_D(lam, chi) @ chi).tocsr()
    inner = np.nonzero(S.level <= S.J2max - 4)[0]           # every neighbour exists
    sub = dc[inner][:, inner]
    scale = abs(comm_D(lam, chi)).max()
    assert sub.nnz == 0 or abs(sub).max() < 1e-12 * scale, abs(sub).max()
    assert scale > 1.0


def test_window_companion_equals_the_edge_crossing_sum():
    """On the window W the pure-gauge eps^2 heat coefficient equals, to rounding,
    -sum_{i in W, j not in W} |chi_ij|^2 f'(lambda_i)(lambda_i - lambda_j), f(x) = exp(-t x^2),
    at fixed tau in {1, 4} and n_max in {2, ..., 5}; the sum is nonzero.

    Wrong answers rejected: an interior (non-crossing) sum, a companion that is zero on
    the window, and an edge weight evaluated at lambda_j instead of lambda_i.
    """
    for tau in (1.0, 4.0):
        for n in (2, 3, 4, 5):
            t = tau / (n + 0.5) ** 2
            _, kpg = window_k(n, t)
            S, a, b = sphere(2 * n + 2)
            chi = (1j * (b - b.getH())).tocoo()
            inW = np.abs(S.lam) <= n + 0.5
            cross = inW[chi.row] & ~inW[chi.col]
            i, j = chi.row[cross], chi.col[cross]
            li, lj = S.lam[i], S.lam[j]
            fprime = -2 * t * li * np.exp(-t * li ** 2)
            edge = -float(np.sum(np.abs(chi.data[cross]) ** 2 * fprime * (li - lj)))
            assert abs(edge) > 1e-3, (tau, n, edge)
            assert abs(kpg - edge) < 1e-10 * max(1.0, abs(edge)), (tau, n, kpg, edge)


# record YC, scored values (tau, n_max) -> rho_pg, as the scorer stored them (yc_score.json,
# 7 significant figures).  The record's disposition table prints 14.94 for (4, 3); the
# stored value is 14.934717, so 14.93 (erratum appended to the record 2026-10-31).
YC_RHO = {(1.0, 2): 4.344522, (1.0, 3): 2.814859, (1.0, 4): 2.412316, (1.0, 5): 2.231772,
          (4.0, 2): 1.912111, (4.0, 3): 14.93472, (4.0, 4): 6.018829, (4.0, 5): 3.372803}


def test_companion_dominates_at_fixed_tau_and_dies_at_fixed_t():
    """At fixed tau the companion ratio rho_pg = |k_pg / k_A| reproduces the record's scored
    values (1.9 to 15) and exceeds 1 at every point, so the gauge functional is dominated by
    the gauge-breaking term at the window edge.  At fixed t = 0.1 it is small and falls
    strictly with n_max (record: 1.2e-2 at n_max = 10 to 1.4e-6 at 14).

    Wrong answers rejected: a companion that is negligible at the window edge (the
    paragraph's 'dominates' would be false), the record table drifting from the
    instrument, and a defect that does not fall when the cutoff moves inside the window.
    """
    for (tau, n), rho_rec in YC_RHO.items():
        kA, kpg = window_k(n, tau / (n + 0.5) ** 2)
        rho = abs(kpg / kA)
        assert rho > 1.0, (tau, n, rho)
        assert abs(rho - rho_rec) < 1e-6 * rho_rec, (tau, n, rho, rho_rec)
    small = []
    for n in (10, 12, 14):
        kA, kpg = window_k(n, 0.1)
        small.append(abs(kpg / kA))
    assert small[0] < 0.02 and small[-1] < 2e-6, small
    assert small[0] > small[1] > small[2], small
