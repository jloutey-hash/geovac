"""Paper 51 G4-1/G4-2 replica-method closed forms (symbolic pins).

S_BH = A Lambda^2 / (12 pi) from the conical-tip action via the replica
derivative, and G_N = 3 pi / Lambda^2 from the area law -- the algebra
behind the paper's eq:S_BH and the G4-2 Bekenstein--Hawking block.

Added 2026-07-04 (group5 certifying /qa run): the panel found the
replica algebra had no test anywhere in the suite.  The former
'factor of 2' between the G7 and G4-2 routes was the rank-2 G2 trace and
is withdrawn (2026-11-07); test_action_route_G_equals_entropy_route_G now
pins the equality from both spectra.  These pins are pure sympy; the
Sommerfeld--Cheeger tip coefficient (1/12)(1/alpha - alpha) is the
classical INPUT (cited, not derived here); what is pinned is the
paper's algebra downstream of it.
"""

import sympy as sp


def test_replica_sbh_closed_form():
    """eq:S_BH: S_BH = -dI_E/dalpha|_1 = r_h^2 Lambda^2/3 = A Lambda^2/(12 pi)."""
    alpha, rh, Lam = sp.symbols("alpha r_h Lambda", positive=True)
    I_conical = (rh ** 2 * Lam ** 2 / 6) * (1 / alpha - alpha)
    S_BH = -sp.diff(I_conical, alpha).subs(alpha, 1)
    assert sp.simplify(S_BH - rh ** 2 * Lam ** 2 / 3) == 0
    A = 4 * sp.pi * rh ** 2
    assert sp.simplify(S_BH - A * Lam ** 2 / (12 * sp.pi)) == 0


def test_replica_tip_term_vanishes_smooth():
    """The tip contribution vanishes at alpha = 1 (smooth disk)."""
    alpha, rh, Lam = sp.symbols("alpha r_h Lambda", positive=True)
    I_conical = (rh ** 2 * Lam ** 2 / 6) * (1 / alpha - alpha)
    assert sp.simplify(I_conical.subs(alpha, 1)) == 0


def test_newton_constant_from_area_law():
    """G4-2: A/(4 G_N) = S_BH = A Lambda^2/(12 pi)  =>  G_N = 3 pi/Lambda^2."""
    rh, Lam, G = sp.symbols("r_h Lambda G_N", positive=True)
    A = 4 * sp.pi * rh ** 2
    S_BH = A * Lam ** 2 / (12 * sp.pi)
    sol = sp.solve(sp.Eq(A / (4 * G), S_BH), G)
    assert len(sol) == 1
    assert sp.simplify(sol[0] - 3 * sp.pi / Lam ** 2) == 0


def test_action_route_G_equals_entropy_route_G():
    """G7 = G4-2 under the rank-4 trace (Paper 51 remark on the former
    factor of 2; rebuilt 2026-11-07, group5 guard pass).

    Both routes are computed from spectra, Gaussian cutoff:
      action route  -- Tr exp(-D^2/Lambda^2) per unit time on S^3_R x R,
                       CH spectrum, rank-4 bundle, matched to
                       -(1/16 pi G) int R_scalar  =>  G = -3 pi / (4 B);
      entropy route -- S_BH = -(d/dalpha) of (1/12)(1/alpha - alpha) times
                       the S^2 Dirac heat trace's leading coefficient, read
                       off the spectrum sum_m 4 m exp(-m^2 t / r_h^2);
                       A / (4 G) = S_BH.
    Rejects the withdrawn bookkeeping ratio 2 (the rank-2 G2 trace).
    Inputs, not tested here: the rank-4 bundle of sec:g2, and the cone
    coefficient (1/12)(1/alpha - alpha) as the TOTAL for the 2D Dirac spinor
    (c = 1; Paper 51's remark on the former factor of 2).  A per-component
    reading would restore the factor 2; the equality rests on it.
    """
    import numpy as np
    from scipy import integrate

    Lam = 8.0
    # action route
    s1 = np.array([[0, 1], [1, 0]], dtype=complex)
    s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    D2 = (1.1 * s1 + 0.4 * s2) @ (1.1 * s1 + 0.4 * s2)
    rank = int(round(np.trace(D2).real / (1.1 ** 2 + 0.4 ** 2)))

    def density(R):
        tot = 0.0
        for n in range(4000):
            l2 = ((n + 1.5) / R) ** 2
            if l2 / Lam ** 2 > 60.0:
                break
            tot += 2 * (n + 1) * (n + 2) * np.exp(-l2 / Lam ** 2) * Lam / (2.0 * np.sqrt(np.pi))
        return rank * tot

    R1, R2 = 1.0, 2.0
    _, B = np.linalg.solve(np.array([[R1 ** 3, R1], [R2 ** 3, R2]]),
                           np.array([density(R1), density(R2)]))
    G_action = -3.0 * np.pi / (4.0 * B)

    # entropy route: leading coefficient c of K_{S^2}(t) ~ c r_h^2 / t
    r_h = 2.0
    ts = np.array([1e-3, 2e-3])
    K = [np.sum(4.0 * np.arange(1, 4000) * np.exp(-np.arange(1, 4000) ** 2 * t / r_h ** 2))
         for t in ts]
    # K = c r_h^2/t + d: two-point solve
    c_rh2, _d = np.linalg.solve(np.array([[1 / ts[0], 1.0], [1 / ts[1], 1.0]]), np.array(K))
    # Gaussian cutoff: t^-1 -> Lambda^2 phi(1) = Lambda^2
    I_coeff = (1.0 / 12.0) * c_rh2 * Lam ** 2      # I(alpha) = I_coeff (1/alpha - alpha)
    S_BH = 2.0 * I_coeff                            # -(d/dalpha)(1/alpha - alpha) at 1 = 2
    A = 4.0 * np.pi * r_h ** 2
    G_entropy = A / (4.0 * S_BH)

    assert abs(G_action / G_entropy - 1.0) < 1e-6, (G_action, G_entropy)
    assert abs(G_action * Lam ** 2 - 3.0 * np.pi) < 1e-5
