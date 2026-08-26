"""Corner asymptotics for T2's outer integral at (s,t)->(0,0).

Matched asymptotic derivation (see debug/sprint_hp_evaluator_memo.md addendum):
with s=rho*alpha, t=rho*(1-alpha), rho=s+t->0 at FIXED alpha, substituting
k=q/sqrt(rho) makes Delta_s=sqrt(alpha q^2+1), Delta_t=sqrt((1-alpha)q^2+1)
BOTH independent of rho, and j0(k b)->1 (b=rho, correction O(rho)). This gives
the LEADING singular behaviour

    J(s,t) ~ S(s,t) = s*t*W(alpha) / sqrt(s+t),   alpha = s/(s+t)

    W(alpha) = int_0^infty K(q;alpha) dq
    K(q;alpha) = e^{-Ds-Dt} (1/Ds^3+3/Ds^4+3/Ds^5)(1/Dt^3+3/Dt^4+3/Dt^5)
    Ds = sqrt(alpha q^2+1), Dt = sqrt((1-alpha) q^2+1)

S(s,t) is homogeneous of degree 3/2 in (s,t) (an s^{3/2}-type corner term, NOT
literally s ln s -- confirmed numerically, see the memo addendum: J(eps,v*eps)/
eps^1.5 converges to a finite v-dependent limit with an O(eps) [not O(eps ln
eps)] correction). This single non-integer-power term already explains the
slow 2D quadrature (rho^{3/2} is C^1 but not C^2 at the corner).

The corner-square [0,delta]^2 integral of S is CLOSED FORM in the radial
variable (elementary power) times a smooth 1D integral in alpha:
    int_{[0,delta]^2} S ds dt = int_0^delta rho^{5/2} drho . int_0^1 H(alpha) dalpha
                               = (2/7) delta^{7/2} . int_0^1 H(alpha) dalpha
    H(alpha) = alpha(1-alpha) W(alpha)
(ds dt = rho drho dalpha; S = rho^{3/2} H(alpha)).
"""
from __future__ import annotations
import sys
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_hp_evaluator import std_gl_nodes  # noqa: E402


def q_grid_sinh_paneled(theta_max, panel_width, n_panel_deg, c0):
    """Same construction as k_grid_sinh_paneled but for the q-integral, with
    q = (2/sqrt(c0)) sinh(theta) so Delta0=sqrt(c0 q^2+1)=cosh(theta) exactly
    at the REFERENCE scale c0 (branch-point-regularizing substitution)."""
    std_nodes = std_gl_nodes(n_panel_deg, mp.mp.prec)
    scale = 2 / mp.sqrt(c0)
    nodes = []
    theta_lo = mp.mpf(0)
    while theta_lo < theta_max - mp.mpf('1e-30'):
        theta_hi = min(theta_lo + panel_width, theta_max)
        half = (theta_hi - theta_lo) / 2
        mid = (theta_hi + theta_lo) / 2
        for x, w in std_nodes:
            theta = mid + half * x
            q = scale * mp.sinh(theta)
            jac = half * scale * mp.cosh(theta)
            nodes.append((q, w * jac))
        theta_lo = theta_hi
    return nodes


def W_of_alpha(alpha, theta_max=None, panel_width=None, n_panel_deg=6):
    """W(alpha) = int_0^inf K(q;alpha) dq via the q = c-scaled sinh grid,
    referenced to c0 = max(alpha, 1-alpha) (the SLOWER-decaying factor sets
    the required q-reach; using the smaller c as reference keeps both
    exponentials resolved)."""
    a = alpha
    b = 1 - alpha
    c0 = min(a, b) if min(a, b) > mp.mpf('1e-30') else max(a, b)
    if theta_max is None:
        theta_max = mp.mpf(16)
    if panel_width is None:
        panel_width = mp.mpf(2)
    nodes = q_grid_sinh_paneled(theta_max, panel_width, n_panel_deg, c0)
    tot = mp.mpf(0)
    for q, w in nodes:
        Ds = mp.sqrt(a * q * q + 1)
        Dt = mp.sqrt(b * q * q + 1)
        fs = 1 / Ds ** 3 + 3 / Ds ** 4 + 3 / Ds ** 5
        ft = 1 / Dt ** 3 + 3 / Dt ** 4 + 3 / Dt ** 5
        tot += w * mp.e ** (-(Ds + Dt)) * fs * ft
    return tot


def H_of_alpha(alpha, **kw):
    if alpha <= 0 or alpha >= 1:
        return mp.mpf(0)
    return alpha * (1 - alpha) * W_of_alpha(alpha, **kw)


def S_corner(s, t, **kw):
    """Leading corner-singular term S(s,t) = s t W(alpha) / sqrt(s+t)."""
    if s <= 0 or t <= 0:
        return mp.mpf(0)
    rho = s + t
    alpha = s / rho
    return s * t * W_of_alpha(alpha, **kw) / mp.sqrt(rho)


def corner_square_integral_of_S(delta, alpha_deg=6, **kw):
    """int_{[0,delta]^2} S ds dt = (2/7) delta^{7/2} * int_0^1 H(alpha) dalpha,
    the alpha-integral done by plain GL (H is smooth/bounded on (0,1),
    vanishing at both ends)."""
    std_nodes = std_gl_nodes(alpha_deg, mp.mp.prec)
    tot = mp.mpf(0)
    for x, w in std_nodes:
        alpha = (x + 1) / 2
        wA = w / 2
        tot += wA * H_of_alpha(alpha, **kw)
    radial = (mp.mpf(2) / 7) * delta ** mp.mpf('3.5')
    return radial * tot


if __name__ == '__main__':
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    mp.mp.dps = dps
    print("Sanity: S(s,t) vs J(s,t) for small rho, several alpha")
    from routeC_probe9 import J_fixed
    from routeC_probe8 import k_grid_sinh_paneled
    knodes = k_grid_sinh_paneled(mp.mpf(20), mp.mpf(2), 6)
    for alpha_str in ['0.5', '0.2', '0.05']:
        alpha = mp.mpf(alpha_str)
        for rho_str in ['0.01', '0.001', '0.0001']:
            rho = mp.mpf(rho_str)
            s = rho * alpha
            t = rho * (1 - alpha)
            Jv = J_fixed(s, t, knodes)
            Sv = S_corner(s, t, n_panel_deg=6)
            print(f"  alpha={alpha_str} rho={rho_str}  J={mp.nstr(Jv,10)}  S={mp.nstr(Sv,10)}  "
                  f"rel_diff={mp.nstr((Jv-Sv)/Jv,6)}")
