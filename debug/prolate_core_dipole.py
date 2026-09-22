"""Phase 1 (relaxable core, POLARIZATION) increment: the Li 1s^2 core's DIPOLE
Hartree screening of the valence, as an anisotropic one-body potential in the
prolate basis.  Phase 0 proved radial breathing is inert; the v5.15.2 mechanism
says the core must change SHAPE (lean toward the bond), i.e. acquire a dipole.

Polarizable core orbital:  phi_core = 1s(zc) + lambda * 2p_z(zc)   (2p_z along +z,
the bond axis, pointing at H).  The 2-electron density's DIPOLAR (l=1) part is the
1s * 2p_z cross term:

    rho_1(r, theta) = u(r) cos(theta),   u(r) = 4 lambda (zc^4/pi) r e^{-2 zc r}
    (normalized 1s = zc^{3/2}/sqrt(pi) e^{-zc r}; 2p_z = zc^{5/2}/sqrt(pi) z e^{-zc r};
     factor 4 = 2 electrons * 2 cross terms).

Its Hartree potential is a pure l=1 field  V(r,theta) = V_H_dip(r) cos(theta), with
the standard multipole radial solution (l=1):

    V_H_dip(r) = (4 pi / 3) [ (1/r^2) INT_0^r u(r') r'^3 dr'  +  r INT_r^inf u(r') dr' ]

both integrals elementary (incomplete-gamma).  Limits verified analytically:
  r->inf : V_H_dip -> d / r^2   with the dipole  d = 4 lambda / zc   (pure dipole tail);
  r->0   : V_H_dip -> 0  (~ r).

On the prolate grid the field point sits at r_A = (R/2)(xi+eta) from Li, and the
angle from the bond axis (measured at the Li focus) is

    cos(theta_A) = z_A / r_A = (xi eta + 1) / (xi + eta)

(=+1 toward H at eta=+1; =-1 away at eta=-1, xi>1; 0/0 at the Li focus xi=1,eta=-1).

Model core polarizability (parameter-free, from zc alone): two hydrogenic core
electrons of effective charge zc, alpha(Z) = 9/(2 Z^4) each ->
    alpha_c = 9 / zc^4  = 0.1725 a.u. at zc=2.6875   (true Li+ = 0.1925; close).
Dipole per unit lambda:  d = 4 / zc.

VALIDATION (run this file): closed-form V_H_dip(r_A) vs a direct 3D integral of
rho_1(r')/|r-r'| with the field point on the +z axis (where cos theta = 1).

Run from root:  python debug/prolate_core_dipole.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from prolate_core_hartree import ZC_LI                        # noqa: E402  (2.6875)

ALPHA_LIP_TRUE = 0.1925      # true Li+ static dipole polarizability (a.u.), for sensitivity


def alpha_c_model(zc=ZC_LI):
    """Parameter-free model core polarizability: 2 hydrogenic electrons, 9/(2 zc^4) each."""
    return 9.0 / zc**4


def d_per_lambda(zc=ZC_LI):
    """Dipole moment of rho_1 per unit admixture lambda:  d = 4 lambda / zc."""
    return 4.0 / zc


def cos_theta_A(xi, eta):
    """Angle from the bond axis at the Li focus:  (xi eta + 1)/(xi + eta)."""
    return (xi * eta + 1.0) / (xi + eta)


def V_H_dip_closed(r_A, zc=ZC_LI, lam=1.0):
    """Closed-form l=1 Hartree potential V_H_dip(r_A) of the 1s*2p_z cross density."""
    a = 2.0 * zc
    K = 4.0 * lam * zc**4 / np.pi
    ex = np.exp(-a * r_A)
    # I_in = K INT_0^r r'^4 e^{-a r'} dr'   (elementary)
    I_in = K * (24.0 / a**5
                - ex * (r_A**4 / a + 4.0 * r_A**3 / a**2 + 12.0 * r_A**2 / a**3
                        + 24.0 * r_A / a**4 + 24.0 / a**5))
    # I_out = K INT_r^inf r' e^{-a r'} dr'  (elementary)
    I_out = K * ex * (r_A / a + 1.0 / a**2)
    return (4.0 * np.pi / 3.0) * (I_in / r_A**2 + r_A * I_out)


def V_H_dip_direct(r_A, zc=ZC_LI, lam=1.0, n=600):
    """Direct 3D integral of rho_1(r')/|r-r'| with field point at r_A on +z axis
    (cos theta = 1 there, so this returns V_H_dip(r_A) itself)."""
    a = 2.0 * zc
    K = 4.0 * lam * zc**4 / np.pi
    rr = np.linspace(1e-4, 30.0 / a, n)
    cc = np.linspace(-1.0, 1.0, n)          # cos(theta')
    RR, CC = np.meshgrid(rr, cc, indexing='ij')
    u = K * RR * np.exp(-a * RR)            # u(r')  (dipolar radial density)
    dist = np.sqrt(r_A**2 + RR**2 - 2.0 * r_A * RR * CC + 1e-30)
    integrand = u * CC / dist * (RR**2)     # rho_1 = u cos(theta') ; d3r' = r'^2 dr' dc dphi
    _trap = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    return 2.0 * np.pi * _trap(_trap(integrand, cc, axis=1), rr)


def validate_closed_form():
    """Validate V_H_dip_closed by three independent handles that DON'T depend on the
    singular 2D reference: (i) the exact dipole tail d=4/zc; (ii) large-r agreement
    with the direct integral (smooth there); (iii) monotone convergence of the
    singular small-r reference toward the closed form as the grid refines.  (Doing
    the angular integral analytically collapses the 'direct' form back INTO the
    closed form, so the singular 2D quadrature is the only independent numeric check,
    and it is singularity-limited at small r_A -- it converges, it does not disagree.)"""
    print("Validation: closed-form V_H_dip(r_A)   (zc=%.4f, lam=1)" % ZC_LI)
    ok = True

    # (i) exact dipole tail
    tail = V_H_dip_closed(30.0) * 30.0**2
    dt = abs(tail - d_per_lambda())
    print("  (i)  dipole tail  V_H_dip(30)*30^2 = %.6f   d=4/zc = %.6f   |diff|=%.1e"
          % (tail, d_per_lambda(), dt))
    ok = ok and dt < 1e-3

    # (ii) large-r agreement with the (smooth-there) direct integral
    print("  (ii) large-r closed vs direct:")
    for r_A in (2.0, 4.0):
        c = V_H_dip_closed(r_A); d = V_H_dip_direct(r_A, n=800)
        rel = abs(c - d) / abs(c)
        print(f"        r_A={r_A}:  closed={c:.6f}  direct={d:.6f}  rel={rel:.2e}")
        ok = ok and rel < 5e-3

    # (iii) small-r: the singular reference must CONVERGE toward closed as n grows
    print("  (iii) small-r convergence of the singular reference toward closed:")
    for r_A in (0.3, 0.6, 1.0):
        c = V_H_dip_closed(r_A)
        errs = [abs(c - V_H_dip_direct(r_A, n=n)) / abs(c) for n in (800, 3200)]
        conv = errs[1] < errs[0]
        print(f"        r_A={r_A}:  closed={c:.6f}  rel(n=800)={errs[0]:.2e} "
              f"-> rel(n=3200)={errs[1]:.2e}  {'converging' if conv else 'DIVERGING'}")
        ok = ok and conv and errs[1] < 8e-3

    print("  alpha_c(model) = %.5f a.u.   (true Li+ = %.4f)" % (alpha_c_model(), ALPHA_LIP_TRUE))
    print("  DIPOLE HARTREE CLOSED FORM", "VALIDATED" if ok else "MISMATCH")
    return ok


if __name__ == '__main__':
    validate_closed_form()
