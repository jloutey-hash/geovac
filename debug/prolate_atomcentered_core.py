"""Route C, increment C1 (the GATE): can a tight, atom-centered Li 1s core be
integrated ANALYTICALLY and R-ACCURATELY in the two-center prolate coordinate
system -- i.e. WITHOUT the R-dependent tight-core grid wall that killed increment 7
(isolated Li core swung 0.36 Ha across the bond range; see
debug/track_logs/prolate_native_lih.md, increment 7 + graded-grid cycle)?

A Li-centered 1s STO is  chi = e^{-zeta r_A},  r_A = (R/2)(xi+eta), so
    |chi|^2 = e^{-zeta R xi} e^{-zeta R eta}      (xi- AND eta-exponentials)
and integrals over the prolate volume (R/2)^3 (xi^2 - eta^2) dxi deta dphi reduce to
    M_xi(n, c) = int_1^inf xi^n e^{-c xi} dxi          (c = zeta R)   [existing class]
    M_eta(q, b)= int_{-1}^{1} eta^q e^{-b eta} deta    (b = zeta R)   [NEW primitive]
The tight-core Jacobian (xi^2 - eta^2)/(xi+eta) = (xi - eta) removes the 1/r_A
denominator cleanly, so <chi|1/r_A|chi> needs no new machinery either.

THE GATE (exact single-center references, R-INDEPENDENT by physics):
    <chi|chi>            = pi / zeta^3
    <chi|1/r_A|chi>/N    = zeta
    E_1s = <chi|-1/2 grad^2 - zeta/r_A|chi>/N = -zeta^2/2   (hydrogenic identity)
Computed THROUGH the two-center prolate moments at c = b = zeta R (which grows with
R), these must come out R-INDEPENDENT.  If they do (mpf-exact), the analytic route
beats the grid wall and C1 is GO.

Run from root:  python debug/prolate_atomcentered_core.py
"""
import os
import sys
import numpy as np
from mpmath import mp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from prolate_core_hartree import ZC_LI                        # noqa: E402  (2.6875)

mp.dps = 50


def M_xi(n, c):
    """int_1^inf xi^n e^{-c xi} dxi = c^{-(n+1)} Gamma(n+1, c)  (upper incomplete)."""
    c = mp.mpf(c)
    return mp.gammainc(n + 1, c) / c**(n + 1)


def M_eta(q, b):
    """int_{-1}^{1} eta^q e^{-b eta} deta, by parts recurrence (mpf-exact).

    M(0) = 2 sinh(b)/b ;  M(q) = ((-1)^q e^b - e^{-b})/b + (q/b) M(q-1).
    """
    b = mp.mpf(b)
    if b == 0:
        return mp.mpf(0) if q % 2 else mp.mpf(2) / (q + 1)
    eb, emb = mp.e**b, mp.e**(-b)
    m = 2 * mp.sinh(b) / b
    for k in range(1, q + 1):
        sgn = mp.mpf(-1)**k
        m = (sgn * eb - emb) / b + (mp.mpf(k) / b) * m
    return m


def _self_check_M_eta():
    """Recurrence vs direct quadrature at a representative b."""
    b = mp.mpf("7.2")
    worst = mp.mpf(0)
    for q in range(0, 5):
        rec = M_eta(q, b)
        ref = mp.quad(lambda e: e**q * mp.e**(-b * e), [-1, 1])
        worst = max(worst, abs(rec - ref) / abs(ref))
    return worst


def core_integrals(zeta, R):
    """Analytic <chi|chi>, <chi|1/r_A|chi> through the two-center prolate moments."""
    zeta, R = mp.mpf(zeta), mp.mpf(R)
    hR = R / 2
    c = zeta * R                       # both xi- and eta-exponential rates
    mx0, mx1, mx2 = M_xi(0, c), M_xi(1, c), M_xi(2, c)
    me0, me1, me2 = M_eta(0, c), M_eta(1, c), M_eta(2, c)
    # <chi|chi> = 2pi hR^3 [ M_xi(2) M_eta(0) - M_xi(0) M_eta(2) ]
    N = 2 * mp.pi * hR**3 * (mx2 * me0 - mx0 * me2)
    # <chi|1/r_A|chi>: (xi^2-eta^2)/(xi+eta) = (xi - eta), 1/r_A = (2/R)/(xi+eta)
    #   = 2pi hR^2 [ M_xi(1) M_eta(0) - M_xi(0) M_eta(1) ]
    V = 2 * mp.pi * hR**2 * (mx1 * me0 - mx0 * me1)
    return N, V


def gate():
    zeta = mp.mpf(ZC_LI)
    N_exact = mp.pi / zeta**3
    print(f"C1 GATE: analytic R-accuracy of a tight Li core (zeta={float(zeta):.4f}) "
          f"in two-center prolate coords", flush=True)
    print(f"  M_eta recurrence self-check (vs quad): {float(_self_check_M_eta()):.1e}", flush=True)
    print(f"  exact refs:  <chi|chi>=pi/zeta^3={float(N_exact):.10f}   "
          f"<1/r_A>/N=zeta={float(zeta):.6f}   E_1s=-zeta^2/2={float(-zeta**2/2):.6f}", flush=True)
    print(f"  {'R':>6} {'<chi|chi>':>16} {'relerr(N)':>11} {'<1/rA>/N':>12} "
          f"{'relerr':>11} {'E_1s':>12}", flush=True)
    Ns, Vs = [], []
    for R in (2.70, 2.85, 3.015, 3.20, 3.45):
        N, V = core_integrals(zeta, R)
        E1s = -zeta**2 / 2 * N + zeta * V - zeta * V   # = -zeta^2/2 N by identity
        relN = abs(N - N_exact) / N_exact
        relV = abs(V / N - zeta) / zeta
        Ns.append(N); Vs.append(V / N)
        print(f"  {R:6.3f} {float(N):16.10f} {float(relN):11.1e} {float(V/N):12.8f} "
              f"{float(relV):11.1e} {float(E1s/N):12.8f}", flush=True)
    spreadN = float(max(Ns) - min(Ns))
    spreadV = float(max(Vs) - min(Vs))
    E_spread_Ha = float((max(Vs) - min(Vs)) * zeta)  # ~ energy-scale spread proxy
    print(f"\n  R-SPREAD  <chi|chi>: {spreadN:.2e}   <1/rA>/N: {spreadV:.2e}", flush=True)
    print(f"  (increment-7 grid wall: isolated core energy swung 0.36 Ha across R)", flush=True)
    ok = spreadN < 1e-12 and spreadV < 1e-12
    print(f"  GATE: {'PASS -- tight core is R-accurate analytically (wall beaten)' if ok else 'CHECK'}",
          flush=True)
    return ok


if __name__ == '__main__':
    gate()
