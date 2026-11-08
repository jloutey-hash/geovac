"""Fix: transform(a) -> transform(0) at the UV end, and a nonzero constant makes the
truncated integral dominated by the missing tail.  Subtract it (its exact cosine
coefficient is 0 for every integer j>=1) and compare gerade DIRECTLY to the raw
symbol, which avoids the modulating |sin| altogether."""
import numpy as np
GL_N = 48
xg, wg = np.polynomial.legendre.leggauss(GL_N)

def coeff(j, kR, tr, Smax=2000.0):
    e, s = [0.0], 0.0
    while s < Smax:
        s += min(np.pi/2, np.pi/max(2*j*kR/(kR*kR+s*s), 1e-300)); e.append(min(s, Smax))
    e = np.asarray(e); lo, hi = e[:-1], e[1:]
    mid, half = 0.5*(lo+hi), 0.5*(hi-lo)
    s = (mid[:, None] + half[:, None]*xg[None, :]).ravel()
    w = (half[:, None]*wg[None, :]).ravel()
    a = np.sinc(s/np.pi)
    f = tr(a) - tr(np.zeros_like(a))          # kill the UV constant (exact c_j=0, j>=1)
    return float(np.dot(w, np.cos(2*j*np.arctan2(kR, s)) * f * 2*kR/(kR*kR+s*s))/np.pi)

kR = 2.0
print("  GERADE (1+a)^-1/2 vs raw symbol a   [predicted ratio -1/2 exactly]")
print(f"  {'j':>7} {'c_j raw':>14} {'c_j gerade':>14}   ratio     Smax-stable")
for j in (256, 1024, 4096, 16384, 65536):
    r  = coeff(j, kR, lambda a: a)
    g  = coeff(j, kR, lambda a: 1/np.sqrt(1+a))
    g2 = coeff(j, kR, lambda a: 1/np.sqrt(1+a), Smax=4000.0)
    print(f"  {j:7d} {r:+14.6e} {g:+14.6e}   {g/r:+7.4f}   {abs(g-g2):.1e}")

print("\n  UNGERADE (1-a)^-1/2 : 1-a ~ (kR)^2 (pi-chi)^2/24 as chi->pi, so the symbol")
print("  ~ 1/(pi-chi) there.  int |symbol| dchi diverges logarithmically:")
for eps in (1e-2, 1e-4, 1e-6, 1e-8):
    chi = np.linspace(eps, np.pi-1e-12, 4_000_00)
    a = np.sinc(kR/np.tan(chi/2)/np.pi)
    m = a < 1-1e-15
    print(f"    int_0^pi |(1-a)^-1/2| dchi  with (pi-chi) cut at {eps:.0e}"
          f"  =  {np.trapezoid(1/np.sqrt(1-a[m]), chi[m]):.3f}")

# ---------------------------------------------------------------------------
# Result (2026-09-11).  Envelope law for the UV pole of the Paper 60 cross-block
# symbol a(chi) = j0(kR cot(chi/2)), companion to eq:sigma_law at the IR pole:
#
#     |c_j| = (2 pi)^(-1/2) 2^(-3/4) (kR)^(-1/4) j^(-5/4) |sin(2 sqrt(2 kR j) + pi/4)|
#
# Parameter-free (stationary phase on the chi->0 chirp); no fitted constants.
# Verified by three independent quadrature routes over j = 64..65536,
# kR in {1, 2, 5}; 1-3% except at zeros of the modulating sine.
# ---------------------------------------------------------------------------
