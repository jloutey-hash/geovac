"""
Sprint GD-2: t-robustness audit of the Moebius substrate-level mechanism.

Route C (single t=1.0) found soft_IR_frac(alpha) -> 1/(2 alpha), giving
F = 1/(2(1 - soft_IR_frac)) = alpha/(2alpha-1). A continuum-Weyl estimate
suggests soft_IR_frac ~ 2 sqrt(t/pi)/alpha (t-DEPENDENT), which would clash
with a clean 1/(2alpha). Audit: sweep t, check robustness.
"""
import numpy as np
from geovac.gravity.warped_dirac import DiscreteWedgeDiracSpectral

N_rho, a, N_0 = 100, 0.05, 60
alphas = [2.0, 3.0, 5.0]
ts = [0.25, 0.5, 1.0, 2.0, 4.0]

def radial_evals_for_alpha(alpha):
    """Eigenvalues per azimuthal mode (t-independent). Returns list of (m_eff, evals)."""
    N_phi = int(round(alpha * N_0)); actual = N_phi / N_0
    w = DiscreteWedgeDiracSpectral(N_rho=N_rho, a=a, N_phi=N_phi, alpha=actual)
    out = []
    for k_idx in range(N_phi):
        k = k_idx if k_idx <= N_phi//2 else k_idx - N_phi
        m_eff = (k + 0.5)/actual
        ev = np.linalg.eigvalsh(w._hermitian_radial_laplacian(m_eff))
        out.append((m_eff, ev))
    return actual, out

print("="*76)
print("t-robustness of soft_IR_frac and (1-soft_IR_frac)*F  [should be 1/2 if robust]")
print("="*76)
print(f"{'alpha':>6}{'t':>7}{'soft_IR_frac':>14}{'1/(2a)':>9}{'(1-X)*F':>10}{'F_meas/F_moeb':>14}")
rows = []
for alpha in alphas:
    actual, modes = radial_evals_for_alpha(alpha)
    F_moeb = actual/(2*actual-1)
    thr = 1.0/actual + 1e-6
    for t in ts:
        Km = [(m_eff, 2.0*np.sum(np.exp(-ev*t))) for m_eff, ev in modes]
        K_tot = sum(k for _, k in Km)
        K_soft = sum(k for m_eff, k in Km if abs(m_eff) <= thr)
        X = K_soft / K_tot
        # F implied by the substrate relation F = 1/(2(1-X)):
        F_meas = 1.0/(2*(1-X))
        rows.append((actual, t, X, 1/(2*actual), (1-X)*F_moeb, F_meas/F_moeb))
        print(f"{actual:>6.2f}{t:>7.2f}{X:>14.5f}{1/(2*actual):>9.4f}"
              f"{(1-X)*F_moeb:>10.4f}{F_meas/F_moeb:>14.4f}")
    print()

print("="*76)
print("VERDICT")
print("="*76)
# robustness: does (1-X)*F stay ~0.5 across t at each alpha?
for alpha in alphas:
    vals = [r[4] for r in rows if abs(r[0]-int(round(alpha*N_0))/N_0)<1e-9]
    spread = max(vals)-min(vals)
    print(f"  alpha={alpha}: (1-X)*F over t in {ts} = "
          f"[{min(vals):.4f}, {max(vals):.4f}], spread={spread:.4f}, "
          f"{'ROBUST' if spread<0.02 else 'T-DEPENDENT'}")
