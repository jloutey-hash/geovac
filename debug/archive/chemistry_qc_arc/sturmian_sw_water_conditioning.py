"""
Does the H2+ gerade lever (Paper 60 sec:molecular/sec:resource) get BETTER for a
higher-symmetry molecule (water, C2v), or is it the same wall in a new costume?

H2+ mechanism: the two centers are EQUIVALENT (inversion), so gerade/ungerade
adaptation isolates a perfectly-conditioned, basis-flat gerade sector (cond~2)
holding the sigma_g ground state -> metric penalty does not grow with basis.

Water: O + 2H, C2v.  For an s-only shared-scale (k=1) Coulomb-Sturmian basis the
only symmetry action is the H1<->H2 swap (O sits on the C2 axis, fixed).  So:
   A1 block (holds the ground state) = {O_n}  U  {(H1_n+H2_n)/sqrt2}
   B2 block                          = {(H1_n-H2_n)/sqrt2}
The A1 block STILL contains the O<->H coupling between symmetry-INEQUIVALENT
centers -- the near-linear-dependence symmetry cannot remove.  Question: does
cond(A1) stay flat like H2+ gerade, or grow?

Two-center s-s SW / L2 integrals depend only on inter-center distance (isotropy),
via the exact momentum-space 1D chi-reduction (from debug/sturmian_sw_momentum.py,
validated there vs 2D position grid + direct 3D p-quad).  Grid-free.  Diagnostic only.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from numpy.linalg import cond, eigvalsh

k = 1.0
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
_M = 600001
_chi = np.linspace(1e-8, np.pi, _M)
_cot = 1.0/np.tan(_chi/2.0)
_1mcos = 1.0 - np.cos(_chi)
_sinj = {}
def _sin(j):
    if j not in _sinj: _sinj[j] = np.sin(j*_chi)
    return _sinj[j]
def _sinc(x):
    out = np.ones_like(x); nz = x != 0.0; out[nz] = np.sin(x[nz])/x[nz]; return out

def block(d, nmax, kind="S"):
    """nmax x nmax two-center s-s block at inter-center distance d (d=0 => intra)."""
    sfac = np.ones_like(_chi) if d == 0.0 else _sinc(k*d*_cot)
    wfac = sfac if kind == "S" else _1mcos*sfac
    B = np.zeros((nmax, nmax))
    for a in range(1, nmax+1):
        for b in range(a, nmax+1):
            v = (2.0/np.pi)*_trapz(_sin(a)*_sin(b)*wfac, _chi)
            B[a-1, b-1] = B[b-1, a-1] = v
    return B

# ---- water geometry (equilibrium): R_OH = 1.809 bohr, angle = 104.5 deg ----
R_OH = 1.809
ang  = np.deg2rad(104.5)
R_HH = 2.0*R_OH*np.sin(ang/2.0)     # 2.861 bohr

def water_raw(nmax, kind="S"):
    """Full 3-center matrix, basis order [O(nmax), H1(nmax), H2(nmax)]."""
    I  = block(0.0,  nmax, kind)
    P  = block(R_OH, nmax, kind)     # O-H
    Q  = block(R_HH, nmax, kind)     # H-H
    Z  = np.zeros((nmax, nmax))
    row0 = np.hstack([I, P, P])
    row1 = np.hstack([P.T, I, Q])
    row2 = np.hstack([P.T, Q.T, I])
    return np.vstack([row0, row1, row2])

def c2v_transform(nmax):
    """Orthogonal T mapping raw [O,H1,H2] -> [A1 (2nmax); B2 (nmax)]."""
    N = 3*nmax
    T = np.zeros((2*nmax + nmax, N))
    r = 0
    for n in range(nmax):                 # A1: O functions
        T[r, n] = 1.0; r += 1
    for n in range(nmax):                 # A1: (H1+H2)/sqrt2
        T[r, nmax+n] = T[r, 2*nmax+n] = 1/np.sqrt(2); r += 1
    for n in range(nmax):                 # B2: (H1-H2)/sqrt2
        T[r, nmax+n] = 1/np.sqrt(2); T[r, 2*nmax+n] = -1/np.sqrt(2); r += 1
    return T

def water_blocks(nmax, kind="S"):
    S = water_raw(nmax, kind); T = c2v_transform(nmax)
    Sr = T @ S @ T.T
    a1 = Sr[:2*nmax, :2*nmax]; b2 = Sr[2*nmax:, 2*nmax:]
    offblk = np.abs(Sr[:2*nmax, 2*nmax:]).max()
    return S, a1, b2, offblk

# H2+ reference (2 equivalent centers), gerade sector, R=2 bohr
def h2p_gerade(nmax, R=2.0, kind="S"):
    I = block(0.0, nmax, kind); Pi = block(R, nmax, kind)
    g = I + Pi; u = I - Pi
    return cond(I+Pi if False else np.block([[I,Pi],[Pi,I]])), cond(g), cond(u)

if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    print("="*80)
    print(f"WATER  (O + 2H, C2v)   R_OH={R_OH} bohr  R_HH={R_HH:.3f} bohr  shared CS scale k=1")
    print("="*80)

    print("\n[VALIDATION] C2v transform is orthogonal & block-diagonalizes S exactly:")
    for nmax in (3, 6):
        S, a1, b2, off = water_blocks(nmax, "S")
        full_ev = np.sort(eigvalsh(S)); blk_ev = np.sort(np.concatenate([eigvalsh(a1), eigvalsh(b2)]))
        print(f"   nmax={nmax}: A1|B2 off-block max={off:.2e}   "
              f"max|spectrum(full)-spectrum(A1(+)B2)|={np.abs(full_ev-blk_ev).max():.2e}")

    print("\n[MAIN] Conditioning vs basis size (per center nmax; total dim 3*nmax):")
    print("  Does the A1 ground-state block stay FLAT (like H2+ gerade ~2) or GROW?")
    print(f"  {'nmax':>4} {'dim':>4} | {'cond(rawS)':>10} {'cond(rawL2)':>11} | "
          f"{'cond A1(gs)':>11} {'cond B2':>9} | {'H2+ gerade':>10} {'H2+ ungr':>9}")
    print("  " + "-"*84)
    nmaxes = list(range(2, 13))
    ca1 = []; craw = []
    for nmax in nmaxes:
        S, a1, b2, off = water_blocks(nmax, "S")
        L2 = water_raw(nmax, "m")
        _, cg, cu = h2p_gerade(nmax, 2.0, "S")
        cA1 = cond(a1); cB2 = cond(b2); cR = cond(S)
        ca1.append(cA1); craw.append(cR)
        print(f"  {nmax:>4} {3*nmax:>4} | {cR:>10.1f} {cond(L2):>11.1f} | "
              f"{cA1:>11.1f} {cB2:>9.2f} | {cg:>10.2f} {cu:>9.1f}")

    N = np.array([3*nm for nm in nmaxes], float)
    def fit(c):
        c = np.array(c, float); p = np.polyfit(np.log(N), np.log(c), 1)
        r2 = 1 - np.sum((np.log(c)-np.polyval(p,np.log(N)))**2)/np.sum((np.log(c)-np.log(c).mean())**2)
        return p[0], r2
    pa, ra = fit(ca1); pr, rr = fit(craw)
    print(f"\n  power-law fits:  cond(A1 ground-state block) ~ N^{pa:.2f} (R2={ra:.3f})   "
          f"cond(raw S) ~ N^{pr:.2f} (R2={rr:.3f})")
    print(f"  H2+ reference: gerade sector cond ~2.1 FLAT (basis-independent) -- the lever that worked.")

    print("\n[VERDICT DATA] spectral driver of A1 growth (lam_min -> 0 = near-linear-dependence):")
    print("     nmax  lam_min(A1)  lam_max(A1)   cond(A1)   | is O<->H coupling the driver?")
    for nmax in (3, 6, 9, 12):
        S, a1, b2, off = water_blocks(nmax, "S")
        ev = eigvalsh(a1)
        # A1 with O-H coupling zeroed (test whether O<->H is the culprit): keep only I (O) + (I+Q) (H+)
        nb = nmax
        a1_noOH = a1.copy(); a1_noOH[:nb, nb:] = 0; a1_noOH[nb:, :nb] = 0
        print(f"     {nmax:>3}   {ev.min():.5f}     {ev.max():.5f}    {cond(a1):>8.1f}   | "
              f"cond(A1 w/ O-H coupling zeroed)={cond(a1_noOH):.2f}")
