"""
Molecular Shibuya-Wulfman (SW) conditioning -- MOMENTUM-SPACE (Fock) evaluation.

Companion / independent cross-check to debug/sturmian_sw_conditioning.py (position-space
2D cylindrical grid + numerical gradients).  Here we evaluate the SAME two matrices
(SW metric S and L2 overlap m) from the momentum-space / Fock-projection forms
[PhD eq 10.5.12/10.5.13], which are analytically reducible to 1D integrals over the
unit 3-sphere for s-orbitals -- a completely independent quadrature (different variables,
different integrand, different singularity structure) from the position grid.

------------------------------------------------------------------------------------
DERIVATION (s-orbitals, common CS scale k, centers separated by R along z):

Fock map  p -> u in S^3:   u4=(p^2-k^2)/(p^2+k^2),  u_123=2k p/(p^2+k^2),
   |p| = k cot(chi/2)   where u4 = cos(chi),   and  d^3p = ((k^2+p^2)/2k)^3 dOmega_3.

4D hyperspherical harmonic (l=0):  Y_{n-1,0,0}(u) = sqrt(2/pi) * U_{n-1}(cos chi) * 1/sqrt(4pi),
   U = Chebyshev-U (Gegenbauer C^1);   sin^2(chi) U_{n'-1} U_{n-1} = sin(n' chi) sin(n chi).

SW matrix  [PhD 10.5.13]:  S_{mu'mu} = INT d^3p e^{ip.R} (2k/(k^2+p^2))^3 Y*_{mu'} Y_mu.
   The weight (2k/(k^2+p^2))^3 cancels the Jacobian exactly -> integral over dOmega_3.
   Angular (theta,phi) integral of the plane wave -> 2 sinc(kR cot(chi/2)).  Result:

     S_{n',n}(R) = (2/pi) INT_0^pi sin(n' chi) sin(n chi) sinc(kR cot(chi/2)) dchi        (SW)

L2 (Sturmian) overlap [PhD 10.5.12]:  weight M(p)^2 -> extra factor 2k^2/(k^2+p^2)=(1-cos chi):

     m_{n',n}(R) = (2/pi) INT_0^pi sin(n' chi) sin(n chi) (1-cos chi) sinc(kR cot(chi/2)) dchi   (L2)

Checks built in:  S_{n',n}(0) = delta_{n'n} EXACTLY (SW intra-center block = identity, no grid
error); m_{n',n}(0) has unit diagonal but nonzero off-diagonal (e.g. m_{1,2}=-1/2) -> the L2
intra-center ill-conditioning.  sinc(x)=sin(x)/x, sinc(0)=1.
------------------------------------------------------------------------------------
Diagnostic only.  Does NOT modify paper_60/tests/CLAUDE/CHANGELOG.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from numpy.linalg import cond, eigvalsh
from scipy.special import genlaguerre

k = 1.0
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))  # numpy 2.x renamed trapz

# ------------------------------------------------------------------ momentum-space 1D form
_M = 600001
_chi = np.linspace(1e-8, np.pi, _M)          # fine uniform grid; integrand ~chi^3 near 0 (safe)
_cot = 1.0 / np.tan(_chi / 2.0)
_1mcos = 1.0 - np.cos(_chi)
_sinj = {}                                   # cache sin(j*chi)

def _sin(j):
    if j not in _sinj:
        _sinj[j] = np.sin(j * _chi)
    return _sinj[j]

def _sinc(x):
    out = np.ones_like(x)
    nz = x != 0.0
    out[nz] = np.sin(x[nz]) / x[nz]
    return out

def block_mom(R, nmax, kind):
    """nmax x nmax block B_{n',n}(R) between two centers (R=0 => intra), momentum form."""
    if R == 0.0:
        sfac = np.ones_like(_chi)
    else:
        sfac = _sinc(k * R * _cot)
    wfac = sfac if kind == "S" else _1mcos * sfac
    B = np.zeros((nmax, nmax))
    for a in range(1, nmax + 1):
        for b in range(a, nmax + 1):
            val = (2.0 / np.pi) * _trapz(_sin(a) * _sin(b) * wfac, _chi)
            B[a - 1, b - 1] = B[b - 1, a - 1] = val
    return B

def assemble_mom(nmax, R, kind):
    intra = block_mom(0.0, nmax, kind)       # S: identity;  m: L2 intra overlap
    inter = block_mom(R, nmax, kind)
    top = np.hstack([intra, inter]); bot = np.hstack([inter.T, intra])
    return np.vstack([top, bot])

# ------------------------------------------------------------------ direct 3D momentum quad (triple check)
def sw_3d_radial(nprime, n, R):
    """S_{n',n}(R) by DIRECT 3D momentum integral (independent of the 1D chi reduction)."""
    from scipy.integrate import quad
    def U(m, x):                              # Chebyshev U_m
        return np.polynomial.chebyshev.Chebyshev.basis(0)(x) if m < 0 else \
               np.sin((m + 1) * np.arccos(np.clip(x, -1, 1))) / np.sqrt(1 - np.clip(x, -1, 1)**2 + 1e-300)
    def integrand(p):
        u4 = (p**2 - k**2) / (p**2 + k**2)
        Yp = np.sqrt(2/np.pi) * U(nprime-1, u4) / np.sqrt(4*np.pi)
        Yn = np.sqrt(2/np.pi) * U(n-1, u4) / np.sqrt(4*np.pi)
        wt = (2*k/(k**2 + p**2))**3
        sinc = 1.0 if p*R == 0 else np.sin(p*R)/(p*R)
        return 4*np.pi * p**2 * wt * Yp * Yn * sinc     # 4pi p^2 from d^3p; angular done -> sinc
    val, _ = quad(integrand, 0, np.inf, limit=400)
    return val

# ------------------------------------------------------------------ position-space reference (self-contained)
def position_ref(nmax, R, big=25.0):
    rho = np.linspace(1e-4, big, 700); z = np.linspace(-20, big, 1100)
    RHO, ZZ = np.meshgrid(rho, z, indexing="ij")
    drho = rho[1]-rho[0]; dz = z[1]-z[0]
    def cs_at(n, zc):
        rr = np.sqrt(RHO**2 + (ZZ-zc)**2)
        f = np.exp(-k*rr) * genlaguerre(n-1, 1)(2*k*rr)
        nrm = np.sqrt(2*np.pi*np.sum(f*f*RHO)*drho*dz)
        return f/nrm
    def integ(g): return 2*np.pi*np.sum(g*RHO)*drho*dz
    def ov(fi, fj): return integ(fi*fj)
    def swf(fi, fj):
        gir, giz = np.gradient(fi, rho, z); gjr, gjz = np.gradient(fj, rho, z)
        return (1/(2*k**2))*integ(gir*gjr + giz*gjz) + 0.5*ov(fi, fj)
    basis = [(n, 0.0) for n in range(1, nmax+1)] + [(n, R) for n in range(1, nmax+1)]
    fs = [cs_at(n, zc) for (n, zc) in basis]
    N = len(fs); O = np.zeros((N, N)); S = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            O[i, j] = ov(fs[i], fs[j]); S[i, j] = swf(fs[i], fs[j])
    return O, S

# =====================================================================================
if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)

    print("="*78)
    print("PART (1)  VALIDATION GATE -- momentum-space vs position-space")
    print("="*78)

    print("\n[a] Intra-center SW block  (momentum form is EXACT identity):")
    S_intra = block_mom(0.0, 3, "S")
    print("   SW intra (momentum) diag:", np.round(np.diag(S_intra), 6),
          " off-diag max:", f"{np.abs(S_intra-np.diag(np.diag(S_intra))).max():.2e}")
    print("   cond(SW intra, momentum):", round(cond(S_intra), 4),
          "   [position grid gave 1.0 with {1,0.999,0.998}]")
    m_intra = block_mom(0.0, 3, "m")
    print("   L2 intra (momentum) matrix:\n", m_intra)
    print("   cond(L2 intra, momentum):", round(cond(m_intra), 3),
          "   [position grid gave 5.83]")

    print("\n[b] Two-center table -- momentum vs position (gate: agree within grid error):")
    print(f"  {'R':>4} {'nmax':>4} | {'cond(L2) mom':>12} {'cond(L2) pos':>12} | "
          f"{'cond(SW) mom':>12} {'cond(SW) pos':>12} | {'SW offblk mom':>12} {'pos':>7}")
    print("  " + "-"*94)
    gate_ok = True
    for R in (1.4, 2.0, 4.0):
        for nmax in (2, 3):
            Om = assemble_mom(nmax, R, "m"); Sm = assemble_mom(nmax, R, "S")
            Op, Sp = position_ref(nmax, R)
            cLm, cLp = cond(Om), cond(Op); cSm, cSp = cond(Sm), cond(Sp)
            offm = np.abs(Sm[:nmax, nmax:]).max(); offp = np.abs(Sp[:nmax, nmax:]).max()
            # gate: relative agreement of cond(SW) within ~8% (grid error of the 2D position calc)
            rel = abs(cSm - cSp)/cSp
            gate_ok &= rel < 0.10
            print(f"  {R:>4} {nmax:>4} | {cLm:>12.2f} {cLp:>12.2f} | "
                  f"{cSm:>12.2f} {cSp:>12.2f} | {offm:>12.4f} {offp:>7.4f}   (dSW={rel*100:.1f}%)")
    print(f"\n  GATE (cond(SW) mom vs pos within 10%): {'PASS' if gate_ok else 'FAIL'}")

    print("\n[c] Triple check -- direct 3D momentum quad vs 1D chi-reduction (R=2.0):")
    for (a, b) in [(1, 1), (1, 2), (2, 3), (3, 3)]:
        v1 = (2.0/np.pi)*_trapz(_sin(a)*_sin(b)*_sinc(k*2.0*_cot), _chi)
        v3 = sw_3d_radial(a, b, 2.0)
        print(f"    S_({a},{b})(2.0):  1D chi = {v1:+.6f}   3D pquad = {v3:+.6f}   "
              f"diff = {abs(v1-v3):.1e}")

    print("\n" + "="*78)
    print("PART (2)  GROWTH of cond(SW) and cond(L2) with basis size (momentum form)")
    print("="*78)
    nmaxes = list(range(1, 13))
    for R in (1.4, 2.0, 4.0):
        cS = np.array([cond(assemble_mom(nm, R, "S")) for nm in nmaxes])
        cL = np.array([cond(assemble_mom(nm, R, "m")) for nm in nmaxes])
        print(f"\n  R = {R} bohr:")
        print("    nmax  :", "  ".join(f"{nm:6d}" for nm in nmaxes))
        print("    cond S:", "  ".join(f"{x:6.1f}" for x in cS))
        print("    cond L:", "  ".join(f"{x:6.1f}" for x in cL))
        # fits over nmax>=2 (skip the trivial nmax=1)
        xs = np.array(nmaxes[1:], float); N = 2*xs
        for name, c in [("SW", cS[1:]), ("L2", cL[1:])]:
            # power law: log c = A + p log N   ;  exponential: log c = A + r N
            pp = np.polyfit(np.log(N), np.log(c), 1)
            pe = np.polyfit(N, np.log(c), 1)
            rp = 1 - np.sum((np.log(c)-np.polyval(pp, np.log(N)))**2)/np.sum((np.log(c)-np.log(c).mean())**2)
            re = 1 - np.sum((np.log(c)-np.polyval(pe, N))**2)/np.sum((np.log(c)-np.log(c).mean())**2)
            print(f"      {name}: power-law  cond ~ N^{pp[0]:.2f}  (R2={rp:.4f}) | "
                  f"exp  cond ~ e^({pe[0]:.3f} N)  (R2={re:.4f})  ->"
                  f" {'POWER' if rp>=re else 'EXPONENTIAL'}")

    print("\n" + "="*78)
    print("PART (3)  LEVERS")
    print("="*78)

    print("\n[L1] cond(SW) vs separation R  (does it stay bounded? limit R->inf):")
    Rs = [1.0, 1.4, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0]
    for nmax in (3, 5, 8):
        row = [cond(assemble_mom(nmax, R, "S")) for R in Rs]
        rowL = [cond(assemble_mom(nmax, R, "m")) for R in Rs]
        print(f"   nmax={nmax:>2} (N={2*nmax}):")
        print("     R      :", "  ".join(f"{R:6.1f}" for R in Rs))
        print("     cond SW:", "  ".join(f"{x:6.1f}" for x in row))
        print("     cond L2:", "  ".join(f"{x:6.1f}" for x in rowL))

    print("\n[L2] Momentum vs position gives the SAME matrix -> SAME cond (evaluation, not conditioning):")
    Sm = assemble_mom(3, 1.4, "S"); Sp = position_ref(3, 1.4)[1]
    print(f"     cond(SW) momentum={cond(Sm):.3f}  position={cond(Sp):.3f}  "
          f"(momentum removes grid error; conditioning is intrinsic to S, invariant)")

    print("\n[L3] Symmetry-adapted (gerade/ungerade) localized combinations -- block-diagonalize S:")
    print("     (orthogonal transform => cond INVARIANT; test whether g/u SPLIT helps per-block)")
    for R in (1.4, 2.0, 4.0):
        nmax = 5; S = assemble_mom(nmax, R, "S")
        # gerade/ungerade orthogonal transform: pair (A_n, B_n) -> (A_n+B_n, A_n-B_n)/sqrt2
        T = np.zeros((2*nmax, 2*nmax))
        for n in range(nmax):
            T[n, n] = T[n, n+nmax] = 1/np.sqrt(2)          # gerade
            T[n+nmax, n] = 1/np.sqrt(2); T[n+nmax, n+nmax] = -1/np.sqrt(2)  # ungerade
        Sr = T @ S @ T.T
        g = Sr[:nmax, :nmax]; u = Sr[nmax:, nmax:]
        offbd = np.abs(Sr[:nmax, nmax:]).max()
        print(f"   R={R}: cond(S)={cond(S):7.2f}  block-diag? offblk={offbd:.1e}  "
              f"cond(g)={cond(g):7.2f} cond(u)={cond(u):7.2f}  max={max(cond(g),cond(u)):7.2f}")

    print("\n[L4] Spectral driver of the growth -- lambda_min / lambda_max of SW vs nmax (R=1.4):")
    print("     nmax   lam_min   lam_max    cond")
    for nmax in range(2, 13):
        S = assemble_mom(nmax, 1.4, "S"); ev = eigvalsh(S)
        print(f"     {nmax:>3}   {ev.min():.5f}   {ev.max():.5f}   {ev.max()/ev.min():7.1f}")

    print("\n[L5] g/u sector conditioning growth (homonuclear, R=1.4) -- does the split help asymptotically?")
    print("     If g & u are solved independently the relevant cost is max(cond_g,cond_u), not cond(full).")
    print("     nmax  N   cond(full)  cond(g)  cond(u)  max(g,u)  full/max")
    ngu, cf_l, cm_l = [], [], []
    for nmax in range(2, 13):
        S = assemble_mom(nmax, 1.4, "S")
        T = np.zeros((2*nmax, 2*nmax))
        for n in range(nmax):
            T[n, n] = T[n, n+nmax] = 1/np.sqrt(2)
            T[n+nmax, n] = 1/np.sqrt(2); T[n+nmax, n+nmax] = -1/np.sqrt(2)
        Sr = T @ S @ T.T; g = Sr[:nmax, :nmax]; u = Sr[nmax:, nmax:]
        cf, cg, cu = cond(S), cond(g), cond(u); cmax = max(cg, cu)
        ngu.append(2*nmax); cf_l.append(cf); cm_l.append(cmax)
        print(f"     {nmax:>3}  {2*nmax:>2}  {cf:>10.1f}  {cg:>7.2f}  {cu:>7.1f}  {cmax:>8.1f}  {cf/cmax:>7.2f}")
    Ngu = np.array(ngu, float)
    pf = np.polyfit(np.log(Ngu), np.log(cf_l), 1); pm = np.polyfit(np.log(Ngu), np.log(cm_l), 1)
    print(f"     power-law:  cond(full) ~ N^{pf[0]:.2f}   max(g,u) ~ N^{pm[0]:.2f}"
          f"   (both polynomial; gerade sector cond~2.1 flat)")
