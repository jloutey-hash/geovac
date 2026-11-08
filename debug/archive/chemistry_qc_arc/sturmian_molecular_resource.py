"""
Paper 60 -- MOLECULAR resource estimate.

Turns the measured polynomial Shibuya-Wulfman conditioning (cond(S) ~ N^1.85) and the
gerade/ungerade lever into an actual qubit / T-gate count for a quantum eigenvalue
estimation of the one-electron molecular (H2+) isoenergetic secular equation
    [ W - k S ] C = 0,   E = -k^2/2,
posed against the SAME problem in a Gaussian LCAO metric.

COST MODEL (grounded in the paper's cited prior art):
  * A generalized eigenproblem A x = lambda S x with S>0 is whitened to the STANDARD
    eigenproblem  M~ = S^{-1/2} A S^{-1/2},  eigenvalues lambda (Liang et al. 2112.02554).
  * S^{-1/2} is applied by QSVT with a polynomial approximating x^{-1/2} on
    [1/kappa, 1] (kappa = cond(S)); its degree is  d_inv ~ kappa * ln(kappa/eps)
    (Gilyen-Su-Low-Wiebe QSVT; Childs-Kothari-Somma inversion).  This is the
    metric penalty -- it is what "cond(S) multiplies the block-encoding cost" means.
  * Eigenvalue estimation to precision eps_k in the scale k needs
    Q_QPE ~ (pi/2) * lambda_eff / eps_k  queries to the block-encoding (qubitization,
    Low-Chuang 1610.06546), lambda_eff = subnormalization (spectral norm) of M~ = k_max.
  * Each block-encoding query = one BE(A) + one BE(S^{-1/2}); the latter is d_inv calls
    to BE(S).  So total block-encoding calls  ~ Q_QPE * (1 + d_inv).
  * Qubits: n_sys = ceil(log2 N) (eigenvector), n_phase = ceil(log2(1/eps_k)),
    n_anc = LCU-select ancillas of BE(A)/BE(S) + 1 QSVT ancilla.

The metric-free ATOMIC case (Paper 60 sec:iso) has kappa = 1 (no metric): d_inv = 0.
The molecular case pays d_inv(kappa_S).  The g/u lever replaces kappa_S(full) with
max(kappa_gerade, kappa_ungerade), and kappa_gerade ~ 2 is flat -- and the H2+ ground
state 1s_sigma_g lives in the GERADE sector.

Diagnostic only.  Does NOT modify paper_60/tests/CLAUDE/CHANGELOG.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from numpy.linalg import cond

k = 1.0
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

# =====================================================================================
# (A) SHIBUYA-WULFMAN metric S and L2 overlap m -- exact momentum-space 1D forms
#     (identical to debug/sturmian_sw_momentum.py; s-orbitals, two centers along z).
# =====================================================================================
_M = 400001
_chi = np.linspace(1e-8, np.pi, _M)
_cot = 1.0 / np.tan(_chi / 2.0)
_1mcos = 1.0 - np.cos(_chi)
_sinj = {}


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
    sfac = np.ones_like(_chi) if R == 0.0 else _sinc(k * R * _cot)
    wfac = sfac if kind == "S" else _1mcos * sfac
    B = np.zeros((nmax, nmax))
    for a in range(1, nmax + 1):
        for b in range(a, nmax + 1):
            v = (2.0 / np.pi) * _trapz(_sin(a) * _sin(b) * wfac, _chi)
            B[a - 1, b - 1] = B[b - 1, a - 1] = v
    return B


def assemble_S(nmax, R, kind="S"):
    intra = block_mom(0.0, nmax, kind)
    inter = block_mom(R, nmax, kind)
    return np.block([[intra, inter], [inter.T, intra]])


def gu_split(S, nmax):
    T = np.zeros_like(S)
    for n in range(nmax):
        T[n, n] = T[n, n + nmax] = 1 / np.sqrt(2)
        T[n + nmax, n] = 1 / np.sqrt(2)
        T[n + nmax, n + nmax] = -1 / np.sqrt(2)
    Sr = T @ S @ T.T
    return Sr[:nmax, :nmax], Sr[nmax:, nmax:]


# =====================================================================================
# (B) GAUSSIAN metric for the SAME H2+: even-tempered s-Gaussians on two centers.
#     Normalized s-Gaussian overlap:  S = (4ab/(a+b)^2)^{3/4} exp(-ab/(a+b) |A-B|^2).
# =====================================================================================
def gauss_overlap(a, b, dAB):
    return (4 * a * b / (a + b) ** 2) ** 0.75 * np.exp(-a * b / (a + b) * dAB ** 2)


def gaussian_metric(nmax, R, a0=0.10, ratio=2.0):
    """N=2*nmax even-tempered s-Gaussians (nmax per center), exponents a0*ratio^i."""
    exps = a0 * ratio ** np.arange(nmax)
    centers = [0.0] * nmax + [R] * nmax
    allexp = list(exps) + list(exps)
    N = 2 * nmax
    Sg = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            Sg[i, j] = gauss_overlap(allexp[i], allexp[j], abs(centers[i] - centers[j]))
    return Sg


# =====================================================================================
# (C) RESOURCE MODEL
# =====================================================================================
def qsvt_inv_sqrt_degree(kappa, eps):
    """Degree of the QSVT polynomial approximating x^{-1/2} on [1/kappa,1] to error eps.
    Model: d ~ kappa * ln(kappa/eps)  (up to an O(1) constant; kappa=1 => 0, metric-free)."""
    if kappa <= 1.0 + 1e-9:
        return 0
    return int(np.ceil(kappa * np.log(kappa / eps)))


def qpe_queries(lambda_eff, eps_k):
    """Qubitization/QPE queries to estimate an eigenvalue to precision eps_k."""
    return int(np.ceil((np.pi / 2) * lambda_eff / eps_k))


def resource(N, kappa, lambda_eff, eps_k, eps_inv, L_A, L_S):
    """Full resource tuple for a generalized-eig solve with metric cond number kappa."""
    n_sys = int(np.ceil(np.log2(N)))
    n_phase = int(np.ceil(np.log2(1.0 / eps_k)))
    n_anc = int(np.ceil(np.log2(L_A))) + int(np.ceil(np.log2(max(L_S, 1)))) + 1
    d_inv = qsvt_inv_sqrt_degree(kappa, eps_inv)
    Qqpe = qpe_queries(lambda_eff, eps_k)
    be_calls = Qqpe * (1 + d_inv)  # each QPE query = 1 BE(A) + d_inv BE(S)
    return dict(n_sys=n_sys, n_phase=n_phase, n_anc=n_anc,
                n_qubits=n_sys + n_phase + n_anc, d_inv=d_inv, Qqpe=Qqpe, be_calls=be_calls)


# =====================================================================================
if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    R_eq = 2.0  # H2+ equilibrium bond length (bohr)

    print("=" * 80)
    print("PART 1 -- metric conditioning: Shibuya-Wulfman vs Gaussian, same H2+ (R=2 bohr)")
    print("=" * 80)
    print(f"  {'nmax':>4} {'N':>3} | {'cond(SW full)':>13} {'cond(SW ger)':>12} {'cond(SW ung)':>12}"
          f" | {'cond(Gauss)':>11}")
    print("  " + "-" * 74)
    sw_full, sw_max_gu, gauss_c, Ns = [], [], [], []
    for nmax in range(2, 11):
        S = assemble_S(nmax, R_eq, "S")
        g, u = gu_split(S, nmax)
        Sg = gaussian_metric(nmax, R_eq)
        cf, cg, cu, cG = cond(S), cond(g), cond(u), cond(Sg)
        Ns.append(2 * nmax)
        sw_full.append(cf)
        sw_max_gu.append(max(cg, cu))
        gauss_c.append(cG)
        print(f"  {nmax:>4} {2 * nmax:>3} | {cf:>13.1f} {cg:>12.2f} {cu:>12.1f} | {cG:>11.3e}")
    Ns = np.array(Ns, float)

    def fit(c):
        c = np.array(c, float)
        pp = np.polyfit(np.log(Ns), np.log(c), 1)
        pe = np.polyfit(Ns, np.log(c), 1)
        rp = 1 - np.sum((np.log(c) - np.polyval(pp, np.log(Ns))) ** 2) / np.sum((np.log(c) - np.log(c).mean()) ** 2)
        re = 1 - np.sum((np.log(c) - np.polyval(pe, Ns)) ** 2) / np.sum((np.log(c) - np.log(c).mean()) ** 2)
        return pp[0], rp, pe[0], re

    for name, c in [("SW full", sw_full), ("SW max(g,u)", sw_max_gu), ("Gaussian", gauss_c)]:
        pw, rp, ex, re = fit(c)
        kind = "POWER N^%.2f (R2=%.3f)" % (pw, rp) if rp >= re else "EXP e^(%.3f N) (R2=%.3f)" % (ex, re)
        print(f"    {name:>12}: {kind}   [alt: {'exp' if rp >= re else 'power'} R2="
              f"{re if rp >= re else rp:.3f}]")

    print("\n" + "=" * 80)
    print("PART 2 -- resource counts at chemical accuracy (eps_E = 1.6 mHa), H2+ sigma_g")
    print("=" * 80)
    eps_E = 1.6e-3      # chemical accuracy in the energy (Ha)
    k_g = 1.485         # sigma_g scale: E_elec=-1.1026 Ha at R=2 => k=sqrt(2*1.1026)
    eps_k = eps_E / k_g  # dE = -k dk  => eps_k = eps_E/k
    eps_inv = 1e-3       # QSVT inverse accuracy
    lam_eff = 1.72       # subnormalization of M~ = k_max from the H2+ solve (validated)
    print(f"  eps_E={eps_E} Ha  k_sigma_g={k_g}  eps_k={eps_k:.2e}  n_phase=ceil(log2(1/eps_k))"
          f"={int(np.ceil(np.log2(1 / eps_k)))}")
    print(f"  lambda_eff (subnormalization) = k_max = {lam_eff}  (validated H2+ solve)")
    print()
    print(f"  {'nmax':>4} {'N':>3} | {'metric/method':>16} {'kappa':>8} {'d_inv':>6} {'Q_QPE':>7}"
          f" {'BE calls':>10} | {'qubits':>7}")
    print("  " + "-" * 80)
    for nmax in (4, 6, 8):
        N = 2 * nmax
        S = assemble_S(nmax, R_eq, "S")
        g, u = gu_split(S, nmax)
        Sg = gaussian_metric(nmax, R_eq)
        rows = [
            ("ATOMIC (no metric)", 1.0, N),
            ("SW full", cond(S), N),
            ("SW gerade (sig_g)", cond(g), nmax),  # g/u split halves the register
            ("Gaussian LCAO", cond(Sg), N),
        ]
        for label, kap, Ndim in rows:
            r = resource(max(Ndim, 2), kap, lam_eff, eps_k, eps_inv,
                         max(Ndim, 2) ** 2, max(Ndim, 2) ** 2 if kap > 1 else 1)
            print(f"  {nmax:>4} {N:>3} | {label:>16} {kap:>8.1f} {r['d_inv']:>6} {r['Qqpe']:>7}"
                  f" {r['be_calls']:>10} | {r['n_qubits']:>7}")
        print()

    print("=" * 80)
    print("PART 3 -- HEAD-TO-HEAD block-encoding cost ratio (SW vs Gaussian), same accuracy")
    print("=" * 80)
    print("  Cost proxy = BE calls = Q_QPE*(1+d_inv);  Q_QPE identical across metrics, so")
    print("  ratio = (1+d_inv[Gauss]) / (1+d_inv[SW]).  g/u uses the gerade sector for sigma_g.")
    print(f"  {'nmax':>4} {'N':>3} | {'d_inv SW-full':>12} {'d_inv SW-ger':>12} {'d_inv Gauss':>12}"
          f" | {'Gauss/SWfull':>12} {'Gauss/SWger':>11}")
    print("  " + "-" * 80)
    for nmax in range(2, 11):
        S = assemble_S(nmax, R_eq, "S")
        g, _ = gu_split(S, nmax)
        Sg = gaussian_metric(nmax, R_eq)
        d_sw = qsvt_inv_sqrt_degree(cond(S), eps_inv)
        d_g = qsvt_inv_sqrt_degree(cond(g), eps_inv)
        d_ga = qsvt_inv_sqrt_degree(cond(Sg), eps_inv)
        rf = (1 + d_ga) / (1 + d_sw)
        rg = (1 + d_ga) / (1 + d_g)
        print(f"  {nmax:>4} {2 * nmax:>3} | {d_sw:>12} {d_g:>12} {d_ga:>12} | {rf:>12.1f} {rg:>11.1f}")

    print("\n" + "=" * 80)
    print("PART 4 -- large-R lever (cond -> 1 for SW, Gaussian cannot follow)")
    print("=" * 80)
    print(f"  {'R':>5} | {'cond(SW,N=10)':>13} {'d_inv SW':>9} | {'cond(Gauss,N=10)':>16} {'d_inv G':>8}")
    print("  " + "-" * 66)
    for R in (1.4, 2.0, 3.0, 4.0, 6.0, 10.0):
        S = assemble_S(5, R, "S")
        Sg = gaussian_metric(5, R)
        cS, cG = cond(S), cond(Sg)
        print(f"  {R:>5.1f} | {cS:>13.2f} {qsvt_inv_sqrt_degree(cS, eps_inv):>9} |"
              f" {cG:>16.3e} {qsvt_inv_sqrt_degree(cG, eps_inv):>8}")
    print("\nDONE.")
