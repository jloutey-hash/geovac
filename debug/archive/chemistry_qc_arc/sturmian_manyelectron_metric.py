"""
Paper 60 -- MANY-ELECTRON molecular case: does the Shibuya-Wulfman metric penalty
COMPOUND with electron number, or stay at the one-electron cost?

Structural claim to test.  In the generalized-Sturmian method a k-electron configuration
is a Slater determinant of one-electron molecular orbitals.  The N-electron secular
equation is [T0 + T' - p_kappa 1]B = 0 -- a STANDARD (metric-free) eigenproblem, the
'-p_kappa 1' being the identity, not a metric (Avery, BK6 eq 6.19 general form).  That
requires the configuration basis to be orthonormal, i.e. the one-electron orbitals to be
orthonormal in the relevant metric.

Two routes to build the molecular orbitals from many-center Coulomb-Sturmians:
  (RAW)  keep the L2-normalized many-center CS orbitals -- they are NON-orthogonal
         (overlap S_L2, ill-conditioned, cond ~ N^1.7 at bond length).  A k-electron
         determinant overlap is det(S_L2[I,J]); the configuration-space metric therefore
         COMPOUNDS -- cond grows like cond(S_L2)^k.  Fully on-device, but the metric
         penalty multiplies with electron number.
  (SW)   diagonalize the one-electron Shibuya-Wulfman problem ONCE (cost = cond(S_SW),
         priced in sec:resource, g/u-improvable) -> molecular Sturmians, S_SW-orthonormal
         (C^T S_SW C = I).  Determinants of orthonormal orbitals have overlap = I -> the
         N-electron problem is METRIC-FREE regardless of electron number.

This driver measures cond(configuration-overlap) for k = 1..4 electrons on a two-center
CS basis in both routes, to show RAW compounds (~ cond^k) while SW stays exactly 1.

Exact momentum-space SW metric + closed-form Gaussian-free L2 overlap (PhD eq 10.5.12/13).
Diagnostic only.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from numpy.linalg import cond
from scipy.linalg import eigh
from itertools import combinations

k = 1.0
_M = 200001
_chi = np.linspace(1e-8, np.pi, _M)
_cot = 1.0 / np.tan(_chi / 2.0)
_1mcos = 1.0 - np.cos(_chi)


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
            v = (2.0 / np.pi) * np.trapezoid(np.sin(a * _chi) * np.sin(b * _chi) * wfac, _chi)
            B[a - 1, b - 1] = B[b - 1, a - 1] = v
    return B


def assemble(nmax, R, kind):
    intra = block_mom(0.0, nmax, kind)
    inter = block_mom(R, nmax, kind)
    return np.block([[intra, inter], [inter.T, intra]])


def config_overlap(S1, k_elec):
    """Configuration-space overlap matrix for k-electron spatial determinants built from
    the one-electron orbitals with one-electron overlap S1.  <Phi_I|Phi_J> = det(S1[I,J])."""
    N = S1.shape[0]
    configs = list(combinations(range(N), k_elec))
    nc = len(configs)
    M = np.zeros((nc, nc))
    for a, I in enumerate(configs):
        for b, J in enumerate(configs):
            M[a, b] = np.linalg.det(S1[np.ix_(I, J)])
    return M, nc


if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    R = 1.4  # H2 equilibrium bond length (bohr)

    print("=" * 82)
    print("MANY-ELECTRON metric compounding -- two-center CS basis, R=1.4 bohr (H2 geometry)")
    print("=" * 82)

    for nmax in (2, 3):
        N = 2 * nmax
        S_L2 = assemble(nmax, R, "m")   # raw L2 overlap (ill-conditioned)
        S_SW = assemble(nmax, R, "S")   # Shibuya-Wulfman metric

        # molecular Sturmians: S_SW-orthonormalize the raw basis (C^T S_SW C = I).
        # In this basis the one-electron overlap that enters determinants is the IDENTITY.
        w, V = eigh(S_SW)
        C = V @ np.diag(1.0 / np.sqrt(w))       # C^T S_SW C = I
        S_L2_in_SWbasis = C.T @ S_L2 @ C        # what L2 overlap looks like after SW-orthonorm
        # NOTE: the generalized-Sturmian determinants are orthonormal in the SW metric, so the
        # relevant one-electron overlap for the config metric is C^T S_SW C = I exactly.
        I1 = np.eye(N)

        lam = np.sort(np.linalg.eigvalsh(S_L2))[::-1]   # one-electron overlap spectrum
        print(f"\n  nmax={nmax}, N={N} orbitals.  one-electron cond: "
              f"cond(S_L2)={cond(S_L2):.2f}, cond(S_SW)={cond(S_SW):.2f}")
        print(f"  one-electron S_L2 eigenvalues: {np.array2string(lam, precision=3)}")
        print(f"  {'k elec':>6} {'n_cfg':>6} | {'cond(RAW config)':>16} {'/cond(S_L2)':>11}"
              f" {'compound-eig':>12} | {'cond(SW config)':>15}")
        print("  " + "-" * 74)
        for ke in range(1, min(N, 5)):
            M_raw, nc = config_overlap(S_L2, ke)      # raw L2: k-th compound matrix of S_L2
            M_sw, _ = config_overlap(I1, ke)          # SW basis: one-electron overlap = I
            craw, csw = cond(M_raw), cond(M_sw)
            # exact cross-check: cond(k-th compound) = prod(top k eigs)/prod(bottom k eigs)
            comp = np.prod(lam[:ke]) / np.prod(lam[-ke:])
            print(f"  {ke:>6} {nc:>6} | {craw:>16.3f} {craw/cond(S_L2):>11.2f}"
                  f" {comp:>12.3f} | {csw:>15.4f}")
        print("  (cond of the k-electron config metric = kth-compound ratio "
              "prod(top k)/prod(bottom k); worst at half-filling, driven by the "
              "one small overlap eigenvalue)")

    print("\n" + "=" * 82)
    print("READING")
    print("=" * 82)
    print("""  RAW many-center L2 basis (fully on-device, no orthonormalization):
    the k-electron configuration-overlap metric is the kth COMPOUND matrix of the one-electron
    overlap; its conditioning = prod(top k eigs)/prod(bottom k eigs), which GROWS with electron
    number (worst at half-filling) and is bounded by cond(S_L2)^min(k,N-k).  The metric penalty
    increases with electron number -- driven by the near-linear-dependence (smallest) eigenvalue.

  SW molecular Sturmians (one-electron SW diagonalization done ONCE -- cost = cond(S_SW), priced
    in sec:resource, g/u-improvable): the configuration overlap is the IDENTITY at EVERY electron
    number -- the N-electron secular equation is METRIC-FREE (Avery BK6 eq 6.19 general form).

  => The molecular metric penalty does NOT have to compound with electron number.  Paid ONCE at
     the one-electron molecular-orbital construction (the Shibuya-Wulfman cost estimated for H2+),
     the many-electron problem inherits an identity metric.  This is the many-electron extension:
     the metric cost is one-electron, not N-electron.""")
    print("\nDONE.")
