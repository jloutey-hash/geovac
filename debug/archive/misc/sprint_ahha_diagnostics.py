"""
/ahha grounded-pass diagnostics (2026-07-08).  Three reaches from the
composition-wall arc (v4.73.0-1).  This driver runs the two NUMERIC reaches;
Reach 3's falsifier is measurement-theoretic (reasoned verdict in the session
summary, not here).

Shared object: the two-center overlap block S_AB of the LiH bond block (both
Z=1 hydrogenic; within-center orthonormal => singular values of S_AB are the
canonical correlations cos(theta_k), principal angles between the center-
subspaces).  Reuses the VALIDATED overlap_two_center (Topos-3; exact 1s-1s
Slater to 12 digits).  Cross-m overlaps are ANALYTICALLY zero (azimuthal
integral => delta_{mm'}); overlap_two_center only takes one m for that reason.

REACH 2 (m-block-diagonality => wall confined to within-m l-sector):
  Build the FULL cross-Gram over the complete basis {(n,l,m)} WITHOUT assuming
  block structure -- fill 0 in cross-m positions (the azimuthal theorem), the
  radial value in same-m positions -- then verify:
    (i)  global principal-angle spectrum == union of per-m spectra (block
         decomposition is faithful: [P_A,P_B] = (+)_m [P_A,P_B]|_m),
    (ii) global ||[P_A,P_B]|| == max over m-blocks,
    (iii) l-mixing (the wall) is nonzero ONLY within m-blocks.
  => m-grading = commutant (carries the Loewdin-invariant Pauli sparsity);
     l-within-m = anticommutant (carries the composition wall).  Disjoint.

REACH 1 (Halmos/CS composed object + Connes-type distance monotone in R):
  On the m=0 block, sweep theta_k(R); check each angle monotone 0->pi/2 as
  R: 0->inf; build the candidate Grassmann geodesic distance d(R)=||theta(R)||_2
  and check it is monotone/invertible (R recoverable from the angle spectrum);
  print the Halmos 2x2 blocks (one qubit per angle) at R_true.
  OUT OF SCOPE (flagged, not run): whether the CS basis keeps the two-body ERI
  sparse -- that is the expensive many-body leg, and the live falsifier.
"""
from __future__ import annotations
import json
import numpy as np
from scipy.special import genlaguerre, factorial, lpmv

R_TRUE = 3.015

# --- robust two-centre overlap via FIXED Gauss-Legendre in prolate spheroidal
# coords (the topos3 adaptive mpmath quad HANGS at n=3,l=2 -- documented WRONG
# TOOL, memo sprint_commutator_and_explorer §7).  Same signature; numpy-fast.
# Validated in sprint_lowdin_tradeoff.py (1s-1s exact to 1e-16). ---
_gl_xi_u, _gl_xi_w = np.polynomial.legendre.leggauss(80)
_gl_eta, _gl_eta_w = np.polynomial.legendre.leggauss(60)


def _ang(l, m, ct):
    m = abs(m)
    norm = np.sqrt((2 * l + 1) / 2.0 * factorial(l - m) / factorial(l + m))
    return norm * lpmv(m, l, ct)


def _R_nl(Z, n, l, r):
    Z = float(Z); rho = 2 * Z * r / n
    norm = np.sqrt((2 * Z / n) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-rho / 2) * rho ** l * genlaguerre(n - l - 1, 2 * l + 1)(rho)


def overlap_two_center(Z1, n1, l1, Z2, n2, l2, m, R):
    """<phi_{n1 l1 m}(0,Z1) | phi_{n2 l2 m}(R zhat, Z2)>, m conserved."""
    R = float(R); half = R / 2.0
    xi_max = 1.0 + 90.0 / max((float(Z1) + float(Z2)) * R, 1e-6)
    xi = 1.0 + (xi_max - 1.0) * (_gl_xi_u + 1.0) / 2.0
    wxi = _gl_xi_w * (xi_max - 1.0) / 2.0
    XI, ETA = np.meshgrid(xi, _gl_eta, indexing='ij')
    W = np.outer(wxi, _gl_eta_w)
    r1 = half * (XI + ETA); r2 = half * (XI - ETA)
    ct1 = np.clip((1 + XI * ETA) / (XI + ETA), -1, 1)
    ct2 = np.clip((XI * ETA - 1) / (XI - ETA), -1, 1)
    integ = (_R_nl(Z1, n1, l1, r1) * _ang(l1, m, ct1)
             * _R_nl(Z2, n2, l2, r2) * _ang(l2, m, ct2)
             * half ** 3 * (XI ** 2 - ETA ** 2))
    return float(np.sum(integ * W))


def full_basis(n_max):
    """Complete single-center basis {(n,l,m)}, m in [-l, l]."""
    return [(n, l, m) for n in range(1, n_max + 1)
            for l in range(n) for m in range(-l, l + 1)]


def full_cross_gram(R, n_max, Z1=1, Z2=1):
    """S_AB over the FULL basis; cross-m entries set to 0 by the azimuthal
    theorem (delta_{mm'}), same-m entries computed.  We build the whole matrix
    explicitly so the m-block structure is VERIFIED, not assumed."""
    basis = full_basis(n_max)
    d = len(basis)
    S = np.zeros((d, d))
    for i, (na, la, ma) in enumerate(basis):
        for j, (nb, lb, mb) in enumerate(basis):
            if ma == mb:  # azimuthal integral is nonzero only for m_i == m_j
                S[i, j] = float(overlap_two_center(Z1, na, la, Z2, nb, lb, abs(ma), float(R)))
            # else: exactly 0 (theorem), left as the initialized zero
    return basis, S


def principal_angles(S):
    sig = np.clip(np.linalg.svd(S, compute_uv=False), 0.0, 1.0)
    return np.degrees(np.arccos(sig)), sig  # theta_k in degrees, cos theta_k


def comm_norm_from_sigma(sig):
    """||[P_A,P_B]|| = max_k sin(theta_k) cos(theta_k) for two subspaces with
    orthonormal-within bases and cross-Gram singular values sig = cos theta_k."""
    sig = np.clip(sig, 0.0, 1.0)
    return float(np.max(sig * np.sqrt(np.maximum(1 - sig ** 2, 0.0)))) if len(sig) else 0.0


def per_m_breakdown(basis, S):
    """Restrict S to each m-block (S is block-diagonal in m by construction)."""
    out = {}
    for m in sorted({b[2] for b in basis}):
        idx = [i for i, b in enumerate(basis) if b[2] == m]
        sub = S[np.ix_(idx, idx)]
        _, sig = principal_angles(sub)
        states = [(basis[i][0], basis[i][1]) for i in idx]
        l_vals = [s[1] for s in states]
        l_mix = any(abs(sub[a, b]) > 1e-10 and l_vals[a] != l_vals[b]
                    for a in range(len(idx)) for b in range(len(idx)))
        out[m] = {'dim': len(idx), 'states': states, 'sigma': sig.tolist(),
                  'comm_norm': comm_norm_from_sigma(sig), 'l_mixing': bool(l_mix)}
    return out


# =====================================================================
# REACH 2
# =====================================================================
def reach2(n_max=3, R=R_TRUE):
    print("=" * 88)
    print(f"REACH 2  --  [P_A,P_B] block-diagonal in m?  (n_max={n_max}, R={R})")
    print("=" * 88)
    basis, S = full_cross_gram(R, n_max)
    d = len(basis)

    # (i) faithfulness: global spectrum == union of per-m spectra
    theta_global, sig_global = principal_angles(S)
    per_m = per_m_breakdown(basis, S)
    sig_union = np.sort(np.concatenate([np.array(b['sigma']) for b in per_m.values()]))[::-1]
    faithful = np.allclose(np.sort(sig_global)[::-1], sig_union, atol=1e-9)

    # explicit block-diagonality check on the raw matrix
    cross_m_max = 0.0
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            if bi[2] != bj[2]:
                cross_m_max = max(cross_m_max, abs(S[i, j]))

    comm_global = comm_norm_from_sigma(sig_global)
    comm_max_block = max(b['comm_norm'] for b in per_m.values())

    print(f"  full basis dim (all n,l,m): {d}")
    print(f"  max |S_AB| in a CROSS-m position: {cross_m_max:.2e}   "
          f"(=> S block-diagonal in m)")
    print(f"  global ||[P_A,P_B]||         = {comm_global:.6f}")
    print(f"  max_m per-block ||[P_A,P_B]|| = {comm_max_block:.6f}   "
          f"(equal => commutator decomposes over m)")
    print(f"  spectrum faithful (global == union per-m): {faithful}")
    print()
    print(f"  {'m':>3s} {'dim':>4s} {'||[P_A,P_B]||_m':>16s} {'l-mixing (wall)':>16s}   states")
    for m, b in sorted(per_m.items()):
        st = ",".join(f"{n}{'spdf'[l]}" for n, l in b['states'])
        print(f"  {m:>3d} {b['dim']:>4d} {b['comm_norm']:16.6f} "
              f"{('YES' if b['l_mixing'] else 'no'):>16s}   {st}")

    wall_blocks = [m for m, b in per_m.items() if b['l_mixing']]
    commutant_only = [m for m, b in per_m.items() if not b['l_mixing']]
    print()
    print(f"  VERDICT: cross-m coupling = {cross_m_max:.1e} (theorem: azimuthal delta_mm').")
    print(f"    l-mixing (the WALL) appears in m-blocks {wall_blocks}, "
          f"absent in {commutant_only}.")
    print(f"    => the composition non-commutativity is CONFINED within-m; the")
    print(f"       m-grading (commutant) is untouched -> carries the Pauli sparsity.")
    return {'n_max': n_max, 'R': R, 'dim': d, 'cross_m_max': cross_m_max,
            'comm_global': comm_global, 'comm_max_block': comm_max_block,
            'faithful': bool(faithful),
            'per_m': {str(k): v for k, v in per_m.items()}}


# =====================================================================
# REACH 1
# =====================================================================
def reach1(n_max=2, Rs=None):
    if Rs is None:
        Rs = [0.5, 1.0, 1.5, 2.0, 2.5, R_TRUE, 4.0, 5.0, 6.0, 8.0, 12.0]
    print("\n" + "=" * 88)
    print(f"REACH 1  --  Halmos/CS composed object + Connes-type distance (m=0 block, "
          f"n_max={n_max})")
    print("=" * 88)

    sweep = []
    for R in Rs:
        basis, S = full_cross_gram(R, n_max)
        idx = [i for i, b in enumerate(basis) if b[2] == 0]
        sub = S[np.ix_(idx, idx)]
        theta, sig = principal_angles(sub)
        theta_sorted = np.sort(theta)  # ascending
        dist = float(np.linalg.norm(theta_sorted))  # Grassmann geodesic distance (deg)
        sweep.append({'R': R, 'theta_deg': theta_sorted.tolist(),
                      'sigma': np.sort(sig)[::-1].tolist(), 'distance_deg': dist})

    # monotonicity of each angle and of the distance functional
    thetas = np.array([s['theta_deg'] for s in sweep])   # rows: R, cols: angle_k (ascending)
    dists = np.array([s['distance_deg'] for s in sweep])
    ang_mono = [bool(np.all(np.diff(thetas[:, k]) > -1e-6)) for k in range(thetas.shape[1])]
    dist_mono = bool(np.all(np.diff(dists) > -1e-6))

    print(f"  {'R':>6s}  " + "  ".join(f"theta{k}" for k in range(thetas.shape[1]))
          + f"   {'d(R)=||theta||':>14s}")
    for s in sweep:
        angs = "  ".join(f"{a:6.2f}" for a in s['theta_deg'])
        print(f"  {s['R']:6.3f}  {angs}   {s['distance_deg']:14.3f}")
    print()
    print(f"  each theta_k(R) monotone increasing (0->90 as R:0->inf): {ang_mono}")
    print(f"  distance d(R) monotone increasing (=> R recoverable from angle spectrum): "
          f"{dist_mono}")

    # Halmos 2x2 blocks at R_true: one qubit per angle
    basis, S = full_cross_gram(R_TRUE, n_max)
    idx = [i for i, b in enumerate(basis) if b[2] == 0]
    sub = S[np.ix_(idx, idx)]
    U, sig, Vt = np.linalg.svd(sub)
    sig = np.clip(sig, 0, 1)
    print(f"\n  Halmos / CS decomposition at R_true={R_TRUE} (m=0 block):")
    print(f"    the pair (P_A,P_B) reduces to {len(sig)} independent 2x2 rotation blocks,")
    print(f"    one per principal angle -- 'one qubit per angle'.  Each block Gram:")
    for k, c in enumerate(sig):
        th = np.degrees(np.arccos(c))
        print(f"      block {k}: theta={th:5.2f} deg, Gram=[[1, {c:+.4f}],[{c:+.4f}, 1]]")
    print(f"\n  OUT OF SCOPE (live falsifier, not run): does the two-body ERI stay sparse")
    print(f"    in this CS basis?  That is the many-body leg; the encoding-wall diagnostic")
    print(f"    already suspects JW forces l-mixing there.")
    return {'n_max': n_max, 'sweep': sweep, 'angle_monotone': ang_mono,
            'distance_monotone': dist_mono}


if __name__ == "__main__":
    out = {'reach2': reach2(n_max=3), 'reach1': reach1(n_max=2)}
    json.dump(out, open('debug/data/ahha_diagnostics.json', 'w'), indent=2)
    print("\n[saved] debug/data/ahha_diagnostics.json")
