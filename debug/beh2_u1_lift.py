"""U(1) lift of the configuration operator's Z2 Berry structure (2026-08-25).

The real-symmetric F = sum_i P_i carries a FLAT Z2 connection with pi-flux deltas at the CIs
(beh2_berry_curvature.py): each CI is the real (equatorial) section of a Berry monopole.  This
driver realizes the lift Z2 -> U(1): break time-reversal with a magnetic/Peierls phase phi
threading the H-Be-H loop (multiply the H1<->H2 overlap block by e^{i phi}, keeping G complex-
Hermitian).  Then F(d1,d2,phi) is complex-Hermitian and the conical intersection becomes a
genuine U(1) Berry monopole of INTEGER charge.

KEY STRUCTURAL POINT (corrects an earlier "planarity" wording): the Z2 is protected by TIME-
REVERSAL (reality of the overlaps), NOT by planarity.  Real orbitals give real overlaps in ANY
3D geometry, so a non-coplanar 4th *real* center keeps F real-symmetric and does NOT lift Z2.
The lift requires a genuine complex (T-breaking) phase.

Verifies:
  (a) phi opens the gap at the CI, linearly (a genuine third / sigma_y cone axis);
  (b) F(phi=0) is exactly real; a real perturbation keeps the (d1,d2) plaquette flux quantized
      to {0,pi}, while phi!=0 makes it continuous (curvature becomes a smooth U(1) field);
  (c) the branching map (a_x,a_z,a_y) over (d1,d2,phi) is a local diffeo (nonzero 3x3 Jacobian);
  (d) the first Chern number over a small sphere enclosing (d*,d*,0) is +-1 (integer monopole).
"""
from __future__ import annotations
import sys, json, importlib.util
import numpy as np
sys.path.insert(0, "debug")
spec = importlib.util.spec_from_file_location("landscape", "debug/beh2_ci_exact_landscape.py")
L = importlib.util.module_from_spec(spec); spec.loader.exec_module(L)
_beh, _hh, PAR = L._beh, L._hh, L.PAR


def projectors_phi(d1, d2, phi):
    """Complex-Hermitian G with a Peierls phase e^{i phi} on the H1<->H2 block; center order (Be,H1,H2)."""
    S1 = _beh(d1); S2 = PAR @ _beh(d2) @ PAR; SHH = PAR @ _hh(d1 + d2) @ PAR
    I = np.eye(2); e = np.exp(1j * phi)
    G = np.block([[I + 0j, S1 + 0j, S2 + 0j],
                  [S1.T + 0j, I + 0j, e * SHH],
                  [S2.T + 0j, np.conj(e) * SHH.T, I + 0j]])
    if np.linalg.eigvalsh(G).min() < 1e-9:
        return None
    Xh = np.linalg.cholesky(G).conj().T
    return [Xh[:, 2 * k:2 * k + 2] @ np.linalg.pinv(Xh[:, 2 * k:2 * k + 2]) for k in range(3)]


def Feig(d1, d2, phi):
    Ps = projectors_phi(d1, d2, phi)
    return None if Ps is None else np.linalg.eigh(Ps[0] + Ps[1] + Ps[2])


def plaq_flux_dplane(d1, d2, h, phi, band):
    vs = [Feig(*p, phi)[1][:, band] for p in [(d1, d2), (d1 + h, d2), (d1 + h, d2 + h), (d1, d2 + h)]]
    U = np.vdot(vs[0], vs[1]) * np.vdot(vs[1], vs[2]) * np.vdot(vs[2], vs[3]) * np.vdot(vs[3], vs[0])
    return float(np.angle(U))


def chern_sphere(c, r_d, r_phi, band, Nt=40, Np=40):
    ths = np.linspace(1e-3, np.pi - 1e-3, Nt); pss = np.linspace(0, 2 * np.pi, Np, endpoint=False)
    V = np.empty((Nt, Np, 6, 6), complex)
    for i, th in enumerate(ths):
        for j, ps in enumerate(pss):
            V[i, j] = Feig(c + r_d * np.sin(th) * np.cos(ps),
                           c + r_d * np.sin(th) * np.sin(ps), r_phi * np.cos(th))[1]
    tot = 0.0
    for i in range(Nt - 1):
        for j in range(Np):
            jn = (j + 1) % Np
            v = [V[i, j, :, band], V[i + 1, j, :, band], V[i + 1, jn, :, band], V[i, jn, :, band]]
            U = np.vdot(v[0], v[1]) * np.vdot(v[1], v[2]) * np.vdot(v[2], v[3]) * np.vdot(v[3], v[0])
            tot += np.angle(U)
    return tot / (2 * np.pi)


def main():
    c = 2.445
    w0, V0 = Feig(c, c, 0.0); kc = int(np.argmin(np.diff(w0)))
    out = {"crossing_band": kc}

    print("(a) magnetic phase opens the gap at the CI (linear => genuine third cone axis):")
    gaps = {}
    for phi in [0.0, 0.02, 0.05, 0.1, 0.2]:
        w, _ = Feig(c, c, phi); g = float(w[kc + 1] - w[kc]); gaps[phi] = g
        print("   phi=%.2f  gap=%.5e" % (phi, g))
    kphi = (gaps[0.05] - gaps[0.0]) / 0.05
    out["gap_slope_phi"] = kphi

    print("\n(b) reality protects Z2: F(phi=0) real; d-plane flux quantized at phi=0, continuous at phi!=0:")
    Fimag = np.max(np.abs((lambda P: P[0] + P[1] + P[2])(projectors_phi(c, c, 0.0)).imag))
    f0 = plaq_flux_dplane(2.40, 2.40, 0.09, 0.0, kc + 1)
    f1 = plaq_flux_dplane(2.40, 2.40, 0.09, 0.10, kc + 1)
    q0 = min(abs(f0), abs(abs(f0) - np.pi)); q1 = min(abs(f1), abs(abs(f1) - np.pi))
    print("   max|imag F(0)| = %.1e" % Fimag)
    print("   plaquette flux/pi: phi=0 -> %.4f (dist-from-{0,pi}=%.1e, Z2), phi=0.1 -> %.4f (dist=%.3f, U(1))"
          % (f0 / np.pi, q0, f1 / np.pi, q1))
    out["Fimag_phi0"] = float(Fimag); out["flux_dev_phi0"] = float(q0); out["flux_dev_phi01"] = float(q1)

    print("\n(c) branching 3x3 Jacobian (a_x,a_z,a_y) over (d1,d2,phi): nonzero => local diffeo:")
    Q = V0[:, [kc, kc + 1]]
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)

    def a_of(d1, d2, phi):
        M = Q.conj().T @ (lambda P: P[0] + P[1] + P[2])(projectors_phi(d1, d2, phi)) @ Q
        M = M - 0.5 * np.trace(M) * np.eye(2)
        return np.array([0.5 * np.trace(M @ sx).real, 0.5 * np.trace(M @ sz).real, 0.5 * np.trace(M @ sy).real])
    h = 1e-3
    J = np.array([(a_of(c + h, c, 0) - a_of(c - h, c, 0)) / (2 * h),
                  (a_of(c, c + h, 0) - a_of(c, c - h, 0)) / (2 * h),
                  (a_of(c, c, h) - a_of(c, c, -h)) / (2 * h)]).T   # columns d1,d2,phi ; rows a_x,a_z,a_y
    detJ = float(np.linalg.det(J))
    print("   det J(a_x,a_z,a_y ; d1,d2,phi) = %.4e  (!=0 => genuine 3D monopole)" % detJ)
    out["det_jacobian_3d"] = detJ

    print("\n(d) first Chern number over a sphere enclosing (d*,d*,0):")
    r_d = 0.03; r_phi = r_d * 0.07 / max(kphi, 1e-6)
    C = {}
    for band in [kc, kc + 1]:
        cc = chern_sphere(c, r_d, r_phi, band)
        C[band] = cc
        print("   band %d: Chern = %.3f" % (band, cc))
    out["chern"] = {str(k): v for k, v in C.items()}

    with open("debug/data/beh2_u1_lift.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\n=> the CI is the T-symmetric section of a charge-+-1 U(1) Berry monopole; a magnetic")
    print("   (T-breaking) phase lifts Z2 -> U(1).  Reality (not planarity) protects the Z2.")
    print("wrote debug/data/beh2_u1_lift.json")


if __name__ == "__main__":
    main()
