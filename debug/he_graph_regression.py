"""He graph-native CI regression (2026-09-25, PI-approved re-pricing).

`casimir_ci._build_graph_h1` (used by `build_graph_native_fci`) sets
h1 = diag(-Z^2/2n^2) + (1/16)*A  (A = GeometricLattice adjacency).
In the k = Z hydrogenic basis the physical one-body operator T + V_ne is
EXACTLY diagonal, so the (1/16)*A "hop" is not part of He's Hamiltonian.

This driver measures everything the corpus rests on, WITH the hop (the
published construction) and WITHOUT it (the physical Hamiltonian in the same
basis, same V_ee at k = Z):

  A. He ground state, n_max = 2..7            (headline 0.216%)
  B. Z-scan: H^- (Z=1), Z_c crossing, Z = 10  (Paper 13 / CHANGELOG CUSP-2)
  C. CUSP-2 Schwartz fit intercept            (ledger row CUSP-2)
  D. Paper 26 Table I (n_max=4) + S ~ Z^alpha
  E. Paper 27 Table ep1 (n_max=2,3) + EP-2c table/fits (n_max=3)
  F. candidate definitions of Paper 13's "~86% of correlation"

V_ee config matrices come from he_graph_regression_vee.py (k = 1 cache;
V_ee is exactly linear in k).  Recomposition H(Z) = Z^2 D1 + Z V1 + hop is
validated against build_graph_native_fci and build_fci_matrix before use.

Run:  python debug/he_graph_regression.py > debug/data/he_graph_regression.log
"""
import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import geovac.casimir_ci as C  # noqa: E402
from geovac.lattice import GeometricLattice  # noqa: E402
from debug.archive.misc.entanglement_geometry import (  # noqa: E402
    build_1rdm_from_singlet_ci, compute_entanglement_measures,
)

# Nonrelativistic infinite-mass He-like ground states (Drake / Nakashima-Nakatsuji)
E_EXACT = {
    1: -0.527751016544, 2: -2.903724377034, 3: -7.279913412669,
    4: -13.655566238423, 5: -22.030971580243, 6: -32.406246601898,
    7: -44.781445148772, 8: -59.156595122757, 9: -75.531712363960,
    10: -93.906806515037,
}
# 1/Z expansion coefficients (Baker, Freund, Hill, Morgan 1990): E = sum e_k Z^(2-k)
E_SERIES = [-1.0, 0.625, -0.157666405, 0.008698991, -0.000888587,
            -0.001036372, -0.000612917, -0.000372187, -0.000242004]


def e_series(Z: float) -> float:
    return float(sum(c * Z ** (2 - k) for k, c in enumerate(E_SERIES)))


# ---------------------------------------------------------------------------
# Config-space building blocks (singlet, M_L = 0, lattice.states ordering --
# identical to build_graph_native_fci's default enumeration)
# ---------------------------------------------------------------------------
class Blocks:
    def __init__(self, n_max: int):
        lat = GeometricLattice(max_n=n_max)
        self.orbitals = list(lat.states)
        self.n_spatial = lat.num_states
        A = lat.adjacency
        A = A.toarray() if hasattr(A, "toarray") else np.array(A)
        d1 = np.diag([-1.0 / (2.0 * n * n) for (n, l, m) in self.orbitals])
        hop = np.where(np.abs(A) > 1e-15, (-1.0 / 16.0) * (-A), 0.0)
        np.fill_diagonal(hop, 0.0)
        self.configs = [(i, j) for i in range(self.n_spatial)
                        for j in range(i, self.n_spatial)
                        if self.orbitals[i][2] + self.orbitals[j][2] == 0]
        self.D1 = self._one_body(d1)       # multiply by Z^2
        self.HOP = self._one_body(hop)     # Z-independent
        self.V1 = np.load(f"debug/data/he_graph_regression_V1_n{n_max}.npz")["V1"]
        assert self.V1.shape == self.D1.shape

    def _one_body(self, h: np.ndarray) -> np.ndarray:
        cf = self.configs
        n = len(cf)
        H = np.zeros((n, n))
        for I, (i, j) in enumerate(cf):
            bra = [(i, j)] + ([(j, i)] if i != j else [])
            for J in range(I, n):
                p, q = cf[J]
                ket = [(p, q)] + ([(q, p)] if p != q else [])
                me = 0.0
                for a, b in bra:
                    for c, d in ket:
                        if b == d:
                            me += h[a, c]
                        if a == c:
                            me += h[b, d]
                me /= np.sqrt(len(bra) * len(ket))
                H[I, J] = H[J, I] = me
        return H

    def H(self, Z: float, hop: bool) -> np.ndarray:
        H = Z * Z * self.D1 + Z * self.V1
        return H + self.HOP if hop else H

    def entropy(self, ci: np.ndarray) -> float:
        rho = build_1rdm_from_singlet_ci(ci, self.configs, self.n_spatial)
        return float(compute_entanglement_measures(rho)["von_neumann_entropy"])


def ground(H: np.ndarray):
    w, v = np.linalg.eigh(H)
    return float(w[0]), v[:, 0]


def pct(E: float, Z: int = 2) -> float:
    return 100.0 * (E - E_EXACT[Z]) / abs(E_EXACT[Z])


def fit(xs, ys):
    lx, ly = np.log(np.asarray(xs, float)), np.log(np.asarray(ys, float))
    a, b = np.polyfit(lx, ly, 1)
    r2 = 1 - np.sum((ly - (a * lx + b)) ** 2) / np.sum((ly - ly.mean()) ** 2)
    return float(a), float(np.exp(b)), float(r2)


def main() -> None:
    avail = [n for n in range(2, 8)
             if os.path.exists(f"debug/data/he_graph_regression_V1_n{n}.npz")]
    B = {n: Blocks(n) for n in avail}
    print(f"available n_max: {avail}")

    # ---------------- validation ------------------------------------------
    print("\n=== VALIDATION ===")
    for n in [m for m in avail if m <= 4]:
        H_pub = C.build_graph_native_fci(Z=2, n_max=n)
        d = np.max(np.abs(H_pub - B[n].H(2, True)))
        E_fm = float(np.linalg.eigvalsh(C.build_fci_matrix(Z=2, n_max=n, k_orb=2.0))[0])
        E_ph = ground(B[n].H(2, False))[0]
        print(f"n_max={n}: max|H_pub - recomposed(hop)| = {d:.2e};  "
              f"E_phys(recomposed) - E(build_fci_matrix k=Z) = {E_ph - E_fm:+.2e}")
        assert d < 1e-12 and abs(E_ph - E_fm) < 1e-10

    # ---------------- A. He ladder -----------------------------------------
    print("\n=== A. He (Z=2) ground state, singlet M_L=0 ===")
    print("n_max  dim   E_hop        %      E_phys       %      <psi_hop|H_phys|psi_hop>  %      hop(mHa)")
    rowsA = {}
    for n in avail:
        Eh, ph = ground(B[n].H(2, True))
        Ep, pp = ground(B[n].H(2, False))
        Em = float(ph @ B[n].H(2, False) @ ph)
        rowsA[n] = (Eh, Ep, Em, ph, pp)
        print(f"{n:4d} {B[n].D1.shape[0]:5d}  {Eh:.6f} {pct(Eh):+.3f}  {Ep:.6f} {pct(Ep):+.3f}  "
              f"{Em:.6f} {pct(Em):+.3f}   {1e3*(Eh-Ep):+.1f}")

    # ---------------- B. Z-scan --------------------------------------------
    print("\n=== B. Z-scan (E - E_exact, mHa; + = above exact = variational) ===")
    print(f"series check: e_series(1)-exact = {1e3*(e_series(1)-E_EXACT[1]):+.3f} mHa, "
          f"e_series(2)-exact = {1e3*(e_series(2)-E_EXACT[2]):+.4f} mHa")
    for n in avail:
        line = f"n_max={n}: "
        for Z in (1, 2, 3, 5, 10):
            Eh = ground(B[n].H(Z, True))[0]
            Ep = ground(B[n].H(Z, False))[0]
            line += (f"Z={Z}: hop {1e3*(Eh-E_EXACT[Z]):+8.1f} ({pct(Eh,Z):+6.2f}%) "
                     f"phys {1e3*(Ep-E_EXACT[Z]):+8.1f} ({pct(Ep,Z):+6.2f}%) | ")
        print(line)
    # H^- boundness: E vs E(H) = -0.5
    for n in avail:
        Eh = ground(B[n].H(1, True))[0]
        Ep = ground(B[n].H(1, False))[0]
        print(f"H^- n_max={n}: E_hop={Eh:.5f} (bound vs -0.5? {Eh < -0.5}; "
              f"err {pct(Eh,1):+.1f}%)   E_phys={Ep:.5f} (bound? {Ep < -0.5}; err {pct(Ep,1):+.1f}%)")
    # Z_c: E_hop(Z) = E_exact(Z) via the 1/Z series
    for n in avail:
        f = lambda Z, hop: ground(B[n].H(Z, hop))[0] - e_series(Z)
        out = f"Z_c n_max={n}: "
        for hop in (True, False):
            zs = np.linspace(1.0, 3.0, 41)
            vals = [f(z, hop) for z in zs]
            cross = [0.5 * (zs[i] + zs[i + 1]) for i in range(len(zs) - 1)
                     if np.sign(vals[i]) != np.sign(vals[i + 1])]
            if cross:
                lo, hi = [z for z in zs if z <= cross[0]][-1], [z for z in zs if z >= cross[0]][0]
                for _ in range(40):
                    mid = 0.5 * (lo + hi)
                    if np.sign(f(mid, hop)) == np.sign(f(lo, hop)):
                        lo = mid
                    else:
                        hi = mid
                zc = f"{0.5*(lo+hi):.4f}"
            else:
                zc = f"none in [1,3] (min E-E_ex = {1e3*min(vals):+.1f} mHa)"
            out += f"{'hop' if hop else 'phys'}: {zc}   "
        print(out)

    # ---------------- C. CUSP-2 intercept -----------------------------------
    print("\n=== C. CUSP-2 fit err_abs = m/(l_max+2)^4 + c over n_max>=4 ===")
    fitn = [n for n in avail if n >= 4]
    if len(fitn) >= 2:
        for Z in (2, 10):
            for hop in (True, False):
                x = [1.0 / (n + 1) ** 4 for n in fitn]
                y = [ground(B[n].H(Z, hop))[0] - E_EXACT[Z] for n in fitn]
                m, c = np.polyfit(x, y, 1)
                print(f"Z={Z} {'hop ' if hop else 'phys'}: c = {1e3*c:+.2f} mHa  "
                      f"(err at n_max={fitn[-1]}: {1e3*y[-1]:+.2f} mHa)")

    # ---------------- D. Paper 26 Table I ------------------------------------
    print("\n=== D. Paper 26 Table I (n_max=4): f_E^h1 f_E^Vee f_S^h1 f_S^Vee S ===")
    b = B[4]
    Zs = [2, 3, 4, 5, 6, 8, 10]
    for hop in (True, False):
        S_list = []
        for Z in Zs:
            H = b.H(Z, hop)
            E, ci = ground(H)
            Vfull = Z * b.V1
            Vd = np.diag(np.diag(Vfull))
            Hd = Z * Z * b.D1
            Hoff = b.HOP if hop else 0.0 * b.HOP
            Ecorr = E - H[0, 0]
            fh1 = float(ci @ Hoff @ ci - Hoff[0, 0]) / Ecorr
            fvee = float(ci @ (Vfull - Vd) @ ci) / Ecorr
            S = b.entropy(ci)
            SA = b.entropy(ground(Hd + Vd)[1])
            SB_ = b.entropy(ground(Hd + Hoff + Vd)[1])
            S_list.append(S)
            print(f"{'hop ' if hop else 'phys'} Z={Z:2d}: E={E:.6f} Ecorr={1e3*Ecorr:+.2f} mHa "
                  f"f_E^h1={fh1:.3f} f_E^Vee={fvee:.3f} f_S^h1={(SB_-SA)/S:.4f} "
                  f"f_S^Vee={(S-SB_)/S:.4f} S={S:.4f}")
        a, _, r2 = fit(Zs, S_list)
        print(f"   {'hop ' if hop else 'phys'}: S ~ Z^{a:.3f} (R2={r2:.4f}) over Z={Zs}")

    # ---------------- E. Paper 27 ------------------------------------------
    print("\n=== E1. Paper 27 Table ep1: n_max, E_full, S_full, S_kin ===")
    for n in (2, 3):
        for hop in (True, False):
            E, ci = ground(B[n].H(2, hop))
            Sk = B[n].entropy(ground(4 * B[n].D1 + (B[n].HOP if hop else 0 * B[n].HOP))[1])
            print(f"n_max={n} {'hop ' if hop else 'phys'}: E_full={E:.4f} S_full={B[n].entropy(ci):.5f} S_kin={Sk:.1e}")
    print("\n=== E2. Paper 27 EP-2c (n_max=3): S_B, w_B~, delta_B; fits ===")
    b = B[3]
    Zs = [2, 3, 4, 6, 8, 10]
    for hop in (True, False):
        S_l, w_l, d_l = [], [], []
        for Z in Zs:
            Hk = Z * Z * b.D1 + (b.HOP if hop else 0 * b.HOP)
            V = Z * b.V1
            E, ci = ground(Hk + V)
            eH, U = np.linalg.eigh(Hk)
            Vin = U.T @ V @ U
            voff = np.sqrt(max(0.0, np.sum(Vin ** 2) - np.sum(np.diag(Vin) ** 2)))
            nk = np.linalg.norm(Hk)
            gaps = np.diff(np.unique(np.round(eH, 10)))
            S_l.append(b.entropy(ci)); w_l.append(voff / nk); d_l.append(gaps[0] / nk)
            print(f"{'hop ' if hop else 'phys'} Z={Z:2d}: S_B={S_l[-1]:.3e} w~={w_l[-1]:.4f} delta={d_l[-1]:.4f}")
        g, A, r2 = fit(w_l, S_l)
        aZ, _, r2Z = fit(Zs, S_l)
        print(f"   {'hop ' if hop else 'phys'}: S=A w~^g: g={g:.3f} A={A:.3g} R2={r2:.4f};  S~Z^{aZ:.3f} (R2={r2Z:.4f})")

    # ---------------- F. "86%" candidates ------------------------------------
    print("\n=== F. candidate readings of Paper 13's '~86% of CI correlation' ===")
    for n in avail:
        Eh, Ep, Em, ph, pp = rowsA[n]
        E1s2 = B[n].H(2, True)[0, 0]
        hop_exp = float(ph @ B[n].HOP @ ph)
        print(f"n_max={n}: E(1s^2)={E1s2:.4f}; <hop>_psi_hop={1e3*hop_exp:+.1f} mHa; "
              f"(E_hop-E_phys)/(E_hop-E1s2)={(Eh-Ep)/(Eh-E1s2):.3f}; "
              f"<hop>/(E_hop-E1s2)={hop_exp/(Eh-E1s2):.3f}; "
              f"(E_hop-E_phys)/(E_ex-E_phys)={(Eh-Ep)/(E_EXACT[2]-Ep):.3f}; "
              f"(E_phys-E1s2)/(E_hop-E1s2)={(Ep-E1s2)/(Eh-E1s2):.3f}")

    # ---------------- G. Paper 27 gamma surface -------------------------------
    print("\n=== G. Paper 27 local slopes dlogS/dlogw~ (Z=15->30) and global Z in [2,30] fits ===")
    Zg = [2, 3, 4, 6, 8, 10, 12, 15, 20, 30]
    for n in avail:
        for hop in (True, False):
            Sw = []
            for Z in Zg + [50, 100]:
                Hk = Z * Z * B[n].D1 + (B[n].HOP if hop else 0 * B[n].HOP)
                V = Z * B[n].V1
                _, ci = ground(Hk + V)
                eH, U = np.linalg.eigh(Hk)
                Vin = U.T @ V @ U
                voff = np.sqrt(max(0.0, np.sum(Vin ** 2) - np.sum(np.diag(Vin) ** 2)))
                Sw.append((Z, B[n].entropy(ci), voff / np.linalg.norm(Hk)))
            d = {Z: (S, w) for Z, S, w in Sw}
            loc = lambda z1, z2: np.log(d[z2][0] / d[z1][0]) / np.log(d[z2][1] / d[z1][1])
            Zw = [2, 3, 4, 6, 10, 15, 30]   # Paper 27 tab:gamma_surface window
            gg, _, _ = fit([d[z][1] for z in Zw], [d[z][0] for z in Zw])
            print(f"n_max={n} {'hop ' if hop else 'phys'}: local(15->30)={loc(15,30):.4f}  "
                  f"local(2->3)={loc(2,3):.3f}  local(50->100)={loc(50,100):.4f}  global[2,30]={gg:.4f}")

    # ---------------- H. Paper 27 sec:cusp (edge attribution, Frobenius) ------
    print("\n=== H. Paper 27 sec:cusp: leave-one-edge-out entropy attribution ===")
    for n in (3, 4):
        b = B[n]
        for label, h1, V in (("published (h1 Z=2 + hop, V at k=1)", 4 * b.D1 + b.HOP, b.V1),
                             ("hop, V at k=Z=2", 4 * b.D1 + b.HOP, 2 * b.V1),
                             ("physical He (h1 diag, V at k=Z=2)", 4 * b.D1, 2 * b.V1)):
            _, c0 = ground(h1 + V)
            S0 = b.entropy(c0)
            ref = int(np.argmax(np.abs(c0)))
            rows = []
            for J in range(len(b.configs)):
                if J == ref or abs(V[ref, J]) < 1e-14:
                    continue
                Vc = V.copy()
                Vc[ref, J] = Vc[J, ref] = 0.0
                rows.append((S0 - b.entropy(ground(h1 + Vc)[1]), abs(V[ref, J]), J))
            rows.sort(reverse=True)
            off = V - np.diag(np.diag(V))
            hot = set(np.unravel_index(np.argmax(np.abs(off)), off.shape))
            topJ = rows[0][2]
            lab = lambda k: tuple(b.orbitals[x] for x in b.configs[k])
            # Frobenius diagonality of V in the H1 eigenbasis + commutator
            eH, U = np.linalg.eigh(h1)
            Vin = U.T @ V @ U
            dfr = np.sqrt(np.sum(np.diag(Vin) ** 2)) / np.linalg.norm(Vin)
            comm = np.linalg.norm(h1 @ V - V @ h1) / np.linalg.norm(h1)
            # control: largest V edge NOT incident on ref
            offc = np.abs(off).copy()
            offc[ref, :] = 0.0
            offc[:, ref] = 0.0
            ca, cb = np.unravel_index(np.argmax(offc), offc.shape)
            Vc = V.copy()
            Vc[ca, cb] = Vc[cb, ca] = 0.0
            dS_ctrl = S0 - b.entropy(ground(h1 + Vc)[1])
            print(f"   control edge {lab(ca)}<->{lab(cb)} |V|={abs(V[ca,cb]):.4f} "
                  f"(top |V|={rows[0][1]:.4f}): dS={dS_ctrl:.1e}")
            from scipy.stats import spearmanr
            rho_s = spearmanr([r[0] for r in rows], [r[1] for r in rows])[0]
            print(f"n_max={n} [{label}]: S_full={S0:.3e}; top edge {lab(ref)}<->{lab(topJ)} "
                  f"carries {rows[0][0]/S0:.1%}; is hottest V edge: {hot == {ref, topJ}}; "
                  f"runner-up {lab(rows[1][2])} {rows[1][0]/S0:.1%}; Spearman={rho_s:+.2f}; "
                  f"diag-Frobenius={dfr:.3f}; ||[H1,V]||/||H1||={comm:.3f}")


if __name__ == "__main__":
    main()
