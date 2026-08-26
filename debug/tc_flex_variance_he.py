"""Stage-2 (He core): FULL deterministic variance-optimized flexible Jastrow.

The well-posed objective (zero-variance principle, from the Stage-1 correction):
  Var[E_L] = <Psi_T|(H-E)^2|Psi_T> / <Psi_T|Psi_T>,   Psi_T = e^{J} Phi,  J = sum_{i<j} u(r_ij).
Using (H-E) e^{J} = e^{J}(Htilde - E), everything reduces to finite-basis objects:
  W2      = FCI matrix of the e^{2u(r12)} two-body operator (the VMC weight) -- rebuilt per c
  E_w     = (W2 @ Htilde)[ref,ref] / W2[ref,ref]                         (VMC correlated energy)
  w       = Htilde[:,ref] - E_w * e_ref
  Var(c)  = (w^T W2 w) / W2[ref,ref]                                     (>= 0, zero-variance principle)
The e^{2u} WEIGHT is exactly what the ill-posed Stage-1 naive sigma^2 dropped.

Jastrow: u(r) = sum_mu c_mu * [-(1/2 g_mu) e^{-g_mu r}], cusp <=> sum c_mu = 1, PHYSICALLY BOUNDED
(|c_mu| <= CMAX) to stay in the BCH-valid regime (Stage-1: unbounded escapes to spurious large-J).

Htilde is assembled by the fast polynomial engine (H is quadratic in c); only W2 is rebuilt per c.
After optimizing Var, TRANSPLANT c* to the cheap TC operator and report E_ground(Htilde), kappa_V.
Gate: does the variance-optimal PHYSICAL Jastrow give E_ground within ~1-2 mHa of exact, kappa_V O(1)?
"""
import os, sys, json, itertools, time
import numpy as np
from scipy.optimize import minimize

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
for p in (_ROOT, _HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
import geovac.transcorrelated_sturmian as TC

EXACT = -2.903724
NS, K, Z, NE = 3, 1.9, 2, 2
GAMMAS = np.array([0.4, 0.7, 1.1, 1.7, 2.5])
CMAX = 4.0
NG, NX = 400, 96


# ---------- Htilde: precomputed polynomial-in-c asym pieces ----------
def _asym_mult(P, Kmat):
    ns = P["ns"]; D = P["D"]
    out = np.zeros((ns, ns, ns, ns))
    for i in range(ns):
        for j in range(ns):
            for kk in range(ns):
                for l in range(ns):
                    out[i, j, kk, l] = D[(i, kk)] @ Kmat @ D[(j, l)]
    return TC.asym_from_phys(TC.transform_2(out, P["Xm"]), P["nso"])


def _asym_conv(P, kA, kB):
    ns = P["ns"]; D = P["D"]; Dd = P["Dd"]
    out = np.zeros((ns, ns, ns, ns))
    for i in range(ns):
        for j in range(ns):
            for kk in range(ns):
                for l in range(ns):
                    out[i, j, kk, l] = -(Dd[(i, kk)] @ kA @ D[(j, l)]) + (D[(i, kk)] @ kB @ Dd[(j, l)])
    return TC.asym_from_phys(TC.transform_2(out, P["Xm"]), P["nso"])


def rhf(h1o, eri_L, maxit=300, tol=1e-11):
    """Closed-shell RHF (2e) in the Lowdin-orthonormal basis (S=I).
    eri_L[i,j,k,l] = (ik|jl) Mulliken.  Returns (U=Fock eigvecs, Ehf)."""
    eps, U = np.linalg.eigh(h1o)
    C = U[:, 0]; Eprev = 0.0
    for _ in range(maxit):
        P = 2.0 * np.outer(C, C)
        J = np.einsum('prqs,rs->pq', eri_L, P)     # (pq|rs) = eri_L[p,r,q,s]
        Kx = np.einsum('pqrs,rs->pq', eri_L, P)     # (pr|qs) = eri_L[p,q,r,s]
        F = h1o + J - 0.5 * Kx
        eps, U = np.linalg.eigh(F)
        C = U[:, 0]
        Ehf = 0.5 * float(np.sum(P * (h1o + F)))
        if abs(Ehf - Eprev) < tol:
            break
        Eprev = Ehf
    return U, Ehf


def setup():
    t0 = time.time()
    r, wr = TC.make_grid(K, Ng=NG)
    S, h1s, Rtab, W = TC.build_one_body(NS, r, wr, K, Z)
    Xm0 = TC.lowdin(S)
    # RHF in the Lowdin basis, then rotate to HF orbitals: X = Xm0 @ U_hf
    h1o_L = TC.transform_1(h1s, Xm0)
    Km0 = TC.build_kernels(r, 1.0, nx=NX)
    eri_c0, _, _ = TC.two_body(NS, Rtab, W, Km0)
    eri_L = TC.transform_2(eri_c0, Xm0)
    U_hf, Ehf = rhf(h1o_L, eri_L)
    # stay in the Lowdin basis (transform_1/2 assume the symmetric Lowdin matrix); carry the
    # HF reference as a CI VECTOR |Phi_HF> = |phi_HF^up phi_HF^dn>, phi_HF = sum_p c_p |p>_Lowdin.
    Xm = Xm0
    h1o = TC.transform_1(h1s, Xm)
    nso = 2 * NS
    dets, didx = TC.make_dets(nso, NE)
    hso = TC.h_spin(h1o, nso)
    c_hf = U_hf[:, 0]
    ref_vec = np.zeros(len(dets))
    for di, det in enumerate(dets):
        a, b = det                                   # spin-orbitals, a<b; 2p=up_p, 2p+1=dn_p
        if a % 2 == 0 and b % 2 == 1:                # up_a  dn_b
            ref_vec[di] = c_hf[a // 2] * c_hf[(b - 1) // 2]
        elif a % 2 == 1 and b % 2 == 0:              # dn_a  up_b  -> sign flip
            ref_vec[di] = -c_hf[b // 2] * c_hf[(a - 1) // 2]
    ref_vec /= np.linalg.norm(ref_vec)
    print("[RHF] E_HF = %.6f  (exact -2.903724)" % Ehf)
    RR = {i: Rtab[i + 1][0] for i in range(NS)}
    dRR = {i: Rtab[i + 1][1] for i in range(NS)}
    Dmat = {(i, k): RR[i] * RR[k] * W for i in range(NS) for k in range(NS)}
    Dd = {(i, k): RR[i] * dRR[k] * W for i in range(NS) for k in range(NS)}
    P = dict(r=r, wr=wr, Rtab=Rtab, W=W, Xm=Xm, hso=hso, nso=nso, ns=NS,
             dets=dets, didx=didx, ref_vec=ref_vec, D=Dmat, Dd=Dd)

    M = len(GAMMAS)
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    xs, ws = np.polynomial.legendre.leggauss(NX)
    Kcoul = np.zeros_like(R1)
    Kwlin = [np.zeros_like(R1) for _ in range(M)]
    KkA = [np.zeros_like(R1) for _ in range(M)]
    KkB = [np.zeros_like(R1) for _ in range(M)]
    pairs = list(itertools.combinations_with_replacement(range(M), 2))
    Kwquad = {p: np.zeros_like(R1) for p in pairs}
    # cache exp tables for the per-c W2 (e^{2u}) kernel
    P["_R1R2x"] = []
    for x, wx in zip(xs, ws):
        r12 = np.sqrt(np.maximum(R1 * R1 + R2 * R2 - 2 * R1 * R2 * x, 1e-30))
        half = 0.5 * wx
        Kcoul += half / r12
        e = [np.exp(-g * r12) for g in GAMMAS]
        for mu, g in enumerate(GAMMAS):
            Kwlin[mu] += half * ((g / 2.0) * e[mu] - e[mu] / r12)
            up = 0.5 * e[mu]
            KkA[mu] += half * up * (R1 - R2 * x) / r12
            KkB[mu] += half * up * (R1 * x - R2) / r12
        for (mu, nu) in pairs:
            Kwquad[(mu, nu)] += half * (-0.25) * np.exp(-(GAMMAS[mu] + GAMMAS[nu]) * r12)
        P["_R1R2x"].append((half, e))          # for e^{2u} rebuild: u = sum c_mu(-1/2g)e
    P["asym_coul"] = _asym_mult(P, Kcoul)
    P["asym_wlin"] = [_asym_mult(P, Kwlin[mu]) for mu in range(M)]
    P["asym_wquad"] = {p: _asym_mult(P, Kwquad[p]) for p in pairs}
    P["asym_Kconv"] = [_asym_conv(P, KkA[mu], KkB[mu]) for mu in range(M)]
    P["pairs"] = pairs
    P["amp"] = np.array([-1.0 / (2.0 * g) for g in GAMMAS])   # u_mu amplitude -(1/2g)
    P["setup_s"] = time.time() - t0
    return P


def asym_TC(P, cs):
    A = P["asym_coul"].copy()
    for mu in range(len(cs)):
        A = A + cs[mu] * P["asym_wlin"][mu] + cs[mu] * P["asym_Kconv"][mu]
    for (mu, nu) in P["pairs"]:
        A = A + cs[mu] * cs[nu] * (1.0 if mu == nu else 2.0) * P["asym_wquad"][(mu, nu)]
    return A


def Htilde(P, cs):
    return TC.build_H(P["dets"], P["didx"], P["hso"], asym_TC(P, cs), P["nso"])


def build_W2(P, cs):
    """FCI matrix of the multiplicative operator e^{2 u(r12)} (VMC weight)."""
    amp_c = P["amp"] * cs                    # coefficients of e^{-g r} in u(r)
    Kw2 = np.zeros((len(P["r"]), len(P["r"])))
    for (half, e) in P["_R1R2x"]:
        u = np.zeros_like(e[0])
        for a, em in zip(amp_c, e):
            u += a * em
        Kw2 += half * np.exp(2.0 * u)
    asymW2 = _asym_mult(P, Kw2)
    zero_h = np.zeros_like(P["hso"])
    return TC.build_H(P["dets"], P["didx"], zero_h, asymW2, P["nso"])


def variance(P, cs):
    """Zero-variance objective with a CI-vector reference Phi (= Phi_HF).
    Var = <(H-E)Psi_T|(H-E)Psi_T>/<Psi_T|Psi_T>, Psi_T=e^{J}Phi, (H-E)e^J=e^J(Htilde-E).
    """
    H = Htilde(P, cs)
    W2 = build_W2(P, cs)
    phi = P["ref_vec"]
    Wphi = W2 @ phi
    nrm = float(phi @ Wphi)
    hphi = H @ phi                      # Htilde |Phi>
    Ew = float((phi @ (W2 @ hphi)) / nrm)
    w = hphi - Ew * phi
    var = float((w @ (W2 @ w)) / nrm)
    return var, Ew


def evaluate(P, cs):
    H = Htilde(P, cs)
    E, im = TC.ground(H, hermitian=False)
    mm = TC.measure_matrix(H)
    var, Ew = variance(P, cs)
    return dict(E_ground=E, im=im, dTC_mHa=(E - EXACT) * 1e3, kV=mm["kV_low"],
                nonnorm=mm["nonnormality"], variance=var, E_weighted=Ew,
                dEw_mHa=(Ew - EXACT) * 1e3)


def main():
    P = setup()
    M = len(GAMMAS)
    print("setup %.1fs  Ng=%d nx=%d  exponents=%s  CMAX=%.1f" % (P["setup_s"], NG, NX, GAMMAS, CMAX))

    print("\n=== single-geminal: variance objective sanity (Var>=0, tracks accuracy) ===")
    print("%6s %12s %10s %12s %10s %8s" % ("gamma", "E_ground", "dGrnd", "Var", "dEw_mHa", "kV"))
    for mu, g in enumerate(GAMMAS):
        cs = np.zeros(M); cs[mu] = 1.0
        d = evaluate(P, cs)
        print("%6.2f %12.6f %+10.3f %12.6e %+10.3f %8.3f" %
              (g, d["E_ground"], d["dTC_mHa"], d["variance"], d["dEw_mHa"], d["kV"]))

    # ---- constrained variance minimization (physical bounded Jastrow) ----
    def free_to_c(free):
        return np.append(free, 1.0 - np.sum(free))
    def penalty(cs):
        return 1e3 * np.sum(np.maximum(0.0, np.abs(cs) - CMAX) ** 2)
    def obj(free):
        cs = free_to_c(free)
        v, _ = variance(P, cs)
        return v + penalty(cs)
    best = None
    for mu in range(M):
        c0 = np.zeros(M); c0[mu] = 1.0
        res = minimize(obj, c0[:-1], method="Nelder-Mead",
                       options=dict(xatol=1e-6, fatol=1e-12, maxiter=6000, maxfev=9000))
        if best is None or res.fun < best.fun:
            best = res
    cs = free_to_c(best.x)
    d = evaluate(P, cs)
    print("\n=== variance-minimized PHYSICAL flexible Jastrow (M=%d, |c|<=%.1f) ===" % (M, CMAX))
    print("c = %s  (sum=%.4f, max|c|=%.2f)" % (np.array2string(cs, precision=3), np.sum(cs), np.max(np.abs(cs))))
    print("Var = %.6e   E_weighted(VMC) dEw = %+.3f mHa" % (d["variance"], d["dEw_mHa"]))
    print("TRANSPLANT -> cheap TC operator:  E_ground dTC = %+.3f mHa   kappa_V = %.3f   nonnorm=%.4f"
          % (d["dTC_mHa"], d["kV"], d["nonnorm"]))

    out = dict(system="He", ns=NS, k=K, gammas=GAMMAS.tolist(), CMAX=CMAX, exact=EXACT,
               variance_opt=dict(c=cs.tolist(), **d),
               refs=dict(R12CI_He_best_mHa=0.80, single_geminal_oracle_mHa=1.2))
    os.makedirs(os.path.join(_HERE, "data"), exist_ok=True)
    json.dump(out, open(os.path.join(_HERE, "data", "tc_flex_variance_he.json"), "w"), indent=2, default=float)
    print("\nwrote debug/data/tc_flex_variance_he.json")

    print("\n--- GATE ---")
    accg = abs(d["dTC_mHa"]); accw = abs(d["dEw_mHa"])
    print("transplanted TC-operator accuracy = %.2f mHa ; VMC energy = %.2f mHa ; kappa_V = %.2f"
          % (accg, accw, d["kV"]))
    if accg <= 2.0 and d["kV"] <= 3.0 and np.max(np.abs(cs)) <= CMAX + 1e-6:
        print("GO: physical variance-optimized Jastrow -> cheap TC operator is accurate. Capstone works.")
    elif accg <= 5.0:
        print("BORDERLINE: within a few mHa; try richer exponent grid / relax CMAX slightly.")
    else:
        print("STOP: even the well-posed variance optimum (physical Jastrow) misses -> report the wall.")


if __name__ == "__main__":
    main()
