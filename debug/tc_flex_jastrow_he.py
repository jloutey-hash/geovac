"""Stage-1 PoC: flexible, cusp-preserving, VARIANCE-optimized Jastrow for the native
non-Hermitian TC on He (Coulomb-Sturmian, s-only).  FAST polynomial-in-c formulation.

Hypothesis (v5.0.7 fork + Track C): the single-parameter geminal's fragility (variance-min
over-transcorrelates to the small-gamma edge) is a RIGIDITY artifact.  A FLEXIBLE multi-term
Jastrow, parameters chosen by minimizing the TC reference-energy VARIANCE, should find a
genuine INTERIOR minimum that is ALSO accurate -- and cheap (kappa_V O(1), no S matrix).

Flexible Jastrow (cusp-preserving):  u(r) = sum_mu c_mu * [ -(1/2 g_mu) e^{-g_mu r} ],
  u'(0)=1/2  <=>  sum_mu c_mu = 1  (singlet Kato cusp).  Free: {c_mu} on a fixed grid {g_mu}.

KEY STRUCTURE -- the TC 2-body operator is a degree-2 POLYNOMIAL in c:
  u'(r)   = (1/2) sum_mu c_mu e^{-g_mu r}                              (K: linear in c)
  w(r;c)  = 1/r  +  sum_mu c_mu wlin_mu(r)  +  sum_{mu,nu} c_mu c_nu wquad_{mu,nu}(r)
    wlin_mu(r)      = (g_mu/2) e^{-g_mu r} - e^{-g_mu r}/r    ( = -Laplacian(u_mu) )
    wquad_{mu,nu}(r)= -(1/4) e^{-(g_mu+g_nu) r}               ( = -(u_mu' u_nu') )
So the ERI tensors (hence asym) are precomputed ONCE per exponent / exponent-pair; each
objective eval is a cheap tensor combination.  (M=1,c=1 reduces to the engine's w_kernel -- validated.)

Objective: sigma^2(c) = sum_{D != Phi} |<D|Htilde|Phi>|^2  (reference-energy variance),
  Phi = aufbau reference det; = squared off-ref entries of the H-tilde reference column.
Gate: does the sigma^2 optimum give E_TC within ~1-2 mHa of exact (-2.903724), kappa_V O(1)?
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
GAMMAS = np.array([0.3, 0.5, 0.8, 1.2, 1.8, 2.6])
NG, NX = 500, 128


def _eri_from_kernel(P, Kmat):
    """<ij|Kmat|kl> physicist tensor (multiplicative kernel), then Lowdin+asym."""
    ns = P["ns"]; D = P["D"]
    out = np.zeros((ns, ns, ns, ns))
    for i in range(ns):
        for j in range(ns):
            for kk in range(ns):
                for l in range(ns):
                    out[i, j, kk, l] = D[(i, kk)] @ Kmat @ D[(j, l)]
    return TC.asym_from_phys(TC.transform_2(out, P["Xm"]), P["nso"])


def _eri_conv(P, kA, kB):
    """convective K tensor from vertex kernels kA,kB, then Lowdin+asym."""
    ns = P["ns"]; D = P["D"]; Dd = P["Dd"]
    out = np.zeros((ns, ns, ns, ns))
    for i in range(ns):
        for j in range(ns):
            for kk in range(ns):
                for l in range(ns):
                    out[i, j, kk, l] = -(Dd[(i, kk)] @ kA @ D[(j, l)]) + (D[(i, kk)] @ kB @ Dd[(j, l)])
    return TC.asym_from_phys(TC.transform_2(out, P["Xm"]), P["nso"])


def setup():
    t0 = time.time()
    r, wr = TC.make_grid(K, Ng=NG)
    S, h1s, Rtab, W = TC.build_one_body(NS, r, wr, K, Z)
    Xm = TC.lowdin(S); h1o = TC.transform_1(h1s, Xm)
    nso = 2 * NS
    dets, didx = TC.make_dets(nso, NE)
    hso = TC.h_spin(h1o, nso)
    de = np.diag(h1o); order = np.argsort(de)
    ref = tuple(sorted((2 * order[0], 2 * order[0] + 1)))
    RR = {i: Rtab[i + 1][0] for i in range(NS)}
    dRR = {i: Rtab[i + 1][1] for i in range(NS)}
    Dmat = {(i, k): RR[i] * RR[k] * W for i in range(NS) for k in range(NS)}
    Dd = {(i, k): RR[i] * dRR[k] * W for i in range(NS) for k in range(NS)}
    P = dict(r=r, wr=wr, Rtab=Rtab, W=W, Xm=Xm, hso=hso, nso=nso, ns=NS,
             dets=dets, didx=didx, ref_idx=didx[ref], D=Dmat, Dd=Dd)

    # --- one angular pass: accumulate all per-exponent / per-pair kernels ---
    M = len(GAMMAS)
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    xs, ws = np.polynomial.legendre.leggauss(NX)
    Kcoul = np.zeros_like(R1)
    Kwlin = [np.zeros_like(R1) for _ in range(M)]
    KkA = [np.zeros_like(R1) for _ in range(M)]
    KkB = [np.zeros_like(R1) for _ in range(M)]
    pairs = list(itertools.combinations_with_replacement(range(M), 2))
    Kwquad = {p: np.zeros_like(R1) for p in pairs}
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

    # --- precompute asym tensors (Lowdin basis) ---
    P["asym_coul"] = _eri_from_kernel(P, Kcoul)
    P["asym_wlin"] = [_eri_from_kernel(P, Kwlin[mu]) for mu in range(M)]
    P["asym_wquad"] = {p: _eri_from_kernel(P, Kwquad[p]) for p in pairs}
    P["asym_Kconv"] = [_eri_conv(P, KkA[mu], KkB[mu]) for mu in range(M)]
    P["pairs"] = pairs
    P["setup_s"] = time.time() - t0
    return P


def asym_TC(P, cs):
    """assemble the full non-Herm TC 2-body asym tensor for coefficients cs (poly in c)."""
    A = P["asym_coul"].copy()
    for mu in range(len(cs)):
        A = A + cs[mu] * P["asym_wlin"][mu] + cs[mu] * P["asym_Kconv"][mu]
    for (mu, nu) in P["pairs"]:
        w = cs[mu] * cs[nu] * (1.0 if mu == nu else 2.0)   # symmetric double-count
        A = A + w * P["asym_wquad"][(mu, nu)]
    return A


def Hmat(P, cs):
    return TC.build_H(P["dets"], P["didx"], P["hso"], asym_TC(P, cs), P["nso"])


def sigma2(P, cs):
    H = Hmat(P, cs)
    col = H[:, P["ref_idx"]].copy(); col[P["ref_idx"]] = 0.0
    return float(np.sum(np.abs(col) ** 2))


def full_from_cs(P, cs):
    H = Hmat(P, cs)
    E, im = TC.ground(H, hermitian=False)
    mm = TC.measure_matrix(H)
    col = H[:, P["ref_idx"]].copy(); col[P["ref_idx"]] = 0.0
    return dict(E=E, im=im, dTC_mHa=(E - EXACT) * 1e3, kV=mm["kV_low"],
                nonnorm=mm["nonnormality"], sigma2=float(np.sum(np.abs(col) ** 2)))


def main():
    P = setup()
    M = len(GAMMAS)
    print("setup %.1fs   grid Ng=%d nx=%d   exponents %s" % (P["setup_s"], NG, NX, GAMMAS))

    # single-geminal calibration (validation: M=1 must reproduce the engine)
    print("\n=== single-geminal (M=1) calibration ===")
    print("%6s %12s %10s %10s %8s" % ("gamma", "E_TC", "dTC_mHa", "sigma2", "kV"))
    sg = []
    for mu, g in enumerate(GAMMAS):
        cs = np.zeros(M); cs[mu] = 1.0
        d = full_from_cs(P, cs); sg.append((g, d))
        print("%6.2f %12.6f %+10.3f %10.4f %8.3f" % (g, d["E"], d["dTC_mHa"], d["sigma2"], d["kV"]))
    sg_bestvar = min(sg, key=lambda t: t[1]["sigma2"])
    print("single-geminal min-sigma2 pick: gamma=%.2f -> dTC=%+.3f mHa (Track A's edge)"
          % (sg_bestvar[0], sg_bestvar[1]["dTC_mHa"]))

    # variance minimization over flexible {c}, multi-start
    def free_to_c(free):
        return np.append(free, 1.0 - np.sum(free))
    def obj(free):
        return sigma2(P, free_to_c(free))
    starts = []
    for mu in range(M):                      # start near each single geminal
        c0 = np.zeros(M); c0[mu] = 1.0
        starts.append(c0[:-1])
    best = None
    for s0 in starts:
        res = minimize(obj, s0, method="Nelder-Mead",
                       options=dict(xatol=1e-7, fatol=1e-14, maxiter=8000, maxfev=12000))
        if best is None or res.fun < best.fun:
            best = res
    cs = free_to_c(best.x)
    vo = full_from_cs(P, cs)
    print("\n=== variance-minimized flexible Jastrow (M=%d) ===" % M)
    print("sigma^2 = %.6e  (single-geminal min was %.4f)" % (vo["sigma2"], sg_bestvar[1]["sigma2"]))
    print("c = %s  (sum=%.4f)" % (np.array2string(cs, precision=4), np.sum(cs)))
    print("E_TC = %.6f   dTC = %+.3f mHa   kappa_V = %.3f   nonnorm = %.4f   im0=%.1e"
          % (vo["E"], vo["dTC_mHa"], vo["kV"], vo["nonnorm"], vo["im"]))

    # oracle: min |E - exact| over the same flexible space (how good CAN this basis+Jastrow be)
    def oacc(free):
        return abs(full_from_cs(P, free_to_c(free))["E"] - EXACT)
    ores = minimize(oacc, best.x, method="Nelder-Mead",
                    options=dict(xatol=1e-7, fatol=1e-10, maxiter=8000, maxfev=12000))
    oc = free_to_c(ores.x); od = full_from_cs(P, oc)
    print("\noracle (min|E-exact| over flexible space): dTC=%+.3f mHa  kV=%.3f" % (od["dTC_mHa"], od["kV"]))

    out = dict(system="He", ns=NS, k=K, gammas=GAMMAS.tolist(), exact=EXACT,
               single_geminal=[dict(gamma=g, **{k2: d[k2] for k2 in ("E", "dTC_mHa", "sigma2", "kV")}) for g, d in sg],
               variance_opt=dict(c=cs.tolist(), **vo),
               oracle=dict(c=oc.tolist(), **od),
               refs=dict(R12CI_He_best_mHa=0.80, single_geminal_oracle_mHa=1.2))
    os.makedirs(os.path.join(_HERE, "data"), exist_ok=True)
    json.dump(out, open(os.path.join(_HERE, "data", "tc_flex_jastrow_he.json"), "w"), indent=2, default=float)
    print("\nwrote debug/data/tc_flex_jastrow_he.json")

    print("\n--- GATE ---")
    acc = abs(vo["dTC_mHa"])
    print("variance-opt accuracy = %.2f mHa ; kappa_V = %.2f" % (acc, vo["kV"]))
    if acc <= 2.0 and vo["kV"] <= 3.0:
        print("GO: flexible variance-optimized Jastrow lands ~chemical accuracy AND stays cheap.")
    elif acc <= 5.0 and vo["kV"] <= 3.0:
        print("BORDERLINE: beats the single-geminal edge but not <2 mHa; try more/other exponents.")
    else:
        print("STOP: variance optimum still not accurate (%.1f mHa) -> tension deeper than flexibility." % acc)


if __name__ == "__main__":
    main()
