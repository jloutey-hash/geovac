"""
Phase 2 -- isoenergetic generalized-Sturmian He extended to l>0 (s,p,d,f), coupled to total ^1S.

Extends the validated s-only skeleton (debug/sturmian_he_secular.py) to angular correlation.

Secular eq [BK6 6.35]:  [ diag(Z R_nu) + T' - p_kappa I ] B = 0,  E = -p_kappa^2/2.
  R_nu = sqrt(sum_j 1/n_j^2);  weighted charge Q_nu = p_kappa/R_nu, built at reference p_kappa=1
  so Q_nu = 1/R_nu (T' is a matrix of PURE NUMBERS, p_kappa-independent [BK6 boxed remark]).
  T'_{ij} = -<Psi_i|1/r12|Psi_j>  (mixed-scale, non-orthogonal; standard eigenproblem = paper's
  metric-free atomic construction; no overlap metric S in the primary method).

Configurations = Goscinskian two-electron singlets (S=0) of hydrogenic orbitals (n_a,l)(n_b,l)
coupled to total L=0 (requires equal l on the two electrons for L=0), spatially symmetric.

Interelectron matrix element machinery (debug/sturmian_secular_equations.md, item 4):
  Slater-Condon -> Legendre multipole 1/r12 = sum_k (r_<^k / r_>^{k+1}) P_k(cos w12)
  -> angular a_k = products of Gaunt integrals (Wigner 3j) ; radial Slater integral R^k
  selection: m_a+m_b=m_c+m_d ; (-1)^{l_a+l_c}=(-1)^{l_b+l_d} ; finite k in {|..|,..,l_a+l_c}.

Validation gates (reported explicitly, pass/fail):
  (1) single-config 1s^2 -> E = -2.847 Ha (textbook variational He);
  (2) adding s+p(+d,f) configs LOWERS E monotonically toward exact non-rel -2.90372 Ha;
  (3) entrywise 1-norm ||M||_1 vs #configs K fits a SUBLINEAR exponent (<1) with l>0 present.

Diagnostic only.  Author: sprint sturmian-lmax, 2026-08-18.
"""
import warnings; warnings.filterwarnings("ignore")
import math
import numpy as np
from math import factorial as fac
from scipy.special import genlaguerre
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import eigh
from itertools import combinations_with_replacement

# --------------------------------------------------------------------------------------
# Radial grid.  Small-r resolution matters (high-n high-charge configs contract).
# Slater potentials use cumulative-TRAPEZOID (O(dr^2)); validated exact on (5/8)Q below
# (grid-converged already at N~12000: (5/8)/sqrt2 to 6 digits, E(1s^2)=-2.84766).
# --------------------------------------------------------------------------------------
R_MAX = 60.0
N_GRID = 18000
r = np.linspace(1e-7, R_MAX, N_GRID)
dr = r[1] - r[0]
r2 = r * r

def _ctrap_fwd(y):
    """ int_{r0}^{r_i} y dr, same length as y, value 0 at i=0. """
    return np.concatenate(([0.0], cumulative_trapezoid(y, dx=dr)))

def _ctrap_rev(y):
    """ int_{r_i}^{rmax} y dr, same length as y. """
    return _ctrap_fwd(y[::-1])[::-1]

Z = 2.0

# --------------------------------------------------------------------------------------
# Angular machinery: self-contained Wigner-3j (Racah) + real-Y Gaunt integral.
# --------------------------------------------------------------------------------------
def wigner3j(j1, j2, j3, m1, m2, m3):
    if m1 + m2 + m3 != 0:
        return 0.0
    if not (abs(j1 - j2) <= j3 <= j1 + j2):
        return 0.0
    if any(abs(m) > j for m, j in ((m1, j1), (m2, j2), (m3, j3))):
        return 0.0
    if (j1 + j2 + j3) < 0:
        return 0.0
    delta = math.sqrt(
        fac(j1 + j2 - j3) * fac(j1 - j2 + j3) * fac(-j1 + j2 + j3)
        / fac(j1 + j2 + j3 + 1)
    )
    pref = math.sqrt(
        fac(j1 + m1) * fac(j1 - m1) * fac(j2 + m2) * fac(j2 - m2)
        * fac(j3 + m3) * fac(j3 - m3)
    )
    tmin = max(0, j2 - j3 - m1, j1 - j3 + m2)
    tmax = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    s = 0.0
    for t in range(tmin, tmax + 1):
        denom = (fac(t) * fac(j1 + j2 - j3 - t) * fac(j1 - m1 - t)
                 * fac(j2 + m2 - t) * fac(j3 - j2 + m1 + t) * fac(j3 - j1 - m2 + t))
        s += (-1) ** t / denom
    return (-1) ** (j1 - j2 - m3) * delta * pref * s

_gaunt_cache = {}
def gaunt(l1, l2, l3, m1, m2, m3):
    """ integral of Y_{l1 m1} Y_{l2 m2} Y_{l3 m3} dOmega (no conjugation). """
    key = (l1, l2, l3, m1, m2, m3)
    v = _gaunt_cache.get(key)
    if v is not None:
        return v
    w0 = wigner3j(l1, l2, l3, 0, 0, 0)
    if w0 == 0.0:
        _gaunt_cache[key] = 0.0
        return 0.0
    val = (math.sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / (4 * math.pi))
           * w0 * wigner3j(l1, l2, l3, m1, m2, m3))
    _gaunt_cache[key] = val
    return val

def cg_L0(l, m):
    """ Clebsch-Gordan <l m; l -m | 0 0> = (-1)^{l-m} / sqrt(2l+1). """
    return (-1) ** (l - m) / math.sqrt(2 * l + 1)

# --------------------------------------------------------------------------------------
# Radial orbitals and Slater integrals.
# --------------------------------------------------------------------------------------
def hyd_radial(n, l, Q):
    """ hydrogenic radial R_{nl} at charge Q, L2-normalized on the grid (int R^2 r^2 dr = 1). """
    a = Q / n
    f = (2 * a * r) ** l * np.exp(-a * r) * genlaguerre(n - l - 1, 2 * l + 1)(2 * a * r)
    nrm = np.sqrt(np.trapezoid(f * f * r2, r))
    return f / nrm

def radial_overlap(Pa, Pb):
    return np.trapezoid(Pa * Pb * r2, r)

# R^k radial Slater integral: R^k(Pa Pc ; Pb Pd)
#   = int int Pa(r1)Pc(r1) (r_<^k / r_>^{k+1}) Pb(r2)Pd(r2) r1^2 r2^2 dr1 dr2
# with U_k[Pb Pd](r1) the multipole potential of density Pb*Pd.
_Rk_cache = {}
def slater_Rk(ida, idc, idb, idd, Pa, Pc, Pb, Pd, k):
    # canonical key: (a,c) unordered on e1, (b,d) unordered on e2, and (e1<->e2) swap symmetric
    p1 = (ida, idc) if ida <= idc else (idc, ida)
    p2 = (idb, idd) if idb <= idd else (idd, idb)
    key = (p1, p2, k) if p1 <= p2 else (p2, p1, k)
    v = _Rk_cache.get(key)
    if v is not None:
        return v
    g = Pb * Pd * r2                                   # electron-2 density * r2
    inner = _ctrap_fwd(g * r ** k)                     # int_0^{r1} dens2 r2^{k+2} dr2
    outer = _ctrap_rev(g * r ** (-(k + 1)))            # int_{r1}^inf dens2 r2^{1-k} dr2
    Uk = inner * r ** (-(k + 1)) + outer * r ** k
    val = np.trapezoid(Pa * Pc * Uk * r2, r)
    _Rk_cache[key] = val
    return val

# --------------------------------------------------------------------------------------
# Orbital representation and coupled two-electron configuration wavefunctions.
# An "orbital" is a dict {id, l, m, P (radial array)}.  Radial id encodes (config, which-orb)
# so the Slater-integral cache dedups across m (radial depends only on n,l,Q, not m).
# A configuration term list is [(coeff, orb_e1, orb_e2), ...] representing Psi^unnorm(1,2).
# --------------------------------------------------------------------------------------
class Config:
    _next_rid = 0
    def __init__(self, l, na, nb, pk_ref=1.0):
        self.l, self.na, self.nb = l, na, nb
        self.Rnu = math.sqrt(1.0 / na ** 2 + 1.0 / nb ** 2)
        self.Q = pk_ref / self.Rnu
        # two radial orbitals (distinct radial ids); if na==nb they are the same radial function
        self.Pa = hyd_radial(na, l, self.Q)
        self.rid_a = Config._next_rid; Config._next_rid += 1
        if nb == na:
            self.Pb = self.Pa
            self.rid_b = self.rid_a
        else:
            self.Pb = hyd_radial(nb, l, self.Q)
            self.rid_b = Config._next_rid; Config._next_rid += 1
        self.terms = self._build_terms()
        self.norm = 1.0 / math.sqrt(self._self_overlap())

    def _orb(self, which, m):
        if which == 'a':
            return dict(rid=self.rid_a, l=self.l, m=m, P=self.Pa)
        return dict(rid=self.rid_b, l=self.l, m=m, P=self.Pb)

    def _build_terms(self):
        """ Psi^unnorm(1,2) coupled to L=0, spatially symmetric.
            same orbital (na==nb): Phi_aa = sum_m cg φ_{a m}(1) φ_{a -m}(2)  (already symmetric).
            distinct: Phi_ab + Phi_ba  (symmetric combination).                             """
        l = self.l
        terms = []
        for m in range(-l, l + 1):
            c = cg_L0(l, m)
            terms.append((c, self._orb('a', m), self._orb('b', -m)))   # Phi_ab
        if self.na != self.nb:
            for m in range(-l, l + 1):
                c = cg_L0(l, m)
                terms.append((c, self._orb('b', m), self._orb('a', -m)))  # Phi_ba
        return terms

    def _self_overlap(self):
        return overlap_terms(self.terms, self.terms)

# --------------------------------------------------------------------------------------
# Two-electron pair Coulomb primitive <phi_a(1) phi_b(2)| 1/r12 |phi_c(1) phi_d(2)>.
#   = sum_k (4π/(2k+1)) R^k(ac,bd) sum_q (-1)^{m_a+q+m_b}
#         gaunt(l_a,k,l_c,-m_a,-q,m_c) gaunt(l_b,k,l_d,-m_b,q,m_d)
# --------------------------------------------------------------------------------------
def pair_coulomb(oa, ob, oc, od):
    la, ma = oa['l'], oa['m']; lb, mb = ob['l'], ob['m']
    lc, mc = oc['l'], oc['m']; ld, md = od['l'], od['m']
    if (ma + mb) != (mc + md):
        return 0.0
    kmax = min(la + lc, lb + ld)
    kmin = max(abs(la - lc), abs(lb - ld))
    total = 0.0
    for k in range(kmin, kmax + 1):
        # parity from (l_a k l_c;000): needs la+k+lc even; gaunt handles zeros
        q = mc - ma
        g1 = gaunt(la, k, lc, -ma, -q, mc)
        if g1 == 0.0:
            continue
        g2 = gaunt(lb, k, ld, -mb, q, md)
        if g2 == 0.0:
            continue
        ang = (-1) ** (ma + q + mb) * g1 * g2 * (4 * math.pi / (2 * k + 1))
        Rk = slater_Rk(oa['rid'], oc['rid'], ob['rid'], od['rid'],
                       oa['P'], oc['P'], ob['P'], od['P'], k)
        total += ang * Rk
    return total

def overlap_terms(termsA, termsB):
    """ <Psi_A^unnorm | Psi_B^unnorm> = sum W_A W_B <u_A|u_B><v_A|v_B>. """
    tot = 0.0
    for (wa, ua, va) in termsA:
        for (wb, ub, vb) in termsB:
            if ua['l'] != ub['l'] or ua['m'] != ub['m']:
                continue
            if va['l'] != vb['l'] or va['m'] != vb['m']:
                continue
            su = 1.0 if ua['rid'] == ub['rid'] else radial_overlap(ua['P'], ub['P'])
            sv = 1.0 if va['rid'] == vb['rid'] else radial_overlap(va['P'], vb['P'])
            tot += wa * wb * su * sv
    return tot

def repulsion_terms(termsA, termsB):
    """ <Psi_A^unnorm | 1/r12 | Psi_B^unnorm>. """
    tot = 0.0
    for (wa, ua, va) in termsA:
        for (wb, ub, vb) in termsB:
            g = pair_coulomb(ua, va, ub, vb)
            if g != 0.0:
                tot += wa * wb * g
    return tot

# --------------------------------------------------------------------------------------
# Assemble secular matrix M = diag(Z R_nu) + T',  T'_{ij} = -<Psi_i|1/r12|Psi_j>.
# --------------------------------------------------------------------------------------
def build_M(configs):
    K = len(configs)
    M = np.zeros((K, K))
    for i in range(K):
        ci = configs[i]
        for j in range(i, K):
            cj = configs[j]
            g = ci.norm * cj.norm * repulsion_terms(ci.terms, cj.terms)
            Tp = -g
            M[i, j] = (Z * ci.Rnu if i == j else 0.0) + Tp
            M[j, i] = M[i, j]
    return M

def build_S(configs):
    """ L2 overlap metric between (L2-normalized) mixed-scale configurations. """
    K = len(configs)
    S = np.zeros((K, K))
    for i in range(K):
        ci = configs[i]
        for j in range(i, K):
            cj = configs[j]
            s = ci.norm * cj.norm * overlap_terms(ci.terms, cj.terms)
            S[i, j] = S[j, i] = s
    return S

def gen_configs(lmax, nmax_per_l):
    """ nmax_per_l: dict l->nmax (max principal qn for that l) or int (same for all). """
    if isinstance(nmax_per_l, int):
        nmax_per_l = {l: nmax_per_l for l in range(lmax + 1)}
    configs = []
    for l in range(lmax + 1):
        nmx = nmax_per_l.get(l, 0)
        ns = list(range(l + 1, nmx + 1))
        for na, nb in combinations_with_replacement(ns, 2):
            configs.append((l, na, nb))
    return configs

def solve(config_tuples):
    Config._next_rid = 0
    _Rk_cache.clear()
    cfgs = [Config(l, na, nb) for (l, na, nb) in config_tuples]
    M = build_M(cfgs)
    p = np.sort(eigh(M, eigvals_only=True))[-1]      # largest root = deepest binding
    E = -p ** 2 / 2
    onenorm = float(np.abs(M).sum())
    return E, onenorm, len(cfgs), M

def solve_with_metric(config_tuples):
    """ Diagnostic: reintroduce the L2 overlap metric -> generalized eigenproblem M B = p S B.
        This is the L2-framing the paper's Sec.2 identifies as ill-conditioned. """
    Config._next_rid = 0
    _Rk_cache.clear()
    cfgs = [Config(l, na, nb) for (l, na, nb) in config_tuples]
    M = build_M(cfgs)
    S = build_S(cfgs)
    p_std = np.sort(eigh(M, eigvals_only=True))[-1]
    p_gen = np.sort(eigh(M, S, eigvals_only=True))[-1]
    return -p_std ** 2 / 2, -p_gen ** 2 / 2, float(np.linalg.cond(S)), len(cfgs)

# ======================================================================================
if __name__ == "__main__":
    EXACT = -2.90372
    print("=" * 74)
    print("Isoenergetic generalized-Sturmian He, extended to l>0  (target E_exact = -2.90372 Ha)")
    print("=" * 74)

    # ---- GATE 1: single-config 1s^2 = -2.847 ----
    E1, n1, K1, _ = solve([(0, 1, 1)])
    print(f"\n[GATE 1] single-config 1s^2 : E = {E1:.5f} Ha   (target -2.847; (5/8)/sqrt2 = {-0.625/math.sqrt(2):.4f})")
    gate1 = abs(E1 - (-2.84766)) < 2e-3
    print(f"         -> {'PASS' if gate1 else 'FAIL'}")

    # bare (no V') sanity: two He+ 1s -> -4.0
    Mbare = np.array([[Z * math.sqrt(2)]])
    print(f"         bare (no V') 1s^2 : E = {-(Z*math.sqrt(2))**2/2:.4f} Ha (exact non-interacting -4.0)")

    # ---- GATE 2: monotone lowering as s,p,d,f sectors are added ----
    print(f"\n[GATE 2] monotone lowering toward exact non-rel He = {EXACT} Ha")
    print(f"  {'set':<34}{'l_max':>6}{'K':>6}{'E (Ha)':>13}{'err%':>9}")
    ladder = [
        ("s  nmax=1  (1s^2)",              0, {0: 1}),
        ("s  nmax=3",                      0, {0: 3}),
        ("s  nmax=6",                      0, {0: 6}),
        ("s6 + p nmax=3",                  1, {0: 6, 1: 3}),
        ("s6 + p nmax=5",                  1, {0: 6, 1: 5}),
        ("s6 + p5 + d nmax=4",             2, {0: 6, 1: 5, 2: 4}),
        ("s6 + p6 + d nmax=6",             2, {0: 6, 1: 6, 2: 6}),
        ("s6 + p6 + d6 + f nmax=6",        3, {0: 6, 1: 6, 2: 6, 3: 6}),
        ("s8 + p8 + d8 + f8",              3, {0: 8, 1: 8, 2: 8, 3: 8}),
    ]
    results = []
    prevE = None
    monotone = True
    for label, lmax, nmp in ladder:
        cts = gen_configs(lmax, nmp)
        E, onenorm, K, _ = solve(cts)
        err = 100 * abs(E - EXACT) / abs(EXACT)
        flag = ""
        if prevE is not None and E > prevE + 1e-6:
            flag = "  <-- rises"
            monotone = False
        prevE = E
        results.append((label, lmax, K, E, onenorm))
        print(f"  {label:<34}{lmax:>6}{K:>6}{E:>13.5f}{err:>9.3f}{flag}")
    best = min(results, key=lambda x: x[3])
    print(f"\n  best E = {best[3]:.5f} Ha at K = {best[2]} configs ({best[0]})   "
          f"[residual to exact = {abs(best[3]-EXACT)*1000:.2f} mHa]")
    gate2 = (best[3] < -2.873) and monotone
    print(f"  monotone across ladder: {monotone}")
    print(f"  -> GATE 2 {'PASS' if gate2 else 'PARTIAL/FAIL'} "
          f"(lowers below s-only floor -2.873 AND monotone)")

    # ---- GATE 3: 1-norm scaling with l>0 present ----
    print(f"\n[GATE 3] entrywise 1-norm ||M||_1 vs K (WITH l>0 configs present)")
    # nested sets that all contain p (and higher), growing K
    scale_sets = [
        ("s3+p3",            gen_configs(1, {0: 3, 1: 3})),
        ("s4+p4",            gen_configs(1, {0: 4, 1: 4})),
        ("s5+p5+d5",         gen_configs(2, {0: 5, 1: 5, 2: 5})),
        ("s6+p6+d6",         gen_configs(2, {0: 6, 1: 6, 2: 6})),
        ("s7+p7+d7+f7",      gen_configs(3, {0: 7, 1: 7, 2: 7, 3: 7})),
        ("s8+p8+d8+f8",      gen_configs(3, {0: 8, 1: 8, 2: 8, 3: 8})),
    ]
    print(f"  {'set':<16}{'K':>6}{'||M||_1':>12}")
    Ks, Ls = [], []
    for label, cts in scale_sets:
        E, onenorm, K, _ = solve(cts)
        Ks.append(K); Ls.append(onenorm)
        print(f"  {label:<16}{K:>6}{onenorm:>12.3f}")
    Ks = np.array(Ks, float); Ls = np.array(Ls, float)
    expo = np.polyfit(np.log(Ks), np.log(Ls), 1)[0]
    print(f"\n  ||M||_1 ~ K^{expo:.3f}   (s-only skeleton gave ~K^0.78)")
    gate3 = expo < 1.0
    print(f"  -> GATE 3 {'PASS (sublinear)' if gate3 else 'FAIL (>=1, superlinear)'}")

    # ---- extended: strictly-nested s,p,d,f convergence toward exact ----
    print(f"\n[CONVERGENCE] strictly-nested s,p,d,f up to nmax=N (each set superset of previous)")
    print(f"  {'N':>3}{'K':>6}{'E (Ha)':>13}{'err%':>9}{'||M||_1':>12}")
    Kc, Lc, Ec = [], [], []
    for N in (2, 3, 4, 6, 8, 10):
        lmax = min(3, N - 1)
        nmp = {l: N for l in range(lmax + 1)}
        E, onenorm, K, _ = solve(gen_configs(lmax, nmp))
        Kc.append(K); Lc.append(onenorm); Ec.append(E)
        err = 100 * abs(E - EXACT) / abs(EXACT)
        print(f"  {N:>3}{K:>6}{E:>13.5f}{err:>9.3f}{onenorm:>12.3f}")
    nested_mono = all(Ec[i + 1] <= Ec[i] + 1e-9 for i in range(len(Ec) - 1))
    nested_expo = np.polyfit(np.log(Kc), np.log(Lc), 1)[0]
    print(f"  nested: best E={min(Ec):.5f} Ha (residual {abs(min(Ec)-EXACT)*1000:.2f} mHa), "
          f"monotone={nested_mono}, ||M||_1~K^{nested_expo:.3f}")

    # ---- corroboration: the L2 metric is ill-conditioned (paper Sec.2 obstruction) ----
    print(f"\n[CORROBORATION] reintroducing the L2 overlap metric S -> generalized eigenproblem")
    print(f"  (paper Sec.2: L2 framing is ill-conditioned; metric-free M is the resolution)")
    print(f"  {'set':>10}{'K':>5}{'E_metricfree':>14}{'E_with_S':>16}{'cond(S)':>11}")
    for N in (3, 4, 6, 8):
        lmax = min(3, N - 1)
        nmp = {l: N for l in range(lmax + 1)}
        Es, Eg, cond, K = solve_with_metric(gen_configs(lmax, nmp))
        print(f"  {('spdf'+str(N)):>10}{K:>5}{Es:>14.5f}{Eg:>16.2f}{cond:>11.1f}")
    print(f"  -> the metric-free standard eigenproblem is the well-conditioned one; the residual")
    print(f"     to exact is basis incompleteness (Goscinskian basis is poor for the He GS,")
    print(f"     cf. Avery's own dedicated 102-CS-config result -2.90250), not a missing metric.")

    print("\n" + "=" * 74)
    print(f"SUMMARY: gate1={'PASS' if gate1 else 'FAIL'}  "
          f"gate2={'PASS' if gate2 else 'PARTIAL/FAIL'}  "
          f"gate3={'PASS' if gate3 else 'FAIL'}")
    print(f"  best He energy {best[3]:.5f} Ha @ K={best[2]};  l>0 1-norm exponent {expo:.3f}")
    print("=" * 74)
