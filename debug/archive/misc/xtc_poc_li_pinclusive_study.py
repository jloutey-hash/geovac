"""
Validation gates + the deliverable sparsity study for the p-inclusive xTC PoC.

Deliverable: does the xTC-contracted effective 2-body operator inherit the angular
Gaunt sparsity, or does the 3-body -> 2-body contraction FILL IN the zero blocks?
Measured as angular density (nonzero spatial Gaunt blocks / total), 1-norm (LCU
lambda proxy), and Pauli count -- PLAIN (Coulomb) vs xTC (w + contracted L3).
"""
import os, sys, json, time
import numpy as np
from math import pi

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE); sys.path.insert(0, os.path.abspath(os.path.join(_HERE, '..')))
import xtc_poc_li_pinclusive as P
from xtc_poc_li_pinclusive import (build_states, radial_table, build_overlap,
                                    build_eri_scalar, coulomb_ML, gA, gB,
                                    angular_density, pauli_count, eri_phys_spinorbital)
from xtc_poc_li import make_grid
from tc_threebody_collapse_angular import four_Y
from geovac.sturmian_solver import SturmianCI, _slater_rk, hydrogenic_radial

TRUE_LI = -7.478060


def pure_angular_supports(max_n=2, ref_ang=(0, 1)):
    """EXACT, radial-free angular support (the model-independent fill-in test):
    Coulomb 2-body vs contracted-L3 (correlator leg contracted over the s-reference).
    Reuses Track 1's four_Y + the correct gaunt (gA legs, gB coulomb-e2)."""
    states = build_states(max_n); ns = len(states); Lmax = 2

    def coulomb_ang(a, b, c, d):
        la, ma = states[a][1:]; lb, mb = states[b][1:]
        lc, mc = states[c][1:]; ld, md = states[d][1:]
        if ma + mb != mc + md: return 0.0
        M = mc - ma; tot = 0.0
        for L in range(5):
            a1 = gA(la, ma, L, M, lc, mc)
            if a1 == 0.0: continue
            b1 = gB(lb, mb, L, M, ld, md)
            if b1 == 0.0: continue
            tot += a1 * b1 / (2 * L + 1)
        return tot

    def l3_contracted(a, b, c, d):   # vertex (a,c), leg2 (b,d), leg3 over s-ref
        la, ma = states[a][1:]; lb, mb = states[b][1:]
        lc, mc = states[c][1:]; ld, md = states[d][1:]
        tot = 0.0
        for o in ref_ang:
            lo, mo = states[o][1:]
            for L in range(Lmax + 1):
                for M in range(-L, L + 1):
                    for Lp in range(Lmax + 1):
                        for Mp in range(-Lp, Lp + 1):
                            fy, _ = four_Y(la, ma, L, M, Lp, Mp, lc, mc)
                            if fy == 0.0: continue
                            g2 = gA(lb, mb, L, M, ld, md)
                            if g2 == 0.0: continue
                            g3 = gA(lo, mo, Lp, Mp, lo, mo)
                            if g3 == 0.0: continue
                            tot += fy * g2 * g3
        return tot

    coul = set(); l3 = set(); badm = 0
    for a in range(ns):
        for b in range(ns):
            for c in range(ns):
                for d in range(ns):
                    if abs(coulomb_ang(a, b, c, d)) > 1e-10:
                        coul.add((a, b, c, d))
                    if abs(l3_contracted(a, b, c, d)) > 1e-10:
                        l3.add((a, b, c, d))
                        if states[a][2] + states[b][2] != states[c][2] + states[d][2]:
                            badm += 1
    return dict(n_coulomb=len(coul), n_l3=len(l3),
                n_fill_in=len(l3 - coul), n_coul_not_l3=len(coul - l3),
                m_violations=badm, total=ns ** 4,
                coul_density=len(coul) / ns ** 4, l3_density=len(l3) / ns ** 4)


def threebody_density(o):
    V3 = o['_arr']['V3o']; ns = o['ns']
    nz = int(np.sum(np.abs(V3) > 1e-9))
    return dict(nnz=nz, total=ns ** 6, density=nz / ns ** 6)


# ----------------------------------------------------------------------------
def gate_coulomb_radial(k=1.5):
    """Validate my grid-multipole Coulomb ERI against a HIGH-ACCURACY reference
    (analytic 5k/8 for the 1s block + adaptive-quad radial multipole for the rest).
    NB: the framework's _slater_rk uses a coarse 500-pt linspace and is itself off
    by ~6e-2 on the 1s block (documented in the s-only memo); we validate against
    the accurate value, not against _slater_rk."""
    from scipy.integrate import quad
    states = build_states(2); ns = len(states)
    r, wr = make_grid(k, Ng=1400); W = r * r * wr
    R = radial_table(states, r, k)
    eri_grid = build_eri_scalar(states, R, W, coulomb_ML(r, 4), 4)

    def Rk_accurate(na, la, nc, lc, nb, lb, nd, ld, L):
        f1 = lambda x: hydrogenic_radial(np.array([x]), na, la, na*k)[0] * \
                       hydrogenic_radial(np.array([x]), nc, lc, nc*k)[0]
        f2 = lambda x: hydrogenic_radial(np.array([x]), nb, lb, nb*k)[0] * \
                       hydrogenic_radial(np.array([x]), nd, ld, nd*k)[0]
        rmax = 80.0 / k
        def yk(r1):
            lo = quad(lambda r2: f2(r2)*r2**2 * r2**L, 0, r1, limit=100)[0] / r1**(L+1)
            hi = quad(lambda r2: f2(r2)*r2**2 / r2**(L+1), r1, rmax, limit=100)[0] * r1**L
            return lo + hi
        return quad(lambda r1: f1(r1)*r1**2 * yk(r1), 0, rmax, limit=100)[0]

    def eri_ref(a, b, c, d):
        na, la, ma = states[a]; nb, lb, mb = states[b]
        nc, lc, mc = states[c]; nd, ld, md = states[d]
        if ma + mb != mc + md:
            return 0.0
        M = mc - ma; tot = 0.0
        for L in range(5):
            av = gA(la, ma, L, M, lc, mc)
            if av == 0.0: continue
            bv = gB(lb, mb, L, M, ld, md)
            if bv == 0.0: continue
            Rk = Rk_accurate(na, la, nc, lc, nb, lb, nd, ld, L)
            tot += (4*pi/(2*L+1)) * av * bv * Rk
        return tot
    mx = 0.0
    for t in [(0,0,0,0),(0,1,0,1),(2,2,2,2),(0,2,0,2),(4,2,3,3),(4,2,2,4),(3,3,2,4)]:
        mx = max(mx, abs(eri_grid[t] - eri_ref(*t)))
    return mx


def gate_framework_mconserving(k=1.5):
    """My correct-gaunt Coulomb matches framework production ERI on m-conserving
    blocks; report how many physical m-CHANGING blocks the framework drops."""
    states = build_states(2); ns = len(states)
    r, wr = make_grid(k, Ng=1400); W = r * r * wr
    R = radial_table(states, r, k)
    eri = build_eri_scalar(states, R, W, coulomb_ML(r, 4), 4)
    ci = SturmianCI(Z=3, n_electrons=3, max_n=2); fw = ci._build_eri(k)
    max_mcons_diff = 0.0; dropped = 0
    for a in range(ns):
        for b in range(ns):
            for c in range(ns):
                for d in range(ns):
                    mine = eri[a, b, c, d]
                    theirs = fw.get((a, b, c, d), 0.0)
                    mcons = (states[a][2] + states[b][2] == states[c][2] + states[d][2])
                    dm = (states[a][2] - states[c][2])  # e1 multipole M
                    if abs(theirs) > 1e-12:  # framework nonzero -> must match
                        max_mcons_diff = max(max_mcons_diff, abs(mine - theirs))
                    elif abs(mine) > 1e-9 and mcons:
                        dropped += 1  # physical block framework drops
    return max_mcons_diff, dropped


def gate_v3_mconservation(o):
    states = o['states']; V3 = o['_arr']['V3o']
    ns = len(states)
    bad = 0; nz = 0
    nzidx = np.argwhere(np.abs(V3) > 1e-9)
    for (a, b, c, d, e, f) in nzidx:
        nz += 1
        if (states[a][2] + states[b][2] + states[c][2]
                != states[d][2] + states[e][2] + states[f][2]):
            bad += 1
    return bad, nz


def gate_geminal_zero(k=1.5):
    rows = []
    for g in [4.0, 10.0, 25.0]:
        o = P.assemble(max_n=2, k=k, gamma=g, Ng=1000, nx=160, want_exact3=True)
        rows.append(dict(gamma=g, E_plain=o['E_plain'], E_xTC=o['E_xTC'],
                         diff=abs(o['E_xTC'] - o['E_plain'])))
    return rows


def gate_contraction_fidelity(o):
    return abs(o['E_xTC'] - o['E_exactTC'])


# ----------------------------------------------------------------------------
# DELIVERABLE: angular density + 1-norm + Pauli, plain vs xTC
# ----------------------------------------------------------------------------
def deliverable(o):
    A = o['_arr']; ns = o['ns']; nso = o['nso']
    states = o['states']

    # ---- (i) pure spatial angular support: Coulomb vs contracted-L3 -----------
    # Coulomb spatial support (transformed orthonormal MO 2-body)
    coul_sp = angular_density(A['eri_coul_o'])
    # contracted-L3 spatial physicist 2-body: contract V3o over ref spin-orbitals,
    # reduce to spatial by symmetric spin-average of the (non-antisymmetrized) piece.
    # Simpler + rigorous at operator level: read v2 (antisymmetrized spin-orbital).
    v2 = A['v2']
    coul_asym = A['asym_coul']
    # spin-orbital supports
    supp_coul = set(map(tuple, np.argwhere(np.abs(coul_asym) > 1e-9)))
    supp_v2 = set(map(tuple, np.argwhere(np.abs(v2) > 1e-9)))
    fill_in = supp_v2 - supp_coul
    # spatial fold of the fill-in (which (l,m) blocks)
    def to_spatial(idx):
        return tuple((states[i // 2][1], states[i // 2][2]) for i in idx)
    fill_spatial_lm = sorted({to_spatial(x) for x in fill_in})

    # ---- (ii) operator 1-norm + nnz (spin-orbital, what gets JW'd) ------------
    m_plain = P_op_metrics(A['hso'], A['asym_coul'], nso)
    asym_x = A['asym_w'] + v2
    m_xtc = P_op_metrics(A['hso_x'], asym_x, nso)
    m_v2 = P_op_metrics(np.zeros((nso, nso)), v2, nso)

    # ---- (iii) Pauli counts via openfermion JW (Hermitian ops) ----------------
    pauli = {}
    try:
        # plain: Coulomb physicist spin-orbital
        eri_c_so = eri_phys_spinorbital(A['eri_coul_o'], nso)
        pauli['plain'] = pauli_count(A['hso'], eri_c_so)
        # xTC Hermitian part: w physicist + contracted-L3.  v2 is antisymmetrized;
        # to feed OF we need a physicist 2-body.  Build eff physicist = w-phys plus
        # a physicist realization of v2 (v2 -> symmetric physicist tensor).
        eri_w_so = eri_phys_spinorbital(A['eri_w_o'], nso)
        v2_phys = asym_to_phys(v2, nso)
        pauli['xTC'] = pauli_count(A['hso_x'], eri_w_so + v2_phys)
    except Exception as e:
        pauli['error'] = repr(e)

    return dict(
        coulomb_spatial=dict(nnz=coul_sp['nnz'], total=coul_sp['total'],
                             density=coul_sp['density'], l1=coul_sp['l1']),
        v2_vs_coulomb=dict(n_supp_coul=len(supp_coul), n_supp_v2=len(supp_v2),
                           n_fill_in=len(fill_in),
                           n_v2_within_coul=len(supp_v2 & supp_coul),
                           fill_spatial_lm_blocks=[str(x) for x in fill_spatial_lm]),
        op_metrics=dict(plain=m_plain, xTC=m_xtc, contractedL3_v2=m_v2),
        pauli=pauli,
    )


def P_op_metrics(h, asym, nso):
    return dict(nnz_1body=int(np.sum(np.abs(h) > 1e-10)),
                l1_1body=float(np.sum(np.abs(h))),
                nnz_2body=int(np.sum(np.abs(asym) > 1e-10)),
                l1_2body=float(np.sum(np.abs(asym))),
                l1_total=float(np.sum(np.abs(h)) + 0.25 * np.sum(np.abs(asym))))


def asym_to_phys(asym, nso):
    """Recover a physicist <pq|rs> from <pq||rs> = <pq|rs> - <pq|sr> by symmetric
    halving (valid because our v2 is a genuine 2-body with <pq|rs>=<rs|pq> and the
    only structure is direct-minus-exchange). Use V[p,q,r,s] = 1/2 asym[p,q,r,s]
    which reproduces the SAME antisymmetrized operator (1/4 sum asym a^a^aa)."""
    # (1/4) sum asym[p,q,r,s] a_p^ a_q^ a_s a_r   ==   (1/2) sum Vphys a_p^ a_q^ a_s a_r
    # with Vphys = 1/2 asym gives identical operator after OF antisymmetrization.
    return 0.5 * asym


# ============================================================================
def main():
    t0 = time.time()
    rep = {'true_Li': TRUE_LI}
    print("=== VALIDATION GATES ===")
    g_rad = gate_coulomb_radial()
    print(f"[G1] grid-multipole Coulomb vs _slater_rk : max|d| = {g_rad:.2e}  (want <1e-4)")
    mcons, dropped = gate_framework_mconserving()
    print(f"[G2] vs framework production ERI: max|d| on m-conserving = {mcons:.2e}; "
          f"framework DROPS {dropped} physical m-changing blocks we keep")

    o = P.assemble(max_n=2, k=1.5, gamma=1.0, Ng=1200, nx=200, want_exact3=True)
    badm, nzv3 = gate_v3_mconservation(o)
    print(f"[G3] V3 m-conservation: {badm}/{nzv3} nonzero entries violate total-m "
          f"(must be 0)")
    fid = gate_contraction_fidelity(o)
    print(f"[G4] xTC contraction fidelity |E_xTC - E_exactTC| = {fid:.2e}  "
          f"(E_xTC={o['E_xTC']:.6f}, E_exactTC={o['E_exactTC']:.6f})")
    gz = gate_geminal_zero()
    print("[G5] geminal->0 (E_xTC -> E_plain):")
    for row in gz:
        print(f"     g={row['gamma']:5.1f}: plain={row['E_plain']:.6f} "
              f"xTC={row['E_xTC']:.6f}  |d|={row['diff']:.2e}")
    rep['gates'] = dict(coulomb_radial=g_rad, framework_mcons_diff=mcons,
                        framework_dropped_blocks=dropped,
                        v3_m_violations=badm, v3_nonzero=nzv3,
                        contraction_fidelity=fid, geminal_zero=gz)

    print("\n=== EXACT ANGULAR SUPPORT (radial-free, model-independent) ===")
    pa = pure_angular_supports()
    print(f"  Coulomb 2-body angular blocks : {pa['n_coulomb']}/{pa['total']} "
          f"({pa['coul_density']*100:.2f}%)")
    print(f"  contracted-L3 2-body blocks   : {pa['n_l3']}/{pa['total']} "
          f"({pa['l3_density']*100:.2f}%)  [m-violations {pa['m_violations']}]")
    print(f"  >>> FILL-IN (L3 not in Coulomb): {pa['n_fill_in']}   "
          f"(Coulomb not in L3: {pa['n_coul_not_l3']}) <<<")
    tb = threebody_density(o)
    print(f"  raw 3-body L3 tensor density   : {tb['nnz']}/{tb['total']} "
          f"({tb['density']*100:.2f}%)  [cf. Track 1 ~5% silver lining]")
    rep['pure_angular'] = pa; rep['threebody_density'] = tb

    print("\n=== DELIVERABLE: angular density + 1-norm + Pauli (plain vs xTC) ===")
    D = deliverable(o)
    rep['deliverable'] = D
    rep['energies'] = dict(E_plain=o['E_plain'], E_TC2=o['E_TC2'],
                           E_exactTC=o['E_exactTC'], E_xTC=o['E_xTC'],
                           imag_xTC=o['imag_xTC'], v0=o['xtc_v0'])

    cs = D['coulomb_spatial']
    print(f"  Coulomb spatial 2-body: nnz={cs['nnz']}/{cs['total']} "
          f"(density {cs['density']*100:.2f}%)  1-norm={cs['l1']:.3f}")
    vc = D['v2_vs_coulomb']
    print(f"  contracted-L3 v2 support: {vc['n_supp_v2']} spin-orb blocks; "
          f"Coulomb asym support: {vc['n_supp_coul']}")
    print(f"  >>> FILL-IN (v2 nonzero where Coulomb zero): {vc['n_fill_in']} <<<")
    print(f"  v2 blocks within Coulomb support: {vc['n_v2_within_coul']}")
    if vc['fill_spatial_lm_blocks']:
        print(f"      fill-in (l,m) blocks: {vc['fill_spatial_lm_blocks']}")
    om = D['op_metrics']
    print(f"  1-norm[2body]  plain={om['plain']['l1_2body']:.3f}  "
          f"xTC={om['xTC']['l1_2body']:.3f}  "
          f"(ratio {om['xTC']['l1_2body']/om['plain']['l1_2body']:.3f})  "
          f"| v2 alone={om['contractedL3_v2']['l1_2body']:.3f}")
    print(f"  1-norm[total]  plain={om['plain']['l1_total']:.3f}  "
          f"xTC={om['xTC']['l1_total']:.3f}  "
          f"(ratio {om['xTC']['l1_total']/om['plain']['l1_total']:.3f})")
    if 'error' not in D['pauli']:
        pp, px = D['pauli']['plain'], D['pauli']['xTC']
        print(f"  Pauli terms    plain={pp['n_pauli']}  xTC={px['n_pauli']}  "
              f"(ratio {px['n_pauli']/pp['n_pauli']:.3f})")
        print(f"  Pauli 1-norm   plain={pp['l1_pauli']:.3f}  xTC={px['l1_pauli']:.3f}  "
              f"(ratio {px['l1_pauli']/pp['l1_pauli']:.3f})")
    else:
        print("  Pauli:", D['pauli']['error'])

    print("\n=== ROBUSTNESS across geminal width gamma (fill-in + 1-norm ratio) ===")
    rep['gamma_scan'] = []
    for g in [0.6, 1.0, 1.5]:
        og = o if abs(g - 1.0) < 1e-9 else P.assemble(max_n=2, k=1.5, gamma=g,
                                                       Ng=1200, nx=200, want_exact3=True)
        Dg = deliverable(og)
        row = dict(gamma=g,
                   n_fill_in=Dg['v2_vs_coulomb']['n_fill_in'],
                   l1_2b_ratio=Dg['op_metrics']['xTC']['l1_2body'] /
                               Dg['op_metrics']['plain']['l1_2body'],
                   l1_tot_ratio=Dg['op_metrics']['xTC']['l1_total'] /
                                Dg['op_metrics']['plain']['l1_total'],
                   d3_mHa=(og['E_exactTC'] - og['E_TC2']) * 1e3,
                   xtc_fidelity_mHa=abs(og['E_xTC'] - og['E_exactTC']) * 1e3)
        if 'error' not in Dg['pauli']:
            row['pauli_ratio'] = Dg['pauli']['xTC']['n_pauli'] / Dg['pauli']['plain']['n_pauli']
            row['pauli_l1_ratio'] = Dg['pauli']['xTC']['l1_pauli'] / Dg['pauli']['plain']['l1_pauli']
        rep['gamma_scan'].append(row)
        print(f"  g={g:.2f}: fill-in={row['n_fill_in']}  "
              f"L1[2b] ratio={row['l1_2b_ratio']:.3f}  L1[tot] ratio={row['l1_tot_ratio']:.3f}  "
              f"Pauli ratio={row.get('pauli_ratio', float('nan')):.3f}  "
              f"3body={row['d3_mHa']:+.2f}mHa (xTC fid {row['xtc_fidelity_mHa']:.3f}mHa)")

    rep['walltime_s'] = time.time() - t0
    os.makedirs(os.path.join(_HERE, 'data'), exist_ok=True)
    with open(os.path.join(_HERE, 'data', 'xtc_poc_li_pinclusive_study.json'), 'w') as fh:
        json.dump(rep, fh, indent=2, default=lambda x: list(x) if isinstance(x, set) else float(x))
    print(f"\nwrote debug/data/xtc_poc_li_pinclusive_study.json  ({rep['walltime_s']:.0f}s)")


if __name__ == '__main__':
    main()
