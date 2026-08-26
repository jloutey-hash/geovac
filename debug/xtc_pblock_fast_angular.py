"""
FAST radial-free angular-support / fill-in demonstrator for the p-block xTC test.

Uses an explicit (l,m) angular basis (one radial function per l up to lmax) and
precomputed multipole tables, so we can reach l=2 (d orbitals) where the naive
O(ns^4 * inner) enumerator is intractable.

Support definitions match the engine EXACTLY:
  Coulomb  block (a,b,c,d) nonzero  iff exists (L,M):
      gA(a;L,M;c) != 0   AND   gB(b;L,M;d) != 0        [e1 uses gA, e2 uses gB]
  L3-xTC   block (a,b,c,d) nonzero  iff exists (L,M,Lp,Mp):
      four_Y(a;L,M;Lp,Mp;c) != 0   (vertex, external a<->c)
      AND gA(b;L,M;d) != 0         (leg-2, external b<->d)
      AND refmult(Lp,Mp) != 0      (leg-3 contracted over the reference density)
  refmult(Lp,Mp) = sum_o gA(lo,mo;Lp,Mp;lo,mo)   (diagonal reference-density multipole)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xtc_poc_li_pinclusive import gA, gB
from tc_threebody_collapse_angular import four_Y


def orbs(lmax):
    return [(l, m) for l in range(lmax + 1) for m in range(-l, l + 1)]


def fast_angular(lmax, ref_ang, Lmax=None, tol=1e-10):
    if Lmax is None:
        Lmax = 2 * lmax
    O = orbs(lmax)
    ns = len(O)

    # reference-density multipoles (Lp,Mp) -> coeff (summed over occupied ref orbitals)
    refmult = {}
    for (lo, mo) in ref_ang:
        for Lp in range(Lmax + 1):
            for Mp in range(-Lp, Lp + 1):
                v = gA(lo, mo, Lp, Mp, lo, mo)
                if abs(v) > tol:
                    refmult[(Lp, Mp)] = refmult.get((Lp, Mp), 0.0) + v
    refmult = {k: v for k, v in refmult.items() if abs(v) > tol}

    # per-pair multipole supports
    ac_gA = {}   # (a,c) -> set of (L,M) with gA(a;L,M;c)!=0
    bd_gB = {}   # (b,d) -> set of (L,M) with gB(b;L,M;d)!=0   (Coulomb e2)
    for i, (la, ma) in enumerate(O):
        for j, (lc, mc) in enumerate(O):
            sA = set(); sB = set()
            for L in range(Lmax + 1):
                for M in range(-L, L + 1):
                    if abs(gA(la, ma, L, M, lc, mc)) > tol:
                        sA.add((L, M))
                    if abs(gB(la, ma, L, M, lc, mc)) > tol:
                        sB.add((L, M))
            ac_gA[(i, j)] = sA
            bd_gB[(i, j)] = sB

    # vertex reachability: for (a,c), which leg-2 multipoles (L,M) are reachable
    # given the reference provides some (Lp,Mp)?
    reach = {}
    for i, (la, ma) in enumerate(O):
        for j, (lc, mc) in enumerate(O):
            s = set()
            for (Lp, Mp) in refmult:
                for L in range(Lmax + 1):
                    for M in range(-L, L + 1):
                        v, _ = four_Y(la, ma, L, M, Lp, Mp, lc, mc)
                        if abs(v) > tol:
                            s.add((L, M))
            reach[(i, j)] = s

    coul = 0; l3 = 0; fill = 0; total = ns ** 4
    fill_examples = []
    for a in range(ns):
        for b in range(ns):
            for c in range(ns):
                for d in range(ns):
                    # m-conservation
                    if O[a][1] + O[b][1] != O[c][1] + O[d][1]:
                        continue
                    is_coul = bool(ac_gA[(a, c)] & bd_gB[(b, d)])
                    is_l3 = bool(reach[(a, c)] & ac_gA[(b, d)])
                    if is_coul:
                        coul += 1
                    if is_l3:
                        l3 += 1
                        if not is_coul:
                            fill += 1
                            if len(fill_examples) < 10:
                                fill_examples.append(
                                    (O[a], O[b], O[c], O[d]))
    return dict(ns=ns, total=total, n_coulomb=coul, n_l3=l3, n_fill_in=fill,
                coul_density=coul / total, l3_density=l3 / total,
                fill_frac=fill / total,
                refmult=sorted(refmult.keys()),
                fill_examples=[str(x) for x in fill_examples])


# reference angular occupations (l,m) with multiplicity
REF_ANG = {
    'Be_s2 (s-ref)':      [(0, 0), (0, 0), (0, 0), (0, 0)],          # 1s^2 2s^2
    'C_2p2_Hund (p-ref)': [(0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (1, 1)],
    'C_2p0sq (p-ref)':    [(0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (1, 0)],
    'O_2p4 (p-ref)':      [(0, 0), (0, 0), (0, 0), (0, 0), (1, -1), (1, -1), (1, 0), (1, 0)],
}

if __name__ == '__main__':
    import json
    out = {}
    for lmax in (1, 2):
        tag = {1: 's+p', 2: 's+p+d'}[lmax]
        print(f"\n===== basis lmax={lmax} ({tag}), {len(orbs(lmax))} angular orbitals =====")
        out[tag] = {}
        for name, ra in REF_ANG.items():
            r = fast_angular(lmax, ra)
            out[tag][name] = r
            print(f"  {name:22s} refmult={r['refmult']}")
            print(f"      Coulomb {r['n_coulomb']}/{r['total']} ({r['coul_density']*100:.1f}%)"
                  f"   L3 {r['n_l3']} ({r['l3_density']*100:.1f}%)"
                  f"   FILL-IN {r['n_fill_in']} ({r['fill_frac']*100:.2f}%)")
            if r['n_fill_in']:
                print(f"      examples (a,b,c,d)=(l,m): {r['fill_examples'][:4]}")
    os.makedirs('debug/data', exist_ok=True)
    with open('debug/data/xtc_pblock_fast_angular.json', 'w') as f:
        json.dump(out, f, indent=2)
    print("\nwrote debug/data/xtc_pblock_fast_angular.json")
