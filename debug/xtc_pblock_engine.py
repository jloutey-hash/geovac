"""
p-BLOCK reference xTC engine helpers (the DECISIVE sparsity test).
===================================================================

The atomic p-INCLUSIVE run (`xtc_poc_li_pinclusive*`, memo
`sprint_xtc_pinclusive_memo.md`) proved the xTC-contracted 2-body inherits
GeoVac's angular Gaunt sparsity with **0 fill-in** -- but ONLY because Li's
reference (1s^2 2s) is all-s.  Contracting a correlator leg over an s orbital
(l=0) forces its multipole L'=0 (monopole -> Coulomb support), so the shared
vertex collapses to a single multipole and cannot fill in the angular zeros.

A **p-BLOCK reference** (occupied 2p, l=1) keeps L' in {0,2}: the leg no longer
collapses to a monopole, so the 3-body -> 2-body contraction CAN fill in the
angular zero-blocks.  This is where "does the sparsity survive xTC in general"
is genuinely at risk.

Deep subtlety (why the task says OPEN p-shells C 2p^2 / O 2p^4):
  The collapse is really about the reference DENSITY being spherically symmetric.
  An s orbital and a CLOSED 2p^6 shell are BOTH spherical (sum_m |Y_1m|^2 = const)
  -> both collapse leg-3 to L'=0 -> no fill-in.  Only an OPEN p-shell breaks the
  symmetry (2p^2 ~ cos^2 theta etc.), switching the L'=2 leg ON.  So carbon /
  oxygen are exactly the references that put the sparsity win at risk.

This module only adds:
  (1) reference-occupation builders for p-block atoms, and
  (2) a generalized, radial-free (model-independent) angular-support test that
      contracts leg-3 over an ARBITRARY reference angular occupation.
Everything else (radial, contraction, PN-projected FCI, 1-norm, Pauli) is reused
verbatim from `xtc_poc_li_pinclusive.assemble`.
"""
import os, sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
_ROOT = os.path.abspath(os.path.join(_HERE, '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from xtc_poc_li_pinclusive import build_states, gA, gB
from tc_threebody_collapse_angular import four_Y

# spatial-orbital indexing in the max_n=2 s+p basis:
#   0=1s  1=2s  2=2p(-1)  3=2p(0)  4=2p(+1)
# spin-orbital p = 2*spatial + spin  (spin 0 = up, 1 = down)
P1S_UP, P1S_DN = 0, 1
P2S_UP, P2S_DN = 2, 3
P2PM_UP, P2PM_DN = 4, 5     # 2p(-1)
P2P0_UP, P2P0_DN = 6, 7     # 2p(0)
P2PP_UP, P2PP_DN = 8, 9     # 2p(+1)


# ----------------------------------------------------------------------------
# Reference occupations (single-determinant reference 1-RDMs, diagonal).
# The xTC contraction sums leg-3 over these occupied spin-orbitals.
# ----------------------------------------------------------------------------
REFERENCES = {
    # Carbon 1s^2 2s^2 2p^2, Hund 3P (M_L=+1, M_S=+1): 2p0 up + 2p+1 up.
    # OPEN p-shell -> non-spherical density -> the decisive fill-in test.
    'C_3P':   dict(Z=6, n_elec=6, ref_occ=(P1S_UP, P1S_DN, P2S_UP, P2S_DN,
                                           P2P0_UP, P2PP_UP)),
    # Carbon 2p0^2 (closed m=0 sub-block): still OPEN shell (density ~ cos^2 th),
    # a second open-p reference to show the fill-in is not a single-config artifact.
    'C_2p0sq': dict(Z=6, n_elec=6, ref_occ=(P1S_UP, P1S_DN, P2S_UP, P2S_DN,
                                            P2P0_UP, P2P0_DN)),
    # Oxygen 1s^2 2s^2 2p^4, 3P (two holes; particle-hole partner of C 2p^2):
    # 2p-1 (up+dn) 2p0 (up+dn) 2p+1 up + one more; here 2p full minus (2p+1 dn):
    # occupy 2p-1^2 2p0^2 2p+1^1 (= 2p^5?)  -> use genuine 2p^4: 2p-1 ud, 2p0 ud.
    'O_3P':   dict(Z=8, n_elec=8, ref_occ=(P1S_UP, P1S_DN, P2S_UP, P2S_DN,
                                           P2PM_UP, P2PM_DN, P2P0_UP, P2P0_DN)),
    # s-reference CONTROL rebuilt in THIS basis (Be-like 1s^2 2s^2, all-s):
    # must reproduce Li's 0 fill-in (spherical reference).
    'Be_1S':  dict(Z=4, n_elec=4, ref_occ=(P1S_UP, P1S_DN, P2S_UP, P2S_DN)),
}


def ref_angular(ref_occ, states):
    """Map occupied spin-orbitals -> list of (l,m) with occupation multiplicity.
    (Leg-3 is contracted over these; multiplicity is the diagonal 1-RDM weight.)"""
    out = []
    for p in ref_occ:
        sp = p >> 1
        out.append((states[sp][1], states[sp][2]))
    return out


# ----------------------------------------------------------------------------
# Radial-free, model-independent angular support test.
#   Coulomb 2-body angular support  vs  contracted-L3 2-body angular support,
#   where leg-3 is summed over the reference angular occupation (with L'=0 collapse
#   switched OFF whenever a p orbital is occupied).
# ----------------------------------------------------------------------------
def pure_angular_supports(max_n, ref_ang, Lmax=2, tol=1e-10):
    states = build_states(max_n)
    ns = len(states)

    def coulomb_ang(a, b, c, d):
        la, ma = states[a][1:]; lb, mb = states[b][1:]
        lc, mc = states[c][1:]; ld, md = states[d][1:]
        if ma + mb != mc + md:
            return 0.0
        M = mc - ma; tot = 0.0
        for L in range(2 * max_n + 1):
            a1 = gA(la, ma, L, M, lc, mc)
            if a1 == 0.0:
                continue
            b1 = gB(lb, mb, L, M, ld, md)
            if b1 == 0.0:
                continue
            tot += a1 * b1 / (2 * L + 1)
        return tot

    def l3_contracted(a, b, c, d):
        # vertex (a,c), leg-2 (b,d), leg-3 contracted over the reference (lo,mo)
        la, ma = states[a][1:]; lb, mb = states[b][1:]
        lc, mc = states[c][1:]; ld, md = states[d][1:]
        tot = 0.0
        for (lo, mo) in ref_ang:
            for L in range(Lmax + 1):
                for M in range(-L, L + 1):
                    for Lp in range(Lmax + 1):
                        for Mp in range(-Lp, Lp + 1):
                            fy, _ = four_Y(la, ma, L, M, Lp, Mp, lc, mc)
                            if fy == 0.0:
                                continue
                            g2 = gA(lb, mb, L, M, ld, md)
                            if g2 == 0.0:
                                continue
                            g3 = gA(lo, mo, Lp, Mp, lo, mo)   # leg-3 on the reference
                            if g3 == 0.0:
                                continue
                            tot += fy * g2 * g3
        return tot

    coul = set(); l3 = set(); badm = 0
    for a in range(ns):
        for b in range(ns):
            for c in range(ns):
                for d in range(ns):
                    if abs(coulomb_ang(a, b, c, d)) > tol:
                        coul.add((a, b, c, d))
                    if abs(l3_contracted(a, b, c, d)) > tol:
                        l3.add((a, b, c, d))
                        if states[a][2] + states[b][2] != states[c][2] + states[d][2]:
                            badm += 1
    fill = l3 - coul
    # which (l,m) content the fill-in touches (on the leg pair b,d = the L'=2 leg)
    fill_lm = sorted({((states[a][1], states[a][2]), (states[b][1], states[b][2]),
                       (states[c][1], states[c][2]), (states[d][1], states[d][2]))
                      for (a, b, c, d) in fill})
    return dict(n_coulomb=len(coul), n_l3=len(l3),
                n_fill_in=len(fill), n_coul_not_l3=len(coul - l3),
                m_violations=badm, total=ns ** 4,
                coul_density=len(coul) / ns ** 4, l3_density=len(l3) / ns ** 4,
                fill_frac=len(fill) / ns ** 4,
                fill_examples=[str(x) for x in fill_lm[:12]])


if __name__ == '__main__':
    states = build_states(2)
    for name in ('Be_1S', 'C_3P', 'C_2p0sq', 'O_3P'):
        spec = REFERENCES[name]
        ra = ref_angular(spec['ref_occ'], states)
        pa = pure_angular_supports(2, ra)
        print(f"{name:8s} ref_ang={ra}")
        print(f"   Coulomb {pa['n_coulomb']}/{pa['total']} ({pa['coul_density']*100:.1f}%)"
              f"  L3 {pa['n_l3']}/{pa['total']} ({pa['l3_density']*100:.1f}%)"
              f"  FILL-IN {pa['n_fill_in']} ({pa['fill_frac']*100:.1f}%)  mviol {pa['m_violations']}")
