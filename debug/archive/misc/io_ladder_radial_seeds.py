"""I/O ladder accounting -- Rung 2 (radial-seed compactness, the drill).

DIAGNOSTIC ONLY. Sharpens Rung 1 (`debug/io_ladder_accounting.py`, shell-quartet
granularity) to the granularity Rung 1's own section 8 named as the next step:
the count of DISTINCT closed-form radial-seed INSTANCES the two-centre engine
(`geovac/two_center_eri.py`) actually emits, vs the naive per-entry count. Also
prices the three-centre elliptic-Bessel-moment class (Paper 59 / routeC) onto
the same ledger, as Rung 1 section 8 named.

Axis under test (Avery's device-I/O thesis): for a Coulomb-Sturmian / hydrogenic
encoding, the ANGULAR structure is generated on-device from (n,l,m) labels via
Gaunt/3j; only 1D RADIAL factors need classical tabulation + shipping. This rung
asks: how big is that radial table, REALLY -- not at shell-quartet granularity,
but at the granularity of the actual scalar arguments fed to {E1, ln, exp, gamma}?

NOT the same axis as:
  - N3b (tensor-density sparsity dies past one centre) -- this counts distinct
    SEED VALUES, not tensor nonzeros.
  - QC-1 (contracted Gaussian matches Slater at matched M; no qubit/Pauli win)
    -- this counts classical radial I/O, not qubits or Pauli terms.

======================================================================
PART A -- TWO-CENTRE ENGINE (geovac/two_center_eri.py)
======================================================================

The engine's own structure (multipole_decomposition / radial_product) makes the
mechanism explicit: the combined radial DECAY RATE of a same-centre orbital pair
(n1,l1,m1)-(n2,l2,m2) is

    b = Z/n1 + Z/n2

-- a function of (Z, n1, n2) ONLY. l1, l2, m1, m2 do not appear (confirmed by
reading `radial_product`: it sums `radial_poly(Z,n,l)`'s decay rate Z/n over the
pair, and that rate has no l-dependence). So every transcendental argument the
engine ever builds -- E1(mu*R), ln(rate ratio), the plain exp(-lambda*R)
factors -- is a function of a HANDFUL of principal-quantum-number-indexed rates,
reused across every (l,m) that shares the same n's. This is the SAME mechanism
Rung 1 found for the angular sector, one level down: the RADIAL content itself
further collapses from "one seed per shell-quartet" to "one seed per rate-tuple",
because many shell quartets (differing only in l, not n) share a rate-tuple.

Four classes, by how the four orbital indices split across the two centres A, B
(one-centre-only quartets, both orbitals same atom, are OUT OF SCOPE for this
engine -- they are ordinary same-centre Slater R_k integrals, a different, already
understood machinery, not part of `two_center_eri.py`):

  AA|BB     bra=(A,A) ket=(B,B).  E1-FREE (structurally, see file header) --
            seed = (bA, bB), bA=Z_A(1/n1+1/n2), bB=Z_B(1/n3+1/n4).
  HYBRID    3 orbitals on one centre (the multipole pair PLUS a third one-centre
            orbital), 1 on the other.  seed = (mu, ad), mu=Z_trio*sum(1/n_i) over
            the THREE trio orbitals, ad=Z_lone/n_lone.  Emits E1(mu R),
            E1((mu+ad) R), ln((mu+ad)/mu) -- Increment 2 docstring, "EQ2".
  EXCHANGE  2+2 split, one A + one B orbital in EACH of bra and ket.  seed =
            (alpha_bra, beta_bra, alpha_ket, beta_ket), alpha=Z_A/n_A-orbital,
            beta=Z_B/n_B-orbital (raw single-orbital rates -- see
            `two_center_spheroidal_product`, alpha/beta come straight from
            `radial_poly`, no l-dependence).  Emits {E1, ln, gamma} (Increment 3c).

Verified directly against the closed-form code (PART 0 below) before trusting the
structural counting: two quartets sharing (n1,n2) but differing l -> identical
rate value; the actual `aabb_closed_form`/`hybrid_closed_form` output's exp/E1/log
ARGUMENTS match the predicted rate exactly.

======================================================================
PART B -- THREE-CENTRE ENGINE (routeC elliptic Bessel moment, Paper 59)
======================================================================

T1 = (XX|YZ): the doubled centre sits inside ONE density -> reduces to the
3-centre ONE-ELECTRON integral, CLOSED (v4.81.0) at weight-1, gamma-FREE:
{exp, E1, ln}. Seed = (mu_X, alpha_Y, beta_Z), mu_X=Z_X(1/na+1/nc) [the X-pair],
alpha_Y=Z_Y/nb, beta_Z=Z_Z/nd (single-orbital rates on Y, Z) -- structurally the
hybrid+exchange rate pattern one level up, SAME (Z,n)-only dependence.

T2 = (XY|XZ): the genuine wall. NOT reducible to a single scalar seed the way
2-centre classes are. Momentum space gives

    y^2 = (c1 k^2 + 1)(c2 k^2 + 1),   c1 = s(1-s),  c2 = t(1-t)

an elliptic curve fibred over Feynman parameters (s,t) in [0,1]^2. c1, c2 are
UNIVERSAL functions of (s,t) -- the SAME continuous 2-parameter family for every
quartet in every molecule (this part of the "seed" is shared/reusable, like the
E1/ln special-function library itself). What varies per quartet is the ELLIPTIC
MOMENT's remaining data: the geometric distances D1=|X-Y|, D2=|X-Z| (fixed once
the shared-centre X and the molecule's geometry are fixed) and the orbital "mass"
content (za, zb, zc, zd) built from the same Z/n rates as the 2-centre engine.
So the DISCRETE parameter count compresses the same way (Z,n-only, not l,m) --
but the OBJECT indexed by that discrete tuple is not a lookup-able scalar the way
E1(mu R) is: it is an entire 2D-Feynman elliptic-period integral that still needs
NUMERICAL quadrature over (s,t) per instance (no closed form / elliptic-polylog
evaluator exists yet -- Paper 59's named frontier). This is stated as a NEGATIVE /
caveat below, not folded into the compression ratio as if it were free.

Run:  python debug/io_ladder_radial_seeds.py
"""
from __future__ import annotations

import sys
from fractions import Fraction
from itertools import product as iproduct
from pathlib import Path

import numpy as np
import sympy as sp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac import two_center_eri as E2  # noqa: E402
from geovac import noci_engine as GE  # noqa: E402


# =============================================================================
# PART 0 -- spot verification: rate arguments depend on n, NOT on l, m
# =============================================================================

def _exp_rates(expr):
    """Distinct decay rates b appearing as exp(-b*R) atoms."""
    out = set()
    for f in expr.atoms(sp.exp):
        out.add(sp.nsimplify(-sp.diff(f.args[0], E2.R_s)))
    return sorted(out)


def _e1_rates(expr):
    """Distinct rates mu appearing as E1(mu*R) atoms (E1 == expint(1, .))."""
    out = set()
    for f in expr.atoms(sp.Function):
        if f.func is sp.expint and f.args[0] == 1:
            out.add(sp.nsimplify(sp.diff(f.args[1], E2.R_s)))
    return sorted(out)


def part0_spot_verify() -> bool:
    print("=" * 78)
    print("PART 0 -- spot-verify: closed-form transcendental arguments are l,m-blind")
    print("=" * 78)
    ok = True

    # --- AA|BB: self-density (2,l,0)-(2,l,0) at l=0 vs l=1 should give the SAME
    # exp decay rate bB = Z(1/2+1/2) = Z, regardless of l.
    Z = Fraction(1)
    expr_s = E2.aabb_closed_form(Z, (1, 0, 0), (1, 0, 0), Z, (2, 0, 0), (2, 0, 0),
                                  simplify=False)
    expr_p = E2.aabb_closed_form(Z, (1, 0, 0), (1, 0, 0), Z, (2, 1, 0), (2, 1, 0),
                                  simplify=False)
    rates_s, rates_p = _exp_rates(expr_s), _exp_rates(expr_p)
    same = rates_s == rates_p
    ok &= same
    print(f"  AA|BB self-density (2s vs 2p, same n): exp-rate set")
    print(f"    2s|2s: {rates_s}")
    print(f"    2p|2p: {rates_p}")
    print(f"    identical -> {'OK (l,m-blind)' if same else 'FAIL'}")
    print(f"    predicted bB = Z*(1/2+1/2) = {Z}  -> "
          f"{'matches' if Z in rates_p else 'MISMATCH'}")

    # --- HYBRID: needs an l>0 pair to reach {E1, ln} at all (an s-type pair
    # dispatches to the elementary _hybrid_direct branch, per the file's own
    # class docstring). trio = (2p_m, 2p_m, 1s) on A for m=0 vs m=+1, lone 1s
    # on B. predicted mu = Z_A*(1/2+1/2+1/1) = 2*Z_A; ad = Z_B/1.
    ZA, ZB = Fraction(3), Fraction(1)
    h_0 = E2.hybrid_closed_form(ZA, (2, 1, 0), (2, 1, 0), (1, 0, 0), ZB, (1, 0, 0))
    h_1 = E2.hybrid_closed_form(ZA, (2, 1, 1), (2, 1, 1), (1, 0, 0), ZB, (1, 0, 0))
    r0, r1 = _e1_rates(h_0), _e1_rates(h_1)
    same2 = r0 == r1
    ok &= same2
    mu_pred = ZA * (Fraction(1, 2) + Fraction(1, 2) + Fraction(1, 1))
    print(f"\n  HYBRID trio (2p_m,2p_m,1s) on A, m=0 vs m=+1, lone 1s on B: "
          f"E1-rate set")
    print(f"    m=0 : {r0}")
    print(f"    m=+1: {r1}")
    print(f"    identical -> {'OK (l,m-blind)' if same2 else 'FAIL'}")
    print(f"    predicted generating rates mu={mu_pred}, ad={ZB}, mu+ad={mu_pred + ZB}"
          f"  -> {'mu+ad present' if mu_pred + ZB in r0 else 'CHECK'} "
          f"(internal expansion carries {len(r0)} derived E1 rates per instance; "
          f"the (mu,ad) PAIR is what varies across instances -- see driver seed key)")

    print(f"\n  PART 0 VERDICT: {'PASS' if ok else 'FAIL'} -- radial seed content is "
          f"(Z, n)-indexed only, confirmed against the live closed-form code.")
    return ok


# =============================================================================
# PART A -- two-centre engine: naive entry count vs distinct rate-seed count
# =============================================================================

def hydrogenic_catalog(n_max: int, center: int):
    out = []
    for n in range(1, n_max + 1):
        for l in range(0, n):
            for m in range(-l, l + 1):
                out.append((center, n, l, m))
    return out


def classify_and_seed(orb, ZA: Fraction, ZB: Fraction, i, j, k, l):
    """Return (class_name, seed_tuple) for one unique quartet (i,j,k,l), or
    None if the quartet is a same-atom (out-of-scope) quartet."""
    ci, cj, ck, cl = orb[i][0], orb[j][0], orb[k][0], orb[l][0]
    ni, nj, nk, nl = orb[i][1], orb[j][1], orb[k][1], orb[l][1]
    centers = (ci, cj, ck, cl)
    ns = (ni, nj, nk, nl)
    Zc = {0: ZA, 1: ZB}
    cset = set(centers)
    if len(cset) == 1:
        return None  # same-atom quartet: a different engine (Slater R_k), OOS

    count0 = centers.count(0)
    if count0 in (1, 3):
        trio_center = 0 if count0 == 3 else 1
        lone_center = 1 - trio_center
        trio_ns = [ns[t] for t in range(4) if centers[t] == trio_center]
        lone_n = [ns[t] for t in range(4) if centers[t] == lone_center][0]
        mu = Zc[trio_center] * sum(Fraction(1, n) for n in trio_ns)
        ad = Zc[lone_center] * Fraction(1, lone_n)
        seed = ("hybrid", trio_center, mu, ad)
        return "hybrid", seed

    # count0 == 2: either AA|BB (bra, ket each single-centre) or EXCHANGE (2+2
    # split within both bra and ket)
    bra_c, ket_c = (ci, cj), (ck, cl)
    if bra_c[0] == bra_c[1] and ket_c[0] == ket_c[1]:
        bra_center, ket_center = bra_c[0], ket_c[0]
        bA = Zc[bra_center] * (Fraction(1, ns[0]) + Fraction(1, ns[1]))
        bB = Zc[ket_center] * (Fraction(1, ns[2]) + Fraction(1, ns[3]))
        seed = ("aabb", tuple(sorted((bra_center, ket_center))),
                tuple(sorted((bA, bB))))
        return "aabb", seed
    else:
        # exchange: bra has one A-orbital + one B-orbital; same for ket.
        bra_A_idx = 0 if ci == 0 else 1
        bra_B_idx = 1 - bra_A_idx
        ket_A_idx = 2 if ck == 0 else 3
        ket_B_idx = 5 - ket_A_idx  # 2+3=5
        alpha_bra = ZA * Fraction(1, ns[bra_A_idx])
        beta_bra = ZB * Fraction(1, ns[bra_B_idx])
        alpha_ket = ZA * Fraction(1, ns[ket_A_idx])
        beta_ket = ZB * Fraction(1, ns[ket_B_idx])
        # (ab|cd) = (cd|ab): bra<->ket swap is a physical symmetry, so
        # canonicalize by sorting the two (alpha,beta) pairs -- otherwise the
        # raw (unreduced) loop double-counts each physical seed under both
        # index orderings.
        seed = ("exchange", tuple(sorted([(alpha_bra, beta_bra),
                                           (alpha_ket, beta_ket)])))
        return "exchange", seed


def part_a_two_center(ZA_val: int, ZB_val: int, n_max: int):
    """RAW (unreduced) quartet count -- the naive M^4 dense-tensor convention
    already used elsewhere in this repo (p58 census 'dense', Poly-0/2 '2401'),
    so the numbers here are directly comparable to established figures. No
    8-fold bra/ket symmetry reduction is applied on the ENTRY side; the SEED
    side counts distinct rate-tuples regardless of how many raw entries share
    one, which is exactly the compression this rung measures."""
    ZA, ZB = Fraction(ZA_val), Fraction(ZB_val)
    orb = hydrogenic_catalog(n_max, 0) + hydrogenic_catalog(n_max, 1)
    M = len(orb)
    print(f"\n  system Z_A={ZA_val}, Z_B={ZB_val}, n_max={n_max}: "
          f"M = {M} orbitals ({M // 2} per centre), M^4 = {M**4} (dense/naive)")

    entries = {"aabb": 0, "hybrid": 0, "exchange": 0}
    seeds = {"aabb": set(), "hybrid": set(), "exchange": set()}
    oos_same_atom = 0

    for i in range(M):
        for j in range(M):
            for k in range(M):
                for l in range(M):
                    res = classify_and_seed(orb, ZA, ZB, i, j, k, l)
                    if res is None:
                        oos_same_atom += 1
                        continue
                    cls, seed = res
                    entries[cls] += 1
                    seeds[cls].add(seed)

    print(f"  (same-atom quartets excluded, out of this engine's scope: {oos_same_atom})")
    rows = []
    tot_entries = tot_seeds = 0
    for cls in ("aabb", "hybrid", "exchange"):
        ne, ns_ = entries[cls], len(seeds[cls])
        ratio = ne / ns_ if ns_ else float("nan")
        rows.append((cls, ne, ns_, ratio))
        tot_entries += ne
        tot_seeds += ns_
        print(f"    {cls:<9} entries = {ne:7d}   distinct seeds = {ns_:5d}   "
              f"ratio = {ratio:8.2f}x")
    tot_ratio = tot_entries / tot_seeds if tot_seeds else float("nan")
    print(f"    {'TOTAL':<9} entries = {tot_entries:7d}   distinct seeds = {tot_seeds:5d}   "
          f"ratio = {tot_ratio:8.2f}x")
    return {"Z_A": ZA_val, "Z_B": ZB_val, "n_max": n_max, "M": M,
            "rows": rows, "total_entries": tot_entries, "total_seeds": tot_seeds,
            "total_ratio": tot_ratio}


# =============================================================================
# PART B -- three-centre engine: T1 (closed) and T2 (elliptic, open) seed counts
# =============================================================================

def _shapes():
    out = {}
    for kind, (l, n_r) in (("1s", (0, 1)), ("2s", (0, 2)), ("2p", (1, 2))):
        a, d, q = GE.fit_sto_shape(l, n_r, n_gauss=8)
        out[kind] = (a, d)
    return out


def water_system():
    """Same geometry/basis as debug/poly0_three_center_burden.py, reused verbatim
    (not re-derived) so the 420/2401 3-centre/total figures are directly
    comparable."""
    BOHR_OH, ANG_HOH = 1.809, 104.5
    th = np.radians(ANG_HOH / 2.0)
    O = np.array([0.0, 0.0, 0.0])
    H1 = np.array([BOHR_OH * np.sin(th), 0.0, BOHR_OH * np.cos(th)])
    H2 = np.array([-BOHR_OH * np.sin(th), 0.0, BOHR_OH * np.cos(th)])
    pos = [O, H1, H2]
    # (center_idx, shell_label, zeta)
    orb = [(0, "1s", 7.66), (0, "2s", 2.25),
           (0, "2p", 2.23), (0, "2p", 2.23), (0, "2p", 2.23),
           (1, "1s", 1.0), (2, "1s", 1.0)]
    return pos, orb


def part_b_three_center():
    """RAW (unreduced) M^4 quartet count -- matches Poly-0/Poly-2's own
    convention exactly (2401 = 7^4 total, 140 T1 + 280 T2 = 420 three-centre),
    so the T1/T2 split is a direct, citable cross-check against established
    repo numbers, not a re-derivation on a different convention."""
    pos, orb = water_system()
    M = len(orb)
    print(f"\n  H2O (bent, C2v), STO-shape minimal basis: M = {M} orbitals, "
          f"M^4 = {M**4} (matches poly0/poly2's 2401 total tensor entries)")

    T1_entries, T2_entries = 0, 0
    T1_seeds, T2_seeds = set(), set()
    other = 0

    def dist(a, b):
        return round(float(np.linalg.norm(pos[a] - pos[b])), 6)

    for i in range(M):
        for j in range(M):
            for k in range(M):
                for l in range(M):
                    cen = (orb[i][0], orb[j][0], orb[k][0], orb[l][0])
                    if len(set(cen)) != 3:
                        other += 1
                        continue
                    bra, ket = (cen[0], cen[1]), (cen[2], cen[3])
                    zeta = (orb[i][2], orb[j][2], orb[k][2], orb[l][2])
                    if bra[0] == bra[1] or ket[0] == ket[1]:
                        # T1: (XX|YZ) -- one density one-centre, other two-centre
                        T1_entries += 1
                        if bra[0] == bra[1]:
                            X, (za, zc) = bra[0], (zeta[0], zeta[1])
                            Y, Z = ket
                            zb, zd = zeta[2], zeta[3]
                        else:
                            X, (za, zc) = ket[0], (zeta[2], zeta[3])
                            Y, Z = bra
                            zb, zd = zeta[0], zeta[1]
                        mu_X = tuple(sorted((round(za, 6), round(zc, 6))))
                        # (ab|cd)=(cd|ab): the (Y,Z) two-centre density's pair
                        # is symmetric, so canonicalize by sorting the two
                        # (centre, zeta, distance-to-X) sides TOGETHER --
                        # otherwise a raw (unreduced) loop can visit the same
                        # physical instance via (Y,Z) and (Z,Y) index order
                        # and record it as two different seeds.
                        side_Y = (Y, round(zb, 6), dist(X, Y))
                        side_Z = (Z, round(zd, 6), dist(X, Z))
                        seed = ("T1", X, mu_X, tuple(sorted((side_Y, side_Z))))
                        T1_seeds.add(seed)
                    else:
                        # T2: (XY|XZ) -- shared centre split across bra and ket
                        T2_entries += 1
                        # shared centre = the one appearing in both bra and ket
                        X = (set(bra) & set(ket)).pop()
                        Y = [c for c in bra if c != X][0]
                        Z = [c for c in ket if c != X][0]
                        a_idx = 0 if bra[0] == X else 1
                        za = zeta[a_idx]
                        zb = zeta[1 - a_idx]
                        c_idx = 2 if ket[0] == X else 3
                        zc = zeta[c_idx]
                        zd = zeta[5 - c_idx]
                        # (ab|cd)=(cd|ab): bra<->ket swap exchanges the
                        # (Y,za,zb,D(X,Y)) side with the (Z,zc,zd,D(X,Z)) side.
                        # Canonicalize by sorting the two sides together.
                        side_bra = (Y, round(za, 6), round(zb, 6), dist(X, Y))
                        side_ket = (Z, round(zc, 6), round(zd, 6), dist(X, Z))
                        seed = ("T2", X, tuple(sorted((side_bra, side_ket))))
                        T2_seeds.add(seed)

    print(f"  (1-centre / 2-centre-only entries excluded, out of 3-centre scope: {other})")
    print(f"  cross-check vs Poly-0/Poly-2 (140 T1 + 280 T2 = 420 total, water): "
          f"{'OK' if (T1_entries, T2_entries) == (140, 280) else 'MISMATCH'} "
          f"(got T1={T1_entries}, T2={T2_entries})")

    r1 = T1_entries / len(T1_seeds) if T1_seeds else float("nan")
    r2 = T2_entries / len(T2_seeds) if T2_seeds else float("nan")
    print(f"\n  T1 (XX|YZ, CLOSED weight-1 gamma-free, v4.81.0):")
    print(f"    entries = {T1_entries:4d}   distinct (mu_X,alpha_Y,beta_Z,geom) "
          f"seeds = {len(T1_seeds):3d}   ratio = {r1:6.2f}x")
    print(f"  T2 (XY|XZ, elliptic Bessel moment, OPEN -- Paper 59):")
    print(f"    entries = {T2_entries:4d}   distinct (za,zb,zc,zd,geom) "
          f"'mass+geometry' instances = {len(T2_seeds):3d}   ratio = {r2:6.2f}x"
          f"  [NEGATIVE CAVEAT: see below]")

    return {"M": M, "T1_entries": T1_entries, "T1_seeds": len(T1_seeds), "T1_ratio": r1,
            "T2_entries": T2_entries, "T2_seeds": len(T2_seeds), "T2_ratio": r2}


# =============================================================================

def main() -> None:
    ok0 = part0_spot_verify()

    print("\n" + "=" * 78)
    print("PART A -- TWO-CENTRE ENGINE (Li/H, Z_A=3 Z_B=1, matches p58 census system)")
    print("=" * 78)
    res_n2 = part_a_two_center(3, 1, 2)
    res_n3 = part_a_two_center(3, 1, 3)

    print("\n" + "=" * 78)
    print("PART B -- THREE-CENTRE ENGINE (H2O, matches poly0/poly2 basis+geometry)")
    print("=" * 78)
    res_b = part_b_three_center()

    print("\n" + "=" * 78)
    print("SUMMARY -- decision gate: distinct-seed count << entry count, sub-M^2 growth?")
    print("=" * 78)
    print(f"  PART 0 (l,m-blindness of the rate arguments): {'PASS' if ok0 else 'FAIL'}")
    print(f"  2c n_max=2: {res_n2['total_entries']} entries -> "
          f"{res_n2['total_seeds']} seeds ({res_n2['total_ratio']:.1f}x)")
    print(f"  2c n_max=3: {res_n3['total_entries']} entries -> "
          f"{res_n3['total_seeds']} seeds ({res_n3['total_ratio']:.1f}x)")
    m2, m3 = res_n2["M"], res_n3["M"]
    e2, e3 = res_n2["total_entries"], res_n3["total_entries"]
    s2, s3 = res_n2["total_seeds"], res_n3["total_seeds"]
    if m3 > m2 and e2 > 0 and s2 > 0:
        entry_exp = np.log(e3 / e2) / np.log(m3 / m2)
        seed_exp = np.log(s3 / s2) / np.log(m3 / m2)
        print(f"  growth exponents (M: {m2}->{m3}): entries ~ M^{entry_exp:.2f}, "
              f"seeds ~ M^{seed_exp:.2f}")
    print(f"  3c T1 (closed, weight-1 gamma-free): {res_b['T1_entries']} entries -> "
          f"{res_b['T1_seeds']} seeds ({res_b['T1_ratio']:.1f}x)")
    print(f"  3c T2 (elliptic, OPEN): {res_b['T2_entries']} entries -> "
          f"{res_b['T2_seeds']} distinct (mass,geometry) instances "
          f"({res_b['T2_ratio']:.1f}x) -- COUNT compresses, but each instance is "
          f"a 2D-Feynman elliptic-period integral, NOT a lookup-able scalar: no "
          f"closed-form/elliptic-polylog evaluator exists yet (Paper 59 frontier).")


if __name__ == "__main__":
    main()
