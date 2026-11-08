"""
Paper-8 hinge test: is the genuine shared-p0 Coulomb-Sturmian two-center overlap
BLOCK-DIAGONAL in n (Paper 8's orthogonal SO(4) D-matrix, => R-independent
eigenvalues => the corpse), or does it MIX n (=> R-dependent => binding accessible)?

Paper 8 §Eigenfunctions-as-D-matrix asserts the second-pole placement is an SO(4)
rotation g(gamma), so the cross-center block = D^(n)(gamma) is block-diagonal in n
(line 324-330), and "the standard Shibuya-Wulfman formula is exact only within a
single n-manifold" (line 337). The Explorer + the hydrogenic commutator probe both
suggest the genuine overlap mixes n. This tests it on the ACTUAL Sturmian object.

Key identity: an L2-normalized Coulomb-Sturmian with common scale k equals the
L2-normalized hydrogenic radial R_nl_modern evaluated at Z = n*k (both give radial
scale 2k, rho = 2kr). So the shared-p0 Sturmian two-center overlap is
    W[n1 l1, n2 l2](m,R) = overlap_two_center(n1*k, n1, l1, n2*k, n2, l2, m, R),
reusing the VALIDATED Topos-3 machinery (no new radial code).

DECISIVE, normalization-robust test: are the CROSS-n entries (n1 != n2, same m,l
allowed to differ) zero (block-diagonal D-matrix, Paper 8) or nonzero (mixing)?
An SO(4) rotation gives EXACTLY zero cross-n. Any nonzero cross-n refutes the
block-diagonal model for the genuine Sturmian overlap.
"""
from __future__ import annotations
import importlib.util
import json
import mpmath as mp

mp.mp.dps = 10

spec = importlib.util.spec_from_file_location(
    "topos3", "debug/compute_topos3_two_center_meet.py")
topos3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(topos3)
overlap_two_center = topos3.overlap_two_center


def sturmian_overlap(k, n1, l1, n2, l2, m, R):
    """Shared-scale-k Coulomb-Sturmian two-center overlap <chi_{n1l1}(A)|chi_{n2l2}(B)>."""
    return float(overlap_two_center(n1 * k, n1, l1, n2 * k, n2, l2, m, mp.mpf(R)))


def m0_states(n_max):
    return [(n, l) for n in range(1, n_max + 1) for l in range(n) if l == 0 or True and 0 <= 0]


def probe(k, R, n_max=3, m=0):
    # m=0 states across n: (n,l) with l>=m, l<n
    states = [(n, l) for n in range(1, n_max + 1) for l in range(m, n)]
    W = [[sturmian_overlap(k, na, la, nb, lb, m, R) for (nb, lb) in states]
         for (na, la) in states]
    # classify entries: within-n (na==nb) vs cross-n (na!=nb)
    cross_n = [(states[i], states[j], W[i][j])
               for i in range(len(states)) for j in range(len(states))
               if states[i][0] != states[j][0]]
    max_cross = max((abs(v) for *_, v in cross_n), default=0.0)
    max_within = max(abs(W[i][j]) for i in range(len(states))
                     for j in range(len(states)) if states[i][0] == states[j][0])
    return states, W, max_cross, max_within, cross_n


if __name__ == "__main__":
    out = {}
    print("=" * 90)
    print("Paper-8 hinge: does the shared-p0 Coulomb-Sturmian two-center overlap mix n?")
    print("(SO(4) rotation / block-diagonal D-matrix => cross-n == 0 exactly)")
    print("=" * 90)
    for k, R in [(1.0, 3.015), (1.0, 5.0), (1.5, 3.015)]:
        states, W, mx_cross, mx_within, cross = probe(k, R, n_max=3, m=0)
        key = f"k={k},R={R}"
        print(f"\n--- {key}, m=0, states {states} ---")
        hdr = "        " + " ".join(f"({n}{l})".rjust(8) for n, l in states)
        print(hdr)
        for i, (na, la) in enumerate(states):
            print(f"  ({na}{la}) " + " ".join(f"{W[i][j]:+8.4f}" for j in range(len(states))))
        print(f"  max |within-n| = {mx_within:.4f}   max |CROSS-n| = {mx_cross:.4f}"
              f"   ratio cross/within = {mx_cross/mx_within:.3f}")
        # biggest cross-n couplings
        cross_sorted = sorted(cross, key=lambda t: -abs(t[2]))[:4]
        for a, b, v in cross_sorted:
            print(f"     cross-n  <{a[0]}{a[1]} | {b[0]}{b[1]}> = {v:+.4f}")
        out[key] = {'states': states, 'W': W, 'max_cross_n': mx_cross,
                    'max_within_n': mx_within, 'ratio': mx_cross / mx_within}

    json.dump(out, open('debug/data/sturmian_hinge.json', 'w'), indent=2)
    print("\n[saved] debug/data/sturmian_hinge.json")
    print("\n" + "=" * 90)
    print("VERDICT: cross-n entries NONZERO and comparable to within-n => the genuine")
    print("Sturmian second-pole operator is NOT block-diagonal in n, NOT an SO(4)")
    print("rotation => Paper 8's D-matrix model (and its R-independence corpse) does")
    print("not describe the genuine two-center Sturmian overlap.  cross-n ~ 0 => corpse")
    print("stands and the Avery route buys only k(R)/SCF, no fixed-graph escape.")
