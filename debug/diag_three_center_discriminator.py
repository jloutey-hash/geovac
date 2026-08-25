"""DIAGNOSTIC addendum (2026-08-25): separate the two invariants that the Paper 32 remark
slightly conflates --
  (A) commutant dim = 1  (IRREDUCIBLE, the CATEGORY jump 2 reducible / 3 irreducible), and
  (B) nested-commutator norm ||[[P_A,P_B],P_C]||  (the STRENGTH of three-body coupling).

Question: is (A) a bonding phenomenon, or generic for any 3 coupled projections?  Sweep the
apex-H distance on the linear enriched-sigma BeH2 {1s,2s,2p0} and watch both.  Hypothesis:
(A) stays 1 until the residual overlap drops below the numerical tol (a knife-edge / generic
property of >=3 coupled projections), while (B) decays smoothly to 0 -- the physical signal.
"""
from __future__ import annotations
import sys
import numpy as np
sys.path.insert(0, "debug")
from diag_three_center_generality import linear_sigma_G, projectors_from_G
from beh2_algebra_classification import commutant


def nested_max(Ps):
    def nst(A, B, C):
        c = A @ B - B @ A
        return float(np.linalg.norm(c @ C - C @ c, 2))
    return max(nst(Ps[0], Ps[1], Ps[2]), nst(Ps[1], Ps[2], Ps[0]), nst(Ps[0], Ps[2], Ps[1]))


def pairmax(Ps):
    return max(float(np.linalg.norm(Ps[i] @ Ps[j] - Ps[j] @ Ps[i], 2))
               for (i, j) in [(0, 1), (1, 2), (0, 2)])


STATES = [(1, 0), (2, 0), (2, 1)]   # enriched sigma {1s,2s,2p0}

print("Linear BeH2 enriched-sigma {1s,2s,2p0}, Be Z=2 / H Z=1 -- distance sweep")
print(" d(BeH)   overlap~   pair||[P,P]||   3-way||[[P,P],P]||   dim(A')  -> ")
for d in [2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 25.0, 30.0]:
    G, sz = linear_sigma_G(d, 2, 1, STATES)
    Ps, evmin = projectors_from_G(G, sz)
    dAp, _ = commutant(Ps)
    pm, tw = pairmax(Ps), nested_max(Ps)
    verdict = "IRREDUCIBLE M9" if dAp == 1 else ("reducible dimA'=%d" % dAp)
    print("  %5.1f   %8.1e   %10.4f     %12.4f       %3d    %s"
          % (d, np.exp(-d), pm, tw, dAp, verdict))

print("\nReading: dim(A')=1 (category = 'irreducible') is a KNIFE-EDGE -- it holds for any")
print("nonzero coupling and only flips to reducible once residual overlap falls below the")
print("numerical tol.  The BONDING-specific signal is the nested-commutator MAGNITUDE (col 4),")
print("which decays smoothly.  Two different invariants; report both.")
