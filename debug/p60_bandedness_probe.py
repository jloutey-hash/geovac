"""Is the overlap BANDED in n on one centre, and does translation destroy it?

Claim under test: bandedness in n is a PER-CENTRE structure (Laguerre three-term
recurrence), so it should die under translation to a second centre for the same
reason l-selection does -- exactly Paper 58's abelian-residue mechanism, but on
the radial label instead of the angular one.

Measured as the decay of |S_{n,n'}| with |n - n'| within a fixed l, separately
for the same-centre block and the cross-centre block.
"""
import os, sys
from fractions import Fraction
import numpy as np
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from debug.p60_sigma_cond_discriminate import s_element

R, E = 2.0, -1.1026342144
k = float(np.sqrt(-2.0 * E)); kf = Fraction(k).limit_denominator(10 ** 12)
Rs = sp.nsimplify(sp.Float(R, 30))
NM, L = 6, 0

A = [("A", kf * n, n, L) for n in range(1, NM + 1)]
B = [("B", kf * n, n, L) for n in range(1, NM + 1)]
same = np.array([[float(sp.N(s_element(a, b, Rs), 30)) for b in A] for a in A])
cross = np.array([[float(sp.N(s_element(a, b, Rs), 30)) for b in B] for a in A])

print("l=%d, shared-scale Coulomb Sturmians, n=1..%d, R=%.1f" % (L, NM, R))
print("magnitude of S by |n - n'| (max over the diagonal band):")
print("  |n-n'|   same-centre      cross-centre")
for d in range(NM):
    s = max(abs(same[i, i + d]) for i in range(NM - d))
    c = max(abs(cross[i, i + d]) for i in range(NM - d))
    print("     %d     %.3e       %.3e" % (d, s, c))
print()
sb = max(abs(same[i, j]) for i in range(NM) for j in range(NM) if abs(i - j) > 2)
cb = max(abs(cross[i, j]) for i in range(NM) for j in range(NM) if abs(i - j) > 2)
print("max |S| OUTSIDE a bandwidth of 2:")
print("   same-centre  = %.3e   <- banded if ~0" % sb)
print("   cross-centre = %.3e   <- dense if O(1)" % cb)
