"""Test the split-shell mechanism supplied by the external source.

Claim: the ground state is starved because BOTH electrons sit at n=1 and are
therefore forced to the identical exponent Q_nu/n_j; excited states (1s2s) get
two different exponents inside one configuration for free, so they should
converge far better under the SAME locked scale.

Prediction: locked-scale gap(1s2s) << locked-scale gap(1s^2).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import debug.p60_engine as E
import geovac.sturmian_secular as S

GROUND = -2.903724377      # He 1^1S exact non-relativistic
EXC1 = -2.145974046        # He 1s2s 1^1S... (2 1S) exact non-relativistic
lmax = int(sys.argv[1])
print("nmax    K   E(1S ground)   gap_mHa  |  E(2 1S)      gap_mHa")
for nmax in [int(x) for x in sys.argv[2:]]:
    E.set_grid(max(80.0, 5.0 * nmax * nmax), 24000, "grade", 2.0)
    cfgs = S.build_configs(E.family(nmax, lmax))
    M = S.build_M(cfgs, Z=2.0)
    p = np.sort(np.linalg.eigvalsh(M))[::-1]
    e0, e1 = -p[0] ** 2 / 2, -p[1] ** 2 / 2
    print("%4d %4d   %.7f  %8.3f  |  %.7f  %8.3f"
          % (nmax, len(cfgs), e0, (e0 - GROUND) * 1000, e1, (e1 - EXC1) * 1000))
