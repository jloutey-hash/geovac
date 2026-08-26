import numpy as np
from geovac.molecular_spec import beh2_spec
from geovac.composed_qubit import build_composed_hamiltonian

spec = beh2_spec(max_n=2)
res = build_composed_hamiltonian(spec, pk_in_hamiltonian=False, verbose=False)
h1 = res['h1']; M = res['M']
ev = np.linalg.eigvalsh(0.5*(h1+h1.T))
print("COMPOSED (block-diagonal) one-body h1 spectrum at equilibrium, M=%d" % M)
np.set_printoptions(precision=4, suppress=True, linewidth=120)
print("eigs:", ev)
gaps = np.diff(ev)
print("min adjacent gap: %.5f Ha at level %d->%d (eigs %.4f, %.4f)" % (
      gaps.min(), gaps.argmin(), gaps.argmin()+1, ev[gaps.argmin()], ev[gaps.argmin()+1]))
# occupied count = 6 electrons / 2 = 3 lowest spatial (frozen-core + valence); HOMO-LUMO-ish
occ = 3
print("HOMO(#%d)=%.4f  LUMO(#%d)=%.4f  gap=%.4f Ha (%.3f eV)" % (
      occ-1, ev[occ-1], occ, ev[occ], ev[occ]-ev[occ-1], (ev[occ]-ev[occ-1])*27.2114))
