import warnings
warnings.simplefilter("ignore")
from geovac.lattice_index import LatticeIndex
from geovac.qubit_encoding import JordanWignerEncoder
S = 2.0 ** 20
for n in (5,):
    li = LatticeIndex(n_electrons=2, max_n=n, nuclear_charge=2, vee_method='slater_full', h1_method='exact')
    li._H1_spatial = li._H1_spatial * S
    li._eri = {k: v * S for k, v in li._eri.items()}
    enc = JordanWignerEncoder(li); enc._v_nn *= S
    q = enc.build_qubit_operator()
    nonid = sum(1 for k, v in q.terms.items() if k and abs(v) / S > 1e-11)
    print(f"He n={n}: scaled-build non-identity terms above 1e-11: {nonid}", flush=True)
