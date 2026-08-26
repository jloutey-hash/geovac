"""PSLQ fit runner: V (T2, corner-subtracted) against the weight-1/weight-2
Eisenstein/CM-Gamma ring, with a decoy (meaningless, same-magnitude constant)
control. See routeC_pslq_fit.py for the basis construction.
"""
import sys
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_pslq_fit import build_basis, run_pslq


def main():
    V_str = sys.argv[1]
    dps_list = [int(a) for a in sys.argv[2:]] if len(sys.argv) > 2 else [25, 35]
    maxcoeff = 10 ** 6

    for dps in dps_list:
        print(f"=== dps={dps} ===")
        mp.mp.dps = dps + 15
        V = mp.mpf(V_str)
        w1, w2, chk = build_basis(dps + 10)
        print(f"  varpi vs K(1/2) consistency check: {mp.nstr(chk, 4)}")

        mp.mp.dps = dps
        # COMBINED weight-1 U weight-2 basis, RANK-DEFICIENCY FIXED by
        # dropping the redundant generator 'pi' from weight-1 (per the
        # coordinator's Task-2 direction: use an independent generating set,
        # not a weight truncation). The classical Legendre relation at the
        # lemniscatic point, 4*varpi*E12 - 2*varpi^2 = pi, makes {pi} U
        # {varpi^2, varpi*E12} linearly dependent; since pi is RECOVERABLE
        # from the weight-2 products, dropping it (keeping {1,varpi,E12,K2}
        # for weight-1 + all 10 weight-2 products = 14 generators) removes
        # the dependency without truncating the ring V is conjectured to
        # live in.
        combined = {}
        for k, v in w1.items():
            if k in ('1', 'pi'):
                continue
            combined[k] = v
        for k, v in w2.items():
            combined[k] = v

        rel = run_pslq(V, combined, dps, maxcoeff=maxcoeff, label=f"REAL V, dps={dps}")

        # DECOY: a same-magnitude, structurally meaningless constant.
        # sqrt(2)*ln(3)/7 + 1/e  is order ~0.3-0.5, no algebraic/transcendental
        # relation to the weight-1/2 CM-Gamma ring by construction.
        decoy = mp.sqrt(2) * mp.log(3) / 7 + 1 / mp.e
        print(f"  decoy value = {mp.nstr(decoy, 15)}")
        rel_decoy = run_pslq(decoy, combined, dps, maxcoeff=maxcoeff, label=f"DECOY, dps={dps}")


if __name__ == '__main__':
    main()
