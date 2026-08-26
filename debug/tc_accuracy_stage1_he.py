"""Stage-1 capstone diagnostic (He): does the native NON-HERMITIAN TC track the
variational R12-CI accuracy anchor, or is its accuracy a fitted-gamma artifact?

Head-to-head at MATCHED (ns, k), scanning the geminal width gamma:
  - E_plain      : s-only Coulomb FCI (the ~24.7 mHa s-limit plateau)
  - E_TC         : native non-Hermitian TC  (D + K), scipy.linalg.eig ground  [CHEAP-qubit]
  - E_sym        : symmetrized (H+H^dag)/2                                     [rejected route]
  - E_R12CI      : variational explicitly-correlated CI (bounded >= exact)     [ACCURATE anchor]
exact He = -2.903724 ; s-limit -2.879029.

GATE:
  GO       - E_TC lands within ~1-2 mHa of exact AND tracks R12-CI (both move together with gamma),
             so the cheap non-Herm operator inherits the variational geminal accuracy.
  CAUTION  - E_TC's closeness is a lone gamma-crossing (monotone through exact) while R12-CI has a
             genuine variational optimum -> accuracy is a fitted-gamma artifact, not inherited.
"""
import os, sys, json
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import xtc_sym_accuracy_validate as V          # non-Herm TC (prep/plain/tc_energies), He 2-body
import ctf12_r12ci_he as R                     # variational R12-CI he_energy

EXACT = -2.903724
NS, K = 3, 1.9
GAMMAS = [0.5, 0.7, 0.9, 1.1, 1.3, 1.6, 2.0, 2.5, 3.0]
NG, NX = 600, 128


def main():
    P = V.prep(NS, K, Z=2, ne=2, Ng=NG)
    E_plain = V.plain_energy(P, nx=NX)
    rows = []
    print("He  ns=%d k=%.2f   exact=%.6f  (s-limit -2.879029)" % (NS, K, EXACT))
    print("plain(s-only) FCI = %.6f  (%+.2f mHa vs exact)\n" % (E_plain, (E_plain - EXACT) * 1e3))
    print("%5s | %11s %8s | %11s %8s | %11s %8s | %8s" %
          ("gamma", "E_TC", "dTC", "E_sym", "dsym", "E_R12CI", "dR12", "TC-R12"))
    print("-" * 86)
    for g in GAMMAS:
        E_TC, im, E_sym = V.tc_energies(P, g, with_L3=False, nx=NX)
        E_r12 = R.he_energy(NS, K, g, n_gem=1, gem_refs=[1], Ng=NG, nx=NX)[0]
        dTC, dsym, dR12 = (E_TC - EXACT) * 1e3, (E_sym - EXACT) * 1e3, (E_r12 - EXACT) * 1e3
        rows.append(dict(gamma=g, E_plain=E_plain, E_TC=E_TC, im_TC=im, E_sym=E_sym,
                         E_R12CI=E_r12, dTC_mHa=dTC, dsym_mHa=dsym, dR12_mHa=dR12,
                         TC_minus_R12_mHa=(E_TC - E_r12) * 1e3))
        print("%5.2f | %11.6f %+8.2f | %11.6f %+8.2f | %11.6f %+8.2f | %+8.2f" %
              (g, E_TC, dTC, E_sym, dsym, E_r12, dR12, (E_TC - E_r12) * 1e3))

    # summaries
    dTC = np.array([r['dTC_mHa'] for r in rows])
    dR12 = np.array([r['dR12_mHa'] for r in rows])
    gs = np.array(GAMMAS)
    print("\n--- read ---")
    i_tc = int(np.argmin(np.abs(dTC)));  i_r12 = int(np.argmin(np.abs(dR12)))
    print("best |E_TC - exact|   = %.2f mHa @ gamma=%.2f  (side: %s)" %
          (abs(dTC[i_tc]), gs[i_tc], "above" if dTC[i_tc] > 0 else "BELOW (overshoot)"))
    print("best |E_R12CI - exact|= %.2f mHa @ gamma=%.2f  (variational anchor; should be >=0)" %
          (abs(dR12[i_r12]), gs[i_r12]))
    r12_all_above = bool(np.all(dR12 > -0.5))
    tc_crosses = bool(np.any(dTC > 0) and np.any(dTC < 0))
    print("R12-CI variational (all >= exact within 0.5 mHa)? %s" % r12_all_above)
    print("E_TC crosses exact over the gamma range (monotone-through)? %s" % tc_crosses)
    # tracking: correlation of the two curves
    corr = float(np.corrcoef(dTC, dR12)[0, 1])
    print("corr(E_TC(gamma), E_R12CI(gamma)) = %.3f" % corr)

    os.makedirs(os.path.join(_HERE, 'data'), exist_ok=True)
    outp = os.path.join(_HERE, 'data', 'tc_accuracy_stage1_he.json')
    json.dump(dict(ns=NS, k=K, exact=EXACT, E_plain=E_plain, rows=rows,
                   best_TC_mHa=float(abs(dTC[i_tc])), best_TC_gamma=float(gs[i_tc]),
                   best_R12_mHa=float(abs(dR12[i_r12])),
                   r12_variational=r12_all_above, tc_crosses_exact=tc_crosses,
                   corr_TC_R12=corr),
              open(outp, 'w'), indent=2, default=float)
    print("\nwrote", outp)


if __name__ == '__main__':
    main()
