r"""Map the ANALYTIC eigenvectors of en_factors (basis [chi_I, xi, W1 psi, W2 psi, (Fcc psi, Fc1 psi)]) onto
the stage-0 LM column layout (debug/en_factors_stage0_lm.py: Phi_I, F Phi_I, class columns), so that
plain VMC re-samples the variational energy of the FIXED analytic vector (unbiased; gate G3 cross-check).
chi_I = Phi_I / 2, psi = sum d_I Phi_I / 2;  P^{key}_{cls,gamma} = G Psi_LM = 2 G psi;
W2 psi = P_cc0 / 2,  W1 psi = (P_cv0 + 2 P_cc0) / 6,  Fcc psi = P_cc1.5 / 2,  Fc1 psi = (P_cv1.5 + 2 P_cc1.5) / 2.
Usage: python debug/en_factors_vmc_vecs.py 1e-2   -> debug/data/en_factors_vmcvecs_1e-02.npz (keys ph1, ph2)
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
KEYS = (0.6, 1.0, -1.0, -2.0)
GAMMAS = (0.0, 1.5)
CLASSES = ("cc", "cv", "vv")


def main():
    t = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-2
    z1 = np.load(os.path.join(DATA, f"en_factors_{t:.0e}.npz"), allow_pickle=True)
    d = z1['d']; N = len(d)
    gold = json.load(open(os.path.join(os.path.dirname(HERE), "tests", "data", "lih_golden", "golden.json"),
                          encoding="utf-8"))
    items = gold.get("items", gold)
    Fbar = float(items["pi_marriage_exp1.5_1e-3" if N > 20 else "pi_marriage_exp1.5_1e-2"]["values"]
                 ["pieces"]["value"]["Fbar"]["repr"])
    psi_cols = [(r0, c, gm) for r0 in KEYS for c in CLASSES for gm in GAMMAS]
    K = 2 * N + len(psi_cols)
    pidx = {k: 2 * N + j for j, k in enumerate(psi_cols)}

    def to_lm(c):
        v = np.zeros(K)
        v[:N] += 0.5 * c[:N]
        v[N:2 * N] += 0.5 * c[N] * d
        v[:N] += -0.5 * c[N] * Fbar * d
        v[pidx[(-1.0, 'cv', 0.0)]] += c[N + 1] / 6.0
        v[pidx[(-1.0, 'cc', 0.0)]] += 2.0 * c[N + 1] / 6.0 + 0.5 * c[N + 2]
        if len(c) > N + 3:
            v[pidx[(-1.0, 'cc', 1.5)]] += 0.5 * c[N + 3] + c[N + 4]
            v[pidx[(-1.0, 'cv', 1.5)]] += 0.5 * c[N + 4]
        return v
    out = dict(ph1=to_lm(z1['vec_A_W12']), E_ph1=json.loads(str(z1['E']))['A_W12'])
    p2 = os.path.join(DATA, f"en_factors_{t:.0e}_phase2.npz")
    if os.path.exists(p2):
        z2 = np.load(p2, allow_pickle=True)
        out.update(ph2=to_lm(z2['vec_A_W12_G']), E_ph2=json.loads(str(z2['E']))['A_W12_G'])
    np.savez(os.path.join(DATA, f"en_factors_vmcvecs_{t:.0e}.npz"), **out)
    print({k: (v if np.ndim(v) == 0 else v.shape) for k, v in out.items()})


if __name__ == "__main__":
    main()
