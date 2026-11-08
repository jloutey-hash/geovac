"""FP step O: is the Sym^3 function P_h4 = N t^3 int K0(t x) h4(x) dx, h4 = x^6 K0^3 = (c1 - c2)^2 k0
(theta chi = x^2 K1, theta^2 chi - theta chi = -x^3 K0), in the t-module of P?
x-side: sum_{i<=I, j<=11} b_ij T^i theta^j g = sum_{i<=I} q_i T^i h4 (mod p nullspace, q != 0).
Also: the minimal t-order of P_h4 (x-scan from h4) -- expected 4 (F_4's class)."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
import flint

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_A_operator as A  # noqa: E402
import four_density_path_C_tmodule as C  # noqa: E402

OUT = HERE / "data" / "four_density_path_O.json"


def vec_of(basis, elem, mp_):
    E = np.zeros((mp_.n, mp_.W), dtype=np.int64)
    idx = {m: i for i, m in enumerate(basis)}
    for mono, cf in elem.items():
        E[idx[tuple(sorted(mono))], -mp_.lo] = (cf % mp_.p)
    return E


H4 = {("c1", "c1", "k0"): 1, ("c1", "c2", "k0"): -2, ("c2", "c2", "k0"): 1}


def search(r, I, Ih, p=C.P0):
    M = A.path_module()
    basis, A0, A1 = C.module_matrices(M, ("c0", "h0"))
    for m in H4:
        assert tuple(sorted(m)) in basis, m
    mp_ = C.ModP(A0, A1, p, -max(I, Ih) - 2, max(I, Ih) + r + 6)
    cols = C.columns(mp_, r, I)
    keys = [(i, j) for i in range(I + 1) for j in range(r + 1)]
    h = vec_of(basis, H4, mp_)
    hcols = [h]
    for _ in range(Ih):
        hcols.append(mp_.T(hcols[-1]))
    mats = np.stack([cols[k].reshape(-1) for k in keys] + [x.reshape(-1) for x in hcols], axis=1)
    mats = mats[np.nonzero(mats.any(axis=1))[0]]
    Mm = flint.nmod_mat(mats.shape[0], mats.shape[1], [int(x) for x in mats.reshape(-1)], p)
    X, nul = Mm.nullspace()
    nk = len(keys)
    with_q = 0
    for k in range(nul):
        if any(int(X[nk + i, k]) for i in range(Ih + 1)):
            with_q += 1
    # rank of the q-part of the nullspace
    qrows = [[int(X[nk + i, k]) for k in range(nul)] for i in range(Ih + 1)]
    qrank = flint.nmod_mat(Ih + 1, nul, [x for rr in qrows for x in rr], p).rank() if nul else 0
    return nul, qrank


def main() -> None:
    rec = json.loads(OUT.read_text()) if OUT.exists() else {}
    t0 = time.time()
    # minimal order of P_h4 alone
    M = A.path_module()
    basis, A0, A1 = C.module_matrices(M, ("c0", "h0"))
    for r in range(1, 7):
        found = None
        for I in range(0, 40):
            mp_ = C.ModP(A0, A1, C.P0, -I - 2, I + r + 6)
            h = vec_of(basis, H4, mp_)
            th = [h]
            for _ in range(r):
                th.append(mp_.theta(th[-1]))
            cols = []
            for E in th:
                cur = E
                for i in range(I + 1):
                    cols.append(cur.reshape(-1))
                    if i < I:
                        cur = mp_.T(cur)
            mats = np.stack(cols, axis=1)
            mats = mats[np.nonzero(mats.any(axis=1))[0]]
            Mm = flint.nmod_mat(mats.shape[0], mats.shape[1], [int(x) for x in mats.reshape(-1)], C.P0)
            if mats.shape[1] - Mm.rank():
                found = I
                break
        print(f"P_h4: order {r}: {'t-degree %d' % (2*found) if found is not None else 'none to t-degree 78'}", flush=True)
        if found is not None:
            rec["P_h4_min"] = [r, 2 * found]
            break
    for I in (20, 40, 60, 80):
        nul, qrank = search(11, I, I)
        print(f"I={I}: nullity {nul}, rank of the q-part {qrank} ({time.time()-t0:.0f}s)", flush=True)
        rec[f"membership_I{I}"] = [nul, qrank]
        OUT.write_text(json.dumps(rec, indent=1))
        if qrank:
            break


if __name__ == "__main__":
    main()
