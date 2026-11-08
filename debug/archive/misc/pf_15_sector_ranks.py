"""How Y's operator L (order 19) is spread over the six sectors of the off-slice tetrahedron system.

M = the rank-34 D-module of the masters, filtered by sector (debug/pf_14_sector_blocks.py: block triangular).  For a
downward-closed set I of sectors, the image of Y in M/M_I has a minimal operator L_I, a right factor of L, and
ord L_I - ord L_I' (I' = I plus one sector s) is the rank that sector s contributes to the composition factors of L.
Also: the minimal operator of a random vector in each diagonal block (its order is the block's cyclic rank).
Everything mod the sample primes, from the fitted systems in ib_tetra_pf{,_p1,_p2}_op.json.
usage: python pf_15_sector_ranks.py"""
from __future__ import annotations
import json, random, sys, time
from pathlib import Path
from flint import nmod_poly
sys.path.insert(0, str(Path(__file__).parent))
import ib_kite_de as K

DATA = Path(__file__).parent / "data"
FILES = ["ib_tetra_pf_op.json", "ib_tetra_pf_p1_op.json", "ib_tetra_pf_p2_op.json"]
ORDERS = {"A": [(3,), (4,), (3, 5), (4, 5), (3, 4), (3, 4, 5)],
          "B": [(4,), (3,), (3, 4), (4, 5), (3, 5), (3, 4, 5)]}


def sector(m):
    n = eval(m)
    return tuple(i for i in (3, 4, 5) if n[i] >= 0)


def rf(e, P):
    return K.RF(nmod_poly(e[0], P), nmod_poly(e[1], P))


def min_order(c, A, idx):
    cc = [c[i] for i in idx]
    if all(x.is_zero() for x in cc):
        return 0
    AA = [[A[i][j] for j in idx] for i in idx]
    return len(K.cyclic_operator(cc, AA)) - 1


def main():
    out = {}
    for fn in FILES:
        o = json.loads((DATA / fn).read_text())
        P = o["P"]
        K.P = P
        t0 = time.time()
        ms = o["masters"]
        sec = [sector(m) for m in ms]
        c = [rf(e, P) for e in o["c"]]
        A = [[rf(e, P) for e in row] for row in o["A"]]
        res = {"P": P}
        for name, order in ORDERS.items():
            ranks = []
            for k in range(len(order) + 1):
                ideal = set(order[:k])
                idx = [i for i in range(len(ms)) if sec[i] not in ideal]
                ranks.append(min_order(c, A, idx))
            contrib = {str(order[k]): ranks[k] - ranks[k + 1] for k in range(len(order))}
            res[f"order_{name}"] = {"ranks": ranks, "contribution": contrib}
            print(f"{fn} order {name}: ranks of Y mod M_I {ranks}; contribution per sector {contrib} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        rng = random.Random(7)
        blocks = {}
        for s in ORDERS["A"]:
            idx = [i for i in range(len(ms)) if sec[i] == s]
            v = [K.RF(nmod_poly([rng.randrange(1, P)], P)) for _ in idx]
            AA = [[A[i][j] for j in idx] for i in idx]
            blocks[str(s)] = {"size": len(idx), "cyclic_rank_random_vector": len(K.cyclic_operator(v, AA)) - 1}
        res["blocks"] = blocks
        print(f"   diagonal blocks {blocks} ({time.time()-t0:.0f}s)", flush=True)
        out[fn] = res
        (DATA / "pf_15_sector_ranks.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
