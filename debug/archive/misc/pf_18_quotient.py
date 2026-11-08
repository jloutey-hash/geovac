"""Y's operator after the relation the IBP reduction cannot see: the two sunrise sectors are one function.

The sector {g1} (geminal ik alone, no Coulomb) and {g2} (geminal jk alone) factorise: the q-loop is a one-line tadpole,
int d^3q / ((q - k)^2 + g^2) = const * g in D = 3 (dimensional regularisation; index 1 for all eight such masters), times
the same density sunrise.  So M_{g1}(s) = (g1/g2) M_{g2}(s) = M_{g2}(s)/2 for masters with the same density shifts s.
The span K of M_{g1}(s) - lam M_{g2}(s) is a D-submodule for every constant lam (the two blocks are equal, pf_17), and
lies in Y's cyclic module (both sunrise blocks contribute rank 4, pf_15).  In M/K Y's minimal operator has order 19 - 4.
Mod p here: the order, for lam = 1/2 (the tadpole ratio) and the control lam = 1/3; checks that the blocks are equal.
usage: python pf_18_quotient.py"""
from __future__ import annotations
import json, sys
from pathlib import Path
from flint import nmod_poly
sys.path.insert(0, str(Path(__file__).parent))
import ib_kite_de as K

DATA = Path(__file__).parent / "data"


def quotient(ms, c, A, lam_num, lam_den, P):
    """drop the {g1} masters, M_{g1}(s) -> lam M_{g2}(s); returns (idx, c', A')"""
    def key(m):
        n = eval(m)
        return n[:3]
    g1 = {key(m): i for i, m in enumerate(ms) if eval(m)[3] >= 0 and eval(m)[4] < 0 and eval(m)[5] < 0}
    g2 = {key(m): i for i, m in enumerate(ms) if eval(m)[3] < 0 and eval(m)[4] >= 0 and eval(m)[5] < 0}
    assert set(g1) == set(g2) and len(g1) == 4
    pair = {g1[k]: g2[k] for k in g1}
    lam = K.RF(nmod_poly([lam_num * pow(lam_den, -1, P) % P], P))
    keep = [i for i in range(len(ms)) if i not in pair]
    pos = {i: k for k, i in enumerate(keep)}
    cn = [c[i] for i in keep]
    for i, j in pair.items():
        cn[pos[j]] = cn[pos[j]] + lam * c[i]
    An = [[A[i][j] for j in keep] for i in keep]
    for r, i in enumerate(keep):
        for a, b in pair.items():
            if not A[i][a].is_zero():
                An[r][pos[b]] = An[r][pos[b]] + lam * A[i][a]
    eq = all((A[a][a2] - A[pair[a]][pair[a2]]).is_zero() for a in pair for a2 in pair) and \
        all(A[a][j].is_zero() for a in pair for j in range(len(ms)) if j not in pair)
    return keep, cn, An, eq


def main():
    out = {}
    for fn in ["ib_tetra_pf_op.json", "ib_tetra_pf_p1_op.json"]:
        o = json.loads((DATA / fn).read_text())
        P = o["P"]; K.P = P
        c = [K.RF(nmod_poly(e[0], P), nmod_poly(e[1], P)) for e in o["c"]]
        A = [[K.RF(nmod_poly(e[0], P), nmod_poly(e[1], P)) for e in row] for row in o["A"]]
        for ln, ld in ((1, 2), (1, 3)):
            keep, cn, An, eq = quotient(o["masters"], c, A, ln, ld, P)
            a = K.cyclic_operator(cn, An)
            degs = [x.n.degree() for x in a]
            out[f"{fn} lam={ln}/{ld}"] = {"order": len(a) - 1, "blocks_equal": eq}
            print(f"{fn} lam = {ln}/{ld}: blocks equal {eq}; Y's order in M/K = {len(a)-1}", flush=True)
    (DATA / "pf_18_quotient.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
