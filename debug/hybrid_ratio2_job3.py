"""v3 of hybrid_ratio2_job.py: identical pieces, faster (see hybrid_ratio2_fastrep.py: flint-backed collector
field, term-by-term integration + collection).

One piece of a parametric hybrid closed form, collected exactly and written to JSON.

Same algebra as `hybrid_ratio_engine.shell_parts` (which is a line-for-line copy of geovac's
`hybrid_closed_form_shell`), restricted to ONE multipole term (index `ti` of `T.hybrid_terms`), ONE branch
(x < R -> `lt`, x > R -> `gt`) and, optionally, one chunk of the inner double's terms (round-robin over the
Add-args of the fixed-x double integral, taken before the weight is applied; `_integrate_with_e1` is a term-by-term sum, so the chunks sum to the piece).
The piece 2 pi / R * coeff * integral is then passed through `hybrid_ratio_rep.rep` (the existing exact
collector, unchanged) and serialised: key -> (numerator terms, denominator terms) over Q[u, R].

The full form is the sum over all (term, branch, chunk) pieces (`hybrid_ratio2_combine.py`).
Nothing here edits the existing hybrid_ratio_* files; they are imported only.

Usage: python debug/hybrid_ratio2_job3.py <a,b,c,d> <region> <ti> <lt|gt> <chunk> <nchunks> <outfile>
       python debug/hybrid_ratio2_job3.py <a,b,c,d> <region> count          (prints the number of terms)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "debug"))
sys.path.insert(0, str(ROOT))

import hybrid_ratio_engine as E  # noqa: E402  (applies the symbolic-Z helpers to geovac.two_center_eri)
import hybrid_ratio_rep as HR  # noqa: E402
import hybrid_ratio2_fastrep as FR  # noqa: E402

FR.patch_field()

T = E.T


def ser_rep(rep):
    out = []
    for (c, kind, kap, consts), v in rep.items():
        num = [[[int(x) for x in m], str(q)] for m, q in v.n.to_dict().items()]
        den = [[[int(x) for x in m], str(q)] for m, q in v.d.to_dict().items()]
        out.append({"c": str(c), "kind": kind, "kap": str(kap),
                    "consts": sorted([n, str(e)] for n, e in consts), "num": num, "den": den})
    return out


def piece(case, region, ti, branch, chunk, nchunks, log):
    a, b, c, d = (E.ORB[x] for x in case)
    ZA, ZB, _ap, _ac, _ad = E.region_Z(region, a, b, c, d)
    t0 = time.time()
    terms = list(T.hybrid_terms(ZA, a, b, c, ZB, d))
    coeff, L, Lp, Mp, ld, md = terms[ti]
    w_s = sp.Symbol("w", positive=True)
    x_s = T.x_s
    rad1, b1 = None, None
    for LL, _M, _g, rad, bb in T.multipole_decomposition(ZA, *a, ZA, *b):
        if LL == L:
            rad1, b1 = rad, bb
            break
    weight = (4 * sp.pi / (2 * L + 1) * sum(cc * x_s ** (k + 2) for k, cc in rad1.items()) * sp.exp(-b1 * x_s))
    e_br = T.shell_inner_double(ZA, a, b, c, ZB, d, L, Lp, Mp, ld, md, x_s, E.R, branch=branch)
    if branch == "gt":
        assert not e_br.atoms(sp.Ei)
    log(f"term {ti} L={L} L'={Lp} {branch}: inner double {time.time() - t0:.1f}s")
    # v2: chunk the inner double's Add-args BEFORE multiplying by the weight (every later step is linear and
    # term-wise; the x-power assert below still runs on every product term).  v1 chunked after the x -> R + w
    # expansion, which repeated that expansion in every chunk.
    eargs = sp.Add.make_args(e_br)
    e_br = sp.Add(*[t for i, t in enumerate(eargs) if i % nchunks == chunk])
    log(f"  inner double {len(eargs)} terms, chunk {chunk}/{nchunks}: {len(sp.Add.make_args(e_br))}")
    ex = E.expand_fixed(weight * e_br)
    for term in sp.Add.make_args(ex):
        pw = sum(int(f.exp) if (f.is_Pow and f.base == x_s) else (1 if f == x_s else 0)
                 for f in sp.Mul.make_args(term))
        assert pw >= 0
    if branch == "gt":
        ex = E.expand_fixed(ex.subs(x_s, E.R + w_s))
        var, X = w_s, sp.oo
    else:
        var, X = x_s, E.R
    args = sp.Add.make_args(ex)
    mine = ex
    log(f"  integrand {len(args)} terms {time.time() - t0:.1f}s")
    rep = FR.integrate_collect(T, mine, var, X, log=log)
    log(f"  integral+collect {time.time() - t0:.1f}s")
    rep = HR.mul(HR.rep(2 * sp.pi / E.R * coeff), rep)
    log(f"  collected {time.time() - t0:.1f}s, {len(rep)} keys")
    return rep, dict(L=L, Lp=Lp, n_integrand_terms=len(args), seconds=round(time.time() - t0, 1))


def main():
    case = tuple(sys.argv[1].split(","))
    region = sys.argv[2]
    lg = lambda m: print(m, flush=True)  # noqa: E731
    if sys.argv[3] == "count":
        a, b, c, d = (E.ORB[x] for x in case)
        ZA, ZB, *_ = E.region_Z(region, a, b, c, d)
        print(len(list(T.hybrid_terms(ZA, a, b, c, ZB, d))))
        return
    ti, branch, chunk, nchunks, outfile = int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), int(sys.argv[6]), sys.argv[7]
    rep, meta = piece(case, region, ti, branch, chunk, nchunks, lg)
    Path(outfile).write_text(json.dumps({"case": sys.argv[1], "region": region, "ti": ti, "branch": branch,
                                         "chunk": chunk, "nchunks": nchunks, "meta": meta, "rep": ser_rep(rep)}))
    lg("done")


if __name__ == "__main__":
    main()
