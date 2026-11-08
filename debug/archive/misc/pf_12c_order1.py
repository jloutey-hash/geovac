"""PF1, PF3, PF4, GR1: does F_8 (the order-12 factor) have an order-one right or left factor?  The exponent-bounded
hyperexponential search of record g32 (correlated_triangle_g32_lib4: every exponent choice allowed by the local data
and the Fuchs bound, then a polynomial solve), on F_8 (right factors) and on its adjoint (left factors).  Control: the
same search on chain_R2 = F_2 o F_1, which has the order-one right factor F_1 by construction.
usage: python pf_12c_order1.py"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib4 as S
import pf_lib  # noqa: F401  (sets the eight cut points for the searches; without it the control fails)

DATA = Path(__file__).parent / "data"
d = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
F8 = G.parse_op(d["factor_8"])
out = {}
for nm, op in (("control chain_R2", G.parse_op(d["chain_R2"])), ("F_8", F8), ("F_8*", G.adjoint(F8))):
    t0 = time.time()
    c = S.candidates(op, "sol")
    hits = [(str(f0), {str(k): str(v) for k, v in cm.items() if v}, N) for f0, cm, N in c if S.test_candidate(op, f0, cm, N, "sol")]
    out[nm] = {"admissible_choices": len(c), "hits": hits, "seconds": time.time() - t0}
    print(nm, "admissible:", len(c), "hits:", hits, f"({time.time()-t0:.0f}s)", flush=True)
(DATA / "pf_12c_order1.json").write_text(json.dumps(out, indent=1))
print("saved", flush=True)
