"""Build tests/data/paper59_kite_p2_2026_10_27.json for the PK / NP / ND / OE guards (records
2026-10-26_parametric_kite_operator, 2026-10-26_kite_ninth_point, 2026-10-27_kite_ninth_point_origin,
2026-10-26_operator_evaluator). Copies stored data only; computes nothing.
usage: python debug/build_kite_guard_data_2026_10_27.py"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "debug" / "data"
op = json.loads((D / "ib_pko_p2_op_Q.json").read_text())
grid = json.loads((D / "pf_grid_d60_h0.03125.json").read_text())
sA = [json.loads((D / f"np_stageA_{p}.json").read_text()) for p in ("p0", "p1")]
sC = json.loads((D / "np_stageC_dense_samples.json").read_text())
out = {
    "note": "Kite guards 2026-10-27. p2_op: the kite's order-11 IBP operator at (a,b,g1,g2) = (1,3/2,7/3,8/3), "
            "sum_k polys[k](z) (d/dz)^k, ascending integer coefficients (debug/data/ib_pko_p2_op_Q.json). "
            "grid_nodes: two nodes of the 60-digit Y grid at (1,2;3/2,3), w and the kite phi(w^2) by quadrature "
            "(debug/data/pf_grid_d60_h0.03125.json, indices 116 and 123). stageA: z9 measured mod P on the family "
            "(1,2,3/2,g2), s = g2^2 (np_stageA_p0/p1). stageC: z9 measured mod P at 53 random (b^2, g1^2, g2^2) with "
            "a^2 = 1 (np_stageC_dense_samples).",
    "p2_op": {"polys": op["polys"], "masters": op["masters"], "masses": op["masses"]},
    "grid_nodes": [{"index": i, "w": grid["w"][i], "phi": grid["phi"][i]} for i in (116, 123)],
    "stageA": [{"P": d["P"], "s": d["s"], "z9": d["z9"]} for d in sA],
    "stageC": {"P": sC["P"], "pts": sC["pts"], "vals": sC["vals"]},
}
(ROOT / "tests" / "data" / "paper59_kite_p2_2026_10_27.json").write_text(json.dumps(out))
print("written")
