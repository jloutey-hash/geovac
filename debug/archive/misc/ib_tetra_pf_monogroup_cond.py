"""diagnostic: conditioning of the scaled transport matrices U (base x_b -> loop point of each target)."""
import sys
sys.argv = [sys.argv[0], "g", "60", "1"]
import json
import mpmath as mp
import ib_tetra_pf_monogroup as M

M.tasks()
res = json.loads((M.DATA / "ib_tetra_pf_monogroup_transfer_d60.json").read_text())
R = M.R
Db = mp.diag([M.XB ** (-k) for k in range(R)])
U = mp.eye(R)
for lab, x0, r in M.TARGETS:
    for q in range(M.CHUNKS[lab]):
        U = M.to_mat(res[f"S {lab} {q}"]) * U
    Ut = M.to_mat(res[f"D {lab}"]) * U if r != M.H else U
    Dr = mp.diag([r ** k for k in range(R)])
    sv = mp.svd_c(Dr * Ut * Db, compute_uv=False)
    a = [abs(v) for v in sv]
    print(lab, "scaled U: largest sv", mp.nstr(max(a), 4), "smallest", mp.nstr(min(a), 4),
          "cond 10^%.1f" % float(mp.log10(max(a) / min(a))), flush=True)
