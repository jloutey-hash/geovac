"""FP step B2: the grid code with H -> chi (leaf d removed; H_0 = chi exactly) reproduces the chain
C(t) = (8/pi^3) t^3 int chi^2 K0(t x) dx of density_loop_B_ode.C_pos (step B's line used chi instead of chi^2:
corrected here)."""
import json, sys
from pathlib import Path
import mpmath as mp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_num as PN  # noqa: E402
import density_loop_B_ode as B  # noqa: E402
G = PN.Grid(30, 32)
c_grid = G.thetas("0.7", 0, which=lambda x, gx, cx, H: cx * cx)[0] * (mp.pi / 2)
mp.mp.dps = 25
c_ref = B.C_pos("0.7")
rec = {"grid_H_to_chi": mp.nstr(c_grid, 22), "C_pos": mp.nstr(c_ref, 22), "rel": mp.nstr(abs(c_grid / c_ref - 1), 3)}
print(rec)
d = json.loads((HERE / "data" / "four_density_path_B.json").read_text())
d["chain_leaf_removed"] = rec
d["chain_leaf_removed_note"] = "step B's first line integrated chi, not chi^2 (bug); corrected by four_density_path_B2_leaf.py"
(HERE / "data" / "four_density_path_B.json").write_text(json.dumps(d, indent=1))
