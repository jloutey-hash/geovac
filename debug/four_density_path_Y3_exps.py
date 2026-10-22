"""FP step Y3: local exponents of Z2 and Z2' at t = 1, 3 (and 0, infinity)."""
import json, sys
from pathlib import Path
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_L_exps as LE  # noqa: E402
X = json.loads((HERE / "data" / "four_density_path_X.json").read_text())
rec = {k: LE.exps(X[k], k, pts=(1, 3)) for k in ["Z2p_Q", "Z2_Q"]}
(HERE / "data" / "four_density_path_Y3.json").write_text(json.dumps(rec, indent=1))
