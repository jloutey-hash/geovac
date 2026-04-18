"""Merge sqrt_n fit into main hp_operator.json."""
import json
from pathlib import Path

ROOT = Path(__file__).parent.parent

main_path = ROOT / "debug" / "data" / "hp_operator.json"
fit_path = ROOT / "debug" / "data" / "hp_operator_sqrt_n_fit.json"

with open(main_path, "r", encoding="utf-8") as f:
    main = json.load(f)
with open(fit_path, "r", encoding="utf-8") as f:
    fit = json.load(f)

main["weyl_law_fits"] = fit

with open(main_path, "w", encoding="utf-8") as f:
    json.dump(main, f, indent=2)
print(f"Merged sqrt_n fits into {main_path}")
