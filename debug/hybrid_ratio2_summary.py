"""Summary + label checks for the second-pipeline runs: applies `hybrid_ratio_summary.check` (unchanged, imported)
to every debug/data/hybrid_ratio2_runs/*.json, attaches the cross-check results, writes
debug/data/hybrid_ratio2_results.json and prints one line per run.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "debug"))

import hybrid_ratio_summary as HS  # noqa: E402


def main():
    out = {}
    for fn in sorted((ROOT / "debug/data/hybrid_ratio2_runs").glob("*.json")):
        run = json.loads(fn.read_text())
        res, form = run["result"], run["form"]
        chk = HS.check(res, form)
        shifts_c = list(res.get("shifts", {}))
        rates = sorted(set(res.get("exp_rates_of_f", [])) | set(shifts_c))
        cfn = ROOT / "debug/data/hybrid_ratio2_checks" / fn.name
        cross = json.loads(cfn.read_text())["points"] if cfn.exists() else None
        out[run["key"]] = {"order": res.get("order"), "d_E": res.get("d_E"), "k": res.get("k"),
                           "rates_of_f": rates, "n_rates": len(rates),
                           "order_minus_rates": (res.get("order") or 0) - len(rates),
                           "Lf_rates": {r: g["zero_factor"] for r, g in res.get("Lf_rates", {}).items()},
                           "E_zero_factors": {c: g["q_first_zero_factor"] for c, g in res.get("shifts", {}).items()},
                           "E_ratios": {c: g["ratio_to_first"] for c, g in res.get("shifts", {}).items()},
                           "constants": res.get("constants"), "issues": chk["issues"],
                           "cross_check": cross, "info": res.get("info")}
        mx = max((float(p["rel_spread"]) for p in cross.values()), default=None) if cross else None
        print(run["key"], "order", res.get("order"), "rates", len(rates), "dE", res.get("d_E"), "k", res.get("k"),
              "issues", chk["issues"], "cross pts", len(cross) if cross else 0, "max spread", mx)
    (ROOT / "debug/data/hybrid_ratio2_results.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
