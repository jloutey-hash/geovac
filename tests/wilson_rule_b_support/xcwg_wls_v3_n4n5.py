"""Continue v3 sprint at n_max=4 and n_max=5 only (v3 was killed during n_max=4)."""
import os, sys, json

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))  # repo root (tests/wilson_rule_b_support/ -> two levels up)
sys.path.insert(0, _ROOT); sys.path.insert(0, _HERE)
sys.stdout.reconfigure(line_buffering=True)

from xcwg_wls_v3 import run_n_max

out = {"continuation": "n_max=4 and n_max=5 only"}

# n_max=4: cap L=6 simple at 300k (it streams about 30k simple cycles per sec on this graph)
# cap L=8 simple at 300k as well
out["n_max_4"] = run_n_max(
    4, [4, 6, 8],
    caps_full={4: 10_000_000, 6: 10_000_000, 8: 600_000},
    caps_simple={4: 10_000_000, 6: 300_000, 8: 300_000},
)
out["n_max_5"] = run_n_max(
    5, [4, 6],
    caps_full={4: 500_000, 6: 400_000},
    caps_simple={4: 500_000, 6: 300_000},
)
out["n_max_5"]["L8_skipped"] = "Compute-intractable at L=8 on n_max=5"

with open(os.path.join(_HERE, 'data', 'xcwg_wls_v3_n4n5.json'), 'w') as f:
    json.dump(out, f, indent=2, default=str)

# Compact table
print()
print("="*78)
print("PARTIAL TABLE (v3 n_max=4,5)")
print("="*78)
for k in ["n_max_4", "n_max_5"]:
    f = out[k]["scaling_fit_full"]
    s = out[k]["scaling_fit_simple"]
    nm = k.replace("n_max_", "")
    print(f"  n_max={nm}: FULL alpha={f['alpha']:.4f} ({f['verdict']});  "
          f"SIMPLE alpha={s['alpha']:.4f} ({s['verdict']})")
