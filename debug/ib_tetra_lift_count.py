"""count how many system entries lift from 2 and from 3 primes (diagnostic)."""
import json, sys
from pathlib import Path
from ib_lift_system import lift_rf
ops = [json.loads(Path(a).read_text()) for a in sys.argv[1:]]
nm = len(ops[0]["c"])
ents = [("c", i, None) for i in range(nm)] + [("A", i, j) for i in range(nm) for j in range(nm)]
def get(o, e):
    return o["c"][e[1]] if e[0] == "c" else o["A"][e[1]][e[2]]
for k in (2, 3):
    use = ops[:k]; mods = [o["P"] for o in use]
    fail = [e for e in ents if lift_rf([get(o, e) for o in use], mods) is None]
    print(k, "primes: failed", len(fail), "of", len(ents), fail[:12], flush=True)
