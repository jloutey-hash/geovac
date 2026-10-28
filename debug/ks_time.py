"""KS timing driver (record docs/predictions/2026-10-28_kite_speed_test.md). One child process at a time.
C6: route Q at h = 1/8 and 1/16 (T12), three repeats each, median ratio.
Main: per target, five interleaved repeats O, Q, O, Q, ... at the scored rungs (O: 60 digits; Q: h = 1/8),
up to ten if max/min > 1.3 (C7). Writes data/ks_timing.json."""
from __future__ import annotations
import json, subprocess, sys, statistics
from pathlib import Path
HERE = Path(__file__).parent
RUNG_O = {"T12": 60, "T9": 60}
RUNG_Q = {"T12": "1/8", "T9": "1/8"}


def child(args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True, cwd=HERE)
    if r.returncode:
        raise RuntimeError(r.stderr[-2000:])
    return json.loads(r.stdout.strip().splitlines()[-1])


def qcost(d):
    return d["cpu"]["stage1"] + d["cpu"]["stage2"]


def ocost(d):
    c = d["cpu"]
    return c["load"] + c["input"] + c["walk"] + c["stage2"]


def main():
    out = {"C6": {}, "runs": {}}
    for h in ("1/8", "1/16"):
        out["C6"][h] = [child(["ks_route_q.py", "T12", h, "correct"]) for _ in range(3)]
    m8 = statistics.median(qcost(d) for d in out["C6"]["1/8"])
    m16 = statistics.median(qcost(d) for d in out["C6"]["1/16"])
    out["C6_ratio"] = m16 / m8
    print("C6 ratio", out["C6_ratio"], flush=True)
    for tg in ("T12", "T9"):
        O, Q = [], []
        n = 0
        while True:
            O.append(child(["ks_route_o.py", tg, str(RUNG_O[tg]), "ok"]))
            Q.append(child(["ks_route_q.py", tg, RUNG_Q[tg], "correct"]))
            n += 1
            if n >= 5:
                so = max(map(ocost, O)) / min(map(ocost, O))
                sq = max(map(qcost, Q)) / min(map(qcost, Q))
                if (so <= 1.3 and sq <= 1.3) or n >= 10:
                    break
        out["runs"][tg] = {"O": O, "Q": Q}
        mo, mq = statistics.median(map(ocost, O)), statistics.median(map(qcost, Q))
        print(tg, "repeats", n, "O median", mo, "spread", max(map(ocost, O)) / min(map(ocost, O)),
              "Q median", mq, "spread", max(map(qcost, Q)) / min(map(qcost, Q)),
              "R_marg", mo / mq, "R_tot", (mo + 9200) / mq, flush=True)
        (HERE / "data" / "ks_timing.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
