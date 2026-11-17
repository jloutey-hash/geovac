"""KS supplementary (disclosed deviation): Windows process_time ticks at 15.625 ms, so route Q's ~30 ms is only
2-3 ticks. Here route Q (stage 1 + stage 2) is repeated 50 times inside one child and the process time divided by 50.
Five children per setting. Writes data/ks_timing_qloop.json. usage: python ks_time_q_loop.py [child <T> <h>]"""
from __future__ import annotations
import json, subprocess, sys, statistics, time
from pathlib import Path
HERE = Path(__file__).parent
REP = 50

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "child":
        import ks_route_q as Q
        tg, h = sys.argv[2], eval(sys.argv[3])
        Q.run(tg, h)                       # warm: outer-node file read once
        t0 = time.process_time(); p0 = time.perf_counter()
        for _ in range(REP):
            Q.run(tg, h)
        print(json.dumps({"cpu_per": (time.process_time() - t0) / REP, "wall_per": (time.perf_counter() - p0) / REP}))
    else:
        out = {}
        for tg, h in (("T12", "1/8"), ("T12", "1/16"), ("T9", "1/8")):
            rs = []
            for _ in range(5):
                r = subprocess.run([sys.executable, __file__, "child", tg, h], capture_output=True, text=True, cwd=HERE)
                rs.append(json.loads(r.stdout.strip().splitlines()[-1]))
            c = [x["cpu_per"] for x in rs]
            out[f"{tg} {h}"] = {"runs": rs, "median_cpu": statistics.median(c), "spread": max(c) / min(c)}
            print(tg, h, out[f"{tg} {h}"]["median_cpu"], out[f"{tg} {h}"]["spread"], flush=True)
        (HERE / "data" / "ks_timing_qloop.json").write_text(json.dumps(out, indent=1))
