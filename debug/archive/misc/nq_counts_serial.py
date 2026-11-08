"""Direct counts, one process (record 2026-10-27_kite_ninth_point_origin, NM0).
usage: python nq_counts_serial.py <tag> <p> <z:eps> ...   -> debug/data/nq_counts_<tag>.jsonl"""
import json
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import nq_engine as E

if __name__ == "__main__":
    tag, p = sys.argv[1], int(sys.argv[2])
    out = Path(__file__).parent / "data" / f"nq_counts_{tag}.jsonl"
    with open(out, "a") as fh:
        for spec in sys.argv[3:]:
            zs, es = spec.split(":")
            t0 = time.process_time()
            cnt, _ = E.sample(p, E.res(Fr(zs), p), E.res(Fr(es), p), want_traces=False)
            r = {"z": zs, "eps": es, "p": p, "count": cnt, "cpu": round(time.process_time() - t0, 1)}
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            print(r, flush=True)
