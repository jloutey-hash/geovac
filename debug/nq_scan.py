"""Phase-2 driver (record 2026-10-27_kite_ninth_point_origin).
  python nq_scan.py scan <tag> <p> <eps> <nsamp> <seed> [workers]   random z residues at fixed eps -> jsonl
  python nq_scan.py counts <tag> <p> <z:eps> [<z:eps> ...]           direct counts at rational (z, eps)
Output: debug/data/nq_scan_<tag>.jsonl (one line per sample, appended)."""
import json
import random
import sys
import time
from fractions import Fraction as Fr
from multiprocessing import Pool
from pathlib import Path

import nq_engine as E

DATA = Path(__file__).parent / "data"


def _work(args):
    p, zr, er, traces, label = args
    t0 = time.process_time()
    cnt, tr = E.sample(p, zr, er, want_traces=traces)
    return {"label": label, "p": p, "z": zr, "eps": er, "count": cnt, "traces": tr,
            "cpu": round(time.process_time() - t0, 1)}


if __name__ == "__main__":
    mode, tag, p = sys.argv[1], sys.argv[2], int(sys.argv[3])
    out = DATA / f"nq_scan_{tag}.jsonl"
    if mode == "scan":
        eps, n, seed = Fr(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
        workers = int(sys.argv[7]) if len(sys.argv) > 7 else 6
        rng = random.Random(seed)
        er = E.res(eps, p)
        jobs = [(p, rng.randrange(1, p), er, True, f"eps={eps}") for _ in range(n)]
    else:
        workers = 6
        jobs = []
        for spec in sys.argv[4:]:
            zs, es = spec.split(":")
            jobs.append((p, E.res(Fr(zs), p), E.res(Fr(es), p), False, spec))
    t0 = time.time()
    with Pool(min(workers, len(jobs))) as pool, open(out, "a") as fh:
        for r in pool.imap_unordered(_work, jobs):
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            print(r["label"], "count", r["count"], "cpu", r["cpu"], flush=True)
    print("wall", round(time.time() - t0, 1), "s", flush=True)
