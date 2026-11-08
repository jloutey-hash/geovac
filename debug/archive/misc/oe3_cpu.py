"""Measure the processor time of one IB4 prime run (p5), sample by sample.

Each worker does one task (maxtasksperchild=1), so time.process_time() at the end of the task is that worker
process's whole CPU, imports included. The main process's own CPU (rational rebuild) is measured separately by
rebuilding the system the same way ib_kite_de.main does, minus writing any file.
"""
import json
import os
import random
import sys
import time
from multiprocessing import Pool

DEBUG = r"C:\Users\jlout\Desktop\Project_Geometric\debug"
sys.path.insert(0, DEBUG)
os.chdir(DEBUG)
PRIMES = open(os.path.join(DEBUG, "data", "ib_primes.txt")).read().split()
os.environ["IB_PRIME"] = PRIMES[5]

import ib_kite_de as K  # noqa: E402
from fractions import Fraction as F  # noqa: E402


def timed(args):
    r = K.sample(args)
    return r, time.process_time()


def main():
    masses = [F(1), F(2), F(3, 2), F(3)]
    nsamp = 80
    rng = random.Random(12345)
    zs = [rng.randrange(10**6, K.P - 1) for _ in range(nsamp)]
    t0 = time.time()
    with Pool(6, maxtasksperchild=1) as pool:
        out = pool.map(timed, [(z, masses, 5, 1, False, F(1)) for z in zs], chunksize=1)
    wall = time.time() - t0
    cpu = [c for _, c in out]
    res = {"prime_index": 5, "prime": PRIMES[5], "nsamp": nsamp, "pool_wall_s": wall,
           "worker_cpu_total_s": sum(cpu), "worker_cpu_min_s": min(cpu), "worker_cpu_max_s": max(cpu),
           "main_cpu_s_after_pool": time.process_time()}
    print(json.dumps(res, indent=1), flush=True)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "oe3_cpu_p5.json"), "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
