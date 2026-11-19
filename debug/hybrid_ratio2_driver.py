"""Schedule the piece jobs of several (case, region) runs on a process pool, then combine each run.

Usage: python debug/hybrid_ratio2_driver.py <workers> <nch_lt> <nch_gt> <case|region> [<case|region> ...]
Each job: one subprocess `hybrid_ratio2_job.py`, 90-minute cap, log in debug/data/hybrid_ratio2_logs/.
Pieces already on disk are skipped (restartable). A run is combined as soon as all its pieces exist.
"""
from __future__ import annotations

import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PIECES = ROOT / "debug/data/hybrid_ratio2_pieces"
LOGS = ROOT / "debug/data/hybrid_ratio2_logs"
PY = sys.executable
CAP = 5400


def stem(case, region):
    return case.replace(",", "_") + "__" + region


def run(cmd, logf, cap=CAP):
    t0 = time.time()
    with open(logf, "w") as fh:
        try:
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, timeout=cap, cwd=ROOT).returncode
        except subprocess.TimeoutExpired:
            rc = "TIMEOUT"
    return rc, round(time.time() - t0, 1)


def main():
    workers, nlt, ngt = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    specs = [s.split("|") for s in sys.argv[4:]]
    PIECES.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    jobs, per_run = [], {}
    for case, region in specs:
        n = int(subprocess.run([PY, "debug/hybrid_ratio2_job.py", case, region, "count"], capture_output=True,
                               text=True, cwd=ROOT).stdout.strip())
        per_run[(case, region)] = 0
        for ti in range(n):
            for br, nch in (("gt", ngt), ("lt", nlt)):
                for ch in range(nch):
                    tag = f"{stem(case, region)}__t{ti}_{br}_{ch}of{nch}"
                    out = PIECES / (tag + ".json")
                    if out.exists():
                        continue
                    jobs.append(((case, region), [PY, "debug/hybrid_ratio2_job.py", case, region, str(ti), br,
                                                  str(ch), str(nch), str(out)], LOGS / (tag + ".log"), tag))
                    per_run[(case, region)] += 1
    print(f"{len(jobs)} jobs", flush=True)
    t0 = time.time()
    remaining = dict(per_run)
    failed = set()

    def combine(cr):
        case, region = cr
        rc, sec = run([PY, "debug/hybrid_ratio2_combine.py", case, region],
                      LOGS / f"combine_{stem(case, region)}.log")
        print(f"[{time.time() - t0:.0f}s] COMBINE {case}|{region} rc={rc} {sec}s", flush=True)

    cex = ThreadPoolExecutor(3)          # combines get their own pool (not queued behind the piece jobs)
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(run, cmd, logf): (cr, tag) for cr, cmd, logf, tag in jobs}
        for cr in [cr for cr, n in remaining.items() if n == 0]:
            cex.submit(combine, cr)
        for f in as_completed(futs):
            cr, tag = futs[f]
            rc, sec = f.result()
            print(f"[{time.time() - t0:.0f}s] {tag} rc={rc} {sec}s", flush=True)
            if rc != 0:
                failed.add(cr)
            remaining[cr] -= 1
            if remaining[cr] == 0 and cr not in failed:
                cex.submit(combine, cr)
    cex.shutdown(wait=True)
    print("failed runs:", sorted(failed), flush=True)


if __name__ == "__main__":
    main()
