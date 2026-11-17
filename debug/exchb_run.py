"""Record EXb driver (docs/predictions/2026-11-16_exchange_unequal_zeta_b.md).

Record EX was VOID at (a, b) = (1, 3/2): CRT + rational reconstruction of the operator failed with four primes
(M ~ 2^103).  This driver changes no frozen file.  It runs the frozen `debug/exch_run.py run 1 3/2 p` on more
primes, then the frozen `heldout 1 3/2` and `numeric 1 3/2` on them, all under the new suffix EXB.

Primes: the four EXS primes first (67108859, 33554393, 67108837, 67108819), then the primes below 2^26 in
descending order, skipping the third prime P3 = 16777213 and duplicates.  A candidate is ACCEPTED when its run gives
an operator with order 9, degree 16 and per-coefficient degrees equal to the EXS run on 67108859 (screens primes
that drop a leading coefficient), annihilating all complete rows, valuation check passed.  The first NACC = 32
accepted are used; rejected candidates are listed (cap 48 candidates).  Candidate runs are written under suffix
EXBC; accepted ones are copied to suffix EXB, which is all `heldout` reads.

Instrument checks recorded in debug/data/exchb_summary.json:
  - determinism: on the four EXS primes the EXBC operator mod p equals the EXS one coefficient for coefficient
    (a crash is recorded as an error string, not raised);
  - stability ladder: CRT + reconstruction on the first k accepted, k = 18, 20, ..., NACC; k0 = the smallest k
    that succeeds; the reconstruction at k0 + 4 (k0 + 4 <= accepted) must succeed and be equal (bits of the
    largest reconstructed numerator / denominator reported).
Then `heldout 1 3/2 --suffix EXB` and `numeric 1 3/2 --suffix EXB` (frozen code).

Usage: python debug/exchb_run.py   (refuses to start if any exch_run_1_3_2_*_EXB.json exists)
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time

import sympy as sp

sys.path.insert(0, "debug")
import exch_run as ER  # noqa: E402  (frozen; only crt_operator and P3 are used)

A, B = "1", "3/2"
TAG = "1_3_2"
EXS = [67108859, 33554393, 67108837, 67108819]
NACC, CAP, KSTEP = 32, 48, 2
OUT = "debug/data"
REF_DEG = json.load(open(os.path.join(OUT, f"exch_run_{TAG}_67108859_EXS.json")))["minimal_operator"]["coeff_degrees"]


def candidates():
    seen = list(EXS)
    p = sp.prevprime(2 ** 26)
    while len(seen) < CAP:
        if p not in seen and p != ER.P3:
            seen.append(int(p))
        p = sp.prevprime(p)
    return seen


def run_json(p, suf):
    return os.path.join(OUT, f"exch_run_{TAG}_{p}_{suf}.json")


def accept(p):
    if not os.path.exists(run_json(p, "EXBC")):
        return False, "no output"
    d = json.load(open(run_json(p, "EXBC")))
    m = d.get("minimal_operator")
    if (m and m["order"] == 9 and m["degree"] == 16 and m["coeff_degrees"] == REF_DEG
            and m["annihilates_all_complete_rows"] and d["generator"]["valuation_check"]):
        return True, ""
    shape = (m["order"], m["degree"], m["coeff_degrees"]) if m else None
    return False, f"shape {shape}, val {d['generator']['valuation_check']}"


def main():
    t0 = time.time()
    stale = [f for f in os.listdir(OUT) if f.startswith(f"exch_run_{TAG}_") and f.endswith("_EXB.json")]
    assert not stale, f"stale EXB files present: {stale}"
    summ = {"pair": [A, B], "accepted": [], "rejected": [], "NACC": NACC, "CAP": CAP, "ref_coeff_degrees": REF_DEG}
    for p in candidates():
        if len(summ["accepted"]) >= NACC:
            break
        rc = subprocess.call([sys.executable, "debug/exch_run.py", "run", A, B, str(p), "--suffix", "EXBC"])
        ok, why = accept(p) if rc == 0 else (False, f"rc={rc}")
        if ok:
            shutil.copyfile(run_json(p, "EXBC"), run_json(p, "EXB"))
            summ["accepted"].append(p)
        else:
            summ["rejected"].append({"p": p, "why": why})
        print(p, "accepted" if ok else "rejected " + why, flush=True)
    # determinism on the EXS primes
    det = {}
    for p in EXS:
        try:
            a = json.load(open(run_json(p, "EXS")))["minimal_operator"]["coeffs_mod_p"]
            b = json.load(open(run_json(p, "EXBC")))["minimal_operator"]["coeffs_mod_p"]
            det[str(p)] = a == b
        except Exception as e:
            det[str(p)] = f"error: {type(e).__name__}: {e}"
    summ["determinism_vs_EXS"] = det
    # stability ladder of the reconstruction
    runs = [json.load(open(run_json(p, "EXB"))) for p in summ["accepted"]]
    ladder, ops = {}, {}
    for k in range(18, len(runs) + 1, KSTEP):
        op, err = ER.crt_operator(runs[:k])
        ladder[str(k)] = {"ok": err is None, "error": err}
        if op is not None:
            ladder[str(k)].update({"max_num_bits": max(abs(x.numerator).bit_length() for c in op for x in c),
                                   "max_den_bits": max(x.denominator.bit_length() for c in op for x in c)})
            ops[k] = op
    k0 = min(ops) if ops else None
    stab = {"ladder": ladder, "k0": k0}
    stab["k0_plus_4_ok"] = k0 is not None and (k0 + 4) in ops
    stab["agree_k0_k0plus4"] = stab["k0_plus_4_ok"] and ops[k0] == ops[k0 + 4]
    stab["agree_all_successful"] = bool(ops) and all(ops[k] == ops[k0] for k in ops)
    summ["stability"] = stab
    summ["phase_a_time_s"] = round(time.time() - t0, 1)
    json.dump(summ, open(os.path.join(OUT, "exchb_summary.json"), "w"), indent=1)
    for mode in ("heldout", "numeric"):
        rc = subprocess.call([sys.executable, "debug/exch_run.py", mode, A, B, "--suffix", "EXB"])
        summ[f"{mode}_rc"] = rc
    summ["total_time_s"] = round(time.time() - t0, 1)
    json.dump(summ, open(os.path.join(OUT, "exchb_summary.json"), "w"), indent=1)
    print("EXB_DONE", flush=True)


if __name__ == "__main__":
    main()
