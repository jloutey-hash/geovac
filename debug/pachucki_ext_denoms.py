"""Reviewer item 1: every intermediate of the Q(t) run, not only K, for one case (a = 1, b = t).

Uniqueness by decay is used at every multi-index, so the proof needs, for every alpha in the downset:
  * the coefficient denominators of f_alpha, F_alpha and d_alpha p^2 to factor into t, t - 1, t + 1 only;
  * every function key to decay for EVERY t > 1 (rate lam, or lam + mu for Ei keys, linear in t), and no Ei
    argument mu to vanish at a t > 1.
The reviewer's debug/pachucki_ext_review_qt.py did this for 1s and 2psig; this driver does any case and saves it.

Engine: tests/_pachucki_extension.py (the guard's copy of debug/pachucki_ext_core.py, same code).  One Q(t) solve
serves three purposes:
  1. the audit above  -> debug/data/pachucki_ext_denoms_<CASE>.json  (verdict PASS / FAIL);
  2. on PASS, tests/data/paper58_pachucki_ext_2026_11_18.json is rebuilt (debug/pachucki_ext_guard_data.py);
  3. the body of the guard's slow test C3 (test_Qt_run_recomputed_residuals_audit_and_closed_form) is then called
     on the same in-process solution: residuals identically 0 by the test's own formula, audit = frozen record,
     K coefficients = the author's frozen Q(t) coefficients.  Its outcome is recorded as "C3_test_body".
Exit code 0 on PASS, 3 on FAIL (the caller stops).
Usage: python debug/pachucki_ext_denoms.py CASE        (CASE = 1s | 2psig | 3dsig | ppsig | pppi0)
"""
import json
import subprocess
import sys
import time
import traceback

sys.path.insert(0, "tests")
import _pachucki_extension as X  # noqa: E402
import test_paper58_pachucki_extension as T  # noqa: E402

case = sys.argv[1]
t0 = time.time()
print(case, "multi-indices", len(X.case_downset(case)), flush=True)
ctx, D, calc, fS, Fs, P2, nfree, K = T.solved_qt(case)
tsolve = time.time() - t0
print("solved over Q(t); free unknowns", nfree, "%.0fs" % tsolve, flush=True)

au = X.audit_intermediates(ctx, D, fS, Fs, P2)
dens = au["denominator_factors_max_multiplicity"]
print("denominator factors (max multiplicity):", dens, flush=True)
ALLOWED = {"t", "t - 1", "t + 1"}
dens_ok = all(set(v) <= ALLOWED for v in dens.values())
rows, bad = [], []
for r in au["keys"]:
    v = X.key_verdict(r)
    row = dict(r, rate=v["rate"], decays_for_all_t_gt_1=v["decays"], ei_argument_vanishes_at_t_gt_1=v["ei_zero"])
    rows.append(row)
    if not v["decays"] or v["ei_zero"]:
        bad.append(row)
    print("  key", row, flush=True)
verdict = "PASS" if (dens_ok and not bad and nfree == 0) else "FAIL"
print("VERDICT", verdict, "| denominators ok:", dens_ok, "| bad keys:", bad, flush=True)
out = {"case": case, "engine": "tests/_pachucki_extension.py", "multi_indices": len(D), "free_unknowns": nfree,
       "denominator_factors_max_multiplicity": dens, "denominators_within_t_tm1_tp1": dens_ok,
       "keys": rows, "bad_keys": bad, "verdict": verdict, "solve_time_s": round(tsolve, 1)}
fn = f"debug/data/pachucki_ext_denoms_{case}.json"
json.dump(out, open(fn, "w"), indent=1)
if verdict != "PASS":
    print("FAILED; stopping")
    sys.exit(3)

# 2. rebuild the guard's frozen data, 3. run the body of C3 on this solution
subprocess.check_call([sys.executable, "debug/pachucki_ext_guard_data.py"])
T.DATA = json.loads(open("tests/data/paper58_pachucki_ext_2026_11_18.json").read())
t1 = time.time()
try:
    if case == "1s":
        T.test_1s_fifteen_residuals_vanish_identically_over_Qt()
        T.test_1s_every_intermediate_has_denominators_t_tm1_tp1_and_decays_for_all_t()
    else:
        T.test_Qt_run_recomputed_residuals_audit_and_closed_form.__wrapped__(case) \
            if hasattr(T.test_Qt_run_recomputed_residuals_audit_and_closed_form, "__wrapped__") \
            else T.test_Qt_run_recomputed_residuals_audit_and_closed_form(case)
    c3 = "passed"
except Exception:                                   # noqa: BLE001
    c3 = "FAILED: " + traceback.format_exc()[-1500:]
print("C3 test body:", c3, "%.0fs" % (time.time() - t1), flush=True)
out["C3_test_body"] = c3
# cases with the same downset share the solve (3dsig, ppsig, pppi0): re-run their C3 body on this process's solution
same = {}
for other in X.CASES:
    if other != case and other != "1s" and X.case_downset(other) == D and other in T.DATA["denoms"]:
        try:
            T.test_Qt_run_recomputed_residuals_audit_and_closed_form(other)
            same[other] = "passed"
        except Exception:                           # noqa: BLE001
            same[other] = "FAILED: " + traceback.format_exc()[-800:]
        print("C3 test body for", other, "(same downset):", same[other], flush=True)
out["C3_test_body_same_downset"] = same
out["C3_time_s"] = round(time.time() - t1, 1)
out["time_s"] = round(time.time() - t0, 1)
json.dump(out, open(fn, "w"), indent=1)
print("done %.0fs" % (time.time() - t0))
sys.exit(0 if c3 == "passed" else 4)
