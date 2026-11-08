"""Re-run TIMEOUT files with longer per-file timeout, fewer parallel workers."""
from __future__ import annotations
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS_PATH = ROOT / "debug" / "data" / "sprint_test_cleanup_timeout_rerun.json"

TIMEOUT_FILES = [
    "tests/test_1rdm_exchange.py",
    "tests/test_ab_initio_pk.py",
    "tests/test_ab_initio_pk_v2.py",
    "tests/test_composed_diatomic.py",
    "tests/test_composed_qubit.py",
    "tests/test_gh_convergence.py",
    "tests/test_gh_convergence_tensor.py",
    "tests/test_ihara_zeta_dirac.py",
    "tests/test_l_dependent_pk.py",
    "tests/test_level4_multichannel.py",
    "tests/test_lorentzian_lichnerowicz.py",
    "tests/test_lorentzian_propinquity.py",
    "tests/test_mo_fci.py",
    "tests/test_n_electron_2d.py",
    "tests/test_prolate_heteronuclear_scf.py",
    "tests/test_prolate_relaxed_ci.py",
    "tests/test_prolate_stress.py",
    "tests/test_r25_l3_lipschitz_bound.py",
]


def run_pytest_on_file(test_file: str) -> dict:
    cmd = [
        sys.executable, "-m", "pytest", test_file,
        "--no-header", "-q", "--tb=no",
        "-p", "no:cacheprovider",
    ]
    try:
        proc = subprocess.run(
            cmd, cwd=str(ROOT),
            capture_output=True, text=True,
            timeout=900,  # 15 minutes per file
        )
        out = proc.stdout + proc.stderr
        passed = failed = errors = skipped = 0
        for line in out.splitlines():
            ln = line.strip()
            if "passed" in ln or "failed" in ln or "error" in ln:
                parts = ln.replace(",", "").split()
                for i, tok in enumerate(parts):
                    if tok.isdigit() and i + 1 < len(parts):
                        kind = parts[i + 1]
                        n = int(tok)
                        if kind.startswith("passed"):
                            passed = n
                        elif kind.startswith("failed"):
                            failed = n
                        elif kind.startswith("error"):
                            errors = n
                        elif kind.startswith("skipped"):
                            skipped = n
        status = "GREEN" if (failed == 0 and errors == 0) else "BROKEN"
        if passed == 0 and failed == 0 and errors == 0 and proc.returncode != 0:
            status = "COLLECTION_ERROR"
        return {
            "file": test_file, "status": status, "returncode": proc.returncode,
            "passed": passed, "failed": failed, "errors": errors, "skipped": skipped,
            "tail": "\n".join(out.splitlines()[-10:]),
        }
    except subprocess.TimeoutExpired:
        return {
            "file": test_file, "status": "TIMEOUT_15min", "returncode": -1,
            "passed": 0, "failed": 0, "errors": 0, "skipped": 0, "tail": "TIMEOUT 900s",
        }


def main():
    print(f"Re-running {len(TIMEOUT_FILES)} TIMEOUT files with 15-min per-file timeout...")
    results = []
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(run_pytest_on_file, f): f for f in TIMEOUT_FILES}
        done = 0
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            done += 1
            print(f"[{done}/{len(TIMEOUT_FILES)}] {r['status']:18} {r['file']:60} F={r['failed']} E={r['errors']} P={r['passed']} S={r['skipped']}")
    results.sort(key=lambda r: r["file"])
    RESULTS_PATH.write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
