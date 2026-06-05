"""Parallel per-file pytest scanner for sprint test cleanup."""
from __future__ import annotations
import json
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests"
RESULTS_PATH = ROOT / "debug" / "data" / "sprint_test_cleanup_scan.json"


def run_pytest_on_file(test_file: str) -> dict:
    """Run pytest on a single file, return parsed status."""
    rel = test_file
    cmd = [
        sys.executable, "-m", "pytest", rel,
        "--no-header", "-q", "--tb=no",
        "-p", "no:cacheprovider",
    ]
    try:
        proc = subprocess.run(
            cmd, cwd=str(ROOT),
            capture_output=True, text=True,
            timeout=180,
        )
        out = proc.stdout + proc.stderr
        # Parse the summary line
        passed = failed = errors = skipped = 0
        for line in out.splitlines():
            ln = line.strip()
            if "passed" in ln or "failed" in ln or "error" in ln:
                # Parse like "5 failed, 12 passed in 1.5s" or "27 passed, 15 skipped in 8.46s"
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
        # Detect collection errors (no tests ran)
        if passed == 0 and failed == 0 and errors == 0 and proc.returncode != 0:
            status = "COLLECTION_ERROR"
        return {
            "file": rel, "status": status, "returncode": proc.returncode,
            "passed": passed, "failed": failed, "errors": errors, "skipped": skipped,
            "tail": "\n".join(out.splitlines()[-5:]),
        }
    except subprocess.TimeoutExpired:
        return {
            "file": rel, "status": "TIMEOUT", "returncode": -1,
            "passed": 0, "failed": 0, "errors": 0, "skipped": 0, "tail": "TIMEOUT 180s",
        }


def main():
    test_files = sorted(
        str(p.relative_to(ROOT)).replace("\\", "/")
        for p in TESTS.glob("test_*.py")
    )
    print(f"Scanning {len(test_files)} test files in parallel...")
    results = []
    with ProcessPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(run_pytest_on_file, f): f for f in test_files}
        for i, fut in enumerate(as_completed(futs)):
            r = fut.result()
            results.append(r)
            status_short = r["status"]
            if status_short == "GREEN":
                marker = "."
            elif status_short == "TIMEOUT":
                marker = "T"
            elif status_short == "COLLECTION_ERROR":
                marker = "C"
            else:
                marker = "F"
            sys.stdout.write(marker)
            sys.stdout.flush()
            if (i + 1) % 60 == 0:
                sys.stdout.write(f" [{i+1}/{len(test_files)}]\n")
                sys.stdout.flush()
    print()
    # Sort results
    results.sort(key=lambda r: r["file"])
    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULTS_PATH.write_text(json.dumps(results, indent=2))
    # Summary
    broken = [r for r in results if r["status"] != "GREEN"]
    print(f"\n--- SUMMARY ---")
    print(f"Total files: {len(results)}")
    print(f"GREEN: {sum(1 for r in results if r['status'] == 'GREEN')}")
    print(f"BROKEN: {sum(1 for r in results if r['status'] == 'BROKEN')}")
    print(f"COLLECTION_ERROR: {sum(1 for r in results if r['status'] == 'COLLECTION_ERROR')}")
    print(f"TIMEOUT: {sum(1 for r in results if r['status'] == 'TIMEOUT')}")
    print(f"\nFailing files:")
    for r in broken:
        print(f"  {r['status']:18} {r['file']:60} F={r['failed']} E={r['errors']} P={r['passed']} S={r['skipped']}")


if __name__ == "__main__":
    main()
