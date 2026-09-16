"""Fire-test: prove the `p12-cusp-as-the-h2-gap` registry entry can FAIL.

The entry reported "clean (exempt/withdrawn-flagged: 0)" on the corrected
corpus -- which is the right answer, but is indistinguishable from a pattern
that matches nothing ever.  Sec. 9: a guard whose rejected answer cannot be
named is not a guard.  So name it: plant each withdrawn formulation, require
the gate to FAIL, remove the plant, require the gate to PASS again.

Plants are removed by exact-string deletion, not by git, because the files
carry uncommitted work.

Run:  python debug/firetest_p12_retracted.py
"""
import io
import subprocess
import sys

GATE = [sys.executable, "debug/qa/check_retracted_terms.py"]

PLANTS = [
    ("papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
     "\n%% FIRE-TEST PLANT: The electron-electron cusp requires the "
     "non-analytic terms that no polynomial prolate spheroidal basis can "
     "represent.\n",
     "cusp requires the non-analytic"),
    ("papers/synthesis/group2_quantum_chemistry_synthesis.tex",
     "\n%% FIRE-TEST PLANT: This exceeds the earlier figure, confirming the "
     "cusp-resolution advantage of the hyperspherical geometry.\n",
     "confirming the cusp-resolution advantage"),
    ("papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
     "\n%% FIRE-TEST PLANT: The diagnosis identifies the next natural "
     "geometry for the coalescence problem.\n",
     "diagnosis identifies the next natural geometry"),
]


def gate_status():
    r = subprocess.run(GATE, capture_output=True, text=True)
    tail = [ln for ln in r.stdout.splitlines() if ln.startswith("RESULT")]
    return r.returncode, (tail[-1] if tail else "(no RESULT line)")


def entry_line():
    r = subprocess.run(GATE, capture_output=True, text=True)
    for ln in r.stdout.splitlines():
        if "p12-cusp-as-the-h2-gap" in ln:
            return ln.strip()
    return "(entry not listed!)"


def main():
    rc0, res0 = gate_status()
    print("baseline : rc=%d  %s" % (rc0, res0))
    print("           %s" % entry_line())
    if rc0 != 0:
        print("ABORT: baseline is already failing; fix that first.")
        return 1

    failures = 0
    for path, plant, label in PLANTS:
        with io.open(path, encoding="utf-8") as fh:
            original = fh.read()
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write(original + plant)
        rc, res = gate_status()
        line = entry_line()
        with io.open(path, encoding="utf-8") as fh:
            now = fh.read()
        assert plant in now, "plant vanished"
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write(now.replace(plant, "", 1))

        fired = rc != 0
        print("plant    : %-45s -> %s" % (label, "FIRED" if fired else "*** DID NOT FIRE ***"))
        print("           %s" % line)
        if not fired:
            failures += 1

    rc1, res1 = gate_status()
    print("restored : rc=%d  %s" % (rc1, res1))
    if rc1 != 0:
        print("*** plants were not fully removed ***")
        failures += 1

    print()
    if failures:
        print("FIRE TEST FAILED: %d plant(s) did not fire." % failures)
        return 1
    print("FIRE TEST PASSED: every withdrawn formulation is caught, and the "
          "corrected corpus is clean.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
