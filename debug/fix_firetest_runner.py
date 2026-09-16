"""Case J plants in the TEST file (d_exp lives there), not the module, but the
runner passed --plant-in MOD for every case.  Give each case an explicit
target.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P = "debug/firetest_p12_azimuthal.py"

with io.open(P, encoding="utf-8") as fh:
    s = fh.read()

s = s.replace(
    '    ("J  Gaussian control vs. removal of the d shell",\n'
    '     ["--plant-in-test", "    d_exp = [0.55, 1.60]=>    d_exp = []"],\n'
    '     "gaussian", True),',
    '    ("J  Gaussian control vs. removal of the d shell",\n'
    '     ["--plant", "    d_exp = [0.55, 1.60]=>    d_exp = []"],\n'
    '     "gaussian", True, TEST),', 1)

# every other case targets the module; give them an explicit target too
s = s.replace(
    'def main():\n'
    '    bad = 0\n'
    '    for label, plants, selector, want_fire in CASES:\n'
    '        cmd = ([sys.executable, "debug/qa/fire_test.py", TEST,\n'
    '                "--plant-in", MOD] + plants + ["-k", selector])',
    'def main():\n'
    '    bad = 0\n'
    '    for case in CASES:\n'
    '        label, plants, selector, want_fire = case[:4]\n'
    '        target = case[4] if len(case) > 4 else MOD\n'
    '        cmd = ([sys.executable, "debug/qa/fire_test.py", TEST,\n'
    '                "--plant-in", target] + plants + ["-k", selector])', 1)

with io.open(P, "w", encoding="utf-8") as fh:
    fh.write(s)

print("runner now supports a per-case plant target")
