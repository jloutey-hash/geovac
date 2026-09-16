"""Fire tests for tests/test_paper12_azimuthal_channels.py.

Every guard in that file, run against the specific wrong answer it claims to
exclude.  Recorded here so the checks are reproducible rather than a one-off
transcript, per the Sec. 9 rule that a guard whose rejected answer cannot be
named is not a guard.

Result on 2026-09-14 (all as recorded in the test docstrings):

  A  kill the m != 0 azimuthal coupling        -> FIRED
  B  revert to Paper 12's printed prefactor    -> FIRED
  C  remove the eta selection rule ALONE       -> did not fire  (expected)
  D  remove the l cap ALONE                    -> did not fire  (expected)
  E  remove BOTH                               -> FIRED

C and D are the informative ones:  the two exactness mechanisms are
deliberately redundant, so no single-point mutation breaks the property and
the guard is right to stay green.  The docstring in the test was corrected to
claim only the conjunction after this run.

Run:  python debug/firetest_p12_azimuthal.py
"""
import subprocess
import sys

TEST = "tests/test_paper12_azimuthal_channels.py"
MOD = "geovac/prolate_general_m.py"

ETA_RULE = ("                    if l > Qq + 2 * s - m or (Qq + l - m) % 2 != 0:"
            "=>                    if False:")
L_CAP = ("    l_neumann = min(l_neumann, q_max + 2 * max(s_set))"
         "=>    l_neumann = l_neumann")

CASES = [
    ("A  kill the m != 0 azimuthal coupling",
     ["--plant", "            tot += math.pi=>            tot += 0.0"],
     "close_the_gap or channel_effect", True),
    ("B  revert to Paper 12's printed prefactor",
     ["--plant",
      "    return (-1.0) ** m * (2 * l + 1) * r * r=>    return r"],
     "kernel_reproduces_coulomb", True),
    ("C  remove the eta selection rule ALONE",
     ["--plant", ETA_RULE], "truncation_is_exact", False),
    ("D  remove the l cap ALONE",
     ["--plant", L_CAP], "truncation_is_exact", False),
    ("E  remove BOTH exactness mechanisms",
     ["--plant", L_CAP, "--plant", ETA_RULE],
     "truncation_is_exact", True),
    # Added 2026-09-14 after the /qa DELTA run found the 99.09% HEADLINE
    # unpinned (the fast test pins the (2,2) row).  F and G prove the new
    # @slow headline guard covers both of the things its docstring claims.
    ("F  headline guard vs. killed azimuthal coupling",
     ["--plant", "            tot += math.pi=>            tot += 0.0"],
     "headline", True),
    ("G  headline guard vs. disabled canonical orthogonalisation",
     ["--plant", "    keep = w > thresh * w[-1]=>    keep = w > -1e300"],
     "headline", True),
    # Case H closes the one real gap the 2026-09-14 DELTA code review found:
    # its plant P10 (tighten the discard threshold) left all seven tests green
    # while the headline configuration returned -4.10 Ha.  The alpha-robust
    # headline guard now catches it.
    ("H  headline guard vs. tightened discard threshold",
     ["--plant",
      "def solve_generalized(H, S, thresh: float = 1e-11):"
      "=>def solve_generalized(H, S, thresh: float = 1e-15):"],
     "headline", True),
    # DELTA #2 closed two more gaps its code review found: the headline guard
    # had no UPPER bound (loosening the threshold a decade left it green while
    # the certified value left the paper's envelope), and the Gaussian control
    # measured a d-less basis, i.e. a different calculation from the one the
    # paper quotes.
    ("I  headline guard vs. LOOSENED discard threshold",
     ["--plant",
      "def solve_generalized(H, S, thresh: float = 1e-11):"
      "=>def solve_generalized(H, S, thresh: float = 1e-10):"],
     "headline", True),
    ("J  Gaussian control vs. removal of the d shell",
     ["--plant", "    d_exp = [0.55, 1.60]=>    d_exp = []"],
     "gaussian", True, TEST),
    # DELTA #3.  Cases H and I plant on the SIGNATURE LITERAL, so they fire
    # through the inspect.signature pin and say nothing about whether the
    # envelope bound does any work.  K loosens the discard BEHAVIOUR by two
    # decades while leaving the signature at 1e-11: the pin sees nothing, and
    # only the 99.0 < pct < 99.2 envelope can catch it.  This is the case that
    # shows the envelope assertion is load-bearing rather than decorative.
    ("K  headline envelope vs. loosened discard BEHAVIOUR (signature intact)",
     ["--plant",
      "    keep = w > thresh * w[-1]=>    keep = w > 1e2 * thresh * w[-1]"],
     "headline", True),
    # DELTA #3.  The Gaussian control certified 99.42% (the ALL-m value)
    # instead of the paper's 99.10%, because adding the d shell turned
    # fci([True]*len(orbs)) into the |m| = 2 calculation.  L plants the old
    # behaviour back -- folding the delta sector into the |m| <= 1 space --
    # and the new two-sided bound must reject it.
    ("L  Gaussian control vs. folding |m| = 2 into the control",
     ["--plant",
      "    m1_idx = [i for i in range(n) if i not in drop]"
      "=>    m1_idx = [i for i in range(n)]"],
     "gaussian", True, TEST),
]


def main():
    bad = 0
    for case in CASES:
        label, plants, selector, want_fire = case[:4]
        target = case[4] if len(case) > 4 else MOD
        cmd = ([sys.executable, "debug/qa/fire_test.py", TEST,
                "--plant-in", target] + plants + ["-k", selector])
        r = subprocess.run(cmd, capture_output=True, text=True)
        fired = r.returncode == 0          # fire_test exits 0 when the plant FAILS the test
        verdict = "FIRED" if fired else "did not fire"
        ok = (fired == want_fire)
        print("  %-42s -> %-13s %s"
              % (label, verdict, "ok" if ok else "*** UNEXPECTED ***"))
        if not ok:
            bad += 1
            print(r.stdout.strip()[-400:])
    print()
    if bad:
        print("FIRE TEST SUITE FAILED: %d case(s) behaved unexpectedly." % bad)
        return 1
    print("FIRE TEST SUITE PASSED: every guard discriminates as documented.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
