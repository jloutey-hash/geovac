r"""Fold the round-3 code/test-backing review into the v5.11.20 CHANGELOG entry.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

CL = "CHANGELOG.md"

ANCHOR = """### Guards

Eight retired claims registered in C16"""

INSERT = """### The code/test-backing reviewer found two more, and one of them was mine from the same day

- **The independent Gaussian control was certifying the wrong number.** `test_independent_gaussian_route_agrees` computed the *all-m* value, 99.42%, not the 99.10% the paper quotes. With the old 8s3p basis every function had |m| <= 1, so taking every orbital WAS the |m| <= 1 calculation; round 2 added the d shell to fix a different defect and thereby turned that same line into the |m| = 2 calculation. Neither assertion could tell 99.10 from 99.42, and the gain assertion cited "the paper quotes 12.36" -- a number in no paper and no registry, whose only home was a `debug/` memo. **The paper's numbers are right**: 92.34 and 99.10 both reproduce exactly once the true m = 0 (dim 30) and |m| <= 1 (dim 50) spaces are built. The test was wrong, and it is the third repetition of the wrong-evaluation-object class in this one test. Those spaces are contractions, not index masks -- xx+yy is m = 0 while xx-yy is |m| = 2 -- so the control now builds them explicitly, bounds both sides, and asserts that the full space sits measurably *above* the |m| <= 1 space, which is the discriminator the old test lacked.
- **"The weakest variational value anywhere on the full grid is 95.5%" is refuted, and it is a sentence written earlier the same day.** Re-measured independently: 95.5% is the cell at alpha = 1.10, threshold 1e-8, and is the grid minimum **only if alpha is restricted to >= 1.10** -- excluding exactly the low-alpha half the same paragraph had just called pathological. Over the declared range alpha in [0.90, 1.30] the weakest variational value is **76.11%**, and at that cell the sigma-only solve returns **91.06%**, so opening the azimuthal channels there makes the answer *worse by fifteen points*. "Robust across the grid" is false at a grid point. A restricted-evaluation artifact: a clean floor produced by deleting the part of the object that breaks it.

  The paragraph now says what was measured -- the value ranges 76.1-99.1% over the declared grid -- and states what the result actually rests on: the tight-threshold points, which all lie in 99.03-99.14%, and the well-conditioned Gaussian route at 99.10%, which needs no discard threshold at all. Not a property of the grid, most of which is conditioning noise.

**This is the fourth consecutive round whose largest defect sat in the previous round's remediation, and the first in which the previous round was the same session.** The withdraw-do-not-replace rule stopped the narrative defects; it did not stop an unverified *number* being carried forward into a correction. The number was never re-measured, only re-framed.

### Also from the code review

- The `[92.0, 92.6]` band written this round to stop hiding the dropped-digit denominator **still admitted it** -- under the wrong denominator the value reads 92.1087, inside the band, and a plant reverting the literal alone passed. Narrowed to `[92.15, 92.35]`.
- The sigma band `(92.2, 92.5)` did not reject the d-less basis it claimed to: 8s3p gives 92.22, which is inside it. Only the other leg was doing the work.
- The `cond(S)` test ran at alpha = 1.0 while the registry declares its two literals at alpha = 1.05, so it could not pin the numbers it was credited with, and its one-sided floors sat about 1.5 decades low. Moved to the declared convention and bounded both sides.
- Two tests in the file were in no claim-matrix row -- and the unregistered one was the false positive. Both now registered.
- Two new fire cases. **K** loosens the discard *behaviour* by two decades while leaving the signature literal untouched: the `inspect.signature` pin is blind to it and only the envelope bound catches it, which is what makes the envelope load-bearing rather than decorative. **L** folds |m| = 2 back into the Gaussian control and must be rejected.

### On the test budget

`test_level4_multichannel.py` justified leaving Paper 15's 96.0% headline untested as "beyond a tractable CI budget" at ~754 s. This same round promoted an ~850 s test into `tests/` and called that cost justified. The justification is withdrawn; the real reason is that nobody has written it. Raised as PI item 4 below.

### Guards

Eight retired claims registered in C16"""

with io.open(CL, encoding="utf-8") as fh:
    t = fh.read()

if ANCHOR not in t:
    print("FAILED: anchor not found")
    sys.exit(1)

with io.open(CL, "w", encoding="utf-8") as fh:
    fh.write(t.replace(ANCHOR, INSERT, 1))

# the registry count moved from 8 to 9
t2 = io.open(CL, encoding="utf-8").read()
t2 = t2.replace(
    "Eight retired claims registered in C16 with declared dependents, and each "
    "proved to discriminate in **both** directions before being trusted -- 11 cases",
    "Nine retired claims registered in C16 with declared dependents, and each "
    "proved to discriminate in **both** directions before being trusted -- 13 cases",
    1)
io.open(CL, "w", encoding="utf-8").write(t2)
print("  + CHANGELOG v5.11.20 gains the code-review section")
