"""DELTA #2 code findings LARGE-1 and LARGE-2.

LARGE-1  The control test I promoted uses 8s3p -- no d shell -- and therefore
         computes 92.22 / 98.49, while the paper quotes 92.34 / 99.10, which
         are the 8s3p2d numbers.  So the test gave a permanent home to a
         DIFFERENT calculation, and the paper's quoted control is still backed
         only by debug/ -- the exact condition the promotion was meant to
         remove.  I chose the smaller basis for runtime and then let the
         docstring claim it validated the paper's figures.
         Fix: add the d shell so the test computes what the paper prints.

LARGE-2  The headline guard asserts only `> 98.9` with no upper bound, so
         loosening the discard threshold by a decade in either direction
         leaves it green while the certified value moves to 98.99% or 99.20%
         -- outside the paper's own stated envelope.
         Fix: bound the value inside the declared envelope, and pin the
         threshold the envelope is stated at.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

T = "tests/test_paper12_azimuthal_channels.py"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ---------------- LARGE-1: give the test the basis the paper quotes
edit(
    """    s_exp = [0.0347, 0.0925, 0.2469, 0.6584, 1.7557, 4.6819, 12.485, 33.293]
    p_exp = [0.25, 0.75, 2.25]

    orbs, is_sigma = [], []
    for c in centers:
        for a in s_exp:
            orbs.append(BasisFn(c, (0, 0, 0), np.array([a]), np.array([1.0])))
            is_sigma.append(True)
        for a in p_exp:
            for lmn, sig in (((0, 0, 1), True), ((1, 0, 0), False),
                             ((0, 1, 0), False)):
                orbs.append(BasisFn(c, lmn, np.array([a]), np.array([1.0])))
                is_sigma.append(sig)""",
    """    # 8s3p2d -- the basis the PAPER quotes.  An earlier version of this test
    # used 8s3p, which computes 92.22 / 98.49 and therefore validated a
    # different calculation from the 92.34 / 99.10 the paper prints.  The d
    # shell costs runtime (this test is @slow for that reason) and buys the
    # only thing that matters here: that the test pins the published control.
    s_exp = [0.0347, 0.0925, 0.2469, 0.6584, 1.7557, 4.6819, 12.485, 33.293]
    p_exp = [0.25, 0.75, 2.25]
    d_exp = [0.55, 1.60]

    orbs, is_sigma = [], []
    for c in centers:
        for a in s_exp:
            orbs.append(BasisFn(c, (0, 0, 0), np.array([a]), np.array([1.0])))
            is_sigma.append(True)
        for a in p_exp:
            for lmn, sig in (((0, 0, 1), True), ((1, 0, 0), False),
                             ((0, 1, 0), False)):
                orbs.append(BasisFn(c, lmn, np.array([a]), np.array([1.0])))
                is_sigma.append(sig)
        for a in d_exp:
            # zz and xx+yy are m = 0; xz, yz are |m| = 1; xy, xx-yy are |m| = 2.
            # Cartesian xx and yy each mix m = 0 with |m| = 2, so neither is
            # sigma-pure; they are excluded from the sigma set, which makes the
            # sigma restriction conservative (it can only UNDER-state the
            # sigma ceiling, never inflate it).
            for lmn, sig in (((0, 0, 2), True), ((2, 0, 0), False),
                             ((0, 2, 0), False), ((1, 0, 1), False),
                             ((0, 1, 1), False), ((1, 1, 0), False)):
                orbs.append(BasisFn(c, lmn, np.array([a]), np.array([1.0])))
                is_sigma.append(sig)""",
    "LARGE-1: control test now uses the basis the paper quotes")

edit(
    """    assert 91.5 < pct_sigma < 93.0, (
        f"Gaussian sigma-only ceiling is {pct_sigma:.2f}%, but Paper 12's "
        f"prolate sigma-only value is 92.4% -- two unrelated bases should "
        f"agree here, and that agreement is the whole point of this control"
    )
    assert pct_all > 98.0, (
        f"releasing |m| = 1 reaches only {pct_all:.2f}%; the gap does not "
        f"close in the independent basis"
    )
    assert 1000.0 * (e_sigma - e_all) > 9.0, (
        f"the azimuthal channels are worth only "
        f"{1000.0*(e_sigma-e_all):.2f} mHa in the Gaussian basis"
    )""",
    """    # The paper prints 92.34 and 99.10 for this basis; bands are tight enough
    # that substituting a smaller basis (8s3p gives 92.22 / 98.49) fails.
    assert 92.2 < pct_sigma < 92.5, (
        f"Gaussian sigma-only ceiling is {pct_sigma:.2f}%; the paper quotes "
        f"92.34%, within 0.2 mHa of its prolate sigma-only value, and that "
        f"agreement across unrelated bases is the whole point of this control"
    )
    assert pct_all > 98.9, (
        f"releasing the azimuthal channels reaches only {pct_all:.2f}%; the "
        f"paper quotes 99.10% for this basis"
    )
    assert 1000.0 * (e_sigma - e_all) > 11.5, (
        f"the azimuthal channels are worth only "
        f"{1000.0*(e_sigma-e_all):.2f} mHa here; the paper quotes 12.36"
    )""",
    "LARGE-1b: bands tightened to the published numbers")

# ---------------- LARGE-2: bound the envelope and pin the threshold
edit(
    """    assert _de_pct(e_pi) > 98.9, (
        f"headline is {_de_pct(e_pi):.2f}% of D_e, paper quotes 99.1% "
        f"(envelope 99.0-99.1%)"
    )""",
    """    # Bounded on BOTH sides, and the threshold pinned.  With only a lower
    # bound, loosening the discard threshold one decade either way left this
    # green while the certified value moved to 98.99% or 99.20% -- outside the
    # envelope the paper states.  A guard for an envelope has to be an
    # envelope.
    import inspect
    default_thresh = inspect.signature(solve_generalized).parameters["thresh"].default
    assert default_thresh == 1e-11, (
        f"the discard threshold is {default_thresh:g}; the paper's stability "
        f"envelope (99.0-99.1%) is stated at 1e-11, so moving the default "
        f"invalidates the published envelope rather than just the number"
    )
    assert 99.0 < _de_pct(e_pi) < 99.2, (
        f"headline is {_de_pct(e_pi):.2f}% of D_e, outside the paper's stated "
        f"stability envelope of 99.0-99.1%"
    )""",
    "LARGE-2: envelope bounded both sides; threshold pinned")

with io.open(T, encoding="utf-8") as fh:
    t = fh.read()

applied, failed = [], []
for old, new, label in EDITS:
    if old in t:
        t = t.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

with io.open(T, "w", encoding="utf-8") as fh:
    fh.write(t)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
