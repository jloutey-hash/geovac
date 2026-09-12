"""Repair test 7: the bracket-direction claim does NOT transfer between sectors.

Written first on the cheap s-only ladder, on the assumption that "the windowed
fit approaches from below" is a property of the extrapolator.  It is not.  On
the s-only ladder the fitted floors FALL:

    4.3098 -> 4.3059 -> 4.3035

so the test failed, correctly.  On the spdf ground-state ladder they RISE, and
reproduce the stored driver output exactly:

    6.3887 -> 6.4192 -> 6.4388      Shanks(last 3) = 6.7165, above all of them

That is a real fact worth keeping rather than a nuisance:  the bracket
construction the paper quotes is specific to the spdf ground-state ladder, not a
general property of the two extrapolators.  A cheap proxy in the wrong sector
would have "backed" the claim while testing something that behaves oppositely --
the exact failure mode a guard is supposed to exclude, reached by writing one.

Cost of the correct version: the spdf ladder K=74..244, 266s measured.  The two
endpoint VALUES [6.47, 6.62] still need K=452 (+7 min) and stay driver-backed.
"""
import io

P = "tests/test_paper60_resource_ladder.py"
s = io.open(P, encoding="utf-8").read()

# --- module docstring: it claims the s-only ladder is used for the bracket.
OLD_DOC = """Partially closed, and declared rather than papered over:  the floor BRACKET
``[6.47, 6.62]`` / ``[1.647, 1.676]``.  Its two endpoints come from the full
spdf ladder K=74..452, which costs ~20 minutes to rebuild -- more than a test
should.  What is backed here is the CLAIM FORM that makes it a bracket at all:
that the windowed three-parameter fit approaches the floor FROM BELOW while a
model-free Shanks extrapolation descends FROM ABOVE, so the two straddle.  That
is tested on the s-only ladder, which is cheap.  The two spdf endpoint VALUES
remain driver-backed and are recorded as such in docs/claim_test_matrix.md."""
NEW_DOC = """Partially closed, and declared rather than papered over:  the floor BRACKET
``[6.47, 6.62]`` / ``[1.647, 1.676]``.  Its two endpoints need the full spdf
ladder K=74..452, ~20 minutes to rebuild -- more than a test should cost.  What
is backed here is the CLAIM FORM that makes it a bracket at all: that the
windowed three-parameter fit approaches the floor FROM BELOW while a model-free
Shanks extrapolation descends FROM ABOVE, so the two straddle.  Tested on the
spdf ladder truncated at K=244 (266s).

That truncation is the cheapest HONEST version.  This test was first written on
the s-only ladder, which is far cheaper -- and it failed, because on that ladder
the fitted floors FALL (4.3098, 4.3059, 4.3035) instead of rising.  The approach
direction is a property of the SECTOR, not of the extrapolator, so the bracket
construction the paper quotes is specific to the spdf ground-state ladder.  A
cheap proxy in the wrong sector would have reported the claim "backed" while
testing something that behaves oppositely.  The two spdf endpoint VALUES remain
driver-backed and are recorded as such in docs/claim_test_matrix.md."""
assert OLD_DOC in s
s = s.replace(OLD_DOC, NEW_DOC, 1)

# --- replace the test body.
i = s.index("@pytest.mark.slow\ndef test_floor_bracket_directions_are_opposite():")
NEW_TEST = '''@pytest.mark.slow
def test_floor_bracket_directions_are_opposite_on_the_spdf_ladder():
    """Why the two extrapolators BRACKET: they approach from opposite sides.

    The paper quotes the ground-state floor as a bracket rather than a single
    value, on the grounds that the windowed three-parameter fit
    ``dE(K) = c + b K^-q`` approaches from BELOW (c rises as the window moves
    out) while a model-free Shanks extrapolation descends from ABOVE.

    Measured here on the spdf ground-state ladder K = 74..244:

        windowed fitted floors : 6.3887 -> 6.4192 -> 6.4388   (rising)
        Shanks on the last 3   : 6.7165                       (above all)

    reproducing the driver's stored values exactly.

    SECTOR-SPECIFICITY, which is why this test costs 266s instead of 60s.  The
    direction is NOT a property of the extrapolator.  On the s-only ladder the
    same windowed fit FALLS (4.3098, 4.3059, 4.3035), so a cheaper s-only
    version of this test fails -- and would have been worse than useless had it
    passed, since it would have certified the claim while measuring a sector
    that behaves oppositely.

    SCOPE.  The two endpoint VALUES [6.47, 6.62] need the full ladder out to
    K=452 (a further ~7 minutes) and are NOT asserted here;  they remain
    driver-backed, recorded in docs/claim_test_matrix.md.  What is backed is the
    claim form that makes "bracket" the right word.

    WRONG ANSWER THIS EXCLUDES.  "The fit and Shanks are two estimates of the
    same quantity, so average them" -- sound only if they straddle.  The
    assertions below establish that they do, in the stated directions, and the
    strict inequality between max(fit) and Shanks rules out their coinciding.
    """
    K, G = [], []
    for nmax in (7, 8, 9, 10, 11, 12):
        _grid_for(nmax)
        tuples = SV.family(nmax, 3)
        e_iso, _m1, _p, _M = SS.solve(tuples, Z=Z)
        K.append(len(tuples))
        G.append((e_iso - EXACT) * 1000.0)
    assert K == [74, 100, 130, 164, 202, 244], K
    K = np.array(K, float)
    G = np.array(G, float)

    def fit_floor(k, y):
        best = None
        for q in np.linspace(0.05, 8.0, 4000):
            A = np.column_stack([np.ones_like(k), k ** (-q)])
            sol, *_ = np.linalg.lstsq(A, y, rcond=None)
            r = float(np.sqrt(np.mean((y - A @ sol) ** 2)))
            if best is None or r < best[0]:
                best = (r, float(sol[0]))
        return best[1]

    def shanks(y):
        a, b, c = y[-3], y[-2], y[-1]
        return float(c - (c - b) ** 2 / ((c - b) - (b - a)))

    c4 = [fit_floor(K[i:i + 4], G[i:i + 4]) for i in range(len(K) - 3)]
    # the fit approaches FROM BELOW
    assert all(c4[i] < c4[i + 1] for i in range(len(c4) - 1)), c4
    for got, want in zip(c4, (6.3887, 6.4192, 6.4388)):
        assert abs(got - want) < 2e-3, (got, want)

    sh = shanks(G)
    assert abs(sh - 6.7165) < 2e-3, sh
    # ...Shanks descends from ABOVE: strictly above every windowed fit
    assert sh > max(c4) + 0.1, (sh, max(c4))
    # the last measured ladder point also sits above the fitted floor
    assert G[-1] > max(c4), (G[-1], max(c4))
'''
s = s[:i] + NEW_TEST
io.open(P, "w", encoding="utf-8").write(s)
print("test 7 rewritten on the spdf ladder; module docstring corrected")
