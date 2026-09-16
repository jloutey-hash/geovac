r"""Round-3 code-review MATERIAL-SMALL and NIT items.

NIT-cond   test_basis_is_linearly_dependent_at_the_largest_truncation runs at
           alpha = 1.0 with floors 1e13 / 1e15, but the registry declares the
           2.6e14 and 2.0e16 literals at alpha = 1.05 (at 1.0 the values are
           3.04e14 / 2.89e16).  The test therefore cannot pin the literals it
           is credited with, and its floors sit about 1.5 decades low -- a
           one-sided guard that would survive the conditioning improving by a
           factor of thirty.  Moved to the declared convention and bounded
           both sides.

SMALL-matrix  Two tests in the file appear in no claim-matrix row, and the
           unregistered one turned out to be the false positive.  Added, along
           with the corrected Gaussian row.  Row 272 also still said "the same
           algebraic V_ee", which the paper stopped saying this round.

SMALL-budget  test_level4_multichannel.py justifies leaving Paper 15's 96.0%
           headline untested as "beyond a tractable CI budget" at ~754 s.  This
           same round promoted an ~850 s @slow Gaussian control into tests/ and
           called that cost justified.  The justification is not tenable as
           written;  corrected to say what is actually true (nobody has written
           it) and the coverage gap is raised to the PI rather than papered
           over with a budget argument.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

TAZ = "tests/test_paper12_azimuthal_channels.py"
TL4 = "tests/test_level4_multichannel.py"
MX = "docs/claim_test_matrix.md"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# ---- cond(S): declared convention, bounded both sides ----------------------
edit(TAZ,
     '''    Paper 12 states cond(S) = 2.6e14 at its headline N = 72 basis, rising past
    double precision once mu = 1 doubles it.  If this ever reads as
    well-conditioned, either the basis changed or the overlap is being built
    wrongly -- and the paper's caution about its own sixth decimal would be
    unfounded.
    """
    alpha = 1.0
    for (j, l, mu, floor) in ((3, 3, 0, 1e13), (3, 3, 1, 1e15)):
        basis = generate_basis(j, l, mu, alpha)
        mom = Moments(2.0 * alpha, 6 * max(j, l) + 6 * (mu + 2) + 20)
        s, _h = one_body(basis, R_DEFAULT, 1.0, mom)
        w = np.linalg.eigvalsh(s)
        cond = w[-1] / w[0]
        assert cond > floor, (
            f"cond(S) = {cond:.2e} at (j,l,mu) = ({j},{l},{mu}), below the "
            f"{floor:.0e} the paper's conditioning caution rests on"
        )''',
     '''    Paper 12 states cond(S) = 2.6e14 at its headline N = 72 basis, rising past
    double precision once mu = 1 doubles it.  If this ever reads as
    well-conditioned, either the basis changed or the overlap is being built
    wrongly -- and the paper's caution about its own sixth decimal would be
    unfounded.

    alpha = 1.05 is the registry's declared convention for these two literals
    (numeric_registry.py: cond_s_33) -- the alpha at which 2.6e14, 2.0e16 and
    the -79 Ha direct-eigensolve figure all land together.  An earlier version
    of this test ran at alpha = 1.0, where the values are 3.04e14 / 2.89e16, so
    it could not pin the numbers it was credited with;  its floors were also
    one-sided and about 1.5 decades low, which would have tolerated the
    conditioning improving thirtyfold without failing.
    """
    alpha = 1.05
    for (j, l, mu, lo, hi) in ((3, 3, 0, 1.0e14, 1.0e15),
                               (3, 3, 1, 5.0e15, 1.0e17)):
        basis = generate_basis(j, l, mu, alpha)
        mom = Moments(2.0 * alpha, 6 * max(j, l) + 6 * (mu + 2) + 20)
        s, _h = one_body(basis, R_DEFAULT, 1.0, mom)
        w = np.linalg.eigvalsh(s)
        cond = w[-1] / w[0]
        assert lo < cond < hi, (
            f"cond(S) = {cond:.2e} at (j,l,mu) = ({j},{l},{mu}), alpha = "
            f"{alpha}, outside [{lo:.0e}, {hi:.0e}].  The paper's caution "
            f"about its own sixth decimal rests on this magnitude; a value "
            f"below the band means the basis or the overlap changed, one "
            f"above means the solve is further past double precision than "
            f"the paper admits"
        )''',
     "NIT-cond: declared alpha, two-sided band")

# ---- the budget justification that this round's own precedent refutes -------
edit(TL4,
     """        NO-TEST (documented, not faked): the converged 96.0% / l_max=6 / 61-ch
        2D+cusp value is intentionally NOT asserted -- a single l_max=6 2D
        solve is ~754 s (Paper 15 Sec. VII.E), beyond a tractable CI budget.
        This l_max=4 case validates the identical machinery at a feasible cost.""",
     """        NO-TEST (documented, not faked): the converged 96.0% / l_max=6 / 61-ch
        2D+cusp value is NOT asserted anywhere.  This l_max=4 case validates
        the identical machinery at a lower cost, but it does not back the
        headline.

        The earlier justification -- that a single l_max=6 2D solve is ~754 s
        (Paper 15 Sec. VII.E) and therefore "beyond a tractable CI budget" --
        is withdrawn (2026-09-14).  The 2026-09-14 delta run promoted an
        ~850 s @slow Gaussian control into tests/ and judged that cost
        justified for a claim of the same standing, so 754 s cannot be the
        reason.  The real reason is that the test has not been written.  The
        gap is raised to the PI: 96.0% is a headline at eight loci (this
        paper's abstract, its summary table, and CLAUDE.md's best-results
        table) with no backing test and no inline tier tag.""",
     "SMALL-budget: the untenable budget justification withdrawn")

# ---- claim-matrix rows -----------------------------------------------------
edit(MX,
     """| 12 | **azimuthal channels close the gap**: sigma-only basis spans only m1=m2=0; restoring |m|<=1 in the same basis with the same algebraic V_ee gives 99.1% of D_e (+11.64 mHa) vs 0.34 mHa from tripling the sigma basis |""",
     """| 12 | **azimuthal channels close the gap**: sigma-only basis spans only m1=m2=0; restoring |m|<=1 in the same basis with the same Neumann kernel gives 99.1% of D_e (+11.64 mHa) vs 0.34 mHa from near-tripling the sigma basis (the quadrature-free property holds in the sigma sector only; mu>0 uses spectral quadrature) |""",
     "matrix: 'algebraic V_ee' -> 'Neumann kernel', quadrature scope added")

edit(MX,
     """| 13 | He 2D-var 0.022% raw (l=7) / 0.004% cusp (l=4) |""",
     """| 12 | **independent Gaussian route**: 8s3p2d two-centre FCI reproduces 92.34% (m=0) and 99.10% (|m|<=1), agreeing with the prolate values across an unrelated basis | `test_paper12_azimuthal_channels.py`: `test_independent_gaussian_route_agrees` (slow) | **BACKED-SOUND (repaired 2026-09-14)** | was a FALSE POSITIVE: it certified 99.42%, the ALL-m value, because adding the d shell turned `fci([True]*len(orbs))` into the |m|=2 calculation. Now builds the true m=0 (dim 30) and |m|<=1 (dim 50) contractions and bounds both sides |
| 12 | **the general-m code path reproduces the original sigma pipeline** at every basis size (161/1.6/1.0/6.7/58 uHa) | `test_paper12_azimuthal_channels.py`: `test_mu0_reproduces_exact_neumann_machinery` | BACKED-SOUND | registered 2026-09-14; previously in no matrix row |
| 13 | He 2D-var 0.022% raw (l=7) / 0.004% cusp (l=4) |""",
     "matrix: the two unregistered tests now have rows")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

applied, failed = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    for old, new, label in items:
        if old in t:
            t = t.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append("%s   [%s]" % (label, path))
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)

print("applied %d of %d" % (len(applied), len(EDITS)))
for a in applied:
    print("  +", a)
if failed:
    print("")
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
