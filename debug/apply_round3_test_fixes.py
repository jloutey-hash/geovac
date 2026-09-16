r"""Round-3 code-review test fixes.

LARGE-1  test_independent_gaussian_route_agrees certified 99.42%, not the
         99.10% the paper quotes.  With 8s3p every function had |m| <= 1, so
         fci([True]*len(orbs)) WAS the |m| <= 1 calculation.  Round 2 added the
         d shell to fix a different defect and thereby made that same line the
         |m| = 2 calculation -- the third repetition of the wrong-evaluation-
         object class in this one test.  Neither assertion could tell 99.10
         from 99.42, and "the paper quotes 12.36" is a number that appears in
         no paper and in no registry;  its only home was a debug/ memo.
         The PAPER's numbers are right -- 92.34 and 99.10 both reproduce
         exactly once the right spaces are built.  The test was wrong.

         Cartesian d shell:  zz and xx+yy are m = 0;  xz, yz are |m| = 1;
         xy and xx-yy are |m| = 2.  So the true m = 0 space adds the xx+yy
         combination to the sigma set (dim 30), and the true |m| <= 1 space is
         everything except xy and xx-yy (dim 50).  Those are contractions, not
         index masks, which is why fci() now takes a coefficient matrix.

S10-bis  The [92.0, 92.6] band written on 2026-09-14 to stop hiding the
         dropped-digit denominator still admits it:  under the wrong
         E_exact the value reads 92.1087, which is inside the band.  The
         comment claimed the narrowing was what made the guard non-blind;  it
         was not.  Reviewer's plant P4 (revert the literal only) passed.

SIGMA    "bands are tight enough that substituting a smaller basis (8s3p gives
         92.22 / 98.49) fails" -- 92.22 is INSIDE (92.2, 92.5).  The plural
         was wrong;  only the pct_all leg rejected 8s3p.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

TAZ = "tests/test_paper12_azimuthal_channels.py"
TNV = "tests/test_neumann_vee.py"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# ---- track the indices the contractions need ------------------------------
edit(TAZ,
     """    orbs, is_sigma = [], []""",
     """    orbs, is_sigma = [], []
    # Indices needed to build the AZIMUTHAL spaces.  A boolean mask cannot
    # express them: xx+yy is m = 0 and xx-yy is |m| = 2, so the two sectors
    # are separated by a rotation of the (xx, yy) pair, not by selection.
    xx_yy_pairs, xy_idx = [], []""",
     "LARGE-1a: index bookkeeping for the azimuthal contractions")

edit(TAZ,
     """            for lmn, sig in (((0, 0, 2), True), ((2, 0, 0), False),
                             ((0, 2, 0), False), ((1, 0, 1), False),
                             ((0, 1, 1), False), ((1, 1, 0), False)):
                orbs.append(BasisFn(c, lmn, np.array([a]), np.array([1.0])))
                is_sigma.append(sig)""",
     """            here = {}
            for lmn, sig in (((0, 0, 2), True), ((2, 0, 0), False),
                             ((0, 2, 0), False), ((1, 0, 1), False),
                             ((0, 1, 1), False), ((1, 1, 0), False)):
                here[lmn] = len(orbs)
                orbs.append(BasisFn(c, lmn, np.array([a]), np.array([1.0])))
                is_sigma.append(sig)
            xx_yy_pairs.append((here[(2, 0, 0)], here[(0, 2, 0)]))
            xy_idx.append(here[(1, 1, 0)])""",
     "LARGE-1b: record the xx/yy pair and xy per d shell")

# ---- build the real spaces and assert against them -------------------------
edit(TAZ,
     """    def fci(mask):
        idx = [i for i, keep in enumerate(mask) if keep]
        c = np.eye(len(orbs))[:, idx]
        s2, h2 = c.T @ s @ c, c.T @ h @ c
        g2 = np.einsum("pi,qj,rk,sl,pqrs->ijkl", c, c, c, c, g, optimize=True)
        x = lowdin_orbitals(s2)
        ht, gt = transform_integrals(x, h2, g2)
        return fci_ground(ht, gt, 2) + e_nuc

    e_sigma = fci(is_sigma)
    e_all = fci([True] * len(orbs))

    pct_sigma = 100.0 * (-1.0 - e_sigma) / DE_EXACT
    pct_all = 100.0 * (-1.0 - e_all) / DE_EXACT

    assert e_all > E_EXACT and e_sigma > E_EXACT, "non-variational"
    # The paper prints 92.34 and 99.10 for this basis; bands are tight enough
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
     """    def fci(cmat):
        s2, h2 = cmat.T @ s @ cmat, cmat.T @ h @ cmat
        g2 = np.einsum("pi,qj,rk,sl,pqrs->ijkl", cmat, cmat, cmat, cmat, g,
                       optimize=True)
        x = lowdin_orbitals(s2)
        ht, gt = transform_integrals(x, h2, g2)
        return fci_ground(ht, gt, 2) + e_nuc

    n = len(orbs)
    eye = np.eye(n)

    # The m = 0 half of each Cartesian (xx, yy) pair.
    plus = np.zeros((n, len(xx_yy_pairs)))
    for k, (i_xx, i_yy) in enumerate(xx_yy_pairs):
        plus[i_xx, k] = plus[i_yy, k] = 1.0 / np.sqrt(2.0)

    # sigma: s, p_z, zz, and xx+yy.  Adding xx+yy is what makes this the TRUE
    # m = 0 space rather than the conservative subset -- it is worth 0.015 mHa,
    # and it is the difference between the paper's 92.34 and 92.33.
    sigma_idx = [i for i, sg in enumerate(is_sigma) if sg]
    c_sigma = np.hstack([eye[:, sigma_idx], plus])

    # |m| <= 1: everything except the |m| = 2 directions xy and xx-yy.  Dropping
    # the raw xx and yy columns removes BOTH, so xx+yy is added back.
    drop = set(xy_idx) | {i for pair in xx_yy_pairs for i in pair}
    m1_idx = [i for i in range(n) if i not in drop]
    c_m1 = np.hstack([eye[:, m1_idx], plus])

    e_sigma = fci(c_sigma)
    e_m1 = fci(c_m1)
    e_all = fci(eye)

    pct_sigma = 100.0 * (-1.0 - e_sigma) / DE_EXACT
    pct_m1 = 100.0 * (-1.0 - e_m1) / DE_EXACT
    pct_all = 100.0 * (-1.0 - e_all) / DE_EXACT

    assert c_sigma.shape[1] == 30 and c_m1.shape[1] == 50 and n == 58, (
        f"azimuthal spaces are the wrong size: sigma={c_sigma.shape[1]} "
        f"(expect 30), |m|<=1={c_m1.shape[1]} (expect 50), all={n} (expect 58)"
    )
    assert e_all > E_EXACT and e_m1 > E_EXACT and e_sigma > E_EXACT, \\
        "non-variational"

    # Bounded BOTH sides.  An earlier version asserted only `pct_sigma > 92.2`
    # and `pct_all > 98.9`, which could not reject 8s3p's 92.22 on the sigma
    # leg, and could not tell the |m| <= 1 value (99.10) from the all-m value
    # (99.42) on the other.
    assert 92.30 < pct_sigma < 92.40, (
        f"Gaussian sigma-only ceiling is {pct_sigma:.2f}%; the paper quotes "
        f"92.34%, within 0.2 mHa of its prolate sigma-only value, and that "
        f"agreement across unrelated bases is the whole point of this control"
    )
    assert 99.0 < pct_m1 < 99.2, (
        f"releasing the azimuthal channels reaches {pct_m1:.2f}%; the paper "
        f"quotes 99.10% for this basis, inside the 99.0-99.1 envelope"
    )
    # The discriminator the old test lacked: |m| = 2 is a DIFFERENT sector and
    # must not be silently folded into the control.  If this stops holding,
    # c_m1 has been built wrong.
    assert pct_all > pct_m1 + 0.15, (
        f"the full space ({pct_all:.2f}%) is not measurably above the "
        f"|m|<=1 space ({pct_m1:.2f}%), so the delta sector is being "
        f"counted inside the control -- exactly the LARGE-1 defect"
    )
    gain = 1000.0 * (e_sigma - e_m1)
    assert 11.5 < gain < 12.1, (
        f"the azimuthal channels are worth {gain:.2f} mHa here; the paper's "
        f"92.34 -> 99.10 implies 11.8 mHa, against 11.64 in the prolate basis"
    )""",
     "LARGE-1c: the control now measures |m|<=1, bounded both sides")

# ---- S10-bis: the band that still accepted the dropped digit ---------------
edit(TNV,
     """        # Band tightened 2026-09-14.  [90, 94] was wide enough to hide a
        # 0.16% error in D_e_exact's denominator for as long as it stood;
        # correcting the literal without narrowing the band would have left
        # the guard exactly as blind as before.
        assert 92.0 < pct < 92.6, \\
            f"H2 D_e fraction {pct:.2f}% outside headline band [92.0, 92.6]\"""",
     """        # Band tightened twice.  [90, 94] was wide enough to hide a 0.16%
        # error in D_e_exact's denominator for as long as it stood.  The first
        # narrowing, to [92.0, 92.6], did NOT fix that: under the wrong
        # denominator the value reads 92.1087, which is still inside it -- a
        # plant reverting the literal alone passed.  The correct value is
        # 92.2539, so the lower edge has to sit above 92.11 to exclude it.
        assert 92.15 < pct < 92.35, \\
            f"H2 D_e fraction {pct:.2f}% outside headline band [92.15, 92.35]\"""",
     "S10-bis: band now actually excludes the dropped-digit value")

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
