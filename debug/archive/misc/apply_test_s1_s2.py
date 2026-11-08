"""DELTA code-review S1 and S2.

S1 -- `test_neumann_truncation_is_exact_not_merely_converged` varies nothing.
      The internal cap collapses l_neumann = 12, 16, 20 to a single
      computation (the cap at (2,2,|m|<=1) is 10 and the true cutoff is 8), so
      the measured spread is EXACTLY 0.0 by construction.  The only regime
      where truncation is real -- l_neumann = 6, 7 -- is untested, and l = 6
      would violate the assertion (it shifts the energy 2.07e-8 Ha).
      Fix: assert invariance ABOVE the cutoff and sensitivity BELOW it.  A
      guard that cannot distinguish "converged" from "nothing varied" is the
      failure this corpus keeps finding.

S2 -- the independent Gaussian cross-check is load-bearing in Paper 12
      ("Second, an independent route agrees") and lives only in
      debug/p12_m_channel_validate.py, which Sec. 9 prunes by design.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

T = "tests/test_paper12_azimuthal_channels.py"

OLD = '''    energies = [_energy(2, 2, 1, l_neumann=lm)[0] for lm in (8, 12, 16, 20)]
    spread = max(energies) - min(energies)
    assert spread < 1e-10, (
        f"energy moves by {spread:.2e} Ha across l_neumann = 8..20; the "
        f"selection rule is not being imposed exactly"
    )
    assert all(e > E_EXACT for e in energies), (
        f"a truncation produced a non-variational energy: {energies}"
    )'''

NEW = '''    above = [_energy(2, 2, 1, l_neumann=lm)[0] for lm in (8, 12, 16, 20)]
    spread = max(above) - min(above)
    assert spread < 1e-10, (
        f"energy moves by {spread:.2e} Ha across l_neumann = 8..20; the "
        f"selection rule is not being imposed exactly"
    )
    assert all(e > E_EXACT for e in above), (
        f"a truncation produced a non-variational energy: {above}"
    )

    # The invariance above is necessary but NOT sufficient: the internal cap
    # (l <= Q + 2s - m) collapses 12, 16 and 20 to one computation, so a
    # spread of exactly 0.0 is also what a test that varies nothing returns.
    # Below the true cutoff -- 8 for this basis -- truncation is real, and the
    # guard must be able to SEE it, or it is measuring its own cap.
    below = _energy(2, 2, 1, l_neumann=6)[0]
    truncation_effect = abs(below - above[0])
    assert truncation_effect > 1e-9, (
        f"dropping to l_neumann = 6 moved the energy by only "
        f"{truncation_effect:.2e} Ha; this test cannot distinguish an exact "
        f"selection rule from an internal cap that makes every tested point "
        f"the same computation"
    )'''

with io.open(T, encoding="utf-8") as fh:
    t = fh.read()

if OLD not in t:
    print("FAILED: S1 body not matched")
    sys.exit(1)
t = t.replace(OLD, NEW, 1)

GAUSS = '''

@pytest.mark.slow
def test_independent_gaussian_route_agrees():
    """REJECTS: the diagnosis being an artifact of the prolate basis or of the
    Neumann kernel.

    Paper 12's "Second, an independent route agrees" is load-bearing: it is
    what turns "our recomputation disagrees with our earlier reading" into
    "the gap is a property of the CONFIGURATION SPACE".  It lived only in
    debug/, which Sec. 9 prunes by design, so it had no permanent home.

    Cartesian Gaussians through the corpus's own McMurchie-Davidson engine --
    different functions, different integrals, different code, and a
    well-conditioned basis, so none of the prolate machinery's linear
    dependence is in play.  Restricting to m = 0 orbitals reproduces the
    sigma-only ceiling; releasing |m| = 1 closes the gap.
    """
    import numpy as np
    from geovac.noci_engine import (
        BasisFn, integral_set_md, lowdin_orbitals, transform_integrals,
        fci_ground,
    )

    r = 1.4011
    centers = [np.array([0.0, 0.0, -r / 2]), np.array([0.0, 0.0, r / 2])]
    s_exp = [0.0347, 0.0925, 0.2469, 0.6584, 1.7557, 4.6819, 12.485, 33.293]
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
                is_sigma.append(sig)

    s, h, g = integral_set_md(orbs, [(c, 1.0) for c in centers])
    e_nuc = 1.0 / r

    def fci(mask):
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
    assert 91.5 < pct_sigma < 93.0, (
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
    )
'''

anchor = "# ======================================================================\n# 7. The conditioning caveat the paper states"
if anchor not in t:
    print("FAILED: Gaussian anchor not found")
    sys.exit(1)
t = t.replace(anchor, GAUSS.strip("\n") + "\n\n\n" + anchor, 1)

with io.open(T, "w", encoding="utf-8") as fh:
    fh.write(t)
print("  + S1: truncation guard now asserts sensitivity below the cutoff")
print("  + S2: Gaussian cross-check promoted to tests/ (@slow)")
