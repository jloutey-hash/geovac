"""Backing tests for Paper 58's NOCI legs, via geovac.noci_engine.

The engine was promoted from an exploratory driver to ``geovac/noci_engine.py``
for this paper, so these tests import a PERMANENT module rather than reaching
into the transient debug/ tree (cite-permanent-records policy).

Covers:
  test_paper58_sto_fit_quality      paper's <fit|STO> = 1.000000 claim
  test_paper58_hydrogenic_limits    analytical-limit check on the fitted shapes
  test_paper58_span_identity        orthogonalization is HARMLESS on the
                                    complete space (the precondition without
                                    which the truncation observation is empty)
  test_paper58_truncation_damage    Obs. "orthogonalization damages only
                                    truncated spaces", fixed-geometry form

Scope, stated: the paper's Table 2 reports NaH BINDING percentages, which need
the dissociated limit and the full zeta-selection pipeline.  What is asserted
here is the fixed-geometry form of the mechanism on LiH -- energy error against
complete-space FCI -- which is the load-bearing content of the observation.
The NaH Table 2 numbers themselves remain OWED.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy", reason="engine needs scipy.special.hyp1f1")

from geovac import noci_engine as E  # noqa: E402

UP, DN = 0, 1


# ---------------------------------------------------------------------------
# Shape fitting
# ---------------------------------------------------------------------------

def _shapes():
    out = {}
    for kind, (l, n_r) in (("1s", (0, 1)), ("2s", (0, 2)),
                           ("2p", (1, 2)), ("3s", (0, 3))):
        a, d, q = E.fit_sto_shape(l, n_r)
        out[kind] = (a, d, q)
    return out


def test_paper58_sto_fit_quality():
    """Paper: 6-Gaussian STO shape fits reach <fit|STO> = 1.000000.

    Guards the claim as stated (six digits), not a loose tolerance -- the
    paper uses this to argue NaH's D_e deficit is basis incompleteness rather
    than a fitting artifact, so a sloppy fit would undercut that reasoning.
    """
    sh = _shapes()
    for kind in ("2p", "3s"):
        q = sh[kind][2]
        assert q >= 0.9999995, (
            f"{kind} fit quality {q:.7f} < 1.000000 to six digits; the paper's "
            "'fits were essentially perfect, so the D_e deficit is basis "
            "incompleteness' argument weakens"
        )
    assert sh["1s"][2] >= 0.99999, f"1s fit quality {sh['1s'][2]:.7f}"


def test_paper58_hydrogenic_limits():
    """Analytical-limit check: fitted shapes reproduce exact hydrogen energies.

    E(1s) = -1/2 and E(2p, zeta=1/2) = -1/8 exactly.  This validates the
    fitted basis against known values rather than against itself.
    """
    sh = {k: (v[0], v[1]) for k, v in _shapes().items()}
    orig = np.array([0.0, 0.0, 0.0])
    nuc = [(orig, 1.0)]

    chi = E.sto_shape_basis(orig, "1s", 1.0, sh, (0, 0, 0))
    s, h, _g = E.integral_set_md([chi], nuc)
    e_1s = h[0, 0] / s[0, 0]
    assert abs(e_1s + 0.5) < 1e-3, f"E(H 1s) = {e_1s:.6f}, exact -0.5"

    chi = E.sto_shape_basis(orig, "2p", 0.5, sh, (0, 0, 1))
    s, h, _g = E.integral_set_md([chi], nuc)
    e_2p = h[0, 0] / s[0, 0]
    assert abs(e_2p + 0.125) < 1e-3, f"E(H 2p) = {e_2p:.6f}, exact -0.125"


# ---------------------------------------------------------------------------
# LiH ladder: span identity and truncation damage
# ---------------------------------------------------------------------------

def _lih_system(R: float = 3.25):
    """Li(1s, 2s) + H(1s), engine-native fitted shapes.  4 electrons."""
    sh = {k: (v[0], v[1]) for k, v in _shapes().items()}
    pos_li = np.array([0.0, 0.0, 0.0])
    pos_h = np.array([0.0, 0.0, R])
    orbs = [
        E.sto_shape_basis(pos_li, "1s", 2.69, sh, (0, 0, 0)),
        E.sto_shape_basis(pos_li, "2s", 0.65, sh, (0, 0, 0)),
        E.sto_shape_basis(pos_h, "1s", 1.00, sh, (0, 0, 0)),
    ]
    nuclei = [(pos_li, 3.0), (pos_h, 1.0)]
    s, h, g = E.integral_set_md(orbs, nuclei)
    return s, h, g


def _ladder_dets():
    core = [(0, UP), (0, DN)]
    cov_a = core + [(1, UP), (2, DN)]
    cov_b = core + [(2, UP), (1, DN)]
    ion_h = core + [(2, UP), (2, DN)]
    return cov_a, cov_b, ion_h


def test_paper58_span_identity():
    """Orthogonalization is HARMLESS on the complete space.

    Without this the truncation observation would be vacuous -- the damage has
    to be attributable to TRUNCATION, not to orthogonalization corrupting the
    integrals.  Non-orthogonal complete-space FCI must equal Loewdin-basis
    bitstring FCI to machine precision.
    """
    s, h, g = _lih_system()
    x = E.lowdin_orbitals(s)
    h_o, g_o = E.transform_integrals(x, h, g)

    e_orth = E.fci_ground(h_o, g_o, 4)

    # Same space, expressed non-orthogonally: all determinants of 3 spatial
    # orbitals with 4 electrons, evaluated by non-orthogonal Slater-Condon.
    from itertools import combinations
    spin_orbs = [(p, sig) for p in range(h.shape[0]) for sig in (UP, DN)]
    dets = [list(c) for c in combinations(spin_orbs, 4)]
    e_nonorth, _ = E.noci_ground_gensc(dets, s, h, g)

    assert abs(e_orth - e_nonorth) < 1e-9, (
        f"complete-space identity broken: Loewdin FCI {e_orth:.12f} vs "
        f"non-orthogonal {e_nonorth:.12f}; the truncation observation cannot "
        "be attributed to truncation if the complete spaces disagree"
    )


def test_paper58_truncation_damage():
    """Obs.: orthogonalization damages TRUNCATED spaces (fixed-geometry form).

    Three fragment-native determinants land close to complete-space FCI; the
    SAME COUNT of determinants after Loewdin orthogonalization is far worse.
    Asserted as a ratio of energy errors so it guards the mechanism rather
    than a fitted magnitude.
    """
    s, h, g = _lih_system()
    cov_a, cov_b, ion_h = _ladder_dets()
    n_spatial = h.shape[0]

    _h_o, _g_o = E.transform_integrals(E.lowdin_orbitals(s), h, g)
    e_fci = E.fci_ground(_h_o, _g_o, 4)

    e_noci3, _ = E.noci_ground_gensc([cov_a, cov_b, ion_h], s, h, g)

    # Same three determinants, built on Loewdin-orthogonalized orbitals.
    x = E.lowdin_orbitals(s)
    h_o, g_o = E.transform_integrals(x, h, g)
    s_id = np.eye(n_spatial)
    e_low3, _ = E.noci_ground_gensc([cov_a, cov_b, ion_h], s_id, h_o, g_o)

    err_noci = e_noci3 - e_fci
    err_low = e_low3 - e_fci

    assert err_noci > -1e-9, (
        f"NOCI-3 ({e_noci3:.9f}) below complete-space FCI ({e_fci:.9f}); "
        "a truncated variational energy cannot be lower"
    )
    assert err_low > -1e-9, "Loewdin-3 below complete-space FCI"
    assert err_low > err_noci, (
        f"Loewdin truncation was not worse: err_noci={err_noci:.6f} Ha vs "
        f"err_low={err_low:.6f} Ha.  Paper 58 Obs. 'orthogonalization damages "
        "only truncated spaces' would need revisiting."
    )
    # The paper's claim is a large gap, not a marginal one (LiH: 97.6% vs 15%).
    assert err_low > 5.0 * err_noci, (
        f"gap too small to support the paper's framing: err_noci="
        f"{err_noci:.6f} Ha, err_low={err_low:.6f} Ha "
        f"(ratio {err_low / err_noci:.2f}x, expected >5x)"
    )
