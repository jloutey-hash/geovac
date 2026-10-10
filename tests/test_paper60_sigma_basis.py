"""Backing for Paper 60's full-sigma-basis claims (added 2026-10-09).

sec:molecular, paragraph "The full sigma basis": on n <= N, l < n, m = 0 per centre the
Shibuya-Wulfman metric's conditioning is exponential in the number of shells, the gerade sector
included, the Bessel deficit of a displaced 1s falls without plateau, and the near-null directions
multiply -- so the s-sector levers (KMS law, flat gerade constant, tri(1,2,1) preconditioner,
one-direction rotation) do not transfer as stated.

sec:obstruction, fixed-S paragraph: on the same sigma basis, the fixed-R Galerkin posing
(eq:sw) beats the fixed-S Shibuya-Wulfman posing by orders of magnitude in H2+ energy error.

Builder: geovac/sturmian_sigma_basis.py (position-space prolate quadrature).  Its s-only path is
anchored to the momentum-space sine-basis plateau of tests/test_paper60_one_direction.py.
"""
import numpy as np
import pytest

from geovac.sturmian_sigma_basis import (
    bessel_deficit, build, galerkin_energy, gerade_ungerade, kernel_energy)

H2P_EXACT_R2 = -1.102634214   # H2+ electronic energy at R = 2 bohr


def test_s_sector_path_reproduces_the_recorded_plateau():
    """Anchor: the s-only path of the new builder reproduces the 0.380/0.696/0.907 plateau that
    tests/test_paper60_one_direction.py measures by an independent momentum-space route.

    WRONG ANSWER REJECTED: a builder bug (normalisation, centre placement, metric definition)
    that would make every sigma-basis number below meaningless.  Three kR values, 0.006
    tolerance, against values computed elsewhere.
    """
    for kR, ref in ((1.0, 0.380), (2.0, 0.696), (4.0, 0.907)):
        assert abs(bessel_deficit(12, kR, sonly=True) - ref) < 0.006, kR


def _avery_sw_closed_forms():
    """Two-centre Coulomb-Sturmian SW integrals <chi_mu'(x-X1)|(-lap+k^2)/(2k^2)|chi_mu(x-X2)>,
    S = k(X2 - X1) along z, after Avery & Avery (2017), table of Shibuya-Wulfman integrals.
    Two entries are written in this test's convention (checked 2026-10-09):
      * the 3s-2p_z entry is <3s(X1)|2p_z(X2)>; <2p_z(X1)|3s(X2)> is its exact negative, as
        inversion symmetry requires;
      * the 2p_z-2p_z entry is normalised to +1 at S = 0 (same-centre orthonormality).
    Keys are ((n', l') on X1, (n, l) on X2), m = 0."""
    e = np.exp
    return {
        ((1, 0), (1, 0)): lambda S: e(-S) * (1 + S),
        ((1, 0), (2, 0)): lambda S: -2 / 3 * e(-S) * S ** 2,
        ((1, 0), (2, 1)): lambda S: -2 / 3 * e(-S) * (1 + S) * S,
        ((1, 0), (3, 0)): lambda S: 1 / 3 * e(-S) * (-2 + S) * S ** 2,
        ((2, 0), (2, 0)): lambda S: 1 / 3 * e(-S) * (3 + 3 * S - 2 * S ** 2 + S ** 3),
        ((2, 0), (2, 1)): lambda S: 1 / 3 * e(-S) * (-1 - S + S ** 2) * S,
        ((2, 0), (3, 0)): lambda S: 2 / 15 * e(-S) * S ** 2 * (-10 + 5 * S - S ** 2),
        ((2, 1), (2, 1)): lambda S: -1 / 3 * e(-S) * (-3 - 3 * S + S ** 3),   # +1 at S = 0
        ((3, 0), (2, 1)): lambda S: 2 / 15 * e(-S) * (4 - S) * S ** 3,        # <3s(X1)|2p_z(X2)>
        ((3, 0), (3, 0)): lambda S: 1 / 45 * e(-S) * (45 + 45 * S - 60 * S ** 2 + 60 * S ** 3
                                                      - 18 * S ** 4 + 2 * S ** 5),
    }


def test_sw_matrix_matches_avery_closed_forms():
    """Independent-route oracle: the position-space quadrature reproduces Avery & Avery's
    closed-form two-centre SW integrals (momentum-space / 4D-harmonic derivation) for all ten
    tabulated entries at three separations.

    WRONG ANSWER REJECTED: any error in the radial or angular normalisation, the sign convention
    of l > 0 functions, the centre placement, or the metric definition -- each moves at least
    one entry by far more than the 1e-9 tolerance (measured agreement ~1e-12).
    """
    forms = _avery_sw_closed_forms()
    for s in (0.5, 1.5, 3.0):
        S, _, labels = build(3, 1.0, s)
        m = len(labels)
        idx = {lab: i for i, lab in enumerate(labels)}
        for (bra, ket), f in forms.items():
            assert abs(S[idx[bra], m + idx[ket]] - f(s)) < 1e-9, (s, bra, ket)


def test_intra_centre_block_is_identity_on_the_sigma_basis():
    """The SW metric's intra-centre block is the identity for l > 0 too (potential-weighted
    orthonormality), so the A-set is S-orthonormal and the Bessel deficit is meaningful.

    WRONG ANSWER REJECTED: a mis-normalised l > 0 radial function or angular factor, which shows
    up as an intra-block deviation of order 1e-2 or worse.
    """
    S, _, labels = build(4, 1.3, 2.0)
    m = len(labels)
    assert np.max(np.abs(S[:m, :m] - np.eye(m))) < 1e-9
    assert np.max(np.abs(S[m:, m:] - np.eye(m))) < 1e-9


def test_sigma_deficit_falls_without_plateau():
    """sec:molecular: with l > 0 the deficit keeps falling, so the s-sector plateau is the l = 0
    projection and not incompleteness of the Coulomb-Sturmian set.

    WRONG ANSWER REJECTED: a plateau (the withdrawn general reading).  Asserted as strict
    monotone decrease over N = 2..6, a last value below a quarter of the s-sector plateau, and a
    last-step decrease of at least 15% (a flat tail would fail it).
    """
    vals = [bessel_deficit(N, 1.0) for N in (2, 3, 4, 5, 6)]
    assert all(b < a for a, b in zip(vals, vals[1:])), vals
    assert vals[-1] < 0.25 * 0.380, vals
    assert vals[-1] < 0.85 * vals[-2], vals


def test_sigma_conditioning_is_exponential_not_polynomial():
    """sec:molecular: cond(S) grows exponentially in the number of shells at kR = 3.

    WRONG ANSWER REJECTED: polynomial growth (the s-sector law).  Going from 3 to 10 functions
    per centre, even an N^4 law would grow cond by about 120x; the measured growth is ~7e5x.
    Asserted as > 1e4 overall AND every shell step multiplying cond by > 100.
    """
    conds = []
    for N in (2, 3, 4):
        S, _, _ = build(N, 1.0, 3.0)
        ev = np.linalg.eigvalsh(S)
        conds.append(ev[-1] / ev[0])
    assert conds[2] / conds[0] > 1e4, conds
    assert all(b / a > 100 for a, b in zip(conds, conds[1:])), conds


def test_sigma_gerade_sector_is_not_flat():
    """sec:molecular: the gerade sector is no refuge on the sigma basis.

    WRONG ANSWER REJECTED: a flat gerade sector (the s-sector constant 2.555).  At kR = 3 the
    s-sector gerade cond is 2.0-2.4; the sigma gerade cond at N = 4 must exceed 1e3 and grow by
    > 5x per shell (a flat sector steps by ~1.1x; measured 9x then 350x), while the gerade/ungerade split itself stays exact (leak < 1e-10, so the
    number is not an artefact of a wrong parity assignment).
    """
    conds = []
    for N in (2, 3, 4):
        S, _, labels = build(N, 1.0, 3.0)
        Sg, _ = gerade_ungerade(S, labels)
        m = len(labels)
        pi = np.diag([(-1.0) ** l for (_, l) in labels])
        U = np.vstack([np.hstack([np.eye(m), np.eye(m)]), np.hstack([pi, -pi])]) / np.sqrt(2)
        assert np.max(np.abs((U.T @ S @ U)[:m, m:])) < 1e-10
        conds.append(np.linalg.cond(Sg))
    s_only_ref, _, s_labels = build(8, 1.0, 3.0, sonly=True)
    assert np.linalg.cond(gerade_ungerade(s_only_ref, s_labels)[0]) < 3.0
    assert conds[2] > 1e3, conds
    assert all(b / a > 5 for a, b in zip(conds, conds[1:])), conds


def test_fixed_R_galerkin_beats_fixed_S_kernel_on_the_same_basis():
    """sec:obstruction: same sigma basis (6 functions per centre), H2+ at R = 2 bohr.

    WRONG ANSWER REJECTED: (i) the two posings being comparable -- asserted as Galerkin error
    < 2e-4 Ha and kernel error > 1e-2 Ha, a > 50x gap; (ii) either posing landing below the
    exact energy (both converge from above here); (iii) a builder that is simply inaccurate --
    the Galerkin bound is tight enough that a broken metric or potential cannot meet it.
    """
    eg = galerkin_energy(3, 2.0) - H2P_EXACT_R2
    ek = kernel_energy(3, 2.0) - H2P_EXACT_R2
    assert 0.0 < eg < 2e-4, eg
    assert ek > 1e-2, ek
    assert ek / eg > 50, (eg, ek)
