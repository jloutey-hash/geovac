"""Backing tests for Paper 12, Sec. "Restoring the Azimuthal Channels" --
the re-conditioning result (:mod:`geovac.prolate_recondition`).

The claim under test: Paper 12's 99.1% headline is a CONDITIONING artifact of the
monomial radial basis xi^j, not a structural ceiling.  Re-basing the *same span*
onto an orthogonal-polynomial family (Laguerre x Legendre, or the mu-adapted
associated-Laguerre x Gegenbauer) is an exact change of basis that leaves the
exact-arithmetic energy unchanged but conditions the generalized eigenproblem, so
the energy climbs monotonically and variationally past the monomial wall to
chemical accuracy.

Each test names the wrong answer it rejects, per the Sec. 9 guard rule.

WHICH TESTS ARE @slow, AND WHY THAT CHANGED.  The @slow set is the one that
computes a re-based energy at a LARGE truncation.  It used to be every
energy-computing test, on the grounds that "the mpf V_ee build is inherently slow
(~2-10 min)" -- no longer true twice over: v5.13.4 replaced the V_ee B-seed
quadrature with a closed form (vee_mp 139 s -> 3.7 s at (3,3,1)), and v5.13.8
made `engine="direct"` the default, which drops the mpf one-body build and one of
the two changes of basis (whole pipeline 2.7x at (3,3,1), 3.2x at (4,4,2)).  So
two energy-computing tests now run unmarked at (2,2,1), and the phase that
actually costs is V_ee -- its mpf build plus its one surviving re-basing -- not
the one-body half.

Fire-tested by ``debug/firetest_p12_recondition.py``.
"""

from __future__ import annotations

import numpy as np
import mpmath as mp
import pytest

from geovac import prolate_general_m as pg
from geovac import prolate_recondition as pr

E_EXACT = pr.E_EXACT          # -1.174475 (Kolos & Wolniewicz)
DE_EXACT = pr.DE_EXACT        # 0.174475


# ======================================================================
# FAST: algebraic invariants (no mpf V_ee build)
# ======================================================================

def test_orthogonal_families_reduce_to_the_mu0_reference():
    """REJECTS: a mu-adapted family whose mu = 0 sector is NOT the plain
    Laguerre x Legendre reference.

    The Gegenbauer weight parameter lam = mu + 1/2 gives lam = 1/2 = Legendre at
    mu = 0, and the generalized Laguerre L_n^{(mu)} gives beta = 0 = plain
    Laguerre.  If either reduction broke, the "gegenbauer" basis would silently
    span a different mu = 0 sector from the reference and the same-span guarantee
    below would be vacuous.
    """
    with mp.workdps(40):
        for n in range(6):
            g = pr.gegenbauer_coeffs(n, mp.mpf('0.5'), 8)
            le = pr.legendre_coeffs(n, 8)
            assert max(abs(g[i] - le[i]) for i in range(8)) < mp.mpf(10) ** -30
            a = pr.assoc_laguerre_coeffs(n, 0, 1.0, 8)
            b = pr.laguerre_coeffs(n, 1.0, 8)
            assert max(abs(a[i] - b[i]) for i in range(8)) < mp.mpf(10) ** -30


def test_factored_change_of_basis_equals_the_dense_one():
    """REJECTS: a wrong Kronecker contraction / transpose in the factored change
    of basis.

    The whole re-conditioning rests on C M C^T being computed EXACTLY as the
    dense congruence, only faster (the memo's 4e-36 factored-vs-dense check).  C
    is block-diagonal in mu with each block = kron(Tr_mu, Ta_mu); this asserts
    the fast :func:`_factored_cob` reproduces that dense product on a random
    symmetric matrix.  A single mis-placed axis in the tensordots would fail here
    but leave the energies looking plausible.
    """
    with mp.workdps(30):
        j_max, l_max, mu_max = 2, 2, 1
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1)
                  if (l + m) % 2 == 0])
        Nmu = mu_max + 1
        N = Nmu * Nr * Na
        Tr_list, Ta_list = pr._transforms_per_mu("gegenbauer", j_max, l_max,
                                                  mu_max, 1.0)
        rng = np.random.default_rng(0)
        M = np.empty((N, N), object)
        for i in range(N):
            for j in range(i, N):
                v = mp.mpf(float(rng.standard_normal()))
                M[i, j] = M[j, i] = v

        # dense C = blockdiag_mu( kron(Tr_mu, Ta_mu) )
        C = np.zeros((N, N), object)
        C[:] = mp.mpf(0)
        for mu in range(Nmu):
            sl = slice(mu * Nr * Na, (mu + 1) * Nr * Na)
            C[sl, sl] = np.kron(Tr_list[mu], Ta_list[mu])
        dense = C @ M @ C.T
        fast = pr._factored_cob(M, Nmu, Nr, Na, Tr_list, Ta_list)

        err = max(abs(dense[i, j] - fast[i, j]) for i in range(N) for j in range(N))
        assert err < mp.mpf(10) ** -25, (
            f"factored change of basis differs from the dense one by {float(err):.2e}; "
            f"the Kronecker contraction is wrong"
        )


def _naive_canonical_solve(S_o, H_o, thresh=1e-11):
    """A plain float64 canonical orthogonalization on the UN-normalized S_o --
    the wrong answer the normalized solve is meant to beat.  Discards directions
    whose raw S eigenvalue is below `thresh` * max, which large norm spread makes
    misjudge."""
    n = S_o.shape[0]
    Sf = np.array([[float(S_o[i, j]) for j in range(n)] for i in range(n)])
    Hf = np.array([[float(H_o[i, j]) for j in range(n)] for i in range(n)])
    w, U = np.linalg.eigh(0.5 * (Sf + Sf.T))
    keep = w > thresh * w[-1]
    X = U[:, keep] / np.sqrt(w[keep])
    return float(np.linalg.eigvalsh(X.T @ (0.5 * (Hf + Hf.T)) @ X)[0])


def test_normalized_solve_recovers_what_norm_spread_hides():
    """REJECTS: dropping the unit-norm rescaling from the solve.

    The re-based functions have wildly different norms, so cond(S_o) is inflated
    by pure norm spread.  A plain canonical orthogonalization then DISCARDS
    directions that carry real variational content -- the reported energy comes
    out too HIGH.  Rescaling each function to unit norm is a diagonal congruence
    that leaves the true generalized eigenvalues unchanged but exposes the real
    near-dependence, recovering the correct lowest energy.

    Constructed with a well-conditioned (S_hat, H_hat) whose true lowest
    eigenvalue is known, then scaled by a huge diagonal congruence D (spread
    1e8): the eigenvalues are unchanged, but the plain solve loses content and
    the normalized solve does not.  A guard asserting mere invariance would be
    vacuous -- generalized eigenvalues are congruence-invariant for ANY sound
    solve -- so this asserts the sharper claim that the normalized solve BEATS
    the plain one on the norm-spread pathology.
    """
    rng = np.random.default_rng(1)
    n = 12
    B = rng.standard_normal((n, n))
    Shat = B @ B.T + n * np.eye(n)                    # well-conditioned SPD
    Dh = np.sqrt(np.diag(Shat))
    Shat = (Shat / Dh[:, None]) / Dh[None, :]         # unit-diagonal correlation
    A = rng.standard_normal((n, n))
    Hhat = 0.5 * (A + A.T)
    from scipy.linalg import eigh as _eigh
    E_true = float(_eigh(Hhat, Shat, eigvals_only=True)[0])

    d = np.geomspace(1.0, 1e8, n)                     # huge norm spread
    S_o = np.empty((n, n), object)
    H_o = np.empty((n, n), object)
    for i in range(n):
        for j in range(n):
            S_o[i, j] = mp.mpf(float(Shat[i, j] * d[i] * d[j]))
            H_o[i, j] = mp.mpf(float(Hhat[i, j] * d[i] * d[j]))

    E_norm, _c, _k, _s = pr._normalized_solve(S_o, H_o)
    E_naive = _naive_canonical_solve(S_o, H_o)

    assert abs(E_norm - E_true) < 1e-8, (
        f"normalized solve gives {E_norm:.8f}, true lowest is {E_true:.8f}; the "
        f"norm-spread rescaling failed to recover the eigenvalue"
    )
    assert E_naive > E_true + 1e-3, (
        f"the plain (un-normalized) solve gave {E_naive:.6f}, not measurably "
        f"above the true {E_true:.6f}; then the norm-spread pathology is absent "
        f"from this fixture and the test does not discriminate the normalization"
    )


# ======================================================================
# SLOW: the physics claims (mpf V_ee build)
# ======================================================================

@pytest.mark.slow
def test_rebasing_preserves_the_energy_where_the_monomial_is_sound():
    """REJECTS: a change of basis that is not span-preserving.

    At (2,2), mu <= 1 the monomial basis is still sound (cond(S) = 6.8e10,
    variational, 98.95% of D_e).  Re-basing the SAME span must land on the SAME
    energy to many digits -- that is what makes it a legitimate re-conditioning
    and not a different calculation.  A per-function-exponent contamination or a
    dropped term would change the energy here.
    """
    alpha = 1.0
    basis = pg.generate_basis(2, 2, 1, alpha)
    mom = pg.Moments(2.0 * alpha, 6 * 2 + 6 * 3 + 20)
    grid = pg.XiGrid(alpha)
    s, h1 = pg.one_body(basis, pr.R_DEFAULT, 1.0, mom)
    v = pg.vee_matrix(basis, pr.R_DEFAULT, grid, l_neumann=14)
    h = h1 + v + (1.0 / pr.R_DEFAULT) * s
    e_mono, _nk, _tot = pg.solve_generalized(h, s)
    assert e_mono > E_EXACT, "monomial (2,2,1) should be variational"

    res = pr.recondition_energy(2, 2, 1, alpha=alpha)
    assert res.variational
    assert abs(res.energy - e_mono) < 1e-5, (
        f"re-based energy {res.energy:.7f} differs from the sound monomial "
        f"{e_mono:.7f} by {1e3 * abs(res.energy - e_mono):.3f} mHa; the change "
        f"of basis is not span-preserving"
    )


@pytest.mark.slow
def test_rebasing_breaks_the_conditioning_wall():
    """REJECTS: the reading that 99.1% is a STRUCTURAL ceiling.

    Same span, same alpha = 1.0, at (3,3), mu <= 1 (N = 144):
      - the monomial solve is NON-variational -- E = -64.3 Ha, cond(S) = 2.9e16;
      - the re-based solve is variational and reaches ~99.2% of D_e.
    The contrast is the whole claim: the monomial wall is conditioning, and
    re-basing removes it.  If the re-based solve were also non-variational, or
    stalled at the monomial 99.09%, this claim would be unsupported.
    """
    alpha = 1.0
    basis = pg.generate_basis(3, 3, 1, alpha)
    mom = pg.Moments(2.0 * alpha, 6 * 3 + 6 * 3 + 20)
    grid = pg.XiGrid(alpha)
    s, h1 = pg.one_body(basis, pr.R_DEFAULT, 1.0, mom)
    v = pg.vee_matrix(basis, pr.R_DEFAULT, grid, l_neumann=14)
    h = h1 + v + (1.0 / pr.R_DEFAULT) * s
    e_mono, _nk, _tot = pg.solve_generalized(h, s)
    assert e_mono < E_EXACT, (
        f"monomial (3,3,1) at alpha=1.0 returned E = {e_mono:.4f}, expected the "
        f"documented non-variational break (cond(S) = 2.9e16); the wall this "
        f"test contrasts against is gone"
    )

    res = pr.recondition_energy(3, 3, 1, alpha=alpha)
    assert res.variational, (
        f"re-based (3,3,1) is non-variational (E = {res.energy:.4f}); the "
        f"re-conditioning failed"
    )
    assert 99.0 < res.de_pct < 99.4, (
        f"re-based (3,3,1) reaches {res.de_pct:.3f}% of D_e; expected ~99.22%, "
        f"decisively past the monomial 99.09% cap"
    )
    assert res.cond_norm < 1e9, (
        f"normalized cond = {res.cond_norm:.1e}; the re-based basis is supposed "
        f"to be well-conditioned (~1e5-1e6), not near the monomial 1e16"
    )


@pytest.mark.slow
def test_gegenbauer_is_the_same_span_and_better_conditioned():
    """REJECTS: (a) a mu-adapted basis that silently changes the span/energy
    (the per-function-exponent trap that kills completeness -- failed-approaches
    ledger 2026-08-26), and (b) a "mu-adapted" family that does NOT condition
    better than plain Laguerre x Legendre (then it is not actually adapted).

    Both families re-base the identical monomial span, so at equal (j,l,mu) they
    MUST give the same energy; the Gegenbauer family, matched to the mu weight,
    must give a smaller normalized condition number.
    """
    ll = pr.recondition_energy(3, 3, 1, alpha=1.0, basis="laguerre_legendre")
    gg = pr.recondition_energy(3, 3, 1, alpha=1.0, basis="gegenbauer")
    assert ll.variational and gg.variational
    assert abs(ll.energy - gg.energy) < 1e-6, (
        f"the two families disagree by {1e6 * abs(ll.energy - gg.energy):.2f} uHa "
        f"at equal (3,3,1); they should span the identical space -- a difference "
        f"means one family changed the exponent, not just the polynomial basis"
    )
    assert gg.cond_norm < ll.cond_norm, (
        f"gegenbauer cond {gg.cond_norm:.2e} is not below laguerre_legendre "
        f"{ll.cond_norm:.2e}; the mu-adaptation bought no conditioning"
    )


def test_direct_engine_matches_the_mpf_reference_one_body():
    """REJECTS: a direct orthogonal-basis one-body build that is not the same
    operator as the mpf monomial build re-based.

    `engine="direct"` never forms the monomial one-body matrices, so it bypasses
    the change of basis entirely rather than doing it more cheaply.  That is only
    legitimate if it lands on the SAME matrices the mpf route produces, which is
    the definition of correctness here.  Checked at (2,2) mu<=1 in BOTH families,
    since the mu-adapted family is a separate code path (different radial AND
    angular polynomials) and `recondition_energy`'s default family is the other
    one.

    The bar is SCALE-relative (max|dX| / max|X|), not entry-relative.  Entry-
    relative is not a soundness criterion for H1: its entries span many orders
    and some pass near zero by T/V_ne cancellation, so a negligible absolute
    deviation over a near-zero denominator yields a large relative number with
    nothing wrong (measured: 9.2e-13 entry-relative at (5,5)+delta comes from a
    single element 10 orders below the matrix scale, whose absolute deviation is
    ~1e-22 of that scale, while scale-relative sits at 3.7e-16).
    """
    alpha = 1.0
    for basis in ("laguerre_legendre", "gegenbauer"):
        Sd, Hd = pr.build_one_body_direct(2, 2, 1, alpha, basis=basis)
        with mp.workdps(pr.DEFAULT_DPS):
            idx = pr._product_index(2, 2, 1)
            fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
            A = pr.ngm._mono_moments(2.0 * alpha, 6 * 2 + 6 * 3 + 20)
            S, H1 = pr.one_body_mp(fns, alpha, pr.R_DEFAULT, A)
            Nr = (2 + 1) ** 2
            Na = len([1 for l in range(3) for m in range(3) if (l + m) % 2 == 0])
            Tr, Ta = pr._transforms_per_mu(basis, 2, 2, 1, alpha)
            Sg = pr._to_f64(pr._factored_cob(S, 2, Nr, Na, Tr, Ta))
            Hg = pr._to_f64(pr._factored_cob(H1, 2, Nr, Na, Tr, Ta))
        for name, D, G in (("S", Sd, Sg), ("H1", Hd, Hg)):
            scale = float(np.abs(D - G).max() / np.abs(G).max())
            assert scale < 1e-13, (
                f"{basis} {name}: direct build differs from the mpf re-based "
                f"reference by {scale:.2e} scale-relative; the direct engine is "
                f"not computing the same operator"
            )


def test_direct_and_mpf_engines_agree_end_to_end():
    """REJECTS: an assembly error that cancels in the matrices but moves the
    energy -- a wrong nuclear-repulsion shift, or a V_ee re-basing applied to the
    wrong engine's matrices.

    The direct path assembles H_o = H1_o + cob(V) + S_o/R in float64, relying on
    linearity of the change of basis; the mpf path re-bases H = H1 + V + S/R as
    one object. Those are different orderings of the same sum, so they must give
    the same eigenvalue. Asserted on the ENERGY (and on the surviving dimension,
    so a silently-truncated solve cannot pass) at a truncation small enough to
    run unmarked.
    """
    a = pr.recondition_energy(2, 2, 1, alpha=1.0, engine="mpf")
    b = pr.recondition_energy(2, 2, 1, alpha=1.0, engine="direct")
    assert a.engine == "mpf" and b.engine == "direct", (
        "the result does not report the engine it used, so this test cannot "
        "tell the two paths apart"
    )
    assert abs(a.energy - b.energy) < 1e-9, (
        f"engines disagree by {1e6 * abs(a.energy - b.energy):.4f} uHa "
        f"(mpf {a.energy:.10f} vs direct {b.energy:.10f}); the float64 "
        f"H_o = H1_o + cob(V) + S_o/R assembly is not the mpf H re-based"
    )
    assert a.n_kept == b.n_kept, (
        f"engines keep different dimensions ({a.n_kept} vs {b.n_kept}); the "
        f"energies could agree while one path silently discards content"
    )


def test_alpha_default_is_the_measured_optimum():
    """REJECTS: the module default drifting away from the exponent the paper's
    headline is quoted at.

    Paper 12's headline (99.81% of D_e, 0.32 mHa) is the (5,5)+delta result at the
    variational optimum alpha = 1.40, registered as `p12_rebased_alpha_opt`.  If
    this default moved, the module would compute one number while the paper
    printed another -- which is exactly the state this test was written to end:
    the re-conditioned table was published at alpha = 1.0 while the paper's own
    text said alpha is "optimized variationally".

    Fast by construction (reads a signature; runs no physics).
    """
    import inspect
    default = inspect.signature(pr.recondition_energy).parameters['alpha'].default
    assert abs(default - 1.40) < 1e-12, (
        f"recondition_energy's alpha default is {default}, not the measured "
        f"optimum 1.40 that Paper 12's headline is quoted at"
    )


@pytest.mark.slow
def test_alpha_optimum_beats_the_fixed_alpha_ladder():
    """REJECTS: the claim that alpha = 1.0 is this basis's best variational point.

    Paper 12 prints a monotone ladder at a FIXED alpha = 1.0 -- correctly, since
    "each point contains its predecessor and lies below it" only holds at a
    consistent exponent -- and its headline was long quoted from that ladder's
    endpoint.  But alpha is a variational parameter, so the endpoint is not the
    basis's best point.  Measured at (5,5)+delta: 0.406 / 0.346 / 0.324 / 0.334
    mHa at alpha = 1.00 / 1.20 / 1.40 / 1.50, a bracketed minimum at 1.40 worth
    0.082 mHa -- a fifth of the residual.

    Checked here at (4,4)+delta instead, which is ~365 s rather than ~25 min.
    That choice is deliberate: the (5,5)+delta headline previously had no
    regression guard because it was "confirmed once, too slow for CI", and that
    gap is what let the v5.13.6 seed regression through.  A cheaper basis that
    exercises the same mechanism is worth more than an un-run exact one.

    Asserts the improvement is REAL and not a trade: lower energy, still
    variational, no functions lost, and conditioning no worse -- because a "gain"
    bought by degrading the metric would not support the headline.
    """
    lo = pr.recondition_energy(4, 4, 2, alpha=1.00, basis="laguerre_legendre")
    hi = pr.recondition_energy(4, 4, 2, alpha=1.40, basis="laguerre_legendre")

    assert lo.variational and hi.variational, "both points must be variational"
    assert hi.energy < lo.energy, (
        f"alpha=1.40 gave {hi.energy:.7f} vs alpha=1.0's {lo.energy:.7f}; the "
        f"optimum does not beat the ladder endpoint, so the headline's alpha is "
        f"not justified"
    )
    gain_mha = abs(lo.err_mha) - abs(hi.err_mha)
    assert gain_mha > 0.03, (
        f"alpha optimization bought only {gain_mha:.4f} mHa (measured 0.057); "
        f"below this the improvement is not distinguishable from noise and the "
        f"headline should stay at the ladder endpoint"
    )
    assert hi.n_kept == hi.n_basis, (
        f"alpha=1.40 kept {hi.n_kept}/{hi.n_basis} functions; a gain that "
        f"discards basis functions is not the same calculation"
    )
    assert hi.cond_norm <= lo.cond_norm, (
        f"conditioning got WORSE at the optimum ({hi.cond_norm:.2e} vs "
        f"{lo.cond_norm:.2e}); then the gain is bought from the metric and the "
        f"downcast solve is less trustworthy, not more"
    )


@pytest.mark.slow
def test_climb_reaches_chemical_accuracy_past_the_cap():
    """REJECTS: a climb that plateaus at the monomial cap, or a non-variational
    "improvement".

    The load-bearing headline: adding the delta channels on the re-based basis
    takes the energy to chemical accuracy (error < 1.6 mHa), variationally and
    monotonically below the (3,3) pi point -- decisively past the 99.1% cap that
    the monomial basis cannot exceed.  (4,4)+delta reaches 99.71% (0.505 mHa);
    the full (5,5)+delta reaches 99.767% and is confirmed once but too slow for
    CI.
    """
    e_pi = pr.recondition_energy(3, 3, 1, alpha=1.0)
    e_delta = pr.recondition_energy(4, 4, 2, alpha=1.0)
    assert e_pi.variational and e_delta.variational
    assert e_delta.energy < e_pi.energy, (
        f"adding delta at a larger basis did not lower the energy "
        f"({e_delta.energy:.7f} vs {e_pi.energy:.7f}); the climb is not monotone"
    )
    # abs() is load-bearing, and its absence made this assertion DEAD.
    # `err_mha` is SIGNED -- (E_exact - E)*1000 -- so every variational result is
    # NEGATIVE and `err_mha < 1.6` passed unconditionally: measured, it passes at
    # E = -1.00, i.e. no binding at all (0.00% of D_e, err = -174.475 mHa).  The
    # only case it rejected was E BELOW exact, which the `variational` assertion
    # above already covers.  So the "inside chemical accuracy" claim had zero
    # enforcement until this line took the magnitude.
    assert abs(e_delta.err_mha) < 1.6, (
        f"(4,4)+delta error is {abs(e_delta.err_mha):.3f} mHa, not inside "
        f"chemical accuracy (1.6 mHa); the re-based climb did not reach it"
    )
    assert e_delta.de_pct > 99.5, (
        f"(4,4)+delta reaches only {e_delta.de_pct:.3f}% of D_e; expected "
        f"~99.71%, decisively past the 99.1% monomial cap"
    )
