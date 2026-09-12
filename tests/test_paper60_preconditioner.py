"""Backing for Paper 60's third conditioning lever (sec:resource, added 2026-09-12):
a band-Toeplitz preconditioner in the sense of Serra, Math. Comp. 66, 651 (1997).

Written as a SEPARATE pass from the paper edit it protects (CLAUDE.md Sec. 9),
against the TRACKED engine `geovac.sturmian_sigma_law`, not against a private
reimplementation.  Each guard names the wrong answer it excludes; each is
fire-tested against that answer via debug/qa/fire_test.py.

Claims backed:
  1. The matching polynomial g = 2 + 2 cos chi has a Hankel part that vanishes
     identically in this basis, so P is EXACTLY tridiagonal(1, 2, 1).
  2. P is diagonalized exactly by the DST-I, with closed-form eigenvalues.
  3. cond(G) = cond(P^-1/2 (I-C) P^-1/2) is BOUNDED while cond(I-C) grows ~ n^2.
  4. The whitening X = P^-1/2 G^-1/2 satisfies X^T (I-C) X = I, so the
     generalized spectrum is unchanged -- the lever is not bought by changing
     the problem.
  5. The locality gain is REAL at 1e-2 and ERODES at 1e-3 (algebraic decay),
     and the clean discriminator is the profile exponent: n-stable for G^-1/2,
     drifting for S^-1/2.  This is the two-pole statement: preconditioning
     cures the chi = pi zero exactly and cannot touch the chi -> 0 chirp.
"""
import numpy as np
import pytest

from geovac.sturmian_sigma_law import sw_cross_block

S_KR = 2.0          # reduced separation k*R used throughout
M_QUAD = 100_001


def tridiag(n, lo=1.0, mid=2.0, hi=1.0):
    return (np.diag(mid * np.ones(n))
            + np.diag(hi * np.ones(n - 1), 1)
            + np.diag(lo * np.ones(n - 1), -1))


def inv_sqrt(M):
    ev, U = np.linalg.eigh(M)
    assert ev.min() > 0, "matrix is not positive definite"
    return U @ np.diag(ev ** -0.5) @ U.T


def bandwidth_for(M, tol):
    n = M.shape[0]
    total = np.linalg.norm(M)
    off = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    for b in range(n):
        if np.linalg.norm(M * (off > b)) <= tol * total:
            return b
    return n


def profile_exponent(M):
    n = M.shape[0]
    d = np.arange(1, n // 2)
    prof = np.array([np.abs(np.diag(M, k)).mean() for k in d])
    m = (d >= 3) & (prof > 0)
    return float(np.polyfit(np.log(d[m]), np.log(prof[m]), 1)[0])


def _G(n):
    C = sw_cross_block(S_KR, n, M=M_QUAD)
    P_is = inv_sqrt(tridiag(n))
    return C, P_is @ (np.eye(n) - C) @ P_is


# ------------------------------------------------------------ the preconditioner
def test_matching_polynomial_gives_an_exactly_tridiagonal_preconditioner():
    """g = 2 + 2 cos chi has g_0 = 2, g_1 = 1, g_{j>=2} = 0.

    In this basis the entry is g_{|a-b|} - g_{a+b}; since a, b >= 1 the Hankel
    index a+b is always >= 2, where g vanishes.  So the Hankel part contributes
    NOTHING and P is exactly tridiagonal -- which is what makes the
    preconditioner local and DST-diagonalizable.

    WRONG ANSWER REJECTED: that the Hankel part matters here (it would spoil
    both the tridiagonality and the closed-form spectrum).  The test builds the
    full Toeplitz-minus-Hankel matrix from the coefficients and compares it to
    an independently constructed tri(1,2,1), so a nonzero Hankel contribution
    fails the comparison.
    """
    g = {0: 2.0, 1: 1.0}
    for n in (6, 17, 40):
        built = np.array([[g.get(abs(a - b), 0.0) - g.get(a + b, 0.0)
                           for b in range(1, n + 1)] for a in range(1, n + 1)])
        assert np.allclose(built, tridiag(n), atol=0, rtol=0)
        assert np.count_nonzero(np.triu(built, 2)) == 0


def test_preconditioner_is_exactly_dst_diagonalized():
    """P V = V diag(2 + 2 cos(k pi/(n+1))) on the DST-I basis, in closed form.

    WRONG ANSWER REJECTED: the sign-flipped spectrum 2 - 2 cos(k pi/(n+1)),
    i.e. mistaking the zero's location for chi = 0 instead of chi = pi.  That
    is the natural error (the standard discrete Laplacian has -1 off-diagonals)
    and it would put the preconditioner's zero at the wrong end of the symbol,
    where it would not cancel anything.
    """
    for n in (16, 64):
        k = np.arange(1, n + 1)
        V = np.sqrt(2 / (n + 1)) * np.sin(np.outer(k, k) * np.pi / (n + 1))
        lam = 2 + 2 * np.cos(k * np.pi / (n + 1))
        assert np.abs(V @ V.T - np.eye(n)).max() < 1e-12
        assert np.abs(tridiag(n) @ V - V * lam).max() < 1e-11


# -------------------------------------------------------------- the conditioning
def test_preconditioned_conditioning_is_bounded_while_raw_grows():
    """cond(I-C) ~ n^2; cond(G) bounded and its increments -> 0.

    WRONG ANSWER REJECTED: that preconditioning merely reduces a constant while
    leaving the growth intact.  Asserted two ways that a constant-factor gain
    cannot satisfy: the raw condition number must QUADRUPLE per doubling (the
    n^2 law), while cond(G)'s successive increments must SHRINK by at least 2x
    per doubling and stay under 2.3 throughout.
    """
    ns = (20, 40, 80, 160)
    raw, pre = [], []
    for n in ns:
        C, G = _G(n)
        raw.append(np.linalg.cond(np.eye(n) - C))
        pre.append(np.linalg.cond(G))

    for a, b in zip(raw, raw[1:]):
        assert 3.5 < b / a < 4.5, f"raw growth {b/a:.2f} is not the n^2 law"
    assert all(p < 2.3 for p in pre), f"cond(G) not bounded: {pre}"
    incs = [b - a for a, b in zip(pre, pre[1:])]
    assert all(i > 0 for i in incs), "expected monotone approach from below"
    for a, b in zip(incs, incs[1:]):
        assert b < a / 2, f"cond(G) increments not collapsing: {incs}"


def test_whitening_preserves_the_generalized_spectrum():
    """X = P^-1/2 G^-1/2 satisfies X^T (I-C) X = I, so eigenvalues are unchanged.

    This is what makes the lever legitimate rather than a change of problem.

    WRONG ANSWER REJECTED: using P^-1/2 ALONE as the whitening (the tempting
    shortcut, since that is the factor carrying the DST).  It does not whiten --
    P^-1/2 (I-C) P^-1/2 = G != I -- and the generalized eigenvalues it returns
    are wrong.  Both halves are asserted.
    """
    n = 64
    C, G = _G(n)
    A = np.eye(n) - C
    rng = np.random.default_rng(60_12)
    H = rng.normal(size=(n, n))
    H = 0.5 * (H + H.T)

    X = inv_sqrt(tridiag(n)) @ inv_sqrt(G)
    assert np.abs(X.T @ A @ X - np.eye(n)).max() < 1e-8

    import scipy.linalg as sla
    ref = np.sort(sla.eigh(H, A, eigvals_only=True))
    via = np.sort(np.linalg.eigvalsh(X.T @ H @ X))
    assert np.allclose(ref, via, rtol=1e-6, atol=1e-8)


# ------------------------------------------------------------------- locality
def test_locality_gain_is_real_at_one_percent_and_erodes_at_one_permille():
    """The honest two-sided claim: big gain at 1e-2, eroding at 1e-3.

    WRONG ANSWER REJECTED: the overclaim that G^-1/2 has an n-INDEPENDENT
    bandwidth at every tolerance.  The second block asserts that the 1e-3
    bandwidth GROWS by more than 2x across n=32..128, which an n-independent
    bandwidth cannot do.  (The first block equally excludes the opposite
    overclaim, that preconditioning buys no locality at all.)
    """
    ns = (32, 64, 128)
    b_S_1, b_G_1, b_G_3 = [], [], []
    for n in ns:
        C, G = _G(n)
        b_S_1.append(bandwidth_for(inv_sqrt(np.eye(n) - C), 1e-2))
        b_G_1.append(bandwidth_for(inv_sqrt(G), 1e-2))
        b_G_3.append(bandwidth_for(inv_sqrt(G), 1e-3))

    # at 1e-2 the raw factor needs a fixed FRACTION of the matrix; G does not
    for n, b in zip(ns, b_S_1):
        assert b / n > 0.5, f"raw bandwidth fraction {b/n:.2f} at n={n}"
    assert b_G_1[-1] < 0.15 * ns[-1]
    assert b_G_1[-1] < 2 * b_G_1[0], f"1e-2 bandwidth for G grew: {b_G_1}"

    # at 1e-3 the advantage erodes -- do not let the paper overclaim
    assert b_G_3[-1] > 2 * b_G_3[0], f"1e-3 bandwidth did NOT grow: {b_G_3}"


def test_profile_exponent_is_n_stable_for_G_and_drifts_for_raw():
    """The clean discriminator, and the two-pole statement in one number.

    G^-1/2 keeps the chirp's own exponent (~ -5/4), n-independently, because
    preconditioning removed the chi = pi singularity and cannot touch chi -> 0.
    S^-1/2 keeps that singularity, so its profile exponent DRIFTS as n grows.

    WRONG ANSWER REJECTED: that both are n-stable (i.e. that preconditioning
    changed nothing qualitatively), and that both drift (i.e. that it fixed
    nothing).  Asserted as a stability bound on one and a drift bound on the
    other, so either collapse fails.
    """
    ns = (64, 128, 256)
    pG, pS = [], []
    for n in ns:
        C, G = _G(n)
        pG.append(profile_exponent(inv_sqrt(G)))
        pS.append(profile_exponent(inv_sqrt(np.eye(n) - C)))

    assert max(pG) - min(pG) < 0.08, f"G exponent not n-stable: {pG}"
    assert all(-1.35 < p < -1.05 for p in pG), f"G exponent off the chirp: {pG}"
    assert pS[-1] - pS[0] > 0.08, f"raw exponent did not drift: {pS}"
    assert all(p > q for p, q in zip(pS, pG)), "raw must decay SLOWER than G"


# ------------------------------------------------- transfer to a polyatomic block
R_OH, R_HH = 1.809, 2.862          # water, ~104.5 degrees, k = 1


def _water_A1(n):
    """A_1 block of water's three-center SW metric under C_2v.

    Only the H1 <-> H2 swap acts (O is on the axis), so A_1 holds BOTH the O
    functions and the symmetric H combination -- and with them the O <-> H
    coupling between symmetry-INEQUIVALENT centers that the group cannot remove.
    """
    C_OH = sw_cross_block(R_OH, n, M=M_QUAD)
    C_HH = sw_cross_block(R_HH, n, M=M_QUAD)
    I = np.eye(n)
    return np.block([[I, np.sqrt(2) * C_OH],
                     [np.sqrt(2) * C_OH.T, I + C_HH]])


def _null_direction_rotation():
    """At chi = pi every block symbol -> j0(0) = 1 whatever the separation, so the
    matrix symbol degenerates to rank-one all-ones.  For A_1 that is
    [[1, sqrt2], [sqrt2, 2]]: singular, trace 3, null direction fixed."""
    A_pi = np.array([[1.0, np.sqrt(2)], [np.sqrt(2), 2.0]])
    w, V = np.linalg.eigh(A_pi)
    assert abs(w[0]) < 1e-12 and abs(w[1] - 3.0) < 1e-12, f"symbol limit wrong: {w}"
    return V


def test_lever_transfers_to_water_A1_block():
    """The breach reaches the case the gerade lever fails: inequivalent centers.

    WRONG ANSWER REJECTED: that the diatomic result is a homonuclear accident and
    the polyatomic block still grows.  The raw column must show the n^2 law
    (quadrupling per doubling) while the preconditioned column stays bounded with
    increments collapsing -- a constant-factor gain satisfies neither.
    """
    ns = (6, 12, 24, 48)
    raw, pre = [], []
    for n in ns:
        A = _water_A1(n)
        Q = np.kron(_null_direction_rotation(), np.eye(n))
        P = np.zeros_like(A)
        P[:n, :n] = tridiag(n)
        P[n:, n:] = np.eye(n)
        P_is = inv_sqrt(P)
        raw.append(np.linalg.cond(A))
        pre.append(np.linalg.cond(P_is @ (Q.T @ A @ Q) @ P_is))

    for a, b in zip(raw, raw[1:]):
        assert 3.5 < b / a < 4.5, f"raw growth {b/a:.2f} is not the n^2 law"
    assert all(p < 50.0 for p in pre), f"preconditioned not bounded: {pre}"
    incs = [b - a for a, b in zip(pre, pre[1:])]
    for a, b in zip(incs, incs[1:]):
        assert b < a / 2, f"increments not collapsing: {incs}"


def test_water_needs_the_null_direction_rotation():
    """The control, and the load-bearing half: it is the ROTATION, not the
    preconditioning as such, that removes the growth.

    WRONG ANSWER REJECTED: "any band preconditioner fixes it."  Applying
    blockdiag(P, P) in the unrotated frame must leave the growth intact -- if it
    did not, the aligned result would prove nothing about the mechanism.
    """
    naive = []
    for n in (12, 24, 48):
        A = _water_A1(n)
        P = np.zeros_like(A)
        P[:n, :n] = tridiag(n)
        P[n:, n:] = tridiag(n)
        P_is = inv_sqrt(P)
        naive.append(np.linalg.cond(P_is @ A @ P_is))

    for a, b in zip(naive, naive[1:]):
        assert b / a > 3.0, f"unrotated control did NOT keep growing: {naive}"
    assert naive[-1] > 1e3, f"unrotated control is unexpectedly small: {naive}"


# ------------------------------------------------------- end-to-end resource pricing
def test_amplitude_floor_is_factorization_invariant():
    """eq:amplitude_floor. Any X with X^T A X = I has ||X|| = ||A^-1/2|| EXACTLY.

    Reason: X^T A X = I forces X = A^-1/2 U for some unitary U.  So the
    subnormalization floor of a block-encoding of the whitening cannot be lowered
    by ANY factorization, and the untreated route already attains it.

    WRONG ANSWER REJECTED: the hope that a cleverer factorization buys amplitude
    as well as depth.  Tested constructively against THREE genuinely different
    whitenings -- the symmetric inverse square root, the preconditioned
    P^-1/2 G^-1/2, and an inverse Cholesky factor (which is triangular, not
    symmetric, so it is not a disguised copy of the first) -- all of which must
    agree to near machine precision.
    """
    for n in (20, 40, 80):
        A = np.eye(n) - sw_cross_block(S_KR, n, M=M_QUAD)

        X_sym = inv_sqrt(A)
        P_is = inv_sqrt(tridiag(n))
        X_pre = P_is @ inv_sqrt(P_is @ A @ P_is)
        X_chol = np.linalg.inv(np.linalg.cholesky(A)).T

        for X in (X_sym, X_pre, X_chol):
            assert np.abs(X.T @ A @ X - np.eye(n)).max() < 1e-8, "not a whitening"

        norms = [np.linalg.norm(X, 2) for X in (X_sym, X_pre, X_chol)]
        assert max(norms) - min(norms) < 1e-8 * max(norms), (
            f"amplitude floor is NOT invariant at n={n}: {norms}")


def test_the_lever_buys_depth_and_costs_amplitude():
    """The honest pricing: depth goes flat, composed amplitude gets WORSE.

    Untreated: alpha ~ n, d_inv ~ n^2  (product ~ n^3).
    Preconditioned with G obtained by COMPOSING P^-1/2 with (I-C): d_inv flat,
    but alpha inherits ||P^-1/2||^2 ~ n^2  (product ~ n^2).

    WRONG ANSWER REJECTED: that preconditioning is a pure win.  The test asserts
    that the composed amplitude exponent is ~2, i.e. STRICTLY WORSE than the
    untreated ~1 -- so a reading in which the lever costs nothing fails here.
    It equally rejects the opposite overclaim by requiring the depth exponent to
    collapse to ~0.
    """
    ns = np.array([20, 40, 80, 160])
    a_naive, d_naive, a_comp, d_pre = [], [], [], []
    for n in ns:
        n = int(n)
        A = np.eye(n) - sw_cross_block(S_KR, n, M=M_QUAD)
        P_is = inv_sqrt(tridiag(n))
        G = P_is @ A @ P_is
        kA, kG = np.linalg.cond(A), np.linalg.cond(G)
        eps = 1.6e-3
        a_naive.append(np.linalg.norm(inv_sqrt(A), 2))
        d_naive.append(kA * np.log(kA / eps))
        a_comp.append(np.linalg.norm(P_is, 2) ** 2 * np.linalg.norm(A, 2))
        d_pre.append(kG * np.log(kG / eps))

    slope = lambda y: float(np.polyfit(np.log(ns), np.log(y), 1)[0])
    assert 0.9 < slope(a_naive) < 1.1, f"untreated alpha exponent {slope(a_naive):.2f}"
    assert 1.9 < slope(d_naive) < 2.3, f"untreated depth exponent {slope(d_naive):.2f}"
    assert 1.9 < slope(a_comp) < 2.1, f"composed alpha exponent {slope(a_comp):.2f}"
    assert abs(slope(d_pre)) < 0.05, f"preconditioned depth not flat: {slope(d_pre):.3f}"
    # the paid-for-nothing factor the open item would recover
    assert a_comp[-1] / a_naive[-1] > 20, "composition penalty vanished unexpectedly"


@pytest.mark.slow
def test_preconditioned_conditioning_stays_flat_at_large_basis():
    """One cutoff past the data the paper quotes -- the guard-asymptotics rule.

    WRONG ANSWER REJECTED: a bound that holds only on the n <= 160 window the
    paper tabulates.  n = 320 is outside it.
    """
    _, G = _G(320)
    assert np.linalg.cond(G) < 2.3
