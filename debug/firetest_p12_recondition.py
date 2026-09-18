"""Fire test for tests/test_paper12_recondition.py guards (Sec. 9).

Each check plants the SPECIFIC wrong answer a guard is meant to reject and
asserts the guard's numeric condition FAILS on it.  A guard whose planted wrong
answer cannot be named, or that passes it, is not a guard.

Cheap by construction: no mpf V_ee build.  The two @slow physics guards
(breaks-the-wall, climb) are fired here on their CHEAP half -- the monomial wall
must actually exist -- and their expensive half (re-based climbs past it) is
pinned by the module run recorded in the claim matrix; the guard would fail if
the re-based path returned the monomial's -64 Ha instead of 99.2%.
"""
from __future__ import annotations

import numpy as np
import mpmath as mp

from geovac import prolate_general_m as pg
from geovac import prolate_recondition as pr

mp.mp.dps = 30
fired = []


def check(name: str, condition: bool, detail: str) -> None:
    tag = "FIRES" if condition else "DID NOT FIRE"
    fired.append((name, condition, detail))
    print(f"  [{tag:12s}] {name}: {detail}")


# --- Guard: test_factored_change_of_basis_equals_the_dense_one (thresh 1e-25).
#     Wrong answer: computing C M C  (missing the transpose on the right).
def fire_factored_cob() -> None:
    j_max, l_max, mu_max = 2, 2, 1
    Nr = (j_max + 1) ** 2
    Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1) if (l + m) % 2 == 0])
    Nmu = mu_max + 1
    N = Nmu * Nr * Na
    Tr_list, Ta_list = pr._transforms_per_mu("gegenbauer", j_max, l_max, mu_max, 1.0)
    rng = np.random.default_rng(0)
    M = np.empty((N, N), object)
    for i in range(N):
        for j in range(i, N):
            v = mp.mpf(float(rng.standard_normal()))
            M[i, j] = M[j, i] = v
    C = np.zeros((N, N), object)
    C[:] = mp.mpf(0)
    for mu in range(Nmu):
        sl = slice(mu * Nr * Na, (mu + 1) * Nr * Na)
        C[sl, sl] = np.kron(Tr_list[mu], Ta_list[mu])
    correct = C @ M @ C.T
    wrong = C @ M @ C                       # missing transpose
    err = max(abs(correct[i, j] - wrong[i, j]) for i in range(N) for j in range(N))
    # the guard's threshold is 1e-25; the planted transpose bug must exceed it
    check("factored_cob==dense", err > mp.mpf(10) ** -25,
          f"missing-transpose deviation = {float(err):.2e} (guard thresh 1e-25)")


# --- Guard: test_normalized_solve_recovers_what_norm_spread_hides.
#     Wrong answer: a solve that does NOT rescale to unit norm, so it discards
#     variational content hidden by norm spread and returns E far above the true
#     lowest eigenvalue.  The guard asserts E_norm ~ E_true; the un-normalized
#     solve violates that.
def fire_norm_spread() -> None:
    from scipy.linalg import eigh as _eigh
    rng = np.random.default_rng(1)
    n = 12
    B = rng.standard_normal((n, n))
    Shat = B @ B.T + n * np.eye(n)
    Dh = np.sqrt(np.diag(Shat)); Shat = (Shat / Dh[:, None]) / Dh[None, :]
    A = rng.standard_normal((n, n)); Hhat = 0.5 * (A + A.T)
    E_true = float(_eigh(Hhat, Shat, eigvals_only=True)[0])
    d = np.geomspace(1.0, 1e8, n)
    S_o = np.empty((n, n), object)
    for i in range(n):
        for j in range(n):
            S_o[i, j] = mp.mpf(float(Shat[i, j] * d[i] * d[j]))
    # the WRONG solve: plain float64 canonical orth on the un-normalized S_o
    Sf = np.array([[float(S_o[i, j]) for j in range(n)] for i in range(n)])
    Hf = np.array([[float(Hhat[i, j] * d[i] * d[j]) for j in range(n)] for i in range(n)])
    w, U = np.linalg.eigh(0.5 * (Sf + Sf.T))
    keep = w > 1e-11 * w[-1]
    X = U[:, keep] / np.sqrt(w[keep])
    E_naive = float(np.linalg.eigvalsh(X.T @ (0.5 * (Hf + Hf.T)) @ X)[0])
    # the guard requires |E_norm - E_true| < 1e-8; the un-normalized solve is far
    check("normalized_solve beats norm-spread", abs(E_naive - E_true) > 1e-8,
          f"un-normalized solve gives {E_naive:.6f} vs true {E_true:.6f} "
          f"(off by {abs(E_naive - E_true):.2e}); normalization is load-bearing")


# --- Guard: test_orthogonal_families_reduce_to_the_mu0_reference.
#     Wrong answer: Gegenbauer with lam = mu (not mu + 1/2) -> at mu=0, lam=0,
#     which is NOT Legendre.
def fire_reduction() -> None:
    bad = pr.gegenbauer_coeffs(2, mp.mpf(0), 4)      # lam = mu = 0 (wrong)
    legendre = pr.legendre_coeffs(2, 4)
    dev = max(abs(bad[i] - legendre[i]) for i in range(4))
    check("gegenbauer(lam=mu+1/2)==Legendre@mu0", dev > mp.mpf(10) ** -30,
          f"lam=mu (wrong) deviates from Legendre by {float(dev):.2e}")


# --- Guard: test_rebasing_breaks_the_conditioning_wall (cheap half).
#     The wall must EXIST: monomial (3,3,1) at alpha=1.0 must be non-variational.
#     If it were variational, the "re-basing fixes a broken solve" contrast is
#     vacuous.  (Fast grid engine.)
def fire_wall_exists() -> None:
    alpha = 1.0
    basis = pg.generate_basis(3, 3, 1, alpha)
    mom = pg.Moments(2.0 * alpha, 6 * 3 + 6 * 3 + 20)
    grid = pg.XiGrid(alpha)
    s, h1 = pg.one_body(basis, pg.R_DEFAULT, 1.0, mom)
    v = pg.vee_matrix(basis, pg.R_DEFAULT, grid, l_neumann=14)
    h = h1 + v + (1.0 / pg.R_DEFAULT) * s
    e_mono, _nk, _tot = pg.solve_generalized(h, s)
    # guard asserts e_mono < E_EXACT (broken); it fires (the wall is real) iff so
    check("monomial-wall-exists @ (3,3,1)", e_mono < pg.E_EXACT,
          f"monomial E = {e_mono:.3f} Ha (< exact {pg.E_EXACT}); the wall is real")


# --- Guard: test_climb_reaches_chemical_accuracy_past_the_cap, the abs() half.
#     Wrong answer: E = -1.00, i.e. NO binding at all (0% of D_e).  `err_mha` is
#     SIGNED -- (E_exact - E)*1000 -- so every variational result is negative and
#     the original `err_mha < 1.6` accepted this; the repaired `abs(...) < 1.6`
#     must reject it.  Asserts BOTH halves, so it also records that the guard was
#     genuinely dead rather than merely loose.
def fire_abs_err_mha() -> None:
    E_bad = -1.00
    err_bad = (pr.E_EXACT - E_bad) * 1000.0
    old_accepts = err_bad < 1.6
    new_rejects = not (abs(err_bad) < 1.6)
    check("abs(err_mha) enforces chemical accuracy", old_accepts and new_rejects,
          f"E=-1.00 (0% of D_e, err={err_bad:+.3f} mHa): unsigned guard accepted "
          f"it ({old_accepts}), abs() guard rejects it ({new_rejects})")


# --- Guard: test_direct_engine_matches_the_mpf_reference_one_body.
#     Wrong answer: a one-body engine hardcoded to ONE basis family (which the
#     debug predecessor was -- gegenbauer only), compared against the other
#     family's mpf reference.  recondition_energy's DEFAULT family is the other
#     one, so this is the mismatch that would ship silently.
def fire_one_body_family() -> None:
    alpha = 1.0
    _Sd, Hd = pr.build_one_body_direct(2, 2, 1, alpha, basis="laguerre_legendre")
    with mp.workdps(pr.DEFAULT_DPS):
        idx = pr._product_index(2, 2, 1)
        fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        A = pr.ngm._mono_moments(2.0 * alpha, 6 * 2 + 6 * 3 + 20)
        _S, H1 = pr.one_body_mp(fns, alpha, pr.R_DEFAULT, A)
        Na = len([1 for l in range(3) for m in range(3) if (l + m) % 2 == 0])
        Tr, Ta = pr._transforms_per_mu("gegenbauer", 2, 2, 1, alpha)
        Hg_wrong = pr._to_f64(pr._factored_cob(H1, 2, 9, Na, Tr, Ta))
    scale = float(np.abs(Hd - Hg_wrong).max() / np.abs(Hg_wrong).max())
    check("one-body basis family", scale > 1e-13,
          f"wrong-family reference deviates {scale:.2e} scale-relative "
          f"(guard bar 1e-13)")


# --- Guard: test_direct_and_mpf_engines_agree_end_to_end.
#     Wrong answer: drop the nuclear-repulsion shift S_o/R from the direct
#     assembly H_o = H1_o + cob(V) + S_o/R.  The direct path adds the three
#     pieces separately (relying on linearity of the change of basis) where the
#     mpf path re-bases their sum, so a dropped term is the characteristic error.
def fire_missing_nuclear_shift() -> None:
    alpha = 1.0
    S_o, H1_o = pr.build_one_body_direct(2, 2, 1, alpha, basis="laguerre_legendre")
    Tr, Ta = pr._transforms_per_mu("laguerre_legendre", 2, 2, 1, alpha)
    idx = pr._product_index(2, 2, 1)
    fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
    with mp.workdps(pr.DEFAULT_DPS):
        V = pr.vee_mp(fns, alpha, pr.R_DEFAULT, 2 * 2 + 4 * 1 + 10)
        V_o = pr._to_f64(pr._factored_cob(V, 2, 9, 5, Tr, Ta))
    E_good = pr._normalized_solve(S_o, H1_o + V_o + (1.0 / pr.R_DEFAULT) * S_o)[0]
    E_bad = pr._normalized_solve(S_o, H1_o + V_o)[0]
    check("Sf*S_o nuclear shift", abs(E_good - E_bad) > 1e-9,
          f"dropping S_o/R moves E by {abs(E_good - E_bad):.4f} Ha "
          f"({E_good:.7f} -> {E_bad:.7f}); guard bar 1e-9")


if __name__ == "__main__":
    print("=== fire test: tests/test_paper12_recondition.py guards ===")
    fire_factored_cob()
    fire_norm_spread()
    fire_reduction()
    fire_wall_exists()
    fire_abs_err_mha()
    fire_one_body_family()
    fire_missing_nuclear_shift()
    n_fire = sum(1 for _, c, _ in fired if c)
    print(f"\n{n_fire}/{len(fired)} guards FIRE on their planted wrong answer.")
    if n_fire != len(fired):
        raise SystemExit("A guard did NOT fire; it does not discriminate its wrong answer.")
