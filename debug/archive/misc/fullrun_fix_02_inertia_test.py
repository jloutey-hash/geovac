"""FULL-run remediation, cert-blocker 2: give C8.16 (the Sylvester-inertia
variational bound + root-by-root correspondence) a STANDALONE backing test.

The FULL run's completeness dimension FAILed partly on this: the [INTERNAL
THEOREM] on which every excited-state number in section 4 rests was verified
1.4e-17 only in the DoD, in no `tests/` file (matrix row 568, declared OPEN,
raised to PI). The code reviewer confirmed its CONSEQUENCES (root-by-root k=0,1)
and INGREDIENTS (T=I-S/2, M=ZW-G, W diagonal) are individually tested, so the
identity is "implied by tested legs" -- but the inertia-COUNTING step itself,
the load-bearing bridge, had no direct test. This adds one.

Written as its own pass (guard-writing rule), reviewed by naming the wrong
answer it rejects rather than asking whether it passes.

THE THEOREM (C8.16). With T = I - S/2 and M = Z*W - G (W diagonal, C8.15):
    H(lam) + (1/2) lam^2 S  =  lam ( lam*I - M )      [algebraic identity]
so, S being positive definite, Sylvester's law of inertia gives, for every lam:
    #{ generalized roots of H(lam)C = E S C that lie below -lam^2/2 }
      ==  #{ eigenvalues of M that exceed lam }.
At lam = p_k (the k-th largest eigenvalue of M) this forces the k-th pencil
root to equal E_iso(k) = -p_k^2/2 -- the root-by-root correspondence.

WRONG ANSWERS REJECTED, named:
  (1) "E_iso is asserted to be the lowest root but the pencil-to-standard-
      eigenproblem reduction does not actually hold" -- caught by the algebraic
      identity leg (fires if T != I - S/2 or M != ZW - G).
  (2) "the root-by-root correspondence is assumed, not forced" -- caught by the
      inertia-count leg: if the count below -lam^2/2 did NOT track #{eig(M)>lam},
      the correspondence would be a coincidence. Fire-tested by perturbing M:
      a shifted M gives a count mismatch.

Idempotent. Appends one test to tests/test_paper60_scale_lock.py.
"""
from __future__ import annotations

import sys

T = "tests/test_paper60_scale_lock.py"

NEW_TEST = '''

def test_c5_inertia_bound_and_root_by_root():
    """C8.16 -- the variational bound and root-by-root correspondence are FORCED
    by Sylvester inertia, not assumed. Standalone backing for matrix row 568.

    Three legs. The identity and the inertia count are the content; the
    root-by-root check is their consequence, pinned here so a future edit that
    breaks the bridge cannot pass on the consequence alone.

    WRONG ANSWER REJECTED (1): the pencil H(lam)C = E S C is claimed to reduce to
    the standard eigenproblem lam(lam*I - M), but does not. Leg A fires if
    T != I - S/2 (the kinetic identity) or if the standard-eigenproblem M does
    not equal Z*W - G.
    WRONG ANSWER REJECTED (2): the root-by-root correspondence E_iso(k) = -p_k^2/2
    is assumed rather than forced. Leg B fires unless, for every lam tested, the
    number of pencil roots below -lam^2/2 equals #{eig(M) > lam}. Fire-tested by
    shifting M: the count then mismatches (asserted inline below).
    """
    c = case(6, 1)
    n = c.M.shape[0]
    I = np.eye(n)

    # ---- Leg A: the algebraic identity, to machine precision ---------------
    # T = I - S/2  (the kinetic identity), and M = Z*W - G  (W diagonal).
    assert np.linalg.norm(c.T - (I - 0.5 * c.S)) < 1e-9, (
        f"kinetic identity T = I - S/2 broken: "
        f"||T - (I - S/2)|| = {np.linalg.norm(c.T - (I - 0.5 * c.S)):.2e}")
    M_from_WG = Z * c.W - c.G
    assert np.linalg.norm(M_from_WG - c.M) < 1e-6, (
        f"standard-eigenproblem M != Z*W - G (route cross-check): "
        f"||ZW - G - M|| = {np.linalg.norm(M_from_WG - c.M):.2e}")
    # hence H(lam) + (1/2) lam^2 S == lam(lam*I - M), for arbitrary lam:
    for lam in (1.0, 2.386, 5.0):
        H = lam ** 2 * c.T + lam * (-Z * c.W + c.G)
        lhs = H + 0.5 * lam ** 2 * c.S
        rhs = lam * (lam * I - c.M)
        assert np.linalg.norm(lhs - rhs) / (n * lam ** 2) < 1e-9, (
            f"identity H(lam) + lam^2 S/2 = lam(lam I - M) broken at lam={lam}")

    # ---- Leg B: the inertia count (the load-bearing bridge) ----------------
    X = _whiten(c.S)                       # S = (X X^T)^{-1}; whitened pencil
    pM = np.sort(np.linalg.eigvalsh(c.M))[::-1]   # descending eigenvalues of M

    def pencil_roots(lam: float) -> np.ndarray:
        H = lam ** 2 * c.T + lam * (-Z * c.W + c.G)
        return np.linalg.eigvalsh(X.T @ H @ X)

    def inertia_ok(Mmat: np.ndarray) -> list:
        eM = np.sort(np.linalg.eigvalsh(Mmat))[::-1]
        out = []
        for lam in (0.8, 1.5, 2.386, 3.0, 4.5):
            H = lam ** 2 * c.T + lam * (-Z * c.W + c.G)
            roots = np.linalg.eigvalsh(X.T @ H @ X)
            below = int(np.sum(roots < -0.5 * lam ** 2 - 1e-9))
            above = int(np.sum(eM > lam + 1e-9))
            out.append((lam, below, above))
        return out

    for lam, below, above in inertia_ok(c.M):
        assert below == above, (
            f"inertia count broken at lam={lam}: #{{roots < -lam^2/2}}={below} "
            f"but #{{eig(M) > lam}}={above}")

    # fire-test the inertia leg: a SHIFTED M must break the count (proves the
    # count is tracking M, not passing vacuously)
    mism = [(lam, b, a) for lam, b, a in inertia_ok(c.M + 0.5 * I) if b != a]
    assert mism, ("inertia leg is vacuous: shifting M by 0.5*I did not change "
                  "any count -- the test would pass for the wrong M")

    # ---- Leg C: the consequence -- level k of H(p_k) equals E_iso(k) -------
    for k in (0, 1, 2, 3):
        lam_k = math.sqrt(-2 * c.E_iso_k(k))
        root_k = float(c.levels(lam_k, k + 1)[k])
        assert abs(root_k - c.E_iso_k(k)) < 1e-9, (
            f"root-by-root broken at k={k}: pencil root {root_k:.9f} "
            f"vs E_iso(k) {c.E_iso_k(k):.9f}")
'''


def main() -> int:
    with open(T, encoding="utf-8") as fh:
        t = fh.read()
    if "def test_c5_inertia_bound_and_root_by_root" in t:
        print("already applied")
        return 0
    # need _whiten in scope
    if "_whiten" not in t:
        # import it alongside the existing var_* imports
        anchor = "from geovac.sturmian_variational import (                             # noqa: E402\n    build, var_energy, var_levels)"
        if t.count(anchor) != 1:
            print(f"  MISS import anchor count={t.count(anchor)}")
            return 3
        t = t.replace(anchor,
                      "from geovac.sturmian_variational import (                             # noqa: E402\n"
                      "    build, var_energy, var_levels, _whiten)")
    with open(T, "w", encoding="utf-8") as fh:
        fh.write(t + NEW_TEST)
    print("  ok    test_c5_inertia_bound_and_root_by_root appended (+ _whiten import)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
