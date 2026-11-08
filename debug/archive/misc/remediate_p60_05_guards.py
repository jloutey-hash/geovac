"""REMEDIATION 5/5 -- the two guards, written as their OWN pass.

CLAUDE.md Sec. 9: guard-writing is a separate, separately-reviewed activity.
Remediations 1-4 completed the content fixes; this pass only touches guards, and
each is reviewed by asking WHAT WRONG ANSWER WOULD THIS ACCEPT.

GUARD 1 -- test_water_needs_the_null_direction_rotation (code-A MATERIAL-1).
  The offered control is blockdiag(T, T) = I2 (x) T, which COMMUTES exactly with
  the rotation Q = V (x) I.  PM-verified: ||PQ - QP|| = 0 and
  ||Q^T P^-1/2 Q - P^-1/2|| = 4e-15, so rotated and unrotated give IDENTICAL
  spectra.  The control's discriminating power over the variable it is NAMED
  after is exactly zero, and planting the rotation into it does not fire.
  Replacement asserts two things instead:
    (i) the commutation, so the blindness is documented and cannot be
        reintroduced by someone "restoring" the old control;
    (ii) the genuinely discriminating counterfactual -- the SELECTIVE
        preconditioner blockdiag(T, I) applied in the UNROTATED frame.
        PM-measured: 1689 -> 20301 -> 290329 -> 4437303 over n = 12..96, i.e. it
        grows FASTER than untreated (699 -> 41700).  Applying the right
        preconditioner in the wrong frame is actively harmful, which is what
        makes the frame load-bearing.

GUARD 2 -- the (b) leg of the K=452 state-dependence test (code-B M4).
  Current: norm_all = |M|.sum(); then assert |M|.sum() == norm_all on an
  unmodified M.  That is x == x; the loop body discards its own result.  Nothing
  is excluded and no wrong answer is reachable.
  Replacement builds the per-state alternative the docstring names -- the SAME
  configuration family rebuilt at a state-dependent reference scale, via
  Config's pk_ref, which no test has ever varied (code-B M5) -- and asserts that
  it DOES move the 1-norm while the pipeline's single M does not.  That makes
  "a per-state matrix would break the equality" a tested statement rather than
  an asserted one, and closes M5's frozen-parameter gap in the same guard.

Idempotent.
"""
from __future__ import annotations

import sys

T_PRE = "tests/test_paper60_preconditioner.py"
T_LAD = "tests/test_paper60_resource_ladder.py"
MARKER1 = "commutes with the rotation"
MARKER2 = "per-state alternative is REACHABLE"

OLD1 = '''def test_water_needs_the_null_direction_rotation():
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
'''

NEW1 = '''def test_uniform_band_preconditioner_commutes_with_the_rotation():
    """Why the OLD control could not work, pinned so it cannot come back.

    The uniform control blockdiag(T, T) is I2 (x) T and the rotation is
    V (x) I, so they COMMUTE: rotated and unrotated frames give identical
    spectra, and a control built that way is blind to the variable it is named
    after.  /qa paper_60 FULL 2026-09-12 found the old guard asserting exactly
    that comparison; planting the rotation into it did not fire.

    WRONG ANSWER REJECTED: "blockdiag(T, T) unrotated is a valid control for the
    rotation."  If this assertion ever fails, the uniform control has become
    frame-sensitive and the old guard could be revived; it is here so that
    cannot happen silently.
    """
    n = 24
    A = _water_A1(n)
    P = np.zeros_like(A)
    P[:n, :n] = tridiag(n)
    P[n:, n:] = tridiag(n)
    Q = np.kron(_null_direction_rotation(), np.eye(n))
    P_is = inv_sqrt(P)
    assert np.linalg.norm(P @ Q - Q @ P) < 1e-10, "uniform P must commute with Q"
    assert np.linalg.norm(Q.T @ P_is @ Q - P_is) < 1e-10
    c_un = np.linalg.cond(P_is @ A @ P_is)
    c_rot = np.linalg.cond(P_is @ (Q.T @ A @ Q) @ P_is)
    assert abs(c_un - c_rot) / c_un < 1e-9, (
        f"uniform control IS frame-sensitive ({c_un} vs {c_rot}) -- the old "
        f"control would then have had discriminating power after all")


def test_water_needs_the_null_direction_rotation():
    """The load-bearing half, with a control that can actually fail.

    WRONG ANSWER REJECTED: "the frame does not matter -- any selective
    preconditioner fixes it."  Applying the SELECTIVE preconditioner
    blockdiag(T, I) in the UNROTATED frame must not merely fail to help: it must
    be worse than doing nothing, because it sharpens the wrong direction.  That
    is what makes the alignment, rather than the preconditioning, load-bearing.
    """
    raw, naive_sel = [], []
    for n in (12, 24, 48, 96):
        A = _water_A1(n)
        P = np.zeros_like(A)
        P[:n, :n] = tridiag(n)
        P[n:, n:] = np.eye(n)
        P_is = inv_sqrt(P)
        raw.append(np.linalg.cond(A))
        naive_sel.append(np.linalg.cond(P_is @ A @ P_is))

    for a, b in zip(naive_sel, naive_sel[1:]):
        assert b / a > 4.5, f"unrotated selective control did not blow up: {naive_sel}"
    assert naive_sel[-1] > 20.0 * raw[-1], (
        f"unrotated selective preconditioning must be WORSE than untreated: "
        f"{naive_sel[-1]:.3e} vs raw {raw[-1]:.3e}")
'''

OLD2 = '''    # (b) one matrix serves both roots -- the 1-norm cannot be per-state
    norm_all = float(np.abs(M).sum())
    assert norm_all > 0
    for k in (0, 1):
        _ = -p[k] ** 2 / 2.0
        assert float(np.abs(M).sum()) == norm_all
'''

NEW2 = '''    # (b) one matrix serves both roots -- the 1-norm cannot be per-state.
    # The old form here asserted |M|.sum() == |M|.sum() on an unmodified M,
    # i.e. x == x, and excluded nothing (/qa paper_60 FULL 2026-09-12).  The
    # exclusion only means something if the per-state alternative is REACHABLE,
    # so build it: the same family at a state-dependent reference scale, via
    # Config's pk_ref -- a parameter no test had ever varied.
    norm_all = float(np.abs(M).sum())
    assert norm_all > 0
    per_state = []
    for k in (0, 1):
        cfgs_k = [SS.Config(c.l, c.na, c.nb, pk_ref=float(p[k]))
                  for c in cfgs[:40]]
        per_state.append(float(np.abs(SS.build_M(cfgs_k, Z=Z)).sum()))
    assert abs(per_state[0] - per_state[1]) / per_state[0] > 0.05, (
        f"a per-state rebuild must MOVE the 1-norm, else this leg excludes "
        f"nothing: {per_state}")
    sub = float(np.abs(SS.build_M(list(cfgs[:40]), Z=Z)).sum())
    for k in (0, 1):
        assert float(np.abs(M).sum()) == norm_all
    assert abs(sub - per_state[0]) / sub > 0.05, (
        "the pipeline's own matrix must differ from the per-state rebuild")
'''


def main() -> int:
    with open(T_PRE, encoding="utf-8") as fh:
        a = fh.read()
    with open(T_LAD, encoding="utf-8") as fh:
        b = fh.read()
    done = 0
    if MARKER1 in a:
        print("  guard 1 ALREADY APPLIED")
        done += 1
    elif a.count(OLD1) != 1:
        print(f"  guard 1 anchor count={a.count(OLD1)}; ABORT")
        return 2
    if MARKER2 in b:
        print("  guard 2 ALREADY APPLIED")
        done += 1
    elif b.count(OLD2) != 1:
        print(f"  guard 2 anchor count={b.count(OLD2)}; ABORT")
        return 3
    if done == 2:
        return 1
    if MARKER1 not in a:
        with open(T_PRE, "w", encoding="utf-8") as fh:
            fh.write(a.replace(OLD1, NEW1))
        print("  ok  guard 1 (water control) replaced + commutation pin added")
    if MARKER2 not in b:
        with open(T_LAD, "w", encoding="utf-8") as fh:
            fh.write(b.replace(OLD2, NEW2))
        print("  ok  guard 2 (tautological 1-norm leg) replaced")
    return 0


if __name__ == "__main__":
    sys.exit(main())
