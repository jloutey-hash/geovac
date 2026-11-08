"""Retry: '### Gates' is not unique across the CHANGELOG. Anchor on the v5.11.0
entry's own Gates paragraph. The renumber already applied; this finishes the
caveat fix and the water subsection."""
from pathlib import Path

C = Path("CHANGELOG.md")
s = C.read_text(encoding="utf-8")

OLD = ("**Adding a QA criterion is a gate change and therefore a candidate minor "
       "(v5.11.0) — PI call; bumped as a patch by default per Sec. 9.**")
NEW = ("**Adding a QA criterion is a gate change, so this ships as a minor (v5.11.0) "
       "rather than a patch — PI-confirmed 2026-09-12.** Under the Sec. 9 rule a moved "
       "second number should tell a reader that something corpus-significant happened "
       "without their having to read the entry; a new QA criterion and a breached wall "
       "both qualify.")
if OLD in s:
    s = s.replace(OLD, NEW)
    print("  version caveat resolved")
else:
    assert NEW in s, "neither caveat form present"
    print("  version caveat already resolved")

WATER = """### The breach reaches the polyatomic case — and the reason is structural

The scope caveat this entry was about to ship with is now measured, and it went the good way. Water's `A_1` block is the case Paper 60 records as defeating the gerade lever: `O` sits on the `C_2` axis so the only symmetry action is the `H <-> H` swap, leaving the `O <-> H` coupling between symmetry-**inequivalent** centers inside the totally-symmetric block where the ground state lives.

It works, because **the degeneracy's direction is geometry-independent**. At `chi = pi` every block symbol tends to `j0(0) = 1` whatever the separation, so for `M` centers the matrix symbol degenerates to the rank-one all-ones matrix and its null space has dimension `M-1` — a fixed subspace, not one that moves with the geometry or the basis. For `A_1` that limit is `[[1, sqrt2], [sqrt2, 2]]`: singular, trace 3, null direction `v ~ (sqrt2, -1)`. Rotating the block space by `v` and applying `tri(1,2,1)` to that component alone:

| `N` | `cond(A_1)` raw | preconditioned |
|--:|--:|--:|
| 12 | 183.0 | 38.45 |
| 48 | 2696.1 | 43.62 |
| 192 | 41699.7 | **44.06** |

The raw column independently reproduces the paper's `N^1.97` (`N^1.96` here). The preconditioned column is bounded with increments collapsing `3.94, 1.24, 0.35, 0.09`. The constant is larger than the diatomic `2.23`, but the *growth* — the thing that makes the metric penalty basis-dependent — is gone.

**The control is the load-bearing half.** Naive `blockdiag(P, P)` without the rotation leaves the growth intact (`2766 -> 42008` over the same range), so it is the alignment onto the symbol's null direction doing the work, not preconditioning as such. Both are asserted in the backing test, and the control's fire test is worth recording: the first plant tried (swapping one block for the identity) **did not fire**, correctly — it left the frame unrotated, so it never tested the guard's subject. `fire_test.py` reported exactly that, and the real plant (making the control secretly the aligned construction) fires. A fire test that fails is the tool working.

"""

A = "### Gates\n\nC10 / C21 / C16 / C22 / C14 / latex-escapes / internal-titles / inline-arxiv PASS"
assert s.count(A) == 1, f"anchor matched {s.count(A)} times"
s = s.replace(A, WATER + A)
C.write_text(s, encoding="utf-8")
print("  water subsection inserted")

# CLAUDE.md Sec.2 bullet
M = Path("CLAUDE.md")
t = M.read_text(encoding="utf-8")
OLDB = ("band-Toeplitz preconditioner, cond -> 2.23 flat. Locality capped by the chirp. "
        "New C23 gate.")
NEWB = ("band-Toeplitz preconditioner bounds cond, diatomic AND water A_1. Locality capped "
        "by the chirp. New C23 gate.")
if OLDB in t:
    M.write_text(t.replace(OLDB, NEWB), encoding="utf-8")
    print("  CLAUDE.md Sec.2 bullet updated")
else:
    print("  CLAUDE.md Sec.2 bullet already updated")
