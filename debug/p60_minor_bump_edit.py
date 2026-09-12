"""PI-approved minor bump v5.10.19 -> v5.11.0 (C23 is a gate change), plus the
water-transfer subsection folded into the release entry."""
from pathlib import Path

for path in ("CHANGELOG.md", "CLAUDE.md", "docs/walls/register.md"):
    P = Path(path)
    s = P.read_text(encoding="utf-8")
    assert "v5.10.19" in s, path
    P.write_text(s.replace("v5.10.19", "v5.11.0"), encoding="utf-8")
    print(f"  renumbered  {path}")

# the C23 paragraph's version caveat is now resolved
C = Path("CHANGELOG.md")
s = C.read_text(encoding="utf-8")
OLD = ("**Adding a QA criterion is a gate change and therefore a candidate minor "
       "(v5.11.0) — PI call; bumped as a patch by default per Sec. 9.**")
NEW = ("**Adding a QA criterion is a gate change, so this ships as a minor (v5.11.0) "
       "rather than a patch — PI-confirmed 2026-09-12.** Under the Sec. 9 rule a moved "
       "second number should tell a reader that something corpus-significant happened "
       "without their having to read the entry; a new QA criterion and a breached wall "
       "both qualify.")
assert s.count(OLD) == 1
s = s.replace(OLD, NEW)

# ---- water transfer subsection, inserted before the Gates section -------------
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
A = "### Gates\n"
assert s.count(A) == 1
s = s.replace(A, WATER + A)
C.write_text(s, encoding="utf-8")
print("  CHANGELOG: water subsection + version caveat resolved")

# ---- CLAUDE.md Sec. 2 bullet gains the transfer -------------------------------
M = Path("CLAUDE.md")
s = M.read_text(encoding="utf-8")
OLDB = ("- **Composition wall is 3 axes; conditioning BREACHED (2026-09-12, v5.11.0):** "
        "band-Toeplitz preconditioner, cond -> 2.23 flat. Locality capped by the chirp. "
        "New C23 gate. See CHANGELOG v5.11.0.\n")
NEWB = ("- **Composition wall is 3 axes; conditioning BREACHED (2026-09-12, v5.11.0):** "
        "band-Toeplitz preconditioner bounds cond, diatomic AND water A_1. Locality capped "
        "by the chirp. New C23 gate. See CHANGELOG v5.11.0.\n")
assert s.count(OLDB) == 1
M.write_text(s.replace(OLDB, NEWB), encoding="utf-8")
print(f"  CLAUDE.md: Sec.2 bullet updated ({len(NEWB.split())} words)")
