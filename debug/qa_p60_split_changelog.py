"""Split the coverage closure into its own v5.10.17 entry.

v5.10.16 is already committed and TAGGED (c2e8423 / tag v5.10.16).  Renaming its
heading to v5.10.17 -- which the previous script did -- would leave that tag
pointing at a released state with no CHANGELOG entry of its own, and would
rewrite history that has already been stamped.  The honest shape is two entries:
v5.10.16 as it was cut, including the "Owed" section that was TRUE at that
commit, and v5.10.17 for the pass that discharges it.
"""
import io

C = "CHANGELOG.md"
c = io.open(C, encoding="utf-8").read()

# 1. put the v5.10.16 heading back.
assert "## [v5.10.17] - 2026-09-11" in c
c = c.replace("## [v5.10.17] - 2026-09-11", "## [v5.10.16] - 2026-09-11", 1)

# 2. pull the closure text back out of the v5.10.16 entry and restore "Owed",
#    which was an accurate statement of that commit's state.
start = c.index("### The coverage LARGE, closed as its own pass")
end = c.index("### Also", start)
closure = c[start:end]
OWED = """### Owed

**The second LARGE is a coverage gap, and it is owed as its own pass:** six abstract-level `[MEASURED]` families (the span-deficit pair, the free-scale matched set, the K=452 state pair, posing-cost roots 2–3, the state-prep overlap, the floor brackets) have driver-only backing in the prunable `debug/` tree. The reviewer independently reproduced every one, so the exposure is regression protection rather than correctness — and §9 requires guard-writing to be separate, separately-reviewed work, so it is not bundled here.

"""
c = c[:start] + OWED + c[end:]

# 3. add v5.10.17 above v5.10.16.
ENTRY = ("## [v5.10.17] - 2026-09-11\n\n"
         "**The coverage LARGE from the v5.10.16 DELTA, discharged as its own pass.** "
         "Six abstract-level `[MEASURED]` families move from prunable-driver backing to "
         "tracked `geovac/` recomputation. New file `tests/test_paper60_resource_ladder.py` "
         "(7 tests, all `@pytest.mark.slow`, 14 min).\n\n"
         + closure.replace("### The coverage LARGE, closed as its own pass\n\n", "", 1)
         + "\n")
anchor = "## [v5.10.16] - 2026-09-11"
c = c.replace(anchor, ENTRY + anchor, 1)

io.open(C, "w", encoding="utf-8").write(c)
print("CHANGELOG: v5.10.16 restored with its Owed section; v5.10.17 added above it")
