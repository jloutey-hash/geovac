r"""v5.11.21 -> v5.12.0, the standing qa-skill minor-bump rule, and the
Lorentzian archive folded into the same entry.

PI direction 2026-09-14: "It should be minor. Really any time we touch the qa
skill that's gonna be a minor bump."  That converts a per-event PI call into a
standing rule, so it is recorded in Sec. 9 rather than left to judgement.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

CL = "CHANGELOG.md"
CM = "CLAUDE.md"

RULE_OLD = """So: patch (x.y.Z) is the standing default for everything — bug fixes, documentation, completed diagnostic arcs, paper updates, benchmark results, sprint closes. Minor (x.Y.0) and major (X.0.0) are reserved for the PI to call, and mark corpus-significant events: a retraction that moves published numbers, a change to the QA gate or the agent protocol, an architectural change, a reorganization of the paper series. If a sprint feels like it warrants more than a patch, say so in the session summary and let the PI decide — do not bump it unilaterally."""

RULE_NEW = RULE_OLD + """

**Standing exception (added 2026-09-14, PI direction): any change to the `/qa` skill is a MINOR bump, automatically.** This one does not wait for a per-event PI call. `.claude/commands/qa.md` defines what the certification gate does, so a change there changes the meaning of every subsequent PASS — a reader seeing the second number move should be able to infer that the instrument itself moved, without reading the entry. Adding, removing or rescoping a deterministic check (C-criterion) counts; fixing a typo in the file does not."""

ARCHIVE_SECTION = """
### The Lorentzian tail is archived (same release, PI direction)

Papers 46, 47, 48 and 49 moved to `papers/archive/` with `git mv`, so history is preserved. This is the first use of the policy above, and it went in the order the policy prescribes: measure, decide, move, record, stamp.

**Scope, and the boundary.** The four archived are the supplanted tail: descoped or partial, with dependents that are mostly other descoped papers. Three Lorentzian-arc papers were deliberately **not** archived, and the reasons are recorded so the boundary is not guesswork:

| Paper | Disposition | Why |
|---|---|---|
| 43, 50 | kept | ACTIVE, with seven and three external citers |
| 45 | kept | DESCOPED but load-bearing: five healthy dependents, including Paper 38, the WH1 keystone |
| 52, 53 | kept | DRAFT with no dependents, which is *unfinished* rather than supplanted — §9 says zero dependents is not grounds |

**All four are CLOSED-VALID, not CLOSED.** Each retains content that is not refuted, and the register names it so a future sprint reuses rather than rebuilds: Lemma 3.2's degeneracy diagnosis (46), the norm-resolvent arrow and three-carrier identification (47), the bridge's categorical design (48), the cocycle-deficit / TICI algebra (49).

**What moved with them.** Four `INDEX.md` tombstones, four register rows with trigger terms, and archive notes in the group1 synthesis and the field guide — both of which narrate the arc at length and were already honest about the descope, so only the *location* changed and no claim was touched. C24 re-run: 11 archived papers, all registered, register intact.

**Nothing was deleted, and nothing was retracted.** The four remain in the repo, keep their Zenodo DOIs, and keep their own abstracts' descope statements.

"""

applied, failed = [], []

# ---- CLAUDE.md: standing rule + version -----------------------------------
with io.open(CM, encoding="utf-8") as fh:
    c = fh.read()

if "Standing exception (added 2026-09-14, PI direction)" in c:
    applied.append("standing rule already present")
elif RULE_OLD in c:
    c = c.replace(RULE_OLD, RULE_NEW, 1)
    applied.append("CLAUDE.md Sec. 9: qa-skill changes are automatically minor")
else:
    failed.append("version-rule anchor")

if "**Version:** v5.11.21 (September 14, 2026)" in c:
    c = c.replace("**Version:** v5.11.21 (September 14, 2026)",
                  "**Version:** v5.12.0 (September 14, 2026)", 1)
    applied.append("CLAUDE.md version -> v5.12.0")
else:
    failed.append("version anchor")

c = c.replace(
    "- **Paper-retirement policy + C24 (2026-09-14, v5.11.21):** four "
    "supplanted papers sat live because nothing owned the decision; "
    "register + dependency-profile gate added. Two metrics failed first. "
    "See CHANGELOG v5.11.21.",
    "- **Paper-retirement policy + C24; Lorentzian tail archived (2026-09-14, "
    "v5.12.0):** nothing owned the retire decision, so four supplanted papers "
    "stayed live. Policy + gate + register; 46-49 archived. See CHANGELOG "
    "v5.12.0.", 1)
if "v5.12.0):**" in c:
    applied.append("CLAUDE.md Sec. 2 one-liner updated")
else:
    failed.append("Sec 2 bullet")

with io.open(CM, "w", encoding="utf-8") as fh:
    fh.write(c)

# ---- CHANGELOG: renumber and extend ---------------------------------------
with io.open(CL, encoding="utf-8") as fh:
    t = fh.read()

if "## [v5.12.0] - 2026-09-14" in t:
    applied.append("changelog already renumbered")
elif "## [v5.11.21] - 2026-09-14" in t:
    t = t.replace("## [v5.11.21] - 2026-09-14", "## [v5.12.0] - 2026-09-14", 1)
    t = t.replace("See CHANGELOG v5.11.21", "See CHANGELOG v5.12.0")
    applied.append("CHANGELOG v5.11.21 -> v5.12.0")
else:
    failed.append("changelog version heading")

# fold the archive in, and replace the now-obsolete PI note
old_note = """**PI note.** Bumped as a patch, but this is arguably a minor: Sec. 9 gains a standing policy and the `/qa` deterministic list gains a gate, both of which Sec. 9 lists as corpus-significant. Flagging rather than bumping unilaterally, per the version rule. **I also wired C24 into `.claude/commands/qa.md` so it actually runs** -- gate changes have historically been PI-directed, and a gate that never runs is the failure mode this corpus has hit four times, so I chose running over inert. Revert that one line if you would rather it stayed out.

**Nothing was retired.** The Lorentzian decision is yours and remains open."""

new_note = """**Version.** Minor, per PI direction, and the rule is now standing: **any change to the `/qa` skill is automatically a minor bump** (Sec. 9). C24 is wired into `.claude/commands/qa.md` so it actually runs -- a gate that never runs is the failure mode this corpus has hit four times.
""" + ARCHIVE_SECTION

if old_note in t:
    t = t.replace(old_note, new_note, 1)
    applied.append("CHANGELOG: PI note replaced; archive section folded in")
else:
    failed.append("changelog PI note")

with io.open(CL, "w", encoding="utf-8") as fh:
    fh.write(t)

print("applied %d" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
