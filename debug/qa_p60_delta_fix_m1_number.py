"""CODE-M1: the one genuinely wrong number in Paper 60, fixed and registered.

The paper's s-sector span-deficit comparison reads

    "0.15 mHa ... at K=136 ... against 4.43 mHa locked"

but 4.43 is the K=105 row.  From the driver's own output
(debug/data/p60_freescale_l0.json, the s-only ladder):

    K=105   free 0.22190   locked 4.43413
    K=136   free 0.14662   locked 4.40376

So the free value is quoted at K=136 and the locked value at K=105 -- a
mismatched pair, in the sentence that carries the floor's re-attribution to the
scale lock.  No result moves ("the span is not the limitation; the lock is" is
unaffected either way), but this paper's entire recent history is exactly this
class, and the pair is the evidence for the keystone.

WHY IT SURVIVED, and the durable half of the fix:  4.43 was an UNREGISTERED
literal, so C21 was structurally blind to it.  Worse, the declared-debt table
written into docs/qa/paper_60.done.md this session asserts "Every one is
currently CORRECT" of its seven listed literals -- and this one is not listed and
is not correct.  A debt declaration that enumerates the wrong set launders the
exposure it was written to disclose.

Both values are therefore registered from the measured driver output and the
paper loci annotated, so the pair cannot drift apart again.
"""
import io, json

# ---- verify against the driver output before touching anything (Sec.15 rule 3:
#      never register a value you have not measured).
rows = {r["K"]: r for r in json.load(open("debug/data/p60_freescale_l0.json"))}
assert abs(rows[136]["gap_iso"] - 4.40376) < 1e-4, rows[136]["gap_iso"]
assert abs(rows[136]["gap_var"] - 0.14662) < 1e-4, rows[136]["gap_var"]
assert abs(rows[105]["gap_iso"] - 4.43413) < 1e-4, rows[105]["gap_iso"]
print("driver verified: K=136 locked %.5f free %.5f  |  K=105 locked %.5f"
      % (rows[136]["gap_iso"], rows[136]["gap_var"], rows[105]["gap_iso"]))

# ---- 1. register both halves of the pair.
R = "debug/qa/numeric_registry.py"
r = io.open(R, encoding="utf-8").read()

OLD_PROV = ('                   "counterpart (against the known exact s-limit -2.879029 Ha) is "\n'
            '                   "0.15 mHa at K=136 vs 4.43 locked. Pipeline unit-tested at K=1, "\n')
NEW_PROV = ('                   "counterpart (against the known exact s-limit -2.879029 Ha) is "\n'
            '                   "p60_span_deficit_sonly_free vs p60_span_deficit_sonly_locked, "\n'
            '                   "both at K=136 (this prose said \'4.43 locked\' until 2026-09-11, "\n'
            '                   "which is the K=105 row -- a mismatched pair). Pipeline unit-tested at K=1, "\n')
assert OLD_PROV in r, "provenance locus not found"
r = r.replace(OLD_PROV, NEW_PROV, 1)

ANCHOR = '    "p60_posing_cost_ground": dict('
NEW_KEYS = '''    "p60_span_deficit_sonly_locked": dict(
        value=4.4038, convention="constant: mHa above the exact He s-limit "
                                 "(-2.879029 Ha) reached by the LOCKED metric-free "
                                 "isoenergetic posing, s-only, K=136 -- the partner "
                                 "of p60_span_deficit_sonly_free at the SAME K",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_freescale_resource.py; read back "
                   "from debug/data/p60_freescale_l0.json 2026-09-11. Registered "
                   "because the paper carried 4.43 here -- the K=105 row -- for "
                   "three days as an unregistered literal C21 could not see. The "
                   "two halves of this comparison are a MATCHED PAIR and must move "
                   "together: quoting one at K=136 and the other at K=105 is the "
                   "defect this key exists to prevent.",
        aliases={4.4341: "K=105", 4.4805: "K=78", 4.5565: "K=55"}),
    "p60_span_deficit_sonly_free": dict(
        value=0.1466, convention="constant: mHa above the exact He s-limit "
                                 "(-2.879029 Ha) reached by a VARIATIONAL CI over "
                                 "the identical s-only Goscinskian span with the "
                                 "global scale lambda optimized, K=136",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_freescale_resource.py; read back "
                   "from debug/data/p60_freescale_l0.json 2026-09-11. Partner of "
                   "p60_span_deficit_sonly_locked at the same K.",
        aliases={0.2219: "K=105", 0.3555: "K=78", 0.6135: "K=55"}),
'''
assert ANCHOR in r
r = r.replace(ANCHOR, NEW_KEYS + ANCHOR, 1)
io.open(R, "w", encoding="utf-8").write(r)
print("registered p60_span_deficit_sonly_locked = 4.4038 / _free = 0.1466")

# ---- 2. fix + annotate both paper loci.
P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
t = io.open(P, encoding="utf-8").read()

pairs = [
 ("same comparison is $0.15$~mHa against $4.43$~mHa at $K=136$.",
  "same comparison is $\\gvq{p60_span_deficit_sonly_free}{0.15}$~mHa against "
  "$\\gvq{p60_span_deficit_sonly_locked}{4.40}$~mHa at $K=136$."),
 ("($s$-only, where $-2.879029$~Ha is known independently) against $4.43$~mHa\nlocked",
  "($s$-only, where $-2.879029$~Ha is known independently) against "
  "$\\gvq{p60_span_deficit_sonly_locked}{4.40}$~mHa\nlocked"),
]
for old, new in pairs:
    assert old in t, "paper locus not found: %.60s" % old
    t = t.replace(old, new, 1)
io.open(P, "w", encoding="utf-8").write(t)
print("paper: 4.43 -> 4.40 at 2 loci, both annotated to the registry")

# ---- 3. correct the declared-debt table's false universal.
D = "docs/qa/paper_60.done.md"
d = io.open(D, encoding="utf-8").read()
OLD_DEBT = "**Every one is currently CORRECT** — none is a retired value."
NEW_DEBT = ("**Every one listed above was verified correct on 2026-09-11** — none is a "
            "retired value.  *But the table was not exhaustive, and that is the "
            "lesson:* the DELTA run found an eighth un-delegated literal, the s-only "
            "locked span deficit, written as `4.43` when the measured K=136 value is "
            "`4.40` (the 4.43 is the K=105 row).  It was not in this table, so the "
            "table's own reassurance did not cover it.  An enumeration offered as "
            "complete is a stronger claim than the literals it lists;  it is now "
            "registered (`p60_span_deficit_sonly_locked` / `_free`) and this table "
            "asserts only what it enumerates.")
assert OLD_DEBT in d, "debt-table locus not found"
d = d.replace(OLD_DEBT, NEW_DEBT, 1)
io.open(D, "w", encoding="utf-8").write(d)
print("done.md: declared-debt table no longer claims a completeness it did not have")
