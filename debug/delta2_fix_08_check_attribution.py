"""DELTA #2 / N8 -- `generate_table --check` could not see the field the last
defect lived in.

DELTA #1's LARGE was a prior-art credit present in the generated markdown and
absent from the generator, so the distributed JSON lacked it and the next
regeneration would have erased the human-readable copy.  That was fixed.  But
the artifact's own drift check compares ONLY the `value` field:

    old = {e["id"]: e["value"] for e in stored["entries"]}

so deleting the Kac-Murdock-Szego credit from `method` (or from the new
`provenance` field) leaves `--check` reporting green.  The guard added after the
defect does not cover the defect.

Widened to compare `method` and `provenance` as well, reported as ATTRIBUTION
DRIFT and failing the check.

Written via the Write tool rather than a heredoc: the replacement text contains
f-string newline escapes, and bash heredocs halve backslashes -- the standing
project rule, violated once already while making this very edit.

Idempotent.
"""
from __future__ import annotations

import sys

G = "benchmarks/certified_reference/generate_table.py"

OLD = '''        bad = [k for k in new if k in old and old[k] != new[k]]
        missing = sorted(set(old) - set(new))
        print(f"{len(new)} entries; {len(bad)} value mismatches; "
              f"{len(missing)} missing")
        for k in bad:
            print(f"  MISMATCH {k}\\n    stored {old[k]}\\n    fresh  {new[k]}")
        return 1 if (bad or missing) else 0
'''

NEW = '''        bad = [k for k in new if k in old and old[k] != new[k]]
        missing = sorted(set(old) - set(new))
        # Added 2026-09-13 (/qa paper_60 DELTA #2).  This check compared ONLY
        # `value`, so the prior-art credit DELTA #1 found missing from the
        # generator could be deleted again from `method`/`provenance` with the
        # gate still green.  A check that cannot see the field the last defect
        # lived in is not guarding that defect.
        attr_fields = ("method", "provenance")
        old_attr = {(e["id"], f): e.get(f, "") for e in stored["entries"]
                    for f in attr_fields}
        new_attr = {(r["id"], f): r.get(f, "") for r in rows
                    for f in attr_fields}
        drift = [k for k in new_attr
                 if k in old_attr and old_attr[k] != new_attr[k]]
        print(f"{len(new)} entries; {len(bad)} value mismatches; "
              f"{len(missing)} missing; {len(drift)} attribution drift(s)")
        for k in bad:
            print(f"  MISMATCH {k}\\n    stored {old[k]}\\n    fresh  {new[k]}")
        for eid, f in drift:
            print(f"  ATTRIBUTION DRIFT {eid}.{f}\\n"
                  f"    stored {old_attr[(eid, f)][:160]}\\n"
                  f"    fresh  {new_attr[(eid, f)][:160]}")
        return 1 if (bad or missing or drift) else 0
'''


def main() -> int:
    with open(G, encoding="utf-8") as fh:
        t = fh.read()
    if "ATTRIBUTION DRIFT" in t:
        print("already applied")
        return 0
    n = t.count(OLD)
    if n != 1:
        print(f"  MISS: anchor count={n}")
        return 3
    with open(G, "w", encoding="utf-8") as fh:
        fh.write(t.replace(OLD, NEW))
    print("  ok    --check now sees method + provenance")
    return 0


if __name__ == "__main__":
    sys.exit(main())
