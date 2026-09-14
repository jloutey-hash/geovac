"""OWED ITEM 3/4 -- the summary-surface reading rule, into CLAUDE.md Sec. 9.

PI-directed 2026-09-12 after the /qa paper_60 FULL run.  Sec. 9 is PM-editable
(the Sec. 13.5 access table: "7-9 Code/Coding/Workflow | Yes").

WHY A READING RULE RATHER THAN ANOTHER PATTERN.  The corpus already has three
mechanisms for stale claims -- the C16 phrase registry, `cited_by` dependents,
and the `rests on:` edges -- and all three are DOCUMENT-granular.  Every defect
in the FULL run was LOCUS-granular: the abstract contradicting the body, the
conclusion contradicting the section, inside ONE file.  The dependents rule had
nothing to say because there was no second document.  And a phrase registry
cannot catch a paraphrase, which is what a summary always is.

The measured token argument, which is why this is affordable:

    abstract + conclusion of Paper 60   ~6,200 tokens
    whole paper                        ~29,500 tokens
    the /qa FULL run that found these   ~1,500,000 tokens

Reading the paper costs ~2% of reviewing it.  We were spending fifty times more
to FIND these than it costs to PREVENT them.

Idempotent.
"""
from __future__ import annotations

import sys

C = "CLAUDE.md"
MARKER = "### Summary-Surface Reading Rule"
ANCHOR = "### Benchmarking Rule\n"

NEW = """### Summary-Surface Reading Rule (added 2026-09-12, PI direction)

**When a claim changes, reread the places that summarize it — in the same edit.**
Not a grep. A read.

1. **Every claim change:** reread the paper's **abstract**, **conclusion**, any
   paragraph labelled **Scope**, and the **Acknowledgments**. For a paper the
   size of Paper 60 that is ~6k tokens.
2. **Every session that touches claims:** read the **whole paper** once. ~30k
   tokens for Paper 60.
3. **A paper's synthesis moves with the paper.** If the claim is summarized in
   `papers/synthesis/`, that file is part of the change, not a separate errand.

*Why a reading rule and not another pattern.* The corpus has three mechanisms
for stale claims — the C16 phrase registry, `cited_by` dependents, and the
§13.8 `rests on:` edges — and **all three are document-granular**. The
2026-09-12 `/qa paper_60` FULL run found every defect **inside one file**: the
abstract claimed a 1953 theorem as derived here while the body carried
`[PRIOR ART]`; the conclusion confined a lever to the case the body had
breached; a retired energy survived as the abstract's validation endpoint under
a paragraph certifying the section converged. `cited_by` correctly reported no
dependent documents, because there were none. And no phrase registry can catch a
paraphrase — which is what a summary is, by construction.

*Why it is affordable.* Measured on Paper 60: abstract + conclusion ≈ 6.2k
tokens, whole paper ≈ 29.5k, against ≈1.5M for the review pass that found these
defects. **Reading the paper costs about 2% of reviewing it.** Cost is not a
reason to skip it.

*What it does not cover, stated so nobody over-trusts it.* Roughly half the
findings of that run were in summary surfaces; the rest were in body text, and
most of those were sentences contradicting themselves rather than stale copies.
Self-contradiction is a different failure and is caught by the reviewers'
internal-consistency mandate, not by this rule.

*The phrase registries stay.* They cost nothing to run and occasionally catch
something. They are a backstop, not the mechanism.

"""


def main() -> int:
    with open(C, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(ANCHOR) != 1:
        print(f"anchor count={t.count(ANCHOR)}; ABORT")
        return 2
    with open(C, "w", encoding="utf-8") as fh:
        fh.write(t.replace(ANCHOR, NEW + ANCHOR))
    print("applied: Summary-Surface Reading Rule added to CLAUDE.md Sec. 9")
    return 0


if __name__ == "__main__":
    sys.exit(main())
