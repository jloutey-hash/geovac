---
description: DEPRECATED — renamed to /checkpoint. Use /checkpoint (or /checkpoint push).
---

**This command was renamed to `/checkpoint` on 2026-08-26. Run `/checkpoint` instead.**

Do not follow this file's former protocol. Two things were wrong with it:

1. **The name was a misnomer.** It never created a GitHub Release. Releases are made manually
   on GitHub and mint a Zenodo DOI via the webhook — a PI action. This command only ever did
   bump + commit + tag + push.
2. **It bundled push into the same atomic step**, so any sprint the PI didn't want pushed got
   skipped entirely — and lost its tag too. v5.1.1, v5.1.2 and v5.1.3 were all committed with
   no tag for exactly that reason.
3. **Its push step targeted the wrong ref.** It said `git push origin main`, which contradicts
   the standing policy that merge-to-main is PI-only (CLAUDE.md §2), and would have pushed a
   stale local `main` from a working branch.

`/checkpoint` fixes all three: push is opt-in and off by default, it pushes the *current
branch* and never `main`, and it refuses to tag without a successful commit on a clean tree.

If the PI typed `/release`, treat it as `/checkpoint` (no push) and say so.
