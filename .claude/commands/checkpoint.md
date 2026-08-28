---
description: Version bump + commit + tag (push is opt-in, never automatic). NOT a GitHub Release.
---

Cut a checkpoint: bump the version, commit, tag. **Stop there unless told to push.**

> **This does NOT create a GitHub Release.** Releases are created manually on GitHub, and
> publishing one mints a Zenodo DOI via the webhook. That is a PI action, always. This command
> only produces a local commit + tag (and optionally pushes the working branch).
>
> *(Renamed from `/release` on 2026-08-26. The old name was a misnomer — it never made a
> Release — and it bundled push into the same atomic step, so any sprint the PI didn't want
> pushed got no tag either. v5.1.1, v5.1.2 and v5.1.3 were all committed with no tag for
> exactly that reason.)*

**Usage.**
- `/checkpoint` — bump, commit, tag, stop. **Default. Does not touch the remote.**
- `/checkpoint push` — the same, then push the **current branch** and the tag.

---

**Preconditions (check before any git action; stop and ask if one is unmet).**

1. `git status` is meaningful — what's staged, modified, untracked. Quote it back to the PI.
   **Call out untracked files explicitly.** A `git commit -a` habit stages modifications but
   silently skips new files; that is how 585 files — including the CERTIFIED Paper 60 — sat
   outside the repo for weeks (2026-08-26).
2. CLAUDE.md §1 version string bumped to the new version.
3. CLAUDE.md §2 has a one-liner for the change being cut (per `/sprint-close`).
4. CHANGELOG.md has an entry under the new version heading.
5. Papers edited → they compile three-pass clean. Confirm.
6. Production code edited → relevant tests pass (`/regression touched`). Confirm.
7. Hard-prohibition check (§13.5): nothing in the staged diff violates them.
8. Repo health gate: `python debug/repo_health_check.py`. On WARN (CLAUDE.md > 150 KB,
   debug/ top-level > 600 files, MEMORY.md > 24 KB) report it alongside the checkpoint —
   it doesn't block, but it must be surfaced so bloat never silently regrows.

**Version-bump policy (revised 2026-08-22, PI direction).**
Bump the **last number only** (x.y.Z → x.y.Z+1) by default. Minor (x.Y.0) and major (X.0.0)
are **PI calls, made explicitly** — they mark corpus-significant events (a retraction that
moves published numbers, a change to the QA gate or agent protocol, an architectural change,
a reorganization of the paper series). If a sprint feels bigger than a patch, say so and let
the PI decide; do not bump it unilaterally. A diagnostic arc testing 10 hypotheses and finding
9 negatives is **one** entry, not 10.

---

**Steps.**

1. **Stage.** Prefer explicit `git add` of named files. If you use a directory-level add, note
   that it stages **modifications as well as new files**, and say so in the commit message —
   do not let the message understate what is in the commit.
2. **Commit.** HEREDOC message, title `vX.Y.Z: short description`, body with Added / Changed /
   Closed sections, trailing:
   ```
   Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
   Claude-Session: <session URL>
   ```
   (Keep the model name in sync with the active model.)
3. **Verify the commit succeeded, then tag.** `git tag -a vX.Y.Z -m "vX.Y.Z: short description"`.
   - **Never tag a dirty tree.** `git status --porcelain` must be empty first.
   - **Never tag without a commit in this invocation.** The tag applies to the HEAD the commit
     just produced. This ordering is what guarantees "a commit every time we tag."
4. **Push — ONLY if the invocation said `push`.**
   - `git push origin "$(git rev-parse --abbrev-ref HEAD)"` — the **current branch**.
   - `git push origin vX.Y.Z`.
   - **NEVER `git push origin main`.** Work happens on a working branch; merge-to-main is
     PI-only (see CLAUDE.md §2).
5. **Report** the commit SHA, the tag, and whether anything was pushed.

---

**Hard prohibitions.**
- NEVER push to `main`. NEVER `git push --force`, anywhere.
- NEVER skip hooks (`--no-verify`) or bypass signing.
- NEVER `git reset --hard` without explicit PI direction.
- NEVER stage `.env`, `credentials.json`, or anything that looks like a secret — flag and stop.
- NEVER create a GitHub Release, or suggest doing so as a follow-up step. PI only.
