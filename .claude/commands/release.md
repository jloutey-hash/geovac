---
description: Version bump + commit + tag + push (mechanical release wrapper)
---

Run the release protocol. Stop and ask if any precondition is unmet.

**Preconditions (check before any git action).**

1. `git status` is meaningful — what's staged, what's modified, what's untracked. Quote it back to the PI.
2. CLAUDE.md §1 version string has been bumped to the new version.
3. CLAUDE.md §2 has at least one one-liner entry for the change being released (per `/sprint-close`).
4. CHANGELOG.md has an entry under the new version heading.
5. If papers were edited:\ they compile three-pass clean. Confirm.
6. If production code was edited:\ relevant tests pass (`pytest tests/<paths>`). Confirm.
7. Hard-prohibition check (§13.5):\ nothing in the staged diff violates the prohibitions.

**Version-bump policy.**
- Patch (x.y.Z) for documentation / paper edits / dead-end recordings.
- Minor (x.Y.0) for new features, completed diagnostic arcs, new paper additions, operational policy changes.
- Major (X.0.0) for architectural changes.

A diagnostic arc that tests 10 hypotheses and finds 9 negative results is **one** minor version, not 10 patches (per §9 Changelog Protocol).

**Release steps.**

1. Stage the right files. **Prefer explicit `git add` of named files.** Avoid `git add -A` or `git add .` — they pick up secrets, large debug data, untracked stray files.
2. `git commit` with a HEREDOC message following the project commit convention (multi-line, leading title in `vX.Y.Z: short description` form, body with Added / Changed / Closed sections, trailing `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>`).
3. `git tag -a vX.Y.Z -m "vX.Y.Z: short description"`.
4. `git push origin main` then `git push origin vX.Y.Z`.
5. Quote the commit SHA and tag back to the PI as confirmation.

**Hard prohibitions.**
- NEVER `git push --force` (especially to main).
- NEVER skip hooks (`--no-verify`).
- NEVER `git reset --hard` without explicit PI direction.
- NEVER stage `.env`, `credentials.json`, secret files. If something looks sensitive, flag and stop.
