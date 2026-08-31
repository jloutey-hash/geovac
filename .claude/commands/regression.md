---
description: Run a regression-test slice of tests/ — derived from the current diff by default, with explicit fast/full/topo scopes available.
---

Run pytest on a slice of `tests/`. The slice depends on the `scope` argument.

**Scopes.**

- `touched` (default) — derive the test selection from the current git diff.
  - Run `git diff --name-only HEAD -- 'geovac/**/*.py'` to enumerate touched production modules.
  - For each touched module `geovac/X.py`, find consumer test files via `grep -l "from geovac\.X\b\|import geovac\.X\b" tests/ -r`.
  - Union those test files with the topological-integrity baseline (`tests/test_fock_projection.py`, `tests/test_fock_laplacian.py`).
  - Plus a small random tail-risk sample from `tests/_durations.json` (5 fast tests, picked uniformly) — catches orthogonal regressions that the diff-derived selection would miss by construction.
  - If no production files are touched, run the topological baseline only.
  - Target wall: 30s–2min.

- `fast` — read `tests/_durations.json` and select every test whose last-recorded duration is `< 0.5s`. Run those. For "I just want a smoke check across the suite." Target wall: <1 min.
  **CURRENTLY UNAVAILABLE (2026-08-31): `tests/_durations.json` is an empty
  2-byte file**, dated 2026-06-07. It selects nothing, and the `touched`
  scope's random tail-risk sample is dead for the same reason. Nothing failed
  when it emptied because its only consumer is *this prompt* rather than code
  — the same shape as any derived artifact with no consumer, no refresh
  trigger, and no failure mode when stale. Do **not** silently fall back;
  surface it (see the bootstrap block below).

- `full` — `pytest tests/ --ignore=tests/_archive --tb=line -n auto`. Everything except the archive.
  **MEASURED 2026-08-31: ~6.7 h serial, ~2.4 h with `-n auto`.** This scope is
  **NOT a sprint-close gate** — it is a *scheduled baseline*. Use `touched` at
  sprint close and run `full` on its own schedule.
  *(The previous "target wall: 10–15 min" was wrong by roughly two orders of
  magnitude and is the reason this scope kept being invoked at close, where it
  was then killed by a timeout or — worse — reported green off a piped exit
  code. See `debug/qa/test_suite_cost_memo.md`.)*
  **Always pass `-n auto`**: 2.77× on the heavy slice, 1.94× on a light one,
  with pass/skip/xfail counts **bit-identical** to serial in every arm, so the
  parallelism is free of ordering or shared-state effects. Thread-pinning
  (`OMP_NUM_THREADS=1` etc.) was tested and does **nothing** (4%, noise) — the
  ceiling is not BLAS oversubscription. Untested lead: `--dist loadfile`, since
  the default `--dist load` re-runs module/session fixtures once *per worker*
  and this suite's fixtures build Hamiltonians (one measured at 71 s).

- `topo` — just the 18 symbolic S³ proofs (`tests/test_fock_projection.py tests/test_fock_laplacian.py`). 5–10s. Use when you've only touched papers.

**Default invocation.** `/regression` with no argument is equivalent to `/regression touched`.

**Bootstrap requirement (one-time).** The `fast` and `touched` scopes both read `tests/_durations.json` (the latter to pick the random tail-risk sample). If that file does not exist, generate it once via:

```
pytest tests/ --ignore=tests/_archive --collect-only -q --no-header 2>/dev/null  # warm collection
pytest tests/ --ignore=tests/_archive --durations=0 --tb=no -q --no-header > tests/_durations_raw.log 2>&1
python -c "
import re, json
durations = []
for line in open('tests/_durations_raw.log'):
    m = re.match(r'^(\\d+\\.\\d+)s call\\s+(.+)$', line.strip())
    if m:
        durations.append({'duration_s': float(m.group(1)), 'nodeid': m.group(2)})
json.dump(sorted(durations, key=lambda d: d['duration_s']), open('tests/_durations.json', 'w'), indent=2)
print(f'wrote tests/_durations.json with {len(durations)} entries')
"
```

If the file is missing, do NOT silently fall back to a broader scope — surface the missing-baseline state to the PI and ask whether to bootstrap now (15–20 min wall) or run a fallback (`topo` baseline + the diff-derived consumer set, skipping the random sample). The fallback is acceptable; running `full` instead is not, because it defeats the speed promise of `/regression`.

**Refresh cadence.** Re-bootstrap `tests/_durations.json` quarterly, or whenever the suite materially grows. Tests whose runtime drifts above 0.5s naturally fall out of the `fast` set on the next baseline.

**Implementation notes.**

- Use `python -m pytest <files> --no-header -q --tb=line -p no:cacheprovider` for all scopes.
- For `touched`, deduplicate the union of consumer files before passing to pytest.
- For the random tail-risk sample in `touched` mode: use a deterministic seed derived from the current git SHA (`hash(commit_sha) % 2**32`) so re-running the same commit reproduces the same sample. This makes the result auditable; "the random sample missed it" is a reproducible state, not a one-shot.
- Report at the end: which scope ran, which test files were selected, pass/fail count, wall time. If `touched` selected zero consumer files (only topo + sample), say so explicitly.

**What this does NOT do.**

- Does not gate any other command. `/checkpoint` continues to ship even if `/regression` would have shown failures — by PI policy, paper-progress sprints should not be blocked by chemistry-test rot. `/regression` is a cheap voluntary check, not enforcement.
- Does not auto-archive failing tests. If `/regression` surfaces a regression, the PM triages per CLAUDE.md §14 (fix-in-place / redirect / archive). The skill reports the failure; the response is human-decided.
- Does not refresh the durations baseline automatically. Stale-baseline drift is a real risk; the quarterly refresh cadence is the protection.

**When to use.**

- After any non-trivial code edit, before declaring sprint-close ready: `/regression` (touched).
- When you've done broad refactors that span more than 2–3 modules: `/regression full`.
- When you only touched papers or memos: `/regression topo`.
- When you want a quick "is anything obviously broken" smoke: `/regression fast`.
