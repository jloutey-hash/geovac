# Paper 59 integrated-class support (durable test backing)

Backing modules for `tests/test_paper59_integrated_class.py`, the guard on Paper 59
sec:scope's integrated-class claim (record `docs/predictions/2026-11-01_integrated_class.md`, IC).

The instrument that scored IC lives in `debug/` (transient by design, CLAUDE.md section 9), so tests
must not import from it. These are verbatim copies, made 2026-11-02:

| File | Source | sha256 |
|:-----|:-------|:-------|
| `intclass_lib.py` | `debug/intclass_lib.py` (frozen with record IC, hash in commit 40093e7c) | `b1655d6450fec0c1ed5f51d29e89e4d0cdb222da7eadf3d570286e1bff7467a6` |
| `intclass_score.py` | `debug/intclass_score.py` (frozen with record IC, hash in commit 40093e7c) | `aaf700e482f53c0a204939c08e6e12c9786d544e48c77e176503f21e5add91ca` |
| `na_modguess.py` | `debug/na_modguess.py` (record NA; unchanged since AO3) | `84691c87a4a280b6823438af574c0fff6af9bf3f7ba7b76d93f2bbe866847d7d` |
| `na_modanalyze.py` | `debug/na_modanalyze.py` (record NA) | `de92dc3acdcc7082001e32b50ef44db8c3bd75d990f615603b5f97f6ea178005` |

The guard checks these hashes after normalizing CRLF to LF (git may rewrite line endings on checkout; `intclass_lib.py`'s LF-normalized sha256 is `2a142e5b5667b432c7de4be885f5b6767daf3d209c86261076dc6ae5e4f641d4`), so an edit here fails it. Do not edit; if the instrument changes, that
is a new record.
