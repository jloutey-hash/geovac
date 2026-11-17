# Paper 59 integrated-class support (durable test backing)

Backing modules for `tests/test_paper59_integrated_class.py`, (added 2026-11-16, record IN) `tests/test_paper59_na_class.py`, (added 2026-11-16, record AK) `tests/test_paper59_aa_bc_class.py` and (added 2026-11-17, record MO) `tests/test_paper59_aa_bc_minimal_operator.py`, the guards on Paper 59
sec:scope's integrated-class claim (record `docs/predictions/2026-11-01_integrated_class.md`, IC).

The instrument that scored IC lives in `debug/` (transient by design, CLAUDE.md section 9), so tests
must not import from it. These are verbatim copies, made 2026-11-02:

| File | Source | sha256 |
|:-----|:-------|:-------|
| `intclass_lib.py` | `debug/intclass_lib.py` (frozen with record IC, hash in commit 40093e7c) | `b1655d6450fec0c1ed5f51d29e89e4d0cdb222da7eadf3d570286e1bff7467a6` |
| `intclass_score.py` | `debug/intclass_score.py` (frozen with record IC, hash in commit 40093e7c) | `aaf700e482f53c0a204939c08e6e12c9786d544e48c77e176503f21e5add91ca` |
| `na_modguess.py` | `debug/na_modguess.py` (record NA; unchanged since AO3) | `84691c87a4a280b6823438af574c0fff6af9bf3f7ba7b76d93f2bbe866847d7d` |
| `na_modanalyze.py` | `debug/na_modanalyze.py` (record NA) | `de92dc3acdcc7082001e32b50ef44db8c3bd75d990f615603b5f97f6ea178005` |
| `na_exact.py` | `debug/na_exact.py` (record NA; frozen with record IN, commit 247836bc) | `85fc39f18bfe9231f47925436b72258dca7afddc6f9619257d9fafc708004afb` (LF-normalized) |
| `intclass_in_check.py` | `debug/intclass_in_check.py` (frozen with record IN, commit 247836bc) | `50fb31871454014a1f2dee653985c778a8e276afbe478a66eecb6e2212fbf4c4` (LF-normalized) |
| `s1_controls.py` | `debug/s1_controls.py` (record MO's pipeline, frozen in commit 56f10c90; copied 2026-11-17 for `tests/test_paper59_aa_bc_minimal_operator.py`) | `1c50016620c5f6291ac735028d9c87a658d8a9c226f28f08a1c667eba9436687` (LF-normalized `a150322d746cb445b3b47e8264cb51a53cfe1dd81b21f8f1ab00eb8f21049453`) |
| `aabc_exact.py` | `debug/aabc_exact.py` (record AB's engine, unchanged since commit 3a2fbfcf; copied 2026-11-16 for record AK's guard `tests/test_paper59_aa_bc_class.py`) | `218b7283f525d25ed8a87df61f0749ca6f576b20abcc4bde1272de32436b2d30` (LF-normalized) |

The guard checks these hashes after normalizing CRLF to LF (git may rewrite line endings on checkout; `intclass_lib.py`'s LF-normalized sha256 is `2a142e5b5667b432c7de4be885f5b6767daf3d209c86261076dc6ae5e4f641d4`), so an edit here fails it. Do not edit; if the instrument changes, that
is a new record.
