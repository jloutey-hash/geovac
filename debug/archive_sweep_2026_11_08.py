"""Repo-hygiene sweep (2026-11-08): second run of the 2026-06-10 sweep with a later cutoff.

Same rules and arc classification as archive_sweep_2026_06_10.py (imported, not copied), with two changes:
  - CUTOFF moved to 2026-10-28, so debug/ top level drops under the 600-file health budget;
  - files referenced from geovac/, debug/qa/, docs/qa/ or .claude/ are kept as well as those referenced
    from tests/ (gates and skills read some top-level debug/ files).

Nothing is deleted. A manifest (old -> new) is written to debug/archive/sweep_manifest_2026_11_08.json.
The 2026-11-08 run was followed by an import-closure pass: 43 moved files that kept code imports by bare
module name (sys.path) or opens by bare filename were moved back and dropped from the manifest. The reference
regex above sees only literal 'debug/<name>' strings, so a later run needs that pass too.
Run with --dry-run to print counts without moving anything.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import archive_sweep_2026_06_10 as first  # noqa: E402

ROOT = first.ROOT
DEBUG = first.DEBUG
ARCHIVE = first.ARCHIVE
CUTOFF = datetime(2026, 10, 28)
PINNED = first.PINNED | {"archive_sweep_2026_11_08.py", "repo_health_check.py"}
REF_DIRS = [ROOT / "tests", ROOT / "geovac", DEBUG / "qa", ROOT / "docs" / "qa", ROOT / ".claude"]


def referenced_basenames() -> set[str]:
    refs: set[str] = set()
    for base in REF_DIRS:
        for p in base.rglob("*"):
            if not p.is_file() or p.suffix not in {".py", ".md", ".json", ".sh", ".txt"}:
                continue
            try:
                text = p.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            for m in first.REF_RE.finditer(text):
                name = m.group(1).rstrip(".")
                if name:
                    refs.add(name)
    return refs


def main(dry_run: bool) -> int:
    refs = referenced_basenames()
    ref_stems = {Path(r).stem for r in refs}
    manifest: dict[str, str] = {}
    kept_recent = kept_ref = kept_pin = 0
    moved: dict[str, int] = {}
    for f in sorted(DEBUG.iterdir()):
        if not f.is_file():
            continue
        name = f.name
        if name in PINNED:
            kept_pin += 1
            continue
        if name in refs or f.stem in ref_stems:
            kept_ref += 1
            continue
        if datetime.fromtimestamp(f.stat().st_mtime) >= CUTOFF:
            kept_recent += 1
            continue
        arc = first.classify(name)
        dest_dir = ARCHIVE / arc
        dest = dest_dir / name
        if dest.exists():
            dest = dest_dir / (f.stem + "__dup2" + f.suffix)
        if not dry_run:
            dest_dir.mkdir(parents=True, exist_ok=True)
            f.rename(dest)
        manifest[f"debug/{name}"] = str(dest.relative_to(ROOT)).replace("\\", "/")
        moved[arc] = moved.get(arc, 0) + 1
    if not dry_run:
        (ARCHIVE / "sweep_manifest_2026_11_08.json").write_text(
            json.dumps(manifest, indent=1, sort_keys=True), encoding="utf-8")
    print(f"{'would move' if dry_run else 'moved'}: {sum(moved.values())} files into {len(moved)} arcs")
    for arc in sorted(moved):
        print(f"  {arc:24s} {moved[arc]}")
    print(f"kept (on/after {CUTOFF.date()}): {kept_recent}")
    print(f"kept (referenced): {kept_ref}")
    print(f"kept (pinned): {kept_pin}")
    print(f"top level after: {kept_recent + kept_ref + kept_pin}")
    return 0


if __name__ == "__main__":
    sys.exit(main("--dry-run" in sys.argv))
