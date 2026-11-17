"""Repo-hygiene sweep (2026-11-17): third run of the 2026-06-10 sweep, cutoff 2026-11-06.

Same rules and arc classification as archive_sweep_2026_06_10.py and the reference set of
archive_sweep_2026_11_08.py (both imported, not copied), with one addition: the import-closure pass that the
2026-11-08 run had to do by hand is automatic here. A top-level .py that is kept (recent, referenced or pinned)
keeps every top-level debug/ module it imports by bare name (`import X`, `from X import ...`), iterated to a
fixed point, so frozen record drivers keep running from debug/.

Nothing is deleted. Manifest (old -> new): debug/archive/sweep_manifest_2026_11_17.json.
Run with --dry-run to print counts without moving anything.
"""
from __future__ import annotations

import json
import re
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import archive_sweep_2026_06_10 as first  # noqa: E402
import archive_sweep_2026_11_08 as second  # noqa: E402

ROOT = first.ROOT
DEBUG = first.DEBUG
ARCHIVE = first.ARCHIVE
CUTOFF = datetime(2026, 11, 6)
PINNED = second.PINNED | {"archive_sweep_2026_11_17.py"}
IMPORT_RE = re.compile(r"^\s*(?:from\s+([A-Za-z_][A-Za-z0-9_]*)\s+import|import\s+([A-Za-z_][A-Za-z0-9_, ]*))", re.M)


def imported_names(path: Path) -> set[str]:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return set()
    out: set[str] = set()
    for m in IMPORT_RE.finditer(text):
        if m.group(1):
            out.add(m.group(1))
        else:
            for part in m.group(2).split(","):
                out.add(part.strip().split(" ")[0])
    return out


def main(dry_run: bool) -> int:
    refs = second.referenced_basenames()
    ref_stems = {Path(r).stem for r in refs}
    files = {f.name: f for f in DEBUG.iterdir() if f.is_file()}
    keep: dict[str, str] = {}
    for name, f in files.items():
        if name in PINNED:
            keep[name] = "pinned"
        elif name in refs or f.stem in ref_stems:
            keep[name] = "referenced"
        elif datetime.fromtimestamp(f.stat().st_mtime) >= CUTOFF:
            keep[name] = "recent"
    changed = True
    while changed:
        changed = False
        for name in list(keep):
            if not name.endswith(".py"):
                continue
            for mod in imported_names(files[name]):
                cand = mod + ".py"
                if cand in files and cand not in keep:
                    keep[cand] = "import-closure"
                    changed = True
    manifest: dict[str, str] = {}
    moved: dict[str, int] = {}
    for name, f in sorted(files.items()):
        if name in keep:
            continue
        arc = first.classify(name)
        dest_dir = ARCHIVE / arc
        dest = dest_dir / name
        if dest.exists():
            dest = dest_dir / (f.stem + "__dup3" + f.suffix)
        if not dry_run:
            dest_dir.mkdir(parents=True, exist_ok=True)
            f.rename(dest)
        manifest[f"debug/{name}"] = str(dest.relative_to(ROOT)).replace("\\", "/")
        moved[arc] = moved.get(arc, 0) + 1
    if not dry_run:
        (ARCHIVE / "sweep_manifest_2026_11_17.json").write_text(
            json.dumps(manifest, indent=1, sort_keys=True), encoding="utf-8")
    print(f"{'would move' if dry_run else 'moved'}: {sum(moved.values())} files into {len(moved)} arcs")
    for arc in sorted(moved):
        print(f"  {arc:24s} {moved[arc]}")
    counts: dict[str, int] = {}
    for why in keep.values():
        counts[why] = counts.get(why, 0) + 1
    print("kept:", counts, "top level after:", len(keep))
    return 0


if __name__ == "__main__":
    sys.exit(main("--dry-run" in sys.argv))
