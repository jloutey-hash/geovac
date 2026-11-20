"""Repo-hygiene sweep (2026-11-20): fourth run of the 2026-06-10 sweep.

Same rules as archive_sweep_2026_11_17.py (imported, not copied): keep what is pinned, referenced, recent, or
imported by something kept; move the rest to debug/archive/<arc>/. The cutoff and the manifest name differ, and
one pass is added: a file named as `debug/<name>` anywhere in CLAUDE.md, README.md, docs/ (outside docs/archive),
tests/, geovac/, papers/, .claude/ or debug/qa/ is kept. The earlier reference set missed docs/predictions/ and
docs/claim_test_matrix.md, whose records name their drivers (76 files on this run).

Nothing is deleted. Manifest (old -> new): debug/archive/sweep_manifest_2026_11_20.json.
Run with --dry-run to print counts without moving anything.
"""
from __future__ import annotations

import json
import re
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import archive_sweep_2026_11_17 as third  # noqa: E402

first = third.first
second = third.second
ROOT, DEBUG, ARCHIVE = third.ROOT, third.DEBUG, third.ARCHIVE
CUTOFF = datetime(2026, 11, 18)
PINNED = third.PINNED | {"archive_sweep_2026_11_20.py"}
TAG = "2026_11_20"


def path_referenced() -> set[str]:
    out: set[str] = set()
    pat = re.compile(r"debug/([A-Za-z0-9_.\-]+)")
    roots = [ROOT / "CLAUDE.md", ROOT / "README.md", ROOT / "docs", ROOT / "tests", ROOT / "geovac",
             ROOT / "papers", ROOT / ".claude", DEBUG / "qa"]
    for r in roots:
        paths = [r] if r.is_file() else [p for p in r.rglob("*") if p.is_file()
                                         and p.suffix in {".md", ".py", ".tex", ".json", ".txt"}]
        for p in paths:
            if "docs/archive" in p.as_posix() or "_archive" in p.as_posix():
                continue
            try:
                out.update(pat.findall(p.read_text(encoding="utf-8", errors="ignore")))
            except OSError:
                pass
    return out


def main(dry_run: bool) -> int:
    refs = second.referenced_basenames() | path_referenced()
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
            for mod in third.imported_names(files[name]):
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
            dest = dest_dir / (f.stem + "__dup4" + f.suffix)
        assert not dest.exists(), dest
        if not dry_run:
            dest_dir.mkdir(parents=True, exist_ok=True)
            f.rename(dest)
        manifest[f"debug/{name}"] = str(dest.relative_to(ROOT)).replace("\\", "/")
        moved[arc] = moved.get(arc, 0) + 1
    if not dry_run:
        (ARCHIVE / f"sweep_manifest_{TAG}.json").write_text(
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
