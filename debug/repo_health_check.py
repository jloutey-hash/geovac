# -*- coding: utf-8 -*-
"""Repo health gate (run as /release precondition 8; added 2026-06-10, v3.110.1).

Warns when the always-loaded / always-paid surfaces regrow past budget:
  CLAUDE.md          > 150 KB   (loaded into every PM session and sub-agent dispatch)
  debug/ top level   > 600 files (re-run debug/archive_sweep with a later CUTOFF)
  MEMORY.md          > 24 KB    (harness silently truncates ~24.4 KB)

Exit 0 = all green; exit 1 = at least one budget exceeded (warn, don't block).
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
MEMORY = Path.home() / ".claude" / "projects" / "C--Users-jlout-Desktop-Project-Geometric" / "memory" / "MEMORY.md"

checks = []

cm_kb = (ROOT / "CLAUDE.md").stat().st_size / 1024
checks.append(("CLAUDE.md size", f"{cm_kb:.0f} KB", cm_kb <= 150, "<= 150 KB"))

n_debug = sum(1 for p in (ROOT / "debug").iterdir() if p.is_file())
checks.append(("debug/ top-level files", str(n_debug), n_debug <= 600,
               "<= 600 (re-run archive sweep with later CUTOFF)"))

if MEMORY.exists():
    mem_kb = MEMORY.stat().st_size / 1024
    checks.append(("MEMORY.md size", f"{mem_kb:.1f} KB", mem_kb <= 24, "<= 24 KB (truncates ~24.4)"))

ok = True
for name, val, passed, budget in checks:
    flag = "OK  " if passed else "WARN"
    print(f"[{flag}] {name}: {val}  (budget {budget})")
    ok = ok and passed

sys.exit(0 if ok else 1)
