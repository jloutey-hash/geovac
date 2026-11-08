"""Second pass of the debug/-refs sweep: family citations.

Some papers cite debug/ artifacts as families — brace or glob notation
like  debug/mr\\_a\\_dirac\\_propinquity\\_*  or
debug/w3\\_pmns\\_recheck.\\{py,\\_memo.md\\}. The first-pass parser
captures only the prefix, and no single file matches it. This pass takes
those prefixes, finds every path ever added under debug/ in HEAD history
that starts with one, and resurrects each member (same mechanism as
sweep_debug_refs.py). Over-restoring a sibling file is harmless in the
frozen repo; a citation whose family has zero members is reported.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Set

ROOT = Path(__file__).resolve().parents[1]

FRAGMENTS = [
    "debug/ads_track_a_s5_",
    "debug/compute_q5p_tc1",
    "debug/data/energy_graph_nmax",
    "debug/data/ep2",
    "debug/data/mr_",
    "debug/data/sprint_q5p_tc1",
    "debug/l2_",
    "debug/ls",
    "debug/mr_",
    "debug/mr_a_dirac_propinquity_",
    "debug/mr_b_spectral_action_rate_",
    "debug/mr_c_l2_subleading_",
    "debug/sprint_q5p_tc1",
    "debug/w3_lepton_mass_recheck",
    "debug/w3_mechanical_basis",
    "debug/w3_pmns_recheck",
    "debug/z_cliff_probe_",
]


def all_history_paths() -> List[str]:
    out = subprocess.run(
        ["git", "log", "--pretty=format:", "--name-only",
         "--diff-filter=A", "--", "debug/"],
        cwd=ROOT, capture_output=True, text=True, encoding="utf-8",
        errors="replace")
    return sorted({p for p in out.stdout.splitlines()
                   if p.startswith("debug/")})


def restore(path: str) -> bool:
    commit = subprocess.run(
        ["git", "rev-list", "-n", "1", "HEAD", "--", path],
        cwd=ROOT, capture_output=True, text=True).stdout.strip()
    if not commit:
        return False
    for ref in (f"{commit}^:{path}", f"{commit}:{path}"):
        out = subprocess.run(["git", "show", ref], cwd=ROOT,
                             capture_output=True)
        if out.returncode == 0:
            target = ROOT / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(out.stdout)
            return True
    return False


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    history = all_history_paths()
    print(f"{len(history)} debug/ paths ever added in HEAD history")
    empty: List[str] = []
    restored = 0
    for frag in FRAGMENTS:
        members: Set[str] = {p for p in history if p.startswith(frag)}
        live = [m for m in members if (ROOT / m).exists()]
        todo = [m for m in members if not (ROOT / m).exists()]
        done = [m for m in todo if restore(m)]
        restored += len(done)
        status = f"{len(live)} live, {len(done)}/{len(todo)} restored"
        if not members:
            empty.append(frag)
            status = "NO MEMBERS IN HISTORY"
        print(f"  {frag:45s} {status}")
    print(f"\nfamily members restored: {restored}")
    if empty:
        print("citations with zero family members (repoint by hand):")
        for f in empty:
            print(f"  {f}")


if __name__ == "__main__":
    main()
