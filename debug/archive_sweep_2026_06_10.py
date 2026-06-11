"""Repo-hygiene sprint (2026-06-10): sweep debug/ top-level files into debug/archive/<arc>/.

Rules:
  KEEP at top level if any of:
    - modified on/after 2026-06-01 (current frontier era, v3.92+)
    - referenced by any file in tests/ (frozen falsifiers + docstring provenance)
    - explicitly pinned (followon_register.md, compile_3pass.sh, this script, README.md)
  Otherwise MOVE to debug/archive/<arc>/ chosen by ordered filename-prefix rules.

Nothing is deleted. debug/data/, debug/plots/, debug/track_logs/ and all other
subdirectories are untouched. A manifest (old -> new) is written to
debug/archive/sweep_manifest_2026_06_10.json so every historical pointer in
CHANGELOG.md / CLAUDE.md remains mechanically resolvable.
"""
from __future__ import annotations

import json
import re
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEBUG = ROOT / "debug"
ARCHIVE = DEBUG / "archive"
TESTS = ROOT / "tests"
CUTOFF = datetime(2026, 6, 1)

PINNED = {
    "followon_register.md",
    "compile_3pass.sh",
    "archive_sweep_2026_06_10.py",
    "sweep_manifest_2026_06_10.json",
    "README.md",
}

# Ordered (first match wins), applied to lowercase filename.
ARC_RULES: list[tuple[str, str]] = [
    ("q5p_tannakian", r"q5p|tannakian|pro_system|prosystem|jlo|kleinschmidt|hain|^hb_|cosmic"),
    ("spectral_triple_arc", r"^r25|^r3[._-]|wh1|berezin|connes|real_structure|almost_commut"
                            r"|operator_system|circulant|gh_conv|fejer|^ts_|track_ts|supertrace"
                            r"|^st_|^na1|chirality|spinor|spectral_triple|propinquity"),
    ("lorentzian_arc", r"lorentz|krein|^l3[abc]|^l[012][_-]|wick|modular|tomita|^bw_|oslpls"
                       r"|mondino|bridge|four_witness|unruh|wedge|carrier|^p4[5-9]|paper4[5-9]"
                       r"|hypertopology|twin_paradox|tici"),
    ("gravity_arc", r"^g[1-8][_-]|^g4_|gravity|graviton|cigar|warp|^bh_|moebius|mobius"
                    r"|sommerfeld|cheeger|replica|newton|fierz|lichnerowicz|^b[1-4]_"
                    r"|disk|cone|seeley|black_hole|^l6_|^gd[0-9]"),
    ("rh_arc", r"^rh_|ihara|ramanujan|alon|smin|galois|riemann|functional_equation"
               r"|spectral_zero|spectral_chi|hp_operator|hashimoto|gue|zeta_zero"),
    ("nuclear_arc", r"nuclear|deuteron|minnesota|he4|moshinsky|^track_n[a-z]|shell_model|^n[efhijk]_"),
    ("precision_arc", r"autopsy|hfs|21cm|^muh|mu_h|zemach|roothaan|hylleraas|polarizab"
                      r"|precision|cesium|^cs_|rydberg|isotope|drake|eides|^rz_|r_z|lamb"
                      r"|fine_structure|hyperfine|^ps_|positronium"),
    ("qed_arc", r"qed|^ls[0-9_-]|ls8|uehling|bethe|^gn[_-]|^vq[_-]|vector|furry|vertex"
                r"|self_energy|vacuum_polar|schwinger|two_loop|three_loop|bound_state"
                r"|photon|wilson|gauge"),
    ("alpha_arc", r"^kcc|alpha|phase4|sm_running|^bf_d|b_f_delta"),
    ("chemistry_qc_arc", r"^w1|nah|lih|beh|mgh|h2o|h2s|hcl|sih|geh|ash|hbr|^kh_|cah|srh|bah"
                         r"|chem|composed|balanced|^pk_|pseudopot|frozen|sturmian|cusp|fci"
                         r"|casimir|^tc_|jastrow|mvs|tapering|dmrg|vqe|uccsd|scf|^hf_|sto"
                         r"|gaussian|pauli|qubit|ecosystem|trotter|molecul|diatomic|hydride"
                         r"|slater|zeff|z_eff|screen|core|valence|orbital|spec|hamiltonian"
                         r"|level[0-9]|hypersph|prolate|laguerre|gegenbauer|he_|h2_|track_[a-d][a-z]"),
    ("audits", r"audit|cleanup|citation|confidence|literature|verify|hallucin|erratum"),
    ("tracks_misc", r"^track_"),
]
DEFAULT_ARC = "misc"

REF_RE = re.compile(r"debug[/\\]([\w\.\-]+)")


def test_referenced_basenames() -> set[str]:
    refs: set[str] = set()
    for p in TESTS.rglob("*.py"):
        try:
            text = p.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for m in REF_RE.finditer(text):
            name = m.group(1).rstrip(".")
            if name:
                refs.add(name)
    return refs


def classify(name: str) -> str:
    low = name.lower()
    for arc, pattern in ARC_RULES:
        if re.search(pattern, low):
            return arc
    return DEFAULT_ARC


def main() -> int:
    keep_refs = test_referenced_basenames()
    manifest: dict[str, str] = {}
    kept_recent, kept_ref, kept_pin = [], [], []
    moved: dict[str, int] = {}

    for f in sorted(DEBUG.iterdir()):
        if not f.is_file():
            continue
        name = f.name
        if name in PINNED:
            kept_pin.append(name)
            continue
        if name in keep_refs or f.stem in {Path(r).stem for r in keep_refs}:
            kept_ref.append(name)
            continue
        if datetime.fromtimestamp(f.stat().st_mtime) >= CUTOFF:
            kept_recent.append(name)
            continue
        arc = classify(name)
        dest_dir = ARCHIVE / arc
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / name
        if dest.exists():
            dest = dest_dir / (f.stem + "__dup" + f.suffix)
        f.rename(dest)
        manifest[f"debug/{name}"] = str(dest.relative_to(ROOT)).replace("\\", "/")
        moved[arc] = moved.get(arc, 0) + 1

    ARCHIVE.mkdir(exist_ok=True)
    (ARCHIVE / "sweep_manifest_2026_06_10.json").write_text(
        json.dumps(manifest, indent=1, sort_keys=True), encoding="utf-8"
    )

    print(f"moved: {sum(moved.values())} files into {len(moved)} arcs")
    for arc in sorted(moved):
        print(f"  {arc:24s} {moved[arc]}")
    print(f"kept (June frontier): {len(kept_recent)}")
    print(f"kept (test-referenced): {len(kept_ref)}")
    print(f"kept (pinned): {len(kept_pin)}")
    remaining = sum(1 for x in DEBUG.iterdir() if x.is_file())
    print(f"debug/ top-level files remaining: {remaining}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
