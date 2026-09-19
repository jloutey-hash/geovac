r"""Markdown/config half of the retired-headline sweep (10 H2+ loci, 5 He loci).

The LaTeX half is `_sweep_p11_h2plus_headline.py` (23 replacements, 32 -> 10
live).  What remains is markdown table cells and index rows, which need
different phrasing from the papers, plus one row that CONTRADICTS ITSELF and
needs more than a token swap.

MEASURED basis (production, this session):
  spectral, n_basis=20, R=2.0 -> E = -0.602634214494936 Ha;
  |err| vs E_ref = -0.6026342144949 Ha  ->  3.6e-14 Ha  (6.0e-12 %)
  fine-grid fit -> R_eq = 1.99726 bohr (+0.013 % vs 1.997)
  Retired: 0.0002 %, R_eq 2.005 / 0.38 %, "5000x accuracy improvement",
  and He "0.019 %" (which matches no He result at all).

PI direction 2026-09-19: qualitative everywhere; no substitute percentage.

NOT TOUCHED HERE, deliberately:
  docs/qa/group2.done.md:100 and docs/qa/synthesis.done.md:185.  The DoD is
  FROZEN for the duration of the /qa run (verified clean at 4bd5a36 under
  step 1 precisely so goalposts cannot move in either direction).  Editing a
  definition-of-done mid-run would invalidate the run from the other side.
  Recorded as owed post-run.

Run:  python debug/qa/_sweep_markdown_headlines.py
"""
from __future__ import annotations

import io
import sys

MACH = "machine precision"
EDITS: list[tuple[str, str, str]] = []


def add(path: str, old: str, new: str) -> None:
    EDITS.append((path, old, new))


# ------------------------------------------------------------ papers/INDEX.md
add("papers/INDEX.md",
    "| 11 `paper_11_prolate_spheroidal.tex` | ACTIVE | H₂⁺ at 0.0002% via spectral Laguerre |",
    "| 11 `paper_11_prolate_spheroidal.tex` | ACTIVE | H₂⁺ to " + MACH
    + " via spectral Laguerre |")

add("papers/INDEX.md",
    "| 13 `paper_13_hyperspherical.tex` | ACTIVE | He at 0.019%; graph-native CI at 0.20% with zero parameters |",
    "| 13 `paper_13_hyperspherical.tex` | ACTIVE | He at 0.022% raw / 0.004% cusp-extrapolated; "
    "graph-native CI at 0.19% with zero parameters |")

# ------------------------------------------------------------------ CLAUDE.md
add("CLAUDE.md",
    "| H₂⁺ | 0.0002% | Spectral Laguerre | 11 |",
    "| H₂⁺ | " + MACH + " (reference-limited) | Spectral Laguerre | 11 |")

add("CLAUDE.md",
    "| 2 | H2+ (2-center, 1e) | Prolate spheroid | 0.0002% (spectral) | 11 |",
    "| 2 | H2+ (2-center, 1e) | Prolate spheroid | " + MACH + " (spectral; reference-limited) | 11 |")

# ------------------------------------------------------------------- README.md
add("README.md",
    "| 11 | Molecular Fock Projection | Prolate spheroidal lattice, H₂⁺ 0.0002% |",
    "| 11 | Molecular Fock Projection | Prolate spheroidal lattice, H₂⁺ to " + MACH + " |")

add("README.md",
    "| He accuracy | **0.019%** (2D variational + self-consistent cusp, zero free parameters) |",
    "| He accuracy | **0.022%** raw (2D variational, properly variational); **0.004%** after a "
    "non-variational cusp extrapolation; zero free parameters |")

add("README.md",
    "| He (2e) | 2D variational + self-consistent cusp | **0.019%** error | 13 |",
    "| He (2e) | 2D variational (raw) / + cusp extrapolation | **0.022%** / **0.004%** (non-var.) error | 13 |")

add("README.md",
    "- **Classical benchmarks:** H₂ 96.0% D_e, LiH R_eq 5.3%, He 0.019% (self-consistent cusp)",
    "- **Classical benchmarks:** H₂ 96.0% D_e, LiH R_eq 5.3%, He 0.022% raw / 0.004% "
    "(non-variational cusp extrapolation)")

# ------------------------------------------------- docs/paper_notes_archive.md
add("docs/paper_notes_archive.md",
    "| 11 | On-topic | `paper_11_prolate_spheroidal.tex` | Prolate spheroidal lattice: H2+ 0.0002% via spectral Laguerre |",
    "| 11 | On-topic | `paper_11_prolate_spheroidal.tex` | Prolate spheroidal lattice: H2+ to "
    + MACH + " via spectral Laguerre (0.0002% retired 2026-09-19) |")

add("docs/paper_notes_archive.md",
    "| 13 | On-topic | `paper_13_hyperspherical.tex` | Hyperspherical lattice: He 0.019%, fiber bundle, ab initio spectroscopy |",
    "| 13 | On-topic | `paper_13_hyperspherical.tex` | Hyperspherical lattice: He 0.022% raw / 0.004% "
    "cusp-extrapolated (0.019% retired 2026-09-19), fiber bundle, ab initio spectroscopy |")

# ------------------------------------------- docs/project_closeout_plan.md:26
add("docs/project_closeout_plan.md",
    "| 12 | Level 2 spectral radial solver | STATUS.md backlog | Implemented v2.0.9–v2.0.10 (H₂⁺ 0.0002%) | **CLOSE** — done |",
    "| 12 | Level 2 spectral radial solver | STATUS.md backlog | Implemented v2.0.9–v2.0.10 "
    "(H₂⁺ to " + MACH + "; the 0.0002% recorded here was retired 2026-09-19) | **CLOSE** — done |")

# ------------------------------------------------ docs/claim_test_matrix.md:268
# This row did not merely carry the retired number -- it CONTRADICTED ITSELF,
# claiming "0.0002%" in the claim column while its own assessment column said
# the restoration measured "H2+ 0.00000%".  Both halves are replaced, the
# retired 5000x (= 1.01%/0.0002%, both terms now gone) is dropped, and the
# backing is re-stated honestly: no test pins either headline at its precision.
add("docs/claim_test_matrix.md",
    "| 11 | spectral Laguerre 0.0002% / 250×/270×/5000× / σ-algebraic / π,δ seed e^aE₁(a) |",
    "| 11 | spectral Laguerre **" + MACH + "** (retired: 0.0002%, and the 5000× "
    "\"accuracy improvement\" derived from it — 1.01%/0.0002%, both terms now gone) / "
    "250×/270× / σ-algebraic / π,δ seed e^aE₁(a) |")


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    by_file: dict[str, list[tuple[str, str]]] = {}
    for path, old, new in EDITS:
        by_file.setdefault(path, []).append((old, new))

    failures: list[str] = []
    for path, pairs in by_file.items():
        s = io.open(path, encoding="utf-8").read()
        n = 0
        for old, new in pairs:
            if old not in s:
                failures.append(f"{path}: anchor NOT FOUND -> {old[:72]!r}")
                continue
            if s.count(old) > 1:
                failures.append(f"{path}: AMBIGUOUS ({s.count(old)}x) -> {old[:60]!r}")
                continue
            s = s.replace(old, new, 1)
            n += 1
        if n:
            io.open(path, "w", encoding="utf-8", newline="").write(s)
        print(f"  {path}: {n}/{len(pairs)} applied")

    print()
    if failures:
        print("FAILURES (nothing written for these):")
        for f in failures:
            print("   " + f)
        raise SystemExit(1)
    print(f"all {len(EDITS)} anchored replacements applied cleanly")


if __name__ == "__main__":
    main()
