"""CLAUDE.md Sec. 5: record that Level 4's placement for H2 rested on the
withdrawn cusp diagnosis, and what survives independently.

Sec. 5 is a Sec. 13.5 hard-prohibition section; this edit is made under
explicit PI direction (2026-09-14, "You can touch section 5 also").  It adds a
note and changes no level, no coordinate system and no table row -- whether H2
should be re-seated is left open as a PI question rather than decided here.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

PATH = "CLAUDE.md"

OLD = """| 4 | H2 (2-center, 2e) | Mol-frame hyperspherical | 96.0% D_e | 15 |
| 4N | LiH (2-center, 4e) | Full mol-frame hypersp. (SO(12)) | R_eq 63.5% (l_max=2, 2D variational; unbound D_e) | 17 |"""

NEW = """| 4 | H2 (2-center, 2e) | Mol-frame hyperspherical | 96.0% D_e | 15 |

*Level 4 note (added 2026-09-14, PI-directed):* **Level 4's placement for H2 was motivated in part by a claim that is now withdrawn, and on accuracy it is no longer the leading geometry for this system.** Paper 12 read its 92.4% prolate-spheroidal residual as the electron-electron cusp being a coordinate singularity in (xi, eta) that demands a coordinate change; Papers 13 and 15 inherited that motivation. The residual was instead the sigma-only restriction — a phi-independent basis spans only m1 = m2 = 0, while a 1Sigma_g+ state constrains only the total M. Restoring the azimuthal channels in the *same* prolate-spheroidal basis reaches **99.09% of D_e**, against Level 4's 96.0% at l_max = 6 with a Schwartz cusp correction; an independent Gaussian route gives 99.10%, and a published grid-based prolate-spheroidal calculation reaches 99.97% (Tao-McCurdy-Rescigno, *Phys. Rev. A* **82**, 023423 (2010)). What Level 4 retains, independently of the withdrawn claim, is *structural*: it is the exact N-electron generalization (SO(3N), S_N antisymmetry) that Level 4N and the composed Level 5 build on, and its cusp sits at a geometry-independent location. Those are architecture properties, not an accuracy advantage, and the hierarchy should not be read as asserting one here. **Open, and deliberately not decided by the PM: whether H2 belongs at Level 2 rather than Level 4.** See CHANGELOG v5.11.18, Paper 12 Sec. "Restoring the Azimuthal Channels".

| 4N | LiH (2-center, 4e) | Full mol-frame hypersp. (SO(12)) | R_eq 63.5% (l_max=2, 2D variational; unbound D_e) | 17 |"""

with io.open(PATH, encoding="utf-8") as fh:
    text = fh.read()

if OLD not in text:
    print("FAILED: Sec. 5 Level-4 rows not matched")
    sys.exit(1)

with io.open(PATH, "w", encoding="utf-8") as fh:
    fh.write(text.replace(OLD, NEW, 1))

print("Sec. 5: Level-4 note added (no level, coordinate system or row changed)")
