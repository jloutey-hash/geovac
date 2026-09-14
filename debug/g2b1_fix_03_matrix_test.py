"""group2 Batch-1 remediation: stale claim-matrix rows, the failing test, and
the Herbst citation scope.

matrix-267 (P8) -- Structural Theorem row still says NO-TEST; the backing was
   backfilled (test_paper8_sturmian_structural.py 10 pass + test_molecular_sturmian
   10 pass, verified this run). -> BACKED-SOUND.
matrix-271 (P12) -- full V_ee row still says NO-TEST/BACKED-WEAK with the
   tautological-B_l note; superseded by genuine independent anchors + two-route
   V_ee exactness (verified this run). -> BACKED-SOUND, KW literal note.
test-fix (P8) -- test_overlap_conserves_m NaNs on an m=1 self-overlap; the m!=0
   overlap quadrature is broken near the coordinate corner, but the paper uses
   only m=0 overlaps. Rewrite the sanity leg onto the m=0 path (which works) and
   document the m!=0 limitation. (The test never built a cross-m element anyway.)
herbst (P8) -- HerbstAveryDreuw2019 (PRA 99, 012512) is atomic Coulomb-Sturmian
   HF, not molecular binding; soften the citation clause.

Write-first; raw strings; idempotent.
"""
from __future__ import annotations
import sys

MTX = "docs/claim_test_matrix.md"
TEST = "tests/test_paper8_overlap_cross_n.py"
P8 = "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex"

EDITS = [
    (MTX, "matrix-267", "verified 2026-09-13 (group2 baseline): R-independent",
     "| — (only constituent E₀=−p₀²/2 tested) | **NO-TEST** | coverage gap — backfill a matrix-level H∝S + R-independence test |",
     "| `test_paper8_sturmian_structural.py` (10) + `test_molecular_sturmian.py` (10) | **BACKED-SOUND** | verified 2026-09-13 (group2 baseline): R-independent eigenvalues bit-exact across R∈{1.5,3,5,8}, congruence pullback UᵀHU=H₀, falsifiability control (`w_break`) + D-variation guard. Scoped to the single-n D-matrix approximation; the genuine cross-n overlap escapes it (`test_paper8_overlap_cross_n`) |"),

    (MTX, "matrix-271", "upgraded 2026-09-13 (group2 baseline)",
     '| only 10%-tol vs inaccurate grid; B_l test tautological (quad-vs-quad) | **NO-TEST / BACKED-WEAK** | gap — machine-precision V_ee-matrix test vs independent 4D quad (§13.4a) + 92.4% recompute test; B_l independent check. ("no integration of any kind" reframed v-this-sprint: B_l uses scipy.quad) |',
     "| `test_neumann_vee.py` (32): `test_vee_element_machine_precision` (<1e-8 vs independent Legendre assembly) + `test_vee_literal_4d_quad_normalization` (<2% vs raw 4D Coulomb) + `test_h2_de_headline_924pct` | **BACKED-SOUND** | upgraded 2026-09-13 (group2 baseline): B_l backing no longer tautological — three independent algebraic anchors (recurrence, E₁/γ closed form, Q₁=ξQ₀−1); full V_ee exact by two independent routes; 92.4% reproduces to 6 digits. KW reference literal corrected −1.17475→−1.174475 (the D_e% column already used the right denominator) |"),

    (TEST, "test-m0-fix", "m=0 self-overlap should be nonzero",
     '    same_m = _sturmian_overlap(1.0, 2, 1, 2, 1, 1, 3.015)    # <chi_2p+1|chi_2p+1>, m=1\n'
     '    assert abs(same_m) > 1e-6, "m=1 self-overlap should be nonzero (sanity)"',
     '    # m-conservation is STRUCTURAL: _sturmian_overlap takes ONE shared m, so it\n'
     '    # couples only equal m by construction (this is why Loewdin stays m-block).\n'
     '    # Sanity on the m=0 path the paper actually uses (1s/2s/2p0 overlaps).\n'
     '    # (The m!=0 overlap quadrature NaNs near the xi->1, eta->+-1 corner -- a\n'
     '    # known limitation of _overlap_two_center for m!=0; the paper Remark uses\n'
     '    # only m=0 overlaps, so no claim depends on it.  Fixed 2026-09-13.)\n'
     '    same_m = _sturmian_overlap(1.0, 1, 0, 1, 0, 0, 3.015)    # <chi_1s|chi_1s>, m=0\n'
     '    assert abs(same_m) > 1e-6, "m=0 self-overlap should be nonzero (sanity)"'),

    (P8, "herbst-scope", "consistent with published Coulomb--Sturmian",
     "consistent with the\npublished binding of molecular Coulomb--Sturmian calculations at\nHartree--Fock and beyond~\\cite{HerbstAveryDreuw2019}.",
     "consistent with published Coulomb--Sturmian Hartree--Fock\ncalculations~\\cite{HerbstAveryDreuw2019} (there for atoms; the molecular\ncross-$n$ binding is the escape this Remark establishes)."),
]


def main() -> int:
    loaded: dict[str, str] = {}
    applied, skipped, missed = [], [], []
    for path, name, marker, old, new in EDITS:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
            skipped.append(name); continue
        if t.count(old) != 1:
            missed.append((name, t.count(old))); continue
        loaded[path] = t.replace(old, new); applied.append(name)
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    for n in applied: print(f"  ok    {n}")
    for n in skipped: print(f"  skip  {n} (already applied)")
    for n, c in missed: print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
