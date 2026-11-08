"""Record group2 baseline FULL run - Batch 3 (Paper 19 + FCI-atoms + FCI-molecules).
v5.11.16. Idempotent."""
from __future__ import annotations
import sys, pathlib

CL = "CHANGELOG.md"
ANCHOR = "## [v5.11.15] - 2026-09-13\n"
MARKER = "## [v5.11.16]"
ENTRY = """## [v5.11.16] - 2026-09-13

**`/qa group2` baseline FULL run - Batch 3 (Paper 19 balanced coupled + FCI-atoms + FCI-molecules) = FAIL, remediated.** Five reviewers, tree frozen. One genuine LARGE (stale post-ERI-fix energies with three RED backing tests) and one MATERIAL (a pair-diagonal zombie in Paper 19's prose contradicting its own tables); the rest SMALL. Every corrected number was re-measured this session, not copied from a reviewer.

### The LARGE - FCI-atoms energies stale since the 2026-08-29 ERI fix

Table I, the convergence-detail table, the abstract and the summary sentence all carried PRE-fix (wrong-sign-q) energies; three backing tests (`test_direct_ci` Li n_max=4/5, Be n_max=4) were pinned to the old values and had been RED since the fix (only the He pin was updated at the time). **Re-measured every row** (He hybrid-h1 via direct CI at the paper's own 5995-determinant convention; Li/Be exact-h1). N_SD is unchanged (basis size); energies are lower (restored m-changing correlation, all still above the variational bound); NNZ is ~3x higher (the same census that moved cross-block ERIs 130->214):

- He n_max=5: -2.8936 -> **-2.8963** Ha (0.35% -> **0.26%**)
- Li n_max=4: -7.3959 -> **-7.3987** (1.10% -> **1.06%**); n_max=5: -7.3978 -> **-7.4007** (1.07% -> **1.03%**)
- Be n_max=3: -14.531 -> **-14.558** (0.93% -> **0.75%**); n_max=4: -14.536 -> **-14.563** (0.90% -> **0.71%**)

**Be n_max=4 shifted 0.027 Ha - 10x the others - so it was independently cross-checked before re-pinning** (independent-route rule): the production grid-quad R_k and the exact analytical hypergeometric R_k agree to **1.5e-4**, and the exact route reproduces F^0(1s,1s)=5Z/8=2.500000 for Be exactly. The Be convergence is monotone (n3 -> n4: -14.558 -> -14.563), both above exact -14.6674. **The shift is physical (restored correlation), not an evaluator artifact.** All three RED tests re-pinned to this session's measured values and verified green (the two slow ones pass by construction - pinned to the identical production path's output).

### MATERIAL - Paper 19 pair-diagonal zombie (prose vs the paper's own tables)

Paper 19's re-measured tables (2026-09-01, exact global-M_L) carried the current figures, but the PROSE was never updated and cited retired pair-diagonal values that directly contradicted those tables - the trunk-level pair-diagonal C16 entry never reached this paper-specific prose. Corrected against the (measured) tables: **n_max=3 balanced 19,959 -> 127,855 Pauli / 448.9 -> 436.9 Ha 1-norm**; Pauli scaling exponent **3.03 -> 3.74**, 1-norm exponent 1.75 -> 1.71 (recomputed from the table values); LiH balanced/composed ratio **2.53x~2.63x -> 3.25x**; the polyatomic ratio sequence **2.63/4.77/7.45x -> 3.25/6.35/10.10x** (= the census table's own ratio column, which the prose sat directly below); ERI census **130/195 -> 214/321** (the GREEN `test_cross_block_eri_count` pin). A focused C16 entry (`pairdiag-p19-balanced-prose-values`) was added and **fire-tested** (FIRE on all four retired wordings, SILENT on all four corrected).

### Also remediated

- **DirectCI4e live-validation re-enabled.** The quarantine premise - "DirectCI4e's closed-form same-spin block assumes 8-fold ERI symmetry and computes a wrong energy on the exact-rule 4-fold tensor" - was tested and is **FALSE**: DirectCI4e(faithful=False) == coupled_fci_energy at **0.0000 mHa**. The energy is a scalar Slater-Condon contraction needing only the two physical symmetries (particle-exchange + hermiticity); the broken single-swap symmetry is the wrong-sign-q artifact and never enters. Leg restored, test green.
- **Citations - the NIST_ASD misattribution.** The exact non-relativistic total energies E_He/E_Li/E_Be were cited to the NIST Atomic Spectra Database, a spectroscopic (ionization-energy) compilation that does not tabulate total non-relativistic energies. Values correct; source wrong. Repointed to the verified primaries: **Pekeris (1958)** for He and **Chakravorty-Gwaltney-Davidson-Parpia-Fischer, Phys. Rev. A 47, 3649 (1993)** for Li and Be (covers 3-18 electrons). NIST_ASD was cited only here -> removed as an orphan. All three papers compile with zero undefined citations.
- **R_eq-drift consistency (M2).** The n_max=2->3 balanced R_eq drift was stated as +0.057 (abstract, Step-4-adjacent) and +0.053 (Step-4, and L2019 which is authoritative: 3.227 -> 3.280); and Step-4 framed it "+0.053 per step" (constant) against the abstract's "decelerating." Reconciled everywhere to **+0.053 then +0.023 bohr, decelerating**.
- **FCI-molecules (SMALL):** "the only sound discrete framework for heteronuclear FCI" scoped to "among the alternatives surveyed here" (M3); "the correct weighting" -> "a qualitatively correct weighting" (M4). The guardrail-negative is CLEAN (no summary drift; MolecularLatticeIndex archive-restored, 48 tests pass per the reviewer).
- **Coupled (CB) 1-norm test re-pin:** 91.65 -> **89.80** Ha (measured; matches the table's 80.5 non-identity + identity). Claim-matrix balanced-Pauli 878 -> 2726 (exact-rule; 878 retired).

### Carried to the group review (declared NITs)

herbst2018 orphan bibitem (P19); the uncited exact-LiH denominator -8.071 Ha (P19); the owed He 0.19%@n_max=7 graph-native NO-TEST/extrapolation footnote (shared with Paper 13's owed 0.19% footnote from Batch 2); FCI-molecules virial-T caveat and unpinned Wigner-D^2 literals; the claim-matrix composed-Pauli-334/333 row (retired-rule; exact-rule composed is 838, tapering-variant ambiguous); FCI-atoms convergence-table timing column (hardware-dependent, left as-is).

### Verdict

**Batch 3: FAIL, remediated** - 1 LARGE (FCI-atoms stale energies + 3 RED tests, all re-measured and Be independently cross-checked), 1 MATERIAL (P19 pair-diagonal prose), the rest SMALL. Deterministic layer green on group2 (C10 compile 3/3, C13/C14/C16/C17/C19/C21/C22 PASS); the C16 guard fire-tested; the re-enabled and re-pinned tests verified green. **Batch 4 (Paper 58 + group2 synthesis C9 + completeness-critic) remains to finish the group2 baseline.**

"""


def main() -> int:
    p = pathlib.Path(CL); t = p.read_text(encoding="utf-8")
    if MARKER in t:
        print("skip"); return 0
    if t.count(ANCHOR) != 1:
        print(f"MISS anchor count={t.count(ANCHOR)}"); return 3
    p.write_text(t.replace(ANCHOR, ENTRY + ANCHOR), encoding="utf-8")
    print("ok CHANGELOG v5.11.16")
    c = pathlib.Path("CLAUDE.md"); ct = c.read_text(encoding="utf-8")
    ct = ct.replace("**Version:** v5.11.15 (September 13, 2026)",
                    "**Version:** v5.11.16 (September 13, 2026)", 1)
    c.write_text(ct, encoding="utf-8"); print("CLAUDE.md v5.11.16")
    return 0


if __name__ == "__main__":
    sys.exit(main())
