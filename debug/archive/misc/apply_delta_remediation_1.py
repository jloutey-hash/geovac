"""DELTA remediation, pass 1: the loci outside the two still-running reviewers'
scope -- docs, README, the QA record, code comments, and the cited_by amendment.

Every item below was verified against primary text by the PM before editing.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# ---- 1. claim_test_matrix: a zombie ledger row instructing a future agent to
#         build backing for a WITHDRAWN claim.
edit(
    "docs/claim_test_matrix.md",
    '| 15 | "exceeds Paper 12" σ-only/σ+π | tests pass via *adiabatic* solver the paper disavows (σ-only asserts opposite) | **FALSE-POSITIVE** | fix tests to use the 2D variational solver |',
    '| 15 | ~~"exceeds Paper 12" σ-only/σ+π~~ | — | **WITHDRAWN 2026-09-14** | The claim itself is withdrawn, not under-tested: the comparison put Level-4 σ+π against Paper 12 σ-only. Do NOT build backing for it. Paper 15 now claims no coordinate-system advantage; the numerical band 92.4 < pct < 100 survives in `test_level4_multichannel.py` as a solver guard only |',
    "claim_test_matrix: zombie row re-tiered to WITHDRAWN")

# ---- 2. The frozen DoD watch-note would RE-CERTIFY the retired claim.
edit(
    "docs/qa/group2.done.md",
    """  - **Paper 12:** H$_2$ Neumann $V_{ee}$ recovers **92.4%** of $D_e$ vs 80.1% numerical;
    the **7.6%** gap is the cusp (a diagnosed limitation, stated as such). Algebraic
    recurrence = *exact*; the surviving transcendental seed is named.""",
    """  - **Paper 12:** H$_2$ Neumann $V_{ee}$ recovers **92.4%** of $D_e$ vs 80.1% numerical.
    **AMENDED 2026-09-14:** the 92.4% is the **σ-only** ceiling, and the 7.6% gap is
    the absent $m \\neq 0$ configurations, **not** the cusp — restoring them in the same
    basis reaches **99.09%**. The old watch-note ("the 7.6% gap is the cusp") is
    retired; grading a corrected paper against it would mark the correction as a
    defect. Algebraic recurrence = *exact* **for σ**; μ>0 uses spectral quadrature.
    The surviving transcendental seed is named.""",
    "group2.done.md: watch-note amended so it cannot re-certify the retired claim")

# ---- 3. README front door: understates Paper 12 (mirror direction).
edit(
    "README.md",
    "| **12** | Algebraic V_ee | Neumann expansion, H₂ 92.4% D_e |",
    "| **12** | Algebraic V_ee | Neumann expansion; H₂ 99.09% D_e (\\|m\\|≤1), 92.4% σ-only |",
    "README: Paper 12 headline no longer understated")

# ---- 4/5. Index + archive staleness: the diagnosis is no longer the cusp.
edit(
    "docs/topic_to_paper_lookup.md",
    "| Cusp diagnosis (7.6% gap) | 12 | Sec VII",
    "| Azimuthal-channel diagnosis of the 7.6% gap (the cusp reading is withdrawn) | 12 | Sec VII",
    "topic lookup: cusp-diagnosis label corrected")

edit(
    "docs/paper_notes_archive.md",
    "Neumann V_ee: H2 92.4% D_e, cusp diagnosis (7.6% gap)",
    "Neumann V_ee: H2 92.4% D_e σ-only / 99.09% at |m|<=1; the 7.6% gap is the absent m != 0 channels (the cusp diagnosis is withdrawn, 2026-09-14)",
    "paper notes archive: corrected + no longer understated")

# ---- 6/7. Production code printing a not-like-for-like superiority verdict.
for mod in ("geovac/level4_multichannel.py", "geovac/level4_sigma_channel.py"):
    edit(
        mod,
        'print(f"  ** IMPROVES on Paper 12 Neumann V_ee (92.4%) **")',
        'print(f"  ** above Paper 12 Neumann V_ee sigma-only (92.4%) -- NOT '
        'like-for-like: P12 |m|<=1 reaches 99.09% **")',
        "%s: superiority print scoped" % mod)

edit(
    "geovac/cusp_factor.py",
    "Paper 12, Section VII (cusp diagnosis: 7.6% D_e gap)",
    "Paper 12, Section VII (the 7.6% D_e gap; its cusp diagnosis was withdrawn "
    "2026-09-14 -- the gap is the absent m != 0 channels)",
    "cusp_factor.py: docstring pointer corrected")

# ---- 8. cited_by amendment: Paper 18 and the surfaces this run found.
edit(
    "debug/qa/check_retracted_terms.py",
    '''        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "reviewed 2026-09-14",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14",
            "docs/validation_benchmarks.md": "reviewed 2026-09-14",
            "docs/claim_test_matrix.md": "reviewed 2026-09-14",
            "tests/test_level4_multichannel.py": "reviewed 2026-09-14",
        },''',
    '''        # AMENDED after the 2026-09-14 DELTA run.  Paper 18 was MISSING from
        # the original list and is the largest dependent of all: its
        # `sec:mu_level4` builds a TAXONOMY classification on this claim.  Two
        # independent reviewers found it; no pattern could, because Paper 18
        # restates the claim in its own words and its own numbers (92.5%,
        # 7.5%).  Left unstamped deliberately -- the correction is a PI call.
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "reviewed 2026-09-14",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14",
            "papers/group3_foundations/paper_18_exchange_constants.tex": None,
            "docs/validation_benchmarks.md": "reviewed 2026-09-14",
            "docs/claim_test_matrix.md": "reviewed 2026-09-14",
            "docs/qa/group2.done.md": "reviewed 2026-09-14",
            "README.md": "reviewed 2026-09-14",
            "docs/topic_to_paper_lookup.md": "reviewed 2026-09-14",
            "docs/paper_notes_archive.md": "reviewed 2026-09-14",
            "tests/test_level4_multichannel.py": "reviewed 2026-09-14",
            "viz/public/papers/paper_12_algebraic_vee.html": None,
            "viz/public/papers/paper_13_hyperspherical.html": None,
            "viz/public/papers/paper_15_level4_geometry.html": None,
        },''',
    "cited_by: Paper 18 + viz pages added (UNSTAMPED -- PI calls)")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

applied, failed = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        text = fh.read()
    ok = True
    for old, new, label in items:
        if old in text:
            text = text.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append("%s  [%s]" % (label, path))
            ok = False
    if ok:
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write(text)

if failed:
    print("FAILED TO MATCH:")
    for f in failed:
        print("  -", f)
    sys.exit(1)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
