r"""Register this round's withdrawals in the C16 retracted-claims gate, and add
the standardized withdrawal token at the loci that quote the retired wording.

HARD RULE (/qa, C16): whenever a run retires a claim, its phrase goes in the
REGISTRY and its `cited_by` dependents are declared.  Eight claims were retired
or scoped this round and none was registered -- which is exactly how the
p12-cusp-as-the-h2-gap dependent (Paper 18) went unnoticed for two rounds.

REGISTRY DISCRIMINATION RULE: every pattern below is proved in BOTH directions
by debug/firetest_round3_c16.py -- it must FIRE on the retired wording and stay
SILENT on the corrected wording.  An entry that fires on nothing is worse than
no entry, because the gate then reports PASS for a class nobody is guarding.

New entries use the standardized `[retracted YYYY-MM-DD: entry-id]` token rather
than hand-authored exemption vocabulary;  the registry's own note records that
the hand-authored class produced three false-clean entries on 2026-09-03.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

CHK = "debug/qa/check_retracted_terms.py"
P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
P13 = "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex"
P15 = "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex"

SY = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
MX = "docs/claim_test_matrix.md"
P18 = "papers/group3_foundations/paper_18_exchange_constants.tex"

NEW_ENTRIES = r'''    {
        "id": "p15-sigma-pi-decoupled",
        "note": "2026-09-14 (round-3 DELTA).  Paper 15 stated that the sigma "
                "(m=0) and pi (|m|=1) sectors are 'completely decoupled in the "
                "angular eigenvalue problem', on the grounds that the nuclear "
                "coupling is diagonal in m and the e-e coupling conserves "
                "M = m1 + m2.  WITHDRAWN: the reason is a non-sequitur -- "
                "(0,0) and (+1,-1) both carry M = 0, so conserving M is "
                "precisely what puts them in the SAME block, and the e-e "
                "multipole expansion couples them (Paper 12 says so explicitly: "
                "its m != 0 terms 'are what couple different mu').  The "
                "measurement agrees: decoupled blocks would make the ground "
                "state a minimum over blocks, so adding pi channels could not "
                "lower it by a near-constant amount, which is what is observed. "
                "The empirical offset survives as an observation; the "
                "decoupling explanation does not.",
        "pattern": r"completely decoupled in the angular"
                   r"|decoupling in the angular eigenvalue problem",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "reviewed 2026-09-14 -- owner, withdrawn at both loci",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": "reviewed 2026-09-14 -- its m != 0 coupling result is the contradicting evidence",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14 -- does not restate the decoupling claim",
            "docs/claim_test_matrix.md": "reviewed 2026-09-14",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/synthesis/geovac_field_guide.tex",
            "docs/claim_test_matrix.md",
        ],
    },
    {
        "id": "p13-hyperspherical-cusp-advantage",
        "note": "2026-09-14 (round-3 DELTA).  Paper 13 claimed at three loci "
                "that placing the e-e cusp on a coalescence manifold is a "
                "'crucial structural advantage' of hyperspherical coordinates "
                "and 'the key advantage over single-electron coordinate "
                "systems'.  WITHDRAWN as an ADVANTAGE claim: it is a structural "
                "DIFFERENCE.  Paper 15 records that it does not translate into "
                "an accuracy advantage, and Paper 12 reaches 99.1% of H2's D_e "
                "in prolate spheroidal coordinates, where the same locus is not "
                "a coordinate surface at all.  Papers 12 and 15 were corrected "
                "in round 2; Paper 13 was in the same review set and was not "
                "reached until round 3 -- the owner-corrected / citer-stale "
                "class again.",
        "pattern": r"crucial structural advantage of these coordinates"
                   r"|key advantage over single[- ]electron\s+coordinate systems"
                   r"|cusp at a boundary\s+condition rather than a coordinate singularity",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14 -- owner, all three loci corrected",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "reviewed 2026-09-14 -- carries the corrected framing already",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": "reviewed 2026-09-14",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/synthesis/geovac_field_guide.tex",
        ],
    },
    {
        "id": "p13-sparsity-causes-005",
        "note": "2026-09-14 (round-3 DELTA).  Paper 13 said the SO(6) angular "
                "sparsity 'is ultimately responsible for the 0.05% accuracy' of "
                "the single-channel adiabatic approximation.  WITHDRAWN: the "
                "same paper's closing section attributes that number to "
                "fortuitous cancellation between the adiabatic approximation "
                "error and the finite-difference discretization error.  The two "
                "readings are incompatible -- one makes the number a structural "
                "property of the framework, the other an artifact -- and the "
                "cancellation reading is the one the evidence supports.",
        "pattern": r"ultimately responsible for the 0\.05",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14 -- owner",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14 -- states the non-variational caveat correctly",
            "docs/validation_benchmarks.md": "reviewed 2026-09-14",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/validation_benchmarks.md",
        ],
    },
    {
        "id": "p13-graph-reproduces-hydrogenic",
        "note": "2026-09-14 (round-3 DELTA).  Paper 13 said 'The discrete graph "
                "Laplacian reproduces hydrogenic energies to < 0.1%'.  "
                "WITHDRAWN: this states as a measured accuracy what is a "
                "property of the construction.  E_0 = kappa * lambda_max holds "
                "BY CONSTRUCTION (CLAUDE.md Sec. 5, Level 1), so any such "
                "figure is a bound on the spectral truncation deficit, not an "
                "accuracy against experiment.  Same class as the 'reproduces "
                "the Coulomb degeneracy' paraphrase retired from Paper 1 in "
                "DELTA #10 (v5.10.1).",
        "pattern": r"reproduces\s+hydrogenic energies to",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis trunk",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14 -- owner",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
        ],
    },
    {
        "id": "p18-mu-rho-r-algebraic",
        "note": "2026-09-14 (round-3 DELTA).  Paper 18's sec:algebraic_curve "
                "bullet said the adiabatic eigenvalues mu(R) AND mu(rho,R) are "
                "algebraic functions whose R-dependence is algebraic at every "
                "truncation.  SCOPED: the construction given is a LINEAR MATRIX "
                "PENCIL H_0 + R*V_C, which is what supplies the single global "
                "P(R,mu) = 0.  The Level-4 sweep is not of that form, and the "
                "same section's own later text says so.  What the Level-4 "
                "eigenvalues are instead is left open rather than reclassified "
                "-- the round-1 and round-2 defects in this subsection were "
                "both replacement readings written in place of a withdrawal.",
        "pattern": r"mu\(\\?rho\s*,\s*R\)[^.]{0,70}are\s+algebraic",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group3",
        "cited_by": {
            "papers/group3_foundations/paper_18_exchange_constants.tex": "reviewed 2026-09-14 -- owner",
            "docs/algebraic_registry.md": "reviewed 2026-09-14",
        },
        "files": [
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "docs/algebraic_registry.md",
        ],
    },
    {
        "id": "p18-mu-needed-for-sub-01",
        "note": "2026-09-14 (round-3 DELTA).  Paper 18 said 'the mu(R) "
                "parameterization is needed to achieve sub-0.1% accuracy'.  "
                "WITHDRAWN -- false in both directions: the adiabatic route "
                "that CARRIES mu(R) floors at 0.19-0.20% and never reaches "
                "sub-0.1%, while sub-0.1% IS reached without it by the 2D "
                "variational solver (0.022% raw at l_max = 7), which treats R "
                "and alpha simultaneously.  The Class-S to Class-C reading in "
                "the following sentence does not depend on the necessity claim "
                "and stands.",
        "pattern": r"parameterization is needed\s+to achieve sub-0\.1",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group3",
        "cited_by": {
            "papers/group3_foundations/paper_18_exchange_constants.tex": "reviewed 2026-09-14 -- owner",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14 -- source of both floor figures",
        },
        "files": [
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
        ],
    },
    {
        "id": "p12-envelope-insensitive",
        "note": "2026-09-14 (round-3 DELTA).  Paper 12's stability-envelope "
                "paragraph -- itself written as the round-2 remediation -- said "
                "the conclusion that the azimuthal channels close essentially "
                "all of the 7.6% gap 'is insensitive to the choice', and then "
                "quoted the number that refutes it: the weakest variational "
                "value on the full (alpha, threshold) grid is 95.5%, which "
                "closes about 41% of the gap, not essentially all.  WITHDRAWN. "
                "The qualitative effect is robust; the quantitative value ranges "
                "95.5-99.1% across the grid.  Recorded because this is the third "
                "consecutive round in which the largest defect was in the "
                "previous round's own remediation.",
        "pattern": r"is insensitive to the choice",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": "reviewed 2026-09-14 -- owner",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14 -- carries the envelope, not the insensitivity claim",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
        ],
    },
    {
        "id": "p15-delta-cbs-above-97",
        "note": "2026-09-14 (round-3 DELTA).  Paper 15's remaining-gap budget "
                "said 'CBS extrapolation including delta would push the limit "
                "above 97%'.  WITHDRAWN: the +0.65 pp delta gain it rests on "
                "was measured against SIGMA-ONLY at l_max = 4, not against the "
                "sigma+pi baseline the budget uses, so the transfer is invalid. "
                "The row it comes from is itself unreconciled -- its N_ch = 37 "
                "implies sigma+pi+delta, which cannot return 87.6% when its own "
                "sigma+pi subset returns 93.6%.",
        "pattern": r"push the limit above 97",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "reviewed 2026-09-14 -- owner",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14 -- does not carry the >97% figure",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
        ],
    },
'''

# --- tokens at the loci that QUOTE the retired wording inside a withdrawal ----
TOKENS = [
    (P15,
     r"""\textbf{[WITHDRAWN 2026-09-14]} Earlier versions of this paragraph
stated that the $\sigma$ ($m = 0$) and $\pi$ ($|m| = 1$) sectors are""",
     r"""\textbf{[WITHDRAWN 2026-09-14]}
[retracted 2026-09-14: p15-sigma-pi-decoupled]  Earlier versions of this
paragraph stated that the $\sigma$ ($m = 0$) and $\pi$ ($|m| = 1$) sectors are""",
     "token: p15-sigma-pi-decoupled (body)"),

    (P15,
     r"""\textbf{[WITHDRAWN 2026-09-14]} Earlier versions read this as
consistent with complete $\sigma$--$\pi$ decoupling;\ the two sectors""",
     r"""\textbf{[WITHDRAWN 2026-09-14]}
[retracted 2026-09-14: p15-sigma-pi-decoupled]  Earlier versions read this
as consistent with complete $\sigma$--$\pi$ decoupling;\ the two sectors""",
     "token: p15-sigma-pi-decoupled (conclusion)"),

    (P15,
     r"""\textbf{[WITHDRAWN 2026-09-14]} Earlier versions carried that gain
    into this budget and concluded that CBS extrapolation including""",
     r"""\textbf{[WITHDRAWN 2026-09-14]}
    [retracted 2026-09-14: p15-delta-cbs-above-97]  Earlier versions carried
    that gain into this budget and concluded that CBS extrapolation including""",
     "token: p15-delta-cbs-above-97"),

    (P13,
     r"""\textbf{[SCOPE 2026-09-14]} Earlier versions called this a ``crucial
structural advantage''.""",
     r"""\textbf{[SCOPE 2026-09-14]}
[retracted 2026-09-14: p13-hyperspherical-cusp-advantage]  Earlier versions
called this a ``crucial structural advantage''.""",
     "token: p13-hyperspherical-cusp-advantage"),

    (P13,
     r"""\textbf{[WITHDRAWN 2026-09-14]} Earlier versions made the sparsity
responsible for the $0.05\%$ figure specifically.""",
     r"""\textbf{[WITHDRAWN 2026-09-14]}
[retracted 2026-09-14: p13-sparsity-causes-005]  Earlier versions made the
sparsity responsible for the $0.05\%$ figure specifically.""",
     "token: p13-sparsity-causes-005"),

    (P13,
     r"""\textbf{[SCOPE 2026-09-14]} Earlier versions said
    here that the discrete graph Laplacian ``reproduces hydrogenic
    energies to $<0.1\%$''.""",
     r"""\textbf{[SCOPE 2026-09-14]}
    [retracted 2026-09-14: p13-graph-reproduces-hydrogenic]  Earlier versions
    said here that the discrete graph Laplacian ``reproduces hydrogenic
    energies to $<0.1\%$''.""",
     "token: p13-graph-reproduces-hydrogenic"),

    (P12,
     r"""\textbf{[SCOPE 2026-09-14]} An earlier version of this sentence called
the conclusion ``insensitive to the choice'' and then quoted the figure""",
     r"""\textbf{[SCOPE 2026-09-14]}
[retracted 2026-09-14: p12-envelope-insensitive]  An earlier version of this
sentence called the conclusion ``insensitive to the choice'' and then quoted
the figure""",
     "token: p12-envelope-insensitive"),
]

applied, failed = [], []

# ---- 1. insert the registry entries ---------------------------------------
with io.open(CHK, encoding="utf-8") as fh:
    src = fh.read()

ANCHOR = '    {\n        "id": "p12-cusp-as-the-h2-gap",'
if "p15-sigma-pi-decoupled" in src:
    applied.append("registry entries already present -- skipped")
elif ANCHOR in src:
    src = src.replace(ANCHOR, NEW_ENTRIES + ANCHOR, 1)
    with io.open(CHK, "w", encoding="utf-8") as fh:
        fh.write(src)
    applied.append("8 registry entries inserted into check_retracted_terms.py")
else:
    failed.append("registry anchor not found")

# ---- 2. tokens -------------------------------------------------------------
by_path = {}
for path, old, new, label in TOKENS:
    by_path.setdefault(path, []).append((old, new, label))

for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    for old, new, label in items:
        if old in t:
            t = t.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append("%s   [%s]" % (label, path))
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)

print("applied %d" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("")
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
