"""The GeoVac claim registry -- load-bearing CLAIMS as linked objects.

WHY THIS EXISTS (2026-09-30, PI-directed after /qa trunk FULL #11)
------------------------------------------------------------------
The numeric registry (``numeric_registry.py`` + C21) made number drift nearly
extinct by turning numbers into linked objects with one owner.  Claims never
got the same treatment, and the measured cost is large: across DELTAs #13-#23
and FULLs #9-#11 the dominant LARGE class was a single shape --

    a claim is corrected (re-tiered, scoped, withdrawn) in its OWNER, and a
    restatement elsewhere keeps the old strength or drops the new qualifier.

C16 (``check_retracted_terms.py``) is a BLOCKLIST: it names known-bad
phrasings and refuses to let them re-surface.  It cannot see the MIRROR
direction -- a restatement that is merely missing the owner's scope
qualifier ("Paper 38 ... (unconditional)" with no "commutative factor";
"case-exhaustion ... INTERNAL THEOREM" when the owner says CONDITIONAL).
That is what this registry adds: each claim declares its owner, its
current tier, the QUALIFIER every restatement must carry, and the
OVERCLAIM wording no restatement may use.

WHAT THE GATE (C25, ``check_claim_registry.py``) DOES
-----------------------------------------------------
  1. Owner check -- the owner still carries the anchor and its tier tag
     (a registry that drifted from its owner is itself a defect).
  2. Front-door surfaces (``FRONT_DOOR``: the INDEX and the claims register,
     which newcomers are told to read first) -- every line that mentions a
     claim must carry its qualifier and must not carry its overclaim.  FAIL.
  3. Annotated loci -- ``\\gvc{<id>}{<text>}`` in LaTeX (renders <text>, like
     ``\\gvq``) and ``<!-- gvc:<id> -->`` in Markdown are explicit citations
     of a claim; the annotated text is held to the same rule.  FAIL.
  4. ``--mentions <id>`` -- every locus, corpus-wide (papers, syntheses,
     front-door docs, tests, geovac docstrings), that mentions the claim,
     with whether its qualifier is nearby.  This is the READING LIST a
     claim-search DELTA reviewer works from (qa.md, DELTA shape).
  5. ``--changed-since YYYY-MM-DD`` -- the claims whose status moved; the
     DELTA's claim-impact set is nominated from here.

MAINTENANCE RULE
----------------
When a claim's tier, scope or status changes: edit its entry HERE FIRST
(``tier``, ``qualifier``, ``overclaim``, ``changed``), run C25, fix every
front-door / annotated locus it names, then hand ``--mentions <id>`` to the
DELTA.  Never register a claim whose owner text you have not read; the
``owner_anchor`` + ``owner_tier`` pair is checked every run so the registry
cannot silently disagree with the paper it points at.

Patterns are Python regexes, matched case-insensitively against a LINE
(Markdown front-door rows are single lines) or against the annotated span.
"""
from __future__ import annotations

FRONT_DOOR = [
    "papers/INDEX.md",
    "docs/claims_register.md",
    # added 2026-10-02 (/qa trunk DELTA #25): README:201 carried a section-13.5 K-tier
    # violation and forced_free_seam.md (cited by Paper 57 as "the full table") carried
    # retired rows, and neither was on the front door, so C25 never read them.
    "README.md",
    "docs/forced_free_seam.md",
]

P32 = "papers/group1_operator_algebras/paper_32_spectral_triple.tex"
P38 = "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex"
P7 = "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex"
P0 = "papers/group3_foundations/Paper_0_Geometric_Packing.tex"
P2 = "papers/group5_qed_gauge/paper_2_alpha.tex"
P55 = "papers/group3_foundations/paper_55_periods_of_geovac.tex"

CLAIMS = {
    "pi-source-case-exhaustion": dict(
        statement="Every pi in a chain of the first fifteen Paper 34 projections engages M1/M2/M3: "
                  "(i)<=>(ii) empirical (Paper 35) + (ii)=>(iii) case by case; the converse "
                  "(iii)=>(i) is withdrawn as NOT ESTABLISHED.",
        owner=P32, owner_anchor=r"label\{thm:pi_source_case_exhaustion\}",
        owner_tier=r"\[CONDITIONAL\]", tier="CONDITIONAL",
        mention=r"case[- ]exhaustion|pi-source|π-source",
        qualifier=r"\bconditional\b|first fifteen|first 15|III\.1",
        overclaim=r"INTERNAL THEOREM|following are equivalent|classification is complete|"
                  r"stated and proven|does not force|every appearance of .{0,10}pi.{0,10} in any|"
                  r"complete (period-theoretic |period-ring )?stratification|classif(ies|y) the transcendentals|"
                  r"which transcendentals can appear|identifies as the \$?M_?\{?1",
        changed="2026-09-30"),
    "wh1-p38-gh-convergence": dict(
        statement="Truncations of the TRUE Camporesi-Higuchi window converge to the round-S^3 triple "
                  "in state-space GH (translation seminorm): COMMUTATIVE FACTOR ONLY; the almost-"
                  "commutative extension and Dirac-seminorm convergence are OPEN.",
        owner=P38, owner_anchor=r"label\{thm:main_unconditional\}",
        owner_tier=r"\[INTERNAL THEOREM", tier="INTERNAL THEOREM", owner_window=3,
        mention=r"Paper[~ ]?38\b.{0,60}convergence|discrete.?→.?continuum convergence|"
                r"state-space G(romov|H).{0,40}convergence",
        qualifier=r"(?<!almost-)(?<!almost )(?<!non-)(?<!non )\bcommutative\b",
        overclaim=r"almost[- ]commutative[^.;]{0,40}\b(is|are|been)\s+(also\s+)?(proved|proven)|"
                  r"almost[- ]commutative[^.;]{0,60}\b(proved|proven)\s+(likewise|too|as well)|"
                  r"\b(proves|proved)\s+(the\s+)?(convergence\s+(of|for)\s+the\s+)?(almost[- ]commutative|non-?commutative gauge)|"
                  r"graph.{0,30}converges? to|\\le\s*C_\{?3\}?\s*(\\cdot\s*)?\\gamma",
        changed="2026-09-29"),
    "p38-rate-constant": dict(
        statement="The displayed rate (4/pi) log n/n is [INTERNAL THEOREM + MEASURED] in the dual-"
                  "Coxeter normalisation (2/pi unit S^3; 2 sqrt2/pi Kac).",
        owner=P38, owner_anchor=r"4/\\pi",
        owner_tier=r"INTERNAL THEOREM \+ MEASURED", tier="INTERNAL THEOREM + MEASURED",
        mention=r"4/π.{0,40}log|\(4/\\pi\)|4/π\)·log|rate .{0,15}4/π|c = 4/π",
        qualifier=r"MEASURED",
        overclaim=r"4/π (is )?canonical",
        changed="2026-09-03"),
    "graph-s3-wall": dict(
        statement="The lattice edge Laplacian D-A does NOT converge to Delta_{S^3} (n_max disconnected "
                  "l-blocks; WALL); what is exact is the (n,l,m) label set + diagonal SO(4) Casimir.",
        owner=P7, owner_anchor=r"\[WALL\]",
        owner_tier=r"\[WALL\]", tier="WALL",
        mention=r"graph Laplacian|edge Laplacian|D\s*[−-]\s*A\b",
        qualifier=r"not|WALL|disconnected|refuted|wall|band edge|block",
        overclaim=r"converges to (the )?(S|Δ|\\Delta)|reproduces the (hydrogen|Rydberg)",
        changed="2026-09-25"),
    "k-combination-rule": dict(
        statement="K = pi(B + F - Delta) ~ alpha^-1 is an OBSERVATION (numerical coincidence); never "
                  "conjecture, derivation, prediction or theorem (CLAUDE.md 13.5).",
        owner=P2, owner_anchor=r"label\{obs:cross_sector\}",
        owner_tier=r"numerical observation", tier="OBSERVATION", owner_window=12,
        mention=r"π\(B ?\+ ?F|\\pi\(B ?\+ ?F|combination rule",
        qualifier=r"Observation|OBSERVATION|coincidence",
        overclaim=r"\b(conjecture|prediction)\b|(?-i:\b(?:rule|K)\b)[^.;]{0,40}\b(is|as)\s+(a\s+)?theorem(?![-\w])|is derived|derives K|"
                  r"is now a theorem|(?-i:\b(?:rule|K)\b)[^.;]{0,60}\b(is\s+)?proved\s+from",
        changed="2026-06-14"),
    "sm-gauge-group-upgrade-b": dict(
        statement="U(1) x SU(2) x SU(3) (and the inner algebra C+H+M3(C)) is FORCED ONLY RELATIVE TO "
                  "the adopted Upgrade B axiom, via an organizing reading -- not a derivation.",
        owner=P32, owner_anchor=r"admitted but not forced",
        owner_tier=r"adopted axiom", tier="ORGANIZING READING + ADOPTED AXIOM",
        mention=r"U\(1\)\s*(?:[×x]|\\times)\s*\\?(?:mathrm\{)?SU\(2\)|inner algebra|"
                r"\\mathbb\{C\}\s*\\oplus\s*\\mathbb\{H\}|Upgrade~?B",
        qualifier=r"Upgrade|organizing reading|adopted|not a derivation|admitted",
        overclaim=r"mostly forced|derived from the packing|construction does reach|autonomously generat|"
                  r"relative theorem of necessity|disjointness forbids",
        changed="2026-09-30"),
    "central-smoothing-cb-norm": dict(
        statement="The central smoothing maps are UCP, so their cb-norm is 1; 2/(n_max+1) (and "
                  "4/((n_a+1)(n_b+1))) are Plancherel MASS maxima, not cb-norms.",
        owner=P38, owner_anchor=r"UCP map has cb-norm",
        owner_tier=r"exactly \$1\$", tier="INTERNAL THEOREM", owner_window=2,
        mention=r"cb-norm",
        qualifier=r"cb-norm[^.;]{0,30}(=|is|of)\s*\$?1\b(?!\))|UCP|unital completely positive|mass max|legacy",
        overclaim=r"cb-norm\s*(=|is)?\s*2\s*/|cb-norm .{0,25}(decay|O\(1/n)|cb-norm[^.;]{0,40}2\s*/\s*\(n",
        changed="2026-09-30"),
    "m2-ring-integer-s": dict(
        statement="M2's pure-Tate ring (+)_k pi^{2k} Q holds for Seeley-DeWitt coefficients and "
                  "spectral-zeta values at INTEGER s; half-integer s >= 5/2 carries odd zeta.",
        owner=P55, owner_anchor=r"label\{thm:m2_mixed_tate\}",
        owner_tier=r"(?<!half-)(?<!half )\binteger\b", tier="THEOREM (integer s)",
        mention=r"pure-Tate even-weight|⊕.?π\^\{?2k|\\bigoplus_k \\pi\^\{2k\}",
        qualifier=r"(?<!half-)(?<!half )\binteger\b|Seeley|SD coeff",
        overclaim=r"at (any|every) s\b|half-integer s.{0,20}pure-Tate|pure-Tate[^.;]{0,40}(at|for) (every|all|any) half-integer|non-integer s",
        changed="2026-09-30"),
    "paper7-symbolic-checks": dict(
        statement="Paper 7's sympy suite is 18 CHECKS of Fock's continuum chain (ten limiting-case, "
                  "definitional or non-discriminating), not 18 independent proofs.",
        owner=P7, owner_anchor=r"18 symbolic checks",
        owner_tier=r"ten of", tier="SYMBOLIC (checks)",
        mention=r"\b18\b.{0,20}symbolic|symbolic.{0,20}\b18\b",
        qualifier=r"check|limiting|definitional",
        overclaim=r"18 (independent )?(symbolic )?proofs|18/18 symbolic proofs|independent proofs",
        changed="2026-09-30"),
    "p7-parseval-f0": dict(
        statement="The single-integral F^0 density-overlap formula is the standard Fourier-Parseval form "
                  "of the Coulomb integral written in the stereographic variable t = q/(2Z) [PRIOR ART]; "
                  "the single-integral structure is Parseval's, not a consequence of S^3 geometry.",
        owner=P7, owner_anchor=r"Fourier--Parseval",
        owner_tier=r"\[PRIOR ART", tier="PRIOR ART",
        mention=r"density[- ]overlap|single integral over one|F\^0.{0,40}S\^3|Slater integrals on \$S\^3",
        qualifier=r"Parseval|PRIOR ART|prior art|stereographic variable",
        overclaim=r"key structural finding|Fock projection converts|consequence of the \$?S\^3",
        changed="2026-09-30"),
    "he-bracket": dict(
        statement="The He fixed-k CI (1.88% above exact) and the cusp-corrected 0.004% (a Schwartz "
                  "extrapolation marginally BELOW exact, non-variational) bracket the exact energy; the "
                  "properly variational 2D figure is 0.022% raw.",
        owner="papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
        owner_anchor=r"marginally below", owner_tier=r"non-variational", tier="MEASURED", owner_all=True, owner_window=1,
        mention=r"0\.004\s*\\?%",
        qualifier=r"extrapolat|non-variational|below the exact|below exact",
        overclaim=r"do not bracket|variational 0\.004|0\.004\s*\\?%[^.;]{0,30}variational bound",
        changed="2026-09-30"),
    "p1-degree-not-monotone-l": dict(
        statement="On the binary label lattice node degree depends on m and on which shell neighbours "
                  "exist, not monotonically on l (2p m=0 -> 3 > 2s = 2 = 2p m=+-1) [MEASURED].",
        owner="papers/group3_foundations/paper_1_spectrum.tex",
        owner_anchor=r"emph\{not\} monotone in \$l\$", owner_tier=r"\[MEASURED\]", tier="MEASURED",
        mention=r"(higher|high)[- ](angular momentum|\$?l\$?)[^.;]{0,40}(degree|connect)|more connections",
        qualifier=r"not monotone|m = 0|m=0|not\} monotone",
        overclaim=r"more connections|occupy higher-degree",
        changed="2026-09-30"),
    "p38-uniform-bound": dict(
        statement="gamma_n <= 6 log n / n is a finite-range monotonicity check [MEASURED], not proved for "
                  "all n_max.",
        owner=P38, owner_anchor=r"6\s*(\\,)?\s*\\log", owner_tier=r"\[MEASURED\]", tier="MEASURED",
        mention=r"6\s*(\\,)?\s*\\log\s*\\?n|6 log n",
        qualifier=r"MEASURED|finite[- ]range|not proved",
        overclaim=r"PANEL-VERIFIED|holds for all|valid for every|proved for (every|all) n",
        changed="2026-09-30"),
    "k-match-precision": dict(
        statement="K = pi(B + F - Delta) itself matches 1/alpha to 4.77e-7 (raw, CODATA 2018); the "
                  "8.8e-8 figure is the root of Paper 2's cubic 1/alpha = K - alpha^2, not K.",
        owner=P2, owner_anchor=r"label\{obs:cross_sector\}", owner_tier=r"cubic's root", tier="OBSERVATION",
        mention=r"8\.8\s*(×|\\times)\s*10|8\.8e-8",
        qualifier=r"cubic|root|\\alpha\^2|α²|K\s*[−-]\s*\\?alpha",
        overclaim=r"\bK\s*(=\s*(π|\\pi)\s*\(B[^)]*\)\s*)?(match(es)?|agrees? with)\s*(1/α|1/\\alpha|\$1/\\alpha\$)?\s*to\s*\$?8\.8|cubic is irrelevant",
        changed="2026-09-30"),
    "p40-general-g-conditional": dict(
        statement="Paper 40 generalises Paper 38 to compact connected semisimple Lie groups "
                  "unconditionally at rank 1 and, for general G, conditionally on the per-group "
                  "spin-window decomposition; scalar compact-group convergence is Gaudillot-Estrada--"
                  "van Suijlekom's.",
        owner="papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
        owner_anchor=r"(?<!un)conditionally on the per-group", owner_tier=r"for general \$G\$", tier="CONDITIONAL",
        mention=r"Paper~?\s?40\b[^.;]{0,120}(generali[sz]|compact|semisimple|universal)|all\s+(simple\s+)?compact\s+(connected\s+)?(semisimple\s+)?Lie\s+groups",
        qualifier=r"condition|rank[~ ]?1|rank 1|spin-window|named gap",
        overclaim=r"resolved by\s+Paper[~\s]?40|closed it on the universal|settled by the present author|"
                  r"machinery to all simple compact|closes it positively|was closed in\s+Paper[~\s]?40",
        changed="2026-10-02"),
    "p1-ladder-exact": dict(
        statement="Paper 1's ladder operators reproduce n on INTERIOR shells and l(l+1) on every shell "
                  "exactly [SYMBOLIC PROOF] (a two-line commutator identity, w_{n-1}^2 - w_n^2 = -n/2); not "
                  "on the truncated top shell; not [PANEL-VERIFIED]. Upgraded from [MEASURED] 2026-10-02.",
        owner="papers/group3_foundations/paper_1_spectrum.tex",
        owner_anchor=r"They do, exactly, on interior shells", owner_tier=r"\[SYMBOLIC PROOF\]",
        tier="SYMBOLIC PROOF", owner_window=1,
        mention=r"(reproduce|return)[^.;]{0,40}quantum numbers|ladder operators[^.;]{0,60}(exact|reproduce)",
        qualifier=r"interior|SYMBOLIC|commutator|n\s*<\s*n_",
        overclaim=r"PANEL-VERIFIED|on every shell[^.;]{0,20}\bn\b|including the top shell[^.;]{0,30}\bn\b",
        changed="2026-10-02"),
    "three-homes-negative": dict(
        statement="That no single spectral construction on S^3 generates B, F and Delta is an empirical "
                  "record (twelve mechanisms eliminated, none exhibited) -- an OBSERVATION, not part of "
                  "Paper 2's Three-homes theorem and not an exhaustion.",
        owner=P2, owner_anchor=r"No single spectral construction on", owner_tier=r"\[OBSERVATION\]",
        tier="OBSERVATION", owner_window=1,
        mention=r"(generat\w*|share)[^.;]{0,40}all three|single[- ](spectral )?(construction|mechanism|generator)",
        qualifier=r"OBSERVATION|empirical|eliminated|exhibited|record",
        overclaim=r"no spectral construction on \$?S\^3\$? (simultaneously )?generates|complete enumeration|"
                  r"three-homes'?\s+theorem\s+(establishes|proves)[^.;]{0,80}no\s+(single|spectral)",
        changed="2026-10-02"),
    "m1-gh-rate-reading": dict(
        statement="The GH-rate constant 4/pi shares its value with the M1 (k = 0) volume slot of the master "
                  "Mellin engine -- a READING, not a derivation (Paper 38); M1 does not produce, govern or "
                  "own the GH rate, and the case-exhaustion does not cover GH rates.",
        owner=P38, owner_anchor=r"appearance of this constant in", owner_tier=r"\[OBSERVATION\]",
        tier="OBSERVATION", owner_window=2,
        mention=r"(M_?1|M\$_\{1\}\$|k\s*=\s*0)[^.;]{0,80}(\bGH\b|Gromov|rate constant|propinquity rate)|"
                r"(\bGH\b|Gromov|rate constant)[^.;]{0,80}(\bM_?1\b|M\$_\{1\}|k\s*=\s*0\s*(volume|sub))",
        qualifier=r"reading|OBSERVATION|shares (its|the) value|membership|not a derivation",
        overclaim=r"mechanism produces the|natural domain of (state-space )?G|M1 owns|pure M1 expression|"
                  r"in \$?M_?1\$? alone|identifies (the state-space GH rate )?with the|"
                  r"exactly what the M1 mechanism predicts|closes it positively",
        changed="2026-10-02"),
    "p55-census-scope": dict(
        statement="Paper 55 is a CENSUS of S^3 spectral-zeta and heat-trace data: its ring memberships are "
                  "theorems, its slot indexing a reading; chemistry seeds, spectral determinants and Paper 61's "
                  "elliptic period lie outside it. Not a classification of all GeoVac transcendentals.",
        owner=P55, owner_anchor=r"Over the census, the answer is affirmative", owner_tier=r"census",
        tier="INTERNAL THEOREM (ring memberships) + OBSERVATION (scope)", owner_window=2,
        mention=r"(all|every)[^.;]{0,30}transcendental[^.;]{0,80}(mixed-Tate|period)|"
                r"Periods of GeoVac|period-theoretic class",
        qualifier=r"census|over the cases treated|spectral(-zeta)? data|reading|identified",
        overclaim=r"We show that the answer is affirmative|across the entire GeoVac corpus|full classification|"
                  r"engine classification|classification of the transcendentals|"
                  r"all of GeoVac'?s\s+transcendentals sit",
        changed="2026-10-02"),
}

TIER_TAGS = ["SYMBOLIC PROOF", "MEASURED", "PANEL-VERIFIED", "INTERNAL THEOREM", "CONDITIONAL",
             "OBSERVATION", "CONJECTURE", "RETRACTED", "PROOF-BY-ARGUMENT", "PRIOR ART", "WALL", "OPEN"]
