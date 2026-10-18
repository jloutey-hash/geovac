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
        statement="4/pi is the asymptotic constant of the FEJER amplitude's moment, gamma_n = (4/pi) "
                  "log n/n [INTERNAL THEOREM + MEASURED] (dual-Coxeter; 2/pi unit S^3; 2 sqrt2/pi Kac). It is "
                  "NOT the convergence rate: the lifted-state bound holds for any real central amplitude, "
                  "a Cesaro taper gives d_GH = O(1/n) [MEASURED], proof and optimal constant OPEN "
                  "(re-priced 2026-10-03, FULL #12).",
        owner=P38, owner_anchor=r"4/\\pi",
        owner_tier=r"INTERNAL THEOREM \+ MEASURED", tier="INTERNAL THEOREM + MEASURED",
        mention=r"4/π.{0,40}log|\(4/\\pi\)|4/π\)·log|rate .{0,15}4/π|c = 4/π|4/\\pi\$?\s*(rate|asymptotic constant)",
        qualifier=r"Fej",
        overclaim=r"4/(π|\\pi)\$? (is )?canonical|rigorous\s+(?:at\s+)?rank|rigorous\s+rank-1|rigorous\s+\$?4/|"
                  r"(convergence|GH) rate\s+(is|of|=|reads)\s+[^.;]{0,15}4/(π|\\pi)|"
                  r"converges? at (the )?rate[^.;]{0,10}4/(π|\\pi)|"
                  r"rate at which[^.;]{0,60}converge[^.;]{0,40}carries|"
                  r"4/(π|\\pi)\$?[^.;]{0,25}\b(sharp|optimal|forced)\b",
        changed="2026-10-03"),
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
                  "van Suijlekom's. The 4/pi universality is reported for SIMPLE G only (derived at "
                  "rank 1, numerically supported at rank >= 2); for several simple factors the constant "
                  "is open (semisimple corollary withdrawn 2026-10-02).",
        owner="papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
        owner_anchor=r"(?<!un)conditionally on the per-group", owner_tier=r"for general \$G\$", tier="CONDITIONAL",
        mention=r"Paper~?\s?40\b[^.;]{0,120}(generali[sz]|compact|semisimple|universal)|all\s+(simple\s+)?compact\s+(connected\s+)?(semisimple\s+)?Lie\s+groups",
        qualifier=r"condition|rank[~ ]?1|rank 1|spin-window|named gap",
        overclaim=r"resolved by\s+Paper[~\s]?40|closed it on the universal|settled by the present author|"
                  r"machinery to all simple compact|closes it positively|was closed in\s+Paper[~\s]?40|"
                  r"universal\s+(?:across|over)\s+(?:the\s+)?compact\s+(?:connected\s+)?(?:\*\*)?semisimple|"
                  r"rigorous\s+(?:at\s+)?rank",
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
        statement="The Fejer-moment constant (2/pi unit S^3; 4/pi dual-Coxeter) shares its value with the M1 "
                  "(k = 0) volume slot of the master Mellin engine -- a READING, not a derivation (Paper 38). "
                  "It is a kernel constant, not the GH convergence rate (a tapered amplitude gives O(1/n), "
                  "2026-10-03); M1 does not produce, govern or own the GH rate.",
        owner=P38, owner_anchor=r"volume slot of the master Mellin engine", owner_tier=r"\[OBSERVATION\]",
        tier="OBSERVATION", owner_window=5,
        mention=r"(M_?1|M\$_\{1\}\$|k\s*=\s*0)[^.;]{0,80}(\bGH\b|Gromov|rate constant|propinquity rate)|"
                r"(\bGH\b|Gromov|rate constant)[^.;]{0,80}(\bM_?1\b|M\$_\{1\}|k\s*=\s*0\s*(volume|sub))",
        qualifier=r"reading|OBSERVATION|shares (its|the) value|membership|not a derivation",
        overclaim=r"mechanism produces the|natural domain of (state-space )?G|M1 owns|pure M1 expression|"
                  r"in \$?M_?1\$? alone|identifies (the state-space GH rate )?with the|"
                  r"exactly what the M1 mechanism predicts|closes it positively|"
                  r"rate at which[^.;]{0,60}converge[^.;]{0,40}carries|"
                  r"(GH|convergence) rate[^.;]{0,30}(sits|lies|is read) in",
        changed="2026-10-03"),
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
    "lipschitz-c3-compression": dict(
        statement="The Lipschitz comparison constant is C_3 = 1 at every cutoff, by compression: single factor "
                  "(Paper 38), every rank (Paper 40), and jointly on the Clifford-doubled S^3 x S^3 module "
                  "(Paper 39, C_3^(2) = 1; k-fold sketch C_3^(k) = 1). It is a constant of the abandoned "
                  "(B,P) route, not a rate constant of the lifted-state theorems; the per-summand "
                  "Dirac-triangle / PRV analysis is not a check of it; per-irrep C_3 values are withdrawn.",
        owner="papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex",
        owner_anchor=r"The sharp statement is by compression on the", owner_tier=r"\[INTERNAL THEOREM\]",
        tier="INTERNAL THEOREM", owner_window=1,
        mention=r"C_\{?3\}?\^\{\(2\)\}|C_\{?3\}?\^\{\(k\)\}|C_3\(G\)|Lipschitz comparison constant",
        qualifier=r"compression|=\s*1|abandoned|withdrawn|not a rate constant|record",
        overclaim=r"C_3\^\{\(2\)\}\s*\\le\s*\\sqrt\s*\{?2\}?\s*(\$|\(Lemma|\\eqref\{eq:C3_2_full\}\) belongs)|"
                  r"asymptotic-tight \$?C_3|C_3\s*\\to\s*1\^-|\\nearrow\s*\\sqrt|independent per-summand check|"
                  r"genuine \$?\\sqrt\{?k",
        changed="2026-10-02"),
    "p39-state-space-not-propinquity": dict(
        statement="Paper 39 proves van Suijlekom state-space GH convergence of the truncated S^3 x S^3 "
                  "Camporesi-Higuchi tensor triple by the lifted-state route, d_GH <= gamma^(ab) <= "
                  "gamma_a + gamma_b with constant 1; Latremoliere propinquity is not claimed and the "
                  "propinquity question stays open.",
        owner="papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex",
        owner_anchor=r"the object established is van~Suijlekom", owner_tier=r"state-space",
        tier="INTERNAL THEOREM", owner_window=1,
        mention=r"Paper~?\s?39[^.;]{0,120}(converg|tensor|closes|closed|propinquity)|"
                r"tensor-product (state-space|propinquity)",
        qualifier=r"state-space|Gromov|\bGH\b|not claimed|abandoned|stays open|lifted-state",
        overclaim=r"closes that gap|propinquity of Paper~?39|Open prior to this work|"
                  r"Tensor-product propinquity convergence for|Paper 39 closes the same-manifold",
        changed="2026-10-02"),
    "p45-rates-spatial": dict(
        statement="Paper 45 proves a degeneracy theorem (the K^+-compressed Lorentzian seminorm vanishes), "
                  "not a propinquity convergence; its surviving spatial statement is d_GH <= gamma_SU "
                  "(constant 1); the product carrier converges at the additive rate gamma_SU + gamma_U; "
                  "the circle Fejer rate is O(T log N_t / N_t); its panel column is the spatial rate.",
        owner="papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex",
        owner_anchor=r"Our main theorem", owner_tier=r"fails structurally",
        tier="INTERNAL THEOREM (degeneracy) + spatial statement", owner_window=1,
        mention=r"Paper~?\s?45[^.;]{0,120}(rate|converg|propinquity|panel|proof)",
        qualifier=r"degenera|spatial|additive|state-space|descoped|withdrawn|retracted|not a",
        overclaim=r"Paper~?\s?45 propinquity proof|O\(T/\\Nt\)|SU\(2\) factor dominates the joint rate|"
                  r"standard Latr\S*\s+propinquity",
        changed="2026-10-02"),
    "p59-four-centre-genus-one": dict(
        statement="The general four-centre Slater ERI (XY|ZW) lives on a genus-one curve: a fourth "
                  "distinct centre adds no genus (no genus two) [MEASURED; agrees with an independent "
                  "McMurchie-Davidson evaluation]. Pre-registered 2026-10-13.",
        owner="papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
        owner_anchor=r"A fourth distinct center adds no genus", owner_tier=r"\[MEASURED\]",
        tier="MEASURED", owner_window=4,
        mention=r"four[- ]cent(er|re)[^.;]{0,80}(genus|elliptic|period|ERI)|fourth (distinct )?cent(er|re)",
        qualifier=r"MEASURED|measured|genus[- ]one|no genus two|adds no genus|elliptic",
        overclaim=r"four[- ]cent(er|re)[^.;]{0,60}\b(theorem|proved|proven)\b|"
                  r"genus[- ]two (curve|surface|period)s? (at|for) four",
        changed="2026-10-15"),
    "p59-three-electron-k3": dict(
        statement="The three-electron chain integral lives on a K3 surface whose holomorphic period is "
                  "a sum of two products of complete elliptic integrals ('elliptic squared'; BBBG 2008 eq. 33, from an identity of Bailey) [MEASURED / classical]; the TRUE fibre collapses to K3 "
                  "periods at equal scales (evidence); at generic scales a PSLQ search at one mass triple finds no "
                  "relation in a K3 + elliptic + Tate basis (evidence, not proof); Calabi-Yau n-folds for longer "
                  "trees are OPEN.",
        owner="papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
        owner_anchor=r"A three-electron chain lives on a K3", owner_tier=r"\[MEASURED / classical\]",
        tier="MEASURED / classical", owner_window=20,
        mention=r"Calabi.{0,3}Yau[^.;|]{0,15}folds?|generic (chain )?(fibre|fiber)|longer (trees|chains)",
        qualifier=r"OPEN|\bopen\b|conjectur|hunch|untested|not tested|no relation|evidence",
        overclaim=r"generic (fibre|fiber)[^.;]{0,40}\b(is|are|reduces to|lies in)\b[^.;]{0,20}\bK3\b|"
                  r"(proves?|proved|shows?|established?)[^.;]{0,30}Calabi.?Yau|"
                  r"Calabi.?Yau[^.;]{0,40}\b(proved|proven|established)\b",
        changed="2026-10-15"),
    "p35-window-circle-adopted": dict(
        statement="Paper 35 reads the observer's finite real-time window as a periodic Euclidean circle; "
                  "this identification is ADOPTED, not derived: a compact interval without periodic "
                  "identification keeps a continuous spectrum, and the KMS circle is Euclidean while the "
                  "window is real-time (Kempf 2010 is the published cap).",
        owner="papers/group6_precision_observations/paper_35_time_as_projection.tex",
        owner_anchor=r"identification the paper adopts", owner_tier=r"not one it derives",
        tier="ADOPTED IDENTIFICATION", owner_window=1,
        mention=r"periodic Euclidean (circle|time)|(observer|real-time|temporal|lifetime|finite)[^.;]{0,40}window[^.;]{0,60}\bcircle\b|"
                r"compactified time[^.;]{0,30}discrete",
        qualifier=r"adopt|not (one it )?deriv|rather than deriv|identif|interval|caveat|Kempf",
        overclaim=r"compactified time is (automatically )?discrete|"
                  r"window[^.;]{0,40}\b(forces|derives|implies)\b[^.;]{0,20}(periodic|circle|discrete)"
                  r"(?![^.;]{0,40}\b(open|unsettled|untested)\b)",
        changed="2026-10-16"),
    "seam-j4-j7-forced": dict(
        statement="Seam rows J4 (fourth centre: genus one) and J7 (three-electron chain: K3, elliptic "
                  "squared) are F (forced) as of 2026-10-15, on Paper 59 sec:beyond; J2 (front) and J8 "
                  "(r12 bridge, a validated reduction) and J9 (the rung rule, re-tiered F -> Obs 2026-10-16) are "
                  "Obs; J3/J6 (closed forms of the three-centre radial core) are W. J10 (four linked densities) is F as of 2026-10-16; J11 (a ring of three densities) is F as of "
                  "2026-10-17. Totals 41 F / 4 Obs / 1 A / 25 C / 3 W = 74.",
        owner="docs/forced_free_seam.md",
        owner_anchor=r"^\| J4 \|", owner_tier=r"\*\*F\*\*",
        tier="F (seam verdict)", owner_window=0,
        mention=r"\bJ[47]\b",
        qualifier=r"\bF\b|forced|FORCED|promot|Paper~?\s?59|genus|K3",
        overclaim=r"\bJ[47]\b[^.;]{0,15}\b(is|are|stays?|remains?)\s+(an?\s+)?\**Obs\b|"
                  r"\bJ[47]\b[^.;]{0,15}\**Obs\**\s*\|",
        changed="2026-10-15"),
    "p59-coulomb-ring-elliptic": dict(
        statement="A Coulomb ring of three linked densities does not climb: at the 3-D fingerprint fibre "
                  "the triangle's operator factors exactly as M o (theta - r) with an algebraic piece and M "
                  "equivalent to Legendre's equation, so it sits on the elliptic rung, below the chain's K3; "
                  "at one centre it carries dilogarithms where the chain is rational [SYMBOLIC + MEASURED]. "
                  "Scope: Coulomb kernels; the drop is the exact (unique) integration of the loop, and the "
                  "correlated (geminal) triangle does not drop.",
        owner="papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
        owner_anchor=r"A fourth finds that a Coulomb ring of three densities does not climb", owner_tier=r"\[SYMBOLIC \+ MEASURED\]",
        tier="SYMBOLIC + MEASURED", owner_window=1,
        mention=r"ring of three|rings? (does|do) not climb|Coulomb (ring|triangle)|"
                r"triangle of three\s+densit|loop (lowers|integrates out)|ring behaves|"
                r"closing [^.;]{0,30}into a ring|loop trades geometry|"
                r"(rung|class) is set by how many densities|kernel sharing a vertex lifts",
        qualifier=r"Coulomb|uniqueness|unique|integrates out|fibre|fiber|correlated|geminal",
        overclaim=r"\b(rings|loops) (never|cannot|can ?not) climb|"
                  r"\b(any|every) (ring|loop)[^.;]{0,40}(lowers|elliptic rung|drops)|"
                  r"(?<!Coulomb )\bring of three densities\b[^.;]{0,30}\bsits\b|"
                  r"(?<!Coulomb )\btriangle of three densities sits|"
                  r"\ba ring does not climb\b|"
                  r"loop lowers the rung when it integrates out exactly|"
                  r"whose (massive )?loop does not( integrate out( exactly)?)?\b|"
                  r"(rung|class) is set by how many densities one integral links",
        changed="2026-10-17"),
    "p59-correlated-triangle-k3": dict(
        statement="The correlated triangle (exponential geminals e^(-gamma r) on the two edges at one "
                  "density) sits with the three-density chain on the K3 rung exactly at gamma = 1 "
                  "(annihilator R3 o R1, R3 equivalent to F_3's operator) [SYMBOLIC + MEASURED]; at "
                  "gamma = 3/2 it reaches the K3 rung by local data mod p plus Levelt rigidity, and whether it climbs further is open. "
                  "One geminal shape; the cusp-correct and linear-exponential geminals are untested. "
                  "The registered elliptic prediction (QG2) was falsified.",
        owner="papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
        owner_anchor=r"the correlated triangle does not drop", owner_tier=r"\[SYMBOLIC \+ MEASURED\]",
        tier="SYMBOLIC + MEASURED", owner_window=1,
        mention=r"correlated (\(geminal\) )?triangle|geminal triangle|geminals on (the )?two edges|"
                r"exponential[- ]geminal triangle|geminal kite|stays K3",
        qualifier=r"gamma|\\gamma|exponential|one geminal|K3|chain|rigidity",
        overclaim=r"correlated (\(geminal\) )?triangle[^.;]{0,40}\b(sits|lies|stays) on the elliptic|"
                  r"\b(all|every) geminals?\b[^.;]{0,40}\bK3\b|"
                  r"\bK3\b[^.;]{0,40}\b(for|at) (all|every) (\$?\\?gamma|geminal)|"
                  r"local data (mod \$?p\$? )?at \$?\\?gamma ?= ?3/2\$?(?![^.]{0,120}(not excluded|open))",
        changed="2026-10-17"),
    "p59-ising-thresholds": dict(
        statement="At t = 1 the tree ladder is the Ising class c_{n,0} = int K0^n: for n <= 8 the finite "
                  "singular points are t^2 = (n-1-2j)^2 and no others, so t = 1 is singular exactly for "
                  "even n [SYMBOLIC + MEASURED]. The values (c_{3,0} CM elliptic, c_{4,0} = 4 pi^2 L(f,2), "
                  "c_{6,0} = 108 pi^2 L(g,4)) are PRIOR ART (BBBG; Lim-Tu-Yu); c_{8,0} fits no single "
                  "L-value tried, so the searches give parity no support as a rule for values; c_{5,0}, c_{7,0} fit no tested "
                  "lower basis at a stated bound. What the ladder adds is the threshold reading.",
        owner="papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
        owner_anchor=r"At \$t=1\$ the ladder is the Ising class", owner_tier=r"\[SYMBOLIC \+ MEASURED\]",
        tier="SYMBOLIC + MEASURED (thresholds); PRIOR ART (values)", owner_window=1,
        mention=r"Ising[- ]class|c_\{?[3-8],0\}?|even moments|parity (rule|heuristic)",
        qualifier=r"\\le ?8|<= ?8|up to \$?8|\$?n ?\\le|prior art|PRIOR ART|heuristic|not the rule",
        overclaim=r"(?<!heuristic that )\beven moments (always )?drop(?![^.;]{0,40}\b(is not|not the rule|heuristic)\b)|"
                  r"parity (is|gives) the rule for values|"
                  r"(we|this paper) (prove|establish|identify|find)s?[^.;]{0,40}\bL\(|"
                  r"for (every|all) \$?n\$?(?![^.;]{0,15}(up to|\\le|<=))[^.;]{0,40}(singular|threshold)|"
                  r"c_\{?8,0\}?[^.;]{0,40}\b(equals|is)\b[^.;]{0,20}\bL\(",
        changed="2026-10-17"),
    "p59-tree-banana-hadamard": dict(
        statement="The tree ladder's operators are the Hadamard products of the Legendre operator with "
                  "the equal-mass banana operators with n-1 propagators, in y = 1/(16 t^2), checked in "
                  "exact arithmetic for n = 2 to 5 [SYMBOLIC]; n = 4 is AESZ #58 and n = 5 the operator "
                  "printed with #39 (prior art: Almkvist-Zudilin; the tables). F_n is irreducible for "
                  "n = 2..5 by local exponents and Fuchs' relation.",
        owner="papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
        owner_anchor=r"So the tree ladder is tied to the banana ladder", owner_tier=r"\[SYMBOLIC\]",
        tier="SYMBOLIC (n = 2..5) + PRIOR ART", owner_window=0,
        mention=r"banana ladder|Hadamard product",
        qualifier=r"\$?n ?= ?2\$? to \$?5|n ?= ?2\.\.5|2 to 5|checked|Almkvist|prior art|PRIOR ART|AESZ",
        overclaim=r"(proved|proven|theorem)[^.;]{0,40}(banana|Hadamard)|"
                  r"(banana|Hadamard)[^.;]{0,40}\bfor (all|every) \$?n\b|"
                  r"\bwhole (tree )?ladder\b|"
                  r"five densities[^.;|]{0,60}\b(open|evidence|not settled)\b|relation to four is open",
        changed="2026-10-17"),
}

TIER_TAGS = ["SYMBOLIC PROOF", "MEASURED", "PANEL-VERIFIED", "INTERNAL THEOREM", "CONDITIONAL",
             "OBSERVATION", "CONJECTURE", "RETRACTED", "PROOF-BY-ARGUMENT", "PRIOR ART", "WALL", "OPEN"]
