# Sprint A4 Scoping memo (2026-06-03)

## TL;DR

**Verdict:\ CLOSED-BY-PRIOR with a one-line sharpening.** A4
(inner-factor mixed-Tate extension on the AC extension
$\mathcal{A}_{\mathrm{GV}} \otimes \mathcal{A}_F$) was scoped as
"multi-week" in Sprint Mixed-Tate Test (v3.45.3) BEFORE the late-2026-06
sub-sprints landed. Three subsequent results — the
$\eta$-trivialization theorem (Paper 18 §IV.6, Sprint H1 close), the
combined-spectral-action factorization theorem (Paper 18
`thm:ac_factorization`), and the Sprint Yukawa-PSLQ empirical
sweep — already fully determine A4's verdict at the structural level.
The combined Seeley--Dewitt coefficients on the AC extension factor
into (outer SD coefficients) $\times$ (inner Mellin output), where the
outer factor is pure-Tate over $\mathbb{Q}$ in
$\bigoplus_k \pi^{2k}\cdot\mathbb{Q}$ (Sprint Mixed-Tate Test
POSITIVE-with-refinement) and the inner factor is a parameter-tied
Dirichlet ring $\mathbb{Q}[y_i^{-2s}]$ whose generators are
empirically not in low-coefficient pure-Tate (Sprint Yukawa-PSLQ).
**The inner-factor SD coefficients are NOT mixed-Tate over
$\mathbb{Q}$** because the Yukawa eigenvalues that generate the
Dirichlet ring are themselves outside the mixed-Tate period ring.

A4 reduces to "one-paragraph closure statement plus Paper 55 §7.4
update plus register marker." Sprint-scale: 1--2 hours, not
multi-week.

## Three-option analysis

### Option (a) — CLOSED-BY-PRIOR — selected

**What the prior sprints already establish.** Three theorem-grade
results jointly determine A4:

1.  **Theorem `thm:eta_trivialization` (Paper 18 §IV.6).** On any
    Krajewski-class even spectral triple, $\{\gamma_F, D_F\} = 0$
    forces
    $\mathrm{Tr}(D_F^k\, e^{-t D_F^2}) \equiv 0$
    for every odd $k$. In particular the M3 sub-mechanism vanishes
    identically on the inner side. **Consequence:\ the inner side of
    the master Mellin engine reaches only the $k \in \{0, 2\}$ slots.**

2.  **Theorem `thm:ac_factorization` (Paper 18 §IV.6).** For the
    combined Dirac $D = D_{\mathrm{GV}} \otimes 1_F + \gamma_{\mathrm{GV}}
    \otimes D_F$, the cross-term $\{D_{\mathrm{GV}} \otimes 1,
    \gamma_{\mathrm{GV}} \otimes D_F\}$ vanishes by outer chirality
    anticommutation, so
    $D^2 = D_{\mathrm{GV}}^2 \otimes 1_F + 1_{\mathrm{GV}} \otimes D_F^2$
    and the combined heat trace factorizes:
    $\mathrm{Tr}\,e^{-tD^2} = \mathrm{Tr}\,e^{-tD_{\mathrm{GV}}^2}
    \cdot \mathrm{Tr}\,e^{-tD_F^2}$.
    **Consequence:\ the AC-extension Mellin output sits in (outer
    Mellin ring) $\otimes$ (inner Dirichlet ring), and the two factors
    do not mix.**

3.  **Empirical Sprint Yukawa-PSLQ (2026-06-03,
    `debug/sprint_yukawa_pslq_memo.md`).** 162-cell PSLQ sweep of the
    nine measured SM Yukawa values against the M1 $\cup$ M2 pure-Tate
    basis $\{1, \pi, \pi^2, \pi^4, \pi^6, \pi^8, 1/\pi, 1/\pi^2,
    1/\pi^4\}$ at coefficient ceiling $M \le 1000$, at both
    $\overline{\rm MS}\ M_Z$ and $\overline{\rm MS}\ 2\times 10^{16}$
    GeV, in three transforms ($y_f,\, y_f^2,\, \log y_f$), returned
    **zero hits**. **Consequence:\ the measured Yukawa values do not
    sit in low-coefficient mixed-Tate.**

**The structural verdict.** Per Eq.~(\ref{eq:inner_dirichlet}) of
Paper 18 §IV.6, the inner Mellin output is

$$
\mathcal{M}\bigl[\mathrm{Tr}(D_F^k\, e^{-tD_F^2})\bigr](s)
\;=\; \Gamma(s) \cdot \sum_i m_i \cdot y_i^{\,k - 2s},
\qquad k \in \{0, 2\}
$$

with rational multiplicities $m_i$. The output ring is
$\mathbb{Q}[y_1^{-2s}, y_2^{-2s}, \ldots]$ — a parameter-tied
Dirichlet ring. Two structural cases:

-   **Case 1 (the structural case).** If the Yukawa eigenvalues
    $y_i$ are algebraic over $\mathbb{Q}$ and the algebraic closure
    inside the mixed-Tate period ring captures the negative-power
    Dirichlet series, the inner factor inherits mixed-Tate. The
    framework provides NO such selection mechanism for $y_i$ (Sprint
    H1 Yukawa non-selection theorem; Forced-Count Theorem
    `thm:forced_count` confirms the moduli dimension is forced but
    the moduli point is not).
-   **Case 2 (the empirical case).** Sprint Yukawa-PSLQ confirms the
    measured Yukawa values are NOT in low-coefficient pure-Tate at
    $M \le 1000$. At the precision available from PDG (2--8 sig
    digits, charged leptons load-bearing at $M = 10$), no honest
    identification with $\mathbb{Q}[\pi, \pi^{-1}]$ exists.

**Consequence:\ the inner-factor SD coefficients are NOT mixed-Tate
over $\mathbb{Q}$ in any framework-derived sense, and are not
empirically mixed-Tate at the available precision.** The combined
AC-extension SD coefficients factor as
(pure-Tate outer) $\times$ (Yukawa-Dirichlet inner), and the inner
factor lives in a categorically disjoint ring.

### Option (b) — NEEDS-DISTINCTION — partially relevant

There IS a structural stratification, but it doesn't reopen A4 as a
multi-week sprint. The stratification is one-line:

-   **Gauge-sector inner-factor content** ($D_F = 0$, pure $\mathcal{A}_F$
    structure with chirality and $J$ but no Yukawa data). On this
    slice the inner Mellin output collapses to the $k = 0$ trace
    $\mathrm{Tr}\,e^{-t \cdot 0} = \dim \mathcal{H}_F$, a finite
    rational integer. The SD coefficients become
    $\dim \mathcal{H}_F$ times the outer pure-Tate values, hence
    still pure-Tate. **This slice IS mixed-Tate over $\mathbb{Q}$,
    trivially.** (And it's the Marcolli--vS-without-Higgs reading,
    empirically falsified by the sphaleron transition per Paper 32
    §VIII.C.)
-   **Yukawa-dependent inner-factor content** ($D_F \neq 0$). On this
    slice the inner Mellin output is $\sum_i m_i \cdot y_i^{-2s}$
    (and $y_i^{2-2s}$ for $k = 2$). **NOT mixed-Tate** unless $y_i$
    are themselves mixed-Tate, which they empirically are not.

The stratification is **forced/free at the Mellin engine level**:
the gauge-sector slice is pure-Tate by inheritance from the outer
factor; the Yukawa-dependent slice is calibration-tied and
categorically disjoint. This is one line of structural framing, not
a sprint-scale question.

### Option (c) — MULTI-WEEK STANDS — rejected

The original A4 framing (Sprint Mixed-Tate Test, open question 4)
assumed that the inner-factor SD coefficients had not yet been
explicitly computed and that an explicit CC spectral-action
computation on $\mathcal{A}_{\mathrm{GV}} \otimes \mathcal{A}_F$ with
the Higgs/Dirac internal operator $D_F$ was needed. After
`thm:ac_factorization` (Paper 18) was proved (Sprint TD Track 3 close,
late May 2026), the inner trace IS explicitly given by the
Dirichlet sum in $y_i^{-2s}$; no further computation is needed. The
question reduces to "are the Yukawas in mixed-Tate?", which Sprint
Yukawa-PSLQ answered empirically. Hence (c) is not the right scope.

The residual genuinely-open content is the **two-loop / one-loop
mixing with the outer factor** (does the M2 / inner-Dirichlet
product produce any new identifiable closed forms?), but that is
NOT A4 as scoped. It's a separate question, named below as open
follow-on.

## Structural argument (the load-bearing observation)

Putting the three theorem-grade ingredients together gives the
A4 verdict in one paragraph.

**Lemma (consequence of Theorems `thm:eta_trivialization` +
`thm:ac_factorization`).** Let
$\mathcal{T}_{\mathrm{AC}} = \mathcal{T}_{\mathrm{GV}} \otimes
\mathcal{T}_F$ be a Connes--Chamseddine almost-commutative extension
with $\mathcal{T}_F = (\mathcal{A}_F, \mathcal{H}_F, D_F, J_F,
\gamma_F)$ Krajewski-class. The Seeley--Dewitt coefficients of the
combined spectral action
$\mathrm{Tr}(f(D / \Lambda))$ factor multiplicatively:
$a_k^{\mathrm{AC}} = \sum_{j} a_{k-j}^{\mathrm{outer}} \cdot
b_j^{\mathrm{inner}}$
where $a_k^{\mathrm{outer}} \in \bigoplus_m \pi^{2m} \cdot
\mathbb{Q}$ (Sprint Mixed-Tate Test) and
$b_j^{\mathrm{inner}} \in \mathbb{Q}[y_1^2, \ldots, y_n^2]$
(traces of even powers of $D_F$; the M3 / $k$-odd traces vanish by
$\eta$-trivialization). The combined SD ring is
$\bigl(\bigoplus_m \pi^{2m} \cdot \mathbb{Q}\bigr) \cdot \mathbb{Q}[y_1^2,
\ldots, y_n^2]$.

**Verdict on mixed-Tate-ness.** This combined ring is mixed-Tate
over $\mathbb{Q}$ iff $\mathbb{Q}[y_1^2, \ldots, y_n^2]$ is
mixed-Tate over $\mathbb{Q}$. The inner ring is a polynomial ring
in the Yukawa-squared values:
-   If $y_i^2 \in \overline{\mathbb{Q}}$ (algebraic), the inner
    ring is contained in $\overline{\mathbb{Q}}$, which is NOT
    inside the mixed-Tate period ring over $\mathbb{Q}$ (the period
    ring contains rationals and rational multiples of MZVs, not
    arbitrary algebraics).
-   If $y_i^2$ are transcendental but not periods of mixed Tate
    motives, the inner ring is outside.
-   The framework provides no mechanism to constrain $y_i$ to either
    case, and Sprint Yukawa-PSLQ confirms no low-coefficient
    pure-Tate identification at the available precision.

**Conclusion.** The AC-extension SD coefficients are mixed-Tate over
$\mathbb{Q}$ if and only if $y_i^2 \in \mathbb{Q}$. The framework
neither forces nor verifies this. **Generically the AC-extension SD
coefficients are NOT mixed-Tate over $\mathbb{Q}$, and the obstruction
sits entirely on the inner side, in the parameter-tied Dirichlet
ring of Paper 18 §IV.6 sixth tier ("inner-factor input data").**

The Sprint Mixed-Tate Test verdict (POSITIVE-with-refinement that
the OUTER M2 sector is pure-Tate $\bigoplus_k \pi^{2k} \cdot
\mathbb{Q}$) is preserved. The combined-triple verdict adds:\ the
inner factor escapes via Yukawa content, but the escape is
classified — it lives in the sixth tier of Paper 18 §IV.6
(parameter-tied Dirichlet, categorically disjoint from M1/M2/M3).

## Closure statement (proposed Paper 55 §7.4 update)

The current Paper 55 §7.4 text (lines 1215--1225):

> The almost-commutative extension $\mathcal{A}_{\mathrm{GV}}
> \otimes (\mathbb{C} \oplus \mathbb{H})$ studied in Paper 32
> §VIII.C (Sprint H1 closure) admits a Higgs structure but does
> not autonomously select the Yukawa coupling matrix. The
> associated inner-factor input data (Yukawa Dirichlet ring,
> Paper 18 §IV tier 6) is classified empirically as Dirichlet-$L$-class,
> but its period-theoretic placement inside the M1/M2/M3 trinity
> (or as a genuinely new sub-mechanism) is open.

Recommended replacement (sprint-scale 1-paragraph update; PI applies
after review):

> The almost-commutative extension $\mathcal{A}_{\mathrm{GV}} \otimes
> \mathcal{A}_F$ studied in Paper 32 §VIII.C (Sprint H1 closure)
> admits a Higgs structure but does not autonomously select the
> Yukawa coupling matrix. By Theorems
> `thm:eta_trivialization` and `thm:ac_factorization` of
> Paper 18 §IV.6, the combined spectral-action Seeley--Dewitt
> coefficients factorize as (outer pure-Tate $\bigoplus_k \pi^{2k}
> \cdot \mathbb{Q}$) $\times$ (inner Dirichlet ring
> $\mathbb{Q}[y_i^{-2s}]$). The inner factor is generically NOT
> mixed-Tate over $\mathbb{Q}$:\ it is mixed-Tate iff $y_i^2 \in
> \mathbb{Q}$, and Sprint Yukawa-PSLQ (2026-06-03,
> `debug/sprint_yukawa_pslq_memo.md`) returned zero hits on a
> 162-cell PSLQ sweep against the low-coefficient pure-Tate basis at
> both $\overline{\rm MS}\ M_Z$ and unification scales. The
> period-theoretic placement is therefore:\ the AC-extension SD
> coefficients split into a pure-Tate outer factor (M2 sub-mechanism
> on $S^3$) and a parameter-tied Dirichlet inner factor (Paper 18
> §IV.6 sixth tier), with no mixing. The "genuinely new
> sub-mechanism" reading of the open question is closed in the
> negative:\ no new sub-mechanism is required, only the
> categorically disjoint sixth tier of Paper 18, which is structurally
> compatible with — but not derived from — the master Mellin engine
> M1/M2/M3.

## Register marker (CLAUDE.md §1.7 WH2)

The current WH2 status reads:

> three of four axis-quadrants filled (Phases 4B-4I, Tiers 2-3,
> Sprint 4 QG/RH-J, Sprint 2 RH-M). [...] **Sprint Mixed-Tate Test
> (2026-06-03):** M2 sub-mechanism on $S^3$ inherits the
> Fathizadeh--Marcolli mixed-Tate classification (arXiv:1611.01815)
> at the standard volume-normalized convention, with a pure-Tate
> sharpening [...]. M3 sub-mechanism (cyclotomic mixed-Tate at
> level 4 over $\mathbb{Z}[i]$) remains an open mini-sprint for
> parallel classification.

Recommended addition (after the M3 mini-sprint line):

> **Sprint A4 scoping (2026-06-03, this session):** A4 (inner-factor
> mixed-Tate extension on the AC extension) closed by prior
> sprints. The combined-spectral-action SD coefficients factor as
> (outer pure-Tate) $\times$ (inner Dirichlet ring). The inner ring
> is generically NOT mixed-Tate over $\mathbb{Q}$ because Yukawa
> values are empirically Class 1 calibration data (Sprint
> Yukawa-PSLQ negative). The inner-factor content lives in Paper
> 18 §IV.6 tier 6 (parameter-tied Dirichlet), categorically
> disjoint from M1/M2/M3. No new sub-mechanism is required.

## Open follow-on (NOT A4 as scoped, but adjacent)

The mixed-Tate test on the **product of outer M2 with inner
Dirichlet** is structurally trivial because the rings are
multiplicative and disjoint. However, there's an adjacent
genuinely-open question:

-   **Joint outer-inner SD coefficients at higher heat-kernel order.**
    The M2 sub-mechanism at orders $k \ge 4$ (where the outer SD
    expansion terminates at $k = 1$ for Dirac on $S^3$ but doesn't
    terminate on $S^5$ Bargmann-Segal per Paper 51) produces $\pi^{2k}
    \cdot \mathbb{Q}$ outer coefficients. When multiplied by the
    inner $\sum_i m_i y_i^{2 - 2s}$, the result at integer $s$ is
    $\pi^{2k} \cdot \mathbb{Q}[y_i^2]$. The empirical test "is
    $\pi^{2k} \cdot$ inner trace in mixed-Tate?" is structurally
    identical to the Yukawa-PSLQ test at higher orders. No new
    content; sprint-scale 1--2 days if needed.

-   **Non-commutative-mixed-Tate extension (paper-grade NCG question).**
    Whether the AC-extension period theory admits a "noncommutative
    mixed-Tate" classification distinct from the commutative
    Brown--Deligne--Goncharov / Fathizadeh--Marcolli ones is a deep
    NCG question. Marcolli's program (Connes--Marcolli 2008) treats
    noncommutative motives for finite spectral triples, but the
    explicit period-ring classification analogous to the
    $\mathbb{Q}[\pi, \pi^{-1}]$ / MZV / cyclotomic-mixed-Tate
    hierarchy on the inner factor is not yet published to my
    knowledge. **Multi-year frontier**, not sprint-scale, and not the
    A4 question as scoped.

-   **Cyclotomic mixed-Tate test on M3 / vertex-parity sector
    (Sprint M3, separate from A4).** Paper 28 vertex-parity content
    ($D_{\mathrm{even}}(s) - D_{\mathrm{odd}}(s) = 2^{s-1}(\beta(s) -
    \beta(s-2))$) sits in mixed-Tate over $\mathbb{Z}[i]$ at level 4
    per the Deligne--Glanois descent of Paper 55 §6. A parallel
    mixed-Tate test on M3 would verify the cyclotomic placement
    directly. **Distinct from A4** (A4 is about the AC extension; M3
    is about the outer vertex-parity sector). Sprint Mixed-Tate Test
    open question 3 is the entry point; Paper 55 §6 outline is the
    placeholder.

## Prerequisite chain (for option (a) closure)

To apply the closure statement above, the following are required and
all already in place:

1.  **Theorem `thm:eta_trivialization` (Paper 18 §IV.6, in current
    repo).** $\checkmark$
2.  **Theorem `thm:ac_factorization` (Paper 18 §IV.6, in current
    repo).** $\checkmark$
3.  **Sprint Yukawa-PSLQ empirical sweep
    (`debug/sprint_yukawa_pslq_memo.md`, 2026-06-03).** $\checkmark$
4.  **Forced-Count Theorem `thm:forced_count` (Paper 32 §VIII,
    Sprint Forced-Count Synthesis 2026-06-03,
    `debug/sprint_forced_count_synthesis_memo.md`).** $\checkmark$
    (provides the formal forced/free seam at the $D_F$ moduli level)
5.  **Sprint Mixed-Tate Test memo
    (`debug/sprint_mixed_tate_test_memo.md`, 2026-06-03).**
    $\checkmark$ (provides the OUTER pure-Tate verdict, $\bigoplus_k
    \pi^{2k} \cdot \mathbb{Q}$)

All prerequisites in place. Closure is mechanical:\ Paper 55 §7.4
text update + register marker.

## Honest scope

-   **Theorem-grade content used here, not produced.** The
    $\eta$-trivialization and AC-factorization theorems were proven
    in Sprint H1 / Sprint TD Track 3 / Sprint Mixed-Tate Test
    pre-cursors and live in Paper 18 §IV.6. Sprint Yukawa-PSLQ
    is the empirical confirmation layer.
-   **What this memo produces:** the *naming* of A4 as
    CLOSED-BY-PRIOR plus the one-paragraph closure statement
    proposing the Paper 55 §7.4 update + register marker. The
    structural verdict is a CONSEQUENCE of the three prior
    theorems plus the empirical sweep; the consequence is what's new
    here.
-   **What's not closed:** the joint M2 $\times$ inner-Dirichlet
    higher-order mixed-Tate test (sprint-scale follow-on, named); the
    noncommutative mixed-Tate extension (paper-grade NCG, multi-year);
    Sprint M3 cyclotomic mixed-Tate (separate from A4, named).
-   **No paper edits applied in this memo.** PI applies the Paper 55
    §7.4 paragraph replacement + CLAUDE.md §1.7 WH2 register marker
    after review, per CLAUDE.md §13.5 / §13.8 paper-edit policy
    (sprint-scoping memos do not auto-edit; structural verdicts of
    "closed by prior" benefit from PI confirmation since A4 was
    explicitly carved out as a follow-on in v3.45.3).
-   **Precision-limited gap.** The Sprint Yukawa-PSLQ negative is at
    $M \le 1000$ coefficient ceiling and PDG precision (2--8 sig
    digits). High-coefficient identities ($M > 1000$) cannot be
    ruled out at the current precision and would require Yukawa
    measurements to 30+ digits. This gap is structural (external
    measurement limit), not a sprint failure. The framework's
    prediction does not specify a coefficient bound.

## Files used

### Memos (read)
-   `debug/sprint_yukawa_pslq_memo.md` (load-bearing empirical
    sweep, 2026-06-03)
-   `debug/sprint_mixed_tate_test_memo.md` (where A4 was originally
    scoped as open question 4; 2026-06-03)
-   `debug/sprint_forced_count_synthesis_memo.md` (Bertrand-analogue
    forced/free seam consolidation; 2026-06-03)
-   `debug/h1_ac_extension_memo.md` (Sprint H1 closure, the
    structural foundation; 2026-05-06)

### Papers (read)
-   `papers/group3_foundations/paper_18_exchange_constants.tex`
    (lines 1914--2089; §IV.6 inner-factor input data tier;
    `thm:eta_trivialization` lines 1927--1949;
    `thm:ac_factorization` lines 1951--1972; tier-6 statement lines
    1990--2000; sphaleron empirical falsification lines 2061--2089)
-   `papers/group1_operator_algebras/paper_32_spectral_triple.tex`
    (lines 2766--2965 §VIII.C Sprint H1 closure; lines 3890--3945
    `thm:forced_count` + `cor:bertrand_analogy` +
    `rem:forced_count_with_seam`)
-   `papers/group3_foundations/paper_55_periods_of_geovac.tex`
    (lines 1215--1225 §7.4 current text, target for update)

### No new files / no new computation
This is a structural scoping memo. No PSLQ runs, no new closed
forms, no code changes. The verdict crystallizes the consequence of
prior theorems.

## Memo location

`debug/sprint_a4_scoping_memo.md` (this file).
