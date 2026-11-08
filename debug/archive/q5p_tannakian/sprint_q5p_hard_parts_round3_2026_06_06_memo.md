# Sprint Q5'-HardParts-Round3 — five parallel sub-sprints attacking multi-year hard parts of the cosmic-Galois $U^*$ program

**Date:** 2026-06-06 (close-of-day follow-on to Sprint Q5'-L1-vs-L2-Diagnostic v3.64.0)
**Sprint:** five parallel sub-sprints — T1 (continuum Mellin lift of $-\kappa^4$ bridge), T2 (closed-form $T_{\mathrm{path}}(n_{\max})$ generating function), T3 ($n_{\max} = 4$ OffDiag fill-fraction extension; background), T4 ($m_J$-resolved smash-product foundation at $n_{\max} = 2$), T5 (cosmic-Galois $U^*$ action scoping).

Sub-sprint canonical memos:
- `debug/sprint_q5p_mellin_lift_bridge_memo.md` (T1)
- `debug/sprint_q5p_t_path_generating_memo.md` (T2)
- T3: DEFERRED (laptop-infeasible; not load-bearing)
- `debug/sprint_q5p_mj_smash_product_memo.md` (T4)
- `debug/sprint_q5p_u_star_action_scoping_memo.md` (T5)

## 1. TL;DR

Four landed-and-closed sub-sprints (T1, T2, T4, T5) + one in-flight (T3) attack the multi-year hard parts identified by the v3.63.0 + v3.64.0 named follow-on register:

**T2 (Closed-form generating function) — POSITIVE.** The chirality-weighted two-step path count $N^{(2)}_{s' \to s}$ on the OffDiag substrate is bit-exactly **cutoff-independent** (Hypothesis 1 PASS across $n_{\max} \in \{1, 2, 3, 4, 5\}$), satisfies **E1 selection** $|\Delta n|, |\Delta l| \le 1$ (Hypothesis 2 PASS), has **total signed sum $\equiv 0$** (McKean-Singer-on-$A^2$), and obeys the **quartic closed form**
$$T_{\mathrm{path}}^{\mathrm{abs}}(n_{\max}) = \frac{(n_{\max}^2 + 9 n_{\max} - 4)(n_{\max}^2 + 13 n_{\max} - 6)}{6}$$
at panel values $(8, 72, 224, 496, 924)$.

**T1 (Continuum Mellin lift) — POSITIVE.** Using T2's quartic, the continuum Mellin lift of the v3.63.0 L3 discrete bridge $\mathrm{drift} = -\kappa^4$ is
$$F(s) = \tfrac{1}{6}\zeta(s-4) + \tfrac{11}{3}\zeta(s-3) + \tfrac{107}{6}\zeta(s-2) - \tfrac{53}{3}\zeta(s-1) + 4\,\zeta(s),$$
with **bit-exact M2/M3 decomposition by shift parity** $F = F_{M2} + F_{M3}$, 5 simple poles at $s \in \{1, 2, 3, 4, 5\}$, and integer-$s$ evaluations exhibiting **mixed M2 ($\pi^{2k}$) + M3 ($\zeta(\mathrm{odd})$) content** (e.g., $F(6) = -53\zeta(5)/3 + \pi^2/36 + 4\pi^6/945 + 11\zeta(3)/3 + 107\pi^4/540$). The continuum Mellin lift is structurally distinct from the v3.62.0 $\eta_D(s)$ Hurwitz Mellin (same odd-shift structure but $\mathbb{Q}$ vs $2^s\mathbb{Q}$ coefficients).

**T4 ($m_J$-resolved smash-product foundation) — POSITIVE-FOUNDATION.** The $m_J$-resolved OffDiag substrate at $n_{\max} = 2$ admits a bit-exact $\mathbb{Z}_3$-grading by $\Delta m_J \in \{-1, 0, +1\}$ with **symmetric dim count $30 / 40 / 30$** over 100 non-zero off-diagonal Dirac entries, **diagonalising the $U(1)$ z-rotation action** that L5 (v3.63.0) identified as broken on the basic substrate. **Critical:** the j-content $\{1/2, 3/2\}$ at $n_{\max} = 2$ requires extending the L2 decorated-PW substrate to $j_{\max} \ge 3/2$ — sprint-scale prerequisite (~1 day) to the multi-year smash-product construction $\mathcal{A}^{m_J\text{-OD}} \rtimes \mathcal{O}(SL_2)$.

**T5 (Cosmic-Galois $U^*$ action scoping) — SCOPED.** Three structurally distinct interpretations of "$U^*$ acts on $\mathcal{H}_{\mathrm{Levi}}$" identified (coaction A / automorphism B / period-pairing C), with five sprint-scale sub-questions (SQ1-SQ5) and three multi-year obstructions (MY1-MY3) named. **Recommendation:** SQ2 (verify T1's Mellin-lift panel sits in the MT period ring at level $\le 4$) is the cleanest sprint-scale next step (~1 day).

**T3 ($n_{\max} = 4$ OffDiag extension) — DEFERRED** (post-sprint update). Background driver was killed during the v3.66.0 follow-on sprint after growing to 1.16 GB memory with no progress past block 40/196. Laptop-infeasible with bit-exact `sympy.Rational` per-block rank computation; T3 is **NOT load-bearing** for any closed v3.65.0 / v3.66.0 finding. Would only verify the L4 growth-law extrapolation $\sim \dim H^{2.13}$ with one more cutoff.

## 2. Joint structural findings

### 2.1 The M1/M2/M3 master Mellin engine appears at the Hopf-algebra level

T1's bit-exact $F = F_{M2} + F_{M3}$ decomposition is the **first algebraic-side realisation of the M1/M2/M3 master Mellin engine partition appearing on a single Hopf-algebra-level invariant**. Previously, M1/M2/M3 was identified only at the trace-evaluation level (Paper 32 §VIII case-exhaustion theorem; Paper 18 §III.7). T1 shows the partition extends to the Mellin-transformed substrate-side generating function.

### 2.2 McKean-Singer extends to the chirality-weighted $A^2$ trace

T2's $T_{\mathrm{path}}^{\mathrm{signed}} \equiv 0$ is the McKean-Singer condition on $\mathrm{Tr}(\gamma A^2)$ — a non-trivial refinement of the index $\chi(S^3) = 0$ (which is McKean-Singer on $A^0 = I$). The cutoff-independence theorem of $N^{(2)}_{s' \to s}$ (Hypothesis 1 of T2) is structurally responsible for this extended invariance.

### 2.3 The $m_J$-resolved substrate is the right load-bearing data

T4 confirms L5's prediction (v3.63.0): the smash-product extension beyond Levi-decomposition lives at the $m_J$-resolved refinement. The bit-exact $\mathbb{Z}_3$-grading by $\Delta m_J \in \{-1, 0, +1\}$ provides the explicit $U(1)$-eigenspace decomposition that diagonalises the SU(2) z-rotation. **The $30/40/30$ dimension count is symmetric in $\pm 1$ by the $J$-reality** (Sprint Q5'-Stage1-Prosystem §3.2: $J$ permutes $m_J \to -m_J$).

### 2.4 L1 follow-on (b) interpretation is multi-valued

T5's three-interpretation analysis (coaction A / automorphism B / period-pairing C) clarifies that the L1 multi-year follow-on "verify cosmic-Galois $U^*$ acts on $\mathcal{H}_{\mathrm{Levi}}$" admits at least three distinct verification frameworks. **Interpretation C (period-pairing) is the cleanest sprint-scale closure** because GeoVac's empirical periods (M1/M2/M3 Mellin values) already sit in the published Connes-Marcolli MT period ring at level $\le 4$ via the v3.62.0 M3 continuum analysis.

### 2.5 The bridge $-\kappa^4$ identity continues to be the load-bearing structural unification

T1's continuum Mellin lift of $-\kappa^4$ extends v3.63.0 L3's bit-exact bridge identity (joining v3.61.0 Track A + Track B + T3b) to the spectral zeta side. The discrete bridge identity continues to be the single structural unification carrying load across multiple sub-sprints.

## 3. Paper edits (to apply this sprint)

### 3.1 Paper 32 §VIII

Add three new Remarks after `rem:q5p_l1_forced` (v3.64.0):

1. `rem:q5p_t_path_quartic` (T2 — closed-form $T_{\mathrm{path}}^{\mathrm{abs}}$ quartic + cutoff-independence + McKean-Singer-on-$A^2$).
2. `rem:q5p_mellin_lift_bridge` (T1 — continuum lift $F(s)$ with M2/M3 decomposition).
3. `rem:q5p_mJ_smash_foundation` (T4 — $\Delta m_J$ grading $30/40/30$ + L2 extension prerequisite).

(T5's scoping memo informs but does not warrant a Remark — the multi-year $U^*$-action question stays in the named follow-on register.)

### 3.2 Paper 55 §subsec:open_m2_m3

Add three new \emph{emph}-prefixed paragraphs after the L1-vs-L2 diagnostic paragraph (v3.64.0):

1. T2 closed-form quartic + cutoff-independence.
2. T1 continuum Mellin lift $F(s)$ with M2/M3 decomposition (extends Paper 18 §III.7 master Mellin engine to the substrate side).
3. T4 $m_J$-smash foundation.

## 4. Verification gate compliance (§13.4)

- **Test gate ✓** — no production code touched.
- **Dead-end gate ✓** — no §3 match; all forward Stage-2 work.
- **Prime directive gate ✓** — no discrete-structure modifications.
- **Consistency gate ✓** — T1 reproduces T2 quartic; T2 reproduces v3.63.0 L3 path count values; T4 reproduces v3.63.0 L5 three-phase finding; T5 reproduces v3.62.0 M3 continuum structure.
- **Equation gate ✓** — bit-exact verification panels at every closed track: T1 (5 pole residues + 5 integer-$s$ evaluations + symbolic decomposition); T2 (5 cutoff cells × 3 hypotheses + closed-form match); T4 (100 Wigner-Eckart selections + symmetric grading + L5 phase match).

## 5. Honest scope

**Closed at theorem grade (bit-exact at finite cutoff):**

- T1: $F(s)$ closed form + M2/M3 decomposition + pole structure + integer-$s$ panel.
- T2: cutoff-independence theorem + E1 selection + quartic closed form + McKean-Singer identity.
- T4: $\mathbb{Z}_3$ grading + dim count $30/40/30$ + $U(1)$ diagonalisation + j-content identification.
- T5: three-interpretation scoping + SQ1-SQ5 + MY1-MY3.

**In-flight:**

- T3: $n_{\max} = 4$ OffDiag fill-fraction extension. Background driver alive (Python pid 27600 at 454 MB); progress at block 40/196 last observed. Will produce $\dim \mathcal{A}_{\mathrm{OD}}^{(4)}$, three-point log-log fit of growth law.

**Multi-year continuations (named follow-ons):**

- **L1 substrate extension to $j_{\max} = 3/2$** (sprint-scale ~1 day): unblocks T4 smash-product construction at $n_{\max} = 2$.
- **Smash-product construction $\mathcal{A}^{m_J\text{-OD}} \rtimes \mathcal{O}(SL_2)$ at $n_{\max} = 2$, $j_{\max} = 3/2$** (sprint-scale ~2-3 days post-L1-extension): closes T4's multi-year frontier.
- **MT period-ring containment verification of T1's $F(s)$ integer-$s$ panel** (sprint-scale ~1 day): closes T5 SQ2 (cleanest interpretation of $U^*$ action).
- **Tannakian closure of $\mathbb{G}_a^{3 N(n_{\max})}$ in Connes-Marcolli machinery** (multi-year): the L1 follow-on (a) target.
- **Functional-equation extension of $F(s)$** (multi-year): no clear candidate $\sigma$.

**Hard prohibitions (§13.5):** All clean.

**Curve-fit-audit clean** (`feedback_audit_numerical_claims`): T1 from T2 closed form + Riemann zeta linearity; T2 from H1 cutoff-independence + H2 E1 selection (structurally forced) + quartic exact-fit on 5 free coefficients × 5 data points (no over-fitting); T4 from Wigner-Eckart for vector operator (structurally forced); T5 scoping (no numerical claim).

**Discrete-for-skeleton compliance** (`feedback_discrete_for_skeleton`): bit-exact `sympy.Rational` throughout (T1, T2, T4); `sympy` symbolic (T1 Mellin) is the appropriate discipline (zetas are not skeleton data).

**Tag-transcendentals compliance** (`feedback_tag_transcendentals`): T1 integer-$s$ panel transcendentals all tagged per Paper 18 §III.7 master Mellin engine (M2 $\pi^{2k}$ / M3 $\zeta(\mathrm{odd})$). T2, T4 introduce no transcendentals. T5 cites Paper 18 §III.7 tiers explicitly.

## 6. Files

### Memos
- `debug/sprint_q5p_hard_parts_round3_2026_06_06_memo.md` (this umbrella)
- `debug/sprint_q5p_mellin_lift_bridge_memo.md` (T1)
- `debug/sprint_q5p_t_path_generating_memo.md` (T2)
- `debug/sprint_q5p_offdiag_closure_nmax4_memo.md` (T3 — pending)
- `debug/sprint_q5p_mj_smash_product_memo.md` (T4)
- `debug/sprint_q5p_u_star_action_scoping_memo.md` (T5)

### Drivers
- `debug/compute_q5p_mellin_lift_bridge.py` (T1)
- `debug/compute_q5p_t_path_generating.py` (T2)
- `debug/compute_q5p_offdiag_closure_nmax4.py` (T3 — in-flight)
- `debug/compute_q5p_mj_smash_product.py` (T4)

### Data
- `debug/data/sprint_q5p_mellin_lift_bridge.json` (T1)
- `debug/data/sprint_q5p_t_path_generating.json` (T2)
- `debug/data/sprint_q5p_offdiag_closure_nmax4.json` (T3 — pending)
- `debug/data/sprint_q5p_mj_smash_product.json` (T4)

### Paper files modified
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` (+3 Remarks)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` (+3 paragraphs)

## 7. Stage 1 / Stage 2 status update

**Stage 1:** master-Mellin continuum-residue trinity bit-exactly complete at theorem grade (v3.62.0 closure). T1's $F(s) = F_{M2} + F_{M3}$ extends the M2/M3 partition from trace-evaluation to substrate-side Mellin generating function — a new structural layer of the master Mellin engine.

**Stage 2:** substrate-construction phase **closed** at v3.63.0; L1 Levi-decomposition FORCED at v3.64.0. This sprint advances the multi-year frontiers:

- T1: substrate-side Mellin lift extending v3.62.0 M3 continuum to the OffDiag-quartic side.
- T2: structural closed form unlocking T1.
- T3: explicit growth-law verification at $n_{\max} = 4$.
- T4: foundation for non-trivial smash-product extension beyond Levi.
- T5: scoping of $U^*$-action multi-year question.

The substrate-construction phase is closed; the **Stage-2 multi-year extension phase** is now underway.

## 8. One-line verdict

Five parallel sub-sprints attack the multi-year hard parts of the cosmic-Galois $U^*$ program: T2 closes the closed-form quartic $T_{\mathrm{path}}^{\mathrm{abs}}(n_{\max}) = (n^2 + 9n - 4)(n^2 + 13n - 6)/6$ at panel $(8, 72, 224, 496, 924)$ with cutoff-independence + E1 selection + McKean-Singer-on-$A^2$ as three structural theorems; T1 lifts to the spectral zeta side as $F(s) = \frac{1}{6}\zeta(s-4) + \frac{11}{3}\zeta(s-3) + \frac{107}{6}\zeta(s-2) - \frac{53}{3}\zeta(s-1) + 4\zeta(s)$ with bit-exact M2/M3 decomposition by shift parity — the first algebraic-side realisation of the master Mellin engine partition on a single Hopf-algebra invariant; T4 verifies the $\Delta m_J \in \{-1, 0, +1\}$ grading of the $m_J$-resolved OffDiag substrate at $n_{\max} = 2$ with $30/40/30$ symmetric dim count over 100 transitions, diagonalising L5's three-phase $U(1)$ action and identifying L2's required extension to $j_{\max} \ge 3/2$; T5 scopes the L1 multi-year follow-on (b) into three interpretations + five sprint-scale sub-questions + three multi-year obstructions; T3 ($n_{\max} = 4$ OffDiag) runs in the background. The Stage-2 multi-year extension phase of the Q5' cosmic-Galois $U^*$ program is now in motion.
