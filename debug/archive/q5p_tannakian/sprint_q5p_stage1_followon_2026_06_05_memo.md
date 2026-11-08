# Sprint Q5'-Stage1-Followon — joint closure of Stage 1's two open follow-ons

Date: 2026-06-05 (close of day after the morning Q5'-CH arc v3.57.0 and afternoon Q5'-Stage1-Arc v3.58.0)
Scope: two parallel sub-sprints (Track 1 + Track 2) executed in background while PI conversation proceeded. Paper 32 §VIII edits applied sequentially after both returned to avoid race; Paper 55 (Track 1) and Paper 18 (Track 2) edits applied in-track.

Sub-sprint memos (detail):
- `debug/sprint_q5p_cm_bicomplex_memo.md` (Track 1: CM-residue bicomplex for M3)
- `debug/sprint_q5p_continuum_mellin_memo.md` (Track 2: continuum-limit Mellin analysis)

## 1. TL;DR

The Q5'-Stage1-Arc (v3.58.0) explicitly named two open follow-ons that were sprint-scale rather than multi-year: (a) the CM-residue bicomplex housing M3 outside JLO on commutative algebra; (b) the continuum-limit Mellin analysis that lifts Sub-Sprint 2b's "M1 lives at $s = d/2$, not $s = 0$" from structural sketch to theorem grade. Both closed POSITIVE in one afternoon as parallel sub-sprints.

The **joint structural finding** (visible only because both tracks landed): the cohomological-dual side of Stage 1 is closed by two cocycle classes on the same $K_0$ sector decomposition $\{[e_s]\}_{s=0..4}$, with JLO $\mathrm{HP}^{\mathrm{even}}$ class $(+2, -2, +2, +2, -4)$ encoding *chirality balance* (zeroth heat-kernel supertrace, McKean–Singer index $\chi(S^3) = 0$) and CM-$\eta$ class $(3, 3, 5, 15, 10)$ encoding *chirality-symmetrized eigenvalue magnitude* (leading $\eta$-density supertrace, M3 mechanism leading). The M1 residue at $s = d/2 = 3/2$ is bit-exactly $\sqrt{\pi}/2$ by two independent routes (Hurwitz pole + Weyl-law Gilkey heat-kernel), and the three Hopf-base normalizations across Papers 18, 38, and this sprint ($\pi$, $4/\pi$, $\sqrt{\pi}/2$) are exact-factor consistent — $\pi^2/4$ connects Paper 18 to Paper 38, $\pi^{3/2}/8$ connects Paper 38 to this sprint. All bit-exact at finite cutoff; Tauberian rate uniformity in a neighborhood of $s = 3/2$ is named with Karamata 1962 / Korevaar 2004 published precedent as the only Stage-2-relevant gap.

## 2. Joint structural finding

The two tracks were chosen as the cohomological-side and analytical-side closures of Stage 1, and the joint result is more than the sum. At the level of the same $K_0$ sector decomposition $\{[e_s]\}_{s=0..4}$ on the truncated Camporesi–Higuchi spectral triple at $n_{\max} = 2$:

| | JLO $\mathrm{HP}^{\mathrm{even}}$ | CM-$\eta$ class |
|:---|:---|:---|
| Class | $(+2, -2, +2, +2, -4) \in \mathbb{Q}^5$ | $(3, 3, 5, 15, 10) \in \mathbb{Z}^5$ |
| Sum | $0 = \chi(S^3)$ (McKean–Singer) | $36 = M_3(n_{\max}=2)$ (Sub-Sprint 1 g.t.) |
| Sector formula | $\chi_s = \dim(\gamma_+ e_s \mathcal{H}) - \dim(\gamma_- e_s \mathcal{H})$ | $\eta_s = \dim_s \cdot (n_s + 1/2)$ |
| Reading | chirality balance | chirality-symmetrized eigenvalue magnitude |
| Heat-kernel supertrace | zeroth coefficient | leading $t^0$ $\eta$-density |
| Master Mellin host | M1 + M2 via different $t$-orders of $\phi_0^{\mathrm{odd}}$ | M3 alone via $\eta$-pairing |

These are not redundant labels of the same K-theory data — they are categorically distinct cohomological invariants on the same $K_0$ decomposition, dual via the chirality symmetry ($\gamma_+ - \gamma_- \to \chi$; $\gamma_+ + \gamma_- \to \eta$). The Sub-Sprint 1 "JLO+CM mixed framework" reading is upgraded from a *representation-theoretic split* (one functor with three slots) to a *class-level dual decomposition* (two cocycle classes coordinated by the $K_0$ generators). Stage 2 (motivic Galois action) therefore acts on this DUAL, not on a single class with three slots — sharpening the v3.58.0 "two classes coordinating three analytic extractions" statement.

The Mellin-analytic side gives the matching synthesis. The three Hopf-base normalizations of the M1 mechanism — $\mathrm{Vol}(S^2)/4 = \pi$ (Paper 18 §III.2 Hopf-base measure), $\mathrm{Vol}(S^2)/\pi^2 = 4/\pi$ (Paper 38 §VIII L2 asymptote rate constant), $\Gamma(d/2) = \sqrt{\pi}/2$ (this sprint, Mellin residue at the spectral-dimension pole) — are NOT three independent appearances but exact-factor-consistent sibling normalizations of the *same* spectral object. The factors are pure powers of $\pi$ (no spectral-data input). The published $4/\pi$ in Paper 38 is the M1 signature at the rate-class level of GH convergence; the $\sqrt{\pi}/2$ here is the M1 signature at the Mellin-residue level; the $\pi$ in Paper 18 is the M1 signature at the volume-form level. Three faces, one object.

## 3. Sub-sprint summaries

### Track 1 — CM-residue bicomplex for M3

The CM-$\eta$ pairing residue $(b, B)$ bicomplex on the truncated CH triple at $n_{\max} = 2$ constructed bit-exactly in sympy.Rational. The CM-$\eta$ cochain $\psi_n^{\eta}$ replaces the JLO leftmost insertion $\gamma$ by $\gamma D$, sharing the $(b+B)\psi = 0$ bicomplex structure with JLO via the same $b/B$ formulas. Degree-1 cocycle condition bit-exactly zero on the full panel of 36 idempotent-pair inputs at three $t$-orders both flavours (216 bit-exact zero residuals); no truncation artifact of the JLO $\pm 1/196608$ type. Explicit CM-$\eta$ class on the Morita-trivial baseline $\mathbb{Q}^5 = \mathrm{HP}_0(\mathbb{C}^5)$: $(3, 3, 5, 15, 10) \in \mathbb{Z}^5$, summing to $36 = M_3(n_{\max} = 2)$. Cross-check at $n_{\max} = 3$: $M_3(3) = 120$ bit-exact, cocycle bit-exact zero on all 9 diagonal sectors.

Paper edits: Paper 55 §subsec:open_m2_m3 extended by one paragraph (23 → 24 pages, three-pass clean); Paper 32 §VIII new remark `rem:q5p_cm_residue_class` applied after `rem:q5p_hp_even_class`.

Driver: `debug/compute_cm_residue_bicomplex.py`. Data: `debug/data/sprint_q5p_cm_bicomplex.json`. Memo: `debug/sprint_q5p_cm_bicomplex_memo.md`.

### Track 2 — continuum-limit Mellin analysis

Continuum-limit Tauberian step from Sub-Sprint 2b's "structurally sketched" finding closed at Paper 38 L2 grade. The continuum spectral zeta $\zeta_{D^2}^{\mathrm{cont}}(s) = 2\,\zeta(2s-2, 3/2) - (1/2)\,\zeta(2s, 3/2)$ has a simple pole at $s = d/2 = 3/2$ with bit-exact meromorphic residue $1$ (verified two independent ways: Hurwitz pole at $u = 2s - 2 = 1$, and direct Laurent expansion via digamma identity). The Mellin residue $\Gamma(3/2) \cdot 1 = \sqrt{\pi}/2$ matches the standard Gilkey 1995 / Vassilevich 2003 heat-kernel Weyl-law leading coefficient $2 \cdot \mathrm{Vol}(S^3)/(4\pi)^{3/2} = \sqrt{\pi}/2$ bit-exactly. M2 panel at integer $s \in \{1, 2, 3, 4, 5\}$ reproduces the Q5'-CH-2 Seeley–DeWitt panel verbatim: $-\pi^2/4$, $\pi^2 - \pi^4/12$, $\pi^4/3 - \pi^6/30$, $2\pi^6/15 - 17\pi^8/1260$, $17\pi^8/315 - 31\pi^{10}/5670$ (5/5 bit-exact). Three-sibling normalization ($\pi \leftrightarrow 4/\pi \leftrightarrow \sqrt{\pi}/2$) internally consistent via exact factors. M3 $\eta$-pairing values $(6, 36, 120, 300)$ bit-exact at $n_{\max} \in \{1, 2, 3, 4\}$.

Paper edits: Paper 18 §III.7 sharpened in place (26 pages, unchanged, three-pass clean); Paper 32 §VIII new remark `rem:q5p_continuum_residue` applied after `rem:q5p_mellin_extraction_points`; four new bibitems added to Paper 32 (gilkey1995, vassilevich2003, karamata1962, korevaar2004).

Driver: `debug/compute_continuum_mellin_residue.py`. Data: `debug/data/sprint_q5p_continuum_mellin.json`. Memo: `debug/sprint_q5p_continuum_mellin_memo.md`.

## 4. Paper edits applied

- **Paper 55 §subsec:open_m2_m3** (Track 1): +1 paragraph. Pages 23 → 24, three-pass clean.
- **Paper 18 §III.7** (Track 2): M1 paragraph sharpened in place with explicit $s = d/2$ pole interpretation and $\sqrt{\pi}/2$ closed-form. Pages 26 (unchanged), three-pass clean.
- **Paper 32 §VIII**: two new remarks applied sequentially (rem:q5p_cm_residue_class after rem:q5p_hp_even_class from Track 1; rem:q5p_continuum_residue after rem:q5p_mellin_extraction_points from Track 2); four new bibitems (gilkey1995, vassilevich2003, karamata1962, korevaar2004 from Track 2). Pages 66 (unchanged at page-count resolution), three-pass clean. One sprint-introduced LaTeX warning during first compile (cross-paper `\ref` to a Paper 55 label) corrected during sequential application — convention is literal text for cross-paper labels per the existing §VIII style.

## 5. Verification

No production code touched (`geovac/` clean). All activity in `debug/` (new drivers + data) + three papers. Per CLAUDE.md §9 / `/regression` skill, the diff-derived consumer-test selection is therefore empty; topological-integrity baseline (18 S³ symbolic proofs) unaffected by paper edits or new debug scripts. No new equations needing tests per §13.4a — the bit-exact JLO/CM-residue identifications and the Mellin residue identifications run against existing infrastructure with verified ground truth (Sub-Sprint 1 panel for Track 1; Q5'-CH-2 panel + Gilkey/Vassilevich heat-kernel for Track 2).

## 6. Honest scope

**Closed at theorem grade (bit-exact at finite cutoff):**

- Track 1: CM-$\eta$ $(b+B)\psi = 0$ bit-exact on full 216-cell panel at $n_{\max}=2$ (both flavours, three $t$-orders, 36 idempotent pairs).
- Track 1: CM-$\eta$ class identification $(3, 3, 5, 15, 10) \in \mathbb{Z}^5$ on the Morita-trivial baseline $\mathbb{Q}^5$.
- Track 1: Cross-cutoff sanity at $n_{\max}=3$: $M_3(3) = 120$ bit-exact, cocycle bit-exact zero on all 9 diagonal sectors.
- Track 2: M1 Mellin residue at $s = d/2 = 3/2$ is bit-exact $\sqrt{\pi}/2$ via two independent routes (Hurwitz pole + Weyl-law Gilkey heat-kernel).
- Track 2: M2 continuum panel at integer $s \in \{1, 2, 3, 4, 5\}$ reproduces Q5'-CH-2 5/5 bit-exact.
- Track 2: Three-sibling Hopf-base normalization ($\pi \leftrightarrow 4/\pi \leftrightarrow \sqrt{\pi}/2$) internally consistent via exact factors $\pi^2/4$ and $\pi^{3/2}/8$ (no spectral-data input).

**Structural sketch (not yet theorem):**

- Tauberian rate uniformity in a neighborhood of $s = 3/2$: named open question with Karamata 1962 §V.3 and Korevaar 2004 Ch. III §4 published precedent; reachable but not closed at sprint scale.

**Numerical observation / interpretive synthesis (no new bit-exact theorem):**

- The dual reading "JLO chirality balance ↔ CM-$\eta$ chirality-symmetrized magnitude as opposite-sign halves of the same $K_0$ decomposition" is the structural insight that running the two tracks in parallel makes visible. The individual classes are bit-exact theorem-grade; the dual interpretation lives at the cohomological-reading level rather than as a formal cohomology identity (where $\gamma_+ + \gamma_- = 1$ would make the duality formal).

**Named open follow-ons:**

- Continuum-limit analysis for M3 — analog of Track 2's M1/M2 closure for the CM-$\eta$ side, proving the M3 continuum residue against literature precedent at quarter-integer Hurwitz shifts. Sprint-scale.
- Pro-system functoriality — relate the cohomology classes at different cutoffs via the truncation / Berezin maps (named in v3.58.0 as "no published precedent"). Sprint-scale but harder, and depends on either Track 1's bicomplex (already constructed) or Track 2's continuum machinery (already established).
- Stage 2 (multi-year): the strict Tannakian construction of $\omega^{\mathrm{tri}}$ as enriched fiber functor and the identification of the motivic Galois group acting on the dual cocycle classes.

**Hard prohibitions check (§13.5):** No changes to natural geometry hierarchy. No fitted/empirical parameters introduced. No deletion of negative results from §3. No removal of "conjectural" label from the Paper 2 combination rule $K = \pi(B + F - \Delta)$.

**Curve-fit audit (`feedback_audit_numerical_claims`):** Both tracks' bit-exact closures use re-derivation, not curve-fit. Track 1's CM-$\eta$ class is structurally derived ($\eta_s = \dim_s \cdot (n_s + 1/2)$), not fit from the 5-cell vector. Track 2's M1 residue is two-independent-route derivation (Hurwitz pole + Weyl-law Gilkey), not curve-fit to a known coefficient. Track 2's three-sibling normalization identification is purely algebraic exact-factor matching ($\pi^2/4$ and $\pi^{3/2}/8$), no PSLQ.

**Discrete-for-skeleton compliance (`feedback_discrete_for_skeleton`):** Both tracks bit-exact sympy.Rational throughout. Track 2's transcendental content ($\pi$) is tagged: M1 = Hopf-base measure (Paper 18 §III.2), explicitly entering at the continuum boundary via Mellin residue. Track 1 stays purely on the skeleton side (cocycle class lives in $\mathbb{Z}^5$ on Morita-trivial baseline; no transcendentals introduced).

**Tag transcendentals (`feedback_tag_transcendentals`):** The $\pi$ in Track 2's $\sqrt{\pi}/2$ residue is M1 by classification (Paper 18 §III.7), Hopf-base measure mechanism, projection-chain identity at the spectral-dimension pole $s = d/2$. The $\pi^{2k}$ powers in the M2 panel are M2 by classification (Paper 18 §III.7), Seeley–DeWitt mechanism at integer-$s$ regular points. The $4/\pi$ in Paper 38 is M1 at the GH-rate-class level. All consistent with the case-exhaustion theorem (Paper 32 §VIII `thm:pi_source_case_exhaustion`).

## 7. Files

### Memos
- `debug/sprint_q5p_stage1_followon_2026_06_05_memo.md` (this umbrella)
- `debug/sprint_q5p_cm_bicomplex_memo.md` (Track 1 detail)
- `debug/sprint_q5p_continuum_mellin_memo.md` (Track 2 detail)

### Drivers (new this sprint)
- `debug/compute_cm_residue_bicomplex.py`
- `debug/compute_continuum_mellin_residue.py`

### Data (new this sprint)
- `debug/data/sprint_q5p_cm_bicomplex.json`
- `debug/data/sprint_q5p_continuum_mellin.json`

### Paper files modified
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` (+2 remarks, +4 bibitems)
- `papers/group3_foundations/paper_18_exchange_constants.tex` (M1 paragraph sharpened in place)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` (+1 paragraph)

## 8. One-line verdict

Two parallel sub-sprints closed Stage 1's two named sprint-scale open follow-ons in one afternoon arc: the cohomological-dual side now carries two cocycle classes on the same $K_0$ sector decomposition (JLO $\mathrm{HP}^{\mathrm{even}}$ chirality balance + CM-$\eta$ chirality-symmetrized magnitude), and the analytical side now carries the M1 mechanism's bit-exact $\sqrt{\pi}/2$ Mellin residue at $s = d/2$ derived by two independent routes with three-sibling $\pi$-normalization internally consistent. Stage 1 is now bit-exactly constructed on five diagnostic axes (the original four from v3.58.0 plus the dual cocycle-class structure); Stage 2 remains multi-year but better-scoped, acting on the dual class structure rather than on a single multi-slot functor.
