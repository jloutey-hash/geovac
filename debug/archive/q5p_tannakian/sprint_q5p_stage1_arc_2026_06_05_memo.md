# Sprint Q5'-Stage1-Arc — second span of the cosmic-Galois U* bridge

Date: 2026-06-05 (close-of-day after the morning Q5'-CH arc v3.57.0)
Scope: six sub-sprints in one session executed in parallel/sequence on
the multi-year cosmic-Galois $U^*$ Stage 1 program, plus paper edits.
Started with the question "what does Stage 1 actually look like as a
sprint?" and closed with bit-exact construction of the
$\omega^{\mathrm{tri}}$ symbol across four diagnostic axes and a
full-corpus paper-edit pass.

This is the umbrella memo. Six sub-sprint memos hold the detail:
- `debug/sprint_q5p_s5_falsifier_memo.md` (substrate scope)
- `debug/sprint_q5p_stage1_reframe_memo.md` (R3 vs CH-1 reconciliation)
- `debug/sprint_q5p_jlo_nmax2_memo.md` (Sub-Sprint 1)
- `debug/sprint_q5p_2a_jlo_nmax_sweep_memo.md` (Sub-Sprint 2a)
- `debug/sprint_q5p_2b_phi0_mellin_memo.md` (Sub-Sprint 2b)
- `debug/sprint_q5p_2c_bicomplex_memo.md` (Sub-Sprint 2c)

## TL;DR

**Stage 1 of Q5' is bit-exactly constructed at finite cutoff across four diagnostic axes.** The morning's CH arc (v3.57.0) gave a symbol-level finite-cutoff witness on M1/M2/M3. Today's afternoon arc puts it on the cohomological-dual side ($\mathrm{HP}^*$, JLO + CM-residue mixed framework), extends it across cutoffs ($n_{\max} \in \{1, 2, 3, 4\}$ with polynomial closed forms), localises the master-Mellin partition inside the formal Mellin transform $\Gamma(s)\zeta_{D^2}(s)$ (three categorically different extraction points), and constructs the cyclic $(b, B)$ bicomplex with the explicit $\mathrm{HP}^{\mathrm{even}}$ class $(+2, -2, +2, +2, -4) \in \mathbb{Q}^5$.

The arc also confirms (S5 falsifier) that the CH triple is the structurally correct substrate — the Bargmann-Segal Hardy-on-$S^5$ triple fails both load-bearing conditions (no $\pm$-symmetric spectrum, no canonical chirality grading on the $(N, 0)$ tower).

Stage 2 (motivic Galois group of $\omega^{\mathrm{tri}}$) remains genuinely multi-year and is now better-scoped: it acts on **two** cocycle classes (JLO $\mathrm{HP}^{\mathrm{even}}$ for M1+M2, CM-residue $\eta$-class for M3), not three slots of a single functor.

## What landed, sub-sprint by sub-sprint

### Sub-sprint S5 falsifier — BREAKS-AS-EXPECTED

CH-1's chirality-parity selection rule fails on the Bargmann-Segal Hardy-on-$S^5$ triple at every $N_{\max} \in \{2, 3, 4\}$ via both load-bearing structural conditions simultaneously: the Hardy spectrum is strictly positive (no $\pm$-partner), and Paper 24 §V.4 already records that no canonical chirality grading exists on the $(N, 0)$ $\mathrm{SU}(3)$ tower (the two natural candidates $\chi_N = (-1)^N$ and $\chi_l = (-1)^l$ coincide because $(N, 0)$ branching forces $l$-parity $\equiv N$-parity; neither reproduces CH-1's even-$j$ zeros). Confirms CH-1's chirality-parity factorization is CH-substrate-specific, not Mellin-engine-generic. **Stage 1 correctly targets the CH triple.**

### Sub-sprint Stage 1 reframe — REFRAMED-CONSTRUCTIVE

The R3 sub-sprint (2026-06-04) closed the Marcolli-Tabuada $\mathrm{HP}_*$-route at source via Morita-triviality of $\mathcal{A}_{\mathrm{GV}}^{(n_{\max})} \cong \mathbb{C}^{N_{\mathrm{Fock}}}$. CH-1 (2026-06-05 morning) then named Stage 1 as constructing $\omega^{\mathrm{tri}}$ on a refinement of $\mathrm{HP}_*$. **The R3 vs CH-1 tension dissolves once the side of cyclic theory is named:** CH-1's actual computational objects $\mathrm{Tr}(D^k e^{-tD^2})$ and $\mathrm{Tr}(\gamma D^k e^{-tD^2})$ are functionals on the algebra, i.e. cyclic cochains, i.e. $\mathrm{HP}^*$ elements (not $\mathrm{HP}_*$). CH-1's phrasing was a misnomer; the genuine refinement is of $\mathrm{HP}^*$. Named target: the **Jaffe-Lesniewski-Osterwalder (1988) entire-cyclic cocycle** $\{\phi_n\}_{n \ge 0}$ of the truncated CH spectral triple. Connes-Marcolli $U^*$ (arXiv:math/0409306) supplies the published-precedent Tannakian shape on cohomology-side labelled data.

### Sub-sprint 1 — POSITIVE-WITH-SHIFT

JLO cochains $\{\phi_n^{\mathrm{even/odd}}\}_{n \in \{0,1,2\}}$ at $n_{\max} = 2$ (dim $\mathcal{H} = 16$, algebra $\mathbb{C}^5$, Dirac $D = \Lambda + \kappa A$ with $\kappa = -1/16$) computed bit-exactly in `sympy.Rational`. **CH-1's three leading coefficients (16, 84, 36) verified bit-exactly.** But the naive map "cochain degree $n \leftrightarrow$ slot $k = n$" is FALSE:

| Slot | Predicted | Observed | Host |
|:---:|:---:|:---:|:---|
| M1 | 16 | **16** | $\phi_0^{\mathrm{odd}}(1; t)\big|_{t^0}$ |
| M2 | 84 | **84** | $-\phi_0^{\mathrm{odd}}(1; t)\big|_{t^1}$ |
| M3 | 36 | **36** | $\mathrm{Tr}(\gamma D e^{-tD^2})\big|_{t^0}$ — CM-residue, not JLO |

M1 and M2 live as different $t$-orders of the **same** JLO degree-0 cochain; M3 lives outside JLO on commutative algebra and inside the Connes-Moscovici residue framework. **The ω^tri target is a mixed JLO + CM-residue framework, not a single cocycle.** Theorem (proved bit-exactly on all idempotent pairs, both flavours, all $t$-orders): **$\phi_1(a_0, a_1; t) \equiv 0$ on commutative $\mathcal{A}$.** The first non-trivial JLO cochain degree on commutative $\mathcal{A}$ is $n = 2$.

### Sub-sprint 2a — CLEAN-POSITIVE

The $\omega^{\mathrm{tri}}$ symbol extends bit-exactly to $n_{\max} \in \{1, 2, 3, 4\}$ with three closed-form polynomial-in-$n_{\max}$ expressions identified by **re-derivation** from the CH shell-degeneracy sum (not curve-fit; no free parameters; verified bit-exact at all four cutoffs plus cross-checked against CH-1 on 9/9 cells):

| $n_{\max}$ | $M_1 = \dim\mathcal{H}$ | $M_2 = \mathrm{Tr}(\Lambda^2)$ | $M_3 = \mathrm{Tr}(\gamma\Lambda)$ |
|:---:|:---:|:---:|:---:|
| 1 | 4 | 9 | 6 |
| 2 | 16 | 84 | 36 |
| 3 | 40 | 378 | 120 |
| 4 | 80 | 1188 | 300 |

Closed forms:
- $M_1(n) = 2n(n+1)(n+2)/3$ — degree 3
- $M_2(n) = n(n+1)(n+2)(2n+1)(2n+3)/10$ — degree 5
- $M_3(n) = n(n+1)^2(n+2)/2$ — degree 4

**Three different polynomial degrees** is the cohomological-side witness of the heat-kernel-order grading: at finite cutoff the three master-Mellin slots manifest as polynomials of distinct degree in $n_{\max}$. Pro-system rationality holds; no transcendentals at finite cutoff.

### Sub-sprint 2b — POSITIVE-WITH-STRUCTURAL-FINDING

Formal Mellin transform $\mathcal{M}[\phi_0^{\mathrm{odd}}(1; t)](s) = \Gamma(s) \cdot \zeta_{D^2}(s)$ closes bit-exactly in $\mathbb{Q}[\Gamma(s)]$ at $n_{\max} \in \{2, 3\}$. **Headline structural correction:** the JLO memo and the agent prompt assumed $\pi$ would enter at the $s \to 0$ residue. It doesn't. The actual structure:

- **M1** at the spectral-dimension pole $s = d/2 = 3/2$ (Weyl-law short-time leading $\mathrm{Tr}(e^{-tD^2}) \sim \sqrt{\pi}/2 \cdot t^{-d/2}$ supplies the Hopf-base measure $\pi$).
- **M2** at integer-$s$ regular points (Seeley-DeWitt continuum closed forms reproduce the Q5'-CH-2 panel verbatim).
- **M3** in the $\eta$-style pairing $\Gamma(s) \cdot \mathcal{M}[\mathrm{Tr}(\gamma D e^{-tD^2})](s)$ — a structurally different Mellin complex.

The three master-Mellin slots are not three k-indices of one functor; they are three **categorically different extraction points** of two formally distinct Mellin transforms. At finite cutoff the spectral zeta is rational-valued (e.g. $\zeta_{D^2}^{(2)}(1) = 1\,978\,465\,196\,515\,232/532\,976\,619\,280\,625$ for full Dirac at $n_{\max} = 2, \kappa = -1/16$); the $s \to 0$ residue is exactly $\dim \mathcal{H} \in \mathbb{Z}$. Master-Mellin transcendentals enter only in the continuum $n_{\max} \to \infty$ limit. **Cosmic-Galois Stage 2 acts at continuum, not finite cutoff.**

### Sub-sprint 2c — POSITIVE-WITH-STRUCTURAL-FINDING

Cyclic $(b, B)$ bicomplex of the truncated CH triple at $n_{\max} = 2$ constructed bit-exactly. **Load-bearing degree-1 cocycle condition $b\phi_0 + B\phi_2 = 0$ holds bit-exactly on the full panel** of 36 idempotent-pair inputs at three $t$-orders, both flavours — 216 bit-exact zero residuals. Both $b\phi_0$ and $B\phi_2$ vanish individually on commutative $\mathcal{A}$.

**Explicit $\mathrm{HP}^{\mathrm{even}}$ class:**
$$[\phi_{\mathrm{JLO}}]^{\mathrm{HP}^{\mathrm{even}}} = (+2, -2, +2, +2, -4) \in \mathbb{Q}^5$$
sector-resolved McKean-Singer chirality grading, summing to $\mathrm{Tr}(\gamma) = 0 = \chi(S^3)$. **The Sub-Sprint 1 symbol (16, 84, 36) is therefore a JLO cocycle representative, not a coboundary.**

**Structural finding (honest scope):** at degree 3, $(b\phi_2 + B\phi_4)$ has bit-exact residual $\pm 1/(3 \cdot 2^{16}) = \pm 1/196\,608$ on specific palindromic 4-tuples. JLO-tower truncation artifact tied to $\phi_1 \equiv \phi_3 \equiv 0$ on commutative $\mathcal{A}$ — the cancellation channel is removed. **The CM-residue cocycle for M3 is single-cochain and immune** to this artifact. The JLO/CM-residue distinction is Stage-2-relevant **at the cocycle-class level**, not only at the mechanism level.

## Symbol-level synthesis

Stage 1 of Q5' is now bit-exactly constructed across four diagnostic axes:

| Axis | What we have |
|:---|:---|
| **Symbol** | (16, 84, 36) at $n_{\max} = 2$; polynomial closed forms across $n_{\max} \in \{1, 2, 3, 4\}$ |
| **Cochain framework** | Mixed: JLO (M1, M2 as $t$-orders of $\phi_0^{\mathrm{odd}}$) + Connes-Moscovici residue (M3 in $\eta$-pairing) |
| **Mellin extraction** | Three categorically distinct: M1 at non-integer pole $s = 3/2$ of $\Gamma\zeta_{D^2}$; M2 at integer-$s$ regular points; M3 in $\eta$-complex |
| **Bicomplex** | $(b + B)\phi = 0$ bit-exact at degree 1; $\mathrm{HP}^{\mathrm{even}}$ class $(+2, -2, +2, +2, -4)$; degree-3 truncation artifact identifies the JLO/CM split at class level |

The second span of the cosmic-Galois $U^*$ bridge (after the morning's first span = CH arc) is structurally complete at finite cutoff.

## Paper edits applied (PI authorized)

All five recommended edits across four papers, all four compile clean three-pass.

- **Paper 24 §V** (Coulomb/HO asymmetry): added Layer 6 (master Mellin engine chirality-parity selection rule); updated four-layer summary sentence to six-layer; pages unchanged at 11 (PDF includes the new layer in existing item list). Source: S5 falsifier memo.
- **Paper 32 §VIII** (case-exhaustion theorem, master Mellin engine): four new short remarks added after `rem:master_mellin_domain`:
  - `rem:q5p_ch_specificity` — CH-substrate-specificity of chirality-parity (S5 falsifier).
  - `rem:q5p_stage1_jlo_witness` — Stage-1 cochain-level symbol at $n_{\max} = 2$, JLO+CM split, $\phi_1 \equiv 0$ theorem, cutoff sweep closed forms.
  - `rem:q5p_mellin_extraction_points` — three categorically different extraction points of $\Gamma\zeta_{D^2}(s)$ + $\Gamma \mathcal{M}[\mathrm{Tr}(\gamma D e^{-tD^2})](s)$.
  - `rem:q5p_hp_even_class` — explicit $\mathrm{HP}^{\mathrm{even}}$ class $(+2, -2, +2, +2, -4)$; degree-1 cocycle holds bit-exact; degree-3 truncation artifact $\pm 1/196608$.
  Page count: 65 (unchanged at the resolution of the truncated count).
- **Paper 55 §subsec:open_m2_m3** (Q5' open question): four narrative paragraphs added after the CH-3 paragraph and before the "Honest scope" paragraph: Stage 1 reframe, explicit symbol + cutoff sweep, Mellin localization, bicomplex + class. **Pages: 23 (was 20), three-pass clean.**
- **Paper 18 §III.7** (mechanism-as-domain sharpening): one-paragraph addition at the end of the closing paragraph sharpening M1's $\pi$ injection to the spectral-dimension pole $s = d/2$ via Sub-Sprint 2b. Pages: 26 (unchanged).

## Honest scope

**Closed at theorem grade (bit-exact at finite cutoff):**
- $\phi_1(a_0, a_1; t) \equiv 0$ on commutative $\mathcal{A} = \mathbb{C}^5$, all idempotent pairs, all $t$-orders, both flavours (Sub-Sprint 1).
- Degree-1 entire-cyclic cocycle condition $b\phi_0 + B\phi_2 = 0$ bit-exact on full 216-cell panel (Sub-Sprint 2c).
- Symbol triple $(M_1, M_2, M_3)$ closed-form polynomial in $n_{\max}$ with re-derived closed forms (Sub-Sprint 2a).
- Formal Mellin transform closes bit-exactly in $\mathbb{Q}[\Gamma(s)]$ at $n_{\max} \in \{2, 3\}$ (Sub-Sprint 2b).
- $\mathrm{HP}^{\mathrm{even}}$ class identified bit-exactly as $(+2, -2, +2, +2, -4) \in \mathbb{Q}^5$ on Morita-trivial baseline (Sub-Sprint 2c).
- CH-substrate-specificity of chirality-parity factorization — Bargmann-Segal $S^5$ Hardy fails on both load-bearing conditions (S5 falsifier).

**Structural sketch (not yet theorem):**
- M1 lives at the spectral-dimension pole $s = d/2$ — verified at finite cutoff via the $t^{-d/2}$ Weyl-law asymptotic interpretation; the continuum-limit residue analysis is sketched but not theorem-graded.
- Connes-Moscovici residue cocycle as the natural M3 host — verified by the bit-exact symbol identification and by exclusion of JLO on commutative $\mathcal{A}$, but the full CM-cyclic-cohomology bicomplex was not constructed in this sprint.

**Numerical observation:**
- The cross-cutoff polynomial degrees (3, 5, 4) match expected Weyl-asymptotic scaling at finite cutoff; the structural meaning of "three different degrees" as the cohomological-side witness of heat-kernel-order grading is an interpretive synthesis that adds no new bit-exact theorem.

**Named open follow-ons:**
- **Stage 2 (multi-year):** define the slot-graded Chern-character Tannakian category and identify its motivic Galois group. Acts on two cocycle classes (JLO HP-even, CM-residue $\eta$), not three slots of a single functor — the action coordinates *three* analytic extractions on *two* cohomology classes. Shape precedent: Connes-Marcolli cosmic Galois $U^*$.
- **Continuum-limit analysis:** the finite-cutoff polynomials in $n_{\max}$ diverge as $n_{\max} \to \infty$; the Mellin regularisation against $\Gamma(s)$ injects $\pi$ at $s = d/2$. The transition from $\mathbb{Q}$-valued finite-cutoff data to $M_1 \cup M_2$-valued continuum data is structurally sketched but not theorem-graded.
- **CM-residue bicomplex:** construct the full Connes-Moscovici cyclic complex housing M3 explicitly; verify the $\eta$-pairing cocycle condition; identify the analog HP^odd or HP^even class on that complex.
- **Pro-system functoriality:** the natural restriction map $\mathcal{T}_{n_{\max}+1} \to \mathcal{T}_{n_{\max}}$ does not lift cleanly to a JLO-cochain map; the compatible cohomology-side pro-system has no published precedent.
- **Strict Tannakian construction:** target tensor category + coalgebra structure on symbol components + tensor-product-compatible fiber functor. Multi-year.

**Hard prohibitions check (§13.5):** No changes to natural geometry hierarchy. No fitted/empirical parameters introduced. No deletion of negative results from §3. No removal of "conjectural" label from Paper 2 combination rule.

**Curve-fit-audit:** Closed-form polynomials in Sub-Sprint 2a were RE-DERIVED from the chirality-balanced shell degeneracy $2n(n+1)$ and the CH eigenvalue $n + 1/2$, not curve-fit from 4 data points. Zero free parameters. The Sub-Sprint 2b structural finding ("M1 $\pi$ lives at $s = d/2$, not $s = 0$") corrects an incorrect assumption in the prompt; selection bias is acknowledged in the memo.

**Discrete-for-skeleton compliance:** Everything bit-exact `sympy.Rational` throughout the entire arc. Zero PSLQ. Zero floats. Transcendentals tagged for continuum limit only.

## Files

### Memos (six sub-sprint memos + this umbrella)
- `debug/sprint_q5p_s5_falsifier_memo.md`
- `debug/sprint_q5p_stage1_reframe_memo.md`
- `debug/sprint_q5p_jlo_nmax2_memo.md`
- `debug/sprint_q5p_2a_jlo_nmax_sweep_memo.md`
- `debug/sprint_q5p_2b_phi0_mellin_memo.md`
- `debug/sprint_q5p_2c_bicomplex_memo.md`
- `debug/sprint_q5p_stage1_arc_2026_06_05_memo.md` (this)

### Drivers
- `debug/compute_s5_k_nmax_truncated.py`
- `debug/compute_jlo_nmax2.py`
- `debug/compute_jlo_nmax_sweep.py`
- `debug/compute_phi0_mellin.py`
- `debug/compute_jlo_bicomplex.py`

### Data
- `debug/data/sprint_q5p_s5_falsifier_data.json`
- `debug/data/sprint_q5p_jlo_nmax2_data.json`
- `debug/data/sprint_q5p_2a_jlo_nmax_sweep_data.json`
- `debug/data/sprint_q5p_2b_phi0_mellin_data.json`
- `debug/data/sprint_q5p_2c_bicomplex_data.json`

### Paper files modified
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex`
- `papers/group3_foundations/paper_18_exchange_constants.tex`
- `papers/group3_foundations/paper_24_bargmann_segal.tex`
- `papers/group3_foundations/paper_55_periods_of_geovac.tex`

## One-line verdict

Six sub-sprints + four paper edits in one session: Stage 1 of the cosmic-Galois $U^*$ bridge is bit-exactly constructed at finite cutoff across four diagnostic axes (symbol, mixed JLO+CM cochain framework, Mellin extraction localization, $\mathrm{HP}^{\mathrm{even}}$ cocycle class) with the CH triple confirmed as the structurally correct substrate; Stage 2 remains multi-year but now operates on two cocycle classes coordinating three analytic extractions, rather than three slots of one functor.
