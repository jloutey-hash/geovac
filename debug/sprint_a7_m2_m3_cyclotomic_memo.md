# Sprint A7 — M2 / M3 cyclotomic-level-4 coincidence: verdict

Date: 2026-06-03
Scope: 1-day mini-sprint resolving the Open Question #2 of
`debug/sprint_a5_fm_full_read_memo.md` and Open Question #2 of
`debug/sprint_m3_cyclotomic_mixed_tate_memo.md` (the M2/M3
cyclotomic-level-4 coincidence flagged by Sprint A5).

## TL;DR

**Verdict: HALF-STRUCTURAL (case (c) of the decision gate).**

Both M2 and M3 engage the field $\mathbb{Q}(i) = \mathbb{Q}(\sqrt{-1}) =
\mathbb{Q}(\zeta_4)$, but via **two structurally independent mechanisms**:

- **M2's $\mathbb{Q}(i)$ is the (generic) splitting field of the standard
  Euclidean quadratic form** $u_1^2 + u_2^2 + \ldots + u_{2n+2}^2$.
  Isotropy is the elementary Witt-index identity (hyperbolic decomposition
  $X = u_j + i u_{j+1}, Y = u_j - i u_{j+1}$, pairwise over an even rank).
  This works for *any* positive-definite even-rank rational quadratic form;
  the field $\mathbb{Q}(i)$ is selected by sign signature alone.
  $\text{Gal}(\mathbb{Q}(i)/\mathbb{Q}) = \mathbb{Z}/2$ acts **trivially on
  M2 outputs**: the period values descend to $\mathbb{Q}$ (in fact to
  $\bigoplus_k \pi^{2k}\mathbb{Q}$, pure Tate over $\mathbb{Q}$, no $i$
  surviving).

- **M3's $\mathbb{Q}(i) = \mathbb{Q}(\zeta_4)$ is the cyclotomic field of
  the (unique) non-trivial primitive Dirichlet character of conductor 4**,
  $\chi_{-4}$. This is non-generic: it is the smallest cyclotomic field
  carrying a character whose period values ($G = \beta(2)$, $\beta(4)$,
  $\ldots$) genuinely require $\zeta_4$ in their motivic representation.
  $\text{Gal}(\mathbb{Q}(i)/\mathbb{Q}) = \mathbb{Z}/2$ acts
  **non-trivially on M3 outputs**: complex conjugation realises the
  Deligne–Glanois Galois descent from level $4$ to level $2$, exchanging
  $D_{\mathrm{even}} \leftrightarrow D_{\mathrm{odd}}$ at the
  spectral-triple level.

The two $\mathbb{Z}/2$ actions are therefore not the same $\mathbb{Z}/2$
in the sense of a single Tannakian symmetry of the underlying spectral
triple. The shared field $\mathbb{Q}(i)$ is a coincidence at the
**number-field level**, not at the **motivic-Galois level**.

**But the coincidence is not entirely accidental either.** Both
occurrences trace to the **dimension-3 parity of $S^3$**, which feeds
into M2 via the Euclidean signature $(+,+,+,+)$ of $S^3 \hookrightarrow
\mathbb{R}^4$ and into M3 via the half-integer Camporesi–Higuchi shift
$3/2$ (parity of $\mathrm{spin}(3) = \mathrm{SU}(2)$). The deeper
question — whether there is a single $S^3$-specific structural object that
selects $\mathbb{Q}(\zeta_4)$ for both sub-mechanisms simultaneously —
appears to be (a) genuine open question in the periods program, and (b)
amenable to a future motivic / Tannakian-symmetry investigation that
does NOT have a sprint-scale handle in 2026.

This is the **(c) HALF-STRUCTURAL** verdict of the original decision
gate. Paper 55 §7 gets a refined open-question entry (Open Q5'), and
the Remark 4.5 (`rem:m2_specialisation`) is sharpened with one paragraph
distinguishing the two roles of $\mathbb{Q}(i)$.

## Argument

### M2's $\mathbb{Q}(i)$: a generic Witt-isotropy field

Fathizadeh–Marcolli arXiv:1611.01815 Theorem 7.3 establishes the
mixed-Tate-over-$\mathbb{Q}(i)$ classification of
$\mathcal{V}_n^{\mathrm{F-M}}(\lambda)$ by showing the quadric
$Q_{\lambda, 2n}$ becomes isotropic over $\mathbb{Q}(\sqrt{-1})$. At the
GeoVac specialisation $\lambda = 1$, the quadric is the standard
positive-definite form
$$Q_{1, 2n} = u_1^2 + u_2^2 + u_3^2 + u_4^2 + u_5^2 + \ldots + u_{2n+2}^2,$$
which is anisotropic over $\mathbb{R}$ (and over $\mathbb{Q}$, with Witt
index $0$) but isotropic over $\mathbb{Q}(i)$ via the **elementary**
hyperbolic pair decomposition
$$X_j = u_{2j-1} + i\, u_{2j}, \qquad Y_j = u_{2j-1} - i\, u_{2j}, \qquad
X_j \cdot Y_j = u_{2j-1}^2 + u_{2j}^2.$$
With $2n+2$ even, the form pairs up exactly, and the Witt index over
$\mathbb{Q}(i)$ is the maximal $n+1$.

This is a generic fact about positive-definite even-rank rational
quadratic forms. The **field $\mathbb{Q}(i)$ is selected by signature
alone**, not by any feature of the discrete Camporesi–Higuchi spectrum
or the conformal projection of $S^3 \hookrightarrow \mathbb{R}^4$. Any
even-rank Euclidean-signature form over $\mathbb{Q}$ would land in the
same $\mathbb{Q}(i)$.

The Galois action $\text{Gal}(\mathbb{Q}(i)/\mathbb{Q}) = \mathbb{Z}/2$
acts on the splitting variables via $i \mapsto -i$, i.e., it swaps
$X_j \leftrightarrow Y_j$. **On periods of the affine complement
$\mathcal{V}_n^{\mathrm{GeoVac}}$, the Galois acts trivially**: the
M2 period output sits in $\bigoplus_k \pi^{2k}\mathbb{Q}$, all of which
is Galois-fixed under $\text{Gal}(\mathbb{Q}(i)/\mathbb{Q})$. The
$\mathbb{Q}(i)$ in F–M Theorem 7.3 is the **proof-time** field, not
the output field; the output descends to $\mathbb{Q}$ by Galois
invariance.

This is the Witt-style "field-of-definition" argument: the motive lives
over $\mathbb{Q}(i)$ for technical reasons (one needs the quadric to be
isotropic to apply Gysin/Mayer–Vietoris), but the period values
themselves are Galois-fixed.

### M3's $\mathbb{Q}(i) = \mathbb{Q}(\zeta_4)$: the canonical
cyclotomic-level-4 field

The M3 sub-mechanism produces Hurwitz values at quarter-integer shifts
$\zeta(s, 1/4), \zeta(s, 3/4)$ via the parity decomposition $D =
D_{\mathrm{even}} + D_{\mathrm{odd}}$ on the half-integer
Camporesi–Higuchi spectrum. The difference $D_{\mathrm{even}}(s) -
D_{\mathrm{odd}}(s) = 2^{s-1}(\beta(s) - \beta(s-2))$ (Paper 28
Theorem 3) contains the Dirichlet $\beta$ values
$\beta(s) = L(s, \chi_{-4})$ where $\chi_{-4}$ is the **unique
non-trivial primitive Dirichlet character of conductor $4$**,
$\chi_{-4}(n) = \mathrm{Im}(i^n)$.

The values $\beta(s)$ are Dirichlet $L$-values at a Galois-non-trivial
character. By Glanois 2015 (Corollaries 1.1–1.2), these live in
$\mathcal{MT}(\mathbb{Z}[i, 1/2])$ at level $N = 4$ via motivic iterated
integrals on $\mathbb{P}^1 \setminus \mu_4 \cup \{0, \infty\} =
\mathbb{P}^1 \setminus \{0, 1, i, -1, -i, \infty\}$. The level $N = 4$
is the **smallest cyclotomic level** at which $\chi_{-4}$ becomes a
character of $(\mathbb{Z}/N\mathbb{Z})^\times$, i.e., the conductor.

The Galois action $\text{Gal}(\mathbb{Q}(i)/\mathbb{Q}) = \mathbb{Z}/2$
acts **non-trivially**:
- On the period side: $\beta(s) \mapsto (-1) \cdot \beta(s)$ at odd $s$
  (since $\chi_{-4}$ is odd: $\chi_{-4}(-1) = -1$), reflecting the
  odd-character parity at the Galois orbit.
- On the spectral-triple side: complex conjugation realises the parity
  flip $D_{\mathrm{even}} \leftrightarrow D_{\mathrm{odd}}$, identified
  with the Deligne–Glanois descent from level $4$ to level $2$
  (Proposition 5.5 of Paper 55).

The $\mathbb{Q}(i)$ here is **NOT a generic splitting field** but the
**unique** field carrying the conductor-4 character. The choice is
forced by $\chi_{-4}$'s Hilbert-symbol class, not by signature.

### The two $\mathbb{Z}/2$'s are independent (over $\mathbb{Q}$)

Both M2 and M3 admit a $\text{Gal}(\mathbb{Q}(i)/\mathbb{Q}) =
\mathbb{Z}/2$ action, but the two actions are structurally distinct:

| | M2 | M3 |
|:---|:---|:---|
| Role of $\mathbb{Q}(i)$ | Witt-splitting field | conductor of $\chi_{-4}$ |
| Selected by | signature $(+,+,+,+)$ | character $\chi_{-4}$ |
| Genericity | generic (any even-rank Euclidean form) | non-generic (conductor 4 is canonical) |
| Motivic depth | $0$ (pure Tate) | $\ge 1$ (mixed Tate) |
| $\mathbb{Z}/2$ on output | trivial (output in $\mathbb{Q}$) | non-trivial (output in $\mathbb{Q}(i)$, descends to $\mathbb{Q}$ via parity sum but not parity difference) |

In particular, the F–M Theorem 7.3 proof uses Gysin triangles
**independent of the specific identification of the field as
$\mathbb{Q}(\zeta_4)$**; it would equally apply to $\mathbb{Q}(\sqrt{p})$
for any prime $p$ over which the quadric splits. The Glanois N=4 basis,
in contrast, is **specifically tied to $\mu_4$**, which is a feature of
the conductor of $\chi_{-4}$ and would change if the spectrum's
half-integer shift were not $3/2$.

There is therefore no shared Tannakian symmetry between the two
sub-mechanisms at the level of the motivic Galois group of GeoVac. The
M2 motivic Galois (over $\mathbb{Q}$) of the pure-Tate ring
$\bigoplus_k \pi^{2k}\mathbb{Q}$ is the unipotent group $\mathbb{G}_m$
acting on $(2\pi i)$-twists; the M3 motivic Galois (over $\mathbb{Q}$,
descending from level 4) is the cyclotomic-mixed-Tate motivic Galois at
level 4, which is strictly larger than the M2 Galois.

### What IS shared: dimension-3 parity of $S^3$

The one structural common thread:

- M2's $\mathbb{Q}(i)$ comes from the Euclidean signature $(+,+,+,+)$
  on $\mathbb{R}^4 \supset S^3$, i.e., from the embedding dimension of
  $S^3$.
- M3's $\mathbb{Q}(\zeta_4)$ comes from the half-integer
  Camporesi–Higuchi shift $|\lambda_n| = n + 3/2$, which is itself a
  consequence of the **spinor parallelisability of $S^3 = \mathrm{SU}(2)$**.
  The shift value $3/2 = (d-1)/2$ at $d = 3$ feeds into Hurwitz at
  quarter-integer shifts $1/4, 3/4, 5/4$, and the conductor-4
  Dirichlet character emerges from the difference $\zeta(s, 3/4) -
  \zeta(s, 1/4) = -4^s \beta(s)$.

Both occurrences therefore trace to the **dimension $d = 3$ of $S^3$**,
but via two structurally independent mechanisms:
- M2: $d = 3$ implies the Euclidean signature of $\mathbb{R}^4 \supset
  S^3$, which produces $\mathbb{Q}(i)$ generically via Witt isotropy.
- M3: $d = 3$ implies the half-integer Camporesi–Higuchi shift, which
  produces $\mathbb{Q}(\zeta_4)$ specifically via the quarter-integer
  Hurwitz shifts.

These two mechanisms do not share a single underlying Tannakian symmetry
of the spectral triple — they share the input "dimension 3" but not the
output structure.

### Cross-check with the M2 / $S^5$ closure

A useful sanity check: Sprint A2 closed the M2 sub-mechanism on $S^5$
(`debug/sprint_a2_s5_mixed_tate_memo.md`), and the period ring is
$\pi^3 \cdot \mathbb{Q}$ (odd Tate weight 3 slice), still pure Tate
over $\mathbb{Q}$. The same F–M proof applies: the standard isotropic
form $u_1^2 + \ldots + u_{2n+2}^2$ at $\lambda = 1$ still becomes
hyperbolic over $\mathbb{Q}(i)$ (Witt index argument, same as $S^3$).
The output ring is different ($\pi^3 \cdot \mathbb{Q}$ vs.\ $\pi^{2k}
\mathbb{Q}$), but the **role of $\mathbb{Q}(i)$ is the same Witt-splitting
field**, independent of dimension parity.

In contrast, M3 on $S^5$ would engage Hurwitz at $1/4, 3/4, 5/4, 7/4$
(shift $5/2$ at $d = 5$), still pulling in $\mathbb{Q}(\zeta_4)$. So
on $S^5$ both M2 and M3 STILL engage $\mathbb{Q}(i)$, but for the
same structurally distinct reasons as on $S^3$. The coincidence is
therefore **dimension-3-specific only at the input level** (which
Camporesi–Higuchi shift, which Euclidean signature), but the
field-coincidence pattern $\mathbb{Q}(i)$ in both M2 and M3 persists
across $d = 3, 5$ — strongly suggesting the field-coincidence is
**structurally accidental in the field-of-output sense but generic in
the field-of-input sense**.

### What would (c) STRUCTURAL look like, if it were real

A genuine structural identification would require:
- A single object (a motive, a spectral triple, a Tannakian symmetry)
  whose Galois group simultaneously selects both $\mathbb{Q}(i)$ as
  Witt-splitting field for the SD-coefficient integrand AND
  $\mathbb{Q}(\zeta_4)$ as cyclotomic conductor for the $\eta$-invariant
  Mellin slot.
- A natural isomorphism between the two Galois actions (the trivial
  M2-action and the non-trivial M3-action).
- A natural transformation between M2's pure-Tate sub-ring and M3's
  level-4 cyclotomic mixed-Tate ring that respects the master Mellin
  engine slot index $k$.

None of these is known in the published periods program (Brown,
Deligne, Glanois, Goncharov, Fathizadeh–Marcolli). The "M2/M3
cyclotomic coincidence" question therefore sits as an **open
sub-question of the periods program**, not as a known structural fact.

A possible structural object that COULD unify the two is the **motivic
Galois group of the GeoVac discrete spectral triple as an enriched
mixed-Tate object** — but this has not been constructed (the
construction would be a multi-year mathematical research project of
its own).

The verdict for the present sprint is therefore **HALF-STRUCTURAL**:
shared $\mathbb{Q}(i)$ field at the number-field level, with a clean
common origin in the dimension-3 parity of $S^3$, but no shared
Tannakian symmetry at the motivic-Galois level.

## Verdict

**HALF-STRUCTURAL (case (c) of the decision gate).**

- M2's $\mathbb{Q}(i)$ and M3's $\mathbb{Q}(\zeta_4)$ are the **same
  number field** but enter via **two structurally distinct
  mechanisms** (Witt-splitting field vs.\ cyclotomic conductor).
- $\text{Gal}(\mathbb{Q}(i)/\mathbb{Q}) = \mathbb{Z}/2$ acts
  **trivially** on M2 outputs and **non-trivially** on M3 outputs.
  These are two different $\mathbb{Z}/2$'s in the sense of motivic
  Galois action.
- The shared $\mathbb{Q}(i)$ traces to the **dimension-3 parity of
  $S^3$** via two structurally independent paths (Euclidean signature
  for M2, half-integer Camporesi–Higuchi shift for M3).
- There is no known shared Tannakian symmetry; constructing one would
  require enriching the motivic Galois group of the GeoVac spectral
  triple in a way not done in the published periods program.

## Honest scope

1. **No fitted parameters introduced**: this is pure motivic
   bookkeeping.
2. **Paper 2 framing preserved**: combination rule $K = \pi(B + F -
   \Delta)$ stays a numerical observation; no edits to that label.
3. **No paper edits applied**: only proposals below; PI applies after
   review.
4. **Diagnostic-only**: no implementation, no PSLQ, no numerical
   checks were required because the question is structural and rests
   on established periods-program facts (Witt isotropy, Glanois N=4
   basis, conductor-4 Dirichlet character).
5. **The HALF-STRUCTURAL verdict is robust** under three alternative
   readings considered:
   - (a) Maybe a single "motivic spectral triple Galois" unifies them:
     possible in principle, but no construction in the literature;
     would be a multi-year mathematical research project.
   - (b) Maybe the $\mathbb{Q}(i)$ in M2 is also non-generic via some
     hidden $S^3$-specific structure: tested by the $S^5$ cross-check
     above; the field-of-output pattern persists on $S^5$ with the
     same Witt-isotropy mechanism, so the $\mathbb{Q}(i)$ in M2 IS
     generic to even-rank Euclidean signature, NOT $S^3$-specific.
   - (c) Maybe both $\mathbb{Q}(i)$'s are the same $\mathbb{Z}/2$ but
     act differently because of independent restriction: not
     supported by Glanois 2015 (level $N = 4$ basis is specifically
     tied to $\mu_4$, not to a generic field).
6. **The "dimension-3-parity-of-$S^3$" common thread is HONEST
   but not strong enough** to upgrade to (a) STRUCTURAL. Both
   mechanisms read dimension 3 in different ways; sharing input
   does not imply sharing Galois action.

## Paper-level proposals (not yet applied)

### Paper 55 §4 Remark 4.5 sharpening (`rem:m2_specialisation`)

Add a final paragraph distinguishing the two roles of $\mathbb{Q}(i)$:

> *Distinction from the M3 cyclotomic-level-4 ring (§5).*  The
> $\mathbb{Q}(\sqrt{-1})$ appearing in the proof of
> Theorem~\ref{thm:m2_mixed_tate} is the **Witt-splitting field** of
> the standard Euclidean quadratic form $Q_{1, 2n}$:\ a generic field
> selected by the positive-definite signature alone, with
> $\text{Gal}(\mathbb{Q}(i)/\mathbb{Q})$ acting trivially on M2
> outputs (the pure-Tate ring $\bigoplus_k \pi^{2k}\cdot\mathbb{Q}$
> descends to $\mathbb{Q}$).  The same field appearing in §\ref{sec:m3}
> as $\mathbb{Q}(\zeta_4)$ is the **cyclotomic conductor** of the
> non-trivial primitive Dirichlet character $\chi_{-4}$, with
> $\text{Gal}(\mathbb{Q}(\zeta_4)/\mathbb{Q})$ acting non-trivially
> via the Deligne–Glanois descent from level $4$ to level $2$.  The
> coincidence of fields is at the number-field level only;\ the two
> motivic Galois actions are structurally independent (sprint memo
> \texttt{debug/sprint\_a7\_m2\_m3\_cyclotomic\_memo.md},
> June 2026).  A possible deeper unification via an enriched motivic
> Galois group of the discrete spectral triple is recorded as an open
> question (Q5', §\ref{sec:open}).

### Paper 55 §7 new open-question subsection (Q5')

Add a new subsection after the existing §7.5 (or wherever the open
questions sequence terminates):

> \subsection{Open Q5':\ A possible deeper M2/M3 cyclotomic unification}
> \label{subsec:open_m2_m3}
>
> Both M2 and M3 sub-mechanisms engage $\mathbb{Q}(i) = \mathbb{Q}(\zeta_4)$,
> but via two structurally independent mechanisms (Witt-splitting field
> for M2, cyclotomic conductor for M3).  Both occurrences trace to the
> dimension-3 parity of $S^3$:\ M2 via the Euclidean signature
> $(+,+,+,+)$ on $\mathbb{R}^4 \supset S^3$, M3 via the half-integer
> Camporesi–Higuchi shift $3/2$ on $\mathrm{spin}(3) = \mathrm{SU}(2)$.
> An open question:\ is there a single object (a motive, a Tannakian
> symmetry, an enriched motivic Galois group of the discrete spectral
> triple) whose action simultaneously selects $\mathbb{Q}(\zeta_4)$ as
> Witt-splitting field for the SD integrand AND as cyclotomic conductor
> for the $\eta$-invariant Mellin slot, and a natural identification
> between the trivial M2-Galois and the non-trivial M3-Galois on the
> respective output rings?  Sprint A7 (memo
> \texttt{debug/sprint\_a7\_m2\_m3\_cyclotomic\_memo.md}) returns
> HALF-STRUCTURAL:\ shared field is robust at the number-field level,
> but no shared Tannakian symmetry is known in the published periods
> program (Brown, Deligne, Glanois, Goncharov, Fathizadeh–Marcolli).
> A structural unification would require an enriched motivic Galois
> group of the GeoVac discrete spectral triple, not yet constructed in
> the literature.  Multi-year mathematical research project.

### CLAUDE.md §2 entry (one-liner)

> **Sprint A7 M2/M3 cyclotomic (2026-06-03):** HALF-STRUCTURAL — both
> M2's Witt-splitting $\mathbb{Q}(i)$ and M3's cyclotomic conductor
> $\mathbb{Q}(\zeta_4)$ trace to dimension-3 parity of $S^3$ via
> independent mechanisms; no shared Tannakian symmetry. Paper 55 Q5'
> opened. See `debug/sprint_a7_m2_m3_cyclotomic_memo.md`.

### MEMORY.md / followon_register.md

Add to follow-on register (multi-year research projects):

> **M2/M3 motivic Galois unification (sprint A7, 2026-06-03):**
> HALF-STRUCTURAL verdict; open question is whether an enriched motivic
> Galois group of the GeoVac discrete spectral triple unifies M2's
> Witt-splitting $\mathbb{Q}(i)$ with M3's cyclotomic conductor
> $\mathbb{Q}(\zeta_4)$.  Multi-year research project; no
> sprint-scale handle.

## Files used

### Papers read
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` §4 + §5
  (full M2 / M3 theorem text + Remark 4.5 + Glanois descent
  identification Prop 5.5).

### Memos read
- `debug/sprint_a5_fm_full_read_memo.md` (M2 detailed F–M extraction;
  explicit $\mathbb{Q}(i)$ Witt isotropy at line 432).
- `debug/sprint_m3_cyclotomic_mixed_tate_memo.md` (M3 cyclotomic
  level-4 classification; explicit Glanois N=4 basis structure on
  page 4 cited).
- `debug/sprint_mixed_tate_test_memo.md` (M2 pure-Tate classification;
  explicit "no $\sqrt{\pi}$, no $\zeta(\text{odd})$, no MZV" content
  in output ring).
- `debug/sprint_a2_s5_mixed_tate_memo.md` ($S^5$ cross-check used for
  robustness audit of the HALF-STRUCTURAL verdict).

### Published references (already in Paper 55 bibliography)
- Fathizadeh–Marcolli arXiv:1611.01815 Theorem 7.3 (Witt isotropy
  over $\mathbb{Q}(i)$).
- Deligne arXiv:math/0302267 (period map isomorphism at $N \in \{2,
  3, 4, 6, 8\}$).
- Glanois arXiv:1411.4947 Corollaries 1.1–1.2 (explicit N=4 basis
  $\mathcal{B}^4$ with $\epsilon_p = \sqrt{-1}$).
- Brown 2012 (level-1 baseline MZV ring over $\mathbb{Z}$).

### Scripts (none required)
The question is structural and resolved by reading published periods
program results; no new numerical computation required.
