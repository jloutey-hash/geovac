# Sprint Q5'-CH-1 — Master Mellin engine sources on the truncated CH triple

Date: 2026-06-05
Scope: 1-day scoping sprint of the Q5' multi-year cosmic-Galois $U^*$
program. First stone of the bridge: compute the master Mellin engine
sources $G_k(t) = \mathrm{Tr}(D^k\,e^{-tD^2})$ at $k \in \{0, 1, 2\}$ on the
truncated Camporesi--Higuchi spectral triple at
$n_{\max} \in \{2, 3, 4\}$, then check whether the bit-exact rational
structure aligns with the master Mellin partition M1/M2/M3
(Paper 32 §VIII Thm `thm:pi_source_case_exhaustion` +
Rem `rem:master_mellin_domain`).

Driver: `debug/compute_ch_k_nmax_truncated.py`. Data:
`debug/data/sprint_q5p_ch1_data.json`. Wall time: 0.9 s total over the
three cutoffs. No floating point; no PSLQ at this stage.

## TL;DR

**Verdict: CLEAN POSITIVE — the master Mellin k-slot factors bit-exactly
at every finite $n_{\max} \in \{2, 3, 4\}$ into a chirality-parity Z/2
component times a heat-kernel-order component.**

Specifically: the chirality-balance + spectral-symmetry of the CH triple
forces the bit-exact selection rule

| $k$ | plain $\mathrm{Tr}(D^k\,e^{-tD^2})$ | supertrace $\mathrm{Tr}(\gamma\,D^k\,e^{-tD^2})$ | Mellin slot |
|:---:|:-----------------------------------:|:------------------------------------------------:|:------------|
| 0   | nonzero, all small-$t$ coefficients | **identically zero** (exact rational 0)          | **M1**      |
| 1   | **identically zero** (exact rational 0) | nonzero, all small-$t$ coefficients          | **M3**      |
| 2   | nonzero, all small-$t$ coefficients | **identically zero** (exact rational 0)          | **M2**      |

Holds at every panel cell with the diagonal CH Lambda *and* with the full
Dirac $D = \Lambda + \kappa\,A$ at $\kappa = -1/16$ (the Paper 0 topological
constant), $j$-moment range $j \in \{0, \ldots, 10\}$ for $\Lambda$ and
$j \in \{0, \ldots, 6\}$ for full $D$ on $n_{\max} = 4$ (dim$\,= 80$). The
parity selection rule is therefore not an artifact of the diagonal-only
case — the $E_1$ adjacency $A$ preserves it, because $A$ flips $l$ by
$\pm 1$ and thus respects the chirality-balance structure of the CH
sector.

**Structural sharpening of the k-slot Tannakian memo (Sprint Q5'-k-slot,
2026-06-04, `debug/sprint_q5p_k_slot_tannakian_memo.md`).** The memo's
BORDERLINE verdict was that the $\mathrm{Z}/3$ k-slot is TANNAKIAN-INVISIBLE
on the standard $\mathrm{HP}_*$ fiber functor (which is $\mathrm{Z}/2$-graded,
not $\mathrm{Z}/3$). The bit-exact finite-$n_{\max}$ data now shows
**exactly which part of the k-slot $\mathrm{HP}_*$ sees and which part it
misses:**

- The **chirality dimension** of the k-slot (M3 separated from M1+M2) IS
  visible to $\mathrm{HP}_*$, because the chirality grading $\gamma$ is
  exactly the standard $\mathrm{Z}/2$ split of cyclic homology
  (plain trace = HP$_0$ sector, supertrace = HP$_1$ sector).
- The **heat-kernel-order dimension** within HP$_0$ (M1 vs M2) is NOT
  visible to $\mathrm{HP}_*$ — both M1 and M2 outputs sit in the same
  chirality sector and are distinguished only by the operator order
  $D^0$ vs $D^2$.

Q5' therefore sharpens from "is the Z/3 k-slot Tannakian-relevant?" to
the more focused

> **Q5'-sharp:** Is the heat-kernel-asymptotic-order distinction within
> a single $\mathrm{HP}_*$ chirality sector recoverable as Tannakian
> content via the candidate enrichment $\omega^{\mathrm{tri}}$?

This is the *operator-order grading* discussed in Paper 32 §VIII
Rem `rem:master_mellin_domain` ("each $k$ classifies both the
sub-mechanism and the natural domain of observable"). The finite-$n_{\max}$
data places it concretely at the bit-exact level: M1 and M2 differ by
heat-kernel ORDER ($t^0$ vs $t^2$ in the small-$t$ expansion of the
spectral-zeta-times-$\Gamma(s)$ pairing), not by $\mathrm{HP}_*$ sector.

## Computational summary

Three sample data points from the panel:

**$n_{\max} = 2$, $\dim\,\mathcal{H} = 16$, $\chi$-balance: 8/8.**

| moment | $j = 0$ | 1 | 2 | 3 | 4 | 5 | 6 |
|:-------|:-------:|:-:|:-:|:-:|:-:|:-:|:-:|
| $\mathrm{Tr}(\Lambda^j)$ | 16 | 0 | 84 | 0 | 489 | 0 | 11901/4 |
| $\mathrm{Tr}(\gamma\Lambda^j)$ | 0 | 36 | 0 | 201 | 0 | 4809/4 | 0 |
| $\mathrm{Tr}(D^j)$ (full) | 16 | 0 | 5401/64 | 0 | 8142371/16384 | 0 | 6449682761/2097152 |
| $\mathrm{Tr}(\gamma D^j)$ | 0 | 36 | 0 | 25947/128 | 0 | 40306423/32768 | 0 |

**$n_{\max} = 3$, $\dim\,\mathcal{H} = 40$, $\chi$-balance: 20/20.**

| moment | $j = 0$ | 1 | 2 | 3 | 4 | 5 | 6 |
|:-------|:-------:|:-:|:-:|:-:|:-:|:-:|:-:|
| $\mathrm{Tr}(\Lambda^j)$ | 40 | 0 | 378 | 0 | 8181/2 | 0 | 376749/8 |
| $\mathrm{Tr}(\gamma\Lambda^j)$ | 0 | 120 | 0 | 1230 | 0 | 27615/2 | 0 |

**$n_{\max} = 4$, $\dim\,\mathcal{H} = 80$, $\chi$-balance: 40/40.**

| moment | $j = 0$ | 1 | 2 | 3 | 4 | 5 | 6 |
|:-------|:-------:|:-:|:-:|:-:|:-:|:-:|:-:|
| $\mathrm{Tr}(\Lambda^j)$ | 80 | 0 | 1188 | 0 | 20493 | 0 | 1516977/4 |
| $\mathrm{Tr}(\gamma\Lambda^j)$ | 0 | 300 | 0 | 4875 | 0 | 350475/4 | 0 |

The parity zeros are exact rational 0, not numerical small. Every nonzero
entry is a rational $p/2^m$ for some small $m$ (the only denominators come
from $\lambda = n+1/2$ raised to a positive power, so $2^j$ at worst).

## Structural argument (why the parity is exact)

Two facts force the bit-exact selection rule:

1. **Spectral symmetry of $\Lambda$:** the CH spectrum is
   $\{\pm(n + 1/2) : n = 1, \ldots, n_{\max}\}$ with equal multiplicity at
   $\pm$ values. Therefore
   $\mathrm{Tr}(\Lambda^j) = \sum_i \lambda_i^j$ vanishes whenever $j$
   is odd (each $\lambda > 0$ contribution is cancelled by its negative
   partner).
2. **Chirality-balance + diagonal alignment:** on the CH triple,
   $\gamma_{ii} = \chi_i$ and $\Lambda_{ii} = \chi_i |\lambda_i|$, so
   $(\gamma \Lambda^j)_{ii} = \chi_i^{j+1} |\lambda_i|^j$. For $j$ even,
   $\chi^{j+1} = \chi$ and the sum is the signed sum
   $\sum_i \chi_i |\lambda_i|^j$, which vanishes by chirality-balance. For
   $j$ odd, $\chi^{j+1} = 1$ and the sum is the unsigned sum
   $\sum_i |\lambda_i|^j > 0$.

Both facts are *structural* properties of the CH Dirac on $S^3$ at any
$n_{\max}$, and both survive the $\kappa A$ perturbation because $A$ is
the parity-respecting $E_1$ dipole adjacency (verified numerically above).

The parity selection rule is the operator-system finite-cutoff witness of
the standard fact that on a closed even-dimensional spin manifold, the
Dirac spectrum $\{\pm |\lambda|\}$ is symmetric and the supertrace of
even powers of $D$ vanishes. On an *odd*-dim manifold like $S^3$, the
γ does not anticommute with $D$ as it would on an even manifold — but it
still satisfies the diagonal-alignment fact above, and that's enough to
force the parity rule on Tr/supertrace of $\Lambda^j$ at every finite
cutoff.

## Implications for the next sprint

Q5' is now structured as a forced two-stage research program:

**Stage 1 (could-be-sprint-scale, ~2-4 weeks):** Construct the candidate
enrichment $\omega^{\mathrm{tri}}: \mathrm{dg}(\mathcal{T}_{n_{\max}}) \to
\mathrm{Vec}_{\mathbb{Q}} \otimes \mathrm{IndexCat}(\{0, 1, 2\})$
explicitly on the truncated CH triple. The Q5'-CH-1 data above is the
finite-cutoff witness: $\omega^{\mathrm{tri}}$ on $G_0(t)$ lands in the
M1 sector with leading coefficient $\dim \mathcal{H}$, on $G_2(t)$ in
the M2 sector with leading coefficient $\mathrm{Tr}(\Lambda^2)$, on
$S_1(t)$ in the M3 sector with leading coefficient $\mathrm{Tr}(\gamma\Lambda)$.
The structural symbol level (which $G$ vs $S$, which power of $t$) is
already explicit; what remains is to lift it to a $\mathrm{Z}/2 \times \mathrm{Z}$-graded
fiber functor on the dg-category.

**Stage 2 (genuine multi-year, the original Q5' target):** Define the
candidate enriched Tannakian category and ask whether its motivic Galois
group has a non-trivial quotient acting on the operator-order grading
within HP$_0$. This is what no published precedent does
(Marcolli-Tabuada 2016 uses HP$_*$, Z/2-graded only; the enrichment
beyond is what would be new).

The bit-exact finite-cutoff data narrows Stage 1 to something
sprint-scale, which is the genuine first stone of the bridge — not the
whole bridge. Stage 2 remains the multi-year frontier.

## Caveats and honest scope

1. **No periods appear at finite $n_{\max}$.** Every quantity computed
   here is bit-exact rational. The M1/M2/M3 period content (π, $\sqrt\pi$,
   Catalan $G$) enters only in the Mellin transform via $\Gamma(s)$ and
   in the continuum limit $n_{\max} \to \infty$. This sprint did not
   compute the continuum limit; it showed the FINITE-CUTOFF structural
   witness of the partition.
2. **The parity selection rule is automatic on any spin manifold.** It
   is not a GeoVac-specific finding; it is a property of the CH Dirac on
   a closed manifold. The GeoVac-specific finding is that the
   *truncated* CH triple at any finite $n_{\max}$ also exhibits it
   bit-exactly, and that this is exactly the bit of the master Mellin
   k-slot that $\mathrm{HP}_*$ sees.
3. **No Tannakian category has been constructed.** This sprint is
   scoping data, not the enrichment construction. Stage 1 (above) is the
   forward sprint.
4. **No PSLQ identification.** At finite $n_{\max}$ everything is
   rational; PSLQ at this stage would identify the rational against
   itself. PSLQ would enter at the continuum-limit Mellin transform
   stage, which is a separate computation.
5. **Curve-fit-audit (per `feedback_audit_numerical_claims`):** the
   claim "the k-slot factors as Z/2 chirality × heat-kernel order"
   has zero free parameters (no fitting), the selection-bias is
   minimal (the master Mellin engine partition was already published
   in Paper 32 §VIII; the chirality-parity is automatic on CH).
   Alternative: the parity selection rule might be specific to the
   CH spectral triple and break on other discrete spectral triples
   (e.g., on a non-spherical manifold with non-symmetric spectrum). The
   robustness test would be to repeat on the Bargmann-Segal Hardy-on-$S^5$
   triple (Paper 24) — where the spectrum is *not* symmetric (only
   positive eigenvalues) and the parity rule should *break*. This is a
   clean follow-on test.
6. **Discrete-for-skeleton compliance (per
   `feedback_discrete_for_skeleton`):** the computation is exact sympy
   `Rational` throughout; the parity selection rule is a Layer 1
   skeleton observation, properly stated in bit-exact form, no PSLQ +
   $n_{\max}$ sweep needed.

## Files produced

- `debug/compute_ch_k_nmax_truncated.py` — driver (~220 lines, fast path
  for diagonal matrices, full-D moments at $j \le 6$).
- `debug/data/sprint_q5p_ch1_data.json` — exact rational data per
  $n_{\max}$.
- `debug/sprint_q5p_ch1_memo.md` — this memo.

## Recommended paper edit (PI to apply, decline, or modify)

### Paper 55 §subsec:open_m2_m3 (Q5') — refine the sharpening from Q5'-k-slot

The 2026-06-04 Q5'-k-slot recommended a one-paragraph sharpening of Q5'.
This sprint refines that further. Suggested replacement paragraph:

> *Bit-exact finite-cutoff witness of the k-slot partition (Sprint
> Q5'-CH-1, June 2026; memo
> `debug/sprint_q5p_ch1_memo.md`).* The master Mellin source
> $\mathrm{Tr}(D^k\,e^{-tD^2})$ on the truncated Camporesi--Higuchi
> spectral triple at $n_{\max} \in \{2, 3, 4\}$ exhibits a bit-exact
> chirality-parity selection rule:\ plain trace nonzero for $k \in
> \{0, 2\}$ and identically zero for $k = 1$; supertrace identically
> zero for $k \in \{0, 2\}$ and nonzero for $k = 1$. This factors the
> Z/3 master Mellin slot index into a Z/2 chirality component (M3
> separated from M1+M2) times a heat-kernel-order component within
> HP$_0$ (M1 vs M2). The Z/2 chirality component is exactly the
> standard $\mathrm{HP}_*$ grading, hence Tannakian-visible via
> Marcolli--Tabuada 2016. The heat-kernel-order component within
> HP$_0$ — distinguishing M1 ($k = 0$) from M2 ($k = 2$) at the
> level of the small-$t$ expansion of $\mathrm{Tr}(D^k\,e^{-tD^2})$ —
> is the genuinely Tannakian-invisible component on the standard
> fiber functor. Q5' sharpens to:\ does the candidate enrichment
> $\omega^{\mathrm{tri}}$ make the heat-kernel order Tannakian-visible
> via a $\mathrm{Z}/2 \times \mathbb{Z}$-graded refinement of
> $\mathrm{HP}_*$?

(Recommendation only; no Paper 55 edits applied. PI direction sought.)

## One-line verdict

CLEAN POSITIVE — bit-exact chirality-parity factorization of the master
Mellin k-slot at every finite $n_{\max}$, sharpening Q5' from
"Z/3 k-slot vs Z/2 HP$_*$" to "heat-kernel-order grading within a
single HP$_*$ chirality sector"; first concrete stone of the cosmic-Galois
$U^*$ multi-year bridge.
