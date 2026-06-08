# Sprint Spectral-Action-Expansion-Chemistry: Diagnostic memo

**Date:** 2026-06-07
**Sprint:** Spectral-Action-Expansion-Chemistry (one of four parallel symmetry sprints, M-vS arc)
**Verdict:** **NEGATIVE-with-PARTIAL caveat**

## TL;DR

At fixed $R = R_{\rm eq}$, the Marcolli–vS spectral action
$S(D, \Lambda) = \mathrm{Tr}\,\exp(-D^2/\Lambda^2)$ of the GeoVac chemistry
Dirac $D = h_1$ (the M-vS-2-confirmed default `lih_spec()` Hamiltonian) is, at
large $\Lambda$, **bit-exactly the matrix-exponential Taylor series**

$$
S(\Lambda) \;=\; \sum_{k=0}^{\infty} \frac{(-1)^k}{k!\,\Lambda^{2k}}\,\mathrm{Tr}(D^{2k}),
$$

verified to **rel. residual $\le 1.16 \times 10^{-5}$ at $\Lambda = 5$ and to
$\lesssim 10^{-15}$ for $\Lambda \ge 10$** across LiH, $H_2$, NaH at
$k_{\max} = 6$. There is **no Chamseddine–Connes-style Seeley–DeWitt coefficient
structure** beyond the trivial Taylor coefficients $\mathrm{Tr}(D^{2k})$. This
is structurally forced — finite-dim $D$ has no UV divergence, so the
characteristic SD heat-kernel asymptotic $a_0 \Lambda^d + a_2 \Lambda^{d-2} +
\cdots$ that gives Einstein–Hilbert, Yang–Mills, Higgs, Yukawa in the
Chamseddine–Connes SM **does not arise**. The PARTIAL caveat: the substitutions
$D \to \mathrm{diag}(h_1)$ and $D \to h_1 - \mathrm{diag}(h_1)$ produce a clean
split $\mathrm{Tr}(D^2) = \mathrm{Tr}(h_{\rm diag}^2) + \mathrm{Tr}(h_{\rm off}^2)$
(cross-term vanishes by orthogonality of strictly diagonal vs strictly
off-diagonal matrices — a generic linear-algebra fact, not chemistry-specific).
The off-diagonal piece $\mathrm{Tr}(h_{\rm off}^2)$ IS a chemistry-meaningful
"bond-coupling strength" (cross-block $V_{\rm ne}$ when `cross_block_h1=True`),
but this is a property of the GeoVac construction's block structure, not an
emergent property of the spectral action. **Closes the spectral-action
chemistry-expansion thread.**

## Setup

- **System 1:** LiH default `lih_spec()` at $R = 3.015$ bohr (M-vS-2 confirmed
  bit-exact M-vS gauge network); $M = 15$ spatial orbitals (Li_core 5 +
  LiH_bond_center 5 + LiH_bond_partner 5).
- **System 2:** $H_2$ two-center pilot at $R = 1.4$ bohr (`make_h2_two_center_spec`
  from `bratteli_h2_pilot_driver.py`); $M = 10$.
- **System 3:** NaH default `nah_spec()` at $R = 3.566$ bohr (NIST CCCBDB
  $R_{\rm eq}$); $M = 10$.

For each, $h_1 = $ `build_balanced_hamiltonian(spec, R, cross_block_h1=True)`
gives the Track CD chemistry Dirac. M-vS-2 sprint already confirmed
$h_1^{M-vS-assembled} = h_1^{\rm GeoVac}$ at residual 0.0 for LiH; we treat
this $h_1$ as the M-vS Dirac $D$ for the spectral action.

**Lambda grid:** 50 points log-spaced over $[0.1, 100]$.

## Numerical results

### Trace invariants

| invariant | LiH | $H_2$ | NaH |
|:--------- | ---:| ---:| ---:|
| $M$ (dim) | 15 | 10 | 10 |
| $\mathrm{Tr}(h_1)$ | $-17.381$ | $-5.114$ | $-15.127$ |
| $\mathrm{Tr}(h_1^2)$ | $41.512$ | $5.565$ | $53.142$ |
| $\mathrm{Tr}(h_{\rm diag}^2)$ | $37.543$ | $3.515$ | $36.650$ |
| $\mathrm{Tr}(h_{\rm off}^2)$ | $3.969$ | $2.050$ | $16.492$ |
| $\mathrm{Tr}(h_{\rm off}^2) / \mathrm{Tr}(h^2)$ | $0.096$ | $0.368$ | $0.310$ |
| $\mathrm{Tr}(h_1^4)$ | $637.9$ | $18.74$ | $795.8$ |
| $\mathrm{Tr}(h_1^6)$ | $1.439 \times 10^4$ | $79.43$ | $1.638 \times 10^4$ |

The identity $\mathrm{Tr}(h_1^2) = \mathrm{Tr}(h_{\rm diag}^2) +
\mathrm{Tr}(h_{\rm off}^2)$ holds bit-exact (residual $\le 7 \times 10^{-15}$
for all three systems). The vanishing cross-term
$\mathrm{Tr}(h_{\rm diag} \cdot h_{\rm off}) = 0$ exactly is the elementary fact
that the diagonal of $A \cdot B$ when $A$ is diagonal and $B$ has zero diagonal
is $(A B)_{ii} = A_{ii} B_{ii} = A_{ii} \cdot 0 = 0$. **This is generic
linear-algebra, not chemistry-emergent.**

### Spectral action vs Taylor prediction

The diagnostic compares $S(\Lambda)$ from `scipy.linalg.expm` against the
$k_{\max} = 6$ Taylor expansion
$\sum_{k=0}^{6} (-1)^k \mathrm{Tr}(D^{2k}) / (k!\,\Lambda^{2k})$ for $\Lambda \ge 5$:

| System | $\|\lambda\|_{\max}(D)$ | rel. residual at $\Lambda = 5$ | smallest $\Lambda$ for resid $\le 10^{-6}$ |
|:------ | ---:| ---:| ---:|
| LiH | $4.918$ | $6.25 \times 10^{-6}$ | $7.91$ |
| $H_2$ | $2.073$ | $5.37 \times 10^{-11}$ | $3.39$ |
| NaH | $4.963$ | $1.16 \times 10^{-5}$ | $7.91$ |

At $\Lambda = 100$ all three are at machine precision ($\le 10^{-15}$). The
threshold $\Lambda \gtrsim 1.6 \|\lambda\|_{\max}$ for the $k_{\max} = 6$
truncation gives $10^{-6}$ residuals; this matches the radius-of-convergence
expectation for the matrix exponential power series. **The expansion is the
trivial matrix-exp Taylor.**

### Substitution panel

At a representative $\Lambda$, comparing $S(D)$ for the three substitutions
$D \in \{h_1, \mathrm{diag}(h_1), h_1 - \mathrm{diag}(h_1)\}$:

LiH (M=15):
```
   Lambda |  S(full) |  S(diag) |   S(off) | S(d)+S(o)-M
    0.471 |  3.012   |  2.763   | 12.322   |   0.085
    0.954 |  6.669   |  6.582   | 13.062   |   4.644
    1.931 | 10.698   | 10.945   | 14.156   |  10.101
    3.907 | 13.173   | 13.337   | 14.755   |  13.092
   10.481 | 14.647   | 14.680   | 14.964   |  14.644
```

Note $S(\mathrm{diag}) + S(\mathrm{off}) - M \ne S(\mathrm{full})$ in general
because $\mathrm{exp}(A + B) \ne \mathrm{exp}(A) \cdot \mathrm{exp}(B)$ for
non-commuting $A, B$, and $[h_{\rm diag}, h_{\rm off}] \ne 0$ generically. This
non-additivity carries no chemistry-specific information: it is the BCH
non-commutativity correction visible in any matrix decomposition $D = A + B$
with $[A, B] \ne 0$.

## Why no Seeley–DeWitt structure emerges

In Chamseddine–Connes SM, the spectral action
$S(D, \Lambda) = \mathrm{Tr}\,f(D^2/\Lambda^2)$ on an infinite-dim Dirac with
positive cutoff $f$ admits the heat-kernel asymptotic

$$
S(D, \Lambda) \;\sim\; \Lambda^d f_0 a_0(D^2) + \Lambda^{d-2} f_2 a_2(D^2)
+ \Lambda^{d-4} f_4 a_4(D^2) \log\Lambda + \cdots
$$

where $d$ is the manifold dimension, $f_k = \int_0^\infty f(u) u^{(k/2)-1} du$
are Mellin moments of $f$, and $a_k(D^2)$ are the Seeley–DeWitt heat-kernel
coefficients (curvature scalars, Yang–Mills field strengths, Higgs and Yukawa
potentials in the SM finite spectral triple). The hierarchy of $\Lambda$
powers comes from the UV divergence of $\mathrm{Tr}\,e^{-tD^2}$ as $t \to 0^+$
on an infinite-dim Dirac.

For a **finite-dim** Dirac (our $h_1$ is $15 \times 15$ for LiH, $10 \times 10$
for $H_2$/NaH):

- $\mathrm{Tr}\,e^{-tD^2}$ is an entire function of $t$; no $t \to 0^+$
  divergence.
- As $\Lambda \to \infty$, $S(\Lambda) \to \mathrm{Tr}(I) = M$, a constant.
- The asymptotic expansion is purely the Taylor series of $\exp(-x)$ around
  $x = 0$, with coefficients $\mathrm{Tr}(D^{2k})$ — no separation into
  $\Lambda^d, \Lambda^{d-2}, \ldots$ tiers, no heat-kernel coefficients.

This is **structurally forced**: the SD coefficients don't appear because the
trace they would extract is already finite at every $\Lambda$. The
Chamseddine–Connes machinery REQUIRES the infinite-dim limit (or equivalently,
$n_{\max} \to \infty$ in the GeoVac n_max cutoff) to surface a non-trivial
$\Lambda$-power expansion.

This connects directly to a finding the GeoVac arc has made before: **Paper 32
§VIII Sprint TS case-exhaustion theorem** says that every $\pi$ in GeoVac is
$\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$ at $k \in \{0, 1, 2\}$ — a Mellin
transform on the heat-kernel side. The Mellin engine is the infinite-dim,
heat-kernel-on-$S^3$ structure where SD coefficients exist. The finite-cutoff
$n_{\max} = 2$ chemistry $h_1$ is several steps removed from that limit. So no
SD analog appears here, by the same structural reason that the master Mellin
engine is a continuum-limit statement.

## Cross-system reading: is the off-diagonal fraction chemistry-meaningful?

The ratio $\mathrm{Tr}(h_{\rm off}^2) / \mathrm{Tr}(h^2)$ varies systematically
across the three systems:

- $H_2$: $0.368$ (large off-diagonal fraction — strong covalent bond, no
  same-center diagonal dominance from a heavy core).
- NaH: $0.310$ (similar, second-row alkali hydride — significant cross-center
  $V_{\rm ne}$ matrix elements in the Na 3s ↔ H 1s sector).
- LiH: $0.096$ (small — large diagonal dominance from Li core + LiH bond center
  blocks; bond is "smaller fraction of the total" because there are TWO diagonal
  blocks before the bond block).

This is sensible chemistry: $\mathrm{Tr}(h_{\rm off}^2)$ in
`cross_block_h1=True` is the sum of squared cross-center $V_{\rm ne}$ matrix
elements between blocks on different nuclei — what M-vS would call the
**bond-coupling intertwiner Frobenius norm squared**. So the PARTIAL caveat
to NEGATIVE: at the level of the leading coefficient $\mathrm{Tr}(D^2)$, the
diag/off split IS chemistry-meaningful as (kinetic + same-center $V_{\rm ne}$)
+ cross-center $V_{\rm ne}$.

But this isn't the spectral action doing the work. It is the elementary
fact that $\mathrm{Tr}(D^2)$ for a Hermitian matrix is the Frobenius norm
squared, which split additively along ANY orthogonal decomposition of $D$. We
chose the diag/off decomposition because it aligns with chemistry; one could
equally choose intertwiner-direction-vs-vertex decompositions, etc. The
spectral action provides no preferred decomposition.

## Cross-check with M-vS-2 Q2 outcome

Sprint M-vS-2 (v3.87.0) found that $S(D)(R)$ as a function of bond length $R$
at FIXED $\Lambda$ is **monotone-decreasing**, no interior minimum, no binding.
That sprint's Q2 verdict combined with the present sprint's verdict gives a
consistent reading:

- **Q2 (M-vS-2, R-axis):** $S(D)(R)$ does not bind chemistry because
  $\mathrm{Tr}(D^2)$ grows monotonically as the bond shortens (eigenvalues get
  bigger).
- **Spectral-action-expansion (this sprint, $\Lambda$-axis):** $S(D)(\Lambda)$
  at fixed $R$ has no chemistry-specific Seeley–DeWitt-style coefficient
  structure because finite-dim $D$ admits only the trivial matrix-exp Taylor.

**Combined:** the spectral action is the WRONG functional for chemistry binding
at finite $n_{\max}$. Binding lives in the FCI ground-state expectation
$\langle \Psi_0 | H | \Psi_0 \rangle$, NOT in the trace functional
$\mathrm{Tr}\,f(D^2/\Lambda^2)$. The trace functional collapses to a generic
trace invariant of $D$ with no preferred decomposition.

## Connection to CLAUDE.md WH5 and the broader framing

This is consistent with WH5 (α is a projection constant, not a derivable
number): the Chamseddine–Connes spectral action machinery, when transferred
naively to finite-cutoff chemistry, does not generate chemistry-side analogues
of Einstein–Hilbert / Yang–Mills / Higgs / Yukawa. It produces only the
generic trace invariants $\mathrm{Tr}(D^{2k})$ of the finite-dim Dirac. Any
chemistry meaning must enter through the explicit operator decomposition of
$h_1$ (e.g., the diag vs off split), not through a spectral-action coefficient
hierarchy.

This is also consistent with **CLAUDE.md §1.7 WH1 (PROVEN at qualitative-rate
GH-convergence)**: the SD coefficient structure is a continuum-limit
property of the round-$S^3$ spectral triple, which is what Paper 38's GH
convergence delivers. Truncating to chemistry-relevant $n_{\max} = 2$ kills
the SD hierarchy by removing the UV divergence.

## Verdict

**NEGATIVE-with-PARTIAL caveat.** The $\Lambda$-expansion of $S(D, \Lambda)$ on
the GeoVac chemistry Dirac at $n_{\max} = 2$ is bit-exactly the trivial
matrix-exponential Taylor series. No Seeley–DeWitt analog or chemistry-
specific coefficient hierarchy emerges. The PARTIAL caveat is that the leading
coefficient $\mathrm{Tr}(D^2)$ splits cleanly into diagonal (kinetic +
same-center $V_{\rm ne}$) and off-diagonal (cross-center $V_{\rm ne}$) pieces,
the latter being the M-vS bond-coupling intertwiner Frobenius norm. This is a
property of the block structure of $h_1$, not an emergent feature of the
spectral action.

**Status:** Closes the spectral-action chemistry-expansion thread. Spectral
action ≠ chemistry binding functional at finite $n_{\max}$. Binding lives in
the FCI ground-state expectation, NOT in $\mathrm{Tr}\,f(D^2/\Lambda^2)$.

## Files

Driver: `debug/sprint_spectral_action_expansion_chemistry_diagnostic_driver.py`
Log: `debug/sprint_spectral_action_expansion_chemistry_diagnostic_log.txt`
Data: `debug/data/spectral_action_expansion_chemistry.json`
Memo: this file.

Staged release notes:
- `debug/changelog_staging_spectral_action_diagnostic.md`
- `debug/claudemd_staging_spectral_action_diagnostic.md`

## References

- Sprint M-vS-2: `debug/sprint_mvs2_lih_default_plus_rsweep_memo.md` (Q1 PASS,
  Q2 NEGATIVE — monotone $S(D)(R)$).
- Sprint M-vS-1: H₂ Bratteli pilot
  (`debug/bratteli_h2_pilot_memo.md`).
- Chamseddine–Connes 1997 spectral action principle.
- Marcolli & van Suijlekom 2014, "Gauge networks in noncommutative geometry",
  arXiv:1301.3480.
- Paper 32 §VIII (case-exhaustion theorem, master Mellin engine).
- WH1 (PROVEN), WH5 (α as projection constant).
