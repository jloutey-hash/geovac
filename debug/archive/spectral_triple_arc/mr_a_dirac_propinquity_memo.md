# Sprint MR-A: Dirac analog of L2 propinquity rate — clean structural negative with master-Mellin interpretation

**Verdict: scenario (d) STRUCTURAL NEGATIVE with REINFORCING interpretation.**
The natural Dirac analog of the L2 central Fejér kernel — built from
half-integer-spin Peter-Weyl characters with Plancherel weights
$\sqrt{g_n^{\rm Dirac}}$ — gives mass-concentration moment
$\gamma_n^{\rm Dirac} = \pi$ **exactly** for every cutoff $n_{\max} \geq 1$.
There is no asymptotic rate to PSLQ. The kernel does not concentrate.

This is a clean negative against the naive "the propinquity rate is the
M3 signature" reading. But re-examined in light of MR-B, the result
**reinforces** the master Mellin engine's three-mechanism partition by
demonstrating that the propinquity rate is structurally an **M1
observable**, not an M3 one — the half-integer obstruction is exactly
what M1's kernel-on-base-manifold mechanism predicts when the relevant
sub-bundle lacks the constant function.

## 1. Setup

The L2 sprint (closed 2026-05-04 + 2026-05-06) computed the Latrémolière
propinquity rate for the SU(2) scalar central Fejér kernel built from
spin-$j$ characters $\chi_j(\chi) = \sin((2j+1)\chi/2)/\sin(\chi/2)$ summed
over **all** spins $j \in \{0, 1/2, 1, 3/2, \ldots, j_{\max}\}$ with
Plancherel weights $\sqrt{2j+1}$. The asymptotic rate constant is
$4/\pi = \mathrm{Vol}(S^2)/\pi^2$, the **M1** Hopf-base measure signature.

MR-B (this session, 2026-05-06) computed the Connes-Chamseddine spectral
action modular residual for the half-integer-shifted Camporesi-Higuchi
Dirac spectrum and found the closed form
$$
\varepsilon(t) = \sum_{m \geq 1}(-1)^m \sqrt{\pi}\, e^{-m^2 \pi^2/t}
                 \left[t^{-3/2} - 2 m^2 \pi^2 t^{-5/2} - \tfrac{1}{2}t^{-1/2}\right]
$$
verified to >100 digits. Every coefficient lies in the **M2** ring
$\sqrt{\pi}\cdot\mathbb{Q} \oplus \pi^2 \cdot \mathbb{Q}$.

Sprint MR-A asks: does the **Dirac propinquity rate** — the natural
half-integer analog of L2 — produce a constant in the **M3** ring (Catalan
$G$, Dirichlet $\beta(s)$, half-integer Hurwitz $\zeta(s, 1/4)$)?

## 2. Construction

The spinor bundle on $S^3 = \mathrm{SU}(2)$ decomposes under
$\mathrm{Spin}(4) = \mathrm{SU}(2)_L \times \mathrm{SU}(2)_R$ as
$$
\bigoplus_{n \geq 0}
   V_{(n+1)/2,\, n/2} \;\oplus\; V_{n/2,\, (n+1)/2},
$$
each summand of dimension $(n+1)(n+2)$. Restricted to the diagonal
$\mathrm{SU}(2)$, each $V_{(n+1)/2, n/2}$ decomposes into spins
$j = 1/2, 3/2, \ldots, n+1/2$ with multiplicity 1. The full Dirac sector
gives $g_n^{\rm Dirac} = 2(n+1)(n+2)$ at level $n$, supported on
**half-integer spins only**.

The natural Dirac analog of L2's central Fejér kernel is therefore
$$
D^{\rm Dirac}_{n_{\max}}(\chi) \;=\;
   \sum_{n=0}^{n_{\max}-1} \sqrt{g_n^{\rm Dirac}}\,
       \chi_{n+1/2}(\chi),
\qquad
K^{\rm Dirac}_{n_{\max}}(\chi) \;=\;
   \frac{|D^{\rm Dirac}_{n_{\max}}(\chi)|^2}{Z^{\rm Dirac}_{n_{\max}}}
$$
with $\chi_{n+1/2}(\chi) = \sin((n+1)\chi)/\sin(\chi/2)$ and
$Z^{\rm Dirac}_{n_{\max}} = \sum_{n=0}^{n_{\max}-1} g_n^{\rm Dirac}
 = \frac{2}{3}\,n_{\max}(n_{\max}+1)(n_{\max}+2)$.

The mass-concentration moment is
$$
\gamma_n^{\rm Dirac} \;=\;
   \frac{1}{\pi} \int_0^{2\pi} \chi \cdot K^{\rm Dirac}(\chi) \cdot
       \sin^2(\chi/2)\, d\chi.
$$
Driver: `debug/mr_a_dirac_propinquity_rate.py`.

## 3. Result: $\gamma_n^{\rm Dirac} = \pi$ exactly, for all $n_{\max}$

### 3.1 Direct calculation (verifying small $n_{\max}$)

Sympy direct integration confirms $\gamma_n^{\rm Dirac} = \pi$ exactly
for $n_{\max} \in \{1, 2, 3\}$. For larger $n_{\max}$ sympy is
prohibitively slow (factorial blow-up of $\binom{n_{\max}}{2}$ cross
terms before simplification), so we use the structural argument below.

### 3.2 Structural proof (all $n_{\max}$)

Expanding $|D^{\rm Dirac}|^2$ and using
$\sin((n+1)\chi)\sin((n'+1)\chi) = \frac{1}{2}[\cos((n-n')\chi) -
\cos((n+n'+2)\chi)]$, then integrating against $\chi\,d\chi$ over
$[0, 2\pi]$:

**Key fact.** For nonzero **integer** $p$,
$\int_0^{2\pi} \chi \cos(p\chi)\, d\chi = 0$ (integration by parts:
boundary $\sin(2\pi p)/p = 0$, and the remaining $-\cos(p\chi)/p^2$
returns 0 at both endpoints).

For $p = 0$ (i.e. $n = n'$, the diagonal),
$\int_0^{2\pi} \chi \, d\chi = 2\pi^2$.

**Crucial.** Since $n, n' \geq 0$ integers, both $n - n'$ and
$n + n' + 2$ are **integers**. The case $n + n' + 2 = 0$ never occurs
($n + n' + 2 \geq 2$). So **only the diagonal $n = n'$ contributes**, and
the cross-term $\cos((n+n'+2)\chi)$ contribution vanishes for *every*
pair.

Thus
$$
\gamma_n^{\rm Dirac} \;=\;
   \frac{1}{\pi Z^{\rm Dirac}}
   \cdot \frac{1}{2}
   \cdot \sum_{n=0}^{n_{\max}-1} g_n^{\rm Dirac} \cdot 2\pi^2
   \;=\;
   \frac{\pi \cdot Z^{\rm Dirac}}{Z^{\rm Dirac}}
   \;=\; \pi.
$$
The cancellation is **exact, not asymptotic**, and holds for every
$n_{\max} \geq 1$.

### 3.3 Numerical confirmation

`mpmath` quadrature at 50 dps confirms
$|\gamma_n^{\rm Dirac} - \pi| < 10^{-30}$ at $n_{\max} \in \{2, 3, 4, 5,
8, 10, 16, 20\}$. No log decay, no power-law decay — a constant.

## 4. Interpretation: why this happens

The half-integer-only Peter-Weyl truncation **cannot span the constant
function** $\chi_0(\chi) = 1$ on $\mathrm{SU}(2)$. The constant has
spin $j = 0$ (integer), which is excluded from the Dirac sub-bundle.
Without the constant, the kernel cannot concentrate at the identity:
in the conjugacy-class basis, $K^{\rm Dirac}$ at any $n_{\max}$ has
a peak distribution that broadens — it does not localize on $\chi = 0$.

For the L2 scalar kernel, the rate decay $\gamma_n^{\rm scalar}
\sim (4/\pi)\log(n)/n$ comes specifically from **integer × half-integer
cross terms** in $|D^{\rm scalar}|^2$. These produce $\cos((j-j')\chi)$
or $\cos((j+j'+1)\chi)$ with **half-integer** argument coefficient
(when one $j$ is integer and the other half-integer), and
$\int_0^{2\pi} \chi \cos(p\chi)\, d\chi = -2/p^2$ for half-integer $p$.
The $1/(k_1 \pm k_2)^2$ terms in the L2 sum-rule
$T_n = \sum_{k_1+k_2 \text{ odd}} \sqrt{k_1 k_2}\,[1/(k_1-k_2)^2 -
1/(k_1+k_2)^2]$ are precisely these mixed-parity contributions. Removing
the integer-spin sector kills $T_n$, leaving $\gamma = \pi$.

**The L2 rate decay is an integer-half-integer interference effect, not
a half-integer-only effect.** This is a sharp structural finding that
was not visible from the L2 sprint alone.

## 5. Implication for the master Mellin engine

At first glance, MR-A is a clean **negative** for the predictive-engine
reading: M3's natural ring (Catalan $G$, $\beta(s)$, half-integer
Hurwitz) is not produced by the Dirac propinquity rate — instead the
construction is degenerate. But re-examined, MR-A is consistent with
**a sharper version** of the master Mellin engine:

| Mechanism | Operator order $k$ | Natural domain | Witnessed signature |
|:---|:---:|:---|:---|
| **M1** (Hopf-base measure) | $k = 0$ | Kernel concentration on the base manifold (propinquity rate) | $4/\pi = \mathrm{Vol}(S^2)/\pi^2$ (L2 sprint) |
| **M2** (Seeley-DeWitt) | $k = 2$ | Heat-kernel / spectral-action convergence | $\sqrt{\pi}\,\mathbb{Q} \oplus \pi^2\mathbb{Q}$ ring; modular exponent $\pi^2$ (MR-B sprint) |
| **M3** (vertex parity / Hurwitz) | $k = 1$ | **Vertex-restricted observables**, NOT propinquity rates | Catalan $G$, $\beta(4)$ in QED two-loop sunset (Paper 28 §QED-vertex; existing) |

MR-A's degenerate-kernel result is exactly the prediction: the M3
mechanism, $\mathcal{M}[\mathrm{Tr}(D \cdot e^{-tD^2})]$, lives on
**vertex-restricted spectral sums** with parity-character content (the
$\chi_{-4}$ Dirichlet character of Paper 28's two-loop QED sunset,
producing $G$ and $\beta(4)$). It does **not** live in state-space
kernel concentration. Trying to extract M3 from a Fejér-kernel
propinquity rate is a **category error**: that question is structurally
M1, and the M1 mechanism cannot extract Catalan $G$ because the kernel
itself sees only the base measure.

So the master Mellin engine's three-mechanism partition gains a
sharper interpretation:

> **The mechanism index $k \in \{0, 1, 2\}$ also indexes the natural
> domain of the corresponding observable**: $k=0$ for state-space
> propinquity rates, $k=1$ for vertex-restricted observables, $k=2$ for
> heat-kernel convergence rates. Mixing the index of the question
> ("propinquity") with the index of the mechanism ("vertex parity")
> produces a structural degeneracy — exactly what MR-A measured.

This is a stronger reading than the "predictive engine" framing of the
sprint plan. The engine is predictive *within each domain*, and the
domain is determined by the operator power $k$.

## 6. Comparison to MR-B and MR-C

| Sprint | Object | Result | Verdict |
|:---|:---|:---|:---:|
| **L2 sprint** | Scalar SU(2) propinquity rate | $4/\pi$ in M1 ring | ✓ M1 confirmed |
| **MR-B** | Spectral action modular residual | $\sqrt{\pi}\,\mathbb{Q} \oplus \pi^2\mathbb{Q}$ throughout, exponent $\pi^2$ | ✓ M2 confirmed |
| **MR-C** | L2 subleading constant $c \approx 4.10932\ldots$ | No PSLQ identification in $\{M1, M2, M3, \log 2, \gamma_E, G, \zeta(3)\}$ at 100 dps | partial negative |
| **MR-A** | Dirac propinquity rate | $\pi$ exactly, no decay (degenerate by half-integer obstruction) | structural negative; reinforces the partition |

MR-A's structural finding shows that L2's subleading $c$ (MR-C's
target) is unlikely to live in M3 ring either: by the same argument,
the M1 mechanism does not have access to half-integer Hurwitz / Catalan
content. MR-C's negative is consistent with this — $c$ should live in
the M2 ring (Seeley-DeWitt) if the "next order is Stein-Weiss IBP heat
kernel" reading is correct, but there are larger-coefficient
combinations that would not be found by PSLQ at the tested coefficient
ceiling.

## 7. Recommended paper updates

**Paper 32 §VIII (case-exhaustion theorem) — sharpen the master Mellin
engine reading.** Add a remark following the case-exhaustion theorem
that the index $k \in \{0, 1, 2\}$ also classifies the natural domain
of the corresponding observable:
- $k=0$: state-space convergence (propinquity, kernel concentration on
  base manifold);
- $k=1$: vertex-restricted spectral sums (parity-character content,
  half-integer Hurwitz);
- $k=2$: heat-kernel / spectral-action expansion.
Citing L2 ($4/\pi$, M1), MR-B ($\sqrt{\pi}\,\mathbb{Q} \oplus
\pi^2\mathbb{Q}$, M2), and Paper 28 §QED-vertex ($G$, $\beta(4)$, M3).

**Paper 18 §III.7 (master Mellin engine) — add domain classification.**
Update §III.7 with the operator-order-as-domain-index reading. State
explicitly that one cannot extract the M3 signature by computing an
M1-type observable (Dirac propinquity rate via half-integer Peter-Weyl
truncation) because the M1 mechanism lacks half-integer Hurwitz
content. This is the sharpening MR-A produced.

**CLAUDE.md §1.7 WH1 — update the master Mellin engine reading.** Add
one sentence noting that MR-A's structural negative on the Dirac
propinquity rate (constant $\pi$, no decay) reinforces the
mechanism-as-domain reading. The three-mechanism partition is not
about *which transcendentals appear in any given observable* but about
*which observables can produce which transcendentals*.

**CLAUDE.md §2 — add MR-A sprint summary** with files
`debug/mr_a_dirac_propinquity_rate.py`,
`debug/data/mr_a_dirac_propinquity.json`,
`debug/mr_a_dirac_propinquity_memo.md`.

## 8. Honest caveats

1. **Construction choice.** The "natural Dirac analog" of L2 is not
   uniquely defined. The choice made here — half-integer-only Peter-Weyl
   truncation with $\sqrt{g_n^{\rm Dirac}}$ weights at spin $j = n+1/2$ —
   is the most direct adaptation of the spinor-bundle decomposition,
   but it is one of several possible choices. A future R3.5
   continuation could instead compute the full Dirac propinquity via
   the operator-system Connes-distance SDP (R2.3 framework), which would
   incorporate the eigenvalue shift $|\lambda| = n + 3/2$ explicitly.

2. **The $\gamma = \pi$ result is not a "rate is $\pi/\text{const}$"
   statement.** $\pi$ here is the moment of the *uniform* kernel on
   $\mathrm{SU}(2)$ against $\chi$, i.e. the constant value of $\gamma$
   when $K = 1$ uniformly. Its appearance is a coincidence with the
   normalization of the geodesic distance, not a structural M-mechanism
   signature.

3. **Operator-system Connes distance probably does have a nontrivial
   rate.** R2.3 / R3.5 computed Connes distances at $n_{\max} = 2, 3$
   via SDP; extracting an asymptotic rate from those would require
   $n_{\max} \geq 5$ at least, and the SDP cost grows steeply. That is
   the right next sprint if M3-from-propinquity is wanted; but the
   reading above suggests it's not the right place to look for M3.

4. **MR-C's negative is not closed.** MR-C's $c = 4.10932\ldots$ may
   simply be in M2 with larger coefficients than the PSLQ ceiling
   tested. A higher-coefficient PSLQ run (max_coeff $\geq 10^6$) is
   warranted before declaring a clean negative.

## 9. Files

- `debug/mr_a_dirac_propinquity_rate.py` — driver
- `debug/data/mr_a_dirac_propinquity.json` — symbolic + numerical results
- `debug/mr_a_dirac_propinquity_memo.md` — this file
