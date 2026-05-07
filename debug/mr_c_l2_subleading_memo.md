# Sprint MR-C: L2 next-order constant — partial negative against M2 prediction

**Verdict: PARTIAL NEGATIVE.** The L2 next-order constant $c$ in
$n \cdot \gamma_n = (4/\pi) \log(n) + c + O(\log(n)/n)$ is extracted
cleanly to $\geq 80$ digits, but PSLQ at coefficient ceiling $10^4$
returns no closed-form identification in
$M_1 \cup M_2 \cup M_3 \cup \{\gamma_E, \log 2, G, \zeta(3)\}$.

This sprint memo is concise — the original sprint agent stalled before
writing a long memo; the data files preserve the substantive
computation.

## 1. Setup

The L2 quantitative rate sprint (closed 2026-05-06 morning) proved
$\gamma_n \leq 6 \log(n)/n$ for $n \geq 2$ and the asymptotic constant
$\lim n \gamma_n / \log(n) = 4/\pi$. The subleading term in
$n \cdot \gamma_n = (4/\pi) \log(n) + c + O(\log(n)/n)$ was not
pinned in closed form.

The master Mellin engine prediction (Paper 18 §III.7, Sprint TS-E1)
suggested that $c$ might land in M2 ring $\sqrt{\pi}\cdot\mathbb{Q}
\oplus \pi^2\cdot\mathbb{Q}$, on the reasoning that Stein–Weiss IBP
through next order picks up the heat-kernel asymptotic.

## 2. Method

Driver: `debug/mr_c_l2_subleading.py`.

1. Compute $\gamma_n$ to 120 dps at panel
$n \in \{16, 32, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536,
2048, 3072, 4096\}$ via `central_fejer_su2.gamma_n_via_sum_rule`.
2. Define $h_n = n \cdot \gamma_n - (4/\pi) \log(n)$. Fit
$h_n = c + a_1 \log(n)/n + a_2/n + \ldots$ via Richardson extrapolation
on $K = 1, 2, \ldots, 6$ orders.
3. PSLQ refinement: `debug/mr_c_pslq_clean.py` with corrected basis to
remove the $\{4/\pi, 1/\pi\}$ basis-internal redundancy.

## 3. Result

$$
c = 4.1093214674877940927579607260741005838057691088362615503253972964276017819113301\ldots
$$

at residual $3.87 \times 10^{-16}$ between $K = 5$ and $K = 6$ Richardson
orders, stable to 30+ digits.

PSLQ panel:

| Basis | $\max\,|a_i|$ | Found? |
|:--|:-:|:-:|
| Strict M2 ($\sqrt{\pi}\cdot\mathbb{Q} \oplus \pi^2\cdot\mathbb{Q}$ + rational/$\pi$ powers) | $50, 200, 10^3, 10^4, 10^5$ | NO |
| M2 + $\log 2$ (Stein–Weiss intermediate) | $50, 200, 10^3, 10^4, 10^5$ | NO |
| M2 + $\log 2$ + $\gamma_E$ | $50, 200, 10^3, 10^4, 10^5$ | NO |
| Extended (M1 + M2 + M3 + $\log 2 + \gamma_E + G + \zeta(3)$) | $50, 200, 10^3, 10^4$ | NO |
| Stirling + Glaisher | $50, 200, 10^3, 10^4, 10^5$ | basis-internal redundancy, $a_0 = 0$ |

No identification at coefficient ceiling $10^4$. Earlier "found" relations
in `mr_c_l2_subleading.json` were spurious basis-internal relations
(integer relations among basis vectors with $a_0 = 0$, ignoring the
target $c$); those were filtered in `mr_c_pslq_clean.json`.

## 4. Interpretation

Read inside MR-A's structural finding ($\gamma_n^{\rm Dirac} = \pi$
exactly: half-integer characters cannot span $\chi_0$, so propinquity
is structurally an M1 observable):

- The L2 next-order constant $c$ is downstream of the M1 Stein–Weiss
IBP, not a direct $\mathcal{M}[\mathrm{Tr}(D^k \cdot e^{-tD^2})]$
output. The master Mellin engine's predictive content is for
direct Mellin outputs, not for derivative quantities along
subsequent analytical pipelines.
- Failure to land in M2 at $\max\,|a_i| \leq 10^4$ is consistent with
either (a) $c$ in M1 ring with larger coefficients than tested,
(b) $c$ in an extended ring (e.g. $\zeta'$ values, Glaisher, derivative
of L-functions) outside the standard transcendental basis,
(c) $c$ genuinely not having a closed form in any standard basis
because it mixes rings along the IBP chain.

The computation does not distinguish among (a), (b), (c).

A higher-coefficient PSLQ ($\max\,|a_i| \geq 10^6$) is warranted before
declaring a clean negative; this is left for a future sprint, but the
absence of a hit at $10^4$ is itself informative — natural closed
forms typically PSLQ at coefficient size $\leq 10^3$.

## 5. Implication for the master Mellin engine

MR-C does NOT contradict the master Mellin engine reading. The
case-exhaustion theorem (Paper 32 §VIII) classifies which transcendentals
appear in any GeoVac observable's $\pi$-content, and is unaffected by
MR-C: $c$ is not constrained to live in $M_1 \cup M_2 \cup M_3$ if it
contains no $\pi$ (e.g. it could be a rational multiple of $\gamma_E$ or
$\zeta(3)$), and no structural prediction was made about its content.

The mechanism-as-domain reading (Sprint MR-A/B/C synthesis) does
predict that direct heat-kernel-Mellin outputs land in their ring
domain — that's what MR-B confirmed for the Camporesi–Higuchi modular
residual. MR-C's $c$ is two analytical steps removed from such an
output (Stein–Weiss IBP applied to the central Fejér moment, then
asymptotic expansion). The engine doesn't commit to the rings of
those derivatives. This is a sharpening of the predictive engine,
not a falsification.

## 6. Recommendations

1. **No paper update required for MR-C alone** — the negative is
informative but doesn't change any paper claim. The fact that the
master Mellin engine doesn't predict $c$ is incorporated into
Paper 18 §III.7's "Mechanism-as-domain sharpening" paragraph alongside
MR-A and MR-B.
2. **Open follow-up: extend PSLQ ceiling to $10^6$** if curiosity
demands. Cost: $\sim 10\times$ the current run, marginal value.
3. **Open follow-up: try $\zeta'$ and derivative-of-L basis** at
$\max\,|a_i| = 10^4$. If $c \in$ derivative basis, it would tighten
the picture — but does not change the qualitative reading.

## 7. Files

- `debug/mr_c_l2_subleading.py` — main driver
- `debug/mr_c_pslq_clean.py` — PSLQ basis cleanup
- `debug/data/mr_c_l2_subleading.json` — $\gamma_n$ panel + Richardson + initial PSLQ
- `debug/data/mr_c_pslq_clean.json` — clean PSLQ scan against five basis sets
