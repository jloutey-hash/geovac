# Sprint: Lorentzian-propinquity rescue probe — Toeplitz under K⁺ compression (2026-06-18)

**Trigger.** PI question during /qa group1 Bite B-3 ("all these rescues has me
wondering — can we rescue the Lorentzian propinquity?"), after the p39 (C₃→√2)
and p53 (disk→plane) rescues. PI: "run the probe now; decide recording after
seeing the result."

## Question

Paper 45 was descoped on the **K⁺ annihilation theorem**: the Krein-self-adjoint
Lorentzian Dirac `D_L = i(J_s⊗D_t + D_GV⊗I_t)` has `{J, D_GV⊗I_t}=0`, so the
positive-compression `P₊ D_L P₊` annihilates the spatial Dirac and the restricted
Lipschitz seminorm `‖[P₊ D_L P₊, P₊ a P₊]‖` was found `= 0` for every multiplier
`a` — kernel = entire operator system. But P45 tested only the **spatial** GeoVac
multipliers. WH7 (2026-06-10) showed the *parallel* P45 "time is Lipschitz-invisible"
result was an artifact of the **momentum-diagonal** algebra `g(D_t)`, repaired by
the **Toeplitz** algebra `S_q = P_K M_{e_q} P_K` (`L(S_q) = ω_q` exactly).

**Decisive probe:** does the K⁺ compression annihilate the *Toeplitz temporal*
multipliers too, or do they survive (→ a nondegenerate seminorm on the compressed
Krein system)?

## Result (verified on the real P45 machinery, cells (2,3) and (3,5))

`krein = CompactTemporalKreinSpace`, `D_L = lorentzian_dirac_compact_matrix`,
`P₊ = ½(I+J)`, multiplier `a = I_spatial ⊗ S_q`:

| Multiplier class | `s_full` | `s_restr` = ‖[P₊D_L P₊, P₊aP₊]‖ |
|:--|:--|:--|
| spatial (P45) | 0.23 / 7.20 | **0.00** (annihilation reproduced) |
| Toeplitz temporal `I_s⊗S_q` | ω_q | **ω_q = 2πq/T, bit-exact** (1.000; 2.000) |

So the K⁺ compression does **NOT** annihilate the Toeplitz temporal seminorm —
it survives equal to the continuum value. P45's annihilation is
**spatial-multiplier-only**, not system-wide (third instance of the
"annihilation = tested-multiplier-selection artifact" pattern, after WH7's
momentum-diagonal and p53's wrong-carrier).

## Structural reason (hand-derivation, confirmed numerically)

`J = γ⁰_spatial ⊗ I_{N_t}` acts **only on space**, identity on time ⇒
`P₊ = P₊_s ⊗ I_t`. Hence:
- spatial: `P₊_s D_GV P₊_s = 0` (since `{J_s, D_GV}=0`) → annihilated;
- temporal: `P₊ D_L P₊ = i·P₊_s ⊗ D_t` (≠0, since `[J,D_t]=0`), and
  `P₊(I_s⊗S_q)P₊ = P₊_s⊗S_q`, giving
  `s_restr = ‖i P₊_s ⊗ [D_t,S_q]‖ = ω_q·‖P₊_s‖·‖S_q‖ = ω_q`.

The K⁺ compression kills **exactly the sector where the Lorentzian signature acts**
(spatial, `{J,D}=0`) and spares **exactly the sector where it does not** (temporal,
`[J,D]=0`).

## v2 probe (temporal Krein involution — PI: "scope it now")

To test whether the signature-triviality is a v1 artifact (J = γ⁰_s ⊗ I_t is
identity on time), re-ran with a genuine temporal Krein involution
**J = γ⁰_s ⊗ sign(D_t)** — the positive/negative-**frequency** (energy) splitting,
i.e. the actual Wick/Feynman Krein structure for time (sign(0)=+1 for the
rest mode). This was chosen over time-reversal `R` (k↔−k), which makes the
temporal Dirac ALSO anticommute with J (`{J, J_s⊗D_t}=0`) and would annihilate
everything.

Legitimacy checks (cells (2,3),(3,5),(2,7)): `J_t²=I` exact; spatial annihilation
preserved (`{J,D_GV}=0` exact); and `[J,D_L]≠0` AND `{J,D_L}≠0` (5/2, 7/4, 5/6) —
so v2 J is a genuine Krein structure that acts nontrivially on time and is NOT
v1's commuting case.

**Result: s_restr_v2(S_q) = ω_q bit-exactly at every cell and every q
(ratio s_restr/ω = 1.0000).** Even though sign(D_t) does not commute with the
Toeplitz shift S_q (the shift crosses the frequency-sign boundary, so the
compressed multiplier P₊ a P₊ is genuinely rearranged), the surviving seminorm
is unchanged from the Euclidean value.

So the indefinite signature is **metrically invisible to the surviving Lipschitz
seminorm regardless of whether the Krein involution acts on time** — the
signature-triviality is ROBUST, not a v1 artifact.

## Verdict (the rescue-vs-convention fork, resolved — robustly)

- **Rescuable from annihilation: YES.** A nondegenerate, convergent Lipschitz
  seminorm exists on the K⁺-compressed Krein operator system (Toeplitz temporal,
  `= ω_q`). The seminorm is not identically zero; P45's annihilation theorem is
  spatial-multiplier-only.
- **Genuinely *indefinite-signature* geometry: NO — and robustly so.** The
  surviving metric is the **Euclidean temporal** Lipschitz constant. v2 shows
  this is not because J happened to be identity on time:\ a genuine temporal Wick
  involution `J_t = sign(D_t)` leaves the surviving seminorm bit-exactly `ω_q`
  too. No choice of temporal Krein involution injects signature content into the
  propinquity Lipschitz seminorm. **The Lorentzian signature is convention
  relative to this metric.**
- **WH7 impact:** the visibility leg's "weakens-to-convention" lean is confirmed
  AND sharpened — the indefinite signature itself is metrically invisible, robust
  under the temporal Krein structure. (Does not touch WH7 leg (i):\ compactness →
  discreteness + π-injection, which is the load-bearing content.)
- **Remaining genuine open:** whether a *non-Lipschitz* Krein metric device (a
  J-indefinite "length" rather than the operator-norm Lipschitz seminorm) could
  carry signature content. That is a different metric primitive, not the
  propinquity seminorm; a separate research question, not a rescue of the
  propinquity object.

## Honest scope

- The numbers are bit-exact at the two existing P45 panel cells; the structural
  derivation is general (`J = γ⁰_s ⊗ I_t` is the construction's defining choice).
- This **sharpens** P45's annihilation theorem (spatial-only) and **confirms**
  WH7's primary-falsifier "weakens-to-convention" lean with an exact
  space/time sector decomposition. It does **not** revive a Lorentzian-metric
  convergence claim for the papers.
- Disposition (paper Status notes / CLAUDE.md §1.7 WH7 / Paper 45 / 47) deferred
  to PI per "decide after seeing the result".

Driver (reproducible inline; freeze as `tests/test_lorentzian_toeplitz_kplus.py`
if promoted): builds `temporal_shift(N_t,q)` in the `fourier_momentum_grid` basis,
injects `I_spatial ⊗ S_q`, compares `s_full` / `s_restr` to `ω_q = 2πq/T`.
