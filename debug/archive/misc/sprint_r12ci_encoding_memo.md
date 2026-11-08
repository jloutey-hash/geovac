# Sprint memo: R12-CI encoding cost — overlap conditioning kappa(S) of the geminal NOQE

**Date:** 2026-08-23  **Scope:** debug/ only; no production/paper/CLAUDE.md edits.
**Question:** R12-CI reaches He 0.45–0.80 mHa at a qubit-proxy of ~6–8 (compact, accurate;
`sprint_ctf12_poc_memo.md`). But its geminal `G = e^{-g r12} R_ref(r1) R_ref(r2)` is a
*non-orthogonal correlated basis function*, so encoding the R12-CI eigenproblem as a quantum
algorithm is a **non-orthogonal quantum eigensolver (NOQE) / generalized eigenproblem** whose
cost is governed by the overlap-matrix conditioning `kappa(S)`. Does the function-count
compactness survive once `kappa(S)` is priced in?

## Method
Reused the existing R12-CI integral machinery **READ-ONLY** (imported `debug/ctf12_r12ci_he.py`
as an integral library; not modified). Same fixed diffuse grid as its `run_sweep` (Ng=600,
nx=128, r_max=42/1.30) so energies reproduce the recorded best-`(k,gamma)` operating points
bit-for-bit (e.g. ns=3,ngem=1 → E=−2.902926, 0.798 mHa, matching the memo). Assembled the
overlap matrix `M` over the basis of `n_s` s-Sturmian orbital-pairs + `n_gem` geminals, and
computed the **normalized** conditioning `kappa2(S)` of the unit-diagonal correlation matrix
`S_ij = M_ij / sqrt(M_ii M_jj)` — the physically-relevant object because NOQE prepares
*normalized* quantum states (`S_ii = 1`). (Raw unnormalized kappa also reported; ~3.5–4.6×
larger, but carries the arbitrary per-function norms, so normalized is primary.)

Driver: `debug/r12ci_encoding_cost.py` → `debug/data/r12ci_encoding_cost.json`.

## Result 1 — kappa(S) at the accurate operating point (the headline table)

| n_s | n_gem | nbf | q=2·n_s | E (Ha) | err (mHa) | **kappa_norm(S)** | kappa_raw | lambda_min(S) |
|----:|------:|----:|-------:|--------:|----------:|------------------:|----------:|-------------:|
| 3 | 1 | 7 | 6 | −2.902926 | 0.798 | **177** | 631 | 0.0157 |
| 3 | 2 | 8 | 6 | −2.903084 | 0.640 | **251** | 942 | 0.0129 |
| 4 | 1 | 11 | 8 | −2.903250 | 0.474 | **189** | 725 | 0.0159 |
| 4 | 2 | 12 | 8 | −2.903275 | 0.449 | **236** | 1086 | 0.0140 |

**kappa_norm(S) ≈ 180–250 at every accurate operating point — an order of magnitude past the
O(1–10) "cheap NOQE" GO threshold.** The smallest overlap eigenvalue is ~0.013–0.016 (moderate
near-linear-dependence, not yet catastrophic collapse, but well away from orthonormal).

## Result 2 — where kappa comes from: two compounding sources

**(a) The shared-scale Sturmian pair block is *itself* ill-conditioned** (n_gem=0, no geminal):

| n_s | nbf | kappa_norm | lambda_min |
|----:|----:|-----------:|-----------:|
| 1 | 1 | 1.0 | 1.00 |
| 2 | 3 | 9.1 | 0.222 |
| 3 | 6 | 34.4 | 0.076 |
| 4 | 10 | 91.0 | 0.033 |
| 5 | 15 | 197 | 0.016 |
| 6 | 21 | 375 | 0.009 |

The shared-`k` Coulomb-Sturmians are strongly non-orthogonal (guardrail Paper 8 records
⟨χ_1s|χ_2s⟩=−0.47), and their pair products compound it: kappa grows steeply with basis size
(9→34→91→197→375). This is the corpus's already-documented **Sturmian/Löwdin conditioning wall**
(§3 "sparsity-destroying option" S non-PSD; Sturmian-CI 2.8–4.5× 1-norm inflation).

**(b) The geminal is 93% collinear with the reference orbital pair.** Normalized overlaps of
`G` with each orbital pair at ns=3,ngem=1,k=1.7,gamma=0.7:

  `⟨G|(1,1)⟩ = +0.931`,  ⟨G|(1,2)⟩ = −0.374,  ⟨G|(1,3)⟩ = −0.148,  ⟨G|(2,2)⟩ = +0.128, …

Adding that single geminal jumps kappa from **34.4 → 176.9 (≈5×)** at ns=3 — the geminal's
near-linear-dependence on the (1,1) pair is the dominant multiplier. A second geminal adds
another ~1.4× (177→251); a third saturates (251→261) and buys no energy (0.640→0.640 mHa).

## Result 3 — the honest tension: accuracy and conditioning are ANTI-correlated

gamma scan at ns=3, n_gem=1, k=1.7 (the variational optimum sits at the *small*-gamma end):

| gamma | E (Ha) | err (mHa) | kappa_norm |
|------:|--------:|----------:|-----------:|
| 0.40 | −2.903408 | **0.316** | **402** |
| 0.50 | −2.903357 | 0.367 | 283 |
| 0.70 | −2.902926 | 0.798 | 177 |
| 1.10 | −2.901211 | 2.514 | 109 |
| 2.00 | −2.896040 | 7.685 | 73 |
| 3.00 | −2.891018 | 12.707 | **60** |

The variational optimum **pulls toward ill-conditioning**: the more accurate the geminal (small
gamma → milder modification → more collinear with the reference pair), the worse kappa(S). You
cannot buy your way to a well-conditioned overlap: pushing gamma up to tame kappa to ~60 costs an
order of magnitude in accuracy (0.32 → 12.7 mHa, back near the plain-FCI floor). Low-kappa and
high-accuracy are mutually exclusive in this basis.

## Result 4 — thresholding gives NO relief (the NOQE cost is unavoidable here)

Canonical-orthogonalization thresholding (discard overlap eigenvalues below tau) is the standard
NOQE conditioning lever (Baek et al. 2608.12830: eigenvalue error ∝ kappa of the *retained*
block). But at every accurate operating point the smallest normalized eigenvalue is ~0.013–0.016,
so **for any tau in [1e-4, 1e-10] all functions are retained and kappa_retained = kappa_full**
(≈177–251). The ill-conditioning is not a few droppable near-null directions — the small-eigenvalue
direction *is* the correlation/cusp direction the geminal contributes; discarding it throws away the
very correlation energy that took the error from 26 mHa (orbital-only) to <1 mHa. There is no
cheap thresholded sub-block.

## NOQE cost implication

- **Fault-tolerant GEVP route** (Liang et al. 2112.02554): runtime carries `kappa_B` (metric
  condition number) as a *direct multiplicative factor*. Here `kappa_B ≈ 180–250` at the accurate
  point → a ~180–250× overhead vs an orthonormal encoding of the same size, and it *grows* toward
  the accuracy target (gamma→small pushes kappa past 400).
- **Near-term NOQE route** (Baek et al. 2205.09039; improved-scaling 2608.12830): measurement /
  shot cost ∝ kappa of the retained overlap block. Thresholding buys nothing (Result 4), so the
  effective cost is the full kappa ≈ 180–250. The 2608.12830 improvement makes this *linear* in
  kappa rather than worse — real but not catastrophic on a 7–12-function basis in absolute terms.
- **Net:** the compactness advantage that motivated the pivot (R12-CI ≈18× fewer functions than
  plain FCI at matched ~1 mHa; qubit-proxy 6–8) is **roughly offset by the ~200× overlap-conditioning
  multiplier** the non-orthogonal geminal imposes. The function-count win does NOT translate cleanly
  into a quantum-resource win once conditioning is priced in.

## Verdict — **CAUTION** (leaning STOP on the "cheap NOQE" reading)

kappa(S) ≈ 180–250 ≫ 100 at the accurate operating point, driven by (a) the shared-scale Sturmian
pair block (already-documented conditioning wall) and (b) a geminal that is 93% collinear with its
reference pair — and the variational optimum actively *pulls toward* the ill-conditioned corner.
This is the same non-orthogonality/conditioning tension the corpus documents for the Löwdin retrofit
and the Shibuya–Wulfman metric: the geminal's compactness and its non-orthogonality are two faces of
one object. Not a categorical wall — kappa≈200 is a constant (not exponentially growing) on a tiny
basis, so a small NOQE is not impossible — but the pivot's headline selling point (compactness →
cheap encoding) is **not established**: the conditioning cost eats the compactness. If pursued, the
honest framing is "trade ~18× fewer functions for a ~200× conditioned metric," not "compact and
cheap," and any real advantage would have to come from a *better-conditioned* correlated basis
(orthogonalized/optimized geminals, or the isoenergetic reformulation of Paper 60 that removes the
metric for atoms) — not from R12-CI as-is.

## Files
- `debug/r12ci_encoding_cost.py` — kappa(S) probe (imports ctf12_r12ci_he.py read-only).
- `debug/data/r12ci_encoding_cost.json` — all numbers (operating points, retained-kappa vs tau,
  orbital-only trend, geminal trend, gamma scan, collinearity probe).
