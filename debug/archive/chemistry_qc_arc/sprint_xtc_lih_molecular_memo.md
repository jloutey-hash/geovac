# Sprint memo — does the xTC angular-sparsity survival extend from an ATOM to a MOLECULE? (LiH, 4e)

**Date:** 2026-08-23  **Verdict: GO (qualified).** 0 fill-in, 1-norm ≤ plain — the s-reference
collapse *transfers* — **but** the headline atomic l-sparsity is a single-center property that the
two-center Coulomb baseline already lacks, so the molecular "win" is the weaker "no densification +
lower 1-norm," not "a rich Gaunt sparsity preserved."

**Question the atomic run left open.** The p-inclusive atomic GO (`sprint_xtc_pinclusive_memo.md`)
showed the xTC-contracted 2-body inherits GeoVac's angular Gaunt sparsity EXACTLY (0 fill-in) — but
only because contracting a correlator leg over a single-center **s** reference forces its leg
multipole L′=0 (monopole), collapsing the four_Y vertex to Coulomb's single-multipole support. A
molecule breaks that premise: the reference σ MO is a two-center object with its own multipole tower
about any single center. Does the collapse survive, or does the two-center reference fill in the zeros?

**Build (all in `debug/`).** Minimal two-center LiH: Li@A s(2.7),s(0.65),p₋₁,₀,₊₁(0.65); H@B s(1.0)
= 6 spatial / 12 spin-orb / C(12,4)=495 dets. Complex Yₗₘ ⇒ every function has a definite **m** about
the molecular z-axis (both centers on z), so m-conservation is the exact angular selection rule and l
is NOT (the point). Integrals from the validated `TwoCenterLM` engine (`two_center_grid_lm.py`,
~1e-6): S, h₁ (correct Z_A=3, Z_B=1), Coulomb ERI, geminal-w ERI (general-kernel multipole assembly),
and the 3-body L3 (grid-moment reformulation of Track 1's four_Y machinery — correlator expanded about
the common origin A, orbital angular content taken as numeric grid moments; radial = scalar multipole
of u′, angular selection exact, same model as the atomic run). RHF closed-shell reference (2 occ σ
MOs); xTC = replace Coulomb by w + contract L3 against the σ-MO 1-RDM (rotate to MO, contract, rotate
v2 back to the orthonormal AO basis). Files: `xtc_lih_molecular.py`, `xtc_lih_molecular_study.py`,
`xtc_lih_molecular_diag.py`, `xtc_lih_geminal_validate.py`, `debug/data/xtc_lih_molecular_*.json`.

## Validation gates (all PASS)
- **Two-center geminal ERI vs DIRECT r12 double-grid quadrature:** cross-density cases |Δ| 6e-7–4e-6;
  localized same-center-pair cases direct-grid-limited ~1e-3 (the multipole method is the accurate one
  — it reproduces the engine's native Coulomb to 1e-6). `xtc_lih_geminal_validate.py`.
- **geminal→0 (γ→∞):** shift E_xTC−E_plain → 0 monotonically (−1.33 @γ0.6 → −0.0011 @γ20); 1-norm
  ratio → 1.000. xTC → plain FCI. PASS.
- **Non-Hermitian GS real:** imag(E_xTC)=0 at every γ, R. PASS.
- **m-conservation:** 0 violations in the contracted-L3 tensor. PASS.
- **E_tot sanity:** −7.96 Ha (E_elec −8.96 + E_nn 1.0), reasonable minimal-basis LiH.

## Deliverable — fill-in + 1-norm, PLAIN vs xTC (R=3.0, γ=1.0)

| quantity | plain (Coulomb) | xTC (w + contracted-L3) | ratio |
|---|---|---|---|
| spatial 2-body density | 454/1296 = **35.0%** | 454/1296 = 35.0% | — |
| **fill-in** (xTC nonzero where Coulomb zero) | — | **0** | — |
| spatial 2-body 1-norm | 14.96 | 13.35 | **0.893** |
| spin-orbital ⟨pq‖rs⟩ nnz (qubit H) | 2368 | 2368 | — |
| **spin-orbital fill-in** | — | **0** | — |
| spin-orbital 1-norm (LCU-λ proxy) | 94.46 | 88.43 | **0.936** |
| AO-diagonal reference control 1-norm ratio | — | — | 0.940 |

**Robustness:** fill-in = 0 and 1-norm ratio ≤ 1 across **γ ∈ [0.6, 20]**, **R ∈ [2.5, 4.0]**, and with
**+H p_z** added. 1-norm ratio 0.914 (γ=0.6) → 1.000 (γ=20).

## The key contrast (why this is qualified, not a clean GO)
The impressive part of the atomic result was preserving a **rich l+m Gaunt sparsity**. That sparsity is
a *single-center* property and it is **gone for two centers before xTC even enters**:

| | Coulomb support | as % of its **m-conservation limit** | reading |
|---|---|---|---|
| **atom** (5-fn s+p) | 107/625 = 17.1% | **55%** of m-limit (195) | REAL l-sparsity: 45% of m-blocks Gaunt-zeroed |
| **molecule** (LiH, 6-fn) | 454/1296 = 35.0% | **100%** of m-limit (454) | NO l-sparsity: every m-allowed block filled |

Raw two-center Coulomb is already 80.6% of the m-limit; Löwdin orthonormalization (required for the
qubit Hamiltonian, and it mixes l across the cross-center overlaps) fills it to 100%. Count is
tol-stable (min |eri|=5e-5, nothing near the 1e-9 cutoff). So xTC has only m-conservation to preserve,
which it does (0 fill-in), and it also *lowers* the 1-norm — but there is no rich structure left to
preserve.

## Mechanism (reference-density multipoles about center A)
| MO | frac of multipole power in L=0 (monopole) |
|---|---|
| Li **core** MO (1s²) | **1.000** — a pure single-center s; the atomic collapse applies exactly |
| **bond** 2σ MO | **0.128** — collapse FULLY broken (87% higher-L) |
| total reference | 0.977 — core-power-dominated |

So the collapse "survives" for LiH not because the σ reference is monopole (the bond MO is emphatically
not), but because the tightly-bound Li **core dominates the density power**. The bond MO breaks the
collapse — yet still causes **0 fill-in**, because the two-center Coulomb baseline is already m-complete
(there are no l-zeros for the reactivated higher-L channels to fill). The contracted-L3 v2 is a real,
substantial operator (nnz=454, ‖v2‖₁=0.33) that co-supports *exactly* with Coulomb.

## Verdict
**GO (qualified).** For a molecule, xTC does **not** densify (0 fill-in at every γ/R/basis) and
**lowers** the 1-norm (0.89–0.94× at γ=1, →1 as γ→0) — the cusp is handled with no sparsity penalty,
so the s-reference collapse transfers in the sense the gate asks. **But** the atomic headline —
preserving a rich l+m Gaunt sparsity — does **not** transfer, because that l-sparsity is a
single-center property the two-center Coulomb baseline lacks (its support = the m-conservation limit
exactly). The molecular statement is therefore: *xTC inherits whatever angular sparsity the two-center
Coulomb has (m-conservation), with 0 fill-in and a lower 1-norm — it just has much less to inherit.*

**Caveats:** minimal mixed-exponent STO basis (not converged accuracy); model-level L3 radial
magnitudes (angular selection exact, same model as the atomic run); single spin-independent geminal;
classical non-Hermitian PoC; LiH-specific — collapse-survival is driven by the heavy Li core's monopole
dominance, a coreless molecule would stress the bond-MO collapse-break harder (though the 0-fill-in
conclusion is basis-geometric, not reference-dependent).
