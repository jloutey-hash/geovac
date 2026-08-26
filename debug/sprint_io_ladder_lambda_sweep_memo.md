# I/O ladder — LG-3 closure: λ (LCU 1-norm) vs basis richness

**Date:** 2026-08-17 · **Type:** diagnostic (no geovac/ / papers / CHANGELOG / version edits)
**Driver:** `debug/io_ladder_lambda_sweep.py` · **Data:** `debug/data/io_ladder_lambda_sweep.json`
**Closes:** the LG-3 leg of the block-encoding crossover (`debug/sprint_io_ladder_costmodel_memo.md`).

## Question
Rung-3 cost model: on-device *generate* beats *load* iff (G+S ≪ L) **and** λ is preserved
(LG-3). Rungs 1–2 cleared the I/O legs. LG-3 was the one unmeasured leg: does the LCU
1-norm λ = Σ|c_i| of the GeoVac qubit Hamiltonian **inflate with basis richness** (the
plane-wave penalty), or stay comparable to a matched Gaussian basis?

## Method
JW Pauli 1-norm, λ = Σ|c_i|, computed through ONE identical pipeline for both families
(so the ratio is convention-independent). Reported λ_excl (identity/constant removed = the
true block-encoding cost, the lead metric) and λ_incl (matches the corpus `one_norm`).
GeoVac: `ecosystem_export.hamiltonian(name, max_n).​_qubit_op`, max_n swept.
Gaussian (no pyscf available): He STO-3G/cc-pVDZ (`gaussian_reference`, hardcoded),
H2 STO-3G/6-31G + LiH STO-3G (openfermion cached MolecularData). Alignment by qubit count Q.

## Results (λ_excl, Q)
| family | points (Q → λ_excl) | exponent p (λ~Q^p) |
|---|---|---|
| GeoVac He | 2→1.69, 10→10.77, 28→64.92, 60→214.3 | **1.43** |
| GeoVac H2 | 10→6.20, 28→33.15, 60→107.0 | **1.59** |
| GeoVac LiH | 30→26.61, 84→156.75 | **1.72** |
| Gauss He (STO-3G,cc-pVDZ) | 2→1.69, 10→33.12 | **1.85** |
| Gauss H2 (STO-3G,6-31G) | 4→1.27, 8→11.45 | **3.17** (2-pt, small-Q) |
| Gauss LiH (STO-3G) | 12→12.37 | — (single point) |

**Matched-Q ratios (GeoVac/Gaussian, λ_excl):**
- He Q=2: 1.688 / 1.688 = **1.00×** — GeoVac He max_n=1 *is* STO-3G He (single 1s); pipeline check.
- He Q=10: 10.77 / 33.12 (cc-pVDZ) = **0.33×** — GeoVac λ is a third of the matched Gaussian.
- LiH Q30 (indicative, extrapolating STO-3G by p=1.85–3.17): **0.12–0.39×**.

## Verdict — LG-3 FAVORABLE
GeoVac λ does **not** inflate with basis richness. Its λ-vs-Q exponent (1.43/1.59/1.72) is
**below** the Gaussian exponent (1.85/3.17) in every case, and at matched qubit count GeoVac
λ is **lower** (He Q10: 0.33×). No plane-wave penalty. Mechanism: GeoVac keeps a compact
Coulomb-Sturmian-like basis rather than a delocalized plane-wave one, so the generate-friendly
representation does not pay λ to buy low I/O. With Rungs 1–2 (I/O legs) already cleared, **all
three legs of the Rung-3 crossover are now on the favorable side.**

## Caveats (do not overstate)
1. Cleanest evidence is **He** (matched-Q at both ends + full exponent). H2/LiH rest on GeoVac
   exponents (3/2 points) vs thin Gaussian data (H2 2-point steep; LiH single point → LiH ratio
   is *indicative* extrapolation, not measured — pyscf missing blocks real LiH multi-basis).
2. Matched-**size**, not matched-**space**: GeoVac hydrogenic shells and contracted Gaussians
   span different 5-orbital subspaces at Q10.
3. H2 R not perfectly matched across sources (minor; He atom has no such issue).
4. Small-molecule regime only (consistent with the whole ladder).

## Side-flag for PI (separate from the LG-3 verdict) — audit of the §1.5 "0.95× vs STO-3G"
This pipeline **reproduces the corpus GeoVac number exactly**: LiH max_n=2 λ_incl = **32.59**
(§1.5 "live 32.6"). But STO-3G LiH from the repo's own openfermion cache gives λ_incl = **16.46**
(λ_excl 12.37), whereas §1.5 quotes **34.3** for the STO-3G comparator — a **2.08× discrepancy**.
It is not an identity-term effect (λ_incl−λ_excl is only ~4). If STO-3G LiH is really ~16.5, the
"1-norm 0.95× vs STO-3G" line (already a mismatched-Q comparison, Q30 vs Q12) would instead read
~2× *above* STO-3G. The 34.3 source is not the openfermion cache (which `lih_sto3g_from_cache`
uses); likely a pyscf-era or different-convention/geometry computation. **Needs reconciliation
before the 0.95× line is relied on.** Does not touch the LG-3 verdict (internally matched pipeline).

## Cross-refs
Ladder: `sprint_io_ladder_accounting_memo.md` (Rung 1), `sprint_io_ladder_radial_seeds_memo.md`
(Rung 2), `sprint_io_ladder_costmodel_memo.md` (Rung 3 model + LG-1/2/3 definitions).
