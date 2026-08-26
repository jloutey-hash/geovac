# Sprint memo — p-inclusive xTC PoC on Li (does the contracted 2-body inherit angular sparsity?)

**Date:** 2026-08-23  **Verdict: GO** (with an honest s-reference scope caveat).

**Question the s-only run could not answer.** The validated s-only Li xTC engine
(`xtc_poc_li*`, memo `sprint_xtc_poc_li_memo.md`) showed the 3-body TC operator L3
contracts to an effective 2-body operator exactly (1.7e-18 on ≤2-exc), recovers
93–98% of the genuine 3-body shift, and does NOT densify (1-norm ~0.90×). But
s-only has ONE angular block, so it could not test the real prize: does the
xTC-contracted 2-body inherit Track 1's angular Gaunt sparsity, or does the
3-body→2-body contraction FILL IN the zero angular blocks?

**Build.** Extended the engine to a minimal s+p Coulomb-Sturmian basis
(`max_n=2` → 1s, 2s, 2p₋₁, 2p₀, 2p₊₁ = 2 s + one p-shell, 5 spatial / 10 spin-orb,
120 dets). L3 built with **Track 1's validated four-harmonic W machinery**
(`four_Y`, reused directly) for the vertex + `G_leg` correlator legs; radial via the
scalar multipole of u′=½e^{−γr₁₂} per line. Coulomb + w(D) built on one grid via
the same correct complex-Y Gaunt. xTC contraction / Wick-to-bare / PN-projected FCI
plumbing reused verbatim from `xtc_poc_li.py`. Files:
`debug/xtc_poc_li_pinclusive.py`, `debug/xtc_poc_li_pinclusive_study.py`,
`debug/data/xtc_poc_li_pinclusive_study.json`.

**Convention note (framework aside, not fixed per scope).** GeoVac's production
`SturmianCI._ck_coefficient` uses a q-sign convention that silently **drops 42
physical m-changing 2-body Coulomb multipoles** at this basis (e.g.
⟨2p₊₁2p₋₁|1/r₁₂|2p₀2p₀⟩ = −0.0342 → returned 0). I build both Coulomb and L3 with
the correct gaunt (matches the framework EXACTLY on every m-conserving block), so the
plain baseline is physically complete and the comparison is internally consistent.

## Validation gates (all PASS)
- **G1 radial:** grid-multipole Coulomb vs high-accuracy adaptive quad, max|Δ|=**4.4e-6**
  (1s block = 0.937504 vs analytic 5k/8 = 0.9375). The framework `_slater_rk` is itself
  ~6e-2 off here — coarse 500-pt linspace, as the s-only memo already noted.
- **G3 m-conservation:** 0 / 2065 nonzero L3 entries violate total-m (rotational invariance holds).
- **G4 contraction fidelity:** |E_xTC − E_exactTC| = **1.7e-5** (E_xTC −7.508418 vs exact-3body −7.508402).
- **G5 sanity control (geminal→0):** E_xTC → E_plain monotonically, |Δ| 7.5e-2 (γ=4) → 2.9e-4 (γ=25).
  Non-Hermitian ground state stays **real** (imag = 0).

## Deliverable — angular density + 1-norm + Pauli, PLAIN vs xTC (γ=1.0)

| quantity | plain (Coulomb) | xTC (w + contracted-L3) | ratio |
|---|---|---|---|
| 2-body angular blocks (spatial, exact/radial-free) | **107/625 (17.1%)** | **107/625 (17.1%)** | — |
| **fill-in** (xTC nonzero where Coulomb zero) | — | **0** | — |
| 2-body 1-norm (LCU λ proxy, spin-orb) | 129.82 | 109.21 | **0.841** |
| total 1-norm | 52.71 | 47.95 | **0.910** |
| Pauli terms (openfermion JW) | 279 | 279 | **1.000** |
| Pauli 1-norm | 14.98 | 13.63 | **0.909** |

Raw 3-body L3 tensor density = 13.2% (2065/15625), consistent with Track 1's
l_ext=1/L_corr=2 regime (11.5%); Track 1's "~5%" was its full l_max=2 case. The point
holds either way: L3 is **angularly sparse** (low double-digit %, not the dense/collapsed
plane-wave-TC object), and the contraction preserves it.

**Robustness across geminal width** (fill-in and Pauli count are γ-invariant):

| γ | fill-in | L1[2b] ratio | L1[tot] ratio | Pauli ratio | 3-body shift | xTC fidelity |
|---|---|---|---|---|---|---|
| 0.60 | 0 | 0.760 | 0.874 | 1.000 | +60.9 mHa | 0.029 mHa |
| 1.00 | 0 | 0.841 | 0.910 | 1.000 | +22.5 mHa | 0.017 mHa |
| 1.50 | 0 | 0.907 | 0.945 | 1.000 | +7.5 mHa | 0.004 mHa |

## Why (structural, and the scope of the win)
Contracting a correlator line over an **s** reference orbital o=(l=0,m=0) forces
`G_leg(0,0,L′,M′,0,0)` → **L′=0**, which collapses the non-abelian four_Y vertex
`four_Y(a;L,M,0,0;d)` to a **single** multipole L (= g_ext). So the effective 2-body
reduces to the same single-multipole (a,d)×(b,e) coupling as Coulomb → identical angular
support, provably no fill-in. Li's aufbau reference (1s²2s) is all-s, so the collapse is
exact here. **Caveat:** this is a property of the *s-only reference*, not of xTC in
general — a p-block reference (occupied 2p) keeps L′≠0 and could fill in. Also: still Li,
single-common-k Sturmian (caps absolute accuracy, hits plain and xTC equally), single
spin-independent geminal, classical non-Hermitian PoC (a quantum algorithm needs QEVE/QITE).

## Verdict
**GO.** The xTC-contracted effective 2-body operator inherits the angular Gaunt
sparsity **exactly** (0 fill-in at every γ, m-conservation intact), keeps the **same
Pauli count** (279), and **lowers** the 1-norm (0.76–0.91× for 2-body, 0.91× Pauli-λ) —
a genuine qubit/sparsity win with the cusp handled, on the first basis that could test it.
Track 1's silver lining is **confirmed** for s-referenced systems.

**Molecular test warranted?** Yes — H₂ (σ from 1s, s-reference) and LiH (Li 1s2s + H 1s,
s-reference) should inherit the collapse; they are the natural next PoC. The sharper,
higher-value test is a **p-block reference** (e.g. C/N/O, occupied 2p), where the L′=0
collapse no longer applies and fill-in becomes possible — that is where "does the sparsity
survive" is genuinely at risk and worth measuring.
