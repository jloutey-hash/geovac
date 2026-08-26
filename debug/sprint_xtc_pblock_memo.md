# Sprint memo — p-BLOCK reference xTC: the decisive sparsity test

**Date:** 2026-08-23  **Verdict: GO** (with an honest not-bit-exact caveat).

**The question the s-reference runs could not answer.** The atomic p-inclusive run
(`sprint_xtc_pinclusive_memo.md`) got **0 fill-in** ONLY because Li's reference
(1s²2s) is all-**s**: contracting a correlator leg over an s orbital forces its
multipole L′=0 (monopole → Coulomb support), collapsing the shared vertex to a
single multipole. A **p-block reference** (occupied 2p, l=1) keeps L′∈{0,2}, so the
3-body→2-body contraction CAN fill the angular zero-blocks. This is where "does the
sparsity survive xTC in general" is genuinely at risk.

**Build.** Reused the validated atomic s+p xTC engine (`xtc_poc_li_pinclusive.py`:
four-harmonic L3 via Track 1's `four_Y`, correct complex-Y Gaunt Coulomb, xTC
contraction, PN-projected non-Hermitian FCI) **unchanged** — only swapped the
reference occupation to a p-block atom. References: **C 1s²2s²2p²** (Hund ³P:
2p₀↑2p₊₁↑), **O 1s²2s²2p⁴**, plus a **Be 1s²2s² s-reference control**. Two probes:
(A) full radial-weighted s+p engine (gates + density/1-norm/Pauli); (B) a fast,
radial-free EXACT angular-support test extended to s+p+d and s+p+d+f (angular
selection is radial-independent, so support fill-in is exact there).
Files: `debug/xtc_pblock_{engine,fast_angular,study}.py`, `debug/data/xtc_pblock_{study,fast_angular}.json`.

## Gates (all PASS with a p-block reference)
- **m-conservation:** 0/2065 nonzero L3 entries violate total-M for Be, C, O
  (the diagonal reference contraction forces M′=0, so rotational invariance survives
  even with a non-spherical open-p reference).
- **Contraction fidelity |E_xTC−E_exactTC|:** C 1.6e-4, Be 4.0e-4, O 5.1e-3 Ha
  (≪ the ~1.4 Ha 3-body shift; O larger only because a single-det 2p⁴ 1-RDM is a
  poorer reference — a reference-quality artifact, not a sparsity issue).
- **Non-Hermitian ground state real:** imag = 0 in every case.
- **Geminal→0 (C):** E_xTC → E_plain monotonically, |Δ| 2.1e-1 (γ=4) → 1.1e-3 (γ=25).

## Deliverable — angular density / fill-in / 1-norm / Pauli, PLAIN vs xTC (s+p, γ=1)

| ref | fill-in (v2 vs Coulomb) | Coulomb density | 1-norm[2b] ratio | Pauli terms | Pauli ratio |
|---|---|---|---|---|---|
| **Be (s-ref)** | **0** | 107/625 = 17.1% | 0.796 | 279 → 279 | 1.000 |
| **C (p-ref)** | **0** | 107/625 = 17.1% | 0.841 | 279 → 279 | 1.000 |
| **O (p-ref)** | **0** | 107/625 = 17.1% | 0.885 | 279 → 279 | 1.000 |

At the **s+p** basis the p-reference gives **0 fill-in** and **identical Pauli count**,
exactly like Li's s-reference — xTC still LOWERS the 2-body 1-norm (0.80–0.89×).

## The decisive nuance — EXACT angular fill-in vs basis richness (radial-free)

| basis | s-ref fill-in | **p-ref fill-in (C, O)** | Coulomb density | p-ref L3 density |
|---|---|---|---|---|
| s+p     | 0 | **0** (0.000%)   | 14.84% | 14.84% |
| s+p+d   | 0 | **4** (0.061%)   | 8.52%  | 8.58%  |
| s+p+d+f | 0 | **96** (0.146%)  | 6.06%  | 6.21%  |

- **s-reference: bit-exact 0 fill-in at EVERY basis** (reference-density multipole =
  {(0,0)} monopole only; the collapse is basis-independent). Li's silver lining is a
  property of *any* s-referenced system at *any* basis.
- **p-reference: reference-density multipole = {(0,0),(2,0)}** — the L′=2 leg is ON.
  **0 fill-in at s+p** (l≤1 is protected: every same-parity orbital pair contains the
  monopole, so Coulomb never misses a same-parity block). Fill-in APPEARS once d
  orbitals are present and grows slowly with l_max (0→4→96) — but the density stays
  **single-digit %** and tracks Coulomb within ≤2.5% relative (6.06→6.21% at s+p+d+f).

**Mechanism.** The leg-3 multipole is a diagonal density multipole ⟨o|Y_{L′M′}|o⟩ =
δ_{M′0}·(nonzero only for EVEN L′). A spherical density (s orbital, or a *closed*
shell) has only L′=0 → vertex resultant Λ = leg-2 multipole L → identical to Coulomb's
single-multipole support → 0 fill-in, exactly. An **open** p orbital carries L′=2 too:
the vertex can then bridge an even-multipole gap Δ=2 between the (a,c) and (b,d)
external pairs, filling blocks Coulomb misses — but only "one even step away" blocks
(e.g. (s,s|s,d): a monopole pair meeting a quadrupole-only pair). The fill-in reach is
set by the **reference's** angular content (p → {0,2}), not the basis; s+p has no
gapped blocks, d orbitals create the first ((s,d)) ones.

## Verdict
**GO.** The xTC-contracted effective 2-body operator **stays sparse with a p-block
reference**: density remains single-digit % (6–9%), the Pauli count is unchanged
(279→279), and the 2-body 1-norm is lowered (0.80–0.89×). All gates pass. So xTC
preserves GeoVac's angular Gaunt sparsity **generally**, not only for s-references —
the qubit/Pauli edge with the cusp handled survives into the p-block.

**Honest caveat (the difference from the s-reference).** It is no longer *bit-exact*.
Unlike the s-reference's 0 fill-in at any basis, an open-p reference produces a small,
structured, slowly-growing fill-in (0.06% at s+p+d, 0.15% at s+p+d+f; ≤2.5% relative
density inflation) as soon as l≥2 orbitals are in the basis, driven by the open p
shell's L′=2 density multipole. It does not densify the operator across the physically
relevant s/p/d/f range, but the fill-in fraction is not asymptotically zero — it rises
with l_max, so a very-high-l basis would warrant a re-measurement. Also still: single
common-k Sturmian (caps absolute accuracy — hits plain and xTC equally, irrelevant to
sparsity), single spin-independent geminal, classical non-Hermitian PoC.
