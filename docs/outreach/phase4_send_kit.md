# Phase 4 Send Kit (prepared 2026-06-10 — NOTHING SENT)

Everything below is staged. Sends require the checklist at the bottom, and each send is an explicit PI action.

## Send-gate status — RESOLVED (2026-06-10)

The S1-vs-S2 decision is **moot**: the G1/G2 gap-closure sprint landed same-day (v3.109.0, translation-seminorm reframing; see `debug/sprint_p38_g1g2_phaseA_memo.md`), so the strongest posture — Option S2's "close the gaps first" — was achieved at no schedule cost. Paper 38's theorem is **unconditional**; N1 and the drafts below state it as such. The remaining gate is simply: **release precedes any send** — the notes point at repository content, which must be the corrected corpus.

## Recipient ladder

### N1 (math.OA note)
| Rung | Recipient | Rationale | Wait before next rung |
|---|---|---|---|
| 1 | Walter D. van Suijlekom (Radboud) | The framework is his; CvS 2021 defers the convergence question; GE–vS is his student's extension | 3 weeks |
| 2 | Frédéric Latrémolière (Denver) | Propinquity originator; Q2/Q3 are framework-fit questions for him | 3 weeks |
| 3 | Malte Leimbach; Eva-Maria Hekkelman; Yvann Gaudillot-Estrada | Early-career, higher reply rates; the theorem is the spinor extension of GE–vS's scalar compact-group result | — |

### N2 (periods note)
| Rung | Recipient | Rationale |
|---|---|---|
| 1 | Francis Brown (Oxford/IHÉS) and/or Axel Kleinschmidt (AEI Potsdam) | The sector-to-ring question is squarely in the cyclotomic-MZV/periods program; identities are 10-minute-checkable |

Recommended sequencing: N1 rung 1 first; N2 ~2 weeks later (any content engagement from the math.OA side strengthens the periods contact; Brown is a one-shot contact).

**Affiliations above must be re-verified at send time** (people move; the project has a measured failure mode of stale facts).

## Email drafts (≤6 sentences each; note attached as PDF; plain text)

### E1 — van Suijlekom (N1, rung 1)
Subject: A degeneracy obstruction for Lorentzian spectral-truncation metrics (3-page note)

> Dear Professor van Suijlekom,
>
> I am an independent researcher working on spectral truncations of the Dirac triple on S³ = SU(2) in your state-space Gromov–Hausdorff framework. The attached three-page note contains (i) a convergence theorem on the chirality-doubled spinor substrate — the spinor extension of the Gaudillot-Estrada–van Suijlekom compact-group result — with explicit rate (4/π)·log n/n, obtained by metrizing the truncations with the translation seminorm (the Dirac-commutator seminorm degenerates at finite cutoff) and an exact-fit spinor lifted state; and (ii) an elementary but, to my knowledge, unrecorded degeneracy theorem: compressing a Krein-self-adjoint product Dirac to the Krein-positive subspace annihilates the spatial Dirac, so the natural route to a "Lorentzian propinquity" is structurally closed. My three questions are in §5 of the note; the most important is whether either observation is already known. Full write-ups, code, and bit-exact verification scripts are at https://github.com/jloutey-hash/geovac. A one-line "this is known, see X" would be exactly as valuable to me as any other reply.
>
> With thanks for your time,
> Josh Loutey

### E2 — Latrémolière (N1, rung 2)
Subject: Is this Krein-compression degeneracy known? (3-page note on spectral truncations)

> Dear Professor Latrémolière,
>
> The attached three-page note records two results on spectral truncations of the Dirac triple on SU(2): a state-space GH convergence theorem at explicit rate (4/π)·log n/n, and a degeneracy theorem showing that the obvious device for a "Lorentzian propinquity" — compressing a Krein-self-adjoint Dirac to the Krein-positive subspace — annihilates the spatial Dirac, so no quantum-metric structure survives. I would particularly value your view on Question 3: whether the repaired problem (a Toeplitz temporal multiplier algebra, and a Krein-compatible Lipschitz seminorm that does not compress the Dirac) is well-posed within, or adjacent to, the propinquity framework, and whether anyone is already working on it. Full materials at https://github.com/jloutey-hash/geovac. Any reply, including "this is known" or "this is ill-posed because X," is valuable.
>
> With thanks,
> Josh Loutey

### E3 — Leimbach / Hekkelman / Gaudillot-Estrada (N1, rung 3; personalize the first sentence per recipient)
Subject: A spinor extension of compact-group truncation convergence (3-page note)

> Dear Dr. [NAME],
>
> The attached three-page note proves a state-space GH convergence theorem for Dirac-triple truncations on SU(2) — in effect the chirality-doubled spinor extension of [your compact-metric-group convergence theorem / the GE–vS compact-group result], with the dual direction supplied by an exact-fit spinor lifted state rather than a transference estimate (the Camporesi–Higuchi shells turn out to be exactly the V_j⊗V_{j±1/2} window blocks). The note also contains a short degeneracy theorem for the Lorentzian extension that may interest you independently. Question 2 of the note is the concrete ask — is the construction known, or is there a slicker route? — and I would value any pointer, including "this is routine" or "this fails because X." Full materials at https://github.com/jloutey-hash/geovac.
>
> Best regards,
> Josh Loutey

### E4 — Brown / Kleinschmidt (N2)
Subject: Three 10-minute identities on S³ spectral zetas, and one motivic question

> Dear Professor [NAME],
>
> The attached two-page note states three exact identities for Dirac spectral zetas on the round S³ — a parity-sector identity collapsing to the Dirichlet-beta line, the vanishing of the spectral zeta at non-positive integers by a Bernoulli mechanism, and an all-orders closed-form heat trace — each checkable in about ten minutes from the definitions given. The question (§ end) is whether the observed sector-to-ring assignment — untwisted spectral sums landing in ⊕ₖ π^{2k}·Q, parity-twisted sums in level-4 cyclotomic mixed Tate, with the even/odd difference realizing what looks like a Galois descent — is an instance of a known functoriality in the cyclotomic mixed-Tate category, and if so what else it predicts. If the descent reading is naive, I would genuinely value knowing that too. Full materials at https://github.com/jloutey-hash/geovac.
>
> With thanks for your time,
> Josh Loutey

## Pre-registered outcomes (logged like sprint results, per plan §4d)

| Outcome | Classification | Action |
|---|---|---|
| Content-engaged reply, incl. refutation ("wrong because X") | **SUCCESS** | Log verbatim; open a gap-closure or correction sprint |
| "This is known, see X" | **SUCCESS** (calibration) | Log; cite X; reposition the result |
| "Interesting, but [objection to a proof step]" | PARTIAL | Open a verification sprint on the named step; reply with findings |
| Silence ≥ 3 weeks | NULL | Move to next rung; log |
| Request for collaboration / more material | SUCCESS+ | PI decision before any commitment |

## Send checklist (each item checked by the PI at send time)

1. [ ] Release done; the public repo contains the corrected corpus (v3.106+ — degeneracy framing, claims register, corrected `.zenodo.json`).
2. [x] `https://github.com/jloutey-hash/geovac` filled in both notes (footnotes) and all four drafts.
3. [ ] Recipient affiliations/emails re-verified (web, day-of).
4. [x] ~~PI decision on Option S1 vs S2~~ — moot (gaps closed 2026-06-10, v3.109.0; see Send-gate status above).
5. [ ] PI has read N1/N2 end-to-end as final sign-off.
6. [ ] One send per rung; outcomes logged in this file + CHANGELOG.
