# Sprint memo: Paper 45 adversarial hardening (Phase 1a/1b of corpus accessibility plan)

**Date:** 2026-06-09. **Status:** verdicts in; PI decision required before any paper edits.
**Plan reference:** `docs/corpus_accessibility_plan.md` Phase 1 (gate: "if 1b finds a real gap → fix or descope honestly before any outreach").
**Agent memos:** `debug/p45_adversarial_L2_memo.md`, `debug/p45_adversarial_L5_memo.md`, `debug/p45_hypothesis_checklist_memo.md`.
**Independent verification artifact:** `debug/p45_kplus_seminorm_check.py` (+ run output, below).

---

## 1. Headline

The adversarial refutation pass on Paper 45's two highest-risk lemmas found that **the main theorem does not survive as stated**. Three independent failures compound:

1. **The K⁺-restricted Lipschitz seminorm is identically zero on the entire operator system** — the quantity bounded by `thm:main` is degenerate in the cited framework (not a quantum compact metric space; the Monge–Kantorovich metric is undefined/infinite). This follows from the paper's own construction by short algebra (§3 below) and is consistent with four previously-celebrated "features" (the L3 structural identity, the bit-exact N_t=1 recovery, Λ^strong = Λ^P45 "free upgrade", and panel values equal to Paper 38's spatial rate at every cell) — all four are symptoms of the same fact: the computed quantity never sees the temporal factor or the restriction.
2. **The cited "Latrémolière Thm 5.5 / Def 3.4 / Def 3.5" do not exist.** Both candidate sources were fetched and checked: arXiv:1811.10843 (Adv. Math. 404 (2022)) has no Theorem 5.5 (its §5 ends at Remark 5.3); its items 3.4/3.5 are unrelated. The actual device used is van Suijlekom's C¹-approximation state-space GH framework (arXiv:2005.08544; Leimbach–vS 2024), dressed in Latrémolière vocabulary. The fabricated reference was inserted by the 2026-05-24 "Phase 1B-A closure" sprint — the sprint that claimed to close the L5 bookkeeping gap "via explicit Thm 5.5 transcription."
3. **Lemma L2's headline equation conflates two objects.** The "Plancherel symbol" (2j+1)/Z_n is the kernel's mass distribution, not the convolution-multiplier symbol; the true Fejér symbol decreases from 1 (forced by normalization: P38's own L2(a) ∫K = 1 implies symbol = 1 at the trivial rep, contradicting the claimed 1/Z_n) and the smoothing map is UCP with cb-norm exactly 1, not 2/(n_max+1). The Bożejko–Fendler 1991 citation (a discrete-groups theorem) and "Pisier 2001 Thm 8.10" do not support the claim; the fact actually needed is elementary. Downstream, the reach_P argument uses the wrong invariant (a forward cb-norm cannot bound a partial inverse; the inverse is controlled by the symbol minimum).

**What survives** (also a Phase-1 finding): Paper 38's spatial analytical content — γ_n closed form and the 4/π rate (kernel moments, independent of the symbol confusion), L3 Lipschitz bound C₃ = 1, L4 Berezin properties under the convolution reading, the J_L-equivariance algebra, and the L3 identity as bit-exact algebra (reinterpreted: it is the *diagnosis* of temporal invisibility, not a feature). The headline bound Λ ≤ C₃γ survives **numerically** because the false constants never enter the panel formulas — but its *meaning* requires reframing in the vS framework, and reach_P requires a genuinely new lemma (antiderivative-transference, Leimbach–vS Lemma 3.4/3.5-style; the independent Gaudillot-Estrada–vS arXiv:2310.14733 corroborates the conclusion for compact metric groups).

---

## 2. Verdict table

| Target | Verdict | Most damaging specific |
|---|---|---|
| P45 L2 object identification | WRONG (counterexample) | True SU(2) symbol at n_max=2 is (1, √2/3, 2/9), decreasing; claimed (1/3, 2/3, 0). UCP map ⟹ cb-norm = 1. |
| P45 L2 citations (BF 1991, Pisier Thm 8.10) | WRONG citation / UNVERIFIED-IMPORT | BF is discrete-groups uniformly-bounded-reps; needed fact is elementary block-scalar cb-norm. |
| reach_P (uses L2) | SERIOUS-GAP | Partial inverse controlled by symbol min (O(n²) on full envelope), not max; needs new transference lemma. |
| U(1) factor / N_t-independence | FIXABLE-GAP | Displayed symbol mismatches displayed kernel; sup = 1 survives; "2/(n+1), N_t-independent" exists only by mixing conventions across factors. |
| P38 L2(b),(c) antecedent | WRONG (inherited) | Same confusion; P40 L2(b),(c) repeats it for all compact groups. |
| K⁺ degeneracy (seminorm kernel) | WRONG — counterexample | Seminorm ≡ 0 on entire restricted system (§3); thm:main bounds a degenerate quantity. |
| "Latrémolière Thm 5.5/Def 3.4/3.5" | WRONG — nonexistent | No Thm 5.5 in either candidate paper; framework is tunnels/quantum isometries, never "UCP tunneling pairs." |
| Operator-system vs C*-hypotheses | FAILS 4/8 rows | Not an algebra; Leibniz not statable; kernel ≠ scalars; MK doesn't metrize weak-*. |
| K⁺ bookkeeping (rem:Kplus_tunnel_valid) | SERIOUS-GAP | Constituents computed against an extraneous classical-gradient ball, not the triples' own Lip-balls (suprema there are +∞). |
| P38 §L5 antecedent | FIXABLE-GAP | Same (B,P) device but honestly labeled "analog"; no fabricated citations; needs vS reframing + kernel check. |
| Citation sweep (36 imports) | 22 VERIFIED / 10 MISSTATED / 2 NOT-FOUND / 1 UNCHECKED / 1 UNVERIFIED | NOT-FOUND: farsi_latremoliere2024 (no such paper); hekkelman_mcdonald2024 (arXiv:2403.18619 resolves to an unrelated OpenMP paper). Plus 2 orphan \cites with no bibitem (compile-level defect). |

---

## 3. Independent main-session verification (verify-the-verifier)

Adversarial agents can over-claim; the two most damaging findings were re-derived in main session before acceptance.

**(a) Panel values are closed-form formulas, not operator computations.** `geovac/lorentzian_propinquity_compact_temporal.py` `LorentzianPropinquityBound.build()` computes `gamma_rate_su2(n_max) + gamma_rate_circle(N_t, T)` and Λ = C₃·γ. No operator-level propinquity is ever evaluated; `reduces_to_paper38_at_N_t_1()` checks a formula identity. Confirmed by direct code read.

**(b) The K⁺ compression annihilates the seminorm — algebra from the code's own construction.** `geovac/lorentzian_dirac_compact.py:144-188` builds D_L = i·(J_s ⊗ D_t + D_GV ⊗ I_t) with D_t anti-Hermitian Fourier-diagonal. The time part is Hermitian and J-commuting. The space part i·(D_GV ⊗ I) is anti-Hermitian (D_GV Hermitian), so Krein-self-adjointness D_L^× = D_L — asserted as L1′(b) — *forces* {J, D_GV ⊗ I} = 0. Then P₊(D_GV ⊗ I)P₊ = ¼(D + JDJ) = ¼(D − D) = 0: **the compression annihilates the spatial Dirac exactly.** The surviving compressed Dirac is momentum-diagonal; the temporal multipliers are momentum-polynomial diagonals (`operator_system_compact_temporal.py:98-144` — functions of momentum, not of time), so every commutator vanishes. The restricted seminorm is identically zero on the whole system. Numerical bit-exact confirmation: `debug/p45_kplus_seminorm_check.py`, run 2026-06-09:

```
Cell (2,3) dim_K=48,  42 multipliers:  ||{J, D_GV⊗I}|| = 0.0 exactly;
  ||P+ (i·D_GV⊗I) P+|| = 0.0 exactly (uncompressed 1.587e+01);
  max_a ||[P+ D_L P+, a+]|| = 0.0 exactly, kernel 42/42;
  max_a ||[D_L, a]|| = 2.251e-01, kernel 30/42.
Cell (3,5) dim_K=200, 275 multipliers:  same pattern bit-exact;
  restricted kernel 275/275; full-space kernel 130/275.
VERDICT: CONFIRMED at all cells.
```

Note the C4 side-finding: even on the FULL (unrestricted) Krein space the seminorm kernel is 30/42 resp. 130/275 — far larger than the scalars. The kernel-condition failure is not created by the K⁺ restriction; the restriction merely completes it (everything → 0). Both the restricted and unrestricted objects fail the compact-quantum-metric-space hypothesis.

**(c) The L2 trivial-rep contradiction is logically self-contained.** P38 L2(a) (∫K = 1, symbolically verified in project code) forces convolution symbol 1 at the trivial rep; the claimed symbol gives 1/Z ≠ 1. No external input needed.

Item 2 of §1 (nonexistence of Thm 5.5) rests on the agent's fetch of both source PDFs with specific structural descriptions; it will be re-verified once more during the rewrite, but the burden has shifted.

---

## 4. Downstream exposure

- **Papers 46, 47, 48, 49** consume P45's Λ quantity, panel, and framing ("free upgrade" Λ^strong = Λ^P45 is now explained as two degenerate quantities agreeing on a formula). Paper 40 L2(b),(c) repeats the symbol error for all compact groups.
- **Paper 38**: L2(b),(c) and the Berezin display need correction; §L5 needs vS reframing + an explicit kernel-condition check for the spatial system; reach_P needs the new transference lemma. The spatial theorem very likely survives repair.
- **CLAUDE.md / MEMORY / CHANGELOG**: "first published Lorentzian propinquity convergence theorem," "L5 closed," "strong-form closed," WH1-adjacent Lorentzian claims — all need honest revision after the PI picks a repair strategy.
- **Process ledger:** the fabricated Thm 5.5 entered via the 2026-05-24 closure sprint; "closure-by-citation-transcription" is now a *measured* logic-level failure mode of the workflow (distinct from string-level citation hallucination). Goes to the claims register and the §3 dead-end table once repairs are decided.

---

## 5. Repair options (PI decision)

**A. Descope + erratum (fast, ~2–4 sprint-days).** Reframe Papers 38/45/46 in the van Suijlekom state-space GH framework actually used; fix L2 to the correct elementary statement; add the reach_P transference lemma (adapt Leimbach–vS; cite Gaudillot-Estrada–vS); restate P45 as: GH convergence of the spatial truncations, carried J_L-equivariantly, with the temporal factor a metrically invisible spectator — explicitly *retracting* the "first Lorentzian propinquity theorem" claim; propagate to P46–49 + corpus records.

**B. Repair the Lorentzian content (research, weeks, uncertain).** Replace momentum-diagonal temporal multipliers with genuine Toeplitz compressions of multiplication operators e^{iqt} (Connes–vS S¹ style) so time becomes metrically visible (costs the L3 identity), and rework the K⁺ device (the compressed-Dirac definition self-annihilates; restriction must act on states, not compress D). Honest target after repair: "GH convergence of spectral truncations of a product geometry equivariant under a BBB Krein structure" — whether any statement deserving the word *Lorentzian* survives is precisely the open research question the paper claimed to close.

**C. A then B (recommended).** Erratum now — honesty first, and it unblocks outreach on the repaired Paper 38 spatial theorem (still plausibly novel at the SU(2)/CH-spinor level); the product/Krein repair becomes the named research target, replacing the freeze-period exception slot.

**N1 outreach note:** on hold either way until the repair lands; under option C it pivots to the repaired Paper 38 statement.

---

## 5a. Execution record (Option C, approved and executed 2026-06-09)

PI selected **Option C** (descope + erratum now; product/Krein repair as named research target). Executed same day:

- **Paper 45 → v2.0**: retitled; abstract replaced; `\section*{Erratum and descope}` (3 items); main theorem replaced by `thm:kplus_annihilation` + `cor:annihilation_numerics` + withdrawn-assembly record (`sec:withdrawn_assembly`, retains `thm:main` label on the explicitly-withdrawn statement) + `prop:spatial_conditional`; L2 corrected (UCP/cb-norm-1 + mass/symbol split); Berezin redefined (convolution form; `rem:berezin_inequivalence` documents the Toeplitz-vs-momentum-diagonal target mismatch); U(1) symbol normalization fixed; appendix Plancherel lemma rewritten with the correct CG symbol (worked n_max=2 values (1, √2/3, 2/9)); G1 reopened with sharpened spec, G2 metric claim withdrawn (norm-resolvent arrow survives); panel + "load-bearing falsifier" reinterpreted honestly; bibliography: farsi_latremoliere2024 (nonexistent) removed, hekkelman_mcdonald2024 (wrong arXiv ID) → hekkelman2021 (verified 2111.13865), latremoliere2018 title/volume fixed, minguzzi title fixed, Toyota initial fixed, vs2021_jgp + gaudillot_vs2023 added (titles web-verified this session). **Three-pass clean, 0 errors.**
- **Paper 38 → v2.0**: retitled (state-space GH); 5-item erratum section; L2(b),(c) + proof corrected; `eq:B_def` → convolution form with σ(N) weights on envelope N ≤ 2n_max−1; L5 reframed in vS framework, reach_P = named gap; thm:main + limit identification conditional; "first non-abelian case" repositioned per GE–vS; "'elsewhere' three occasions" footnote removed (not reproducible); fabricated hekkelman cites removed; microtype disabled (P45's established environmental fix); 4 pre-existing double-subscript idioms repaired. **Three-pass clean, 0 errors, 18 pp.** Note: this paper's own pre-submission flag (3) in `memory/paper_38_drafted.md` — "pin which exact propinquity is used" — anticipated the framework failure and was never executed; recorded as process evidence.
- **Papers 40/46/47/48/49 + Paper 32** (sub-agent, log at `debug/p45_descope_siblings_log.md`): erratum blocks/remarks applied per template; all three-pass clean; several pre-existing compile defects in P32/47/49 repaired in passing.
- **Falsifier frozen**: `tests/test_p45_kplus_degeneracy.py` (3 tests, pass in 1.3 s) — anticommutation + annihilation; restricted-kernel = 100%; truthful-vs-offdiag spatial kernel counts (10/14 vs 1/14).
- **Records**: CHANGELOG v3.106.0; CLAUDE.md version + §2 one-liner + §3 two dead-end rows + §1.7 WH1 PROVEN → PROVEN-conditional + §6 flags on 9 paper entries; memory files wh1_proven / paper_45_drafted / paper_46_drafted / paper_38_drafted / g2_metric_closed updated + MEMORY.md index lines replaced (no additions; index over budget).
- **Verified-citation discipline held**: every newly added citation (vs2021_jgp, gaudillot_vs2023, hekkelman2021, farsi_latremoliere2024-in-P38, F–L 2404.00240) was web-verified before insertion; one self-caught violation (a from-memory journal reference for Hekkelman) was stripped to arXiv-only data before compile.

**Outreach status:** N1 remains blocked. It unblocks when Paper 38's two named gaps close (reach_P transference adapted to the spinor substrate; either an offdiag-Dirac L3 re-derivation or a proved truthful→offdiag comparison). The option-B target (genuine Lorentzian/Krein quantum metric: Toeplitz temporal algebra, non-compressed Krein seminorm, kernel condition) is the standing freeze-period exception.

## 5b. De-versioning follow-up (2026-06-10, v3.106.1)

PI directive: no "v2 papers" — papers are a single source of truth corrected in place; git/Zenodo is the version record. All 8 papers stripped of erratum/v1/v2 scaffolding same-day (corrections fully retained; one History-and-retraction remark each in P38/P45; "Status and scope" notes in the siblings). All GATE: PASS via the new `debug/compile_3pass.sh`; zero leftover version vocabulary by grep. Companion PI directive on subagent token burn produced the standing rule `memory/feedback_subagent_token_budget.md`; the 6-paper sibling pass re-ran under it (sonnet, hard caps, scripted gate) at 141.7k tokens vs 366k for the equivalent pass the day before. Records synced in CHANGELOG v3.106.1.

## 6. What this validates about the plan

Phase 1 fired exactly as designed: one day of adversarial review prevented sending a refutable flagship to van Suijlekom and Latrémolière — the two people who would have found these errors in an afternoon, at the cost of the project's one shot at engagement. The finding *is* the success mode. It also recalibrates the assessment ranking: the most-likely-to-survive external contact is now the repaired Paper 38, not Paper 45.
