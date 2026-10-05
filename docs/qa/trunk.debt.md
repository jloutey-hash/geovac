# Trunk SMALL-debt register

SMALL findings that could not be fixed on sight (criteria.md RESEARCH-GRADE).

## SMALL sweep (edit bar) -- SWEPT 2026-10-03 (v5.20.2)

The FULL #12 and FULL #13 SMALL lists were fixed in one scheduled sweep (four
file-disjoint agents: Paper 32; Paper 38 + Fejer modules; Papers 0/1/7 + tests;
syntheses + non-trunk papers + misc), backstopped by C16 entry
`trunk-small-sweep-2026-10-03`.  The swept units show CHANGED and are reviewed by
the next run.  Item-by-item record: CHANGELOG v5.20.2.  Residue (not fixed):

- **Needs a new guard (separate guard pass).** Paper 38 L4 convolution form B = P(K*f)P on the true window; good-kernel estimate (iii).
- **API renames deferred.** geovac/modular_hamiltonian.py `restrict_K_alpha_to_wedge` (a relabelling); central_fejer_su2 `asymptotic_rate_constant` / `gamma_rate` / `quantitative_rate_*`; geovac `UNIVERSAL_KINETIC_SCALE` constant name.
- **[CLOSED 2026-10-04, v5.20.4] Infrastructure.** The genuinely heavy tests (true-window odd_no_grading[3] 11.6 GB; trunk_qa_splitting n_max 25/30 1.7-4 GB; full_dirac propagation 2.2 GB) are in `xdist_group("heavy_memory")`; tests/conftest.py switches -n runs to --dist loadgroup. The other reported offenders were collateral (they failed beside the heavy ones).
- **Chronicle (do-not-edit).** Korevaar Remark 3 pointer in docs/chronicles/paper_32_chronicle.tex unverified.
- **History notes outside swept units.** Further dated / "earlier version" notes in Papers 39 and 40 (e.g. P39 L221, L507, L1229-1273; P40 L1057, L1931-1960), several holding C16 retraction tokens -- belongs to the chronicle-relocation pass (item 5). Paper 7 Sec. 6 Li nu=2 paragraph has no inline tier tag.
- **[CLOSED 2026-10-04, v5.20.4] Paper 32 audit item 4.** The obstruction is forced by J^2 = -1 (both m_j = +-1/2 partners hop to the same neighbour with the same real value), not by the (-1)^l phase, which only selects which hops fail; residual = 2 max(alpha, coupling). Paper 32 item 4 rewritten; tests/test_real_structure.py gains two independent, fire-tested checks.
- **[CLOSED 2026-10-04, v5.20.4] Pre-existing test failures.** Two Paper 27 tests: the 2026-08-29 sign fix (b32581a) gauge-flips V within degenerate H1 eigenspaces, so the 'diagonal fraction in the H1 eigenbasis' is basis-dependent -> pinned the basis-free eigenspace-block fraction instead; the Be EP-2N occupations were stale and broke 1S symmetry -> Paper 27 updated (S_full 0.794 -> 0.821, reduction 42% -> 40%, EP-2j FCI S 0.902 -> 0.916). Paper 39 status test asserted a reworded phrase -> asserts substance. All fire-tested. Open wording item: Paper 27 should define the diagonal fraction basis-free (degenerate H1).
- **Not exercised in FULL #13 (owed before PASS).** C23 inverse citation (C_sat closed form; gearing 2 sqrt6/3; taper O(1/n) vs Jackson/Fejer-Riesz; c^2(n,l)). Dedicated arithmetic-audit pass. Stamp leak (P32 ::Conclusion) -> double the AUDIT fraction next run.

## Open

- 2026-10-03 (DELTA #26) — Paper 40 PRV upgrade owed: (INT) holds on every PRV summand at all ranks (Cas(sigma) = |sigma|^2 + 2<sigma,rho> >= |lambda + w lambda'^*|^2 >= |lambda - lambda'|^2); sound by inspection, needs an exact-arithmetic test before it replaces the v5.19.0 'rests on panels' caveat.
- 2026-10-03 (DELTA #26 code review) — `tests/test_paper39_lifted_state.py` T2 hardcodes sigma_J rather than computing it from `central_fejer_kernel_su2`; a guard improvement, separate pass.
- 2026-10-03 (DELTA #26 self-check #10) — Paper 28's 'Sigma(1) carries Catalan G' rests on an unrecorded-precision PSLQ at 50 dps with no regression test (now labelled as such in the text); coverage gap.
- 2026-10-03 (DELTA #26) — five of the ten new C16 entries use vocabulary exemptions (e.g. 'corrected', 'withdrawn') rather than the per-entry `[retracted ...]` token the registry now prefers; convert when next touched.
- 2026-10-03 (DELTA #26) — `viz/public/papers/*.html` (incl. group5 synthesis L29) still carries the pre-DELTA-#26 wording (standing PI-gated generated-layer item).
- 2026-10-02 (DELTA #25, qa.md rule 3) — `README.md` carries its headline figures as bare values (~30: Pauli coefficients, equal-qubit ratios, LiH 838 / 907 / 276, 1-norms, nuclear 688 / 828 / 710, H2 99.97%, LiH +0.2% / −8.0490, He 1.88%, H2O 19.4%, ...). Under the new rule they should point at numeric-registry keys or the owner paper (Markdown `<!-- gvq:key -->` or a pointer); every value was re-sourced against its owner in this round, so this is a conversion, not a correctness debt. Blocked on extending C21 to read Markdown annotations (it reads `.tex` only), so the pointers can be checked rather than merely written.

- 2026-10-02 (DELTA #25 self-check #15, out of scope) — Paper 39 states the elementary triangle bound for the displayed $D_{a,b}$ that the same paper calls unrealisable ($\gamma_a = 0$); the bound is trivially true there, but the 'tight at top-shell harmonics' test context (the code's graded multipliers) should be named.
- **[CLOSED 2026-10-03, DELTA #26]** 2026-10-02 (self-check #15, out of scope) — Paper 55's M1 table source cell still reads '(Hopf base)' (label retired corpus-wide at DELTA #19). Now '(quarter unit-sphere area)'.
- 2026-10-02 (self-check #15, out of scope) — Paper 34 ~L1856 keeps H2 96.0% as the Level-4 anchor (legitimate for Level 4; not the corpus-best H2 figure, which is Level 2).
- 2026-10-02 (DELTA #25 loop scope note) — rounds #12-#14 read whole summary surfaces (group1/3/4 syntheses, field guide, README, one-pager, claims register) and found pre-existing debt each time, all fixed. **Not yet read in full under the new cite-don't-paraphrase rule:** the group2, group5 and group6 syntheses and docs/geovac_positioning.md. Owed to the next FULL run.
- 2026-10-02 (standing, PI-gated) — `viz/public/papers/*.html` is a stale generated layer rebuilt from the Zenodo manifest (declared gap since paper_60 FULL); every DELTA #25 correction is absent from it.

_The five items below were closed 2026-10-02 (v5.18.1 debt pass); kept as the record._

## Closed

- **[CLOSED 2026-10-02]** 2026-09-29 (DELTA #15) — `paper_34_projection_taxonomy.tex` ~L3947–4769: an 820-line `[h]` float table sits at TeX's maximum dimension ("Dimension too large" at `\end{table}` once ~2 lines are added). Blocked on a structural edit (convert to `longtable`); not a content defect.
- **[CLOSED 2026-10-02]** 2026-09-30 (FULL #11) — `paper_32_spectral_triple.tex`: 19 inline `Author~Year` attributions resolve to no `\bibitem` (C20 reports them; they sit inside its ratchet baseline, so the gate passes). Blocked on primary-source verification of each before a bibitem is added — adding unverified metadata would launder a guess into a citation (criteria C4). List: `python debug/qa/check_inline_attributions.py --gate trunk`.
- **[CLOSED 2026-10-02]** 2026-10-02 (DELTA #24) — group3 synthesis "Kac form" (Cas(ad) = 2h∨ under (θ|θ) = 2): the citation (Kac, *Infinite Dimensional Lie Algebras*, 3rd ed.) is confirmed but the specific normalisation statement was not fetched from the book (UNVERIFIABLE). Blocked on access to Kac Ch. 6.
- **[CLOSED 2026-10-02]** 2026-10-02 (DELTA #24) — Paper 7 New Contributions item 3 claims the scale-invariance framing "has not, to our knowledge, been stated explicitly in the literature" [OBSERVATION]; the paper itself concedes (§Discussion) it is implicit in Fock and the hyperspherical literature. Blocked on a primary-source check of Bander–Itzykson 1966 and Avery's Sturmian treatment before the novelty claim stands or is withdrawn.
- **[CLOSED 2026-10-02]** 2026-10-02 (DELTA #24) — Paper 32 tags ζ_{D²}(5/2) = 14ζ(3) − 31/2 ζ(5) under the M2 sub-mechanism (~L2441) while §~3528/3613 place odd ζ in M3. Both can hold (an M2 output whose value lies in the M3 ring), but the paper does not say so. Blocked on an owner ruling (Paper 55) on whether slot-indexing is by mechanism or by ring.

### How each was closed (2026-10-02)
- **P34 table:** converted to `longtable` (caption + label at the head, repeated header); compiles at 162 pages with no "Dimension too large". Its He 0.004% row was scoped to "extrapolation, non-variational" in the same edit.
- **P32 inline attributions (19):** a citation-reviewer verified each against primary sources (17 metadata CONFIRMED, CCM year WRONG 2010 -> 2007, "Clarified 2026" a false positive). 16 bibitems added, 3 resolved to existing keys (Krajewski 1998, Glanois 2015, Marcolli--van Suijlekom 2014); sentence fixes: Adams covers rungs n >= 5 only (S^7 is an H-space), Hawkins defines a Berezin-type limit quantization (not "the Berezin map"), Hochschild "Ch. VII" dropped, CM04's pro-unipotent factor not "conjectural" and U^* a group, Sylvetsky-Martin labelled a general (not NaH-specific) estimate. P32 now has 0 unresolved attributions. **Instrument defect found on the way:** C20 baselined only 5 of 57 papers and skipped any paper missing from the baseline, so ~100 attributions in 29 papers could never fail; unbaselined papers now FAIL (`tests/test_inline_attributions_gate.py`, fire-tested) and the baseline covers all 57 papers (147 tags recorded as corpus ratchet debt; the trunk's six documents are at 0).
- **Kac normalisation:** no longer rests on the book. `tests/test_kac_normalisation_casimir.py` computes Cas(ad) = 2h^vee under (theta|theta) = 2 for su(2), su(3) (and 1 under the Killing form; fire-tested); the synthesis cites Kac only for the convention's name.
- **P7 novelty claim:** withdrawn and re-tagged [PRIOR ART]. Bander--Itzykson (Rev. Mod. Phys. 38, 330 and 346, 1966; Part II read in the open-access copy) rescale momentum by sqrt(2 mu |E|) onto the unit sphere / unit hyperboloid; Part II bibitem added.
- **P32 M2/M3 odd zeta:** no owner ruling needed. Paper 55 (M2 section) already indexes slots by mechanism and states ring membership separately; P32 now says so explicitly (zeta_{D^2}(5/2) is an M2 output whose value lies in the mixed-Tate ring M3's values also occupy).
