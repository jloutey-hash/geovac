# /qa claim-search DELTA on v5.20.21 (frozen at 31549dd8) -- PM ledger

Verified against primary text by the PM before acceptance.

## Gates
- All deterministic gates green (trunk/group1/2/3/6/synthesis at checkpoint; paper_59 scope this run).
- C25 FAIL x1 (new entry p59-coulomb-ring-elliptic): papers/INDEX.md:104 "a ring of three densities sits on the elliptic rung, below the chain" -- no Coulomb scope. SURFACE-LARGE (front door; correlated ring sits on K3).

## Verified findings
| # | Locus | Family | Grade | Note |
|---|---|---|---|---|
| G1 | INDEX.md:104 | coulomb-ring | SURFACE-LARGE | C25 caught |
| T1 | P59 L1222-1224 "a loop lowers the rung when it integrates out exactly" | correlated-triangle / coulomb-ring | OWNER-LARGE | contradicted by body L1199-1203: massive loop integrates exactly (solid angle) and does not lower; mechanism is uniqueness; also a reading under a SYMBOLIC tag |
| T2 | P59 conclusion L1363-1364 "whose loop does not [integrate out exactly]" + P57 L159 same | same | SURFACE (P59 concl in owner -> with T1) LARGE | false contrast, same family as T1 |
| T3 | P59 abstract L95-98 "drop comes from the exact integration of the Coulomb loop ... correlated triangle ... sits with the chain, on the K3 rung" | same | SMALL (rides with T1 fix) | no gamma scope |
| T4 | P57 L172 (J11) "exactly at gamma=1 and by local data at gamma=3/2" | correlated-triangle | SURFACE-LARGE | drops "factor above K3 not excluded" |
| T5 | P59 L1213 "no operator of order three annihilates Y" unbounded (memo: up to t^2-degree 99); CHANGELOG L15 | correlated-triangle | SMALL | not load-bearing for the rung |
| T6 | P59 L1215 "So the correlated triangle sits ... K3" no gamma | correlated | SMALL | rides with T1 |
| T7 | P57 L170 J9, seam L157, atlas L43 no gamma | correlated | SMALL | |
| T8 | P59 L1216 Levelt / beukers_heckman1989 "from memory, NOT VERIFIED" in memo | citation | route to citation reviewer | |
| T9 | P59 L1211 "at most 2e-26" vs memo 2.3e-26 | | NIT | |
| T10 | P59 abstract L92 "a ring does not climb" generic | coulomb-ring | SMALL | |
| B1 | P59 abstract L90-91 "the whole tree ladder is the banana ladder" vs body n=2..5 | banana | OWNER-LARGE | general-n series argument in memo (verify) |
| B2 | P57 L171 J10 "whole tree ladder"; "operators are prior art" | banana | SMALL (rides with B1 family close) | |
| B3 | atlas L44 door 20 unscoped + "five densities, evidence" stale | banana | SMALL | |
| B4 | seam L158 "relation to four is open" vs same entry "n=5 no new type" | banana | SMALL (self-contradiction) -> consider LARGE? front door | seam is FRONT_DOOR in C25 |
| B5 | density_loop.md L235, CHANGELOG L14, matrix L167 unscoped | banana | SMALL | |
| I1 | P59 L733-735 [OBSERVATION] conductor-4 absent from Bessel-moment corpus "levels 3,6,8,15" vs Lim-Tu-Yu 48.3.e.a = 12.3.c.a x chi_-4 for c_{3,0}, 24.6.a.a | ising / P59 placement | OWNER-LARGE (doubt) | PI item: does the unclaimed-slot conclusion survive? |
| I2 | P59 L1050-1051 "registered cases hold ... every integer-relation search controlled"; CHANGELOG "IS1-IS3 hold as registered" | ising | SMALL | 4 void homogeneous runs; IS3 conifold attachment untested |
| I3 | P59 L1052 "a statement at a stated coefficient bound" -- no bound stated | ising | SMALL | |
| I4 | record ising_class.md L114 "n = 1..8 ... 286 digits" vs n=1 route B 93-101 | ising | SMALL | |
| I5 | P59 L1046 "it is the solution without the logarithm" uniqueness overstated; Observation under SYMBOLIC tag | ising | SMALL | |
| I6 | P59 L1061 n=5,7 determinants; omits n=8 | ising | NIT | |

## Instrument
- p59-ising-thresholds overclaim "for (every|all) n ... singular" matches owner L1042-1043 across a line break: tighten with up-to/le-8 lookahead.
- p59-ising-thresholds qualifier too permissive (threshold|Lim).
- p59-tree-banana-hadamard overclaim misses "whole (tree )?ladder".
- New patterns from reviewers to add before close.

## Wave 2 (ring claim search + P57/seam/atlas owner read), PM-verified
| # | Locus | Grade | Note |
|---|---|---|---|
| S1 | P57 abstract l.70 "40 forced entries ... " (sums 73) vs body 41/74 at l.185/685/694; omits J10/J11 classes | OWNER-LARGE (headline count wrong at the catalogue owner's front door; introduced by v5.20.21 edit omission) | test_atlas_register checks only seam .md |
| S2 | P57 l.159 "the one-centre limit of every tree integral below is rational" (new in v5.20.21) | SURFACE-LARGE | P59 shows chain only; four-density one-centre: "closed forms" (Fromm-Hill, Padhy) |
| S3 | atlas l.17 "A ring does not climb: the triangle of three densities sits on the elliptic rung" no Coulomb; "one centre, rationals and pi" no ring exception | SURFACE-LARGE (G1 family) | |
| S4 | field guide l.482-486 and l.734-737 "a second kernel sharing a vertex lifts the curve to a surface ... observed pattern across the tested cases" -- J11 is a tested counter-case; no trees scope | SURFACE-LARGE | summary surface; ring never mentioned |
| S5 | P57 l.291 + seam l.177 "every forced entry except C9 at MF=1" vs l.185 four exceptions (C9, J7, J10, J11) | SMALL | stale since J7 |
| S6 | seam l.158 J10 "relation to four is open" vs same row "n=5 no new type" (=B4) | SMALL | |
| S7 | INDEX l.69 "60-entry" catalogue | SMALL | pre-existing |
| S8 | group2 synthesis l.808-822 omits ring/correlated triangle | SMALL | |
| S9 | density_loop.md caps l.184 "geminal triangle untested" (QG tested it) | SMALL | |
| S10 | P59 l.1355-1356 "the rung is set by how many densities one integral links"; l.917 "a ring behaves differently" | SMALL | |
| S11 | field guide l.747 "Two rows were pre-registered" | SMALL | |
| U1 | seam J11 l.159 "within the searched maps ... does not split" UNDERCLAIM (owner: proved) | SMALL (two-way) | |
| N  | P57 l.394 36/60 vs 34/60 labellings; matrix l.164 rests-on; seam F definition vs J rows | NIT / PI | |

## Wave 3 (Paper 59 owner read in full), PM-verified
- Converges with T1 (OWNER-LARGE mechanism), B1 (abstract "whole"), T3 (gamma scope; reviewer grades abstract+conclusion LARGE -> accept SURFACE-LARGE), S10 (conclusion count rule l.1355 -> upgrade to SURFACE-LARGE).
| # | Locus | Grade |
|---|---|---|
| P1 | abstract l.92 "a ring does not climb" (only the triangle; QR STOP; gamma=3/2 climb not excluded) | SURFACE-LARGE (merges T10) |
| P2 | conclusion l.1355 "the rung is set by how many densities one integral links" | SURFACE-LARGE (was S10) |
| P3 | abstract/conclusion correlated triangle no gamma scope | SURFACE-LARGE (merges T3) |
| P4 | abstract l.93 "splits" vs factorization | SMALL |
| P5 | l.1233 "a ring of three does not sit with the chain ... but lower" no Coulomb | SMALL (rides with P1 family) |
| P6 | l.1057 c80 bounds absent; n=8 determinant dropped (=I6) | SMALL |
| P7 | l.1051 bound not stated (=I3); c60 precision absent | SMALL |
| P8 | l.1046 "the solution" (=I5) | SMALL |
| P9 | new transcendentals (pi^3, Li2, ln, L(f,2), L(g,4), Gamma(1/3)) untagged vs Paper 18/34 scheme | SMALL |
| P10 | l.834-838 fifth registration omits ring of four + its STOP | SMALL |
| NITs | l.1099 1e-12 vs 1.2e-12; l.1080 uncited; l.975 MUM; l.1199 falsified not said; remiddi1991 no title; Bourjaily 2(L-1) "bound" vs conjecture (pre-existing, v5.20.20) |

## Wave 4 (citations), PM-checked
- No fabricated IDs, no LARGE. 17 keys checked against primary text.
| C1 | lim_tu_yu2026 credit: LTY fn.3 credit AMNT26 + Zagier for the c30/c40 L-value identifications; only c60 numerically theirs | SMALL (credit given away, but to the wrong party) |
| C2 | bonisch2021 title "integrals" not "amplitudes" (JHEP version) | SMALL |
| C3 | vanhove2014: W_N period statement not found in pp. 22-28 | UNVERIFIABLE -> SMALL (narrow or re-source) |
| C4 | cydb "same Yukawa coupling" | UNVERIFIABLE (DB query failed) |
| C5 | beukers_heckman1989 metadata CONFIRMED (secondary); Levelt needs 3 genuine points: memo l.161 shows N3 singular only at t=0,+-2,inf + apparent -> argument holds; paper should state it | SMALL |
| C6 | ppr2004 journal ref missing | NIT |

## Wave 5 (code), PM-checked
- No invalid backing. Both files pass --slow (14 in 81.7 s; 16 in 192.5 s). Independent re-derivations all positive (LMFDB a_p via API; c40 from the eta product to 1.7e-44; n = 6 thresholds by recurrence guessing; ring factorisation to O(t^39); chain at new t; no-split lemma recomputed). Fire tests 16/16 meaningful (a_59 plant insensitive: NIT).
| K1 | P59 l.1213 "no operator of order three annihilates Y" true (split search: rank 48 < 49) but unguarded | SMALL |
| K2 | P59 l.1154 "pieces ... elementary, algebraic and elliptic" + abstract "splits": L_G has NO split (rank 7 < 8); G+1 is an elliptic period integrated against an algebraic kernel | SMALL (abstract word rides with P1/P3 fix) |
| K3 | stated precisions stronger than guarded (1e-28, 1e-12, 18 digits, 2e-26, 5e-261); several claims untested | SMALL |
| K4 | P59 l.1057 c80 conclusion rests on untested bounded negative | SMALL (= P6) |
| UP | Fuchs excludes every factor for F_6, F_7, F_8 -> the n <= 8 threshold statement is about the minimal operator | upgrade candidate |

## Close
- C25 entries extended with the reviewers' patterns; discrimination: fire on every verified locus, silent on corrected wording (0 false +); selftest PASS. New hit: field guide l.760 "relation to four is open" (SMALL, B4 family).
- C25 now FAILs on 4 loci (owner abstract l.92, INDEX:104 x2, seam l.158) until remediation; tests/test_claim_registry.py::test_live_corpus_passes red until then.

## Verdict
DEFECTS. RESEARCH-GRADE not earned. OWNER-LARGE 2 (T1 mechanism sentence; I1 conductor-4 observation). SURFACE-LARGE 10 (G1 INDEX ring; B1 abstract whole ladder; S1 P57 abstract count; T4 P57 J11 gamma=3/2 cap; S2 P57 every-tree rational; S3 atlas preamble ring; S4 field guide rung rule; P1 abstract ring; P2 conclusion count rule; P3 abstract/conclusion gamma scope). SMALL ~35.

## Remediation (2026-10-18)
- All 12 LARGE fixed gate-first (8 C16 entries, 30 loci). Self-check loop: six fresh Opus rounds, LARGE 4, 5, 3, 2, 1, 0 (16 further C16 entries, 29 loci). SMALLs: docs/qa/v5.20.21_claim_delta.debt.md. CHANGELOG v5.20.22.
