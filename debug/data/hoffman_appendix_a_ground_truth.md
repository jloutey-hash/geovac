# PM ground truth: Hoffman Appendix A vs GeoVac pinned t-value identifications

Independent verification artifact (PM main session, 2026-06-11), extracted
directly from the primary source PDF (`debug/data/hoffman_1612.05232.pdf`,
arXiv:1612.05232v5) BEFORE reading the T2 literature-grounding agent's memo.
Purpose: hallucination-proof diff base for the agent's citations.

## Source verification (fetched from arxiv.org directly)

- Hoffman, "An Odd Variant of Multiple Zeta Values", Commun. Number Theory
  Phys. 13 (2019) 529-567; arXiv:1612.05232 [math.NT] (v5 dated 13 Oct 2020;
  paper dated January 29, 2019). REAL, title/venue exact.
- Bluemlein-Broadhurst-Vermaseren, "The Multiple Zeta Value Data Mine",
  Comput. Phys. Commun. 181 (2010) 582-625; arXiv:0907.2557. REAL.
  Euler sums (level 2) explicitly reduced through weight 12.
- In Hoffman's own bibliography: [17] Kaneko-Tasaka Math. Ann. 357? (printed
  "307 (2013), 1091-1118"); [18] Kaneko-Tsumura "On a variant of multiple
  zeta values of level two" (preprint at his writing); [23] Nakamura-Tasaka
  "Remarks on double zeta values of level 2", J. Number Theory 133 (2013)
  48-54; [26] B. Saha, "A conjecture about multiple t-values",
  arXiv:1712.06325; [5] Broadhurst-Borwein-Bradley, Electron. J. Combin.
  4(2) (1997) art. 5.

## Convention (his eq. (1.1), page 2)

t(i1,...,ik) = sum over n1 > n2 > ... > nk >= 1, all n odd, of
1/(n1^i1 ... nk^ik); i1 > 1; FIRST index on LARGEST variable; NO
normalization factor; t(i) = (1-2^-i) zeta(i). IDENTICAL to GeoVac's
t2/t3 convention (= our lambda for depth 1).

## Key structural results

- Thm 3.2: symmetric-sum theorem (the stuffle machinery our R-relations
  instantiate; published ancestor).
- Cor 4.1 + eq (4.4) + eq (4.5): exact conversion of any depth-k t-value
  to 2^(k-1) alternating MZVs; for doubles
  t(a,b) = (1/2 - 2^-(a+b)) zeta(a,b) + (1/2) zeta(abar,bbar).
  => combined with BBV (alternating reductions PROVEN through w=12), every
  GeoVac double-t reduction at w <= 11 is DERIVABLE-grade even where no
  explicit Hoffman table row exists.
- Tables stop at weight 7 (Appendix A, pages 33-36). No parity theorem
  for t-values is stated as such in this paper; the odd-weight depth-2
  reductions through w7 are explicit table rows instead.

## Appendix A rows vs our pins (manual conversion + machine check at 35 dps)

| Object | Hoffman row (his notation) | Our pinned form | Residual |
|---|---|---|---|
| t(2,1) | -1/2 t(3) + t(2) log2 | -(7/16)z3 + (1/8)pi^2 ln2 | exact conversion |
| t(2,3) | -1/2 t(5) + (4/7) t(2)t(3) | -(31/64)z5 + (1/16)pi^2 z3 | exact conversion |
| t(4,1) | -1/2 t(5) - (1/7) t(2)t(3) + t(4) log2 | (1/96)pi^4 ln2 - (1/64)pi^2 z3 - (31/64)z5 | exact conversion |
| t(4,3) | -1/2 t(7) + (6/7) t(3)t(4) - (10/31) t(2)t(5) | (1/128)pi^4 z3 - (5/128)pi^2 z5 - (127/256)z7 | 7.5e-37 |
| t(5,1) | -(73/84)t(6) + (17/98)t(3)^2 - (1/2)zeta(5bar,1) + t(5)log2 | direct sum t2(5,1) | 1.3e-36 |
| t(4,1,1) | (45/112)t(6) - (31/196)t(3)^2 + (1/4)zeta(5bar,1) - (1/2)t(5)log2 - (1/7)t(2)t(3)log2 + (1/2)t(4)log^2 2 | our t3(4,1,1) pin incl. -(1/2)t2(5,1) | 1.0e-36 |
| t(3,2,1) | (37/112)t(6) - (33/98)t(3)^2 + (1/4)zeta(5bar,1) - (1/2)t(5)log2 + (3/7)t(2)t(3)log2 | our t3(3,2,1) pin incl. -(1/2)t2(5,1) | 1.1e-36 |
| t(2,1,1) | (11/60)t(4) + (1/4)zeta(3bar,1) - (1/2)t(3)log2 + (1/2)t(2)log^2 2 | (1/2)Li4(1/2) + (1/48)ln^4 2 + (1/24)pi^2 ln^2 2 - (19/5760)pi^4 | 0.0 |
| a6 = t(4,2)-t(2,4) | rows t(4,2) = -(1/7)t(6)+(1/7)t(3)^2, t(2,4) = (11/28)t(6)-(1/7)t(3)^2 | -pi^6/1792 + (7/32)z3^2 | exact conversion |
| t(2,2) | (1/4) t(4) | stuffle (lam2^2-lam4)/2 = pi^4/384 | exact conversion |

Also explicit in Appendix A (w7 doubles): t(6,1), t(5,2), t(3,4), t(2,5)
=> our stage-A t2(6,1), t2(2,5) identifications are table rows (MATCHED).

## Verdicts pre-staged for the T2 diff

1. ALL GeoVac w<=7 pinned identifications are verbatim rows of Hoffman
   Appendix A (after eliminating zeta(5bar,1) via his own t(5,1) row):
   MATCHED grade. Our PSLQ runs rediscovered the published table.
2. t2(5,1) generator status: SUPPORTED (Hoffman represents it via the
   alternating MZV zeta(5bar,1); no product reduction exists in his table).
3. w8-w11 items (t3(3,4,1), t2(6,2), t2(2,6), w9/w11 odd reductions,
   t2(5,3)/t2(7,1) as generators): beyond Hoffman's tables; grade
   DERIVABLE via Cor 4.1 conversion + BBV data mine (proven w<=12).
4. No "parity theorem for multiple t-values" exists IN THIS PAPER; any
   such citation by the agent must point elsewhere (and be verified).

Machine checks: driver inline (see session log 2026-06-11); residuals
column above at mpmath 35 dps.
