# Trunk SMALL-debt register

SMALL findings that could not be fixed on sight (criteria.md RESEARCH-GRADE).

- 2026-09-29 (DELTA #15) — `paper_34_projection_taxonomy.tex` ~L3947–4769: an 820-line `[h]` float table sits at TeX's maximum dimension ("Dimension too large" at `\end{table}` once ~2 lines are added). Blocked on a structural edit (convert to `longtable`); not a content defect.
- 2026-09-30 (FULL #11) — `paper_32_spectral_triple.tex`: 19 inline `Author~Year` attributions resolve to no `\bibitem` (C20 reports them; they sit inside its ratchet baseline, so the gate passes). Blocked on primary-source verification of each before a bibitem is added — adding unverified metadata would launder a guess into a citation (criteria C4). List: `python debug/qa/check_inline_attributions.py --gate trunk`.
