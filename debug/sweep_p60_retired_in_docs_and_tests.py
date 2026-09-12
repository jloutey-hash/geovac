"""Sweep retired Paper-60 values out of the claim registry and the test prose.

Found while reviewing the commit: geovac/sturmian_secular.py's docstring carried
TWO retired values, which prompted a wider look.  docs/claim_test_matrix.md --
the corpus's own claim->artifact registry, and a goalpost the next certifying
run measures against -- still ratifies four of them as LIVE:

  * K^0.84                    -> retired; converged window fit is K^0.82
  * "local slope RISES ... 0.906"  -> retired; converged slope FALLS to 0.766
  * "T' is SUPERlinear (K^1.05)"   -> retired; T' full is K^0.88, SUBLINEAR
  * "cond(S) 4 -> 3673" as evidence -> retired; box artifact, converged 16.0

tests/test_sturmian_secular.py's ASSERTIONS were correctly repaired on
2026-09-07 (test_L2_metric_conditioning_is_box_dependent documents the artifact
properly), but its module docstring and two comment blocks were never swept.
No assertion is touched here -- prose only.
"""
import io

MTX = "docs/claim_test_matrix.md"
TST = "tests/test_sturmian_secular.py"

# ------------------------------------------------------ claim_test_matrix.md
m = io.open(MTX, encoding="utf-8").read()
subs = [
 # header block
 ("**K^0.84** (`test_sturmian_secular.py`, measured 0.842 over K≤100)",
  "**K^0.82** (`test_sturmian_secular.py` asserts a sublinear band; the exact "
  "window fit is 0.82 on a CONVERGED radial box — 0.842/K^0.84 is RETIRED, "
  "measured on an unconverged 60-bohr box)"),
 ("The L²-divergence\n> (cond(S) 4.1→3672.6), eq:secular",
  "The cond(S) sequence\n> (4.1→3672.6 — RETIRED as a radial-box artifact 2026-09-07; converged "
  "cond(S)=16.0, growing ~0.12·K), eq:secular"),
 # eq:sublinear row
 ("full basis K^0.84 fitted on K=74..164)",
  "full basis K^0.82 fitted on K=74..164, converged box)"),
 ("**Not an asymptotic regime** (/qa 2026-09-07): the local slope rises monotonically past the fitted window — 0.850, 0.868, 0.882, 0.906 at K=202, 244, 290, 340 — so the exponent trends toward 1.",
  "**Not an asymptotic regime**: the local slope FALLS monotonically on a converged box — 0.827→0.766 across K=100..514 — so no window fit is stable. (The 2026-09-07 reading that it *rises* to 0.906 is RETIRED: that rise was the growth of a 60-bohr box's error across the fit range, not the matrix.)"),
 # eq:sublinear_split row
 ("eq:sublinear_split — the sublinearity is carried by the NUCLEAR DIAGONAL T⁰=Z R_ν (K^0.70, stable), not by the pure-number block T′, which is SUPERlinear (K^1.05); the total exponent drifts up because T′'s share of the norm grows (64.6%→52.2%)",
  "eq:sublinear_split — the sublinearity is carried by the NUCLEAR DIAGONAL T⁰=Z R_ν, which is not a power law at all but Z√(2K)·ln(K/2), asymptotic exponent 1/2. **RETIRED 2026-09-08:** the companion claim that T′ is SUPERlinear at K^1.05 — on a converged box T′ is K^0.88 full / K^0.94 off-diagonal, i.e. SUBLINEAR, and K^0.70 for T⁰ is merely what a log-log fit returns inside K=74..164"),
 # l>0 row
 ("(full s+p+d+f K^0.84; converges to −2.897, residual ~7 mHa = basis incompleteness)",
  "(full s+p+d+f K^0.82; converges to −2.897, residual ~7 mHa — **NOT basis incompleteness**: a variational CI over the identical span reaches 1.28 mHa, so the residual is the scale lock, Paper 60 eq:scale_lock)"),
 ("(exact K^0.84 over the K=74..164 window — NOT an asymptote, /qa 2026-09-07; cond(S) 4→3673 collapse)",
  "(exact K^0.82 over the K=74..164 window on a converged box — NOT an asymptote; the cond(S) 4→3673 collapse is RETIRED as a box artifact)"),
 ("consistent with the paper's 0.78 s-only → 0.84 windowed",
  "consistent with the paper's 0.73 s-only → 0.82 windowed"),
 ("measures p=**0.842** over K={9..100} (28 s), via tracked `geovac/sturmian_secular.py`",
  "measures p=**0.842** over K={9..100} (28 s) on the module's default box, via tracked "
  "`geovac/sturmian_secular.py`. That value is box-limited; the converged window fit is **0.82**"),
 # obstruction row -- its whole verdict rested on the artifact
 ("| 60 | §obstruction — restoring the L² metric (generalized eigenproblem MB=pSB) on the converged spdf He diverges as cond(S) climbs 4→3673 |",
  "| 60 | §obstruction — **WITHDRAWN 2026-09-07/09-08.** The L² metric does NOT diverge: cond(S) 4→3673 is a radial-box artifact (converged 16.0, growing ~0.12·K, an ordinary Gram matrix). The real obstruction to freeing the metric is eq:scale_lock, not conditioning |"),
 ("**LIVE (engine sprint 2026-08-18)** | BACKED-SOUND: cond(S) sequence 4.1→5.7→9.8→**3672.6** (N=3,4,6,8) reproduced + monotone-growth assert; corroborates that the metric-free posing is the only well-conditioned one.",
  "**WITHDRAWN** | The test was REPLACED 2026-09-07 by `test_L2_metric_conditioning_is_box_dependent`, which pins the artifact as an artifact and asserts the converged value; the old test asserted `2000 < cond < 6000` and would have FAILED on repair. The retired verdict — “the metric-free posing is the only well-conditioned one” — does not survive: it is the only *metric-free* one, which is a different and now-derived statement (eq:scale_lock)."),
]
for old, new in subs:
    assert old in m, "matrix locus not found: %.70s" % old
    m = m.replace(old, new, 1)
io.open(MTX, "w", encoding="utf-8").write(m)
print("claim_test_matrix.md: %d retired-value loci swept" % len(subs))

# --------------------------------------------------- tests/ prose (NO asserts)
t = io.open(TST, encoding="utf-8").read()
tsubs = [
 ('        approaching ``~K^0.84`` for the full s+p+d+f basis (Paper 60 ``eq:sublinear``);',
  '        approaching ``~K^0.82`` over the window K = 74..164 on a CONVERGED radial\n'
  '        box (Paper 60 ``eq:sublinear``).  ``K^0.84`` is RETIRED -- a 60-bohr-box\n'
  '        artifact.  This file measures 0.842 on the module default box;  that is\n'
  '        box-limited, which is why the assertion below is a BAND, not a value;'),
 ('        is ill-conditioned -- ``cond(S)`` grows from ~4 into the thousands',
  '        is NOT ill-conditioned -- the "~4 into the thousands" growth is a\n'
  '        radial-box artifact (converged cond(S) = 16.0, growing ~0.12*K)'),
 ('# rises with basis size from ~0.78-0.80 at small K toward the full-K asymptote ~0.84\n'
  '# (K~164, s8p8d8f8 + N=10 in the driver).  We assert a band bracketing the measured\n'
  '# value, not the asymptote itself.',
  '# is box-dependent.  There is NO asymptote: on a converged box the local slope FALLS\n'
  '# monotonically (0.827 -> 0.766 across K=100..514), so no window fit is stable.  We\n'
  '# assert a BAND bracketing the measured value, never an asymptote.'),
 ('    # Defensible band bracketing the measured value (~0.84 at K=100; asymptote ~0.84).',
  '    # Defensible band bracketing the measured value (0.842 on this box; converged\n'
  '    # window fit 0.82).  Deliberately a band: there is no asymptote to pin.'),
]
for old, new in tsubs:
    assert old in t, "test-prose locus not found: %.60s" % old
    t = t.replace(old, new, 1)
io.open(TST, "w", encoding="utf-8").write(t)
print("tests/test_sturmian_secular.py: %d prose loci swept (0 assertions touched)" % len(tsubs))
