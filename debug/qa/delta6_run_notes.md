# Delta-6 verification run (2026-08-30)

**Verdict: DEFECTS.** The full certifying run stays blocked.

## Calibration: SPLIT

| dim | seeds | result |
|:--|:--|:--|
| claims | S1 reinstated QWC direction; S2 unbacked convergence; S5 synthesis "essentially unchanged" | **3/3 CAUGHT — calibrated** |
| code | S3 vacuous guard; S4 reverted physics fix | **1/2 — DE-CALIBRATED** |

The code agent caught S4 thoroughly -- three independent verifications, and
it monkeypatched the fix to confirm the expected values -- but **missed S3**,
a vacuously-true assertion in the very file it was auditing, despite an
explicit instruction to sweep for tautologies. Per the tiering rule that
de-calibrates the dimension: its "clean" on the tautology axis carries no
information. **Re-dispatched at Opus tier** with a seed set weighted toward
the axis that failed (a vacuous comparison, a tolerance that admits the
known-wrong value, and a weakened structural assertion), and with the
reverted-fix seed deliberately NOT re-planted, since it dominated attention
last time.

Both of the code dimension's MATERIAL findings were S4 or downstream of it;
the real corpus was verified to carry the corrected loop throughout.

## The pattern, third run running: propagation, not arithmetic

Every recomputed number checked out again. The claims reviewer independently
re-derived `eq:eri_union` (|P_k| = 7/12/9, 59+48 = 107), the pure-d
sum-of-squares (85), the `eq:alpha_exact` identity, and the four-point
lambda fit. **What keeps failing is reach**, and this run found the sharpest
instance yet.

### The best catch: a table inconsistent with the exponent fitted from it

`tab:onenorm`'s lambda column still held retired values at Q=60 and Q=110
(261.57, 657.07 against the canonical 275.718, 790.007). The proof is the
elegant part: **fitting the printed values gives slope 1.6944 -- the retired
exponent exactly -- while the canonical four give 1.7737, the 1.774 the
paper states two pages later.** The table and the fit derived from it had
come apart, and the arithmetic identified which side was stale. After the
fix the printed values re-fit to 1.7737, matching the stated value.

### The most consequential: a false load-bearing mechanism, asserted 3x

The synthesis said, in its abstract, its discipline bullet, and its "What is
robust" paragraph: *"the correction changes ERI values, not which entries
are nonzero."* It changed exactly which entries are nonzero, 65 -> 107 per
block. If only values had moved, no count could have moved -- yet every
count moved, as the same document says elsewhere.

The underlying confusion is worth naming, because it explains how a false
sentence became load-bearing. **Two defects were corrected together and they
behave oppositely:** the q-sign error changed the SUPPORT (65 -> 107, which
is why every count and exponent moved), while the factor-order error leaves
support exactly invariant and flips signs only (measured on the relativistic
block: identical 2,676 entries, identical sum|ERI|, 40.4% of signs
reversed). "Values, not entries" is true of the second and false of the
first; the synthesis generalised it to the whole correction. The two are now
separated explicitly.

The same paragraph also listed `1-norm O(Q^1.69) (retired-rule vintage)` and
`QWC O(Q^3.36)` under **"What is robust"** while stating two sentences later
that the exponents had moved -- a paragraph contradicting itself about its
own subject.

### The backing gap, worse than the reviewer could see

`tests/test_paper14_scaling.py` still asserted the retired 3.15 / 1.69 /
3.36, with a QWC gate of [3.1, 3.6] against a live 4.013. The reviewer
flagged it as "either failing or exercising a retired path". It was neither:
**every test in the file is slow-marked, so it ran "4 skipped" and was
silent.** The corpus's headline exponents had no active backing at all --
precisely what the corrected criteria name as MATERIAL. Re-centred on the
measured values, and a structural assertion added that QWC exceeds the Pauli
exponent, since the "no measurement-group advantage" claim rests on that
ordering.

## Remediation applied

- `tab:onenorm` lambda cells (Q=60, Q=110) + the prose lambda at Q=110
- `tab:pauli` He rows, which disagreed with `tab:qwc` on the same systems
  (14,211/252,095 vs the measured 14,078/250,402), plus one propagated locus
- the n_max=5 count, which existed in **three mutually exclusive states**
  (227,338 live at three loci, "pending re-measurement" at a fourth, 2.4
  million at two more) -- unified to the measured 2,434,441
- `tab:spinor_resource` observations 2 and 3, which cited retired lambda and
  QWC values *and sourced them to the corrected table*, and which kept alive
  the matched-Q comparator the caption says was removed
- P20's abstract "one to three orders of magnitude", contradicted two
  sentences later by its own "near parity"; replaced with the honest span
  (0.28x to 317x, depending on the comparison)
- the groups-per-term range and the inference drawn from it: the ratio
  *rises* with basis, so grouping gets less efficient, not more
- the M^-0.85 ERI-density family at four loci (live: M^-0.49), including a
  figure caption that plotted an M^-0.49 fit and described it as M^-0.85
- two inverted provenance tags, one understated survived advantage
  (6.8x -> the live 7.2x), and the commutator-bound subtraction (0.22 was
  1.69 - 1.47; against the live 1.774 it is 0.30)
- the "favorable accuracy-per-Pauli-term ratio" conclusion, which the
  paper's own sparsity section now contradicts

## Upgrade accepted

The reviewer improved one of my own arguments. I had justified the 55 direct
Pauli terms' convention-independence by the k=0 monopole coefficient, which
does not by itself establish non-vanishing -- the Slater sum carries both
signs at k>0. The stronger and simpler argument is positivity:
<ab|ab> = int int |phi_a|^2 |phi_b|^2 / r_12 > 0 strictly, for any basis and
any convention. To be applied.

## Verification

Gates 7/7 on group4 and group6. P14 30 pp, P20 12 pp, synthesis 5 pp, all
compiling with zero undefined refs/citations/control sequences.

---

## Delta-6b: code dimension re-dispatched at Opus — 3/3, RE-CALIBRATED

The Sonnet pass missed a vacuously-true assertion in the file it was
auditing. Re-dispatched at Opus with the seed set weighted toward that axis
(and the reverted-fix seed deliberately withheld, since it had dominated
attention).

| seed | class | result |
|:--|:--|:--|
| T1 | vacuous `agree <= total` | **CAUGHT** |
| T2 | tolerances admitting the retired ordering | **CAUGHT**, and quantified |
| T3 | weakened structural assertion | **CAUGHT** |

**3/3 — the tier upgrade was the fix, and the calibration net is what
detected the need for it.** The agent enumerated 17 assertions in one file
and 7 in the other, stating the counts, rather than sampling.

### Genuine findings against my own work

**M3 -- the discriminator is SCALE-BLIND, and its docstring over-credits
it.** `_compare` least-squares-fits `scale` and measures the residual only
AFTER rescaling, so a uniform rescale of X_k passes everything. The agent
demonstrated it by mutating `jj_angular_Xk -> X/2`: signs 168/168, scale
4.0000, residual 0.0000, at every multipole. The k=0 control cannot catch it
either, because **X_0 is exactly the identity matrix**, so k=0 pins no 3j or
reduced-matrix-element normalization at all.

This does not weaken the ordering verdict -- sign patterns are what
discriminate order, and a positive rescale cannot flip a sign -- but the
header credited this test with the machine-precision agreement, an
attribution the assertions did not support once `scale` was free. Fixed by
asserting `scale == 1` tightly at every multipole (free: the corrected
ordering already gives 1.0000) and by stating the scope honestly: this is a
**cross-basis consistency test**; the physical anchoring of the scalar
convention lives outside the file. Verified the new guard bites: with
`X -> X/2` it fires at scale 4.0000.

**M4 -- a band that admitted the value named in its own docstring.** The
1-norm test asserted `1.65 <= alpha <= 1.95` while its docstring said
"Retired pair-diagonal value 1.69" -- which is *inside the band*. Zero
discrimination: a regression restoring the q-sign error would put lambda
back at 1.69 and this test would pass while both siblings fired. Tightened
to [1.72, 1.90], which excludes 1.69 and admits the live 1.774/1.792.

**M5 -- an unreachable assertion labelled "the key fault-tolerant claim".**
`assert alpha < 2.0` cannot fail given the band above it. Replaced with the
comparison the argument actually rests on: lambda grows strictly slower than
the term count.

**NIT with teeth -- `casimir_ci`'s own docstring stated the RETIRED factor
order**, contradicting its corrected code 48 lines below. That is the module
the discriminator builds its reference from, so a reader checking the
reference against the docstring would have concluded the reference was
wrong. Fixed.

### A self-correction the new fast pin forced

Adding a default-run pin surfaced that `tab:pauli` uses the
INCLUDING-identity convention (Q=10 row reads 288 = 287 + 1). My delta-6 fix
had re-priced its Q=28 and Q=60 rows to NON-identity values, leaving the
table in a mixed convention. The retired 14,211 / 252,095 were still
genuinely stale -- neither convention -- so the finding stood; only my
replacement values were wrong. Corrected to 14,079 / 250,403. **The test
disagreeing with my number is what caught it**, which is the argument for
having written it.

### Coverage restored

The re-centred exponent tests PASS under `--slow` (4 passed, 70 min). But 70
minutes means they run only under `--slow`, leaving the headline exponents
with no default-run protection -- the same shape as the propinquity gap
found earlier in this arc. Added `test_small_basis_points_are_pinned`: 5.8 s,
runs by default, pins the two cheapest MEASURED points the fit is built from
(Q=10: 287 / 11.175; Q=28: 14,078 / 74.207) and thereby locks out a silent
revert to the retired rule.

### Verification

Gates 7/7 on group4 and group6. 15 passed in the discriminator + eri-rule
set; 13 passed + 4 slow-skipped in the ordering + scaling set.

Seeds were worktree-only throughout; the real corpus was verified clean of
all six delta-6/6b seed strings. The worktree's git registration is removed
and its contents deleted (0 files); the empty directory itself is held by an
OS handle and could not be unlinked -- it contains nothing.

### Remaining NITs cleared, and two of them had teeth

The Opus agent re-notified with one correction to its own report: the
`--slow` run **did** complete (4 passed, 4144 s), not blocked -- its earlier
"did not complete" was its own timeout. That independently reproduces the
69-minute figure measured here (4175 s). Everything else in its report
stood.

Six NITs cleared, two of which were more than cosmetic:

**The mask support was never asserted.** `agree == total` is the decisive
statistic, but nothing pinned `total`. A *fully* collapsed mask raises
loudly, so that case was safe -- but a **partial** collapse is silent, and
demonstrated: restricting to the largest entries drops the support 168 -> 48
and `agree == total` still reads True. The discriminator would report
perfect agreement while comparing a quarter of the tensor. Now
`EXPECTED_SUPPORT = {(2p,2): 168, (3d,2): 1004, (3d,4): 788}` is asserted
first, and it fires on the collapse.

**The production-source guard was cwd-relative.** `open("geovac/...")` can
be defeated by running pytest from another directory; anchored to
`R.__file__`.

Also: the module docstring claimed "~20 s" against a measured 69 min (the
estimate predated the correction, which multiplied the term counts); all six
n_max=1,2 rows of `tab:spinor_resource` are now pinned individually rather
than LiH alone, so each independently locks out a revert; the k=0 control's
docstring no longer claims that its passing means "the comparison machinery
is sound" -- at k=0 the X matrix is exactly the identity, so it is a null
control and nothing more; a redundant bound and two stray escape sequences
removed.

The ordering + scaling suite goes 12 -> **18 passed**.

### Delta-6 close

Gates 7/7 on group4 and group6. **80 passed, 5 skipped** across the
discriminator, scaling, eri-rule, heavy-hydride, spin-ful and topological
sets -- including the 18 symbolic S^3 proofs.

Verdict stands at **DEFECTS** for delta-6 (the run found real defects and
they were remediated). Calibration: claims 3/3; code 1/2 at Sonnet, then
**3/3 on Opus re-dispatch**. The tiering rule did exactly what it exists for
-- the miss was detected by the seed, not by luck, and the re-dispatch
converted a de-calibrated dimension into a calibrated one.

A further delta is owed over THIS remediation before any certifying run.
