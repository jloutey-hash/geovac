"""Items 1-3 docs: amend the v5.11.0 Prop-A claim, update the walls register,
add the matrix rows, close the C23 owed table, and write v5.11.2."""
from pathlib import Path

# ---- item 1: amend the refuted claim in place (Sec. 13.11 rule 9) -----------
C = Path("CHANGELOG.md")
s = C.read_text(encoding="utf-8")

OLD = ("- **The overcompleteness wall is a theorem, and an elementary one.** If `g != 0` lies in "
       "the closed span of `{f_i}`, then `lam_min(G_N) <= dist(g, V_N)^2 -> 0`. Completeness of "
       "the one-centre set alone forces it; \"translate\" is incidental. Ron-Shen fiberization "
       "gives the operator form: Riesz sequence iff `ess inf (1 - |sigma|) > 0` iff "
       "`||sigma||_inf < 1`, and ours is exactly 1, attained.")
NEW = ("- **The overcompleteness wall is a theorem** — via Ron-Shen fiberization: Riesz sequence "
       "iff `ess inf (1 - |sigma|) > 0` iff `||sigma||_inf < 1`, and ours is exactly 1, attained. "
       "~~The elementary route — if `g != 0` lies in the closed span of `{f_i}` then "
       "`lam_min(G_N) <= dist(g, V_N)^2 -> 0`, so completeness of the one-centre set alone forces "
       "it~~ **[CORRECTED 2026-09-12, v5.11.2: that hypothesis is FALSE for this basis.** The "
       "Bessel deficit `eps_N^2 = 1 - sum_i <chi^A_i, chi^B_1>^2` plateaus at 0.380 / 0.696 / "
       "0.907 for `kR = 1/2/4`, flat over `N = 16..256` — the one-centre set is far from complete "
       "in the MOLECULAR metric, so the bound holds only vacuously and the degeneracy is not "
       "\"one span already contains the other\". Ron-Shen is the mechanism; it says something "
       "narrower and more useful, namely that the near-dependence is ONE DIRECTION. See v5.11.2.]")
assert s.count(OLD) == 1, "v5.11.0 Prop-A passage not found"
s = s.replace(OLD, NEW)

# ---- item 1+2+3: the v5.11.2 entry -----------------------------------------
ENTRY = """## [v5.11.2] - 2026-09-12

**The overcompleteness mechanism was wrong, and the corrected version is more useful: the near-dependence is ONE DIRECTION, not a property of the basis.** Plus the C23 owed-citation table closed against primaries. Probe `debug/p60_completeness_hypothesis_probe.py`; backing `tests/test_paper60_one_direction.py`.

### The measurement that inverts the reading

The v5.11.0 entry adopted a frames-theoretic mechanism: completeness of the ONE-CENTRE set forces `lam_min -> 0`, with the second centre incidental. That hypothesis is directly measurable here — the SW intra-centre block is exactly the identity, so the `A`-set is orthonormal in that metric and Bessel's inequality gives the deficit of a displaced Sturmian against the one-centre span as `eps_N^2 = 1 - sum_{i<=N} <chi^A_i, chi^B_1>^2`.

**It does not go to zero.** It plateaus at `0.380 / 0.696 / 0.907` for `kR = 1/2/4`, flat from `N = 16` to `N = 256`. Between 38% and 91% of a displaced basis function lies outside the one-centre span. So the hypothesis fails, the bound holds only vacuously (`1e-4 <= 0.70` says nothing), and **"overcompleteness is the price of one-centre completeness" is withdrawn**. It never reached the paper — only the v5.11.0 entry, now amended in place.

**Ron-Shen was always the better route and it survives**: `||sigma||_inf = 1` because `j0(0) = 1`. And it says something narrower and more useful. `sigma_max -> 1` asserts that *some combination* in the `B`-span is captured by the `A`-span — not that individual functions are, and measurably they are not. **The near-dependence is one direction, not a diffuse property of the basis.** That reconciles every structural fact the arc turned up: the rank-one (rank `M-1`) degeneracy, the fixed geometry-independent null direction, the small participation ratio, and why one rotation plus a tridiagonal preconditioner fixes it so cleanly. Captured as a new `[MEASURED]` paragraph in `sec:molecular`; the paired guard asserts that `1 - sigma_max` collapses as `N^-2` **while** `eps_N^2` does not, which no single-sided bug satisfies.

### Two escapes checked and closed

Asked whether the degeneracy is eliminable by angular surgery:

- **Adding `l>0` cannot help.** The original cross block is a submatrix of the enlarged one and `sigma_max` is non-decreasing under submatrix extension, so the `s`-sector degeneracy survives inside the larger problem. (A 2026-09-10 guess that higher `l` would be better-conditioned was backwards.)
- **Excluding `l=0` does not help either.** `l>=1` functions vanish *at* `p=0`, but a normalized combination can still concentrate *near* it, which is all `sigma_max -> 1` requires. Reasoning, not measurement — flagged as such.

### A caveat on the basis's motivation

The plateau has a second consequence worth stating. Paper 60 motivates Coulomb Sturmians as evading the Gaussian coverage-versus-linear-dependence trade-off "because they are a complete set at one scale". That completeness is in the **atomic** metric; the molecular problem is posed in `V_0`, where the one-centre set is measurably *not* complete (70% deficit at `kR = 2`). The property that sells the basis is not the property the molecular problem uses. Nothing measured contradicts anything, but the motivating sentence is now qualified in place rather than left to be discovered later.

### C23 run #1: the owed table is closed

The WebSearch budget was raised (project `env`, 200 -> 500), which took effect immediately. Ten primaries verified this session and cited; each was checked at source, not taken from the scan:

| now cited | verified how |
|:--|:--|
| Jordan, *Bull. Soc. Math. France* **3**, 103 (1875) | NUMDAM record; principal-angle priority confirmed |
| Björck & Golub, *Math. Comp.* **27**, 579 (1973) | DOI 10.1090/S0025-5718-1973-0348991-3 |
| Eijkhout & Vassilevski, *SIAM Rev.* **33**, 405 (1991) | title/volume/pages |
| Hartman & Wintner, *Amer. J. Math.* **76**, 867 (1954) | title + the content used (self-adjoint Toeplitz spectrum = convex hull of the essential range) |
| Löwdin, *J. Chem. Phys.* **18**, 365 (1950) | read in Slater–Koster's own footnote 12 |
| Slater & Koster, *Phys. Rev.* **94**, 1498 (1954) | **PRIMARY READ** — Sec. II, p. 1500 introduces Löwdin orthogonalization for exactly this multi-centre non-orthogonality problem |
| Rokob, Szabados & Surján | existence + abstract (the symmetry operator's transformation matrix must be unitary; fails for Cartesian `d`/`f`); the ELTE PDF threw a certificate error, so the primary is unread |
| Jaffard, *Ann. IHP C* **7**, 461 (1990) | title/volume/pages |
| Gröchenig & Leinert, *TAMS* **358**, 2695 (2006) | title/journal/volume |
| Driscoll & Fornberg, *Comput. Math. Appl.* **43**, 413 (2002) | title/volume/pages; coined "flat limit" |

**One scan claim was dropped rather than cited**: that Slater–Koster's *Appendix* states the symmetry theorem. What was read is their Sec. II use of Löwdin, which is what the citation now carries — the Appendix-specific assertion is unverified and not relied on. **Still named in prose with no bibitem**: the "Jordan–Wielandt" label (Stewart–Sun and Horn–Johnson are books, unopened).

### C23 scope revised

Run #1's two best catches were its two *newest* claims — both within 48 hours of authorship — while four older ones had survived three DELTA runs and a FULL run. C23 now has **two triggers**: at FULL certification, and **at authorship** for any new claim matching a priority signature (a clean closed-form constant; an external field entered sideways; a derivation under a page). Per-claim, not per-paper. PI-approved.

### Gates

C10 / C21 / C16 / C22 / C14 / escapes / titles / arxiv / duration PASS in scope `paper_60`. New guards fire-tested three ways, including both halves of the paired claim (make the gap not collapse; make the deficit collapse too).

"""
A = "## [v5.11.1] - 2026-09-12"
assert s.count(A) == 1
s = s.replace(A, ENTRY + A)
C.write_text(s, encoding="utf-8")
print("CHANGELOG: v5.11.0 amended in place + v5.11.2 inserted")

# ---- walls register --------------------------------------------------------
W = Path("docs/walls/register.md")
t = W.read_text(encoding="utf-8")
OLDW = ("**Falsifier for the split.** A congruence that is simultaneously (i) `l`-block diagonal "
        "and (ii) orthogonalizing, on a metric with nonzero inter-center coupling")
NEWW = ("**Mechanism correction (2026-09-12, v5.11.2).** The conditioning axis's obstruction is "
        "`||sigma||_inf = 1` attained at `p = 0` (Ron-Shen), **not** completeness of the "
        "one-centre set. The Bessel deficit against the one-centre span plateaus at 0.380 / 0.696 "
        "/ 0.907 for `kR = 1/2/4`, flat over `N = 16..256`, so that set is far from complete in "
        "the molecular metric. The corrected reading is stronger operationally: the "
        "near-dependence is **one direction**, which is why a fixed rank-`M-1` rotation removes "
        "it. Two escapes checked and closed: adding `l>0` cannot help (`sigma_max` is "
        "non-decreasing under submatrix extension), and excluding `l=0` does not either (an "
        "`l>=1` combination can still concentrate near `p=0`) -- the second is reasoning, not "
        "measurement. Backing `tests/test_paper60_one_direction.py`.\n\n"
        "**Falsifier for the split.** A congruence that is simultaneously (i) `l`-block diagonal "
        "and (ii) orthogonalizing, on a metric with nonzero inter-center coupling")
assert t.count(OLDW) == 1, "walls falsifier anchor not found"
W.write_text(t.replace(OLDW, NEWW), encoding="utf-8")
print("walls register: mechanism correction added")

# ---- matrix row ------------------------------------------------------------
M = Path("docs/claim_test_matrix.md")
m = M.read_text(encoding="utf-8")
ROW = ("| 60 | §molecular [MEASURED] — the degeneracy is ONE DIRECTION, not one-centre "
       "completeness: the Bessel deficit `eps_N^2 = 1 - sum_i <chi^A_i,chi^B_1>^2` PLATEAUS at "
       "0.380/0.696/0.907 for kR=1/2/4 (flat N=16..256), so the one-centre set is far from "
       "complete in the molecular metric and `sigma_max -> 1` asserts only that some COMBINATION "
       "is captured. **Withdraws the 2026-09-11 frames reading** (completeness of the one-centre "
       "set forces lam_min -> 0) as inapplicable here | "
       "`tests/test_paper60_one_direction.py`"
       "``::test_one_centre_set_is_not_complete_in_the_molecular_metric`` (3 kR) + "
       "``::test_the_degeneracy_is_one_direction_not_the_whole_basis`` | tracked "
       "`geovac/sturmian_sigma_law.py` | **NEW 2026-09-12** | BACKED-SOUND. The load-bearing form "
       "is PAIRED: `1 - sigma_max` must collapse as N^-2 WHILE `eps_N^2` must not, so neither "
       "single-sided failure passes. Plateau asserted three ways a decaying sequence cannot meet "
       "(value, ratio, fitted slope). Fire-tested three ways: summing the whole cross block; "
       "making the gap not collapse; making the deficit collapse too. rests on: eq:sigma_law |")
A2 = "| 60 | §molecular [SYMBOLIC + MEASURED] eq:chirp_decay"
i = m.index(A2)
M.write_text(m[:i] + ROW + "\n" + m[i:], encoding="utf-8")
print("claim_test_matrix: +1 one-direction row")

# ---- C23 record: close the owed table --------------------------------------
R = Path("docs/qa/c23_run_001_paper_60.md")
r = R.read_text(encoding="utf-8")
OLDR = ("**To close this record:** raise the WebSearch budget, verify the table above, convert "
        "the prose names to bibitems, and re-run C23's verdict on A1/A3/A4/B3.")
NEWR = ("**CLOSED 2026-09-12 (v5.11.2).** The WebSearch budget was raised (project `env`, "
        "200 -> 500) and took effect immediately. Ten of the eleven owed primaries were verified "
        "at source and are now cited: Jordan 1875, Björck–Golub 1973, Eijkhout–Vassilevski 1991, "
        "Hartman–Wintner 1954, Löwdin 1950, Slater–Koster 1954 (**primary read** — Sec. II, "
        "p. 1500), Rokob–Szabados–Surján (existence + abstract; PDF cert-blocked), Jaffard 1990, "
        "Gröchenig–Leinert 2006, Driscoll–Fornberg 2002. **Two items were deliberately not "
        "closed:** the \"Jordan–Wielandt\" label stays prose-only (Stewart–Sun and Horn–Johnson "
        "are books, unopened), and the scan's claim that Slater–Koster's *Appendix* states the "
        "symmetry theorem is **dropped, not cited** — what was read is their Sec. II use of "
        "Löwdin, which is what the citation carries. Böttcher–Spitkovsky remains cited-but-unread "
        "and flagged. A1/A3/A4/B3 verdicts stand as PRIOR ART, now with primaries.")
assert r.count(OLDR) == 1
R.write_text(r.replace(OLDR, NEWR), encoding="utf-8")
print("C23 record: owed table closed")

# ---- CLAUDE.md -------------------------------------------------------------
L = Path("CLAUDE.md")
l = L.read_text(encoding="utf-8")
old = "**Version:** v5.11.1 (September 12, 2026)"
assert l.count(old) == 1
l = l.replace(old, "**Version:** v5.11.2 (September 12, 2026)")
bullet = ("- **Overcompleteness is ONE DIRECTION (2026-09-12, v5.11.2):** Bessel deficit plateaus "
          "0.38-0.91, so the one-centre set is NOT complete in the molecular metric; the frames "
          "reading is withdrawn. C23 citations closed. See CHANGELOG v5.11.2.\n")
A3 = "- **C23 run #1: 6 of 8 claims already known (2026-09-12, v5.11.1):**"
assert l.count(A3) == 1
L.write_text(l.replace(A3, bullet + A3), encoding="utf-8")
print(f"CLAUDE.md: bumped + Sec.2 one-liner ({len(bullet.split())} words)")
