"""CHANGELOG v5.10.18 + CLAUDE.md Sec.2 one-liner and version bump."""
from pathlib import Path

ENTRY = """## [v5.10.18] - 2026-09-11

**`eq:sigma_law` is Kac-Murdock-Szego (1953): the corpus had rediscovered a known asymptotic, exponent *and* constant, and claimed it.** Also: the `l`-selection loss is separated from conditioning and comes out *stronger*; two citation defects fixed. PI-directed conversational thread, not a `/qa` run; three parallel literature scans. Memos `debug/lit_scan/{toeplitz_finite_section,frames_riesz_overcompleteness,sturmian_conditioning_prior_art}_memo.md`, probe `debug/p60_symbol_pole_decay_note.md`.

### The attribution

Paper 60 derives `1 - sigma_max = (kR)^2 pi^2 / (24 n^2)` for the two-centre Shibuya-Wulfman metric and presents it as its own. It is the **Kac-Murdock-Szego extreme-eigenvalue asymptotic**: in the normal form `|1-t|^{2a} b(t)` used by Boettcher & Widom (arXiv:math/0412269), `lam_min ~ (c_a / n^{2a}) b(1)` with `c_1 = pi^2` due to Kac, Murdock & Szego (*J. Rational Mech. Anal.* **2**, 767 (1953)). Our symbol is the `a = 1` case with curvature `b(1) = (kR)^2/24`; the product reproduces the printed constant exactly. **The paper carried 26 bibitems and zero Toeplitz-family references.**

Both legs were re-verified here rather than taken on the scan's report: `c_1 = pi^2` from the standard `(2,-1)` tridiagonal, whose spectrum `4 sin^2(k pi / 2(n+1))` is exact (`lam_min (n+1)^2 = 9.869604` at `n = 10^4`), and `b(1)` as an exact sympy series coefficient. **What survives as ours is the identification**, and it is worth having: that the SW metric in the sine basis *is* such a finite section -- Toeplitz minus Hankel, `<n|a|m> = c_{n-m} - c_{n+m}` -- with symbol `j_0(kR cot(chi/2))`; equivalently that the SW operator is multiplication by the translation phase `e^{ip.R}` on the Fock sphere, whose angular average goes trivial at `p = 0`. The conditioning exponent 2 is then just the order of the symbol's maximum.

Added with it, because the mechanism is now legible: `1 - sigma_max = (1/6)(R/L_max)^2` with `L_max = 2n/(pi k)` the longest wavelength the truncated basis carries. **The degeneracy switches on exactly when the basis starts carrying wavelengths longer than the bond -- which a complete basis must eventually do.** Overcompleteness is the price of completeness, not a defect of the basis.

### An honest residue, and a second pole

Our symbol does **not** satisfy the smoothness hypothesis under which Boettcher-Widom prove the constant, and the constant holds anyway. At the opposite end `chi -> 0` the symbol is a chirp (amplitude `~chi`, phase `~2kR/chi`) whose Fourier coefficients decay only as `|c_j| ~ j^{-5/4}`, so `sum_j j|c_j|` diverges. That exponent is new here, parameter-free from stationary phase:

    |c_j| = (2pi)^{-1/2} 2^{-3/4} (kR)^{-1/4} j^{-5/4} |sin(2 sqrt(2 kR j) + pi/4)|

verified over `j = 64..65536` and `kR` in {1,2,5}, including the **sign pattern** -- a phase prediction, not a fit. Two adaptive-quadrature routes fail intermittently past `j ~ 512` (they disagree and return values that *grow* with `j`, impossible for a continuous symbol); a deterministic Gauss-Legendre route on phase-resolved panels is stable to `1e-11` and is the one to trust. The Toeplitz scan independently measured `k^{-1.25}` for the same object.

Consequence: **the banded-Loewdin lever is closed, negative.** The gerade symbol `(1+a)^{-1/2}` is bounded and smooth at the IR pole -- `cond(I+C) -> 2.555041`, flat -- but inherits the chirp linearly, coefficient ratio measured `-0.4988/-0.5018/-0.5002` at `j = 4096/16384/65536` against the predicted `-1/2`. `l1` band-truncation error therefore falls only as `b^{-1/4}`: bandwidth `~4e5` for `1e-2`, `~4e9` for `1e-3`. The ungerade symbol `(1-a)^{-1/2} ~ 1/(pi-chi)` is not even in `L^1`. **The gerade sector is perfectly conditioned and still not local** -- conditioning fails at the IR pole, locality at the UV pole, and no single operation reaches both. That is a mechanism for the cost-conservation pattern the walls register carries as an observation.

### Proposition D: the `l`-selection loss is not a conditioning effect

If `X` is invertible and block diagonal w.r.t. `H = ⊕_l H_l` and `X†SX` is block diagonal, then `S = X^{-†}(X†SX)X^{-1}` is block diagonal too. Contrapositively, if `S` is not `l`-block diagonal then **no** block-diagonal congruence -- Loewdin, canonical, Cholesky -- orthogonalizes it. The two-centre metric couples `l` while preserving `m`, so `m`-selection survives and within-`m` `l`-selection cannot, **at every `cond(S) > 1`**, and does not relax as `cond(S) -> 1+`.

So the sparsity cost is *independent* of `eq:sigma_law`, not a functional of the sigma spectrum. Paper 60's "the two walls are functionals of one object" is correct for `cond(S)` and `||[P_A,P_B]||`; the `l`-block loss was riding along with them and is a third thing. The wall is **stronger** than the paper stated, not weaker. (One correction to the scan that surfaced this: it phrased the claim as failing "at any condition number, even 1". At `cond = 1` exactly the coupling vanishes, `S = I`, and selection is perfect -- the accurate statement is *every* `cond > 1`, a discontinuity at zero coupling rather than a large-`kappa` effect.)

### Citation defects

- **Halmos (1969) does not state the commutator norm.** "Two subspaces" gives the canonical form; `||[P_A,P_B]|| = max_k sigma_k sqrt(1-sigma_k^2)` is a one-line consequence of it (in the `2x2` block at principal angle `theta_k` the commutator has norm `sin theta_k cos theta_k`). Reworded to cite the canonical form and derive the norm inline. The identity itself is correct.
- `bottcher_spitkovsky2010` rescoped with it; `west_ruedenberg2013` remains over-characterised and **unread from both directions** (HTTP 403 on the full text), logged as owed.
- The paper's "`~1%` at `n=160`" is the asymptotic's own `O(1/n)` term, not scatter: the relative residue halves under each doubling, `0.134 -> 0.010` across `n = 10..160` at `kR = 2`.

### What the scans settled, and what they did not

- **Prior art for the conditioning analysis: ABSENT**, and the repo's two standing claims survive. Aquilanti/Cavalli/Coletti/Calderini, the Avery canon, Shibuya-Wulfman and successors are about completeness, closed-form integral evaluation and *energy*-convergence -- never the metric's spectrum. Herbst-Avery-Dreuw (PRA **99**, 012512) was fetched and searched directly: zero hits. **Weakest link, flagged:** the two Avery books could not be read in full, so that leg is search-index absence rather than a verified read.
- **The overcompleteness wall is a theorem, and an elementary one.** If `g != 0` lies in the closed span of `{f_i}`, then `lam_min(G_N) <= dist(g, V_N)^2 -> 0`. Completeness of the one-centre set alone forces it; "translate" is incidental. Ron-Shen fiberization gives the operator form: Riesz sequence iff `ess inf (1 - |sigma|) > 0` iff `||sigma||_inf < 1`, and ours is exactly 1, attained.
- **Balian-Low does NOT transfer** and must not be cited: `ab = 1` is essential (at redundancy > 1 the obstruction disappears) and the mechanism is topological, needing a lattice we do not have. Beurling density / Ramanathan-Steger likewise. BCHL is technically available but buys a weaker conclusion at the cost of an `l1`-localization hypothesis.
- **Not claimed, left open:** whether a weaker hypothesis (Serra-Capizzano, *LAA* **270** (1998) is the likely home) covers our symbol class. The citation was found by search but not verified to primary-source standard, so it is *not* in the paper.
- **A lattice correction worth keeping.** `(n,l,m)` *is* a lattice -- the SO(4)/SU(2) weight lattice, with the S^3 harmonics as Peter-Weyl matrix elements. What is missing is a lattice in the *translation* direction (`{0,R}` is two points, not a subgroup; make it one and you have a crystal). And the fibration Ron-Shen wants already exists here with `n` dual to `chi`, so the right analogy is band theory with the Fock angle as quasi-momentum: `sigma(chi)` is the band function, `1 +/- sigma` the two branches, and the lower band touches **zero** at `chi = pi`. Localized Wannier/Loewdin functions need that band bounded off zero. The *topological* Wannier no-gos still do not apply -- no Bloch bundle.

### Backing and gates

New `tests/test_paper60_kms_attribution.py` (9 tests, 13 s), written as a **separate pass** from the edits it protects per Sec. 9 and fire-tested against the specific wrong answer each guard names: `c_1 -> pi^2/2`, `b(1) -> (kR)^2/6`, the KMS product halved, envelope `-5/4 -> -3/2`, `L_max` prefactor `1/6 -> 1/24`, the chirp replaced by a smooth symbol, and -- for Proposition D -- the block-diagonal `X` swapped for the eigenvector matrix. **All seven fire.** Proposition D's contrapositive is deliberately run at `eps = 1e-6` (`cond = 1 + 2e-6`) as well as at `0.3`, because a guard testing only an ill-conditioned `S` would accept exactly the wrong reading it exists to exclude.

**C22 caught the author.** The new matrix row said "closed-form ... eigenvalues", which trips the Paper 34 zombie guard (that paper's *negative* result about eigenvector closed forms). Unrelated claim, genuinely ambiguous phrase -- reworded to "exact tridiagonal spectrum" rather than weakening the guard.

Claim-impact sweep: four live dependents restated the law as derived (`certified_reference_values.md` `anchor.collapse_pi2_24`, `topic_to_paper_lookup.md`, the original `claim_test_matrix` row, and `docs/qa/paper_60.done.md`, which was ratifying the attribution -- the ".done.md as re-infection vector" class, given a supersession note). `development_frontier_archive.md` deliberately **not** edited: it is a verbatim historical record (Sec. 13.11 rule 10).

Gates: C10 / C21 / C16 / C14 / C22 / latex-escapes / headline-numbers / inline-attributions / inline-arxiv / internal-titles / duration-language all PASS in scope `paper_60`.

### Process note for the PI

**`/qa` could not have caught this.** C11 verifies that a *cited* source says what the paper claims; nothing asks the inverse -- *is this uncited derived result already a named theorem?* Three DELTA runs and a FULL run walked past `eq:sigma_law`. That looks like a genuine criteria gap rather than an execution miss, and whether it becomes a criterion is a PI call.

"""

P = Path("CHANGELOG.md")
s = P.read_text(encoding="utf-8")
anchor = "## [v5.10.17] - 2026-09-11"
assert s.count(anchor) == 1
s = s.replace(anchor, ENTRY + anchor)
P.write_text(s, encoding="utf-8")
print("CHANGELOG: v5.10.18 inserted")

# ---- CLAUDE.md: version bump + Sec.2 one-liner -----------------------------
C = Path("CLAUDE.md")
s = C.read_text(encoding="utf-8")

old_v = "**Version:** v5.10.17 (September 11, 2026)"
assert s.count(old_v) == 1, "version line not found"
s = s.replace(old_v, "**Version:** v5.10.18 (September 11, 2026)")

bullet = ("- **eq:sigma_law is Kac-Murdock-Szego (2026-09-11, v5.10.18):** a 1953 theorem, "
          "claimed as ours; the identification survives. l-selection loss separated from "
          "conditioning. See CHANGELOG v5.10.18.\n")
anchor2 = "- **/qa paper_60 DELTA = DEFECTS, remediated (2026-09-11, v5.10.17):**"
assert s.count(anchor2) == 1
s = s.replace(anchor2, bullet + anchor2)
C.write_text(s, encoding="utf-8")
print("CLAUDE.md: version bumped, Sec.2 one-liner added ({} words)".format(
    len(bullet.split())))
