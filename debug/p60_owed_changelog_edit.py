"""Matrix rows + CHANGELOG v5.10.19 + CLAUDE.md for the owed-items pass."""
from pathlib import Path

T = "tests/test_paper60_preconditioner.py"
ROWS = [
    ("| 60 | sec:resource (third lever) — the ill-conditioning is a symbol ZERO of known order "
     "and location, so a band-Toeplitz preconditioner removes it (Serra, Math. Comp. 66, 651 "
     "(1997)): `g = 2+2cos(chi)` gives `P` EXACTLY tridiagonal(1,2,1) (its Hankel part vanishes "
     "identically since `a+b>=2` and `g_j=0` for `j>=2`), DST-I diagonalizable in closed form; "
     f"`cond(G)` -> 2.23 FLAT in n against 19127 at n=160 | `{T}`"
     "``::test_matching_polynomial_gives_an_exactly_tridiagonal_preconditioner`` + "
     "``::test_preconditioner_is_exactly_dst_diagonalized`` + "
     "``::test_preconditioned_conditioning_is_bounded_while_raw_grows`` | tracked "
     "`geovac/sturmian_sigma_law.py` (`sw_cross_block`) | **NEW 2026-09-12** | BACKED-SOUND. "
     "The conditioning test excludes a CONSTANT-FACTOR reading two ways a constant gain cannot "
     "satisfy: raw must quadruple per doubling (the n^2 law) while `cond(G)`'s increments must "
     "halve. Fire-tested: removing the preconditioner FIRES; flipping the zero to `chi=0` "
     "(`2-2cos`, the natural error — the standard discrete Laplacian) FIRES. Cross-checked "
     "against an independent Gauss-Legendre panel route to every printed digit. "
     "rests on: eq:sigma_law (the zero's order and location ARE the KMS input) |"),

    ("| 60 | sec:resource (third lever, legitimacy) — the lever does not change the problem: any "
     "`X` with `X^T S X = I` preserves the generalized spectrum, and `X = P^-1/2 G^-1/2` is such "
     f"an `X`, so the QSVT degree is set by `cond(G)` not `cond(S)` | `{T}`"
     "``::test_whitening_preserves_the_generalized_spectrum`` | same | **NEW 2026-09-12** | "
     "BACKED-SOUND. Asserts BOTH `X^T(I-C)X = I` and equality of the generalized eigenvalues "
     "against `scipy.linalg.eigh(H, A)`. Fire-tested with the tempting shortcut `X = P^-1/2` "
     "alone (the factor carrying the DST): FIRES on both halves |"),

    ("| 60 | sec:resource (third lever, LIMIT) — preconditioning cures the `chi=pi` pole exactly "
     "and CANNOT touch the `chi->0` chirp, so locality is capped: `G^-1/2` bandwidth 11->17 at "
     "1e-2 (vs a fixed 0.72n for `S^-1/2`) but 27->141 at 1e-3; the clean discriminator is the "
     "profile exponent, n-STABLE at -1.19 for `G^-1/2` and DRIFTING -0.90->-0.78 for `S^-1/2` | "
     f"`{T}``::test_locality_gain_is_real_at_one_percent_and_erodes_at_one_permille` + "
     "``::test_profile_exponent_is_n_stable_for_G_and_drifts_for_raw`` | same | "
     "**NEW 2026-09-12** | BACKED-SOUND, and written specifically to block the OVERCLAIM. The "
     "source scan reported `G` has 'n-independent decay'; that is right about the exponent and "
     "wrong about the bandwidth, so the test asserts the 1e-3 bandwidth GROWS by >2x — which an "
     "n-independent bandwidth cannot do — as well as the 1e-2 gain. The exponent test rejects "
     "BOTH collapses (all-stable = preconditioning changed nothing; all-drifting = it fixed "
     "nothing). Fire-tested both. rests on: the j^-5/4 chirp law (test_paper60_kms_attribution) |"),
]

P = Path("docs/claim_test_matrix.md")
s = P.read_text(encoding="utf-8")
A = "| 60 | §molecular [SYMBOLIC] — Proposition D:"
i = s.index(A); end = s.index("\n", i) + 1
P.write_text(s[:end] + "\n".join(ROWS) + "\n" + s[end:], encoding="utf-8")
print(f"claim_test_matrix: +{len(ROWS)} rows")

# ------------------------------------------------------------------ CHANGELOG
ENTRY = """## [v5.10.19] - 2026-09-12

**The owed items, and the conditioning half of the composition wall turns out to be breachable.** Follow-up to v5.10.18. Probes `debug/p60_{preconditioner,locality}_probe.py`; backing `tests/test_paper60_preconditioner.py`.

### The third conditioning lever, verified independently

The Toeplitz scan reported that Serra's theorem (*Math. Comp.* **66**, 651 (1997)) applies here. It does, and the construction is cleaner than expected. Because the ill-conditioning is a symbol **zero of known order and location**, the matching trigonometric polynomial is `g = 2 + 2cos(chi)` — quadratic zero at `chi = pi`, exactly where `1 - sigma`'s is — and in this basis its Hankel part vanishes identically (`a+b >= 2` while `g_j = 0` for `j >= 2`), so the preconditioner is **exactly** `tri(1,2,1)`, which the DST-I diagonalizes in closed form (verified against closed-form eigenpairs to `9e-15`).

Measured on the tracked `geovac.sturmian_sigma_law.sw_cross_block`, cross-checked against an independent Gauss-Legendre panel route to every printed digit:

| `n` | `cond(I+C)` | `cond(I-C)` | `cond(G)` |
|--:|--:|--:|--:|
| 10 | 2.383 | 81.9 | 2.096 |
| 40 | 2.537 | 1225.7 | 2.219 |
| 160 | 2.554 | 19126.5 | **2.229** |

The lever is legitimate rather than a change of problem: any `X` with `X^T S X = I` preserves the generalized spectrum, and `X = P^-1/2 G^-1/2` is such an `X`. On Paper 60's own resource model that is `d_inv ~ 16`, **flat in basis size**, against `3e5` untreated. It escapes the Bernstein `Theta(kappa)` floor the paper quotes rather than contradicting it — that floor constrains polynomial approximation of `x^-1/2` on `[1/kappa, 1]`, and preconditioning changes the *operator*, so `x^-1/2` is never approximated on the bad interval. **Paper 60's "mitigated, not dissolved" is therefore too pessimistic on the conditioning axis**, and now says so.

### And it stops exactly at the other pole

Preconditioning cures `chi = pi` and cannot touch the `chi -> 0` chirp. The first locality measurement here was **wrong and was discarded**: it fitted `exp(-d/L)` to a profile the symbol analysis already says is algebraic, so the returned `L` tracked the fit window and appeared to grow like `0.14n` for everything. Measured operationally instead —

- bandwidth for fixed relative accuracy at `1e-2`: `S^-1/2` needs a fixed **fraction** (`b/n = 0.72`, flat over `n = 32..256`); `G^-1/2` needs `b = 11 -> 17`.
- at `1e-3` the advantage **erodes**: `G^-1/2` needs `b = 27 -> 141`.
- clean discriminator, the profile exponent: **n-stable** at `-1.19` for `G^-1/2` (the chirp's own `-5/4`), **drifting** `-0.90 -> -0.78` for `S^-1/2` as the `chi = pi` singularity sharpens.

So the scan's "G has n-independent decay" is right about the exponent and wrong about the bandwidth; the backing test asserts the erosion explicitly so the paper cannot drift into the overclaim.

### `west_ruedenberg2013` removed

Cited alongside Amos–Hall (1961) and King *et al.* (1967) for "the cosines of the principal angles". Its full text is unreachable (HTTP 403 from two directions); its abstract describes a web of localizing orbital transformations and fast localization methods for quasi-atomic and split-localized orbitals — **no principal angles, no SVD, no corresponding orbitals**. It cannot be verified to support the attribution, so it is dropped from the citation and the bibliography. Amos–Hall and King are the verified lineage and remain. Not a claim that the paper is wrong; a claim that we could not check it, which is sufficient reason not to lean on it.

### A near-miss worth recording

Drafting the Serra bibitem, the applier carried a **fabricated** reference (a *BIT* 1994 entry assembled from memory) while the scan's actual source was *Math. Comp.* 66, 651 (1997). The correction was routed through a bash heredoc, which halves backslashes (`memory/feedback_no_heredoc_backslashes.md`); the guard assertion failed, the correction did not apply, **and the applier then ran with the unverified reference anyway**. Caught, verified against AMS (vol 66, no 218, April 1997, pp 651–665 — and the abstract matches the measurement: circulant preconditioners fail when the generating function has zeros, band-Toeplitz from trigonometric polynomials gives a condition number bounded independent of `n`), and replaced via a script file. Two lessons, both already written down: the heredoc rule is not optional, and *a citation found by a search summary is not a verified citation*. Now C23's second hard rule.

### Walls register: the composition wall is three axes, not one

The register carried `||[P_A,P_B]|| = 0.50` as a measurement; it is the saturation value of an exact formula over the same singular spectrum as `cond(S)`. Splitting what was one row: **conditioning = BREACHED** (above), **locality = STANDING but capped** (the chirp, `j^-5/4`), **`l`-block structure = STANDING, HARD and strengthened** (Proposition D — no block-diagonal congruence orthogonalizes a non-block-diagonal metric, at every `cond > 1`). Dispatch consequence: a conditioning-only proposal should no longer be rejected on Wall-B grounds, and conversely a conditioning gain is no longer evidence of progress toward sparsity. Scope: the breach is measured on the homonuclear two-center `s`-sector; whether it reaches water's `A_1` block — the case the gerade lever already fails — is **untested and is the next probe**.

### C23: the inverse-citation criterion

`/qa` could not have caught the KMS rediscovery, and no reviewer was at fault: C11 and the citation dimension ask whether a *cited* source says what we claim, and an uncited claim is invisible to them by construction. **C23** asks the inverse over the `[SYMBOLIC]`-tier claims a paper presents as its own, prioritised by the signatures that predict prior art (a clean closed-form constant; a well-developed external field entered sideways; a derivation under a page). Verdicts `PRIOR ART` / `ABSENT` / `UNVERIFIABLE`, with the hard rules that the identification must be re-verified locally before editing, that no unverified citation may be added while fixing a citation defect, and that prior art re-tiers attribution rather than retracting truth. Runs on FULL runs only. **Adding a QA criterion is a gate change and therefore a candidate minor (v5.11.0) — PI call; bumped as a patch by default per Sec. 9.**

### Gates

C10 / C21 / C16 / C22 / C14 / latex-escapes / internal-titles / inline-arxiv PASS in scope `paper_60`. New tests 6 passed + 1 slow-skipped in 2.1 s; all six guards fire-tested against the wrong answer each names (including "remove the preconditioner" and "use the `P^-1/2` shortcut").

"""

P = Path("CHANGELOG.md")
s = P.read_text(encoding="utf-8")
A = "## [v5.10.18] - 2026-09-11"
assert s.count(A) == 1
P.write_text(s.replace(A, ENTRY + A), encoding="utf-8")
print("CHANGELOG: v5.10.19 inserted")

C = Path("CLAUDE.md")
s = C.read_text(encoding="utf-8")
old = "**Version:** v5.10.18 (September 11, 2026)"
assert s.count(old) == 1
s = s.replace(old, "**Version:** v5.10.19 (September 12, 2026)")
bullet = ("- **Composition wall is 3 axes; conditioning BREACHED (2026-09-12, v5.10.19):** "
          "band-Toeplitz preconditioner, cond -> 2.23 flat. Locality capped by the chirp. "
          "New C23 gate. See CHANGELOG v5.10.19.\n")
A2 = "- **eq:sigma_law is Kac-Murdock-Szego (2026-09-11, v5.10.18):**"
assert s.count(A2) == 1
C.write_text(s.replace(A2, bullet + A2), encoding="utf-8")
print(f"CLAUDE.md: bumped + Sec.2 one-liner ({len(bullet.split())} words)")
