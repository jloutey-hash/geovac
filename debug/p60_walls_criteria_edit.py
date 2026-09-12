"""Owed items 3 and 4 (2026-09-12, PI-directed):
  * walls register -- split the composition wall into its three independent axes
  * docs/qa/criteria.md -- C23, the inverse-citation criterion

Scope note: this is the targeted consolidation the PI authorized, NOT a full
`/walls` re-verification sweep of the cluster (that remains PI-invoked).
"""
from pathlib import Path

# --------------------------------------------------------------- walls register
W = Path("docs/walls/register.md")
s = W.read_text(encoding="utf-8")

DELTA = """
### Delta 2026-09-12 (PI-directed consolidation, v5.10.19)

**The composition wall was one row; it is three independent axes, and they now have different statuses.** Driver: the Paper 60 KMS/preconditioner arc (CHANGELOG v5.10.18-19), three literature scans, `tests/test_paper60_{kms_attribution,preconditioner}.py`.

The register carried `||[P_A,P_B]|| = 0.50` as a measurement. It is the *saturation value* of an exact formula, `max_k sigma_k sqrt(1 - sigma_k^2)`, over the **same** singular spectrum that carries `cond(S) = (1 + sigma_max)/(1 - sigma_max)`. So the commutator and the conditioning are one object, as Paper 60 already said. What was riding along with them, and should not have been, is the `l`-block-structure loss. Splitting the three:

| Axis | Status | Mechanism | Changed? |
|:--|:--|:--|:--|
| **Conditioning** (`cond(S) ~ n^2`) | **BREACHED** | The symbol's zero at `chi = pi` has known order and location, so a band-Toeplitz preconditioner in the sense of Serra (*Math. Comp.* **66**, 651 (1997)) removes it: `P = tri(1,2,1)` exactly, DST-I diagonalizable in closed form, `cond -> 2.23` **flat in n** against 19 127 at `n = 160`. Whitening-invariant, so the spectrum is unchanged. | **YES** -- this axis was thought hard and is not |
| **Locality** (`S^{-1/2}` dense) | **STANDING** (capped, not removable) | The *other* pole. Preconditioning cures `chi = pi` and cannot touch the `chi -> 0` chirp, which fixes the off-diagonal envelope at `j^{-5/4}`. Measured: `G^{-1/2}` bandwidth 11->17 at 1e-2 (vs a fixed 0.72`n` for `S^{-1/2}`) but 27->141 at 1e-3. Profile exponent n-STABLE for `G^{-1/2}` (-1.19), DRIFTING for `S^{-1/2}` (-0.90 -> -0.78). | scoped |
| **`l`-block structure** | **STANDING, HARD** | **Proposition D**: if `S` is not `l`-block diagonal, no block-diagonal congruence orthogonalizes it -- Loewdin, canonical or Cholesky. Holds at **every** `cond(S) > 1` and does not relax as `cond(S) -> 1+`. Not a functional of the sigma spectrum at all. | **STRENGTHENED** |

**Consequence for the cluster's dispatch rule.** Rule (B) ("molecular chemical accuracy via a better basis/integral -> STOP unless the proposal breaches Wall B") is unchanged in outcome but its *reason* is now sharper: a proposal that improves conditioning no longer counts as progress toward sparsity, because Proposition D makes those independent. Conversely a conditioning-only proposal should no longer be rejected on Wall-B grounds -- that axis is open, and the quantum-resource lane (rule C) is where it pays.

**Falsifier for the split.** A congruence that is simultaneously (i) `l`-block diagonal and (ii) orthogonalizing, on a metric with nonzero inter-center coupling -- which Proposition D forbids outright; or a locality repair that reaches the `chi -> 0` chirp, which would have to change the *symbol*, not the matrix.

**Honest scope.** The breach is measured on the homonuclear two-center `s`-sector symbol, where the parity blocks are `I +- C`. Whether the construction reaches a polyatomic block with symmetry-inequivalent centers (water's `A_1`, `cond ~ N^1.97`) is **untested** -- and that is the case the gerade lever already fails, so it is the one that matters. Next probe.

---
"""

anchor = "---\n\n## Promotion candidates (PI-gated)"
assert s.count(anchor) == 1
s = s.replace(anchor, DELTA + "\n## Promotion candidates (PI-gated)")
W.write_text(s, encoding="utf-8")
print("walls register: 3-axis split appended")

# ------------------------------------------------------------------- criteria
C = Path("docs/qa/criteria.md")
s = C.read_text(encoding="utf-8")

C23 = """
## C23 -- inverse citation: is an UNcited result already known? (added 2026-09-12, PI direction)

**The gap this closes, measured.** On 2026-09-11 a literature scan found that Paper 60's `eq:sigma_law` -- derived, tested, tiered `[SYMBOLIC + MEASURED]`, and presented as the paper's own -- is the Kac-Murdock-Szego extreme-eigenvalue asymptotic (1953), exponent **and** constant. The paper carried 26 bibitems and zero Toeplitz-family references. **Three DELTA runs and a FULL run walked past it**, and none of them was at fault: every citation criterion in this document runs in one direction only.

- **C11** and the `citation-reviewer` dimension ask: *does the source we cited say what we claim?* They can only examine citations that exist.
- Nothing asks the inverse: *is this result, which we cite nobody for, already a named theorem?*

An uncited claim is invisible to a citation gate by construction. That is a criteria gap, not an execution miss, and it is the exact blind spot in which a rediscovery survives certification.

**The criterion.** For each target, enumerate the claims the paper presents as **its own derivations** -- the `[SYMBOLIC]` and `[SYMBOLIC + MEASURED]` tier, plus anything phrased as "we derive", "is in fact derived", "the law is", "we show that" -- and for each, ask whether it is a named result elsewhere. Priority order, because the criterion is not free:

1. **Closed-form asymptotics and scaling laws** with a clean constant (`pi^2`, `1/24`, an exponent that comes out an integer or a small rational). A named constant is the signature of a named theorem.
2. **Results in a well-developed external field** the corpus has entered sideways -- Toeplitz/finite-section, frame theory, harmonic analysis, spectral geometry. GeoVac reaches these through physics, so its authors are not embedded in their literatures.
3. **Anything whose derivation took under a page.** Short derivations of clean results are the ones most likely to have been done before.

**Verdict vocabulary.** `PRIOR ART` (found -- the paper must cite, and may claim only what remains: usually the identification, which is often the real contribution) / `ABSENT` (searched, not found -- record the search so the claim's novelty is defensible) / `UNVERIFIABLE` (could not reach the primary source -- say so in the record, and do not let an inaccessible source become an implied clearance).

**Hard rules.**

- **Re-verify the identification yourself.** When a scan reports prior art, reproduce the load-bearing leg independently before editing the paper. For the KMS case both legs were re-derived here (`c_1 = pi^2` from the exact tridiagonal spectrum; `b(1)` as a symbolic series coefficient) rather than taken on report.
- **Never add an unverified citation while fixing a citation defect.** Two references were dropped from that remediation for exactly this reason, and one fabricated reference was caught mid-edit and replaced with the verified one. A citation found by a search summary is not verified.
- **Finding prior art is not a retraction.** It re-tiers attribution, not truth. The measured content stands; what changes is who is credited, and the corpus's own contribution usually survives in sharper form. Record it as a re-attribution and sweep the dependents (this is the owner-corrected/citer-stale shape of Sec. 9).

**Scope.** Runs on FULL certification runs, not on every DELTA -- the cost is a literature scan per keystone claim. A DELTA run inherits the last FULL run's C23 verdicts unless the claim itself changed.

"""

anchor2 = "## C22 -- test-claim backing integrity (added 2026-08-31, PI direction)"
assert s.count(anchor2) == 1
i = s.index(anchor2)
nxt = s.find("\n## ", i + 10)
if nxt == -1:
    s = s + "\n" + C23
else:
    s = s[:nxt] + "\n" + C23 + s[nxt:]
C.write_text(s, encoding="utf-8")
print("criteria.md: C23 added after the C22 block")
