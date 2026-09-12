"""C23 insertion, retry: the C22 heading uses an em-dash, not '--'."""
from pathlib import Path
import re

C = Path("docs/qa/criteria.md")
s = C.read_text(encoding="utf-8")
assert "## C23" not in s, "C23 already present"

C23 = """
## C23 — inverse citation: is an UNcited result already known? (added 2026-09-12, PI direction)

**The gap this closes, measured.** On 2026-09-11 a literature scan found that Paper 60's `eq:sigma_law` — derived, tested, tiered `[SYMBOLIC + MEASURED]`, and presented as the paper's own — is the Kac–Murdock–Szegő extreme-eigenvalue asymptotic (1953), exponent **and** constant. The paper carried 26 bibitems and zero Toeplitz-family references. **Three DELTA runs and a FULL run walked past it**, and none of them was at fault: every citation criterion in this document runs in one direction only.

- **C11** and the `citation-reviewer` dimension ask: *does the source we cited say what we claim?* They can only examine citations that exist.
- Nothing asks the inverse: *is this result, which we cite nobody for, already a named theorem?*

An uncited claim is invisible to a citation gate by construction. That is a criteria gap, not an execution miss, and it is the exact blind spot in which a rediscovery survives certification.

**The criterion.** For each target, enumerate the claims the paper presents as **its own derivations** — the `[SYMBOLIC]` and `[SYMBOLIC + MEASURED]` tier, plus anything phrased as “we derive”, “is in fact derived”, “the law is”, “we show that” — and for each, ask whether it is a named result elsewhere. Priority order, because the criterion is not free:

1. **Closed-form asymptotics and scaling laws** with a clean constant (`pi^2`, `1/24`, an exponent that comes out an integer or a small rational). A named constant is the signature of a named theorem.
2. **Results in a well-developed external field** the corpus has entered sideways — Toeplitz/finite-section, frame theory, harmonic analysis, spectral geometry. GeoVac reaches these through physics, so its authors are not embedded in their literatures.
3. **Anything whose derivation took under a page.** Short derivations of clean results are the ones most likely to have been done before.

**Verdict vocabulary.** `PRIOR ART` (found — the paper must cite, and may claim only what remains: usually the identification, which is often the real contribution) / `ABSENT` (searched, not found — record the search so the claim's novelty is defensible) / `UNVERIFIABLE` (could not reach the primary source — say so in the record, and do not let an inaccessible source become an implied clearance).

**Hard rules.**

- **Re-verify the identification yourself.** When a scan reports prior art, reproduce the load-bearing leg independently before editing the paper. For the KMS case both legs were re-derived here (`c_1 = pi^2` from the exact tridiagonal spectrum; `b(1)` as a symbolic series coefficient) rather than taken on report.
- **Never add an unverified citation while fixing a citation defect.** Two references were dropped from that remediation for exactly this reason, and one fabricated reference was caught mid-edit and replaced with the verified one. A citation found by a search summary is not verified.
- **Finding prior art is not a retraction.** It re-tiers attribution, not truth. The measured content stands; what changes is who is credited, and the corpus's own contribution usually survives in sharper form. Record it as a re-attribution and sweep the dependents (this is the owner-corrected/citer-stale shape of Sec. 9).

**Scope.** Runs on FULL certification runs, not on every DELTA — the cost is a literature scan per keystone claim. A DELTA run inherits the last FULL run's C23 verdicts unless the claim itself changed.

"""

m = re.search(r"^## C22 .*?test-claim backing integrity.*?$", s, re.M)
assert m, "C22 heading not found"
nxt = s.find("\n## ", m.end())
s = (s + "\n" + C23) if nxt == -1 else (s[:nxt] + "\n" + C23 + s[nxt:])
C.write_text(s, encoding="utf-8")
print("criteria.md: C23 inserted after the C22 block")
