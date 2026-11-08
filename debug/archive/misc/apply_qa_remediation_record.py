"""Record the /qa paper_60 FULL run and its remediation. Idempotent."""
from __future__ import annotations

import sys

CL = "CHANGELOG.md"
DOD = "docs/qa/paper_60.done.md"
CL_ANCHOR = "## [v5.11.4] - 2026-09-12\n"
CL_MARKER = "## [v5.11.5]"

ENTRY = """## [v5.11.5] - 2026-09-12

**`/qa paper_60` FULL, PI-invoked = FAIL, remediated. NOT certified.** Six dimensions, unseeded, against criteria frozen with PI approval after extending the C8 headline list to cover v5.11.0-v5.11.4 (thirteen new headlines, 15-27; adding claims to be checked, no goalpost relaxed). Deterministic layer 14/14. **Zero mathematical defects in any dimension.**

### The asymmetry, again, and sharper than usual

Both code reviewers rebuilt the pipelines on independent routes. **Every published converged constant reproduced** -- the atomic exponents to 4-5 decimals including all seven falling local slopes, the two encoding exponents at three box sizes, the conditioning columns, the ratio-symbol sup, the water columns, the M-centre orders (by a symbolic Newton polygon, *stronger* than the test asserts), the chirp constant and phase by an independent mpmath route, the gerade constant to 12 digits. All 46 bibitems resolve at source. Twenty honest-scope sentences checked for direction; none inverted.

**Every defect was staleness, scoping, or a guard that does not discriminate** -- and one root cause: each re-tiering of the last two days reached the body paragraph owning the claim and missed the abstract, the conclusion, the scope paragraph, the module docstring and the synthesis. The Sec. 9 dependents rule cannot fire on this, because the stale citers are *inside the same file* as the corrected owner.

### The five that mattered

1. **A wrong number in the abstract.** The He chain endpoint read `-2.897` at K=164; the registry's own alias gives `-2.8964`. Confirmed three ways -- the alias, and the paper's own `eq:no_selection` interlacing applied to its own K=244 value, which forces `E(164) >= -2.896667` and which `-2.897` violates. A retired 60-bohr value, unregistered and therefore invisible to C21, sitting under a paragraph certifying the section converged. Fixed at both loci and **registered** (`p60_he_chain_spdf_k164`). **OWED:** the two earlier rungs are unverified and are candidates for the same defect.
2. **Abstract and conclusion contradicted the body** on the newest headline, still confining the molecular lever to the case v5.11.0 breached, and still calling the gerade sector "perfectly conditioned" against `tab:resource`'s own row. Rewritten; the closing verdict now reads removable-on-conditioning, capped-on-locality, structurally-closed-on-sparsity.
3. **The pre-registered load-bearing control was blind.** `blockdiag(T,T)` is `I2 (x) T` and the rotation is `V (x) I`: they COMMUTE, so rotated and unrotated frames give identical spectra (`||PQ-QP|| = 0`, verified). Planting the rotation into it did not fire. Replaced by the genuinely discriminating counterfactual -- the SELECTIVE preconditioner in the unrotated frame, which is *worse than untreated* (4.4e6 vs 4.2e4 at n=96) -- plus a commutation pin so the blind version cannot be restored silently.
4. **The tau prior-art surrender of v5.11.4 was too generous, and it was ours.** The paper glossed the DST-I algebra as "Toeplitz minus Hankel, the structure of this section". Form is not membership: tau needs the coefficient sequence to TERMINATE, and the chirp's does not. Measured: the DST-I leaves 4.1%/2.0%/1.1% off-diagonal on the cross block at n=16/64/160 against 2e-14 for a genuinely-tau matrix. The claim-matrix row written the same day records the gap flatly while the paper asserted the opposite. The corrected form is **stronger**: the identity covers the tau idealisation, and the departure is exactly the residue -- the grid-sampled prediction reproduces `pi^2/24` to 9e-6 while the true `sigma_max` departs by 0.26%.
5. **A guard that was literally `x == x`**, backing an abstract headline: it computed a 1-norm then asserted that norm equalled itself on an unmodified matrix, with the loop body discarding its own result. Replaced by building the per-state alternative the docstring names -- via `pk_ref`, a parameter **no test had ever varied** -- and asserting it moves the 1-norm while the pipeline's matrix does not.

### Two reviewer findings the PM overturned

Verification is between the reviewer and the record, and it earned its place twice. A prior-art attribution was reported as a work that could not be located after three targeted searches, with a recommendation to drop the count and re-price the paper's novelty concession from three to two; **the work exists** (Wulfman & Takahata, *J. Chem. Phys.* **47**, 488 (1967)) and its abstract names the Lie algebras of E4, R5 and O(4,1) exactly as attributed -- found in two queries. The concession stands at three counts; the defect is a missing bibitem, not a phantom. Separately, a wrong-basis-point finding was reported against the abstract as well as the conclusion; the abstract states it correctly.

### Instrument findings

- **Three gates reported PASS on scope `trunk` rather than `paper_60`** on the first invocation, because they take a different flag. Caught and re-run. This is the gate-self-audit class: a gate that scopes its verdict away is indistinguishable from a working one.
- **C16 reports clean on a file it scopes** while a retired reading lives in it, twice -- the rising-slope sequence inside `eq:sublinear`'s own backing test, and the withdrawn ill-conditioning mechanism in `sturmian_l2_encoding.py`'s docstring. Patterns match wording; the zombies are written in different words.
- **A claim-matrix row cited a test deleted on 2026-09-07** and was marked BACKED-SOUND.
- **The staleness banner measures one file of the two in scope** for single-paper certs: its glob takes the paper and not the synthesis, though the trunk and group branches both add theirs. C9 is GATING, and the group2 synthesis has five commits since the certified date.
- **C21 examines nothing in the molecular half** (141 unannotated decimal literals past the last `\\gvq`), and three registry keys added 2026-09-12 are cited from no `.tex` at all.

### Remediated

Six passes, content before guards per Sec. 9: the wrong number + registration; the abstract; the conclusion and Acknowledgments; eight body corrections; the two guards as their own reviewed pass; the stale-text sweep across the synthesis, two `geovac/` docstrings, a test comment and the claim matrix. **Four new guard assertions, all fire-tested** (restore the blind control; rotate the control's frame; break the commutation pin; make `pk_ref` inert). Deterministic layer 14/14 green in the correct scope; 141 + 64 tests pass.

**NOT certified.** The verdict stands at FAIL until a delta-verification run over this remediation comes back clean.

"""

DOD_ANCHOR = "## Change log\n"
DOD_MARKER = "2026-09-12 — **FULL run"
DOD_NEW = """## Change log
- 2026-09-12 — **FULL run (PI-invoked) = FAIL, remediated. NOT certified.** Six dimensions,
  unseeded. C8 extended to headlines 15–27 first (PI-confirmed) so v5.11.0–v5.11.4 were
  pre-registered rather than unmeasured. **Zero mathematical defects**; every published
  converged constant independently reproduced by both code reviewers; 46/46 bibitems resolve.
  Five load-bearing findings: a retired `-2.897` He chain endpoint in the abstract (violating
  the paper's own interlacing bound, now `-2.896` and registered); abstract + conclusion
  pre-breach against `sec:resource`; the water control provably blind (it commutes with the
  rotation it names); the τ prior-art surrender over-scoped (our matrices are not in that
  algebra — measured 1.1% off-diagonal at n=160); and an `x == x` guard behind an abstract
  headline. Two reviewer findings OVERTURNED by PM verification (a "missing" 1967 reference
  exists; a wrong-basis-point finding did not apply to the abstract). Instrument findings:
  three gates scoped to `trunk` on first invocation; C16 clean on files it scopes while
  zombies live in them; a claim-matrix row citing a deleted test; the staleness banner
  measuring the paper and not the synthesis on a GATING dimension; C21 blind to the molecular
  half. Remediated in six passes, content before guards; four new guard assertions all
  fire-tested. **Next: a delta-verification run over the remediation.**
"""


def main() -> int:
    n = 0
    with open(CL, encoding="utf-8") as fh:
        t = fh.read()
    if CL_MARKER not in t and t.count(CL_ANCHOR) == 1:
        with open(CL, "w", encoding="utf-8") as fh:
            fh.write(t.replace(CL_ANCHOR, ENTRY + CL_ANCHOR))
        n += 1
        print("  ok   CHANGELOG.md")
    with open(DOD, encoding="utf-8") as fh:
        d = fh.read()
    if DOD_MARKER not in d and d.count(DOD_ANCHOR) == 1:
        with open(DOD, "w", encoding="utf-8") as fh:
            fh.write(d.replace(DOD_ANCHOR, DOD_NEW))
        n += 1
        print("  ok   paper_60.done.md")
    print(f"applied {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
