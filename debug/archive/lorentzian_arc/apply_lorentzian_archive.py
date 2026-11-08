r"""Archive the supplanted Lorentzian tail: Papers 46, 47, 48, 49.

PI direction 2026-09-14. Scope is the tail identified by C24, not the whole
Lorentzian arc:

  ARCHIVED  46, 47, 48, 49 -- descoped/partial, and their dependents are mostly
            other descoped papers (C24 check B: unhealthy counts 4, 6, 2, 1
            against healthy 2, 2, 1, 1).
  KEPT      45 -- DESCOPED but load-bearing: five healthy dependents including
            Paper 38, the WH1 keystone. Partially rebuilt in the action-seminorm
            framework.
  KEPT      43, 50 -- ACTIVE, seven and three external citers.
  KEPT      52, 53 -- DRAFT with no dependents, which is UNFINISHED rather than
            supplanted. CLAUDE.md Sec. 9 says explicitly that zero dependents is
            not grounds; these need finishing or dropping, a different decision.

All four are classed CLOSED-VALID, not CLOSED. Each has surviving content that
the INDEX already records, and a future sprint should reuse it rather than
rebuild it:
  46  Lemma 3.2, the degeneracy diagnosis
  47  the norm-resolvent arrow and the three-carrier identification
  48  the bridge's categorical design
  49  the cocycle-deficit / TICI algebra

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

INDEX = "papers/INDEX.md"
REG = "docs/retired_papers.md"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# ---- INDEX: remove the four live rows ------------------------------------
edit(INDEX,
     """| 46 `paper_46_strong_form_lorentzian_propinquity.tex` | **DESCOPED** | Strong-form claims pending product repair; Lemma 3.2 = the degeneracy diagnosis |
| 47 `paper_47_two_rate_hybrid_convergence.tex` | PARTIAL | Norm-resolvent arrow + three-carrier identification stand; propinquity arrow descoped |
| 48 `paper_48_krein_ms_bridge.tex` | PARTIAL | Bridge design conditional on repair; T3/T6 descoped |
| 49 `paper_49_oslpls_strong_form_bridge.tex` | PARTIAL | Λ inheritance descoped; cocycle-deficit / TICI algebra survives |
| 50 `paper_50_cft3_partition_function.tex` | ACTIVE | Bit-exact F-theorem match on S³ and S⁵ (Klebanov–Pufu–Safdi) |""",
     """| 46--49 | **ARCHIVED 2026-09-14** | The strong-form Lorentzian tail. Moved to `archive/`; see the Archive table below and `docs/retired_papers.md`. Paper 45 stays: it is descoped but load-bearing (five healthy dependents, including Paper 38). |
| 50 `paper_50_cft3_partition_function.tex` | ACTIVE | Bit-exact F-theorem match on S³ and S⁵ (Klebanov–Pufu–Safdi) |""",
     "INDEX: group1 rows replaced by an archived pointer")

# ---- INDEX: add the tombstones -------------------------------------------
edit(INDEX,
     """| 21 `paper_21_geometric_vacuum_synthesis.tex` | Superseded by the two group syntheses |""",
     """| 21 `paper_21_geometric_vacuum_synthesis.tex` | Superseded by the two group syntheses |
| 46 `paper_46_strong_form_lorentzian_propinquity.tex` | **Archived 2026-09-14.** Strong-form construction descoped (degenerate seminorm). *Survives:* Lemma 3.2, the degeneracy diagnosis. |
| 47 `paper_47_two_rate_hybrid_convergence.tex` | **Archived 2026-09-14.** Propinquity arrow descoped. *Survives:* the norm-resolvent arrow and the three-carrier identification. |
| 48 `paper_48_krein_ms_bridge.tex` | **Archived 2026-09-14.** Metric-level theorems open pending a repair the register records as closed. *Survives:* the bridge's categorical design. |
| 49 `paper_49_oslpls_strong_form_bridge.tex` | **Archived 2026-09-14.** Λ inheritance descoped. *Survives:* the cocycle-deficit / TICI algebra. |""",
     "INDEX: four archive tombstones added")

# ---- register rows -------------------------------------------------------
edit(REG,
     """| `paper_21_geometric_vacuum_synthesis.tex` | SUCCESSOR-COVERED | pre-2026-05-22 | Superseded by the two group syntheses. | — |""",
     """| `paper_21_geometric_vacuum_synthesis.tex` | SUCCESSOR-COVERED | pre-2026-05-22 | Superseded by the two group syntheses. | — |
| `paper_46_strong_form_lorentzian_propinquity.tex` | CLOSED-VALID | 2026-09-14 | Strong-form construction descoped: the full-Krein operator-norm Lipschitz seminorm is degenerate, so the closed-form bound evaluates a rate formula rather than a distance. **Survives: Lemma 3.2, the degeneracy diagnosis** — reuse it, do not re-derive it. | strong-form Lorentzian propinquity; full-Krein operator-norm Lipschitz seminorm; strong-form Lipschitz seminorm |
| `paper_47_two_rate_hybrid_convergence.tex` | CLOSED-VALID | 2026-09-14 | Propinquity arrow descoped with the Paper 45 substrate. **Survives: the norm-resolvent arrow and the three-carrier identification.** | two-rate hybrid convergence; de-compactification gap; compact-to-non-compact resolvent comparison |
| `paper_48_krein_ms_bridge.tex` | CLOSED-VALID | 2026-09-14 | Metric-level theorems open pending a repair whose own verdict is that it weakens the claim to convention. **Survives: the bridge's categorical design** (the contravariant functor and the F2-mismatch resolution via thermal time). | Krein-MS bridge; pinned proper quantum metric space; covered Lorentzian pre-length space; F2 mismatch |
| `paper_49_oslpls_strong_form_bridge.tex` | CLOSED-VALID | 2026-09-14 | Λ inheritance descoped. **Survives: the cocycle-deficit / TICI algebra.** | OSLPLS; operator-system Lorentzian pre-length space; cocycle deficit; TICI algebra |""",
     "register: four CLOSED-VALID rows with surviving content named")

# ---- register: record the arc-level verdict once -------------------------
edit(REG,
     """## Adding a row""",
     """## Note: the Lorentzian tail (Papers 46--49, archived 2026-09-14)

These four were archived together as one decision. The model they rest on --
a Lorentzian quantum metric via K⁺ compression of the Krein Dirac -- was
withdrawn by the Paper 45 annihilation theorem, and the identified repair path
(Toeplitz temporal compressions) carries its own recorded verdict that it
weakens the claim to convention. The register lists them as CLOSED-VALID
rather than CLOSED because **each has surviving content that is not refuted**;
the trigger terms exist so that a future sprint reuses those parts instead of
rebuilding the whole arc to reach them.

**Papers 43, 45 and 50 were deliberately NOT archived.** Paper 45 is descoped
but load-bearing, with five healthy dependents including Paper 38, the WH1
keystone. Papers 43 and 50 are ACTIVE with seven and three external citers.
**Papers 52 and 53 were also not archived**: they are DRAFT with no dependents,
which is unfinished rather than supplanted, and §9 says explicitly that zero
dependents is not grounds.

---

## Adding a row""",
     "register: the arc-level decision and its boundary recorded")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

applied, failed = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    for old, new, label in items:
        if old in t:
            t = t.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append("%s   [%s]" % (label, path))
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)

print("applied %d of %d" % (len(applied), len(EDITS)))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
