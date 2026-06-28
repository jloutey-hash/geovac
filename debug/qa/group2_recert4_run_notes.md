# /qa group2 — re-cert run #4 (2026-06-28): FAIL → remediated (thin residual layer)

Fourth `/qa group2` invocation (confirmation re-cert of the v4.50.2 corpus, HEAD 56e0d03).
Runs #1–3 FAIL→remediated; run #4 peeled a **thin** residual layer (2 loci) on a
**perfectly calibrated** panel.

## Verdict: FAIL (trustworthy) — but the closest yet

Calibration **PERFECT: sensitivity 11/11, specificity 6/6.** Every seed caught (5 code
incl. **cs3/P17 via the fast-test carry-forward** — the run-#3 P17 over-run is fixed; and
**cs5/FCI-M rs8-class caught a 3rd time**), every control clean. The panel is fully
calibrated per-agent across all five dimensions this run.

### Per-dimension
| Dimension | Exercised | Calibrated | Genuine material | Verdict |
|:--|:--|:--|:--|:--|
| Code C1/C2 | 9/9 (no over-run) | 5/5 | none (NITs + documented NO-TEST) | clean |
| Claims C3/5/6/8 | yes | 3/3 | 1 (fci_molecules C6) | FAIL |
| Citations C4 | yes | 2/2 | none (NITs only) | clean |
| Synthesis C9 | yes | 1/1 | 1 (cusp-variational lag) | FAIL |
| Deterministic C10–16 | yes | n/a | clean | PASS |

## Genuine material findings (verified, non-seed) — 2, both thin framing-syncs
1. **Synthesis cusp 0.004% labeled "variational"** at 3 loci (abstract:73-74, body:337-339,
   conclusion:776) + Table-I He row (152, NIT) — contradicts Paper 13's corrected
   "non-variational" (run-#2 fix). Synthesis-lag my run-#3 sweep scoped to paper_13 but
   missed in the synthesis. [LARGE]
2. **fci_molecules "graph Laplacian eigenvalues ε_n=−Z²/(2n²)"** (138-139, 169-170) — C6
   discrete-vs-continuum misstatement (the graph Laplacian is positive-semidefinite; the
   negative hydrogenic energies are its continuum image under the Fock projection, not
   graph-Laplacian eigenvalues). Long-standing; fci_atoms already says "hydrogenic
   eigenvalues" correctly. [SMALL]

Code-dimension non-seed "MATERIAL" flags reconciled to fix-on-sight NITs / documented
NO-TEST coverage gaps (P12 +0.01 loose guard on a claim true by 0.098 Ha; P17 structural
cross-check; P19/FCI-M n_max=3 magnitudes — §3-established / n_max=2-backed). None
cert-blocking.

## Remediation — DONE + VERIFIED (released v4.50.3)
Full-class sweep (avoiding patch-and-rerun leaks): all 3 synthesis cusp-variational prose
loci + the Table-I He row → labeled non-variational cusp extrapolation; both fci_molecules
C6 loci → "hydrogenic eigenvalues … continuum image under the Fock projection, not
graph-Laplacian eigenvalues". Verified: re-grep shows the class clean (remaining
"0.004+variational" hits are the corrected non-variational labels); deterministic
C11/C14/C16 PASS; synthesis + fci_molecules compile content_err=0, undef=0. Worktree torn
down; leak-scan clean.

## Honest ceiling / next
- The other "graph Laplacian eigenvalue" mentions corpus-wide (Paper 7, P18, P11:90, archive)
  are CORRECT usage (the L=D−A spectrum genuinely has eigenvalues) — NOT touched.
- Pattern across 4 runs: restorations (#1) → framing/§3-reassertion+P8-theorem (#2) →
  cross-paper second-locus propagation (#3) → 2 thin residual framing-syncs (#4). The layers
  are thinning monotonically; **run #5 is the certified-PASS attempt** and has a strong shot.
