# Paper Update Change Log — Task 12

Date: 2026-03-27
Source: Diagnostic arc report (Tasks 1–11), diagnostic_arc_report.md

## Paper 17: Composed Natural Geometries (paper_17_composed_geometries_updated.tex)

### Change 1: Expanded l_max divergence section (Sec V.A, Outlook)
**Original lines:** 697–730
**Action:** Replaced the existing 34-line divergence discussion with an ~95-line expanded version.

**What was changed:**
- Extended the l-dependent PK divergence table from l_max=2–4 to l_max=2–7, adding N_ch and ΔR_eq columns
- Added `\paragraph{Mechanism: differential angular correlation.}` explaining the root cause with quantitative evidence (w(0,0) = 0.9998 e-e only vs 0.243 nuclear only)
- Added `\paragraph{Eigenvalue convergence at fixed geometry.}` noting β ≈ 0.42 exponential convergence rate
- Added `\paragraph{R-dependent pseudopotential scaling.}` with Eq.~(\ref{eq:rdep_pk}), comparison table of l-dependent vs R-dependent vs self-consistent at l_max=4 and 7, and discussion of the asymptotic drift rate convergence
- Replaced the "path forward" paragraph with updated text reflecting what was tested: l-dependent (2.5x improvement), R-dependent (~2x additional at moderate l_max, converges asymptotically), self-consistent (50% drift reduction), and the clean resolution (full 4-electron Level 4)

**New equation:** `\label{eq:rdep_pk}` — R-dependent PK scaling formula

### Change 2: Updated conclusion (Sec VI)
**Original lines:** 1052–1056 (after \item 5)
**Action:** Added a new enumerated item (item 6) summarizing the l_max divergence findings.

**Content:** The l_max divergence mechanism (differential angular correlation, 100% nuclear origin), R-dependent PK result (13.5% → 2.0% at l_max=4), and the pointer to full 4-electron Level 4 as the clean resolution. References Eq.~(\ref{eq:rdep_pk}).

### Summary of additions to Paper 17
- ~60 net lines added
- 1 new equation (eq:rdep_pk)
- 1 new table (l-dep vs R-dep vs self-consistent comparison)
- 1 existing table expanded (l_max=2–4 → l_max=2–7)
- 1 existing conclusion item added
- No content deleted (existing text replaced/expanded in place)

---

## Paper 15: Level 4 Geometry (paper_15_level4_geometry_updated.tex)

### Change 1: Heteronuclear convergence paragraph (Sec VII, Discussion)
**Original line:** After line 850 (end of "Heteronuclear systems and the charge-center origin" subsection)
**Action:** Added `\paragraph{Heteronuclear convergence and the composed-geometry boundary.}` (2 paragraphs, ~15 lines)

**Content:** Exponential convergence of μ₀ at fixed (R, R_e) with β ≈ 0.42 and R² = 0.997. The nuclear-only origin of wavefunction spreading (w(0,0) = 0.9998 vs 0.243). Clarification that the R_eq divergence is a composed-geometry (PK) property, not a Level 4 limitation. References Paper 17 via \cite{paper17}.

### Change 2: Channel weight distribution paragraph (Sec V.D, Heteronuclear extension)
**Original line:** After line 541 (end of charge-center origin paragraph, before "Dissociation limit")
**Action:** Added `\paragraph{Channel weight distribution.}` (~15 lines)

**Content:** Quantitative channel weights for H₂ (77% in (0,0), 3 channels for 90%) vs 6:1 asymmetry (26% in (0,0), 8 channels for 90%). The dipolar channels carrying 36% of the norm. The nuclear-only origin of spreading. Practical implication: l_max ≥ 4 needed for asymmetric bonds vs l_max = 2 for homonuclear.

### Change 3: Added \bibitem{paper17} to bibliography
**Original line:** After line 1016 (after \bibitem{paper13})
**Action:** Added Paper 17 bibliography entry for the new \cite{paper17} cross-reference.

### Summary of additions to Paper 15
- ~35 net lines added
- 2 new paragraphs in Discussion
- 1 new paragraph in Heteronuclear extension
- 1 new bibliography entry
- No content deleted

---

## Companion files
- `papers/core/diagnostic_arc_report.md` — comprehensive technical report from Task 11
- `debug/data/task11_higher_lmax.json` — Task 11 Part A numerical results
- `debug/plots/task11_req_vs_lmax.png` — R_eq vs l_max plot (all PK modes)
- `debug/plots/task11_drift_increments.png` — drift increment plot
- `debug/plots/task11_pes_rdependent.png` — PES curves (R-dependent PK)
