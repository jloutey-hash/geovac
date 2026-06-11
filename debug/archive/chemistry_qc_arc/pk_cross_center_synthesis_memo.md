# Phillips-Kleinman Cross-Center Sprint — Synthesis Memo

**Date:** 2026-05-08
**Sprint:** PK cross-center, named engineering closure for the W1b-residual
orthogonality wall (Phase C-W1c memo §6, Track 2 chemistry re-test
synthesis §"Recommended next sprint")
**Worker fork directive:** Implement Phillips-Kleinman-class projection on
the H valence orbital against the [Ne] core in cross-center geometry, and
test whether NaH binding recovers under balanced coupled FCI.

---

## 1. Verdict — short form

**Phillips-Kleinman cross-center barrier reduces the W1c-residual NaH
overattraction by 14.6% (0.357 Ha → 0.305 Ha) at standard physical
parameters but does NOT close the wall: NaH PES is still monotonically
descending with no internal equilibrium minimum.**

The W1b-residual orthogonality wall named in Phase C-W1c memo §6 is
**NOT the dominant residual** at NaH. PK orthogonality of the H valence
basis against the [Ne]-core orbitals on Na is necessary but not
sufficient for second-row chemistry binding. The residual ~290 mHa
overattraction (vs experimental D_e ≈ 75 mHa) sits in mechanisms beyond
W1c + W1b: candidate sources include (a) inner-core penetration
physics that the screened multipole expansion handles correctly in form
but may overestimate in magnitude near the unscreened nucleus,
(b) the hydrogenic Z_eff=1 basis on the Na valence side being
qualitatively wrong (it is not the actual Na 3s orbital), (c) the
cross-block ERI treatment between Na valence and H valence at
mismatched effective charges.

This is an HONEST NEGATIVE RESULT for the named engineering closure.
The wall is structurally deeper than W1b orthogonality.

## 2. What was implemented

New module `geovac/phillips_kleinman_cross_center.py` (~280 lines, 21
tests) provides the standard PK operator-form approximation of strict
Schmidt orthogonalization:

$$
\Delta H_{pq}^{\rm PK} = \sum_{c \in {\rm core}} (E_v - E_c)\,S_{pc}\,S_{cq}
$$

where $p, q$ index valence orbitals on the off-center side (e.g. the H
block in NaH), $c$ indexes the frozen-core orbitals on the heavy atom
(1s, 2s, 2p$_{\pm 1, 0}$ of [Ne] for $Z = 11{-}18$), $E_c$ is the
Clementi-Roetti HF orbital eigenvalue, $E_v$ is the valence reference
energy (default 0 → "absolute PK", purely repulsive since
$E_c < 0$), and $S_{pc}$ is the cross-center overlap.

**Cross-center overlaps**: computed via 2D Gauss-Legendre × radial-grid
numerical quadrature in the valence-orbital-centered frame:

$$
S_{pc}(R_{AB}) = \int_0^\infty r^2\,R_p(r)\,\Bigl[
\int_{-1}^{1} R_c(\rho(r,u))\,P_{l_p}^{|m|}(u)\,P_{l_c}^{|m|}(u_B(r,u))\,du
\Bigr]\,dr \cdot N_{\rm angular}
$$

with $\rho = \sqrt{r^2 + R^2 - 2rRu}$, $u_B = (ru - R)/\rho$, real-SH
azimuthal selection $m_p^{\rm signed} = m_c^{\rm signed}$ enforced
exactly. Direction-based rotation reuses
`shibuya_wulfman._build_block_rotation_matrix` for non-collinear nuclei.

**Hydrogenic radial wavefunctions at arbitrary points** computed
via the Bethe-Salpeter analytical normalization
$N_{nl} = \sqrt{(2Z/n)^3 (n-l-1)!/[2n(n+l)!]}$, consistent with
`composed_qubit._radial_wf_grid` to numerical precision.

**Core orbital eigenvalues**: hardcoded Clementi-Roetti HF values for
$Z = 11{-}18$ ([Ne] core). For $Z = 11$ (Na): $E_{1s} = -40.479$ Ha,
$E_{2s} = -2.797$ Ha, $E_{2p} = -1.518$ Ha. Other cores ([Ar], [Kr],
[Xe]) currently return empty list (PK barrier = 0 → wall documented but
not addressed for $Z \geq 19$ in this sprint).

**Wiring**: `geovac/balanced_coupled.build_balanced_hamiltonian` gained
two new kwargs: `pk_cross_center: bool = False` and
`pk_E_valence_ref: float = 0.0`. Bit-exactly backward compatible when
`pk_cross_center = False` (verified: `h1_pk_cross` is exact zero matrix
when off, and N_pauli / 1-norm match Sprint 7 NaH baseline at R=3.5).

**Tests**: `tests/test_phillips_kleinman_cross_center.py` 21/21 passing.
Coverage includes:
- Real-SH normalization at $l \in \{0, 1\}$, $m \in \{0, \pm 1\}$
- Hydrogenic radial normalization (1s, 2p)
- Self-overlap at $R \to 0$ → 1
- Decay at large $R$
- Azimuthal selection (different $m$ → 0)
- Decay with increasing $Z_{\rm orb}$
- PK barrier symmetry, positive-semidefinite at $E_v = 0$
- First-row $Z_{\rm nuc}$ → barrier identically zero
- Direction-based rotation = z-axis at direction (0,0,1)

## 3. NaH PES results

Driver: `debug/pk_cross_center_nah.py`. Data: `debug/data/pk_cross_center_nah.json`.

R-grid: 7 points from 2.5 to 10.0 bohr. Configurations:
1. **screened+PK** (W1c + W1b PK)
2. **screened only** (W1c only — Track 2 baseline)
3. **PK only** (no screening, PK barrier on bare V_ne)
4. **bare** (Sprint 7 baseline, no screening, no PK)

| Config | $E_{\min}$ (Ha) | $R_{\min}$ (bohr) | Descent depth (Ha) | Internal min? |
|:-------|----------------:|------------------:|-------------------:|:-------------:|
| screened+PK    | -163.264 | 2.5 | 0.305 | **NO** |
| screened only  | -163.316 | 2.5 | 0.357 | NO |
| PK only        | -170.949 | 2.5 | 6.097 | NO |
| bare (Sprint 7)| -171.097 | 2.5 | 6.244 | NO |

**Reductions**:
- Bare → screened: 17.5× (W1c does the heavy lifting)
- Screened → screened+PK: **1.17×** (W1b PK adds 14.6% on top of W1c)
- Bare → screened+PK: 20.5×

**Comparison to physical NaH**: experimental $R_{eq} = 3.566$ bohr,
$D_e \approx 0.075$ Ha. Our screened+PK descent of 0.305 Ha is still
**4.1× deeper than physical $D_e$**. The PES does not have a turning
point in $R \in [2.5, 10]$ bohr, so $R_{eq}$ is not extracted.

PK trace per R-point (E_v=0):
| R (bohr) | PK trace (Ha) |
|---------:|--------------:|
| 2.5  | 0.106  |
| 3.0  | 0.081  |
| 3.5  | 0.065  |
| 4.0  | 0.052  |
| 5.0  | 0.031  |
| 7.0  | 8.5e-3 |
| 10.0 | 8.2e-4 |

Trace decays exponentially with $R$ as expected (overlaps decay with
core orbital extent).

## 4. Strong-PK probe (E_v scan)

Tested at the four key R points (2.5, 3.5, 5.0, 10.0 bohr):

| $E_v$ (Ha) | Descent depth (Ha) | Note |
|:----------:|-------------------:|:-----|
| 0          | 0.305 | Standard absolute PK |
| 100        | 0.148 | Hint at internal min: $E$(R=3.5) < $E$(R=2.5) |
| 1000       | 0.253 | Overcompensates at small R |
| 10000      | 0.321 | Strongly overcompensates |

The minimum-depth signature at $E_v = 100$ Ha is suggestive but
**unphysical as a parameter choice**: standard PK theory uses $E_v$ at
the actual valence-state eigenvalue (~$-0.4$ Ha for Na 3s, $-0.5$ Ha
for H 1s in molecule), giving $(E_v - E_c) \approx 40$ Ha for Na 1s,
nearly identical to my $E_v = 0$ default. The "sweet spot" near
$E_v = 100$ is a tuning artifact, not a derivation.

The strict-Schmidt limit is $E_v \to \infty$ with proper renormalization
of the orthogonalized basis. Naive scaling of $E_v$ in the PK operator
form (without renormalization) does not converge to strict Schmidt;
it overshoots and ruins the PES at large $E_v$. A proper Schmidt
implementation would require recomputing all valence-block matrix
elements in the orthogonalized basis and renormalizing — substantially
more work than the PK operator approximation.

## 5. Why doesn't PK close the wall? — diagnostic

The Phase C-W1c memo §6 listed three candidate sources for the
W1c-residual:
1. **Inner-core penetration of the H orbital** — a real H 1s tail
   extending into the unscreened Na nucleus region experiences strong
   bare-Coulomb attraction. The screened multipole expansion captures
   this in form (Z_eff(ρ) approaches Z=11 as ρ → 0) but the framework
   may overestimate the magnitude because the H orbital amplitude in
   the inner-core region is not properly suppressed by Pauli
   exclusion (which is what PK is supposed to enforce).
2. **Missing kinetic-energy repulsion (Pauli on H 1s vs [Ne])**
   — W1b territory. THIS sprint addresses this; result is 14.6%
   reduction.
3. **W1a cross-register two-body operator** — recoil / quantum
   nuclear motion. Already partially closed by Phase C-W1a; the
   leading-order recoil at NaH is too small (~10⁻⁴ Ha) to close 0.3 Ha.

The fact that strict-Schmidt-equivalent PK ($E_v \to$ "sweet spot
near 100 Ha") does NOT close the wall either suggests sources beyond
this list. Plausible additional candidates:

4. **Hydrogenic Z_eff=1 basis on Na is qualitatively wrong**. The Na
   side of the bond block uses $Z_{\rm orb} = 1$ hydrogenic orbitals.
   These are NOT the actual Na 3s/3p valence orbitals — they have
   wrong nodal structure, wrong scale, wrong polarization response.
   The cross-block ERI between H 1s (proper) and Na "valence" (wrong)
   carries this error. A proper basis would use either explicit
   atomic-orbital functions or radial Schrödinger eigenstates of the
   FrozenCore Z_eff(r) potential.
5. **Cross-block ERI mass scale**. The balanced builder's cross-block
   ERI infrastructure was developed for first-row LiH where both
   centers have small Z. At Na with frozen core, the cross-block
   integration may not have the same accuracy.
6. **Multipole truncation at $L_{\max} = 2 l_{\max} = 4$**. For
   $l_{\max} = 2$ orbitals, $L_{\max} = 4$ multipoles is exact for
   the cross-V_ne by Gaunt selection. But at very small $R$ and
   strongly compressed orbitals, the multipole expansion convergence
   should be re-examined (probably not the issue here, included for
   completeness).

Without a more refined basis on the Na valence side and/or a more
careful cross-V_ee treatment, the second-row chemistry-solver wall
remains structurally open after this sprint.

## 6. Comparison: same-center PK vs cross-center PK

LiH explicitly carries a Li-core block in the balanced architecture,
giving 4 electrons in 2 blocks (Li_core + LiH_bond) with proper FCI
between core and valence. The same-center PK (`ab_initio_pk.py`) is
applied within the Li_core block and produces a clean ~5% R_eq error
at $l_{\max} = 2$ (the Track CD baseline).

In NaH the Na [Ne] core is FROZEN — there is no Na_core block in the
FCI sector. The frozen-core treatment relies on the Hartree screening
$Z_{\rm eff}(r)$ to give the correct effective potential to the
valence electrons, and on the core electrons being well-separated from
valence in energy (Koopmans). The cross-center PK barrier is the
analog of LiH's same-center PK, but applied across the bond — the H
valence orbital orthogonalized against orbitals living on a different
nucleus.

The KEY structural difference is that LiH's same-center PK lives in a
finite n_max basis where all orbitals (core + valence) are sampled
simultaneously, so antisymmetric Slater determinants enforce Pauli
automatically. NaH's cross-center PK has only the valence orbitals in
the FCI sector; the [Ne] core is external to the FCI Hilbert space and
its "presence" is encoded only via Z_eff(r) screening (a one-body
local potential) and now via the PK barrier (a one-body non-local
projector). One-body operators alone cannot fully encode the Pauli
exclusion that an explicit core block provides.

This may be the structural limitation: the frozen-core approximation
trades off Hilbert-space size against Pauli enforcement, and what the
W1b-residual is telling us is that the trade-off is not benign for
short bond distances where valence and core orbitals overlap
substantially.

## 7. Sprint outcome and recommendation

**Sprint deliverable status:**
- ✅ Module: `geovac/phillips_kleinman_cross_center.py` (~280 lines, clean)
- ✅ Tests: `tests/test_phillips_kleinman_cross_center.py` (21/21 passing)
- ✅ Wiring: `balanced_coupled.py` extended with `pk_cross_center` kwarg, bit-exact backward compat
- ✅ NaH PES: `debug/pk_cross_center_nah.py` + JSON
- ❌ Binding recovery: NaH still unbound under PK + W1c (negative engineering result)
- ➖ MgH₂ extension: not run (per directive, only run if NaH binds)

**Sprint verdict**: PK cross-center is necessary and clean as an
architectural addition (the framework now properly enforces
Pauli-exclusion-class orthogonality between valence and frozen-core
orbitals via a standard operator approximation), but it is NOT
sufficient to close the second-row chemistry-solver wall.

**Recommended next sprint** (in priority order):
1. **Refined Na valence basis**: replace the $Z_{\rm orb} = 1$ hydrogenic
   orbitals on the Na center with radial Schrödinger eigenstates of
   the FrozenCore Z_eff(r) potential. The infrastructure already
   exists in `neon_core._solve_screened_radial`. This addresses
   diagnostic candidate #4 (the basis on Na is qualitatively wrong)
   directly.
2. **Strict Schmidt orthogonalization**: re-implement the PK
   correction by computing ALL matrix elements (h1, ERI) in the
   orthogonalized valence basis after explicit Schmidt-of-core
   subtraction. Substantially more work than PK barrier but
   mathematically exact in the strict-Schmidt sense.
3. **Cross-block ERI audit**: re-examine the cross-block V_ee infrastructure
   for the Na-side / H-side asymmetric-Z case. If the cross-block ERI is
   evaluating Z_orb-1 hydrogenic orbitals on Na (consistent with the
   valence basis) but the Na-side electron physically should see screened
   nuclear attraction during electron-electron interaction, the ERI
   evaluation may have a hidden inconsistency.

(1) is the cheapest and most likely to move the needle; (2) is
mathematically definitive but expensive; (3) is exploratory.

Each of these sprints would be a separate engineering effort. The
total program of "make second-row chemistry-solver bind" requires
sequential progress on these wall components.

## 8. CLAUDE.md update

A short paragraph for §2 active frontier and a row in §3 failed-approaches
table. The §3 row is honest as the Phillips-Kleinman cross-center
barrier did not close the wall as the named engineering closure for W1b.

## 9. Files changed

| Path | Change |
|:-----|:-------|
| `geovac/phillips_kleinman_cross_center.py` | New (~280 lines): cross-center overlap + PK barrier with real-SH azimuthal selection, direction rotation, Clementi-Roetti [Ne] core eigenvalues. |
| `geovac/balanced_coupled.py` | Added `pk_cross_center: bool = False` and `pk_E_valence_ref: float = 0.0` kwargs to `build_balanced_hamiltonian`. Bit-exact backward compat when off. PK details exported in result dict (`h1_pk_cross`, `pk_cross_center_count`, `pk_cross_center_details`). |
| `tests/test_phillips_kleinman_cross_center.py` | New: 21 unit tests (all passing). |
| `debug/pk_cross_center_nah.py` | New: 7-point NaH PES driver, 4 configs (bare / screened / PK only / both). |
| `debug/data/pk_cross_center_nah.json` | New: PES results in machine-readable form. |
| `debug/pk_cross_center_synthesis_memo.md` | This memo. |
| `papers/group2_quantum_chemistry/paper_17_composed_geometries.tex` | Architectural clarification subsection (Track 2 retrospective: drift lives in `composed_diatomic.py`, balanced coupled is a different solver). |
| `papers/group2_quantum_chemistry/paper_19_coupled_composition.tex` | New subsection on W1c-residual orthogonality wall + PK cross-center engineering closure (negative result). |
| `CLAUDE.md` | §2 sprint outcome paragraph (post-Track-2 PK cross-center sprint). §3 row for negative engineering result. |

---

**End of synthesis memo.**
