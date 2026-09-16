"""CHANGELOG entry + CLAUDE.md version cursor and Sec. 2 one-liner.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

ENTRY = """## [v5.11.18] - 2026-09-14

**Paper 12's 7.6% H2 residual is the sigma-only restriction, not the electron-electron cusp -- and restoring the azimuthal channels in its own basis reaches 99.09% of D_e.** Follow-up to the 2026-09-13 accuracy-path scan, which surfaced Tao-McCurdy-Rescigno as the one published chemical-accuracy result in GeoVac's own Level-2 coordinate system. Canonical memo `debug/sprint_tmr_method_memo.md`; new module `geovac/prolate_general_m.py`; backing `tests/test_paper12_azimuthal_channels.py`; fire tests `debug/firetest_p12_azimuthal.py`, `debug/firetest_p12_retracted.py`.

### The defect

Paper 12's basis carries no azimuthal dependence, so it spans only `m1 = m2 = 0`, and its kernel is projected onto the `m = 0` Neumann component. The paper justified that as "`1Sigma_g+` states have m = 0 symmetry". A `1Sigma_g+` state constrains the **total** `M = m1 + m2`, not each `m_i`: every `pi^2` and `delta^2` configuration has `M = 0` and is `1Sigma_g+`, and those carry the **angular** part of the correlation. The two restrictions are mutually consistent, which is why the calculation converged cleanly to a wrong answer and no test caught it for a decade.

### The measurement

| (j,l) | N | sigma only | N | \\|m\\| <= 1 |
|:--|--:|--:|--:|--:|
| (2,2) | 27 | 92.25% | 54 | 98.96% |
| (3,2) | 46 | 92.37% | 92 | 99.00% |
| (3,3) | 72 | 92.42% | 144 | **99.09%** |

`+11.64 mHa` at the largest basis. The sigma column reproduces Paper 12's published values to 161 / 1.6 / 1.0 / 6.7 / 58 uHa.

**Not a basis-count effect** -- the control is Paper 12's own convergence table: `N = 27 -> 46 -> 72` buys **0.34 mHa**, while opening the azimuthal axis at fixed `(j,l)` buys **11.6 mHa**. Asserted as the discriminating test.

### Three independent validations

- mu = 0 V_ee reproduces `geovac.neumann_vee` (the recurrence-based, quadrature-free path) **elementwise to 1.1e-9**.
- An independent Cartesian-Gaussian full CI -- different functions, different integrals, different code -- gives a sigma-only ceiling of **92.34%** (within **0.2 mHa** of Paper 12's value in an unrelated basis) and **99.10%** at `|m| <= 1`, against 99.09% here. Controls: H atom -0.499888, H2 at R=20 -0.999265, RHF -1.133270 vs the -1.13363 limit.
- The general-m kernel reproduces `1/|r1-r2|` **pointwise to 2e-6**.

### Two further defects in Paper 12, found by making m != 0 load-bearing

- **Eq. `neumann_full` as printed was wrong**: missing `(-1)^m` and `(2l+1)`, factorial ratio unsquared. It diverges pointwise. Invisible for a decade because every calculation used only the `m = 0` specialisation, where all three discrepancies vanish. Corrected; the test asserts the old form FAILS, so it discriminates.
- **Eq. `r12_prolate` omitted `cos(phi1 - phi2)`**, holding only at `phi1 = phi2`. Corrected.

### Two numerical traps, both recorded in the module docstring

- The Neumann sum terminates exactly only if the termination is **imposed**: the eta moment is identically zero for `l > Q + 2s - m`, but left to floating point the ~1e10 Legendre-derivative coefficients leave a residue the radial integral amplifies -- **E = -3.2e8 Ha** at `l_neumann = 18`. Now enforced twice (selection rule + l cap). The fire test showed the two are *redundant*: removing either alone does not break exactness, removing both does; the guard's docstring was corrected to claim only the conjunction.
- **Paper 12's headline basis is linearly dependent**: `cond(S) = 2.6e14` at N = 72, `2.0e16` at N = 144. A direct `eigh(H,S)` returned **-79 Ha**. Canonical orthogonalisation throughout. Paper 12's published `-1.161304` sits **58 uHa below** the conditioned value, so the last two of its six decimals are linear dependence, not physics. No conclusion turns on them.

### Dependents swept (Sec. 9 retraction -> dependents)

- **Paper 12**: abstract, intro, Sec. `neumann`, Sec. `gap` (retitled, cusp reading withdrawn), the hierarchy list (its fourth level -- "hyperspherical for the coalescence cusp" -- removed, since it was presented as a consequence of this residual), conclusion, plus a new Sec. "Restoring the Azimuthal Channels".
- **Paper 15**: its comparison of Level-4 `sigma+pi` (94.1%) against Paper 12 `sigma`-only (92.4%) was read as a cusp-resolution advantage for hyperspherical coordinates. Not like-for-like, and **the matched comparison reverses it**: 99.09% prolate vs 96.0% hyperspherical at l_max=6 with a cusp correction, and 99.97% for TMR. Withdrawn at 6 loci incl. abstract and conclusion.
- **Paper 13**: its motivation cited Paper 12's diagnosis. Restated on independent grounds (helium's coalescence is genuinely three-body).
- **group2 synthesis** (2 loci), `docs/validation_benchmarks.md`, `docs/claim_test_matrix.md` (3 new rows), the Paper-15 figure README, and `tests/test_level4_multichannel.py` (docstrings only -- **no assertion changed**, the 92.4 < pct < 100 band remains valid).
- New `check_retracted_terms.py` entry `p12-cusp-as-the-h2-gap` with 6 declared `cited_by`, **fire-tested on three withdrawn formulations**.
- 5 numeric-registry entries; the registry's own `test_every_convention_string_parses` caught their convention strings, so a new `accuracy` kind was added.

### Scope, stated plainly

`|m| = 2` was **not** obtained (the `d^4 Q_l` cancellations near `xi = 1` defeat this quadrature), so the delta contribution (~0.5 mHa) is known only from the Gaussian route. And the `mu > 0` radial integrals use a graded-panel spectral quadrature, **not** the `A_l/B_l/X_l` recurrences -- so the quadrature-free property does **not** yet extend to `mu > 0`. Generalising those three tables to associated Legendre functions is the open piece, now with a validated numerical reference to check against.

### Gates

C10 (compile, group2, 13 papers) / C14 / C21 / latex-escapes / retracted-terms (72/72 entries) all PASS. Regression: topo baseline + consumers + bearing-on set, 82 passed / 4 skipped. `tests/_durations.json` is still the known 2-byte stub, so the random tail-risk sample was skipped per the skill's documented fallback rather than silently widening scope.

**PI note:** shipped as a patch per the Sec. 9 default, but this is a retraction that moves published claims across four documents and adds a QA registry entry -- arguably a minor. PI call.

"""

CL = "CHANGELOG.md"
with io.open(CL, encoding="utf-8") as fh:
    text = fh.read()
anchor = "## [v5.11.17] - 2026-09-13"
if anchor not in text:
    print("FAILED: CHANGELOG anchor missing")
    sys.exit(1)
text = text.replace(anchor, ENTRY + anchor, 1)
with io.open(CL, "w", encoding="utf-8") as fh:
    fh.write(text)
print("  + CHANGELOG v5.11.18")

# ---------------------------------------------------------------- CLAUDE.md
CM = "CLAUDE.md"
with io.open(CM, encoding="utf-8") as fh:
    cm = fh.read()

if "**Version:** v5.11.17 (September 13, 2026)" not in cm:
    print("FAILED: version cursor missing")
    sys.exit(1)
cm = cm.replace("**Version:** v5.11.17 (September 13, 2026)",
                "**Version:** v5.11.18 (September 14, 2026)", 1)

bullet_anchor = "- **/qa group2 baseline FULL run COMPLETE"
if bullet_anchor not in cm:
    print("FAILED: Sec. 2 anchor missing")
    sys.exit(1)
new_bullet = (
    "- **Paper 12's H2 gap is the sigma-only restriction, not the cusp "
    "(2026-09-14, v5.11.18):** azimuthal channels restored in its own basis "
    "give 99.09% of D_e (+11.64 mHa vs 0.34 from tripling the sigma basis); "
    "2 more printed-equation defects; P13/P15/synthesis swept. "
    "See CHANGELOG v5.11.18.\n"
)
cm = cm.replace(bullet_anchor, new_bullet + bullet_anchor, 1)

with io.open(CM, "w", encoding="utf-8") as fh:
    fh.write(cm)
print("  + CLAUDE.md version cursor + Sec. 2 one-liner")
print("done")
