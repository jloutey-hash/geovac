"""Register the azimuthal-channel numbers, and add the claim_test_matrix rows.

Sec. 15: a value appearing at more than one locus is a linked object, not free
text.  99.09% now appears in Paper 12, Paper 15, the group2 synthesis and
docs/validation_benchmarks.md.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

REG = "debug/qa/numeric_registry.py"
MATRIX = "docs/claim_test_matrix.md"

ANCHOR = '''    "lih_composed_pauli": dict('''

NEW_ENTRIES = '''    # ---- Paper 12 azimuthal channels (2026-09-14) ------------------------
    # The sigma-only restriction, not the cusp, is Paper 12's 7.6% residual.
    # Memo: debug/sprint_tmr_method_memo.md.  Backing:
    # tests/test_paper12_azimuthal_channels.py.
    "p12_sigma_only_de_pct": dict(
        value=92.42, convention="% of D_e, H2 R=1.4011, (j,l)=(3,3) sigma "
                                "only, canonical orthogonalisation",
        provenance="MEASURED 2026-09-14; Paper 12's own published 92.45 "
                   "(E=-1.161304) is 58 uHa lower and reflects cond(S)=2.6e14 "
                   "linear dependence rather than physics",
        aliases={92.45: "Paper 12's unconditioned N=72 value",
                 92.25: "(j,l)=(2,2), N=27"}),
    "p12_azimuthal_de_pct": dict(
        value=99.09, convention="% of D_e, H2 R=1.4011, (j,l)=(3,3), |m|<=1",
        provenance="MEASURED 2026-09-14 in Paper 12's own basis; reproduced "
                   "to 99.10 by an independent Cartesian-Gaussian full CI",
        aliases={98.96: "(j,l)=(2,2)", 99.00: "(j,l)=(3,2)",
                 99.10: "independent Gaussian-basis route"}),
    "p12_azimuthal_gain_mha": dict(
        value=11.64, convention="mHa gained by opening |m|<=1 at (3,3)",
        provenance="MEASURED 2026-09-14; the sigma-axis control over the same "
                   "enlargement is 0.34 mHa"),
    "p12_sigma_growth_mha": dict(
        value=0.34, convention="mHa gained by N=27 -> 72 along the sigma axis",
        provenance="MEASURED; read from Paper 12's own Table tab:convergence "
                   "(92.2 -> 92.4 -> 92.4 % of D_e)"),
    "p12_cond_s_33": dict(
        value=2.6e14, convention="cond(S), (j,l)=(3,3) sigma only, N=72",
        provenance="MEASURED 2026-09-14; rises to 2.0e16 at |m|<=1, N=144"),

'''

MATRIX_ROWS = '''| 12 | **azimuthal channels close the gap**: sigma-only basis spans only m1=m2=0; restoring |m|<=1 in the same basis with the same algebraic V_ee gives 99.09% of D_e (+11.64 mHa) vs 0.34 mHa from tripling the sigma basis | `test_paper12_azimuthal_channels.py`: `test_azimuthal_channels_close_the_gap`, `test_gain_is_a_channel_effect_not_a_count_effect` (the discriminator), `test_sigma_only_reproduces_paper12` | **BACKED-SOUND** | new 2026-09-14. Fire-tested (`debug/firetest_p12_azimuthal.py`): killing the m!=0 coupling FIRES. rests on: the general-m Neumann kernel, itself checked pointwise against 1/r12 |
| 12 | **general-m Neumann kernel** Eq. `neumann_full` (the printed form was missing (-1)^m and (2l+1) and unsquared the factorial ratio) | `test_paper12_azimuthal_channels.py`: `test_general_m_kernel_reproduces_coulomb` (against the Cartesian distance, and asserts the OLD form fails) | **BACKED-SOUND** | new 2026-09-14. The m=0 tests could never see this: all three discrepancies vanish at m=0 |
| 12 | **exact termination + conditioning caveat**: the Neumann sum terminates by selection rule (l > Q+2s-m), and the basis needs canonical orthogonalisation (cond(S)=2.6e14 at N=72) | `test_paper12_azimuthal_channels.py`: `test_neumann_truncation_is_exact_not_merely_converged`, `test_basis_is_linearly_dependent_at_the_largest_truncation` (@slow) | **BACKED-SOUND** | new 2026-09-14. Fire test showed the two exactness mechanisms are redundant; the guard discriminates the conjunction and its docstring says so |
'''

# ---------------------------------------------------------------- registry
with io.open(REG, encoding="utf-8") as fh:
    reg = fh.read()
if ANCHOR not in reg:
    print("FAILED: registry anchor not found")
    sys.exit(1)
reg = reg.replace(ANCHOR, NEW_ENTRIES + ANCHOR, 1)
with io.open(REG, "w", encoding="utf-8") as fh:
    fh.write(reg)
print("  + numeric_registry: 5 azimuthal entries")

# ------------------------------------------------------------ claim matrix
with io.open(MATRIX, encoding="utf-8") as fh:
    mat = fh.read()
marker = "| 12 | full V_ee \"exact\""
idx = mat.find(marker)
if idx < 0:
    print("FAILED: claim_test_matrix P12 row not found")
    sys.exit(1)
line_end = mat.index("\n", idx) + 1
mat = mat[:line_end] + MATRIX_ROWS + mat[line_end:]
with io.open(MATRIX, "w", encoding="utf-8") as fh:
    fh.write(mat)
print("  + claim_test_matrix: 3 rows for the azimuthal claims")
print("done")
