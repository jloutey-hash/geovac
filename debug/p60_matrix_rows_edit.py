"""Add claim_test_matrix rows for the Paper 60 KMS-attribution pass (2026-09-11)."""
from pathlib import Path

P = Path("docs/claim_test_matrix.md")
src = P.read_text(encoding="utf-8")

T = "tests/test_paper60_kms_attribution.py"
ROWS = [
    ("| 60 | §molecular [PRIOR ART] — eq:sigma_law IS the Kac-Murdock-Szego extreme-eigenvalue "
     "asymptotic, not an independent derivation: c_1 = pi^2 (KMS 1953, via Boettcher-Widom "
     "arXiv:math/0412269) x curvature b(1) = (kR)^2/24 reproduces the printed constant. What is "
     f"ours is the IDENTIFICATION of the SW metric as such a finite section | `{T}`"
     "``::test_kms_constant_c1_is_pi_squared`` + ``::test_symbol_curvature_b1_is_kR_squared_over_24`` "
     "+ ``::test_kms_product_reproduces_sigma_law`` | self-contained (closed-form tridiagonal "
     "eigenvalues + sympy series) | **NEW 2026-09-11** | BACKED-SOUND. c_1 is checked in CLOSED FORM "
     "(the standard (2,-1) tridiagonal, eigenvalues 4sin^2(k pi/2(n+1))), so the constant does not "
     "rest on the literature scan that found it. b(1) is an exact symbolic coefficient. Fire-tested: "
     "planting c_1 -> pi^2/2 and b(1) -> (kR)^2/6 each FIRE. rests on: eq:sigma_law |"),

    ("| 60 | §molecular [MEASURED] — the Boettcher-Widom smoothness hypothesis FAILS on our symbol "
     "and the constant holds anyway: the chi->0 chirp gives |c_j| ~ j^(-5/4), so sum_j j|c_j| "
     f"diverges | `{T}``::test_chirp_envelope_exponent_is_five_fourths` (3 kR values) + "
     "``::test_bottcher_widom_smoothness_hypothesis_fails`` | self-contained (deterministic "
     "Gauss-Legendre on phase-resolved panels) | **NEW 2026-09-11** | BACKED-SOUND. The envelope is "
     "taken as a running max over half-decade windows, so the predicted sin(2sqrt(2kRj)+pi/4) "
     "modulation cannot fake the exponent. Tolerance 0.08 excludes 1, 3/2 and 2. Divergence is "
     "asserted as SUSTAINED partial-sum growth (>1.35 across three doublings), which no convergent "
     "series can meet. Fire-tested: exponent -> -3/2 FIRES; replacing the chirp with a smooth symbol "
     "FIRES the divergence test. Independently reproduced by the 2026-09-11 Toeplitz literature scan "
     "(|b_k| ~ k^-1.25) |"),

    ("| 60 | §molecular [PRIOR ART] — physical restatement: 1 - sigma_max = (1/6)(R/L_max)^2 with "
     "L_max = 2n/(pi k) the longest wavelength the truncated basis carries, so the degeneracy switches "
     f"on exactly when L_max exceeds the bond length | `{T}``::test_Lmax_restatement_matches_sigma_law` "
     "| self-contained (exact symbolic identity) | **NEW 2026-09-11** | BACKED-SOUND. Exact symbolic "
     "difference against eq:sigma_law, so a wrong prefactor or a wrong L_max convention cannot pass. "
     "Fire-tested: 1/6 -> 1/24 FIRES. rests on: eq:sigma_law |"),

    ("| 60 | §molecular [SYMBOLIC] — Proposition D: a block-diagonal congruence cannot orthogonalize "
     "a metric that is not itself block diagonal, so within-m l-selection is lost at EVERY cond(S)>1 "
     "and the loss does not relax as cond(S)->1+. The l-sparsity cost is therefore INDEPENDENT of "
     f"eq:sigma_law, not a functional of the sigma spectrum | `{T}`"
     "``::test_prop_d_block_diagonal_congruence_preserves_block_structure`` | self-contained | "
     "**NEW 2026-09-11** | BACKED-SOUND. Tests BOTH directions: S block diagonal => S^-1/2 block "
     "diagonal (the m-selection half, which survives), and the contrapositive sampled over random "
     "block-diagonal X. Critically the contrapositive is run at eps=1e-6 (cond(S)=1+2e-6, essentially "
     "perfect conditioning) as well as at eps=0.3 — a guard testing only an ill-conditioned S would "
     "accept the wrong reading that this is a conditioning effect. Fire-tested: swapping the "
     "block-diagonal X for the eigenvector matrix FIRES. rests on: the composition-wall reading in "
     "Paper 58 / walls register (which this SHARPENS — the wall is stronger than stated there) |"),
]

ANCHOR = ("| 60 | sec:resource (validation) — the one-electron molecular isoenergetic")
i = src.index(ANCHOR)
end = src.index("\n", i) + 1
src = src[:end] + "\n".join(ROWS) + "\n" + src[end:]

P.write_text(src, encoding="utf-8")
print(f"inserted {len(ROWS)} rows after the sec:resource (validation) row")
