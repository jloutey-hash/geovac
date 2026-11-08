"""claim_test_matrix rows for the end-to-end resource pricing."""
from pathlib import Path

T = "tests/test_paper60_preconditioner.py"
ROWS = [
    ("| 60 | eq:amplitude_floor — any `X` with `X^T S X = I` satisfies `||X|| = ||S^-1/2||` EXACTLY "
     "(since `X = S^-1/2 U`), so the block-encoding subnormalization floor is factorization-"
     "INVARIANT and the untreated route already attains it; preconditioning buys depth and only "
     f"depth | `{T}``::test_amplitude_floor_is_factorization_invariant` | tracked "
     "`geovac/sturmian_sigma_law.py` | **NEW 2026-09-12** | BACKED-SOUND. Tested against THREE "
     "genuinely different whitenings — symmetric `A^-1/2`, the preconditioned `P^-1/2 G^-1/2`, and "
     "an inverse-Cholesky factor (triangular, so not a disguised copy of the first) — all agreeing "
     "to 1e-8 relative, and each separately verified to BE a whitening. Fire-tested by scaling one "
     "factorization off the constraint: FIRES. This is the claim that caps the lever, so it is "
     "written to resist the hopeful reading rather than to confirm it |"),

    ("| 60 | sec:resource (pricing) — in exponents: untreated `alpha ~ n`, `d_inv ~ n^2` (product "
     "`n^3`); preconditioned with `G` COMPOSED, `d_inv` flat but `alpha` inherits `||P^-1/2||^2 ~ "
     "n^2` (product `n^2`); with a DIRECT encoding of `G` at the floor it would be `n`. So the "
     "lever is worth ONE power of n as priced, TWO if a direct block-encoding of `G` is found — "
     f"the open item, sized by `||G|| = 0.372` against a composed `alpha = 3196` at n=160 | `{T}`"
     "``::test_the_lever_buys_depth_and_costs_amplitude`` | same | **NEW 2026-09-12** | "
     "BACKED-SOUND. Asserts all four exponents from one ladder (n=20..160), including that the "
     "composed amplitude exponent is ~2 i.e. STRICTLY WORSE than the untreated ~1 — so the "
     "'preconditioning is a pure win' reading fails here, and the opposite overclaim fails via the "
     "flat-depth assertion. Fire-tested by pretending composition is free: FIRES. "
     "rests on: eq:amplitude_floor |"),
]

P = Path("docs/claim_test_matrix.md")
s = P.read_text(encoding="utf-8")
A = "| 60 | sec:resource (third lever, TRANSFER) —"
i = s.index(A)
P.write_text(s[:i] + "\n".join(ROWS) + "\n" + s[i:], encoding="utf-8")
print(f"claim_test_matrix: +{len(ROWS)} pricing rows")
