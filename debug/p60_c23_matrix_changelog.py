"""C23 run #1: matrix row updates + CHANGELOG v5.11.1 + CLAUDE.md."""
from pathlib import Path

# ---- matrix: the chirp claim was upgraded from an exponent to a closed form ---
M = Path("docs/claim_test_matrix.md")
s = M.read_text(encoding="utf-8")

OLD = ("| 60 | §molecular [MEASURED] — the Boettcher-Widom smoothness hypothesis FAILS on our symbol "
       "and the constant holds anyway: the chi->0 chirp gives |c_j| ~ j^(-5/4), so sum_j j|c_j| "
       "diverges |")
NEW = ("| 60 | §molecular [SYMBOLIC + MEASURED] eq:chirp_decay — the Boettcher-Widom smoothness "
       "hypothesis FAILS on our symbol and the constant holds anyway; the chi->0 chirp gives the "
       "FULL closed form `|c_j| = (2pi)^-1/2 2^-3/4 (kR)^-1/4 j^-5/4 |sin(2 sqrt(2 kR j) + pi/4)|` "
       "(**re-attributed 2026-09-12, C23 run #1: this is a Bessel asymptotic, DLMF 10.32.10 at "
       "nu=2 plus 10.40.2 -- no stationary-phase argument needed, and the pi/4 is the (pi/2z)^1/2 "
       "BRANCH phase, not a stationary-phase signature**), so sum_j j|c_j| diverges while "
       "sum_j |c_j| CONVERGES (Wiener algebra, 5/4 > 1) |")
assert s.count(OLD) == 1
s = s.replace(OLD, NEW)

EXTRA = (
    "| 60 | §molecular eq:chirp_decay (constant + phase) — the closed form's leading constant "
    "`(2pi)^-1/2 2^-3/4` and its `pi/4` branch phase, not merely the `-5/4` envelope | "
    "`tests/test_paper60_kms_attribution.py``::test_chirp_closed_form_constant_and_phase` "
    "(3 kR values) | self-contained | **NEW 2026-09-12** | BACKED-SOUND. The PHASE is pinned "
    "TOLERANCE-FREE by sign agreement across every sampled j where the modulation is not near a "
    "zero -- so the guard tests the right reason, not just the right number (the pi/4 would have "
    "come out the same from the wrong mechanism). Fire-tested three ways: dropping the pi/4, "
    "dropping 2^-3/4 from the constant, and flipping the kR exponent -- all FIRE. "
    "rests on: DLMF 10.32.10 / 10.40.2 |\n"
    "| 60 | §molecular (flat limit) — the rank-`M-1` all-ones degeneracy at `chi=pi` is the FLAT "
    "LIMIT of the RBF/kernel literature (Barthelme-Usevich 2021), not new; what is ours is that "
    "the block symbols realise it at a SYMBOL POINT rather than a shape-parameter limit, which is "
    "why a fixed rank-`M-1` rotation removes it | `tests/test_paper60_preconditioner.py`"
    "``::test_lever_transfers_to_water_A1_block`` + ``::test_water_needs_the_null_direction_rotation`` "
    "| tracked `geovac/sturmian_sigma_law.py` | **re-attributed 2026-09-12 (C23 run #1)** | "
    "BACKED-SOUND; the tests are unchanged, only the attribution moved. See "
    "`docs/qa/c23_run_001_paper_60.md` |\n")
A = "| 60 | sec:resource (third lever) —"
i = s.index(A)
M.write_text(s[:i] + EXTRA + s[i:], encoding="utf-8")
print("claim_test_matrix: chirp row re-attributed, +2 rows")

# ---------------------------------------------------------------- CHANGELOG
ENTRY = """## [v5.11.1] - 2026-09-12

**The two owed items: the resource lever priced honestly (it shrinks), and C23's first run (six prior-art catches, two of them on claims written this week).** Probes `debug/p60_resource_pricing_probe.py`; scans `debug/lit_scan/c23_paper60_{linalg,analysis}_memo.md`; record `docs/qa/c23_run_001_paper_60.md`.

### Pricing: the lever buys depth, and nothing else

An invariance settles it. Any `X` with `X^T S X = I` satisfies `X = S^{-1/2} U`, so `||X|| = ||S^{-1/2}||` **exactly, independent of the factorization** — verified to `6e-13` across three genuinely different whitenings (symmetric, preconditioned, inverse-Cholesky). **No factorization can lower the block-encoding subnormalization**, and the untreated route already attains it. New `eq:amplitude_floor`.

At `n = 160`, `kR = 2`, chemical accuracy:

| route | `alpha` | `d_inv` | product |
|:--|--:|--:|--:|
| untreated | 125.4 | 3.1e5 | 3.9e7 |
| preconditioned, `G` composed | 3196 | 16.1 | 5.1e4 |
| preconditioned, `G` direct | 125.4 | 16.1 | 2.0e3 |

In exponents: untreated `alpha ~ n`, `d_inv ~ n^2`, product `n^3`; preconditioned with `G` obtained by *composing*, the degree goes flat but `alpha` inherits `||P^{-1/2}||^2 ~ n^2`, giving `n^2`; with a direct encoding of `G` at the floor it would be `n`. **So the lever is worth one power of `n` as priced, two if a direct block-encoding of `G` is found.** The v5.11.0 entry's `19000x` was depth-only; the honest end-to-end figure is `758x` at `n=160`, growing like `n`. The prize is now sized rather than named: `||G|| = 0.372` against a composed `alpha = 3196`, a factor `8.6e3` paid for nothing but the order of operations. `P^{-1/2}` itself costs no block-encoding calls — DST-I has an `O(log^2 N)` circuit (Klappenecker–Rötteler, verified).

### C23 run #1: six of eight audited claims were already known

The criterion was written because `eq:sigma_law` turned out to be Kac–Murdock–Szegő. On its first run it caught five more, **including two claims written this week**. Nothing was false; C23 re-tiers attribution, not truth, and every identity was re-verified numerically here before any edit.

- **A1** — the `spec{1±σ}` / `cond = (1+σ)/(1−σ)` / principal-angles chain is classical (Jordan; Jordan–Wielandt; the two-block CBS constant). Now used rather than derived.
- **A2** — Halmos's primary re-read: Theorem 2 is the canonical form **only**, no norm anywhere. Loring 2014 added for the identity; the inline derivation stands.
- **A3 — the catch.** "Proposition D", written 2026-09-11, is the Löwdin symmetry-preservation property specialized to the `l` grading — known in this paper's own field since Slater–Koster (1954), and in operator terms the statement that block-diagonal matrices are a commutant and therefore inverse-closed. **Demoted** from Proposition; the `l`-vs-`m` application is what the paper now claims.
- **A4** — the constant `2.555041…` is unnamed, but its mechanism is textbook (a self-adjoint Toeplitz spectrum is the convex hull of the symbol's essential range). Now named.
- **B1** — the `j^{-5/4}` law **needed no deriving**: the model integral is DLMF 10.32.10 at `nu = 2`, and 10.40.2 delivers constant *and* phase. The claim is therefore **upgraded**, from an exponent to `eq:chirp_decay` in full.
- **B4** — the rank-`M-1` all-ones degeneracy, written *yesterday*, is the **flat limit** of the RBF/kernel literature (Barthelmé–Usevich 2021). What survives as ours is that the block symbols realize it at a *symbol point* rather than a shape-parameter limit — which is exactly why a fixed rotation removes it.
- **B2** ABSENT (nobody composes `j0` with a cotangent). **B3** PRIOR ART *and our mechanism reading is correct*: with decay `5/4 > 1` the metric sits in Jaffard's class, where inversion preserves decay **given bounded invertibility** — so the escape is the spectrum touching zero, not the localization class.

### Two corrections to our own reasoning

1. **The `pi/4` is a branch phase**, from the `(pi/2z)^{1/2}` prefactor, *not* the stationary-phase `sign(phi'')·pi/4`. The same number for the wrong reason. The new guard pins it by **sign agreement**, tolerance-free, so the mechanism is tested rather than the value.
2. **`sum |c_j|` CONVERGES** — `5/4 > 1`, so the symbol is in the Wiener algebra. What diverges is `sum j|c_j|`, Böttcher–Widom's hypothesis, a different condition. The paper said the right thing; the distinction was not drawn, and it is load-bearing for B3.

### What was deliberately not done

The WebSearch budget (200) ran out partway. `WebFetch` survived, so anything with a known URL was reachable and anything needing a *search* was not. Cited here: **only** DLMF 10.32.10/10.40.2 (quoted verbatim with phase conditions), Loring 2014, Barthelmé–Usevich 2021. Eleven further primaries — Jordan 1875, Slater–Koster, Hartman–Wintner, Jaffard 1990, Gröchenig–Leinert TAMS 358, Driscoll–Fornberg, Björck–Golub among them — are **named in prose and given no bibitem**, per C23's own second hard rule. The over-claims are gone now; the attributions are owed, and `docs/qa/c23_run_001_paper_60.md` lists every one with what it is for. Also recorded there: `bottcher_spitkovsky2010` is cited by this paper for content nobody here has read (paywalled, no preprint); its title matches its use, so it is retained and flagged rather than dropped.

### Backing and gates

`tests/test_paper60_kms_attribution.py` gains the closed-form constant-and-phase guard (fire-tested three ways: drop the `pi/4`, drop `2^{-3/4}`, flip the `kR` exponent — all fire). `tests/test_paper60_preconditioner.py` gains the amplitude-floor and pricing guards, both written against the *hopeful* reading: the floor is tested across three whitenings, and the pricing guard asserts the composed amplitude exponent is ~2, i.e. strictly worse than untreated, so "preconditioning is a pure win" cannot pass. All gates PASS in scope `paper_60`.

**Process finding worth acting on:** C23's two best catches were its two newest claims. That argues for running it close to authorship rather than only at certification — a scope change to the criterion, and a PI call.

"""

C = Path("CHANGELOG.md")
s = C.read_text(encoding="utf-8")
A = "## [v5.11.0] - 2026-09-12"
assert s.count(A) == 1
C.write_text(s.replace(A, ENTRY + A), encoding="utf-8")
print("CHANGELOG: v5.11.1 inserted")

L = Path("CLAUDE.md")
s = L.read_text(encoding="utf-8")
old = "**Version:** v5.11.0 (September 12, 2026)"
assert s.count(old) == 1
s = s.replace(old, "**Version:** v5.11.1 (September 12, 2026)")
bullet = ("- **C23 run #1 = 6 of 8 claims already known (2026-09-12, v5.11.1):** two written this "
          "week; Prop D is Slater-Koster 1954, the chirp law is DLMF Bessel. Lever priced: one "
          "power of n. See CHANGELOG v5.11.1.\n")
A2 = "- **Composition wall is 3 axes; conditioning BREACHED (2026-09-12, v5.11.0):**"
assert s.count(A2) == 1
L.write_text(s.replace(A2, bullet + A2), encoding="utf-8")
print(f"CLAUDE.md: bumped + Sec.2 one-liner ({len(bullet.split())} words)")
