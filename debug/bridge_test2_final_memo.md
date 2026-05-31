# Bridge Test 2 (final): Is GeoVac atomic correlation entropy MODULAR?

**Date:** 2026-05-30
**Final driver:** `debug/bridge_test2_final_v3.py` (calibrated; supersedes
`bridge_test2_final.py` (v1) and `bridge_test2_final_v2.py` (v2))
**Data:** `debug/data/bridge_test2_final.json`
**Tier:** debug-only. No `geovac/` or paper edits.

## The question

Sprint TD Track 5 (`debug/data/sprint_td_track5.json`) PSLQ'd the GeoVac atomic
correlation entropy `S_full(GS)` against a **spectral-side** mechanical basis
(π-powers, ζ-values, logs, Catalan/β) — a 12,312-form basis frozen before
testing — and returned a clean **NULL** for all three targets, concluding the
correlation entropy is "ring-orphan."

Open challenge: did Track 5 test the **wrong dictionary**? The wedge KMS
(thermal) entropy lives in **modular** arithmetic (`q = e^{-2π}`, `τ = i`
lemniscatic point) — `debug/bridge_test1_wedge_kms_entropy.py` shows the wedge
KMS entropy has the form `S = log(integer) + q-series`. If the atomic
correlation entropy were secretly modular, it would come back NULL against
Track 5's spectral basis **by construction**, and the "ring-orphan" conclusion
would be premature. This test settles it by re-running PSLQ against a basis that
ADDS the modular ring.

## Targets (full ~150-digit strings from sprint_td_track5.json)

| System | S_full(GS) (lead digits) |
|---|---|
| He n_max=3  | 0.040811051366… |
| Li⁺ n_max=3 | 0.011211717937… |
| He n_max=4  | 0.041879430089… |

## Basis (frozen before reading targets, W3 protocol), n = 20

Kept **Q-linearly clean** to avoid the documented v2-calibrated pathology
(multiplicatively-related transcendentals together → near-dependencies that
destabilize PSLQ). Only Q-linearly-independent generators are included.

- **Spectral group (14):** `1, π², π⁴, π⁶, ζ(3), ζ(5), log2, log3, log5,
  Catalan(=β(2)), β(4), √2, √3, √5`. (Deliberately NOT ζ(2)=π²/6, ζ(4)=π⁴/90,
  1/ζ(n) — Q-dependent on the π-powers; they would break cleanliness.)
- **Modular group ADDED (6), the τ=i / q=e^{-2π} ring:**
  `q = e^{-2π}`, `e^{-π}`, `q² = e^{-4π}`, `q³ = e^{-6π}` (genuine q-series
  terms — squaring is nonlinear, so these are linearly independent over ℚ),
  plus the single independent lemniscatic generator `log Γ(1/4)` and `log π`.
  - `log ϖ` (lemniscate const) and `log η(i)` (Dedekind eta) are BOTH ℚ-linear
    combinations of `{log Γ(1/4), log2, log π}` and were EXCLUDED to keep the
    basis clean — the lemniscatic ring is fully represented by `log Γ(1/4)`
    plus the spectral logs.

## Calibration (the load-bearing methodological point)

An intermediate run (v2) at `dps=320` flagged the 20-element basis as
"DEGENERATE" because PSLQ at `maxcoeff=1e6` manufactured a **spurious** internal
relation with coefficients ~10⁵. The 20 transcendentals are genuinely
ℚ-linearly independent; the "relation" is over-large-lattice noise: a 20-element
integer lattice with `maxcoeff=1e6` at 320 dps is dense enough that a fake
relation with ~10⁵ coefficients can close to 1e-90. The same noise produced
huge-coefficient pseudo-"hits" on the targets (resid ~1e-88 with coeffs ~10⁵) —
**not** identifications.

The final run (v3) fixes this with the curve-fit-audit discipline:
- **dps = 500** (margin: 20·log₁₀(1000) ≈ 60 ≪ 500).
- **Trustworthy-hit rule:** a PSLQ relation counts as a HIT only if residual
  `< 1e-120` **AND** max|coeff| ≤ 1000. A genuine entropy-in-the-ring relation
  would have small coefficients; large-coeff relations are lattice noise and are
  rejected.
- **Cleanliness check at maxcoeff 1e4** (where a genuine small-int dependency
  would show but spurious large-coeff ones cannot form).

This makes the instrument honest.

## Machinery validity (decision gate) — ALL PASS

- **Basis cleanliness:** CLEAN — PSLQ returns no small-coefficient internal
  relation at maxcoeff 1e4. (The v2 "DEGENERATE" flag was the precision/lattice
  artifact above, now eliminated.)
- **Positive controls (6, structurally different) — ALL HIT** at residual
  `< 1e-120` with small coefficients (≤7):
  - PC1 (spectral): `3·π² − 2·ζ(3) + 5·log2` ✓
  - PC2 (modular): `7·q − 3·logΓ(1/4) + 2·logπ` ✓
  - PC3 (mixed): `−ζ(5) + 4·Catalan − 6·e^{-π} + logπ` ✓
  - PC4 (cross): `2 − π⁴ + 3·logΓ(1/4) − 5·√3` ✓
  - **PC5 (modular-entropy-shaped, load-bearing):** `log3 + 5·q − 2·q²` —
    **exactly the wedge-KMS form** (`log(integer) + q-series`, per
    `bridge_test1`). Recovered exactly → the basis CAN detect a genuine modular
    entropy, so a NULL on the targets is a real statement.
  - PC6 (deep q-series, exercises q³): `q − 3·q² + 4·q³ − β(4)` ✓
- **Null control** (arbitrary structureless decimal): **NULL** — no false
  positive.

**MACHINERY VALID = True.**

## Sanity: spectral-only reproduces Track 5

Against the spectral group ALONE, all three targets return **NULL** —
reproducing Track 5's headline with a cleaner basis. The instrument is
consistent with the prior work.

## Target verdict

Against the **full spectral + modular** basis:

| Target | maxcoeff 1e3 | maxcoeff 1e4 | maxcoeff 1e6 |
|---|---|---|---|
| He n_max=3  | no relation | no relation | noise only (coeffs ~9.6×10⁵) |
| Li⁺ n_max=3 | no relation | no relation | noise only (coeffs ~8.3×10⁵) |
| He n_max=4  | no relation | no relation | noise only (coeffs ~7.9×10⁵) |

**No trustworthy hit anywhere.** At maxcoeff 1e3 and 1e4 PSLQ returns nothing at
all. At maxcoeff 1e6 the only relations found have coefficients ~10⁵–10⁶ at
residual ~1e-118 — over-large-lattice noise (rejected by the trust cap, exactly
as the same lattice produces for the v2 spurious internal relation). Cross-system
ratios (He/Li⁺, He_n3/He_n4, He_n4/Li⁺) also NULL.

## Conclusion

**Headline: `headline-dead-Track-5-hardened`.**

With **valid machinery** (all six positive controls HIT — including a
modular-entropy-shaped control that proves the basis can detect the wedge-KMS
form — null control NULL, basis genuinely clean), the GeoVac atomic correlation
entropy `S_full(GS)` returns **NULL against the modular ring as well as the
spectral ring.**

Track 5 did **not** test the wrong dictionary. The correlation entropy is not
spectral AND not modular: it is a genuine von Neumann entropy of the algebraic
1-RDM occupations, lying outside any closed-form transcendental ring tested —
spectral (master Mellin engine M1/M2/M3) or modular (q-series / lemniscatic).
This is consistent with the structural reading recorded in CLAUDE.md §2 (Sprint
TD Track 5): "GeoVac correlation entropy is genuinely von Neumann-only," joining
`K = π(B+F−Δ)`, the Sprint MR-C L2 constant `c`, and the Wolfenstein parameters
as a natural-but-ring-orphan constant.

The "Bridge Test 2" challenge is **resolved in the negative**: the bridge from
the wedge-KMS thermal entropy (which IS modular) to the atomic correlation
entropy does **not** go through. The bridge stops at correlation entropy.

### Honest scope

- The modular basis tested is the natural τ=i / q=e^{-2π} ring (q-powers +
  lemniscatic log Γ(1/4)). A genuinely different modular structure (e.g. built
  on a *different* CM point τ, or involving raw Γ(1/4)-powers with their
  near-dependency hazards) is not ruled out by this test — but there is no
  structural reason from the wedge-KMS side to expect one, and `bridge_test1`'s
  wedge entropy uses exactly q = e^{-2π}.
- The NULL is at trustworthy maxcoeff (≤1e4, where PSLQ returns nothing) at
  500 dps; the huge-coeff relations at 1e6 are demonstrably lattice noise (same
  signature as the spurious basis-internal relation), not identifications.
- This is a **real negative result**, stated plainly. It strengthens (does not
  weaken) Track 5.

### File trail
- `debug/bridge_test2_final.py` — v1 (minimal modular extension, n=14): VALID
  machinery, targets NULL. Superseded by v3 (richer modular ring).
- `debug/bridge_test2_final_v2.py` — v2 (n=20, dps=320, maxcoeff 1e6 uncapped):
  flagged DEGENERATE due to precision/lattice artifact. Diagnostic value:
  exhibited the over-large-lattice noise that v3's calibration removes.
- `debug/bridge_test2_final_v3.py` — **final** (n=20, dps=500, small-coeff trust
  cap): VALID machinery, basis CLEAN, targets NULL.
- `debug/data/bridge_test2_final.json` — final v3 output.
