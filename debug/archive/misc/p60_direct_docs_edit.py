"""Direct-encoding docs: amend the superseded v5.11.1 verdict, sweep its live
dependents, and write v5.11.3."""
from pathlib import Path

# ---- claim_test_matrix: the pricing row is superseded, + a new row ----------
M = Path("docs/claim_test_matrix.md")
m = M.read_text(encoding="utf-8")
OLD = ("n^2` (product `n^2`); with a DIRECT encoding of `G` at the floor it would be `n`. So the "
       "lever is worth ONE power of n as priced, TWO if a direct block-encoding of `G` is found — "
       "the open item, sized by `||G|| = 0.372` against a composed `alpha = 3196` at n=160")
NEW = ("n^2` (product `n^2`); with a DIRECT encoding of `G` at the floor it is `n`. "
       "**SUPERSEDED 2026-09-12 (same day): the direct encoding is CONSTRUCTED (eq:ratio_symbol), "
       "so the lever is worth TWO powers of n, not one** — see the row below")
assert m.count(OLD) == 1, "pricing row anchor not found"
m = m.replace(OLD, NEW)

ROW = ("| 60 | sec:resource eq:ratio_symbol — the DIRECT block-encoding exists: `G`'s symbol is a "
       "ratio of two symbols with matching quadratic zeros, `ratio(s) = (1-j0(s))(s^2+(kR)^2)/"
       "(4s^2)`, whose sup is `0.3716 = ||G||` EXACTLY, so a circulant-embedded Toeplitz-minus-"
       "Hankel encoding carries `alpha = 0.372` against the composed 3196. `B != G` (11.5-11.8% "
       "in operator norm, NOT growing) but that is harmless: `cond(B^-1/2 G B^-1/2) -> 1.234` and "
       "`||P^-1/2 B^-1/2||` reaches the amplitude floor to 0.13%. Metric penalty `n^3 -> n` | "
       "`tests/test_paper60_direct_encoding.py``::test_ratio_symbol_sup_equals_the_norm_of_G` + "
       "``::test_direct_object_whitens_to_bounded_residual_conditioning`` + "
       "``::test_direct_whitening_attains_the_amplitude_floor`` (+ a slow n=160 row) | "
       "tracked `geovac/sturmian_sigma_law.py` + self-contained ratio quadrature | "
       "**NEW 2026-09-12** | BACKED-SOUND. The middle test asserts BOTH that the 12% discrepancy "
       "is real AND that it is harmless, so neither misreading passes; the floor test requires the "
       "ratio to TIGHTEN toward 1, which a fixed offset fails. Fire-tested three ways: a wrong "
       "ratio symbol, the composed whitening substituted for the direct one, and the sup compared "
       "against the composed alpha. rests on: eq:amplitude_floor (the spectrum-preservation "
       "argument is what makes `B != G` acceptable) |")
A = "| 60 | sec:resource (third lever, TRANSFER) —"
i = m.index(A)
M.write_text(m[:i] + ROW + "\n" + m[i:], encoding="utf-8")
print("claim_test_matrix: pricing row superseded, +1 direct-encoding row")

# ---- CHANGELOG: amend v5.11.1's verdict in place, then add v5.11.3 ----------
C = Path("CHANGELOG.md")
s = C.read_text(encoding="utf-8")
OLDC = ("**So the lever is worth one power of `n` as priced, two if a direct block-encoding of `G` "
        "is found.**")
NEWC = ("~~**So the lever is worth one power of `n` as priced, two if a direct block-encoding of "
        "`G` is found.**~~ **[SUPERSEDED 2026-09-12, v5.11.3: the direct encoding was constructed "
        "the same day — `G`'s symbol is a bounded ratio whose sup is exactly `||G||`, so the lever "
        "is worth TWO powers of `n` and the metric penalty goes `n^3 -> n`.]**")
assert s.count(OLDC) == 1, "v5.11.1 verdict sentence not found"
s = s.replace(OLDC, NEWC)

ENTRY = """## [v5.11.3] - 2026-09-12

**The open item closed the same day: `G` has a direct block-encoding, so the metric penalty goes `n^3 -> n`.** Probe `debug/p60_direct_encoding_probe.py`; backing `tests/test_paper60_direct_encoding.py`.

### The construction

v5.11.1 priced the preconditioner lever at one power of `n` and named the obstacle precisely: composing `P^-1/2` with `(I-C)` makes the subnormalization inherit `||P^-1/2||^2 ~ n^2`, wasting a factor `8.6e3` at `n = 160` on operation order alone, against `||G|| = 0.372`.

The waste is recoverable because **`G`'s symbol is a ratio of two symbols that vanish to the same order.** Both `1 - sigma` and `g = 2 + 2cos(chi)` have quadratic zeros at `chi = pi`, so the quotient is finite there (`(kR)^2/24`) and tends to `1/4` at `chi -> 0`. In `s = kR cot(chi/2)`:

    ratio(s) = (1 - j0(s)) (s^2 + (kR)^2) / (4 s^2),     ||ratio||_inf = 0.3716 = ||G||

The sup **coincides with `||G||`**, so the Toeplitz-minus-Hankel matrix `B` built from this symbol's own cosine coefficients is directly constructible and a circulant-embedded encoding of it carries `alpha = 0.372` — exactly the factor the composition threw away. New `eq:ratio_symbol`.

### Why `B != G` does not matter

`B` differs from `G` by the finite-section commutator: **11.5–11.8% in operator norm, and not growing with `n`.** That is real and is asserted as such. It is also harmless, because `B` enters only as a whitening, and by `eq:amplitude_floor` any `X` with `X^T S X = I` preserves the spectrum. Taking `X = P^-1/2 B^-1/2`:

| `n` | `cond(B)` | `\\|G-B\\|/\\|G\\|` | `cond(B^-1/2 G B^-1/2)` | `\\|X\\|`/floor |
|--:|--:|--:|--:|--:|
| 20 | 2.198 | 0.115 | 1.222 | 1.0097 |
| 80 | 2.228 | 0.118 | 1.232 | 1.0025 |
| 160 | 2.229 | 0.118 | **1.234** | **1.0013** |

The residual conditioning settles near `1.23`, so one further `O(1)`-degree transformation absorbs it; and `||X||` lands on the amplitude floor, tightening toward it rather than sitting at a fixed offset.

### What the metric factor now is

Three explicitly-known pieces, none growing with basis size: a **DST-I** (`O(log^2 N)` circuit, Klappenecker–Rötteler), a **diagonal** computed from the index, and a **degree-≈24 QSVT** on a directly-encoded Toeplitz-minus-Hankel matrix with `O(1)` subnormalization. On Table `tab:resource`'s model at `n = 160`: `alpha = 125.5`, `d_inv = 24.3`, product `3.1e3` against the untreated `3.9e7` — and the *scaling* is now `n` against `n^3`, so the ratio grows as `n^2`.

### Honest limits, restated

- **The circuit is cited, not compiled.** The sine transform is referenced and the circulant embedding is standard but not laid out, so these remain a resource model rather than a gate count. Said so in the paper's scope paragraph.
- **`s`-sector shared-scale bases**, `M = 2` and `M = 3`.
- **None of this recovers `l`-selection.** Proposition D is untouched: the metric is now cheap to *apply*, and the sparsity it destroys stays destroyed. The two open programs identified in v5.11.2 are unchanged in status — run cost moved, sparsity competitiveness is still closed by theorem.

### Dependents swept

The v5.11.1 verdict sentence is amended in place (Sec. 13.11 rule 9) and the `claim_test_matrix` pricing row marked SUPERSEDED, both pointing here; CLAUDE.md Sec. 2's "one power of n" corrected. Gates: C10 / C21 / C16 / C22 / C14 / escapes / titles / arxiv / duration PASS in scope `paper_60`. New guards fire-tested three ways.

"""
A2 = "## [v5.11.2] - 2026-09-12"
assert s.count(A2) == 1
C.write_text(s.replace(A2, ENTRY + A2), encoding="utf-8")
print("CHANGELOG: v5.11.1 verdict amended + v5.11.3 inserted")

# ---- walls register --------------------------------------------------------
W = Path("docs/walls/register.md")
w = W.read_text(encoding="utf-8")
OLDW = ("Remaining scope: `s`-sector "
        "shared-scale bases at `M = 2, 3`; the end-to-end resource claim still needs the "
        "sine-transform circuit and a block-encoding of `G` priced.")
NEWW = ("Remaining scope: `s`-sector "
        "shared-scale bases at `M = 2, 3`. **Resource claim closed 2026-09-12 (v5.11.3):** `G` "
        "has a DIRECT block-encoding — its symbol is a bounded ratio, `||ratio||_inf = 0.3716 = "
        "||G||`, so a circulant-embedded Toeplitz-minus-Hankel encoding carries `alpha = O(1)` "
        "instead of the composed `O(n^2)`; the resulting whitening reaches the amplitude floor to "
        "0.13% with residual conditioning 1.234. The metric penalty scales as `n` against the "
        "untreated `n^3`. What remains uncompiled is the circuit, not the construction.")
assert w.count(OLDW) == 1, "walls scope anchor not found"
W.write_text(w.replace(OLDW, NEWW), encoding="utf-8")
print("walls register: resource claim closed")

# ---- CLAUDE.md -------------------------------------------------------------
L = Path("CLAUDE.md")
l = L.read_text(encoding="utf-8")
old = "**Version:** v5.11.2 (September 12, 2026)"
assert l.count(old) == 1
l = l.replace(old, "**Version:** v5.11.3 (September 12, 2026)")

oldb = "chirp law = DLMF Bessel. Lever = one power of n. See CHANGELOG v5.11.1."
newb = "chirp law = DLMF Bessel. Lever priced. See CHANGELOG v5.11.1."
assert l.count(oldb) == 1
l = l.replace(oldb, newb)

bullet = ("- **Metric penalty n^3 -> n (2026-09-12, v5.11.3):** G's symbol is a bounded ratio, "
          "sup = ||G|| exactly, so a DIRECT encoding reaches the amplitude floor. See CHANGELOG "
          "v5.11.3.\n")
A3 = "- **Overcompleteness is ONE DIRECTION (2026-09-12, v5.11.2):**"
assert l.count(A3) == 1
L.write_text(l.replace(A3, bullet + A3), encoding="utf-8")
print(f"CLAUDE.md: bumped, v5.11.1 bullet corrected, +1 bullet ({len(bullet.split())} words)")
