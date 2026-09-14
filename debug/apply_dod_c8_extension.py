"""Extend the Paper 60 DoD C8 headline list to cover v5.11.0..v5.11.4.

PI-confirmed 2026-09-12 before any reviewer was dispatched. This COMPLETES the
pre-registration for claims that did not exist when the list was last written
(2026-09-11, v5.10.18); it does not relax any existing goalpost. Adding claims
to be checked makes the gate harder.

Per the block's own rule, values are NOT written here -- each headline names the
claim, its tier, and the gate that owns its number.

Idempotent.
"""
from __future__ import annotations

import sys

D = "docs/qa/paper_60.done.md"
MARKER = "15. **The conditioning lever"

ANCHOR = """## Seeding plan (worktree only; never touches the real corpus)
"""

NEW = """15. **The conditioning lever: a symbol zero of known order and location [SYMBOLIC +
    MEASURED].** The ill-conditioning is a zero of known order and location, so a
    band-Toeplitz preconditioner in Serra's sense removes it; the matching polynomial is
    EXACTLY tridiagonal (its Hankel part vanishes identically) and DST-I diagonalizable in
    closed form; `cond(G)` is FLAT in `n` against the raw `n^2` growth. Backing
    `tests/test_paper60_preconditioner.py`. **Legitimacy leg, which must stay attached:**
    any `X` with `X^T S X = I` preserves the generalized spectrum, so this is not a change
    of problem. **Scope that must stay attached:** `s`-sector shared-scale, `M = 2, 3`; and
    it does NOT recover `l`-selection — Proposition D is untouched. Prose implying the lever
    restores sparsity = MATERIAL.
16. **The lever reaches water's `A_1`, and the rotation is what does the work [MEASURED].**
    The symmetry-inequivalent-centre case that defeats the gerade lever (C8#9); raw growth
    reproduces `N^1.97` independently while the preconditioned column is bounded. **The
    CONTROL is load-bearing:** naive block-diagonal preconditioning WITHOUT the
    null-direction rotation leaves the growth intact. Reporting the gain without the control
    = MATERIAL. Backing `tests/test_paper60_preconditioner.py`.
17. **The M-centre null space is geometry-independent; its RATES are not [MEASURED,
    SCOPE].** At the symbol point every block symbol tends to `j0(0)`, so the null space is
    the constants' orthogonal complement for every arrangement — but the order at which each
    direction opens is governed by `rank(P D2 P)`, which is ONE for collinear centres, so a
    linear polyatomic opens at orders 2, 4, ..., 2(M-1). The lever is established for `M = 2`
    and NON-COLLINEAR `M = 3`; **claiming it for a linear polyatomic = MATERIAL.** Backing
    `tests/test_paper60_mcentre_orders.py`. External: Batenkov-Demanet-Goldman-Yomdin.
18. **The amplitude floor [SYMBOLIC].** `eq:amplitude_floor` — any `X` with `X^T S X = I`
    satisfies `||X|| = ||S^{-1/2}||` EXACTLY, independent of the factorization, so no
    whitening can lower the block-encoding subnormalization. Backing
    `tests/test_paper60_preconditioner.py`. Prose implying a factorization buys amplitude =
    MATERIAL.
19. **The direct block-encoding of `G` [SYMBOLIC + MEASURED].** `eq:ratio_symbol` — `G`'s
    symbol is a bounded ratio of two symbols with zeros of the same order, and its sup
    COINCIDES with `||G||`, so a circulant-embedded Toeplitz-minus-Hankel encoding carries
    `O(1)` subnormalization and the metric penalty scales as `n` rather than `n^3`. Backing
    `tests/test_paper60_direct_encoding.py`. **Two honest limits that must stay attached:**
    the circuit is CITED, not compiled (a resource model, not a gate count); and `B != G` by
    a stated operator-norm fraction that does not grow with `n`, harmless only because `B`
    enters as a whitening and the amplitude floor (C8#18) makes any such `X`
    spectrum-preserving. Presenting this as a compiled circuit = MATERIAL.
20. **Overcompleteness is ONE DIRECTION, not a property of the basis [MEASURED].** The
    Bessel deficit of a displaced Sturmian against the one-centre span PLATEAUS and does not
    tend to zero, so the one-centre set is measurably far from complete IN THE MOLECULAR
    METRIC. **The 2026-09-11 frames reading — "overcompleteness is the price of one-centre
    completeness" — is WITHDRAWN;** re-asserting it = MATERIAL. Ron-Shen is the surviving
    mechanism. The paired claim form is load-bearing: the gap must collapse WHILE the deficit
    does not, which no single-sided bug satisfies. Backing
    `tests/test_paper60_one_direction.py`. **Consequent caveat that must stay attached:** the
    basis's motivating completeness is in the ATOMIC metric, not the molecular one.
21. **`eq:sigma_law` is Kac-Murdock-Szego [PRIOR ART].** Not derived here; the corpus claims
    only the IDENTIFICATION of the metric as such a finite section, Toeplitz minus Hankel
    with the stated symbol. Re-claiming the asymptotic = MATERIAL. Backing
    `tests/test_paper60_kms_attribution.py`. **Residue leg (2026-09-12):** the ~1% figure is
    DOMINATED by the `n -> n+1` grid convention, with a genuine `O(1/n)` term surviving both
    conventions; calling it purely an asymptotic tail = MATERIAL-SMALL.
22. **The `l`-selection loss is NOT a conditioning effect [SYMBOLIC].** Proposition D — a
    block-diagonal congruence cannot orthogonalize a metric that is not block diagonal, at
    EVERY `cond(S) > 1`, and it does not relax as `cond(S) -> 1+`. **Re-attributed 2026-09-12
    (C23 run #1): this is Löwdin symmetry-preservation specialized to the `l` grading, known
    since Slater-Koster (1954); what the paper claims is the `l`-vs-`m` application.**
    Presenting Proposition D as a new proposition = MATERIAL. Backing
    `tests/test_paper60_kms_attribution.py`.
23. **`eq:chirp_decay` is a Bessel asymptotic [SYMBOLIC + PRIOR ART].** Constant AND phase
    from DLMF, no stationary-phase argument needed; the `pi/4` is the BRANCH phase of the
    square-root prefactor, not a stationary-phase signature; `sum|c_j|` CONVERGES while
    `sum j|c_j|` diverges, which is the Böttcher-Widom hypothesis and a different condition.
    Backing `tests/test_paper60_kms_attribution.py`. Asserting a stationary-phase origin, or
    conflating the two sums, = MATERIAL-SMALL.
24. **Transcendental tagging of this section [SYMBOLIC + MEASURED].** Both constants are
    calibration-tier M2 and the tagging is PROVENANCE ONLY. **Prior art that must be
    credited (C23 run #3, 2026-09-12):** the Dirichlet-eigenvalue reading of the constant,
    the extremal Wirtinger-Sobolev problem behind it, and the independence of the constant
    from the rest of the symbol are ALL Böttcher-Widom's, in a source the paper already
    cites. **The Bessel-free measurement is a change of REPRESENTATION, not an independent
    route** — presenting it as independent corroboration = MATERIAL. **The removability
    corollary is WITHDRAWN** ("truncation-side prices are matrix-level and reachable;
    continuum-side prices are symbol-level and untouchable"): both halves are false, because
    the preconditioner is built FROM the symbol and preconditioning IS a congruence.
    Re-asserting it = MATERIAL. Backing `tests/test_paper60_contraction_window.py`.
25. **The minimiser carries the antipodal parity [MEASURED].** The band minimiser is the
    Dirichlet ground state in the band index ONLY with the alternating factor; without it
    the two are EXACTLY ORTHOGONAL, so the bare statement is not an approximation of the
    right one. Stating it bare = MATERIAL-SMALL. Backing
    `tests/test_paper60_contraction_window.py`.
26. **The law is carried by the TRANSLATION, not the metric [MEASURED + PRIOR ART].** The
    generalized symbol is the quotient, so a smooth positive radial weight cancels; the
    vanishing-weight CONTROL is the load-bearing half and reporting the agreement without it
    = MATERIAL. **Prior art (2026-09-12): Ahmad et al. give this EXACTLY for the tau/DST-I
    algebra — the structure this paper works in — so the measurement confirms a theorem
    rather than establishing one;** claiming it as novel = MATERIAL. **Scope that must stay
    attached:** this is NOT `V_0`-independence, since a position-space-local `V_0` acts by
    convolution and leaves the class. Backing `tests/test_paper60_contraction_window.py`.
27. **The translation identification is NOT ours [PRIOR ART].** Shibuya-Wulfman's own
    abstract builds the molecular operator from one unitary per nucleus; the explicit
    group formulation and the Coulomb-Sturmian translation-operator reading are both later
    published work. **What survives as ours is the SYMBOL.** Re-claiming the translation
    reading = MATERIAL. **Companion resolution that must stay intact:** Monkhorst-Jeziorski's
    "no linear dependence" and this paper's measured conditioning are the SAME pencil and
    both true — they never INVERT the overlap, and GeoVac inverts because a block-encoding
    wants a standard Hermitian eigenproblem, so the exposure belongs to the ENCODING
    REQUIREMENT. **Provenance cap:** that paper's two-page body is UNREAD (closed, no
    repository copy); the mechanism is reconstructed from the lineage and the paper says so.
    Dropping that cap = MATERIAL-SMALL.

> **Pre-registration completed 2026-09-12 (PI-confirmed) for v5.11.0–v5.11.4.** Headlines
> 15–27 cover the preconditioner lever, the water `A_1` transfer, the M-centre rate scoping,
> the amplitude floor, the direct block-encoding, the one-direction correction, and the
> five prior-art re-tierings of 2026-09-12. This ADDS claims to be checked; no existing
> goalpost was relaxed. **Standing caution carried forward:** three of these (20, 24, 26)
> record a WITHDRAWN reading, and a withdrawn reading re-surfacing is the corpus's most
> frequent defect class — C16 entries are owed for 24's removability corollary and 20's
> frames reading.

"""


def main() -> int:
    with open(D, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(ANCHOR) != 1:
        print(f"anchor count={t.count(ANCHOR)}")
        return 2
    with open(D, "w", encoding="utf-8") as fh:
        fh.write(t.replace(ANCHOR, NEW + ANCHOR))
    print("applied: C8 headlines 15-27 pre-registered")
    return 0


if __name__ == "__main__":
    sys.exit(main())
