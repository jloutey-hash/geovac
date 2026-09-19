r"""Rewrite the stale C8 criteria in the two group2/synthesis definition-of-done
files (PI direction 2026-09-19), self-declaring the mid-run change.

WHY THIS IS SENSITIVE, stated in the artifact it modifies.  qa.md step 1 has the
PM verify the DoD frozen before any review, "so the goalposts are fixed before
any review begins" -- in EITHER direction.  This run verified group2.done.md
clean at commit 4bd5a36, then ran the CODE dimension on Papers 11/12/13 and
found a DoD-listed headline wrong by ~7.6 orders.  Claims, citations, synthesis
and the completeness-critic have NOT run.

So this edit changes the criteria BETWEEN dimensions of one run.  The PI
directed it; the cost is that a later verdict is graded against different
criteria than the CODE dimension was.  The mitigation is to make that
reconstructible rather than invisible: every changed line carries the date, the
measurement, and a marker; and the DoD's own Change log records the split.  A
moved goalpost that declares itself is auditable; one that does not is what
makes a PASS meaningless.

WRITTEN AS A FILE, not a heredoc: both target lines carry LaTeX/Unicode
(`H$_2^+$`, `$n_{\rm basis}=20$`, `H2+` as U+2082/U+207A), and four heredoc
escape incidents today ended with the CHANGELOG paragraph documenting the first
three being corrupted by the fourth.

Run:  python debug/qa/_rewrite_dod_criteria.py
"""
from __future__ import annotations

import io
import sys

MARK11 = "[retracted 2026-09-19: p11-h2plus-0002pct-retired]"

G2 = "docs/qa/group2.done.md"
SY = "docs/qa/synthesis.done.md"

# --- group2.done.md: the Paper 11 C8 watch-note (line ~100) ----------------
G2_OLD_P11 = (
    "  - **Paper 11:** H$_2^+$ **0.0002%** energy (spectral Laguerre, $n_{\\rm basis}=20$);\n"
    "    the FD 1.01% is an *artifact* (must be flagged, not a competing result).\n"
)
G2_NEW_P11 = (
    "  - **Paper 11 — REWRITTEN 2026-09-19 (PI direction), mid-run:** H$_2^+$ energy is\n"
    "    reproduced to **machine precision** (spectral Laguerre, $n_{\\rm basis}=20$):\n"
    "    measured $|E - E_{\\rm ref}| = 3.6\\times10^{-14}$~Ha against\n"
    "    $E_{\\rm ref} = -0.6026342144949$~Ha, i.e. at or below the precision to which\n"
    "    that reference is conventionally quoted, so **no percentage is a valid\n"
    "    criterion here** and none is stated. The FD 1.01% is an *artifact* (must be\n"
    "    flagged, not a competing result).\n"
    "    " + MARK11 + " The criterion this line carried until 2026-09-19 was\n"
    "    **0.0002%** ($=1.21\\times10^{-6}$~Ha), which understates the method by ~7.6\n"
    "    orders: even $n_{\\rm basis}=5$ (3.5e-6 %) is 57x better, and nine\n"
    "    quantity-x-reference combinations failed to reproduce it (nearest 4.64e-4 %,\n"
    "    a coarse-grid PES fit, off by 2.3x). Also retired with it: R_eq 2.005 bohr /\n"
    "    0.38% (a coarse-grid *fit* artifact; fine-grid fit gives 1.99726, +0.013%)\n"
    "    and the \"5000x accuracy improvement\" derived as 1.01%/0.0002%.\n"
    "    Registry: `p11_h2plus_err_ha`, `p11_h2plus_req_bohr`; C17 family\n"
    "    `p11-h2plus-0002pct-retired`.\n"
    "    **Grading note:** the CODE dimension of the 2026-09-19 run was graded against\n"
    "    the OLD criterion and reported its falsification; any later dimension is\n"
    "    graded against this one.\n"
)

# --- group2.done.md: the Paper 12 watch-note, stale TWICE -------------------
G2_OLD_P12_HEAD = "  - **Paper 12:** H$_2$ Neumann $V_{ee}$ recovers **92.4%** of $D_e$ vs 80.1% numerical.\n"
G2_NEW_P12_HEAD = (
    "  - **Paper 12 — watch-note refreshed 2026-09-19 (PI direction):** the live\n"
    "    headline is **99.81% of $D_e$ / 0.32 mHa** at the variational optimum\n"
    "    $\\alpha=1.40$, with **99.767% / 0.41 mHa** the fixed-$\\alpha=1.0$ ladder\n"
    "    endpoint (registry `p12_rebased_de_pct_aopt`, `p12_rebased_de_pct`; both real,\n"
    "    neither superseding the other). The 92.4% below is the $\\sigma$-only monomial\n"
    "    ceiling and 99.1% the azimuthal-restored monomial value — both historical\n"
    "    rungs, not the current criterion. **$\\mu>0$ is quadrature-FREE** since\n"
    "    v5.13.4 (closed-form $\\{E_1,\\gamma,\\ln\\}$ B-seeds, enforced by\n"
    "    `test_B_table_uses_closed_form_not_quadrature`); the \"spectral quadrature\"\n"
    "    clause below is superseded. **The 80.1% numerical comparator has NO backing\n"
    "    test** and its only guard tolerates a ~31 mHa wrong-direction swing — raised\n"
    "    to the PI 2026-09-19, unresolved.\n"
    "    H$_2$ Neumann $V_{ee}$ recovers **92.4%** of $D_e$ vs 80.1% numerical.\n"
)

# --- synthesis.done.md: the hierarchy-table criterion (line ~185) -----------
SY_OLD = (
    "- **Natural-geometry hierarchy** table: He 0.004% cusp / 0.022% raw / 0.19% CI; H₂⁺\n"
    "  0.0002%; H₂ 96.0% D_e; LiH R_eq 5.3%; etc. (Paper 13/11/15/17 certified values).\n"
)
SY_NEW = (
    "- **Natural-geometry hierarchy** table: He 0.004% cusp / 0.022% raw / 0.19% CI;\n"
    "  H₂⁺ **machine precision** (reference-limited; REWRITTEN 2026-09-19, PI\n"
    "  direction — " + MARK11 + " the 0.0002% this line carried\n"
    "  understated the method by ~7.6 orders, measured 3.6e-14 Ha at\n"
    "  $n_{\\rm basis}=20$); H₂ 96.0% D_e; LiH R_eq 5.3%; etc.\n"
    "  (Paper 13/11/15/17 values; H₂⁺ re-measured 2026-09-19).\n"
)

# --- group2.done.md Change log entry ---------------------------------------
G2_LOG_ANCHOR = "## Change log\n"
G2_LOG_ENTRY = (
    "## Change log\n"
    "- 2026-09-19 — **CRITERIA REWRITTEN MID-RUN (PI direction).** The `/qa group2 full`\n"
    "  run of 2026-09-19 (PI-scoped to Papers 11/12/13) verified this file frozen at\n"
    "  commit `4bd5a36`, ran the **CODE** dimension, and found the Paper-11 C8 headline\n"
    "  **0.0002%** falsified by ~7.6 orders (measured 3.6e-14 Ha at $n_{\\rm basis}=20$).\n"
    "  The Paper-11 and Paper-12 watch-notes above were rewritten at PI direction\n"
    "  **after** that dimension and **before** claims / citations / synthesis /\n"
    "  completeness-critic, which had not run (session rate limit). **Consequence, stated\n"
    "  so a later verdict stays auditable:** the CODE dimension was graded against the\n"
    "  OLD criteria and reported their falsification; any subsequent dimension is graded\n"
    "  against the NEW ones. The 2026-09-19 run's status is **INCONCLUSIVE** (four gating\n"
    "  dimensions unexercised) and is not a certification under either set.\n"
)


def apply(path: str, pairs: list[tuple[str, str]]) -> None:
    s = io.open(path, encoding="utf-8").read()
    n = 0
    for old, new in pairs:
        if old not in s:
            print(f"  {path}: ANCHOR NOT FOUND -> {old[:70]!r}")
            continue
        if s.count(old) != 1:
            print(f"  {path}: AMBIGUOUS ({s.count(old)}x) -> {old[:60]!r}")
            continue
        s = s.replace(old, new, 1)
        n += 1
    io.open(path, "w", encoding="utf-8", newline="").write(s)
    print(f"  {path}: {n}/{len(pairs)} applied")


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    apply(G2, [(G2_OLD_P11, G2_NEW_P11),
               (G2_OLD_P12_HEAD, G2_NEW_P12_HEAD),
               (G2_LOG_ANCHOR, G2_LOG_ENTRY)])
    apply(SY, [(SY_OLD, SY_NEW)])

    # Verify C17 now sees nothing live in either file.
    import importlib.util
    import re
    spec = importlib.util.spec_from_file_location(
        "chn", "debug/qa/check_headline_numbers.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    e = [x for x in mod.REGISTRY
         if x["id"] == "p11-h2plus-0002pct-retired"][0]
    pat, req, exm = (re.compile(e["pattern"]), re.compile(e["require_nearby"]),
                     re.compile(e["exempt_if_nearby"]))
    win = getattr(mod, "WINDOW", 3)
    for f in (G2, SY):
        lines = io.open(f, encoding="utf-8").read().splitlines()
        live = []
        for k, ln in enumerate(lines):
            if not pat.search(ln):
                continue
            ctx = "\n".join(lines[max(0, k - win):k + win + 1])
            if req.search(ctx) and not exm.search(ctx):
                live.append(k + 1)
        print(f"  {f}: live occurrences now -> {live or 'none'}")


if __name__ == "__main__":
    main()
