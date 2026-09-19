r"""Extend C21's `_KINDS` with a `residual` kind, rather than mis-labelling.

WHY EXTEND INSTEAD OF REWORD.  `tests/test_numeric_registry.py::
test_every_convention_string_parses` failed after I added
`p11_h2plus_err_ha`, naming three keys: my new one plus the pre-existing
`p12_rebased_err_mha` and `p12_rebased_err_mha_aopt`.

The recognised accuracy-side needles are "% of d_e", "mha gained" and
"cond(s)".  None of them describes an ABSOLUTE ENERGY RESIDUAL -- neither my
`|E - E_ref|` in Ha, nor the two p12 keys' "mHa above exact".  Rewording a
convention string to match a needle it does not mean would quiet the gate by
registering something false, which is exactly what `provenance is not
decoration; a registry of guesses launders them into authority` forbids.

So the honest fix is a new kind.  Like `density`, `constant` and `accuracy`,
a residual carries no identity-in/out convention -- what its convention
string MUST carry is the reference it is measured against and the basis /
truncation, because an absolute residual quoted against a reference of fewer
digits than the residual itself is meaningless (the Paper-11 defect that
prompted this: a 1e-14 Ha residual quoted against a 13-digit reference).

Written as a file, not a heredoc: this edit contains regex-free text but the
same session already shipped two doubled-escape patterns through heredocs,
and the standing rule is that edits to gate internals go through Write.

Run:  python debug/qa/_fix_c21_residual_kind.py
"""
from __future__ import annotations

import io

PATH = "debug/qa/numeric_registry.py"

ANCHOR = '    ("accuracy", ("% of d_e", "mha gained", "cond(s)")),\n'

ADDITION = '''    # Absolute energy RESIDUALS (Paper 11 H2+, Paper 12 re-based H2;
    # registered 2026-09-19).  |E - E_ref| in Ha, or mHa above exact.  No
    # identity-in/out convention applies.  What the convention string must
    # carry instead is (a) the REFERENCE the residual is measured against and
    # (b) the basis/truncation -- because a residual smaller than the
    # reference's own quoted precision is not a measurement of the method.
    # That is precisely the Paper-11 defect this kind was added for: a
    # 3.6e-14 Ha residual reported as "0.0002%" against a reference given to
    # 13 digits.  Added rather than reworded: the three keys it covers
    # (p11_h2plus_err_ha, p12_rebased_err_mha, p12_rebased_err_mha_aopt) are
    # genuinely residuals, and relabelling them "% of d_e" to satisfy the
    # parser would register a false convention to quiet a gate.
    ("residual", ("in ha, against", "mha above exact", "residual")),
'''


def main() -> None:
    s = io.open(PATH, encoding="utf-8").read()
    if '("residual",' in s:
        print("residual kind already present; nothing to do")
        return
    if ANCHOR not in s:
        raise SystemExit("anchor line not found -- _KINDS layout changed")
    s = s.replace(ANCHOR, ANCHOR + ADDITION, 1)
    io.open(PATH, "w", encoding="utf-8", newline="").write(s)

    # Prove it parses the three keys it was written for, and that it did not
    # swallow anything else.
    import importlib.util
    spec = importlib.util.spec_from_file_location("nr", PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for key in ("p11_h2plus_err_ha", "p12_rebased_err_mha",
                "p12_rebased_err_mha_aopt"):
        if key in mod.MEASURED:
            print(f"  {key}: parses={mod.parses(key)} kind={mod.family(key)[0]}")
    unparsed = [k for k in list(mod.MEASURED) + list(mod.CITED)
                if not mod.parses(k)]
    print(f"  remaining unparsed conventions: {unparsed}")


if __name__ == "__main__":
    main()
