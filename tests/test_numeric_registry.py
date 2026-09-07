"""Self-test for the numeric registry (C20).

The registry is only worth its maintenance cost if it actually discriminates.
A gate that fires on nothing reports PASS for a class it never checked --
this project has been bitten by exactly that twice, once by a C17 pattern
with doubled escapes and once by a thin-space regex mangled in a heredoc,
which left C20 blind to half the corpus's large numbers on its first run.

So these tests check the machinery, not the corpus:

  * every DERIVED expression resolves and recomputes
  * every RETIRED entry points at a registered replacement
  * the matcher sees BOTH LaTeX thin-space conventions
  * the gate FIRES on a planted retired value and stays SILENT on the
    corrected one (the two-way discrimination rule)
  * provenance is present on every entry, since a registry of unsourced
    numbers would launder guesses into authority
"""
from __future__ import annotations

import pathlib
import re
import sys

import pytest

QA = pathlib.Path(__file__).resolve().parents[1] / "debug" / "qa"
sys.path.insert(0, str(QA))

import numeric_registry as REG          # noqa: E402
import check_numeric_consistency as C21  # noqa: E402


def test_every_derived_expression_recomputes():
    """A DERIVED entry whose inputs are missing is dead weight."""
    for name, (expr, tol, desc) in REG.DERIVED.items():
        val = REG.evaluate(expr)
        assert isinstance(val, float)
        assert val == val, f"{name} evaluated to NaN"


def test_retired_values_point_at_registered_replacements():
    """A retired entry must be able to say what the locus SHOULD read."""
    for value, (key, need, forbid) in REG.RETIRED.items():
        REG.resolve(key)          # raises if unregistered
        assert need, f"retired {value} has no required context -- it would "\
                     f"fire on every bare occurrence of the numeral"
        re.compile(need)
        if forbid:
            re.compile(forbid)


def test_matcher_sees_both_thin_space_conventions():
    """The corpus writes 1{,}413 and 1\\,413; the gate must read both.

    Regression guard: a heredoc once halved this pattern's backslashes,
    silently reducing it to an escaped comma.  The gate then reported PASS
    while skipping every \\,-formatted number in the corpus.
    """
    brace = "1{,}413"
    thin = "1" + chr(92) + ",413"
    for text in (brace, thin):
        found = C21.NUM.findall(text)
        assert found, f"matcher missed {text!r}"
        assert C21._norm(found[0]) == 1413.0, \
            f"{text!r} normalised to {C21._norm(found[0])}, expected 1413"


def test_gate_fires_on_a_planted_retired_value(tmp_path):
    """Two-way discrimination: fires on the retired form, silent on the fixed."""
    # The gate matches required context against a +/-8 line WINDOW, not the
    # row alone -- which is how a bare table row fires at all, since the
    # "Pauli" context lives in the header above.  Model that faithfully: a
    # test that reads context differently from the gate proves nothing
    # about the gate.
    header = r"Molecule & $Q$ & GeoVac $N_{\mathrm{Pauli}}$ & Chawla RaH-18q \\"
    retired = r"LiH & 30 & 1\,413 & 12\,556 & 0.113$\times$ \\"
    corrected = r"LiH & 30 & 1\,501 & 12\,556 & 0.120$\times$ \\"

    def hits(row: str) -> int:
        ctx = header + "\n" + row
        n = 0
        for m in C21.NUM.finditer(row):
            v = C21._norm(m.group(1))
            if v is None or v not in REG.RETIRED:
                continue
            key, need, forbid = REG.RETIRED[v]
            if need and not re.search(need, ctx, re.I):
                continue
            if forbid and re.search(forbid, ctx, re.I):
                continue
            n += 1
        return n

    assert hits(retired) >= 1, "gate did not fire on the retired value"
    assert hits(corrected) == 0, "gate fired on the corrected value"


def test_every_entry_carries_provenance():
    """A value with no stated source cannot be audited."""
    for key, d in REG.MEASURED.items():
        assert d.get("provenance"), f"{key} has no provenance"
        assert d.get("convention"), f"{key} has no declared convention"
    for key, d in REG.CITED.items():
        assert d.get("source"), f"{key} has no source"


def test_aliases_are_distinct_from_their_base_value():
    """An alias records a DIFFERENT convention, not a duplicate."""
    for key, d in REG.MEASURED.items():
        for alias in (d.get("aliases") or {}):
            assert alias != d["value"], \
                f"{key}: alias {alias} equals the base value"


@pytest.mark.parametrize("scope", ["group4", "group6"])
def test_gate_passes_on_current_corpus(scope):
    """The corpus is consistent with the registry as it stands."""
    assert C21.check_derivations() == 0
    assert C21.check_retired(scope) == 0, \
        f"{scope}: live retired value(s) -- run "\
        f"python debug/qa/check_numeric_consistency.py --gate {scope}"


# ---------------------------------------------------------------------------
# The annotation layer (\gvq foreign keys)
# ---------------------------------------------------------------------------

def test_gvq_pattern_handles_braced_numerals():
    """The corpus writes large numbers as 21{,}607.

    Regression guard: the first version of this pattern used a [^}]* value
    argument, which truncated such a value mid-brace.  The annotation pass
    caught it via its own no-op invariant before it reached the corpus, but
    the CHECKER carried the same bug and would have mis-read every braced
    citation.
    """
    text = r"\gvq{he_ccpvtz_pauli}{21{,}607}"
    m = C21.GVQ.search(text)
    assert m, "pattern missed a braced value"
    assert m.group(1) == "he_ccpvtz_pauli"
    assert C21._norm(m.group(2)) == 21607.0


def test_annotation_check_accepts_rounding_but_not_drift():
    """A paper may quote fewer digits than the value was measured to.

    The check compares at the precision of the LITERAL, so correct display
    rounding passes while a retired value fails -- 11.18 is a rounding of
    11.175; 11.29 (the retired 1-norm) is not, at any precision.
    """
    canon = REG.resolve("he_n2_lambda")

    def renders(lit: str) -> bool:
        dec = len(lit.split(".")[1]) if "." in lit else 0
        return abs(float(lit) - round(canon, dec)) <= 10 ** (-dec) / 2 + 1e-9

    assert renders("11.175")
    assert renders("11.18"), "correct rounding must be accepted"
    assert renders("11.2")
    assert not renders("11.29"), "the retired value must still be rejected"


def test_every_annotation_key_is_registered():
    """An annotation naming a key the registry lacks is a dangling edge."""
    for scope in ("group4", "group6"):
        for path in C21._files(scope):
            txt = path.read_text(encoding="utf-8", errors="ignore")
            for m in C21.GVQ.finditer(txt):
                REG.resolve(m.group(1))   # raises on an unknown key


def test_annotations_are_a_rendering_no_op():
    r"""Stripping \gvq must leave valid LaTeX with the literal in place.

    Guards the invariant the annotation pass was built around: the macro
    renders its second argument, so a citation and a bare numeral are the
    same document.  If a future edit nests braces so that stripping changes
    the text, the PDFs move silently.
    """
    # Reuse the checker's own pattern rather than redefining it: a test that
    # parses \gvq differently from the gate cannot guard the gate.
    for scope in ("group4", "group6"):
        for path in C21._files(scope):
            txt = path.read_text(encoding="utf-8", errors="ignore")
            for m in C21.GVQ.finditer(txt):
                assert C21.GVQ.sub(lambda x: x.group(2), m.group(0)) \
                    == m.group(2), \
                    f"{path.name}: stripping {m.group(0)!r} changes the text"
                assert C21._norm(m.group(2)) is not None, \
                    f"{path.name}: {m.group(0)!r} literal does not parse"


# ---------------------------------------------------------------------------
# Stage 3: convention checking
# ---------------------------------------------------------------------------

def test_every_convention_string_parses():
    """A convention the parser cannot read drops silently out of check E.

    The `convention` field is written for a human, so nothing stops a new
    entry being phrased in a way the parser misses -- and it would then be
    exempt from convention checking without anyone noticing.  This test is
    what makes the free-text field safe.
    """
    unparsed = [k for k in list(REG.MEASURED) + list(REG.CITED)
                if not REG.parses(k)]
    assert not unparsed, (
        f"unrecognised convention on: {unparsed}. Add the phrasing to "
        f"numeric_registry._KINDS or reword the convention string."
    )


def test_measured_pauli_and_lambda_declare_an_identity_convention():
    """Our own counts must say whether the identity term is included.

    This is the pairing that value-checking cannot see: both 287 and 288 are
    correct, and only the mixture is wrong.
    """
    for key in REG.MEASURED:
        kind, ident = REG.family(key)
        if kind in ("pauli", "lambda"):
            assert ident in ("in", "out"), \
                f"{key} does not declare an identity convention"


def test_cited_baselines_do_not_assert_an_unknown_convention():
    """External numbers are convention-unknown unless the source says.

    Regression guard: the parser assigns "identity-included" whenever
    "non-identity" is absent, so an unstated convention silently became an
    asserted one -- and check E then reported a conflict it could not
    actually establish.  Registry rule 3: never register what you have not
    measured or cited.
    """
    for key in REG.CITED:
        kind, ident = REG.family(key)
        if kind in ("pauli", "lambda"):
            assert ident is None, (
                f"{key} asserts identity={ident!r}, but it is an external "
                f"baseline -- set identity=None unless the source states it"
            )


def test_check_e_flags_a_genuinely_mixed_table():
    """Two of OUR quantities in opposite conventions must be caught."""
    ours = [k for k in REG.MEASURED if REG.family(k)[0] == "pauli"]
    outs = [k for k in ours if REG.family(k)[1] == "out"]
    assert outs, "no non-identity Pauli quantity registered"
    # a table citing the same kind under both conventions is the defect;
    # verify the family() call that check E rests on distinguishes them
    k = outs[0]
    assert REG.family(k)[1] == "out"
    aliases = REG.MEASURED[k].get("aliases") or {}
    assert aliases, (
        f"{k} has no identity-included alias, so the gate cannot tell "
        f"'other convention' from 'wrong value' for it"
    )

def test_exemption_window_is_tighter_than_the_context_window():
    """A disclosure about one number must not exempt its neighbours.

    C21 used one +-8 window for BOTH finding require-context and accepting
    disclosure markers.  The second use let a single vintage note blanket
    its neighbourhood: Paper 14's tab:tc_composed carried three retired
    cells (334 / 556 / 778) reported clean because a note sat five lines
    below the table, and a "the previously printed 354.9" sentence exempted
    306.4, 373.4 and 66.0 alongside it.

    The split is only worth anything if it cannot be widened back without
    a test noticing.
    """
    assert C21.EXEMPT_WINDOW < C21.WINDOW, (
        "the exemption window must stay strictly tighter than the context "
        "window; equalising them restores the blanket-exemption failure"
    )

    retired = sorted(REG.RETIRED)[0]
    key, need, _forbid = REG.RETIRED[retired]

    far = "\n".join(
        ["% doc", f"An earlier value was retired here ({need.split('|')[0]})."]
        + ["filler"] * (C21.EXEMPT_WINDOW + 2)
        + [f"cell & {retired} & {need.split('|')[0]}"])
    near = f"The retired {retired} ({need.split('|')[0]}) was superseded."

    def live(text):
        lines = text.split("\n")
        hits = 0
        for i, line in enumerate(lines):
            if line.lstrip().startswith("%"):
                continue
            ex = "\n".join(lines[max(0, i - C21.EXEMPT_WINDOW):
                                  i + C21.EXEMPT_WINDOW + 1])
            if C21.EXEMPT.search(ex):
                continue
            ctx = "\n".join(lines[max(0, i - C21.WINDOW):i + C21.WINDOW + 1])
            for m in C21.NUM.finditer(line):
                if C21._norm(m.group(1)) != retired:
                    continue
                if need and not re.search(need, ctx, re.I):
                    continue
                hits += 1
        return hits

    assert live(far) > 0, (
        "a retired value far from the disclosure marker must FIRE -- "
        "otherwise the gate reports PASS on a class it cannot see"
    )
    assert live(near) == 0, (
        "a retired value on the disclosure marker's own line must stay "
        "silent -- otherwise every honest erratum becomes a defect"
    )


def test_no_duplicate_literal_keys_in_retired():
    """A retired-value dict literal must not name the same float twice.

    Python collapses `1.69` and `1.690` to one key, so the later entry
    silently overwrites the earlier one and its detection anchors vanish --
    the gate keeps reporting PASS on a class it can no longer see. That is
    the "guard that cannot fail" failure mode arriving through a dict
    literal rather than through an assertion.

    Found 2026-09-06: `exp_lambda_4pt` was written as both 1.69 and 1.690;
    the surviving entry had lost its `alpha|lambda` anchors.

    The check must read the SOURCE, not the parsed dict -- by the time the
    module is imported the duplicate is already gone.
    """
    import pathlib
    import re
    from collections import Counter

    src = pathlib.Path(__file__).resolve().parents[1] / "debug" / "qa" / "numeric_registry.py"
    text = src.read_text(encoding="utf-8")
    block = text[text.index("RETIRED"):]
    keys = re.findall(r"^\s{4}([0-9][0-9_.eE+-]*):\s*\(", block, re.M)
    assert keys, "could not locate any RETIRED literal keys -- parser drifted"

    counts = Counter(float(k) for k in keys)
    dupes = {v: [k for k in keys if float(k) == v]
             for v, n in counts.items() if n > 1}
    assert not dupes, (
        f"duplicate float keys in RETIRED silently discard detection "
        f"anchors: {dupes}"
    )
