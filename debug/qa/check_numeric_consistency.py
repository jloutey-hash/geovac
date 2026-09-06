"""C21 — numeric consistency gate.

Verifies the corpus against ``debug/qa/numeric_registry.py``.  Four checks,
each aimed at a defect class that survived repeated judgment review:

  **A. Derivations recompute.**  Every DERIVED entry is evaluated from its
  registry inputs.  This is the class no reviewer catches reliably, because
  it requires arithmetic rather than reading: lambda/Q, N^2, a ratio against
  a fixed baseline, an exponent fitted from a printed column.  When a base
  value moves, its derivatives silently do not.

  **B. Annotations agree with the registry.**  Load-bearing quantities are
  cited as ``\\gvq{key}{literal}`` -- the value AND a foreign key into the
  registry.  The literal renders (so the PDF stays self-contained and the
  .tex stays hand-editable), while the key supplies the referential
  integrity the corpus otherwise lacks.  This is what makes the registry a
  linked structure rather than a parallel list of numbers.

  **C. Retired values are not live**, reported *together with the registry
  key that replaces them* -- so the message says what the locus should say,
  not merely that it is wrong.

  **D. Salience report (advisory).**  Numerals appearing in multiple
  documents but NOT registered.  The maintenance mechanism: how the registry
  learns what it is missing, and why it cannot quietly fall behind.

WHY ANNOTATE RATHER THAN RESOLVE
--------------------------------
An earlier design would have papers carry ``\\gvq{key}`` alone and resolve
the value at build time.  Rejected for two reasons that outweigh the
convenience:

  * **The prose would decouple silently.**  A resolved macro updates the
    numeral and leaves the sentence around it untouched.  This corpus has
    live examples -- "far below" was true at 5.8x and false at 0.92x; "the
    d-block is the sparsest" reversed outright.  Resolution would have
    updated the number and left the claim standing, with no diff to catch
    it.  Our actual failure mode is value-to-prose coupling, and resolution
    makes it unobservable.
  * **The PDFs are DOI-stamped.**  A paper whose numbers resolve from a live
    registry is not self-contained; reproducing an archived artifact would
    need the registry state of that date.

Annotation keeps the literal in the file (visible in every diff, present in
the PDF) and adds the edge.  When a value moves, ``--index`` names every
locus that cites it, so the surrounding prose gets re-read rather than
silently updated.

Usage:  python debug/qa/check_numeric_consistency.py [--gate <branch>]
                                                     [--salience]
                                                     [--index <key>]
Exit 0 = PASS.  Self-test: tests/test_numeric_registry.py
"""
from __future__ import annotations

import math
import pathlib
import re
import sys
from collections import defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numeric_registry as REG  # noqa: E402
import os

# Shared --gate scope resolution (see debug/qa/qa_scopes.py): named
# scopes resolve to an explicit file list and every RESULT line carries
# the file count, so a gate can never report PASS on an empty scope.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qa_scopes  # noqa: E402


ROOT = pathlib.Path(__file__).resolve().parents[2]

SCOPES = {
    "group4": ["papers/group4_quantum_computing/paper_14_qubit_encoding.tex",
               "papers/group4_quantum_computing/paper_20_resource_benchmarks.tex",
               "papers/synthesis/group4_quantum_computing_synthesis.tex"],
    "group6": ["papers/group6_precision_observations/paper_26_entanglement.tex"],
    # Added 2026-08-30: the angular-ERI-density family (global-M_L vs
    # pair-diagonal) is owned here, and the M9 trace found the retired
    # convention live in both of these.
    "group3": ["papers/group3_foundations/paper_22_angular_sparsity.tex",
               "papers/synthesis/group3_foundations_synthesis.tex"],
}

EXEMPT = re.compile(
    r"retired|corrected|Corrected 20|vintage|artifact|previously|earlier|"
    r"withdrawn|superseded|dissolve|DISSOLVED|former|was measured|"
    r"pair-diagonal rule gave|rule gave|deliberately excluded|"
    r"does not survive|not the four-point|pair-diagonal composed|"
    r"pair-diagonal count|\\to |\\rightarrow|moving .{0,40} from", re.I)

# 8 rather than 3: a table's vintage note sits in its caption, which can be
# many rows above the stale cell.  Narrower windows reported disclosed rows.
WINDOW = 8

# ...but the SAME width must not be used to accept a disclosure marker.
# At +-8 a single retirement note blankets its whole neighbourhood: the
# pre-fix BeH2 paragraph's "the previously printed 354.9" exempted 306.4,
# 373.4 and 66.0 alongside it, and a vintage note five lines below
# tab:tc_composed exempted that table's three retired cells outright.
# A disclosure about one number is not a disclosure about its neighbours.
#
# Swept 8/4/3/2/1 against both directions (probe: does it fire on the
# pre-fix text; does it stay silent on the corrected corpus).  +-3 is the
# widest setting that gains sensitivity at zero false-positive cost.
EXEMPT_WINDOW = 3

# The corpus uses BOTH thin-space conventions -- 1{,}413 and 1\,413.  An
# earlier version matched only the first and was therefore blind to roughly
# half the large numbers in the corpus, including every row of tab:sunaga.
NUM = re.compile(
    r"(?<![\w.])(\d{1,3}(?:(?:\{,\}|\\,)\d{3})+|\d+\.\d+|\d{2,})(?![\w.])")

# \gvq{registry-key}{literal-value}
GVQ = re.compile(
    # One level of nested braces: the corpus writes large numbers as
    # 21{,}607, and a [^}]* value argument truncates them mid-brace.
    r"\\gvq\{([A-Za-z0-9_]+)\}\{((?:[^{}]|\{[^{}]*\})*)\}")


def _files(gate):
    """Files this gate examines for a target.

    SCOPE WIDENED 2026-08-31 (pre-flight for the re-certification sweep).
    C21 used to accept only the three hand-curated targets above -- the
    papers that happened to own registry entries -- and to hard-exit on
    every other target name.  Two consequences, both bad for a sweep:
    group1/2/5/synthesis/trunk and the single-paper targets had NO numeric
    guarding at all (the Stage-4 memo's open follow-on 3), and the
    unregistered-multi-document-numeral census could only ever count the
    corner of the corpus already being watched.

    Scope now resolves through debug/qa/qa_scopes.py, so C21 examines the
    target's full pre-registered document set.  SCOPES is retained below as
    the record of which papers OWN registry entries; it no longer bounds
    what gets checked.
    """
    files, warnings = qa_scopes.resolve(gate or "")
    qa_scopes.emit_warnings(warnings)
    if gate and not files:
        raise SystemExit(
            f"unknown gate {gate!r}; known scopes: "
            f"{', '.join(sorted(qa_scopes.SCOPES))}")
    return [pathlib.Path(f) for f in files]


def _norm(tok):
    """Parse a gvq literal to a float.

    Plain decimals parse directly.  For the trunk's symbolic constants
    (symbolic-literal support added 2026-09-06, FULL #8 follow-up) also
    reduce the common LaTeX forms -- \\frac{a}{b}, \\pi, \\pi^2, inline a/b,
    and a single coefficient variable Z (e.g. the Slater F^0 = 5Z/8, whose
    registered quantity is the coefficient 5/8) -- to a number, so C21 can
    verify symbolic-fraction annotations instead of reporting them UNPARSED.
    Returns None if the token cannot be reduced to a pure-arithmetic
    expression (then the caller reports UNPARSED, as before).
    """
    t = (tok.replace("{,}", "").replace("\\,", "")
         .replace("$", "").replace("~", "").strip())
    try:
        return float(t)
    except ValueError:
        pass
    s = t
    # \frac / \dfrac / \tfrac {A}{B} -> ((A)/(B))
    s = re.sub(r"\\[dt]?frac\{([^{}]*)\}\{([^{}]*)\}", r"((\1)/(\2))", s)
    s = (s.replace("\\pi^{2}", f"({math.pi ** 2})")
         .replace("\\pi^2", f"({math.pi ** 2})")
         .replace("\\pi", f"({math.pi})"))
    s = s.replace("Z", "")          # coefficient variable: 5Z/8 -> 5/8
    s = (s.replace("\\cdot", "*").replace("\\times", "*")
         .replace("^", "**").replace("{", "(").replace("}", ")").strip())
    if not re.fullmatch(r"[0-9.\s()*/+\-]+", s):
        return None
    try:
        return float(eval(s, {"__builtins__": {}}, {}))
    except Exception:
        return None


# ---------------------------------------------------------------------------

def check_derivations(verbose=True):
    """A. Every derived quantity recomputes from its declared inputs."""
    if verbose:
        print("A. derivations recompute from their inputs")
    bad = 0
    for name, (expr, tol, desc) in sorted(REG.DERIVED.items()):
        try:
            val = REG.evaluate(expr)
        except KeyError as e:
            print(f"   [BROKEN] {name}: {e}")
            bad += 1
            continue
        if verbose:
            print(f"   [ok] {name:24s} = {val:12.4f}   ({desc})")
    return bad


def check_annotations(gate, verbose=True):
    """B. Every \\gvq annotation's literal matches its registry key."""
    if verbose:
        print("\nB. annotations agree with the registry")
    bad = 0
    seen = 0
    for path in _files(gate):
        lines = path.read_text(encoding="utf-8", errors="ignore").split("\n")
        for i, line in enumerate(lines):
            for m in GVQ.finditer(line):
                key, literal = m.group(1), m.group(2)
                seen += 1
                try:
                    canon = REG.resolve(key)
                except KeyError:
                    print(f"   [UNKNOWN KEY] {path.name}:{i+1}  "
                          f"\\gvq{{{key}}} -- not in the registry")
                    bad += 1
                    continue
                lit = _norm(literal)
                if lit is None:
                    print(f"   [UNPARSED] {path.name}:{i+1}  "
                          f"\\gvq{{{key}}}{{{literal}}}")
                    bad += 1
                    continue
                aliases = (REG.MEASURED.get(key, {}).get("aliases") or {})
                cands = [canon] + [float(a) for a in aliases]
                # Two disjoint acceptance rules, chosen by literal KIND:
                #   * a SYMBOLIC literal (\frac, /, ^ ...) is EXACT -- it
                #     evaluates to the value, not a rounded display of it, so
                #     it must TIGHTLY equal a candidate.  (DELTA #12 fix: the
                #     earlier version OR'd a loose dec=0 display window ahead
                #     of this, which subsumed it -- \frac{4}{\pi} would have
                #     passed on a 2/pi key.  Gating by kind makes the tight
                #     match load-bearing for symbolic annotations.)
                #   * a plain NUMBER (int or decimal) may be a rounded display
                #     of a longer measured value; compare at its own precision.
                literal_clean = literal.rstrip("$~ ")
                is_symbolic = bool(re.search(r"[\\/^]", literal_clean))
                if is_symbolic:
                    ok = any(abs(lit - c) <= 1e-6 * max(1.0, abs(c))
                             for c in cands)
                else:
                    dec = len(literal_clean.split(".")[1]) \
                        if "." in literal_clean else 0
                    ok = any(abs(lit - round(c, dec)) <= 10 ** (-dec) / 2 + 1e-9
                             for c in cands)
                if not ok:
                    print(f"   [MISMATCH] {path.name}:{i+1}  \\gvq{{{key}}} "
                          f"says {lit}, registry says {canon}")
                    bad += 1
    if verbose:
        print(f"   [ok] {seen} annotation(s) checked, {bad} mismatched")
    return bad


def check_retired(gate, verbose=True):
    """C. No retired value live in a gated document."""
    if verbose:
        print("\nC. retired values are not live")
    live = 0
    for path in _files(gate):
        lines = path.read_text(encoding="utf-8", errors="ignore").split("\n")
        for i, line in enumerate(lines):
            if line.lstrip().startswith("%"):
                continue
            ctx = "\n".join(lines[max(0, i - WINDOW):i + WINDOW + 1])
            ex_ctx = "\n".join(
                lines[max(0, i - EXEMPT_WINDOW):i + EXEMPT_WINDOW + 1])
            if EXEMPT.search(ex_ctx):
                continue
            for m in NUM.finditer(line):
                v = _norm(m.group(1))
                if v is None or v not in REG.RETIRED:
                    continue
                key, need, forbid = REG.RETIRED[v]
                if need and not re.search(need, ctx, re.I):
                    continue
                if forbid and re.search(forbid, ctx, re.I):
                    continue
                try:
                    should = REG.resolve(key)
                except KeyError:
                    should = "?"
                print(f"   [LIVE] {path.name}:{i+1}  {m.group(1)} "
                      f"-> should be {should} ({key})")
                print(f"          {line.strip()[:96]}")
                live += 1
    if verbose and not live:
        print("   [ok] no live retired value in scope")
    return live


def check_salience(gate, show=False):
    """D. Unregistered but multi-document numerals (advisory)."""
    canon = REG.all_canonical_values()
    retired = set(REG.RETIRED)
    occ = defaultdict(lambda: defaultdict(int))
    for path in _files(gate):
        txt = re.sub(r"%.*", "", path.read_text(encoding="utf-8",
                                                errors="ignore"))
        # arXiv IDs (YYMM.NNNNN) and DOIs parse as decimals -- arXiv:
        # 2401.03705 was reported as the quantity "2401.04" living in three
        # documents.  MEASURED effect: 27 of 728 entries (3.7%).  The corpus
        # holds 388 arXiv IDs, but one only reaches this list if it appears
        # in TWO OR MORE documents, so most never enter it; the leak is the
        # small set of IDs shared across papers.  Cheap to strip, so kept --
        # but it is a minor cleanup, not a structural finding.  (An earlier
        # version of this comment claimed "every bibliography was feeding the
        # worklist"; that was one instance generalised without measurement.)
        txt = re.sub(r"arXiv:\s*\d{4}\.\d{4,5}(v\d+)?", " ", txt,
                     flags=re.I)
        txt = re.sub(r"10\.\d{4,9}/[^\s{}]+", " ", txt)
        for m in NUM.finditer(txt):
            v = _norm(m.group(1))
            if v is None or v < 10:
                continue
            occ[v][path.name] += 1
    unreg = {v: d for v, d in occ.items()
             if round(v, 6) not in canon and v not in retired and len(d) > 1}

    # The raw count is not a worklist.  Its only filters are value >= 10
    # and "appears in more than one document", which across 62 papers
    # admits every year (2026 leads with ~1,460 occurrences), every paper
    # and section number, and every qubit count.  Registering those would
    # violate registry rule 3 -- never register a value you have not
    # measured or cited.  So report the raw count, then the subset that is
    # actually MEASUREMENT-SHAPED, which is the part worth registering.
    shaped = {v: d for v, d in unreg.items() if _measurement_shaped(v)}

    print(f"\nD. salience report (advisory): "
          f"{len(unreg)} unregistered multi-document numerals, "
          f"of which {len(shaped)} are measurement-shaped")
    if show:
        print("   -- measurement-shaped (the actionable worklist) --")
        for v, d in sorted(shaped.items(),
                           key=lambda kv: -sum(kv[1].values()))[:40]:
            where = " ".join(f"{k[:9]}x{n}" for k, n in sorted(d.items()))
            print(f"   {v:>12g}  {sum(d.values()):3d}  [{where}]")
    else:
        print("   (run with --salience to list the measurement-shaped "
              "subset; registering those is the maintenance path, not a "
              "defect)")
    return len(unreg)


# Years, paper/section numbers and round counts are not quantities.  A numeral
# is treated as measurement-shaped only if it carries the fingerprints of a
# measured value: a decimal part, or a magnitude large enough that it is not a
# structural label.  Deliberately conservative -- a false negative here costs
# one unregistered value, a false positive costs a junk registry entry.
def _measurement_shaped(v: float) -> bool:
    if v != v or v in (float("inf"), float("-inf")):
        return False
    # Years and the project's own version-adjacent integers.
    if 1900 <= v <= 2100 and float(v).is_integer():
        return False
    # A decimal part means somebody measured it.
    if not float(v).is_integer():
        return True
    # Bare integers below 100 are overwhelmingly labels (paper numbers,
    # section numbers, qubit counts, quantum numbers) in this corpus.
    return v >= 100

def check_table_conventions(gate, verbose=True):
    """E. One identity convention per kind, within a table.

    Catches the class where every value is individually correct and only the
    pairing is wrong -- a table quoting one row's Pauli count including the
    identity term and another's excluding it.  Pure value-checking cannot
    see this, because both numbers are right.
    """
    if verbose:
        print("\nE. one identity convention per kind, within a table")
    bad = 0
    tables = 0
    for path in _files(gate):
        txt = path.read_text(encoding="utf-8", errors="ignore")
        for tm in re.finditer(r"\\begin\{table\*?\}(.*?)\\end\{table\*?\}",
                              txt, re.S):
            body = tm.group(1)
            keys = [m.group(1) for m in GVQ.finditer(body)]
            if not keys:
                continue
            tables += 1
            seen = {}
            for k in keys:
                try:
                    kind, ident = REG.family(k)
                except KeyError:
                    continue
                if kind is None or ident is None:
                    continue
                seen.setdefault(kind, {}).setdefault(ident, []).append(k)
            label = re.search(r"\\label\{([^}]*)\}", body)
            name = label.group(1) if label else "(unlabelled)"
            for kind, byident in seen.items():
                if len(byident) > 1:
                    detail = "; ".join(
                        f"{i}: {', '.join(sorted(set(v)))}"
                        for i, v in sorted(byident.items()))
                    print(f"   [MIXED] {path.name}  {name}  "
                          f"{kind} quoted in two conventions -- {detail}")
                    bad += 1
    if verbose:
        print(f"   [ok] {tables} annotated table(s) checked, {bad} mixed")
    return bad


def reverse_index(key):
    """Every locus citing a quantity -- the foreign key, queried.

    This is the payoff of annotating.  When a measured value moves, this
    names the prose that has to be RE-READ, not merely the numerals that
    have to be changed: the sentence around a number is a claim keyed to
    that number's magnitude, and it does not update itself.
    """
    try:
        canon = REG.resolve(key)
    except KeyError:
        print(f"unregistered quantity: {key}")
        return 1
    print(f"{key} = {canon}")
    d = REG.MEASURED.get(key) or REG.CITED.get(key) or {}
    for field in ("convention", "provenance", "source"):
        if d.get(field):
            print(f"  {field:11s}: {d[field]}")
    print("\n  cited at:")
    n = 0
    for path in _files(None):
        lines = path.read_text(encoding="utf-8", errors="ignore").split("\n")
        for i, line in enumerate(lines):
            if any(m.group(1) == key for m in GVQ.finditer(line)):
                print(f"    {path.name}:{i+1}  {line.strip()[:88]}")
                n += 1
    if not n:
        print("    (no annotated citation yet -- annotate to link it)")
    dep = [nm for nm, (expr, _, _) in REG.DERIVED.items()
           if re.search(rf"\b{re.escape(key)}\b", expr)]
    if dep:
        print(f"\n  feeds derivations: {', '.join(sorted(dep))}")
    return 0


def main():
    if "--index" in sys.argv:
        return reverse_index(sys.argv[sys.argv.index("--index") + 1])

    gate = None
    if "--gate" in sys.argv:
        gate = sys.argv[sys.argv.index("--gate") + 1]
    show = "--salience" in sys.argv
    scope = f"scope '{gate}'" if gate else "ALL scopes"
    print(f"numeric consistency gate (C21)   [{scope}]\n")

    broken = check_derivations()
    mismatched = check_annotations(gate)
    live = check_retired(gate)
    mixed = check_table_conventions(gate)
    check_salience(gate, show)

    print()
    if broken or mismatched or live or mixed:
        print(f"RESULT: FAIL ({broken} broken derivation(s), "
              f"{mismatched} annotation mismatch(es), "
              f"{live} live retired value(s), "
              f"{mixed} mixed-convention table(s) in {scope})")
        return 1
    print(f"RESULT: PASS (derivations recompute; annotations agree; "
          f"no live retired value; conventions consistent in {scope})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
