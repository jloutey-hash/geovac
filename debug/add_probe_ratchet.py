r"""Ratchet C24's re-derivation probe against a recorded baseline.

Without this the probe reports the same seven hits every run -- all legitimate,
because the live papers that narrate the Lorentzian arc genuinely discuss those
topics. A report that never changes is a report nobody reads, and then the
mechanism built to answer "has this been attempted before?" is decoration.

C22 already uses this pattern, with the rule attached: a ratchet that hides its
own size is how debt becomes permanent. So the baseline size is printed every
run, and only NEW hits are surfaced.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P = "debug/qa/check_paper_retirement.py"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


edit(
    'INDEX = os.path.join(ROOT, "papers", "INDEX.md")',
    'INDEX = os.path.join(ROOT, "papers", "INDEX.md")\n'
    'PROBE_BASELINE = os.path.join(os.path.dirname(os.path.abspath(__file__)),\n'
    '                              "retirement_probe_baseline.json")',
    "baseline path")

edit(
    "import argparse\nimport collections",
    "import argparse\nimport collections\nimport json",
    "json import")

edit(
    '''    print()
    print("C. re-derivation probe (ADVISORY)")
    hits = check_c(rows, paths)
    if not hits:
        print("   [ok] no retired-paper trigger term is live in the corpus")
    for fname, term, where in hits:
        broad = len(where) > 8
        print("   %-42s %-38s %d live doc(s)%s"
              % (fname, '"%s"' % term[:36], len(where),
                 "  <- term too broad to be useful" if broad else ""))''',
    '''    print()
    print("C. re-derivation probe (ADVISORY, ratcheted)")
    hits = check_c(rows, paths)

    baseline = {}
    if os.path.exists(PROBE_BASELINE):
        try:
            baseline = json.load(io.open(PROBE_BASELINE, encoding="utf-8"))
        except ValueError:
            baseline = {}

    known, fresh = [], []
    for fname, term, where in hits:
        key = "%s|%s" % (fname, term)
        was = set(baseline.get(key, []))
        new_docs = sorted(set(where) - was)
        (fresh if new_docs else known).append((fname, term, where, new_docs))

    # The baseline size is ALWAYS printed. A ratchet that hides its own size is
    # how debt becomes permanent (the C22 rule).
    print("   [baseline] %d known hit(s) across %d recorded term(s) -- these are"
          % (len(known), len(baseline)))
    print("              the live papers that narrate an archived arc, not"
          " re-derivations")
    if not hits:
        print("   [ok] no retired-paper trigger term is live in the corpus")
    if not fresh:
        print("   [ok] no NEW locus has taken up a retired paper's topic")
    for fname, term, where, new_docs in fresh:
        print("   [NEW] %-38s %-34s now also in: %s"
              % (fname[:38], '"%s"' % term[:32], ", ".join(new_docs[:4])))
        print("         ^ that topic was attempted before. Read the archived"
              " paper before rebuilding it.")
    for fname, term, where, _ in known:
        if len(where) > 8:
            print("   [broad] %-38s %-34s %d docs -- term too generic to be useful"
                  % (fname[:38], '"%s"' % term[:32], len(where)))

    if args.update_baseline:
        snap = {"%s|%s" % (f, t): sorted(w) for f, t, w in hits}
        io.open(PROBE_BASELINE, "w", encoding="utf-8").write(
            json.dumps(snap, indent=2, sort_keys=True) + "\\n")
        print("   [baseline updated] %d term(s) recorded" % len(snap))
    hits = fresh''',
    "ratcheted probe reporting")

edit(
    '    ap.add_argument("--gate", default=None, help="accepted for symmetry; C24 is corpus-wide")',
    '    ap.add_argument("--gate", default=None, help="accepted for symmetry; C24 is corpus-wide")\n'
    '    ap.add_argument("--update-baseline", action="store_true",\n'
    '                    help="record the current probe hits as the baseline")',
    "--update-baseline flag")

edit(
    '''    print("RESULT: PASS (register intact; %d retirement candidate(s) and %d "
          "re-derivation hit(s) reported as advisory)" % (len(flagged), len(hits)))''',
    '''    print("RESULT: PASS (register intact; %d retirement candidate(s) and %d "
          "NEW re-derivation hit(s) reported as advisory)" % (len(flagged), len(hits)))''',
    "result line counts NEW hits")

with io.open(P, encoding="utf-8") as fh:
    t = fh.read()

applied, failed = [], []
for old, new, label in EDITS:
    if old in t:
        t = t.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

with io.open(P, "w", encoding="utf-8") as fh:
    fh.write(t)

print("applied %d of %d" % (len(applied), len(EDITS)))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
