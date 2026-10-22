"""C25 -- claim-registry gate (claims as linked objects).  See claim_registry.py.

Usage:
  python debug/qa/check_claim_registry.py                 # gate: FAIL on owner drift, front-door
                                                          #   or annotated-locus violations
  python debug/qa/check_claim_registry.py --mentions ID   # DELTA reading list for one claim
  python debug/qa/check_claim_registry.py --changed-since 2026-09-30
  python debug/qa/check_claim_registry.py --selftest      # proves every check FIRES
Exit code 0 = PASS, 1 = FAIL.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import claim_registry as R  # noqa: E402

GVC_OPEN = re.compile(r"\\gvc\{([a-z0-9-]+)\}\{")
GVC_MD = re.compile(r"<!--\s*gvc:([a-z0-9-]+)\s*-->")
# A negation exempts an overclaim only when it is ADJACENT to it ("not a conjecture",
# "never derived", "rather than a theorem").  The earlier window rule ("not" anywhere in
# the 24 chars before) let "not only checked but INTERNAL THEOREM" and "no: K is
# derived" through (DELTA #24 code review).
_LIST = r"(?:[\w-]+(?:,|\s+or|\s+nor)\s+){0,4}"
# A comma-separated negated list exempts an overclaim only if the list CLOSES with or/nor --
# before the overclaim ("not a conjecture or theorem") or after it ("never a derivation,
# conjecture, prediction or theorem").  "is not conditional, INTERNAL THEOREM" is not a
# list at all and must fire (DELTA #25 code review: the v5.18.1 list rule let it through).
_LIST_CLOSED_BEFORE = re.compile(r"\b(?:or|nor)\s+" + "[\"'`\u201c\u2018]*" + r"$", re.I)
_LIST_CLOSED_AFTER = re.compile(r"^[\"'`\u201d\u2019]*\s*(?:,\s*[\w-]+(?:\s+[\w-]+)?\s*){0,4},?\s*(?:or|nor)\s+\w", re.I)
_ADJ_NEG = re.compile(r"\b(?:not|never|nor|neither)\s+[\"'`\u201c\u2018]*(?:(?:a|an|the|be|been|is|are|was|were|presented\s+as)\s+){0,2}[\"'`\u201c\u2018]*$"
                      r"|\b(?:withdrawn|downgraded|retired|demoted)\s+(?:as|from)\s+(?:(?:a|an|the)\s+)?[\"'`\u201c\u2018]*$"
                      r"|\bno\s+longer\s+(?:(?:a|an|the)\s+)?[\"'`\u201c\u2018]*$"
                      r"|\b(?:rather\s+than|instead\s+of)\s+(?:(?:a|an|the)\s+)?[\"'`\u201c\u2018]*$", re.I)
_Q = "[\"'`\u201c\u2018]*"
NEGATION = re.compile(r"(\b(?:not|never|nor)\s+(?:(?:a|an|the|be|been|is|are|was|were)\s+)?" + _LIST + _Q +
                      r"|\bno\s+longer\s+(?:(?:a|an|the)\s+)?" + _Q +
                      r"|\b(?:rather\s+than|instead\s+of)\s+(?:(?:a|an|the)\s+)?" + _LIST + _Q + r")$", re.I)
ARROW = re.compile(r"\s*(→|->|-->|\\to\b|\$\\to\$)\s*")
OWNER_WINDOW = 8


def gvc_spans(t: str):
    """Yield (id, span, start) for every \\gvc{id}{span}, braces balanced to any depth.
    An unterminated \\gvc yields span None so the caller can report it."""
    for m in GVC_OPEN.finditer(t):
        depth, i = 1, m.end()
        while i < len(t) and depth:
            if t[i] == "\\" and i + 1 < len(t):
                i += 2
                continue
            depth += {"{": 1, "}": -1}.get(t[i], 0)
            i += 1
        yield m.group(1), (t[m.end():i - 1] if depth == 0 else None), m.start()

MENTION_GLOBS = ["papers/group*/*.tex", "papers/synthesis/*.tex", "papers/INDEX.md",
                 "docs/claims_register.md", "docs/claim_test_matrix.md",
                 "tests/*.py", "geovac/*.py"]


def _read(rel: str) -> str:
    with open(os.path.join(ROOT, rel), encoding="utf-8", errors="replace") as f:
        return f.read()


def _c(p: str) -> re.Pattern:
    return re.compile(p, re.I)


def _norm(span: str) -> str:
    """Strip LaTeX formatting wrappers and ties so "\\emph{not}~a theorem" reads
    "not a theorem" and "mostly \\textbf{forced}" reads "mostly forced"."""
    span = re.sub(r"\\(?:emph|textbf|textit|textrm|textsc|textsf|texttt|textup|underline|mathrm|text)\{", "", span)
    span = re.sub(r"(?<!\\)[{}]", "", span)
    return span.replace("~", " ")


def check_text_span(cid: str, span: str) -> list[str]:
    """Hold one restatement span to the claim's rule; return problems."""
    c = R.CLAIMS[cid]
    probs = []
    # strip LaTeX formatting wrappers so "mostly \textbf{forced}" reads "mostly forced"
    span = _norm(span)
    if not _c(c["qualifier"]).search(span):
        probs.append(f"missing qualifier /{c['qualifier']}/")
    m = first_live_overclaim(cid, span)
    if m:
        probs.append(f"overclaim '{m}'")
    return probs


def first_live_overclaim(cid: str, span: str):
    c = R.CLAIMS[cid]
    for m in _c(c["overclaim"]).finditer(span):
        # an ADJACENT negation ("not a conjecture", "never derived") is the corrected form;
        # a negated LIST counts only if the list closes with or/nor (before or after)
        pre = span[max(0, m.start() - 40):m.start()]
        if _ADJ_NEG.search(pre):
            continue
        if NEGATION.search(pre) and (_LIST_CLOSED_BEFORE.search(pre) or _LIST_CLOSED_AFTER.match(span[m.end():m.end() + 80])):
            continue
        # a tier TRANSITION record is history only when the arrow points AT the current
        # tier ("conjecture→observation"); "INTERNAL THEOREM -> now proven" is not
        a = ARROW.match(span, m.end())
        if a and _c(c["qualifier"]).match(span, a.end()):
            continue
        return m.group(0)
    return None


def owner_checks() -> list[str]:
    out = []
    for cid, c in R.CLAIMS.items():
        try:
            t = _read(c["owner"])
        except OSError:
            out.append(f"[{cid}] owner file missing: {c['owner']}")
            continue
        lines = t.split("\n")
        hits = [i for i, l in enumerate(lines) if _c(c["owner_anchor"]).search(l)]
        if not hits:
            out.append(f"[{cid}] owner anchor /{c['owner_anchor']}/ not found in {c['owner']}")
            continue
        w = c.get("owner_window", OWNER_WINDOW)
        wins = [_norm("\n".join(lines[max(0, i - w):i + w + 1])) for i in hits]
        # owner_all: EVERY anchor locus must carry the tier (one locus cannot vouch for another)
        quant = all if c.get("owner_all") else any
        if not quant(_c(c["owner_tier"]).search(x) for x in wins):
            out.append(f"[{cid}] owner tier /{c['owner_tier']}/ not within {w} lines of "
                       f"{'every' if c.get('owner_all') else 'the'} anchor in {c['owner']}")
        # the owner must not itself carry the overclaim next to its anchor (a registry that
        # only checks a tier word is present cannot see the owner drifting STRONGER)
        for x in wins:
            bad = first_live_overclaim(cid, x)
            if bad:
                out.append(f"[{cid}] owner carries its own overclaim '{bad}' near the anchor in {c['owner']}")
                break
    return out


def front_door_checks() -> list[str]:
    out = []
    for rel in R.FRONT_DOOR:
        try:
            t = _read(rel)
        except OSError:
            continue
        for n, line in enumerate(t.split("\n"), 1):
            for cid, c in R.CLAIMS.items():
                if _c(c["mention"]).search(line):
                    for p in check_text_span(cid, line):
                        out.append(f"[{cid}] {rel}:{n}: {p} :: {line.strip()[:140]}")
    return out


def annotation_checks() -> list[str]:
    out = []
    for g in MENTION_GLOBS[:5]:
        for f in sorted(glob.glob(os.path.join(ROOT, g))):
            rel = os.path.relpath(f, ROOT).replace("\\", "/")
            t = _read(rel)
            for cid, span, start in gvc_spans(t):
                ln = t.count("\n", 0, start) + 1
                if cid not in R.CLAIMS:
                    out.append(f"[unknown id '{cid}'] {rel}:{ln}")
                    continue
                if span is None:
                    out.append(f"[{cid}] {rel}:{ln}: unterminated \\gvc annotation")
                    continue
                for p in check_text_span(cid, span):
                    out.append(f"[{cid}] {rel}:{ln}: annotated span: {p}")
            for n, line in enumerate(t.split("\n"), 1):
                for m in GVC_MD.finditer(line):
                    cid = m.group(1)
                    if cid not in R.CLAIMS:
                        out.append(f"[unknown id '{cid}'] {rel}:{n}")
                        continue
                    for p in check_text_span(cid, GVC_MD.sub("", line)):
                        out.append(f"[{cid}] {rel}:{n}: annotated row: {p}")
    return out


def mentions(cid: str) -> list[tuple[str, int, bool, str]]:
    c = R.CLAIMS[cid]
    res = []
    for g in MENTION_GLOBS:
        for f in sorted(glob.glob(os.path.join(ROOT, g))):
            rel = os.path.relpath(f, ROOT).replace("\\", "/")
            lines = _read(rel).split("\n")
            for i, l in enumerate(lines):
                if _c(c["mention"]).search(l):
                    win = "\n".join(lines[max(0, i - 2):i + 3])
                    res.append((rel, i + 1, bool(_c(c["qualifier"]).search(win)), l.strip()[:120]))
    return res


def selftest() -> int:
    ok = True
    # each claim: a span with its overclaim and without qualifier must FAIL; a span with the
    # qualifier and no overclaim must PASS.
    probes = {
        "pi-source-case-exhaustion": ("case-exhaustion theorem | INTERNAL THEOREM",
                                      "case-exhaustion theorem [CONDITIONAL], first fifteen projections"),
        "wh1-p38-gh-convergence": ("Paper 38 -- the discrete→continuum convergence theorem (unconditional)",
                                   "Paper 38 -- convergence of the true-window spectral truncations (commutative factor)"),
        "sm-gauge-group-upgrade-b": ("The inner algebra structure is mostly forced",
                                     "inner algebra, an organizing reading plus the adopted Upgrade B axiom"),
        "paper7-symbolic-checks": ("established with 18 independent proofs",
                                   "18 symbolic checks (ten limiting-case or definitional)"),
        "central-smoothing-cb-norm": ("the cb-norm = 2/(n_max+1) decays", "the cb-norm is 1 (UCP)"),
        "k-combination-rule": ("K = pi(B + F - Delta), a CONJECTURE of Paper 2",
                               "K = π(B + F − Δ): an Observation (not a conjecture or derivation)"),
        "p38-rate-constant": ("the convergence rate is (4/π) log n/n | INTERNAL THEOREM + MEASURED",
                              "the Fejér moment is (4/π) log n/n [INTERNAL THEOREM + MEASURED]; the distance is O(1/n) [MEASURED]"),
        "graph-s3-wall": ("the graph Laplacian converges to Δ on S^3",
                          "the graph Laplacian does not converge to Δ (WALL; disconnected l-blocks)"),
        "p7-parseval-f0": ("the density-overlap formula is a key structural finding of S^3",
                           "the density-overlap formula is the Fourier--Parseval form [PRIOR ART]"),
        "he-bracket": ("He 0.004% is a variational 0.004% bound",
                       "He 0.004% (a non-variational extrapolation, below exact)"),
        "p1-degree-not-monotone-l": ("higher-l nodes have more connections",
                                     "degree is not monotone in l; only the m = 0 nodes order D_2p > D_2s"),
        "p38-uniform-bound": ("gamma_n <= 6 log n/n holds for all n",
                              "gamma_n <= 6 log n/n [MEASURED] on a finite range, not proved for all n"),
        "k-match-precision": ("K = π(B+F−Δ) matches 1/α to 8.8×10^-8",
                              "K matches 1/α to 4.8×10^-7; the cubic's root to 8.8×10^-8"),
        "p40-general-g-conditional": ("Paper 40 generalises it to all compact Lie groups, closes it positively",
                                      "Paper 40 generalises it to compact semisimple groups (conditionally on the per-group spin-window decomposition)"),
        "m2-ring-integer-s": ("the pure-Tate even-weight ring holds at every s",
                              "the pure-Tate even-weight ring holds at integer s (Seeley-DeWitt coefficients)"),
        "three-homes-negative": ("no spectral construction on S^3 simultaneously generates all three",
                                 "no single construction generating all three has been exhibited (twelve eliminated; an empirical record)"),
        "p1-ladder-exact": ("[PANEL-VERIFIED] the ladder operators reproduce the quantum numbers",
                            "[SYMBOLIC PROOF] the ladder operators reproduce the quantum numbers on interior shells"),
        "m1-gh-rate-reading": ("the M1 mechanism produces the 4/π GH rate constant",
                               "the GH rate constant shares its value with M1 (a reading, not a derivation)"),
        "p55-census-scope": ("all of GeoVac's transcendentals sit in mixed-Tate periods: the engine classification",
                             "every transcendental in the census of S^3 spectral data sits in mixed-Tate periods"),
        "lipschitz-c3-compression": ("the joint Lipschitz comparison constant C_3^{(2)} \\le \\sqrt2 $ of the route",
                                     "the joint Lipschitz comparison constant C_3^{(2)} = 1, by compression"),
        "p39-state-space-not-propinquity": ("Paper~39 closes that gap in the tensor-product case",
                                            "Paper~39 proves the tensor-product case at the state-space GH level"),
        "p45-rates-spatial": ("Paper~45 propinquity proof gives the rate O(T/\\Nt)",
                              "Paper~45 proves a degeneracy theorem; its spatial rate is the surviving statement"),
        "p59-four-centre-genus-one": ("The four-centre ERI is proved to be genus one in every case.",
                                      "A fourth distinct centre adds no genus [MEASURED]: the four-centre ERI stays on a genus-one curve."),
        "p59-three-electron-k3": ("Longer trees give Calabi-Yau n-folds, as Paper 59 proves.",
                                  "Whether longer trees give Calabi-Yau n-folds is OPEN."),
        "p35-window-circle-adopted": ("The real-time window is a periodic Euclidean circle.",
                                      "Whether a real-time window forces the circle stays open; Paper 35 adopts the identification rather than derives it."),
        "seam-j4-j7-forced": ("J4 stays Obs until a guard exists.",
                              "J4 and J7 promoted from Obs 2026-10-15 (Paper 59 sec:beyond)."),
        "p59-coulomb-ring-elliptic": ("Rings never climb the period ladder.",
                                      "A Coulomb ring of three densities does not climb: its loop integrates out exactly."),
        "p59-correlated-triangle-k3": ("The correlated triangle also sits on the elliptic rung.",
                                       "The correlated triangle with exponential geminals sits with the chain on the K3 rung at gamma = 1."),
        "p59-ising-thresholds": ("In the Ising class the even moments drop a rung.",
                                 "In the Ising class t = 1 is a threshold exactly for even n up to 8, and parity is not the rule for values."),
        "p59-tree-banana-hadamard": ("The banana ladder relation holds for all $n$.",
                                     "The tree ladder is tied to the banana ladder by a Hadamard product, checked for n = 2 to 5."),
        "p59-slater-slice-membership": ("On the Slater slice, P G = 0 follows from the measured annihilation by W.",
                                        "On the Slater slice both annihilations are derived, with the Mellin contour shift as the one cap."),
    }
    missing = set(R.CLAIMS) - set(probes) - {"pi-source-case-exhaustion"} | (set(probes) - set(R.CLAIMS))
    if missing - {"pi-source-case-exhaustion"}:
        print(f"SELFTEST FAIL: claims without a fire+silent probe (or stale probes): {sorted(missing)}")
        ok = False
    # The exact strings that ESCAPED the first version (DELTA #24 code review): each must fire.
    escapes = [
        ("pi-source-case-exhaustion", "case-exhaustion theorem (unconditional): every pi engages M1/M2/M3"),
        ("wh1-p38-gh-convergence", "Paper 38 proves convergence for the almost-commutative extension (unconditional)"),
        ("m2-ring-integer-s", "pure-Tate even-weight ring holds at every half-integer s"),
        ("pi-source-case-exhaustion", "case-exhaustion (first fifteen) -- not only checked but INTERNAL THEOREM"),
        ("k-combination-rule", "the combination rule is an observation? no: K is derived from the spectral action"),
        ("pi-source-case-exhaustion", "case-exhaustion [CONDITIONAL] upgraded: INTERNAL THEOREM -> now proven"),
        ("pi-source-case-exhaustion", "the case-exhaustion is not conditional, INTERNAL THEOREM"),
        ("pi-source-case-exhaustion", "the case-exhaustion (first fifteen) is not conditional, INTERNAL THEOREM"),
        ("p40-general-g-conditional", "has been resolved by Paper~40 (unconditionally at rank~1; conditionally on the per-group spin-window decomposition)"),
        ("p40-general-g-conditional", "resolved by Paper 40, conditionally on the spin-window decomposition"),
        ("k-combination-rule", "the combination rule, an Observation, \\textsc{conjecture}"),
        ("k-combination-rule", "combination rule, once an Observation, is now a theorem"),
        ("k-combination-rule", "combination rule (an Observation) is proved from the spectral action"),
        ("wh1-p38-gh-convergence", "Paper 38 convergence (commutative factor), and the almost commutative extension is proved too"),
        ("wh1-p38-gh-convergence", "Paper 38 proves convergence for the non-commutative gauge extension"),
        ("m2-ring-integer-s", "pure-Tate even-weight ring at non-integer s"),
        ("central-smoothing-cb-norm", "the cb-norm of the smoothing, 2/(n_max+1), vanishes"),
        ("p38-uniform-bound", "6 log n/n [MEASURED], valid for every n_max (proved)"),
        ("k-match-precision", "K agrees with 1/alpha to 8.8e-8 (the cubic is irrelevant)"),
        # FULL #12 (2026-10-03): 4/pi presented as the convergence rate, with "Fejer" present elsewhere
        ("p38-rate-constant", "the Fejér kernel aside, the GH rate is (4/π) log n/n"),
        ("p38-rate-constant", "the rate at which the truncated triples converge (Fejér) carries the transcendental signature"),
        ("m1-gh-rate-reading", "the state-space GH rate sits in the M1 slot (a reading)"),
    ]
    for cid, bad in escapes:
        if not check_text_span(cid, bad):
            print(f"SELFTEST FAIL: {cid} escaped: {bad}")
            ok = False
    # history records that must stay silent
    for cid, good in [("k-combination-rule", "combination rule: an Observation (conjecture→observation downgrade 2026-06-14)"),
                      ("k-combination-rule", "the combination rule is an Observation, never a derivation, conjecture, prediction or theorem"),
                      ("k-combination-rule", "combination rule: an Observation; \\emph{not}~a conjecture"),
                      ("k-combination-rule", "an Observation; the Mellin slot k is a theorem-grade index elsewhere"),
                      ("wh1-p38-gh-convergence", "Paper 38: commutative-factor convergence leg proved as of 2026-06-10, the almost-commutative structure's convergence open"),
                      ("k-combination-rule", "K = π(B+F−Δ) is an Observation rather than a theorem")]:
        if check_text_span(cid, good):
            print(f"SELFTEST FAIL: {cid} fired on corrected text: {good} -> {check_text_span(cid, good)}")
            ok = False
    for cid, (bad, good) in probes.items():
        if not check_text_span(cid, bad):
            print(f"SELFTEST FAIL: {cid} did not fire on: {bad}")
            ok = False
        if check_text_span(cid, good):
            print(f"SELFTEST FAIL: {cid} fired on corrected text: {good} -> {check_text_span(cid, good)}")
            ok = False
    # annotation parser: nested braces are parsed, not silently skipped
    spans = list(gvc_spans(r"x \gvc{sm-gauge-group-upgrade-b}{the \emph{mostly \textbf{forced}} inner algebra} y"))
    if not (spans and spans[0][1] and "forced" in spans[0][1] and check_text_span(spans[0][0], spans[0][1])):
        print("SELFTEST FAIL: \\gvc nested-brace parser"); ok = False
    if list(gvc_spans(r"\gvc{k-combination-rule}{unterminated"))[0][1] is not None:
        print("SELFTEST FAIL: unterminated \\gvc not reported"); ok = False
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mentions")
    ap.add_argument("--changed-since")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.mentions:
        if a.mentions not in R.CLAIMS:
            print("unknown claim id"); return 1
        rows = mentions(a.mentions)
        print(f"# {a.mentions}: {R.CLAIMS[a.mentions]['statement']}")
        for rel, ln, q, txt in rows:
            print(f"{'  ' if q else '!!'} {rel}:{ln}: {txt}")
        print(f"# {len(rows)} mention(s); '!!' = no qualifier within +-2 lines (read these first)")
        return 0
    if a.changed_since:
        for cid, c in R.CLAIMS.items():
            if c["changed"] >= a.changed_since:
                print(f"{cid}  (changed {c['changed']})  owner {c['owner']}")
        return 0
    problems = owner_checks() + front_door_checks() + annotation_checks()
    for p in problems:
        print(p)
    n = len(R.CLAIMS)
    if problems:
        print(f"RESULT: FAIL ({len(problems)} violation(s) over {n} registered claims; front door = {', '.join(R.FRONT_DOOR)})")
        return 1
    print(f"RESULT: PASS ({n} registered claims; owners anchored; front door + annotated loci clean)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
