r"""Per-unit certification stamps: review what changed, not everything again.

WHY THIS EXISTS (2026-10-03, PI-directed after /qa trunk FULL #12)
-----------------------------------------------------------------
Every FULL run re-reviewed every dimension of every paper from scratch.  By
FULL #9-#12 the code dimension had stopped finding owner defects (every keystone
re-derived clean, run after run) and the citation dimension was converging the
same way, yet both were re-run in full each time.  The prose dimension had the
opposite problem: an 8,000-line paper is read by a different reviewer each run,
each one samples a different region, and each finds debt nobody had read
before -- so a FULL run measured "how much unread prose is left", not "did the
last remediation hold".

So a certification is recorded per UNIT, keyed to a content hash of that unit's
inputs, and a later run reviews only the units whose inputs changed, whose
claims moved, or that fall in a random audit sample of the still-current ones.

UNITS
-----
  prose:<paper>::<section / subsection title>
      The text of one \section or \subsection chunk (the part before the first
      \section -- title, abstract, front matter -- is the unit "front").
      Reviewed by the claims dimension (synthesis files: the synthesis
      dimension).
  cite:<paper>
      The bibliography block plus every line carrying a \cite.  Reviewed by the
      citation dimension.
  code:<paper>
      Every test the paper backs itself with (cited inline; test_paper{N}_*;
      matrix rows naming the paper) plus every geovac module those tests import.
      Reviewed by the code dimension.

Hashes are taken on text normalised to LF line endings with trailing
whitespace stripped, so a line-ending or whitespace-only rewrite does not
invalidate a stamp.

A UNIT NEEDS REVIEW WHEN
------------------------
  NEW          it has no stamp;
  CHANGED      its hash differs from the stamp;
  CLAIM-MOVED  its hash is unchanged but a C25 claim it mentions (by the claim's
               `mention` pattern) changed status after the stamp, or a C16 entry
               dated after the stamp names its paper as a `cited_by` dependent --
               the text did not move, the claim under it did;
  AUDIT        a deterministic random sample of the CURRENT units (10%, min 2 per
               dimension).  The audit measures stamp decay: a LARGE found in an
               audited unit means stamps are leaking and the next run widens.

USAGE
-----
  python debug/qa/cert_stamps.py --target trunk --status
  python debug/qa/cert_stamps.py --target trunk --plan [--seed 2026-10-03]
  python debug/qa/cert_stamps.py --target trunk --stamp --run "FULL #13" \
         [--dims prose,cite,code] [--exclude SUBSTR ...]
  python debug/qa/cert_stamps.py --selftest

`--stamp` records the CURRENT hashes as reviewed-clean.  Run it only after the
run's findings in those units are remediated AND re-reviewed (the self-check
loop) -- remediated text is not clean text.  `--exclude` keeps units with open
findings unstamped (substring match on the unit id).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import importlib.util
import json
import random
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STAMPS = ROOT / "docs" / "qa" / "stamps"
TESTS = ROOT / "tests"
MATRIX = ROOT / "docs" / "claim_test_matrix.md"
DATE = re.compile(r"20\d{2}-\d{2}-\d{2}")


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def target_papers(target: str) -> list[str]:
    """Repo-relative .tex files a target covers (same scope as the staleness audit)."""
    staleness = _load("_stale", "debug/qa/check_cert_staleness.py")
    out: list[str] = []
    for p in staleness.paths_for(target):
        q = ROOT / p
        if q.is_dir():
            out += [str(f.relative_to(ROOT)).replace("\\", "/") for f in sorted(q.glob("*.tex"))]
        elif q.suffix == ".tex":
            out.append(p)
    return sorted(set(out))


# --------------------------------------------------------------------------- units
def _norm(text: str) -> str:
    return "\n".join(line.rstrip() for line in text.replace("\r\n", "\n").split("\n"))


def _h(text: str) -> str:
    return hashlib.sha256(_norm(text).encode("utf-8")).hexdigest()[:16]


SECTION = re.compile(r"^\s*\\(section|subsection)\*?(?:\[[^\]]*\])?\{(.*)$")


def _title(raw: str) -> str:
    depth, out = 1, []
    for ch in raw:
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                break
        out.append(ch)
    t = re.sub(r"\\[A-Za-z]+\*?|[{}$~\\]", " ", "".join(out))
    return re.sub(r"\s+", " ", t).strip()[:70] or "(untitled)"


def prose_units(rel: str, text: str) -> dict[str, str]:
    """{unit_id: chunk text}; chunks split at \\section and \\subsection."""
    units: dict[str, str] = {}
    seen: dict[str, int] = {}
    cur_sec, cur_id, buf = "", f"prose:{rel}::front", []

    def flush():
        if buf:
            units[cur_id] = units.get(cur_id, "") + "\n".join(buf)

    for line in text.replace("\r\n", "\n").split("\n"):
        m = SECTION.match(line)
        if m:
            flush()
            buf = []
            t = _title(m.group(2))
            if m.group(1) == "section":
                cur_sec, name = t, t
            else:
                name = f"{cur_sec} / {t}" if cur_sec else t
            seen[name] = seen.get(name, 0) + 1
            if seen[name] > 1:
                name = f"{name} #{seen[name]}"
            cur_id = f"prose:{rel}::{name}"
        buf.append(line)
    flush()
    return units


BIB = re.compile(r"\\begin\{thebibliography\}.*?\\end\{thebibliography\}", re.S)


def cite_text(text: str) -> str:
    bib = "\n".join(BIB.findall(text))
    body = BIB.sub("", text)
    cites = [ln for ln in body.replace("\r\n", "\n").split("\n") if "\\cite" in ln]
    return bib + "\n----\n" + "\n".join(cites)


TEST_REF = re.compile(r"tests/test(?:\\?_[A-Za-z0-9*]+)+(?:\.py)?|test(?:\\?_[A-Za-z0-9*]+)+\.py")
GEOVAC_IMPORT = re.compile(r"^\s*(?:from\s+(geovac(?:\.[A-Za-z0-9_]+)*)\s+import|import\s+(geovac(?:\.[A-Za-z0-9_]+)*))", re.M)


def _paper_number(rel: str) -> str | None:
    m = re.search(r"[Pp]aper_(\d+)_", rel)
    return m.group(1) if m else None


def _matrix_tests(n: str) -> set[str]:
    out: set[str] = set()
    if not MATRIX.exists():
        return out
    for line in MATRIX.read_text(encoding="utf-8", errors="replace").splitlines():
        cells = [c.strip() for c in line.split("|")]
        if len(cells) < 4:
            continue
        papers = re.findall(r"\d+", cells[1])
        if n in papers:
            out |= set(re.findall(r"\btest_[A-Za-z0-9_]+(?=\.py|\b)", line))
    return out


def code_inputs(rel: str, text: str) -> list[Path]:
    """Test files backing the paper + the geovac modules they import."""
    names: set[str] = set()
    for ref in TEST_REF.findall(text):
        s = ref.replace("\\_", "_").replace("\\", "")
        s = s.split("/", 1)[1] if "/" in s else s
        names.add(s[:-3] if s.endswith(".py") else s)
    n = _paper_number(rel)
    if n:
        names |= _matrix_tests(n)
    files: set[Path] = set()
    for name in names:
        files |= {p for p in TESTS.glob(f"{name}.py")}
    if n:
        files |= set(TESTS.glob(f"test_[Pp]aper{n}_*.py")) | set(TESTS.glob(f"test_[Pp]aper_{n}_*.py"))
    mods: set[Path] = set()
    for f in files:
        for a, b in GEOVAC_IMPORT.findall(f.read_text(encoding="utf-8", errors="replace")):
            dotted = a or b
            p = ROOT / Path(*dotted.split("."))
            for cand in (p.with_suffix(".py"), p / "__init__.py"):
                if cand.exists():
                    mods.add(cand)
                    break
    return sorted(files | mods)


def units_for(target: str) -> dict[str, dict]:
    """{unit_id: {"hash", "dim", "paper", "text"}} for every unit of the target."""
    out: dict[str, dict] = {}
    for rel in target_papers(target):
        text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
        pdim = "synthesis" if "/synthesis/" in rel else "prose"
        # the bibliography belongs to the cite unit; leaving it in the last prose
        # chunk made every bibitem edit invalidate "Acknowledgments"
        for uid, chunk in prose_units(rel, BIB.sub("", text)).items():
            out[uid] = dict(hash=_h(chunk), dim=pdim, paper=rel, text=chunk)
        ct = cite_text(text)
        out[f"cite:{rel}"] = dict(hash=_h(ct), dim="cite", paper=rel, text=ct)
        files = code_inputs(rel, text)
        if files:
            blob = "\n".join(f"== {f.relative_to(ROOT).as_posix()}\n" + f.read_text(encoding="utf-8", errors="replace")
                             for f in files)
            out[f"code:{rel}"] = dict(hash=_h(blob), dim="code", paper=rel, text="",
                                      inputs=[f.relative_to(ROOT).as_posix() for f in files])
    return out


# --------------------------------------------------------------------------- stamps
def load_stamps(target: str) -> dict:
    f = STAMPS / f"{target}.json"
    if f.exists():
        return json.loads(f.read_text(encoding="utf-8"))
    return {"target": target, "units": {}}


def save_stamps(target: str, data: dict) -> None:
    STAMPS.mkdir(parents=True, exist_ok=True)
    (STAMPS / f"{target}.json").write_text(json.dumps(data, indent=1, sort_keys=True) + "\n", encoding="utf-8")


def _claim_moves() -> list[tuple[str, str, re.Pattern]]:
    """(claim id, changed date, mention pattern) for every C25 claim."""
    try:
        reg = _load("_c25", "debug/qa/claim_registry.py")
    except Exception:
        return []
    out = []
    for cid, c in getattr(reg, "CLAIMS", {}).items():
        if c.get("mention") and c.get("changed"):
            out.append((cid, c["changed"], re.compile(c["mention"], re.I)))
    return out


def _c16_dependents() -> list[tuple[str, str, set[str]]]:
    """(entry id, date, cited_by paper set) for every C16 entry with dependents."""
    try:
        reg = _load("_c16", "debug/qa/check_retracted_terms.py")
    except Exception:
        return []
    out = []
    for e in getattr(reg, "REGISTRY", []):
        deps = set((e.get("cited_by") or {}).keys())
        dates = DATE.findall(str(e.get("note", "")) + " " + str(e.get("id", "")))
        if deps and dates:
            out.append((str(e.get("id")), max(dates), {d.replace("\\", "/") for d in deps}))
    return out


def classify(units: dict, stamps: dict, claims=None, c16=None) -> dict[str, tuple[str, str]]:
    """{unit_id: (state, reason)} with state in NEW / CHANGED / CLAIM-MOVED / CURRENT."""
    claims = _claim_moves() if claims is None else claims
    c16 = _c16_dependents() if c16 is None else c16
    out = {}
    for uid, u in units.items():
        st = stamps["units"].get(uid)
        if not st:
            out[uid] = ("NEW", "")
            continue
        if st["hash"] != u["hash"]:
            out[uid] = ("CHANGED", "")
            continue
        moved = [cid for cid, d, pat in claims if d > st["date"] and u["text"] and pat.search(u["text"])]
        if u["dim"] in ("prose", "synthesis"):
            moved += [eid for eid, d, deps in c16 if d > st["date"] and u["paper"] in deps
                      and uid.endswith("::front")]
        out[uid] = ("CLAIM-MOVED", ", ".join(moved)) if moved else ("CURRENT", "")
    return out


def audit_sample(states: dict, units: dict, seed: str, frac: float = 0.10, floor: int = 2) -> list[str]:
    rng = random.Random(seed)
    picks: list[str] = []
    by_dim: dict[str, list[str]] = {}
    for uid, (s, _) in states.items():
        if s == "CURRENT":
            by_dim.setdefault(units[uid]["dim"], []).append(uid)
    for dim in sorted(by_dim):
        pool = sorted(by_dim[dim])
        k = min(len(pool), max(floor, round(frac * len(pool))))
        picks += rng.sample(pool, k)
    return sorted(picks)


# --------------------------------------------------------------------------- CLI
def _lines(units: dict, uid: str) -> int:
    return units[uid]["text"].count("\n") + 1 if units[uid]["text"] else 0


def cmd_status(target: str, plan: bool, seed: str) -> int:
    units = units_for(target)
    stamps = load_stamps(target)
    states = classify(units, stamps)
    review = {u for u, (s, _) in states.items() if s != "CURRENT"}
    audit = audit_sample(states, units, seed) if plan else []
    counts: dict[str, dict[str, int]] = {}
    for uid, (s, _) in states.items():
        counts.setdefault(units[uid]["dim"], {}).setdefault(s, 0)
        counts[units[uid]["dim"]][s] += 1
    print(f"target {target}: {len(units)} units, {len(stamps['units'])} stamped")
    for dim in sorted(counts):
        print(f"  {dim:<10} " + "  ".join(f"{k}={v}" for k, v in sorted(counts[dim].items())))
    if plan:
        tot = sum(_lines(units, u) for u in units if units[u]["dim"] in ("prose", "synthesis"))
        rev = sum(_lines(units, u) for u in review | set(audit) if units[u]["dim"] in ("prose", "synthesis"))
        print(f"\nREVIEW PLAN (seed {seed}): prose lines to read {rev} of {tot}")
        for uid in sorted(review | set(audit)):
            s, why = states[uid]
            tag = "AUDIT" if uid in audit else s
            extra = f"  [{why}]" if why else ""
            n = f" ({_lines(units, uid)} lines)" if units[uid]["text"] else ""
            print(f"  {tag:<12} {uid}{n}{extra}")
    else:
        for uid in sorted(review):
            s, why = states[uid]
            print(f"  {s:<12} {uid}" + (f"  [{why}]" if why else ""))
    return 0


def _units_at(target: str, rev: str) -> dict[str, str]:
    """{unit_id: hash} for the prose and cite units as they stood at git revision `rev`."""
    import subprocess
    out: dict[str, str] = {}
    for rel in target_papers(target):
        text = subprocess.run(["git", "show", f"{rev}:{rel}"], capture_output=True, text=True,
                              encoding="utf-8", errors="replace", cwd=ROOT).stdout
        if not text:
            continue
        for uid, chunk in prose_units(rel, BIB.sub("", text)).items():
            out[uid] = _h(chunk)
        out[f"cite:{rel}"] = _h(cite_text(text))
    return out


def cmd_stamp(target: str, run: str, dims: set[str], exclude: list[str],
              unchanged_since: str | None = None) -> int:
    """Stamp current hashes.  With `unchanged_since`, stamp only units whose text is identical
    at that revision (the review happened there; remediated units wait for the next run).
    Code units are never stamped this way -- their inputs are not reconstructed at `rev`."""
    units = units_for(target)
    data = load_stamps(target)
    today = _dt.date.today().isoformat()
    then = _units_at(target, unchanged_since) if unchanged_since else None
    n = 0
    for uid, u in units.items():
        if u["dim"] not in dims or any(x in uid for x in exclude):
            continue
        if then is not None and then.get(uid) != u["hash"]:
            continue
        data["units"][uid] = dict(hash=u["hash"], run=run, date=today)
        n += 1
    live = set(units)
    stale = [k for k in data["units"] if k not in live]
    for k in stale:
        del data["units"][k]
    save_stamps(target, data)
    print(f"stamped {n} units ({', '.join(sorted(dims))}) for {target} as '{run}' on {today}; "
          f"dropped {len(stale)} stamps for units that no longer exist")
    return 0


def _selftest() -> int:
    ok = True

    def check(name, cond):
        nonlocal ok
        print(("[ok]   " if cond else "[FAIL] ") + name)
        ok &= bool(cond)

    tex = ("\\title{T}\n\\begin{abstract}A.\\end{abstract}\n\\section{Intro}\nx \\cite{a}\n"
           "\\subsection{Sub}\ny\n\\section{Intro}\nz\n"
           "\\begin{thebibliography}{9}\\bibitem{a} A.\\end{thebibliography}\n")
    pu = prose_units("p.tex", tex)
    check("units split at section and subsection, duplicates numbered",
          set(pu) == {"prose:p.tex::front", "prose:p.tex::Intro", "prose:p.tex::Intro / Sub",
                      "prose:p.tex::Intro #2"})
    check("CRLF / trailing-space rewrite keeps the hash", _h("a \r\nb") == _h("a\nb"))
    check("a one-word edit changes the hash", _h("a\nb") != _h("a\nc"))
    units = {uid: dict(hash=_h(t), dim="prose", paper="p.tex", text=t) for uid, t in pu.items()}
    stamps = {"units": {uid: dict(hash=u["hash"], date="2026-10-01", run="r") for uid, u in units.items()}}
    st = classify(units, stamps, claims=[], c16=[])
    check("unchanged + no claim moves -> all CURRENT", all(s == "CURRENT" for s, _ in st.values()))
    edited = dict(units)
    edited["prose:p.tex::Intro / Sub"] = dict(units["prose:p.tex::Intro / Sub"], hash=_h("y2"))
    st = classify(edited, stamps, claims=[], c16=[])
    check("an edited unit is CHANGED and only it",
          [u for u, (s, _) in st.items() if s != "CURRENT"] == ["prose:p.tex::Intro / Sub"])
    claims = [("c1", "2026-10-02", re.compile("z")), ("c0", "2026-09-01", re.compile("x"))]
    st = classify(units, stamps, claims=claims, c16=[])
    check("a claim moved after the stamp flags exactly the units that mention it",
          {u for u, (s, _) in st.items() if s == "CLAIM-MOVED"} == {"prose:p.tex::Intro #2"})
    c16 = [("e1", "2026-10-02", {"p.tex"})]
    st = classify(units, stamps, claims=[], c16=c16)
    check("a later C16 entry citing the paper flags its front matter",
          st["prose:p.tex::front"][0] == "CLAIM-MOVED")
    st = classify(units, {"units": {}}, claims=[], c16=[])
    check("no stamps -> all NEW", all(s == "NEW" for s, _ in st.values()))
    big = {f"prose:p::{i}": dict(dim="prose") for i in range(50)}
    sts = {k: ("CURRENT", "") for k in big}
    a1, a2 = audit_sample(sts, big, "s"), audit_sample(sts, big, "s")
    check("audit sample is deterministic and 10%", a1 == a2 and len(a1) == 5)
    check("audit sample floor is 2", len(audit_sample({"x": ("CURRENT", "")}, {"x": dict(dim="cite")}, "s")) == 1
          and len(audit_sample({k: sts[k] for k in list(sts)[:5]}, big, "s")) == 2)
    ct = cite_text(tex)
    check("cite unit holds the bibliography and the citing lines", "bibitem{a}" in ct and "x \\cite{a}" in ct)
    print("\nSELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--seed", default=_dt.date.today().isoformat())
    ap.add_argument("--stamp", action="store_true")
    ap.add_argument("--run")
    ap.add_argument("--dims", default="prose,synthesis,cite,code")
    ap.add_argument("--exclude", action="append", default=[])
    ap.add_argument("--unchanged-since", help="stamp only units whose text is identical at this git revision")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return _selftest()
    if not a.target:
        ap.error("--target required")
    if a.stamp:
        if not a.run:
            ap.error("--stamp needs --run")
        return cmd_stamp(a.target, a.run, set(a.dims.split(",")), a.exclude, a.unchanged_since)
    return cmd_status(a.target, a.plan, a.seed)


if __name__ == "__main__":
    sys.exit(main())
