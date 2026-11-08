"""Build the Zenodo per-paper upload manifest.

Walks papers/group{1..6}_*/ and papers/synthesis/, extracts \\title and the
abstract environment from each .tex, pairs it with its tracked PDF, and emits
debug/data/zenodo_manifest.json — one entry per record for zenodo_upload.py.

Titles/abstracts are lightly de-TeXed for the Zenodo description field; inline
math is left as $...$ (readable, and honest about content). Entries that fail
extraction are flagged, not dropped — fix them by hand in the .tex or accept
the flagged manifest entry after review.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path(__file__).resolve().parents[1]
PAPER_DIRS = sorted(ROOT.glob("papers/group*_*")) + [ROOT / "papers" / "synthesis"]
OUT_PATH = ROOT / "debug" / "data" / "zenodo_manifest.json"

BASE_KEYWORDS = [
    "GeoVac", "spectral graph theory", "graph Laplacian", "quantum chemistry",
    "Fock projection", "hyperspherical harmonics",
]
GROUP_KEYWORDS: Dict[str, List[str]] = {
    "group1_operator_algebras": [
        "noncommutative geometry", "spectral triple", "operator systems",
        "Gromov-Hausdorff propinquity", "spectral truncation",
    ],
    "group2_quantum_chemistry": [
        "electronic structure", "Sturmian basis", "molecular Hamiltonians",
        "pseudopotential", "configuration interaction",
    ],
    "group3_foundations": [
        "mathematical physics", "periods", "mixed Tate motives",
        "transcendental number theory", "exchange constants",
    ],
    "group4_quantum_computing": [
        "quantum simulation", "qubit Hamiltonian", "VQE", "Pauli decomposition",
        "quantum resource estimation",
    ],
    "group5_qed_gauge": [
        "quantum electrodynamics", "gauge theory", "spectral action",
        "fine structure constant", "lattice QED",
    ],
    "group6_precision_observations": [
        "precision spectroscopy", "hyperfine structure", "Lamb shift",
        "atomic physics", "QED tests",
    ],
    "synthesis": [
        "research synthesis", "noncommutative geometry", "quantum simulation",
    ],
}


def strip_comments(tex: str) -> str:
    """Remove LaTeX comments (% to EOL, keeping escaped \\%)."""
    return re.sub(r"(?<!\\)%.*", "", tex)


def find_braced_arg(tex: str, command: str) -> Optional[str]:
    """Return the brace-balanced argument of the first \\command{...}."""
    m = re.search(r"\\" + command + r"\s*\{", tex)
    if not m:
        return None
    depth, start = 1, m.end()
    for i in range(start, len(tex)):
        if tex[i] == "{" and tex[i - 1] != "\\":
            depth += 1
        elif tex[i] == "}" and tex[i - 1] != "\\":
            depth -= 1
            if depth == 0:
                return tex[start:i]
    return None


# TeX accent → unicode (the forms that occur in this corpus)
ACCENTS = [
    (r'\"o', "ö"), (r'\"a', "ä"), (r'\"u', "ü"), (r'\"O', "Ö"),
    (r"\'e", "é"), (r"\`e", "è"), (r"\'a", "á"), (r"\'o", "ó"),
]
# project/custom macros → readable unicode (word-boundary matched)
MACROS = [  # NB: (?![a-zA-Z]) not \b — macros may abut _, digits, unicode.
    # Symbol macros FIRST: replacing a letter-producing macro (\Uone → U(1))
    # before an adjacent symbol macro (\times) would poison its lookahead.
    (r"\\rtimes(?![a-zA-Z])", "⋊"), (r"\\otimes(?![a-zA-Z])", "⊗"),
    (r"\\times(?![a-zA-Z])", "×"), (r"\\pi(?![a-zA-Z])", "π"),
    (r"\\alpha(?![a-zA-Z])", "α"), (r"\\zeta(?![a-zA-Z])", "ζ"),
    (r"\\sthree(?![a-zA-Z])", "S³"), (r"\\stwo(?![a-zA-Z])", "S²"),
    (r"\\sfive(?![a-zA-Z])", "S⁵"),
    (r"\\nmax(?![a-zA-Z])", "n_max"), (r"\\lmax(?![a-zA-Z])", "l_max"),
    (r"\\SU(?![a-zA-Z])", "SU"), (r"\\Uone(?![a-zA-Z])", "U(1)"),
    (r"\\SL(?![a-zA-Z])", "SL"), (r"\\Ga(?![a-zA-Z])", "G_a"),
    (r"\\fibfun(?![a-zA-Z])", "ω"), (r"\\Aut(?![a-zA-Z])", "Aut"),
    (r"\\mathbb\{R\}", "ℝ"), (r"\\R(?![a-zA-Z])", "ℝ"),
]


def detex(s: str, title_mode: bool = False) -> str:
    """Lightly flatten LaTeX prose for Zenodo title/description fields.

    title_mode additionally strips math delimiters and sub/superscript
    braces so record titles are plain searchable text; abstracts keep
    their $...$ inline math (readable, and Zenodo renders MathJax).
    """
    s = re.sub(r"\\texorpdfstring\{([^{}]*)\}\{[^{}]*\}", r"\1", s)
    s = re.sub(r"\\thanks\{[^{}]*\}", "", s)
    # unwrap one level of styling macros, repeatedly for nesting
    for _ in range(4):
        s = re.sub(
            r"\\(?:emph|textit|textbf|texttt|textrm|textsc|mbox|text)\{([^{}]*)\}",
            r"\1", s)
    s = re.sub(r"\\(?:cite[pt]?|citealp)\*?(?:\[[^\]]*\])?\{[^{}]*\}", "[ref]", s)
    s = re.sub(r"\\(?:ref|eqref|autoref|cref|Cref)\{[^{}]*\}", "(ref)", s)
    s = re.sub(r"\\footnote\{[^{}]*\}", "", s)
    # flatten list/theorem-like environments that occur inside abstracts
    s = re.sub(r"\\begin\{(?:itemize|enumerate|prediction)\*?\}", " ", s)
    s = re.sub(r"\\end\{(?:itemize|enumerate|prediction)\*?\}", " ", s)
    s = re.sub(r"\\item\b", " • ", s)
    s = re.sub(r"\\(?:begin|end)\{equation\*?\}", " $$ ", s)  # display math
    s = re.sub(r"\\\\(\[[^\]]*\])?", " ", s)  # forced line breaks → space
    s = re.sub(r"\\(?:large|Large|LARGE|small|footnotesize|normalsize)\b", "", s)
    for tex, uni in ACCENTS:
        s = s.replace(tex + "{", uni + "{").replace(tex, uni)
        s = s.replace(uni + "{" + "}", uni)  # collapse empty brace remnant
    for pat, uni in MACROS:
        s = re.sub(pat, uni, s)
    s = s.replace(r"\%", "%").replace(r"\&", "&").replace(r"\_", "_")
    s = s.replace(r"\,", " ").replace("~", " ").replace(r"\ ", " ")
    s = re.sub(r"(?<!\\)\\(?:noindent|par|smallskip|medskip|bigskip)\b", " ", s)
    s = s.replace("---", "—").replace("--", "–")
    if title_mode:
        s = s.replace("$", "")
        # _{X}/^{X} → _X / ^(long) so titles stay plain text
        s = re.sub(r"([_^])\{(\w)\}", r"\1\2", s)
        s = re.sub(r"([_^])\{([^{}]+)\}", r"\1(\2)", s)
        s = s.replace("{", "").replace("}", "")
    s = re.sub(r"\s+", " ", s)
    return s.strip(" :")


def extract(tex_path: Path) -> Dict[str, object]:
    tex = strip_comments(tex_path.read_text(encoding="utf-8", errors="replace"))
    title = find_braced_arg(tex, "title")
    m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.DOTALL)
    abstract = m.group(1) if m else None
    return {
        "title": detex(title, title_mode=True) if title else None,
        "abstract": detex(abstract) if abstract else None,
    }


def paper_number(stem: str) -> Optional[str]:
    m = re.match(r"paper_(\w+?)_", stem)
    return m.group(1) if m else None


def main() -> None:
    entries: List[Dict[str, object]] = []
    for d in PAPER_DIRS:
        group = d.name
        for tex_path in sorted(d.glob("*.tex")):
            stem = tex_path.stem
            pdf_path = tex_path.with_suffix(".pdf")
            meta = extract(tex_path)
            flags = []
            if not meta["title"]:
                flags.append("NO_TITLE")
            if not meta["abstract"]:
                flags.append("NO_ABSTRACT")
            if not pdf_path.exists():
                flags.append("NO_PDF")
            entries.append({
                "id": stem,
                "group": group,
                "paper_no": paper_number(stem),
                "title": meta["title"],
                "abstract": meta["abstract"],
                "pdf": str(pdf_path.relative_to(ROOT)).replace("\\", "/"),
                "keywords": BASE_KEYWORDS + GROUP_KEYWORDS.get(group, []),
                "flags": flags,
            })
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(entries, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    flagged = [e for e in entries if e["flags"]]
    print(f"{len(entries)} entries -> {OUT_PATH.relative_to(ROOT)}")
    print(f"{len(flagged)} flagged:")
    for e in flagged:
        print(f"  {e['id']}: {','.join(e['flags'])}")


if __name__ == "__main__":
    main()
