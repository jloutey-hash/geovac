"""Generate static per-paper abstract pages for the GitHub Pages site.

Consumes debug/data/zenodo_manifest.json (built by zenodo_manifest_build.py)
and writes crawlable static HTML into viz/public/ (Vite copies public/
verbatim into dist/, so these deploy through the existing deploy-viz.yml
workflow on the next push):

    viz/public/papers/<id>.html   one page per paper (title, abstract, links)
    viz/public/papers/index.html  grouped listing
    viz/public/sitemap.xml        absolute URLs for crawlers
    viz/public/robots.txt         allow-all + sitemap pointer

If debug/data/zenodo_upload_results.json exists (written by zenodo_upload.py
--publish), each paper's minted DOI is linked on its page — re-run this script
after publishing to enrich the pages.

Rationale (docs/project_closeout_plan.md §B2): the viz app is a React SPA and
serves an empty <div> to crawlers; these static pages are the surface an AI
or search engine actually reads.
"""
from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "debug" / "data" / "zenodo_manifest.json"
RESULTS = ROOT / "debug" / "data" / "zenodo_upload_results.json"
OUT_DIR = ROOT / "viz" / "public" / "papers"
PUBLIC_DIR = ROOT / "viz" / "public"

SITE_BASE = "https://jloutey-hash.github.io/geovac"
REPO_URL = "https://github.com/jloutey-hash/geovac"
CORPUS_DOI = "10.5281/zenodo.20482394"
AUTHOR = "Josh Loutey"

GROUP_LABEL = {
    "group1_operator_algebras": "Operator Algebras / Noncommutative Geometry",
    "group2_quantum_chemistry": "Quantum Chemistry",
    "group3_foundations": "Mathematical Foundations",
    "group4_quantum_computing": "Quantum Computing",
    "group5_qed_gauge": "QED and Gauge Theory",
    "group6_precision_observations": "Precision Observations",
    "synthesis": "Synthesis / Reading Guides",
}
GROUP_ORDER = list(GROUP_LABEL)

CSS = """
:root { color-scheme: light dark; }
body { max-width: 46rem; margin: 2rem auto; padding: 0 1rem;
       font: 1rem/1.6 system-ui, sans-serif; }
a { color: #2a78d6; }
h1 { line-height: 1.3; font-size: 1.5rem; }
.meta { color: #777; font-size: 0.9rem; }
.abstract { text-align: justify; }
nav { margin-bottom: 1.5rem; font-size: 0.9rem; }
footer { margin-top: 2.5rem; border-top: 1px solid #7774;
         padding-top: 1rem; font-size: 0.85rem; color: #777; }
li { margin: 0.4rem 0; }
"""

MATHJAX = ('<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/'
           'tex-chtml.js"></script>')


def page(title: str, description: str, body: str,
         with_math: bool = False) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<meta name="description" content="{html.escape(description)}">
<title>{html.escape(title)}</title>
<style>{CSS}</style>
{MATHJAX if with_math else ''}
</head>
<body>
{body}
<footer>
<p>GeoVac — an independent research project. Code, tests, and all papers:
<a href="{REPO_URL}">{REPO_URL}</a> · Corpus archive DOI:
<a href="https://doi.org/{CORPUS_DOI}">{CORPUS_DOI}</a> ·
<a href="../index.html">Interactive visualizations</a></p>
</footer>
</body>
</html>
"""


def paper_page(e: Dict[str, Any], doi: Optional[str]) -> str:
    paper_no = e.get("paper_no")
    label = f"GeoVac Paper {paper_no}" if paper_no else "GeoVac synthesis document"
    group = GROUP_LABEL.get(e["group"], e["group"])
    pdf_url = f"{REPO_URL}/blob/main/{e['pdf']}"
    doi_line = (f' · DOI: <a href="https://doi.org/{doi}">{doi}</a>'
                if doi else "")
    desc = e["abstract"][:300].rsplit(" ", 1)[0] + "…"
    body = f"""<nav><a href="index.html">← All GeoVac papers</a></nav>
<h1>{html.escape(e['title'])}</h1>
<p class="meta">{AUTHOR} · {html.escape(label)} · {html.escape(group)}{doi_line}
 · <a href="{pdf_url}">PDF</a></p>
<h2>Abstract</h2>
<p class="abstract">{html.escape(e['abstract'])}</p>
"""
    return page(f"{e['title']} — GeoVac", desc, body, with_math=True)


def index_page(entries: List[Dict[str, Any]]) -> str:
    sections = []
    for g in GROUP_ORDER:
        items = [e for e in entries if e["group"] == g]
        if not items:
            continue
        lis = "\n".join(
            f'<li><a href="{e["id"]}.html">{html.escape(e["title"])}</a></li>'
            for e in items)
        sections.append(f"<h2>{html.escape(GROUP_LABEL[g])}</h2>\n<ul>{lis}</ul>")
    body = (
        "<h1>The GeoVac Paper Series</h1>\n"
        "<p>GeoVac is a spectral-graph-theory framework for computational "
        "quantum chemistry and mathematical physics: the discrete graph "
        "Laplacian on S³ is conformally equivalent to the Schrödinger "
        "equation (Fock 1935), and the framework exploits this to build "
        "structurally sparse qubit Hamiltonians, a discrete almost-"
        "commutative spectral triple with proven Gromov–Hausdorff "
        "convergence, and a classification of where transcendental "
        "constants enter physical observables. Abstracts below; PDFs and "
        "code in the repository.</p>\n" + "\n".join(sections))
    return page("The GeoVac Paper Series — abstracts and PDFs",
                "Abstracts and PDFs for the ~60-paper GeoVac series: sparse "
                "qubit Hamiltonians from spectral graph theory, the discrete "
                "spectral triple, periods, QED and precision observables.",
                body)


def main() -> None:
    entries: List[Dict[str, Any]] = json.loads(
        MANIFEST.read_text(encoding="utf-8"))
    dois: Dict[str, str] = {}
    if RESULTS.exists():
        results = json.loads(RESULTS.read_text(encoding="utf-8"))
        dois = {k: v["doi"] for k, v in results.items() if v.get("doi")}

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for e in entries:
        (OUT_DIR / f"{e['id']}.html").write_text(
            paper_page(e, dois.get(e["id"])), encoding="utf-8")
    (OUT_DIR / "index.html").write_text(index_page(entries), encoding="utf-8")

    urls = [f"{SITE_BASE}/", f"{SITE_BASE}/papers/index.html"] + [
        f"{SITE_BASE}/papers/{e['id']}.html" for e in entries]
    sitemap = ('<?xml version="1.0" encoding="UTF-8"?>\n'
               '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
               + "\n".join(f"  <url><loc>{u}</loc></url>" for u in urls)
               + "\n</urlset>\n")
    (PUBLIC_DIR / "sitemap.xml").write_text(sitemap, encoding="utf-8")
    (PUBLIC_DIR / "robots.txt").write_text(
        f"User-agent: *\nAllow: /\nSitemap: {SITE_BASE}/sitemap.xml\n",
        encoding="utf-8")

    print(f"{len(entries)} paper pages + index -> {OUT_DIR.relative_to(ROOT)}")
    print(f"sitemap ({len(urls)} URLs) + robots.txt -> "
          f"{PUBLIC_DIR.relative_to(ROOT)}")
    print(f"DOI links: {len(dois)} (re-run after zenodo_upload.py --publish)")


if __name__ == "__main__":
    main()
