"""Build the site's Downloads page (viz/public/downloads/) and the library LiH files.

1. Exports the library's composed LiH Hamiltonian, `geovac.hamiltonian('LiH')` at its defaults (the same object the
   site's QC dashboard shows), as an FCIDUMP, a gzipped Pauli-term list and a metadata JSON.
2. Writes viz/public/downloads/index.html from every *.json metadata file in that directory. No physics numeral is
   typed into the page: each one is read from a metadata JSON written by a generator that checks it
   (this script for LiH; benchmarks/export_paper62_water_hamiltonian.py for the Paper 62 water basis).

Run the water exporter first, then this:
    python benchmarks/export_paper62_water_hamiltonian.py
    python benchmarks/build_site_downloads.py
"""
from __future__ import annotations

import gzip
import html
import json
import os
import sys
from typing import Any, Dict, List

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "viz", "public", "downloads")
SITE_BASE = "https://jloutey-hash.github.io/geovac"
REPO_URL = "https://github.com/jloutey-hash/geovac"

# registry lih_composed_pauli (838 identity included), lih_composed_lambda (identity excluded)
LIH_PAULI_WITH_IDENTITY = 838
LIH_LAMBDA = 27.715


def export_lih() -> None:
    from geovac.ecosystem_export import hamiltonian, read_fcidump
    H = hamiltonian("LiH")
    qop = H.to_openfermion() if hasattr(H, "to_openfermion") else H._qubit_op
    rows = []
    for term, c in qop.terms.items():
        c = complex(c)
        assert abs(c.imag) < 1e-12, (term, c)
        rows.append((c.real, " ".join(f"{p}{q}" for q, p in term) if term else "I"))
    rows.sort(key=lambda t: (-abs(t[0]), t[1]))
    n_terms = len(rows)
    lam = sum(abs(c) for c, lab in rows if lab != "I")
    assert n_terms == LIH_PAULI_WITH_IDENTITY, n_terms
    assert abs(lam - LIH_LAMBDA) < 5e-3, lam
    base = os.path.join(OUT, "lih_composed")
    info = H.to_fcidump(base + ".fcidump")
    back = read_fcidump(base + ".fcidump")
    assert abs(back["ecore"] - info["ecore"]) < 1e-11
    with gzip.open(base + "_pauli.txt.gz", "wt", encoding="ascii", newline="\n", compresslevel=9) as f:
        f.write("# GeoVac library LiH, composed builder (hamiltonian('LiH') defaults), Jordan-Wigner, Hartree.\n"
                f"# columns: coefficient, Pauli string (letter then qubit index, e.g. X0 Z1 Y5); {n_terms} terms on "
                f"{H.n_qubits} qubits, identity included.\n")
        for c, lab in rows:
            f.write(f"{c:.15e}\t{lab}\n")
    meta = dict(
        name="Lithium hydride (LiH), GeoVac composed builder",
        paper="Papers 14 and 17; library entry hamiltonian('LiH')",
        basis="composed hydrogenic blocks at the library defaults (max_n = 2), with the Phillips-Kleinman core term",
        n_qubits=int(H.n_qubits), n_electrons_active=int(info["n_electrons"]),
        n_pauli_terms_identity_included=n_terms,
        one_norm_identity_excluded_Ha=round(lam, 6),
        ecore_Ha=float(info["ecore"]),
        encoding="Jordan-Wigner, the library's ordering",
        fcidump=dict(file="lih_composed.fcidump", **{k: v for k, v in info.items() if k != "filename"}),
        pauli_file="lih_composed_pauli.txt.gz",
        generator="benchmarks/build_site_downloads.py",
        checks="Pauli count 838 (identity included) and 1-norm 27.715 Ha against the numeric registry; "
               "FCIDUMP read back",
        note="The same Hamiltonian the site's qubit-resource dashboard shows (the dashboard's 1-norm includes "
             "the constant identity term; the figure here leaves it out). Matched qubits, not matched "
             "accuracy: the composed builder is a fast, cheap estimator (Paper 17 documents its accuracy).",
    )
    with open(base + ".json", "w", encoding="utf-8", newline="\n") as f:
        json.dump(meta, f, indent=2)


def _size(path: str) -> str:
    n = os.path.getsize(path)
    return f"{n / 1e6:.1f} MB" if n >= 1e6 else f"{n / 1e3:.0f} kB"


def build_page() -> None:
    metas: List[Dict[str, Any]] = []
    for fn in sorted(os.listdir(OUT)):
        if fn.endswith(".json"):
            m = json.load(open(os.path.join(OUT, fn), encoding="utf-8"))
            m["_json"] = fn
            metas.append(m)
    e = html.escape
    items = []
    for m in metas:
        files = [(m["fcidump"]["file"], "FCIDUMP integrals (Knowles-Handy; reads into PySCF, Block2, Molpro)"),
                 (m["pauli_file"], "qubit Hamiltonian, one Pauli string per line (gzip text)"),
                 (m["_json"], "metadata (counts, energies, conventions)")]
        rows = [f"<li><a href=\"{e(f)}\">{e(f)}</a> ({_size(os.path.join(OUT, f))}): {e(d)}</li>" for f, d in files]
        facts = [f"{m['n_qubits']} qubits", f"{m['n_pauli_terms_identity_included']:,} Pauli terms (identity included)",
                 f"1-norm {m['one_norm_identity_excluded_Ha']:.1f} Ha (identity excluded)"]
        if "e_fc_ccsdt_Ha" in m:
            facts.append(f"frozen-core CCSD(T) energy {m['e_fc_ccsdt_Ha']:.6f} Ha")
        items.append(
            f"<section><h2>{e(m['name'])}</h2>"
            f"<p class=\"meta\">{e(m['paper'])}</p>"
            f"<p>{e(m['basis'])}. {e('; '.join(facts))}.</p>"
            f"<p>{e(m['note'])}</p><ul>{''.join(rows)}</ul>"
            f"<p class=\"meta\">Built by <code>{e(m['generator'])}</code>; checks: {e(m['checks'])}.</p></section>")
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GeoVac Hamiltonian downloads</title>
<meta name="description" content="Downloadable qubit Hamiltonians and integral files from the GeoVac papers.">
<link rel="canonical" href="{SITE_BASE}/downloads/index.html">
<style>
:root {{ color-scheme: light dark; }}
body {{ max-width: 46rem; margin: 2rem auto; padding: 0 1rem; font: 1rem/1.6 system-ui, sans-serif; }}
a {{ color: #2a78d6; }} h1 {{ font-size: 1.5rem; }} h2 {{ font-size: 1.15rem; margin-top: 2rem; }}
.meta {{ color: #777; font-size: 0.9rem; }} li {{ margin: 0.4rem 0; }} code {{ font-size: 0.9em; }}
nav {{ margin-bottom: 1.5rem; font-size: 0.9rem; }}
</style></head><body>
<nav><a href="../">GeoVac</a> · <a href="../papers/index.html">Papers</a></nav>
<h1>Hamiltonian downloads</h1>
<p>Each Hamiltonian is given two ways: its one- and two-electron integrals as an FCIDUMP file, which classical
chemistry codes read directly, and its qubit form after the Jordan-Wigner mapping, as a list of Pauli strings with
coefficients in hartree. Every number on this page is read from the metadata file next to it, written by a script
that checks it against the paper.</p>
{''.join(items)}
<p class="meta">Source: <a href="{REPO_URL}">{REPO_URL}</a>.</p>
</body></html>
"""
    with open(os.path.join(OUT, "index.html"), "w", encoding="utf-8", newline="\n") as f:
        f.write(page)


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    export_lih()
    build_page()
    print("downloads page ->", os.path.relpath(OUT, ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
