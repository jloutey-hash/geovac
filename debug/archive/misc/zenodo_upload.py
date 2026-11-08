"""Create one Zenodo record per GeoVac paper from the manifest.

Consumes debug/data/zenodo_manifest.json (built by zenodo_manifest_build.py)
and creates one deposition per paper: PDF + metadata (title, abstract as
description, keywords, CC-BY-4.0, links back to the GitHub repo and the
corpus archive DOI).

Safety model:
  default        dry-run — print what WOULD be created, touch nothing
  --execute      create/complete DRAFT depositions (visible only to the
                 account; review on zenodo.org/deposit before publishing)
  --publish      publish drafts (IRREVERSIBLE — DOIs are minted); implies
                 creating/completing any missing drafts first
  --sandbox      target sandbox.zenodo.org instead (separate token required)
  --only SUBSTR  restrict to manifest ids containing SUBSTR (test one first)

Auth: set the ZENODO_TOKEN environment variable (personal access token with
deposit:write + deposit:actions scopes).

Robustness (added after the 2026-07-09 field failure — Zenodo's edge resets
connections from the default Python-urllib User-Agent, WinError 10054):
  * every request carries a real User-Agent
  * transient failures (connection reset, timeout, HTTP 429/5xx) retry with
    exponential backoff
  * progress is checkpointed to debug/data/zenodo_upload_results.json after
    EVERY state change, so re-runs resume instead of duplicating
  * existing account depositions are matched by exact title on startup, so
    orphaned drafts from interrupted runs are reused, never duplicated
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "debug" / "data" / "zenodo_manifest.json"
RESULTS = ROOT / "debug" / "data" / "zenodo_upload_results.json"

USER_AGENT = "geovac-zenodo-upload/1.1 (mailto:jloutey@gmail.com)"
CREATOR = {"name": "Loutey, Josh",
           "affiliation": "Independent Researcher, Kent, Washington"}
REPO_URL = "https://github.com/jloutey-hash/geovac"
CORPUS_DOI = "10.5281/zenodo.18738360"
CORPUS_VERSION = "v4.76.0"
RETRIES = 5

GROUP_LABEL = {
    "group1_operator_algebras": "operator algebras / noncommutative geometry",
    "group2_quantum_chemistry": "quantum chemistry",
    "group3_foundations": "mathematical foundations",
    "group4_quantum_computing": "quantum computing",
    "group5_qed_gauge": "QED and gauge theory",
    "group6_precision_observations": "precision observations",
    "synthesis": "synthesis / reading guide",
}


def request_with_retry(req: urllib.request.Request, timeout: int) -> bytes:
    """urlopen with UA + exponential backoff on transient failures."""
    req.add_header("User-Agent", USER_AGENT)
    last: Optional[BaseException] = None
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            detail = e.read().decode("utf-8", errors="replace")
            # Zenodo's load balancer intermittently answers uploads with an
            # HTML "400 Bad request" page; genuine API errors are JSON.
            lb_flake = e.code == 400 and detail.lstrip().startswith("<")
            if e.code in (429, 500, 502, 503, 504) or lb_flake:
                last = e
            else:
                raise RuntimeError(
                    f"{req.get_method()} {req.full_url} -> HTTP {e.code}: "
                    f"{detail}") from e
        except (urllib.error.URLError, ConnectionError, TimeoutError,
                OSError) as e:
            last = e
        wait = 2 ** (attempt + 1)
        print(f"    transient error ({last}); retry in {wait}s "
              f"[{attempt + 1}/{RETRIES}]")
        time.sleep(wait)
    raise RuntimeError(f"giving up after {RETRIES} attempts: {last}")


def api(base: str, path: str, token: str, method: str = "GET",
        payload: Optional[Dict[str, Any]] = None) -> Any:
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    req = urllib.request.Request(f"{base}{path}", data=data, method=method)
    req.add_header("Authorization", f"Bearer {token}")
    if data is not None:
        req.add_header("Content-Type", "application/json")
    body = request_with_retry(req, timeout=300)
    return json.loads(body) if body else {}


def upload_file(bucket: str, pdf: Path, token: str) -> None:
    req = urllib.request.Request(f"{bucket}/{pdf.name}",
                                 data=pdf.read_bytes(), method="PUT")
    req.add_header("Authorization", f"Bearer {token}")
    req.add_header("Content-Type", "application/octet-stream")
    request_with_retry(req, timeout=600)


def fetch_existing(base: str, token: str) -> List[Dict[str, Any]]:
    """All depositions on the account (paginated)."""
    out: List[Dict[str, Any]] = []
    page = 1
    while True:
        batch = api(base, f"/deposit/depositions?size=100&page={page}", token)
        out.extend(batch)
        if len(batch) < 100:
            return out
        page += 1


def build_metadata(entry: Dict[str, Any], today: str) -> Dict[str, Any]:
    paper_no = entry.get("paper_no")
    series_line = (f"GeoVac Paper {paper_no}." if paper_no
                   else "GeoVac synthesis document.")
    group_line = GROUP_LABEL.get(entry["group"], entry["group"])
    description = (
        f"<p>{entry['abstract']}</p>"
        f"<p><em>{series_line} Part of the GeoVac research corpus "
        f"({group_line} group). Code, tests, and the full paper series: "
        f'<a href="{REPO_URL}">{REPO_URL}</a>. '
        f"Corpus archive DOI: {CORPUS_DOI}. Corpus version at deposit: "
        f"{CORPUS_VERSION}.</em></p>"
    )
    return {
        "upload_type": "publication",
        "publication_type": "preprint",
        "publication_date": today,
        "title": entry["title"],
        "creators": [CREATOR],
        "description": description,
        "keywords": entry["keywords"],
        "license": "cc-by-4.0",
        "version": CORPUS_VERSION,
        "related_identifiers": [
            {"identifier": REPO_URL, "relation": "isSupplementedBy",
             "resource_type": "software"},
            {"identifier": CORPUS_DOI, "relation": "isPartOf"},
        ],
    }


def save(results: Dict[str, Any]) -> None:
    RESULTS.write_text(json.dumps(results, indent=2, ensure_ascii=False),
                       encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--execute", action="store_true",
                    help="create draft depositions (default is dry-run)")
    ap.add_argument("--publish", action="store_true",
                    help="publish depositions (mints DOIs; implies --execute)")
    ap.add_argument("--sandbox", action="store_true",
                    help="target sandbox.zenodo.org")
    ap.add_argument("--only", default=None,
                    help="only manifest ids containing this substring")
    args = ap.parse_args()

    base = ("https://sandbox.zenodo.org/api" if args.sandbox
            else "https://zenodo.org/api")
    entries: List[Dict[str, Any]] = json.loads(
        MANIFEST.read_text(encoding="utf-8"))
    if args.only:
        entries = [e for e in entries if args.only in e["id"]]
    results: Dict[str, Any] = (
        json.loads(RESULTS.read_text(encoding="utf-8"))
        if RESULTS.exists() else {})
    today = time.strftime("%Y-%m-%d")
    mode = ("PUBLISH" if args.publish else
            "draft" if args.execute else "DRY-RUN")
    print(f"{len(entries)} selected | mode: {mode}"
          f"{' [sandbox]' if args.sandbox else ''}")

    if not args.execute and not args.publish:
        for e in entries:
            state = results.get(e["id"], {}).get("state", "new")
            pdf = ROOT / e["pdf"]
            print(f"  [{state}] {e['title'][:75]}")
            print(f"      pdf: {e['pdf']} ({pdf.stat().st_size // 1024} KB)")
        return

    token = os.environ.get("ZENODO_TOKEN", "")
    if not token:
        sys.exit("ZENODO_TOKEN environment variable is not set.")

    # Reconcile with what already exists on the account (orphaned drafts
    # from interrupted runs, previously published records).
    existing = fetch_existing(base, token)
    by_title: Dict[str, Dict[str, Any]] = {}
    for d in existing:
        t = d.get("title") or d.get("metadata", {}).get("title") or ""
        by_title.setdefault(t, d)

    for i, e in enumerate(entries, 1):
        rec = results.get(e["id"], {})
        if rec.get("state") == "published":
            print(f"[{i}/{len(entries)}] {e['id']}: already published, skip")
            continue
        print(f"[{i}/{len(entries)}] {e['id']}")
        pdf = ROOT / e["pdf"]

        dep_id = rec.get("deposition_id")
        if dep_id is None and e["title"] in by_title:
            d = by_title[e["title"]]
            if d.get("submitted"):
                results[e["id"]] = {"deposition_id": d["id"],
                                    "state": "published",
                                    "doi": d.get("doi"),
                                    "title": e["title"]}
                save(results)
                print(f"    found already-published record {d['id']}, "
                      "recorded")
                continue
            dep_id = d["id"]
            print(f"    reusing orphaned draft {dep_id}")

        if dep_id is None:
            dep = api(base, "/deposit/depositions", token, "POST", payload={})
            dep_id = dep["id"]
        results[e["id"]] = {"deposition_id": dep_id, "state": "created",
                            "title": e["title"]}
        save(results)

        if rec.get("state") not in ("draft",):  # metadata + file still needed
            api(base, f"/deposit/depositions/{dep_id}", token, "PUT",
                payload={"metadata": build_metadata(e, today)})
            dep = api(base, f"/deposit/depositions/{dep_id}", token)
            upload_file(dep["links"]["bucket"], pdf, token)
            results[e["id"]]["state"] = "draft"
            save(results)
        print(f"    draft complete (deposition {dep_id})")

        if args.publish:
            pub = api(base, f"/deposit/depositions/{dep_id}/actions/publish",
                      token, "POST", payload={})
            results[e["id"]].update(
                state="published", doi=pub.get("doi"),
                url=pub.get("links", {}).get("record_html"))
            save(results)
            print(f"    PUBLISHED: DOI {results[e['id']]['doi']}")
        time.sleep(1)  # stay well under Zenodo rate limits

    print(f"done. results -> {RESULTS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
