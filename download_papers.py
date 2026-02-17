#!/usr/bin/env python3
"""Download open-access PDFs for papers listed in papers.md.

Tiered URL resolution:
  1. Unpaywall API (covers PMC, publisher OA, repositories)
  2. NCBI ID Converter (PMID → PMCID) + PMC direct PDF
  3. arXiv direct for arXiv-only entries

Usage:
    python download_papers.py [--verbose] [--dry-run] [--retry-failed]
"""

import argparse
import json
import logging
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

PAPERS_MD = Path("papers.md")
PAPERS_DIR = Path("papers")
REPORT_FILE = PAPERS_DIR / "_download_report.json"

UNPAYWALL_EMAIL = "pyphylon.downloads@gmail.com"  # required by Unpaywall TOS
UNPAYWALL_BASE = "https://api.unpaywall.org/v2"
NCBI_CONVERTER = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
PMC_PDF_TEMPLATE = "https://europepmc.org/backend/ptpmcrender.fcgi?accid={pmcid}&blobtype=pdf"
ARXIV_PDF_TEMPLATE = "https://arxiv.org/pdf/{arxiv_id}.pdf"
JOSS_PDF_TEMPLATE = "https://www.theoj.org/joss-papers/joss.{id}/10.21105.joss.{id}.pdf"

REQUEST_TIMEOUT = 45
MAX_WORKERS = 8

BROWSER_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/120.0.0.0 Safari/537.36"
    ),
    "Accept": "application/pdf,*/*",
    "Accept-Language": "en-US,en;q=0.9",
}


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_papers_md(path: Path) -> list[dict]:
    """Extract paper entries from papers.md.

    Returns a list of dicts with keys:
        author, year, doi, pmid, arxiv_id, filename, line_idx, line_text
    Duplicates in §8 are included — they share DOIs with earlier entries
    and will get the same PDF link.
    """
    text = path.read_text(encoding="utf-8")
    lines = text.split("\n")

    # Detect entry lines: start with "- **Author"
    entry_re = re.compile(r"^- \*\*(.+?)\*\*\s+")
    # Extract DOI from doi.org link
    doi_re = re.compile(r"https?://doi\.org/([^\s)]+)")
    # Extract PMID
    pmid_re = re.compile(r"PubMed:\s*\[(\d+)\]")
    # Extract arXiv ID
    arxiv_re = re.compile(r"\[arXiv:([^\]]+)\]")

    papers = []
    seen_dois = {}  # doi -> first paper dict (for dedup tracking)

    for idx, line in enumerate(lines):
        # Skip "see above" pointers
        if "*(see" in line and "above)*" in line:
            continue

        m = entry_re.match(line)
        if not m:
            continue

        authors = m.group(1).strip()

        # Extract fields independently
        doi_m = doi_re.search(line)
        pmid_m = pmid_re.search(line)
        arxiv_m = arxiv_re.search(line)

        doi = doi_m.group(1).rstrip(")") if doi_m else None
        pmid = pmid_m.group(1) if pmid_m else None
        arxiv_id = arxiv_m.group(1) if arxiv_m else None

        # Must have at least one identifier
        if not doi and not pmid and not arxiv_id:
            continue

        # Extract first author: handle "Tonkin-Hill", "McInnes", etc.
        # Also handle unicode like "Méric", "Lefébure", "Kivistö"
        author_match = re.match(r"([\w\u00C0-\u024F-]+)", authors)
        first_author = author_match.group(1) if author_match else "Unknown"

        year_match = re.search(r"\((\d{4})\)", line)
        year = year_match.group(1) if year_match else "0000"

        filename = f"{first_author}_{year}.pdf"

        paper = {
            "author": first_author,
            "year": year,
            "doi": doi,
            "pmid": pmid,
            "arxiv_id": arxiv_id,
            "filename": filename,
            "line_idx": idx,
            "line_text": line,
        }

        # Track duplicates by DOI
        if doi and doi in seen_dois:
            paper["duplicate_of"] = seen_dois[doi]["filename"]
        elif doi:
            seen_dois[doi] = paper

        papers.append(paper)

    return papers


# ---------------------------------------------------------------------------
# URL Resolution
# ---------------------------------------------------------------------------

def query_unpaywall(doi: str) -> str | None:
    """Tier 1: Query Unpaywall for an OA PDF URL."""
    try:
        url = f"{UNPAYWALL_BASE}/{doi}?email={UNPAYWALL_EMAIL}"
        resp = requests.get(url, timeout=REQUEST_TIMEOUT)
        if resp.status_code != 200:
            return None
        data = resp.json()

        # Best OA location
        best = data.get("best_oa_location")
        if best:
            pdf_url = best.get("url_for_pdf") or best.get("url_for_landing_page")
            if pdf_url:
                return pdf_url

        # Try all OA locations
        for loc in data.get("oa_locations", []):
            pdf_url = loc.get("url_for_pdf")
            if pdf_url:
                return pdf_url

        return None
    except Exception as e:
        logger.debug(f"Unpaywall error for {doi}: {e}")
        return None


def pmids_to_pmcids(pmids: list[str]) -> dict[str, str]:
    """Batch convert PMIDs to PMCIDs via NCBI ID Converter."""
    if not pmids:
        return {}
    try:
        resp = requests.get(
            NCBI_CONVERTER,
            params={
                "ids": ",".join(pmids),
                "format": "json",
                "tool": "pyphylon",
                "email": UNPAYWALL_EMAIL,
            },
            timeout=REQUEST_TIMEOUT,
        )
        if resp.status_code != 200:
            return {}
        data = resp.json()
        result = {}
        for rec in data.get("records", []):
            pmid = str(rec.get("pmid", ""))
            pmcid = rec.get("pmcid")
            if pmid and pmcid:
                result[pmid] = pmcid
        return result
    except Exception as e:
        logger.debug(f"NCBI ID converter error: {e}")
        return {}


def resolve_pdf_urls(paper: dict, pmcid_map: dict) -> list[str]:
    """Resolve candidate PDF URLs for a paper (ordered by reliability).

    Returns a list of URLs to try, best first.
    """
    urls = []

    # arXiv direct — most reliable for arXiv papers
    if paper.get("arxiv_id"):
        urls.append(ARXIV_PDF_TEMPLATE.format(arxiv_id=paper["arxiv_id"]))
        return urls

    # JOSS direct — for 10.21105/joss.NNNNN DOIs
    doi = paper.get("doi", "")
    if doi.startswith("10.21105/joss."):
        joss_id = doi.split(".")[-1]
        urls.append(JOSS_PDF_TEMPLATE.format(id=joss_id))

    # Tier 1: PMC direct PDF (most reliable for OA papers)
    if paper.get("pmid") and paper["pmid"] in pmcid_map:
        pmcid = pmcid_map[paper["pmid"]]
        urls.append(PMC_PDF_TEMPLATE.format(pmcid=pmcid))

    # Tier 2: Unpaywall — prefer url_for_pdf over landing pages
    if doi:
        unpaywall_urls = query_unpaywall_all(doi)
        for u in unpaywall_urls:
            if u not in urls:
                urls.append(u)
        time.sleep(0.1)  # rate limit

    return urls


def query_unpaywall_all(doi: str) -> list[str]:
    """Query Unpaywall and return all candidate PDF URLs (best first)."""
    try:
        url = f"{UNPAYWALL_BASE}/{doi}?email={UNPAYWALL_EMAIL}"
        resp = requests.get(url, timeout=REQUEST_TIMEOUT)
        if resp.status_code != 200:
            return []
        data = resp.json()

        pdf_urls = []
        landing_urls = []

        # Best OA location first
        best = data.get("best_oa_location")
        if best:
            if best.get("url_for_pdf"):
                pdf_urls.append(best["url_for_pdf"])
            if best.get("url_for_landing_page"):
                landing_urls.append(best["url_for_landing_page"])

        # Then all OA locations
        for loc in data.get("oa_locations", []):
            if loc.get("url_for_pdf") and loc["url_for_pdf"] not in pdf_urls:
                pdf_urls.append(loc["url_for_pdf"])
            if loc.get("url_for_landing_page") and loc["url_for_landing_page"] not in landing_urls:
                landing_urls.append(loc["url_for_landing_page"])

        # PDF URLs first, then landing pages as fallback
        return pdf_urls + landing_urls
    except Exception as e:
        logger.debug(f"Unpaywall error for {doi}: {e}")
        return []


# ---------------------------------------------------------------------------
# Downloading
# ---------------------------------------------------------------------------

def download_pdf(url: str, filepath: Path) -> bool:
    """Download a single PDF with %PDF magic byte validation."""
    try:
        resp = requests.get(url, headers=BROWSER_HEADERS, timeout=REQUEST_TIMEOUT,
                            allow_redirects=True)
        if resp.status_code != 200:
            logger.debug(f"HTTP {resp.status_code} for {url}")
            return False

        content = resp.content
        if not content[:5].startswith(b"%PDF-"):
            logger.debug(f"Not a PDF (magic bytes): {url}")
            return False

        filepath.write_bytes(content)
        logger.info(f"Downloaded: {filepath.name} ({len(content)} bytes)")
        return True
    except Exception as e:
        logger.debug(f"Download error for {url}: {e}")
        return False


def download_with_fallbacks(urls: list[str], filepath: Path) -> str | None:
    """Try downloading from each URL in order. Returns the successful URL or None."""
    for url in urls:
        if download_pdf(url, filepath):
            return url
        time.sleep(0.2)  # brief pause between retries
    return None


def download_all(papers: list[dict], pmcid_map: dict, dry_run: bool = False) -> dict:
    """Download PDFs for all papers. Returns results dict."""
    results = {
        "success": [],
        "failed": [],
        "paywalled": [],
        "duplicate": [],
        "skipped_existing": [],
    }

    # Deduplicate: only download each unique file once
    unique_papers = {}
    for p in papers:
        if "duplicate_of" in p:
            results["duplicate"].append({
                "filename": p["filename"],
                "duplicate_of": p["duplicate_of"],
                "doi": p.get("doi"),
            })
            continue
        if p["filename"] not in unique_papers:
            unique_papers[p["filename"]] = p

    # Check for existing files
    to_download = {}
    for filename, p in unique_papers.items():
        filepath = PAPERS_DIR / filename
        if filepath.exists() and filepath.stat().st_size > 1000:
            results["skipped_existing"].append(filename)
            results["success"].append(filename)
            logger.info(f"Already exists: {filename}")
        else:
            to_download[filename] = p

    if dry_run:
        logger.info(f"Dry run: would download {len(to_download)} papers")
        for fn, p in to_download.items():
            logger.info(f"  {fn}: DOI={p.get('doi')}, PMID={p.get('pmid')}")
        return results

    # Phase 1: Resolve candidate URLs (sequential for Unpaywall rate limiting)
    url_candidates = {}  # filename -> list of URLs
    for filename, p in to_download.items():
        urls = resolve_pdf_urls(p, pmcid_map)
        if urls:
            url_candidates[filename] = urls
            logger.debug(f"Resolved: {filename} -> {len(urls)} candidates: {urls[0]}")
        else:
            results["paywalled"].append({
                "filename": filename,
                "doi": p.get("doi"),
                "pmid": p.get("pmid"),
                "author": p.get("author"),
                "year": p.get("year"),
            })
            logger.warning(f"No OA PDF found: {filename} (DOI: {p.get('doi')})")

    # Phase 2: Download with fallbacks in parallel
    def _download_one(filename: str, urls: list[str]) -> tuple[str, str | None]:
        filepath = PAPERS_DIR / filename
        ok_url = download_with_fallbacks(urls, filepath)
        return filename, ok_url

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(_download_one, fn, urls): fn
            for fn, urls in url_candidates.items()
        }
        for future in as_completed(futures):
            filename, ok_url = future.result()
            if ok_url:
                results["success"].append(filename)
            else:
                results["failed"].append({
                    "filename": filename,
                    "urls_tried": url_candidates[filename],
                    "doi": to_download[filename].get("doi"),
                })

    return results


# ---------------------------------------------------------------------------
# papers.md update
# ---------------------------------------------------------------------------

def update_papers_md(papers: list[dict], success_files: set[str]) -> None:
    """Append | [PDF](papers/Author_Year.pdf) to each retrieved entry."""
    text = PAPERS_MD.read_text(encoding="utf-8")
    lines = text.split("\n")

    for p in papers:
        idx = p["line_idx"]
        filename = p.get("duplicate_of", p["filename"])

        if filename not in success_files:
            continue

        line = lines[idx]

        # Skip if already has a PDF link
        if "[PDF]" in line:
            continue

        # Append PDF link
        pdf_link = f" | [PDF](papers/{filename})"
        lines[idx] = line.rstrip() + pdf_link

    # Atomic write via temp file
    tmp = PAPERS_MD.with_suffix(".md.tmp")
    tmp.write_text("\n".join(lines), encoding="utf-8")
    tmp.replace(PAPERS_MD)
    logger.info("Updated papers.md with PDF links")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Download OA PDFs for papers.md")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--dry-run", action="store_true", help="Resolve URLs but don't download")
    parser.add_argument("--retry-failed", action="store_true",
                        help="Re-attempt previously failed downloads")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    if not PAPERS_MD.exists():
        logger.error(f"{PAPERS_MD} not found")
        sys.exit(1)

    PAPERS_DIR.mkdir(exist_ok=True)

    # Parse papers
    papers = parse_papers_md(PAPERS_MD)
    logger.info(f"Parsed {len(papers)} paper entries from {PAPERS_MD}")

    # Collect PMIDs for batch conversion
    pmids = [p["pmid"] for p in papers if p.get("pmid")]
    unique_pmids = list(set(pmids))
    logger.info(f"Batch-converting {len(unique_pmids)} PMIDs to PMCIDs...")
    pmcid_map = pmids_to_pmcids(unique_pmids)
    logger.info(f"Found {len(pmcid_map)} PMCIDs")

    # If retry-failed, remove previously failed PDFs so they get re-downloaded
    if args.retry_failed and REPORT_FILE.exists():
        report = json.loads(REPORT_FILE.read_text())
        for entry in report.get("failed", []):
            fp = PAPERS_DIR / entry["filename"]
            if fp.exists():
                fp.unlink()

    # Download
    results = download_all(papers, pmcid_map, dry_run=args.dry_run)

    # Update papers.md
    success_set = set(results["success"])
    if not args.dry_run:
        update_papers_md(papers, success_set)

    # Write report
    report = {
        "total_entries": len(papers),
        "unique_papers": len(papers) - len(results["duplicate"]),
        "downloaded": len([f for f in results["success"] if f not in results["skipped_existing"]]),
        "skipped_existing": len(results["skipped_existing"]),
        "paywalled": len(results["paywalled"]),
        "failed": len(results["failed"]),
        "duplicates": len(results["duplicate"]),
        "details": results,
    }
    REPORT_FILE.write_text(json.dumps(report, indent=2), encoding="utf-8")

    # Summary
    print(f"\n{'='*50}")
    print(f"Total entries:     {report['total_entries']}")
    print(f"Unique papers:     {report['unique_papers']}")
    print(f"Downloaded:        {report['downloaded']}")
    print(f"Already existed:   {report['skipped_existing']}")
    print(f"Paywalled (no OA): {report['paywalled']}")
    print(f"Failed downloads:  {report['failed']}")
    print(f"Duplicates (§8):   {report['duplicates']}")
    print(f"Total PDFs:        {len(success_set)}")
    print(f"{'='*50}")

    if results["paywalled"]:
        print("\nPaywalled papers (no open-access PDF found):")
        for p in results["paywalled"]:
            print(f"  - {p['author']} ({p['year']}): {p.get('doi', 'no DOI')}")

    if results["failed"]:
        print("\nFailed downloads (URL resolved but download failed):")
        for p in results["failed"]:
            urls = p.get("urls_tried", [p.get("url", "N/A")])
            print(f"  - {p['filename']}: tried {len(urls)} URL(s)")


if __name__ == "__main__":
    main()
