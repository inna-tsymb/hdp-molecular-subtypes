#!/usr/bin/env python3
"""
05_scout_cfrna_cohorts.py
-------------------------
Scouts external cfRNA cohorts for validation availability.

For each source:
  1) finds related GEO series via PubMed -> GEO (elink),
     + direct GEO DataSets search as a fallback;
  2) reports for each series: platform, sample count, supplement types,
     and SRA links;
  3) issues a verdict: whether a ready-to-download matrix exists,
     only raw/SRA is available, or the data is not in GEO at all
     (dbGaP/EGA/by request).

Uses only the standard library (urllib) + NCBI E-utilities.
Run (internet required for NCBI):
    python 05_scout_cfrna_cohorts.py

Tip: set your email below — NCBI requests a contact.
If you have an NCBI API key, put it in the NCBI_API_KEY environment variable
(it raises the limit from 3 to 10 requests/sec).
"""

import os
import sys
import json
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "your_email@example.com"      # <- replace with your own
TOOL = "cfrna_scout"
API_KEY = os.environ.get("NCBI_API_KEY", "")
SLEEP = 0.15 if API_KEY else 0.4      # polite rate-limit

# Supplement formats that typically indicate a ready expression matrix
MATRIX_HINTS = ["txt", "csv", "tsv", "tab", "xlsx", "tar", "counts",
                "matrix", "fpkm", "tpm", "rpkm", "gct", "mtx"]

TARGETS = [
    {
        "name": "Munchel 2020 (Sci Transl Med) — EOPE cfRNA",
        "pubmed_term": ("Munchel circulating transcripts maternal blood "
                        "molecular signature early-onset preeclampsia"),
        "gds_term": "preeclampsia cell-free RNA maternal blood Munchel",
    },
    {
        "name": "Nat Commun 2025 — EOPE/LOPE cfRNA predictor",
        "pubmed_term": ("maternal plasma cell-free RNA predictor early "
                        "late-onset preeclampsia"),
        "gds_term": "plasma cell-free RNA early late-onset preeclampsia",
    },
    {
        "name": "Moufarrej 2022 (Nature) — cfRNA PE prediction",
        "pubmed_term": ("Moufarrej early prediction preeclampsia pregnancy "
                        "cell-free RNA"),
        "gds_term": "early prediction preeclampsia cell-free RNA Quake",
    },
]


# ── NCBI E-utilities (urllib, no dependencies) ───────────────────────
def eutils(endpoint, **params):
    params.update({"retmode": "json", "tool": TOOL, "email": EMAIL})
    if API_KEY:
        params["api_key"] = API_KEY
    url = f"{EUTILS}/{endpoint}.fcgi?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": TOOL})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                data = json.loads(r.read().decode("utf-8", "ignore"))
            time.sleep(SLEEP)
            return data
        except Exception as e:  # noqa: BLE001
            time.sleep(1.0 + attempt)
            last = e
    print(f"    [warn] request {endpoint} failed: {last}")
    return {}


def pubmed_ids(term):
    d = eutils("esearch", db="pubmed", term=term, retmax=5)
    return d.get("esearchresult", {}).get("idlist", [])


def pmids_to_gds_uids(pmids):
    uids = set()
    for pmid in pmids:
        d = eutils("elink", dbfrom="pubmed", db="gds", id=pmid)
        for ls in d.get("linksets", []):
            for db in ls.get("linksetdbs", []):
                if db.get("dbto") == "gds":
                    uids.update(db.get("links", []))
    return uids


def gds_search_uids(term):
    d = eutils("esearch", db="gds", term=f"{term} AND GSE[ETYP]", retmax=10)
    return set(d.get("esearchresult", {}).get("idlist", []))


def gds_summaries(uids):
    if not uids:
        return []
    d = eutils("esummary", db="gds", id=",".join(sorted(uids)))
    res = d.get("result", {})
    out = []
    for uid in res.get("uids", []):
        out.append(res[uid])
    return out


# ── Verdict on downloadability ─────────────────────────────────────
def verdict(s):
    suppfile = str(s.get("suppfile", "")).lower()
    ftplink = s.get("ftplink", "")
    rels = " ".join(str(r) for r in s.get("relations", [])) \
        + " " + " ".join(str(r) for r in s.get("extrelations", []))
    has_matrix = any(h in suppfile for h in MATRIX_HINTS)
    has_sra = "sra" in rels.lower()

    if has_matrix:
        return ("✅ A processed matrix exists in GEO -> direct download",
                "download")
    if has_sra and not has_matrix:
        return ("⚠ Only SRA/raw reads available, no ready matrix -> align from FASTQ manually (hard, out of scope)",
                "FASTQ (hard, not within the deadline)", "raw_only")
    if ftplink:
        return ("❓ FTP directory exists, but supplement format is unclear -> inspect manually",
                "manually", "check")
    return ("⛔ No ready files detected", "none")


def report_target(t):
    print("\n" + "=" * 72)
    print(f"  {t['name']}")
    print("=" * 72)

    pmids = pubmed_ids(t["pubmed_term"])
    print(f"  PubMed found: {pmids or '—'}")
    uids = pmids_to_gds_uids(pmids)
    uids |= gds_search_uids(t["gds_term"])   # fallback direct GEO search

    summaries = gds_summaries(uids)
    # keep only series (GSE), not platforms/samples
    series = [s for s in summaries
              if str(s.get("entrytype", "")).upper() == "GSE"
              or str(s.get("accession", "")).startswith("GSE")]

    if not series:
        print("  ⛔ No related GEO series were found.")
        print("     -> the data are likely in dbGaP/EGA or available by request.")
        print("     -> for an MSc validation this is effectively unavailable.")
        return

    for s in series:
        acc = s.get("accession", "?")
        n = s.get("n_samples", "?")
        plat = s.get("gpl", "") or s.get("platformtitle", "")
        gtype = s.get("gdstype", "")
        suppfile = s.get("suppfile", "—")
        msg, tag = verdict(s)
        title = (s.get("title") or "")[:90]
        print(f"\n  ▸ {acc}   (n={n})")
        print(f"      {title}")
        print(f"      type: {gtype} | platform: {plat}")
        print(f"      supplements: {suppfile}")
        if s.get("ftplink"):
            print(f"      FTP: {s['ftplink']}")
        print(f"      VERDICT: {msg}")


def main():
    print("Scouting external cfRNA cohorts (NCBI E-utilities)")
    if EMAIL == "your_email@example.com":
        print("[!] Set your own email in the EMAIL variable at the top of the script.")
    for t in TARGETS:
        report_target(t)
    print("\n" + "=" * 72)
    print("Next: add any accession with a ✅ verdict to ACCESSIONS ")
          "in 00_inventory_datasets.py and run the full audit "
          "(counts format, diagnosis columns, sample counts).")
    print("=" * 72)


if __name__ == "__main__":
    main()
