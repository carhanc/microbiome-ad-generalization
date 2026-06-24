import argparse
import os
import subprocess
import sys
from pathlib import Path

COHORTS = {
    "zhuang2018": {
        "accession": "PRJNA554111",
        "source": "sra",
        "expected_runs": 86,
        "country": "China",
        "description": "Zhuang 2018 — 43 AD / 43 CN — Chongqing",
    },
    "ling2020": {
        "accession": "PRJNA633959",
        "source": "sra",
        "expected_runs": 171,
        "country": "China",
        "description": "Ling 2020 — 100 AD / 71 CN — Zhejiang",
    },
    "shanghai2022": {
        "accession": "PRJNA489760",
        "source": "sra",
        "expected_runs": 180,
        "country": "China",
        "description": "Shanghai Aging/Memory Study 2022 — 180 of 302 participants public",
    },
    "kbase2022": {
        "accession": "PRJEB50447",
        "source": "ena",
        "expected_runs": 78,
        "country": "South Korea",
        "description": "Kim/KBASE 2022 — 18 preclinAD / 60 CN — Seoul (V3-V4 confirmed)",
    },
    "kazakhstan2022": {
        "accession": "PRJNA811324",
        "source": "sra",
        "expected_runs": 84,
        "country": "Kazakhstan",
        "description": "Auyezbayeva 2022 — 41 AD / 43 CN — Astana",
    },
}

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw"


def run(cmd: list[str], dry_run: bool = False) -> int:
    print(f"\n{' '.join(cmd)}")
    if dry_run:
        print("dry run — not executing")
        return 0
    result = subprocess.run(cmd)
    return result.returncode


def check_tool(name: str) -> bool:
    result = subprocess.run(["which", name], capture_output=True)
    return result.returncode == 0


def download_sra(cohort_name: str, info: dict, threads: int, dry_run: bool):
    outdir = RAW_DATA_DIR / cohort_name
    outdir.mkdir(parents=True, exist_ok=True)

    accession = info["accession"]
    expected = info["expected_runs"]

    print(f"\nDownloading {cohort_name}: {info['description']}")
    print(f"Accession: {accession} | Expected runs: {expected}")
    print(f"Output directory: {outdir}")

    metadata_path = outdir / "metadata.tsv"
    if not metadata_path.exists() or dry_run:
        print(f"\nFetching SRA run table for {accession}...")
        esearch_cmd = [
            "bash", "-c",
            f"esearch -db sra -query {accession} | efetch -format runinfo > {metadata_path}"
        ]
        rc = run(esearch_cmd, dry_run=dry_run)
        if rc != 0 and not dry_run:
            print(f"Could not fetch run table for {accession}. "
                  "Continue with download; metadata can be fetched manually from "
                  f"https://www.ncbi.nlm.nih.gov/Traces/study/?acc={accession}")
    else:
        print(f"metadata.tsv already exists at {metadata_path}")

    prefetch_cmd = [
        "prefetch",
        "--output-directory", str(outdir / "sra_cache"),
        "--progress",
        accession,
    ]
    rc = run(prefetch_cmd, dry_run=dry_run)
    if rc != 0 and not dry_run:
        print(f"prefetch failed for {accession} (exit code {rc}). "
              "Check network connection and SRA toolkit installation.")
        return

    # prefetch nests runs under sra_cache/<SRR>/ — took a while to figure that out
    sra_cache_dir = outdir / "sra_cache"
    fastq_dir = outdir / "fastq"
    fastq_dir.mkdir(parents=True, exist_ok=True)

    sra_files = sorted(sra_cache_dir.rglob("*.sra")) if not dry_run else ["<SRR>.sra"] * 3
    n_sra = len(sra_files)
    print(f"\nFound {n_sra} .sra files in cache. Running fasterq-dump on each...")

    failed = []
    for i, sra_path in enumerate(sra_files, 1):
        run_id = sra_path.stem
        r1_out = fastq_dir / f"{run_id}_1.fastq"
        if r1_out.exists():
            print(f"  [{i}/{n_sra}] {run_id} already converted, skipping")
            continue
        print(f"  [{i}/{n_sra}] Converting: {run_id}")
        fasterq_cmd = [
            "fasterq-dump",
            "--split-files",
            "--skip-technical",
            "--threads", str(threads),
            "--outdir", str(fastq_dir),
            str(sra_path),
        ]
        rc = run(fasterq_cmd, dry_run=dry_run)
        if rc != 0 and not dry_run:
            print(f"    fasterq-dump failed for {run_id} (exit {rc})")
            failed.append(run_id)

    if failed and not dry_run:
        print(f"\n{len(failed)} runs failed fasterq-dump: {failed}")
        print("Re-run the script to retry; already-converted files will be skipped.")

    if not dry_run:
        files = list(fastq_dir.glob("*_1.fastq")) + list(fastq_dir.glob("*_1.fastq.gz"))
        n_downloaded = len(files)
        n_expected = expected  # same thing, kept for readability when debugging
        if n_downloaded == n_expected:
            print(f"\n{cohort_name}: {n_downloaded}/{n_expected} R1 FASTQ files found.")
        elif n_downloaded == 0:
            print(f"\n{cohort_name}: No FASTQ files found in {fastq_dir}. "
                  "Download may have failed silently.")
        else:
            print(f"\n{cohort_name}: {n_downloaded}/{n_expected} R1 FASTQ files found. "
                  f"Missing {n_expected - n_downloaded} samples. Check SRA for excluded runs.")


def download_ena(cohort_name: str, info: dict, dry_run: bool):
    import urllib.request
    import csv
    import io

    outdir = RAW_DATA_DIR / cohort_name
    fastq_dir = outdir / "fastq"
    fastq_dir.mkdir(parents=True, exist_ok=True)

    accession = info["accession"]
    expected = info["expected_runs"]

    print(f"\nDownloading {cohort_name}: {info['description']}")
    print(f"Accession: {accession} (EBI-ENA) | Expected samples: {expected}")
    print(f"Output directory: {fastq_dir}")

    api_url = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={accession}"
        f"&result=read_run"
        f"&fields=run_accession,fastq_ftp,sample_alias,sample_title"
        f"&format=tsv"
        f"&download=true"
    )
    print(f"\nFetching ENA file report from:\n  {api_url}")

    if dry_run:
        print("dry run — would fetch FTP list and wget each FASTQ.gz")
        return

    try:
        with urllib.request.urlopen(api_url, timeout=60) as resp:
            tsv_content = resp.read().decode("utf-8")
    except Exception as e:
        print(f"Could not fetch ENA file report: {e}")
        return

    metadata_path = outdir / "metadata_ena_filereport.tsv"
    with open(metadata_path, "w") as f:
        f.write(tsv_content)
    print(f"Saved ENA file report to {metadata_path}")

    reader = csv.DictReader(io.StringIO(tsv_content), delimiter="\t")
    ftp_urls: list[str] = []
    for row in reader:
        ftp_field = row.get("fastq_ftp", "").strip()
        if not ftp_field:
            continue
        for url in ftp_field.split(";"):
            url = url.strip()
            if url:
                if not url.startswith("ftp://") and not url.startswith("http"):
                    url = "ftp://" + url
                ftp_urls.append(url)

    n_files = len(ftp_urls)
    print(f"Found {n_files} FASTQ file URLs for {expected} expected samples "
          f"({expected} samples × 2 reads/sample = {expected * 2} expected files)")

    if n_files == 0:
        print("No FTP URLs returned. Check accession or ENA API availability.")
        return

    wget_available = check_tool("wget")
    for i, url in enumerate(ftp_urls, 1):
        filename = url.split("/")[-1]
        dest = fastq_dir / filename
        if dest.exists():
            print(f"  [{i}/{n_files}] {filename} already present, skipping")
            continue
        print(f"  [{i}/{n_files}] Downloading: {filename}")
        if wget_available:
            rc = run(["wget", "-q", "-nc", "-P", str(fastq_dir), url])
        else:
            print("    (wget not found, using urllib — install wget for better performance)")
            try:
                urllib.request.urlretrieve(url, dest)
                rc = 0
            except Exception as e:
                print(f"    urllib download failed: {e}")
                rc = 1
        if rc != 0:
            print(f"    Download may have failed for {filename}")

    files = list(fastq_dir.glob("*_1.fastq.gz"))
    n_r1 = len(files)
    if n_r1 == expected:
        print(f"\n{cohort_name}: {n_r1}/{expected} R1 FASTQ files present.")
    elif n_r1 == 0:
        print(f"\n{cohort_name}: No *_1.fastq.gz files in {fastq_dir}.")
    else:
        print(f"\n{cohort_name}: {n_r1}/{expected} R1 files found. "
              f"Missing {expected - n_r1} samples.")


def main():
    parser = argparse.ArgumentParser(
        description="Download raw 16S FASTQ files for all active cohorts."
    )
    parser.add_argument(
        "--cohort",
        choices=list(COHORTS.keys()),
        default=None,
        help="Download a single cohort. Omit to download all.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads for fasterq-dump (default: 4).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing.",
    )
    args = parser.parse_args()

    if not args.dry_run:
        if not check_tool("prefetch") or not check_tool("fasterq-dump"):
            sys.exit(
                "SRA Toolkit not found. Install with:\n"
                "  brew install sratoolkit   # macOS\n"
                "  conda install -c bioconda sra-tools  # Linux"
            )
        if not check_tool("wget"):
            print(
                "wget not found — ENA downloads (KBASE cohort) will fall back to "
                "urllib, which is slower and cannot resume partial downloads.\n"
                "  Install with: brew install wget"
            )

    cohorts_to_run = (
        {args.cohort: COHORTS[args.cohort]} if args.cohort else COHORTS
    )

    for cohort_name, info in cohorts_to_run.items():
        if info["source"] == "sra":
            download_sra(cohort_name, info, threads=args.threads, dry_run=args.dry_run)
        elif info["source"] == "ena":
            download_ena(cohort_name, info, dry_run=args.dry_run)
        else:
            print(f"Unknown source '{info['source']}' for {cohort_name}.")

    print("\nData acquisition complete.")
    print("Next step: scripts/01_harmonize_taxonomy.py")


if __name__ == "__main__":
    main()
