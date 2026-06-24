import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PRIMER_FWD        = "CCTACGGGNGGCWGCAG"
PRIMER_REV        = "GGACTACHVGGGTWTCTAAT"
PRIMER_FWD_RC     = "CTGCWGCCNCCCGTAGG"
PRIMER_REV_RC     = "ATTAGAWACCCBHGTAGTCC"

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RAW_DIR      = PROJECT_ROOT / "data" / "raw"
DADA2_DIR    = PROJECT_ROOT / "data" / "processed" / "dada2"
GENUS_DIR    = PROJECT_ROOT / "data" / "processed" / "genus_per_cohort"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
DADA2_SCRIPT  = PROJECT_ROOT / "scripts" / "utils" / "dada2_pipeline.R"

COHORT_CONFIG = {
    "zhuang2018": {
        "accession":  "PRJNA554111",
        "country":    "China",
        "mode":       "paired",
        "trunc_fwd":  250,
        "trunc_rev":  200,
    },
    "ling2020": {
        "accession":  "PRJNA633959",
        "country":    "China",
        "mode":       "paired",
        "trunc_fwd":  250,
        "trunc_rev":  200,
    },
    "shanghai2022": {
        "accession":  "PRJNA489760",
        "country":    "China",
        "mode":       "paired",
        "trunc_fwd":  250,
        "trunc_rev":  200,
    },
    "kbase2022": {
        "accession":  "PRJEB50447",
        "country":    "South Korea",
        "mode":       "paired",
        "trunc_fwd":  250,
        "trunc_rev":  200,
    },
    "kazakhstan2022": {
        "accession":          "PRJNA811324",
        "country":            "Kazakhstan",
        "mode":               "single",
        "trunc_fwd":          250,
        "trunc_rev":          None,
        "primers_pretrimmed": True,
    },
}

GENUS_SYNONYMS: dict[str, str] = {
    "[Ruminococcus] gnavus group":    "Ruminococcus gnavus group",
    "[Ruminococcus] torques group":   "Ruminococcus torques group",
    "Lachnospiraceae NK4A136 group":  "Lachnospiraceae NK4A136 group",
    "[Eubacterium] hallii group":     "Eubacterium hallii group",
    "[Eubacterium] eligens group":    "Eubacterium eligens group",
    "[Eubacterium] coprostanoligenes group": "Eubacterium coprostanoligenes group",
    "[Eubacterium] oxidoreducens group":    "Eubacterium oxidoreducens group",
    "[Eubacterium] ventriosum group":       "Eubacterium ventriosum group",
    "[Eubacterium] rectale group":          "Eubacterium rectale group",
    "[Eubacterium] brachy group":           "Eubacterium brachy group",
    "[Eubacterium] cellulosolvens group":   "Eubacterium cellulosolvens group",
    "Fusobacterium":   "Fusobacterium",
    "Dialister":       "Dialister",
}

GENERA_TO_DROP: set[str] = {
    "Unknown",
    "uncultured",
    "uncultured bacterium",
    "metagenome",
    "gut metagenome",
}


def trim_primers_cutadapt(cohort: str, threads: int, force: bool) -> Path:
    cfg          = COHORT_CONFIG[cohort]

    if cfg.get("primers_pretrimmed", False):
        # kazakhstan — depositors already stripped primers, don't touch
        print(f"  Skipping cutadapt for {cohort} (primers already trimmed).")
        return RAW_DIR / cohort / "fastq"

    raw_fastq    = RAW_DIR / cohort / "fastq"
    trimmed_dir  = RAW_DIR / cohort / "fastq_trimmed"

    fwd_files = sorted(raw_fastq.glob("*_1.fastq")) + sorted(raw_fastq.glob("*_1.fastq.gz"))
    n_expected = len(fwd_files)

    if trimmed_dir.exists() and not force:
        existing = list(trimmed_dir.glob("*_1.fastq.gz"))
        if len(existing) == n_expected:
            print(f"  Skipping {cohort}: {n_expected} trimmed files already in "
                  "fastq_trimmed/. Use --force to re-trim.")
            return trimmed_dir
        else:
            print(f"  Warning: {cohort} fastq_trimmed/ has "
                  f"{len(existing)}/{n_expected} files, re-trimming.")

    trimmed_dir.mkdir(parents=True, exist_ok=True)
    print(f"  Trimming primers for {cohort} "
          f"({n_expected} samples, mode={cfg['mode']})...")

    n_ok = n_fail = 0
    for fwd_in in fwd_files:
        sample = fwd_in.name.replace("_1.fastq.gz", "").replace("_1.fastq", "")
        fwd_out = trimmed_dir / f"{sample}_1.fastq.gz"

        if cfg["mode"] == "paired":
            rev_in = raw_fastq / fwd_in.name.replace("_1.fastq", "_2.fastq")
            if not rev_in.exists():
                rev_in = raw_fastq / fwd_in.name.replace("_1.fastq.gz", "_2.fastq.gz")
            if not rev_in.exists():
                print(f"    No reverse read for {sample}, skipping.")
                n_fail += 1
                continue

            rev_out = trimmed_dir / f"{sample}_2.fastq.gz"
            cmd = [
                "cutadapt",
                "-g", PRIMER_FWD,
                "-G", PRIMER_REV,
                "-a", PRIMER_REV_RC,
                "-A", PRIMER_FWD_RC,
                "--discard-untrimmed",
                "--minimum-length", "150",
                "--error-rate", "0.15",
                "-j", str(threads),
                "-o", str(fwd_out),
                "-p", str(rev_out),
                str(fwd_in), str(rev_in),
            ]
        else:
            cmd = [
                "cutadapt",
                "-g", PRIMER_FWD,
                "-a", PRIMER_REV_RC,
                "--discard-untrimmed",
                "--minimum-length", "150",
                "--error-rate", "0.15",
                "-j", str(threads),
                "-o", str(fwd_out),
                str(fwd_in),
            ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    cutadapt failed for {sample}:\n{result.stderr[:300]}")
            n_fail += 1
        else:
            n_ok += 1

    print(f"  {cohort}: {n_ok} samples trimmed, {n_fail} failed.")
    if n_fail > 0:
        print(f"  Warning: {n_fail} samples failed primer trimming. "
              "Check that primer sequences match those used in library prep.")
    return trimmed_dir


def run_dada2(cohort: str, silva_dir: Path, threads: int,
              inspect_only: bool, force: bool) -> bool:
    cfg      = COHORT_CONFIG[cohort]
    trimmed_dir = RAW_DIR / cohort / "fastq_trimmed"
    if inspect_only:
        fastq_dir = RAW_DIR / cohort / "fastq"
    elif cfg.get("primers_pretrimmed", False):
        fastq_dir = RAW_DIR / cohort / "fastq"
    elif trimmed_dir.exists() and any(trimmed_dir.iterdir()):
        fastq_dir = trimmed_dir
    else:
        fastq_dir = RAW_DIR / cohort / "fastq"
        print(f"  Warning: {cohort} has no trimmed reads in fastq_trimmed/. "
              "Running DADA2 on raw reads — primer contamination may affect results. "
              "Run primer trimming first for best results.")
    out_dir   = DADA2_DIR / cohort
    genus_out = out_dir / "genus_table.csv"

    if genus_out.exists() and not inspect_only and not force:
        print(f"  Skipping {cohort}: genus_table.csv already exists. "
              "Use --force to re-run DADA2.")
        return True

    silva_path = silva_dir / "silva_nr99_v138.1_train_set.fa.gz"

    cmd = [
        "Rscript", str(DADA2_SCRIPT),
        "--cohort",       cohort,
        "--fastq_dir",    str(fastq_dir),
        "--out_dir",      str(out_dir),
        "--mode",         cfg["mode"],
        "--trunc_fwd",    str(cfg["trunc_fwd"]),
        "--trunc_rev",    str(cfg["trunc_rev"] or 0),
        "--silva_path",   str(silva_path),
        "--threads",      str(threads),
        "--inspect_only", str(inspect_only).upper(),
    ]

    print(f"\n  Starting DADA2 for {cohort} "
          f"({cfg['mode']}, trunc={cfg['trunc_fwd']}"
          + (f"/{cfg['trunc_rev']}" if cfg["mode"] == "paired" else "") + ")")
    print(f"  Running: {' '.join(cmd)}\n")

    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        print(f"  dada2_pipeline.R failed for {cohort} "
              f"(exit code {result.returncode})")
        return False
    return True


def load_genus_table(cohort: str) -> pd.DataFrame:
    path = DADA2_DIR / cohort / "genus_table.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"genus_table.csv not found for {cohort} at {path}\n"
            "Run Stage A (DADA2) first: python 01_harmonize_taxonomy.py --dada2_only"
        )
    df = pd.read_csv(path, index_col=0)
    df.index.name = "sample_id"
    return df


def apply_synonym_map(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns=GENUS_SYNONYMS)
    df = df.T.groupby(level=0).sum().T
    return df


def drop_uninformative_genera(df: pd.DataFrame) -> pd.DataFrame:
    cols_to_drop = [c for c in df.columns if c in GENERA_TO_DROP
                    or c.startswith("uncultured")
                    or c.startswith("metagenome")]
    if cols_to_drop:
        print(f"    Dropping {len(cols_to_drop)} uninformative genera: "
              f"{cols_to_drop[:5]}{'...' if len(cols_to_drop) > 5 else ''}")
    return df.drop(columns=cols_to_drop, errors="ignore")


def apply_prevalence_filter(
    cohort_tables: dict[str, pd.DataFrame],
    min_prevalence: float = 0.20,
    min_cohorts: int = 1,
) -> set[str]:
    genus_pass_counts: dict[str, int] = {}
    for cohort, df in cohort_tables.items():
        if df.empty:
            continue
        prevalence = (df > 0).mean(axis=0)
        passing = set(prevalence[prevalence >= min_prevalence].index)
        for g in passing:
            genus_pass_counts[g] = genus_pass_counts.get(g, 0) + 1

    retained = {g for g, n in genus_pass_counts.items() if n >= min_cohorts}
    print(f"  Prevalence filter (≥{min_prevalence*100:.0f}% in ≥{min_cohorts} "
          f"cohort): {len(retained)} genera retained")
    return retained


def load_sample_metadata(cohort: str) -> pd.DataFrame:
    meta_path = RAW_DIR / cohort / "metadata.tsv"
    if not meta_path.exists():
        print(f"    Warning: metadata.tsv not found for {cohort}. "
              "Run 00_data_acquisition.py first.")
        return pd.DataFrame()

    meta = pd.read_csv(meta_path, on_bad_lines="skip")

    if "Run" not in meta.columns:
        print(f"    Warning: no 'Run' column in {cohort} metadata. "
              "Columns found: " + str(list(meta.columns[:10])))
        return pd.DataFrame()

    meta = meta.set_index("Run")

    cfg = COHORT_CONFIG[cohort]
    meta["cohort"]  = cohort
    meta["country"] = cfg["country"]
    meta["mode"]    = cfg["mode"]

    return meta


def harmonize_all(
    cohorts: list[str],
    min_prevalence: float = 0.20,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    print("\nStage B: cross-cohort taxonomy harmonization")

    raw_tables: dict[str, pd.DataFrame] = {}
    tables = raw_tables  # alias, habit from notebook work
    for cohort in cohorts:
        print(f"\n  Loading {cohort}...")
        try:
            df = load_genus_table(cohort)
        except FileNotFoundError as e:
            print(f"  {e}")
            continue
        df = apply_synonym_map(df)
        df = drop_uninformative_genera(df)
        n_samples, n_genera = df.shape
        print(f"    {n_samples} samples × {n_genera} genera after synonym map")
        tables[cohort] = df

    if not tables:
        raise RuntimeError("No genus tables loaded. Did Stage A complete?")

    print("\n  Applying prevalence filter...")
    retained_genera = apply_prevalence_filter(tables, min_prevalence)

    aligned: dict[str, pd.DataFrame] = {}
    for cohort, df in tables.items():
        present_retained = [g for g in retained_genera if g in df.columns]
        missing          = retained_genera - set(df.columns)
        aligned_df = df[present_retained].copy()
        for g in missing:
            aligned_df[g] = 0
        aligned_df = aligned_df[sorted(retained_genera)]
        aligned[cohort] = aligned_df
        print(f"  {cohort}: {len(present_retained)} of {len(retained_genera)} "
              f"retained genera present ({len(missing)} filled with 0)")

    print("\n  Genus presence across cohorts (top 20 cohort-specific genera):")
    presence_df = pd.DataFrame({
        c: (df > 0).mean() for c, df in aligned.items()
    }).T
    presence_range = presence_df.max() - presence_df.min()
    top_discordant = presence_range.nlargest(20)
    for genus, diff in top_discordant.items():
        row = presence_df[genus]
        print(f"    {genus:<45} range={diff:.2f} | " +
              " | ".join(f"{c}:{row[c]:.2f}" for c in cohorts if c in row.index))

    print("\n  Loading sample metadata...")
    meta_frames = []
    for cohort in cohorts:
        m = load_sample_metadata(cohort)
        if not m.empty:
            meta_frames.append(m)
    metadata_df = pd.concat(meta_frames) if meta_frames else pd.DataFrame()

    all_dfs = []
    for cohort, df in aligned.items():
        df = df.copy()
        df.insert(0, "cohort", cohort)
        all_dfs.append(df)

    unified = pd.concat(all_dfs, axis=0)
    unified.index.name = "sample_id"

    print(f"\n  Unified matrix: {unified.shape[0]} samples × "
          f"{unified.shape[1] - 1} genera (+1 cohort column)")

    return unified, metadata_df


def main():
    parser = argparse.ArgumentParser(
        description="Run DADA2 + taxonomy harmonization for all 5 cohorts."
    )
    parser.add_argument(
        "--silva_dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "silva",
        help="Directory containing SILVA 138.1 .fa.gz files (default: data/silva/)",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Threads for DADA2 (default: 4)",
    )
    parser.add_argument(
        "--cohort",
        choices=list(COHORT_CONFIG.keys()),
        default=None,
        help="Run a single cohort. Omit to run all.",
    )
    parser.add_argument(
        "--dada2_only", action="store_true",
        help="Run Stage A (DADA2) only; skip Stage B harmonization.",
    )
    parser.add_argument(
        "--harmonize_only", action="store_true",
        help="Run Stage B (harmonization) only; skip Stage A DADA2.",
    )
    parser.add_argument(
        "--inspect_only", action="store_true",
        help="Generate DADA2 quality plots only (no DADA2, no harmonization).",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run DADA2 and primer trimming even if outputs already exist.",
    )
    parser.add_argument(
        "--skip_trim", action="store_true",
        help="Skip primer trimming (Stage A0). Use only if reads are pre-trimmed "
             "or for debugging. Not recommended for production runs.",
    )
    parser.add_argument(
        "--min_prevalence", type=float, default=0.20,
        help="Min fraction of samples a genus must appear in (default: 0.20).",
    )
    args = parser.parse_args()

    cohorts = [args.cohort] if args.cohort else list(COHORT_CONFIG.keys())

    if not args.harmonize_only:
        print("\nStage A: DADA2 processing")

        silva_path = args.silva_dir / "silva_nr99_v138.1_train_set.fa.gz"
        if not args.inspect_only and not silva_path.exists():
            sys.exit(
                f"\nSILVA training set not found at {silva_path}\n"
                "Download with:\n"
                "  mkdir -p data/silva && cd data/silva\n"
                "  wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz\n"
                "  wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz\n"
            )

        if not args.inspect_only and not args.skip_trim:
            print("\nStage A0: primer trimming (cutadapt)")
            for cohort in cohorts:
                trim_primers_cutadapt(cohort, args.threads, args.force)

        DADA2_DIR.mkdir(parents=True, exist_ok=True)
        failed_cohorts = []
        for cohort in cohorts:
            success = run_dada2(
                cohort, args.silva_dir, args.threads,
                args.inspect_only, args.force
            )
            if not success:
                failed_cohorts.append(cohort)

        if failed_cohorts:
            print(f"\nWarning: DADA2 failed for: {failed_cohorts}")
            print("Fix the errors above and re-run. Already-completed cohorts will be skipped.")
            if not args.dada2_only:
                print("Continuing to Stage B with completed cohorts...")

        if args.inspect_only:
            print("\nQuality plots written. Review PDFs in data/processed/dada2/<cohort>/")
            print("Then re-run without --inspect_only to run full DADA2.")
            return

        if args.dada2_only:
            print("\nStage A complete.")
            return

    GENUS_DIR.mkdir(parents=True, exist_ok=True)

    unified, metadata_df = harmonize_all(
        cohorts,
        min_prevalence=args.min_prevalence,
    )

    cohort_col = unified["cohort"]
    genera_only = unified.drop(columns=["cohort"])
    for cohort in cohorts:
        mask = cohort_col == cohort
        sub = genera_only[mask]
        rel = sub.div(sub.sum(axis=1), axis=0).fillna(0)
        rel.to_csv(GENUS_DIR / f"{cohort}_genus_relabund.csv")
        print(f"  Wrote {GENUS_DIR}/{cohort}_genus_relabund.csv "
              f"({rel.shape[0]} samples × {rel.shape[1]} genera)")

    out_path = PROCESSED_DIR / "unified_genus_matrix.csv"
    unified.to_csv(out_path)
    print(f"\n  Unified count matrix: {out_path}")

    if not metadata_df.empty:
        meta_out = PROCESSED_DIR / "sample_metadata.csv"
        metadata_df.to_csv(meta_out)
        print(f"  Sample metadata: {meta_out}")
    else:
        print("\n  Warning: no metadata loaded. Sample metadata will need to be "
              "curated manually from paper supplementary tables for each cohort.\n"
              "  Expected: SRA BioSample attributes should contain diagnosis "
              "(AD/MCI/CN) fields.\n"
              "  Action: inspect data/raw/<cohort>/metadata.tsv for each cohort "
              "and confirm the diagnosis column name.")

    summary_rows = []
    for cohort in cohorts:
        stats_path = DADA2_DIR / cohort / "processing_stats.csv"
        if stats_path.exists():
            stats = pd.read_csv(stats_path)
            row = {
                "cohort":        cohort,
                "country":       COHORT_CONFIG[cohort]["country"],
                "mode":          COHORT_CONFIG[cohort]["mode"],
                "n_samples":     len(stats),
                "reads_input":   stats["input"].sum() if "input" in stats.columns else None,
                "reads_nonchim": stats["nonchim"].sum() if "nonchim" in stats.columns else None,
            }
            if row["reads_input"] and row["reads_nonchim"]:
                row["pct_retained"] = round(
                    100 * row["reads_nonchim"] / row["reads_input"], 1
                )
            summary_rows.append(row)

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_out = PROCESSED_DIR / "pipeline_summary.csv"
        summary_df.to_csv(summary_out, index=False)
        print(f"\n  Pipeline summary: {summary_out}")
        print(summary_df.to_string(index=False))

    print("\nStage B complete.")
    print("Next step: python scripts/02_clr_transform.py")


if __name__ == "__main__":
    main()
