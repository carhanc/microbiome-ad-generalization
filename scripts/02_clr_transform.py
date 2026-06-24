import numpy as np
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"


def multiplicative_replacement(counts: np.ndarray, delta: float = 0.65) -> np.ndarray:
    result = counts.astype(float).copy()
    for i in range(result.shape[0]):
        row = result[i]
        n_zeros = (row == 0).sum()
        if n_zeros == 0:
            continue
        total = row.sum()
        if total == 0:
            result[i] = delta / result.shape[1]
            continue
        zero_mask  = row == 0
        nonzero_mask = ~zero_mask
        replacement_total = n_zeros * delta
        if replacement_total >= total:
            delta_adj = total * 0.01 / n_zeros
            replacement_total = n_zeros * delta_adj
            result[i, zero_mask] = delta_adj
        else:
            result[i, zero_mask] = delta
        scale = (total - replacement_total) / row[nonzero_mask].sum()
        result[i, nonzero_mask] = row[nonzero_mask] * scale
    return result


def clr_transform(counts: np.ndarray) -> np.ndarray:
    log_counts = np.log(counts)
    log_geomean = log_counts.mean(axis=1, keepdims=True)
    return log_counts - log_geomean


def main():
    print("02_clr_transform.py — CLR transform")

    in_path = PROCESSED_DIR / "unified_genus_matrix.csv"
    print(f"\nLoading {in_path}")
    matrix = pd.read_csv(in_path, index_col=0)
    print(f"  {matrix.shape[0]} samples × {matrix.shape[1]} columns")

    cohort_col = matrix["cohort"].copy()
    genera = matrix.drop(columns=["cohort"])
    print(f"  {genera.shape[1]} genus columns")
    print(f"  cohorts: {cohort_col.value_counts().to_dict()}")

    n_zeros = (genera == 0).sum().sum()
    total_cells = genera.shape[0] * genera.shape[1]
    pct_zeros = round(100 * n_zeros / total_cells, 1)
    print(f"\n  {n_zeros:,} zeros / {total_cells:,} cells ({pct_zeros}%)")
    zeros_per_sample = (genera == 0).sum(axis=1)
    print(f"  zeros per sample: min={zeros_per_sample.min()}, "
          f"median={zeros_per_sample.median():.0f}, "
          f"max={zeros_per_sample.max()}")

    print(f"\n  multiplicative replacement (δ=0.65) on {n_zeros:,} zeros...")
    counts_arr = genera.to_numpy()
    counts_arr = counts_arr  # already a copy from to_numpy, fine
    counts_replaced = multiplicative_replacement(counts_arr, delta=0.65)

    assert (counts_replaced == 0).sum() == 0, "Zero replacement failed — zeros remain"
    print("  done, no zeros left")

    print("\n  applying CLR...")
    clr_arr = clr_transform(counts_replaced)

    row_sums = clr_arr.sum(axis=1)
    max_row_sum = max(abs(x) for x in row_sums)  # same as np.abs().max(), just clearer to me
    print(f"  row-sum check: max |sum| = {max_row_sum:.2e}")
    assert max_row_sum < 1e-8, f"CLR sanity check failed: max row sum = {max_row_sum}"

    clr_df = pd.DataFrame(clr_arr, index=genera.index, columns=genera.columns)
    clr_df.insert(0, "cohort", cohort_col)

    print("\n  per-cohort CLR summary:")
    summary_rows = []
    for cohort in cohort_col.unique():
        mask = cohort_col == cohort
        sub = clr_arr[mask.values]
        summary_rows.append({
            "cohort":          cohort,
            "n_samples":       sub.shape[0],
            "clr_mean":        round(sub.mean(), 4),
            "clr_std":         round(sub.std(), 4),
            "clr_min":         round(sub.min(), 4),
            "clr_max":         round(sub.max(), 4),
            "pct_positive":    round(100 * (sub > 0).mean(), 1),
        })
    summary_df = pd.DataFrame(summary_rows)
    print(summary_df.to_string(index=False))

    out_matrix  = PROCESSED_DIR / "clr_matrix.csv"
    out_summary = PROCESSED_DIR / "clr_summary.csv"

    clr_df.to_csv(out_matrix)
    summary_df.to_csv(out_summary, index=False)

    print(f"\n  wrote {out_matrix} ({clr_df.shape[0]} × {clr_df.shape[1]})")
    print(f"  wrote {out_summary}")
    print("\nFinished. Next: python scripts/03_within_cohort_baseline.py")


if __name__ == "__main__":
    main()
