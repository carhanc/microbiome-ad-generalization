"""Prevalence-restricted Jaccard null sensitivity (Reviewer 4).

Reviewer 3's statistic-matching issue (the null must average six pairwise
Jaccards among four random sets, matching the observed statistic's own
construction) is already fixed in scripts/08_shap_taxa_comparison.py. This
script addresses a separate issue Reviewer 4 raises: the uniform null draws
random genera from the full 396-genus pool as if every genus were equally
available in every cohort, but genus prevalence varies substantially by
cohort (Kazakhstan alone passes >=20% prevalence for 356/396 genera; the
three MiSeq cohorts each pass for ~115-131). A genus that is essentially
absent in a given cohort could still be drawn into that cohort's random
comparison set under the uniform null, which is not a realistic null model
for "how much feature-availability-adjusted top-N overlap would arise by
chance in this specific cohort."

This script restricts each cohort's random draws (and its own top-N ranking)
to genera that actually reach >=20% prevalence within that cohort's own
AD/CN-labeled samples -- a cohort-specific eligible universe, not the shared
396-genus pool. The four cohorts' eligible universes differ in both content
and size; verified below that every cohort has >=50 eligible genera, so no N
adaptation is needed for N in {10, 20, 50}.
"""
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

SUPERVISED_COHORTS = ["zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022"]
MIN_PREVALENCE = 0.20
N_LIST = [10, 20, 50]
N_PERM = 100_000
RANDOM_STATE = 42


def jaccard(a, b):
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def eligible_universe_per_cohort():
    unified = pd.read_csv(PROCESSED_DIR / "unified_genus_matrix.csv", index_col=0)
    diag = pd.read_csv(PROCESSED_DIR / "diagnosis_labels.csv")
    labeled = diag[diag["diagnosis"].isin(["AD", "CN"])].set_index("run_id")
    genus_cols = [c for c in unified.columns if c != "cohort"]

    universes = {}
    for c in SUPERVISED_COHORTS:
        ids = labeled[labeled["cohort"] == c].index
        sub = unified.loc[unified.index.intersection(ids), genus_cols]
        prevalence = (sub > 0).mean(axis=0)
        eligible = set(prevalence[prevalence >= MIN_PREVALENCE].index)
        universes[c] = sorted(eligible)
        print(f"  {c}: n={len(sub)}, eligible genera (>={MIN_PREVALENCE*100:.0f}% prevalence) = {len(eligible)}")
    return universes


def restricted_top_n(shap_df, model_name, cohort, eligible_set, top_n):
    sub = shap_df[(shap_df["model"] == model_name) & (shap_df["cohort"] == cohort)]
    sub = sub[sub["taxon"].isin(eligible_set)].sort_values("mean_abs_shap", ascending=False)
    return set(sub["taxon"].head(top_n))


def observed_six_pair_mean(top_sets):
    pairs = [(a, b) for i, a in enumerate(SUPERVISED_COHORTS)
             for b in SUPERVISED_COHORTS[i + 1:]]
    assert len(pairs) == 6
    vals = [jaccard(top_sets[a], top_sets[b]) for a, b in pairs]
    return float(np.mean(vals)), vals


def null_replicate(universes, top_n, rng):
    draws = {c: set(rng.choice(universes[c], size=min(top_n, len(universes[c])), replace=False))
             for c in SUPERVISED_COHORTS}
    pairs = [(a, b) for i, a in enumerate(SUPERVISED_COHORTS)
             for b in SUPERVISED_COHORTS[i + 1:]]
    vals = [jaccard(draws[a], draws[b]) for a, b in pairs]
    return float(np.mean(vals))


def main():
    print("15_jaccard_prevalence_null.py -- Reviewer 4 prevalence-restricted sensitivity")

    print("\nComputing per-cohort eligible genus universes:")
    universes = eligible_universe_per_cohort()
    for c in SUPERVISED_COHORTS:
        for n in N_LIST:
            if len(universes[c]) < n:
                print(f"  WARNING: {c} has only {len(universes[c])} eligible genera, "
                      f"fewer than N={n}. This should not occur for N<=50 per the pre-check; "
                      f"if it does, N will be capped to the available count for that cohort "
                      f"and this is documented in the output row.")

    shap_df = pd.read_csv(TABLES_DIR / "shap_within_cohort_importance.csv")

    rows = []
    for model_name in ["logreg", "lgbm"]:
        for top_n in N_LIST:
            print(f"\n  {model_name}, N={top_n}")
            top_sets = {c: restricted_top_n(shap_df, model_name, c, set(universes[c]), top_n)
                       for c in SUPERVISED_COHORTS}
            for c in SUPERVISED_COHORTS:
                print(f"    {c}: restricted top-{top_n} set size = {len(top_sets[c])} "
                      f"(eligible universe = {len(universes[c])})")

            observed_mean, pair_vals = observed_six_pair_mean(top_sets)

            rng = np.random.default_rng(RANDOM_STATE)
            null_vals = np.array([null_replicate(universes, top_n, rng) for _ in range(N_PERM)])
            null_mean = float(null_vals.mean())
            ci_lo, ci_hi = np.percentile(null_vals, [2.5, 97.5])
            p_value = (1 + np.sum(null_vals >= observed_mean)) / (N_PERM + 1)

            where = "ABOVE" if observed_mean > ci_hi else ("BELOW" if observed_mean < ci_lo else "WITHIN")
            print(f"    observed mean-of-six (prevalence-restricted) = {observed_mean:.4f} | "
                  f"null mean = {null_mean:.4f} [{ci_lo:.4f}, {ci_hi:.4f}] | "
                  f"{where} null 95% CI | p = {p_value:.6f}")

            rows.append({
                "model": model_name, "top_n": top_n,
                "observed_mean_jaccard_prevalence_restricted": round(observed_mean, 4),
                "null_mean": round(null_mean, 4),
                "null_ci_lo": round(float(ci_lo), 4),
                "null_ci_hi": round(float(ci_hi), 4),
                "empirical_p": round(p_value, 6),
                "n_perm": N_PERM,
                "eligible_universe_zhuang2018": len(universes["zhuang2018"]),
                "eligible_universe_ling2020": len(universes["ling2020"]),
                "eligible_universe_shanghai2022": len(universes["shanghai2022"]),
                "eligible_universe_kazakhstan2022": len(universes["kazakhstan2022"]),
            })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(TABLES_DIR / "jaccard_prevalence_restricted_null.csv", index=False)
    print(f"\n  wrote {TABLES_DIR}/jaccard_prevalence_restricted_null.csv")
    print(out_df.to_string(index=False))
    print("\nDone.")


if __name__ == "__main__":
    main()
